# ~# This code takes a stack of detector images from a theta rock, and a working ORmatrix from spec file, and generates a reciprocal space volume #~#
# ~# Original version written by J.P. Castellan, Argonne, 2011.  Subsequently modified by Castellan and Ruff, for use at CHESS with Pilatus300K   #~#
# ~# THIS is the CHESS kludge for 6M data, to test against the mighty ANL discovery engine workflow
## Modified to accept a nexus datastack, calibrated/ corrected via pyFAI, 2020/2021

import gc
from nexusformat.nexus import *
from pylab import *
import os

from anglerock import anglerock


def Pil6M_HKLConv_3D(H, K, L, workingdir, ormdir, stacks_list):
    print('H = ', H)
    print('K = ', K)
    print('L = ', L)
    print("n_bins = " + str(len(H) * len(K) * len(L)))

    # create the data storage arrays
    data = np.zeros((len(H) * len(K) * len(L)), dtype=np.float32)
    norm = np.zeros((len(H) * len(K) * len(L)), dtype=np.float32)
    errors = np.zeros((len(H) * len(K) * len(L)), dtype=np.float32)

    # some physical constants
    PC = 4.13566733e-15  # planks constant in si
    c = 2.99792458e8  # speed of light

    workingdir = os.path.normpath(workingdir) + "\\"
    ormdir = os.path.normpath(ormdir) + "\\"

    rotations_list = sorted([file for file in os.listdir(workingdir) if file.startswith("stack")])
    print("\nWorking directory: " + workingdir)
    selection_list = [i-1 for i in stacks_list]
    rotations_list = [rotations_list[i] for i in selection_list]
    print("\nSelected stacks:")
    [print("[" + str(i) + "] " + rotation) for i, rotation in enumerate(rotations_list)]
    num_rotations = len(rotations_list)

    # If file already exists, append 001
    file_name = str(num_rotations) + 'rot_hkli.nxs'

    suffix = 1
    while os.path.exists(workingdir + file_name):
        file_name = file_name[:-4] + "_" + str(suffix) + ".nxs"
        suffix += 1

    outpath = workingdir + file_name

    # where to write the nexus file (the output)
    nxs_name = file_name[:-4]

    print('\nOutput file: ' + str(outpath))

    nxsetmemory(100000)

    # load the orientation info
    ormat = nxload(ormdir + "ormatrix_auto.nxs")

    # ormatrix from file
    U = ormat.ormatrix.U

    # U=reshape(U,(3,3))
    U = U.T  # To match spec convention??
    # U=np.matrix(U)

    for i in range(0, len(rotations_list)):
        stack = nxload(workingdir + rotations_list[i])

        # incident energy and wavelength
        WL = stack.geo.wl  # wavelength
        print('\nwavelength = ' + str(WL))

        success = False
        while not success:
            try:
                print("\nProcess rotation " + str(i + 1) + " of " + str(len(rotations_list)) + "...")
                W2 = anglerock(H, K, L, stack, ormat, WL, U, data, norm, errors, workingdir, i + 1)
                success = True
                print("Stack " + str(i + 1) + " processed successfully.")

            except MemoryError:
                print("Unable to allocate enough memory! Try closing other programs.")
            except OSError:
                print("Unable to allocate enough memory! Try closing other programs.")
            except Exception as e:
                print(e)

            if not success:
                print("Trying again...")
        print("Clearing stack from memory...")
        del stack
        gc.collect()
        print("Memory cleared.")

    dataout = data.clip(0.0) / norm.clip(0.9)

    del norm
    dataout = dataout.reshape(len(H), len(K), len(L))
    H = H.astype('float32')
    K = K.astype('float32')
    L = L.astype('float32')
    H = NXfield(H, name='H', long_name='H')
    K = NXfield(K, name='K', long_name='K')
    L = NXfield(L, name='L', long_name='L')
    dataout = NXfield(dataout, name='counts', long_name='counts')

    G = NXdata(dataout, (H, K, L))

    G.save(outpath)
