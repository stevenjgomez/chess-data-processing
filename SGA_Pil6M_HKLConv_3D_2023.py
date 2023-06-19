#~# This code takes a stack of detector images from a theta rock, and a working ORmatrix from spec file, and generates a reciprocal space volume #~#
#~# Original version written by J.P. Castellan, Argonne, 2011.  Subsequently modified by Castellan and Ruff, for use at CHESS with Pilatus300K   #~#
#~# THIS is the CHESS kludge for 6M data, to test against the mighty ANL discovery engine workflow
## Modified to accept a nexus datastack, calibrated/ corrected via pyFAI, 2020/2021

from PIL import Image
import numpy as np
import ctypes as ct
import os
import sys
import signal
import time
from nexusformat.nexus import *
import hkl
from pylab import *
from scipy import *
import gc

#some physical constants

PC=4.13566733e-15 	# planks constant in si
c=2.99792458e8 		# speed of light

projd="/nfs/chess/id4baux/2023-2/pollak-3572-b/"
#samdi="/conb3s6_31/sample1/"
#samdi="/conb3s6_33/sample2/"
samdi="NdSb2_lowT/xtal1/"
#samdi="/Li6p75LaZrNbO/Li6p75_LDFZ/"
#samdi="/y2mo2o7/sample5/"

nxsetmemory(100000)
stackhklpath=projd+samdi+sys.argv[1]+"/"

stack1=nxload(stackhklpath+"stack1.nxs")
#stack2=nxload(stackhklpath+"stack2.nxs")
#stack3=nxload(stackhklpath+"stack3.nxs")


#load the orientation info
ormat=nxload(stackhklpath+"ormatrix_doubledP150.nxs")

#######################################################
####  Input parameters - things to change go here  ####
## Ei
## data_path1
## UB matix
## fout
#######################################################

#incident energy and wavelength

#WL=stack1.geo.wl*1e10 #wavelength  old style
WL=stack1.geo.wl #wavelength
print(WL)

#define the hkl space to histogram and define/zero data structures

H=np.arange(-10,10,0.05)
K=np.arange(-10,10,0.05)
L=np.arange(-10,10,0.05)

#H=np.arange(-10,10,0.02)
#K=np.arange(-10,10,0.02)
#L=np.arange(-18,18,0.02)

#H=np.arange(-21,21,0.02)
#K=np.arange(-20,21,0.02)
#L=np.arange(-20,21,0.02)

#H=np.arange(-6.5,6.5, 0.01)
#K=np.arange(-6.5,6.5, 0.01)
#L=np.arange(-6.5,4.0, 0.02)

#H=np.arange(-9.1,9.1, 0.01)
#K=np.arange(-9.1,9.1, 0.01)
#L=np.arange(-3.1,3.1, 0.02)

#H=np.arange(-7.1,7.1, 0.02)
#K=np.arange(-7.1,7.1, 0.02)
#L=np.arange(-7.1,7.1, 0.02)

#where to write the nexus file (the output)
nxs_path=stackhklpath

#name of output file
nxs_name='3scan_HKLI'
#nxs_name='checkorm_HKLI'

print(nxs_name)

#######################################################
####  End of input parameters                      ####
#######################################################

print('H = ',H)
print('K = ',K)
print('L = ',L)

print(len(H)*len(K)*len(L))


#create the data storage arrays
data=np.zeros((len(H)*len(K)*len(L)),dtype=np.float32)
norm=np.zeros((len(H)*len(K)*len(L)),dtype=np.float32)
errors=np.zeros((len(H)*len(K)*len(L)),dtype=np.float32)

#ormatric from file
U=ormat.ormatrix.U

#U=reshape(U,(3,3))
U=U.T# To match spec convension??
#U=np.matrix(U)

def anglerock(H,K,L,stack):

    ap=time.time()
    dtype = np.dtype(np.uint16)

    #counter for number of files processed
    count=0

    #import the spec read command ... it reads the spec scan to get the angles and monitor information for the scan
    angs=np.linspace(0,0,6)
    #in this case, we are flyscanning phi. Need to change for other angle scans
    tth=0.0
    eta=stack.psic.eta.nxdata
    chi=stack.psic.chi.nxdata
    phi=stack.data.phi.nxdata

    #print "Phi = ",phi
    nu=0.0
    mu=0.0

    dpsi=ormat.dspi.dpsi.nxdata
    myaz=np.asarray(stack.geo.az.data)
    mypol=np.asarray(stack.geo.pol.data)
    mytth=np.sqrt((myaz*myaz) + (mypol*mypol))
    mypsi=np.asarray(stack.geo.psi.data)
    pol2=mytth*np.cos(mypsi+(np.pi*0.5)+dpsi)
    az2=-mytth*np.sin(mypsi+(np.pi*0.5)+dpsi)

    n=range(0,len(phi))

    print("Loading stack...")
    Iall=stack.data.counts.nxdata
    print("Stack loaded...")


    #process the frames
    for i in range(0,len(phi)):
#    for i in range(0,900):
        print(i)
        if (stack.norm.icnorm[i]>0.00):
            framenorm=1.0/(stack.norm.icnorm[i]*stack.norm.solidangle) #normalize solid angle and ionchamber
            count=count+1
            IN=hkl.Calc_HKL(pol2,az2,eta,mu,chi,phi[i],WL,U)
            hkl.HIST(IN,(Iall[:,:,i]*framenorm).ravel(),1.0,H,K,L,data,norm,errors)

    bp=time.time()
    n=len(n)
    print("It took ",bp-ap," seconds to process ",count," frames!!!")

    Iall=0
    return U





#process rotation 1
print("Process rotation 1...")
W2=anglerock(H,K,L,stack1)
del stack1
gc.collect()

#process rotation 2
print("Process rotation 2...")
W2=anglerock(H,K,L,stack2)
del stack2
gc.collect()

#process rotation 3
print("Process rotation 3...")
W2=anglerock(H,K,L,stack3)
del stack3
gc.collect()



fout=nxs_path+nxs_name+'.nxs'

dataout=data.clip(0.0)/norm.clip(0.9)

del norm
dataout=dataout.reshape(len(H),len(K),len(L))
H=H.astype('float32')
K=K.astype('float32')
L=L.astype('float32')
H=NXfield(H,name='H',long_name='H')
K=NXfield(K,name='K',long_name='K')
L=NXfield(L,name='L',long_name='L')
dataout=NXfield(dataout,name='counts',long_name='counts')

G=NXdata(dataout,(H,K,L))

if not os.path.exists(nxs_path):
    os.makedirs(nxs_path)

while os.path.exists(fout):
    fout=fout[:-4]+"more.nxs"

G.save(fout)
