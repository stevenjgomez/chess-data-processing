import pandas as pd
from nexusformat.nexus import *
from scipy.spatial.transform import Rotation as R
from scipy.optimize import basinhopping
import hkl
import sys, os
import logging
from os.path import exists
import gc
from _checkccode import checkccode

checkccode()

def get_peaklist(projectdir, valmin, valmax, lower_bound, upper_bound):
    projdir = projectdir
    percofmax = float(valmax)
    percofmin = float(valmin)
    stack_file = sorted([f for f in os.listdir(projectdir) if f.startswith("stack")])[0]
    stack = nxload(projectdir + stack_file)
    nxsetmemory(100000)
    print("Loading stack...")
    Iall = stack.data.counts.nxdata
    print("Stack loaded.")
    peaksmax = np.max(Iall)
    print("Max intensity is ", peaksmax)
    print("Finding peaks...")
    peaks = Iall > percofmax * peaksmax
    peaks.shape
    listofpeaks = np.asarray(np.where(peaks)).T

    while len(listofpeaks) < lower_bound:
        print("Number of peaks found: " + str(len(listofpeaks)))
        percofmax = percofmax - 0.1
        print('Not enough peaks found.')
        gc.collect()
        print('Decreasing to: ' + str(percofmax))

        peaks = np.logical_and(Iall < percofmax * peaksmax, Iall > percofmin * peaksmax)
        peaks.shape
        listofpeaks = np.asarray(np.where(peaks)).T

    if len(listofpeaks) > upper_bound:
        while len(listofpeaks) > upper_bound:
            print("Number of peaks found: " + str(len(listofpeaks)))
            percofmax = percofmax + 0.02
            print('Too many peaks found.')
            print('Increasing to: ' + str(percofmax))

            peaks = np.logical_and(Iall < percofmax * peaksmax, Iall > percofmin * peaksmax)
            peaks.shape
            listofpeaks = np.asarray(np.where(peaks)).T
        if len(listofpeaks) < lower_bound:
            while len(listofpeaks) < lower_bound:
                print("Number of peaks found: " + str(len(listofpeaks)))
                percofmax = percofmax - 0.01
                print('Not enough peaks found.')
                gc.collect()
                print('Decreasing to: ' + str(percofmax))

                peaks = np.logical_and(Iall < percofmax * peaksmax, Iall > percofmin * peaksmax)
                peaks.shape
                listofpeaks = np.asarray(np.where(peaks)).T
    print(len(listofpeaks))

    peaksout = projdir + "peaklist1"
    np.save(peaksout, listofpeaks)


def get_length_plist(projectdir):
    peaklist = np.load(projectdir + "peaklist1.npy")
    print('Number of peaks found: ' + str(len(peaklist)))


def find_euler(projectdir):
    proj = projectdir
    peakfile = projectdir + "peaklist1.npy"
    stack_file = sorted([f for f in os.listdir(projectdir) if f.startswith("stack")])[0]
    w = nxload(projectdir + stack_file)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(proj + "ormfinder.log"),
            logging.StreamHandler(sys.stdout)
        ]
    )

    peaklist = np.load(peakfile)
    peaknum = len(peaklist)
    logging.info("Loaded {n} peak points from file.".format(n=peaknum))
    prelimflag = 0
    if (peaknum > 75):
        shortpeaklist = peaklist[:20]
        shortpeaknum = len(shortpeaklist)
        prelimflag = 1

    if exists(projectdir + "cell.txt"):
        df = pd.read_csv(projectdir + "cell.txt", header=None)
        uca = float(df[0])
        ucb = float(df[1])
        ucc = float(df[2])
        ucal = float(df[3])
        ucbe = float(df[4])
        ucga = float(df[5])
    else:
        print("Enter unit cell paramters (units of angstroms & degrees)\n")
        uca = float(input("a = "))
        ucb = float(input("b = "))
        ucc = float(input("c = "))
        ucal = float(input("alpha = "))
        ucbe = float(input("beta = "))
        ucga = float(input("gamma = "))

    logging.info("UCell: {a} {b} {c} {al} {be} {ga}".format(a=uca, b=ucb, c=ucc, al=ucal, be=ucbe, ga=ucga))
    philu = w.data.phi.nxdata

    wl = w.geo.wl
    logging.info("Wavelength  = {wl}".format(wl=wl))
    eta = w.psic.eta
    mu = w.psic.mu
    chi = w.psic.chi
    dpsi = 0.0 * np.pi / 180.0
    myaz = np.asarray(w.geo.az.data)
    mypol = np.asarray(w.geo.pol.data)
    # mytth=np.sqrt((myaz*myaz) + (mypol*mypol))
    # mypsi=np.asarray(w.geo.psi.data)
    # newpol=mytth*np.cos(mypsi+(np.pi*0.5)+dpsi)
    # newaz=-mytth*np.sin(mypsi+(np.pi*0.5)+dpsi)
    newpol = mypol
    newaz = myaz

    def XYtoPOLAZ(X, Y):
        return newpol[X, Y], newaz[X, Y]

    # speed this up by only doing calcB once, passing B to the minfunc
    def calcB(a, b, c, alpha, beta, gamma):
        a1 = a
        a2 = b
        a3 = c
        alpha1 = alpha * np.pi / 180.0
        alpha2 = beta * np.pi / 180.0
        alpha3 = gamma * np.pi / 180.0
        Vt = (1.0 - (np.cos(alpha1) * np.cos(alpha1)) - (np.cos(alpha2) * np.cos(alpha2)) - (
                np.cos(alpha3) * np.cos(alpha3)))
        Vt = Vt + (2.0 * np.cos(alpha1) * np.cos(alpha2) * np.cos(alpha3))
        V = (Vt ** 0.5) * a1 * a2 * a3
        b1 = 2.0 * np.pi * a2 * a3 * np.sin(alpha1) / V
        b2 = 2.0 * np.pi * a3 * a1 * np.sin(alpha2) / V
        b3 = 2.0 * np.pi * a1 * a2 * np.sin(alpha3) / V
        betanum = ((np.cos(alpha2) * np.cos(alpha3)) - np.cos(alpha1))
        betaden = (np.sin(alpha2) * np.sin(alpha3))
        beta1 = np.arccos(betanum / betaden)
        betanum = ((np.cos(alpha1) * np.cos(alpha3)) - np.cos(alpha2))
        betaden = (np.sin(alpha1) * np.sin(alpha3))
        beta2 = np.arccos(betanum / betaden)
        betanum = ((np.cos(alpha1) * np.cos(alpha2)) - np.cos(alpha3))
        betaden = (np.sin(alpha1) * np.sin(alpha2))
        beta3 = np.arccos(betanum / betaden)
        UB = np.zeros(9)
        UB[0] = b1
        UB[1] = b2 * np.cos(beta3)
        UB[2] = b3 * np.cos(beta2)
        UB[3] = 0.0
        UB[4] = -b2 * np.sin(beta3)
        UB[5] = b3 * np.sin(beta2) * np.cos(alpha1)
        UB[6] = 0.0
        UB[7] = 0.0
        UB[8] = -b3
        return UB

    def calcUB(eu1, eu2, eu3, UB):
        r = R.from_euler('zxz', [eu1, eu2, eu3])
        UBR = np.matmul(r.as_matrix(), UB.reshape(3, 3))
        return UBR

    def chisq_U(peaklist, wl, UBR, eta, mu, chi):
        peaknum = len(peaklist)
        chisq = 0.0
        for i in range(0, peaknum):
            pol, az = XYtoPOLAZ(peaklist[i][0], peaklist[i][1])
            phi = philu[peaklist[i][2]]
            IN = hkl.Calc_HKL(np.asarray([pol]), np.asarray([az]), eta, mu, chi, phi, wl, UBR)
            dvec = IN - np.rint(IN)
            chisq = chisq + np.linalg.norm(dvec)
        return chisq / peaknum

    def minfuncU(eu, UB, peaklist, wl):
        UBR = calcUB(eu[0], eu[1], eu[2], UB)
        chisq = chisq_U(peaklist, wl, UBR.T, eta, mu, chi)
        return chisq

    def minfuncWL(fitwl, eu0, eu1, eu2, a, b, c, alpha, beta, gamma, peaklist):
        UBR = calcUB(eu0, eu1, eu2, a, b, c, alpha, beta, gamma)
        chisq = chisq_U(peaklist, fitwl[0], UBR.T, eta, mu, chi)
        return chisq

    myUB = calcB(uca, ucb, ucc, ucal, ucbe, ucga)
    print("Calculated B", myUB)

    if (prelimflag):
        x0 = [0, 0, 0]
        logging.info("Preliminary fitting of first 20 points, starting config A")
        minimizer_kwargs = {"method": "BFGS", "args": (myUB, shortpeaklist, wl)}
        #    pre1 = basinhopping(minfuncU, x0, stepsize=np.pi/30.0, T=0.25, minimizer_kwargs=minimizer_kwargs, niter=100, seed=20)
        pre1 = basinhopping(minfuncU, x0, T=0.002, minimizer_kwargs=minimizer_kwargs)
        logging.info("Eu angles: {angs},  chisq: {chisq}".format(angs=pre1.x, chisq=pre1.fun))
        x0 = [0.5, 1.0, 0.25]
        logging.info("Preliminary fitting of first 20 points, starting config B")
        minimizer_kwargs = {"method": "BFGS", "args": (myUB, shortpeaklist, wl)}
        #    pre2 = basinhopping(minfuncU, x0, stepsize=np.pi/30.0, T=0.25, minimizer_kwargs=minimizer_kwargs, niter=100, seed=21)
        pre2 = basinhopping(minfuncU, x0, T=0.002, minimizer_kwargs=minimizer_kwargs)
        logging.info("Eu angles: {angs},  chisq: {chisq}".format(angs=pre2.x, chisq=pre2.fun))

    # if (pre1.fun<pre2.fun):
    #	x0=pre1.x
    # else:
    #	x0=pre2.x

    x0 = [0, 0, 0]

    logging.info("Fitting all peaks. Starting Eu angles: {angs}".format(angs=x0))

    minimizer_kwargs = {"method": "BFGS", "args": (myUB, peaklist, wl)}
    # ret1 = basinhopping(minfuncU, x0, stepsize=np.pi/30.0, T=0.2, minimizer_kwargs=minimizer_kwargs, niter=100, seed=22)
    ret1 = basinhopping(minfuncU, x0, T=0.002, minimizer_kwargs=minimizer_kwargs)
    logging.info("Eu angles: {angs},  chisq: {chisq}".format(angs=ret1.x, chisq=ret1.fun))

    x0 = ret1.x
    minimizer_kwargs = {"method": "BFGS", "args": (myUB, peaklist, wl)}
    ret2 = basinhopping(minfuncU, x0, minimizer_kwargs=minimizer_kwargs, niter=300, seed=23)
    logging.info("Eu angles: {angs},  chisq: {chisq}".format(angs=ret2.x, chisq=ret2.fun))

    x0 = ret2.x
    minimizer_kwargs = {"method": "BFGS", "args": (myUB, peaklist, wl)}
    ret3 = basinhopping(minfuncU, x0, stepsize=np.pi / 30.0, T=0.005, minimizer_kwargs=minimizer_kwargs, niter=100,
                        seed=24)
    logging.info("Eu angles: {angs},  chisq: {chisq}".format(angs=ret3.x, chisq=ret3.fun))

    return uca, ucb, ucc, ucal, ucbe, ucga, ret3.x


def get_ubr(projectdir, uca, ucb, ucc, ucal, ucbe, ucga, angs):
    stack_file = sorted([f for f in os.listdir(projectdir) if f.startswith("stack")])[0]
    w = nxload(projectdir + stack_file)
    auto_eu = angs
    wl = w.geo.wl

    def calcUB(eu1, eu2, eu3, a, b, c, alpha, beta, gamma):
        a1 = a
        a2 = b
        a3 = c
        alpha1 = alpha * np.pi / 180.0
        alpha2 = beta * np.pi / 180.0
        alpha3 = gamma * np.pi / 180.0
        Vt = (1.0 - (np.cos(alpha1) * np.cos(alpha1)) - (np.cos(alpha2) * np.cos(alpha2)) - (
                np.cos(alpha3) * np.cos(alpha3)))
        Vt = Vt + (2.0 * np.cos(alpha1) * np.cos(alpha2) * np.cos(alpha3))
        V = (Vt ** 0.5) * a1 * a2 * a3
        b1 = 2.0 * np.pi * a2 * a3 * np.sin(alpha1) / V
        b2 = 2.0 * np.pi * a3 * a1 * np.sin(alpha2) / V
        b3 = 2.0 * np.pi * a1 * a2 * np.sin(alpha3) / V
        betanum = ((np.cos(alpha2) * np.cos(alpha3)) - np.cos(alpha1))
        betaden = (np.sin(alpha2) * np.sin(alpha3))
        beta1 = np.arccos(betanum / betaden)
        betanum = ((np.cos(alpha1) * np.cos(alpha3)) - np.cos(alpha2))
        betaden = (np.sin(alpha1) * np.sin(alpha3))
        beta2 = np.arccos(betanum / betaden)
        betanum = ((np.cos(alpha1) * np.cos(alpha2)) - np.cos(alpha3))
        betaden = (np.sin(alpha1) * np.sin(alpha2))
        beta3 = np.arccos(betanum / betaden)
        UB = np.zeros(9)
        UB[0] = b1
        UB[1] = b2 * np.cos(beta3)
        UB[2] = b3 * np.cos(beta2)
        UB[3] = 0.0
        UB[4] = -b2 * np.sin(beta3)
        UB[5] = b3 * np.sin(beta2) * np.cos(alpha1)
        UB[6] = 0.0
        UB[7] = 0.0
        UB[8] = -b3
        r = R.from_euler('zxz', [eu1, eu2, eu3])
        UBR = np.matmul(r.as_matrix(), UB.reshape(3, 3))
        return UBR

    UBRfinal = calcUB(auto_eu[0], auto_eu[1], auto_eu[2], uca, ucb, ucc, ucal, ucbe, ucga)

    return UBRfinal, w, wl


def check_peaklist(projectdir, UBRfinal, w, wl):
    peakfile = projectdir + "peaklist1.npy"
    philu = w.data.phi.nxdata
    # print(wl)
    eta = w.psic.eta
    mu = w.psic.mu
    chi = w.psic.chi
    # dpsi=0.09*np.pi/180.0
    dpsi = 0.0 * np.pi / 180.0
    myaz = np.asarray(w.geo.az.data)
    mypol = np.asarray(w.geo.pol.data)
    mytth = np.sqrt((myaz * myaz) + (mypol * mypol))
    mypsi = np.asarray(w.geo.psi.data)
    newpol = mytth * np.cos(mypsi + (np.pi * 0.5) + dpsi)
    newaz = -mytth * np.sin(mypsi + (np.pi * 0.5) + dpsi)
    newpol = mypol
    newaz = myaz

    def XYtoPOLAZ(X, Y):
        return newpol[X, Y], newaz[X, Y]

    peaklist = np.load(peakfile)
    peaknum = len(peaklist)

    def printconvpeaks(peaklist, wl, UBR, eta, mu, chi):
        peaknum = len(peaklist)
        chisq = 0.0
        for i in range(0, peaknum):
            pol, az = XYtoPOLAZ(peaklist[i][0], peaklist[i][1])
            #        phi=peaklist[i][2]
            phi = philu[peaklist[i][2]]
            IN = hkl.Calc_HKL(np.asarray([pol]), np.asarray([az]), eta, mu, chi, phi, wl, UBR)
            print(IN)

    printconvpeaks(peaklist, wl, UBRfinal.T, eta, mu, chi)


def save_orm(projectdir, UBRfinal, uca, ucb, ucc, ucal, ucbe, ucga):
    ormout = projectdir + "ormatrix_auto.nxs"
    dpsi = 0.0

    print(UBRfinal.ravel().tolist())
    np.save(ormout, UBRfinal)

    ormnex = NXroot()
    ormnex.unitcell = NXentry()
    ormnex.unitcell.a = NXfield(uca, name='a')
    ormnex.unitcell.b = NXfield(ucb, name='b')
    ormnex.unitcell.c = NXfield(ucc, name='c')
    ormnex.unitcell.alpha = NXfield(ucal, name='alpha')
    ormnex.unitcell.beta = NXfield(ucbe, name='beta')
    ormnex.unitcell.gamma = NXfield(ucga, name='gamma')

    ormnex.ormatrix = NXentry()
    ormnex.ormatrix.U = NXfield(UBRfinal, name='Orientation_Matrix')
    ormnex.dspi = NXentry()
    ormnex.dspi.dpsi = NXfield(dpsi, name='detector psi offset')
    ormnex.save(ormout)


def auto_ormfinder(ormdir):
    # Check if peaklist exists
    if exists(ormdir + "peaklist1.npy"):
        print('Peaklist found.')
        ans = input("Generate new peaklist? (y/n): ")
    else:
        ans = 'y'

    # Find peaks if needed/desired
    if ans == 'y':
        find_peaks = True
        while find_peaks:
            default = input("Use default peak finding values? (y/n): ")
            if default == 'y':
                lower_bound = 50
                upper_bound = 150
                valmin = 0.9
                valmax = 1.0
            else:
                lower_bound = int(input("Set minimum number of peaks: "))
                upper_bound = int(input("Set maximum number of peaks: "))
                valmin = float(input("Enter lower bound for peak intensity (0.0 <-> 1.0): "))
                valmax = float(input("Enter upper bound for peak intensity (0.0 <-> 1.0): "))

            get_peaklist(ormdir, valmin, valmax, lower_bound, upper_bound)

            get_length_plist(ormdir)

            ans = input("Proceed (y) or retry (n): ")
            if ans == 'n':
                pass
            if ans == 'y':
                find_peaks = False

    uca, ucb, ucc, ucal, ucbe, ucga, angs = find_euler(ormdir)
    UBRfinal, w, wl = get_ubr(ormdir, uca, ucb, ucc, ucal, ucbe, ucga, angs)
    check_peaklist(ormdir, UBRfinal, w, wl)
    save_orm(ormdir, UBRfinal, uca, ucb, ucc, ucal, ucbe, ucga)
