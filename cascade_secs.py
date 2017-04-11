#Produces solution when secondary neutrinos are produced by tau lepton decay
#Not as elegant-looking, but it works.
#nuFATE Aaron Vincent 2017
import numpy as np
import scipy as sp
from numpy import linalg as LA


def get_RHS_matrices(energy_nodes, sigma_fname, sig3fname, dxs_fname, secname,
                     regenname):
    NumNodes = energy_nodes.shape[0]
    sigma_array1 = np.loadtxt(sigma_fname)
    sigma_array2 = np.loadtxt(sig3fname)
    dsigmady = np.loadtxt(dxs_fname)
    emuregen = np.loadtxt(secname)
    tauregen = np.loadtxt(regenname)
    DeltaE = np.diff(np.log(energy_nodes))
    RHSMatrix1 = np.zeros((NumNodes, NumNodes))
    RHSMatrix2 = np.zeros((NumNodes, NumNodes))
    RHSMatrix3 = np.zeros((NumNodes, NumNodes))
    RHSMatrix4 = np.zeros((NumNodes, NumNodes))
    # fill in diagonal terms

    #matrix 1: nue or numu NC:
    for i in range(NumNodes):  #E_out
        for j in range(i + 1, NumNodes):  #E_in
            RHSMatrix1[i][j] = DeltaE[j - 1] * dsigmady[j][i] * energy_nodes[
                j]**-1 * energy_nodes[i]**2

    RHSMatrix1 = -np.diag(sigma_array1) + RHSMatrix1

    # matrix 2: nue/mu production
    for i in range(NumNodes):  #E_out
        for j in range(i + 1, NumNodes):  #E_in
            RHSMatrix2[i][j] = DeltaE[j - 1] * emuregen[j][i] * energy_nodes[
                j]**-1 * energy_nodes[i]**2

    #matrix 4 (3 is zero) tau regen + tau NC
    for i in range(NumNodes):  #E_out
        for j in range(i + 1, NumNodes):  #E_in
            RHSMatrix4[i][j] = DeltaE[j - 1] * (
                dsigmady[j][i] + tauregen[j][i]
            ) * energy_nodes[j]**-1 * energy_nodes[i]**2
    RHSMatrix4 = -np.diag(sigma_array2) + RHSMatrix4

    RHSMatrix = np.vstack((np.hstack((RHSMatrix1, RHSMatrix2)), np.hstack(
        (RHSMatrix3, RHSMatrix4))))

    return RHSMatrix


def get_eigs(flavor, gamma=2., logemin=3, logemax=10, NumNodes=200):

    if flavor == -1:
        sigma_fname = "data/nuebarxs.dat"
    elif flavor == -2:
        sigma_fname = "data/numubarxs.dat"
    elif flavor == 1:
        sigma_fname = "data/nuexs.dat"
    elif flavor == 2:
        sigma_fname = "data/numuxs.dat"
    else:
        print "OH NO! You need to specify a flavor that is not tau (= 3 or -3). About to crash because of this. "

    if flavor > 0:
        dxs_fname = "data/dxsnu.dat"
        sig3fname = "data/nutauxs.dat"
        secname = "data/secfull.dat"
        regenname = "data/tfull.dat"
    else:
        dxs_fname = "data/dxsnubar.dat"
        sig3fname = "data/nutaubarxs.dat"
        secname = "data/secbarfull.dat"
        regenname = "data/tbarfull.dat"

    #Note that the solution is scaled by E^2; if you want to modify the incoming spectrum a lot, you may need to change this here, as well as in the definition of RHS.
    energy_nodes = np.logspace(logemin, logemax, NumNodes)
    RHSMatrix = get_RHS_matrices(energy_nodes, sigma_fname, sig3fname,
                                 dxs_fname, secname, regenname)

    phi_0 = np.hstack((energy_nodes**(2 - gamma), energy_nodes**(2 - gamma)))
    w, v = LA.eig(RHSMatrix)
    #if v is singular, then nothing is really happening. This shouldn't occur for standard model
    #    if LA.det(v)<1.e-40:
    #        w = -sigma_array
    #        v = np.identity(NumNodes)
    #        ci = 1.
    #    else:
    ci = LA.solve(v, phi_0)  # alternative to lstsq solution
    #    ci = LA.lstsq(v,phi_0)[0]

    return w, v, ci, energy_nodes, phi_0
