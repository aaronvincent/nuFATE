""" Funcitons to compute the right hand side matrix of neutrino evolution and diagonalized it.
"""

import numpy as np
import scipy as sp
from numpy import linalg as LA


def get_RHS_matrices(energy_nodes, sigma_fname, dxs_fname):
    """ Returns the right hand side (RHS) matrices.

    Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        sigma_fname: cross section table filename.
        dxs_fname: differential cross section filename.

    Returns:
        RHSMatrix: matrix of size n_nodes*n_nodes containing the E^2 weighted differential
                   cross sections in units of cm^2 GeV.
        sigma_array: one dimensional numpy array containing the total cross section
                     per energy node in cm^2.
    """
    NumNodes = energy_nodes.shape[0]
    sigma_array = np.loadtxt(sigma_fname)
    dsigmady = np.loadtxt(dxs_fname)
    DeltaE = np.diff(np.log(energy_nodes))
    RHSMatrix = np.zeros((NumNodes, NumNodes))
    # fill in diagonal terms

    for i in range(NumNodes):  #E_out
        for j in range(i + 1, NumNodes):  #E_in
            RHSMatrix[i][j] = DeltaE[j - 1] * dsigmady[j][i] * energy_nodes[
                j]**-1 * energy_nodes[i]**2
    return RHSMatrix, sigma_array


def get_eigs(flavor, gamma=2., logemin=3, logemax=10, NumNodes=200):
    """ Returns the eigenvalues for a given flavor, spectral index, and energy range.

    Args:.
        flavor: specifidies the neutrino flavor of interest. Negative numbers for
                antineutrinos, positive numbers for neutrinos.
                1: electron neutrino,
                2: muon neutrino,
                and 3: tau neutrino.
        gamma: spectral index of the initial flux, E^-gamma.
        logemin: log10 of the lower enegy of interest. Energy in GeV.
        logemax: log10 of the maximum energy of interest. Energy in GeV.
        NumNodes: number of energy nodes to use

    Returns:
        w: right hand side matrix eigenvalues in unit of cm^2 GeV.
        v: right hand side matrix normalized eigenvectors.
        ci: coordinates of the input spectrum in the eigensystem basis.
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        phi_0: input spectrum.
    """
    #    print flavor
    if flavor == -1:
        sigma_fname = "data/nuebarxs.dat"
    elif flavor == -2:
        sigma_fname = "data/numubarxs.dat"
    elif flavor == -3:
        sigma_fname = "data/nutaubarxs.dat"
    elif flavor == 1:
        sigma_fname = "data/nuexs.dat"
    elif flavor == 2:
        sigma_fname = "data/numuxs.dat"
    elif flavor == 3:
        sigma_fname = "data/nutauxs.dat"

    if flavor > 0:
        dxs_fname = "data/dxsnu.dat"
    else:
        dxs_fname = "data/dxsnubar.dat"

    #Note that the solution is scaled by E^2; if you want to modify the incoming 
    #spectrum a lot, you'll need to change this here, as well as in the definition of RHS.
    energy_nodes = np.logspace(logemin, logemax, NumNodes)
    RHSMatrix, sigma_array = get_RHS_matrices(energy_nodes, sigma_fname,
                                              dxs_fname)

    #tau regenration
    if flavor == -3:
        RHregen, s1 = get_RHS_matrices(energy_nodes, sigma_fname,
                                       "data/tbarfull.dat")
        RHSMatrix = RHSMatrix + RHregen

    elif flavor == 3:
        RHregen, s1 = get_RHS_matrices(energy_nodes, sigma_fname,
                                       "data/tfull.dat")
        RHSMatrix = RHSMatrix + RHregen

    phi_0 = energy_nodes**(2 - gamma)
    w, v = LA.eig(-np.diag(sigma_array) + RHSMatrix)
    ci = LA.solve(v, phi_0)  # alternative to lstsq solution
    #    ci = LA.lstsq(v,phi_0)[0]

    return w, v, ci, energy_nodes, phi_0
