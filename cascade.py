""" Funcitons to compute the right hand side matrix of neutrino evolution and diagonalized it.
"""

import tables
import numpy as np
import scipy as sp
from numpy import linalg as LA

def get_RHS_matrices(energy_nodes, sigma_array, dxs_array):
    """ Returns the right hand side (RHS) matrices.

    Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        sigma_fname: one dimensional numpy array with total cross sections in cm^2.
        dxs_fname: two dimensional numpy array with the differential cross section cm^2 GeV^-1.

    Returns:
        RHSMatrix: matrix of size n_nodes*n_nodes containing the E^2 weighted differential
                   cross sections in units of cm^2 GeV.
        sigma_array: one dimensional numpy array containing the total cross section
                     per energy node in cm^2.
    """
    NumNodes = len(energy_nodes)
    DeltaE = np.diff(np.log(energy_nodes))
    RHSMatrix = np.zeros((NumNodes, NumNodes))
    # fill in diagonal terms

    for i in range(NumNodes):  #E_out
        for j in range(i + 1, NumNodes):  #E_in
            RHSMatrix[i][j] = DeltaE[j - 1] * dsigmady[j][i] * energy_nodes[
                j]**-1 * energy_nodes[i]**2
    return RHSMatrix, sigma_array


def get_eigs(flavor, gamma, h5_filename):
    """ Returns the eigenvalues for a given flavor, spectral index, and energy range.

    Args:.
        flavor: specifidies the neutrino flavor of interest. Negative numbers for
                antineutrinos, positive numbers for neutrinos.
                1: electron neutrino,
                2: muon neutrino,
                and 3: tau neutrino.
        gamma: spectral index of the initial flux, E^-gamma.
        h5_filename: complete path and filename of the h5 object that contains the cross sections.

    Returns:
        w: right hand side matrix eigenvalues in unit of cm^2.
        v: right hand side matrix normalized eigenvectors.
        ci: coordinates of the input spectrum in the eigensystem basis.
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        phi_0: input spectrum.
    """

    xsh5 = tables.open(h5_filename,"r")

    if flavor == -1:
        sigma_array = xsh5.root.total_cross_section.nuebarxs[:]
    elif flavor == -2:
        sigma_array = xsh5.root.total_cross_section.numubarxs[:]
    elif flavor == -3:
        sigma_array = xsh5.root.total_cross_section.nutaubarxs[:]
    elif flavor == 1:
        sigma_array = xsh5.root.total_cross_section.nuexs[:]
    elif flavor == 2:
        sigma_array = xsh5.root.total_cross_section.numuxs[:]
    elif flavor == 3:
        sigma_array = xsh5.root.total_cross_section.nutauxs[:]

    if flavor > 0:
        dxs_array = xsh5.root.differential_cross_sections.dxsnu[:]
    else:
        dxs_array = xsh5.root.differential_cross_sections.dxsnubar[:]

    logemax = np.log10(root.total_cross_sections._v_attrs.max_energy)
    logemin = np.log10(root.total_cross_sections._v_attrs.min_energy)
    NumNodes = root.total_cross_sections._v_attrs.number_energy_nodes
    energy_nodes = np.logspace(logemin, logemax, NumNodes)

    #Note that the solution is scaled by E^2; if you want to modify the incoming
    #spectrum a lot, you'll need to change this here, as well as in the definition of RHS.
    RHSMatrix, sigma_array = get_RHS_matrices(energy_nodes, sigma_array,
                                              dxs_array)

    #tau regenration
    if flavor == -3:
        RHregen, s1 = get_RHS_matrices(energy_nodes, sigma_fname,
                                       xsh5.root.tau_decay_spectrum.tbarfull[:])
        RHSMatrix = RHSMatrix + RHregen

    elif flavor == 3:
        RHregen, s1 = get_RHS_matrices(energy_nodes, sigma_fname,
                                       xsh5.root.tau_decay_spectrum.tfull[:])
        RHSMatrix = RHSMatrix + RHregen

    phi_0 = energy_nodes**(2 - gamma)
    w, v = LA.eig(-np.diag(sigma_array) + RHSMatrix)
    ci = LA.solve(v, phi_0)  # alternative to lstsq solution
    #    ci = LA.lstsq(v,phi_0)[0]

    xsh5.close()

    return w, v, ci, energy_nodes, phi_0
