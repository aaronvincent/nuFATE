""" Funcitons to compute the right hand side matrix of neutrino secondaries evolution and diagonalized it.
"""

import tables
import numpy as np
import scipy as sp
from numpy import linalg as LA


def get_RHS_matrices(energy_nodes, sigma_array, sig3_array, dxs_array, sec_array,
                     regen_array):
    """ Returns the right hand side (RHS) matrices including secondaries.

    Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        sigma_fname: one dimensional numpy array with total cross sections in cm^2.
        sig3fname: tau neutrino cross section array.
        dxs_fname: differential cross section array.
        dxs_fname: two dimensional numpy array with the differential cross section cm^2 GeV^-1.
        secname: tau regeneration secondaries into nue or numu tables.
        regenname: tau regeneration secondaries into nutau tables.

    Returns:
        RHSMatrix: matrix of size n_nodes*n_nodes containing the E^2 weighted differential
                   cross sections in units of cm^2 GeV.
    """

    NumNodes = len(energy_nodes)
    sigma_array1 = sigma_array
    sigma_array2 = sig3_array
    dsigmady = dxs_array
    emuregen = sec_array
    tauregen = regen_array
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


def get_eigs(flavor, gamma, h5_filename):
    """ Returns the eigenvalues for a given flavor, spectral index, and energy range.

    Args:.
        flavor: specifidies the neutrino flavor of interest. Negative numbers for
                antineutrinos, positive numbers for neutrinos.
                1: electron neutrino,
                2: muon neutrino.
                The specify flavor cannot be tau, i.e. 3 or -3.
        gamma: spectral index of the initial flux, E^-gamma.
        h5_filename: complete path and filename of the h5 object that contains the cross sections.

    Returns:
        w: right hand side matrix eigenvalues in unit of cm^2.
        v: right hand side matrix normalized eigenvectors.
        ci: coordinates of the input spectrum in the eigensystem basis.
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        phi_0: input spectrum.
    """
    xsh5 = tables.open_file(h5_filename,"r")

    if flavor == -1:
        sigma_array = xsh5.root.total_cross_sections.nuebarxs[:]
    elif flavor == -2:
        sigma_array = xsh5.root.total_cross_sections.numubarxs[:]
    elif flavor == 1:
        sigma_array = xsh5.root.total_cross_sections.nuexs[:]
    elif flavor == 2:
        sigma_array = xsh5.root.total_cross_sections.numuxs[:]
    else:
        raise ValueError("OH NO! You need to specify a flavor that is not tau (= 3 or -3). About to crash because of this.")

    if flavor > 0:
        dxs_array = xsh5.root.differential_cross_sections.dxsnu[:]
        sig3_array = xsh5.root.total_cross_sections.nutauxs[:]
        sec_array = xsh5.root.tau_decay_spectrum.secfull[:]
        regen_array = xsh5.root.tau_decay_spectrum.tfull[:]
    else:
        dxs_array = xsh5.root.differential_cross_sections.dxsnubar[:]
        sig3_array = xsh5.root.total_cross_sections.nutaubarxs[:]
        sec_array = xsh5.root.tau_decay_spectrum.secbarfull[:]
        regen_array = xsh5.root.tau_decay_spectrum.tbarfull[:]

    logemax = np.log10(xsh5.root.total_cross_sections._v_attrs.max_energy)
    logemin = np.log10(xsh5.root.total_cross_sections._v_attrs.min_energy)
    NumNodes = xsh5.root.total_cross_sections._v_attrs.number_energy_nodes
    energy_nodes = np.logspace(logemin, logemax, NumNodes)

    # Note that the solution is scaled by E^2; if you want to modify the incoming spectrum a lot,
    # you may need to change this here, as well as in the definition of RHS.
    RHSMatrix = get_RHS_matrices(energy_nodes, sigma_array, sig3_array,
                                 dxs_array, sec_array, regen_array)

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

    xsh5.close()

    return w, v, ci, energy_nodes, phi_0
