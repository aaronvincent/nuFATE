""" Funcitons to compute the right hand side matrix of neutrino evolution and diagonalized it.
"""

import tables
import numpy as np
import scipy as sp
from numpy import linalg as LA

def get_RHS_matrices(energy_nodes, sigma_array, dxs_array, ReverseTime):
    """ Returns the right hand side (RHS) matrices.

    Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        sigma_fname: one dimensional numpy array with total cross sections in cm**2.
        dxs_fname: two dimensional numpy array with the differential cross section cm**2 GeV**-1.

    Returns:
        RHSMatrix: matrix of size n_nodes*n_nodes containing the E**2 weighted differential
                   cross sections in units of cm**2 GeV.
        sigma_array: one dimensional numpy array containing the total cross section
                     per energy node in cm**2.
    """
    NumNodes = len(energy_nodes)
    DeltaE = np.diff(np.log(energy_nodes))
    RHSMatrix = np.zeros((NumNodes, NumNodes))
    # fill in diagonal terms
    if(ReverseTime):
        for i in range(NumNodes):  #E_out
            for j in range(i+1, NumNodes):  #E_in
                RHSMatrix[j][i] = DeltaE[j-1] * dxs_array[j][i] * energy_nodes[
                    j]**-1 * energy_nodes[i]**2
        return RHSMatrix, sigma_array
    else:
        for i in range(NumNodes):  #E_out
            for j in range(i + 1, NumNodes):  #E_in
                RHSMatrix[i][j] = DeltaE[j - 1] * dxs_array[j][i] * energy_nodes[
                    j]**-1 * energy_nodes[i]**2
        return RHSMatrix, sigma_array


def get_eigs(flavor, gamma, h5_filename, ReverseTime, Efinal):
    """ Returns the eigenvalues for a given flavor, spectral index, and energy range.

    Args:.
        flavor: specifidies the neutrino flavor of interest. Negative numbers for
                antineutrinos, positive numbers for neutrinos.
                1: electron neutrino,
                2: muon neutrino,
                and 3: tau neutrino.
        gamma:  If gamma is a string, this is the path and file name of the input spectrum (e.g. an atmospheric flux)
                If gamma is a number, it is used as the spectral index of the initial flux, E**-gamma.
        h5_filename: complete path and filename of the h5 object that contains the cross sections.

    Returns:
        w: right hand side matrix eigenvalues in unit of cm**2.
        v: right hand side matrix normalized eigenvectors.
        ci: coordinates of the input spectrum in the eigensystem basis.
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        phi_0: E^2 * input spectrum.
    """

    xsh5 = tables.open_file(h5_filename,"r")



    if flavor == -1:
        sigma_array = xsh5.root.total_cross_sections.nuebarxs[:]
    elif flavor == -2:
        sigma_array = xsh5.root.total_cross_sections.numubarxs[:]
    elif flavor == -3:
        sigma_array = xsh5.root.total_cross_sections.nutaubarxs[:]
    elif flavor == 1:
        sigma_array = xsh5.root.total_cross_sections.nuexs[:]
    elif flavor == 2:
        sigma_array = xsh5.root.total_cross_sections.numuxs[:]
    elif flavor == 3:
        sigma_array = xsh5.root.total_cross_sections.nutauxs[:]

    if flavor > 0:
        dxs_array = xsh5.root.differential_cross_sections.dxsnu[:]
    else:
        dxs_array = xsh5.root.differential_cross_sections.dxsnubar[:]

    logemax = np.log10(xsh5.root.total_cross_sections._v_attrs.max_energy)
    logemin = np.log10(xsh5.root.total_cross_sections._v_attrs.min_energy)
    NumNodes = xsh5.root.total_cross_sections._v_attrs.number_energy_nodes
    energy_nodes = np.logspace(logemin, logemax, NumNodes)

    #Note that the solution is scaled by E**2; if you want to modify the incoming
    #spectrum a lot, you'll need to change this here, as well as in the definition of RHS.
    RHSMatrix, sigma_array = get_RHS_matrices(energy_nodes, sigma_array,
                                              dxs_array, ReverseTime)

    #tau regenration
    if flavor == -3:
        RHregen, s1 = get_RHS_matrices(energy_nodes, sigma_array,
                                       xsh5.root.tau_decay_spectrum.tbarfull[:], ReverseTime)
        RHSMatrix = RHSMatrix + RHregen

    elif flavor == 3:
        RHregen, s1 = get_RHS_matrices(energy_nodes, sigma_array,
                                       xsh5.root.tau_decay_spectrum.tfull[:], ReverseTime)
        RHSMatrix = RHSMatrix + RHregen
    elif flavor == -1:
        sigma_array = sigma_array + get_glashow_total(energy_nodes)/2
        RHSMatrix = RHSMatrix + get_glashow_partial(energy_nodes)/2

    # Select initial condition: if gamma is string, load file, otherwise use power law E^-gamma
    if type(gamma) == str:
        phi_in = np.loadtxt(gamma)
        if phi_in.size != energy_nodes.size:
            raise Exception('Input spectrum must have the same size as the energy vector (default 200x1).')
        phi_0 = energy_nodes**2*phi_in
    elif ReverseTime:
        phi_0 = time_reversed_phi0(Efinal, energy_nodes)
    else:
        phi_0 = energy_nodes**(2 - gamma)

    w, v = LA.eig(-np.diag(sigma_array) + RHSMatrix)
    ci = LA.solve(v, phi_0)  # alternative to lstsq solution
    #    ci = LA.lstsq(v,phi_0)[0]

    xsh5.close()

    return w, v, ci, energy_nodes, phi_0

def time_reversed_phi0(E, energy_nodes):
    phi_0 = np.zeros(len(energy_nodes))
    for x in range(len(energy_nodes)-1):
        if(energy_nodes[x] <= E <= energy_nodes[x+1]):
            phi_0[x] = 1.
    return phi_0


def get_glashow_total(energy_nodes):
    """ Returns the total nubar e --> W cross section
        
        Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        
        
        Returns:
        total cross section (cm^2)
        """
    Enu = energy_nodes
    GF = 1.16e-5
    hbarc=1.97e-14
    GW = 2.085
    MW = 80.385e0
    mmu=.106e0
    me=511.e-6
    selectron= 2.e0*me*Enu
    sig = 1.e0/3.e0*GF**2*selectron/np.pi*(1.-(mmu**2-me**2)/selectron)**2/((1.-selectron/MW**2)**2+GW**2/MW**2)*.676/.1057 * hbarc**2
    return sig



def get_glashow_partial(energy_nodes):
    """ Returns the partial d sigma/dE nubar e --> nubar e cross section
        
        Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        
        
        Returns:
        differential cross section (cm^2/GeV)
        """
    [Enuin, Enu] = np.meshgrid(energy_nodes,energy_nodes);
    y = 1-Enu/Enuin
    GF = 1.16e-5
    hbarc=1.97e-14
    GW = 2.085
    MW = 80.385
    MZ = 91.18
    me=511.e-6
    s2t = 0.23
    gL =  s2t-0.5
    gR = s2t
    selectron= 2.*me*Enuin
    den = (1-selectron/MW**2)**2 + GW**2/MW**2
    t1 = gR**2/(1.+y*selectron/MZ**2)**2
    t2 = gL/(1.+y*selectron/MZ**2) + (1-selectron/MW**2)/den
    t3 = GW/MW/den
    heaviside = np.piecewise(y, [y < 0, y >= 0], [0., 1.])
    dsig = GF**2*selectron/np.pi*(t1 + (t2**2+t3**2)*(1-y)**2)*hbarc**2*heaviside
    dsig = dsig/Enuin #dy --> dE
    return dsig
