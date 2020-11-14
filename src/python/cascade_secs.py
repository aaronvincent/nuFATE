""" Funcitons to compute the right hand side matrix of neutrino secondaries evolution and diagonalized it.
"""

import tables
import numpy as np
import scipy as sp
from numpy import linalg as LA
import gr_xs


def get_RHS_matrices(energy_nodes, sigma_array, sig3_array, dxs_array, sec_array,
                     regen_array):
    """ Returns the right hand side (RHS) matrices including secondaries.

    Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        sigma_fname: one dimensional numpy array with total cross sections in cm**2.
        sig3fname: tau neutrino cross section array.
        dxs_fname: differential cross section array.
        dxs_fname: two dimensional numpy array with the differential cross section cm**2 GeV**-1.
        secname: tau regeneration secondaries into nue or numu tables.
        regenname: tau regeneration secondaries into nutau tables.

    Returns:
        RHSMatrix: matrix of size n_nodes*n_nodes containing the E**2 weighted differential
                   cross sections in units of cm**2 GeV.
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
            # Comparing with NuFate paper: multiply by E_j (= E_in) to account
            # for log scale, then by E_i^2/E_j^2 to account for variable change
            # phi -> E^2*phi
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


def get_eigs(flavor, gamma, h5_filename, ato=None):
    """ Returns the eigenvalues for a given flavor, spectral index, and energy range.

    Args:.
        flavor: specifidies the neutrino flavor of interest. Negative numbers for
                antineutrinos, positive numbers for neutrinos.
                1: electron neutrino,
                2: muon neutrino.
                The specify flavor cannot be tau, i.e. 3 or -3.
        gamma:  If gamma is a string, this is the path and file name of the input spectrum (e.g. an atmospheric flux)
                If gamma is a number, it is used as the spectral index of the initial flux, E**-gamma.
        h5_filename: complete path and filename of the h5 object that contains the cross sections.
        ato: Atom on which to include electron velocity for GR doppler broadening. Can be one of 'H O Mg Si Ca Fe'
             Defaults to assuming at rest electrons.

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

    # Note that the solution is scaled by E**2; if you want to modify the incoming spectrum a lot,
    # you may need to change this here, as well as in the definition of RHS.
    RHSMatrix = get_RHS_matrices(energy_nodes, sigma_array, sig3_array,
                                 dxs_array, sec_array, regen_array)

    if flavor == -1: #add glashow pieces
        glashow_piece = (-np.diag(get_glashow_total(energy_nodes, ato))+ get_glashow_partial(energy_nodes, ato))/2.
        z = np.zeros((NumNodes, NumNodes))
        bigG = np.vstack((np.hstack((glashow_piece, z)), np.hstack((z, z))))
        RHSMatrix = RHSMatrix+bigG

    # Select initial condition: if gamma is string, load file, otherwise use power law E^-gamma
    if type(gamma) == str:
        phi_in = np.loadtxt(gamma)
        if phi_in.size != energy_nodes.size:
            raise Exception("Input spectrum must have the same size as the energy vector (default 200x1).")
        phi_0 = energy_nodes**2*phi_in
    else:
        phi_0 = energy_nodes**(2 - gamma)


    phi_0 = np.hstack((phi_0, phi_0))
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

def get_glashow_total(energy_nodes, ato=None):
    """ Returns the total nubar e --> W cross section

        Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        ato: Atom on which to include electron velocity for GR doppler broadening. Can be one of 'H O Mg Si Ca Fe'
             Defaults to assuming at rest electrons.


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
    if ato is None:
        return sig
    return sig*gr_xs.sigma_edopp(Enu, ato)/gr_xs.sigma_erest(Enu)



def get_glashow_partial(energy_nodes, ato=None):
    """ Returns the partial d sigma/dE nubar e --> nubar e cross section

        Args:
        energy_nodes: one dimensional numpy array containing the energy nodes in GeV.
        ato: Atom on which to include electron velocity for GR doppler broadening. Can be one of 'H O Mg Si Ca Fe'
             Defaults to assuming at rest electrons.


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
    if ato is None:
        return dsig
    return dsig*gr_xs.sigma_edopp(energy_nodes, ato)/gr_xs.sigma_erest(energy_nodes)
