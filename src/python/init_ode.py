""" Function to initialize the cross section matrix
"""

import tables
import h5py
import numpy as np
import scipy as sp
from numpy import linalg as LA

def init(flavor,gamma,h5_filename,prefactor=1e70):
    """ Returns the RHSMatrix (excluding density factor) for a given flavor and cross sector file

    Args:
        flavor: specifidies the neutrino flavor of interest. Negative numbers for
                antineutrinos, positive numbers for neutrinos. Only sign matters.
        gamma:  If gamma is a string, this is the path and file name of the input spectrum (e.g. an atmospheric flux)
                If gamma is a number, it is used as the spectral index of the initial flux, E**-gamma.
        h5_filename: complete path and filename of the h5 object that contains the cross sections.
        prefactor: a factor added to the initial fluxes to aviod negative solution at very large zenith angles, default 1e70

    Returns:
        phi_0: E^2 * input spectrum.
        RHSMatrix: three dimensional numpy array with the differential cross section cm**2 GeV**-1, prepared by the init_ode module (which contains the right-hand-side matices excluding density factor).
        energy_nodes: one dimensional numpy array containing the energy nodes of neutrinos in GeV.
        energy_tau: one dimensional numpy array containing the energy nodes of tau in GeV.
    """
    
    
    xsh5 = tables.open_file(h5_filename,"r")
    logemax = np.log10(xsh5.root.total_cross_sections._v_attrs.max_energy)
    logemin = np.log10(xsh5.root.total_cross_sections._v_attrs.min_energy)
    NumNodes = xsh5.root.total_cross_sections._v_attrs.number_energy_nodes
    energy_nodes = np.logspace(logemin, logemax, NumNodes)
    NumTau = 200
    dlg = (logemax - logemin)/(NumTau - 1);
    Ei=np.logspace(logemin - dlg/2,logemax + dlg/2,NumTau + 1)
    energy_tau = np.logspace(logemin, logemax, NumTau)    

    # Select initial condition: if gamma is string, load file, otherwise use power law E^-gamma
    if type(gamma) == str:
        phi_in = np.loadtxt(gamma)
        if phi_in.size != 3*energy_nodes.size+energy_tau.size:
            raise Exception('Input spectrum must have the same size as the energy vector (default 800x1).')
        phi_0 = np.concatenate((np.concatenate((energy_nodes,energy_nodes,energy_nodes), axis=0)**2*phi_in[0:3*NumNodes],energy_tau**2*phi_in[3*NumNodes:]),axis=0)
    else:
        phi_nu = energy_nodes**(2 - gamma)
        phi_0 = np.concatenate((phi_nu, phi_nu, phi_nu, np.zeros(NumTau)),axis=0)
    phi_0 = phi_0*prefactor

    RHSMatrix = np.zeros((9, NumNodes, NumNodes))
    if flavor > 0:
        dxs_array = xsh5.root.differential_cross_sections.dxsnu[:]
    elif flavor < 0:
        dxs_array = xsh5.root.differential_cross_sections.dxsnubar[:]
    else:
        raise Exception('Input flavor must be nonzero, positive for neutrinos and negative for antineutrinos.')
    #Note that the solution is scaled by E**2; if you want to modify the incoming
    #spectrum a lot, you'll need to change this here, as well as in the definition of RHS, or modify the prefactor
    RHSMatrix_NC = get_matrices(energy_nodes, energy_nodes, dxs_array)

    for i in range(1,4):
        f = np.sign(flavor)*i
        if f == -1:
            sigma_array = xsh5.root.total_cross_sections.nuebarxs[:]
        elif f == -2:
            sigma_array = xsh5.root.total_cross_sections.numubarxs[:]
        elif f == -3:
            sigma_array = xsh5.root.total_cross_sections.nutaubarxs[:]
        elif f == 1:
            sigma_array = xsh5.root.total_cross_sections.nuexs[:]
        elif f == 2:
            sigma_array = xsh5.root.total_cross_sections.numuxs[:]
        elif f == 3:
            sigma_array = xsh5.root.total_cross_sections.nutauxs[:]
        RHSMatrix[i-1,:,:] = -np.diag(sigma_array) + RHSMatrix_NC

        if f == -1:
            RHSMatrix[i-1,:,:] = RHSMatrix[i-1,:,:] - np.diag(get_glashow_total(energy_nodes)/2) + get_glashow_partial(energy_nodes)/2
    xsh5.close()

    #tau energy loss
    elossfname = "../../resources/eloss.h5"
    xselossh5 = h5py.File(elossfname,"r")
    elements = ['O','Si','Al','Fe','Ca','Na','K','Mg','Ni','S']
    #atomic number
    A = np.array([16,28,27,56,40,23,39,24,58.7,32])
    frac_crust = np.array([0.467,0.277,0.08,0.05,0.03,0.03,0.03,0.02,0,0])
    frac_mantle = np.array([0.448,0.215,0.02,0.06,0.02,0,0,0.228,0,0])
    frac_core = np.array([0,0,0,0.89,0,0,0,0,0.06,0.05])

    #tau energy loss in propogation
    dsigma_crust = np.zeros((NumTau, NumTau))
    dsigma_mantle = np.zeros((NumTau, NumTau))
    dsigma_core = np.zeros((NumTau, NumTau));
    for i in range(len(elements)):
        dsigma_crust = dsigma_crust + xselossh5['epp/'+elements[i]][:,:]*frac_crust[i]/A[i]
        dsigma_crust = dsigma_crust + xselossh5['pn/'+elements[i]][:,:]*frac_crust[i]/A[i]
        dsigma_mantle = dsigma_mantle + xselossh5['epp/'+elements[i]][:,:]*frac_mantle[i]/A[i]
        dsigma_mantle = dsigma_mantle + xselossh5['pn/'+elements[i]][:,:]*frac_mantle[i]/A[i]
        dsigma_core = dsigma_core + xselossh5['epp/'+elements[i]][:,:]*frac_core[i]/A[i];
        dsigma_core = dsigma_core + xselossh5['pn/'+elements[i]][:,:]*frac_core[i]/A[i]

    #remove the diagonal entries which don't actually count
    dsigma_crust = dsigma_crust - np.diag(np.diag(dsigma_crust))
    dsigma_mantle = dsigma_mantle - np.diag(np.diag(dsigma_mantle))
    dsigma_core = dsigma_core - np.diag(np.diag(dsigma_core))
    RHSMatrix_crust = get_matrices_tau(Ei,energy_tau,dsigma_crust)
    RHSMatrix_mantle = get_matrices_tau(Ei,energy_tau,dsigma_mantle)
    RHSMatrix_core = get_matrices_tau(Ei,energy_tau,dsigma_core)
    #total loss cross section
    diffE = np.diff(Ei)
    diffE.shape = (1, diffE.size)
    integrand = diffE[None, :] * dsigma_crust
    integrand = integrand[0,...]
    xsdiag_crust = np.sum(integrand,1)
    RHSMatrix_crust = - np.diag(xsdiag_crust) + RHSMatrix_crust
    integrand = diffE[None, :] * dsigma_mantle
    integrand = integrand[0,...]
    xsdiag_mantle = np.sum(integrand,1)
    RHSMatrix_mantle = - np.diag(xsdiag_mantle) + RHSMatrix_mantle
    integrand = diffE[None, :] * dsigma_core
    integrand = integrand[0,...]
    xsdiag_core = np.sum(integrand,1)
    RHSMatrix_core = - np.diag(xsdiag_core) + RHSMatrix_core
    RHSMatrix[3,:,:] = RHSMatrix_crust
    RHSMatrix[4,:,:] = RHSMatrix_mantle
    RHSMatrix[5,:,:] = RHSMatrix_core

    #tau production from nutau CC interactions
    if flavor > 0:
        dtauCC = np.loadtxt('../../resources/CT14/differential_cross_sections/dxstau.dat')
    else:
        dtauCC = np.loadtxt('../../resources/CT14/differential_cross_sections/dxstaubar.dat')
    RHSMatrix_43 = get_matrices(energy_nodes,energy_tau, dtauCC)
    RHSMatrix[6,:,:] = RHSMatrix_43
    
    #compute tau decay branching ratios
    dndz_tau = np.zeros((NumTau,NumNodes))
    dndz_emu = np.zeros((NumTau,NumNodes))
    for i in range(NumTau):
        Etau = energy_tau[i]
        for j in range(NumNodes):
            Enutau=energy_nodes[j]
            z = Enutau/Etau
            if z < 1:
                if flavor > 0:
                    dndz_tau[i,j] = dntaudE(z,0,-1) + dntaudE(z,1,-1)
                else:
                    dndz_tau[i,j] = dntaudE(z,0,1) + dntaudE(z,1,1)
                dndz_emu[i,j] = dnnottaudE(z)
            dndz_tau[i,j] = dndz_tau[i,j]/Etau**2
            dndz_emu[i,j] = dndz_emu[i,j]/Etau**2
    RHSMatrix_34 = get_matrices(energy_tau,energy_nodes, dndz_tau);
    RHSMatrix_14 = get_matrices(energy_tau,energy_nodes, dndz_emu);
    #RHSMatrix_24 = RHSMatrix_14 #numu from tau decay
    RHSMatrix[7,:,:] = RHSMatrix_34 #NA*c*ttau*rhobar/mtau is to be divided by later
    RHSMatrix[8,:,:] = RHSMatrix_14 #NA*c*ttau*rhobar/mtau is to be divided by later
    
    return phi_0, RHSMatrix, energy_nodes, energy_tau

def get_matrices(energy_in, energy_out, dxs_array):
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
    [NumIN, NumOUT] = dxs_array.shape
    DeltaE = np.diff(np.log(energy_in))
    RHSMatrix = np.zeros((NumIN, NumOUT))

    for i in range(NumOUT):  #E_out
        for j in range(1, NumIN):  #E_in
            RHSMatrix[i][j] = DeltaE[j - 1] * dxs_array[j][i] * energy_in[
                j]**(-1) * energy_out[i]**2
    return RHSMatrix

def  get_matrices_tau(energy_edges, energy_nodes, dxs_array):
    NumNodes = len(energy_nodes)
    DeltaE = np.diff(np.log(energy_edges))
    RHSMatrix = np.zeros((NumNodes, NumNodes))

    for i in range(NumNodes):  #E_out
        for j in range(i + 1, NumNodes):  #E_in
            RHSMatrix[i][j] = DeltaE[j] * dxs_array[j][i] * energy_nodes[
                j]**(-1) * energy_nodes[i]**2
    return RHSMatrix

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
    [Enuin, Enu] = np.meshgrid(energy_nodes,energy_nodes)
    y = 1-Enu/Enuin
    GF = 1.16e-5
    hbarc=1.97e-14
    GW = 2.085
    MW = 80.385
    MZ = 91.18
    me = 511.e-6
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


def dntaudE(z,channel,pol):
    #tau- (from nutau) --> pol = -1; nutaubar --> pol = +1
    if channel == 0:  #electron channel, to be multiplied by 2
       g0 = 5./3. - 3.*z**2 + 4./3.*z**3
       g1 = 1./3. - 3.*z**2 + 8./3.*z**3
       dndz =  0.18*(g0 + pol*g1)
    elif channel == 1: #hadron channel
        dndz = 0.
        #1)pions
        r = 0.1395**2/1.776**2
        if 1. > r + z:
           g0 = 1./(1. - r) #*heaviside(1.e0-r-z)
           g1 = - (2.*z - 1. + r)/(1. - r)**2 #*heaviside(1.e0-r-z)
           dndz = dndz + 0.12*(g0 + pol*g1)
           #2) rho
           r = 0.77**2/1.776**2
           if 1. > r + z:
               g0 = 1./(1. - r) #*heaviside(1.e0-r-z)
               g1 = - (2.*z - 1. + r)/(1. - r)*(1. - 2*r)/(1 + 2*r) #*heaviside(1.d0-r-z)
               dndz = dndz  + 0.26*(g0 + pol*g1);
               #3) a1
               r = 1.26**2/1.776**2
               if 1. > r + z:
                   g0 = 1./(1. - r) #*heaviside(1.e0-r-z)
                   g1 = - (2.*z - 1. + r)/(1. - r)*(1. - 2*r)/(1 + 2*r) #*heaviside(1.d0-r-z)
                   dndz = dndz + 0.13*(g0 + pol*g1)
        #4) X (everything else hadronic)
        if z > 0.3:
            g0 = 1./0.30 #*heaviside(0.3-z)
            dndz = dndz + 0.13*g0
    return dndz

def dnnottaudE(z):
    #nue or numu from tau decay
    dndz = 0.18*(4.0 - 12.*z + 12.*z**2 - 4.*z**3)
    return dndz
