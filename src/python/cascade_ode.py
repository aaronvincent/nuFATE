""" Funcitons to solve the neutrino and tau evolution using ODE solver
"""

import numpy as np
from scipy.integrate import solve_ivp
import earth

def cascade(zenith,phi_0,RHSMatrix,energynodes,energytau,prefactor=1e70,relerr=1e-3): 
    """ Returns the final neutrino and tau fluxes after propogating through the Earth with a specific zenith angle.

    Args:
        zenith: zenith angle in radians.
        phi_0: one dimensional numpy array containing the inital 3 neutrino+tau fluxes.
        RHSMatrix: three dimensional numpy array with the differential cross section cm**2 GeV**-1, prepared by the init_ode module (which contains the right-hand-side matices excluding density factor).        
        energynodes: one dimensional numpy array containing the energy nodes of neutrinos in GeV.
        energytau: one dimensional numpy array containing the energy nodes of tau in GeV.
        prefactor: a factor added to the initial fluxes to aviod negative solution at very large zenith angles, default 1e70           relerr: relative error for the ODE solver, default 1e-3.
    Returns:
        att_nue: one dimensional numpy array containing the attenuation factor of electron neutrino after propogation in Earth.
        att_numu: one dimensional numpy array containing the attenuation factor of muon neutrino after propogation in Earth.
        att_nutau: one dimensional numpy array containing the attenuation factor of tau neutrino after propogation in Earth.
        phi_tau: one dimensional numpy array containing the tau flux after propogation in Earth.
    """
    r_e = np.ones(len(energynodes))
    r_mu = r_e
    r_tau = r_e
    phi_tau = np.zeros(len(energytau))    
    if zenith <= np.pi/2:
        return r_e, r_mu, r_tau, phi_tau
    elif zenith > np.pi:
        raise ValueError('Unphysical zenith angle, must be in between pi/2 and pi.')    
    REarth=6371. #radius of the Earth in km
    energy_nodes = energynodes
    NumNodes = energy_nodes.size
    energy_tau = energytau
    theta = np.pi - zenith
    xmax = 2*REarth*np.cos(theta)*1e5 #in cm
    def odefun(x,phi): #integrand representing the RHSMatrix*phi
        return np.matmul(getdphidx(zenith,x,phi,RHSMatrix,energy_nodes,energy_tau),phi)
    def jac(x,phi): #jacobian which is the RHSMatrix
        return getdphidx(zenith,x,phi,RHSMatrix,energy_nodes,energy_tau)
    
    #one is free to reset the solver with a different method, relative or absolute error, step size, etc.
    sol = solve_ivp(odefun, [0, xmax], phi_0, method = 'BDF', jac = jac, rtol = relerr).y
    phisol = sol[:,-1]
    att_nue = phisol[0:NumNodes]/phi_0[0:NumNodes]
    att_numu = phisol[NumNodes:2*NumNodes]/phi_0[NumNodes:2*NumNodes]
    att_nutau = phisol[2*NumNodes:3*NumNodes]/phi_0[2*NumNodes:3*NumNodes]
    phi_tau = phisol[3*NumNodes:]/energy_tau**2/prefactor  

    return att_nue, att_numu, att_nutau, phi_tau
    
def getdphidx(zenith,x,phi,RHSMatrix,energy_nodes,energy_tau):
    """ Returns dphi/dx at a specific distance x.

    Args:
        zenith: zenith angle in radians.
        x: distance neutrinos have travelled in the Earth in cm.
        phi: one dimensional numpy array containing the 3 neutrino+tau fluxes at x (not used in this function).
        RHSMatrix: three dimensional numpy array with the differential cross section cm**2 GeV**-1, prepared by the init_ode module (which contains the right-hand-side matices excluding density factor).
        energy_nodes: one dimensional numpy array containing the energy nodes of neutrinos in GeV.
        energy_tau: one dimensional numpy array containing the energy nodes of tau in GeV.
    Returns:
        dphidx: numpy matrix dphidx at distance x which gives right-hand-side matrix in ODE
    """       
    NA = 6.02e+23 #Avogadro's constant
    c = 3.e+10 #speed of light in cm/s
    mtau = 1.777 #tau mass in GeV
    ttau = 2.906e-13 #tau lifetime in rest frame in s
    #define crust, mantle and core, it depends on the modeling of the Earth
    REarth = 6371.*1e5 #radius of the Earth in cm
    nodes_earth = np.array([0., 1221., 3480., 5721., 5961., 6347., 6356.,6368., REarth])*1e5 #Earth layers
    Rcrust = nodes_earth[5]
    Rmantle=nodes_earth[2]
    cth = np.cos(np.pi - zenith)
    xtoR = x/REarth
    r = np.sqrt(1 + xtoR**2 - 2.*xtoR*cth)
    radius = REarth*r
    rho = earth.rho_earth(zenith, x/1e5) #in g/cm^3
    
    NumNodes = len(energy_nodes)
    NumTau = len(energy_tau)
    #restore matrix elements
    RHSMatrix_11=RHSMatrix[0,:,:]
    RHSMatrix_22=RHSMatrix[1,:,:]
    RHSMatrix_33=RHSMatrix[2,:,:]
    if radius>=Rcrust: #crust
        RHSMatrix_44 = RHSMatrix[3,:,:]
    elif radius >= Rmantle: #mantle
        RHSMatrix_44 = RHSMatrix[4,:,:]
    else: #core
        RHSMatrix_44 = RHSMatrix[5,:,:]
    RHSMatrix_44 = RHSMatrix_44 - np.diag(mtau/(c*ttau*rho*NA*energy_tau)) #tau decay
    RHSMatrix_43 = RHSMatrix[6,:,:]
    RHSMatrix_34 = RHSMatrix[7,:,:]*mtau/(c*ttau*rho*NA)
    RHSMatrix_14 = RHSMatrix[8,:,:]*mtau/(c*ttau*rho*NA)
    RHSMatrix_24 = RHSMatrix_14
    #construct the huge RHS Matrix
    z0 = np.zeros((NumNodes,NumNodes))
    z0tau = np.zeros((NumTau,NumNodes))
    M = np.block([[RHSMatrix_11, z0, z0, RHSMatrix_14],
       [z0, RHSMatrix_22, z0, RHSMatrix_24],
       [z0, z0, RHSMatrix_33, RHSMatrix_34],
       [z0tau, z0tau, RHSMatrix_43, RHSMatrix_44]])
    dphidx = M*rho*NA
    return dphidx