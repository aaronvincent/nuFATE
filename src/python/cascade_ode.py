""" Funcitons to solve the neutrino and tau evolution using ODE solver
"""

import numpy as np
from scipy.integrate import solve_ivp
import earth

def cascade(theta,phi_0,RHSMatrix,energynodes,energytau,relerr=1e-3): 
    """ regeneration bits are included here and three nus+tau are solved at the 
        same time. For treatment of on-the-spot decay of tau or without
        regenerations, refer to the old code
        input: length of path in a sublayer-x (in cm), average density the 
        sublayer-rhobar (in g*cm^-3), input flux-phi_0, radius of layer-radius, 
        RHS matrix elements-RHSMatrix, energy nodes-energy_nodes
        output: output flux-phi_1, energy nodes-energy_nodes
    """
    REarth=6371. #radius of the Earth in km
    energy_nodes = energynodes
    energy_tau = energytau
    xmax = 2*REarth*np.cos(theta)*1e5 #in cm
    def odefun(x,phi):
        return getdphidx(theta,x,phi,RHSMatrix,energy_nodes,energy_tau)
    def jac(x,phi):
        return jacobian(theta,x,phi,RHSMatrix,energy_nodes,energy_tau)

    sol = solve_ivp(odefun, [0, xmax], phi_0, method = 'BDF', jac = jac, rtol = relerr).y
    phisol = sol[:,-1]

    return phisol
    
def getdphidx(theta,x,phi,RHSMatrix,energy_nodes,energy_tau):
    NA = 6.02e+23 #Avogadro's constant
    c = 3.e+10 #speed of light in cm/s
    mtau = 1.777 #tau mass in GeV
    ttau = 2.906e-13 #tau lifetime in rest frame in s
    #define crust, mantle and core, it depends on the modeling of the Earth
    REarth = 6371.*1e5 #radius of the Earth in cm
    nodes_earth = np.array([0., 1221., 3480., 5721., 5961., 6347., 6356.,6368., REarth])*1e5
    Rcrust = nodes_earth[5]
    Rmantle=nodes_earth[2]
    cth = np.cos(theta)
    xtoR = x/REarth
    r = np.sqrt(1 + xtoR**2 - 2.*xtoR*cth)
    radius = REarth*r
    rho = earth.rho_earth(np.pi - theta, x/1e5) #in g/cm^3
    
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
    return np.matmul(M, phi)*rho*NA

def jacobian(theta,x,phi,RHSMatrix,energy_nodes,energy_tau):
    NA = 6.02e+23 #Avogadro's constant
    c = 3.e+10 #speed of light in cm/s
    mtau = 1.777 #tau mass in GeV
    ttau = 2.906e-13 #tau lifetime in rest frame in s
    #define crust, mantle and core, it depends on the modeling of the Earth
    REarth = 6371.*1e5 #radius of the Earth in cm
    nodes_earth = np.array([0., 1221., 3480., 5721., 5961., 6347., 6356.,6368., REarth])*1e5
    Rcrust = nodes_earth[5]
    Rmantle=nodes_earth[2]
    cth = np.cos(theta)
    xtoR = x/REarth
    r = np.sqrt(1 + xtoR**2 - 2.*xtoR*cth)
    radius = REarth*r;
    rho = earth.rho_earth(np.pi - theta, x/1e5) #in g/cm^3
    
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
    return M*rho*NA
