import numpy as np
import scipy as sp
from numpy import linalg as LA


def get_RHS_matrices(energy_nodes,sigma_fname,dxs_fname):
    NumNodes = energy_nodes.shape[0]
    sigma_array = np.loadtxt(sigma_fname)
    dsigmady = np.loadtxt(dxs_fname)
    DeltaE = np.diff(np.log(energy_nodes))
    RHSMatrix = np.zeros((NumNodes,NumNodes))
    # fill in diagonal terms
    
    for i in range(NumNodes): #E_out
        for j in range(i+1,NumNodes): #E_in
            RHSMatrix[i][j] = DeltaE[j-1]*dsigmady[j][i]*energy_nodes[j]**-1*energy_nodes[i]**2
    return RHSMatrix, sigma_array

def get_eigs(flavor,gamma = 2.,logemin=3,logemax=10,NumNodes=200):
#    print flavor
    if flavor==-1:
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

    if flavor > 0 :
        dxs_fname = "data/dxsnu.dat"
    else:
        dxs_fname = "data/dxsnubar.dat"
    
    #Note that the solution is scaled by E^2; if you want to modify the incoming spectrum a lot, you'll need to change this here, as well as in the definition of RHS.
    energy_nodes = np.logspace(logemin,logemax,NumNodes)
    RHSMatrix, sigma_array = get_RHS_matrices(energy_nodes,sigma_fname,dxs_fname)
    
    
    #tau regenration
    if flavor ==-3:
        RHregen,s1 = get_RHS_matrices(energy_nodes,sigma_fname,"data/tbarfull.dat")
        RHSMatrix = RHSMatrix + RHregen

    elif flavor == 3:
        RHregen,s1 = get_RHS_matrices(energy_nodes,sigma_fname,"data/tfull.dat")
        RHSMatrix = RHSMatrix + RHregen


    phi_0 = energy_nodes**(2-gamma)
    w,v = LA.eig(-np.diag(sigma_array)+RHSMatrix)
    ci = LA.solve(v,phi_0)# alternative to lstsq solution
#    ci = LA.lstsq(v,phi_0)[0]

    return w,v,ci,energy_nodes

