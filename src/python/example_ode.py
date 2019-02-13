""" Brief example of how to use nuFATE ODE solver
"""

import numpy as np
import cascade_ode as code
import init_ode
import matplotlib.pyplot as plt


#Choose the flavor & index you want
flavor = -2  # 1,2,3 = e, mu, tau; negative sign for antiparticles
gamma = 2.  # Power law index of isotropic flux E^-gamma
zenith = np.radians(120.)


def get_att_ode(flavor,g,zenith,E,h5_filename,prefactor=1e70,relerr=1e-4):
    #theta = np.radians(180. - zenith)
    theta = np.pi - zenith
    logE = np.log10(E)
    r_e = np.ones(len(E))
    r_mu = r_e
    r_tau = r_e
    tau = np.zeros(E.shape)
    if zenith <= np.pi/2:
        return r_e, r_mu, r_tau, tau
    elif zenith > np.pi:
        raise ValueError('Unphysical zenith angle, must be in between 90 and 180!')
    [RHSMatrix,energy_nodes,energy_tau] = init_ode.init(flavor,h5_filename)
    NumNodes = len(energy_nodes)
    NumTau = len(energy_tau)
    phi_nu = energy_nodes**(2 - g)
    #need a prefactor to avoid negative solutions, default value works for
    #g<=4, increase it if you want sofer spectrum

    phi_0 = np.concatenate((phi_nu, phi_nu, phi_nu, np.zeros(NumTau)),axis=0)*prefactor
    phi_out = code.cascade(theta,phi_0,RHSMatrix,energy_nodes,energy_tau,relerr)
    att_nue = phi_out[0:NumNodes]/phi_0[0:NumNodes]
    att_numu = phi_out[NumNodes:2*NumNodes]/phi_0[NumNodes:2*NumNodes]
    att_nutau = phi_out[2*NumNodes:3*NumNodes]/phi_0[2*NumNodes:3*NumNodes]
    att_tau = phi_out[3*NumNodes:]/energy_tau**2
    r_e = np.interp(logE,np.log10(energy_nodes),att_nue)
    r_mu = np.interp(logE,np.log10(energy_nodes),att_numu)
    r_tau = np.interp(logE,np.log10(energy_nodes),att_nutau)
    tau = np.interp(logE,np.log10(energy_tau), att_tau)/prefactor
    return r_e, r_mu, r_tau, tau

#prefactor to avoid negative solutions, absent to use default value 1e60
prefactor = 1e70
#relative error of the solutions, absent to use default value 1e-4
relerr = 1e-3
E=np.logspace(3,10,500)
h5_filename = '../../resources/NuFATECrossSections.h5'

[r_e,r_mu,r_tau,tau]=get_att_ode(flavor,gamma,zenith,E,h5_filename,prefactor,relerr)


plt.figure(figsize=(6,5))
plt.semilogx(E,r_e,c='r')
plt.semilogx(E,r_mu,c='g')
plt.semilogx(E,r_tau,c='b')
plt.xlabel(r"Neutrino Energy (GeV)")
plt.ylim(0.,1.1)
plt.ylabel(r"Attenuation")
plt.legend(["NuEBar", "NuMuBar", "NuTauBar"])
plt.grid()
plt.show()