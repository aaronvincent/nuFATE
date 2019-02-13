""" Brief example of how to use nuFATE ODE solver
"""

import numpy as np
import cascade_ode as code
import init_ode
import matplotlib.pyplot as plt


#Choose the flavor & index you want
flavor = -2  # negative sign for antiparticles, only sign matters in ODE solver
gamma = 2  # Power law index of isotropic flux E^-gamma
zenith = np.radians(120.)

#prefactor to avoid negative solutions, absent to use default value 1e60
prefactor = 1e70
#relative error of the solutions, absent to use default value 1e-4
relerr = 1e-3
h5_filename = '../../resources/NuFATECrossSections.h5'

def get_att_ode(flavor,gamma,zenith,h5_filename,prefactor=1e70,relerr=1e-3):

    [phi_0,RHSMatrix,energy_nodes,energy_tau] = init_ode.init(flavor,gamma,h5_filename,prefactor)
    [r_e, r_mu, r_tau, phi_tau] = code.cascade(zenith,phi_0,RHSMatrix,energy_nodes,energy_tau,prefactor,relerr)

    return energy_nodes, energy_tau, r_e, r_mu, r_tau, phi_tau

[energy_nodes, energy_tau, r_e,r_mu,r_tau,phi_tau]=get_att_ode(flavor,gamma,zenith,h5_filename,prefactor,relerr)


plt.figure(figsize=(6,5))
plt.semilogx(energy_nodes,r_e,c='r')
plt.semilogx(energy_nodes,r_mu,c='g')
plt.semilogx(energy_nodes,r_tau,c='b')
plt.xlabel(r"Neutrino Energy (GeV)")
plt.ylim(0.,1.1)
plt.ylabel(r"Attenuation")
plt.legend(["NuEBar", "NuMuBar", "NuTauBar"])
plt.grid()
plt.show()


plt.figure(figsize=(6,5))
plt.loglog(energy_tau,phi_tau,c='b')
plt.xlabel(r"Tau Energy (GeV)")
plt.ylabel(r"Tau Flux")
plt.grid()
plt.show()