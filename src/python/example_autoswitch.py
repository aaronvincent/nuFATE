""" Example of switching betwen nuFATE and ODE solver
"""

import numpy as np
import cascade as cas
import cascade_secs as csx
import cascade_ode as code
import earth
import init_ode
import matplotlib.pyplot as plt

#Choose the flavor & index you want
f = -1  # negative sign for antiparticles, only sign matters
gamma = 2  # Power law index of isotropic flux E^-gamma
zenith = np.radians(100.)

#specify an angle to switch between classical nuFATE and ODE solver, below that angle use ODE solver while above that angle use nuFATE instead to speed up
switchangle = np.radians(110.) #recommended
#whether or not if you want tau flux as output, this is only achievable in ODE solver
tauflag = False
if tauflag:
    switchangle = np.radians(180.)

def get_att_switch(switchangle,f,gamma,zenith,h5_filename):
    if zenith <= switchangle: #use ODE solver
        prefactor = 1e70
        relerr = 1e-3
        [phi_0,RHSMatrix,energy_nodes,energy_tau] = init_ode.init(f,gamma,h5_filename,prefactor)
        [r_e, r_mu, r_tau, phi_tau] = code.cascade(zenith,phi_0,RHSMatrix,energy_nodes,energy_tau,prefactor,relerr)
        return energy_nodes, energy_tau, r_e, r_mu, r_tau, phi_tau
    else: #use classical nuFATE
        if f > 0:
            flavor = [1, 2, 3]
        else:
            flavor = [-1,-2,-3]
        for i in np.abs(flavor):
            if i==3:
                w, v, ci, energy_nodes, phi_0 = cas.get_eigs(flavor[i-1], gamma, h5_filename)
            else:
                w, v, ci, energy_nodes, phi_0 = csx.get_eigs(flavor[i-1], gamma, h5_filename)
            Na = 6.022e23
            t = earth.get_t_earth(zenith) * Na
            phisol = np.dot(v, (ci * np.exp(w * t))) / phi_0  #this is phi/phi_0, i.e. the relative attenuation
            if i == 1:
                r_e = phisol[0:energy_nodes.size]
            elif i == 2:
                r_mu = phisol[0:energy_nodes.size]
            else:
                r_tau = phisol[0:energy_nodes.size]
        return energy_nodes, energy_nodes, r_e, r_mu, r_tau, np.zeros(energy_nodes.size)
    
h5_filename = '../../resources/NuFATECrossSections.h5'
energy_nodes, energy_tau, r_e, r_mu, r_tau, phi_tau = get_att_switch(switchangle,f,gamma,zenith,h5_filename)

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