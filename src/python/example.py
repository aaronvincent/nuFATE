""" Brief example of how to use nuFATE
"""

import numpy as np
import cascade as cas
import earth

#Choose the flavor & index you want
flavor = 2  # 1,2,3 = e, mu, tau; negative sign for antiparticles
gamma = 2.  # Power law index of isotropic flux E^-gamma
ReverseTime = False #You want to go backwards or forward? True for backwards, false for forward in time
Efinal = 0.5e9 #If you're going backwards in time, set the final neutrino energy you're trying to 'unfold'

#gamma = 'data/phiHGextrap.dat' #This is an example Honda Gaisser atmospheric flux. You can use this or add your own file, being careful to follow the energy spacing

#solve the cascade equation once
w, v, ci, energy_nodes, phi_0 = cas.get_eigs(flavor, gamma, "../../resources/NuFATECrossSections.h5", ReverseTime, Efinal)


#this function just interpolates the solution
def get_att_value(w, v, ci, energy_nodes, zenith, E,phi_in):
    Na = 6.0221415e23
    logE = np.log10(E)
    t = earth.get_t_earth(zenith) * Na
    # g/ cm^2
    #    phi = np.dot(v,(ci*np.exp(w*t)))*energy_nodes**(-2) #this is the attenuated flux
    if(ReverseTime):
        t = -1.*t
        phisol = np.dot(v, (ci * np.exp(w * t))) * energy_nodes**(-2)
        #print phisol
    else:
        phisol = np.dot(v, (ci * np.exp(w * t))) * energy_nodes**(-2) / phi_in #this is phi/phi_inital, i.e. the relative attenuation
    return np.interp(logE, np.log10(energy_nodes), phisol)


#specify a zenith angle and energy you are interested in
zenith = np.radians(120.) # zenith angle in radians
E = 100e3  #GeV
att = get_att_value(w, v, ci, energy_nodes, zenith, E,energy_nodes**-gamma)

print "Flux at E  =", E, " GeV , zenith = ", np.degrees(zenith), " degrees will be attenuated by a factor of ", att
#done

# The next section shows how to include secondary electron, muon neutrinos
import cascade_secs as csx
w, v, ci, energy_nodes, phi_0 = csx.get_eigs(
    flavor, gamma, "../../resources/NuFATECrossSections.h5")


def get_att_value_secs(w, v, ci, energy_nodes, zenith, E):
    Na = 6.022e23
    logE = np.log10(E)
    t = earth.get_t_earth(zenith) * Na
    # g/ cm^2
    #    phi = np.dot(v,(ci*np.exp(w*t)))*energy_nodes**(-2) #this is the attenuated flux
    phisol = np.dot(v, (ci * np.exp(w * t))
                   ) / phi_0  #this is phi/phi_0, i.e. the relative attenuation
    phisol = phisol[0:200]  #the non-tau bit.
    return np.interp(logE, np.log10(energy_nodes), phisol)


att = get_att_value_secs(w, v, ci, energy_nodes, zenith, E)

print "If I include secondaries, that factor becomes: ", att

#done done
