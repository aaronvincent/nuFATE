#Quick example: how to get attenuations from nuFATE
# Aaron Vincent 2017
import numpy as np
from cascade import *
from rho_earth import *

#These are made to agree with precomputed cross sections in data/ folder;
#do not change them without recomputing and retabulating the cross sections
logemin = 3
logemax = 10
NumEnodes= 200

#Choose the flavor & index you want
flavor = -2 # 1,2,3 = e, mu, tau; negative sign for antiparticles
gamma = 2.2 # Power law index of isotropic flux E^-gamma

#solve the cascade equation once
w,v,ci,energy_nodes = get_eigs(flavor,gamma,logemin,logemax,NumEnodes)


#this function just interpolates the solution
def get_att_value(w,v,ci,energy_nodes,zenith,E):
    Na = 6.022e23
    theta = (zenith-180)/180*np.pi
    logE = np.log10(E)
    t = get_t(theta)*Na*1000/100**2; # g/ cm^2
#    phi = np.dot(v,(ci*np.exp(w*t)))*energy_nodes**(-2) #this is the attenuated flux
    phisol = np.dot(v,(ci*np.exp(w*t)))*energy_nodes**(gamma-2) #this is phi/phi_0, i.e. the relative attenuation
    return np.interp(logE,np.log10(energy_nodes),phisol)


#specify a zenith angle and energy you are interested in
zenith = 130. #degrees
E = 100e3 #GeV
att = get_att_value(w,v,ci,energy_nodes,zenith,E)

print "Flux at E  =", E, " GeV , zenith = ", zenith, " degrees will be attenuated by a factor of ", att
#done

# The next section shows how to include secondary electron, muon neutrinos
import cascade_secs as csx
w,v,ci,energy_nodes,phi_0 = csx.get_eigs(flavor,gamma,logemin,logemax,NumEnodes) #note the extra argument

def get_att_value_secs(w,v,ci,energy_nodes,zenith,E):
    Na = 6.022e23
    theta = (zenith-180)/180*np.pi
    logE = np.log10(E)
    t = get_t(theta)*Na*1000/100**2; # g/ cm^2
    #    phi = np.dot(v,(ci*np.exp(w*t)))*energy_nodes**(-2) #this is the attenuated flux
    phisol = np.dot(v,(ci*np.exp(w*t)))/phi_0 #this is phi/phi_0, i.e. the relative attenuation
    phisol = phisol[0:200] #the non-tau bit.
    return np.interp(logE,np.log10(energy_nodes),phisol)

att = get_att_value_secs(w,v,ci,energy_nodes,zenith,E)

print "If I include secondaries, that factor becomes: ", att

#done done


