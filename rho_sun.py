import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate

Rsun=695980.0 # km

sun_model=np.genfromtxt("/Users/carlos/Programs/nuSQuIDS/data/astro/bs05_agsop.dat")

sun_radius=sun_model[:,1]
sun_density=sun_model[:,2]
inter_density=interpolate.interp1d(sun_radius,sun_density)

def rho_sun(x,b):
    # returns the Sun density in gr/cm^3
    r = np.sqrt(Rsun**2+x**2-2*x*np.sqrt(Rsun**2-b**2))/Rsun
    if r <= sun_radius[0]:
        return sun_density[0]
    elif r > sun_radius[-1]:
        return 0.
    else:
        return inter_density(r)

def get_t_sun(b):
    if(b<0. or b>Rsun):
        
    # returns the Sun column density in gr/cm^2
    kmTocm=1.0e5
    xmax=2.*np.sqrt(Rsun**2-b**2)
    return integrate.quad(lambda x: rho_sun(x,b),0,xmax,epsrel=1.0e-3 ,epsabs=1.0e-18)[0]*kmTocm