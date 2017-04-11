import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate

RSun = 695980.0  # km

_sun_model = np.genfromtxt(
    "/Users/carlos/Programs/nuSQuIDS/data/astro/bs05_agsop.dat")

_sun_radius = _sun_model[:, 1]
_sun_density = _sun_model[:, 2]
_inter_density = interpolate.interp1d(_sun_radius, _sun_density)

def rho_sun(x, b):
    """ Returns the Sun density in gr/cm^3.

    Args:
        x: position along the trayectory in km.
        b: impact parameter in km.

    Returns:
        rho: density in kg/m^2
    """
    r = np.sqrt(RSun**2 + x**2 - 2 * x * np.sqrt(RSun**2 - b**2)) / RSun
    if r <= _sun_radius[0]:
        return _sun_density[0]
    elif r > _sun_radius[-1]:
        return 0.
    else:
        return _inter_density(r)


def get_t_sun(b):
    """ Returns the Sun column density for a given impact parameter.

    Args:
        b: impact parameter in km.

    Returns:
        rho: density in gr/cm^2
    """

    if (b < 0. or b > Rsun):
        raise NameError(
            "Impact parameter cannot be negative or larger than the Sun radius.")

    # returns the Sun column density in gr/cm^2
    kmTocm = 1.0e5
    xmax = 2. * np.sqrt(RSun**2 - b**2)
    return integrate.quad(
        lambda x: rho_sun(x, b), 0, xmax, epsrel=1.0e-3,
        epsabs=1.0e-18)[0] * kmTocm
