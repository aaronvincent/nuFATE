import numpy as np
import scipy.integrate as integrate
Re = 6371.  # Earth radius


def rho(theta, x):
    #	theta = angle from down vector (0 = direction north pole...if you're at IceCube)
    # piecewise polynomial fit to Reference earth model STW105
    # you could also load a Ref earth model if you want.

    r = np.sqrt(Re**2 + x**2 - 2. * Re * x * np.cos(theta))

    if r < 1221.:
        p1 = -0.0002177
        p2 = -4.265e-06
        p3 = 1.309e+04
    elif r < 3480:
        p1 = -0.0002409
        p2 = 0.1416
        p3 = 1.234e+04
    elif r < 5721:
        p1 = -3.764e-05
        p2 = -0.1876
        p3 = 6664
    elif r < 5961:
        p1 = 0.
        p2 = -1.269
        p3 = 1.131e+04
    elif r < 6347:
        p1 = 0.
        p2 = -.725
        p3 = 7887.
    elif r < 6356:
        p1 = 0
        p2 = 0
        p3 = 2900
    elif r < 6368:
        p1 = 0
        p2 = 0
        p3 = 2600
    else:
        p1 = 0
        p2 = 0
        p3 = 1020

    rho = p1 * r**2 + p2 * r + p3

    return rho  # kg/m^3. Remember units if you integrate along the l.o.s.; r/x is in km


def get_t(theta):
    xmax = 2 * Re * np.cos(theta)
    if xmax <= 0:
        t = 0
    else:
        n = lambda x: rho(theta, x)  #mass density
        t = integrate.quad(
            lambda x: n(xmax - x), 0, xmax, epsrel=1.0e-3,
            epsabs=1.0e-18)[0] * 1000.  #kg/m^2
    return t
