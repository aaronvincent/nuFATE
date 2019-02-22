
import nuFATEpy as nuf
import numpy as np
import sys

flavor_id = 2
gamma_index = 2.2
h5_filename = "../../../resources/NuFATECrossSections.h5"
include_secondaries = False
zenith = 180.* np.pi/180

nufate = nuf.nuFATE(flavor_id, gamma_index, h5_filename, include_secondaries)

print "generate nuFATE"

earth_t = nufate.get_earth_column_density(zenith)
print "earth_thickness[g/cm2] is ", earth_t

# check values.
result = nufate.get_eigensystem()
print "got result"

# debug print CXX function 
index = 10;

eigenvec = (result.eigenvectors())[index]
eigenval = (result.eigenvalues())[index]
coeff = (result.coefficients())[index]

for i in range(10) :
    r = nufate.get_eigensystem()
    e_val = (r.eigenvalues())[index]
    assert (eigenval == e_val), "eigenvalue didn't match, eval1 = %e, eval2 = %e" % (eigenval, e_val)
    ci = (r.coefficients())[index]
    assert (coeff == ci), "coefficient didn't match, c1 = %e, c2 = %e" % (coeff, ci)
    e_vec = (r.eigenvectors())[index]
    assert (eigenvec == e_vec), "eigenvector didn't match"


