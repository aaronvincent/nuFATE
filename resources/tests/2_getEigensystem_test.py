#==========================================
# test 2 : getEigensystem call
# make sure if getEigensystem returns
# same answers for multiple calles 
#==========================================

import nuFATEpy as nuf
import numpy as np
import sys

flavor_id = -1
gamma_index = 2.2
h5_filename = "../NuFATECrossSections.h5"
include_secondaries = False
zenith = 180.* np.pi/180

nufate = nuf.nuFATE(flavor_id, gamma_index, h5_filename, include_secondaries)

print "generate nuFATE"

earth_t = nufate.get_earth_column_density(zenith)
print "earth_thickness[g/cm2] is ", earth_t

# check values.
result = nufate.get_eigensystem()
print "got result"


# test for NuEBar
index = 10
ntrials = 10
eigenvec = (result.eigenvectors())[index]
eigenval = (result.eigenvalues())[index]
coeff = (result.coefficients())[index]

print "start NuEBar test"

for i in range(ntrials) :
    r = nufate.get_eigensystem()
    e_val = (r.eigenvalues())[index]
    assert (eigenval == e_val), "eigenvalue didn't match, eval1 = %e, eval2 = %e" % (eigenval, e_val)
    ci = (r.coefficients())[index]
    assert (coeff == ci), "coefficient didn't match, c1 = %e, c2 = %e" % (coeff, ci)
    e_vec = (r.eigenvectors())[index]
    assert (eigenvec == e_vec), "eigenvector didn't match"

print "passed NuEBar test"

# test for NuTauBar
nufate2 = nuf.nuFATE(-3, gamma_index, h5_filename, include_secondaries)
result2 = nufate2.get_eigensystem()

eigenvec = (result2.eigenvectors())[index]
eigenval = (result2.eigenvalues())[index]
coeff = (result2.coefficients())[index]

print "start NuTauBar test"

for i in range(ntrials) :
    r = nufate2.get_eigensystem()
    e_val = (r.eigenvalues())[index]
    assert (eigenval == e_val), "eigenvalue didn't match, eval1 = %e, eval2 = %e" % (eigenval, e_val)
    ci = (r.coefficients())[index]
    assert (coeff == ci), "coefficient didn't match, c1 = %e, c2 = %e" % (coeff, ci)
    e_vec = (r.eigenvectors())[index]
    assert (eigenvec == e_vec), "eigenvector didn't match"

print "passed NuTauBar test"

print "all tests passed!"
