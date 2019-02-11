
import nuFATEpy as nuf
import numpy as np
import sys

flavor_id = 2
gamma_index = 2.2
h5_filename = "../../resources/NuFATECrossSections.h5"
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
print "*** output from C++ function ***"
result.Print(index);

# debug print pybinding 

print ""
print "*** output from pybinding ***"
print ("eigenvec %d th row is " % (index))
print (result.eigenvectors())[index]
print ""
print ("eigenvalue[%d] is " % (index))
print (result.eigenvalues())[index]
print ""
print ("ci[%d] is " % (index))
print (result.coefficients())[index]
print ""
print ("energy_node[%d] is " % (index))
print (result.energy_nodes())[index]


