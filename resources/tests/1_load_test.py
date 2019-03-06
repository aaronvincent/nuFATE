#===================================================
# test 1 : load test
# Is nuFATE object loading cross sections correctly?
#===================================================

import nuFATEpy as nuf
import numpy as np
import sys
import tables

flavor_id = 2 # NuMu
gamma_index = 2.2
include_secondaries = False
zenith = 180.* np.pi/180
Na = 6.0221415e23

#---------------------------------------
# test 1:
# load test for default constructor
#---------------------------------------

# default h5 cross section 
# Energy range : 10^3 to 10^10 GeV
# number of energy bins : 200
h5_filename = "../NuFATECrossSections.h5"
xsh5 = tables.open_file(h5_filename,"r")
logemax = np.log10(xsh5.root.total_cross_sections._v_attrs.max_energy)
logemin = np.log10(xsh5.root.total_cross_sections._v_attrs.min_energy)
numnodes = xsh5.root.total_cross_sections._v_attrs.number_energy_nodes
default_enodes = np.logspace(logemin, logemax, numnodes)
default_sigma_array = xsh5.root.total_cross_sections.numuxs[:]
default_dxs_array = xsh5.root.differential_cross_sections.dxsnu[:]

nufate1 = nuf.nuFATE(flavor_id, gamma_index, h5_filename, include_secondaries)
enodes = np.array(nufate1.energy_nodes())
np.testing.assert_allclose(enodes, default_enodes, rtol=1e-5, atol=0, err_msg="test1, energy_nodes didn't match")

sigma_array = np.array(nufate1.total_crosssections())
np.testing.assert_allclose(sigma_array, default_sigma_array, rtol=1e-5, atol=0, err_msg="test1, sigma_array didn't match")

dsigma_array = np.array(nufate1.nc_differential_crosssections())
np.testing.assert_allclose(default_dxs_array, dsigma_array, rtol=1e-5, atol=0, err_msg="test1, dxs_array didn't match")

earth_t = nufate1.get_earth_column_density(zenith)
phi_sol1 = np.array(nufate1.get_relative_attenuation(earth_t*Na))

phi_sol1_true = np.loadtxt("load_test_out1.txt")
np.testing.assert_allclose(phi_sol1, phi_sol1_true, rtol=1e-5, atol=0, err_msg="test1, phi_sol didn't match")

#---------------------------------------
# test 2:
# load test for 3rd constructor
#---------------------------------------

# these are nuSQuIDS cross sections (csms)
# Energy range : 10^2 to 10^9 GeV
# number of energy bins : 501 
cc_sigma_file = "../nuSQuIDSCrossSections/nusigma_sigma_CC.dat"
nc_sigma_file = "../nuSQuIDSCrossSections/nusigma_sigma_NC.dat"
nc_dsigma_dE_file = "../nuSQuIDSCrossSections/nusigma_dsde_NC.dat"
include_secondaries = False
logemin = 2
logemax = 9
numnodes = 501

default_enodes = np.logspace(logemin, logemax, numnodes)
ccx = np.loadtxt(cc_sigma_file)
ncx = np.loadtxt(nc_sigma_file)
numu_index = 3 # 0:ene, 1:nue, 2:nuebar, 3:numu
default_sigma_array = ccx[:,numu_index] + ncx[:,numu_index] 
dsde = np.loadtxt(nc_dsigma_dE_file)
numu_index = 4 # 0:enein, 1:eneout, 2:nue, 3:nuebar, 4:numu
default_dxs_array = dsde[:,numu_index]
default_dxs_array = np.reshape(default_dxs_array, (numnodes,numnodes))

nufate2 = nuf.nuFATE(flavor_id, gamma_index, cc_sigma_file, nc_sigma_file, nc_dsigma_dE_file, include_secondaries)

enodes = np.array(nufate2.energy_nodes())
np.testing.assert_allclose(enodes, default_enodes, rtol=1e-5, atol=0, err_msg="test2, energy_nodes didn't match")

sigma_array = np.array(nufate2.total_crosssections())
np.testing.assert_allclose(sigma_array, default_sigma_array, rtol=1e-5, atol=0, err_msg="test1, sigma_array didn't match")

dsigma_array = np.array(nufate2.nc_differential_crosssections())
np.testing.assert_allclose(default_dxs_array, dsigma_array, rtol=1e-5, atol=0, err_msg="test1, dxs_array didn't match")

earth_t = nufate2.get_earth_column_density(zenith)
phi_sol2 = np.array(nufate2.get_relative_attenuation(earth_t*Na))

phi_sol2_true = np.loadtxt("load_test_out2.txt")
np.testing.assert_allclose(phi_sol2, phi_sol2_true, rtol=1e-5, atol=0, err_msg="test2, phi_sol didn't match")

print "all tests passed!"
