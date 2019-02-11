 /******************************************************************************
 *    This program is free software: you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                             *
 ******************************************************************************/

#define H5Gopen_vers 2
#define H5Gcreate_vers 2
#define H5Eset_auto_vers 2
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/list.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/overloads.hpp>
#include <nuFATE/nuFATE.h>
#include <vector>
#include <iostream>

using namespace boost::python;
using namespace nufate;
namespace np = boost::python::numpy;

// converting 1D vector array to python list
template<class T>
struct VecToList
{
   static PyObject* convert(const std::vector<T>& vec) {
      boost::python::list* l = new boost::python::list();
      for (size_t i =0; i < vec.size(); i++) {
         (*l).append(vec[i]);
      }
      return l->ptr();
   }
};

// converting Result::Square_matrix_double object to python list (NxN dimension)
struct VecToList2D 
{
   static PyObject* convert( const Result::Square_matrix_double &vec) {
      unsigned int dim = vec.dim_;
      boost::python::list* l = new boost::python::list();
      for (unsigned int i = 0; i < dim; i++) {
         boost::python::list* ll = new boost::python::list();
         for (unsigned int j = 0; j<dim; j++) { 
            double val = *(vec.evec_.get() + i*dim + j);
            (*ll).append(val);
         }
         (*l).append(*ll);
      }
      return l->ptr();
    }
};


// nuFATEpy module definitions
BOOST_PYTHON_MODULE(nuFATEpy)
{
  to_python_converter< std::vector<double, class std::allocator<double> >, VecToList<double> > ();
  to_python_converter< Result::Square_matrix_double, VecToList2D> ();

  class_<Result>("Result")
    .def(init<>())
    .def("eigenvalues", &Result::get_eigenvalues)
    .def("eigenvectors", &Result::get_eigenvec_matrix)
    .def("coefficients", &Result::get_coefficients)
    .def("energy_nodes", &Result::get_energy_nodes)
    .def("phi_0", &Result::get_phi_0)
    .def("Print", &Result::Print, arg("index"))
  ;

  class_<nuFATE, boost::noncopyable, std::shared_ptr<nuFATE> >("nuFATE",no_init)
    .def(init<int, double, std::string, bool>(args("flavor_id","gamma_index","h5_filename","include_secondaries")))
    .def(init<int, double, std::vector<double>, std::vector<double>, std::vector<std::vector<double> >, bool>(args("flavor_id","gamma_index","energy_nodes","sigma_array","dsigma_dE","include_secondaries")))
    .def("get_eigensystem", &nuFATE::getEigensystem)
    .def("get_earth_column_density", &nuFATE::getEarthColumnDensity, arg("zenith"))
    .def_readonly("flavor", &nuFATE::getFlavor)
    .def_readonly("gamma", &nuFATE::getGamma)
    .def_readonly("numnodes", &nuFATE::getNumNodes)
    .def("set_add_secondaries",&nuFATE::setAddSecondaries, arg("opt"))
  ;

}
