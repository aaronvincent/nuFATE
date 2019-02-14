![nuFATE Logo](/resources/nufate.png)

nuFATE is a code that rapidly calculates the attenuated neutrino flux as a function of energy, neutrino flavor and zenith angle, due to the earth's opacity, for neutrino energies above 1 TeV. The python and matlab version of the the software contains both the classical nuFATE and a new ODE solver which takes care of the tau decay and energy loss in propogation and provides tau flux as a side product. The C++ version of ODE solver is not available yet. The software as implemented employs a user-specified power-law isotropic neutrino flux, the STW105 reference earth model, and neutrino-nucleon cross sections computed with the CT10nlo PDF distribution. The attenuation rates can be used to calculate the upgoing nu and nubar fluxes for high-energy neutrino observatories such as IceCube or ANTARES. Full description is available here: https://arxiv.org/abs/1706.09895

Prerequisites
-------------

The following packages are required to use the library, and
can probably be obtained from your favorite package manager:

* numpy: http://www.numpy.org/
* scipy: http://www.scipy.org/
* tables: https://pypi.python.org/pypi/tables
* h5py: https://www.h5py.org/

Recommended:
* ipython: http://ipython.org/
* jupyter: http://jupyter.readthedocs.io/
* matplotlib: http://matplotlib.org/

For the C++ version, you will need:

* hdf5 with c bindings: http://www.hdfgroup.org/HDF5/
* gsl (>= 1.15): http://www.gnu.org/software/gsl/
* C++ compiler with C++11 support


Compiling
---------

The Python interface can be installed by simply running:

  python setup.py install

Without write permissions, you can install it using:

  python setup.py install --user

The library can be compiled by running:

	make

An example program demonstrating usage and functionality
can be compiled with the command:

	make examples

The resulting example executables can then be found in the
subdirectory of `examples`

Finally the library can be installed using:

	make install


Example
-------

The one thing to remember is that the solution provided by nuFATE is E^2 * phi

The example script is called example.py. Additonal examples of how to use ODE solver can also be found in src/python folder. To run the example script do

python example.py

There is also an iPython notebook (notebook.ipynb) which shows some examples, including plots. This requires the "recommended" packages, above. notebook_ode.ipynb is also available to show new features of ODE solver.

To run the C++ example:

./example.exe

Citation
--------

If you want cite this work, or want to look at further description
please refer to

High-energy neutrino attenuation in the Earth and its associated uncertainties

Aaron C. Vincent, Carlos A. Arguelles, and A. Kheirandish

arXiv:1706.09895

Contributors
------------

- Aaron C. Vincent
- Carlos A. Arguelles
- Ali Kheirandish
- Ibrahim Safa
- Ningqiang Song

