from setuptools import setup, find_packages

with open("README.md") as f:
    long_description = f.read()

setup(
      name='nuFATE',
      version='1.0',
      description='Calculates relative attenuation of a neutrino flux through the Earth',
      long_description=long_description,
      package_dir = {'': 'src/python'},
      py_modules = ['cascade', 'cascade_secs', 'earth', 'sun'],
      install_requires=['numpy', 'tables','scipy'],
      author='Aaron Vincent, Ali Kheirandish, Carlos Arguelles, Ibrahim Safa',
      url='git@github.com:aaronvincent/nuFATE.git'
        )
