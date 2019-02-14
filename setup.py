from setuptools import setup, find_packages

with open("README.md") as f:
    long_description = f.read()

setup(
      name='nuFATE',
      version='1.0',
      description='Calculates relative attenuation of a neutrino flux through the Earth',
      long_description=long_description,
      packages=find_packages('./src/python'),
      install_requires=['numpy', 'tables','scipy','h5py'],
      author='Aaron Vincent, Ali Kheirandish, Carlos Arguelles, Ibrahim Safa','Ningqiang Song',
      url='git@github.com:aaronvincent/nuFATE.git'
        )
