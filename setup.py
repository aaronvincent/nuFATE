from setuptools import setup

with open("README.md") as f:
    long_description = f.read()

setup(
      name='nuFATE',
      version='1.0',
      description='Calculates relative attenuation of a neutrino flux through the Earth',
      long_description=long_description,
      packages=['src/python'],
      install_requires=['numpy', 'tables','scipy'],
      author='Aaron Vincent, Ali Kheirandish, Carlos Arguelles, Ibrahim Safa',
      url='git@github.com:aaronvincent/nuFATE.git'
        )
