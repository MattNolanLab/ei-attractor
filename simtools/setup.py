'''Setup script for simtools.'''
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

all_packages = ['simtools',
                'simtools.plotting',
                'simtools.storage']

setup(name='simtools',
      version='0.1.0',
      description='Python simulation and analysis helpers.',
      author='Lukas Solanka',
      author_email='lsolanka@gmail.com',
      packages=all_packages,
      install_requires=['numpy >= 1.8.1',
                        'matplotlib >= 1.3.1',
                        'configobj >= 5.0.6',
                        'h5py >= 2.4.0',
                        'six >= 1.9.0'],
      )

