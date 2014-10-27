from setuptools import setup, find_packages


setup(
    name='noisefigs',
    version='0.1.0-dev',
    description='Simulations of grid cell networks.',
    author='Lukas Solanka',
    author_email='l.solanka@sms.ed.ac.uk',
    packages=find_packages(),
    install_requires=['numpy>=1.8.0',
                      'scipy>=0.13.3',
                      'matplotlib>=1.3.01',
                      'h5py>=2.3.0']
)
