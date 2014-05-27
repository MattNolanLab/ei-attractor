from distutils.core import setup

setup(
        name='grid_cell_model',
        version='0.1.0dev',
        description='Simulations of grid cell networks.',
        author='Lukas Solanka',
        author_email='l.solanka@sms.ed.ac.uk',
        packages=['grid_cell_model'],
        install_requires=['numpy>=1.8.0',
                          'scipy>=0.13.3',
                          'matplotlib>=1.3.01',
                          'h5py>=2.3.0']
)

