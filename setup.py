from setuptools import setup

packages = [
    'grid_cell_model',
]

entry_points = {
    'console_scripts': [
        'gcm_hello_world = grid_cell_model.entry_points.gcm_hello_world:main',
        'sweepls = grid_cell_model.entry_points.sweepls:main'
    ],
}

setup(
    name='grid_cell_model',
    version='0.1.0-dev',
    description='Simulations of grid cell networks.',
    author='Lukas Solanka',
    author_email='l.solanka@sms.ed.ac.uk',
    packages=packages,
    install_requires=['numpy>=1.8.0',
                      'scipy>=0.13.3',
                      'matplotlib>=1.3.01',
                      'h5py>=2.3.0',
                      'gridcells==0.1.2',
                      'enum34>=1.0.4'],
    entry_points=entry_points,
)
