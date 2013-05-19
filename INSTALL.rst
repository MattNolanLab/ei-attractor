-------------------------------------------------
Installation instructions for the grid cell model
-------------------------------------------------

In fact the model does not need any installation, it is recommended to work
directly from the source directory. However there are a number of
Prerequisites_ that the user must install in order to start simulation.

Prerequisites
=============

- Python          2.6/2.7
- numpy           >= 1.6.2
- scipy           >= 0.11.0
- matplotlib      >= 1.1.1
- brian simulator >= 1.4.0

All of these should be possible to install with your favorite installer (pip,
easy_install, ...)


Running a simulation
====================

All the simulation scripts are in the grid_cell_model/ directory. It is
advisable to inspect the simulation/submitting files before running.

The main files are described below:

*submit_\*.py*
    A script that submits a particular simulation run(s). It can run a single
    simulation or one can define a parameter sweep over a given range of
    parameters.

*simulation_\*.py*
    A script that performs the actual simulation run. It sets up the model and
    the parameters, runs the simulation and extracts all the necessary output
    data, and usually saves the data into a user-specified file format
    (currently a MAT file) or plots them and saves the figures.

*default_params.py*
    Contains default parameters for this particular set of simulations. Each
    specific setup within the simulation directory can then override these
    default parameters.


The submitting scripts currently support submitting to a workstation or to Sun
Grid Engine cluster (using qsub). This can be selected inside the submit_*.py
script. Also, the workstation runs can be submitted in a blocking mode (each
successive simulation waits until the previous one has finished) or in a
non-blocking mode (all simulations are run at once concurrently). See the
blocking parameter in these scripts.
