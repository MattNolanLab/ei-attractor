-------------------------------------------------
Installation instructions for the grid cell model
-------------------------------------------------

In fact the model does not need any installation, it is recommended to work
directly from the source directory. However there are a number of
Prerequisites_ that the user must install in order to start simulation.

Also it is recommended to include the grid_cell_model directory into your
PYTHONPATH environment variables, as all the submitting scripts reference
modules in that directory.


Prerequisites
=============

The model requires some common prerequisites:

- Python 2.6/2.7
- numpy  >= 1.6.2
- scipy  >= 0.11.0
- ideally matplotlib >= 1.1.1
- in some cases the hdf5 library, but this is rather optional. If you don't
  want to use hdf5 in your simulation scripts, you will need to set the
  correct output data format (this is currently in a TODO state).
- Currently, the model also requires the gcc C++ compiler, as some of the
  network setup routines are implemented using C++ code generator from the
  scipy package (weave) to speed up network setup. Also, some data analysis
  tools use this code generation technique.

All of these should be possible to install with your favorite installer (pip,
easy_install, ...)


Brian version
-------------

If you are planning to use brian_ (which is obsolete however), you will need to
install the `brian simulator`_ (>= 1.4.0).

However, the brian version is no longer maintained.


NEST version
------------

If you are planning to use NEST_ for simulations, you need these extra
programs:

- C++ compiler (gcc, mingw, ...)
- `NEST simulator`_ >= 2.0.0


.. _brian: http://briansimulator.org
.. _brian simulator: http://briansimulator.org
.. _NEST: http://www.nest-initiative.org
.. _NEST simulator: http://www.nest-initiative.org


Running a simulation
====================

All the simulation scripts are in the grid_cell_model/simulations directory.
For every separate simulation step, create one subdirectory with a meaningful
name and organize your simulation runs and data into that folder.

There are several types of files:

*default_params.py*
    Contains default parameters for this particular set of simulations. Each
    specific setup within the simulation directory can then override these
    default parameters.

*submit_\*.py*
    A script that submits a particular simulation run(s). It can run a single
    simulation or one can define a parameter sweep over a given range of
    parameters.

    As an example on how to override the default parameters, see for instance
    005_figures/submit_bump_stability.py

    These scripts are **simulator independent**.

*simulation_\*.py*
    A script that performs the actual simulation run. It sets up the model and
    the parameters, runs the simulation and extracts all the necessary output
    data, and usually saves the data into a user-specified file format
    (currently a MAT file).

    These scripts are usually simulator dependent, at least for now. That means
    that a Brian version will not work with a NEST version.

The submitting scripts currently support submitting to a workstation or to Sun
Grid Engine cluster (using qsub). This can be selected inside the submit_*.py
script. Also, the workstation runs can be submitted in a blocking mode (each
successive simulation waits until the previous one has finished) or in a
non-blocking mode (all simulations are run at once concurrently). See the
blocking parameter in these scripts.

