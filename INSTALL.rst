-------------------------------------------------
Installation instructions for the grid cell model
-------------------------------------------------

In fact the model does not need any installation, it is recommended to work
directly from the source directory. However there are a number of
Prerequisites_ that the user must install in order to start simulation.


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
