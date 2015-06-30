===================
E-I attractor model
===================

This repository holds the source code for the Excitatory-Inhibitory attractor
model which was used in [SOLANKA2015]_. It is based partly on Python (numpy,
scipy, matplotlib) and/or C++ (NEST simulator).


Installation
============

For installation information, see the ``INSTALL.rst`` file. This file also
contains information on how to build the documentation only.


Repository content
==================

The repository contains several folders, however the most important one is
*grid_cell_model*. This directory contains the grid cell model and all its
necessary components.

*data/*
    Data files useful in the simulations. Contains rodent tracking data and
    preprocessing scripts.

*doc/*
    Documentation files. Please read ``INSTALL.rst`` on how to build the
    documentation. This does not require installing the whole project.

*grid_cell_model/*
    Grid cell model simulation scripts and source files. This is the main
    directory. If you just want to work with the grid cell model, you can
    simply ignore all the other directories.

*manuscripts/*
    Various manuscripts for publication purposes.

*noisefigs/*
    Python package that contains figure generation code for
    ``grid_cell_model/simulations/007_noise``.

*simtools*
    The simtools_ package (as a git submodule).


.. _simtools: https://github.com/lsolanka/simtools

License
=======

``ei-attractor`` is distributed under the GPL license. See ``LICENSE.txt`` in
the root of the repository.

References
==========

.. [SOLANKA2015] Solanka, L, van Rossum, M.C.W., and Nolan, M.F. (2015). Noise
   promotes independent control of gamma oscillations and grid firing within
   recurrent attractor networks. In Preparation.

