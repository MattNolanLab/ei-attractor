====================
Grid cell repository
====================

This repository holds the source code of the Grid Cell model. It is based
partly on Python (numpy, scipy, matplotlib) and/or C++ (NEST simulator).

Installation
============

For installation information, see the INSTALL.rst file.


Repository content
==================

The repository contains several folders, however the most important one is
*grid_cell_model*. This directory contains the grid cell model and all its
necessary components.

*data/*
    Data files useful in the simulations. Contains rodent tracking data and
    preprocessing scripts.

*doc/*
    Documentation files. Not much there now, however there is ample
    documentation in the form of docstrings.

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
