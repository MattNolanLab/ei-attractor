====================
Grid cell repository
====================

This repository holds the source code of the Grid Cell model. It is based
partly on Python (brian, numpy, scipy, matplotlib) and/or C++ (NEST simulator).

There are two versions therefore, that can (but it is not recommended) be
intermixed:

Brian version
    based on the Brian simulator. This is rather obsolete, Brian has become
    painfully slow with this model and as of current state does not allow for
    parallelisation of the simulation.  It is however very flexible to use,
    because model descriptions can be easily manipulated by simply changing the
    differential equations and a few additional settings.
 
NEST version
    Optimized version that implements the grid cell model as a module for the
    NEST simulator. It is about 10-20 times faster than brian version
    (theoretically, depending on processor type, thread support, etc.) but a
    little harder to comprehend when one needs to change single cell
    properties.

While the idea that a user accesses only a common interface is nice, this was
not possible to achieve completely during this project. Therefore currently,
both models, and especially the simulation scripts (simulation_*.py), are
incompatible.


Installation
============

For installation information, see the INSTALL.rst file.


Repository content
==================

The repository contains several folders, however the most important one is
grid_cell_model. It contains the model and all its necessary component. All
other folders are either experiments not very related to the actual model, but
in the future could be useful.

Cellular_GA
    A very simple MPI implementation of a Cellular Genetic algorithm [KU1995]_.
    This version, or the `Evolving objects`_ can be later used for more
    automatic parameter optimization.

Grid_Cells_ModelDB
    A version of grid cells that will be prepared to ModelDB. This will be
    removed in the future.

cuda_test
    A simple test of GPU computing. Nothing significant.

data
    Data files useful in the simulations. Contain rodent tracking data and
    preprocessing scripts.

data_analysis
    Lots of older MATLAB analysis scripts used in my MSc. thesis. No longer
    maintained, as all these analysis scripts have been ported to Python.

graphs
    Graph theory scripts

grid_cell_model
    Grid cell model simulation scripts and source files

ideas
    Ideas for future modeling that are written down

model_fitting
    Scripts for estimating I-V curve data and subsequently fitting integrate
    and fire models onto the I-V curves. This is no longer maintained actively.

mult_bump_spiking_net
    The first, simplified implementation of the spiking bump attractor network,
    based on [BURAK2009]_. This is quite dead (but will be resurrected when I
    am writing up my thesis - soon)

nest_interface
    A proposal. Possible work in progress to design a C++ NEST interface. God
    knows if it is useful anyhow yet.

oscillations
    Some ancient MATLAB scripts dealing with simulation and analysis of
    gamma oscillations. This was a side-project that was later transferred to
    the theta-gamma oscillatory attractor model (grid_cell_model).

ramp_model
    Another side project that dealt with Hugh Pastoll's ramp gamma model. Due
    to lack of time, this was stopped.



-------------------------------------------------------------------------------

.. _Cellular Genetic algorithm: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.3205

.. _Evolving objects: http://eodev.sourceforge.net

.. [KU1995] K. W. C. Ku, M. W. Mak, and W. C. Siu. A cellular genetic algorithm for
   training recurrent neural networks. In Proceedings of the International
   Conference on Neural Networks and Signal Processing, pages 140--143, 1995.
   <http://citeseer.comp.nus.edu.sg/295805.html>

.. [BURAK2009] Burak, Y., & Fiete, I. R. (2009). Accurate path integration in
   continuous attractor network models of grid cells. PLoS Computational Biology,
   5(2), e1000291. doi:10.1371/journal.pcbi.1000291


