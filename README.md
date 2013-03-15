# Grid cell repository #

This repository holds the source code of the Grid Cell model. It is based
partly on Python (brian, numpy, scipy, matplotlib) and/or C++ (NEST simulator).

There are two versions therefore, that can (but it is not recommended) be
intermixed:

    - **Brian version**: based on the Brian simulator. This is rather obsolete,
      Brian has become painfully slow with this model and as of current state
      does not allow for parallelisation of the simulation.

      It is however very flexible to use, because model descriptions can be
      easily manipulated by simply changing the differential equations and a
      few additional settings.

    - **NEST version**: Optimized version that implements the grid cell model
      as a module for the NEST simulator. It is about 10-20 times faster than
      brian version (theoretically, depending on processor type, thread
      support, etc.) but a little harder to comprehend when one needs to change
      single cell properties.

While the idea that a user accesses only a common interface is nice, this was
not possible to achieve completely during this project. Therefore currently,
both models, and especially the simulation scripts (simulation_\*.py), are
incompatible.

# Installation #

For installation, see the INSTALL file
