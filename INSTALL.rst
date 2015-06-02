-------------------------------------------------
Installation instructions for the grid cell model
-------------------------------------------------

Please follow all the following instructions carefully. Skip through the relevant
sections only if you know what you are doing.


Download and initialize the repository
======================================

You will need git to clone this repository using git. After cloning, you need
to initialize the submodules that the repository depends on by running ``git
submodule init`` and ``git submodule update``.


Install the project
===================

The model requires some prerequisites, in particular, it requires some
non-python packages and some python packages.

Install non-Python packages
---------------------------

- C/C++ compiler.
- `NEST simulator`_ == 2.2.0 (and ideally not using the later versions, since
  there have been API changes).
- The HDF5_ library. Most Linux/Unix distributions come with HDF5 installed.
  However, it might be worth to check the version and update to the latest one
  (>= 1.8.12). On Windows, you need to help yourself as you can...
- SWIG_ >= 3.0. This is not required by the project, but it is required by the
  gridcells_ package and so is also listed here.

Install Python packages
-----------------------

There are a number of requirements, but once you have installed the `non-Python
packages`_. Simply go to the root of the project and do these two steps (it is
recommended to do this in a `virtual environment`_ to avoid clashes with other
python packages):

1. ``pip install -r requirements_first.txt``
2. ``pip install -r requirements.txt``

Make sure that there are no errors during the process and that all the packages
are installed successfully.


Make sure PyNEST is properly installed
--------------------------------------

Try to change directory to somewhere outside of the NEST_ installation root
directory. Then try to run python and import the ``nest`` package: ``import
nest``. If you see something like: 

::

  >>> import nest
  
                -- N E S T --
  
    Copyright (C) 2004 The NEST Initiative
    Version 2.2.2 Jun  2 2015 09:37:25
  
  This program is provided AS IS and comes with
  NO WARRANTY. See the file LICENSE for details.
  
  Problems or suggestions?
    Website     : http://www.nest-initiative.org
    Mailing list: nest_user@nest-initiative.org
  
  Type 'nest.help()' to find out more about NEST.

Then the installation was successful. If you see an import error:

::

  >>> import nest
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  ImportError: No module named nest

then there is a problem. There are several ways how to troubleshoot this. One
option is to set the ``PYTHONPATH`` environment variable to point to the installed
package. However, if you are using the `virtual environment`_ (highly
recommended anyway), the best solution is to do the following steps:

1. Activate your virtual environment

2. After installing NEST_ itself (``make and make install``), change to
   directory ``pynest`` *inside* the NEST_ source package and run ``pip
   install``. This should install all the necessary Python bindings into your
   virtual environment.  After this, change directory to *outside* the NEST_
   source package and repeat the steps above to test whether PyNEST is present.


Install the grid cell NEST_ module
----------------------------------

It is also necessary to compile and install the gridcells NEST_ module, which
is a set of C++ classes describing the behavior of single neurons within the
simulator kernel. The module is present in ``grid_cell_model/nest/gridcells``.
Please follow the instructions in `Writing an extension module`_ on the NEST_
website. There is no need to set up anything related to SLI (unless you want to
use SLI). After successful installation, do not forget to set the
``LD_LIBRARY_PATH`` environment variable to point to the location of the
installed module (on OSX the ``DYLD_LIBRARY_PATH`` variable; see the
instructions).


.. _HDF5: https://www.hdfgroup.org/HDF5/ 
.. _SWIG: http://www.swig.org
.. _NEST: http://www.nest-simulator.org
.. _NEST simulator: http://www.nest-simulator.org
.. _Writing an extension module: http://nest.github.io/nest-simulator/extension_modules
.. _gridcells: https://github.com/lsolanka/gridcells
.. _virtual environment: http://docs.python-guide.org/en/latest/dev/virtualenvs/
.. _non-python packages: `Install non-Python packages`_



Running a demo simulation
=========================

Change directory to ``grid_cell_model/simulations/simulation_demo`` and run
``./submit_test_EI.py -h``. This should print the help and description of the
simulation parameters. The typical usage of the way simulations are submitted
in this project is ``./submit_test_EI.py -v DEBUG --time=1e3 --ntrials=1
workstation output_dir``.

This will run the simulation and save data into an appropriate directory in
``output_dir``. If you get any errors, you need to go back to the previous
steps and make sure they are all completed successfully. If this step is
successfull you will find some HDF5_ file(s) in the output directory and you
are ready to run your own simulations.
