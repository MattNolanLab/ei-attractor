Figure 6 -- Simulations without theta input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These simulations show seizure-like activity and gridness score in networks
where theta frequency inputs are replaced with a constant input with the same
mean amplitude. To generate this figure it is necessary to run simulations of
the stationary attractors, velocity calibration, and simulations of animal
movement. The procedure is very similar to the one for the previous figures,
except that some parameters need to be changed. Therefore, there is a separate
set of simulation scripts that are pertinent to this figure. These scripts have
the ``_no_theta`` suffix in their file names.

.. _no_theta_bumps:

Simulations of stationary attractors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run the simulations, simply run the following command:

::

    $ ./submit_param_sweep_gamma_no_theta.py -v DEBUG cluster output/no_theta/gamma_bump/ --ntrials=5 --rtLimit="03:00:00"

And afterwards perform the analysis step:

::

    $ ./submit_analysis_EI.py --rtLimit="02:00:00" --shape 31 31 --ignoreErrors cluster output/no_theta/gamma_bump/ bump gamma --ns_all

And the 'aggregation' step:

::

    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing output/no_theta/gamma_bump/ --ns=0pA
    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing output/no_theta/gamma_bump/ --ns=150pA
    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing output/no_theta/gamma_bump/ --ns=300pA


Velocity calibration simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step is basically the same as in Section :ref:`grids_main_3noise`, except
that we now have to run slightly different scripts. Here is the slightly
abridged description of what needs to be run.

::

  $ ./submit_param_sweep_velocity_no_theta.py -v DEBUG --ntrials=10 cluster output/no_theta/velocity --rtLimit="12:00:00"

When this step is complete, you need to run data post-processing:


::

    $ ./submit_analysis_EI.py --rtLimit="03:00:00" -v DEBUG --shape 31 31 --ignoreErrors cluster output/no_theta/velocity/ velocity

This script will analyze the data and save the results into the original files.
Next you would want to run the 'aggregation' step for the generated velocity
data:

::

    $ ./aggregate_velocity.py -v DEBUG --ntrials=10 --shape 31 31 even-spacing output/no_theta/velocity --ns=0pA
    $ ./aggregate_velocity.py -v DEBUG --ntrials=10 --shape 31 31 even-spacing output/no_theta/velocity --ns=150pA 
    $ ./aggregate_velocity.py -v DEBUG --ntrials=10 --shape 31 31 even-spacing output/no_theta/velocity --ns=300pA 

Next, you need to update these generated bump slope data in the git
repository itself. Now change to the ``007_noise/bump_slope_data`` directory
and run:

::

   $ ./update_no_theta.sh

This will extract the bump slope data from the ``reductions.h5`` files
described previously, and put this data into special files, named
``bump_slope_no_theta_XXXpA.h5``, where ```XXX``` will be 0, 150 or 300.
Currently there is no way how to have several versions of this data in one
place, so every time you run this script, the old data will be overwritten. As
the last step, do not forget to commit this new data into the repository.


Simulations of animal movement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are again very similar to the full simulations in Section
:ref:`grids_main_3noise`, but you will use a different script and in this case
we run only 3 trials:

::

   $ ./submit_param_sweep_grids_no_theta.py -v DEBUG --rtLimit="32:00:00" --ntrials=3 cluster output/no_theta/grids

Next, perform the analysis step:

::

   $ ./submit_analysis_EI.py -v DEBUG --rtLimit="02:00:00" --shape 31 31 --ignoreErrors cluster output/no_theta/grids/ grids --ns_all

and run the 'aggregation' script:

::

   $ ./aggregate_grids.py -v DEBUG --ntrials=3 --shape 31 31 even-spacing output/main_network/grids/ --ns=0pA

Also note that ``--ntrials`` has to be explicitly stated on the command
line. The system is not sophisticated enough to be able to determine how
many trials have been run.


Figure generation
^^^^^^^^^^^^^^^^^

After you have completed all the simulations and analysis, you are ready to
generate the figures. Change your working directory to
``grid_cell_model/simulations/007_noise/figures/paper/no_theta``. Now you have
two options:

 * The first one is to simply run ``make``. This will generate all the
   necessary figure panels into the ``panels`` directory. However there will
   be many more panels and figures present than what is in the main Figure 6
   (some of them are in the supplementary materials).

 * Or run ``./figure_seizures.py`` *and* ``./figure_grids.py``, which is a
   subset of scripts that are called by ``Make``.

Now the fully assembled figure is in ``ai/figure_grids_main.py``. If you just
want to inspect the panels separately, then the files you should be looking for
are:

  * ``paper_bumps_popMaxFR_sweep*.pdf``

  * ``grids_examples*.pdf``

  * ``grid_sweeps*.pdf``



