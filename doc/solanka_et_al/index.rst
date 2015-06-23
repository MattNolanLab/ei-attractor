.. _solanka_et_al:

=================================================================
Documentation for the Solanka et al., 2015 paper ([SOLANKA2015]_)
=================================================================

Data directory structure in order to regenerate figures
-------------------------------------------------------

**TODO**: Highlight the necessity of keeping the correct directory structure.


General rules of how to run simulations
---------------------------------------

**TODO**:

1. explain why there are several steps necessary.

2. Note that when running the same script on the same directory more than
   once, the simulations will attempt to resume instead of overwriting.


1. Generate data by running the parameter sweep

2. Analyze the data

3. Perform the reduction scripts

4. Run figure generation scripts.


Figure format in the git repository
-----------------------------------

**TODO**


How to reproduce main figures
-----------------------------

Figure 1
~~~~~~~~

There are two data sets necessary to generate this figure. The first one
generates data for the single neuron examples using the
``submit_single_neuron.py`` script:

::

    $ ./submit_single_neuron.py -v DEBUG --ntrials=1 workstation output/main_network

This is a short script which runs only for a few seconds. Next the sample of
connection weights need to be generated, using the
``submit_param_sweep_connections.py`` script:

::

    $ ./submit_param_sweep_connections.py -v DEBUG --ntrials=1 --probabilistic_synapses=0 workstation output/main_network --nCPU=4

This simulation saves connection weight matrices of E->I and I->E connections
for various values of gE and gI. Only a subset is used. the --nCPU parameter
can be changed to speed up the data generation in case more than 4 processors
are available.

The data for figure 1 do not require any analysis or reduction steps. Therefore
to generate the figures, change to
``grid_cell_model/simulations/007_noise/figures/paper`` and run

::

    $ ./figure_model.py

This should create all the necessary panels for Figure 1 in the ``panels``
directory. The associated files are ``network_layers.png``,
``fig_conn_func_E_surr.pdf``, ``figure_connections_examples_*.pdf`` and
``grids_Vm_example_*.pdf``. Alternatively, ``ai/figure_model.ai`` contains the
full figure with these panels.


.. _fig2:

Figure 2
~~~~~~~~

This figure shows how changes in gE and gI influence grid firing fields. When
running everything from scratch, there are 3 important steps in order to
generate the data:

1. Run the velocity calibration simulations (cf. methods section in the
   publication).

2. Generate/update the bump slope data in the repository from step 1.

3. Run the full simulations of animal movement. This step will produce data
   necessary to analyze spatial firing fields.

The data from steps 1. and 2. are already present in the git repository and
take a relatively long time to run. Therefore, if you strictly do not need to
work with this data you can completely skip these steps and only run the
simulations of animal movement.


.. _grids_main_3noise:

Simulations for the three noise levels (0, 150, 300 pA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Velocity calibration simulations**.

   To generate data to calibrate the velocity inputs in the network, run the
   following scripts:

   ::

     $ ./submit_param_sweep_velocity.py -v DEBUG --ntrials=10 cluster output/main_network/velocity --rtLimit="12:00:00"

   Note that you can replace the ``cluster`` parameter with ``workstation`` and
   an appropriate ``--nCPU`` parameter setting if you do not have Sun Grid
   Engine environment which provides the ``qsub`` command. Note, however, that
   these simulations will need to run 961 separate simulation runs (one run for
   each value of gE and gI) and you will therefore need at least 100-1000
   processor ready in order for this step to complete in a reasonable time.

   When this step is complete, you need to run data post-processing, that,
   again, can take up to an hour (or perhaps more) for each of the gE and gI
   values:

   ::

      $ ./submit_analysis_EI.py --rtLimit="03:00:00" -v DEBUG --ns=150 --shape 31 31 --ignoreErrors cluster output/main_network/velocity/ velocity

   This script will analyze the data and save the results into the original
   files.

2. **Update bump slope data**.

   After you have successfully finished step 1., you need to "aggregate" the
   data into a more compact form that will then be used for submission of the
   full simulations in the next step. This is done by the
   ``aggregate_velociy.py`` script:

   ::

      $ ./aggregate_velocity.py -v DEBUG --ns=0pA --ntrials=10 --shape 31 31 even-spacing output/main_network/velocity

   This script extracts all important data from each simulation run and creates
   one file (in the directory that is relevant to each of the noise levels),
   ``reductions.h5``, which can be accessed much faster than the individual
   simulation run files. Note that you will need to repeat the step for each of
   the noise levels that is specified by the ``--ns`` parameter (abbreviation
   for noise_sigma), i.e.  for sigma = 0 pA you use ``--ns=0pA``, for sigma =
   150 pA you would use ``--ns=150pA`` and for sigma = 300 pA you would use
   ``--ns=300pA``.

   Next, you need to update these generated bump slope data in the git
   repository itself. Now change to the ``007_noise/bump_slope_data`` directory
   and run:

   ::

      $ ./update_reductions.sh

   This will extract the bump slope data from the ``reductions.h5`` file
   described previously, and put this data into special files, named
   ``bump_slope_XXXpA.h5``, where ```XXX``` will be 0, 150 or 300. Currently
   there is no way how to have several versions of this data in one place, so
   every time you run this script, the old data will be overwritten. As the
   last step, do not forget to commit this new data into the repository.

3. **Full simulations**.

   After the first steps are complete, you then need to run the simulations of
   an animal moving in an arena. This is accomplished by the
   ``submit_param_sweep_grids.py`` script. Again, this is ideal to run as a
   batch job on a cluster, by issuing the following command:

   ::

      $ ./submit_param_sweep_grids.py -v DEBUG --rtLimit="32:00:00" --ntrials=4 cluster output/main_network/grids

   .. note::

      Note here that the ``--rtLimit`` parameter is quite high. The run time
      for some of the gE and gI parameters can be up to 8h. Some of the
      simulations will not finish at all (as is the case with the velocity
      calibration simulations). In general, for each trial in this
      simulations set you will need at least 8h of run time. If your cluster
      does not allow you to use 32h run time limit, you can perform the
      simulation script 4 times with ``rtLimit="08:00:00" --ntrials=1``
      (waiting for each batch of trials to complete fully). When the same
      script is run with the same output directory, the simulation will try
      to append unfinished trials to the already exisiting ones, instead of
      overwriting the old data.

   As with the velocity calibration simulations, after this step is complete
   you need to run the analysis script:

   ::

      $ ./submit_analysis_EI.py -v DEBUG --rtLimit="01:30:00" --shape 31 31 --ignoreErrors cluster output/main_network/grids/ grids 

   This will perform analysis of firing fields and will save the data for each
   trial into the original files.

   After you are complete with this step, you need to 'aggregate' the data
   again, by running ``aggregate_grids.py`` and using the correct value of the
   ``--ns`` parameter for each noise level:

   ::

      $ ./aggregate_grids.py -v DEBUG --ntrials=4 --shape 31 31 even-spacing output/main_network/grids/ --ns=0pA

   Also note that ``--ntrials`` has to be explicitly stated on the command
   line. The system is not sophisticated enough to be able to determine how
   many trials have been run.


.. _grids_main_detailed_noise:

Simulations for noise levels with finer increase (0 - 300 pA, 10 pA steps)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the panels in Figure 2 contains gridness score of networks as a function
of the noise level. Here the noise level varies from 0 to 300 pA, in 10 pA
steps. The steps to generate the data are similar to the previous section, but
the script have a ``_noise`` suffix in their name.

1. **Velocity calibration simulations**.

   For the velocity calibration you need to run (cf. :ref:`grids_main_3noise`)

   ::

      $ ./submit_param_sweep_velocity_noise.py -v DEBUG --rtLimit="12:00:00" --ntrials=10 --EI_type=EI-1_3 cluster output/main_network/detailed_noise/velocity
      $ ./submit_param_sweep_velocity_noise.py -v DEBUG --rtLimit="12:00:00" --ntrials=10 --EI_type=EI-3_1 cluster output/main_network/detailed_noise/velocity

   The first line is for simulations in which gE = 1 nS and gI = 3 nS, and the
   second line is for simulations simulations in which gE = 3 nS and gI = 1 nS.

   Once the simulations are complete, the next step is to run the data
   analysis:

   ::

      $ ./submit_analysis_detailed_noise.py --rtLimit="03:00:00" -v DEBUG --where=output/main_network/detailed_noise/velocity/ --type=velocity --env=cluster --all-positions --ignoreErrors

   This is only needed to be run once for both of the simulation runs described
   above.

   As in the other simulations, you now need to 'aggregate' some of the data,
   by running:

   ::

      $ ./aggregate_velocity.py -v DEBUG --shape 31 9 detailed-noise output/main_network/detailed_noise/velocity/ --position=EI-1_3
      $ ./aggregate_velocity.py -v DEBUG --shape 31 9 detailed-noise output/main_network/detailed_noise/velocity/ --position=EI-3_1

   **Be very careful** to keep the shape parameter as ``--shape 31 9``,
   otherwise you will not be able to successfully complete the next steps. This
   will produce the ``reductions.h5`` file for each of the directories in
   ``output/main_network/detailed_noise/velocity``.

2. **Update bump slope data**.

   Here you simply change directory to ``007_noise/bump_slope_data`` and run

   ::

      $ ./update_detailed_noise.sh

   Again, this will overwrite the old data and it is also good to commit the
   changes into the repository.

3. **Full simulations**.

   This step generates the data from simulations of animal movement, but in
   this case the noise is varied using much finer steps. You need to run
   separate batches for the different network conditions (gE and gI values):

   ::

      $ ./submit_param_sweep_grids_noise.py -v DEBUG --where=output/main_network/detailed_noise/grids --env workstation --position EI-1_3
      $ ./submit_param_sweep_grids_noise.py -v DEBUG --where=output/main_network/detailed_noise/grids --env workstation --position EI-1_3

   When this is complete, the next step is to run the analysis on these two
   data sets (only the following command is necessary):

   ::

      $ ./submit_analysis_detailed_noise.py --rtLimit="01:30:00" -v DEBUG --where=output/main_network/detailed_noise/grids/ --type=grids --env=cluster --all-positions --ignoreErrors

   And after that 'aggregate' the important data from all the data sets:

   ::

      $ ./aggregate_grids.py -v DEBUG --shape 31 9 detailed-noise output/main_network/detailed_noise/grids/ --position=EI-1_3 --ntrials=1
      $ ./aggregate_grids.py -v DEBUG --shape 31 9 detailed-noise output/main_network/detailed_noise/grids/ --position=EI-3_1 --ntrials=1


Generate the figure
^^^^^^^^^^^^^^^^^^^

After you have successfully completed all the main steps from Sections
:ref:`grids_main_3noise` and :ref:`grids_main_detailed_noise`, you should be
ready to generate all the panels for Figure 2. To do this, change directory to
``grid_cell_model/simulations/007_noise/figures/paper`` and run

::

    $ ./figure_grids.py --grids --examplesFlag --examples_colorbar --detailed_noise --diff_sweep

This will generate PDF files with the ``grids_`` prefix in the ``panels``
directory. The assembled figure is in ``ai/figure_grids.ai``. To properly show
the figure (since the AI file contains only **links** to the figure panels and
these links are absolute) you will need to open it and point the editor to the
correct files that are in *your* ``panels`` directory.


Figures 3, 4 and 5
~~~~~~~~~~~~~~~~~~

All of these figures use data from the common data set which simulates a
stationary bump attractor with velocity and place cell inputs switched off.
Moreover, to generate scatter plots where gridness score appears on the Y axis,
you need to have completed all the steps simulations from section
:ref:`grids_main_3noise` because the generation process requires gridness
scores from this data set as well.

Generate common data of stationary bump attractors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step is more straightforward, because the simulations do not use any
velocity input calibration. Again, change the directory to the root of the
simulation scripts (``grid_cell_model/simulations/007_noise``) and run

::

    $ ./submit_param_sweep_gamma.py -v DEBUG cluster output/main_network/gamma_bump/ --ntrials=5 --rtLimit="03:00:00"

These simulations usually take much shorter and it should also be possible to
run them on a simple workstation in a reasonable time with 32 -- 64 processors.

Once this step is complete, it is necessary to run the analysis script (this
will do the work for all noise levels).

::

    $ ./submit_analysis_EI.py --rtLimit="01:30:00" --shape 31 31 --ignoreErrors cluster output/main_network/gamma_bump/ bump gamma --ns_all

And after this step is done, 'aggregate' the data into a more compact form,
this time using 3 commands:

::

    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing output/main_network/gamma_bump/ --ns=0pA
    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing output/main_network/gamma_bump/ --ns=150pA
    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing output/main_network/gamma_bump/ --ns=300pA


Simulations with finer noise level increase (0 -- 300 pA, 10 pA steps)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we need to do a similar procedure as in the description of :ref:`fig2`,
except that the simulations will take much shorter time. First, run the
simulation scripts, for both values of gE and gI (i.e. ``--position=EI-1_3``
stands for gE = 1 nS; gI = 3 nS):

::

    $ ./submit_param_sweep_gamma_noise.py -v DEBUG --ntrials=5 --rtLimit="03:00:00" cluster output/main_network/detailed_noise/gamma_bump --position=EI-1_3
    $ ./submit_param_sweep_gamma_noise.py -v DEBUG --ntrials=5 --rtLimit="03:00:00" cluster output/main_network/detailed_noise/gamma_bump --position=EI-3_1

Now run the analysis scripts for both cases

::

    $ ./submit_analysis_detailed_noise.py --where=output/main_network/detailed_noise/gamma_bump/ --type gamma --env cluster --all-positions --ignoreErrors --rtLimit="01:30:00"
    $ ./submit_analysis_detailed_noise.py --where=output/main_network/detailed_noise/gamma_bump/ --type bump --env cluster --all-positions --ignoreErrors --rtLimit="01:30:00"

And when finished, 'aggregate' the data into a more compact form:

::

    $ ./aggregate_bumps.py -v DEBUG --shape 31 9 detailed-noise output/main_network/detailed_noise/gamma_bump/ --position=EI-1_3 --ntrials=5 --positions --AC 
    $ ./aggregate_bumps.py -v DEBUG --shape 31 9 detailed-noise output/main_network/detailed_noise/gamma_bump/ --position=EI-3_1 --ntrials=5 --positions --AC 

Again, make sure that the ``--shape 31 9`` parameter is entered exactly as it
is here, since not doing so will produce incorrect data and the figure
generation steps will then fail.

Once this is done, you are ready to generate the figures.


Generate the figures
^^^^^^^^^^^^^^^^^^^^

To generate the figures, change your working directory to
``grid_cell_model/simulations/007_noise/figures/paper`` and follow the next
steps.

 1. **Figure 3 - gamma activity** -- run: ``./figure_gamma.py``. This will
    generate figure panels with the ``gamma_`` prefix [#gamma_fnames]_. The
    fully assembled figure is then in ``ai/figure_gamma.ai``. As with other
    AI files you will need to set the links to the figure panels properly when
    you first open the file (after you have run the figure generation script).

 2. **Figure 4 - bump attractor activity** -- run ``./figure_bumps.py``. This
    will generate figures with the ``bumps_`` prefix [#bumps_fnames]_. The
    fully assembled figure is then in ``ai/figure_bumps.ai``.

 3. **Figure 5 - seizure-like activity** -- run ``./figure_seizures.py``. This
    will generate figures with various (and perhaps a little confusing)
    prefixes in the ``panels`` directory. The file names to look for are the
    following

      * ``bumps_raster*.pdf``

      * ``bumps_rate*.pdf``

      * ``bumps_popMaxFR_sweep*.pdf``

      * ``bumps_seizureProportion_sweep0.pdf``

      * ``maxFR_gridness_scatter_all.pdf``

      * ``PSeizure_gridness_scatter_all.pdf``

    Again, the fully assembled figure is in ``ai/figure_seizures.ai``.

    
.. [#gamma_fnames] Some file names will have a ``gammaFreq_`` prefix. This
                   script also generates panels for Figure 3 -- figure
                   supplement 4. These will have a ``gridness_filt_`` prefix in
                   their file name.

.. [#bumps_fnames] The file names for some of the scatter plots will be
                   ``gamma_scatter_gamma_pbumps_all_exp.pdf`` and
                   ``gamma_scatter_gamma_pbumps_all.pdf``.


Figure 6
~~~~~~~~

Figure 7
~~~~~~~~


How to reproduce Figure supplements
-----------------------------------


References
----------

.. [SOLANKA2015] Solanka, L, van Rossum, M.C.W., and Nolan, M.F. (2015). Noise
   promotes independent control of gamma oscillations and grid firing within
   recurrent attractor networks. In Preparation.
