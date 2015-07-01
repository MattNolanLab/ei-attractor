.. _solanka_et_al:

=======================================================
Documentation for Solanka et al., 2015 ([SOLANKA2015]_)
=======================================================

This documentation provides the necessary information to reproduce simulation
data for all main figures, the figures themselves and some figure supplements
in [SOLANKA2015]_. All the steps described here are necessary and have to
complete successfully in order to reproduce the results obtained in the
publication.


Prerequisite information to run simulations and generate figures
----------------------------------------------------------------

Before running the simulations and generating figures, you need to know the
structure of the repository and associated important initialization steps that
you need to perform. Only skip these if you know exactly what you are doing.

Data directory structure necessary to regenerate the figures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All the scripts necessary to run simulations and generate figures are in the
``grid_cell_model/simulations/007_noise`` directory. The figure generation
scripts are in the ``figures/paper`` sub-directory. These scripts use symbolic
links (symlinks) in order to access one and the same data repository. In
addition, the figure generation scripts require that the data repository has a
specific structure and thus when running the actual simulations you need to
honor this structure, or otherwise you will have to change some low level
configuration file settings. Fortunately, the only thing that is needed is to
create a symlink in the root of the project source code::

    $ ln -s <path_to_data_repository_root> noise_grids_gamma_data

When you then run the simulation scripts you can access your central data
repository in ``007_noise`` via the ``output`` symlink. For instance, you would
use the second positional argument to ``submit_param_sweep_grids.py`` to save
the data into ``output/main_network/grids`` directory (see Section
:ref:`sim_scripts`), since output will then be linked to
``<path_to_data_repository_root>``. The figure generation scripts are also
symlinked to this output directory so if you honor the default directory
structure (shown in the tree below), everything should run smoothly.

As already mentioned, in order for the figure generation scripts to find the
data it is necessary that all the important simulation commands generate their
outputs into correct sub-directories of the data store you are using. The
following is an example of the directory structure that was used to produce all
the figures and supplementary materials for [SOLANKA2015]_ (there is in fact
much more data than necessary) [#1]_. This structure is also followed in the
samples of commands which you will be using to generate data in Section
:ref:`howto_main_figs`.

::

    noise_grids_gamma
    +-- data
        +-- ee_connections/
            +-- gamma_bump/
                +-- 0pA/
                    +-- iterparams.h5
                    +-- job00000_output.h5
                    +-- job00001_output.h5
                    +-- ....
                    +-- jobN_output.h5
                    +-- reductions.h5
                +-- 150pA/
                +-- 300pA/
            +-- grids/
            +-- velocity/
        +-- ee_connections_ei_flat/
            +-- g_EE_total_vs_pEE_sigma/
            +-- g_EE_total_vs_pEE_sigma_AMPA_3060_GABA_1020/
            +-- standard_sweep_g_EE_3060_pEE_sigma_0_0833/
        +-- i_place_cells/
            +-- 10_trials_rate_100_field_std_80/
        +-- i_surround/
            +-- pastoll_et_al/
                +-- gamma_bump/
                +-- grids/
                +-- grids_pc_weight_1/
                +-- grids_pc_weight_3/
                +-- velocity/
        +-- ii_connections/
            +-- gE_vs_gI/
                +-- gamma_bump/
                +-- grids/
                +-- velocity/
        +-- main_network/
            +-- connections/
            +-- const_position/
            +-- detailed_noise/
            +-- gamma_bump/
            +-- grids/
            +-- grids_no_pc_input/
            +-- grids_no_velocity/
            +-- single_neuron/
            +-- velocity/
        +-- no_theta/
            +-- gamma_bump/
            +-- grids/
            +-- velocity/
        +-- probabilistic_connections/
            +-- connections/
            +-- gamma_bump/
            +-- grids/
            +-- velocity/

A few notes here. The sub-directory ``data/ee_connections/gamma_bump`` contains
three further sub-directories -- ``0pA``, ``150pA`` and ``300pA``. Each of
these is a place holder for all the simulation data for each simulated noise
level in the publication. In the tree list only ``ee_connections/gamma_bump``
shows these three noise levels, however almost all other sub-directories will
contain them after the simulations are complete. Also, it is in general not a
good idea to temper with the contents of the ``0/100/300pA`` directories,
however it is safe to move them as a whole, e.g. as part of the whole
``gamma_bump``/``grids``/``velocity`` set.

Finally, an important and useful thing to notice is that the figure generation
scripts mentioned in the next sections use configuration files which can be
used to change where the figure scripts look for their particular data sets, as
well as change the visual appearance of figure panels (e.g. figure sizes,
annotations, X and Y labels, etc.). These configuration files are present in
different places:

1. There is the **default** configuration file in
   ``noisefigs/noisefigs/default_config.py``. This is normally imported in
   the beginning and has to be overriden, otherwise configuration settings
   from this file will be used

2. Each directory in ``grid_cell_model/simulations/007_noise/figures/paper``
   and associated sub-directories will usually contain its own ``config.py``
   file which is used to override the settings in ``default_config.py``. For
   instance, ``figures/paper/ii_connections/config.py`` contains configuration
   values that are specific for generating figures from networks with I-->I
   synapses.

.. [#1] This data set should already be publicly available at the time you are
        reading this. If not, please contact the corresponding author of
        [SOLANKA2015]_.


.. _sim_scripts:

Simulation scripts
~~~~~~~~~~~~~~~~~~

.. highlight:: console

All simulations scripts for this paper are present in the
``grid_cell_model/simulations/007_noise`` directory. They have the ``submit_``
prefix in their file name. These scripts essentially run on top of an
abstraction layer that allows the user to run the same set of simulations
either on a cluster system supporting the ``qsub`` command, or on a standard
multi-core workstation. Each script is an executable that accepts several
parameters. For example, to reproduce some of the data for Figure 2, we would
use the ``submit_param_sweep_grids.py`` script (see below). Every script print
a help text when called with the ``-h`` parameter (here only a part of the help
shown)::

    $ ./submit_param_sweep_grids.py -h
    usage: submit_param_sweep_grids.py [-h] [--all]
                                       [-v {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                                       [--time TIME] --ntrials NTRIALS
                                       [--rtLimit RTLIMIT] [--printout {0,1}]
                                       [--nCPU NCPU] [--dry_run]
                                       [--ns {0,150,300}] [--row ROW] [--col COL]
                                       {workstation,cluster} where
    
    positional arguments:
      {workstation,cluster}
                            How to run the simulations. If `workstation`, run
                            locally on the current machine. If 'cluster', run on
                            the SGE cluster using the qsub command.
      where                 Root directory of output data. This will be passed on
                            to the simulation script.
    
    optional arguments:
      --rtLimit RTLIMIT     Run time limit. Applicable only when submitting the
                            simulation on a cluster using qsub.
      --nCPU NCPU           Number of processors when running on a workstation.
                            This can be used to run several simulations in
                            parallel.
      --dry_run             Do no run anything nor save any meta-data
      --ns {0,150,300}

There are two important **positional arguments**. The first one selects the
environment type (``workstation`` or ``cluster``) and the second one specifies
the output directory of the whole simulation batch.

Using a multi-core workstation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a simpler method, because it does not require any extra environment
settings. Simply set the environment (first) positional argument to
``workstation`` and specify the number of jobs you want to run in parallel with
the ``--nCPU`` parameter. This should in general match the number of cores the
machine has. Note that some of the simulations require 100--1000 cores to
complete in a reasonable time, while some simulations are shorter and might as
well run in a few days when using ~30--50 cores.

.. _sge_info:

Using a Sun Grid Engine (SGE) cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running on an SGE cluster can cut down the simulation time several fold by
submitting several instances of the simulation with different parameter
settings in parallel. However, it takes slightly more work to set up the
environment.


Assuming the current working directory is
``grid_cell_model/simulations/007_noise``, the very first step is to provide
the right settings in the ``cluster_submit.sh`` script. The settings and
environment variables will depend on where you installed the project and what
version of Python you are using. One way to start is to consult a sample
version of the ``cluster_submit.sh`` script in
``grid_cell_model/simulations/simulation_demo``, test the correct values by
running a few short demo simulations, and then update the script in the
``007_noise`` directory.

When the ``cluster_submit.sh`` script is correctly updated, you simply run all
the simulation commands with the first positional argument set to ``cluster``
instead of workstation. Not that all the descriptions of commands for
simulation submission assume that you have the cluster environment set up and
therefore the environment positional argument is set to ``cluster``. To run the
simulations on a multi-core workstation, simply replace ``cluster`` with
``workstation`` and add an appropriate ``--nCPU=XX`` parameter, where ``XX`` is
the number of cores you want to utilize for the simulation run.


General rules about how to run simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The process to create data ready for the figure generation process consists of
the following steps:

1. Generate data by running the parameter sweep, using one of the
   ``submit_param_sweep*.py`` scripts.

   **Important notes**:
     * For some values of gE and gI the simulations might take much longer
       because of increased spiking activity of E or I cells. If you are
       running these simulations on a cluster with a fixed run time limit, they
       might be stopped if not finished in the pre-specified amount of time. If
       this happens, it is possible to run the whole batch of simulations again
       (with **exactly** the same parameter set as before) and the simulation
       scripts will attempt to re-run the missing simulation trials and append
       them to the ones completed previously.

     * Some simulations need more memory than others and if you also have a
       memory limit (usually on the cluster), they might not complete either.
       In this case re-running the simulation batch will not help, but adding a
       parameter that asks for more memory on a worker node might help (see
       Section :ref:`sge_info` or run ``man qsub`` or alternatively consult
       your local cluster administrator on how to do this).

2. Run data analysis on the generated data, using the submit
   ``submit_analysis_EI.py`` script for most of the simulations, or the
   ``submit_analysis_detailed_noise.py`` script for simulations of detailed
   noise levels (Figure 2H and 3H).

3. Perform the reduction scripts, using either ``aggregate_bumps.py``,
   ``aggregate_grids.py`` or ``aggregate_velocity.py``, depending on the
   simulation type.

4. Run the figure generation scripts.

    

Figure format in the git repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The figure generation scripts in fact only generate the figure panels itself.
The workflow I have adopted here is to use Adobe Illustrator (AI; ``*.ai``)
files and create links from withing th AI files to the panels. There are
advantages to this approach, since one can then easily update separate panels
in the assembled figure with only updating the PDF files on disk. However,
there are also caveats with this approach. The main one [#2]_ is that
unfortunately the paths to the linked files saved in the AI file are absolute
(yes!). Therefore, when you make a different copy of the repository and
generate figures, when you open these AI files for the first time, you will be
asked to provide the paths to the files. This can be done in a batch for each
file, for instance by pointing the editor to the first file, which will
automatically be set as the base directory of other files. I recommend to do
this re-linking process separately for each file, otherwise the editor will get
confused and not link the files properly. You can always cross-check with the
publication.

.. [#2] The other one is having to commit binary files into the git repository.
        Ugly, ugly!


.. _howto_main_figs:

How to reproduce the figures
----------------------------

Figure 1 -- model description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Figure 2 -- simulations of grid fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

      $ ./submit_analysis_EI.py --rtLimit="03:00:00" -v DEBUG --shape 31 31 --ignoreErrors cluster output/main_network/velocity/ velocity --ns_all

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

      $ ./submit_analysis_EI.py -v DEBUG --rtLimit="01:30:00" --shape 31 31 --ignoreErrors cluster output/main_network/grids/ grids --ns_all

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


Figures 3, 4 and 5 -- gamma activity, bump attractors and seizure-like activity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of these figures use data from the common data set which simulates a
stationary bump attractor with velocity and place cell inputs switched off.
Moreover, to generate scatter plots where gridness score appears on the Y axis,
you need to have completed all the steps simulations from section
:ref:`grids_main_3noise` because the generation process requires gridness
scores from this data set as well.

.. _bumps_common_3noise:

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



Figure 7 -- I --> I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The protocol to generate the data for simulations of networks with I --> I
synapses is basically similar to the one described in Section
:ref:`grids_main_3noise`, excepts that the scripts use a slightly updated form
of parameters that need to be supplied on command line. To generate all the
parts of the figure, follow all the next sub-sections.

Simulations of stationary attractors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we again have our three usual steps of data generation and analysis. First
run

::

    $ ./submit_param_sweep_gamma_ii_connections.py -v DEBUG --time=10e3 --ntrials=5 \
            --rtLimit="03:00:00" cluster output/ii_connections/gE_vs_gI/gamma_bump  \
            g_AMPA_total g_GABA_total --range1 0 6120 204 --range2 0 6120 204

After this step is finished, analyze the data by running: 

::

    $ ./submit_analysis_EI.py --rtLimit="01:30:00" --shape 31 31 --ignoreErrors \
            cluster output/ii_connections/gE_vs_gI/gamma_bump bump gamma --ns_all

And afterwards run the 'aggregation' step:

::

    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing \
            output/ii_connections/gE_vs_gI/gamma_bump --ns=0pA

    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing \
            output/ii_connections/gE_vs_gI/gamma_bump --ns=150pA

    $ ./aggregate_bumps.py -v DEBUG --ntrials=5 --shape 31 31 even-spacing \
            output/ii_connections/gE_vs_gI/gamma_bump --ns=300pA
 


Velocity calibration simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step is basically the same as in Section :ref:`grids_main_3noise`, except
that we now have to run slightly different scripts, with the
``_ii_connections.py`` suffix. Here is a slightly abridged description of what
needs to be run.

::

    $ ./submit_param_sweep_velocity_ii_connections.py -v DEBUG --time=10e3 \
            --ntrials=10 --rtLimit="12:00:00" \
            cluster output/ii_connections/gE_vs_gI/velocity \
            g_AMPA_total g_GABA_total --range1 0 6120 204 --range2 0 6120 204

When this step is complete, you need to run data post-processing:

::

    $ ./submit_analysis_EI.py --rtLimit="03:00:00" -v DEBUG --shape 31 31 --ignoreErrors \
            cluster output/ii_connections/gE_vs_gI/velocity velocity

This script will analyze the data and save the results into the original files.
Next you would want to run the 'aggregation' step for the generated velocity
data:

::

    $ ./aggregate_velocity.py -v DEBUG --ntrials=10 --shape 31 31 even-spacing \
            output/ii_connections/gE_vs_gI/velocity --ns=0pA

    $ ./aggregate_velocity.py -v DEBUG --ntrials=10 --shape 31 31 even-spacing \
            output/ii_connections/gE_vs_gI/velocity --ns=150pA 

    $ ./aggregate_velocity.py -v DEBUG --ntrials=10 --shape 31 31 even-spacing \
            output/ii_connections/gE_vs_gI/velocity --ns=300pA 

Next, you need to update these generated bump slope data in the git
repository itself. Now change to the ``007_noise/bump_slope_data`` directory
and run:

::

   $ ./update_ii_connections.sh

This will extract the bump slope data from the ``reductions.h5`` files
described previously, and put this data into special files, named
``bump_slope_ii_connections_XXXpA.h5``, where ```XXX``` will be 0, 150 or 300.
Currently there is no way how to have several versions of this data in one
place, so every time you run this script, the old data will be overwritten. As
the last step, do not forget to commit this new data into the repository.


Simulations of animal movement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are again very similar to the full simulations in Section
:ref:`grids_main_3noise`, but you will use a different script and in this case
we in fact run only 1 trial:

::

    $ ./submit_param_sweep_grids_ii_connections.py -v DEBUG --time=600e3 \
            --ntrials=1 --rtLimit="08:00:00" \
            cluster output/ii_connections/gE_vs_gI/grids \
            g_AMPA_total g_GABA_total --range1 0 6120 204 --range2 0 6120 204

Next, perform the analysis step:

::

   $ ./submit_analysis_EI.py -v DEBUG --rtLimit="02:00:00" --shape 31 31 --ignoreErrors \
        cluster output/ii_connections/gE_vs_gI/grids grids --ns_all

and run the 'aggregation' script:

::

   $ ./aggregate_grids.py -v DEBUG --ntrials=1 --shape 31 31 \
        even-spacing output/ii_connections/gE_vs_gI/grids --ns=0pA

   $ ./aggregate_grids.py -v DEBUG --ntrials=1 --shape 31 31 \
        even-spacing output/ii_connections/gE_vs_gI/grids --ns=150pA

   $ ./aggregate_grids.py -v DEBUG --ntrials=1 --shape 31 31 \
        even-spacing output/ii_connections/gE_vs_gI/grids --ns=300pA

Also note that ``--ntrials`` has to be explicitly stated on the command
line. The system is not sophisticated enough to be able to determine how
many trials have been run.


Figure generation
^^^^^^^^^^^^^^^^^

Since the main figure contains panels from different data sets, it is necessary
to run more than one figure generation script. Change your working directory to
``grid_cell_model/simulations/007_noise/figures/paper/ii_connections``. Now run
exactly these scripts, which will generate files into the ``panels``
sub-directory:

  1. ``./figure_grids.py``. This will generate files with a ``grids_`` prefix.

  2. ``./figure_gamma.py``. This will generate files with a ``gamma_`` prefix.

The fully assembled figure is now in ``ai/figure_gamma_grids_mainfig.ai``.
Again, you need to set up links to files in the ``panels`` sub-directory right
after you have opened the file for the first time.

.. note::

    In this directory there are many more files than necessary. If you have the
    original (and hopefully now already published) data set, you could run
    ``make``, which will recreate the full set of figures relevant to I-I
    synapses. These figures are not published in [SOLANKA2015]_, but might
    nevertheless be useful.


How to reproduce Figure supplements
-----------------------------------

There are very similar scripts to run simulations and generate figures for the
supplements as well. Some of them use the data already generated by the scripts
for the main figures, some of them need the whole cycles of velocity
calibration / grid field simulations. You can find the simulation scripts in
``grid_cell_model/simulations/007_noise`` and the figure generation sripts in
the ``figures/paper`` sub-directory.


Figure 1 -- figure supplement 1 -- Connection weights for scaled and probabilistic networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Figure 2 -- figure supplement 1 -- examples of E cell firing fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This figure contains examples of E cell firing fields from data in Section
:ref:`grids_main_3noise`. Once you have finished the full simulations of animal
movement and performed the data analysis and 'aggregation' steps, you can
simply go to ``007_noise/figures/paper`` and run
``suppFigure_grid_examples.py`` in that directory (provided you have saved the
data according to the instructions). Unlike with other figures, the separate
pages of this figure supplement are in the ``panels`` directory and are named
``suppFigure_grid_examples_NN.pdf``, where ``NN`` stands for page number.


Figure 2 -- figure supplement 2 -- Gridness scores of E cells in probabilistic networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Figure 2 -- figure supplement 3 -- Spatial information and sparsity of E and I cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We again assume that the working directory is in ``007_noise/figures/paper``
and that all the steps in Section :ref:`grids_main_3noise`, especially the full
simulations of animal movement, have been completed.  In order to generate the
panels, run::
    
    ./figure_grids.py --sparsity --spatial_info

This will generate PDF files in the ``panels`` sub-directory. You are looking
for ``grids_spatial_info*.pdf`` and ``grids_spatial_sparsity*.pdf``. The
assembled AI figure is in ``ai/figure_grids_info_scores.ai``.


Figure 2 -- figure supplement 4 -- gridness scores of I cells.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we assume the same conditions hold as in the previous figure supplement.
To generate the panels, run::

    ./figure_grids.py --grids

and you are looking for files titled ``grids_sweeps*I.pdf``, i.e. gridness
scores of I cells. The assembled figure is in ``ai/figure_grids_i_fields.ai``.


Figure 2 -- figure supplement 5 -- Uncorrelated spatial inputs to I cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is quite tricky, because the simulations and data analysis are
non-standard.

Prerequisites
^^^^^^^^^^^^^

These simulations require that at least the velocity calibration step from
Section :ref:`grids_main_3noise` is fully finished, including updating the bump
slope data. These bump velocity slope data are required here, because the
network must be properly calibrated before connecting place cells to I cells.
If you have not done these simulations yet, you need to go back to Section
:ref:`grids_main_3noise` and follow the *Velocity calibration simulations* and
*Update bump slope data* steps.


Generate simulation data
^^^^^^^^^^^^^^^^^^^^^^^^

Next, ensure your working directory is
``grid_cell_model/simulations/007_noise`` and then run the grid field
simulation with the ``_ipc`` suffix::

    ./submit_param_sweep_grids_ipc.py -v DEBUG --time=600e3 --ntrials=1 \
            --g_AMPA_total=3060 --g_AMPA_row=15 --g_GABA_total=1020 --g_GABA_col=5 \
            cluster output/i_place_cells/10_trials_rate_100_field_std_80 \
            ipc_weight master_seed \
            --range2 123456 123546 10 --range1 0 6 0.25 --ns=150 \
            --ipc_field_std=80 --ipc_max_rate=100 --ipc_nconn=3 --nrec_spikes_i=510 \
            --rtLimit="08:00:00"

This command will run the full animal movement simulations with various values
of connection weights from place cells to I cells and with various values of
random seeds (random seed here essentially stands for the trial number). This
actually produces more data than needed, since only simulations where
``ipc_weight=4`` are used for plotting the data.

Run analysis
^^^^^^^^^^^^

To analyze the data, run::

    ./submit_analysis_EI.py --shape 25 10 --ns=150 --ignoreErrors --rtLimit="10:00:00" \
            cluster output/i_place_cells/10_trials_rate_100_field_std_80 grids-ipc

This analysis can take several hours for a single job run (i.e. one value of
``ipc_weight`` and ``master_seed``) and **it is important** that the analysis
code **is not interrupted** in the process, because this particular analysis cannot
recover from crashes. Therefore, make sure you have enough of memory and
computing run time limit to do this step.

Once this is finished, there is no need to run the 'aggregation' step, since
during the figure generation process the data will be accessed directly.

.. note::

    Because the analysis script chooses 100 random neurons from each neural
    populations and this particular analysis code does not have a means to
    control random seeds, if you run a new batch of simulations you will
    probably get a different set of 100 neurons that will have been analysed.


Generate the figure
^^^^^^^^^^^^^^^^^^^

Once the data is ready, change your working directory to
``007_noise/figures/paper/i_place_cells`` and run::

    ./figure_grids_trials.py

It might take a while for the script to crunch the data and produce output, but
it should run without errors. This script actually produces a PDF file for each
of the analysed neurons and put the files into
``panels_weight_sparsity_trials``. The histogram plots are the following ones:

    * ``histogram_gridness_150.pdf``

    * ``histogram_info_150.pdf``

    * ``histogram_spasity_150.pdf``

And grid field examples named ``grid_example*.pdf``. As noted earlier, because
of an issue with random seeds here, if you re-generate the data from scratch,
the example firing fields will look different, while the histogram plots should
look very similar. An assembled figure is then in
``ai/ipc_examples_and_histograms.ai``.


Figure 3 -- figure supplements 1-4 -- Figure supplements for gamma activity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Figure 4 -- figure supplement 1 -- Average bump drift
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The bump drift figures are generated from networks that simulate only the
stationary bump attractors. To produce the figures, you need to complete
Section :ref:`bumps_common_3noise` (including the 'aggregation' step). You also
need the gridness score data (for 3 noise levels), generation of which is
described in Section :ref:`grids_main_3noise`.

Once this is done, change your working directory to ``007_noise/figures/paper``
and run::

    $ ./figure_seizures.py --theta
    $ ./figure_drifts.py --bumpDriftSweep


Figure 5 -- figure supplement 1 -- Raster plots of network activity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 6 -- figure supplement 1 -- Gamma activity in networks without theta frequency inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 6 -- figure supplement 2 -- Difference in gridness score in networks without theta frequency inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Figure 6 -- figure supplement 3 -- Mean firing rate of E cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Figure 6 -- figure supplement 4 -- Velocity gain calibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 6 -- figure supplement 5 -- Effectivity of place cell inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.



Figure 7 -- figure supplement 1 -- E cell firing fields in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 7 -- figure supplement 2 -- Bump attractors in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 7 -- figure supplement 3 -- Average bump drift in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 7 -- figure supplement 4 -- Velocity calibration in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 7 -- figure supplement 5 -- Seizure activity in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Figure 7 -- figure supplements 6-9 -- Networks with E-->E synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Figure 7 -- figure supplement 10 -- Networks with E-->E synapses and unstructured E-I connectivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.
