.. _sol_prerequisites:

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



