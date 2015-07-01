.. _fig7_II:

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
