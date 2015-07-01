Supplements for Figure 2
------------------------

Supplement 1 -- Examples of E cell firing fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This figure contains examples of E cell firing fields from data in Section
:ref:`grids_main_3noise`. Once you have finished the full simulations of animal
movement and performed the data analysis and 'aggregation' steps, you can
simply go to ``007_noise/figures/paper`` and run
``suppFigure_grid_examples.py`` in that directory (provided you have saved the
data according to the instructions). Unlike with other figures, the separate
pages of this figure supplement are in the ``panels`` directory and are named
``suppFigure_grid_examples_NN.pdf``, where ``NN`` stands for page number.


Supplement 2 -- Gridness scores of E cells in probabilistic networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Supplement 3 -- Spatial information and sparsity of E and I cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We again assume that the working directory is in ``007_noise/figures/paper``
and that all the steps in Section :ref:`grids_main_3noise`, especially the full
simulations of animal movement, have been completed.  In order to generate the
panels, run::
    
    ./figure_grids.py --sparsity --spatial_info

This will generate PDF files in the ``panels`` sub-directory. You are looking
for ``grids_spatial_info*.pdf`` and ``grids_spatial_sparsity*.pdf``. The
assembled AI figure is in ``ai/figure_grids_info_scores.ai``.


Supplement 4 -- gridness scores of I cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we assume the same conditions hold as in the previous figure supplement.
To generate the panels, run::

    ./figure_grids.py --grids

and you are looking for files titled ``grids_sweeps*I.pdf``, i.e. gridness
scores of I cells. The assembled figure is in ``ai/figure_grids_i_fields.ai``.


Supplement 5 -- Uncorrelated spatial inputs to I cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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



