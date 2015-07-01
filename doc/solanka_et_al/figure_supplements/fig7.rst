Supplements for Figure 7
------------------------

Supplement 1 -- E cell firing fields in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For this figure supplement you need data from the networks with I-I synapses
(Section :ref:`fig7_II`). In particular, you need to run the *Velocity
calibration simulations* and *Simulations of animal movement*. Once the data
has been generated, change your working directory to
``grid_cell_model/simulations/007_noise/figure/paper/ii_connections`` and run::

    $ ./figure_grids.py --examplesFlag --examples_colorbar

The generated panels are in ``panels/`` and are named ``grids_examples*.pdf``
and ``grids_examples_colorbar.pdf``. The assembled figure is in
``ai/figure_grids_examples.ai``.


Supplement 2 -- Bump attractors in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This figure supplement, again, requires data from networks with I-I synapses
(Section :ref:`fig7_II`). In this case all data are required because the
contour plots are the gridness score data generated from simulations of animal
movement. Once you have the data switch yourself to
``grid_cell_model/simulations/007_noise/figures/paper/ii_connections`` and
run::

    $ ./figure_bumps.py

This will generate the panels into the ``panels/`` directory. The relevant
files are:

    * ``bumps_isBumpSnapshotExamples_*.pdf``

    * ``bumps_mainFig_isBumpFracTotal_sweeps_annotated*.pdf``

    * and a colorbar: ``bumps_examples_colorbar.pdf``.

The assembled figure is in the usual location: ``ai/figure_bumps.ai``.


Supplement 3 -- Average bump drift in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have the data from stationary attractors in I-I networks
(:ref:`fig7_II`), switch yourself to 
``grid_cell_model/simulations/007_noise/figures/paper/ii_connections`` and run
the usual command::

    $ ./figure_seizures.py --theta
    $ ./figure_drifts.py

The generated files in ``panels`` are ``bumps_drift_at_time_sweeps*.pdf``. The
assembled figure is in ``ai/suppFigure_bumps_drift.ai``.


Supplement 4 -- Velocity calibration in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The velocity calibration simulations require data from I-I networks (Section
:ref:`fig7_II`; only the velocity calibration and animal movement simulation
data). Once you have the data, go to 
``grid_cell_model/simulations/007_noise/figures/paper/ii_connections`` and run
the usual command::

    $ ./suppFigure_velocity.py

This will generate the panels:

    * ``suppFigure_velocity_err_sweeps*.pdf``

    * ``velocity_slope_examples_*.pdf``

    * ``velocity_slope_sweeps*.pdf``

The assembled figure is in ``ai/suppFigure_velocity.ai``.


Supplement 5 -- Seizure activity in networks with I-->I synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As in other seizure figures, you need all data from Section :ref:`fig7_II`.
Once you have it, go to the I-I figure directory:
``grid_cell_model/simulations/007_noise/figures/paper/ii_connections`` and run::

    $ ./figure_seizures.py

This will generate the panels into ``panels`` sub-directory:

    * ``bumps_raster*.pdf``

    * ``bumps_rate*.pdf``

    * ``bumps_popMaxFR_sweep*.pdf``

    * ``bumps_seizureProportion_sweep0.pdf``

    * ``maxFR_gridness_scatter_all.pdf``

    * ``PSeizure_gridness_scatter_all.pdf``

The fully assembled figure is in ``ai/figure_seizures.ai``


Supplements 6-9 -- Networks with E-->E synapses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.


Supplement 10 -- Networks with E-->E synapses and unstructured E-I connectivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not described yet. Perhaps at some point in the future.

