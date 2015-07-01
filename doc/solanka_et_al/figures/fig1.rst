.. _fig1:

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
