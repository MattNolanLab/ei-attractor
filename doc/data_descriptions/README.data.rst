Noise promotes independent control of gamma oscillations and grid firing within recurrent attractor networks
============================================================================================================

Directory structure of the data associated with the publication. Each
subsection lists and describes the directory in the **data** folder.


ee_connections/
---------------

Simulations equivalent to the ones in main_network but with direct E-->E
connections added to the network. For descriptions of subdirectories, see
`main_network`_.


ee_connections_ei_flat/
-----------------------

These simulations are with networks that contain E-->E connections.
In this setting the E-->I and I-->E connection weights are constant (as
opposed to `main_network`_, where they are distance dependent) and E-->E
connection weights are distance dependent with the synaptic weight function
being a Gaussian function centered at zero, i.e. cells that are close to each
other on the twisted torus have strongest connections.

These data are from simulations of a stationary bump attractor networks in
which velocity and place cell inputs are switched off

    **g_EE_total_vs_pEE_sigma/gamma_bump/**

        Parameter sweep of a stationary bump attractor in which the total
        amount of E-->E synaptic weight (g_EE_total parameter) and the width of
        the profile (pEE_sigma parameter) are varied. g_AMPA_total = 1020 nS,
        g_GABA_total = 110 nS, i.e. gE = 1, gI ~ 0.108.
    
    **g_EE_total_vs_pEE_sigma_AMPA_3060_GABA_1020/**
    
        The same simulation as the one above, but g_AMPA_total = 3060,
        g_GABA_total=1020, i.e. gE = 3, gI = 1.
    
    **standard_sweep_g_EE_3060_pEE_sigma_0_0833/**
    
        Standard 2D sweep of gE and gI in which the g_EE = 3060 nS and
        pEE_sigma = 0.0833, i.e. these two parameteres are now fixed and gE and
        gI vary.


i_place_cells/
--------------

Simulations in which I cells receive uncorrelated spatial input in the form of
place cells connected to them.


    **10_trials_rate_100_field_std_80/**
     
        Simulations in which a 2D parameter sweep of the random seed
        (master_seed parameter; essentially defining a trial) and the synaptic
        weight between place cells and I cells (ipc_weight parameter) is
        performed.
     
        This simulation data is used in the Figure supplement in
        [SOLANKA2015]_.
     
    **10_trials_rate_200_field_std_40/**
     
        The same simulation as above, but the max. firing rate of place cells
        (those that are connected to I cells) is doubled and the width of their
        fields are halved (field_std).
     
    **grids/**
     
        A small 2D parameter sweep of the number of connections each I cell
        receives (ipc_nconn) and the synaptic weight from the place cells to I
        cells (ipc_weight).
     
    **grids_max_rate_100_field_std_80/**
     
        Another small 2D parameter sweep of the number of connections between
        place cells and I cells (ipc_nconn) and their synaptic weights
        (ipc_weight). This time, the max. firing rate of place cells is 100 Hz
        and the width of the field is 80 cm.
     
    **grids_max_rate_150_field_std_80/**
     
        The same as above, but the max. firing rate of place cells is 150 Hz.
     
    **grids_max_rate_200_field_std_80/**
     
        The same as above, but the max. firing rate of place cells is 200 Hz.


i_surround/
-----------

Simulations of networks that are in the I-surround configuration (cf.
[PASTOLL2013]_). All these networks have the AMPA_gaussian parameter set to 1
to enable the I-surround configuration.

    **e_theta_475/gamma_bump**

        Simulations of a stationary bump attractor pameter sweep with gE, gI =
        {0..6} nS and theta modulated input to E cells (Iext_e_theta) increased
        to 375 pA.
    
    **gamma_bump_Iext_e_theta_vs_Iext_e_const/**

        Parameter sweep of a stationary bump attractor in which the theta input
        (Iext_e_theta) and constant background input (Iext_e_const) to E cells
        is varied. g_AMPA_total = 3060 nA (gE); g_GABA_total = 1020 nS (gI).
        Iext_e_theta = <375, 675> pA in steps of 10 pA; Iext_e_const = <300,
        450> pA in steps of 5 pA.
    
    **gamma_bump_ratio_corrected/**

        Parameter sweep of a stationary bump attractor in which the amounts of
        excitation and inhibition (gE and gI respectively) are corrected for
        the changes in the total amount of input a postsynaptic cell would
        receive. The correction factor is 10.5, i.e. g_AMPA_total (gE) range is
        <0, 64260> and g_GABA_total (gI) = <0, 582>. The size of the sweep is
        retained and is 31 x 31 items. All other parameters should be as in
        `main_network`_.
    
    **Iext_e_const_375/gamma_bump**
        
        Standard parameter sweep of gE and gI and a stationary bump attractor,
        with the constant input to E cells increased to 375 pA. All other
        parameters are as in `main_network`_.
    
    **Iext_e_const_450/**
    
        Standard parameter sweep of gE and gI and a stationary bump attractor,
        with the constant input to E cells increased to 450 pA. All other
        parameters are as in `main_network`_.
    
    **Iext_e_const_vs_uni_GABA/gamma_bump**

        Parameter sweep of a stationary bump attractor in which the constant
        background input to E cells (Iext_e_const) and the fraction of uniform
        inhibitory input (uni_GABA_frac) are varied. gE and gI have some
        non-standard values here: g_AMPA_total (gE) = 1400 nS; g_GABA_total
        (gI) = 2160 nS. Iext_e_const = <300, 450> pA in steps of 15 pA;
        uni_GABA_frac = <0, 0.2> in steps of 0.01.
    
    **original_e_surround/**

        This is not an E-surround configuration. These are parameter sweeps (gE
        and gI are varied) of a stationary bump attractors (gamma_bump), animal
        movement simulations (grids) and velocity valibration (velocity). Here
        we use I-surround configuration in which only the parameter
        AMPA_gaussian = 1. Otherwise all the parameters are as in
        `main_network`_.
    
    **pastoll_et_al/**
        
        Simulations of networks in the I-surround configuration where the
        parameters are equivalent to networks in [PASTOLL2013]_ (except that
        the size of the E population is smaller and set to 34x30 neurons).
    
        **gamma_bump/**

            Standard parameter sweeps (gE and gI) of the stationary attractors.
     
        **grids/**

            Standard parameter sweeps (gE and gI) of animal movements to
            analyze grid firing fields.
     
        **grids_pc_weight_1/**

            Standard parameter sweeps (gE and gI) of animal movements to
            analyze grid firing fields, but with the max. input from place
            cells to E cells (pc_conn_weight parameter) increased 2-fold to 1
            nS. These data use the same velocity calibration inputs as in the
            ``grids`` folder here.
     
        **grids_pc_weight_3/**

            Standard parameter sweeps (gE and gI) of animal movements to
            analyze grid firing fields, but with the max. input from place
            cells to E cells (pc_conn_weight parameter) increased 6-fold to 3
            nS. These data use the same velocity calibration inputs as in the
            ``grids`` folder here.
     
        **velocity/**

            Standard parameter sweeps (gE and gI) of the calibration of the
            netowrk velocity inputs.
    
    **pastoll_et_al_pc_max_rate_vs_weight/**

        **grids/**

            Parameter sweep simulations of the place cell (PC) maximum firing
            rate (pc_max_rate) and the maximum weight of the connections
            (pc_conn_weight; here PCs --> E cells). Here g_AMPA_total (gE) =
            4080 and g_GABA_total (gI) = 1020 nS.  The velocity calibration
            inputs have been determined in the ``pastoll_et_al`` part of the
            data. pc_max_rate = {50, 100} Hz and pc_conn_weight = <0.5, 10> in
            steps of 0.5 nS.
     
        **grids_no_velocity/**

            The same as above, but the velocity inputs are switched off to test
            to what extent the place cell input controls the firing rate of E
            cells.


ii_connections/
---------------

Parameter sweeps of networks in which I-->I connections have been added.
These I-->I connections are constant and do not depend on the distance (on
the twisted torus) between the pre- and post-synaptic cells.


    **g_II_total_sweep/gamma_bump**
        
        A 1D parameter sweep of a stationary bump attractor network in which
        the strength of I-->I connections is varied. The parameter is called
        g_II_total and varies in the range <0, 400> nS. However, this is the
        total input conductance to a post-synaptic cell. The actual value
        depends on the size of the I population as well as the sparsity of the
        connections (g_II_uni_density parameter). The E->I and I->E synaptic
        scaling parameters are set to g_AMPA_total (gE) = 1020 nS; g_GABA_total
        = 3060 nS.


    **g_II_total_sweep_high_gE/gamma_bump**

        The same parameter sweep as the one above, but g_AMPA_total = 3060 nS;
        g_GABA_total = 1020 nS.


    **gE_vs_gI/**
        
        Standard parameter sweeps of gE and gI in networks that contain I-->I
        synapses. gE and gI vary in the range of 0 to 6 nS.

            **gamma_bump/**

                Simulations of stationary bump attractor - used for
                bump/gamma/drift simulations.

            **grids/**

                Full simulations of animal movement - to produce grid firing
                fields.

            **velocity/**

                Calibration of the gain of the velocity input.


.. _main_network:

main_network/
-------------

Simulation data for most of the main figures, i.e. Figures 1 - 5.

    **connections/**
    
        Data for figures of connection weights.
    
    **detailed_noise/**
    
        Data for figures where noise is increased in more fine grained steps.
    
        **gamma_bump/**
      
            Simulations of stationary bump attractor - used for
            bump/gamma/drift simulations.
    
        **grids/**
    
            Full simulations of animal movement - to produce grid firing
            fields.
    
        **velocity/**
    
            Calibration of the gain of the velocity input.
    
    **gamma_bump/**
    
        Simulations of stationary bump attractor - used for bump/gamma/drift
        simulations.
    
    **grids/**
    
        Full simulations of animal movement - to produce grid firing fields.
    
    **grids_no_pc_input/**
    
        Full simulations of animal movement for grid firing fields.  However,
        place cell input is not active here.
    
    **grids_no_velocity/**
    
        Full simulations of animal movement for grid firing fields. However,
        velocity input is not active here.
    
    **velocity/**
    
        Calibration of the gain of the velocity input.


no_theta/
---------

Simulations that are equivalent to the simulations `main_network`_ folder, but
with theta input replaced with a constant input with the same mean amplitude.

    **gamma_bump/**
    
        Simulations of stationary bump attractor - used for bump/gamma/drift
        simulations.
    
    **grids/**
    
        Full simulations of animal movement - to produce grid firing fields.
    
    **single_neuron/**
    
        Simple simulations of a single neuron to produce membrane potential
        examples (Figures describing the model).
    
    **velocity/**
    
        Calibration of the gain of the velocity input.


probabilistic_connections/
--------------------------

Simulations that are equivalent to the simulations in the main_network folder.
The difference here is that instead the connection weights being drawn
according to the synaptic profile function, their weight is constant (but still
scaled according to gE and gI) and the probability of connections between a
pair of neurons scales according to the synaptic profile function.

    **connections/**
    
        Data for figures of connection weights.
    
    **gamma_bump/**
    
        Simulations of stationary bump attractor - used for bump/gamma/drift
        simulations.
    
    **grids/**
    
        Full simulations of animal movement - to produce grid firing fields.
    
    **velocity/**
    
        Calibration of the gain of the velocity input.


References
----------

.. [PASTOLL2013] Pastoll, H., Solanka, L., van Rossum, M.C.W., and Nolan, M.F.
   (2013). Feedback inhibition enables theta-nested gamma oscillations and grid
   firing fields. Neuron 77, 141–154. 

.. [SOLANKA2015] Solanka, L, van Rossum, M.C.W., and Nolan, M.F. (2015). Noise
   promotes independent control of gamma oscillations and grid firing within
   recurrent attractor networks. In Preparation.

