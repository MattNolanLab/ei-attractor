#!/bin/sh

##############################################################################
# Run a simple test script that assembles the network of EI neurons and prints
# basic results of the simulation
##############################################################################

export PYTHONPATH="$BASE"
echo $PYTHONPATH

Ne=400
Ni=100

AMPA_density="0.4"
GABA_density="0.4"

Iext_e="900e-12"
Iext_i="250e-12"

theta_freq="8"

taum_e="9.3e-3"
taum_e_spread="3.1e-3"
EL_e="-68.5e-3"
EL_e_spread="2.0e-3"
Vt_e="-45e-3"
Vr_e=$EL_e
Rm_e="44e6"
ad_tau_e_mean="40e-3"
ad_tau_e_std="5e-3"
ad_e_g_inc="1.136e-8"
deltaT_e="0.4e-3"

taum_i="5e-3"
taum_i_spread="4e-3"
EL_i="-60e-3"
EL_i_spread="20e-3"
Vt_i="-45e-3"
Vr_i="$EL_i"
Rm_i="44e6"
ad_tau_i_mean="7.5e-3"
ad_tau_i_std="0.5e-3"  # Unused in the simulation for now
ad_i_g_inc="2.273e-8"
deltaT_i="0.4e-3"

tau_AMPA="1e-3"
g_AMPA_total="3.5e-8"
g_AMPA_std="600e-12"
tau_GABA_rise="1e-3"
tau_GABA_fall="5e-3"
g_GABA_total="4.00e-8"

Vrev_AMPA="0e-3"
Vrev_GABA="-75e-3"

noise_sigma="2e-3"
sigma_init_cond="10e-3"

refrac_abs="0.1e-3"

time=5
sim_dt="0.1e-3"
spike_detect_th="40e-3"
Vclamp="-50e-3"

ntrials=1

output_dir="output_local"
update_interval=10
job_num=1



python2.6 -i test_EI_network.py \
--Ne $Ne \
--Ni $Ni \
--AMPA_density $AMPA_density \
--GABA_density $GABA_density \
--Iext_e $Iext_e \
--Iext_i $Iext_i \
--theta_freq $theta_freq \
--taum_e $taum_e \
--taum_e_spread $taum_e_spread \
--EL_e $EL_e \
--EL_e_spread $EL_e_spread \
--Vt_e $Vt_e \
--Vr_e $Vr_e \
--Rm_e $Rm_e \
--ad_tau_e_mean $ad_tau_e_mean \
--ad_tau_e_std $ad_tau_e_std \
--ad_e_g_inc $ad_e_g_inc \
--deltaT_e $deltaT_e \
--taum_i $taum_i \
--taum_i_spread $taum_i_spread \
--EL_i $EL_i \
--EL_i_spread $EL_i_spread \
--Vt_i $Vt_i \
--Vr_i $Vr_i \
--Rm_i $Rm_i \
--ad_tau_i_mean $ad_tau_i_mean \
--ad_tau_i_std $ad_tau_i_std \
--ad_i_g_inc $ad_i_g_inc \
--deltaT_i $deltaT_i \
--tau_AMPA $tau_AMPA \
--g_AMPA_total $g_AMPA_total \
--g_AMPA_std $g_AMPA_std \
--tau_GABA_rise $tau_GABA_rise \
--tau_GABA_fall $tau_GABA_fall \
--g_GABA_total $g_GABA_total \
--Vrev_AMPA $Vrev_AMPA \
--Vrev_GABA $Vrev_GABA \
--noise_sigma $noise_sigma \
--sigma_init_cond $sigma_init_cond \
--refrac_abs $refrac_abs \
--time $time \
--sim_dt $sim_dt \
--spike_detect_th $spike_detect_th \
--Vclamp $Vclamp \
--output_dir $output_dir \
--update_interval $update_interval \
--job_num $job_num \
--ntrials $ntrials \


