#!/bin/sh

BASE=./
#export PYTHONPATH="/disk/scratch/s0966762/lib/python2.6/site-packages:$BASE"
export PYTHONPATH="$BASE"
echo $PYTHONPATH
EDDIE=0  # if eddie, submit on a cluster using qsub

QSUB_PARAMS="-N mult_bump_spiking_net -P inf_ndtc -cwd -l h_rt=00:15:00 -pe memory-2G 2"

# Network parameters
P_N_i="0   1   2   3   4   5   6   7   8   9  10  11  12  13   14"
P_Ne=(50 100 150 200 250 300 350 400 450 500 600 700 800 900 1000)
P_Ni=(12  25  35  50  60  75  85 100 110 125 150 175 200 225  250)

REPEAT=2


AMPA_density="0.8"
GABA_density="0.8"

Iext_e="475e-12"
Iext_i="100e-12"

taum_e="9.3e-3"
EL_e="-68.5e-3"
Vt_e="-50e-3"
Vr_e=$EL_e
Rm_e="44e6"
ad_tau_e_mean="1e-3"
ad_tau_e_std="0e-3"
ad_e_g_inc="1.136e-8"
deltaT_e="1.5e-3"

taum_i="10e-3"
EL_i="-60e-3"
Vt_i="-50e-3"
Vr_i="$EL_i"
Rm_i="44e6"
ad_tau_i_mean="3e-3"
ad_tau_i_std="0e-3"  # Unused in the simulation for now
ad_i_g_inc="2.27e-8"
deltaT_i="1.5e-3"

tau_AMPA="2e-3"
g_AMPA_total="700e-9"
g_AMPA_std="50e-12"
tau_GABA_rise="1e-3"
tau_GABA_fall="9e-3"
g_GABA_total="400e-9"

Vrev_AMPA="0e-3"
Vrev_GABA="-75e-3"

noise_sigma="0.5e-3"
sigma_init_cond="10e-3"

refrac_abs="0.1e-3"

time=5
sim_dt="0.1e-3"
spike_detect_th="20e-3"

output_dir="output_local/net_size/no_stim"
update_interval=5


job_id=1000
for it in $P_N_i; do
#####################
Ne=${P_Ne[it]}
Ni=${P_Ni[it]}

repeat=1
while [ $repeat -le $REPEAT ]; do
    echo "Ne=$Ne"
    echo

    nice python  test_EI_network.py \
    --Ne $Ne \
    --Ni $Ni \
    --AMPA_density $AMPA_density \
    --GABA_density $GABA_density \
    --Iext_e $Iext_e \
    --Iext_i $Iext_i \
    --taum_e $taum_e \
    --EL_e $EL_e \
    --Vt_e $Vt_e \
    --Vr_e $Vr_e \
    --Rm_e $Rm_e \
    --ad_tau_e_mean $ad_tau_e_mean \
    --ad_tau_e_std $ad_tau_e_std \
    --ad_e_g_inc $ad_e_g_inc \
    --deltaT_e $deltaT_e \
    --taum_i $taum_i \
    --EL_i $EL_i \
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
    --output_dir $output_dir \
    --update_interval $update_interval \
    --job_num $job_id&

    echo


    let job_id=$job_id+1
    let repeat=$repeat+1
done
#####################

done


