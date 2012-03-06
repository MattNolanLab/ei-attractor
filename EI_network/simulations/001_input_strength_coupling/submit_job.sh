#!/bin/sh

########################################
# Submit fiete_path_integration.py jobs
# with different network parameters
# into the cluster
########################################

BASE=../../
#export PYTHONPATH="/disk/scratch/s0966762/lib/python2.6/site-packages:$BASE"
export PYTHONPATH="$BASE"
echo $PYTHONPATH
EDDIE=0  # if eddie, submit on a cluster using qsub

QSUB_PARAMS="-N mult_bump_spiking_net -P inf_ndtc -cwd -l h_rt=00:50:00"


Ne=30
Ni=30

#Iext_coeff="0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6"
#coupling_coeff="0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5"
Iext_coeff="1.6"
coupling_coeff="1.5"

AMPA_density="1.0"
GABA_density="1.0"

Iext_e_1="475*10^-12"
Iext_i_1="150*10^-12"

taum_e="9.3e-3"
EL_e="-68.5e-3"
Vt_e="-50e-3"
Vr_e=$EL_e
Rm_e="44e6"
ad_tau_e_mean="1e-3"
ad_tau_e_std="0e-3"
ad_e_g_inc="0e-8"
deltaT_e="1.5e-3"

taum_i="10e-3"
EL_i="-60e-3"
Vt_i="-50e-3"
Vr_i="$EL_i"
Rm_i="44e6"
ad_tau_i_mean="3e-3"
ad_tau_i_std="0e-3"  # Unused in the simulation for now
ad_i_g_inc="0e-8"
deltaT_i="1.5e-3"

tau_AMPA="2e-3"
g_AMPA_total_1="1500*10^-9"
g_AMPA_std="10e-12"
tau_GABA_rise="1e-3"
tau_GABA_fall="9e-3"
g_GABA_total_1="900*10^-9"

Vrev_AMPA="0e-3"
Vrev_GABA="-75e-3"

noise_sigma="2e-3"
sigma_init_cond="10e-3"

refrac_abs="0.1e-3"

time=30
sim_dt="0.5e-3"
spike_detect_th="20e-3"

ntrials=20

output_dir="output"
update_interval=30


job_num=1

for Iext_c in $Iext_coeff; do
    for coupling_c in $coupling_coeff; do
#####################
    Iext_e=`echo "$Iext_e_1 * $Iext_c" | bc -l`
    Iext_i=`echo "$Iext_i_1 * $Iext_c" | bc -l`

    g_AMPA_total=`echo "$g_AMPA_total_1 * $coupling_c" | bc -l`
    g_GABA_total=`echo "$g_GABA_total_1 * $coupling_c" | bc -l`


    if [ $EDDIE -eq 1 ]
    then
        qsub $QSUB_PARAMS eddie_submit.sh \
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
            --job_num $job_num \
            --ntrials $ntrials
    else
        pwd
        nice python simulation.py \
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
            --job_num $job_num \
            --ntrials $ntrials \

    fi

    echo

    let job_num=$job_num+1
#####################
    done
done
