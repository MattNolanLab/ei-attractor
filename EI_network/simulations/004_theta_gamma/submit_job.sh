#!/bin/sh

########################################
# Submit fiete_path_integration.py jobs
# with different network parameters
# into the cluster
########################################

BASE=../../
export PYTHONPATH="/exports/work/inf_ndtc/s0966762/python-modules/lib/python2.6/site-packages:$BASE"
EDDIE=1  # if eddie, submit on a cluster using qsub

QSUB_PARAMS="-P inf_ndtc -cwd -l h_rt=02:30:00"

net_generations=10

Ne=400
Ni=100

AMPA_density="0.4"
GABA_density="0.4"

Iext_coeff="0.8 0.9 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4"
AMPA_coeff="0.8 0.9 1.0 1.05 1.1 1.15 1.2 1.25 1.4 1.35 1.4"

Iext_e_1="900*10^-12"
Iext_i_1="250*10^-12"

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
g_AMPA_total_1="3.5*10^-8"
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
sim_dt="0.05e-3"
spike_detect_th="40e-3"
Vclamp="-50e-3"

ntrials=2

output_dir="output"
update_interval=10
job_num=1000


for Iext_c in $Iext_coeff; do
    for AMPA_c in $AMPA_coeff; do
#####################
    Iext_e=`echo "$Iext_e_1 * $Iext_c" | bc -l`
    Iext_i=`echo "$Iext_i_1 * $Iext_c" | bc -l`

    g_AMPA_total=`echo "$g_AMPA_total_1 * $AMPA_c" | bc -l`

    echo $Iext_e, $Iext_i, $g_AMPA_total

    if [ $EDDIE -eq 1 ]
    then
        job_num_str=`printf "job%04d" $job_num`
        qsub $QSUB_PARAMS -N $job_num_str  -j y -o $output_dir \
            eddie_submit.sh \
            --net_generations $net_generations \
            --Ne $Ne \
            --Ni $Ni \
            --AMPA_density $AMPA_density \
            --GABA_density $GABA_density \
            --Iext_e $Iext_e \
            --Iext_i $Iext_i \
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
            --ntrials $ntrials
    else
        pwd
        nice python -i simulation.py \
        --net_generations $net_generations \
        --Ne $Ne \
        --Ni $Ni \
        --AMPA_density $AMPA_density \
        --GABA_density $GABA_density \
        --Iext_e $Iext_e \
        --Iext_i $Iext_i \
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

    fi

    echo

    let job_num=$job_num+1
#####################
    done
done
