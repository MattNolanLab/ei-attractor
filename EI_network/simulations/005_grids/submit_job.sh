#!/bin/sh

########################################
# Submit fiete_path_integration.py jobs
# with different network parameters
# into the cluster
########################################

BASE=../../
PYTHON_JUP="/disk/scratch/s0966762/lib/python2.6/site-packages"
export PYTHONPATH="/exports/work/inf_ndtc/s0966762/python-modules/lib/python2.6/site-packages:$BASE:$PYTHON_JUP"
EDDIE=1  # if eddie, submit on a cluster using qsub
dry_run=0

QSUB_PARAMS="-P inf_ndtc -cwd -l h_rt=03:00:00 -pe memory-2G 2"

net_generations=1

P_Ivel_max="40e-12 50e-12 75e-12 100e-12 150e-12"
P_pAMPA_sigma="0.7"

Ne=68
Ni=34

P_prefDirC="4"

Iext_e_coeff="0.75"
Iext_i_coeff="0.9" #"0.4 0.5 0.6 0.7 0.8 0.9"
Iext_e_min_coeff="0.8"
AMPA_coeff="40"
GABA_coeff="54"
extraGABA_coeff="0.7"
adapt_inc_coeff="1.0" #"1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9"
adapt_coeff="1.0" #"0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5"

heterog_e_coeff="0.1"
heterog_i_coeff="0.0"

tau_GABA_rise_coeff="0.1"
tau_GABA_fall_coeff="1.0"

# Percent NMDA of AMPA synapse
P_NMDA_amount="2"

Iext_e_1="900*10^-12"
Iext_i_1="250*10^-12"
Iext_e_min_1="375*10^-12"
Iext_i_min="200e-12"


taum_e="9.3e-3"
taum_e_spread_1="3.1*10^-3"
EL_e="-68.5e-3"
EL_e_spread_1="2.0*10^-3"
Vt_e="-50e-3"
Vr_e="-68.5e-3"
Rm_e="44e6"
deltaT_e="0.4e-3"
Eahp_e="-80e-3"
g_ahp_e="5e-9"
tau_ahp_e="20e-3"

taum_i="10e-3"
taum_i_spread_1="4*10^-3"
EL_i="-60e-3"
EL_i_spread_1="20*10^-3"
Vt_i="-45e-3"
Vr_i="$EL_i"
Rm_i="44e6"
ad_tau_i_mean="7.5e-3"
ad_tau_i_std="0.5e-3"  # Unused in the simulation for now
ad_i_g_inc="2.273e-8"
deltaT_i="0.4e-3"

tau_AMPA="1e-3"
g_AMPA_total_1="3.5*10^-8"
g_AMPA_std_1="600*10^-12"
tau_GABA_rise_1="1*10^-3"
tau_GABA_fall_1="5*10^-3"
g_GABA_total_1="4.00*10^-8"

# Uniform GABA I-->E only
g_extraGABA_total_1="4.00*10^-8"
extraGABA_density="0.4"

Vrev_AMPA="0e-3"
Vrev_GABA="-75e-3"

noise_sigma="2e-3"
sigma_init_cond="10e-3"

refrac_abs="0.1e-3"

time=200
sim_dt="0.1e-3"
spike_detect_th="40e-3"
Vclamp="-50e-3"

ntrials=1

output_dir="output"
readme_file="$output_dir/README_JOBS_`date "+%Y_%m_%dT%H_%M_%S"`"
update_interval=10
job_num=2000



for Iext_e_c in $Iext_e_coeff; do
    for Iext_i_c in $Iext_i_coeff; do
        for AMPA_c in $AMPA_coeff; do
            for GABA_c in $GABA_coeff; do
                for adapt_inc_c in $adapt_inc_coeff; do
                    for adapt_c in $adapt_coeff; do
                        for Ivel_max in $P_Ivel_max; do
                            for pAMPA_sigma in $P_pAMPA_sigma; do
                                for heterog_e_c in $heterog_e_coeff; do
                                    for heterog_i_c in $heterog_i_coeff; do
                                        for tau_GABA_rise_c in $tau_GABA_rise_coeff; do
                                            for tau_GABA_fall_c in $tau_GABA_fall_coeff; do
                                                for extraGABA_c in $extraGABA_coeff; do
                                                    for NMDA_amount in $P_NMDA_amount; do
                                                        for prefDirC in $P_prefDirC; do
                                                            for Iext_e_min_c in $Iext_e_min_coeff; do
net_it=0
while [ $net_it -lt $net_generations ]; do
#####################
    Iext_e=`echo "$Iext_e_1 * $Iext_e_c" | bc -l`
    Iext_i=`echo "$Iext_i_1 * $Iext_i_c" | bc -l`
    Iext_e_min=`echo "$Iext_e_min_1 * $Iext_e_min_c" | bc -l`

    g_AMPA_total=`echo "$g_AMPA_total_1 * $AMPA_c" | bc -l`
    g_AMPA_std=`echo "$g_AMPA_std_1 * $AMPA_c" | bc -l`
    g_GABA_total=`echo "$g_GABA_total_1 * $GABA_c" | bc -l`
    g_extraGABA_total=`echo "$g_extraGABA_total_1 * $extraGABA_c" | bc -l`

    taum_e_spread=`echo "$taum_e_spread_1 * $heterog_e_c" | bc -l`
    EL_e_spread=`echo "$EL_e_spread_1 * $heterog_e_c" | bc -l`

    taum_i_spread=`echo "$taum_i_spread_1 * $heterog_i_c" | bc -l`
    EL_i_spread=`echo "$EL_i_spread_1 * $heterog_i_c" | bc -l`
    
    tau_GABA_rise=`echo "$tau_GABA_rise_1 * $tau_GABA_rise_c" | bc -l`
    tau_GABA_fall=`echo "$tau_GABA_fall_1 * $tau_GABA_fall_c" | bc -l`

    if [ $dry_run -eq 1 ]
    then
        echo "job`printf "%04d" $job_num`" >> $readme_file
        echo "    Iext_e       = `printf "%1.3e" $Iext_e`" >> $readme_file
        echo "    Iext_i       = `printf "%1.3e" $Iext_i`" >> $readme_file
        echo "    g_AMPA_total = `printf "%1.3e" $g_AMPA_total`" >> $readme_file
        echo "    g_AMPA_std   = `printf "%1.3e" $g_AMPA_std`" >> $readme_file
        echo "    g_GABA_total = `printf "%1.3e" $g_GABA_total`" >> $readme_file
        echo
    else
        if [ $EDDIE -eq 1 ]
        then
            job_num_str=`printf "job%04d" $job_num`
            qsub $QSUB_PARAMS -N $job_num_str  -j y -o $output_dir \
                eddie_submit.sh \
                --Ivel_max $Ivel_max \
                --pAMPA_sigma $pAMPA_sigma \
                --Ne $Ne \
                --Ni $Ni \
                --Iext_e $Iext_e \
                --Iext_i $Iext_i \
                --Iext_e_min $Iext_e_min \
                --Iext_i_min $Iext_i_min \
                --taum_e $taum_e \
                --taum_e_spread $taum_e_spread \
                --EL_e $EL_e \
                --EL_e_spread $EL_e_spread \
                --Vt_e $Vt_e \
                --Vr_e $Vr_e \
                --Rm_e $Rm_e \
                --deltaT_e $deltaT_e \
                --Eahp_e $Eahp_e \
                --g_ahp_e $g_ahp_e \
                --tau_ahp_e $tau_ahp_e \
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
                --NMDA_amount $NMDA_amount \
                --tau_GABA_rise $tau_GABA_rise \
                --tau_GABA_fall $tau_GABA_fall \
                --g_GABA_total $g_GABA_total \
                --g_extraGABA_total $g_extraGABA_total \
                --extraGABA_density $extraGABA_density \
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
                --prefDirC $prefDirC \
                --ntrials $ntrials
        else
            pwd
            nice python2.6 simulation.py \
            --Ivel_max $Ivel_max \
            --pAMPA_sigma $pAMPA_sigma \
            --Ne $Ne \
            --Ni $Ni \
            --Iext_e $Iext_e \
            --Iext_i $Iext_i \
            --Iext_e_min $Iext_e_min \
            --Iext_i_min $Iext_i_min \
            --taum_e $taum_e \
            --taum_e_spread $taum_e_spread \
            --EL_e $EL_e \
            --EL_e_spread $EL_e_spread \
            --Vt_e $Vt_e \
            --Vr_e $Vr_e \
            --Rm_e $Rm_e \
            --deltaT_e $deltaT_e \
            --Eahp_e $Eahp_e \
            --g_ahp_e $g_ahp_e \
            --tau_ahp_e $tau_ahp_e \
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
            --NMDA_amount $NMDA_amount \
            --tau_GABA_fall $tau_GABA_fall \
            --g_GABA_total $g_GABA_total \
            --g_extraGABA_total $g_extraGABA_total \
            --extraGABA_density $extraGABA_density \
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
            --prefDirC $prefDirC \
            --ntrials $ntrials&

        fi
    fi

    echo

    let job_num=$job_num+1
    let net_it=$net_it+1
    sleep 0.5
#####################
done
                                                            done
                                                        done
                                                    done
                                                done
                                            done
                                        done
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done
