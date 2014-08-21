#!/bin/bash

dst_dir=~/work/PhdThesis/src/fig/
output_fig_dir='panels'

files="
    $output_fig_dir/grids_detailed_noise_gscore.pdf
    $output_fig_dir/grids_diff_sweep0.pdf
    $output_fig_dir/gamma_scatter_gamma_grids_all.pdf
    $output_fig_dir/bumps_scatter_grids_vs_bumpFracTotal_exp.pdf
    $output_fig_dir/bumps_scatter_grids_vs_bumpFracTotal.pdf
    $output_fig_dir/raster_examples*.pdf
    $output_fig_dir/suppFigure_grid_examples*.pdf
    $output_fig_dir/suppFigure_gamma.pdf
    $output_fig_dir/suppFigure_grids_vs_bumps.pdf
    "

for f in $files
do
    cp -v $f $dst_dir/
done


