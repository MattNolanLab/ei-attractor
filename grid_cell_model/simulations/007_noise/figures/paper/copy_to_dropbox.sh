#!/bin/bash

dst_dir=~/Dropbox/shared_data/noise/figures/
output_fig_dir='output_figures'

files="
    $output_fig_dir/figure1.pdf
    $output_fig_dir/figure2.pdf
    $output_fig_dir/figure3.pdf
    $output_fig_dir/figure4.pdf
    $output_fig_dir/figure_drifts.pdf
    $output_fig_dir/figure_isBump.pdf
    $output_fig_dir/suppFigure_bump_sweeps.pdf
    $output_fig_dir/suppFigure_gamma.pdf
    $output_fig_dir/suppFigure_grids_vs_bumps_exp.pdf
    $output_fig_dir/suppFigure_grids_vs_bumps.pdf
    $output_fig_dir/suppFigure_grids_vs_line_fit_err.pdf
    $output_fig_dir/suppFigure_grids_vs_line_slope.pdf
    $output_fig_dir/suppFigure_line_fit_error_vs_slope.pdf
    $output_fig_dir/suppFigure_seizures.pdf
    $output_fig_dir/suppFigure_velocity.pdf
    "

echo $files
echo
ls $dst_dir
for f in $files
do
    cp -r $f $dst_dir
done

cp -r panels $dst_dir/
