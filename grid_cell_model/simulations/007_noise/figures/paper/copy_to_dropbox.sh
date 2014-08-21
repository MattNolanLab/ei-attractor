#!/bin/bash

dst_dir=~/Dropbox/shared_data/noise/figures/
output_fig_dir='output_figures'
panels_dir='panels'

files="
    $output_fig_dir/figure1.pdf
    $output_fig_dir/figure2.pdf
    $output_fig_dir/figure3.pdf
    $output_fig_dir/figure4.pdf
    $output_fig_dir/figure_drifts.pdf
    $output_fig_dir/figure_isBump.pdf
    $output_fig_dir/suppFigure_bump_sweeps.pdf
    $panels_dir/suppfigure_gamma.pdf
    $panels_dir/suppFigure_gammaF_grids_scatter.pdf
    $output_fig_dir/suppFigure_grids_vs_line_fit_err.pdf
    $output_fig_dir/suppFigure_grids_vs_line_slope.pdf
    $output_fig_dir/suppFigure_line_fit_error_vs_slope.pdf
    $output_fig_dir/suppFigure_seizures.pdf
    $output_fig_dir/suppFigure_velocity.pdf
    $output_fig_dir/suppFigure_gamma_probabilities.pdf
    $output_fig_dir/suppFigure_grids_probability.pdf
    $panels_dir/raster_examples*.pdf
    $panels_dir/gamma_scatter_gamma_pbumps_all.pdf
    $panels_dir/suppFigure_grid_examples.pdf
    $panels_dir/suppFigure_grids_vs_bumps_exp.pdf
    $panels_dir/suppFigure_grids_vs_bumps.pdf
    "

local_files="
    $panels_dir/suppFigure_gamma.pdf
    $panels_dir/suppFigure_gammaF_grids_scatter.pdf
    $panels_dir/suppFigure_grids_vs_bumps_exp.pdf
    $panels_dir/suppFigure_grids_vs_bumps.pdf
"

# Copy to dropbox
echo "Copying to dropbox..."
echo $files
echo
ls $dst_dir
for f in $files
do
    cp -rv $f $dst_dir
done

cp -r panels $dst_dir/
echo "Done"

# Copy local files to output_figures
echo "Copying from panels to output_figures"
for f in $local_files
do
    cp -rv $f $output_fig_dir
done
echo "Done"
