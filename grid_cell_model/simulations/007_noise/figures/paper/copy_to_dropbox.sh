#!/bin/bash

dst_dir=~/Dropbox/shared_data/noise/figures/
output_fig_dir='output_figures'
no_theta_output_fig_dir="no_theta/$output_fig_dir"
panels_dir='panels'

files="
    $output_fig_dir/figure1.pdf
    $output_fig_dir/figure2.pdf
    $output_fig_dir/figure3.pdf
    $output_fig_dir/figure4.pdf
    $output_fig_dir/figureS1.pdf
    $output_fig_dir/figureS2.pdf
    $output_fig_dir/figureS3.pdf
    $output_fig_dir/figureS4.pdf
    $output_fig_dir/figureS5.pdf
    $output_fig_dir/figureS6.pdf
    $output_fig_dir/figureS7.pdf
    $output_fig_dir/figureS8.pdf
    $output_fig_dir/figureS9.pdf
    $output_fig_dir/figureS10.pdf

    $output_fig_dir/suppFigure_grids_vs_line_fit_err.pdf
    $output_fig_dir/suppFigure_grids_vs_line_slope.pdf
    $output_fig_dir/suppFigure_line_fit_error_vs_slope.pdf
    $panels_dir/suppFigure_grids_vs_bumps_exp.pdf
    $panels_dir/suppFigure_grids_vs_bumps.pdf
    "

no_theta_files="
    $no_theta_output_fig_dir/figure2.pdf
    $no_theta_output_fig_dir/figure3.pdf
    $no_theta_output_fig_dir/figure3_details.pdf
    $no_theta_output_fig_dir/figure4.pdf
    $no_theta_output_fig_dir/figureS8.pdf
"

# Copy to dropbox
echo "Copying to main figures to dropbox..."
echo $files
echo
ls $dst_dir
for f in $files
do
    cp -rv $f $dst_dir
done

cp -r panels $dst_dir/


echo "Copying no theta figures to dropbox..."
no_theta_dst_dir=$dst_dir/no_theta
echo $no_theta_files
echo
ls $no_theta_dst_dir
for f in $no_theta_files
do
    cp -rv $f $no_theta_dst_dir
done

cp -r no_theta/panels $no_theta_dst_dir/
echo "Done"
