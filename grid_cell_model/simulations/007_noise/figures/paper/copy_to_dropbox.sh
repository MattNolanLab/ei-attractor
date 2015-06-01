#!/bin/bash

dst_dir=~/Dropbox/Independent_grids_and_gamma/noise/e_life_v2/figures
supp_dst_dir=~/Dropbox/Independent_grids_and_gamma/noise/e_life_v2/figure_supplements
output_fig_dir='output_figures'
panels_dir='panels'

files="
$output_fig_dir/Figure1.pdf
$output_fig_dir/Figure2.pdf
$output_fig_dir/Figure3.pdf
$output_fig_dir/Figure4.pdf
$output_fig_dir/Figure5.pdf
$output_fig_dir/Figure6.pdf
$output_fig_dir/Figure7.pdf
"

supplementary_files="
$output_fig_dir/Figure1_S1.pdf
$output_fig_dir/Figure2_S2.pdf
$output_fig_dir/Figure2_S3.pdf
$output_fig_dir/Figure2_S4.pdf
$output_fig_dir/Figure2_S5.pdf
$output_fig_dir/Figure3_S1.pdf
$output_fig_dir/Figure3_S2.pdf
$output_fig_dir/Figure3_S3.pdf
$output_fig_dir/Figure3_S4.pdf
$output_fig_dir/Figure6_S1.pdf
$output_fig_dir/Figure6_S2.pdf
$output_fig_dir/Figure6_S3.pdf
$output_fig_dir/Figure6_S4.pdf
$output_fig_dir/Figure6_S5.pdf
$output_fig_dir/Figure6_S6.pdf
$output_fig_dir/Figure7_S1.pdf
$output_fig_dir/Figure7_S2.pdf
$output_fig_dir/Figure7_S3.pdf
$output_fig_dir/Figure7_S4.pdf
$output_fig_dir/Figure7_S5.pdf
$output_fig_dir/Figure7_S6.pdf
$output_fig_dir/Figure7_S7.pdf
$output_fig_dir/Figure7_S8.pdf
$output_fig_dir/Figure7_S9.pdf
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

#cp -r panels $dst_dir/

# Copy to dropbox
echo "Copying to supplementary figures to dropbox..."
echo $supplementary_files
echo
ls $supp_dst_dir
for f in $supplementary_files
do
    cp -rv $f $supp_dst_dir
done

