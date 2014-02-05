#!/bin/bash
#
#   copy_to_dropbox.sh
#
#   Copy figures to the shared folder in dropbox
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

dst_dir=~/Dropbox/shared_data/noise/figures/
output_fig_dir='output_figures'

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
    $output_fig_dir/suppFigure_grid_sweeps.pdf
    $output_fig_dir/suppFigure_grids_FR-histograms.pdf
    $output_fig_dir/suppFigure_grids_vs_line_fit_err.pdf
    $output_fig_dir/suppFigure_grids_vs_line_slope.pdf
    $output_fig_dir/suppFigure_bump_sweeps.pdf
    $output_fig_dir/suppFigure_line_fit_error_vs_slope.pdf
    $output_fig_dir/suppFigure_no_theta.pdf
    "

echo $files
echo
ls $dst_dir
for f in $files
do
    cp -r $f $dst_dir
done

cp -r panels $dst_dir/
