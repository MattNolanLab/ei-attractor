#!/bin/bash
#
#   update_reductions.sh
#
#   These are the I-surround simulations in which all other paramters are set
#   up so that they match the original E-surround config (main figures).
#
#   Copy the reduction files as the bump_slopeXXX.h5. The reduction files have
#   to contain the lineFitErr and lineFitSlope data sets (with the correct
#   dimensions).
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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

root_dir="../output/submission/i_surround/pastoll_et_al/velocity/"
target_file="bump_slope_i_surround_pastoll_"

I_vec="0pA 150pA 300pA"
for I in $I_vec; do
    src=${root_dir}${I}/reductions.h5
    dst=${target_file}${I}.h5
    echo "removing $dst"
    rm -f $dst
    echo "src: $src"
    echo "dst: $dst"
    h5copy -i $src -o $dst -s /analysis/lineFitSlope -d /lineFitSlope
    h5diff -r $src $dst /analysis/lineFitSlope /lineFitSlope >log_${I}
done
