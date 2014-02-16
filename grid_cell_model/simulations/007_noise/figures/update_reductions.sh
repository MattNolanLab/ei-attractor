#!/bin/bash
#
#   update_reductions.py
#
#   Update reductions from central data store to the local directory.
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

noise_sigmas="0 150 300"
detailed_dirs="EI-1_3 EI-3_1"

src_dir="../../../../central_data_store/simulation_data/grid_cell_model/007_noise/"
dst_dir="output_local/"
sweep_dir="even_spacing"
detailed_noise="detailed_noise"
reduction_file="reductions.h5"

# no arguments --> update everything
args=$@
if [[ $# -eq 0 || ($# -eq 1 && $1 == "--detailed" ) ]]; then
    args="gamma_bump velocity grids grids_no_velocity const_position"
fi

if [ "$1" == "--detailed" ]; then
    detailed=true
else
    detailed=false
fi

for arg in $args; do
    if [ $arg == "--detailed" ];then
        continue
    fi

    echo
    echo

    if $detailed; then
        echo "Updating the detailed noise reduction files."
        for detailed_dir in $detailed_dirs; do
            if [ $arg == "gamma_bump" ]; then
                subdir="$detailed_noise/gamma_bump/$detailed_dir/$reduction_file"
            elif [ $arg == "velocity" ]; then
                subdir="$detailed_noise/velocity/$detailed_dir/$reduction_file"
            elif [ $arg == "grids" ]; then
                subdir="$detailed_noise/grids/$detailed_dir/$reduction_file"
            else
                echo "Unknown detailed noise option: $arg"
                exit 1
            fi

            resolved_src="$src_dir/$subdir"
            resolved_dst="$dst_dir/$subdir"
            cp -v "$resolved_src" "$resolved_dst"
        done

    else
        echo "Updating the even_spacing reduction files."
        for noise_sigma in $noise_sigmas; do
            if [ $arg == "gamma_bump" ]; then
                subdir="$sweep_dir/gamma_bump/${noise_sigma}pA/$reduction_file"
            elif [ $arg == "velocity" ]; then
                subdir="$sweep_dir/velocity/${noise_sigma}pA/$reduction_file"
            elif [ $arg == "grids" ]; then
                subdir="$sweep_dir/grids/${noise_sigma}pA/$reduction_file"
            elif [ $arg == "grids_no_velocity" ]; then
                subdir="$sweep_dir/grids_no_velocity/${noise_sigma}pA/$reduction_file"
            else
                echo "Non-specific even spacing option detected."
                subdir="$sweep_dir/$arg/${noise_sigma}pA/$reduction_file"
            fi

            echo "$noise_sigma pA..."
            resolved_src="$src_dir/$subdir"
            resolved_dst="$dst_dir/$subdir"
            cp -v "$resolved_src" "$resolved_dst"
        done
    fi
done
