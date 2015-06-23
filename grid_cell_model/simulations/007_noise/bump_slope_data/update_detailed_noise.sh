#!/bin/bash

root_dir="../output/main_network/detailed_noise/velocity/"
target_file="bump_slope_detailed_noise_"

I_vec="EI-1_3 EI-3_1"
for I in $I_vec; do
    src=${root_dir}${I}/reductions.h5
    dst=${target_file}${I}.h5
    echo "removing $dst"
    rm -f $dst
    echo "src: $src"
    echo "dst: $dst"
    h5copy -i $src -o $dst -s /analysis/lineFitSlope -d /lineFitSlope
done
