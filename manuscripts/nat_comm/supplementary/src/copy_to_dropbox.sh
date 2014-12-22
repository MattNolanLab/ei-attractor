#!/bin/bash

dst_dir=~/Dropbox/Independent_grids_and_gamma/noise/

files="
    SI_text.pdf
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


