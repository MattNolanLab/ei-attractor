#! /bin/bash

CWD=`pwd`
DEMO_DIR=gridcells_brian_demo

rm -f *.pdf *.png
cd ..
zip $DEMO_DIR $DEMO_DIR $DEMO_DIR/*.py

cd $CWD
