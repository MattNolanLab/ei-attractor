#! /bin/bash

DEMO_DIR=gridcells_brian_demo
CWD=`pwd`

cd $DEMO_DIR
rm *.pdf *.png
cd $CWD

zip $DEMO_DIR $DEMO_DIR $DEMO_DIR/*.py
