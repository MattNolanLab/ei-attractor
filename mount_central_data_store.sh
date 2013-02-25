#!/bin/bash
mkdir -p central_data_store

sshfs \
s0966762@eddie.ecdf.ed.ac.uk:/exports/work/inf_ndtc/s0966762/central_data_store \
central_data_store -o nonempty
