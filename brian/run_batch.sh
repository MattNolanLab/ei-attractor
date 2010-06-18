#!/bin/bash

size="32 40 64"

for i in $size; do
    python attractor_network_fiete.py -s $i -w
done
