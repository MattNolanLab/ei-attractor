#!/bin/bash

input="0.2 0.25 0.3 0.35"

for i in $input; do
    python fiete_path_integration.py -s 40 -w -t 5 -i $i -u 1
done
