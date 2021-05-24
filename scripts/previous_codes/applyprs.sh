#!/bin/bash

mkdir -p ../files/multiprs/

for prs in ../files/multiprs_ss/*; do
    python3 42_multiprs.py $prs
done
