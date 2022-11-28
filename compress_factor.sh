#!/bin/bash


base_file="$1"
new_file="$2"

T0=`ls -l "$base_file" | awk '{print $5}'`
T1=`ls -l "$new_file" | awk '{print $5}'`

echo "scale=4 ; $T0 / $T1" | bc