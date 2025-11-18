#!/bin/bash

# Pesicek & Ryberg generic Vp volcano model:
# Vp = 0.000101*z^3 - 0.007799*z^2 + 0.241784*z + 4.301992 

# Check if at least two arguments provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <z_start> <z_end> [z_step]"
    echo "Example: $0 -3 30 0.1"
    exit 1
fi

# Set default values
z_start=$1
z_end=$2
z_step=${3:-0.1}

awk -v start="$z_start" -v end="$z_end" -v step="$z_step" 'BEGIN {
    for (z = start; z <= end; z += step) {
	Vp = 0.000101*z^3 - 0.007799*z^2 + 0.241784*z + 4.301992  
      printf "%.2f %.2f\n", Vp, z
    }
}'


