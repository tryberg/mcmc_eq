#!/bin/bash
# Convert a model.inp file (used for plotting) to a model.dat file
# required for mcmc_eq fixed velocity format (see manual).
#
# Usage: ./inp2dat.sh model.inp
#
# Input format: Depth Vp Vp/Vs   (3 columns per line)
# Note: # of layers must equal nz in the mcmc_eq cfg file.

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <model.inp>" >&2
    exit 1
fi

synth_model=$1

if [[ ! -f "$synth_model" ]]; then
    echo "Error: input file '$synth_model' not found" >&2
    exit 1
fi

# Check every line has exactly 3 columns
bad_nc=$(awk 'NF != 3 { print NR": "$0; count++ } END { print count+0 }' "$synth_model")
bad_nc_count=$(echo "$bad_nc" | tail -1)
if [[ "$bad_nc_count" -ne 0 ]]; then
    echo "Error: '$synth_model' has lines without exactly 3 columns:" >&2
    echo "$bad_nc" | head -n -1 >&2
    exit 1
fi

# Check column 3 (Vp/Vs) < column 2 (Vp) on every line
bad_ratio=$(awk '$3 >= $2 { print NR": "$0; count++ } END { print count+0 }' "$synth_model")
bad_ratio_count=$(echo "$bad_ratio" | tail -1)
if [[ "$bad_ratio_count" -ne 0 ]]; then
    echo "Error: column 3 (Vp/Vs) must be less than column 2 (Vp) on every line." >&2
    echo "Violations in '$synth_model':" >&2
    echo "$bad_ratio" | head -n -1 >&2
    exit 1
fi

# Construct output (could also be done with velModTable2trondVelFiles)
#awk '{print "STAN", $1, $2, 0, $3, 0, $2, 0, $3, 0, $2, $3, 0.01}' "$synth_model"
awk '{printf "STAN %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",$1, $2, 0, $3, 0, $2, 0, $3, 0, $2, $3, 0.01}' "$synth_model"
