#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bi>"
    exit 1
fi
# evaluated models threshold:
emt=$1

# run all post processing scripts at once with defaults:
dispe.sh
disp_m_average_sl.sh $emt 100 0
disp_eq.sh picks.mcmc

# optionally
disp_noise.sh
disp_eq_z.sh $emt
disp_eq_z2.sh $emt
disp_eq_evo.sh $emt
disp_msft_dist.sh picks.mcmc

# output results files
outputModels.sh

# compare to input synthetics (or pre-existing model and locations)
disp_compare.sh

# new station correction plots
plotStaCorHists.sh
