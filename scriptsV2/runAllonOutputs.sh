#!/bin/bash

# run all post processing scripts at once with defaults:
dispe.sh
disp_m_average_sl.sh 225000 100 0
disp_eq.sh picks.mcmc

# optionally
disp_noise.sh
disp_eq_z.sh 225000
disp_eq_evo.sh
disp_msft_dist.sh picks.mcmc

# output results files
outputModels.sh

# compare to input synthetics (or pre-existing model and locations)
disp_compare.sh
