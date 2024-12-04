#!/bin/bash

# Example script from Jeremy to run in parallel on HPC cluster using srun, sbatch, slurm
# tuned for HOVENWEEP at USGS ARC: 
# https://www.usgs.gov/advanced-research-computing/usgs-hovenweep-supercomputer
# which has a default time limit of 2 days for jobs
# remember to launch this script as "sbatch srun_mcmc_eq.sh"

#SBATCH --nodes=1
#SBATCH --account=vhp
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jpesicek@usgs.gov
#SBATCH --array=1-100

# for jobs longer than 2 days:
#SBATCH --qos=seven_days_max
#SBATCH --time=5-0

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

id=`echo $SLURM_ARRAY_TASK_ID | awk '{printf "%03d\n",$0}'`

cfg=config_eqx.dat	# configuration file
exe=mcmc_eq 		# executable
mdf=picks.mcmc	 	# pickfile

which $exe
ls $cfg $mdf ./$exe
test -f ./$exe || exit
test -f $cfg || exit
test -f $mdf || exit

srun $exe $cfg rjx-${id}_${SLURM_ARRAY_JOB_ID}.out $mdf
