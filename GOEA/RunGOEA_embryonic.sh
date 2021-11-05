#!/bin/bash

#SBATCH --job-name=runGOEA                          # job name
#SBATCH --partition=128GB                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=1-00:00:00                                  # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=JobLog.%j.out                         # standard output file name
#SBATCH --error=JobLog.%j.err                         # standard error output file name
#SBATCH --mail-user=boxun.li@utsouthwestern.edu         # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

echo "hello world"

module load python/3.7.x-anaconda
conda activate py37_res_diffxpy

python Embryonic_RunGOEA.py
