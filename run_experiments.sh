#!/bin/bash -l	
#SBATCH --job-name="ericBenchmarking"
#SBATCH --output="slurm.log"
#SBATCH --partition=parallel
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
 
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
conda init
conda activate ggrn


for experiment in `ls -1 experiments  | grep -E $1`
do
    echo "Starting ${experiment}"
    echo "Monitor progress:
less experiments/${experiment}/err.txt
less experiments/${experiment}/stdout.txt
"
    pereggrn --experiment_name $experiment --amount_to_do missing_models \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt
done