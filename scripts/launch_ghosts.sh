#!/bin/bash
#SBATCH --job-name=ghosts
#SBATCH --qos=gp_debug
#SBATCH --account=bsc40

#SBATCH --output=../joberrout/%x_%j.out
#SBATCH --error=../joberrout/%x_%j.err

#SBATCH --ntasks=112
#SBATCH --time=02:00:00

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

module load R

mkdir -p ../outputs
python3 01_get_branch_set.py
Rscript 02_simulate_transfers.R
Rscript 03_plotting_shifts.R
Rscript 04_groups_data.R

