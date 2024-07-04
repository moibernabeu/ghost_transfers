#!/bin/bash
#SBATCH --job-name=ghosts
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc40

#SBATCH --output=../joberrout/%x_%j.out
#SBATCH --error=../joberrout/%x_%j.err

#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00

module load R

mkdir ../outputs
Rscript 00_parse_tree.R
python3 01_get_branch_set.py
Rscript 02_simulate_transfers.R
Rscript 03_plotting_shifts.R
Rscript 04_groups_data.R
Rscript 05_heatmap.R

