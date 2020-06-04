#!/bin/bash
#SBATCH --partition owners,normal
#SBATCH --cpus-per-task 4
#SBATCH --time 00:30:00
#SBATCH --export ALL
#SBATCH --mail-type FAIL
#SBATCH --mail-user etchin@stanford.edu

date 

ml R/3.6

Rscript microsims_dynamic_worker_only.R 20 2 1 0.5 ${SLURM_ARRAY_TASK_ID}

date
