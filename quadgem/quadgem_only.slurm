#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:30:00
#SBATCH --job-name=quadgem
#SBATCH --output=quadgem_%A_91rCo_290t3a_%a.out
#SBATCH --array=1-8
# #SBATCH --reservation=AtRiskUnstableEnvironment
##### Uncomment for debugging 
# #SBATCH --partition=debug
# #SBATCH --constraint=skylake

module purge
module restore gfpp
cd /home/zoua1/gem_simulations/quadgem

./gem_part ${SLURM_ARRAY_TASK_ID} 4gemf_290dv
