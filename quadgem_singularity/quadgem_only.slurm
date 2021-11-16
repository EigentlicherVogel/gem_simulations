#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --job-name=quadgem
#SBATCH --output=quadgem_%A_73rCo_350dv_pionb2t_%a.out
#SBATCH --array=1-8
# #SBATCH --reservation=AtRiskUnstableEnvironment
##### Uncomment for debugging 
# #SBATCH --partition=debug
# #SBATCH --constraint=skylake

module purge
cd /home/zoua1/gem_simulations/sim_newver

singularity exec -i /home/zoua1/garfieldpp.20210814.sif ./shellrunjob_g.sh ${SLURM_ARRAY_TASK_ID}
