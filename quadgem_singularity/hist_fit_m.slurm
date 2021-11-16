#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=14:00:00
#SBATCH --job-name=quadgem
#SBATCH --output=quadgem_%A_73ArCo_ROXY_45D_300C_20um_%a.out
#SBATCH --array=1-8
# #SBATCH --reservation=AtRiskUnstableEnvironment
##### Uncomment for debugging 
# #SBATCH --partition=debug
# #SBATCH --constraint=skylake

module purge
cd /home/zoua1/gem_simulations/sim_newver

singularity exec -i /home/zoua1/garfieldpp.20210814.sif ./shellrunjob.sh ${SLURM_ARRAY_TASK_ID} 

