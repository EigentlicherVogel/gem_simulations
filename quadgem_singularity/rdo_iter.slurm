#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=quadgem
#SBATCH --output=rdo_%A_73rCo_350dv_50pion_rdmdt_15bftr2_60it_%a.out
#SBATCH --array=1-8
# #SBATCH --reservation=AtRiskUnstableEnvironment
##### Uncomment for debugging 
# #SBATCH --partition=debug
# #SBATCH --constraint=skylake

module purge
cd /home/zoua1/gem_simulations/sim_newver

singularity exec -i /home/zoua1/garfieldpp.20210814.sif ./shellrunjob.sh ${SLURM_ARRAY_TASK_ID} 

