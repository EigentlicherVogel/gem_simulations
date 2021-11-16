#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=05:30:00
#SBATCH --job-name=quadgem
#SBATCH --output=quadgem_%A_91rCo_300t2a_%a.out
#SBATCH --array=1-8
# #SBATCH --reservation=AtRiskUnstableEnvironment
##### Uncomment for debugging 
# #SBATCH --partition=debug
# #SBATCH --constraint=skylake

module purge
cd /home/zoua1/gem_simulations/quadgem

TEMPVARIA = $(printf "source /home/zoua1/setupgpp.sh && cd gem_simulations/qgem_newver/build && ./quadgem %s 4gemf_300dv" ${SLURM_ARRAY_TASK_ID})
singularity exec -i /home/zoua1/garfieldpp.20210428.sif bash -c TEMPVARIA
