#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=9:00:00
#SBATCH --job-name=quadgem
#SBATCH --output=quadgem_%A_73ArCo_350dv_4rot90_%a.out
#SBATCH --array=1-5
# #SBATCH --reservation=AtRiskUnstableEnvironment
##### Uncomment for debugging 
# #SBATCH --partition=debug
# #SBATCH --constraint=skylake

module purge
cd /home/zoua1/gem_simulations/sim_newver

echo "slurm script check"
singularity exec -i /home/zoua1/garfieldpp.20210814.sif ./shellrunjob_i.sh ${SLURM_ARRAY_TASK_ID} 

