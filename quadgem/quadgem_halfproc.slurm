#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --job-name=quadgem
#SBATCH --output=quadgem_%A_task_%a.out
#SBATCH --array=1
##### Uncomment for debugging 
# #SBATCH --partition=debug
# #SBATCH --constraint=skylake

module purge
module restore gfpp
cd /home/zoua1/gem_simulations/quadgem

./rdo_part 050421
./rdo_part_rect 050421