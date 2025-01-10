#!/bin/bash

#SBATCH --job-name=compare_eigengenes
#SBATCH --output=compare_eigengenes.out
#SBATCH --error=compare_eigengenes.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=bcf

start=$(date +%s)

source /homes/mkessler/no_backup/miniconda3/etc/profile.d/conda.sh
conda activate RanOmics

python ../scripts/compare_eigengenes.py

end=$(date +%s) 
elapsed=$(( end - start )) 

echo "Runtime: $elapsed seconds" >> compare_eigengenes.out
