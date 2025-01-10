#!/bin/bash

#SBATCH --job-name=create_adatas
#SBATCH --output=create_adatas.out
#SBATCH --error=create_adatas.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --partition=bcf

start=$(date +%s)

source /homes/mkessler/no_backup/miniconda3/etc/profile.d/conda.sh
conda activate RanOmics2

python ../scripts/create_adatas.py

end=$(date +%s) 
elapsed=$(( end - start )) 

echo "Runtime: $elapsed seconds" >> create_adatas.out
