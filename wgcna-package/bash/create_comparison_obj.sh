#!/bin/bash

#SBATCH --job-name=WGCNA_comp
#SBATCH --output=WGCNA_comp.out
#SBATCH --error=WGCNA_comp.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=bcf

start=$(date +%s)

source /homes/mkessler/no_backup/miniconda3/etc/profile.d/conda.sh
conda activate RanOmics

python ../scripts/create_comparison_obj.py

end=$(date +%s) 
elapsed=$(( end - start )) 

echo "Runtime: $elapsed seconds" >> WGCNA_comp.out
