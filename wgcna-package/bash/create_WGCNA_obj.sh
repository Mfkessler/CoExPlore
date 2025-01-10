#!/bin/bash

#SBATCH --job-name=WGCNA
#SBATCH --output=WGCNA.out
#SBATCH --error=WGCNA.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --partition=bcf

start=$(date +%s)

source /homes/mkessler/no_backup/miniconda3/etc/profile.d/conda.sh
conda activate RanOmics2

python ../scripts/create_WGCNA_obj.py

end=$(date +%s) 
elapsed=$(( end - start )) 

echo "Runtime: $elapsed seconds" >> WGCNA.out
