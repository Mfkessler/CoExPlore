#!/bin/bash

#SBATCH --job-name=adata_analysis
#SBATCH --output=adata_analysis.out
#SBATCH --error=adata_analysis.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=bcf

start=$(date +%s)

source /homes/mkessler/no_backup/miniconda3/etc/profile.d/conda.sh
conda activate RanOmics2

python ../scripts/analyse_adata.py

end=$(date +%s) 
elapsed=$(( end - start )) 

echo "Runtime: $elapsed seconds" >> adata_analysis.out
