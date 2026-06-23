#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH --job-name=qc_sample_gcta
#SBATCH -c 4
#SBATCH -o logs/o_qc_sample_gcta.txt
#SBATCH -e logs/e_qc_sample_gcta.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# Load GCTA module if needed
# module load gcta

# Generate global Genetic Relationship Matrix (GRM) using multi-threading
gcta --bfile $FILENAME \
     --make-grm \
     --out GRM_matrix \
     --autosome \
     --maf 0.05 \
     --thread-num 4

# Apply relatedness threshold pruning (0.125 removes cousins)
gcta --grm-cutoff 0.125 \
     --grm GRM_matrix \
     --make-grm \
     --out GRM_matrix_0125 \
     --thread-num 4

# Subset the PLINK dataset to keep only unrelated individuals
plink --bfile $FILENAME \
      --keep GRM_matrix_0125.grm.id \
      --make-bed \
      --out ${FILENAME}_relatedness \
      --threads 4

echo "**** Job ends: $(date) ****"