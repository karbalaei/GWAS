#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=15G
#SBATCH --job-name=qc_variant_missing
#SBATCH -c 1
#SBATCH -o logs/o_qc_variant_missing.txt
#SBATCH -e logs/e_qc_variant_missing.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# Basic variant missingness thresholding
plink --bfile $FILENAME \
      --geno 0.05 \
      --make-bed \
      --out ${FILENAME}_geno

# Group patient tracks (MDD=1, BP=1) vs Controls (0) in phenotype mapping prior to testing
plink --bfile ${FILENAME}_geno \
      --test-missing \
      --out missing_snps 

# Prune variants showing differential missingness at P <= 1E-4
awk '{if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt

plink --bfile ${FILENAME}_geno \
      --exclude missing_snps_1E4.txt \
      --make-bed \
      --out ${FILENAME}_missing1

echo "**** Job ends: $(date) ****"