#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=qc_sample_het
#SBATCH -c 1
#SBATCH -o logs/o_qc_sample_het.txt
#SBATCH -e logs/e_qc_sample_het.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# Prune for linkage disequilibrium first to avoid biased clustering
plink --bfile $FILENAME --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile $FILENAME --extract pruning.prune.in --make-bed --out pruned_data

# Calculate heterozygosity statistics
plink --bfile pruned_data --het --out prunedHet

# Extract outliers outside the range [-0.15, 0.15]
awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers_low.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers_high.txt
cat outliers_low.txt outliers_high.txt > HETEROZYGOSITY_OUTLIERS.txt

# Remove outliers from main dataset
cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt
plink --bfile $FILENAME \
      --remove all_outliers.txt \
      --make-bed \
      --out ${FILENAME}_after_heterozyg

echo "**** Job ends: $(date) ****"