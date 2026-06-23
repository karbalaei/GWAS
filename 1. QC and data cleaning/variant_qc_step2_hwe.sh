#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=15G
#SBATCH --job-name=qc_variant_hwe
#SBATCH -c 1
#SBATCH -o logs/o_qc_variant_hwe.txt
#SBATCH -e logs/e_qc_variant_hwe.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# Evaluate haplotype-specific missingness parameters
plink --bfile ${FILENAME}_missing1 \
      --test-mishap \
      --out missing_hap 

# Parse and format the rejected variant identifiers
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\n/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt

# Filter out haplotype outliers
plink --bfile ${FILENAME}_missing1 \
      --exclude missing_haps_1E4_final.txt \
      --make-bed \
      --out ${FILENAME}_missing2

# Calculate HWE P-values across healthy controls exclusively at a P < 1E-4 cutoff
plink --bfile ${FILENAME}_missing2 \
      --filter-controls \
      --hwe 1E-4 \
      --write-snplist

# Extract the validated HWE invariant set and filter on MAF
plink --bfile ${FILENAME}_missing2 \
      --extract plink.snplist \
      --maf 0.01 \
      --make-bed \
      --out ${FILENAME}_final_clean

echo "**** Job ends: $(date) ****"