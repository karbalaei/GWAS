#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=qc_sample_sex
#SBATCH -c 1
#SBATCH -o logs/o_qc_sample_sex.txt
#SBATCH -e logs/e_qc_sample_sex.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# Run sex check on non-PAR regions of Chromosome 23 (X) for hg38
plink --bfile $FILENAME \
      --chr 23 \
      --from-bp 2781479 \
      --to-bp 155701383 \
      --maf 0.05 \
      --geno 0.05 \
      --hwe 1E-5 \
      --check-sex 0.25 0.75 \
      --out gender_check

# Isolate misidentified individuals
grep "PROBLEM" gender_check.sexcheck > GENDER_FAILURES.txt
cut -f 1,2 GENDER_FAILURES.txt > sex_samples_to_remove.txt

# Remove failed samples
plink --bfile $FILENAME \
      --remove sex_samples_to_remove.txt \
      --make-bed \
      --out ${FILENAME}_after_gender

echo "**** Job ends: $(date) ****"