#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=qc_sample_mind
#SBATCH -c 1
#SBATCH -o logs/o_qc_sample_mind.txt
#SBATCH -e logs/e_qc_sample_mind.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# Load PLINK module if required by your cluster infrastructure
# module load plink/1.90b

plink --bfile $FILENAME \
      --mind 0.05 \
      --make-bed \
      --out ${FILENAME}_after_call_rate

# Export list of removed individuals
mv ${FILENAME}_after_call_rate.irem CALL_RATE_OUTLIERS.txt

echo "**** Job ends: $(date) ****"