#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH --job-name=gwas_imputation_filter
#SBATCH -c 4
#SBATCH -o logs/o_gwas_imputation_filter_%a.txt
#SBATCH -e logs/e_gwas_imputation_filter_%a.txt
#SBATCH --array=1-22%4
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Task id (Chromosome): ${SLURM_ARRAY_TASK_ID}"

# Define the target chromosome from the SLURM array index
CHR=${SLURM_ARRAY_TASK_ID}

# -------------------------------------------------------------------------
# CONFIGURATION PARAMETERS
# -------------------------------------------------------------------------
# Set your preferred Rsq (INFO) metric threshold for filtering
# Choose: 0.30 (lenient), 0.60 (standard), or 0.80 (high-stringency/TWAS prep)
RSQ_THRESHOLD=0.80

# Define your directory paths
RAW_IMPUTED_DIR="data/imputed_raw"
FILTERED_OUTPUT_DIR="data/imputed_filtered_rsq${RSQ_THRESHOLD}"
mkdir -p ${FILTERED_OUTPUT_DIR}

# Input file naming conventions
INFO_FILE="${RAW_IMPUTED_DIR}/chr${CHR}.info.gz"
VCF_FILE="${RAW_IMPUTED_DIR}/chr${CHR}.dose.vcf.gz"

# Output file naming conventions
RSQ_PASSED_SNPS="${FILTERED_OUTPUT_DIR}/chr${CHR}_rsq_passed_snplist.txt"
CLEAN_VCF_OUT="${FILTERED_OUTPUT_DIR}/chr${CHR}.dose.filtered.vcf.gz"

# -------------------------------------------------------------------------
# ENVIRONMENT MODULES
# -------------------------------------------------------------------------
# Load required software tools on your cluster ecosystem
module load vcftools
# module load bcftools

# -------------------------------------------------------------------------
# STEP 1: PARSE INFO FILE & EXTRACT VARIANTS MEETING RSQ CRITERIA
# -------------------------------------------------------------------------
echo ">> Step 1: Filtering variants in ${INFO_FILE} with Rsq > ${RSQ_THRESHOLD}..."

# The .info.gz file contains columns: SNP, ALT_Frq, MAF, Rsq, Genotyped, etc.
# This blocks extracts the 'SNP' column (col 1) if the 'Rsq' value (col 7) is greater than the threshold
zcat ${INFO_FILE} | awk -v thresh="${RSQ_THRESHOLD}" '
    NR == 1 {print "Header skipped"}; 
    NR > 1 { if ($7 != "-" && $7 > thresh) print $1 }
' > ${RSQ_PASSED_SNPS}

NUM_PASSED=$(wc -l < ${RSQ_PASSED_SNPS})
echo ">> Total variants passing stringency on Chromosome ${CHR}: ${NUM_PASSED}"

# -------------------------------------------------------------------------
# STEP 2: FILTER DOSAGE VCF FILES BASED ON EXTRACTED INVARIANT SET
# -------------------------------------------------------------------------
echo ">> Step 2: Extracting variants from VCF and compressing output..."

# Filter the dosage VCF using vcftools and output directly to a zipped stream
vcftools --gzvcf ${VCF_FILE} \
         --snps ${RSQ_PASSED_SNPS} \
         --recode \
         --stdout | bgzip -c -@ 4 > ${CLEAN_VCF_OUT}

# -------------------------------------------------------------------------
# STEP 3: INDEX THE NEWLY FILTERED DOSAGE VCF FILE
# -------------------------------------------------------------------------
echo ">> Step 3: Indexing filtered VCF file using tabix..."

# Generate the .tbi index file required by Phase 3 association engines
tabix -p vcf ${CLEAN_VCF_OUT}

echo ">> Process completed successfully for Chromosome ${CHR}."
echo "**** Job ends: $(date) ****"