#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH --job-name=gwas_phase4_meta
#SBATCH -c 4
#SBATCH -o logs/o_gwas_phase4_meta.txt
#SBATCH -e logs/e_gwas_phase4_meta.txt
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# -------------------------------------------------------------------------
# CONFIGURATION PARAMETERS
# -------------------------------------------------------------------------
# Define the specific pairwise contrast context to process 
# Options: MDD_vs_Controls, BP_vs_Controls, or BP_vs_MDD
CONTRAST="MDD_vs_Controls"

COHORT_LIST=("cohort1" "cohort2" "cohort3")
RAW_ASSOC_DIR="data/phase3_association"
HARMONIZED_DIR="data/phase4_meta_input"
mkdir -p ${HARMONIZED_DIR} logs

# Load required cluster dependencies
module load conda_R/4.4
# module load metal

# -------------------------------------------------------------------------
# STEP 1: CONCATENATE CHROMOSOME ASSOC AND INFO FILES PER COHORT
# -------------------------------------------------------------------------
for COHORT in "${COHORT_LIST[@]}"; do
    echo ">> Consolidating genome-wide outputs for: ${COHORT} (${CONTRAST})..."
    
    # Concatenate all regional SingleWald association outputs, stripping internal header rows
    cat ${RAW_ASSOC_DIR}/${COHORT}_${CONTRAST}.chr*.SingleWald.assoc | grep -v 'N_INFORMATIVE' > ${HARMONIZED_DIR}/${COHORT}_${CONTRAST}_allChrs.assoc || true
    
    # Concatenate all regional imputation info files, stripping internal header rows
    cat ${RAW_ASSOC_DIR}/${COHORT}_maf001rsq03minimums_chr*.info | grep -v 'Rsq' > ${HARMONIZED_DIR}/${COHORT}_allChrs.Info || true
done

# -------------------------------------------------------------------------
# STEP 2: HARMONIZE SUMMARY STATISTICS VIA HIGH-SPEED R BACKEND
# -------------------------------------------------------------------------
echo ">> Launching R script to execute allele orientation and data cleaning..."

Rscript - <<EOF
library(data.table)
library(dplyr)

setDTthreads(4)

cohort_list <- c("${COHORT_LIST[@]}")
contrast <- "${CONTRAST}"

for (cohort in cohort_list) {
    message(paste("Processing integration matrices for:", cohort))
    
    info_path <- file.path("${HARMONIZED_DIR}", paste0(cohort, "_allChrs.Info"))
    assoc_path <- file.path("${HARMONIZED_DIR}", paste0(cohort, "_", contrast, "_allChrs.assoc"))
    output_path <- file.path("${HARMONIZED_DIR}", paste0(cohort, "_", contrast, "_harmonized.tab"))
    
    # Check for empty files before processing
    if (file.size(info_path) == 0 || file.size(assoc_path) == 0) {
        stop(paste("Error: Missing concatenated input files for cohort:", cohort))
    }
    
    # Load data via data.table for optimal performance
    infos <- fread(info_path, header = FALSE, col.names = c("SNP", "ALT_Frq", "Rsq"))
    assoc <- fread(assoc_path, header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "N_INFORMATIVE", "Test", "Beta", "SE", "Pvalue"))
    
    # Combine datasets on unique variant keys
    data <- merge(infos, assoc, by.x = "SNP", by.y = "Test", all.y = TRUE)
    
    # Filter out unrealistic effect sizes and missing parameters
    dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(Pvalue))
    
    # Standardize chromosome and variant nomenclature
    dat\$chr <- paste0("chr", dat\$CHROM)
    dat\$markerID <- paste(dat\$chr, dat\$POS, sep = ":")
    
    # Systematically orient effects relative to the Minor Allele
    dat\$minorAllele <- ifelse(dat\$ALT_Frq <= 0.5, as.character(dat\$ALT), as.character(dat\$REF))
    dat\$majorAllele <- ifelse(dat\$ALT_Frq <= 0.5, as.character(dat\$REF), as.character(dat\$ALT))
    dat\$beta        <- ifelse(dat\$ALT_Frq <= 0.5, dat\$Beta, dat\$Beta * -1)
    dat\$se          <- dat\$SE
    dat\$maf         <- ifelse(dat\$ALT_Frq <= 0.5, dat\$ALT_Frq, 1 - dat\$ALT_Frq)
    dat\$P           <- dat\$Pvalue
    dat\$N           <- dat\$N_INFORMATIVE
    
    # Select clean final output schema matrix
    dat0 <- dat[, c("markerID", "minorAllele", "majorAllele", "beta", "se", "maf", "P", "N")]
    
    # Write tab-separated clean files
    fwrite(dat0, file = output_path, quote = FALSE, sep = "\t")
}
message("Allele orientation and structural data cleaning successfully completed.")
EOF

# -------------------------------------------------------------------------
# STEP 3: RUN FIXED-EFFECTS META-ANALYSIS VIA METAL
# -------------------------------------------------------------------------
echo ">> Running METAL meta-analysis..."

# Execute the engine using your custom configuration scheme
metal meta_analysis_scheme.txt

# -------------------------------------------------------------------------
# STEP 4: DIAGNOSTIC SORTING & FILTERING
# -------------------------------------------------------------------------
echo ">> Running post-meta-analysis diagnostic filtering..."

# The system automatically sorts the final results table by the P-value column (Column 10)
META_OUTPUT="gwas_pairwise_contrast_meta_results1.tbl"

if [ -f "${META_OUTPUT}" ]; then
    sort -gk 10 ${META_OUTPUT} > gwas_pairwise_contrast_meta_results_SORTED.tbl
    echo ">> Final sorted meta-analysis file compiled successfully."
else
    echo ">> Warning: Final output file not found. Check your cohort alignment configurations."
fi

echo "**** Job ends: $(date) ****"