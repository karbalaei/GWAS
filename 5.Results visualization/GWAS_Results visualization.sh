#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH --job-name=gwas_phase5_vis
#SBATCH -c 4
#SBATCH -o logs/o_gwas_phase5_vis.txt
#SBATCH -e logs/e_gwas_phase5_vis.txt
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts: $(date) ****"

# -------------------------------------------------------------------------
# CONFIGURATION PARAMETERS
# -------------------------------------------------------------------------
# Set the specific pairwise contrast being analyzed to dynamically adjust sample sizes
# Options: "MDD_vs_Controls", "BP_vs_Controls", "BP_vs_MDD"
CONTRAST="MDD_vs_Controls"

INPUT_GWAS="data/phase4_meta/gwas_${CONTRAST}_meta_results_SORTED.tbl"
OUTPUT_DIR="analysis/plots/${CONTRAST}"
mkdir -p ${OUTPUT_DIR} logs

# Supply sample counts to compute standardized lambda1000 metrics
if [ "${CONTRAST}" == "MDD_vs_Controls" ]; then
    N_CASES=588     # MDD samples
    N_CONTROLS=540  # Control samples as example reference
elif [ "${CONTRAST}" == "BP_vs_Controls" ]; then
    N_CASES=503     # BP samples
    N_CONTROLS=540
elif [ "${CONTRAST}" == "BP_vs_MDD" ]; then
    N_CASES=503
    N_CONTROLS=588
fi

# Load required cluster software dependencies
module load conda_R/4.4

# -------------------------------------------------------------------------
# EXECUTE DIAGNOSTIC RENDERING ENGINE VIA R
# -------------------------------------------------------------------------
echo ">> Processing Diagnostic Visualization for Contrast: ${CONTRAST}..."
echo ">> Active Sample Matrices -> Cases: ${N_CASES} | Controls: ${N_CONTROLS}"

Rscript - <<EOF
library(data.table)
library(ggplot2)
library(dplyr)

# Speed up file read IO
setDTthreads(4)

# Load summary statistics data frame
gwas <- fread("${INPUT_GWAS}")

# Standardize typical input column maps (adjust strings if mapping raw RVTESTS vs METAL)
# For this script, we assume columns named 'Chromosome', 'Position', 'Pvalue' (or 'P-value')
if ("P-value" %in% colnames(gwas)) setnames(gwas, "P-value", "Pvalue")
if ("MarkerName" %in% colnames(gwas)) {
    # If using METAL default output format, extract components
    gwas[, c("CHR_str", "BP_str") := tstrsplit(MarkerName, ":", fixed=TRUE)]
    gwas\$Chromosome <- as.numeric(gsub("chr", "", gwas\$CHR_str))
    gwas\$Position <- as.numeric(gwas\$BP_str)
}

# Remove structural NA or null entries to prevent downstream layout failures
gwas_clean <- gwas[!is.na(Pvalue) & Pvalue > 0 & Pvalue <= 1]

p_vector <- gwas_clean\$Pvalue
n_variants <- length(p_vector)

# -------------------------------------------------------------------------
# CALCULATION: LAMBDA & STANDARDIZED LAMBDA 1000
# -------------------------------------------------------------------------
message(">> Computing Genomic Inflation Factor parameters...")
x2obs <- qchisq(p_vector, df = 1, lower.tail = FALSE)
x2exp <- qchisq((1:n_variants) / n_variants, df = 1, lower.tail = FALSE)

lambda <- median(x2obs, na.rm = TRUE) / median(x2exp, na.rm = TRUE)
n_cases <- ${N_CASES}
n_controls <- ${N_CONTROLS}

# Normalize unbalanced numbers to an idealized cohort of 1000 cases and 1000 controls
lambda1000 <- 1 + (lambda - 1) * (1/n_cases + 1/n_controls) / (1/1000 + 1/1000)

cat(paste0(">>> RESULTS FOR ", "${CONTRAST}", " <<<\n"))
cat(paste0("Observed Lambda (Raw): ", round(lambda, 4), "\n"))
cat(paste0("Standardized Lambda (1000): ", round(lambda1000, 4), "\n"))

# Save values to a quick metrics log file
writeLines(
    c(paste0("Contrast: ", "${CONTRAST}"), paste0("Lambda_Raw: ", lambda), paste0("Lambda_1000: ", lambda1000)),
    con = file.path("${OUTPUT_DIR}", "inflation_metrics.txt")
)

# -------------------------------------------------------------------------
# GENERATION: QUANTILE-QUANTILE (Q-Q) PLOT
# -------------------------------------------------------------------------
message(">> Generating high-resolution Q-Q vector display...")
observed_logP <- -log10(sort(p_vector))
expected_logP <- -log10((1:n_variants) / (n_variants + 1))

qq_df <- data.frame(Expected = expected_logP, Observed = observed_logP)

# Downsample non-significant background variants (P > 0.01) to keep vector image size lightweight
qq_sig <- qq_df %>% filter(Observed >= 2)
qq_nonsig <- qq_df %>% filter(Observed < 2) %>% sample_frac(0.05)
qq_plot_data <- rbind(qq_sig, qq_nonsig)

qq_plot <- ggplot(qq_plot_data, aes(x = Expected, y = Observed)) +
    geom_point(color = "black", alpha = 0.5, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
    labs(
        title = paste0("Q-Q Plot: ", "${CONTRAST}"),
        subtitle = paste0("Raw Lambda = ", round(lambda, 3), " | Lambda1000 = ", round(lambda1000, 3)),
        x = "Expected (-log10 P-value)",
        y = "Observed (-log10 P-value)"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))

ggsave(
    filename = file.path("${OUTPUT_DIR}", "gwas_qqplot.pdf"),
    plot = qq_plot, width = 7, height = 7, units = "in", dpi = 300
)

# -------------------------------------------------------------------------
# GENERATION: MANHATTAN PLOT (FAST BASE-R ARCHITECTURE)
# -------------------------------------------------------------------------
message(">> Generating genome-wide Manhattan display...")
pdf(file.path("${OUTPUT_DIR}", "gwas_manhattan.pdf"), width = 11, height = 6)

# Select a subset of variants for fast Manhattan plotting to avoid RAM bottlenecks
manhattan_data <- gwas_clean %>% 
    filter(Pvalue < 0.05) %>% 
    select(Chromosome, Position, Pvalue) %>% 
    arrange(Chromosome, Position)

with(manhattan_data, {
    # Generate continuous cumulative base-pair indexes across chromosomes
    chr_lengths <- tapply(Position, Chromosome, max)
    chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
    names(chr_offsets) <- names(chr_lengths)
    
    plot_x <- Position + chr_offsets[as.character(Chromosome)]
    plot_y <- -log10(Pvalue)
    
    # Define alternating colors per chromosome
    color_palette <- rep(c("navyblue", "skyblue3"), 12)
    plot_colors <- color_palette[Chromosome]
    
    # Calculate chromosome label tick marks
    tick_positions <- sapply(names(chr_lengths), function(ch) {
        mean(range(plot_x[Chromosome == as.numeric(ch)]))
    })
    
    plot(
        plot_x, plot_y, col = plot_colors, pch = 20, cex = 0.6,
        xlab = "Chromosome", ylab = "-log10(P-value)",
        main = paste0("Genome-Wide Association Signal Matrix: ", "${CONTRAST}"),
        xaxt = "n", las = 1, ylim = c(0, max(c(max(plot_y), 8)))
    )
    
    axis(1, at = tick_positions, labels = names(chr_lengths), cex.axis = 0.8, las = 2)
    
    # Overlay standard genome-wide alpha lines
    abline(h = -log10(5e-8), col = "red", linetype = "solid", lwd = 1.5)  # Suggestive/Significant line
    abline(h = -log10(1e-5), col = "blue", linetype = "dashed", lwd = 1)   # Discovery line
})

dev.off()
message("Downstream diagnostic plots successfully rendered.")
EOF

echo "**** Job ends: $(date) ****"

