# 🛠️ Phase 1: Quality Control & Data Filtering

Upstream data cleaning is the most critical determinant of downstream imputation accuracy. Genotyping arrays naturally introduce systematic errors, probe hybridization failures, and sample contamination that must be statistically pruned prior to reference panel matching. This pipeline utilizes **PLINK (v1.9)** and **GCTA** to execute rigorous sample- and variant-level filtering.
## 👤 1. Sample-Level Filtering & Ancestry Uniformity

Sample-level QC identifies and removes low-quality DNA samples, contaminated specimens, chromosomal abnormalities, and population outliers.

**A. Call Rate Outliers** (--mind)

Samples with low call rates indicate a failed genotyping experiment (e.g., poor DNA quality or low concentration). A strict threshold of $>95\%$ call rate is enforced.
```
plink --bfile $FILENAME \
      --mind 0.05 \
      --make-bed \
      --out ${FILENAME}_after_call_rate
# Export list of removed individuals
mv ${FILENAME}_after_call_rate.irem CALL_RATE_OUTLIERS.txt
```

**B. Heterozygosity Outliers** (--het)

Severe deviations in mean heterozygosity indicate either sample contamination (inflated heterozygosity) or high inbreeding coefficients/poor sample preparation (deflated heterozygosity). Samples outside $\pm 0.15$ of the target inbreeding coefficient ($F$) distribution are removed.

```
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
```

**C. Genetic Sex Discordance Checks** (--check-sex)

Identifies sample swaps or chromosomal anomalies by comparing cataloged sex metadata against calculated X-chromosome homozygosity ($F$-coefficients).

- Note for hg38: Ensure coordinates match the hg38 pseudo-autosomal regions (PAR) boundary limits: 2,781,479 to 155,701,383.

```
# Run sex check on non-PAR regions of Chromosome 23 (X)
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
```

**D. Cryptic Relatedness & Kinship Pruning** (GCTA)

To satisfy the baseline sample-independence assumptions of standard linear/logistic association models, close genetic relatives must be systematically pruned. A genomic relationship matrix (GRM) cutoff of $0.125$ is used to eliminate second-degree cousins or closer.

- LMM Note:* If your downstream configuration relies on a Linear Mixed Model (LMM), you can increase this cutoff to $0.80$ strictly to eliminate duplicate samples while retaining related families, as LMMs account for kinship matrices natively.
```
# Generate global Genetic Relationship Matrix (GRM)
gcta --bfile $FILENAME \
     --make-grm \
     --out GRM_matrix \
     --autosome \
     --maf 0.05 

# Apply relatedness threshold pruning (0.125 removes cousins)
gcta --grm-cutoff 0.125 \
     --grm GRM_matrix \
     --make-grm \
     --out GRM_matrix_0125 

# Subset the PLINK dataset to keep only unrelated individuals
plink --bfile $FILENAME \
      --keep GRM_matrix_0125.grm.id \
      --make-bed \
      --out ${FILENAME}_relatedness
```

**E. Ancestry Outlier Pruning**

Prevents false-positive association signals driven by population stratification. Genotypes are projected against a reference panel (such as HapMap3) using Principal Component Analysis (PCA), and samples falling beyond $\pm 6$ standard deviations from the combined CEU/TSI population mean cluster are excluded.(Executed via the dedicated script inside your local hapmap ancestry check directory).

##🧬 2. Variant-Level Quality Control

Variant-level QC filters out poorly performing genotyping probes, alignment failures, and markers that exhibit skewed inheritance patterns.

**A. Variant Missingness** (--geno)

Filters out individual SNPs with high missing data scores, which indicate bad array probe hybridization chemistry.

```
# Generate global Genetic Relationship Matrix (GRM)
plink --bfile $FILENAME \
      --geno 0.05 \
      --make-bed \
      --out ${FILENAME}_geno
```


**B. Phenotype-Differential Missingness** (--test-missing)

Tests whether variant missingness correlates with phenotypic diagnosis. Significant differences point to systematic technical anomalies during plate processing or call biases rather than biological effects.

- **Three-Group Cohort Adjustment**: Because PLINK’s --test-missing command evaluates binary phenotypes (Case vs. Control), standard three-group factors (Controls, MDD, BP) must be split or combined into binary comparisons. The most common approach is to group both patient categories together as "All Cases" vs. Controls to screen for technical artifacts between healthy and diseased sample tracks.

```
# Step 1: Group both patient tracks (MDD=1, BP=1) vs Controls (0) in your phenotype mapping
plink --bfile $FILENAME \
      --test-missing \
      --out missing_snps 

# Step 2: Prune variants showing differential missingness at P <= 1E-4
awk '{if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt
plink --bfile $FILENAME \
      --exclude missing_snps_1E4.txt \
      --make-bed \
      --out ${FILENAME}_missing1
```

**C. Haplotype-Specific Missingness** (--test-mishap)

Investigates whether the missingness of a variant depends on its surrounding haplotype structure, catching subtle directional processing biases.

```
plink --bfile $FILENAME \
      --test-mishap \
      --out missing_hap 

# Parse and format the rejected variant identifiers
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\n/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt

# Remove failed variants
plink --bfile $FILENAME \
      --exclude missing_haps_1E4_final.txt \
      --make-bed \
      --out ${FILENAME}_missing2
```

**D. Hardy-Weinberg Equilibrium (HWE) in Controls** (--hwe)

Severe deviations from HWE typically signal systemic genotyping errors, duplicating artifacts, or copy number expansions. To prevent true disease risk variants from being aggressively filtered out, this evaluation is explicitly restricted to **healthy controls only**.

```
# Calculate HWE P-values across healthy controls exclusively at a P < 1E-4 cutoff
plink --bfile $FILENAME \
      --filter-controls \
      --hwe 1E-4 \
      --write-snplist

# Extract the validated HWE invariant set
plink --bfile $FILENAME \
      --extract plink.snplist \
      --make-bed \
      --out ${FILENAME}_HWE
```

**E. Minor Allele Frequency** (--maf)

Removes variants with ultra-low frequencies to ensure adequate statistical power during downstream regression modeling.

- Note: This step may be bypassed or set lower if your array contains customized exome content or rare risk variants that you intentionally wish to screen.

```
plink --bfile $FILENAME \
      --maf 0.01 \
      --make-bed \
      --out ${FILENAME}_MAF

```

## 📋 3. Post-QC Imputation Preparation

Before uploading the final datasets to the Michigan Imputation Server, the strand directions, variant names, and coordinates must be accurately checked and converted into chromosome-specific VCF formats.

```
# 1. Download Will Rayner's strand and alignment checking tool
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim.v4.2.5.zip

# 2. Generate a frequency file from your cleaned post-QC dataset
plink --bfile $FILENAME --freq --out $FILENAME

# 3. Check alleles, strands, positions, and frequencies against the HRC panel reference
perl HRC-1000G-check-bim.pl \
     -b $FILENAME.bim \
     -f $FILENAME.frq \
     -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
     -h

# 4. Execute the automated fix script generated by the tool to flip strands/fix alleles
sh Run-plink.sh

# 5. Convert fixed datasets into chromosome-specific VCF format (Chromosomes 1-23)
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}; do
  plink --bfile YOURFILE-updated-chr$chnum \
        --recode vcf \
        --chr $chnum \
        --out YOURFILE_chr$chnum 
done

# 6. Sort and zip the files using vcftools/bgzip for server submission
module load vcftools

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}; do
  vcf-sort YOURFILE_chr$chnum.vcf | bgzip -c > pre_impute_YOURFILE_chr$chnum.vcf.gz
done

# The system is now optimized and ready for upload to the Michigan Imputation Server!
```