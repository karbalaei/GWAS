## 🛠️ Phase 1: Quality Control & Data Filtering

Upstream data cleaning is the most critical determinant of downstream imputation accuracy. Genotyping arrays naturally introduce systematic errors, probe hybridization failures, and sample contamination that must be statistically pruned prior to reference panel matching. This pipeline utilizes **PLINK (v1.9)** and **GCTA** to execute rigorous sample- and variant-level filtering.
### 👤 1. Sample-Level Filtering & Ancestry Uniformity

Sample-level QC identifies and removes low-quality DNA samples, contaminated specimens, chromosomal abnormalities, and population outliers.

**A. Call Rate Outliers ** (--mind)

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


* **Ancestry Homogeneity Tracking:** Infers genetic ancestry by projecting sample genotypes against a reference dataset (such as HapMap3 or 1000 Genomes) using Principal Component Analysis (PCA). Samples falling beyond $\pm 6$ standard deviations from the target cluster mean are excluded to prevent false-positive signals driven by population stratification.
* **Cryptic Relatedness Pruning:** Evaluates a Genetic Relationship Matrix (GRM) to calculate identity-by-descent (IBD) coefficients. A threshold is set at a genomic relationship cutoff of $0.125$ to systematically eliminate first- and second-degree cousins, ensuring sample independence for standard linear regression modeling. Alternatively, the cutoff can be adjusted upwards (e.g., $0.80$) solely to catch duplicate samples or monozygotic twins, enabling the downstream use of a Linear Mixed Model (LMM) that accounts for family structures.

### 🧬 Variant-Level Quality Control Parameters
* **Variant Call Rate (Missingness):** Filters out specific Single Nucleotide Polymorphisms (SNPs) exhibiting high missingness coefficients. High variant missingness indicates poor probe hybridization chemistry or low-quality genomic regions, introducing experimental noise.
* **Differential Missingness by Phenotype Status:** Compares variant call rates between case and control cohorts using a specialized test. Any systematic deviation in missingness indicates technical batch effects or technical artifacts during sample plate processing rather than true biological variance.
* **Haplotype-Specific Missingness:** Tests whether the missingness of a specific variant is dependent on its surrounding flanking haplotype. This catches localized sequencing drops or directional alignment biases.
* **Hardy-Weinberg Equilibrium (HWE) Deviations:** Assesses deviations from HWE within the control cohort exclusively. Severe deviations from historical HWE baseline assumptions signal systematic genotyping errors, copying duplication artifacts, or severe probe misalignments rather than selective evolutionary pressures.
* **Minor Allele Frequency (MAF) Thresholds:** Imposes a strict lower limit on allele frequencies to eliminate ultra-rare variants that lack the statistical power to resolve true linear associations. This filter is applied conditionally based on whether the downstream array design is explicitly optimized for rare-variant discovery.

### 📋 Post-QC Imputation Preparation
Following variant filtration, the cleaned dataset undergoes strand alignment synchronization. Using specialized strand-checking utilities, variant alleles, genomic coordinates, and reference frequencies are matched against the reference assembly to ensure consistency before submission to remote imputation servers.


# 1. QC and data cleaning
The QC and data cleaning is very important prior imputation since the cleaner data you put in there the better and more accurate data you will get out. PLINK and GCTA are very easy program for doing this.

## Sample QC parameters


- Heterozygosity outliers (--het), F cut-off of -0.15 and <- 0.15 for inclusion 

High or low heterozygosity is an indication of failed experiments or contamination of DNA sample

```
plink --bfile $FILENAME --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile $FILENAME --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet
awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > HETEROZYGOSITY_OUTLIERS.txt
cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt
plink --bfile $FILENAME --remove all_outliers.txt --make-bed --out $FILENAME_after_heterozyg
```

- Call rate outliers (--mind), Call rate of >95% is preferred which is --mind 0.05

Low call rate is an indication of a failed experiment

```
plink --bfile $FILENAME --mind 0.05 --make-bed --out $FILENAME_after_call_rate
mv $FILENAME_after_call_rate.irem CALL_RATE_OUTLIERS.txt
```

- Genetic sex fails (--check-sex), F cut-off of 0.25 and 0.75 for inclusion --check-sex 0.25 0.75

Genetic sex fails are indication of a failed experiment or a sample-switch

Note use 2,781,479 and 155,701,383 for hg38

```
plink --bfile $FILENAME --check-sex 0.25 0.75 --maf 0.05 --out gender_check1
plink --bfile $FILENAME --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out gender_check2 
grep "PROBLEM" gender_check1.sexcheck > problems1.txt
grep "PROBLEM" gender_check2.sexcheck > problems2.txt
cat problems1.txt problems2.txt > GENDER_FAILURES.txt
cut -f 1,2 GENDER_FAILURES.txt > samples_to_remove.txt
plink --bfile $FILENAME --remove samples_to_remove.txt --make-bed --out $FILENAME_after_gender
```

## optional steps depending general purpose of data cleaning
- No ancestry outliers -> based on Hapmap3 PCA data, should be near combined CEU/TSI, 6SD+/-

```
see script hapmap ancestry check folder this
```

- No relatedness (--grm-cutoff), Typically this would be set at 0.125 to remove cousins or more related individuals

Note you can also increase the grm-cut-off to e.g. 0.8 and remove duplicate samples. This would allow you keep more related individuals and use a linear mixed model.

```
gcta --bfile $FILENAME --make-grm --out GRM_matrix --autosome --maf 0.05 
gcta --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0125 --make-grm
plink --bfile $FILENAME --keep GRM_matrix_0125.grm.id --make-bed --out $FILENAME_relatedness
```

## Variant QC parameters

- Variant missingness, basically if a variant has high missingness this suggests that the genotyping probe for this variant is not great 

```
plink --bfile $FILENAME --make-bed --out $FILENAME_geno --geno 0.05
```


- Missingness by case control (--test-missing), using P > 1E-4

```
plink --bfile $FILENAME --test-missing --out missing_snps 
awk '{if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt
plink --bfile $FILENAME --exclude missing_snps_1E4.txt --make-bed --out $FILENAME_missing1
```

- Missing by haplotype (--test-mishap), using P > 1E-4
```
plink --bfile $FILENAME --test-mishap --out missing_hap 
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\n/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt
plink --bfile $FILENAME --exclude missing_haps_1E4_final.txt --make-bed --out $FILENAME_missing2
```


- Hardy Weinberg SNP from controls only (--filter-controls --hwe 1E-4)
```
plink --bfile $FILENAME --filter-controls --hwe 1E-4 --write-snplist
plink --bfile $FILENAME --extract plink.snplist --make-bed --out $FILENAME_HWE
```

# optional steps depending general purpose of data cleaning 
- Minor allele frequency (--maf), In some instances this might not be used if specific rare variants are on your array that you do not want to exclude.

```
plink --bfile $FILENAME --maf 0.01 --make-bed --out $FILENAME_MAF
```

## Preparation of data prior to submission to Imputation server

Prepare data for HRC imputation using Will Rayner tool and continue with the plink file that is cleaned with the above described options.


```
# download file to check
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim.v4.2.5.zip

# make your .frq file
plink --bfile $FILENAME --freq --out $FILENAME

perl HRC-1000G-check-bim.pl -b $FILENAME.bim -f $FILENAME.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

# then run to fix your data
sh Run-plink.sh

# then make vcf files

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  plink --bfile YOURFILE-updated-chr$chnum --recode vcf --chr $chnum --out YOURFILE$chnum 
done

## then sort and zip

module load vcftools

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
	vcf-sort YOURFILE$chnum.vcf | bgzip -c >  pre_impute_YOURFILE_$chnum.vcf.gz
done

# and then you are ready to submit to the imputation server
```


