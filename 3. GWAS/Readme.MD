# Phase 3: Association Testing & Covariate Regression

Following genotype imputation, association testing evaluates whether an allele's dosage covaries significantly with the target phenotype across parallel genomic models.  This page outlines the procedures for generating population structure covariates via `FlashPCA` and running the association engine via `RVTESTS`.  

## 👥 1. Latent Covariate Generation

To prevent false-positive signals driven by population stratification, genome-wide association models must systematically control for ancestry using Principal Components (PCs), typically alongside biological covariates like `AGE` and `SEX`.


### A. Long-Range LD Region Inversion Exclusions

Prior to computing PCs, highly polymorphic regions prone to long-range linkage disequilibrium (LD)—such as the Major Histocompatibility Complex (MHC) on chromosome 6 or inversion loops on chromosomes 5 and 8—must be strictly excluded. This ensures the calculated PCs capture true macro-level macro-geographic stratification rather than localized, long-range LD blocks driven by independent selective sweeps. 

- **Recommended hg38 Exclusions**:  

```
Chr1     47534328    51534328    r1
Chr2    133742429   137242430    r2
Chr2    182135273   189135274    r3
Chr3     47458510    49962567    r4
Chr3     83450849    86950850    r5
Chr5     98664296   101164296    r6
Chr5    129664307   132664308    r7
Chr5    136164311   139164311    r8
Chr6     24999772    35032223    r9
Chr6    139678863   142178863    r10
Chr8      7142478    13142491    r11
Chr8    110987771   113987771    r12
Chr11   87789108    90766832    r13
Chr12   109062195   111562196    r14
Chr20   33412194    35912078    r15

```

### B. LD Pruning and PC Calculation
Using the raw dataset, variants are strictly filtered, pruned for pairwise linkage disequilibrium, and processed with `FlashPCA` to extract the top 5–10 eigenvectors.

---

## 🛠️ 2. Phenotype Design for Three-Group Cohorts

Because we are modeling a three-group study design (**Healthy Controls, MDD Patients, and BP Patients**), standard binary case-control regressions are insufficient to capture the full spectrum of risk. To ensure appropriate statistical modeling, we must approach this using one of two strategies:

### Strategy A: Pairwise Direct Group Contrasts (Recommended)

Split the association phase into independent downstream models using the `--peopleIncludeFile` flag to subset our cohorts. This evaluates clear biological boundaries without violating regression assumptions:

1. **Model 1 (MDD vs. Controls):** Isolates unipolar depression mechanisms.
2. **Model 2 (BP vs. Controls):** Isolates bipolar disorder mechanisms.
3. **Model 3 (BP vs. MDD):** Directly maps the genomic differences differentiating the two distinct psychiatric conditions.

### Strategy B: Multi-Categorical or Polytomous Regression

If modeling all three groups within a single continuous framework, we must convert the phenotype column into dummy-coded categorical indicators or leverage an LMM architecture capable of processing multi-allelic multinomial outcomes.

## 💻 Script Description & SLURM Automation

The automation architecture for this step is divided into two separate pipeline scripts:

- `gwas_step3_pca.sh` : (**Ancestry PC Computation**) ; The automated high-performance computing batch scripts to compute covariates and run the association models.

- `gwas_step3_rvtests.sh` : (**Parallelized Association Testing**) ; This script uses a SLURM job array to run association testing across chromosomes 1–22 concurrently, pulling the raw dosage data and testing it using a Wald single-variant model.

### 📊 Parameters Guide

- `--dosage DS`: Standardizes inputs to extract continuous probability dosage data directly from the input VCF arrays rather than hard-called genotypes.

- `--rangeFile`: References a regional positional matrix filter (compiled in Phase 1) to bypass low-frequency or poorly imputed variant coordinates, speeding up processing.

- `--single wald`: Evaluates fixed-effects coefficients and significance levels using standard single-variant Wald testing.

## 📚 References

- **FlashPCA Matrix**: Purcell, S. et al. (2007) 'PLINK: a tool set for whole-genome association and population-based linkage analyses', The American Journal of Human Genetics, 81(3), pp. 559–575. doi:10.1086/519795.

- **RVTESTS Engine**: Zhan, X. et al. (2016) 'RVtests: an efficient and comprehensive tool for rare variant association analysis using sequence data', Bioinformatics, 32(9), pp. 1423–1426. doi:10.1093/bioinformatics/btw079.


