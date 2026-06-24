# Phase 2: Genotype Imputation Protocol

Following variant filtering and strand harmonization in Phase 1, this phase utilizes the remote Michigan Imputation Server framework to statistically estimate unobserved genotypes[cite: 43]. By leveraging dense, deep-sequencing reference panels, imputation massively increases the density of testable genetic markers, allowing for fine-mapped downstream association testing[cite: 43].

---

## 🌐 Imputation Reference Architecture

* **Imputation Engine:** Michigan Imputation Server hosted by the University of Michigan[cite: 43].
* **Haplotype Reference Panel:** Haplotype Reference Consortium (HRC) Release 1.1 (2016), optimized to maximize haplotype coverage across European ancestries[cite: 43].
* **Pre-Phasing Algorithm:** Eagle (v2.3), used to computationally estimate local haplotype phases prior to multi-variant statistical imputation[cite: 43].

---

## ⚖️ Phenotype-Blind Robustness for Three-Group Cohorts

Because genotype imputation calculates alternative allele probabilities based purely on background flanking haplotype matches within reference panels, the imputation process is entirely phenotype-blind[cite: 43]. This makes the protocol highly robust for multi-group study designs (**healthy controls, MDD patients, and BP patients**), as it estimates variant dosages objectively across all chromosomes without introducing diagnostic case-control bias[cite: 43].

---

## 📦 Imputation Server Outputs

Upon completion of the imputation pipeline, the server generates three specific compressed files per chromosome (e.g., for Chromosome 21)[cite: 43]:

* `chr21.dose.vcf.gz`: The primary data file containing imputed genotype configurations stored as continuous probabilistic alternative allele counts, known as dosages[cite: 43].
* `chr21.dose.vcf.gz.tbi`: The tabix index file required by association engines for high-throughput, localized genomic lookups[cite: 43].
* `chr21.info.gz`: A detailed summary file tracking the diagnostic quality, marker names, and alternative frequencies across all imputed or array-genotyped variants[cite: 43].

---

## 📊 Evaluation Metrics & Downstream Filtering Strategy

The `.info.gz` diagnostic file tracks several vital properties per variant to guide downstream quality control[cite: 43]:


```
SNP      REF(0)  ALT(1)  ALT_Frq  MAF      AvgCall  Rsq      Genotyped  LooRsq  EmpR  EmpRsq  Dose0  Dose1
9:11739  G       A       0.00004  0.00004  0.99996  0.00023  Imputed    -       -     -       -      -
9:14665  G       A       0.20555  0.20555  0.82608  0.20321  Imputed    -       -     -       -      -
9:62179  G       A       0.22655  0.22655  0.99994  0.99976  Genotyped  0.998   0.999  0.99717 0.99920 0.00078
```

### 🔑 Key Informational Columns

* `SNP`: Variant identifier formatted as `Chromosome:Basepair` coordinates[cite: 43].
* `REF(0)` / `ALT(1)`: The reference allele (0) and the alternative/effect allele (1)[cite: 43].
* `ALT_Frq` / `MAF`: The calculated alternative allele frequency and the corresponding Minor Allele Frequency within the sample population[cite: 43].
* `Genotyped`: Indicates whether the marker was directly verified by the upstream array chip (`Genotyped`) or computationally inferred (`Imputed`)[cite: 43].
* `Rsq` (**Imputation Quality Score**): Measures the statistical precision of the estimated allele dosages[cite: 43].

---

## 🛑 Downstream Quality Control Filter Thresholds

The `Rsq` score acts as a filter to keep unreliable variants from causing false positives during Phase 3 association mapping[cite: 43]. Depending on the specific polygenic architecture and stringency required by the study design, variants should be pruned using one of the following standard thresholds[cite: 43]:

* `Rsq > 0.30`: Lenient threshold, commonly utilized to maximize genome-wide variant discovery for highly polygenic traits[cite: 43].
* `Rsq > 0.60`: Standard publication-grade threshold balance[cite: 43].
* `Rsq > 0.80`: High-stringency threshold, strictly required if down-stream fine-mapping or transcriptomic weight computation (e.g., TWAS/isoTWAS pipelines) is planned[cite: 43].

---

## 💻 Script Description & SLURM Automation

Since the imputation itself takes place on the remote Michigan Imputation Server web GUI, the downstream filtering script executes the crucial next step: parsing the server's `.info.gz` files, extracting high-quality variants based on your chosen $R^2$ stringency threshold, and filtering the corresponding dosage VCF files so they are optimized for your Phase 3 multi-group association tests[cite: 43].

### 💡 Why the Script is Optimized for This Workflow:

* **SLURM Job Array Architecture (`--array=1-22%4`):** It runs across chromosomes 1 to 22 concurrently, throttling to process a maximum of 4 chromosomes at a time to prevent resource over-allocation on the shared cluster nodes[cite: 43].
* **Phenotype-Blind Alignment:** Since the sample tracks consist of three distinct groups (Controls, MDD, BP), filtering by a unified `Rsq` index across the chromosome ensures that no artificial allele imbalances are introduced between diagnostic groups prior to testing[cite: 43].
* **Automated Clean-Up Verification:** Generates standard `.vcf.gz` and accompanying `.tbi` structures ready to be directly imported into `RVTESTS`, `PLINK`, or complex Linear Mixed Models (LMM)[cite: 43].

---

## 📚 References

* **Michigan Imputation Server:** Das, S. et al. (2016) 'Next-generation genotype imputation service and methods', *Nature Genetics*, 48(10), pp. 1284–1287. doi:10.1038/ng.3656[cite: 43].
* **Haplotype Reference Consortium (HRC):** McCarthy, S. et al. (2016) 'A reference panel of 64,976 haplotypes for genotype imputation', *Nature Genetics*, 48(10), pp. 1279–1283. doi:10.1038/ng.3643[cite: 43].
* **Eagle Phasing Framework:** Loh, P. R. et al. (2016) 'Reference-based genotype phasing across all chromosomes with Eagle2', *Nature Genetics*, 48(11), pp. 1443–1448. doi:10.1038/ng.3679[cite: 43].


