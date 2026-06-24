# Phase 5: Downstream Diagnostic Visualization Protocol

Following the calculation of cross-study meta-analysis summary statistics or single-site regression outputs, global association distributions must be strictly evaluated. Visualizing results via Manhattan and Quantile-Quantile (Q-Q) plots, alongside calculating the Genomic Inflation Factor ($\lambda$), serves as a critical diagnostic gate to verify that observed associations represent true polygenic architecture rather than systematic statistical inflation.

Because a three-group design (**Healthy Controls, MDD Patients, and BP Patients**) is modeled via independent pairwise direct group contrasts (e.g., MDD vs. Controls, BP vs. Controls, or BP vs. MDD), **diagnostic visualization metrics must be executed independently per summary statistic contrast file**.

---

## 📊 1. Diagnostic Metrics and Principles

### A. Genome-Wide Association Distributions (Manhattan Plot)

Manhattan plots display the negative log-transformed association $P$-values ($-\log_{10}(P)$) across chronological chromosomal coordinates. This layout identifies sharp, localized genomic signals ("peaks") rising above the standard genome-wide significance threshold ($P < 5 \times 10^{-8}$).

* **Implementation:** Highly optimized frameworks such as the command-line standalone engine `Manhattan-Plotter` or R-based vector graphics suites are utilized to render publication-ready high-density displays.

### B. Quantile-Quantile (Q-Q) Layouts

Q-Q plots chart the observed distribution of p-values against the uniform distribution expected under a true global null hypothesis. Deviations at the extreme tail reveal true polygenic traits, whereas early systemic upward deviations alert researchers to systemic errors, unmodeled population stratification, or technical batch artifacts.

### C. Genomic Inflation Assessment ($\lambda$ and $\lambda_{1000}$)

* **$\lambda$ (Genomic Inflation Factor):** The ratio of the median observed $\chi^2$ statistic to the median expected $\chi^2$ statistic under the null. A $\lambda \approx 1.0$ flags a clean, well-controlled model; values exceeding $1.05$ require closer evaluation.
* **$\lambda_{1000}$ (Sample-Size Standardization):** Because genomic inflation scales mathematically with sample size, the standard $\lambda$ metric can appear artificially inflated in high-powered discovery designs. The inflation factor is standardized to a hypothetical cohort of 1,000 cases and 1,000 controls using the formula below to cleanly separate true polygenic background from structural artifacts:

$$\lambda_{1000} = 1 + (\lambda - 1) \times \frac{\frac{1}{N_{\text{cases}}} + \frac{1}{N_{\text{controls}}}}{\frac{1}{1000} + \frac{1}{1000}}$$

---

## 💻 2. SLURM Automation Pipeline Script

This robust, production-grade SLURM batch script handles data parsing via fast `data.table` backends, automatically evaluates $\lambda$ and standardized $\lambda_{1000}$ coefficients based on your dynamic case-control allocations, and prints vectorized diagnostic plots per pairwise group context.


---

## 📚 References

* **Manhattan Plotting Standards:** Turner, S. D. (2018) 'qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots', *Journal of Open Source Software*, 3(25), p. 731. doi:10.21105/joss.00731.
* **Genomic Inflation Calculations:** Devlin, B. and Roeder, K. (1999) 'Genomic control for association studies', *Biometrics*, 55(4), pp. 997–1004. doi:10.1111/j.0006-341x.1999.00997.x.
* **Lambda1000 Standardization Metric:** Freedman, M. L. et al. (2004) 'Assessing the heritability of complex traits in stratified populations', *Nature Genetics*, 36(4), pp. 388–393. doi:10.1038/ng1331.