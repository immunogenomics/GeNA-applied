# In this repository
We provide code and analyses supporting the `GeNA` manuscript. We provide scripts to:
- Apply `GeNA` to real single-cell profiling and simulated genotypes to evaluate `GeNA`'s calibration (type I error, see `null` folder) and statsitical power (type II error, see `nonnull_sims` folder)
- Apply `GeNA` to identify cell state abundance QTLs (csaQTLs) in the OneK1K dataset (see `run_gwas` and `leadsnps_perm` folders)
- Test associations to each lead SNP in single-cell objects with cis-genes removed (see `mask_cis_genes`)
- Perform GWAS of cluster-based cell type proportion traits for comparison (see `cluster_gwas`)
- Evaluate the replication by `GeNA` of csaQTLs previously identified using flow cytometry
- Evaluate the replication of csaQTLs from the OneK1K discovery cohort in five replication cohorts
- Examine cell state abundance associations to polygenic risk scores using CNA

# Citation
The `GeNA` manuscript can be found and cited at
[link to preprint]

# To use GeNA
Please refer to the `GeNA` repository at [immnogenomics/GeNA](https://github.com/immunogenomics/GeNA)

# Data availability
All datasets used in these analyses are previously published:
1. Yazar, S. et al. Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease. Science 376, eabf3041 (2022).
2. Perez, R. K. et al. Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. Science (American Association for the Advancement of Science) 376, eabf1970–eabf1970 (2022).
3. Oelen, R. et al. Single-cell RNA-sequencing of peripheral blood mononuclear cells reveals widespread, context-specific gene expression regulation upon pathogenic exposure. Nat Commun 13, 3267 (2022).
4. Randolph, H. E. et al. Genetic ancestry effects on the response to viral infection are pervasive but cell type specific. Science 374, 1127–1133 (2021).

# Contact
Please contact Laurie Rumker (Laurie_Rumker@hms.harvard.edu) with any questions about these analyses.
