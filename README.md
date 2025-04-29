# In this repository
We provide code and analyses supporting the `GeNA` publication. We provide scripts to:
- Apply `GeNA` to real single-cell profiling and simulated genotypes to evaluate `GeNA`'s calibration (`null` folder) and statistical power (`nonnull_sims` folder)
- Apply `GeNA` to identify cell state abundance QTLs (csaQTLs) in the OneK1K dataset (`run_gwas`, `leadsnps_perm`, `suggestive_loci`, `retest_subcohorts`, `molecularQTLs` folders)
- Test associations to each lead SNP in single-cell objects with cis-genes removed (`mask_cis_genes` folder) or suggestive trans-eGenes removed (`mask_trans_eGenes` folder)
- Perform GWAS of cluster-based cell type proportion traits for comparison (`cluster_gwas` folder)
- Evaluate the replication of csaQTLs from the OneK1K discovery cohort in five replication cohorts (`replication` folder)
- Evaluate the replication of csaQTLs previously identified using flow cytometry in our neighborhood-based framework for single-cell data (`published_csaQTLs` folder)
- Examine cell state abundance associations to polygenic risk scores (`prs`)
- Evaluate the sensitivity of our results to various aspects of the primary analysis (`ccg_retained`, `k_sensitivity`, `conditional_testing` folders)
- Apply `GeNA` to a dataset of cells in early neural differentiation (`neural_dset` folder) 

We also provide the `notebooks` used to generate figures and key reported values.

# Citation
If you use `GeNA` in your work, you can cite our paper [here](https://www.nature.com/articles/s41588-024-01909-1) 

# Dependencies
This work was completed with GeNA version v1.0.0, which has the following dependencies:
- Python version 3.8.10
- R version 4.1.1
- PLINK version 2.00a2.3
- CNA version 0.1.6
- Rmpfr version 0.8-7
Scripts in this repo will not run smoothly with later versions of GeNA, which we updated to maintain compatability with new input/output formatting in CNA.

# To use GeNA
Please refer to the `GeNA` repository at [immnogenomics/GeNA](https://github.com/immunogenomics/GeNA)

# Data availability
All datasets used in these analyses are previously published:
1. Yazar, S. et al. Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease. Science 376, eabf3041 (2022).
2. Perez, R. K. et al. Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. Science (American Association for the Advancement of Science) 376, eabf1970–eabf1970 (2022).
3. Oelen, R. et al. Single-cell RNA-sequencing of peripheral blood mononuclear cells reveals widespread, context-specific gene expression regulation upon pathogenic exposure. Nat Commun 13, 3267 (2022).
4. Randolph, H. E. et al. Genetic ancestry effects on the response to viral infection are pervasive but cell type specific. Science 374, 1127–1133 (2021).
5. Jerber, J. et al. Population-scale single-cell RNA-seq profiling across dopaminergic neuron differentiation. Nat. Genet. 53, 304–312 (2021).

# Contact
Please contact Laurie Rumker (Laurie_Rumker AT hms.harvard.edu) with any questions about these analyses.
