#!/bin/bash

# Calls a script to compute and store, per SNP in the SLE PRS: 
# 1) neighborhood-level phenotype values
# 2) sample-level phenotype values per individual

cd /data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA-applied/run_gwas/ # where define_loci.py script is located

celltype="Myeloid"
outfile_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/sle_prs/"
sc_object="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/sle_prs/${celltype}_wgeno.h5ad" # includes genotypes for PRS SNPs in d.samplem
loci="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/sle_prs/sle_prs_snps_sumstats.tsv" #GeNA sumstats for PRS SNPs
covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'

command="python -u define_phenotype.py \
                 --loci $loci \
                 --covs $covs \
                 --outfile_path $outfile_path \
                 --sc_object $sc_object"
echo $command
eval $command
