#!/bin/bash

# Calls a script to compute and store, per suggestive locus lead SNP: 
# 1) neighborhood-level phenotype values
# 2) sample-level phenotype values per individual

# Here we use a different version of define_phenotype.py than in the primary OneK1K GWAS *only* in that this 
# version does not incorporate a batch effect

for celltype in "HipSci"
do
    outfile_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/"
    sc_object="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/${celltype}_wgeno.h5ad" # includes genotypes for lead SNPs in d.samplem
    loci="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/gwas_loci.tsv" #GeNA sumstats for lead SNPs
    covs='gPC1,gPC2,gPC3,gPC4,gPC5'
	
    command="python -u define_phenotype.py \
                 --loci $loci \
                 --covs $covs \
                 --outfile_path $outfile_path \
                 --sc_object $sc_object"
    echo $command
    eval $command
done
