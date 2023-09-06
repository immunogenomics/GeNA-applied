#!/bin/bash

# Calls a script to compute and store, per SNP: 
# 1) neighborhood-level phenotype values
# 2) sample-level phenotype values per individual

for celltype in "Myeloid" "NK"
do
    outfile_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/"
    sc_object="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/${celltype}_wgeno.h5ad"
    loci="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/gwas_loci.tsv" #GeNA sumstats
    covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
	
    command="python -u define_phenotype.py \
                 --loci $loci \
                 --covs $covs \
                 --outfile_path $outfile_path \
                 --sc_object $sc_object"
    echo $command
    eval $command
done
