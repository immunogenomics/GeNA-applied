#!/bin/bash

# This script generates simulated genotypes with the specified MAF and 
# true associations to the input cell state abundance traits. 
# To introduce noise and vary the amount of variance explained in the trait
# by the simulated genotype, for each phenotype we permuted genotype values for 
# 1%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80% and 100% of samples. For each 
# of these ten tiers of sample count to permute, the script will generate 'nsim' 
# genotype permutations (simulates). 

nsim=100
maf=0.25

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/"
for celltype in "T" "NK" "allcells" "Myeloid" "B"
do
    trait_file="${res_folder}traits/${celltype}_traits.txt"
    header_file="${res_folder}vcf_header.txt"
    outfile="${res_folder}sim_genotypes/${celltype}_sim_genotypes.vcf"
    outfile_meta="${res_folder}sim_genotypes/${celltype}_sim_genotypes_meta.tsv"
    
    command="python -u sim_genotypes.py --trait_file $trait_file \
                 --outfile_genos $outfile --outfile_meta $outfile_meta --nsim $nsim  --celltype $celltype --maf $maf"
    echo $command
    eval $command
    
    command="cat ${header_file} ${outfile} > ${outfile}_temp"
    echo $command
    eval $command
    
    command="mv ${outfile}_temp ${outfile}"
    echo $command
    eval $command
    
    command="gzip ${outfile}"
    echo $command
    eval $command
done
