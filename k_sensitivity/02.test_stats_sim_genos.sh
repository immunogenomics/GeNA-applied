#!/bin/bash

# This script runs GeNA in a csaQTL gwas for each cell type, using the simulated genotypes with true associations
# Unlike our primary power analyses, this script applies a modified version of GeNA. The modifications are: 
# 1) storage of more intermediate objects and 2) no multiple testing correction across the total values of k considered
# Further, in contrast to our primary power analyses where default GeNA parameters were used, here we use non-default 
# values of k, forcing GeNA to consider all k between 1 and 50

pwd=$(eval "pwd")
path_to_GeNA="mod_GeNA/"
cd $path_to_GeNA

covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
k_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA-applied/k_sensitivity/large_k.txt"

for celltype in "Myeloid" "NK" "T" "B" "allcells"
do
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
    res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/k_sensitivity/gwas_${celltype}_nonnull_sims/"    

    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/sim_genotypes/${celltype}_sim_genotypes"

    mkdir -p $res_folder
        
    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs -k ${k_file}"
    echo $command
    eval $command
done