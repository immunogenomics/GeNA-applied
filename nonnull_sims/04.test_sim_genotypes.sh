#!/bin/bash

# This script runs GeNA for each genotype simulate within the corresponding single-cell data object

path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/"

for celltype in "Myeloid" "NK" "T" "B" "allcells"
do
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
    gtypes="${res_folder}sim_genotypes/${celltype}_sim_genotypes"

    outfolder="${res_folder}results/${celltype}/"
    mkdir -p $outfolder
    	
    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $outfolder -c $covs"
    echo $command
    eval $command

done