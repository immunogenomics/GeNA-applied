#!/bin/bash

# This script applies GeNA to null (permuted) simulated genotypes to test for csaQTLs in order to quantify
# false association rates. GeNA is applied to these null genotypes for each of five real single-cell objects, 
# each corresponding to a different denominator of included cells (T, B, NK, Myeloid or all cells)

path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

# Run GeNA with default parameters
for celltype in "NK" "T" "B" "allcells" "Myeloid"
do
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
    covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'

    for MAF_group in "MAF5" "contMAF" "MAF1"
    do
	res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/${celltype}/${MAF_group}/"
	mkdir -p $res_folder

	gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_${MAF_group}.perm_${celltype}"
	
	command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs"
	echo $command
	eval $command
    done
done

# Specify value of k in order to assess variation in calibration with k value
for MAF_group in "MAF5" "contMAF"
do
    for celltype in "NK" "Myeloid" "T" "B" "allcells"
    do
	sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
	covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
	res_folder_stem="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/${celltype}/${MAF_group}/"
	gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_${MAF_group}.perm_${celltype}"
	
	k_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/${celltype}/"
	
	for k_type in "small_k" "large_k"
	do
	    ks_file="${k_folder}${k_type}.txt"
	    res_folder="${res_folder_stem}${k_type}/"
	    mkdir -p $res_folder
	    
	    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs -k $ks_file"
	    echo $command
	    eval $command
	done
    done
done