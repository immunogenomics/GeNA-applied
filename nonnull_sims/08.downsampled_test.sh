#!/bin/bash

# This script runs GeNA in a csaQTL gwas for each cell type, using the simulated genotypes with true associations

path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/"

for suffix in "0pt8" "0pt6" "0pt4" "0pt2"
do
    for celltype in "Myeloid" "NK" "B" "T" "allcells"
    do
	sc_object="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/downsampled/${celltype}_${suffix}.h5ad"
	gtypes="${res_folder}sim_genotypes/${celltype}_${suffix}_sim_genotypes"
	
	outfolder="${res_folder}results/${celltype}/downsampled_${suffix}/"
	mkdir -p $outfolder
    	
	command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $outfolder -c $covs"
	echo $command
	eval $command
    done
done