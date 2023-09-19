#!/bin/bash

# This script runs GeNA on the vcf containing permuted versions of genotype values per sample 
# at each lead SNP

path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/leadsnps_perm/results/"

lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
celltypes=( "Myeloid" "NK" "NK" "NK" "NK" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    snp=${lead_snps[i]}
    celltype=${celltypes[i]}
    
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
    gtypes="${res_folder}${snp}.perm"

    outfolder="${res_folder}${celltype}_${snp}"
    mkdir -p $outfolder
    outfolder="${outfolder}/"

    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $outfolder -c $covs"
    echo $command
    eval $command

done