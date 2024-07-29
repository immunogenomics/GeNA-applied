#!/bin/bash

# This script applies GeNA to generate test statistics per NAM-PC for larger values of k than GeNA 
# considers by default for all genome-wide significant loci from the OneK1K csaQTL GWAS
# We again use the modified version of GeNA in order to save the necessary additional intermediate 
# objects and to eliminate multiple testing correction across the 50 values of k

pwd=$(eval "pwd")
path_to_GeNA="mod_GeNA/"
cd $path_to_GeNA

celltypes=( "Myeloid" "NK" "NK" "NK" "NK" )
lead_snps=( "15:80263217:C:T" "11:128070535:A:G" "2:111851212:C:T" "12:10583611:C:T" "19:16441973:G:A" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
    res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/k_sensitivity/gwas_${celltype}/${lead_snp}/"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis"
    covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
    k_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA-applied/k_sensitivity/large_k.txt"

    mkdir -p ${res_folder}
    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs -k ${k_file}"
    echo $command
    eval $command
done