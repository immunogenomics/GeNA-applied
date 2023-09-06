#!/bin/bash

# This script subsets the GeNA summary statistics file to a cis window (+/- 1MB)
# around each lead SNP.

celltypes=( "Myeloid" "NK" "NK" "NK" "NK" "NK")
lead_snps=( "15:80267501:A:G" "2:111851212:C:T" "5:161884695:G:T" "11:128070535:A:G" "12:10583611:C:T" "19:16442019:G:A")
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    # get chromosome and position of lead snp
    IFS=':'
    read -a snp_vals <<< "$lead_snp"
    chr=${snp_vals[0]}
    lead_pos=${snp_vals[1]}
    IFS=' '

    # define 2MB cis window around lead snp
    MB=1000000 # 1 megabase
    start_pos=$(($lead_pos-$MB))
    end_pos=$(($lead_pos+$MB))

    infile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/GeNA_sumstats.txt"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/gwas_res_${lead_snp}_cis.tsv"
    command="awk 'NR==1 || \$1==$chr' $infile | awk 'NR==1 || \$2<$end_pos' | awk 'NR==1 || \$2>$start_pos' > $outfile"
    echo $command
    eval $command
    
done