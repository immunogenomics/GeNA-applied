#!/bin/bash

# This script reformats the per-SNP sample-level phenotype file output from GeNA into
# PLINK(1) format phenotype files per lead SNP. Then it calls PLINK to generate fixed-
# phenotype summary statistics for each csaQTL.

for celltype in "NK" "Myeloid"
do
    infile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/spheno.tsv"
    outfile_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/spheno_"

    command="python -u reformat_spheno.py --infile $infile --outfile_path $outfile_path"
    echo $command
    eval $command
done

module load plink

celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    # get chromosome and position of lead snp                                                                                      
    IFS=':'
    read -a snp_vals <<<"$lead_snp"
    chr=${snp_vals[0]}
    lead_pos=${snp_vals[1]}
    IFS=' '

    trait_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/spheno_${lead_snp}.tsv"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/coloc/gwas_local_${celltype}_${lead_snp}_cis"

    command="plink --bfile ${gtypes} --pheno ${trait_file} --assoc --prune --allow-no-sex --out ${outfile}"
    echo $command    
    eval $command
done