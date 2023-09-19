#!/bin/bash

# This script generates one vcf file per lead SNP containing nperm permutations of the observed
# genotypes values for that lead SNP across all samples

nperm=1000000 #1e-06
n_header=15

lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
celltypes=( "Myeloid" "NK" "NK" "NK" "NK" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    snp=${lead_snps[i]}
    celltype=${celltypes[i]}

    nampc_file="/data/srlab/lrumker/datasets/onek1k/pheno/nampcs/${celltype}_batch_covs.csv"
    
    header_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/leadsnps_perm/results/${snp}_vcf_header.txt"
    infile="/data/srlab/lrumker/MCSC_Project/cna-qtl/leadsnps_perm/results/${snp}.vcf.gz"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/leadsnps_perm/results/${snp}.perm.vcf"
    
    command="python -u perm_dose.py \
                 --outfile $outfile \
                 --nperm $nperm \
                 --n_header $n_header \
                 --nampc_file $nampc_file \
                 --infile $infile"
    echo $command
    eval $command

    command="zcat $infile | head -n ${n_header} > ${header_file}"
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
