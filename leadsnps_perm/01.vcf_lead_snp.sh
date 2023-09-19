#!/bin/bash

# This script generares one vcf file per lead SNP containing only genotype information for that SNP

module load samtools
module load bcftools

lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
chrs=( "15" "2" "11" "12" "19" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    snp=${lead_snps[i]}
    chr=${chrs[i]}

    echo "$snp" > tempfile.txt

    infile="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.chr${chr}.phased.mm3.imputed.dose.snpQC_MAF5.renamed.vcf.gz"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/leadsnps_perm/results/${snp}.vcf.gz" 
    
    command="zcat ${infile} | bcftools view --include ID==@tempfile.txt -Oz -o ${outfile}"                                      
    echo $command                                                                                                                    
    eval $command

    rm tempfile.txt

done