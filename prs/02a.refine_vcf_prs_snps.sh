#!/bin/bash

module load samtools
module load bcftools

for scorename in "RA" "SLE"
do
    in_filename="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.all_chr.phased.mm3.imputed.dose.snpQC_MAF5.renamed.vcf.gz"
    snplist="/data/srlab/lrumker/MCSC_Project/cna-prs/sumstats/snplists/${scorename}_snplist.txt"
    out_filename="/data/srlab/lrumker/MCSC_Project/cna-prs/results/geno/${scorename}.all_chr.vcf.gz"
    
    command="zcat $in_filename | bcftools view --include ID==@${snplist} -Oz -o $out_filename"  
    echo $command
    eval $command
done

