#!/bin/bash

module load samtools
module load bcftools

chr=22
# Generate VCFs for only the selected snps at MAF 5% and all-MAF 
for MAF_group in "MAF5" "contMAF"
do
    snplist="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/sel_snps_${MAF_group}.txt"
    infile="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.chr${chr}.phased.mm3.imputed.dose.snpQC_MAF5.renamed.vcf.gz"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_${MAF_group}.vcf.gz"
    command="zcat ${infile} | bcftools view --include ID==@${snplist} -Oz -o ${outfile}"
    echo $command
    eval $command
done

# Generate VCF for only the selected snps at MAF 1%
snplist="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/sel_snps_MAF1.txt"
infile="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.chr${chr}.phased.mm3.imputed.dose.snpQC.renamed.vcf.gz"
outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_MAF1.vcf.gz"

command="zcat ${infile} | bcftools view --include ID==@${snplist} -Oz -o ${outfile}"
echo $command
eval $command