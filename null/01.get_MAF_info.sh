#!/bin/bash

module load samtools
module load bcftools

chr=22

# Pull MAF information for all candidate snps (already passed MAF 5% threshold) on chr 22 
in_filename="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.chr${chr}.phased.mm3.imputed.dose.snpQC_MAF5.renamed.vcf.gz"
out_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/mafs/1K1K.chr${chr}.info.txt"

command="zcat $in_filename | bcftools query -f '%ID\t%MAF\n' > $out_filename"
echo $command
eval $command

# Pull MAF information for all chr22 snps inclding low-frequency snps (MAF 1% threshold)
in_filename="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.chr${chr}.phased.mm3.imputed.dose.snpQC.renamed.vcf.gz"
out_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/mafs/1K1K.chr${chr}.info_allMAFs.txt"

command="zcat $in_filename | bcftools query -f '%ID\t%MAF\n' > $out_filename"
echo $command
eval $command