#!/bin/bash

module load samtools
module load bcftools

geno_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/vcfs/"
in_filename="${geno_folder}1K1K.all_chr.selsnps.vcf.gz"
out_filename="${geno_folder}1K1K.all_chr.selsnps.DS.vcf.gz"
    
command="zcat $in_filename | bcftools query -f '%ID [\t%DS]\n' | bgzip -c > $out_filename"
echo $command
eval $command
