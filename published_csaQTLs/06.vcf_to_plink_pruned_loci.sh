#!/bin/bash

module load plink

geno_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/vcfs/"
in_filename="${geno_folder}1K1K.all_chr.selsnps.vcf.gz"
out_stem="${geno_folder}1K1K.all_chr.selsnps"

command="plink --vcf $in_filename --make-bed --out $out_stem"
echo $command
eval $command
