#!/bin/bash

module load samtools
module load bcftools

chr_start=1
chr_end=22

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/"
geno_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/vcfs/"

for chr in $(eval echo "{$chr_start..$chr_end}")
do
    in_filename="${geno_folder}1K1K.chr${chr}.selsnps.vcf.gz"
    out_filename="${geno_folder}1K1K.chr${chr}.selsnps.sorted.vcf.gz"
    command="bcftools view -S ${res_folder}sample_list.txt  $in_filename -Oz -o $out_filename"
    echo $command
    eval $command
done

for chr in $(eval echo "{$chr_start..$chr_end}")
do
    in_filename="${geno_folder}1K1K.chr${chr}.selsnps.sorted.vcf.gz"
    command="bcftools index $in_filename"
    echo $command
    eval $command
done

out_filename="${geno_folder}1K1K.all_chr.selsnps.vcf.gz"
command="bcftools concat ${geno_folder}1K1K.chr*.selsnps.sorted.vcf.gz -Oz -o $out_filename"
echo $command
eval $command

