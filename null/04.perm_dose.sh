#!/bin/bash

# Generates VCF files containing permutations (across samples) of the SNPs in the input VCF

nperm=200
n_header=15

for celltype in "Myeloid" "B" "T" "NK" "allcells"
do
    nampc_file="/data/srlab/lrumker/datasets/onek1k/pheno/nampcs/${celltype}_batch_covs.csv"
    
    for MAF_group in "MAF5" "contMAF" "MAF1"
    do
	header_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/vcf_header_${MAF_group}.txt"
	infile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_${MAF_group}.vcf.gz"
	outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_${MAF_group}.perm_${celltype}.vcf"

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
done
