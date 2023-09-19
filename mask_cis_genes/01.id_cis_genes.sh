#!/bin/bash

# This script identifies the genes within a 2MB window centered on the lead SNP for each locus

gtf_path="/data/srlab/lrumker/gencode/gencode.v38lift37.annotation.gtf" 

for celltype in "NK" "Myeloid"
do
    outfolder="/data/srlab/lrumker/MCSC_Project/cna-qtl/mask_cis/results/"
    loci="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/gwas_loci.tsv"
	
    command="python -u id_cis_genes.py \
                 --outfolder $outfolder \
                 --loci $loci \
                 --gtf_path $gtf_path \
                 --celltype $celltype"
    echo $command
    eval $command
done
