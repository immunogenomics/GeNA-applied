#!/bin/bash

# This script identifies the top variable genes *outside* a 2MB window centered on the lead SNP for each locus
# with a suggestive eQTL association from a pseudobulk eQTL model

celltypes=( "Myeloid" "NK" "NK" "NK" "NK" )
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/mask_trans_eGenes/trans_eGenes_${celltype}_${lead_snp}.txt"
	
    command="python -u id_trans_eGenes.py \
                 --outfile $outfile \
                 --lead_snp $lead_snp \
                 --celltype $celltype"
    echo $command
    eval $command
done

# Combine by cell type
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/mask_trans_eGenes/"
command="cat ${res_folder}trans_eGenes_NK_*.txt > ${res_folder}trans_eGenes_NK.txt"
echo $command
eval $command

command="cat ${res_folder}trans_eGenes_Myeloid_*.txt > ${res_folder}trans_eGenes_Myeloid.txt"
echo $command
eval $command