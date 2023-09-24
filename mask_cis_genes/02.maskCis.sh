#!/bin/bash

# Generates a single-cell object lacking any variable genes from the discovery data
# object that lie within a 2MB cis window around the lead SNP
# Defines the phenotype most correlated with the lead SNP in the masked object

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/mask_cis/results/"

sel_ks=( 3 17 17 17 17 )
celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}
    sel_k=${sel_ks[i]}
    
    # get chromosome and position of lead snp                                                                                                                      
    IFS=':'
    read -a snp_vals <<<"$lead_snp"
    chr=${snp_vals[0]}
    lead_pos=${snp_vals[1]}
    IFS=' '

    cis_genes="${res_folder}cis_genes_${celltype}_${lead_snp}.tsv"
    out_sc_object="${res_folder}masked_${celltype}_${lead_snp}.h5ad"
    cna_res="${res_folder}${celltype}_${lead_snp}_maskCis.p"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}_sigsnps.DS.vcf.gz"
    gtype_samples="/data/srlab/lrumker/datasets/onek1k/geno/sample_labels/sample_list_rowsep.txt"
	
    command="python -u test_loci_maskCis.py \
                 --out_sc_object $out_sc_object \
                 --cna_res $cna_res \
                 --cis_genes $cis_genes \
                 --lead_snp $lead_snp \
                 --gtypes $gtypes \
                 --gtype_samples $gtype_samples \
                 --sel_k $sel_k \
                 --celltype $celltype"
    echo $command
    eval $command

done
