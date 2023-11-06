#!/bin/bash

# Integrate published sumstats information with lifted-over positions

sumstats_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/BLUEPRINT/"

csaQTL_celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
n_snps=${#lead_snps[@]}

for i_ref in $(eval echo "{21..38}")                                                                                                              
do
    for j in $(eval echo "{0..$(($n_snps-1))}")
    do
        csaQTL_celltype=${csaQTL_celltypes[j]} # for file naming                                                                                     
        lead_snp=${lead_snps[j]} # for file naming                                                                                                   
	
	infile_hg19_pos="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.hg19.csv"
	infile_sumstats="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.csv" 
	outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.hg19.renamed.csv"
	
	command="python -u reformat_cis_sumstats.py --infile_hg19_pos $infile_hg19_pos --infile_sumstats $infile_sumstats --outfile $outfile"
	echo $command
	eval $command
     done
done