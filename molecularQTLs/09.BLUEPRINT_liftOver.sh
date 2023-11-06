#!/bin/bash

# liftOver published summary statistics for BLUEPRINTS database in csaQTL cis windows to hg19 build

cd /data/srlab/lrumker/lift_over/
chainfile="hg38ToHg19.over.chain.gz"

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
	
	infile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.csv"
	bedfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.bed"
	outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.hg19.csv"
	unmapped="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.unmapped"
	
	sub_command='{print "chr"$2 "\t" $3-1 "\t" $3 "\t" NR-1}'
	command="awk 'NR>1' $infile | awk '${sub_command}' > $bedfile"
	echo $command
	eval $command
	
	command="./liftOver $bedfile $chainfile $outfile $unmapped"
	echo $command
	eval $command
    done
done