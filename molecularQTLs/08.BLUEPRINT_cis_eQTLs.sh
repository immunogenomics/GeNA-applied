#!/bin/bash

# Subset the published summary statistics to cis windows around the csaQTL lead snps
# Note that published sumstats are in GRCh38 so we convert the position of the csaQTL lead snp

sumstats_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/BLUEPRINT/"

csaQTL_celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
chrs=( 15 2 11 12 19 )
GRCh38_pos=( "79970875" "111093635" "128200640" "10431012" "16331162" )
n_snps=${#lead_snps[@]}

for i_ref in $(eval echo "{21..38}")
do
    sumstats_file="QTD0000${i_ref}.cc.tsv.gz"
    sumstats_path="${sumstats_folder}${sumstats_file}"

    for j in $(eval echo "{0..$(($n_snps-1))}")
    do
	csaQTL_celltype=${csaQTL_celltypes[j]} # for file naming 
	lead_snp=${lead_snps[j]} # for file naming
	lead_pos=${GRCh38_pos[j]} # for indexing into sumstats file
	chr=${chrs[j]}
	
        # define 2MB cis window around lead snp in GRCh38                                                                          
	MB=1000000 # 1 megabase                                                                                                                            
	start_pos=$(($lead_pos-$MB))
	end_pos=$(($lead_pos+$MB))
	
	outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD0000${i_ref}_cis.csv"
	
	command="zcat $sumstats_path | awk 'NR==1||\$2==${chr}' | awk 'NR==1||\$3>${start_pos}' | awk 'NR==1||\$3<${end_pos}' | awk 'NR==1||\$9<5e-4' > $outfile"
	echo $command
	eval $command
    done
done