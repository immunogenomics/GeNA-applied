#!/bin/bash

# Refine sumstats file from eQTLgen to reflect associations to BCL2A1 specifically (nominated by survey for colocalizing associations)
# Consider only cis window around csaQTL on chromosome 15

sumstats_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/eQTLgen/"

# Get sumstats for selected cis region, gene and celltype
csaQTL_celltype="Myeloid"
lead_snp="15:80263217:C:T"

sumstats_file="2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
sumstats_path="${sumstats_folder}${sumstats_file}"

# get chromosome and position for lead snp
IFS=":"
read -a snp_vals <<< "$lead_snp"
chr=${snp_vals[0]}
lead_pos=${snp_vals[1]}
IFS=" "

# define 2MB cis window around lead snp (in hg19 build)                                                                                               
MB=1000000 # 1 megabase
                                                                                                                                                      
start_pos=$(($lead_pos-$MB))
end_pos=$(($lead_pos+$MB))
    
awk_str='"BCL2A1"'
outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_cis_BCL2A1.csv"
command="zcat $sumstats_path | awk 'NR==1||\$3==${chr}' | awk 'NR==1||\$4>${start_pos}' | awk 'NR==1||\$4<${end_pos}' | awk 'NR==1||\$9==${awk_str}' > $outfile"
echo $command
eval $command