#!/bin/bash

# Import published eQTL sumstats from eQTLgen database, refine to cis windows around lead csaQTL lead snps            

# Pull sumstats
sumstats_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/eQTLgen/"
cd $sumstats_folder
command="wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
echo $command
eval $command

# Get cis regions for lead snps csaQTLs
csaQTL_celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
n_snps=${#lead_snps[@]}

sumstats_file="2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
sumstats_path="${sumstats_folder}${sumstats_file}"

for j in $(eval echo "{0..$(($n_snps-1))}")
do
    csaQTL_celltype=${csaQTL_celltypes[j]} # for file naming                                                                                           
    lead_snp=${lead_snps[j]} # for file naming                                                                                                        

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
    
    outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_cis.csv"
    command="zcat $sumstats_path | awk 'NR==1||\$3==${chr}' | awk 'NR==1||\$4>${start_pos}' | awk 'NR==1||\$4<${end_pos}' | awk 'NR==1||\$1<5e-4' > $outfile"
    echo $command
    eval $command
done

command="cat ${sumstats_folder}NK_19:16441973:G:A_cis.csv | grep 'KLF2' > ${sumstats_folder}NK_19:16441973:G:A_KLF2_cis.csv"
echo $command
eval $command

command="cat ${sumstats_folder}NK_11:128070535:A:G_cis.csv | grep 'ETS1' > ${sumstats_folder}NK_11:128070535:A:G_ETS1_cis.csv"
echo $command
eval $command

command="cat ${sumstats_folder}Myeloid_15:80263217:C:T_cis.csv | grep 'BCL2A1' > ${sumstats_folder}Myeloid_15:80263217:C:T_BCL2A1_cis.csv"
echo $command
eval $command