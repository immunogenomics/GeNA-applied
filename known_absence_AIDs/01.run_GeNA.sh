#!/bin/bash

# Apply GeNA to re-test the association to the lead SNP, using the masked data object

path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/known_absence_AIDs/"

celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    # get chromosome of lead snp
    IFS=':'
    read -a snp_vals <<<"$lead_snp"
    chr=${snp_vals[0]}
    IFS=' '

    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}_sigsnps"

    sc_object="/data/srlab/lrumker/MCSC_Project/cna-prs/results/sc_objects/${celltype}_PRSqc.h5ad"
    res_subfolder="${res_folder}${lead_snp}"
    mkdir -p $res_subfolder
    res_subfolder="${res_subfolder}/"
    covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'

    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_subfolder -c $covs"
    echo $command
    eval $command
done
