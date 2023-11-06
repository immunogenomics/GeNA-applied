#!/bin/bash

module load plink

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/"
for celltype in "B" "T" "allcells" "Myeloid" # none in NK
do
    trait_file="orru_trait_est_${celltype}_nampcs.csv"
    gtypes="${res_folder}vcfs/1K1K.all_chr.selsnps"
    outfile="${res_folder}plink_nampcs/nampcs"
    
    command="plink --bfile ${gtypes} --pheno ${trait_file} --all-pheno --linear --ci 0.95 \
                 --prune --allow-no-sex --out ${outfile}"
    echo $command
    eval $command
done