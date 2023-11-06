#!/bin/bash

module load plink

traits_start=1
traits_end=40

for n_trait in $(eval echo "{$traits_start..$traits_end}")
do
    trait_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/orru_trait_est_clusters.csv"
    covs_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/covs.csv"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/orru/vcfs/1K1K.all_chr.selsnps"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/plink_clusters/T${n_trait}_chr$chr"
    pheno_name="T${n_trait}"
    
    command="plink --bfile ${gtypes} --pheno ${trait_file} --pheno-name ${pheno_name} --linear --ci 0.95 \
                 --covar ${covs_file} --covar-name sex_M, age, gPC1, gPC2, gPC3, gPC4, gPC5 --hide-covar \
                 --prune --allow-no-sex --out ${outfile}"
    echo $command
    eval $command
done