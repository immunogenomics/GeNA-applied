#!/bin/bash
#debug=1 #un-commenting this line will cause script to execute commands rather than submitting            

# This script applies GeNA to test associations for SNPs with p<5e-8 from the discovery 
# csaQTL GWAS in custom sub-cohorts of the overall OneK1K cohort

pwd=$(eval "pwd")
path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

# Cohorts reflecting samples from individuals with a known absence of autoimmune disease
for celltype in "NK" "Myeloid"
do
    src_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/"
    sc_object="${src_folder}gwas_${celltype}/${celltype}_noAIDs.h5ad"
    gtypes="${src_folder}/geno_munge/sig_snps/${celltype}_sigsnps"
    res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/retest_subcohorts/${celltype}_noAIDs/"
    mkdir -p $res_folder
    covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'
    
    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs"
    echo $command
    eval $command
done

# Other custom NK objects
celltype="NK"
for suffix in "KnownMeta" "noAsthma"
do
    src_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/"
    sc_object="${src_folder}gwas_${celltype}/${celltype}_${suffix}.h5ad"
    gtypes="${src_folder}/geno_munge/sig_snps/${celltype}_sigsnps"
    res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/retest_subcohorts/${celltype}_${suffix}/"
    mkdir -p $res_folder
    covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'

    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs"
    echo $command
    eval $command
done
