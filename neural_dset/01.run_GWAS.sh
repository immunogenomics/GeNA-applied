#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting            

# This script applies GeNA to perform one GWAS for csaQTLs in the HipSci dataset

pwd=$(eval "pwd")
path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

sc_object="/data/srlab/lrumker/datasets/HipSci/pheno/HipSci.h5ad"
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_HipSci/"
gtypes="/data/srlab/lrumker/datasets/HipSci/geno/final/HipSci.all_chr.MAF5"
covs='gPC1,gPC2,gPC3,gPC4,gPC5'

command="./GeNA.sh -s $sc_object -g $gtypes -o $res_folder -c $covs"
echo $command

if [ -z "$debug" ]
then
    mkdir -p ${pwd}/out
    name="run_gwas"
    bsub -J $name -q big \
        -oo ${pwd}/out/${name}.out \
        -eo ${pwd}/out/${name}.err \
        -R 'select[hname!=cn001]' \
        -R 'select[hname!=cn002]' \
        -R 'select[hname!=cn003]' \
        -R 'select[hname!=cn004]' \
        -R 'select[hname!=cn005]' \
        "$command"
else
    LSB_JOBINDEX=1
    eval $command
    exit
fi