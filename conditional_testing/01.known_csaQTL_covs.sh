#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting            

# This script applies GeNA to test the csaQTL with lead SNP rs3003, with additional covariates
# relative to the discovery GWAS that approximate cell state abundance traits found to associate
# with this locus in flow cytometry studies of peripheral blood

# Genotypes into pfile format
in_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/NK_12:10583611:C:T_cis.vcf.gz"
gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/NK_12:10583611:C:T_cis"   
command="plink2 --vcf $in_filename  'dosage=DS' --make-pgen --out ${gtypes}"
echo $command
eval $command

# Apply GeNA
pwd=$(eval "pwd")
path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/GeNA_revisions_other/NK.h5ad"
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_NK_flow_covs/"

covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6,Lymphocyte_Count_INT,Monocyte_Count_INT,CD4_TEM_Pct_T_INT'

command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs"
echo $command
eval $command