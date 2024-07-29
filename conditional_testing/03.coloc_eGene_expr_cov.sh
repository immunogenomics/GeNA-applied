#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting            

# This script applies GeNA to re-test two csaQTLs with additional covariates                   
# relative to the discovery GWAS that reflect mean expression across cells per donor of eGenes
# implicated in eQTLs that colocalize with the csaQTL

# Genotypes into pfile format
in_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/Myeloid_15:80263217:C:T_cis.vcf.gz"
gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/Myeloid_15:80263217:C:T_cis"   
command="plink2 --vcf $in_filename  'dosage=DS' --make-pgen --out ${gtypes}"
echo $command
eval $command

# Apply GeNA
pwd=$(eval "pwd")
path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/GeNA_revisions_other/Myeloid_BCL2A1_covar.h5ad"
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_Myeloid_eGene_covs/"
gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/Myeloid_15:80263217:C:T_cis"
covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6,BCL2A1'

command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs"
echo $command
eval $command

sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/GeNA_revisions_other/NK_KLRC1_covar.h5ad"
res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_NK_eGene_covs/"
gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/NK_12:10583611:C:T_cis" 
covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6,KLRC1'

command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs"
echo $command
eval $command