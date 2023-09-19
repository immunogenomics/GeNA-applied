#!/bin/bash

# This script converts the permutation-filled vcf files per lead SNP to pgen format for use as input to GeNA

module load plink

lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
celltypes=( "Myeloid" "NK" "NK" "NK" "NK" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    snp=${lead_snps[i]}
    celltype=${celltypes[i]}

    in_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/leadsnps_perm/results/${snp}.perm.vcf.gz"
    out_stem="/data/srlab/lrumker/MCSC_Project/cna-qtl/leadsnps_perm/results/${snp}.perm"
    
    command="plink2 --vcf $in_filename 'dosage=DS' --make-pgen --out $out_stem"
    echo $command
    eval $command
done