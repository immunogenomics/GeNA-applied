#!/bin/bash

# This script converts cis region genotypes from VCF to pgen format for all genome-wide sig OneK1K csaQTLs

celltypes=( "Myeloid" "NK" "NK" "NK" "NK" )
lead_snps=( "15:80263217:C:T" "11:128070535:A:G" "2:111851212:C:T" "12:10583611:C:T" "19:16441973:G:A" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    in_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis.vcf.gz"
    out_stem="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis"

    command="plink2 --vcf $in_filename 'dosage=DS' --make-pgen --out $out_stem"
    echo $command
    eval $command
done