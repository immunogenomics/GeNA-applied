#!/bin/bash

# This script identifies the SNPs in each csaQTL GWAS that pass a p<5e-8 threshold

for celltype in "allcells" "Myeloid" "T" "B" "NK"
do
    all_sumstats="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/wCCGs/gwas_${celltype}/GeNA_sumstats.txt"
    sig_snps="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/wCCGs/gwas_${celltype}/sig_snps.txt"
    command="awk 'NR==1 || \$9<5e-8' $all_sumstats > $sig_snps"
    echo $command
    eval $command
done