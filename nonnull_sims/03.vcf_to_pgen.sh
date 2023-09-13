#!/bin/bash

# This script converts our simulated genotypes from VCF to pgen format

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/"
for celltype in "Myeloid" "T" "NK" "B" "allcells"
do
    in_filename="${res_folder}sim_genotypes/${celltype}_sim_genotypes.vcf.gz"
    out_stem="${res_folder}sim_genotypes/${celltype}_sim_genotypes"
    
    command="plink2 --vcf $in_filename 'dosage=DS' --make-pgen --out $out_stem"
    echo $command
    eval $command
done