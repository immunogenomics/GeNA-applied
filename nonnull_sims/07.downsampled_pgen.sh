#!/bin/bash

# This script converts our simulated genotypes from VCF to pgen format

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/"

for suffix in "0pt8" "0pt6" "0pt4" "0pt2"
do
    for celltype in "allcells" "Myeloid" "NK" "B" "T"
    do
	in_filename="${res_folder}sim_genotypes/${celltype}_${suffix}_sim_genotypes.vcf.gz"
	out_stem="${res_folder}sim_genotypes/${celltype}_${suffix}_sim_genotypes"
	
	command="plink2 --vcf $in_filename 'dosage=DS' --make-pgen --out $out_stem"
	echo $command
	eval $command
    done
done