#!/bin/bash

# Converts permuted genotypes to pfile in preparation for GeNA
for celltype in "Myeloid" "T" "NK" "B" "allcells"
do
    for MAF_group in "MAF5" "contMAF" "MAF1"
    do
	in_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_${MAF_group}.perm_${celltype}.vcf.gz"
	out_stem="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/1K1K.selsnps_${MAF_group}.perm_${celltype}"
    
	command="plink2 --vcf $in_filename 'dosage=DS' --make-pgen --out $out_stem"
	echo $command
	eval $command
    done
done