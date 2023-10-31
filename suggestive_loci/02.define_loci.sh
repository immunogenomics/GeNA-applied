#!/bin/bash

# This script calls a function to define independent loci (for suggestive associations with p<1e-6) and their respsective lead SNPs among the
# GeNA csaQTL summary statistics. Then, it calls a function to add genotype information for those SNPs to the single-cell object.
# Note that all genome-wide significant loci are still included here

p_thresh=1e-06
cd /data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA-applied/run_gwas/

for celltype in "B" "T" "NK" "Myeloid"
do
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/gwas_suggestive_loci.tsv"
    gwas_res="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/GeNA_sumstats.txt"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sugg_snps/${celltype}_suggsnps"
    n_gtype_header=40

    command="python -u define_loci.py \
                       --celltype $celltype \
                       --outfile $outfile \
                       --p_thresh ${p_thresh} \
                       --gtypes $gtypes \
                       --n_gtype_header $n_gtype_header \
                       --gwas_res ${gwas_res}"
    echo $command
    eval $command
    
    source_sc_obj="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
    new_sc_obj="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/${celltype}_wgeno.h5ad"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sugg_snps/${celltype}_suggsnps.DS.vcf.gz"
    gtype_samples="/data/srlab/lrumker/datasets/onek1k/geno/sample_labels/sample_list_rowsep.txt"
    
    command="python -u add_gtypes_to_sc.py \
                       --sc_object $source_sc_obj \
                       --outfile $new_sc_obj \
                       --gtypes $gtypes \
                       --gtype_samples $gtype_samples"
    echo $command
    eval $command

done