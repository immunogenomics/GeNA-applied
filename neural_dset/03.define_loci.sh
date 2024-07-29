#!/bin/bash

p_thresh=1e-06

# Identify suggestive loci

for celltype in "HipSci"
do
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/gwas_loci.tsv"
    gwas_res="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/GeNA_sumstats.txt"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}_sigsnps"
    n_gtype_header=40

    command="python -u /data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA-applied/run_gwas/define_loci.py \
                       --celltype $celltype \
                       --outfile $outfile \
                       --p_thresh ${p_thresh} \
                       --gtypes $gtypes \
                       --n_gtype_header $n_gtype_header \
                       --gwas_res ${gwas_res}"
    echo $command
    eval $command
    
    source_sc_obj="/data/srlab/lrumker/datasets/HipSci/pheno/${celltype}.h5ad"
    new_sc_obj="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/${celltype}_wgeno.h5ad"
    gtypes="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}_sigsnps.DS.vcf.gz"
    gtype_samples="/data/srlab/lrumker/datasets/HipSci/geno/sample_list_rowsep.txt"
    
    command="python -u /data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA-applied/run_gwas/add_gtypes_to_sc.py \
                       --sc_object $source_sc_obj \
                       --outfile $new_sc_obj \
                       --gtypes $gtypes \
                       --gtype_samples $gtype_samples"
    echo $command
    eval $command

done