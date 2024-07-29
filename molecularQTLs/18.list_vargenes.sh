#!/bin/bash

# Generates files listing the variable genes used for each single-cell data object

out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/trans_eQTL_candidates/"

for celltype in "NK" "Myeloid"
do
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
    outfile="${out_path}${celltype}_vargenes.txt"

    command="python -u list_vargenes.py --sc_object $sc_object --outfile $outfile"
    echo $command
    eval $command
done