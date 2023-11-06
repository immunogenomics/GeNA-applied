#!/bin/bash

# Generates input objects for single-cell PME eQTL modeling

pbulk_res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/pseudobulk/"
out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/sceQTL/inputs/"
p_thresh=5e-4

for eQTL_celltype in "allcells" "B" "T" "NK" "Myeloid"
do
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${eQTL_celltype}.h5ad"
    expr_object="/data/srlab/lrumker/datasets/onek1k/pheno/${eQTL_celltype}_expr.h5ad"

    command="python -u gather_sceQTL_inputs.py --eQTL_celltype $eQTL_celltype --expr_object $expr_object \
                --sc_object $sc_object --pbulk_res_folder $pbulk_res_folder --p_thresh $p_thresh --out_path $out_path"
    echo $command
    eval $command
done