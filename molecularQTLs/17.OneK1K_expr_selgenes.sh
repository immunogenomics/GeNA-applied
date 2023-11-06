#!/bin/bash

# Gather expression for selected genes (nominated by the csaQTL phenotype annotation process)
# to use as input for differential expression analyses

out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/sceQTL/inputs/custom_"
csaQTL_celltype="Myeloid"
eQTL_celltype="allcells"
lead_snp="15:80263217:C:T"
sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${eQTL_celltype}.h5ad"
expr_object="/data/srlab/lrumker/datasets/onek1k/pheno/${eQTL_celltype}_expr.h5ad"

command="python -u custom_sceQTL_inputs.py --eQTL_celltype $eQTL_celltype --expr_object $expr_object \
                   --sc_object $sc_object --csaQTL_celltype $csaQTL_celltype --out_path $out_path --lead_snp $lead_snp"
echo $command
eval $command