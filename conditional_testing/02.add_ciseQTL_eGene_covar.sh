#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting            

# This script adds mean (normalized) eGene expression per donor to d.samplem as an available covariate

eGene="KLRC1"
celltype="NK"
sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
outfile="/data/srlab/lrumker/datasets/onek1k/pheno/GeNA_revisions_other/${celltype}_${eGene}_covar.h5ad"

command="python -u add_eGene_covar.py --outfile $outfile --sc_object $sc_object --eGene $eGene"
echo $command
eval $command

eGene="BCL2A1"
celltype="Myeloid"
sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}.h5ad"
outfile="/data/srlab/lrumker/datasets/onek1k/pheno/GeNA_revisions_other/${celltype}_${eGene}_covar.h5ad"

command="python -u add_eGene_covar.py --outfile $outfile --sc_object $sc_object --eGene $eGene"
echo $command
eval $command