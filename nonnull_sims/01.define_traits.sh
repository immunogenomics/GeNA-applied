#!/bin/bash

# Calls a script to define the real single-cell state abundance phenotypes
# within each data object

for majortype in "allcells" "T" "B" "NK" "Myeloid"
do
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/${majortype}.h5ad"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/traits/${majortype}_traits.txt"
    command="python -u define_traits.py --majortype ${majortype} --sc_object $sc_object --outfile $outfile"
    echo $command
    eval $command
done