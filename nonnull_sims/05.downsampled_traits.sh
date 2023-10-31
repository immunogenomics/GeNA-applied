#!/bin/bash

# Calls a script to quantify per-sample values for real cell state abundance phenotypes
# within each data object. Here, the data objects are downsampled versions of the full
# OneK1K objects containing progressively fewer total samples.

for suffix in "0pt8" "0pt6" "0pt4" "0pt2"
do
    for majortype in "Myeloid" "allcells" "NK" "B" "Myeloid" "T"
    do
	sc_object="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/downsampled/${majortype}_${suffix}.h5ad"
	outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/nonnull_sims/traits/${majortype}_${suffix}_traits.txt"
	command="python -u define_traits.py --majortype ${majortype} --sc_object $sc_object --outfile $outfile"
	echo $command
	eval $command
    done
done