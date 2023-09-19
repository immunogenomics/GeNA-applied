#!/bin/bash

# Calls a script to define cluster-based cell type proportion trait values per individual in the OneK1K dataset

outfolder="/data/srlab/lrumker/MCSC_Project/cna-qtl/cluster_gwas/"
command="python -u define_cluster_traits.py --outfolder $outfolder"
echo $command
eval $command