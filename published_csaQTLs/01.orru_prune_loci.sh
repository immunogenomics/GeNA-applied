#!/bin/bash

# Prune published loci to set of independent associations

res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/"
published_loci="${res_folder}Orru_loci.txt"

for celltype in "T" "B" "Myeloid" "allcells"
do
    infile_sumstats="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/GeNA_sumstats.txt"

    # prune Orrú loci to set of independent SNPs, gather p-values from our csaQTL assoc. tests per locus 
    outfile_pruned_loci="${res_folder}${celltype}_Orru_replication.csv"
    outfile_snplist="${res_folder}${celltype}_snplist.csv"
    command="python -u orru_prune_loci.py --published_loci $published_loci --outfile_pruned_loci $outfile_pruned_loci \
                     --outfile_snplist $outfile_snplist --infile_sumstats $infile_sumstats"
    echo $command
    eval $command

done

