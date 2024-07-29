#!/bin/bash

# Re-organize the output files to compile the information we need to assess sensitivity of GeNA results to k value

for celltype in "allcells" "NK" "T" "B" "Myeloid"
do
    out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/k_sensitivity/gwas_${celltype}_nonnull_sims/"
    outfile="${out_path}P_k_t_sugg_snps.txt"

    command="awk 'NR>1' ${out_path}GeNA_sumstats.txt > ${out_path}snp_info.txt"
    echo $command
    eval $command

    command="awk 'NR>1' ${out_path}P_k.txt > ${out_path}P_k_noHeader.txt"
    echo $command
    eval $command 

    command="awk 'NR>1' ${out_path}P_k.txt_all_k > ${out_path}P_k_noHeader.txt_all_k"
    echo $command
    eval $command

    command="paste ${out_path}P_k_noHeader.txt ${out_path}snp_info.txt ${out_path}t_per_nampc.txt ${out_path}P_k_noHeader.txt_all_k > $outfile"
    echo $command
    eval $command
done

