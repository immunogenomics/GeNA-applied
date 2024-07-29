#!/bin/bash

# Re-organize the output files to compile the information we need to assess sensitivity of GeNA results to k value

p_thresh=1e-6

celltypes=( "Myeloid" "NK" "NK" "NK" "NK" )
lead_snps=( "15:80263217:C:T" "11:128070535:A:G" "2:111851212:C:T" "12:10583611:C:T" "19:16441973:G:A" )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/k_sensitivity/gwas_${celltype}/${lead_snp}/"
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

    command="paste ${out_path}P_k_noHeader.txt ${out_path}snp_info.txt ${out_path}t_per_nampc.txt ${out_path}P_k_noHeader.txt_all_k | awk '\$1<$p_thresh' > $outfile"
    echo $command
    eval $command
done

