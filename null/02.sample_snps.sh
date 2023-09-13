#!/bin/bash

# Sample 400 SNPs that just pass a MAF 5% threshold, at random
infile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/mafs/1K1K.chr22.info.txt"
outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/sel_snps_MAF5.txt"
sample_sz=400
command="python -u sample_targetMAF.py --outfile $outfile --sample_sz $sample_sz --infile $infile --target_maf 0.05"
echo $command
eval $command

# Sample at random 400 SNPs across all MAF above the MAF 5% threshold, split equally among MAF deciles               
min_MAF=0.05
snp_info="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/mafs/1K1K.chr22.info.txt"
outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/sel_snps_contMAF.txt"
prev_snps="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/sel_snps_MAF5.txt"
tot_snps=400
command="python -u sample_contMAF.py --outfile $outfile --tot_snps $tot_snps --prev_snps $prev_snps --snp_info $snp_info --min_MAF $min_MAF"
echo $command
eval $command

# Sample 100 SNPs that just pass a MAF 1% threshold, at random 
infile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/mafs/1K1K.chr22.info_allMAFs.txt"
outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/null/results/sel_snps_MAF1.txt"
sample_sz=100
command="python -u sample_targetMAF.py --outfile $outfile --sample_sz $sample_sz  --infile $infile --target_maf 0.01"
echo $command
eval $command