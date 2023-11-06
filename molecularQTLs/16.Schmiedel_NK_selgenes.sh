#!/bin/bash

# Refine sumstats file from Schmiedel to reflect associations to ETS1 specifically (nominated by survey for colocalizing associations)         
# Consider only cis window around csaQTL on chromosome 11, convert positions from GRCh38 to hg19 with liftOver

cwd=$(pwd)
sumstats_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/Schmiedel/"

# Get cis region of eQLT sumstats for selected gene
csaQTL_celltype="NK"
lead_snp="11:128070535:A:G"
chr=11
GRCh38_pos="128200640"
gene_id="ENSG00000134954"

sumstats_file="Schmiedel_2018_ge_NK-cell_naive.all.tsv.gz"
sumstats_path="${sumstats_folder}${sumstats_file}"

# define 2MB cis window around lead snp in GRCh38                                                                                               
MB=1000000 # 1 megabase
lead_pos=$GRCh38_pos                                                                                                                                                   
start_pos=$(($lead_pos-$MB))
end_pos=$(($lead_pos+$MB))
    
outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.csv"
awk_str='"ENSG00000134954"'
command="zcat $sumstats_path | awk 'NR==1||\$2==${chr}' | awk 'NR==1||\$3>${start_pos}' | awk 'NR==1||\$3<${end_pos}' | awk 'NR==1||\$1==${awk_str}' > $outfile"
echo $command
eval $command

# liftOver
cd /data/srlab/lrumker/lift_over/
chainfile="hg38ToHg19.over.chain.gz"

infile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.csv"
bedfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.bed"
outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.hg19.csv"
unmapped="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.unmapped"

sub_command='{print "chr"$2 "\t" $3-1 "\t" $3 "\t" NR-1}'
command="awk 'NR>1' $infile | awk '${sub_command}' > $bedfile"
echo $command
eval $command
    
command="./liftOver $bedfile $chainfile $outfile $unmapped"
echo $command
eval $command


# reformat sumstats with hg19 positions
cd $cwd
infile_hg19_pos="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.hg19.csv"
infile_sumstats="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.csv"
outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_NK_cis_${gene_id}.hg19.renamed.csv"
    
command="python -u reformat_cis_sumstats.py --infile_hg19_pos $infile_hg19_pos --infile_sumstats $infile_sumstats --outfile $outfile"
echo $command
eval $command
