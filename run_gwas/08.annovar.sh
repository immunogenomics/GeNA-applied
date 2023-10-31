#!/bin/bash

# This script runs ANNOVAR to annotate all variants in the cis region around each csaQTL.

cd /data/srlab/lrumker/annovar/

celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    # get chromosome and position of lead snp                                                                                        
    IFS=':'
    read -a snp_vals <<<"$lead_snp"
    chr=${snp_vals[0]}
    lead_pos=${snp_vals[1]}
    IFS=' '

    infile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis.vcf.gz"
    outfile_stem="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/annovar/${celltype}_${lead_snp}_cis"

    command="perl table_annovar.pl $infile humandb/ -buildver hg19 -out $outfile_stem -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a  -operation gx,r,f,f,f -nastring . -vcfinput -polish -xreffile example/gene_fullxref.txt"
    echo $command
    eval $command
done
