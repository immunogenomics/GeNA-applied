#!/bin/bash

# Import published eQTL sumstats from Gilchrist et al, refine to cis windows around lead csaQTL lead snps                                                        
# and liftOver published results from GRCh38 to hg19   

# Pull sumstats
cwd=$(pwd)
sumstats_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/Gilchrist/"
cd $sumstats_folder
command="wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000014/QTD000115/QTD000115.all.tsv.gz"
echo $command
eval $command


# Get cis regions for lead snps csaQTLs
csaQTL_celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
chrs=( 15 2 11 12 19 )
GRCh38_pos=( "79970875" "111093635" "128200640" "10431012" "16331162" )
n_snps=${#lead_snps[@]}

i_ref=115
sumstats_file="QTD000${i_ref}.all.tsv.gz"
sumstats_path="${sumstats_folder}${sumstats_file}"

for j in $(eval echo "{0..$(($n_snps-1))}")
do
    csaQTL_celltype=${csaQTL_celltypes[j]} # for file naming                                                                                           
    lead_snp=${lead_snps[j]} # for file naming                                                                                                        
    lead_pos=${GRCh38_pos[j]} # for indexing into sumstats file                                                                                       
    chr=${chrs[j]}
    
    # define 2MB cis window around lead snp in GRCh38                                                                                               
    MB=1000000 # 1 megabase
                                                                                                                                                      
    start_pos=$(($lead_pos-$MB))
    end_pos=$(($lead_pos+$MB))
    
    outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.csv"
    command="zcat $sumstats_path | awk 'NR==1||\$2==${chr}' | awk 'NR==1||\$3>${start_pos}' | awk 'NR==1||\$3<${end_pos}' | awk 'NR==1||\$9<5e-4' > $outfile"
    echo $command
    eval $command
done

# liftOver

cd /data/srlab/lrumker/lift_over/
chainfile="hg38ToHg19.over.chain.gz"

for j in $(eval echo "{0..$(($n_snps-1))}")
do
    csaQTL_celltype=${csaQTL_celltypes[j]} # for file naming                                                                                                      
    lead_snp=${lead_snps[j]} # for file naming                                                                                                                    
    
    infile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.csv"
    bedfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.bed"
    outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.hg19.csv"
    unmapped="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.unmapped"
    
    sub_command='{print "chr"$2 "\t" $3-1 "\t" $3 "\t" NR-1}'
    command="awk 'NR>1' $infile | awk '${sub_command}' > $bedfile"
    echo $command
    eval $command
    
    command="./liftOver $bedfile $chainfile $outfile $unmapped"
    echo $command
    eval $command
done


# reformat sumstats with hg19 positions

cd $cwd
for j in $(eval echo "{0..$(($n_snps-1))}")
do
    csaQTL_celltype=${csaQTL_celltypes[j]} # for file naming                                                                                                      
    lead_snp=${lead_snps[j]} # for file naming                                                                                                                    
    
    infile_hg19_pos="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.hg19.csv"
    infile_sumstats="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.csv"
    outfile="${sumstats_folder}${csaQTL_celltype}_${lead_snp}_QTD000${i_ref}_cis.hg19.renamed.csv"
    
    command="python -u reformat_cis_sumstats.py --infile_hg19_pos $infile_hg19_pos --infile_sumstats $infile_sumstats --outfile $outfile"
    echo $command
    eval $command
done
