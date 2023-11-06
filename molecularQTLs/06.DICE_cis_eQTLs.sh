#!/bin/bash

module load tabix
module load bcftools

# Import published eQTL data from DICE database, subset results to cis windows around csaQTL lead snps

src="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/DICE/"

for celltype in "MONOCYTES" "NK" "B_CELL_NAIVE" "CD4_NAIVE" "CD4_STIM" "CD8_NAIVE" "CD8_STIM" "M2" "TFH" "TH17" "TH1" "TH2" "THSTAR" "TREG_MEM" "TREG_NAIVE"
do
    command="bgzip ${src}${celltype}.vcf && tabix -p vcf ${src}${celltype}.vcf.gz"
    echo $command
    eval $command

    command="bcftools annotate --rename-chrs update_chr_strs.txt ${src}${celltype}.vcf.gz | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${src}${celltype}.renamed.vcf.gz"
    echo $command
    eval $command

    command="tabix -p vcf ${src}${celltype}.renamed.vcf.gz"
    echo $command
    eval $command

    csaQTL_celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
    lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
    n_snps=${#lead_snps[@]}

    for i in $(eval echo "{0..$(($n_snps-1))}")
    do
	csaQTL_celltype=${csaQTL_celltypes[i]}
	lead_snp=${lead_snps[i]}
	
	# get chromosome and lead snp
	IFS=":"
	read -a snp_vals <<< "$lead_snp"
	chr=${snp_vals[0]}
	lead_pos=${snp_vals[1]}
	IFS=' '
	
        # Define 2MB cis window around lead snp
	MB=1000000 # 1 megabase
	start_pos=$(($lead_pos-$MB))
	end_pos=$(($lead_pos+$MB))
	
	infile="${src}${celltype}.renamed.vcf.gz"
	outfile="${src}cis_res/${celltype}.renamed.${csaQTL_celltype}_${lead_snp}_cis.vcf.gz"
	
        # Obtain regional alleles
	command="bcftools view --regions ${chr}:${start_pos}-${end_pos} ${infile} -Oz -o ${outfile}"
	echo $command
	eval $command
    
    done
done