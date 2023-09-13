#!/bin/bash

# This script produces genotype files for each listed (lead) SNP containing information
# about all SNPs with MAF > 1% in a 2MB window centered on the selected SNP.

module load samtools
module load bcftools

celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A")
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype=${celltypes[i]}
    lead_snp=${lead_snps[i]}

    # Get chromosome and position of lead snp
    IFS=':'
    read -a snp_vals <<< "$lead_snp"
    chr=${snp_vals[0]}
    lead_pos=${snp_vals[1]}
    IFS=' '

    # Define 2MB cis window around lead snp
    MB=1000000 # 1 megabase
    start_pos=$(($lead_pos-$MB))
    end_pos=$(($lead_pos+$MB))

    infile="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.chr${chr}.phased.mm3.imputed.dose.snpQC.renamed.vcf.gz"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis.vcf.gz"
    dose_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis.DS.vcf.gz"

    # Obtain regional alleles
    command="bcftools view --regions ${chr}:${start_pos}-${end_pos} ${infile} -Oz -o ${outfile}"
    echo $command
    eval $command
    
    # Get dose file
    command="zcat $outfile | bcftools query -f '%ID [\t%DS]\n' | bgzip -c > $dose_file"
    echo $command
    eval $command

    # Generate bfile format for PLINK
    in_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis.vcf.gz"
    out_stem="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis"

    command="plink --vcf $in_filename --make-bed --out $out_stem"
    echo $command
    eval $command

done