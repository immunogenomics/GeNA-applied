#!/bin/bash

# This script produces a genotype file for each (suggestive) csaQTL (indexed by lead SNP) containing
# information about all SNPs with MAF > 5% within a 2MB window centered on the lead SNP

module load samtools
module load bcftools

celltypes=( "HipSci")
lead_snps=( "18:54746518" )
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

    infile="/data/srlab/lrumker/datasets/HipSci/geno/final/HipSci.chr${chr}.MAF5.vcf.gz"
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