#!/bin/bash

# This script identifies the SNPs passing a 1e-6 threshold for suggestive associations
# ("suggsnps") and generates a set of small genotype data objects storing only information for those SNPs

p_thresh=1e-06 # Retain only SNPs with P < threshold 

# List the SNPs to include
for celltype in "NK" "Myeloid" "T" "B" "allcells"
do
    gwas_res="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/GeNA_sumstats.txt"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sugg_snps/${celltype}.tsv"
    
    command="awk 'NR>1{if(\$9 < $p_thresh) print \$3;next}' $gwas_res > $outfile"
    echo $command
    eval $command   
done

# Generate genotype files with only SNPs that passed the threshold
module load samtools
module load bcftools

for celltype in "T" "NK" "B" "Myeloid" # no sugg snps for allcells csaQTL gwas
do
    snplist="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sugg_snps/${celltype}.tsv"
    infile="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.all_chr.phased.mm3.imputed.dose.snpQC_MAF5.renamed.vcf.gz"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sugg_snps/${celltype}_suggsnps.vcf.gz"
    dosefile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sugg_snps/${celltype}_suggsnps.DS.vcf.gz"
    pfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sugg_snps/${celltype}_suggsnps"

    # VCF of sigsnps
    command="zcat $infile | bcftools view --include ID==@${snplist} -Oz -o $outfile"
    echo $command
    eval $command

    # Genotype dose file of sigsnps
    command="zcat $outfile | bcftools query -f '%ID [\t%DS]\n' | bgzip -c > $dosefile"
    echo $command
    eval $command

    # Sigsnps in pfile format
    command="plink2 --vcf $outfile --make-pgen --out $pfile"
    echo $command
    eval $command
done