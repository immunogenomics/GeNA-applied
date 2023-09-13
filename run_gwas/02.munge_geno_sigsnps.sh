#!/bin/bash

# This script identifies the SNPs passing a 5e-8 threshold for genome-wide significant csaQTLs
# ("sigsnps") and generates a set of small genotype data objects storing only information for those SNPs

# List alleles that passed a genome-wide threshold for association 
p_thresh=5e-08 # Retain only SNPs with P < threshold 

for celltype in "NK" "Myeloid"
do
    gwas_res="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/GeNA_sumstats.txt"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}.tsv"
    
    command="awk 'NR>1{if(\$9 < $p_thresh) print \$3;next}' $gwas_res > $outfile"
    echo $command
    eval $command   
done

# Generate genotype files with only alleles that passed a genome-wide threshold for association
module load samtools
module load bcftools

for celltype in "Myeloid" "NK"
do
    snplist="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}.tsv"
    infile="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.all_chr.phased.mm3.imputed.dose.snpQC_MAF5.renamed.vcf.gz"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}_sigsnps.vcf.gz"
    dosefile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}_sigsnps.DS.vcf.gz"
    pfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}_sigsnps"

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