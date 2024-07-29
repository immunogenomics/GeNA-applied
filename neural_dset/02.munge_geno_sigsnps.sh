#!/bin/bash

p_thresh=1e-06 # Retain only SNPs with P < threshold 

for celltype in "HipSci"
do
    gwas_res="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/gwas_${celltype}/GeNA_sumstats.txt"
    outfile="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}.tsv"
    
    command="awk 'NR>1{if(\$9 < $p_thresh) print \$3;next}' $gwas_res > $outfile"
    echo $command
    eval $command   
done

# Generate genotype files with only the selected alleles
module load samtools
module load bcftools

for celltype in "HipSci"
do
    snplist="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/sig_snps/${celltype}.tsv"
    infile="/data/srlab/lrumker/datasets/HipSci/geno/final/HipSci.all_chr.MAF5.vcf.gz"
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