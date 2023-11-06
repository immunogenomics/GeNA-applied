#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting

module load samtools
module load bcftools

chr_start=1
chr_end=22
for chr in $(eval echo "{$chr_start..$chr_end}")
do
    in_filename="/data/srlab/lrumker/datasets/onek1k/geno/1K1K.phased.chr${chr}.phased.mm3.imputed.dose.snpQC.renamed.vcf.gz"
    snplist="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/replication_snps.txt"
    out_filename="/data/srlab/lrumker/MCSC_Project/cna-qtl/trait_replication/vcfs/1K1K.chr${chr}.selsnps.vcf.gz"
    
    command="zcat $in_filename | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | bcftools view --include ID==@${snplist} -Oz -o $out_filename"
    echo $command
    
    if [ -z "$debug" ]
    then
	mkdir -p out
	name="getsnps_chr${chr}"
	bsub -J $name -q short \
	    -oo out/${name}.out \
	    -eo out/${name}.err \
	    -R 'select[hname!=cn001]' \
	    -R 'select[hname!=cn002]' \
	    -R 'select[hname!=cn003]' \
	    -R 'select[hname!=cn004]' \
	    -R 'select[hname!=cn005]' \
	    "$command"
    else
	LSB_JOBINDEX=1
	eval $command
	#exit
    fi
done

