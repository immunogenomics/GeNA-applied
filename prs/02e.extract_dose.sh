#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting
# NOTE: debug model executes all jobs in loop

module load samtools
module load bcftools

geno_folder="/data/srlab/lrumker/MCSC_Project/cna-prs/results/geno/"

for scorename in "SLE" "RA"
do
    mkdir -p out

    in_filename="${geno_folder}${scorename}.all_chr.vcf.gz"
    out_filename="${geno_folder}${scorename}.all_chr.DS.vcf.gz"

    command="zcat $in_filename | bcftools query -f '%ID [\t%DS]\n' | bgzip -c > $out_filename"

    echo $command
    if [ -z "$debug" ]
    then
	name="extract_dose_${dz}"
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
