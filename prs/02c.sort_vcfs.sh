#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting

module load samtools
module load bcftools

geno_folder="/data/srlab/lrumker/MCSC_Project/cna-prs/results/geno/"

for scorename in "RA" "SLE"
do
    chr_start=1
    chr_end=22
    for chr in $(eval echo "{$chr_start..$chr_end}")
    do
        in_filename="${geno_folder}${scorename}.chr${chr}.vcf.gz"
        out_filename="${geno_folder}${scorename}.chr${chr}.sorted.vcf.gz"
        sample_list="${geno_folder}sample_list.txt"

	if test -f "$in_filename"; then
            command="zcat $in_filename | bcftools view -S ${sample_list} -Oz -o $out_filename"
	    echo $command
	    
	    if [ -z "$debug" ]
	    then
		mkdir -p out
		name="sortvcf_${scorename}_chr${chr}"
		bsub -J $name -q vshort \
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
	fi
    done
done

