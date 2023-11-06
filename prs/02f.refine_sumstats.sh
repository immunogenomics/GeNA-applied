#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting
# Note: debug mode executes all jobs in a loop

geno_folder="/data/srlab/lrumker/MCSC_Project/cna-prs/results/geno/"

for scorename in "SLE" "RA"
do
    mkdir -p out
    geno_file="${geno_folder}${scorename}.all_chr.DS.vcf.gz"

    command="python -u refine_sumstats.py  --geno_file $geno_file --scorename $scorename"

    echo $command
    if [ -z "$debug" ]
    then
	name="refine_sumstats_${dz}"
	bsub -J $name -q big \
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
