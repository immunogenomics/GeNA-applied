#!/bin/bash
# only needs to be run once (not per trait)
debug=1 #un-commenting this line will cause script to execute commands rather than submitting

chr=22
scorename="RA"
geno_folder="/data/srlab/lrumker/MCSC_Project/cna-prs/results/geno/"
command="python -u save_sample_labels.py ${geno_folder}${scorename}.chr${chr}.vcf.gz > ${geno_folder}sample_list.txt"
echo $command
	
if [ -z "$debug" ]
then
    mkdir -p out
    name="list_samples"
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
