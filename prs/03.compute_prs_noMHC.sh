#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting
# NOTE: in debug mode, all jobs in loop are executed

remove_MHC="True"
geno_folder="/data/srlab/lrumker/MCSC_Project/cna-prs/results/geno/"

for scorename in "SLE" "RA"
do
    mkdir -p out

    command="python -u compute_prs.py  --scorename $scorename --remove_MHC $remove_MHC"

    echo $command
    if [ -z "$debug" ]
    then
	name="compute_prs_${dz}"
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
