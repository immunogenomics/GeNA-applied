#!/bin/bash
#debug=1 #un-commenting this line will cause script to execute commands rather than submitting 

# Test replication of sex-associated cell state abundance shifts from OneK1K in the Perez et al. cohorts

pwd=$(eval "pwd")

for celltype in "Myeloid" "B" "NK" "T"
do
    for cohort in "ASI" "EUR"                                                                                                                    
    do
	command="python -u replication_NAM.py --celltype $celltype  --cohort $cohort"
	echo $command

	if [ -z "$debug" ]
	then
            mkdir -p ${pwd}/out
            name="replication_${celltype}_${cohort}"
            bsub -J $name -q big \
		-R 'rusage[mem=100000]' \
		-oo ${pwd}/out/${name}.out \
		-eo ${pwd}/out/${name}.err \
		-R 'select[hname!=cn001]' \
		-R 'select[hname!=cn002]' \
		-R 'select[hname!=cn003]' \
		-R 'select[hname!=cn004]' \
		-R 'select[hname!=cn005]' \
		"$command"
	else
            LSB_JOBINDEX=1
            eval $command
            exit
	fi
    done
done