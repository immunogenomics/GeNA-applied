#!/bin/bash
#debug=1

# Applies PLINK to GWAS for cluster abundance traits
# Considers age, sex, gPC1-6 as covariates
# Excludes individuals without phenotypes available for the trait (--prune)

module load plink

srcfolder="/data/srlab/lrumker/MCSC_Project/cna-qtl/cluster_gwas/"

traits_file="${srcfolder}cluster_traits.tsv"
covs="${srcfolder}covs.tsv"
gtypes="/data/srlab/lrumker/datasets/onek1k/geno/plink/onek1k.QCed_MAF5.all_chr"

n_traits=28
for n_trait in $(eval echo "{1..$n_traits}")
do
    outfile="${srcfolder}gwas_res/T${n_trait}"

    command="plink --bfile ${gtypes} --pheno ${traits_file} --mpheno ${n_trait} --linear hide-covar \
                       --allow-no-sex --covar ${covs} --covar-name age,sex_M,gPC1-gPC6 \
                       --out ${outfile}"

    echo $command    
    if [ -z "$debug" ]
    then
	mkdir -p out
	name="T${n_trait}"
	bsub -J $name -q big \
	    -oo ${srcfolder}out/${name}.out \
	    -eo ${srcfolder}out/${name}.err \
	    -R 'select[hname!=cn001]' \
	    -R 'select[hname!=cn002]' \
	    -R 'select[hname!=cn003]' \
	    -R 'select[hname!=cn004]' \
	    -R 'select[hname!=cn005]' \
	    -R 'select[hname!=cn007]' \
	    "$command"
    else
	LSB_JOBINDEX=1
	eval $command
	exit
    fi
done