#!/bin/bash
#debug=1 #un-commenting this line will cause script to execute commands rather than submitting            

# This script applies GeNA to perform one GWAS for csaQTLs in each of five single-cell objects,
# each corresponding to a different denominator of included cells (T, B, NK, myeloid or all cells)
# This analysis matches our primary csaQTL GWAS except that here we use objects with alternative 
# preprocessing: cell cycle genes are retained during construction of the cell embedding

pwd=$(eval "pwd")
path_to_GeNA="/data/srlab/lrumker/MCSC_Project/cna-qtl/GeNA"
cd $path_to_GeNA

for celltype in "allcells" "Myeloid" "T" "B" "NK"
do
    sc_object="/data/srlab/lrumker/datasets/onek1k/pheno/CCGs_retained_GeNA_revisions/${celltype}_wCCGs.h5ad"
    res_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/wCCGs/gwas_${celltype}/"
    gtypes="/data/srlab/lrumker/datasets/onek1k/geno/plink/onek1k.QCed_MAF5.all_chr"
    covs='age,sex_M,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6'

    command="./GeNA.sh -s $sc_object -b 'True' -g $gtypes -o $res_folder -c $covs"
    echo $command

    if [ -z "$debug" ]
    then
	mkdir -p ${pwd}/out
	name="gwas_${celltype}"
	bsub -J $name -q big \
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