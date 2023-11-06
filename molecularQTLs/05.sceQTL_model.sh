#!/bin/bash

# Use single-cell PME eQTL model to test for eQTLs 
# For a given triple of (csaQTL lead SNP, eGene with association detected in pseudobulk eQTL model, major cell type within which we model expression)  
# We test for an eQTL within the expression of celltype_expr type associating with a SNP that has a csaQTL within the celltype_geno type

#debug=1 #un-commenting this line will cause script to execute commands rather than submitting                                                  

src_filepath="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/sceQTL/inputs/"
res_filepath="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/sceQTL/"
csaQTL_celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
chrs=( 15 2 11 12 19 )
n_snps=${#lead_snps[@]}

p_thresh=1e-6 # if lead SNP does not pass this threshold, wider region is not tested for other eQTLs

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    celltype_geno=${csaQTL_celltypes[i]}
    lead_snp=${lead_snps[i]}
    chr=${chrs[i]}
    geno="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype_geno}_${lead_snp}_cis.DS.vcf.gz"
    geno_ids="/data/srlab/lrumker/datasets/onek1k/geno/sample_labels/sample_list_chr${chr}.txt"

    for celltype_expr in "Myeloid" "NK" "T" "B" "allcells"
    do
    
	command="Rscript sceQTL_model.R --lead_snp $lead_snp --celltype_expr $celltype_expr --celltype_geno $celltype_geno \
                 --src_filepath $src_filepath --res_filepath $res_filepath --geno $geno --geno_ids $geno_ids --p_thresh $p_thresh"
	echo $command
	
	if [ -z "$debug" ]
	then
            mkdir -p out
            name="sceQTL_${celltype_expr}_for_${lead_snp}_${celltype_geno}"
            bsub -J $name -q big \
		-oo out/${name}.out \
		-eo out/${name}.err \
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
            #exit
	fi
    done
done