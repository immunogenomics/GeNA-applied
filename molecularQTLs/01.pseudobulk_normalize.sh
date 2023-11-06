#!/bin/bash

# Generate pseudobulk expression values per gene per individual in the OneK1K cohort

debug=1 #un-commenting this line will cause script to execute commands rather than submitting                                        

for celltype in "B" "T" "allcells" "Myeloid" "NK"
do
    celltype_expr="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}_expr.mtx"
    celltype_meta="/data/srlab/lrumker/datasets/onek1k/pheno/${celltype}_cellmeta.csv"
    genenames="/data/srlab/lrumker/datasets/onek1k/pheno/T_passQC_genenames.csv" # same for all celltypes
    out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/pseudobulk/inputs/${celltype}"
    sample_meta="/data/srlab/lrumker/datasets/onek1k/pheno/sample_meta.csv"
    
    command="Rscript pseudobulk_normalize.R --celltype_expr $celltype_expr --sample_meta $sample_meta \
              --celltype_meta $celltype_meta  --genenames $genenames --out_path $out_path"
    echo $command
    
    if [ -z "$debug" ]
    then
        mkdir -p out
        name="pseudobulk_norm_${celltype}"
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