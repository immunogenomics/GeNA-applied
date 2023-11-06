#!/bin/bash

# Import published eQTL data for BLUEPRINT database (sumstats retrieved via eQTL Catalog)

sumstats_folder="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/sumstats/BLUEPRINT/"
cd $sumstats_folder

for i_ref in 21 26 31
do
    command="wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD0000${i_ref}/QTD0000${i_ref}.all.tsv.gz"
    echo $command
    eval $command
done

for i_ref in 22 23 24 25 27 28 29 30 32 33 34 35
do
    command="wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD0000${i_ref}/QTD0000${i_ref}.cc.tsv.gz"
    echo $command
    eval $command
done

command="cp QTD000021.all.tsv.gz QTD000036.cc.tsv.gz"
echo $command
eval $command

command="cp QTD000026.all.tsv.gz QTD000037.cc.tsv.gz"
echo $command
eval $command

command="cp QTD000031.all.tsv.gz QTD000038.cc.tsv.gz"
echo $command
eval $command