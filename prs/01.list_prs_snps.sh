#!/bin/bash

# Export list of SNPs in the PRS

for scorename in "RA" "SLE"
do
    command="python -u listsnps.py --scorename $scorename"
    echo $command
    eval $command
done
