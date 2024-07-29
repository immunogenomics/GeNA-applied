#!/bin/bash

# Create file with desired values of k (1-50) listed to use as input to GeNA

outfile="large_k.txt"
n_nampcs=50
command="seq 1 ${n_nampcs} > $outfile"
echo $command
eval $command