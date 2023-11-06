#!/bin/bash

# Define values for each OneK1K cohort member with respect to the previously-studied traits
# We include only traits that can be quantified using cell assignments to clusters in this dataset

command="python -u define_and_filter_traits.py" 
echo $command
eval $command