# This script adds genotype information for selected loci to the samplem attribute of a single-cell 
# multianndata object.
#
# Inputs
#   sc_object: filepath to a single-cell data object in multianndata format
#   gtypes: filepath to a gzipped VCF file containing only dose information for the SNPs of interest
#   gtype_samples: list of sample ID corresponding to the rows in gtypes
#   outfile: filepath to save a single-cell object with genotype information added to samplem attribute

import argparse
import pandas as pd
import numpy as np
import cna
np.random.seed(0)

# Parse Arguments                                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("--sc_object",type=str)
parser.add_argument("--outfile",type=str)
parser.add_argument("--gtypes",type=str)
parser.add_argument("--gtype_samples",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

d = cna.read(args.sc_object)
G = pd.read_csv(args.gtypes, sep = "\t", index_col = 0, header = None)
G_samples = pd.read_csv(args.gtype_samples, header = None, sep = "\t")
G.columns = G_samples.iloc[:,0].values
G.index = [G.index[i].strip() for i in np.arange(G.shape[0])]
d.samplem = d.samplem.join(G.T)
d.write(args.outfile)
