import numpy as np
import pandas as pd
import argparse, cna

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sc_object",type=str)
parser.add_argument("--eGene",type=str)
parser.add_argument("--outfile",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

np.random.seed(0)

d = cna.read(args.sc_object)
d.obs[args.eGene] = d.X[:,np.where(d.var.index==args.eGene)[0][0]]
d.obs_to_sample([args.eGene], np.mean)
d.write(args.outfile)
