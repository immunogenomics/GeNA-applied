import numpy as np
import pandas as pd
import cna
import argparse
np.random.seed(0)

# Parse Arguments                                                                                                                                          
parser = argparse.ArgumentParser()
parser.add_argument("--sc_object",type=str)
parser.add_argument("--outfile",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

# Import data object
d = cna.read(args.sc_object)

# Export vargenes
pd.DataFrame({"gene":d.var.index}).to_csv(args.outfile, header = None, index = False)
