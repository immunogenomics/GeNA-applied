#!/usr/bin/env python

import sys
import gzip
import pandas as pd
vcf = sys.argv[1]
with gzip.open(vcf, "rt", "utf_8") as f:
	for line in f:
		line = line.rstrip()
		if line[:len('#CHROM')] == '#CHROM':
			labels = "\n".join(line.split("\t")[9:])
			print(labels)
			quit

