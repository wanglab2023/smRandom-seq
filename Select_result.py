
# Import the necessary packages.
import sys
from pyfastx import Fastq
import gzip
import pandas as pd
import re
import pysam

expr = sys.argv[1]
thres = int(sys.argv[2])
outPrefix = sys.argv[3]

sExpr = outPrefix + '.counts.tsv'

df = pd.read_table(expr, index_col = 0)
df = df.T
obs = pd.DataFrame(index = df.index)
obs['geneCount'] = (df > 0).astype(float).sum(axis=1)
obs.sort_values(by='geneCount', axis=0, ascending=False, inplace=True)
selectCells = obs[0:thres].index
pattern = '|'.join(selectCells)
indexbool = df.index.str.contains(pattern)
sdf = df[indexbool]
sdf = sdf.T
sdf.to_csv(sExpr, sep = '\t')
