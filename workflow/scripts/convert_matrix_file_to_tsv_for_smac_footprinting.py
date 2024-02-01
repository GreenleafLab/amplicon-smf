import sys
import pandas as pd
import os
from os import path
import matplotlib.pyplot as plt

matrix_file = sys.argv[1]
tsv_out_file = sys.argv[2]

df = pd.read_table(matrix_file, index_col=0, header=0, sep='\t')

chrom = df.index.name.strip('#')
col_range = list(map(int, df.columns.tolist()))
start, end = min(col_range), max(col_range)
strand = '+'

# grab only methylatable cols
# might want to avoid the duplicates here?
cols_to_write = df.loc[:,df.mean()>0].columns.tolist()

# iterate through all reads
to_write = []

for i in range(len(df)):
    read = df.iloc[i]
    positions = ','.join(map(str,cols_to_write))
    loglikes = ','.join(map(str,read[cols_to_write].tolist()))
    read_lst = [chrom, start, end, strand, 'read_{}'.format(i), '.', positions, loglikes]
    to_write.append(read_lst)
    
pd.DataFrame(to_write).to_csv(tsv_out_file, sep='\t', header=False, index=False)
