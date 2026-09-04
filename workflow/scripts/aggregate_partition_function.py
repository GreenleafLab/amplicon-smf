# aggregate partition function outputs

import pandas as pd
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns
from sklearn.cluster import KMeans
import pickle
from copy import deepcopy
import gc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collates all partition function model outputs into one big file')
    parser.add_argument("--partition_function_dir", type=str, help="Path to highlevel folder containing multiple partition function runs")
    parser.add_argument("--samples", type=str, nargs='+', default=None, help="Space separated list of samples to aggregate")
    parser.add_argument("--output_dir", type=str, help="Output director for file")
    parser.add_argument("--out_prefix", type=str, help="Output prefix for file")
    parser.add_argument("--partition_model_lib", type=str, help="What we are fitting on (default opjs4)", default='opjs4')

    args = parser.parse_args()

    series_dict = []

    for samp in args.samples:
        p = path.join(args.partition_function_dir, f'{samp}.fit.{args.partition_model_lib}.independent_fit.txt')
        if path.exists(p):
            s = pd.read_table(p, header=None, index_col=0)[1]
            s.name = samp
            series_dict.append(s)

    df = pd.DataFrame(series_dict)
    df.index.name = 'sample'
    
    outpath = path.join(args.output_dir, f'{args.out_prefix}.partition_function_model_params.txt.gz')
    df.to_csv(outpath, sep='\t', index=True, header=True)
