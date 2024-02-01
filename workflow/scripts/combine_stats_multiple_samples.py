#!/usr/bin/env python

import sys
from os import path
import pandas as pd
import numpy as np
import argparse

import mip_tools

def summarize_emul(output_counts, input_counts, output_stats, input_stats):
    input_count_dfs = [pd.read_table(f, header=None, index_col=0, names=[path.basename(f).split('_emul.')[0]]) for f in input_counts]
    merged_counts = pd.concat(input_count_dfs, axis=1).fillna(0)
    merged_counts.to_csv(output_counts, sep='\t')

    input_stats_dfs = [pd.read_table(f, index_col=0) for f in input_stats]
    merged_stats = pd.concat(input_stats_dfs, axis=1).fillna(0)
    merged_stats.to_csv(output_stats, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the count matrix from a Hypr-seq experiment and return useful plots.')
    
    parser.add_argument("-d", "--output_counts", dest="output_counts", type=str, help="Where to store the output count table")
    parser.add_argument("-c", "--input_counts", dest="input_counts", type=str, nargs='+', help="Input count files")
    parser.add_argument("-t", "--output_stats", dest="output_stats", type=str, help="Where to store the output stats table")
    parser.add_argument("-s", "--input_stats", dest="input_stats", type=str, nargs='+', help="Input stats files")
    
    args = parser.parse_args()

    summarize_emul(args.output_counts, args.input_counts, args.output_stats, args.input_stats)