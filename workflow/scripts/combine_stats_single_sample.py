#!/usr/bin/env python

import sys
from os import path
import pandas as pd
import numpy as np
import argparse

import mip_tools

def get_sample_name(output_file_name, experiment_type):
    return path.basename(output_file_name).split(experiment_type)[0]


def combine_stats(input_stats_files, output_stats_file, experiment):
    sample_name = get_sample_name(output_stats_file, experiment)
    input_stats = [pd.read_table(f, header=None, index_col=0, names=[sample_name]) for f in input_stats_files]
    df = input_stats.pop(0)
    appended = df.append(input_stats)
    appended.to_csv(output_stats_file, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the count matrix from a Hypr-seq experiment and return useful plots.')
    
    parser.add_argument("-o", "--output_stats", dest="output_stats", type=str, help="Where to store the stats table")
    parser.add_argument("-i", "--input_stats", dest="input_stats", type=str, nargs='+', help="Input stats files")
    parser.add_argument("-e", "--experiment", dest="experiment", type=str, help="What kind of experiment (as in, the prefix to use to grab the sample name)")
    
    args = parser.parse_args()

    combine_stats(args.input_stats, args.output_stats, args.experiment)