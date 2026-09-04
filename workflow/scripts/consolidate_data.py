# this loads in all individually aggregated data and outputs 3 big files that I should be able to use for all downstream analyses...

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

from common import load_single_molecule_matrix, load_tfbs_positions

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collates all binding model files into one big file')
    parser.add_argument("--promoter_single_molecule", type=str, help="Path to promoter single molecule table")
    parser.add_argument("--binding_single_molecule", type=str, help="Path to binding model single molecule table")
    parser.add_argument("--promoter_amp", type=str, default=None, help="Path to promoter aggregated by amp table")
    parser.add_argument("--binding_amp", type=str, default=None, help="Path to binding model aggregated by amp table")
    parser.add_argument("--promoter_samp", type=str, default=None, help="Path to promoter aggregated by sample table")
    parser.add_argument("--binding_samp", type=str, default=None, help="Path to binding model aggregated by sample table")
    parser.add_argument("--potency", type=str, default=None, help="Potency table")
    parser.add_argument("--partition_function", type=str, default=None, help="Partition function table")
    parser.add_argument("--output_dir", type=str, help="Output director for file")
    parser.add_argument("--out_prefix", type=str, help="Output prefix for file")

    args = parser.parse_args()

    # aggregate single molecule
    promoter_single_molecule = pd.read_table(args.promoter_single_molecule, index_col=0)
    binding_single_molecule = pd.read_table(args.binding_single_molecule, index_col=0)

    merged_single_molecule = promoter_single_molecule.join(binding_single_molecule, rsuffix='_')
    outpath = path.join(args.output_dir, f'{args.out_prefix}.single_molecule.txt.gz')
    merged_single_molecule.to_csv(outpath, sep='\t', index=True, header=True)

    # aggregate amplicon
    promoter_amp_level = pd.read_table(args.promoter_amp).set_index(['sample', 'amplicon']).sort_index()
    binding_amp_level = pd.read_table(args.binding_amp).set_index(['sample', 'amplicon']).sort_index()

    merged_amp_level = promoter_amp_level.join(binding_amp_level, rsuffix='_')
    outpath = path.join(args.output_dir, f'{args.out_prefix}.amplicon_level.txt.gz')
    merged_amp_level.to_csv(outpath, sep='\t', index=True, header=True)

    # aggregate sample
    promoter_samp_level = pd.read_table(args.promoter_samp).set_index('sample')
    binding_samp_level = pd.read_table(args.binding_samp).set_index('sample')
    potency = pd.read_table(args.potency).set_index('sample')
    partition_function = pd.read_table(args.partition_function).set_index('sample')

    merged_samp_level = promoter_samp_level.join([binding_samp_level, potency, partition_function], how='outer')
    outpath = path.join(args.output_dir, f'{args.out_prefix}.sample_level.txt.gz')
    merged_samp_level.to_csv(outpath, sep='\t', index=True, header=True)



