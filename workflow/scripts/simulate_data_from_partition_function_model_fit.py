import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy.stats import pearsonr, entropy
from scipy.optimize import minimize
import os
from os import path
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import fcluster, fclusterdata, linkage, dendrogram 
from copy import copy, deepcopy
import random
import statsmodels.api as sm
import math
from scipy import interpolate, ndimage, fft
import pickle
from scipy.optimize import curve_fit
from scipy.special import comb
import argparse
import warnings

from common import load_single_molecule_matrix, \
                   get_methyl_positions, \
                   load_tfbs_positions, \
                   filter_all_converted_reads, \
                   decorate_single_read_plot2, \
                   plot_single_read, \
                   adjust_gcgs, \
                   fix_missing_data

import fit_partition_function_model_v3

# # BD 231114
# Goal is to simplify the generation of simulated data
# We'll do it with a script, where we pass in the results of the partition function fit

# TO DO:
# - Make this for other models


def load_parameters(output_path, model):
    model_outputs = pd.read_table(output_path, header=None, index_col=0)

    if model == '2param':
        tf_e, nuc_e = model_outputs.loc['E_TF',1], model_outputs.loc['E_nuc',1]
        return (tf_e, nuc_e)
    elif model == '3param_nuc':
        tf_e, nuc_e, delta_e = model_outputs.loc['E_TF',1], model_outputs.loc['E_nuc',1], model_outputs.loc['deltaE_nuc_when_TF_bound',1]
        return (tf_e, nuc_e, delta_e)
    elif model == '3param_tfcoop':
        tf_e, nuc_e, delta_e = model_outputs.loc['E_TF',1], model_outputs.loc['E_nuc',1], model_outputs.loc['deltaE_tf_when_TF_bound',1]
        return (tf_e, nuc_e, delta_e)
    else:
        print('Not a valid model this should not happen')
        assert False

    # return (tf_e, nuc_e, delta_e)

def simulate_data(data, params, amplicons, model, n_samples=10000):
    '''
    Use input parameters to simulate data
    '''

    simulated_data = {}

    for idx, amplicon in enumerate(amplicons):
        molecules, states = data[idx]

        if model == '2param':
            tf_e, nuc_e = params
            probs = fit_partition_function_model_v3.convert_to_probabilities(fit_partition_function_model_v3.assign_energy_two_param(states, tf_e, nuc_e, fit_partition_function_model_v3.fixed_e_bounds))
        elif model == '3param_nuc':
            tf_e, nuc_e, delta_e = params
            probs = fit_partition_function_model_v3.convert_to_probabilities(fit_partition_function_model_v3.assign_energy_three_param(states, tf_e, nuc_e, delta_e, fit_partition_function_model_v3.fixed_e_bounds))
        elif model == '3param_tfcoop':
            tf_e, nuc_e, delta_e = params
            probs = fit_partition_function_model_v3.convert_to_probabilities(fit_partition_function_model_v3.assign_energy_three_param_tfcooperativity(states, tf_e, nuc_e, delta_e, fit_partition_function_model_v3.fixed_e_bounds))
        else:
            print('Not a valid model this should not happen')
            assert False

        samp = np.random.choice(states['idx'], p=probs['prob'].values, size=n_samples)
        simulated_data[amplicon] = states.loc[samp].copy()

    return simulated_data

def compute_summary_stats(simulated_data):
    summary_stats = []

    for amplicon, df in simulated_data.items():
        n_tfbs = int(amplicon.split('_')[1][0]) if 'TetO' in amplicon else 0

        accessibility = (df.filter(like='tfbs_').sum(axis=1)>0).mean()
        avg_tf = df.filter(like='tfbs_').sum(axis=1).mean()
        tf_hist = ','.join(df.filter(like='tfbs').sum(axis=1).value_counts(normalize=True).sort_values('index', ascending=False).reset_index().apply(lambda row: '{:.0f}:{:.4f}'.format(row['index'], row[0]), axis=1).to_list())

        summary_stats.append({'amplicon': amplicon, 'n_tfbs': n_tfbs, 'avg_tf_bound': avg_tf, 'fraction_enhancers_tf_bound': accessibility, 'binding_histogram': tf_hist})

    return pd.DataFrame(summary_stats)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate simulated molecules from the partition function')
    parser.add_argument("--basedir", dest="basedir", type=str, help="Path to location where the classification pipeline was run")
    parser.add_argument("--sample", dest="sample", type=str, help="Sample name on which the binding model was run")
    parser.add_argument("--amplicons", dest="amplicons", type=str, help="Comma-separated list of amplicons to fit on")
    parser.add_argument("--params", dest="param_path", type=str, help="Path to partition function model parameters")
    parser.add_argument("--output", dest="output", type=str, help="Path to outputting the summary stats df")
    parser.add_argument("--model_positions", dest="model_positions", type=str, default="/oak/stanford/groups/wjg/mhinks/projects/smf/220919_P028_doxdose_opJS45/amplicon-smf/amplicon-info/opJS45.positions.long.txt", help="Path to positions file for partition function script")
    parser.add_argument("--model", dest="model", type=str, default="3param_nuc", help="Which flavor of the model to fit (options are '3param_nuc', '2param', '3param_tfcoop', and '4param')")

    args = parser.parse_args()

    # load data
    data, tf_positions = fit_partition_function_model_v3.load_data(args.basedir, args.sample, args.amplicons.split(','), args.model_positions)

    # valid_models = ['3param_nuc', '2param', '3param_tfcoop', '4param']
    valid_models = ['3param_nuc', '2param', '3param_tfcoop']

    if args.model in valid_models:
        # load params
        params = load_parameters(args.param_path, args.model)
    else:
        print('Invalid model string, please choose from: {}'.format(','.join(valid_models)))
        assert False

    # generate simulated data
    simulated_data = simulate_data(data, params, args.amplicons.split(','), args.model)

    # generate summary stats
    summary_stats = compute_summary_stats(simulated_data)
    
    # write to output file
    summary_stats.to_csv(args.output, sep='\t', header=True, index=False)
