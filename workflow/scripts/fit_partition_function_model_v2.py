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


# ### BD 230606
# Finally moving the partition function model into a script
# This is for a few reasons:
# 1. We want to be able to run this reproducibly and avoid any of the weird notebook issues with saving variables
# 2. We want to be able to get other people to run it (which is hard in a notebook)
# 3. We want to make it easy to be able to iterate on model designs, have separate scripts (or functions/arguments) for 2 vs 3 parameter models 
# 4. We want to streamline the procedure and make sure we can do all the downstream stuff (on the simulated molecules) in a principled format 

# Outstanding questions are:
# 1. When we reparametrize the binding model outputs, we'll also need to change some functions here, make that doable
# 2. How are we dealing with nucs "off the edge" in Will speak? Make this an option to customize
# 3. What are the outputs we need from the partition function model itself and what are checks we do after (in a different script)
# 4. How will we deal with the different models where we want to fix certain parameters and fit across multiple conditions and stuff?
# 5. How will we accommodate things like sensitivity analysis? Also how to get error bars on the fits?

fixed_e_bounds = (200,600)

def remove_fake_states(df):
    '''
    I guess it's still an open question whether we want to exclude the molecules that are all 1s... This lets us do that
    '''
    return df.loc[df.idx!=0].copy()

def load_data(basedir, sample, amplicons, model_positions_file, filter_states=False):
    '''
    Reads in single molecule classifications and valid states from the binding model and loads into an object to pass to the partition function model
    Returns a list of (classified single molecules, all possible states) tuples ((df,df)) for all the input amplicons
    TODO: Figure out how to do this across experiments to allow info sharing?
    '''

    # read paths to valid states and sm classifications
    valid_paths = [path.join(basedir,sample,'{}.{}.valid_states_table.txt'.format(sample, amplicon)) for amplicon in amplicons]
    classification_paths = [path.join(basedir,sample,'{}.{}.single_molecule_classification.txt'.format(sample, amplicon)) for amplicon in amplicons]

    # states = possible states, assignments = observed states
    states_list = [pd.read_table(valid_path) for valid_path in valid_paths]
    assignments_list = [pd.read_table(classification_path) for classification_path in classification_paths]

    # optionally remove the all 1s state (and others?)
    if filter_states:
        assignments_list = [remove_fake_states(df) for df in assignments_list]
        states_list = [remove_fake_states(df) for df in states_list]

    # have to mess around with some nas here, I think
    assignments_list_ = []
    for idx in range(len(assignments_list)):
        df = assignments_list[idx].copy()
        cols = df.filter(like='present').columns
        for c in cols:
            df[c] = df[c].fillna(False)
        assignments_list_.append(df)
        
    states_list_ = []
    for idx in range(len(states_list)):
        df = states_list[idx].copy()
        cols = df.filter(like='present').columns
        for c in cols:
            df[c] = df[c].fillna(False)
        states_list_.append(df)

    data = list(zip(assignments_list_,states_list_))

    model_positions = load_tfbs_positions(model_positions_file)

    tf_positions = [model_positions[a] for a in amplicons]

    return (data, tf_positions)

def assign_energy_three_param(mat, tf_e, nuc_e, delta_e, enhancer_bounds):
    '''
    Function for assigning energies to each state in a data frame
    This is our "3 parameter" model, which has a TF-dependent nucleosome energy
    Here, "delta_e" is a discounting factor on the nuc energy (0 is no difference with TF, 1 is NO nuc energy when TF bound, etc)
    '''
    num_tfs = mat.filter(like='tfbs_').sum(axis=1)
    num_nuc = pd.Series(0, index=mat.index)
    
    for idx_nuc in range(len(mat.filter(like='nuc').columns.tolist())//3):
        # iterate through each nucleosome
        tmp_nuc_num = pd.Series(0, index=mat.index)
        tmp_nuc_num = tmp_nuc_num.mask((mat['nuc{}_present'.format(idx_nuc+1)]==True) & (((mat['nuc{}_start'.format(idx_nuc+1)]+mat['nuc{}_end'.format(idx_nuc+1)])/2).between(*enhancer_bounds)), 1)
        num_nuc += tmp_nuc_num.copy()
        
    # return num_tfs * tf_e + num_nuc * nuc_e * (1 - delta_e * (num_tfs>0))
    return num_tfs * tf_e + num_nuc * (nuc_e - delta_e * (num_tfs>0))


def assign_energy_two_param(mat, tf_e, nuc_e, enhancer_bounds):
    '''
    Function for assigning energies to each state in a data frame
    This is our "3 parameter" model, which has a TF-dependent nucleosome energy
    Here, "delta_e" is a discounting factor on the nuc energy (0 is no difference with TF, 1 is NO nuc energy when TF bound, etc)
    '''
    num_tfs = mat.filter(like='tfbs_').sum(axis=1)
    num_nuc = pd.Series(0, index=mat.index)
    
    for idx_nuc in range(len(mat.filter(like='nuc').columns.tolist())//3):
        # iterate through each nucleosome
        tmp_nuc_num = pd.Series(0, index=mat.index)
        tmp_nuc_num = tmp_nuc_num.mask((mat['nuc{}_present'.format(idx_nuc+1)]==True) & (((mat['nuc{}_start'.format(idx_nuc+1)]+mat['nuc{}_end'.format(idx_nuc+1)])/2).between(*enhancer_bounds)), 1)
        num_nuc += tmp_nuc_num.copy()
        
    # return num_tfs * tf_e + num_nuc * nuc_e * (1 - delta_e * (num_tfs>0))
    return num_tfs * tf_e + num_nuc * nuc_e

def convert_to_probabilities(energies, kt=1):
    '''
    Function to convert energies to probabilities using Boltzmann
    '''
    unscaled_probs = np.exp(-1*energies.astype(np.float64)/kt)
    return pd.DataFrame(unscaled_probs / unscaled_probs.sum(), index=energies.index, columns=['prob'])

def neg_log_like_multiple(x, *args):
    '''
    Likelihood function for computing how likely we were to see the data (observed state calls) given a specific energy (nuc, tf)
    x is (tf_e, nuc_e, (delta_e))
    # args is [[(molecules_1, states_1),...]]
    '''
    tf_e = x[0]
    nuc_e = x[1]
    delta_e = x[2]
    
    nll = 0
    
    n_enhancers = len(args[0])
    tf_positions = args[1]
    
    for idx in range(n_enhancers):
        molecules, states = args[0][idx]
        
        pos = tf_positions[idx]
        # pos = tf_positions[amplicons[idx]]
        
        # determine enhancer_bounds
        # basically, for "enhancers" where the TFBS only take up half of the read, do we want to count DNA that is "flanking" the enhancer
        # I literally have no idea what is best to do here...
        if pos:
            enhancer_bounds = (pos[0][0]-50, pos[-1][1]+100)
        else:
            enhancer_bounds = (200,300)
            
#         enhancer_bounds = (240,550)
        enhancer_bounds = fixed_e_bounds
            
        # compute prob per state
        prob_per_state = convert_to_probabilities(assign_energy_three_param(states, tf_e, nuc_e, delta_e, enhancer_bounds))

        # get the combined likelihood for all observed states
        data_probs = prob_per_state.loc[molecules['idx'],'prob']

        # log and sum 
        # make each n_tfbs count the same?
        nll += -np.log(data_probs).sum()
    
    return nll

def neg_log_like_multiple_two_param(x, *args):
    '''
    Likelihood function for computing how likely we were to see the data (observed state calls) given a specific energy (nuc, tf)
    x is (tf_e, nuc_e)
    # args is [[(molecules_1, states_1),...]]
    '''
    tf_e = x[0]
    nuc_e = x[1]
    
    nll = 0
    
    n_enhancers = len(args[0])
    tf_positions = args[1]
    
    for idx in range(n_enhancers):
        molecules, states = args[0][idx]
        
        pos = tf_positions[idx]
        # pos = tf_positions[amplicons[idx]]
        
        # determine enhancer_bounds
        # basically, for "enhancers" where the TFBS only take up half of the read, do we want to count DNA that is "flanking" the enhancer
        # I literally have no idea what is best to do here...
        if pos:
            enhancer_bounds = (pos[0][0]-50, pos[-1][1]+100)
        else:
            enhancer_bounds = (200,300)
            
#         enhancer_bounds = (240,550)
        enhancer_bounds = fixed_e_bounds
            
        # compute prob per state
        prob_per_state = convert_to_probabilities(assign_energy_two_param(states, tf_e, nuc_e, enhancer_bounds))

        # get the combined likelihood for all observed states
        data_probs = prob_per_state.loc[molecules['idx'],'prob']

        # log and sum 
        # make each n_tfbs count the same?
        nll += -np.log(data_probs).sum()
    
    return nll

def fit_model(input_data, num_params):
    '''
    Use scipy's minimization routine to get the energies that minimize the (negative) log-likelihood
    '''
    if num_params == 2:
        starting_values = (-1, -1)
        res = minimize(neg_log_like_multiple_two_param, (1,1), args=input_data, method='Nelder-Mead')
    else:
        starting_values = (-1, -1, 0)
        res = minimize(neg_log_like_multiple, (1,1,0), args=input_data, method='Nelder-Mead')

    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit partition function model on single molecule data')
    parser.add_argument("--basedir", dest="basedir", type=str, help="Path to location where the classification pipeline was run")
    parser.add_argument("--sample", dest="sample", type=str, help="Sample name on which the binding model was run")
    parser.add_argument("--amplicons", dest="amplicons", type=str, help="Comma-separated list of amplicons to fit on")
    parser.add_argument("--output", dest="output", type=str, help="Path to outputting the parameters of the partition function model")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file (optional)")
    parser.add_argument("--model_positions", dest="model_positions", type=str, default="/oak/stanford/groups/wjg/mhinks/projects/smf/220919_P028_doxdose_opJS45/amplicon-smf/amplicon-info/opJS45.positions.long.txt", help="Path to positions file for partition function script")
    parser.add_argument("--num_params", dest="num_params", type=int, default=3, help="Number of paremeters in fit")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots, open(args.output, 'w') as output:
        # load data
        data = load_data(args.basedir, args.sample, args.amplicons.split(','), args.model_positions)

        # fit model
        mle_params = fit_model(data, args.num_params)

        # write to output file
        # extract args
        if args.num_params == 2:
            tf_e, nuc_e = mle_params.x[0], mle_params.x[1]
            output.write('E_TF\t{}\n'.format(tf_e))
            output.write('E_nuc\t{}\n'.format(nuc_e))
        else:
            tf_e, nuc_e, delta_e = mle_params.x[0], mle_params.x[1], mle_params.x[2]
            output.write('E_TF\t{}\n'.format(tf_e))
            output.write('E_nuc\t{}\n'.format(nuc_e))
            output.write('deltaE_nuc_when_TF_bound\t{}\n'.format(delta_e))
