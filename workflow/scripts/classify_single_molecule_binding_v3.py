import pandas as pd
import numpy as np
import argparse
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
from itertools import product, chain, combinations
from scipy.stats import bernoulli
from copy import copy
from os import path
import pickle
from sklearn.linear_model import LogisticRegression
import warnings


# ignore divide by zero error (intended behavior) in the broadcasting function
# np.seterr(divide='ignore')
warnings.filterwarnings('ignore')



from common import load_single_molecule_matrix, \
                   get_methyl_positions, \
                   load_tfbs_positions, \
                   filter_all_converted_reads, \
                   decorate_single_read_plot2, \
                   plot_single_read, \
                   adjust_gcgs, \
                   fix_missing_data

# re-writing the classification model
# this time we will just do a maximum likelihood approach over ALL possible underlying sequences
# use some np vectorization tricks from Ben P. to speed everything up computations
# BD big changes 230329 (from v2):
# first of all, we want to separate the prob of methylation from the prob of conversion, conversion is fixed by experiment, and it's worth having this as an explicitly separate parameter
# what this also allows us to do is then have (state, GpC)-specific _methylation_ probabilities, as a way to both incorporate biology and penalize states we don't like
# for instance, we can penalize getting nuc wrong over the tetos by saying that you'll _never_ methylate at the center of a nucleosome (but you can at edges)
# this might require re-writing how we parametrize nucleosomes and stuff (see below), but this is both theoretically correct and biologically correct
# ideally, once we have parametrized things in a new way, we can also add in an EM-step to estimate the methylation (and conversion?) probs from data
# big change to nuc and TF parametrization is that methylation probability depends on "distance from the center of the motif" (as a logistic regression problem)
# this will presumably also allow us to more easily accomodate the rTetR only (by changing p(methylation))

# to do list
# add plot of nuc (and TF?) LR to gague whether it learned the right thing
# add likelihood of each molecule under the (finalized) model
# allow for filtering based on this likelihood
# potentially make it easier to play around with specific molecules?

def compute_conversion_hyperparameters():
    '''
    currently estimating this from what I expect things to be
    but will build in the capabilities for EM algorithm
    '''

    # conversion probabilities
    # t is the thing we measure as a 1 in the data (i.e. converted = unmethylated = bound = protected)
    # note that these, for now, are fixed across experiments, but ideally we can get experiment specific values from the em-seq
    p_t_given_meth = 0.15 # 0.05
    p_t_given_unmeth = 0.95 # 0.99

    # note, for now we are ignoring the p_methylations (since these will be computed below)

    return (p_t_given_meth, p_t_given_unmeth)

def return_global_amplicon_properties():
    '''
    define the constants for amplicon length and stuff we'll use for the rest of the code
    mostly helpful for translating between real (GpC) and binned (nuc midpoint) space
    '''

    # length of the amplicon from F to R primer (this one probably changes the most, might want command line arg)
    amp_width = 600

    # length of a nucleosome (every underlying nuc is the same length)
    nuc_length = 140

    # bin size resolution for nucleosome midpoints (and therefore starting positions)
    bin_size = 10

    return (amp_width, nuc_length, bin_size)

# global parameters 
amp_width, nuc_length, bin_size = return_global_amplicon_properties()

def recursive_1_placer(lst, spacing, start_idx):
    '''
    places 1s into an empty array of 0s subject to spacing constraints between adjacent 1s
    helper function for below, recursive
    '''
    to_return = [copy(lst)]
    if start_idx < len(lst):
        for idx in range(start_idx, len(lst)):
            new_lst = copy(lst)
            assert new_lst[idx] == 0
            new_lst[idx] = 1
            to_return += recursive_1_placer(copy(new_lst), spacing, idx+spacing)
    return to_return        

def enumerate_nucleosomal_states(amp_width, nuc_length, bin_size):
    '''
    puts nucleosome midpoints (1s) into an array of bins of variable size
    this is how we are going to start parametrizing nucleosomes, for a variety of reasons
    '''
    # figure out nuc_width in bin terms
    nuc_bin_width = nuc_length // bin_size
    
    # figure out length
    input_lst = [0 for _ in range(amp_width // bin_size)]

    nucs = recursive_1_placer(input_lst, nuc_bin_width, 0)
    
    return pd.DataFrame(nucs, columns=list(range(len(input_lst))))

def powerset_generator(n):
    '''
    returns a list of lists of all possible binary strings of length n
    '''
    if n == 0:
        return [[]]
    else:
        return [[0]+x for x in powerset_generator(n-1)] + [[1]+x for x in powerset_generator(n-1)]

def enumerate_tf_states(tfbs_positions):
    '''
    returns a pandas dataframe of all possible TF binding configurations
    '''
    return pd.DataFrame(powerset_generator(len(tfbs_positions)), columns=list(range(len(tfbs_positions))))

# def find_first_matching_sublist(lst, span, start):
#     '''
#     helper function for below
#     basically, given a start position, figure out the first subsequent position that satisfies the run being in the range
#     returns the index of this position, or else returns -1 if it falls off the end
#     '''
#     min_span, max_span = span
    
#     end = start
    
#     while end < len(lst):
#         temp_span = lst[end] - lst[start] + 1 # fencepost
#         if temp_span >= min_span and temp_span <= max_span:
#             return end
#         else:
#             end += 1

#     return -1

# def sublist_finder(lst, span, start_ptr):
#     '''
#     generic function to compute all possible valid runs of a given length from a list of integers
#     input is a lst of positive integers (no duplicates, sorted)
#     also a range (low and high for how long the run can be, in the units of the entries of the list)
#     also a start (since this is a recursive function)
#     each valid configuration is a list of tuples of (start,end) of valid runs (so a single output could be [] or [(1,150)] or [(1,150),(161,300)])
#     output is a list of these lists
#     '''
#     valid_substrings = []

#     start = start_ptr
    
#     while start < len(lst):
#         # from start, find a valid range that works
#         end = find_first_matching_sublist(lst, span, start)
#         if end != -1:
#             # save this tuple
#             temp_sublist = [(lst[start],lst[end])]
#             # recursive call, starting from one past the end of the current
#             # note, can make this end+2 to ensure that the model inserts one accessible base between nucs for spacer, not sure if necessary
#             recursive_valid_sublists = sublist_finder(lst, span, end+1) 
#             # append the empty run (so that we can get the "just this on run" state)
#             # probably not the best way to do it but whatever
#             recursive_valid_sublists += [[]]
#             # add the current run to every recursive one 
#             for sublist in recursive_valid_sublists:
#                 valid_substrings.append(temp_sublist + sublist)
#         # now start at the next element in the sequnce
#         start += 1

#     # only once, append the empty run
#     if start_ptr == 0:
#         valid_substrings += [[]]
        
#     return valid_substrings

# def powerset(iterable):
#     '''
#     computes the powerset of an iterable
#     stolen from SO
#     '''
#     s = list(iterable)
#     return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# def augment_methyl_positions(methyl_positions, step_size=10, overhang=30, nuc_length=130):
#     '''
#     adds extra gpc positions to the end of the list (i.e. not observed ones, imaginary)
#     this is to allow nucs to slide off the end and still be valid
#     can modify how many and by how much a nucleosome has to overlap the molecule
#     220911 modified to allow for nucs to slide off "front" end too, which allows us to do +1 nuc
#     '''
#     start = methyl_positions[0]
#     end = methyl_positions[-1]
#     return list(range(start-nuc_length+overhang,start-step_size,step_size)) + methyl_positions + list(range(end+step_size, end+nuc_length-overhang, step_size))

def get_nuc_start_and_end_from_midpoint(bin_idx, bin_size, nuc_length):
    '''
    Helper function used a couple of times, given which bin the nuc's center falls in, return tuple of nuc start and nuc end (given fixed width nucleosome)
    '''
    center = bin_idx * bin_size + bin_size // 2
    bounds = center - nuc_length // 2, center + nuc_length // 2

    return bounds


def construct_nuc_positions_array(nuc_states, amp_len, nuc_length, bin_size):
    '''
    
    '''

    # figure out nuc_width in bin terms
    nuc_bin_width = nuc_length // bin_size

    to_return = {}

    for idx in nuc_states.index:
        state = nuc_states.loc[idx]

        # initialize empty array
        tmp = np.zeros(amp_len)

        # iterate through nuc_positions
        for i in state.index:
            # check if this state is bound
            if state[i] > 0:
                bounds = get_nuc_start_and_end_from_midpoint(i, bin_size, nuc_length)
                # add all covered bases to list
                for pos in range(*bounds):
                    if pos >= 0 and pos < amp_len:
                        tmp[pos] = 1

        # add to dictionary
        to_return[idx] = tmp
    
    return to_return

def construct_tf_positions_array(tf_states, tfbs_positions, amp_len=600):
    '''

    '''
    to_return = {}

    # iterate through all states
    for idx in tf_states.index:
        state = tf_states.loc[idx]

        # initialize empty array
        tmp = np.zeros(amp_len)

        # iterate through tfbs_positions
        for i in range(len(tfbs_positions)):
            # check if this state is bound
            if state[i] > 0:
                bounds = tfbs_positions[i][0], tfbs_positions[i][1]
                # add all covered bases to list
                for pos in range(*bounds):
                    tmp[pos] = 1

        # add to dictionary
        to_return[idx] = tmp

    return to_return

def pair_is_valid(tf_state, nuc_state):
    '''
    checks whether a pair of nucleosome and tf state is "valid" (in the sense that the TFs and nucs don't overlap)
    '''

    return (tf_state * nuc_state).sum() == 0

def prune_states(tf_states, nuc_states):
    '''
    iterate through all pairs of tf and nuc states, and retain the ones that are valid
    return a pair of dfs (one for nucs and one for tfs) with the satisfactory pairs
    need to figure out the TF adjacent to nuc thing in new implementation
    '''

    valid_states_tf = []
    valid_states_nuc = []

    # construct the mapping between TFs and nucs 
    nuc_dict = construct_nuc_positions_array(nuc_states, amp_len=amp_width, nuc_length=nuc_length, bin_size=bin_size)
    tf_dict = construct_tf_positions_array(tf_states, tfbs_positions, amp_len=amp_width)

    for idx_tf, idx_nuc in product(tf_states.index, nuc_states.index):
        tf_state = tf_states.loc[idx_tf]
        nuc_state = nuc_states.loc[idx_nuc]

        tf_occupancy_state = tf_dict[idx_tf]
        nuc_occupancy_state = nuc_dict[idx_nuc]

        if pair_is_valid(tf_occupancy_state, nuc_occupancy_state):
            valid_states_tf.append(tf_state)
            valid_states_nuc.append(nuc_state)

    return pd.DataFrame(valid_states_tf, index=list(range(len(valid_states_tf)))), pd.DataFrame(valid_states_nuc, index=list(range(len(valid_states_nuc))))

def enumerate_states(methyl_positions, tfbs_positions):
    '''
    generates the list of all possible underlying occupancy states
    first adds in a few "ghost" gpcs at the end of the list to allow nucs to slide off the edge
    then computes all possible nucleosome occupancy events
    then compute all possible tf binding events
    returns the outer product of the two lists
    '''

    # generate all individual possible states
    tf_states = enumerate_tf_states(tfbs_positions)
    nuc_states = enumerate_nucleosomal_states(amp_width, nuc_length, bin_size) # need to feed in here amplicon length and shit

    # trim states
    valid_tf_states, valid_nuc_states = prune_states(tf_states, nuc_states)

    return valid_tf_states, valid_nuc_states



    # augmented_methyl_positions = augment_methyl_positions(methyl_positions)
    # possible_nuc_positions = sublist_finder(augmented_methyl_positions, nuc_span, 0)
    # possible_binding_events = powerset(tfbs_positions)

    # return product(possible_nuc_positions, possible_binding_events)

# def nuc_overlaps_tf_binding(methyl_positions, state):
#     '''
#     check whether TF binding overlaps nucleosome binding
#     almost the same as generate_predicted_protection (below)
#     first fills in all gpcs that are covered by nucs
#     then figures out if any TFs are bound there too
#     also checks promoter state, too
#     '''
#     nuc_positions, binding_sites = state
#     protection = [0 for _ in range(len(methyl_positions))]
    
#     # do nucs
#     for idx, p in enumerate(methyl_positions):
#         for nuc in nuc_positions:
#             if p >= nuc[0] and p <= nuc[1]:
#                 protection[idx] = 1
            
#     # do tfbs
#     for idx, p in enumerate(methyl_positions):
#         for tfbs in binding_sites:
#             if p >= tfbs[0] and p <= tfbs[1]:
#                 if protection[idx] == 1: # checking to see if nucleosome overlaps a bound TF
#                     return True
#                 if (idx > 0 and protection[idx-1] == 1) or (idx+1 < len(methyl_positions) and protection[idx+1] == 1): # also checking if a nucleosome overlaps the _flanks_ of a bound tf, this way a nucleosome cannot abut a bound tf (ambiguous case, plus impossible?)
#                     return True
    
#     return False

# def prune_states(all_states, methyl_positions):
#     '''
#     remove (most of) the states where there is TF and nucleosome bound at the same position
#     '''
#     return [s for s in all_states if not nuc_overlaps_tf_binding(methyl_positions, s)]

# def generate_predicted_protection(methyl_positions, state):
#     '''
#     iterate through each gpc position, if it's bound (by either a nucleosome or a TF) mark the base protected
#     '''
#     nuc_positions, binding_sites = state
#     protection = [0 for _ in range(len(methyl_positions))]
    
#     # do nucs
#     for idx, p in enumerate(methyl_positions):
#         for nuc in nuc_positions:
#             if p >= nuc[0] and p <= nuc[1]:
#                 protection[idx] = 1
            
#     # do tfbs
#     for idx, p in enumerate(methyl_positions):
#         for tfbs in binding_sites:
#             if p >= tfbs[0] and p <= tfbs[1]:
#                 protection[idx] = 1

#     return protection

# def get_flanking_gcs(methyl_positions, tfbs, index=True):
#     '''
#     Returns a list (2 member) of the bases that flank a given TFBS
#     Defaults to returning index of the bases, but can be changed to return the positions themselves
#     '''
#     pos = []
#     for idx, p in enumerate(methyl_positions):
#         if p >= tfbs[0] and p <= tfbs[1]:
#             if index:
#                 pos.append(idx)
#             else:
#                 pos.append(p)
#     assert len(pos) == 2
#     return pos

# def logistic_regression(x, mu, s):
#     '''
#     temporary implementation of LR until I get sklearn's up and running (for actually fitting)
#     p = 1/(1+exp(-(x-mu)/s))
#     Figure out how to parametrize this, b0/b1 or mu/s (depends I guess on sklearn implementation?)
#     '''

#     return 1/(1+np.exp(-(x-mu)/s))


# def generate_unmethylation_prob_nucleosome(nuc_states, methyl_positions, prob_unmeth_given_nuc, mu=60, s=-1, nuc_length=140, bin_size=10):
def generate_unmethylation_prob_nucleosome(nuc_states, methyl_positions, nuc_lr, nuc_length=140, bin_size=10):
    '''
    Generates a matrix of shape (states,gpcs) where entry ij is the probability that given GC j in state i is unmethylated _due to nucleosome protection_
    We are going to construct the (stats,gpcs) matrix by multiplying the states x nucs (0/1) with a nuc x gpc matrix that is the prob of each gpc being protected by each nuc
    Change this to logistic regression (with params we can fit by EM?)  
    BD 230518: originally, thought we might want to memoize the bit where we compute the distance of each base to the nuc center so we don't have to recompute each time when doing EM
    I think this would only marginally speed things up though, since we only have to compute per each nuc, not each state
    '''

    num_gcs = len(methyl_positions)
    num_nucs = len(nuc_states.columns)

    methylation_given_single_nuc_mat = np.zeros((num_nucs, num_gcs))

    # can almost certainly first compute the distance matrix (using helper function below) and then just predict on the whole matrix
    for idx in range(num_nucs):
        center = idx * bin_size + bin_size // 2
        # methylation_given_single_nuc_mat[idx] = 1*(np.abs(methyl_positions-center) < nuc_length / 2)
        # methylation_given_single_nuc_mat[idx] = logistic_regression(np.abs(methyl_positions-center), nuc_lr_mu, nuc_lr_s) 
        methylation_given_single_nuc_mat[idx] = nuc_lr.predict_proba(np.abs(methyl_positions-center).reshape(-1,1))[:,1]

    # by multipling the single nuc matrix by the states by which single nucs they contain matrix, we can quickly generate the full states by probability matrix
    # BD note 230518: check that this never makes and p(unmeth) for a given base >1 (which might be possible under certain LR parameters if the range of two adjacent nucs overlap or something?)
    # or the LR function we use we can set to 0 past some distance manually (but then it's another parameter to use? can potentially do some other operation besides summing?)
    nuc_protections = nuc_states.values @ methylation_given_single_nuc_mat

    # nuc_protections = np.zeros((len(states), len(methyl_positions)))

    # for i, state in enumerate(states):
    #     nuc_positions, binding_sites = state

    #     protection = [0 for _ in range(len(methyl_positions))]

    #     for idx, p in enumerate(methyl_positions):
    #         for nuc in nuc_positions:
    #             if p >= nuc[0] and p <= nuc[1]:
    #                 protection[idx] = 1

    #     nuc_protections[i] = protection

    # return nuc_protections * prob_unmeth_given_nuc
    return nuc_protections

def generate_unmethylation_prob_tf(tf_states, methyl_positions, tfbs_positions, prob_unmeth_given_tf, tf_lr_mu, tf_lr_s):
    '''
    Generates a matrix of shape (states,gpcs) where entry ij is the probability that given GC j in state i is unmethylated _due to TF protection_
    '''

    num_gcs = len(methyl_positions)
    num_tfs = len(tf_states.columns)

    methylation_given_single_tf_mat = np.zeros((num_tfs, num_gcs))

    for idx in range(num_tfs):
        tf_start, tf_end = tfbs_positions[idx][0], tfbs_positions[idx][1]
        center = (tf_start + tf_end) // 2
        width = tf_start - center
        methylation_given_single_tf_mat[idx] = 1 * (methyl_positions < tf_end+2) * (methyl_positions > tf_start-2) * prob_unmeth_given_tf
        # methylation_given_single_tf_mat[idx] = logistic_regression(np.abs(methyl_positions-center), tf_lr_mu, tf_lr_s) 

    tf_protections = tf_states.values @ methylation_given_single_tf_mat 

    # tf_protections = np.zeros((len(states), len(methyl_positions)))

    # for i, state in enumerate(states):
    #     nuc_positions, binding_sites = state

    #     protection = [0 for _ in range(len(methyl_positions))]

    #     for idx, p in enumerate(methyl_positions):
    #         for tfbs in binding_sites:
    #             if p >= tfbs[0] and p <= tfbs[1]:
    #                 protection[idx] = 1

    #     tf_protections[i] = protection

    # return tf_protections * prob_unmeth_given_tf
    return tf_protections

def generate_unmethylation_prob_unbound(states, methyl_positions, promoter_positions, prob_unmeth_given_open, prob_unmeth_given_open_promoter=0.5):
    '''
    Generates a matrix of shape (states,gpcs) where entry ij is the probability that given GC j in state i is unmethylated _with nothing sitting there_
    NOTE: We need to go back and do the promoter stuff here!!! Like just set all the promoter columns to some other value! Easy
    '''
    open_protections = np.ones((len(states), len(methyl_positions))) * prob_unmeth_given_open

    # fix promoters
    promoter_low, promoter_high = promoter_positions
    promoter_indices = [i for i, p in enumerate(methyl_positions) if p >= promoter_low and p <= promoter_high]
    open_protections[:,promoter_indices] = prob_unmeth_given_open_promoter

    return open_protections


def create_bernoulli_logprob_matrices(tf_states, nuc_states, methyl_positions, tfbs_positions, em_hyperparams, p_t_given_meth, p_t_given_unmeth, promoter_positions=(0,250)):
    '''
    generates the p(T | state) and p(C | state) matrices that are used for the likelihood computation

    '''
    n_gpcs = len(methyl_positions)
    n_states = len(tf_states)

    # prob_unmeth_given_nuc = 0.99
    # prob_unmeth_given_tf = 0.99
    # prob_unmeth_given_open = 0.01

    # unpack hyperparams
    # nuc_lr_mu = em_hyperparams['nuc_lr_mu']
    # nuc_lr_s = em_hyperparams['nuc_lr_s']
    nuc_lr = em_hyperparams['nuc_lr']
    prob_unmeth_given_tf = em_hyperparams['prob_unmeth_given_tf']
    # tf_lr_mu = em_hyperparams['tf_lr_mu']
    # tf_lr_s = em_hyperparams['tf_lr_s']
    prob_unmeth_given_open = em_hyperparams['prob_unmeth_given_open']

    # get the unmethylated matrices for nucs, tfs, and other
    prob_unmethylated_nuc = generate_unmethylation_prob_nucleosome(nuc_states, methyl_positions, nuc_lr)
    prob_unmethylated_tf = generate_unmethylation_prob_tf(tf_states, methyl_positions, tfbs_positions, prob_unmeth_given_tf, 0, 0)
    prob_unmethylated_open = generate_unmethylation_prob_unbound(nuc_states, methyl_positions, promoter_positions, prob_unmeth_given_open)

    # join together
    prob_unmethylated = np.max(np.stack([prob_unmethylated_nuc, prob_unmethylated_tf, prob_unmethylated_open], axis=0), axis=0)

    # make complimentary matrix
    prob_methylated = 1 - prob_unmethylated

    # now turn this into a prob of seeing conversion (T = 1 = unmethylated)
    prob_t = prob_unmethylated * p_t_given_unmeth + prob_methylated * p_t_given_meth
    prob_c = prob_unmethylated * (1-p_t_given_unmeth) + prob_methylated * (1-p_t_given_meth)
        
    log_prob_states_t = np.log(prob_t)
    log_prob_states_c = np.log(prob_c)

    return (log_prob_states_t, log_prob_states_c)

def determine_chunk_size(n_mols, n_states, data_size, memory_threshold):
    '''
    Figures out how many times we need to split the data into (along the single molecule axis) to not go above memory limits
    '''
    return int((n_mols * n_states * data_size) // (memory_threshold * 1e9)) + 1

def classify_all_molecules(single_molecules, log_prob_states_t, log_prob_states_c, n_matrix_partitions):
    '''
    Compute Bernoulli likelihood all at once 
    note that log_prob_states_unmeth and meth and not strictly reciprocals of each other, since we need to change unoccupied promoter probs (see above)
    single_molecules is gpcs x molecules
    log_prob_states_(un)meth are states x gpcs
    result of matmul is states x molecules, argmax down columns
    return a list of index of best state per molecule
    BD 221213: rewriting to split this into chunks, for a smaller number of single molecules, so as not to exceed memory limits
    BD 230707: adding in the option to weight here by a prior (to make this bayesian and to provide info about how likely we think certain states are)
    '''
    n_mols = single_molecules.shape[1]
    chunk_size = n_mols // n_matrix_partitions + 1

    to_return = []

    for idx in range(n_matrix_partitions):
        sub_mat = single_molecules[:,idx*chunk_size:(idx+1)*chunk_size]
        res = np.argmax(np.exp(log_prob_states_t @ sub_mat + log_prob_states_c @ (1 - sub_mat)), axis=0)
        to_return += list(res)

    return to_return

def build_distance_to_nuc_matrix(nuc_states, methyl_positions, nuc_length=140, bin_size=10):
    '''
    Generates a matrix of shape (states,gpcs) where entry ij is the distance of GpC j to the closest nucleosome center in state i
    BenP had lots of advice for the np.newaxis trick used to make this possible (and there are more tricks to speed up in slack if needed)
    We will use this a bunch to both fit and apply the LR methylation model, and we can memoize these if needed (I think not?)
    '''

    num_gcs = len(methyl_positions)
    num_nucs = len(nuc_states.columns)

    methylation_given_single_nuc_mat = np.zeros((num_nucs, num_gcs))

    for idx in range(num_nucs):
        center = idx * bin_size + bin_size // 2
        # methylation_given_single_nuc_mat[idx] = 1*(np.abs(methyl_positions-center) < nuc_length / 2)
        methylation_given_single_nuc_mat[idx] = np.abs(methyl_positions - center)

    # BenP trick for matmuls with other operation besides summing over the interior axis
    # basically this takes the binary encoding of states and "multiplies" with the per state distance, but instead of summing over the internal axis we take the min (along the middle shared dimension)
    # this works due to np broadcasting rules
    inverted_nuc_states = (1 / nuc_states.values) - 1 # 1->0 and 0->np.inf

    # print(inverted_nuc_states.values)
    # print(methylation_given_single_nuc_mat)
    # print(inverted_nuc_states.values[:,:,np.newaxis].shape)
    # print(methylation_given_single_nuc_mat[np.newaxis,:,:])

    # basically, we add axes such that they are the same dims (mols,nuc_idxs,gcps), add infinity to any of the distances that not present, and min along the nuc_idx axis
    distance_to_nuc = np.min(inverted_nuc_states[:,:,np.newaxis] + methylation_given_single_nuc_mat[np.newaxis,:,:], axis=1)

    return distance_to_nuc

def build_distance_to_tf_matrix(tf_states, methyl_positions, tfbs_positions):
    '''
    Generates a matrix of shape (states,gpcs) where entry ij is the distance of GpC j to the closest TF center in state i
    For now we are just going to have 1 (if inside some fixed width) and 0 (if outside), but we should probably turn it into the LR thing as above
    '''

    num_gcs = len(methyl_positions)
    num_tfs = len(tf_states.columns)

    methylation_given_single_tf_mat = np.zeros((num_tfs, num_gcs))

    for idx in range(num_tfs):
        tf_start, tf_end = tfbs_positions[idx][0], tfbs_positions[idx][1]
        center = (tf_start + tf_end) // 2
        width = tf_start - center
        methylation_given_single_tf_mat[idx] = 1 * (methyl_positions < tf_end+2) * (methyl_positions > tf_start-2)
        # methylation_given_single_tf_mat[idx] = logistic_regression(np.abs(methyl_positions-center), tf_lr_mu, tf_lr_s) 

    tf_protections = tf_states.values @ methylation_given_single_tf_mat 

    return tf_protections

def get_initial_em_hyperparams():
    '''
    Function that returns initial conditions for EM. If we are doing EM great we can update, if not we'll just use these
    '''
    nuc_lr_mu = 60
    nuc_lr_s = -1
    prob_unmeth_given_tf = 0.9
    tf_lr_mu = 10 
    tf_lr_s = -1
    prob_unmeth_given_open = 0.05

    # tmp: instantiate initial regressions based on fake data
    # redo this
    X_tmp = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90]).reshape(-1,1)
    y_tmp = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0])
    nuc_lr = LogisticRegression().fit(X_tmp, y_tmp)

    em_hyperparams = {'nuc_lr_mu': nuc_lr_mu, 
                      'nuc_lr_s': nuc_lr_s, 
                      'nuc_lr': nuc_lr, 
                      'prob_unmeth_given_tf': prob_unmeth_given_tf, 
                      'tf_lr_mu': tf_lr_mu, 
                      'tf_lr_s': tf_lr_s, 
                      'prob_unmeth_given_open': prob_unmeth_given_open}

    return em_hyperparams

def estimate_latent_params_from_state_assignments(single_molecules, methyl_positions, tfbs_positions, maximum_likelihood_assignments, tf_states, nuc_states, rep, nuc_distance_threshold=140, promoter_positions=(0,250)):
    '''
    this is the E (I think??) of the EM step. basically, what we've done is we have done MLE (M??) to get the thetas (states) for each molecule
    now we need to go back and estimate the latent parameters (e.g. the things we use to compute the p_unmeth_given_state matrices for nuc, tf, open)
    basically, we will first look at all "unbound" GpCs and compute how many of them are methylated (we'll have to think about promoters here, somehow)
    then we can take some width around each nuc center and then _fit_ a LR on all the bases within that window 0/1 as a f(distance) (to get those params)
    then we can do the same for the TFBSes, either LR or just one param
    BD 230523: THINK THIS WON'T WORK UNTIL I ADD A STEP THAT REINDEXES POSSIBLE STATES TO OBSERVED STATES GIVEN THE MLE ASSIGNMENTS AND THE STATE LIST!
    '''

    # first, compute dist of each GpC in each state to closest nuc/TF
    
    # compute distance to nearest nuc
    distance_to_nuc_states = build_distance_to_nuc_matrix(nuc_states, methyl_positions)
    np.savetxt('nuc_states.{}.txt'.format(rep), distance_to_nuc_states)

    # compute distance to TFBS center (or, now, whether base is in TFBS or not)
    distance_to_tf_states = build_distance_to_tf_matrix(tf_states, methyl_positions, tfbs_positions)
    np.savetxt('tf_states.{}.txt'.format(rep), distance_to_tf_states)


    # reindex to get observed states, not all states
    distance_to_nuc = distance_to_nuc_states[maximum_likelihood_assignments,:]
    distance_to_tf = distance_to_tf_states[maximum_likelihood_assignments,:]

    # now, we estimate nuc
    # threshold on distance (so we don't get confused by "far away" guys that are actually bound e.g. by TF) and remove TF bound dudes
    distance_to_nuc_array = distance_to_nuc[(distance_to_nuc < nuc_distance_threshold) & (distance_to_tf == 0)]

    # apply same filter to data
    observed_data_nuc = single_molecules.T[(distance_to_nuc < nuc_distance_threshold) & (distance_to_tf == 0)].astype(np.int64)

    # print(distance_to_nuc_array)
    print('num nuc bases: {}'.format(len(observed_data_nuc)))
    # np.savetxt('nuc_dist.{}.txt'.format(rep), distance_to_nuc_array)
    # np.savetxt('nuc_data.{}.txt'.format(rep), observed_data_nuc)

    # feed into LR classifier
    nuc_logistic_classifier = LogisticRegression(random_state=0).fit(distance_to_nuc_array.reshape(-1,1), observed_data_nuc)

    # now do TFs
    # # only select sites that are covered by TFs
    # distance_to_tf_array = distance_to_tf[distance_to_tf == 1]

    # apply same filter to data (don't think we need to do any nuc filtering here since TFs and nucs can't overlap...)
    observed_data_tf = single_molecules.T[distance_to_tf == 1]
    # print(observed_data_tf)
    np.savetxt('tf_data.{}.txt'.format(rep), observed_data_tf)
    print('num tf bases: {}'.format(len(observed_data_tf)))

    # get mean protection
    prob_unmeth_given_tf = observed_data_tf.mean()

    # estimate open
    # get temporary matrix that encodes whether a site is in the promoter
    gc_in_promoter = np.zeros((len(single_molecules.T), len(methyl_positions)))
    promoter_low, promoter_high = promoter_positions
    promoter_indices = [i for i, p in enumerate(methyl_positions) if p >= promoter_low and p <= promoter_high]
    gc_in_promoter[:,promoter_indices] = 1

    # remove all states/GpCs that are bound by a nuc or under a TF (or in promoter)
    observed_data_open = single_molecules.T[(distance_to_tf == 0) & (distance_to_nuc > 80) & (gc_in_promoter == 0)]
    np.savetxt('open_data.{}.txt'.format(rep), observed_data_open)
    print('num open bases: {}'.format(len(observed_data_open)))

    # get mean protection
    obs_prob_t = observed_data_open.mean()
    # now need to convert this to p_unmeth_given_open from the observed unmeth
    p_t_given_meth, p_t_given_unmeth = compute_conversion_hyperparameters()
    prob_unmeth_given_open = (obs_prob_t - p_t_given_meth) / (p_t_given_unmeth - p_t_given_meth)
    prob_unmeth_given_open = 0.05

    # construct dictionary to return
    params_to_return = {'nuc_lr': nuc_logistic_classifier, 'prob_unmeth_given_tf': prob_unmeth_given_tf, 'prob_unmeth_given_open': prob_unmeth_given_open}

    return params_to_return 

def has_em_converged(old_em_hyperparams, new_em_hyperparams, fractional_tolerance=0.1):
    '''
    function to determine whether the EM algorithm has converged
    compares the hyperparams with some tolerance and decides whether or not it needs another iteration
    '''

    # params_to_check = ['nuc_lr_mu', 'nuc_lr_s', 'prob_unmeth_given_tf', 'tf_lr_mu', 'tf_lr_s', 'prob_unmeth_given_open']
    params_to_check = ['prob_unmeth_given_tf', 'prob_unmeth_given_open']
    # params_to_check = ['nuc_lr_mu', 'nuc_lr_s', 'prob_unmeth_given_tf', 'prob_unmeth_given_open']

    print('hyperparams')
    print(old_em_hyperparams['prob_unmeth_given_tf'])
    print(new_em_hyperparams['prob_unmeth_given_tf'])
    print(old_em_hyperparams['prob_unmeth_given_open'])
    print(new_em_hyperparams['prob_unmeth_given_open'])
    print(old_em_hyperparams['nuc_lr'].coef_)
    print(old_em_hyperparams['nuc_lr'].intercept_)
    print(new_em_hyperparams['nuc_lr'].coef_)
    print(new_em_hyperparams['nuc_lr'].intercept_)

    for p in params_to_check:
        if np.abs(old_em_hyperparams[p] - new_em_hyperparams[p]) / old_em_hyperparams[p] > fractional_tolerance:
            # print(np.abs(old_em_hyperparams[p] - new_em_hyperparams[p]) / old_em_hyperparams[p])
            return False

    return True


def expand_state(tf_state, nuc_state, idx, tfbs_positions, bin_size=10, nuc_length=140):
    '''
    takes a single "state" and returns a legible version that we can plunk in a df and return
    note that this is clunky, partially to retrofit current formatting
    will probably rewrite this, later...
    '''
    # nuc_positions, binding_sites = state

    # print(tf_state)

    expanded_state = {}

    # first do nuc
    # nuc_num = 1
    for idx in nuc_state.index:
        expanded_state['nuc{}_present'.format(idx)] = nuc_state[idx] == 1
        bounds = get_nuc_start_and_end_from_midpoint(idx, bin_size, nuc_length)
        expanded_state['nuc{}_start'.format(idx)] = bounds[0]
        expanded_state['nuc{}_end'.format(idx)] = bounds[1]
        # nuc_num += 1 

    # then do tfbs
    # tfbs_num = 1
    # print(tf_state.index)
    for idx in tf_state.index:
        expanded_state['tfbs_{}'.format(idx)] = tf_state[idx] == 1
        # tfbs_num += 1

    # finally tack on index of state 
    expanded_state['idx'] = idx

    return expanded_state

def summarize_assignments(tf_states, nuc_states, assignments, tfbs_positions):
    '''
    gets legible state assignment for every single molecule
    put into df and return
    optinally fill in nas (in nucleosomes) here?
    '''
    expanded_assignments = []

    for a in assignments:
        # expanded_assignments.append(expand_state(state_list[a], a, tfbs_positions))
        expanded_assignments.append(expand_state(tf_states.loc[a], nuc_states.loc[a], a, tfbs_positions))

    return pd.DataFrame(expanded_assignments)

def write_valid_states(states, out_path, tfbs_positions):
    '''
    we are going to write the valid states to a file
    this will be in parallel to the file we are already writing (e.g. the per molecule file)
    we can use this to assign energies to each state and compare the predicted abundances (boltzman) to the true ones
    we'll also annotate the single molecule classifications with the index of each state, so we can aggregate by it later
    '''
    expanded_assignments = []

    for idx, state in enumerate(states):
        expanded_assignments.append(expand_state(state, idx, tfbs_positions))

    to_write = pd.DataFrame(expanded_assignments)

    to_write.to_csv(out_path, sep='\t', header=True, index=False)

def compute_classifications(single_molecules, methyl_positions, tfbs_positions, valid_states_path=None, precomputed_states_file=None, memory_threshold=1., do_em=False):
    # initialize hyperparams
    p_t_given_meth, p_t_given_unmeth = compute_conversion_hyperparameters()

    # check whether there already exists single-molecule states for this amplicon
    if precomputed_states_file and path.exists(precomputed_states_file):
        # then we can save time and just load the precomputed possible methylation states
        with open(precomputed_states_file, "rb") as in_file:
            tf_states, nuc_states = pickle.load(in_file)
    else:
        # generate possible single molecule states
        print('generating states')
        tf_states, nuc_states = enumerate_states(methyl_positions, tfbs_positions)

        print('states: {}'.format(len(tf_states)))

        # # filter for non-overlapping states
        # valid_states = prune_states(possible_states, methyl_positions)

        # write to file
        if precomputed_states_file:
            with open(precomputed_states_file, "wb") as out_file:
                print('started dump')
                pickle.dump((tf_states, nuc_states), out_file)
                print('ended dump')

    print('valid states: {}'.format(len(tf_states)))

    # # delete me later
    # if valid_states_path:
    #     write_valid_states(valid_states, valid_states_path, tfbs_positions)

    # determine memory usage
    data_size = 8 # can get this down
    n_matrix_partitions = determine_chunk_size(single_molecules.shape[1], len(tf_states), data_size, memory_threshold)
    print('n chunks: {}'.format(n_matrix_partitions))

    # optionally start the EM loop
    if do_em:
        print('starting EM algorithm')
        em_loops = 0
        converged = False

        hyperparams = get_initial_em_hyperparams()

        # then loop
        while not converged:
            print('num loops: {}'.format(em_loops))
            # turn states + hyperparams into log_prob matrices for MLE
            log_prob_states_t, log_prob_states_c = create_bernoulli_logprob_matrices(tf_states, nuc_states, methyl_positions, tfbs_positions, hyperparams, p_t_given_meth, p_t_given_unmeth)

            # then do MLE for states
            maximum_likelihood_assignments = classify_all_molecules(single_molecules, log_prob_states_t, log_prob_states_c, n_matrix_partitions)

            np.savetxt('ml_assignments.{}.txt'.format(em_loops), maximum_likelihood_assignments)

            # then get new hyperparams (E step?)
            new_em_hyperparams = estimate_latent_params_from_state_assignments(single_molecules, methyl_positions, tfbs_positions, maximum_likelihood_assignments, tf_states, nuc_states, em_loops)

            # compare new and old hyperparams to see if we should stop iterating
            converged = has_em_converged(hyperparams, new_em_hyperparams)

            # iterate the loop one more time
            hyperparams = new_em_hyperparams
            em_loops += 1
    else:
        print('classification, no EM algorithm this time')

        hyperparams = get_initial_em_hyperparams()

        # turn states + hyperparams into log_prob matrices for MLE
        log_prob_states_t, log_prob_states_c = create_bernoulli_logprob_matrices(tf_states, nuc_states, methyl_positions, tfbs_positions, hyperparams, p_t_given_meth, p_t_given_unmeth)

        # then do MLE for states
        maximum_likelihood_assignments = classify_all_molecules(single_molecules, log_prob_states_t, log_prob_states_c, n_matrix_partitions)


    # # turn these into (log) prob matrices
    # log_prob_states_t, log_prob_states_c = create_bernoulli_logprob_matrices(tf_states, nuc_states, methyl_positions, tfbs_positions, p_t_given_meth, p_t_given_unmeth)

    # # classify single molecules
    # print('doing classification')
    # maximum_likelihood_assignments = classify_all_molecules(single_molecules, log_prob_states_t, log_prob_states_c, n_matrix_partitions)

    # coerce (new) assignments to output
    print('summarizing output')
    assignments_to_return = summarize_assignments(tf_states, nuc_states, maximum_likelihood_assignments, tfbs_positions)

    # optionally write valid_states to a file, so we can start building the equilibrium model to predict all microstates
    if valid_states_path:
        write_valid_states(valid_states, valid_states_path, tfbs_positions)

    return assignments_to_return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Bayesian classification algorithm on an input single-molecule matrix')
    parser.add_argument("--input", dest="input", type=str, help="Path to single molecule matrix")
    parser.add_argument("--output", dest="output", type=str, help="Path for writing the classification output tables to")
    parser.add_argument("--all_states_output", dest="all_states_output", type=str, default=None, help="Path for writing the valid state matrix")
    parser.add_argument("--precomputed_states_file", dest="precomputed_states_file", type=str, default=None, help="Path to a precomputed methylation file to save computation time, if it doesn't exist it will be generated")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file (optional)")
    parser.add_argument("--positions", dest="positions_file", type=str, help="Path with positions for calling bound/unbound")
    parser.add_argument("--amplicon_name", dest="amplicon_name", type=str, help="Name of amplicon aligned to for grabbing TFBS positions")
    parser.add_argument("--reads_to_use", dest="reads_to_use", type=int, default=0, help="Number of reads to sample for classification (if 0, use all)")
    parser.add_argument("--filter_threshold", dest="filter_threshold", type=float, default=1., help="Fraction of converted Cs above which a molecule will be thrown out (to get rid of the totally 'closed' molecules")
    parser.add_argument("--convert_ambiguous_gcgs", dest="convert_ambiguous_gcgs", type=str, default="53,49,129,122", help="Comma-separated list of (pairs of) positions to apply the ambiguous GCG correction to, first entry is the gcg, second entry is the one to replace with (when methylated)")
    parser.add_argument("--reads_to_plot", dest="reads_to_plot", type=int, default=10, help="Number of classifications to plot, to make sure the algorithm is performing as expected")
    parser.add_argument("--memory_threshold", dest="memory_threshold", type=float, default=1., help="Amount of memory to stay under for the big mat mul (in Gb)")
    parser.add_argument("--individual_reads_to_plot", dest="individual_reads_to_plot", type=str, default="", help="Comma-separated list of molecules that the script has to plot for you")
    parser.add_argument('--do_em', dest='do_em', action='store_true', help="Whether or not to perform expectation maximazation on the hyperparameters")
    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        # read files
        print('loading data')
        mat = load_single_molecule_matrix(args.input)
        methyl_positions = np.array(get_methyl_positions(mat)) 

        pos_file = load_tfbs_positions(args.positions_file)
        tfbs_positions = pos_file[args.amplicon_name]

        # optionally filter and downsample the matrix
        # remove reads with -1s still in them??
        mat = filter_all_converted_reads(mat, args.filter_threshold)
        mat = fix_missing_data(mat)

        # optionally apply the 4/3 and 9/8 corrections (e.g. replace methylated 4/9 with the neighbor)
        if args.convert_ambiguous_gcgs:
            positions_to_fix = list(map(int, args.convert_ambiguous_gcgs.split(",")))
            assert len(positions_to_fix) % 2 == 0
            num_fixes = len(positions_to_fix) // 2
            for idx in range(num_fixes):
                mat = adjust_gcgs(mat, positions_to_fix[2*idx], positions_to_fix[2*idx+1])
            # mat = adjust_gcgs(mat, 53, 49) # 4/3
            # mat = adjust_gcgs(mat, 129, 122) # 9/8

        print('num molecules: {}'.format(len(mat)))

        # transpose the matrix (to gpcs x molecules, and also convert to numpy array)
        single_molecules = mat.T.values

        # get assignments
        print('starting classification')
        assignments = compute_classifications(single_molecules, methyl_positions, tfbs_positions, valid_states_path=args.all_states_output, precomputed_states_file=args.precomputed_states_file, memory_threshold=args.memory_threshold, do_em=args.do_em)
        assignments.index = mat.index

        # write this to a file
        assignments.to_csv(args.output, sep='\t', header=True, index=True)

        # optionally plot some
        if args.reads_to_plot:
            print('plotting')
            num_reads_to_plot = min(args.reads_to_plot, len(mat))
            to_plot = assignments.sample(num_reads_to_plot)
            for idx, assignment in to_plot.iterrows():
                fig, ax = plt.subplots(figsize=(12,8))
                plot_single_read(mat, idx, ax, fillbetween=tfbs_positions)
                decorate_single_read_plot2(assignment, ax, tfbs_positions)
                plots.savefig()
                plt.close()

            # optionally plot specific molecules, for reproducibility during testing
            if args.individual_reads_to_plot:
                indiv_reads = map(int, args.individual_reads_to_plot.split(','))
                for idx in indiv_reads:
                    if idx in mat.index.tolist():
                        fig, ax = plt.subplots(figsize=(12,8))
                        plot_single_read(mat, idx, ax, fillbetween=tfbs_positions)
                        decorate_single_read_plot2(assignments.loc[idx], ax, tfbs_positions)
                        plots.savefig()
                        plt.close()




            