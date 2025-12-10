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

# to do list:
# - decrease the "overhang" parameter, so that for e.g. 6x a nucleosome can extend a bit more off the edge and permit the last guy to be bound
# - - was 40, now making it 30? might make more sense to just require idk at least 3 gpcs at the end to call it but this seems harder to code rn
# - maybe require middle guys to be accessible? but i guess we still can't explain the hemi-methylated guys...
# - there are still some cases (see screenshot on desktop) where i guess the exact positioning of nucleosomes means it's better to mess up two gpcs next to the flanks and get in another nuc, rather than have one fewer nuc and get binding right, have to fix this (upweighting getting the bases around tetos right?)
# - - i think both last two points can potentially be solved by increasing the weighting in the log likelihood of the middle bases, or the guys immediately next to the flanking guys, although sometimes there is still weird behavior. idk it's mostly right but we can talk about this
# - i think we want to decrease the nucleosome length, to what we measure in our data, i don't think a 150bp long streak (from gpc start to gpc end) is reasonable
# - - was (130,150), now trying (100,140), which I think is basically FWHM of the longest protection length distribution 


def compute_hyperparameters():
    '''
    currently estimating this from what I expect things to be
    but in the future can grab from paired EM-seq data? not sure it matters...
    '''
    p_unmeth_given_bound = 0.99 # 0.9
    p_unmeth_given_unbound = 0.01 # 0.1
    # t is the thing we measure as a 1 in the data (i.e. unmethylated = bound = protected whatever)
    p_t_given_meth = 0.05
    p_t_given_unmeth = 0.99
    
    # marginalize
    p_t_given_bound = p_unmeth_given_bound * p_t_given_unmeth + (1 - p_unmeth_given_bound) * p_t_given_meth
    p_t_given_unbound = p_unmeth_given_unbound * p_t_given_unmeth + (1 - p_unmeth_given_unbound) * p_t_given_meth

    return (p_t_given_bound, p_t_given_unbound)

def find_first_matching_sublist(lst, span, start):
    '''
    helper function for below
    basically, given a start position, figure out the first subsequent position that satisfies the run being in the range
    returns the index of this position, or else returns -1 if it falls off the end
    '''
    min_span, max_span = span
    
    end = start
    
    while end < len(lst):
        temp_span = lst[end] - lst[start] + 1 # fencepost
        if temp_span >= min_span and temp_span <= max_span:
            return end
        else:
            end += 1

    return -1

def sublist_finder(lst, span, start_ptr):
    '''
    generic function to compute all possible valid runs of a given length from a list of integers
    input is a lst of positive integers (no duplicates, sorted)
    also a range (low and high for how long the run can be, in the units of the entries of the list)
    also a start (since this is a recursive function)
    each valid configuration is a list of tuples of (start,end) of valid runs (so a single output could be [] or [(1,150)] or [(1,150),(161,300)])
    output is a list of these lists
    '''
    valid_substrings = []

    start = start_ptr
    
    while start < len(lst):
        # from start, find a valid range that works
        end = find_first_matching_sublist(lst, span, start)
        if end != -1:
            # save this tuple
            temp_sublist = [(lst[start],lst[end])]
            # recursive call, starting from one past the end of the current
            # note, can make this end+2 to ensure that the model inserts one accessible base between nucs for spacer, not sure if necessary
            recursive_valid_sublists = sublist_finder(lst, span, end+1) 
            # append the empty run (so that we can get the "just this on run" state)
            # probably not the best way to do it but whatever
            recursive_valid_sublists += [[]]
            # add the current run to every recursive one 
            for sublist in recursive_valid_sublists:
                valid_substrings.append(temp_sublist + sublist)
        # now start at the next element in the sequnce
        start += 1

    # only once, append the empty run
    if start_ptr == 0:
        valid_substrings += [[]]
        
    return valid_substrings

def powerset(iterable):
    '''
    computes the powerset of an iterable
    stolen from SO
    '''
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def augment_methyl_positions(methyl_positions, step_size=10, overhang=30, nuc_length=130):
    '''
    adds extra gpc positions to the end of the list (i.e. not observed ones, imaginary)
    this is to allow nucs to slide off the end and still be valid
    can modify how many and by how much a nucleosome has to overlap the molecule
    220911 modified to allow for nucs to slide off "front" end too, which allows us to do +1 nuc
    '''
    start = methyl_positions[0]
    end = methyl_positions[-1]
    return list(range(start-nuc_length+overhang,start-step_size,step_size)) + methyl_positions + list(range(end+step_size, end+nuc_length-overhang, step_size))

def enumerate_states(methyl_positions, tfbs_positions, nuc_span=(110,140)):
    '''
    generates the list of all possible underlying occupancy states
    first adds in a few "ghost" gpcs at the end of the list to allow nucs to slide off the edge
    then computes all possible nucleosome occupancy events
    then compute all possible tf binding events
    returns the outer product of the two lists
    '''
    augmented_methyl_positions = augment_methyl_positions(methyl_positions)
    possible_nuc_positions = sublist_finder(augmented_methyl_positions, nuc_span, 0)
    possible_binding_events = powerset(tfbs_positions)

    return product(possible_nuc_positions, possible_binding_events)

def nuc_overlaps_tf_binding(methyl_positions, state):
    '''
    check whether TF binding overlaps nucleosome binding
    almost the same as generate_predicted_protection (below)
    first fills in all gpcs that are covered by nucs
    then figures out if any TFs are bound there too
    also checks promoter state, too
    '''
    nuc_positions, binding_sites = state
    protection = [0 for _ in range(len(methyl_positions))]
    
    # do nucs
    for idx, p in enumerate(methyl_positions):
        for nuc in nuc_positions:
            if p >= nuc[0] and p <= nuc[1]:
                protection[idx] = 1
            
    # do tfbs
    for idx, p in enumerate(methyl_positions):
        for tfbs in binding_sites:
            if p >= tfbs[0] and p <= tfbs[1]:
                if protection[idx] == 1: # checking to see if nucleosome overlaps a bound TF
                    return True
                if (idx > 0 and protection[idx-1] == 1) or (idx+1 < len(methyl_positions) and protection[idx+1] == 1): # also checking if a nucleosome overlaps the _flanks_ of a bound tf, this way a nucleosome cannot abut a bound tf (ambiguous case, plus impossible?)
                    return True
    
    return False

def prune_states(all_states, methyl_positions):
    '''
    remove (most of) the states where there is TF and nucleosome bound at the same position
    '''
    return [s for s in all_states if not nuc_overlaps_tf_binding(methyl_positions, s)]

def generate_predicted_protection(methyl_positions, state):
    '''
    iterate through each gpc position, if it's bound (by either a nucleosome or a TF) mark the base protected
    '''
    nuc_positions, binding_sites = state
    protection = [0 for _ in range(len(methyl_positions))]
    
    # do nucs
    for idx, p in enumerate(methyl_positions):
        for nuc in nuc_positions:
            if p >= nuc[0] and p <= nuc[1]:
                protection[idx] = 1
            
    # do tfbs
    for idx, p in enumerate(methyl_positions):
        for tfbs in binding_sites:
            if p >= tfbs[0] and p <= tfbs[1]:
                protection[idx] = 1

    return protection
      
def create_bernoulli_logprob_matrices(predicted_methylations, methyl_positions, p_t_given_bound, p_t_given_unbound, promoter_positions=(75,175)):
    '''
    for each state, convert the predicted methylation profile into two matrices to use for computing bernoulli likelihood
    have to muck with the bases in the promoter (as in, if the promoter is accessible, we shouldn't care if there is protection, since we aren't classifying TATA and stuff yet in this model)
    '''
    n_gpcs = len(methyl_positions)
    n_states = len(predicted_methylations)

    # initialize empty matrices
    prob_states_unmeth = np.zeros((n_states,n_gpcs))
    prob_states_meth = np.zeros((n_states,n_gpcs))

    promoter_low, promoter_high = promoter_positions

    for idx in range(n_states):
        p_m_array = np.array(predicted_methylations[idx])

        # when we see "T" = 1 = protected = unmethylated
        prob_states_unmeth[idx] = p_t_given_bound * p_m_array + p_t_given_unbound * (1 - p_m_array)

        # when we see "C" = 0 = accessible = methylated!
        prob_states_meth[idx] = (1 - p_t_given_bound) * p_m_array + (1 - p_t_given_unbound) * (1 - p_m_array)
        
        # do a little promoter hacking
        # i.e. when the promoter dudes are 0s aka accessible, we should allow for seeing both meth or unmeth!
        # commenting this out since now we are including promoter states
        for i, p in enumerate(methyl_positions):
            if p >= promoter_low and p <= promoter_high and p_m_array[i] == 0:
                prob_states_unmeth[idx,i] = (1 - p_t_given_unbound)
        
    log_prob_states_meth = np.log(prob_states_meth)
    log_prob_states_unmeth = np.log(prob_states_unmeth)

    return (log_prob_states_meth, log_prob_states_unmeth)

def determine_chunk_size(n_mols, n_states, data_size, memory_threshold):
    '''
    Figures out how many times we need to split the data into (along the single molecule axis) to not go above memory limits
    '''

    return int((n_mols * n_states * data_size) // (memory_threshold * 1e9)) + 1
        

def classify_all_molecules(single_molecules, log_prob_states_unmeth, log_prob_states_meth, n_matrix_partitions):
    '''
    Compute Bernoulli likelihood all at once 
    note that log_prob_states_unmeth and meth and not strictly reciprocals of each other, since we need to change unoccupied promoter probs (see above)
    single_molecules is gpcs x molecules
    log_prob_states_(un)meth are states x gpcs
    result of matmul is states x molecules, argmax down columns
    return a list of index of best state per molecule
    BD 221213: rewriting to split this into chunks, for a smaller number of single molecules, so as not to exceed memory limits
    '''
    n_mols = single_molecules.shape[1]
    chunk_size = n_mols // n_matrix_partitions + 1

    to_return = []

    for idx in range(n_matrix_partitions):
        sub_mat = single_molecules[:,idx*chunk_size:(idx+1)*chunk_size]
        res = np.argmax(np.exp(log_prob_states_unmeth @ sub_mat + log_prob_states_meth @ (1 - sub_mat)), axis=0)
        to_return += list(res)

    return to_return
    # Old: return np.argmax(np.exp(log_prob_states_unmeth @ single_molecules + log_prob_states_meth @ (1 - single_molecules)), axis=0)

def expand_state(state, idx, tfbs_positions):
    '''
    takes a single "state" and returns a legible version that we can plunk in a df and return
    note that this is clunky, partially to retrofit current formatting
    will probably rewrite this, later...
    '''
    nuc_positions, binding_sites = state

    expanded_state = {}

    # first do nuc
    nuc_num = 1
    for nuc in nuc_positions:
        expanded_state['nuc{}_present'.format(nuc_num)] = True
        expanded_state['nuc{}_start'.format(nuc_num)] = nuc[0]
        expanded_state['nuc{}_end'.format(nuc_num)] = nuc[1]
        nuc_num += 1 

    # then do tfbs
    tfbs_num = 1
    for tfbs in tfbs_positions:
        expanded_state['tfbs_{}'.format(tfbs_num)] = tfbs in binding_sites
        tfbs_num += 1

    # finally tack on index of state 
    expanded_state['idx'] = idx

    return expanded_state

def summarize_assignments(state_list, assignments, tfbs_positions):
    '''
    gets legible state assignment for every single molecule
    put into df and return
    optinally fill in nas (in nucleosomes) here?
    '''
    expanded_assignments = []

    for a in assignments:
        expanded_assignments.append(expand_state(state_list[a], a, tfbs_positions))

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

def compute_classifications(single_molecules, methyl_positions, tfbs_positions, valid_states_path=None, precomputed_states_file=None, memory_threshold=1.):
    # initialize hyperparams
    p_t_given_bound, p_t_given_unbound = compute_hyperparameters()

    # check whether there already exists single-molecule states for this amplicon
    if precomputed_states_file and path.exists(precomputed_states_file):
        # then we can save time and just load the precomputed possible methylation states
        with open(precomputed_states_file, "rb") as in_file:
            valid_states = pickle.load(in_file)
    else:
        # generate possible single molecule states
        print('generating states')
        possible_states = list(enumerate_states(methyl_positions, tfbs_positions))

        print('total states: {}'.format(len(possible_states)))

        # filter for non-overlapping states
        valid_states = prune_states(possible_states, methyl_positions)


        # write to file
        if precomputed_states_file:
            with open(precomputed_states_file, "wb") as out_file:
                print('started dump')
                pickle.dump(valid_states, out_file)
                print('ended dump')

    print('valid states: {}'.format(len(valid_states)))

    # generate predicted methylation profiles per underlying state
    predicted_methylations = [generate_predicted_protection(methyl_positions, s) for s in valid_states]

    # turn these into matrices
    log_prob_states_meth, log_prob_states_unmeth = create_bernoulli_logprob_matrices(predicted_methylations, methyl_positions, p_t_given_bound, p_t_given_unbound)

    # determine memory usage
    data_size = 8 # can get this down
    n_matrix_partitions = determine_chunk_size(single_molecules.shape[1], len(valid_states), data_size, memory_threshold)
    print('n chunks: {}'.format(n_matrix_partitions))

    # classify single molecules
    print('doing classification')
    maximum_likelihood_assignments = classify_all_molecules(single_molecules, log_prob_states_unmeth, log_prob_states_meth, n_matrix_partitions)

    # coerce (new) assignments to output
    print('summarizing output')
    assignments_to_return = summarize_assignments(valid_states, maximum_likelihood_assignments, tfbs_positions)

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
    parser.add_argument("--reads_to_plot", dest="reads_to_plot", type=int, default=10, help="Number of classifications to plot, to make sure the algorithm is performing as expected")
    parser.add_argument("--memory_threshold", dest="memory_threshold", type=float, default=1., help="Amount of memory to stay under for the big mat mul (in Gb)")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        # read files
        print('loading data')
        mat = load_single_molecule_matrix(args.input)
        methyl_positions = get_methyl_positions(mat) 

        pos_file = load_tfbs_positions(args.positions_file)
        tfbs_positions = pos_file[args.amplicon_name]

        # optionally filter and downsample the matrix
        # remove reads with -1s still in them??
        mat = filter_all_converted_reads(mat, args.filter_threshold)
        mat = fix_missing_data(mat)

        # optionally apply the 4/3 and 9/8 corrections (e.g. replace methylated 4/9 with the neighbor)
        # mat = adjust_gcgs(mat, 53, 49) # 4/3
        # mat = adjust_gcgs(mat, 129, 122) # 9/8

        print('num molecules: {}'.format(len(mat)))

        # transpose the matrix (to gpcs x molecules, and also convert to numpy array)
        single_molecules = mat.T.values

        # get assignments
        print('starting classification')
        assignments = compute_classifications(single_molecules, methyl_positions, tfbs_positions, valid_states_path=args.all_states_output, precomputed_states_file=args.precomputed_states_file, memory_threshold=args.memory_threshold)
        assignments.index = mat.index

        # write this to a file
        assignments.to_csv(args.output, sep='\t', header=True, index=True)

        # optionally plot some
        if args.reads_to_plot:
            print('plotting')
            to_plot = assignments.sample(args.reads_to_plot)
            for idx, assignment in to_plot.iterrows():
                fig, ax = plt.subplots(figsize=(12,8))
                plot_single_read(mat, idx, ax, fillbetween=tfbs_positions)
                decorate_single_read_plot2(assignment, ax, tfbs_positions)
                plots.savefig()
                plt.close()

            
