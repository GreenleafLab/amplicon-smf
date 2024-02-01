import pandas as pd
import numpy as np
from pymc3 import find_MAP, Uniform, Beta, Binomial, sample, Model, NUTS, NegativeBinomial, Normal, HalfStudentT, DiscreteUniform, Mixture, ZeroInflatedPoisson, Bernoulli , Data, Deterministic, set_data, model_to_graphviz, Potential
from pymc3.math import switch, or_, and_
# import arviz as az
import argparse
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')

# turn off the printing from pymc3!
import logging
logger = logging.getLogger('pymc3')
logger.setLevel(logging.ERROR)

import warnings

from common import load_single_molecule_matrix, get_methyl_positions, load_tfbs_positions, filter_all_converted_reads, decorate_single_read_plot, coerce_classification_dict_to_state_assignment, plot_single_read

# ## TODO:
# Add support for removing reads that don't fit well (by examining posterior or posterior checks)
#     e.g. look at Gelman-Rubin stats or autocorr or something!


def convert_posterior_to_assignment(trace, ntfbs, nuc_thresh=0.5, tfbs_thresh=0.4):
    classification = {}

    nuc1_present = np.array(trace.posterior.nuc1_present).mean() > nuc_thresh
    classification['nuc1_present'] = nuc1_present
    nuc1_start = int(np.array(trace.posterior.nuc1_start).mean())
    nuc1_end = nuc1_start + int(np.array(trace.posterior.nuc1_length).mean())
    classification['nuc1_start'] = nuc1_start
    classification['nuc1_end'] = nuc1_end

    nuc2_present = np.array(trace.posterior.nuc2_present).mean() > nuc_thresh
    classification['nuc2_present'] = nuc2_present
    nuc2_start = int(np.array(trace.posterior.nuc2_start).mean())
    nuc2_end = nuc2_start + int(np.array(trace.posterior.nuc2_length).mean())
    classification['nuc2_start'] = nuc2_start
    classification['nuc2_end'] = nuc2_end

    for b in range(1,ntfbs+1):
        tfbs_state = np.array(getattr(trace.posterior, 'b{}'.format(b))).mean() > tfbs_thresh
        classification['tfbs_{}'.format(b)] = tfbs_state

    return classification

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Bayesian classification algorithm on an input single-molecule matrix')
    parser.add_argument("--input", dest="input", type=str, help="Path to single molecule matrix")
    parser.add_argument("--output", dest="output", type=str, help="Path for writing the classification output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file (optional)")
    parser.add_argument("--positions", dest="positions_file", type=str, help="Path with positions for calling bound/unbound")
    parser.add_argument("--amplicon_name", dest="amplicon_name", type=str, help="Name of amplicon aligned to for grabbing TFBS positions")
    parser.add_argument("--reads_to_use", dest="reads_to_use", type=int, default=0, help="Number of reads to sample for classification (if 0, use all)")
    parser.add_argument("--filter_threshold", dest="filter_threshold", type=float, default=1., help="Fraction of converted Cs above which a molecule will be thrown out (to get rid of the totally 'closed' molecules")
    parser.add_argument("--reads_to_plot", dest="reads_to_plot", type=int, default=10, help="Number of classifications to plot, to make sure the algorithm is performing as expected")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        # read files
        mat = load_single_molecule_matrix(args.input)
        methyl_positions = get_methyl_positions(mat) 

        pos_file = load_tfbs_positions(args.positions_file)
        tfbs_positions = pos_file[args.amplicon_name]

        n_tfbs = len(tfbs_positions)

        # optionally filter and downsample the matrix
        # remove reads with -1s still in them??
        mat = filter_all_converted_reads(mat, args.filter_threshold)

        if args.reads_to_use:
            # check whether reads_to_use is > number of molecules, if so replace
            reads_to_sample = min(args.reads_to_use, len(mat)-1)
            mat = mat.sample(n=reads_to_sample)

        # build dataset (kind of weird, need to do for pymc3 I guess)
        observed_data = [mat.loc[idx] for idx in mat.index.tolist()]

        # initialize model
        with Model() as smf:
            # Hyperparameters
            p_unmeth_given_bound = 0.99 # 0.9
            p_unmeth_given_unbound = 0.01 # 0.1
            # t is the thing we measure as a 1 in the data
            p_t_given_meth = 0.05
            p_t_given_unmeth = 0.99
            
            # marginalize
            p_t_given_bound = p_unmeth_given_bound * p_t_given_unmeth + (1 - p_unmeth_given_bound) * p_t_given_meth
            p_t_given_unbound = p_unmeth_given_unbound * p_t_given_unmeth + (1 - p_unmeth_given_unbound) * p_t_given_meth
                
            # nucleosome
            nuc1_present = Bernoulli('nuc1_present', p=0.5)
            nuc2_present = Bernoulli('nuc2_present', p=0.5)
            
            nuc1_start = DiscreteUniform('nuc1_start', lower=100, upper=390)
            nuc1_length = Normal('nuc1_length', mu=147, sigma=7)
            nuc1_end = Deterministic('nuc1_end', nuc1_start+nuc1_length)
            
            nuc2_start = DiscreteUniform('nuc2_start', lower=250, upper=490)
            nuc2_length = Normal('nuc2_length', mu=147, sigma=7)
            nuc2_end = Deterministic('nuc2_end', nuc2_start+nuc2_length)
            
            non_overlapping_potential = Potential('non_overlapping_potential', switch(and_(nuc1_end-10>nuc2_start+10,
                                                                                           and_(nuc1_present>0.5,nuc2_present>0.5)), -np.inf, 0))
            
            order_potential = Potential('order_potential', switch(and_(nuc1_start>nuc2_start,
                                                                       and_(nuc1_present>0.5,nuc2_present>0.5)), -np.inf, 0))
            
            # tfbs
            p_bound = 0.6 # can also modify this to have site-specific prior, EB on the dataset?
            
            # do conditional on nucleosome position to get p_tf_bound_i
            binding_probs = [switch(or_(and_(nuc1_present>0.5,
                                              and_(nuc1_start<=tfbs_positions[idx][1],
                                                   nuc1_end>=[idx][0])),
                                         and_(nuc2_present>0.5,
                                              and_(nuc2_start<=tfbs_positions[idx][1],
                                                   nuc2_end>=[idx][0]))),
                                    0.01, # CHANGE ME!!!
                                    p_bound) 
                             for idx in range(n_tfbs)]

            # make a bernoulli RV from p_bound
            binding_bernoullis = [Bernoulli('b{}'.format(idx+1), p=binding_probs[idx]) for idx in range(n_tfbs)]
            
            # make masks
            masks = [[1 if (p>=tfbs_positions[idx][0] and p<=tfbs_positions[idx][1]) else 0 for p in methyl_positions] for idx in range(n_tfbs)]
            series_masks = [pd.Series(masks[idx], index=methyl_positions) for idx in range(n_tfbs)]

            # convert the states into occupancies at each site
            tf_bound = 0
            for idx in range(n_tfbs):
                tf_bound += binding_bernoullis[idx] * series_masks[idx]

            nuc_bound = switch(or_(and_(nuc1_present>0.5,and_(methyl_positions>=nuc1_start,methyl_positions<=nuc1_end)),and_(nuc2_present>0.5,and_(methyl_positions>=nuc2_start,methyl_positions<=nuc2_end))),1,0)
            bound_bernoulli_p = switch(or_(tf_bound,nuc_bound),p_t_given_bound,p_t_given_unbound)

            # data
            data = Data("data", observed_data[0])
            
            # likelihood
            observed_seq = Bernoulli('observed_seq', p=bound_bernoulli_p, observed=data)

        # generate one trace for each dataset
        # wrap this all in a block to avoid pringint a bunch of runtime warnings to the console..
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            traces = []
            for idx,data_vals in enumerate(observed_data):
                if idx%10 == 0:
                    print('{} molecules processed'.format(idx))
                # need to wrap in try/catch, sometimes certain molecules are badly formatted? idk but need to
                try:
                    with smf:
                        set_data({"data": data_vals})
                        traces.append(sample(1000, return_inferencedata=True, progressbar=False))
                except Exception as e: 
                    print('Molecule {} failed'.format(mat.index.tolist()[idx]))
                    print(e)
                    traces.append(None)



        # inspect traces, compute MAP, classify, and add to list
        # note that we are maybe not necessarily getting the MAP but we are assuming independence and then marginalizing across posterior
        # can do more fancy fitting with the full posterior later?
        classifications = []
        for idx,trace in enumerate(traces):
            if trace: # have to take into account empty traces e.g. errors
                assignment = convert_posterior_to_assignment(trace, n_tfbs)
                index = mat.index.tolist()[idx]
                assignment['index'] = index
                classifications.append(assignment)

        # write to file
        pd.DataFrame(classifications).to_csv(args.output, sep='\t', header=True, index=False)

        # plot
        if args.reads_to_plot:
            for idx in range(args.reads_to_plot):
                read_idx = classifications[idx]['index']
                fig, ax = plt.subplots(figsize=(12,8))
                plot_single_read(mat, read_idx, ax, fillbetween=tfbs_positions)
                state_classification = coerce_classification_dict_to_state_assignment(classifications[idx], n_tfbs)
                decorate_single_read_plot(state_classification, ax, tfbs_positions)
                plots.savefig()
                plt.clf()




