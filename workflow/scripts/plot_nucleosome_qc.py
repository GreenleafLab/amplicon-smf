import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.smoothers_lowess import lowess
import argparse
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
from copy import copy
from sklearn.mixture import GaussianMixture
from Bio import SeqIO
from os import path

from common import plot_bulk_smf_trace_from_matrix, load_single_molecule_matrix

def get_all_protection_streaks(methyls, positions):
    # append a 0 to methyls so we always finish with a nonmethyl (so we don't need edge cases)
    methyls = methyls + [0]
    # also need to do likewise to positions
    positions = positions + [np.inf]
    streak = False
    start = -1
    end = -1
    length = 0
    
    streaks = []

    for idx in range(len(methyls)):
        if streak:
            if methyls[idx] == 1:
                end = positions[idx]
                length = end - start
            else:
                streak = False
                streaks.append((length,start,end))
                start = positions[idx]
                end = positions[idx]
                length = 0
        else:
            if methyls[idx] == 1:
                start = positions[idx]
                end = positions[idx]
                length = 0
                streak = True
            else:
                pass

    return sorted(streaks,key=lambda x: x[0],reverse=True)

def fit_gmm_to_hist(data):
    gm = GaussianMixture(n_components=2, random_state=0).fit(np.array(data).reshape(-1, 1))
    return gm
    

def plot_nuc_len_histogram(lens, plots, title='', gmm=None):
    '''
    
    '''
    fig, ax = plt.subplots()

    ax.hist(lens, bins='auto', density=True)

    if gmm:
        # gmm = fit_gmm_to_hist(lens)
        xs = np.linspace(0,600,600)
        pdf = np.exp(gmm.score_samples(xs.reshape(-1, 1)))
        responsibilities = gmm.predict_proba(xs.reshape(-1, 1))
        pdf_individual = responsibilities * pdf[:, np.newaxis]

        # ax.plot(xs, pdf, '-k')
        ax.plot(xs, pdf_individual, '--k')

        ax.text(0.7*ax.get_xlim()[-1], 0.9*ax.get_ylim()[-1], 'Lower mean: {:.0f}\nUpper mean: {:.0f}'.format(sorted(gmm.means_.flatten())[0],sorted(gmm.means_.flatten())[1]))

    ax.set_xlabel('Longest protection length')
    ax.set_ylabel('Molecules')
    ax.set_title(title)

    plt.tight_layout()
    plots.savefig()
    plt.close()

def plot_nucleosome_qc(input_prefix, amplicon_fa, plots, results, fit_gmm):
    # load amplicon fasta
    amplicon_to_seq_dict = {}
    
    for r in list(SeqIO.parse(amplicon_fa, "fasta")):
        amplicon_to_seq_dict[r.id] = r.seq

    amplicons = amplicon_to_seq_dict.keys()

    amplicon_to_nuc_len_dict = {}

    # iterate through amplicons, load single molecule matrix, filter on columns, grab means, plot
    for amplicon in amplicons:
        matrix_path = '{}.{}.full_unclustered.matrix'.format(input_prefix, amplicon)

        # check whether this matrix was actually written
        if path.exists(matrix_path):
            # load matrix
            mat = load_single_molecule_matrix(matrix_path)
            positions = mat.columns.tolist()
            lengths = []

            for idx in mat.index:
                methyls = mat.loc[idx].tolist()
                protection_streaks = get_all_protection_streaks(methyls, positions)
                if protection_streaks:
                    longest_streak = protection_streaks[0][0]
                    lengths.append(longest_streak)

            amplicon_to_nuc_len_dict[amplicon] = copy(lengths)

    # consolidate all the results into one list for consolidated plotting
    # [item for sublist in l for item in sublist]
    all_lengths = [l for lens in amplicon_to_nuc_len_dict.values() for l in lens]
    amplicon_to_nuc_len_dict['All_Amplicons'] = all_lengths

    # plot
    # do all
    data = amplicon_to_nuc_len_dict['All_Amplicons']
    gmm = fit_gmm_to_hist(data)
    plot_nuc_len_histogram(data, plots, title='All_Amplicons', gmm=gmm if fit_gmm else None)

    # write results
    results.write('{}\t{:.0f}\n'.format('nuc_lower_mean_all', sorted(gmm.means_.flatten())[0]))
    results.write('{}\t{:.0f}\n'.format('nuc_upper_mean_all', sorted(gmm.means_.flatten())[1]))
    results.write('{}\t{:.2f}\n'.format('nuc_frac_in_lower_mode_all', np.mean(gmm.predict(np.array(data).reshape(-1, 1))==(0 if gmm.means_.flatten()[0] < gmm.means_.flatten()[1] else 1))))

    # do no tfbs
    no_tfbs = ['opJS4_0x_TetO_21bp_no_CG', 'opJS5_0xTetO_18bp_b1', 'opJS5_0xTetO_18bp_b2', 'BD24', 'no_TFBS', 'background']
    lower_mean, upper_mean, frac, reads = 0, 0, 0, 0

    for amplicon in no_tfbs:
        if amplicon in amplicon_to_nuc_len_dict:
            data = amplicon_to_nuc_len_dict[amplicon]
            gmm = fit_gmm_to_hist(data)
            plot_nuc_len_histogram(data, plots, title=amplicon, gmm=gmm if fit_gmm else None)
            if len(data) > reads:
                reads = len(data)
                lower_mean = sorted(gmm.means_.flatten())[0]
                upper_mean = sorted(gmm.means_.flatten())[1]
                frac = np.mean(gmm.predict(np.array(data).reshape(-1, 1))==(0 if gmm.means_.flatten()[0] < gmm.means_.flatten()[1] else 1))

    if reads > 0:
        results.write('{}\t{:.0f}\n'.format('nuc_lower_mean_best_0x', lower_mean))
        results.write('{}\t{:.0f}\n'.format('nuc_upper_mean_best_0x', upper_mean))
        results.write('{}\t{:.2f}\n'.format('nuc_frac_in_lower_mode_best_0x', frac))

    for amplicon in amplicons:
        if amplicon in no_tfbs:
            continue
        else:
            if amplicon in amplicon_to_nuc_len_dict:
                plot_nuc_len_histogram(amplicon_to_nuc_len_dict[amplicon], plots, title=amplicon, gmm=None)
            else:
                continue

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Path to the directory with single-molecule matrices")
    parser.add_argument("--amplicon", dest="amplicon_fa", type=str, help="FASTA file containing the amplicons that the reads were aligned to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--results", dest="results_file", type=str, default='', help="Path to output stats file")
    parser.add_argument("--gmm", dest="gmm", action='store_true', help="Whether to fit GMM on the historgram")
    
    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots, open(args.results_file, 'w') as results:
        plot_nucleosome_qc(args.input, args.amplicon_fa, plots, results, args.gmm)


