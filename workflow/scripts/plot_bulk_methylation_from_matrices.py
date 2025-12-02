import pandas as pd
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
from Bio import SeqIO
import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns

from common import plot_bulk_smf_trace_from_matrix, load_single_molecule_matrix


def plot_bulk_methylation(input_prefix, amplicon_fa, plots):
    # load amplicon fasta
    amplicon_to_seq_dict = {}
    
    for r in list(SeqIO.parse(amplicon_fa, "fasta")):
        # print(r)
        amplicon_to_seq_dict[r.id] = r.seq

    amplicons = amplicon_to_seq_dict.keys()

    # iterate through amplicons, load single molecule matrix, filter on columns, grab means, plot
    for amplicon in amplicons:
        matrix_path = '{}.{}.full_unclustered.matrix'.format(input_prefix, amplicon)

        # check whether this matrix was actually written
        if path.exists(matrix_path):
            # load matrix and convert -1s -> NAs for .mean()ing on columns
            mat = load_single_molecule_matrix(matrix_path, every_other=False)
            print(len(mat.loc[:,mat.mean() > 0].mean()))
            mat = mat.replace(-1,np.nan)

            average_methylation = 100 * mat.mean()

            # plot bulk methyl signal
            plot_bulk_smf_trace_from_matrix(average_methylation, plots, title=amplicon, ylim=(0,103))

            # # plot the average signal from GpC and non-GpC
            # all_cpg_df['context'] = 'CpG'
            # all_gpc_df['context'] = 'GpC'
            # all_other_c_df['context'] = 'C'
            # # all_gpc_df['GpC'] = True
            # # all_other_c_df['GpC'] = False

            # fig, ax = plt.subplots()
            # sns.barplot(data=all_other_c_df.append([all_cpg_df,all_gpc_df]), x='context', y='pct')
            # # sns.barplot(data=all_gpc_df.append(all_other_c_df), x='GpC', y='pct')

            # ax.set_ylabel('%Methylation')

            # plt.tight_layout()
            # plots.savefig()
            # plt.close()

            # # return merged table
            # columns = ['chr', 'start', 'SMF', 'total_reads']
            # # output_tables.append(all_gpc_df[columns])
            # output_tables.append(cs_to_plot[columns])
        else:
            print('No data for {}'.format(amplicon))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Path to the directory with single-molecule matrices")
    parser.add_argument("--amplicon", dest="amplicon_fa", type=str, help="FASTA file containing the amplicons that the reads were aligned to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    
    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        plot_bulk_methylation(args.input, args.amplicon_fa, plots)


# # calculate GpC positions (also CpG)
#     gpc_pos_dict = {}
#     cpg_pos_dict = {}
#     for amplicon in amplicons:
#         gpc_pos_dict[amplicon] = [i for i, b in enumerate(amplicon_to_seq_dict[amplicon]) if b=='C' and amplicon_to_seq_dict[amplicon][max(i-1, 0)]=='G']
#         cpg_pos_dict[amplicon] = [i for i, b in enumerate(amplicon_to_seq_dict[amplicon]) if b=='C' and amplicon_to_seq_dict[amplicon][min(i+1, len(amplicon_to_seq_dict[amplicon])-1)]=='G']

#     output_tables = []

# gpcs = gpc_pos_dict[amplicon]
#             cpgs = cpg_pos_dict[amplicon]

            
#             # subset by C type
#             all_gpc_df = all_c_df.loc[all_c_df.start.isin(gpcs)].copy()
#             all_cpg_df = all_c_df.loc[all_c_df.start.isin(cpgs)].copy()
#             all_other_c_df = all_c_df.loc[~all_c_df.start.isin(gpcs+cpgs)].copy()

#             # also make a plot with only GpC and only CpG if also using CpG MTase, to see if everything worked
#             if include_cpg:
#                 plot_bulk_smf_trace(all_gpc_df, plots, title='{}: GpC only'.format(amplicon), ylim=(0,103))
#                 plot_bulk_smf_trace(all_cpg_df, plots, title='{}: CpG only'.format(amplicon), ylim=(0,103))
