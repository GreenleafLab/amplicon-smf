# input two matrices, plot clustered matrices separately
# then do KMeans on them together (and optionally separately) and return the cluster proportions for each
# change to barplot not pie chart
# then plot shapes of the clusters

import pandas as pd
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns
from sklearn.feature_selection import mutual_info_classif

from common import fix_missing_data, load_tfbs_positions, get_methyl_positions


def plot_mutual_info(mat, methyl_positions, output_plot, name, positions_file, amplicon):
    # do some general formatting
    title = amplicon if amplicon else ''
    name = name if name else ''

    # build matrix
    mi_mat = np.zeros((len(methyl_positions),len(methyl_positions)))

    for i in range(len(methyl_positions)):
        for j in range(len(methyl_positions)):
            mi_mat[i,j] = mutual_info_classif(mat[[methyl_positions[i]]], mat[methyl_positions[j]], discrete_features=True)

    # get motif positions
    # pos_dict = None
    # if positions_file:
    #     pos_dict = load_tfbs_positions(positions_file)
    # if pos_dict and amplicon:
    #     # have to write something here to make sure that the amplicon name georgi style and my style are the same...
    #     positions = pos_dict[amplicon]
    #     add_fillbetween(ax, positions)

    # format labels
    # labels_series = pd.Series(by_teto)
    # lut = dict(zip(labels_series.unique(), sns.hls_palette(len(set(labels_series)))))
    # row_colors = labels_series.map(lut).to_numpy()
    row_colors = None

    sns.clustermap(mi_mat, row_cluster=False, col_cluster=False, cmap='vlag', center=0, row_colors=row_colors, col_colors=row_colors)
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.clf()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Path to single molecule matrix")
    # parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--input_name", dest="input_name", type=str, default=None, help="Name of input sample (e.g. CTCF)")
    parser.add_argument("--amplicon_name", dest="amplicon_name", type=str, default=None, help="Name of amplicon aligned to (for plots and also to grab TFBS positions)")

    args = parser.parse_args()

    # read files
    mat = pd.read_table(args.input, index_col=0)
    mat.columns = [int(x) for x in mat.columns]
    mat = fix_missing_data(mat)

    methyl_positions = get_methyl_positions(mat)

    plot_mutual_info(mat, methyl_positions, args.plot_file, args.input_name, args.positions_file, args.amplicon_name)
