# input a matrix (maybe make it do 2, afterwards?)
# do PCA and plot the top N PCs

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns
from sklearn.decomposition import PCA

from common import fix_missing_data


def plot_pca(mat, methyl_positions, plots, name, positions_file, amplicon):
    # do some general formatting
    # colors = cmap([_ for _ in range(n_clusters)])
    title = amplicon if amplicon else ''
    name = name if name else ''

    # fix matrix
    X = fix_missing_data(mat[methyl_positions]).values

    # run_PCA
    pca = PCA(n_components=10)
    pca.fit(X)

    # find number of PCAs that explain 75% variance
    var_frac = 0.75
    exp_var_sum = np.cumsum(pca.explained_variance_ratio_)
    idx = list(exp_var_sum > var_frac).index(True) + 1

    fig, ax = plt.subplots()

    for i,(pc,v) in enumerate(zip(pca.components_, pca.explained_variance_)):
        if i <= idx:
            ax.plot(methyl_positions, v*pc, label=i)

    pos_dict = None
    if positions_file:
        pos_dict = load_tfbs_positions(positions_file)
    if pos_dict and amplicon:
        positions = pos_dict[amplicon]
        add_fillbetween(ax, positions)

    plt.legend()
    plt.xlabel('Position along plasmid')
    plt.ylabel('PC Weight')
    plots.savefig()
    plt.clf()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Path to single molecule matrix")
    # parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--input_name", dest="input_name", type=str, default=None, help="Name of input sample (e.g. -dox)")
    parser.add_argument("--amplicon_name", dest="amplicon_name", type=str, default=None, help="Name of amplicon aligned to (for plots and also to grab TFBS positions)")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        # read files
        mat = pd.read_table(args.input, index_col=0)
        mat.columns = [int(x) for x in mat.columns]

        # grab methyl positions and subset prior to PCA
        methyl_positions = mat.loc[:,mat.mean()>0].columns.tolist()

        plot_pca(mat, methyl_positions, plots, args.input_name, args.positions_file, args.amplicon_name)


