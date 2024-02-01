import pandas as pd
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
from functools import reduce

from common import plot_bulk_smf_trace, load_tfbs_positions, add_fillbetween

def plot_multiple_conditions(smf_tbls_lst, labels, plots, positions_file, ylim, ylabel, col):
    smf_tbls = []

    for idx in range(len(smf_tbls_lst)):
        smf_df = pd.read_table(smf_tbls_lst[idx])
        if labels:
            cond = labels[idx]
        else:
            cond = 'Cond_{}'.format(idx)
        smf_df['condition'] = cond
        smf_tbls.append(smf_df)

    df = pd.concat(smf_tbls, axis=0).reset_index()

    amplicons = df['chr'].unique()

    pos_dict = None
    if positions_file:
        pos_dict = load_tfbs_positions(positions_file)

    if ylim:
        ylim = tuple(map(float, ylim.split(',')))

    for amplicon in amplicons:
        df_subset = df.loc[df['chr'] == amplicon]

        fig, ax = plt.subplots()
        sns.lineplot(data=df_subset, x='start', y=col, hue='condition', ax=ax, ls='-', marker='.', ms=10)

        ax.set_title(amplicon)

        if ylabel:
            ax.set_ylabel(ylabel)

        if ylim:
            ax.set_ylim(ylim)

        if pos_dict:
            positions = pos_dict[amplicon]
            add_fillbetween(ax, positions)
            plt.legend()

        plt.tight_layout()
        plots.savefig()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--inputs", dest="inputs", type=str, nargs='+', help="Space-separated list of SMF tables")
    parser.add_argument("--labels", dest="labels", type=str, nargs='+', default=None, help="[Optional] Space-separated list of labels corresponding to SMF tables")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--ylim", dest="ylim", type=str, default=None, help="Fix bounds of plots for all samples (e.g. -20,50)")
    parser.add_argument("--ylabel", dest="ylabel", type=str, default=None, help="Rename vertical axis")
    parser.add_argument("--col", dest="col", type=str, default='SMF', help="What is the column name to plot the replicates on")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        plot_multiple_conditions(args.inputs, args.labels, plots, args.positions_file, args.ylim, args.ylabel, args.col)
