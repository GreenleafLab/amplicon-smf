import pandas as pd
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
from functools import reduce

from common import plot_bulk_smf_trace, load_tfbs_positions, add_fillbetween

def plot_smf_replicates(smf_tbls_lst, plots, output_tbl, positions_file, ylim, ylabel, col):
    cols_to_use = ['chr', 'start', col]
    smf_tbls = [pd.read_table(t)[cols_to_use] for t in smf_tbls_lst]
    df = reduce(lambda left, right: pd.merge(left, right, on=['chr','start']), smf_tbls)
    df_melted = df.melt(id_vars=['chr', 'start'], value_name='SMF_reps')

    amplicons = df_melted['chr'].unique()

    pos_dict = None
    if positions_file:
        pos_dict = load_tfbs_positions(positions_file)

    ylim = tuple(map(float, ylim.split(',')) if ylim else None)

    for amplicon in amplicons:
        df_subset = df_melted.loc[df_melted['chr'] == amplicon]

        fig, ax = plt.subplots()
        sns.lineplot(data=df_subset, x="start", y="SMF_reps", ls='-', ax=ax, ms=10, marker='.', err_style="bars", ci=95, color='k') 

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

    # still to do is actually calculate the mean + sem manually and write the file back, maybe?
    # smf_df[['chr', 'start', 'delta_SMF']].to_csv(output_tbl, sep='\t', header=True, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--inputs", dest="inputs", type=str, nargs='+', help="Space-separated list of SMF tables")
    parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the summary table to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--ylim", dest="ylim", type=str, default=None, help="Fix bounds of plots for all samples (e.g. -20,50)")
    parser.add_argument("--ylabel", dest="ylabel", type=str, default=None, help="Rename vertical axis")
    parser.add_argument("--col", dest="col", type=str, default='SMF', help="What is the column name to plot the replicates on")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        plot_smf_replicates(args.inputs, plots, args.output, args.positions_file, args.ylim, args.ylabel, args.col)
