import pandas as pd
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
from functools import reduce

from common import plot_bulk_smf_trace, load_tfbs_positions, add_fillbetween

def plot_and_join_smf(smf_tbls_lst, plots, output_tbl, positions_file, ylim, ylabel):
    smf_tbls = [pd.read_table(t) for t in smf_tbls_lst]
    df = pd.concat(smf_tbls, axis=0)

    amplicons = df['chr'].unique()

    pos_dict = None
    if positions_file:
        pos_dict = load_tfbs_positions(positions_file)

    ylim = tuple(map(float, ylim.split(',')) if ylim else None)

    for amplicon in amplicons:
        df_subset = df.loc[df['chr'] == amplicon].sort_values('start')

        positions = pos_dict[amplicon] if pos_dict else None

        plot_bulk_smf_trace(df_subset, plots, title=amplicon, fillbetween=positions, ylim=ylim)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--inputs", dest="inputs", type=str, nargs='+', help="Space-separated list of SMF tables")
    parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the summary table to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--ylim", dest="ylim", type=str, default=None, help="Fix bounds of plots for all samples (e.g. -20,50)")
    parser.add_argument("--ylabel", dest="ylabel", type=str, default=None, help="Rename vertical axis")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        plot_and_join_smf(args.inputs, plots, args.output, args.positions_file, args.ylim, args.ylabel)
