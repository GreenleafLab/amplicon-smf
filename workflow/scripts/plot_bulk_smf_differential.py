import pandas as pd
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
plt.switch_backend('agg')

from common import plot_bulk_smf_trace, load_tfbs_positions

def plot_bulk_differential(smf_tbl1, smf_tbl2, plots, output_tbl, positions_file, ylim):
    smf_df1 = pd.read_table(smf_tbl1)
    smf_df2 = pd.read_table(smf_tbl2)

    smf_df = smf_df1.merge(smf_df2, on=['chr', 'start'])
    smf_df['delta_SMF'] = smf_df['SMF_x'] - smf_df['SMF_y']

    amplicons = smf_df['chr'].unique()

    pos_dict = None
    if positions_file:
        pos_dict = load_tfbs_positions(positions_file)

    ylim = tuple(map(float, ylim.split(',')) if ylim else None)

    for amplicon in amplicons:
        df_subset = smf_df.loc[smf_df['chr'] == amplicon]

        positions = pos_dict[amplicon] if pos_dict else None
        
        plot_bulk_smf_trace(df_subset, plots, smf_col='delta_SMF', title=amplicon, ylabel='delta %SMF', fillbetween=positions, ylim=ylim)

    smf_df[['chr', 'start', 'delta_SMF']].to_csv(output_tbl, sep='\t', header=True, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input1", dest="input1", type=str, help="Prefix for the first SMF table")
    parser.add_argument("--input2", dest="input2", type=str, help="Prefix for the second SMF table")
    parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged/differential output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--ylim", dest="ylim", type=str, default=None, help="Fix bounds of plots for all samples (e.g. -20,50)")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        plot_bulk_differential(args.input1, args.input2, plots, args.output, args.positions_file, args.ylim)
