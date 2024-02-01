import pandas as pd
import numpy as np
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
plt.switch_backend('agg')

# from common import plot_bulk_smf_trace, load_tfbs_positions

def make_l2fc_plot(data, sample1, sample2, plots, title=None):
    # make l2fc df
    gb = data.groupby('sample')['cluster'].value_counts(normalize=True)
    gb.name='fraction'
    gb = gb.reset_index()
    gb = gb.pivot(columns='sample', values='fraction', index='cluster')
    gb = gb.fillna(gb.min().min()) # is this the right thing to do?
    gb['fold_change'] = np.log2(gb[sample2] / gb[sample1])
    gb = gb.sort_values('fold_change')

    # plot
    fig, ax = plt.subplots(figsize=(8,4))

    ax.plot(list(range(len(gb))), gb['fold_change'], '.', c='k', ms=10)
    ax.axhline(0, ls='--', c='gray')

    ax.set_xticks(list(range(len(gb))))
    ax.set_xticklabels(gb.index)

    ax.set_ylabel('log2-fold change cluster abundance')
    ax.set_xlabel('Cluster')

    if title:
        ax.set_title(title)

    plt.tight_layout()
    plots.savefig()
    plt.close()

def plot_promoter_cluster_l2fc(promoter_cluster_tbl, sample1, sample2, plots, output_tbl, error_bars):
    promoter_clusters = pd.read_table(promoter_cluster_tbl)

    samp_clusters = promoter_clusters.loc[promoter_clusters['sample'].isin([sample1,sample2])]

    # start making plots
    make_l2fc_plot(samp_clusters, sample1, sample2, plots, title='All Amplicons')

    amplicons = sorted(promoter_clusters['amplicon'].unique())

    for amplicon in amplicons:
        print(amplicon)

        samp_clusters_amp = samp_clusters.loc[samp_clusters['amplicon'] == amplicon]
        
        make_l2fc_plot(samp_clusters_amp, sample1, sample2, plots, title=amplicon)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Path to promoter clustering matrix")
    parser.add_argument("--sample1", dest="sample1", type=str, help="Sample index for the 'denominator' of the l2fc plots")
    parser.add_argument("--sample2", dest="sample2", type=str, help="Sample index for the 'numerator' of the l2fc plots")
    parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged/differential output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--error_bars", dest="error_bars", action='store_true', help="Whether to split up the inputs and compute the l2fc separately on pairs and then use the separate results for error bars (default False, also not implemented yet!)")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        plot_promoter_cluster_l2fc(args.input, args.sample1, args.sample2, plots, args.output, args.error_bars)
