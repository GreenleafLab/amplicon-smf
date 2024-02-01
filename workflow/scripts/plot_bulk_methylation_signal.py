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

from common import plot_bulk_smf_trace


def plot_bulk_methylation(input_prefix, amplicon_fa, plots, thresh, include_cpg, no_endog_meth, save_individual_png):
    # load amplicon fasta
    amplicon_to_seq_dict = {}
    
    for r in list(SeqIO.parse(amplicon_fa, "fasta")):
        # print(r)
        amplicon_to_seq_dict[r.id] = r.seq

    amplicons = amplicon_to_seq_dict.keys()

    # calculate GpC positions (also CpG)
    gpc_pos_dict = {}
    cpg_pos_dict = {}
    for amplicon in amplicons:
        gpc_pos_dict[amplicon] = [i for i, b in enumerate(amplicon_to_seq_dict[amplicon]) if b=='C' and amplicon_to_seq_dict[amplicon][max(i-1, 0)]=='G']
        cpg_pos_dict[amplicon] = [i for i, b in enumerate(amplicon_to_seq_dict[amplicon]) if b=='C' and amplicon_to_seq_dict[amplicon][min(i+1, len(amplicon_to_seq_dict[amplicon])-1)]=='G']

    # load bedgraphs
    bedgraph_chg_str = input_prefix + '_CHG.bedGraph'
    bedgraph_chh_str = input_prefix + '_CHH.bedGraph'
    bedgraph_cpg_str = input_prefix + '_CpG.bedGraph'

    bedgraph_chg = pd.read_table(bedgraph_chg_str, skiprows=1, header=None, names=['chr','start','end','pct','meth','unmeth'])
    bedgraph_chh = pd.read_table(bedgraph_chh_str, skiprows=1, header=None, names=['chr','start','end','pct','meth','unmeth'])
    bedgraph_cpg = pd.read_table(bedgraph_cpg_str, skiprows=1, header=None, names=['chr','start','end','pct','meth','unmeth'])

    output_tables = []

    # iterate through amplicons, merge the tables, plot, and join tables
    for amplicon in amplicons:
        gpcs = gpc_pos_dict[amplicon]
        cpgs = cpg_pos_dict[amplicon]

        # join the tables
        if include_cpg or no_endog_meth:
            all_c_df = bedgraph_chg.loc[bedgraph_chg['chr']==amplicon].append([bedgraph_chh.loc[bedgraph_chh['chr']==amplicon],bedgraph_cpg.loc[bedgraph_cpg['chr']==amplicon]]).sort_values('start')
        else:
            all_c_df = bedgraph_chg.loc[bedgraph_chg['chr']==amplicon].append(bedgraph_chh.loc[bedgraph_chh['chr']==amplicon]).sort_values('start')

        # process the table
        all_c_df['SMF'] = 100 - all_c_df['pct']

        # filter on low counts
        all_c_df['total_reads'] = all_c_df['meth'] + all_c_df['unmeth']
        all_c_df = all_c_df.loc[all_c_df['total_reads'] > thresh] 

        # check if data 
        if len(all_c_df) > 0:
            # subset by C type
            all_gpc_df = all_c_df.loc[all_c_df.start.isin(gpcs)].copy()
            all_cpg_df = all_c_df.loc[all_c_df.start.isin(cpgs)].copy()
            all_other_c_df = all_c_df.loc[~all_c_df.start.isin(gpcs+cpgs)].copy()

            # now construct the thing we will plot and return, based on whether or not we included CpG MTase
            if include_cpg:
                cs_to_plot = all_gpc_df.append(all_cpg_df).sort_values('start')
            else:
                cs_to_plot = all_gpc_df

            # plot bulk methyl signal
            plot_bulk_smf_trace(cs_to_plot, plots, title=amplicon, ylim=(0,103))

            # also make a plot with only GpC and only CpG if also using CpG MTase, to see if everything worked
            if include_cpg:
                plot_bulk_smf_trace(all_gpc_df, plots, title='{}: GpC only'.format(amplicon), ylim=(0,103))
                plot_bulk_smf_trace(all_cpg_df, plots, title='{}: CpG only'.format(amplicon), ylim=(0,103))

            # plot the average signal from GpC and non-GpC
            all_cpg_df['context'] = 'CpG'
            all_gpc_df['context'] = 'GpC'
            all_other_c_df['context'] = 'C'
            # all_gpc_df['GpC'] = True
            # all_other_c_df['GpC'] = False

            fig, ax = plt.subplots()
            sns.barplot(data=all_other_c_df.append([all_cpg_df,all_gpc_df]), x='context', y='pct')
            # sns.barplot(data=all_gpc_df.append(all_other_c_df), x='GpC', y='pct')

            ax.set_ylabel('%Methylation')

            plt.tight_layout()
            plots.savefig()
            plt.close()

            # return merged table
            columns = ['chr', 'start', 'SMF', 'total_reads']
            # output_tables.append(all_gpc_df[columns])
            output_tables.append(cs_to_plot[columns])
        else:
            print('No data for {}'.format(amplicon))

    pd.concat(output_tables).to_csv('{}.all_GpC.txt'.format(input_prefix), sep='\t', header=True, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Prefix for the bedgraph files from MethylDackel extract (with suffix e.g. _CHG.bedGraph")
    parser.add_argument("--amplicon", dest="amplicon_fa", type=str, help="FASTA file containing the amplicons that the reads were aligned to")
    # parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--thresh", dest="thresh", type=int, default=100, help="How many times a C has to be observed before it is included in the output/plot")
    parser.add_argument('--include_cpg', dest='include_cpg', action='store_true', help="Whether to include CpGs along with GpCs")
    parser.add_argument('--no_endog_meth', dest='no_endog_meth', action='store_true', help="Whether GpCpGs are safe to consider or not")
    parser.add_argument('--save_png', dest='save_individual_png', action='store_true', help="Save each individual plot as separate png")
    parser.set_defaults(save_individual_png=False)

    args = parser.parse_args()

    print(args.include_cpg)

    with PdfPages(args.plot_file) as plots:
        # plot_bulk_methylation(args.input, args.amplicon_fa, args.output, plots, args.thresh, args.save_individual_png)
        plot_bulk_methylation(args.input, args.amplicon_fa, plots, args.thresh, args.include_cpg, args.no_endog_meth, args.save_individual_png)
