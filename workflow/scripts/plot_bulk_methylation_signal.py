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


def load_amplicon_context_positions(amplicon_fa):
    '''
    Parses the amplicon fasta and returns, per amplicon, the 0-based positions of GpC Cs and
    CpG Cs. The two sets overlap at ambiguous GCG sites (a C that is both preceded and followed
    by G) - callers that need to treat those specially should intersect/subtract the sets
    themselves (see plot_background_methylation_qc).
    '''
    amplicon_to_seq_dict = {}

    for r in list(SeqIO.parse(amplicon_fa, "fasta")):
        amplicon_to_seq_dict[r.id] = r.seq

    gpc_pos_dict = {}
    cpg_pos_dict = {}
    for amplicon in amplicon_to_seq_dict:
        gpc_pos_dict[amplicon] = [i for i, b in enumerate(amplicon_to_seq_dict[amplicon]) if b=='C' and amplicon_to_seq_dict[amplicon][max(i-1, 0)]=='G']
        cpg_pos_dict[amplicon] = [i for i, b in enumerate(amplicon_to_seq_dict[amplicon]) if b=='C' and amplicon_to_seq_dict[amplicon][min(i+1, len(amplicon_to_seq_dict[amplicon])-1)]=='G']

    return amplicon_to_seq_dict, gpc_pos_dict, cpg_pos_dict


def load_bedgraphs(input_prefix):
    bedgraph_chg = pd.read_table(input_prefix + '_CHG.bedGraph', skiprows=1, header=None, names=['chr','start','end','pct','meth','unmeth'])
    bedgraph_chh = pd.read_table(input_prefix + '_CHH.bedGraph', skiprows=1, header=None, names=['chr','start','end','pct','meth','unmeth'])
    bedgraph_cpg = pd.read_table(input_prefix + '_CpG.bedGraph', skiprows=1, header=None, names=['chr','start','end','pct','meth','unmeth'])
    return bedgraph_chg, bedgraph_chh, bedgraph_cpg


def plot_bulk_methylation(input_prefix, amplicon_fa, plots, thresh, include_cpg, no_endog_meth, deaminase, save_individual_png):
    amplicon_to_seq_dict, gpc_pos_dict, cpg_pos_dict = load_amplicon_context_positions(amplicon_fa)
    amplicons = amplicon_to_seq_dict.keys()

    bedgraph_chg, bedgraph_chh, bedgraph_cpg = load_bedgraphs(input_prefix)

    output_tables = []

    # iterate through amplicons, merge the tables, plot, and join tables
    for amplicon in amplicons:
        print(amplicon)
        gpcs = gpc_pos_dict[amplicon]
        cpgs = cpg_pos_dict[amplicon]

        # print(bedgraph_chg.loc[bedgraph_chg['chr']==amplicon])

        # join the tables
        # figure out here what kinds of Cs we want! 
        # here, we are just joining the tables, we still have to go in and fish out the GpCs or CpGs or Cs later
        # so here it's just do we want to EXCLUDE CpG or not, which we only do if we JUST want GpC and there might be endog meth
        if include_cpg or no_endog_meth or deaminase: 
            all_c_df = pd.concat([bedgraph_chg.loc[bedgraph_chg['chr']==amplicon], bedgraph_chh.loc[bedgraph_chh['chr']==amplicon], bedgraph_cpg.loc[bedgraph_cpg['chr']==amplicon]]).sort_values('start')
        else:
            all_c_df = pd.concat([bedgraph_chg.loc[bedgraph_chg['chr']==amplicon], bedgraph_chh.loc[bedgraph_chh['chr']==amplicon]]).sort_values('start')

        print(all_c_df)

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

            # now construct the thing we will plot and return, based on whether or not we included CpG MTase or deaminase
            if deaminase:
                # all Cs!
                cs_to_plot = all_c_df
            else:
                # check whether we want to include CpG or not
                if include_cpg:
                    cs_to_plot = pd.concat([all_gpc_df, all_cpg_df]).sort_values('start')
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
            sns.barplot(data=pd.concat([all_other_c_df,all_cpg_df,all_gpc_df]), x='context', y='pct')
            # sns.barplot(data=all_other_c_df.append([all_cpg_df,all_gpc_df]), x='context', y='pct')
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


def plot_background_methylation_qc(input_prefix, amplicon_fa, plots, summary_plot, thresh,
                                    amplicon_stats_path, sample_stats_path):
    '''
    Always-on background methylation QC, independent of include_cpg/no_endog_meth/deaminase:
    per-amplicon traces of (1) "clean" CpG-context methylation (C followed by G, excluding
    ambiguous GCG sites where the C is also preceded by G - see amplicon-smf/CLAUDE.md for why
    those are excluded) and (2) "other C" methylation (not preceded or followed by G, matching
    the definition used to filter reads in mark-nonconverted-reads-and-plot.py), so background/
    endogenous CpG methylation can be confirmed negligible regardless of which enzyme mode a
    sample was run in. Also writes amplicon-level and sample-level summary stats tables.
    '''
    amplicon_to_seq_dict, gpc_pos_dict, cpg_pos_dict = load_amplicon_context_positions(amplicon_fa)
    bedgraph_chg, bedgraph_chh, bedgraph_cpg = load_bedgraphs(input_prefix)

    amplicon_rows = []
    summary_tables = []

    for amplicon in amplicon_to_seq_dict:
        gpcs = set(gpc_pos_dict[amplicon])
        cpgs = set(cpg_pos_dict[amplicon])
        clean_cpgs = cpgs - gpcs

        all_c_df = pd.concat([
            bedgraph_chg.loc[bedgraph_chg['chr'] == amplicon],
            bedgraph_chh.loc[bedgraph_chh['chr'] == amplicon],
            bedgraph_cpg.loc[bedgraph_cpg['chr'] == amplicon],
        ]).sort_values('start')

        all_c_df['total_reads'] = all_c_df['meth'] + all_c_df['unmeth']
        all_c_df = all_c_df.loc[all_c_df['total_reads'] > thresh]

        if len(all_c_df) == 0:
            print('No data for {} (background methylation QC)'.format(amplicon))
            continue

        clean_cpg_df = all_c_df.loc[all_c_df.start.isin(clean_cpgs)].copy()
        other_c_df = all_c_df.loc[~all_c_df.start.isin(gpcs | cpgs)].copy()

        if len(clean_cpg_df) > 0:
            plot_bulk_smf_trace(clean_cpg_df, plots, title='{}: background CpG (non-GCG)'.format(amplicon),
                                 smf_col='pct', ylabel='%methylated', ylim=(0, 103))
        if len(other_c_df) > 0:
            plot_bulk_smf_trace(other_c_df, plots, title='{}: other C (non-GpC, non-CpG)'.format(amplicon),
                                 smf_col='pct', ylabel='%methylated', ylim=(0, 103))

        amplicon_rows.append({
            'amplicon': amplicon,
            'cpg_methylation_pct': clean_cpg_df['pct'].mean() if len(clean_cpg_df) else float('nan'),
            'other_c_methylation_pct': other_c_df['pct'].mean() if len(other_c_df) else float('nan'),
            'n_cpg_positions': len(clean_cpg_df),
            'n_other_c_positions': len(other_c_df),
        })

        if len(clean_cpg_df) > 0:
            summary_tables.append(clean_cpg_df.assign(amplicon=amplicon, context='CpG (non-GCG)'))
        if len(other_c_df) > 0:
            summary_tables.append(other_c_df.assign(amplicon=amplicon, context='other C'))

    # amplicon-level table stays wide/headered (for plotting straight from the file)
    amplicon_df = pd.DataFrame(amplicon_rows)
    amplicon_df.to_csv(amplicon_stats_path, sep='\t', index=False, header=True)

    # sample-level file is tidy metric\tvalue (no header), matching every other simple
    # stats file in stats/ (mapped/unmapped, failure_modes, nuc_len_qc), so a generic
    # collector can grab it without knowing anything about what it means
    with open(sample_stats_path, 'w') as f:
        f.write('cpg_methylation_pct\t{:.4f}\n'.format(amplicon_df['cpg_methylation_pct'].mean()))
        f.write('other_c_methylation_pct\t{:.4f}\n'.format(amplicon_df['other_c_methylation_pct'].mean()))
        f.write('n_amplicons\t{}\n'.format(len(amplicon_df)))

    if summary_tables:
        fig, ax = plt.subplots(figsize=(max(6, 0.5 * amplicon_df['amplicon'].nunique()), 5))
        sns.barplot(data=pd.concat(summary_tables), x='amplicon', y='pct', hue='context', ax=ax)
        ax.set_ylabel('%methylated')
        ax.set_xlabel('')
        plt.xticks(rotation=90)
        plt.tight_layout()
        fig.savefig(summary_plot)
        plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Prefix for the bedgraph files from MethylDackel extract (with suffix e.g. _CHG.bedGraph")
    parser.add_argument("--amplicon", dest="amplicon_fa", type=str, help="FASTA file containing the amplicons that the reads were aligned to")
    # parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--thresh", dest="thresh", type=int, default=100, help="How many times a C has to be observed before it is included in the output/plot")
    parser.add_argument('--include_cpg', dest='include_cpg', action='store_true', help="Whether to include CpGs along with GpCs")
    parser.add_argument('--no_endog_meth', dest='no_endog_meth', action='store_true', help="Whether GpCpGs are safe to consider or not")
    parser.add_argument('--deaminase', dest='deaminase', action='store_true', help="Whether this was a deaminase experiment")
    parser.add_argument('--save_png', dest='save_individual_png', action='store_true', help="Save each individual plot as separate png")
    parser.set_defaults(save_individual_png=False)
    parser.add_argument("--background_plot", dest="background_plot_file", type=str, help="Path to output background CpG/other-C trace plots file (always generated, regardless of include_cpg/no_endog_meth/deaminase)")
    parser.add_argument("--background_summary_plot", dest="background_summary_plot_file", type=str, help="Path to output background methylation per-amplicon summary plot")
    parser.add_argument("--background_amplicon_stats", dest="background_amplicon_stats_file", type=str, help="Path to output amplicon-level background methylation stats table")
    parser.add_argument("--background_sample_stats", dest="background_sample_stats_file", type=str, help="Path to output sample-level (one row) background methylation stats table")

    args = parser.parse_args()

    print(args.include_cpg)

    with PdfPages(args.plot_file) as plots:
        # plot_bulk_methylation(args.input, args.amplicon_fa, args.output, plots, args.thresh, args.save_individual_png)
        plot_bulk_methylation(args.input, args.amplicon_fa, plots, args.thresh, args.include_cpg, args.no_endog_meth, args.deaminase, args.save_individual_png)

    with PdfPages(args.background_plot_file) as background_plots:
        plot_background_methylation_qc(args.input, args.amplicon_fa, background_plots, args.background_summary_plot_file,
                                        args.thresh, args.background_amplicon_stats_file, args.background_sample_stats_file)
