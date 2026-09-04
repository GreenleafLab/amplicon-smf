# this basically loads in all the binding model calls and computes all the metrics that I've done in the past
# it also aggregates the per amplicon results file that I'll need for eg plotting the partition function model

import pandas as pd
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns
from sklearn.cluster import KMeans
import pickle
from copy import deepcopy
import gc

from common import load_single_molecule_matrix, load_tfbs_positions

def compute_remodel_metric(df, positions_list):
    last_teto = max(num for sublist in positions_list for num in sublist[:2]) if positions_list else np.inf
    subset_df = df[[c for c in df.columns if c > last_teto]]
    accessibility = (1-subset_df).mean(axis=1)
    return accessibility

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collates all binding model files into one big file')
    parser.add_argument("--singlemol_dir", type=str, help="Path to highlevel folder containing multiple single molecule matrics")
    parser.add_argument("--classifications_dir", type=str, help="Path to highlevel folder containing multiple single molecule classifications")
    parser.add_argument("--samples", type=str, nargs='+', default=None, help="Space separated list of samples to aggregate")
    parser.add_argument("--amplicons", type=str, nargs='+', default=None, help="Space separated list of amplicons to use")
    parser.add_argument("--output_dir", type=str, help="Output director for file")
    parser.add_argument("--out_prefix", type=str, help="Output prefix for file")
    parser.add_argument("--pos_dicts", type=str, default='/oak/stanford/groups/wjg/bgrd/projects/smf/220829_P026_opJS45/amplicon-info/opJS4.positions.long.txt,/oak/stanford/groups/wjg/bgrd/projects/smf/220829_P026_opJS45/amplicon-info/opJS5.positions.long.txt', help="comma-separated list of positions files")

    args = parser.parse_args()

    # load tfbs positions
    dicts = [load_tfbs_positions(pos_file) for pos_file in args.pos_dicts.split(',')]
    positions_dict = {k: v for d in dicts for k, v in d.items()}

    all_molecule_annotations_list = []

    # load the single molecule classifications
    missing_files = []

    for samp in args.samples:
        tmp_dict = {}
        
        for amplicon in args.amplicons:
            p = path.join(args.classifications_dir, samp,'{}.{}.single_molecule_classification.txt'.format(samp, amplicon))
            if path.exists(p):
                # load df
                df = pd.read_table(p, index_col=0)
                df.index = ['{}-{}-{}'.format(samp,amplicon,i) for i in df.index]

                # annotate df
                df['sample'] = samp
                df['amplicon'] = amplicon
                df['background'] = ('b1' if 'b1' in amplicon else 'b2') if 'opJS5' in amplicon else 'b0'

                # add stats
                df['n_bound'] = df.filter(like='tfbs_').sum(axis=1)
                df['tf_bound'] = df['n_bound']>0
                df['n_tfbs'] = df.apply(lambda row: int(row.amplicon.split('_')[1][0]) if 'TetO' in row.amplicon else np.nan, axis=1)

                # add the remodel metric
                p1 = path.join(args.singlemol_dir, samp, 'matrices', '{}.{}.dedup.full_unclustered.matrix'.format(samp, amplicon))
                df1 = load_single_molecule_matrix(p1, every_other=False)
                df1.index = ['{}-{}-{}'.format(samp,amplicon,i) for i in df1.index]

                peripheral_accessibility = compute_remodel_metric(df1, positions_dict[amplicon])
                df['peripheral_accessibility'] = peripheral_accessibility

                # write df
                all_molecule_annotations_list.append(df)
            else:
                if 'CTCF' in amplicon or 'BD24' in amplicon or '0x_TetO' in amplicon or '0xTetO' in amplicon:
                    continue
                else:
                    missing_files.append(p)


    # join and write
    all_molecule_annotations = pd.concat(all_molecule_annotations_list)

    outpath = path.join(args.output_dir, f'{args.out_prefix}.binding_model.txt.gz')

    all_molecule_annotations.to_csv(outpath, sep='\t', header=True, index=True)

    # also now aggregate to make the file that has the per amplicon stats
    def generate_hist_string(series):
        """
        Custom aggregator to create the 'n_bound:frequency' string.
        """
        counts = series.value_counts(normalize=True).sort_index(ascending=False)
        return ','.join([f"{int(n)}:{freq:.4f}" for n, freq in counts.items()])

    per_amp_df = all_molecule_annotations.groupby(['sample', 'amplicon']).agg(
        avg_tf_bound=('n_bound', 'mean'),
        peripheral_accessibility=('peripheral_accessibility', 'mean'),
        molecules=('n_bound', 'count'),
        binding_histogram=('n_bound', generate_hist_string)
    ).reset_index()

    per_amp_df['background'] = per_amp_df['amplicon'].apply(
        lambda x: ('b1' if 'b1' in x else 'b2') if 'opJS5' in x else 'b0'
    )
    per_amp_df['n_tfbs'] = per_amp_df['amplicon'].apply(
        lambda x: int(x.split('_')[1][0]) if 'TetO' in x else 0
    )

    outpath = path.join(args.output_dir, f'{args.out_prefix}.binding_model_amplicon_level.txt.gz')

    per_amp_df.to_csv(outpath, sep='\t', header=True, index=False)

    # also do the per sample thing, which ideally will be partition function stuff but for now will just be binding metrics for 6x I think...
    # also the remodel metric for 4x?
    amp = 'opJS4_6x_TetO_21bp_no_CG'
    per_samp_tf_bound = per_amp_df.loc[per_amp_df.amplicon==amp,['sample', 'avg_tf_bound']].set_index('sample')

    amp = 'opJS4_4x_TetO_21bp_no_CG'
    per_samp_remodel = per_amp_df.loc[per_amp_df.amplicon==amp,['sample', 'peripheral_accessibility']].set_index('sample')
    
    per_samp_df = per_samp_tf_bound.join(per_samp_remodel)
    outpath = path.join(args.output_dir, f'{args.out_prefix}.binding_model_sample_level.txt.gz')
    per_samp_df.to_csv(outpath, sep='\t', index=True, header=True)
