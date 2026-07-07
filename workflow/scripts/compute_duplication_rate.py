#!/usr/bin/env python3
"""
Reduces a sample's per-amplicon amplicon_stats.txt (amplicon, total_reads, observed_states,
reads_per_state) into a single per-sample duplication-rate summary, written in the same
tidy metric\tvalue (no header) convention as the other simple files in stats/, so it can be
picked up by a generic cross-experiment stats collector alongside everything else.
"""
import argparse
import pandas as pd


def compute_duplication_rate(input_path, output_path):
    df = pd.read_table(input_path)

    # Pooled (count-weighted) rate - sum of totals over sum of uniques - rather than an
    # unweighted mean of each amplicon's reads_per_state. Amplicons with very few reads/
    # states would otherwise get equal say to well-covered ones and add noise to the ratio.
    pooled_reads_per_state = df['total_reads'].sum() / df['observed_states'].sum()

    with open(output_path, 'w') as f:
        f.write('avg_reads_per_state\t{:.4f}\n'.format(pooled_reads_per_state))
        f.write('total_reads\t{}\n'.format(df['total_reads'].sum()))
        f.write('total_observed_states\t{}\n'.format(df['observed_states'].sum()))
        f.write('n_amplicons\t{}\n'.format(len(df)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reduce per-amplicon amplicon_stats.txt into a single per-sample duplication rate')
    parser.add_argument('--input', required=True, help='Path to {sample}.amplicon_stats.txt')
    parser.add_argument('--output', required=True, help='Path to output tidy per-sample stats file')
    args = parser.parse_args()

    compute_duplication_rate(args.input, args.output)
