#!/usr/bin/env python3
"""
Walks one or more amplicon-smf project directories (each containing its own results/
folder) and collects every simple stats/*.txt file - tidy, headerless, two-column
"metric\tvalue" files like {sample}.bwameth.contig_filtered.stats.txt,
{sample}.nuc_len_qc.stats.txt, {sample}.background_cpg_methylation.sample_level.txt,
{sample}.duplication_rate.stats.txt - into one long-format table:

    project_dir, experiment, sample, check, metric, value

Wide/headered per-amplicon breakdown files (e.g. *.amplicon_level.txt) don't fit the
two-column shape and are skipped automatically. There's no hardcoded list of checks or
metric names: any future step that writes a tidy metric\tvalue file into a sample's
stats/ directory is picked up the next time this is run, no changes needed here.

Usage:
    python collect_sample_stats.py --dirs /path/to/project_a /path/to/project_b --output all_stats.tsv
"""
import argparse
import glob
import os
import pandas as pd


def is_tidy_metric_value_file(path):
    '''
    A "tidy" stats file is a headerless TSV where every line has exactly two
    tab-separated fields and the second parses as a number. Wide/headered tables (e.g.
    amplicon-level breakdowns, whose header row's second field is a text column name)
    fail this check and are skipped.
    '''
    with open(path) as f:
        lines = [line.rstrip('\n') for line in f if line.strip()]

    if not lines:
        return False

    for line in lines:
        fields = line.split('\t')
        if len(fields) != 2:
            return False
        try:
            float(fields[1])
        except ValueError:
            return False

    return True


def derive_check_name(stats_path, sample):
    '''
    {sample}.bwameth.contig_filtered.stats.txt -> bwameth.contig_filtered.stats
    {sample}.duplication_rate.stats.txt        -> duplication_rate.stats
    '''
    fname = os.path.basename(stats_path)
    prefix = sample + '.'
    if fname.startswith(prefix):
        fname = fname[len(prefix):]
    if fname.endswith('.txt'):
        fname = fname[:-len('.txt')]
    return fname


def collect_stats(project_dirs):
    rows = []

    for project_dir in project_dirs:
        pattern = os.path.join(project_dir, 'results', '*', '*', 'stats', '*.txt')
        for stats_path in sorted(glob.glob(pattern)):
            if not is_tidy_metric_value_file(stats_path):
                print('Skipping (not a tidy metric/value file): {}'.format(stats_path))
                continue

            # stats_path = .../results/{experiment}/{sample}/stats/{file}
            stats_dir = os.path.dirname(stats_path)
            sample_dir = os.path.dirname(stats_dir)
            experiment_dir = os.path.dirname(sample_dir)
            sample = os.path.basename(sample_dir)
            experiment = os.path.basename(experiment_dir)
            check = derive_check_name(stats_path, sample)

            with open(stats_path) as f:
                for line in f:
                    if not line.strip():
                        continue
                    metric, value = line.rstrip('\n').split('\t')
                    rows.append({
                        'project_dir': project_dir,
                        'experiment': experiment,
                        'sample': sample,
                        'check': check,
                        'metric': metric,
                        'value': float(value),
                    })

    return pd.DataFrame(rows, columns=['project_dir', 'experiment', 'sample', 'check', 'metric', 'value'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect simple per-sample stats files across one or more amplicon-smf project directories into one long-format table')
    parser.add_argument('--dirs', nargs='+', required=True, help='One or more project directories, each containing its own results/ folder')
    parser.add_argument('--output', required=True, help='Path to write the collected long-format table (TSV)')
    args = parser.parse_args()

    stats_df = collect_stats(args.dirs)
    stats_df.to_csv(args.output, sep='\t', index=False, header=True)

    n_files = stats_df[['project_dir', 'experiment', 'sample', 'check']].drop_duplicates().shape[0] if len(stats_df) else 0
    print('Collected {} metric rows from {} stats files across {} directories -> {}'.format(
        len(stats_df), n_files, len(args.dirs), args.output))
