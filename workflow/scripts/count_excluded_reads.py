#!/usr/bin/env python3
"""
Script to take in two sam files, one with included and one with excluded reads
Plot the usage of reads (i.e. how many thrown away)
"""

import sys
import argparse
import pysam
import re
import pandas as pd
from matplotlib import pyplot as plt
from itertools import product

def argparser():
    parser = argparse.ArgumentParser()
    # parser.add_argument("--retained", required=True, help="Input sam file with proper reads")
    # parser.add_argument("--removed", required=True, help="Input sam file with problematic reads")
    parser.add_argument("--stats", required=True, help="Input stats file with mapped and unmapped reads")
    parser.add_argument("--plot", required=True, help="Path to output plot")

    args = parser.parse_args()
    return args

def count_unique_molecules_from_sam_file(sam):
    all_reads = []

    for read in sam:
        all_reads.append(read.query_name)
        
    return len(set(all_reads))

def count_reads():
    args = argparser()

    # retained_sam = pysam.AlignmentFile(args.retained, "r")
    # removed_sam = pysam.AlignmentFile(args.removed, "r")

    # retained_reads = count_unique_molecules_from_sam_file(retained_sam)
    # removed_reads = count_unique_molecules_from_sam_file(removed_sam)

    stats = {}

    with open(args.stats, 'r') as stats_file:
        for line in stats_file:
            if line:
                mapped, val = line.split('\t')
                stats[mapped] = int(val)

    retained_reads, removed_reads = stats['mapped'], stats['unmapped']

    fig, ax = plt.subplots()
    # ax.bar(['Retained', 'Filtered'], [retained_reads, removed_reads])
    ax.bar(['Total', 'Retained'], [retained_reads+removed_reads, retained_reads])
    ax.set_ylabel('Reads')    

    plt.tight_layout()
    plt.savefig(args.plot)

if __name__ == "__main__":
    count_reads()
