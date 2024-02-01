#!/usr/bin/env python3
"""
Script to take a bam file resulting from bwa -a (all alignments) and filter to only alignments that stem from the same contig.
The reason we have to do this is that bwa sometimes gets confused on some our our amplicons (e.g. 0x-5x) that have F and R reads
both ambiguously mapping to >1 amplicon.
We will iterate through all sets of alignments corresponding to the same sequencing cluster, check whether the primary alignment has 
matching contigs, if so write, if not find the alternate pair with matching contigs and the highest alignment score (greater than
some threshold) and write that.
"""

import sys
import argparse
import pysam
import re
import pandas as pd
from matplotlib import pyplot as plt
from itertools import product

from common import parse_fasta

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="Input bam or sam file (must end in .bam or .sam)")
    parser.add_argument("--out", required=True, help="Name for output sam file")
    parser.add_argument("--problem", required=False, help="Name for output file for problematic reads")
    parser.add_argument("--min_len_threshold1", required=False, default=140, type = int, help="Per-read length threshold to call properly aligned read")
    parser.add_argument("--min_len_threshold2", required=False, default=140, type = int, help="Per-read length threshold to call properly aligned read")
    parser.add_argument("--min_as_frac", required=False, default=0.95, type = float, help="Per-read alignment score threshold, as a function of aligned length")
    parser.add_argument('--ignore_bounds', dest='ignore_bounds', action='store_true')
    parser.add_argument("--amplicon", required=False, type = str, help="Amplicon fasta file")

    args = parser.parse_args()
    return args

def get_primary_alignment(aln_list):
    to_ret = []
    for a in aln_list:
        if not a.is_secondary:
            to_ret.append(a)
    return to_ret

def get_contig(aln):
    return aln.reference_name

def has_matching_contigs(aln_pair):
    assert(len(aln_pair)==2)
    # return aln_pair[0].reference_name == aln_pair[1].reference_name
    return get_contig(aln_pair[0]) == get_contig(aln_pair[1])

def get_matching_contigs(aln_list):
    read1s = [a for a in aln_list if a.is_read1]    
    read2s = [a for a in aln_list if a.is_read2]
    
    combinations = product(read1s, read2s)
    
    valid_combos = []
    
    for c in combinations:
        if has_matching_contigs(c):
            valid_combos.append(c)
    
    return valid_combos

def score_single(aln):
    return aln.get_tag('AS')

def score_pair(aln_pair):
    assert(len(aln_pair)==2)
    return score_single(aln_pair[0]) + score_single(aln_pair[1])

def get_read_len(aln):
    return len(aln.query_sequence)

def get_aln_len(aln_pair):
    assert(len(aln_pair)==2)
    return get_read_len(aln_pair[0]) + get_read_len(aln_pair[1])

def choose_best_aln_per_contig(aln_list):
    contig_to_aln_dict = {}

    for pair in aln_list:
        contig = get_contig(pair[0])
        if contig in contig_to_aln_dict.keys():
            contig_to_aln_dict[contig].append(pair)
        else:
            contig_to_aln_dict[contig] = [pair]
    
    to_ret = []
    
    for aln_lst in contig_to_aln_dict.values():
        if len(aln_lst) == 1:
            to_ret.append(aln_lst[0])
        else:
            lens = [get_aln_len(pair) for pair in aln_lst]
            top_len = max(lens)
            longest_aln_pair = aln_lst[lens.index(top_len)]
            to_ret.append(longest_aln_pair)
            
    return to_ret

def get_aln_bounds(aln_pair):
    aln_positions = aln_pair[0].get_reference_positions() + aln_pair[1].get_reference_positions()
    return min(aln_positions), max(aln_positions)

def verify_alignment_bounds(aln_pair, fa_dict, perf_primer_len=30):
    contig = get_contig(aln_pair[0])
    if contig:
        amplicon = fa_dict.get(contig)
        aln_bounds = get_aln_bounds(aln_pair)
        return (aln_bounds[0] < perf_primer_len) and (aln_bounds[1] > len(amplicon) - perf_primer_len)
    else:
        return False

def is_valid_alignment(aln_pair, args, fa_dict):
    '''
    Want to test a few things here
    First, make sure that the query sequence aligns enough
    This is a bit tricky, since we sequence variable amounts, but I think thresholding on like 120 will generally be fine, although we should remember this
    Second, make sure that the AS is a high fraction of the query_len
    This is important because we want to make sure the AS is high and that most of the read actually did align
    We can't use this in isolation, since short partial reads will be included, which is why we need both
    Should also check both reads separately, since a good F might compensate for a bad R
    Finally, we want to check that the pair spans the entirety of the amplicon
    This is a bit tricky, since this is variable within an experiment, between, etc.
    Could just check starting position? I think for now we will just ask whether the start of the read (easy) and the end of the read (hard, rev pos plus len) are within perf primer len of the amplicon
    This will require passing in the fasta but whatever. Not perfect, but gets the job done for now
    '''
    assert(len(aln_pair)==2)

    passes_min_len = (get_read_len(aln_pair[0]) > args.min_len_threshold1) and (get_read_len(aln_pair[1]) > args.min_len_threshold2)
    passes_as_frac = (score_single(aln_pair[0]) / get_read_len(aln_pair[0]) > args.min_as_frac) and (score_single(aln_pair[1]) / get_read_len(aln_pair[1]) > args.min_as_frac)
    passes_total_covered_length = args.ignore_bounds or verify_alignment_bounds(aln_pair, fa_dict) # can override bounds checking

    return passes_min_len and passes_as_frac and passes_total_covered_length

def fix_read(read, mate):
    '''
    Change the properties of the read such that it doesn't get filtered out at the next step(s)
    '''
    read.mapping_quality = max(50, read.mapping_quality)
    read.is_proper_pair = True
    read.is_secondary = False
    read.is_qcfail = False
    read.is_duplicate = False
    read.next_reference_start = mate.reference_start
    read.next_reference_name = '='

def filter_contigs(bam_file, out, problem_reads, args):
    """
    Iterate through bam file, collect all reads corresponding to same cluster, check primary alignment,
    if not good find best alternate alignment and write
    """

    read_types = {'primary': 0, 'secondary': 0, 'neither': 0}

    # read fasta
    fa_dict = parse_fasta(args.amplicon)

    # construct dict of all alignments corresponding to same read
    read_dict = {}

    for read in bam_file:
        cluster = read.query_name
        if cluster in read_dict:
            read_dict[cluster].append(read)
        else:
            read_dict[cluster] = [read]
            
    # iterate through groupings
    for aln_list in read_dict.values():
        matching_contigs = get_matching_contigs(aln_list)
        filtered_matching_contigs = choose_best_aln_per_contig(matching_contigs)
        if len(filtered_matching_contigs) > 0:
            scores = [score_pair(pair) for pair in filtered_matching_contigs]
            # if len(scores) != len(set(scores)):
            #     print(matching_contigs)
            best_score = max(scores)
            best_matching_alignment = filtered_matching_contigs[scores.index(best_score)]
            if is_valid_alignment(best_matching_alignment, args, fa_dict):
            # if best_score > args.thresh:
                read1, read2 = best_matching_alignment
                # futz with mapping quality to make sure it doesn't get filtered out in the next step
                fix_read(read1, read2)
                fix_read(read2, read1)
                out.write(read1)
                out.write(read2)
                read_types['secondary'] += 1
            else:
                read_types['neither'] += 1
                for (read1, read2) in filtered_matching_contigs:
                    problem_reads.write(read1)
                    problem_reads.write(read2)

        # primary_alignment = get_primary_alignment(aln_list)
        # if has_matching_contigs(primary_alignment):
        #     read1, read2 = primary_alignment
        #     out.write(read1)
        #     out.write(read2)
        #     read_types['primary'] += 1
        # else:
        #     matching_contigs = get_matching_contigs(aln_list)
        #     if len(matching_contigs) > 0:
        #         scores = [score_pair(pair) for pair in matching_contigs]
        #         best_score = max(scores)
        #         best_matching_alignment = matching_contigs[scores.index(best_score)]
        #         if best_score > args.thresh:
        #             read1, read2 = best_matching_alignment
        #             # futz with mapping quality to make sure it doesn't get filtered out in the next step
        #             read1.mapping_quality = max(50, read1.mapping_quality)
        #             read2.mapping_quality = max(50, read2.mapping_quality)
        #             read1.is_proper_pair = True
        #             read2.is_proper_pair = True
        #             out.write(read1)
        #             out.write(read2)
        #             read_types['secondary'] += 1
        #         else:
        #             read_types['neither'] += 1
        #     else:
        #         read_types['neither'] += 1


    print(read_types)
    return None


def run_filter():
    args = argparser()

    print(args.min_len_threshold1)
    print(args.min_len_threshold2)
    print(args.ignore_bounds)

    # If a sam/bam file is specified, use it, otherwise use stdin
    if args.bam:
        if args.bam.endswith(".bam"):
            mysam = pysam.AlignmentFile(args.bam, "rb")
        elif args.bam.endswith(".sam"):
            mysam = pysam.AlignmentFile(args.bam, "r")
    else:
        # sys.stdin.isatty() checks that there is some data coming in through stdin
        assert sys.stdin.isatty() is False, "You didn't pipe me any data or specify an input "\
                                            "file! Use 'mark-nonconverted-reads.py -h' for more "\
                                            "usage information."
        mysam = pysam.AlignmentFile("-", "r")

    # If an output bam is specified, write there, otherwise write to stdout
    if args.out:
        if args.out.endswith(".bam"):
            out = pysam.AlignmentFile(args.out, "wb", template = mysam)
        elif args.out.endswith(".sam"):
            out = pysam.AlignmentFile(args.out, "w", template = mysam)
    else:
        out = pysam.AlignmentFile("-", "w", template = mysam)

    # If an output problematic bam file is specified, write there, otherwise write to stdout
    if args.problem:
        if args.problem.endswith(".bam"):
            problem = pysam.AlignmentFile(args.problem, "wb", template = mysam)
        elif args.problem.endswith(".sam"):
            problem = pysam.AlignmentFile(args.problem, "w", template = mysam)
    else:
        problem = pysam.AlignmentFile("-", "w", template = mysam)
    
    filter_contigs(mysam, out, problem, args)
    

if __name__ == "__main__":
    run_filter()
