#!/usr/bin/env python3
"""
Script to take a Sam/Bam file (streamed though stdin, by default) of a bisulfite or EM-seq treated
library, and filter out 'unconverted' reads if they have 3 or more unconverted Cs in a
non-CpG context.

Unconverted reads are flagged as failing platform / vendor quality checks (0x200) and have
the tag "XX:Z:UC" added at the end.

Note: One probably doesn't want to do this for non-human organisms, as they may have
significant levels of non-CpG methylation!

210819 Ben Edits:
Want to do two things
1. Make this record converted vs unconverted reads and output a plot like I'm used to making (i.e. number of reads per amplicon, split converted vs unconverted)
2. Don't flag the reads with the weird thing just throw them out (since I'm not actually sure that our downstream code checks...)
3. OPTIONAL Fix the weird softclipping thing that is maybe a problem or maybe not I can't tell, better solution is probably just to trim the reads manually

Notes: Since we RC everything, I think "top" strand reads are such that they are actually Cs and Ts and this script should work normally. Not sure but I think.

"""

import sys
import argparse
import pysam
import re
import pandas as pd
from matplotlib import pyplot as plt

from common import parse_fasta

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required = False, help = "Reference fasta file")
    parser.add_argument("--bam", required = False, help = "Input bam or sam file (must end in .bam or .sam) [default = stdin]")
    parser.add_argument("--out", required = False, help = "Name for output sam file [default = stdout]")
    parser.add_argument("--plot", required = False, help = "Where to output the per amplicon plots")
    parser.add_argument("--c_count", required = False, default = 3, type = int, help = "Minimum number of nonconverted Cs on a read"\
                                                               " to consider it nonconverted [default = 3]")
    parser.add_argument("--c_frac", required = False, default = 0.85, type = float, help = "Fraction of Cs that must be converted to include a read")
    # parser.add_argument("--ignore_cpg", action='store_true', help="Whether CpG MTase was included (and unconverted Cs in CpG contexts shouldn't count towards the total")
    # parser.add_argument("--flag_reads", required = False, default = False, action="store_true", \
    #                     help = "Set the 'Failing platform / vendor quality check flag")

    args = parser.parse_args()
    return args


# def parse_fasta(reference):
#     """
#     Parses the reference fasta file, and reads it into a dictionary. For the header, will
#     strip anything after the first space. For example:
    
#     >chr1 AC:XXX gi:XXX LN:XXX
    
#     will be saved to the dictionary as:
#     chr1
#     """

#     fasta_dict = {}
#     fasta = open(reference, "r")
#     for line in fasta:

#         if line.startswith(">"): # Header line
#             header = line.strip().split(" ")[0][1:] # Remove '>' and anything after ' '
#             fasta_dict[header] = []

#         else: # Sequence line
#             fasta_dict[header].append(line.strip())
#     fasta.close()    

#     # If it's a multiline fasta, join the individual lines to one sequence
#     for header in fasta_dict:
#         fasta_dict[header] = "".join(fasta_dict[header])

#     return fasta_dict
	
def find_reference(header):
    """
    If a reference fasta isn't provided by the user, looks for the reference in the @PG
    lines of the bam header
    """ 
     
    for command in header['PG']:
        # Assuming bwa-meth was used for aligning
        if command['ID'] == "bwa-meth":
            assert "--reference " in command['CL'], "Couldn't find the reference fasta in "\
                    "the bam header, please provide one with the --reference option"
            
            # A regex might be better, but this will suffice for the time being
            start = command['CL'].rindex("--reference ") + len("--reference ")
            tmp = command['CL'][start:]
            end = tmp.index(" ")
            ref = tmp[:end]
            return ref   

def parse_bam(bam_file, fasta_dict, out, args):
    """
    Reads through each alignment in the bam file, if it's part of a proper pair then checks
    for non-CpG Cs. If there are 3 or more, calls filter_snps() to check if they are
    unconverted Cs, or some snp / misalignment
    """

    # Populate a dictionary to hold nonconverted read counts with 0s
    nonconverted_dict = {}

    for chrom in fasta_dict:
        nonconverted_dict[chrom] = {'converted': 0, 'total': 0, 'unconverted': 0}

    strand = {'top': 0, 'bottom': 0}
    enough_cs = {'enough': 0, 'not_enough': 0}


    # make a list of reads to write and a list of read names that are filtered
    # we need to do this to filter out halves of read pairs if their mate is bad but they pass
    to_write = []
    bad_read_names = []
    converted_fracs = []
        
    for read in bam_file:
        nonconverted_dict[read.reference_name]['total'] += 1
    
        # might want to ignore these guys too idk 
        # Ignore reads not part of a proper pair
        # Ben 210823 doesn't look like any reads are "correct" according to this, so commenting this out.
        # if not read.is_proper_pair or read.is_qcfail or read.is_duplicate or read.is_secondary \
        # or read.is_supplementary:
        #     out.write(read)
        
        if (read.is_reverse and read.is_read2) or (read.mate_is_reverse and read.is_read1): # 'Top' strand
            strand['top'] += 1
            # If there are 3 or more non-CpG Cs in the read, call filter_snps
            # if read.query_alignment_sequence.count("C") - (read.query_alignment_sequence.count("CG") + read.query_alignment_sequence.count("GC")) >= args.c_count:
            if True:
                # nonconverted_dict[read.reference_name]['unconverted'] += filter_snps(read, fasta_dict[read.reference_name], out, args)
                enough_cs['enough'] += 1
                unconverted_pred, converted_fraction = filter_snps(read, fasta_dict[read.reference_name], out, args)
                converted_fracs.append(converted_fraction)
                if unconverted_pred:
                    # don't write read
                    # increment unconverted tally
                    nonconverted_dict[read.reference_name]['unconverted'] += 1
                    # put read name into naughty list
                    bad_read_names.append(read.query_name)
                else:
                    # don't increase tally
                    # add read to writing list
                    to_write.append(read)
            else:
                enough_cs['not_enough'] += 1
                # out.write(read)
                to_write.append(read)
        
        elif (read.is_reverse and read.is_read1) or (read.mate_is_reverse and read.is_read2): # 'Bottom' strand
            strand['bottom'] += 1
            # print('EL BOTTON')
            # If there are 3 or more non-CpG Gs in the read, call filter_snps
            # if read.query_alignment_sequence.count("G") - (read.query_alignment_sequence.count("CG") + read.query_alignment_sequence.count("GC")) >= args.c_count:
            if True:
                # nonconverted_dict[read.reference_name]['unconverted'] += filter_snps(read, fasta_dict[read.reference_name], out, args)
                enough_cs['enough'] += 1
                unconverted_pred, converted_fraction = filter_snps(read, fasta_dict[read.reference_name], out, args)
                converted_fracs.append(converted_fraction)
                if unconverted_pred:
                    # don't write read
                    # increment unconverted tally
                    nonconverted_dict[read.reference_name]['unconverted'] += 1
                    # put read name into naughty list
                    bad_read_names.append(read.query_name)
                else:
                    # don't increase tally
                    # add read to writing list
                    to_write.append(read)
            else:
                enough_cs['not_enough'] += 1
                # out.write(read)
                to_write.append(read)
        else:
            pass

    # now we can do a final pass over the reads, eliminating any single guys that got their mate eliminated
    for read in to_write:
        if read.query_name in bad_read_names:
            # don't write
            # fix nonconverted_dict
            nonconverted_dict[read.reference_name]['unconverted'] += 1
        else:
            out.write(read)

    for chrom in fasta_dict:
        nonconverted_dict[chrom]['converted'] = nonconverted_dict[chrom]['total'] - nonconverted_dict[chrom]['unconverted']
        del nonconverted_dict[chrom]['total']

    print(enough_cs)
    print(strand)

    return nonconverted_dict, converted_fracs

def filter_snps(read, sequence, out, args):
    """
    For reads with 3 or more unconverted Cs, check to be sure they actually align to a C
    on the reference genome (so are not mutations, snps, misalignments, etc.)
    """

    # For all bases in read, gets genomic coordinate of alignment
    coords = read.get_reference_positions(full_length = True)
    my_seq = sequence

    # If the alignment was softclipped at all, call softclip() to remove those indices
    # from the coords list
    if coords[0] is None or coords[-1] is None:
        coords = softclip(coords) # Replacing coords with a subset of itself

    my_unconverted = 0
    my_total_c = 0
    if (read.is_reverse and read.is_read2) or (read.mate_is_reverse and read.is_read1): # 'Top' strand
        # print('top')
        for coord in coords:
            # if the reference base is a C
            if coord != None and my_seq[coord] == 'C':
                # if the C/G is not part of a GpC or CpG
                if (coord+1 < len(my_seq) and my_seq[coord+1] != 'G') and (coord > 0 and my_seq[coord-1] != 'G'): 
                    # check the read
                    if read.query_sequence[coords.index(coord)] == 'C':
                        # increment unconverted
                        my_unconverted += 1
                    # always increment total
                    my_total_c += 1
        # for coord in coords:

        #     # If the base on the read is a C and not part of an indel
        #     if read.query_sequence[coords.index(coord)] == "C" and coord != None:
                
        #         # If the C is not part of a CpG or GpC
        #         if (coords.index(coord) < len(read.query_sequence) - 1 and read.query_sequence[coords.index(coord) + 1] != "G") \
        #             and (coords.index(coord) > 0 and read.query_sequence[coords.index(coord) - 1] != "G"):
                    
        #             # If the C maps to a C on the genome (it's unconverted)
        #             if my_seq[coord] == "C":
        #                 my_unconverted += 1

        #             # then always increment total c
        #             my_total_c += 1
                
                # don't think I need this for now?
                # # If the C is the last base of the read
                # elif coords.index(coord) == len(read.query_sequence) - 1:
                    
                #     # If the C maps to a C on the genome (it's unconverted)
                #     if my_seq[coord] == "C":
                #         my_unconverted += 1

    elif (read.is_reverse and read.is_read1) or (read.mate_is_reverse and read.is_read2): # 'Bottom' strand
        # print('bot')
        for coord in coords:
            # if the reference base is a C (i.e. G)
            if coord != None and my_seq[coord] == 'G':
                # if the C/G is not part of a GpC or CpG
                if (coord+1 < len(my_seq) and my_seq[coord+1] != 'C') and (coord > 0 and my_seq[coord-1] != 'C'): 
                    # check the read
                    if read.query_sequence[coords.index(coord)] == 'G':
                        # increment unconverted
                        my_unconverted += 1
                    # always increment total
                    my_total_c += 1

        # for coord in coords:
            
        #     # If the base on the read is a G and not part of an indel
        #     if read.query_sequence[coords.index(coord)] == "G" and coord != None:

        #         # If the C is not part of a CpG or GpC
        #         if (coords.index(coord) < len(read.query_sequence) - 1 and read.query_sequence[coords.index(coord) + 1] != "C") \
        #             and (coords.index(coord) > 0 and read.query_sequence[coords.index(coord) - 1] != "C"):
                
        #         # # If the G is not part of a CpG
        #         # if coords.index(coord) > 0 and read.query_sequence[coords.index(coord) - 1] != "C":
                    
        #             # If the G maps to a G on the genome (it's unconverted)
        #             if my_seq[coord] == "G":
        #                 my_unconverted += 1

        #             # then always increment total c
        #             my_total_c += 1
                
                # # If the G is the first base of the read
                # elif coords.index(coord) == 0:
                    
                #     # If the G maps to a G on the genome (it's unconverted)
                #     if my_seq[coord] == "G":
                #         my_unconverted += 1

    # If there are c_count or more unconverted Cs on the read, set some flags and return
    # Don't write the read! Instead just filter it out. I don't think I'm smart enough to deal with the flags.
    # print(read.query_sequence)
    # print(my_seq)
    # print(my_unconverted)
    # print(my_total_c)
    converted_fraction = 1 - (my_unconverted / my_total_c)

    # if my_unconverted >= args.c_count:
    if converted_fraction < args.c_frac:
        return 1, converted_fraction
    else:
        # out.write(read)
        return 0, converted_fraction
        # read.set_tags(read.get_tags() + [("XX", "UC")])
        # if args.flag_reads is True:
        #     read.flag += 512
    #     out.write(read)
    #     return 1

    # out.write(read)
    # return 0

def softclip(coord_list):
    """
    If the alignment coords list includes softclipped bases ('None' in the list), cuts them
    from the beginning / end of the list
    """

    start = None
    end = None
    
    x = 0
    # Find first non-'None' index in the list
    while coord_list[x] is None:
        x += 1
    start = x
    
    x = len(coord_list) - 1
    # Find the last non-'None' index in the list
    while coord_list[x] is None:
        x -= 1
    end = x

    # The new coord list goes from the first - last aligned base
    new_coords = coord_list[start:end + 1]
    return new_coords

def run_filter():
    args = argparser()
    
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

    # If a reference fasta is specified, use it, otherwise look in the bam header
    if args.reference:
        reference = args.reference
    else:
        reference = find_reference(mysam.header.to_dict())
            
    fasta_dict = parse_fasta(reference)

    # If an output bam is specified, write there, otherwise write to stdout
    if args.out:
        out = pysam.AlignmentFile(args.out, "wb", template = mysam)
    else:
        out = pysam.AlignmentFile("-", "wb", template = mysam)
    
    nonconverted_counts, converted_fracs = parse_bam(mysam, fasta_dict, out, args)
    print(nonconverted_counts)


    nonconverted_counts_df = pd.DataFrame(nonconverted_counts).T

    nonconverted_counts_df = nonconverted_counts_df / 2 # to correct for 2xreads per molecule

    fig, ax = plt.subplots()

    nonconverted_counts_df.plot.bar(stacked=True, ax=ax)

    plt.tight_layout()
    plt.savefig(args.plot)
    plt.clf()

    fix, ax = plt.subplots()

    ax.hist(converted_fracs, bins=25)
    ax.axvline(args.c_frac, ls='--', c='r')

    plt.tight_layout()
    plt.savefig(args.plot.replace('.png','.fraction.png'))
    plt.clf()


    # # Write each chromosome in the reference and corresponding number of nonconverted reads
    # # to stderr
    # for chrom in mysam.references:
    #     if chrom in nonconverted_counts:
    #         sys.stderr.write("{}\t{}\n".format(chrom, nonconverted_counts[chrom]))
    #     else:
    #         sys.stderr.write("{}\t{}\n".format(chrom, 0))

if __name__ == "__main__":
    run_filter()
