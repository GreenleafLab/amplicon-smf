## BD 230328
## this script takes as input an amplicon fa and a list of motifs
## and outputs a positions file for the binding model
##
## COORDINATE / OFF-BY-ONE WARNING (verified 2026-07-13, see
## workflow/scripts/NOTES_v5_hsmm_status.md and the parent CLAUDE.md GCG/coordinate notes):
##  - The amplicon is reverse-complemented here (matching the pipeline's bottom_strand=TRUE), so
##    the emitted coordinates are in the RC / bottom-strand frame -- the SAME frame as the
##    single-molecule matrix columns (verified offset 0 vs RC-reference GC dinucleotides).
##  - Each motif window is written as HALF-OPEN [start, end): the base at `end` is NOT included.
##    With `end = match.start()+len(motif)+r_offset`, the right edge is right-exclusive, so it is
##    easy to drop the motif's right-flanking GpC by 1 bp. The production opJS45 positions.long.txt
##    did exactly this (left flank C kept, right flank C at the excluded `end` -> dropped). It was
##    harmless only because the v4/v5 classifier adds tf_margin=2 at runtime, which re-includes it.
##  - The offsets therefore DEFAULT TO 0 (changed 2026-08-03; they used to default to 2/2): the file
##    holds the BARE motif span and the classifier's tf_margin does the flank extension
##    symmetrically at runtime. Then ALWAYS verify that each flanking GpC readout C lands in the
##    runtime window [start-tf_margin, end+tf_margin) by matching against the RC amplicon + the
##    matrix GpC columns before trusting a new file.
##  - BACKWARD COMPATIBILITY: positions files generated before this change were built with the old
##    2/2 defaults. Regenerating one now will NOT reproduce it byte-for-byte -- pass
##    --l_offset 2 --r_offset 2 explicitly if you specifically need the old (double-extended)
##    behavior. The v5 HSMM classifier expects the bare-span (0/0) convention.

import re
from Bio import SeqIO, Seq
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert amplicon.fa and motif list to positions for script file')
    parser.add_argument("--input_fa", dest="input_fa", type=str, help="Path to amplicon fa")
    parser.add_argument("--input_motifs", dest="input_motifs", type=str, help="Path to input motifs file")
    parser.add_argument("--output_positions", dest="output_positions", type=str, help="Path for writing the output positions file to")
    # NOTE: default 0/0 -- let the classifier's tf_margin extend the flanks (see header warning).
    # These offsets are baked into the FILE and then ADDED to again by tf_margin at runtime, so a
    # nonzero value here double-extends the motif window.
    parser.add_argument("--l_offset", dest="l_offset", type=int, default=0, help="Extend motif to the left (bp). Default 0 (bare motif span); see header off-by-one warning.")
    parser.add_argument("--r_offset", dest="r_offset", type=int, default=0, help="Extend motif to the right (bp). Default 0 (bare motif span); right edge is half-open/exclusive.")
    
    args = parser.parse_args()

    # grab motifs
    motif_dict = {}

    for record in list(SeqIO.parse(args.input_motifs, "fasta")):
        motif_dict[record.id] = str(record.seq).upper()

    # iterate through amplicon.fa IDing motifs
    with open(args.output_positions, 'w') as output_positions:
        for record in list(SeqIO.parse(args.input_fa, "fasta")):
            output_positions.write('>' + record.id + '\n')
            # for now, since we are RCing the fa, we also need to RC the seqs here, but this will change
            amplicon_seq = str(record.seq.reverse_complement()).upper()
            # iterate through motifs
            for motif in motif_dict.keys():
                motif_seq = motif_dict[motif]
                for match in re.finditer(motif_seq, amplicon_seq):
                    output_positions.write('{},{},{},r\n'.format(match.start()-args.l_offset, match.start()+len(motif_seq)+args.r_offset, motif))
                # also try the reverse complement
                motif_rc = str(Seq.Seq(motif_seq).reverse_complement())
                for match in re.finditer(motif_rc, amplicon_seq):
                    output_positions.write('{},{},{},r\n'.format(match.start()-args.l_offset, match.start()+len(motif_seq)+args.r_offset, motif))

        

