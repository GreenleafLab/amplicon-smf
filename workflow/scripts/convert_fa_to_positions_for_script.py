## BD 230328
## this script takes as input an amplicon fa and a list of motifs
## and outputs a positions file for the binding model

import re
from Bio import SeqIO, Seq
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert amplicon.fa and motif list to positions for script file')
    parser.add_argument("--input_fa", dest="input_fa", type=str, help="Path to amplicon fa")
    parser.add_argument("--input_motifs", dest="input_motifs", type=str, help="Path to input motifs file")
    parser.add_argument("--output_positions", dest="output_positions", type=str, help="Path for writing the output positions file to")
    parser.add_argument("--l_offset", dest="l_offset", type=int, default=2, help="How many positions to extend the motif to the left")
    parser.add_argument("--r_offset", dest="r_offset", type=int, default=2, help="How many positions to extend the motif to the right")
    
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

        

