import sys
import pandas as pd
from Bio import SeqIO


amplicon_fa = sys.argv[1]
peaklist_bed = sys.argv[2]

with open(peaklist_bed, 'w') as out_bed:
    for r in list(SeqIO.parse(amplicon_fa, "fasta")):
        line_to_write = '{}\t0\t{}\t+\n'.format(r.id, len(r.seq))
        out_bed.write(line_to_write)
