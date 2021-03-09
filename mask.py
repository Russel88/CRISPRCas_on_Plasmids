#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[3]
xfile = sys.argv[2]

# Read Xtable
Xtab = pd.read_csv(xfile, sep=",", header=None)

with open(outfile, 'w') as out_file:
        falist = SeqIO.parse(open(infile, 'r'), 'fasta')
        # For each sequence
        for fas in falist:
                name = fas.id
                Xsub = Xtab[Xtab[0].str.contains(name)]
                # Only fastas found in Xtable
                if not Xsub.empty:
                        seq = str(fas.seq)
                        # Each row of Xtable
                        for row in Xsub.itertuples(index=True):
                                # From where to where
                                Xfrom = min(int(getattr(row, '_2')), int(getattr(row, '_3')))
                                Xto = max(int(getattr(row, '_2')), int(getattr(row, '_3')))
                                # New sequence
                                seq1 = seq[:int(Xfrom)-1]
                                seqX = 'N'*(Xto-Xfrom+1)
                                seq2 = seq[int(Xto):]
                                seq = seq1+seqX+seq2
                        fas.seq = Seq(seq)
                # Write sequence
                SeqIO.write(fas, out_file, "fasta")

