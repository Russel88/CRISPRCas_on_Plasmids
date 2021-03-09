#!/usr/bin/env python

# conda activate ccp2

import pandas as pd
from Bio import pairwise2

f = open('CRISPRs/Repeat_sim.csv', 'w')

def repeat_sim(which):
    dat = pd.read_csv('CRISPRs/{}/crisprs_final.tab'.format(which), sep=' ', header=None)
    for k in set(dat[1]):
        tmp = dat[dat[1] == k]
        seq = list(tmp[6])
        for i in range(len(seq)):
            for j in range(len(seq)):
                if i != j:
                    aln = pairwise2.align.globalxx(seq[i], seq[j], score_only=True)
                    f.write(','.join([str(list(tmp[0])[i]), 
                                      str(list(tmp[0])[j]), 
                                      str(aln/max(len(seq[i]), len(seq[j])))]))
                    f.write('\n')

for k in ("Archaea", "Archaea_Hosts", "PLSDB", "Hosts"):
    repeat_sim(k)

f.close()
