#!/usr/bin/env python

# conda activate ccp2

import numpy as np
import pandas as pd
import re

def orf_overlap(which):
    '''
    Check if CRISPRs overlap with an ORF and output a True False csv
    '''

    print(which)
    
    dat = pd.read_csv("CRISPRs/{}/crisprfinder.txt".format(which), header=None, sep="\t")

    genes = pd.read_csv("{}/genepos_sub.txt".format(which), sep=" ", header=None)

    genes['Acc'] = [x[:x.rindex('_')] for x in genes[0]]

    lst_orf = []
    for ind, row in dat.iterrows():
        print(row[1])
        tmp_gene = genes[genes['Acc'] == row[1]]
        overlap = any((tmp_gene[1] <= row[6]) & (row[5] <= tmp_gene[2]))
        lst_orf.append([row[4], overlap])
    df_orf = pd.DataFrame(lst_orf)

    df_orf.to_csv("CRISPRs/{}/orf_overlap.txt".format(which), index=False, header=None)

orf_overlap("PLSDB")
orf_overlap("Archaea")
orf_overlap("Hosts")
orf_overlap("Archaea_Hosts")

