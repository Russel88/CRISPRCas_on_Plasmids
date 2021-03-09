#!/usr/bin/env python

# conda activate ccp2

import numpy as np
import pandas as pd
import re

def orfoverlap(dd, ddg):
    '''
    Check if new arrays overlap with an ORF
    dd is array table
    ddg is table with gene positions
    '''

    lst_orf = []
    for ind, row in dd.iterrows():
        tmp_gene = ddg[ddg['Acc'] == row['Acc']]
        overlap = any((tmp_gene[1] <= row['End']) & (row['Start'] <= tmp_gene[2]))
        lst_orf.append([row['Array_id'], overlap])

    return pd.DataFrame(lst_orf)

def run_oo(which):
    '''
    Run the ORF overlap function a dataset
    '''

    dat = pd.read_csv("CRISPRs/{}/RepeatMatch/arrays.tab".format(which), sep="\t")
    genes = pd.read_csv("{}/genepos_sub.txt".format(which), sep=" ", header=None)
    genes['Acc'] = [x[:x.rindex('_')] for x in genes[0]]
    df_orf = orfoverlap(dat, genes)
    df_orf.to_csv("CRISPRs/{}/RepeatMatch/orf_overlap.txt".format(which), index=False, header=None)

run_oo("PLSDB")
run_oo("Hosts")
run_oo("Archaea")
run_oo("Archaea_Hosts")


