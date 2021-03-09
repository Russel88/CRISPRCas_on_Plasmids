#!/usr/bin/env python

# conda activate ccp2

import numpy as np
import pandas as pd
import statistics as st
from Bio import pairwise2
from Bio import SeqIO
import subprocess
import re
import numbers
import math

# Testing overlap functions
def overlap(x,y):
    return x[0] <= y[1] and y[0] <= x[1]

def overlap_any(x,ll):
    return any([overlap(x,y) for y in ll])

# Distance calculation functions
def dist(x,y):
    return y[0]-x[1] if y[0]>x[1] else x[0]-y[1]

def dist_min(x,ll):
    return min([dist(x,y) for y in ll])

def get_seq(acc, first, second, db):
    return db[acc][first:second+1]

# Functions for identity calculation between sequences
def ident(seq1, seq2):
        align = pairwise2.align.globalxx(seq1, seq2)
        score = align[0][2]/align[0][4]*100
        return score

def loopingAlign(i, j, lst):
    if j > i:
        return(ident(lst[i], lst[j]))

# Cluster repeat matches in arrays
def cluster_matches(name_up, max_dist=60, evalue=0.001):

    # Load and prepare data
    dat = pd.read_csv("CRISPRs/{}/RepeatMatch/blast.m8".format(name_up), sep="\t", header=None)

    dat = dat[dat[10] <= evalue]

    dat['start'] = [min(x,y) for x,y in zip(dat[8],dat[9])]
    dat['end'] = [max(x,y) for x,y in zip(dat[8],dat[9])]
    dat['strand'] = [1 if x < y else 2 for x,y in zip(dat[8],dat[9])]

    dat = dat.sort_values([1, 11, 'strand'], ascending=[True,False,True])

    with open(name_up + "/" + name_up.lower() + "_masked.fna", 'r') as handle:
        database = {}
        for fa in SeqIO.parse(handle, 'fasta'):
            database[fa.id] = str(fa.seq)

    # Remove overlapping matches, keep only first (highest score)
    acc = set(dat[1])
    dat_no_overlap = []
    for i in acc:
        tmp = dat[dat[1] == i]
        pos = tmp[['start','end']].values
        keep = []
        matches_all = []
        # For each match
        for ind, k in enumerate(pos):
            # Keep first match
            if len(matches_all) == 0:
                keep.append(1)
            else:
                # If overlaps with any previous, discard
                if overlap_any(k, matches_all):
                    keep.append(0)
                # No overlap with previous, keep
                else:
                    keep.append(1)
            matches_all.append(k)
        tmp.insert(len(tmp.columns), 'keep', keep)
        dat_no_overlap.append(tmp)

    dat_new = pd.concat(dat_no_overlap)
    dat_new = dat_new[dat_new['keep'] == 1].sort_values('start')

    # Cluster matches if within Xbp
    dat_cluster_list = []
    for i in acc:
        tmp = dat_new[dat_new[1] == i]
        pos = tmp[['start','end']].values
        cluster = 0
        cluster_list = []
        arrays_cluster = []
        # For each match
        for ind, k in enumerate(pos):
            # Keep first match
            if len(arrays_cluster) == 0:
                cluster_list.append(cluster)
                arrays_cluster.append(k)
            else:
                # If match within Xbp of any previous, add match to current cluster
                if dist_min(k, arrays_cluster) <= max_dist:
                    cluster_list.append(cluster)
                    arrays_cluster.append(k)
                # If match > Xbp from previous, initiate new cluster
                else:
                    cluster += 1
                    arrays_cluster = []
                    arrays_cluster.append(k)
                    cluster_list.append(cluster)
        tmp.insert(len(tmp.columns), 'cluster', cluster_list)
        dat_cluster_list.append(tmp)

    dat_cluster = pd.concat(dat_cluster_list)
    dat_cluster.insert(len(dat_cluster.columns), 'cluster_id',[x+"-"+str(y) for x,y in zip(dat_cluster[1],dat_cluster['cluster'])])
    dat_arrays = dat_cluster[[1,'start','end','cluster', 'cluster_id']]

    # Getting sequences
    clusters = set(dat_cluster['cluster_id'])

    final_list = []
    for cid in clusters:
        tmp = dat_arrays[dat_arrays['cluster_id'] == cid]
        acc = list(tmp[1])[0]
        # If only 1 match, no repeat or spacer comparisons
        if len(tmp) == 1:
            first = int(tmp['start'])
            second = int(tmp['end'])
            this = [acc, cid, first, second, get_seq(acc, first, second, database), '', 1, '', '', '']
        else:
            repeats = []
            spacers = []
            ident_spacers = ""
            ident_repeats = ""
            sem_length = ""
            mean_length = ""

            # Gather sequences
            for ind, row in tmp.reset_index().iterrows():
                # First matches, get repeat seq
                if ind == 0:
                    first = int(row['start'])
                    second = int(row['end'])
                    repeats.append(get_seq(acc, first, second, database))
                # Subsequence matches, get repeat seq, and spacer seq using end of match of prior repeat
                else:
                    repeats.append(get_seq(acc, int(row['start']), int(row['end']), database))
                    spacers.append(get_seq(acc, second+1, int(row['start'])-1, database))
                    second = int(row['end'])

            # Repeat ident
            range_repeats = range(len(repeats))
            matches_repeats = [loopingAlign(k,l,repeats) for k in range_repeats for l in range_repeats]
            if len(matches_repeats) > 1:
                ident_repeats = st.mean([x for x in matches_repeats if isinstance(x, numbers.Number)])

            # Spacer ident
            range_spacers = range(len(spacers))
            matches_spacers = [loopingAlign(k,l,spacers) for k in range_spacers for l in range_spacers]
            if len(matches_spacers) > 1:
                ident_spacers = st.mean([x for x in matches_spacers if isinstance(x, numbers.Number)])
                sem_length = st.stdev([len(x) for x in spacers])/math.sqrt(len(spacers))
            if len(matches_spacers) >= 1:
                mean_length = st.mean([len(x) for x in spacers])
            this = [acc, cid, first, second, repeats, spacers, len(tmp), ident_repeats, ident_spacers, sem_length, mean_length]
        final_list.append(this)

    # Save
    final = pd.DataFrame(final_list, columns=['Acc','Array_id','Start','End','Repeats','Spacers','N_repeats','Repeat_Cons','Spacer_Cons',"Spacer_lengthSEM", "Spacer_Length"])
    final.to_csv("CRISPRs/{}/RepeatMatch/arrays.tab".format(name_up), sep="\t", index=False)

# Run
print("PLSDB")
cluster_matches("PLSDB")
print("Hosts")
cluster_matches("Hosts")
print("Archaea")
cluster_matches("Archaea")
print("Archaea_Hosts")
cluster_matches("Archaea_Hosts")
