#!/usr/bin/env python

# Functions
## Align and return identity
def ident(seq1, seq2):
    '''
    Align two sequences and output identity
    '''

    from Bio import pairwise2
    align = pairwise2.align.globalxx(seq1, seq2)
    score = align[0][2]/align[0][4]*100
    return score

## Loop all pairwise, except self
def loopingAlign(i, j, lst):
    '''
    Get identities of all pairs of sequences from a list
    '''
    
    if j > i:
        return(ident(lst[i], lst[j]))

## Get stats for an array
def array_stats(crisp_acc, crisprs, repeats, spacers):
    '''
    Calculate stats on CRISPR arrays
    '''

    from joblib import Parallel, delayed
    import re
    import numbers
    import statistics as st

    parent = 'Parent=' + crisp_acc.split('=')[1]
    name = re.sub('_Crispr', '', list(crisprs[crisprs['name'] == crisp_acc]['id'])[0].split('=')[1])

    ### Repeats
    repeats_sub = [x.split('=')[1] for x in repeats[repeats['parent'] == parent]['sequence']]
    range_repeats = range(len(repeats_sub))
    matches_repeats = Parallel(n_jobs=10)(delayed(loopingAlign)(k, l, repeats_sub) for k in range_repeats for l in range_repeats)
    if len(matches_repeats) > 1:
        ident_repeats = st.mean([x for x in matches_repeats if isinstance(x, numbers.Number)])
    else:
        ident_repeats = "NA"

    ### Spacers
    spacers_sub = [x.split('=')[1] for x in spacers[spacers['parent'] == parent]['sequence']]
    range_spacers = range(len(spacers_sub))
    matches_spacers = Parallel(n_jobs=10)(delayed(loopingAlign)(k, l, spacers_sub) for k in range_spacers for l in range_spacers)
    mean_length = st.mean([len(x) for x in spacers_sub])
    if len(matches_spacers) > 1:
        ident_spacers = st.mean([x for x in matches_spacers if isinstance(x, numbers.Number)])
        sd_length = st.stdev([len(x) for x in spacers_sub])
    else:
        ident_spacers = "NA"
        sd_length = "NA"

    return [name,ident_repeats,ident_spacers,mean_length,sd_length]

## Read gff and output stats
def crispy_stats(acc, which):
    '''
    Run the array_stats on all CRISPRs
    '''

    import pandas as pd
    import re
    acc_noVersion = re.sub(r'\..*', '', acc)
    dat = pd.read_csv("CRISPRs/{}/CRISPRFinder/{}/Here/GFF/{}.gff".format(which, acc, acc_noVersion),
                      sep='\t', comment='#', header=None)
    crisprs = dat[dat[2] == "CRISPR"].copy()
    repeats = dat[dat[2] == "CRISPRdr"].copy()
    spacers = dat[dat[2] == "CRISPRspacer"].copy()
    repeats[['sequence','parent','id']] = repeats[8].str.split(';', expand=True)
    spacers[['sequence','name','parent','id']] = spacers[8].str.split(';', expand=True)
    crisprs[['dr','dr_length','nspacers','name','id','potential_direction','etc']] = crisprs[8].str.split(';', expand=True)
    crispy_lst = []

    ### For each array
    for i in crisprs['name']:
        crispy_lst.append(array_stats(i, crisprs, repeats, spacers))
    return crispy_lst

