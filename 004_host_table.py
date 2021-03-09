#!/usr/bin/env python

# conda activate ccp2

import pandas as pd
import numpy as np
import os

def hosts(which):
    '''
    Get host information from the assembly reports 
    '''

    # For each assembly report
    list_assem = list()
    for file in os.listdir('{}/assembly'.format(which)):

        ## Initialize list
        dat = list()
        dat.append(file)

        df_report = pd.read_csv('{}/assembly/'.format(which) + file, comment='#', sep='\t', header=None)
        arr_assem = df_report.loc[:, (3,4,6)].values

        ## Split in roles
        chrom = [list(x) for x in arr_assem if x[0] == "Chromosome"]
        plas = [list(x) for x in arr_assem if x[0] == "Plasmid"]
        rest = [list(x) for x in arr_assem if x[0] not in ["Plasmid","Chromosome"]]

        ## Refseq acc if possisble, else genbank acc
        chrom_acc = [x[2] if x[2] != "na" else x[1] for x in chrom]
        plas_acc = [x[2] if x[2] != "na" else x[1] for x in plas]
        rest_acc = [x[2] if x[2] != "na" else x[1] for x in rest]

        ## Add to list
        dat.append(chrom_acc)
        dat.append(plas_acc)
        dat.append(rest_acc)

        list_assem.append(dat)

    # To DataFrame and write
    df_assem = pd.DataFrame(list_assem, columns = ("Acc", "Chromosomes", "Plasmids", "Unplaced"))
    df_assem['Chromosomes'] = [','.join(x) for x in df_assem['Chromosomes']]
    df_assem['Plasmids'] = [','.join(x) for x in df_assem['Plasmids']]
    df_assem['Unplaced'] = [','.join(x) for x in df_assem['Unplaced']]
    df_assem.to_csv("{}/{}_host_etc.tsv".format(which, which.lower()), sep="\t", index=False)

hosts("PLSDB")
hosts("Archaea")
