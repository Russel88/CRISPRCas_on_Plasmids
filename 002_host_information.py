#!/usr/bin/env python

# conda activate ccp2

from Bio import SeqIO
import pandas as pd
import numpy as np

def assem_info(which):
    '''
    Extract assembly and sample information for each accession
    '''

    # Get accessions
    acc = []
    f = open("{}/acc_plasmid.txt".format(which), "r")
    for i in f.readlines():
        acc.append(i.rstrip('\n'))

    # Extract database cross-references
    print("Extracting database info")
    dicts = []
    for i in acc:
        for record in SeqIO.parse("{}/gbs/{}.gb".format(which, i), "genbank"):
            subdict = {"Acc": record.id}
            sub = [x.split(":") for x in record.dbxrefs]
            if len(sub) > 0:
                for j in sub:
                    if len(j) == 2:
                        subdict[j[0]] = j[1]
        dicts.append(subdict)
    df = pd.DataFrame(dicts)

    # Load assembly data info
    assem = pd.read_csv("{}/assembly_summary_genbank_all.txt".format(which), sep="\t", header=0, dtype={'relation_to_type_material': str})

    # Merge tables with available Assembly or BioSample info
    ## Merge if Assembly info matches
    df_A = df[df["Assembly"].notna()]
    assem_A = assem[assem["gbrs_paired_asm"].notna()]
    df_merge_A = pd.merge(df_A, assem_A, left_on='Assembly', right_on='gbrs_paired_asm').drop_duplicates(['Acc'])

    ## Merge remaining if BioSample info matches
    df_B = df[df["Assembly"].isna()]
    df_B = df_B[df_B["BioSample"].notna()]
    assem_B = assem[assem["biosample"].notna()]
    df_merge_B = pd.merge(df_B, assem_B, left_on='BioSample', right_on='biosample').drop_duplicates(['Acc'])

    ## Put together
    df_merge = pd.concat([df_merge_A, df_merge_B])

    ## When Refseq acc exists change path to assembly report
    m = df_merge['gbrs_paired_asm'] != "na"
    df_merge.loc[m, 'ftp_path'] = df_merge.loc[m, 'ftp_path'].str.replace("GCA", "GCF")

    ## Write
    df_merge.drop(columns=['biosample', 'bioproject', 'Assembly']).to_csv("{}/{}_assembly_metadata.tsv".format(which, which.lower()), sep="\t", index=False)

assem_info("PLSDB")
assem_info("Archaea")
