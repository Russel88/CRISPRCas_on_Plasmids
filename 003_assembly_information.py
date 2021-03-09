#!/usr/bin/env python

# conda activate ccp2

import pandas as pd
import numpy as np
import wget

def download(which):
    '''
    Dównload assembly reports
    '''

    # Load table
    df_merge = pd.read_csv("{}/{}_assembly_metadata.tsv".format(which, which.lower()), sep="\t", header=0)

    # Download assembly reports
    print("Downloading assembly reports")
    these = list(range(0, len(df_merge), int(len(df_merge)/100)))
    list_assem = []
    for ind, row in df_merge.iterrows():

        if ind in these:
            print(str(round(ind / len(df_merge) * 100, 1))+"%", end="\r")

        acc = row[0]

        ## Download
        url = row['ftp_path'] + "/" + row['ftp_path'].split("/")[-1] + "_assembly_report.txt"

        # Try Refseq assembly first, if not possible change to GenBank assembly
        try:
            try:
                wget.download(url, out = "{}/assembly/{}".format(which, acc))
            except:
                url2 = url.replace("GCF","GCA")
                wget.download(url2, out = "{}/assembly/{}".format(which, acc))
        except:
            print(acc + " failed")

download("PLSDB")
download("Archaea")
