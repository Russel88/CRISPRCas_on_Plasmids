#!/usr/bin/env python

# conda activate ccp2

from crispyStats import crispy_stats, array_stats, ident, loopingAlign
import pandas as pd

# PLSDB
f = open("CRISPRs/PLSDB/acc.txt", "r")
these_plsdb = f.read().split("\n")[:-1]

res_lst = []
for accession in these_plsdb:
    print(accession)
    res_lst = res_lst + crispy_stats(accession, "PLSDB")

res_df = pd.DataFrame(res_lst)

res_df.to_csv("CRISPRs/PLSDB/stats.tsv", sep="\t", header=False, index=False)

# Archaea
f = open("CRISPRs/Archaea/acc.txt", "r")
these_arch = f.read().split("\n")[:-1]

res_lst = []
for accession in these_arch:
    print(accession)
    res_lst = res_lst + crispy_stats(accession, "Archaea")

res_df = pd.DataFrame(res_lst)

res_df.to_csv("CRISPRs/Archaea/stats.tsv", sep="\t", header=False, index=False)

# PLSDB Hosts
f = open("CRISPRs/Hosts/acc.txt", "r")
these_plsdb = f.read().split("\n")[:-1]

res_lst = []
for accession in these_plsdb:
    print(accession)
    res_lst = res_lst + crispy_stats(accession, "Hosts")

res_df = pd.DataFrame(res_lst)

res_df.to_csv("CRISPRs/Hosts/stats.tsv", sep="\t", header=False, index=False)

# Archaea Hosts
f = open("CRISPRs/Archaea_Hosts/acc.txt", "r")
these_arch = f.read().split("\n")[:-1]

res_lst = []
for accession in these_arch:
    print(accession)
    res_lst = res_lst + crispy_stats(accession, "Archaea_Hosts")

res_df = pd.DataFrame(res_lst)

res_df.to_csv("CRISPRs/Archaea_Hosts/stats.tsv", sep="\t", header=False, index=False)

