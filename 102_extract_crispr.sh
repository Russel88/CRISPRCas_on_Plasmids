#!/bin/bash

# PLSDB
for i in CRISPRs/PLSDB/CRISPRFinder/*
do
	tail -n+2 ${i}/Here/TSV/Crisprs_REPORT.tsv
done | awk 'NF>0' > CRISPRs/PLSDB/crisprfinder.txt

tail -n+2 CRISPRs/PLSDB/crisprfinder.txt | awk -F'\t' '{print $2}' | sort | uniq > CRISPRs/PLSDB/acc.txt

# Archaea
for i in CRISPRs/Archaea/CRISPRFinder/*
do
	tail -n+2 ${i}/Here/TSV/Crisprs_REPORT.tsv
done | awk 'NF>0' > CRISPRs/Archaea/crisprfinder.txt

tail -n+2 CRISPRs/Archaea/crisprfinder.txt | awk -F'\t' '{print $2}' | sort | uniq > CRISPRs/Archaea/acc.txt

# PLSDB Hosts
for i in CRISPRs/Hosts/CRISPRFinder/*
do
        tail -n+2 ${i}/Here/TSV/Crisprs_REPORT.tsv
done | awk 'NF>0' > CRISPRs/Hosts/crisprfinder.txt

tail -n+2 CRISPRs/Hosts/crisprfinder.txt | awk -F'\t' '{print $2}' | sort | uniq > CRISPRs/Hosts/acc.txt

# Archaea Hosts
for i in CRISPRs/Archaea_Hosts/CRISPRFinder/*
do
        tail -n+2 ${i}/Here/TSV/Crisprs_REPORT.tsv
done | awk 'NF>0' > CRISPRs/Archaea_Hosts/crisprfinder.txt

tail -n+2 CRISPRs/Archaea_Hosts/crisprfinder.txt | awk -F'\t' '{print $2}' | sort | uniq > CRISPRs/Archaea_Hosts/acc.txt

