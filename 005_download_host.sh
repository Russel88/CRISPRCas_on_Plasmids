#!/bin/bash

# conda activate ccp2

# PLSDB
awk -F'\t' '{print $2}' PLSDB/plsdb_host_etc.tsv \
	| tail -n+2 \
	| sed 's/,/\n/g' \
	| awk 'NF > 0' \
	| sort \
	| uniq > PLSDB/acc_host.txt

mkdir -p Hosts/fastas

for i in $(cat PLSDB/acc_host.txt)
do
	echo $i
	efetch -db sequences -format fasta -id $i > Hosts/fastas/${i}.fna
done

# Archaea
awk -F'\t' '{print $2}' Archaea/archaea_host_etc.tsv \
	| tail -n+2 \
	| sed 's/,/\n/g' \
	| awk 'NF > 0' \
	| sort \
	| uniq > Archaea/acc_host.txt

mkdir -p Archaea_Hosts/fastas

for i in $(cat Archaea/acc_host.txt)
do
	echo $i
	efetch -db sequences -format fasta -id $i > Archaea_Hosts/fastas/${i}.fna
done
