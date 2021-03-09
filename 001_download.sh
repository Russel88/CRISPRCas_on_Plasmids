#!/bin/bash

# conda activate ccp2 

# PLSDB metadata
wget https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip
mv index.html?zip plsdb.zip
unzip plsdb.zip
mkdir PLSDB
mv plsdb* PLSDB
mv README.md PLSDB

# Extract Acc
awk -F'\t' '{print $2}' PLSDB/plsdb.tsv | tail -n+2 > PLSDB/acc_plasmid.txt

# Download GenBank and Fastas
mkdir PLSDB/gbs
mkdir PLSDB/fastas
mkdir PLSDB/assembly

for i in $(cat PLSDB/acc_plasmid.txt)
do
	echo $i
	efetch -db sequences -format gb -id $i > PLSDB/gbs/${i}.gb
	efetch -db sequences -format fasta -id $i > PLSDB/fastas/${i}.fna
done

# Assembly information
wget ftp://ftp.ncbi.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt
wget ftp://ftp.ncbi.nih.gov/genomes/genbank/assembly_summary_genbank.txt
cat <(tail -n+2 assembly_summary_genbank.txt | sed '1 s/# //') <(grep -v '^#' assembly_summary_genbank_historical.txt) > PLSDB/assembly_summary_genbank_all.txt
rm assembly_summary_genbank*

# Download Archaea
mkdir -p Archaea/gbs
mkdir -p Archaea/fastas
mkdir -p Archaea/assembly

ln -s ${PWD}/PLSDB/assembly_summary_genbank_all.txt Archaea/assembly_summary_genbank_all.txt

for i in $(cat Archaea/acc_plasmid.txt)
do
        echo $i
        efetch -db sequences -format gb -id $i > Archaea/gbs/${i}.gb
        efetch -db sequences -format fasta -id $i > Archaea/fastas/${i}.fna
done

