#!/bin/bash

# conda activate ccp2

# PLSDB
mkdir -p PLSDB/prot
mkdir -p PLSDB/gff

doit(){
	nice -n10 prodigal -i PLSDB/fastas/${1}.fna -a PLSDB/prot/${1}.faa -p meta -f gff -o PLSDB/gff/${1}.gff
}

export -f doit

cat PLSDB/acc_plasmid.txt | parallel -j20 doit

# Archaea
mkdir -p Archaea/prot
mkdir -p Archaea/gff

doit(){
        nice -n10 prodigal -i Archaea/fastas/${1}.fna -a Archaea/prot/${1}.faa -p meta -f gff -o Archaea/gff/${1}.gff
}

export -f doit

cat Archaea/acc_plasmid.txt | parallel -j20 doit
