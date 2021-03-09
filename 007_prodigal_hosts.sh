#!/bin/bash

# conda activate ccp2

# PLSDB
mkdir -p Hosts/prot
mkdir -p Hosts/gff

doit(){
	nice -n10 prodigal -i Hosts/fastas/${1}.fna -a Hosts/prot/${1}.faa -f gff -o Hosts/gff/${1}.gff
}

export -f doit

cat PLSDB/acc_host.txt | parallel -j20 doit

# Archaea
mkdir -p Archaea_Hosts/prot
mkdir -p Archaea_Hosts/gff

doit(){
        nice -n10 prodigal -i Archaea_Hosts/fastas/${1}.fna -a Archaea_Hosts/prot/${1}.faa -f gff -o Archaea_Hosts/gff/${1}.gff
}

export -f doit

cat Archaea/acc_host.txt | parallel -j20 doit

