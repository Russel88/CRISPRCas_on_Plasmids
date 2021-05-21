#!/bin/bash

# conda activate mobsuite

# PLSDB
mkdir -p PLSDB/mob

doit(){
	nice -n10 mob_typer --infile PLSDB/fastas/${1}.fna --out_file PLSDB/mob/${1}
}

export -f doit

cat PLSDB/acc_plasmid.txt | parallel -j20 doit

# Archaea
mkdir -p Archaea/mob

doit(){
	nice -n10 mob_typer --infile Archaea/fastas/${1}.fna --out_file Archaea/mob/${1}
}

export -f doit

cat Archaea/acc_plasmid.txt | parallel -j20 doit
