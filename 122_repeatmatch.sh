#!/bin/bash

# conda activate ccp2

matching(){
mkdir -p CRISPRs/${1}/RepeatMatch
echo $1
blastn -query CRISPRs/Repeats/nr_${3}.fasta -db ${1}/${2}_masked -task blastn-short \
	-perc_identity 95 -qcov_hsp_perc 95 \
	-dust no -soft_masking false -num_threads 8 \
	-outfmt=6 > CRISPRs/${1}/RepeatMatch/blast.m8
}

matching "PLSDB" "plsdb" "plsdb"
matching "Hosts" "hosts" "plsdb"
matching "Archaea" "archaea" "arch"
matching "Archaea_Hosts" "archaea_hosts" "arch"

