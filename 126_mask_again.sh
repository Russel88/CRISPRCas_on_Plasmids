#!/bin/bash

# conda activate ccp2

mask(){

    # Masking file
    cut -f 2,3,4 CRISPRs/${1}/RepeatMatch/arrays_sub.tab | sed 's/\t/,/g' > CRISPRs/${1}/RepeatMatch/mask.csv

    # Run masking
    echo "Masking ${1}"
    ./mask.py ${1}/${2}_masked.fna CRISPRs/${1}/RepeatMatch/mask.csv ${1}/${2}_masked2.fna

}

mask "PLSDB" "plsdb"
mask "Hosts" "hosts"
mask "Archaea" "archaea"
mask "Archaea_Hosts" "archaea_hosts"


