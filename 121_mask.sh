#!/bin/bash

# conda activate ccp2

mask(){

    echo "Concat ${1}"
    cat ${1}/fastas/*.fna > ${1}/${2}.fna

    # Run masking
    echo "Masking ${1}"
    ./mask.py ${1}/${2}.fna CRISPRs/${1}/crisprs_mask.csv ${1}/${2}_masked.fna

    # Clean up
    echo "Cleaning up"
    rm ${1}/${2}.fna

    # Blast database
    echo "Creating databases"
    makeblastdb -dbtype nucl -in ${1}/${2}_masked.fna -out ${1}/${2}_masked

}

mask "PLSDB" "plsdb"
mask "Hosts" "hosts"
mask "Archaea" "archaea"
mask "Archaea_Hosts" "archaea_hosts"

