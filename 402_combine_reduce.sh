#!/bin/bash

# conda activate ccp2

for i in PLSDB Hosts Archaea Archaea_Hosts
do
    # Combine
    cat Spacers/${i}.fna Spacers/${i}_newArrays.fna > Spacers/${i}_all.fna 
    
    # Reduce
    cd-hit-est -d 0 -c 1 -s 1 -i Spacers/${i}_all.fna -o Spacers/${i}_nr.fna

    # Table
    ./cdhit_cluster_parse.py Spacers/${i}_nr.fna.clstr > Spacers/${i}_clstr.tab
done

