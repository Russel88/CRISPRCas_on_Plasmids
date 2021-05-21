#!/bin/bash

# conda activate ccp2

for i in PLSDB Hosts Archaea Archaea_Hosts
do
    # Combine
    cat Spacers/${i}.fna Spacers/${i}_newArrays.fna > Spacers/${i}_all.fna 
done

