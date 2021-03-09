#!/bin/bash

# conda activate ccp2

doit(){
    mkdir -p Cas/${2}/logs
    NAME=$(echo $1 | sed 's|.*\/||;s|\.fna||')
    cctyper ${1} Cas/${2}/${NAME} --no_plot -t 1 --prodigal $3 &> Cas/${2}/logs/${NAME}.log
}

export -f doit

echo "Archaea"
parallel -j40 doit ::: Archaea/fastas/* ::: Archaea ::: meta
echo "Archaea Hosts"
parallel -j40 doit ::: Archaea_Hosts/fastas/* ::: Archaea_Hosts ::: single
echo "PLSDB"
parallel -j40 doit ::: PLSDB/fastas/* ::: PLSDB ::: meta
echo "PLSDB Hosts"
parallel -j40 doit ::: Hosts/fastas/* ::: Hosts ::: single



