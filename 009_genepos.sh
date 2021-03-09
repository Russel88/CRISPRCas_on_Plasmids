#!/bin/bash

genepos(){
    echo "${1}"
    for i in ${1}/prot/*.faa; do grep ">" $i | cut -d"#" -f1,2,3 | sed 's|>||;s|# ||g'; done > ${1}/genepos.txt

    # Good quality ORFs
    for i in ${1}/gff/*.gff; do grep -v "^#" $i | cut -f1,9 | cut -d";" -f1,7 | sed 's/;conf=/\t/;s/\tID=1//'; done | awk '$2 >= 90' | cut -f1 > ${1}/good_genes.tab

    join -j 1 <(sort -k1,1 ${1}/genepos.txt) <(sort ${1}/good_genes.tab) > ${1}/genepos_sub.txt
}

genepos "PLSDB"
genepos "Hosts"
genepos "Archaea"
genepos "Archaea_Hosts"


