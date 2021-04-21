#!/bin/bash

sed 's/\(.*\)\(_[0-9]*\)\b/\1 \1\2/' PLSDB/genepos_sub.txt | sort -k1,1 > PLSDB/genepos_sub2.txt
sed 's/\(.*\)\(_[0-9]*\)\b/\1 \1\2/' Archaea/genepos_sub.txt | sort -k1,1 > Archaea/genepos_sub2.txt

for THIS in Archaea PLSDB Archaea_Hosts Hosts
do

    echo $THIS

    # Prepare comparison
    awk 'function min(v,w) {{if( v < w ) { MIN=v } if( w <= v ) { MIN=w } }{ return MIN }}; function max(v,w) {{if( v < w ) { MAX=w } if( w <= v ) { MAX=v } }{ return MAX }}{print $2,$1,min($9,$10),max($9,$10)}' Spacers/${THIS}_plsdb.m8 | sort -k1,1 > Spacers/${THIS}_plsdb2.m8

    if [[ $THIS =~ "Archaea" ]] 
    then
        THAT="Archaea"
    else
        THAT="PLSDB"
    fi

    echo $THAT

    # Only ORF matching
    join -j1 Spacers/${THIS}_plsdb2.m8 ${THAT}/genepos_sub2.txt | awk 'function ol(s1,e1,s2,e2) {{if( s1 <= e2 && e1 >= s2) {OL="TRUE" } else {OL="FALSE"} }{ return OL }}{print $0,ol($3,$4,$6,$7)}' | awk '$8 == "TRUE"' | awk '{print $2,$1,$3,$4,$5}'| sort | uniq > Spacers/${THIS}_plsdb_orf.m8

done
