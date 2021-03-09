#!/bin/bash

mkdir Spacers

for i in PLSDB Hosts Archaea Archaea_Hosts
do
    awk '{print $1}' CRISPRs/${i}/crisprs_final_pos.tab | sort | uniq > CRISPRs/${i}/acc_with_crispr
    
    for k in $(cat "CRISPRs/${i}/acc_with_crispr")
    do
        k2=$(echo $k | sed 's/\.[0-9]*$//')
        join -1 2 -2 1 <(awk '{print $2,$1"_"$3"_"$4}' CRISPRs/${i}/crisprs_final_pos.tab | sort -k2,2) \
            <(grep -h -v "^#" CRISPRs/${i}/CRISPRFinder/${k}/Here/GFF/${k2}.gff | grep CRISPRspacer \
                | cut -f 9 \
                | sed 's/sequence=//;s/;Name=.*;Parent=/\t/;s/;ID=.*//' \
                | awk '{print $2,$1}' \
                | sort -k1,1) \
            | awk '{print $2,$3}' \
            | sed 's/^/>/' \
            | awk 'sub(" ","@"NR"\n",$0)1'
    done > Spacers/${i}.fna
done

