#!/bin/bash

for i in PLSDB Hosts Archaea Archaea_Hosts
do
    HERE="CRISPRs/${i}/RepeatMatch"
    # Merge
    join -t$'\t' -1 2 -2 1 <(sort -k2,2 ${HERE}/arrays.tab) <(sed 's/,/\t/' ${HERE}/orf_overlap.txt | sort -k1,1) > ${HERE}/arrays_merge.tab
    # Subset
    awk -F'\t' '$8 >= 70 && $9 < 50 && $10 < 3 && $12 == "False"' ${HERE}/arrays_merge.tab > ${HERE}/arrays_sub.tab
done




