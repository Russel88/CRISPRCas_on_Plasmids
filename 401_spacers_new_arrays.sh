#!/bin/bash

for i in PLSDB Hosts Archaea Archaea_Hosts
do
    HERE="CRISPRs/${i}/RepeatMatch"
    cut -f 1,6 ${HERE}/arrays_sub.tab | sed "s/\['//;s/'\]//;s/^/>/;s/', '/__/g" \
        | awk 'gsub("__","\n"$1"\n",$0)1' \
        | sed 's/\t/\n/' \
        | awk 'sub(">.*$",$1"@"NR/2+1/2,$1)1' > Spacers/${i}_newArrays.fna
done

