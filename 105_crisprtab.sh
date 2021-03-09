#!/bin/bash

crispyTab () {
	join -t $'\t' -j 1 <(join -t $'\t' -1 5 -2 1 <(sort -t $'\t' -k5,5 ${1}crisprfinder.txt) \
	        <(sort -t $'\t' -k1,1 ${1}stats.tsv)) \
        	<(sed 's|,|\t|' ${1}orf_overlap.txt | sort -t $'\t' -k1,1) \
        	| awk -F'\t' '{print $1,$3,$6,$7,$8,$9,$11,$14,$15,$20,$21,$22,$23,$27,$28,$29,$30,$31,$32}' > ${1}crisprs.tab

}

crispyTab "CRISPRs/PLSDB/"
crispyTab "CRISPRs/Hosts/"
crispyTab "CRISPRs/Archaea/"
crispyTab "CRISPRs/Archaea_Hosts/"
