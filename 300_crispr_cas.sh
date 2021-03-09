#!/bin/bash

mkdir -p CRISPRCas

crisprcas(){

cut -f 1,2,3,4 Cas/${1}_prediction.tab | sort -k1,1 > Cas/${1}_pos.tab

awk '{print $2,$1,$3,$4}' CRISPRs/${1}/crisprs_final.tab > CRISPRs/${1}/crisprs_final_pos.tab
awk '{print $2,$1,$3,$4}' CRISPRs/${1}/RepeatMatch/arrays_sub.tab > CRISPRs/${1}/RepeatMatch/crisprs_final_pos.tab

cat CRISPRs/${1}/crisprs_final_pos.tab CRISPRs/${1}/RepeatMatch/crisprs_final_pos.tab | sort -k1,1 > CRISPRs/${1}_crisprs_all_pos.tab

# Distances
join -j1 <(join -j 1 CRISPRs/${1}_crisprs_all_pos.tab Cas/${1}_pos.tab) <(sort -k1,1 ${1}/size2.tab) \
	| awk '{print $1,$2,$5,$3-$7,$6-$4,$8-$4+$6,$8-$7+$3}' \
	| awk 'function abs(v) {return v < 0 ? -v : v} \
		function min(v,w) {{if( v < w ) { MIN=v } if( w <= v ) { MIN=w } }{ return MIN }} \
		($4 >= 0 || $5 >= 0) {print $1,$2,$3,min(min(abs($4),abs($5)),min(abs($6),abs($7)))} \
		($4 < 0 && $5 < 0) {print $1,$2,$3,0}' > CRISPRCas/${1}_crispr_cas_distance.tab

## Pick nearest
sort -k2,2 -k4,4n CRISPRCas/${1}_crispr_cas_distance.tab | sort -k2,2 -u > CRISPRCas/${1}_crispr_cas_distance_nearest.tab

}

crisprcas "PLSDB"
crisprcas "Archaea"
crisprcas "Hosts"
crisprcas "Archaea_Hosts"
