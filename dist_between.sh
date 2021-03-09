#!/bin/bash

# Find the distance between pairs of genomic elements on contigs. E.g. CRISPR arrays and Cas operons
# Two inputs are needed each containing the positions of the elements. E.g. first could contain position of CRISPR arrays, the other positions of Cas operons.
# Both inputs should be whitespace-delimited with four columns: Contig, Element name, start, end. Start has to be lower than end.
#
#
# Example:
# ./dist_between.sh File1 File2
# 
# File1:
# Contig1 CRISPR_1 10 1200
# Contig1 CRISPR_2 7000 9000
# Contig2 CRISPR_1 120 2310
# Contig2 CRISPR_2 18500 19000
#
# File2:
# Contig1 Cas_1 1400 4000
# Contig2 Cas_1 12000 15000
# Contig2 Cas_2 18000 20000
#
# Output:
# Contig1 CRISPR_1 Cas_1 200
# Contig1 CRISPR_2 Cas_1 3000
# Contig2 CRISPR_1 Cas_1 9690
# Contig2 CRISPR_1 Cas_2 15690
# Contig2 CRISPR_2 Cas_1 3500
# Contig2 CRISPR_2 Cas_2 0

join -j 1 <(sort -k1,1 ${1}) <(sort -k1,1 ${2}) \
    | awk '{print $1,$2,$5,$3-$7,$6-$4}' \
    | awk 'function abs(v) {return v < 0 ? -v : v} \
                            function min(v,w) {{if( v < w ) { MIN=v } if( w < v ) { MIN=w } }{ return MIN }} \
                            ($4 > 0 || $5 > 0) {print $1,$2,$3,min(abs($4),abs($5))} \
                            ($4 < 0 && $5 < 0) {print $1,$2,$3,0}'

