#!/bin/bash

# Find the nearest distance between pairs of genomic elements on contigs. E.g. CRISPR arrays and Cas operons
# Two inputs are needed each containing the positions of the elements. E.g. first could contain position of CRISPR arrays, the other positions of Cas operons.
# Both inputs should be whitespace-delimited with four columns: Contig, Element name, start, end. Start has to be lower than end.
#
# Unique elements in file 1 is kept. For example, if file1 contains CRISPR arrays and file2 contains Cas operons,
# then the output contains the nearest Cas operon to each CRISPR array. Cas operons might therefore be duplicated.
#
# Example:
# ./dist_nearest.sh File1 File2
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
# Contig2 CRISPR_2 Cas_2 0

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

${DIR}/dist_between.sh ${1} ${2} | sort -k1,1 -k2,2 -k4,4g | sort -k1,1 -k2,2 -u
