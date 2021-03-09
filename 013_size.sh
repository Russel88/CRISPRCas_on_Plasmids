#!/bin/bash

for i in Archaea/fastas/*.fna
do
    NAME=$(echo $i | sed 's/.*\///;s/\.fna//')
    printf "${NAME}\t"
    grep -v ">" $i | tr -d '\n' | wc -m
done > Archaea/size.tab


for i in PLSDB/fastas/*.fna
do
    NAME=$(echo $i | sed 's/.*\///;s/\.fna//')
    printf "${NAME}\t"
    grep -v ">" $i | tr -d '\n' | wc -m
done > PLSDB/size.tab


for i in Archaea_Hosts/fastas/*.fna
do
    NAME=$(echo $i | sed 's/.*\///;s/\.fna//')
    printf "${NAME}\t"
    grep -v ">" $i | tr -d '\n' | wc -m
done > Archaea_Hosts/size.tab


for i in Hosts/fastas/*.fna
do
    NAME=$(echo $i | sed 's/.*\///;s/\.fna//')
    printf "${NAME}\t"
    grep -v ">" $i | tr -d '\n' | wc -m
done > Hosts/size.tab


# Fake sizes for linear PLSDB
awk -F'\t' '$5 == "linear" {print $2}' PLSDB/plsdb.tsv > PLSDB/linear
cat <(grep -v -f PLSDB/linear PLSDB/size.tab) <(sed 's/$/\t1000000000/' PLSDB/linear) > PLSDB/size2.tab

for k in Archaea_Hosts Hosts Archaea
do
    ln -f -s ${PWD}/${k}/size.tab ${k}/size2.tab
done
