#!/bin/bash

# Download IMGVR database from homepage

# Get reference phages
awk -F'\t' '$10 == "Reference" {print $1}' IMGVR_all_Sequence_information.tsv > IMG_refs
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' IMG_refs IMGVR_all_nucleotides.fna > IMGVR_refs_nucleotides.fna

# Split spacers in chunks
mkdir PLSDB_spacers
mkdir Hosts_spacers
mkdir Archaea_spacers
mkdir Archaea_Hosts_spacers

awk -v var="200" -v name="Hosts_spacers/chunk" 'BEGIN {n=0;} /^>/ {if(n%var==0){f=sprintf(name"%d.fna",n);} print >> f; n++; next;} { print >> f; }' Spacers/Hosts_all.fna
awk -v var="200" -v name="PLSDB_spacers/chunk" 'BEGIN {n=0;} /^>/ {if(n%var==0){f=sprintf(name"%d.fna",n);} print >> f; n++; next;} { print >> f; }' Spacers/PLSDB_all.fna
awk -v var="200" -v name="Archaea_spacers/chunk" 'BEGIN {n=0;} /^>/ {if(n%var==0){f=sprintf(name"%d.fna",n);} print >> f; n++; next;} { print >> f; }' Spacers/Archaea_all.fna
awk -v var="200" -v name="Archaea_Hosts_spacers/chunk" 'BEGIN {n=0;} /^>/ {if(n%var==0){f=sprintf(name"%d.fna",n);} print >> f; n++; next;} { print >> f; }' Spacers/Archaea_Hosts_all.fna

# Submit spacer matching jobs
for i in /home/projects/cu_10108/data/Generated/Russel/CCP2/PLSDB_spacers/*; do sleep 1s; qsub -v INPUT="${i}" 404_spacer_plsdb.qsh; done
for i in /home/projects/cu_10108/data/Generated/Russel/CCP2/Hosts_spacers/*; do sleep 1s; qsub -v INPUT="${i}" 405_spacer_hosts.qsh; done
for i in /home/projects/cu_10108/data/Generated/Russel/CCP2/Archaea_spacers/*; do sleep 1s; qsub -v INPUT="${i}" 406_spacer_archaea.qsh; done
for i in /home/projects/cu_10108/data/Generated/Russel/CCP2/Archaea_Hosts_spacers/*; do sleep 1s; qsub -v INPUT="${i}" 407_spacer_archaea_hosts.qsh; done
