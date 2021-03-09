#!/bin/bash

wget https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=dr_34.zip
unzip 'DownloadFile?filename=dr_34.zip'
rm 'DownloadFile?filename=dr_34.zip'
mkdir -p CRISPRs/Repeats
mv 20210121_dr_34.fasta CRISPRs/Repeats/
cat CRISPRs/Repeats/20210121_dr_34.fasta | tr '\n' ' ' | sed 's/ >/\n>/g' | awk '$2 ~ /^[ACTG]*$/ {print ">X"NR"\n"$2}' > CRISPRs/Repeats/CCF_DB.fasta
