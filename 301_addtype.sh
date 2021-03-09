#!/bin/bash

# conda activate ccp2

# Download newest model
wget http://mibi.galaxy.bio.ku.dk/russel/repeattyper/2021/repeat_model_2021_03_01.tgz
tar -xvzf repeat_model_2021_03_01.tgz
rm repeat_model_2021_03_01.tgz
rm xgb_report

for i in PLSDB Archaea Hosts Archaea_Hosts
do
    paste CRISPRs/${i}/crisprs_final.tab <(repeatType --db repeat_model <(awk '{print $7}' CRISPRs/${i}/crisprs_final.tab) | cut -f2,3) > CRISPRs/${i}/crisprs_final_typed.tab
done
