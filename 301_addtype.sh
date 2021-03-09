#!/bin/bash

for i in PLSDB Archaea Hosts Archaea_Hosts
do
    paste CRISPRs/${i}/crisprs_final.tab <(repeatType <(awk '{print $7}' CRISPRs/${i}/crisprs_final.tab) | cut -f2,3) > CRISPRs/${i}/crisprs_final_typed.tab
done
