#!/bin/bash

for k in Archaea Archaea_Hosts PLSDB Hosts
do
    for i in Cas/${k}/*;
    do
        tail -n+2 $i/cas_operons.tab
        tail -n+2 $i/cas_operons_putative.tab
    done | awk '$5 != "False"' | awk '$5 != "Putative"' > Cas/${k}_prediction.tab
done

