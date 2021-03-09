#!/bin/bash

for i in PLSDB Hosts Archaea Archaea_Hosts
do
    for j in IMG plsdb
    do
        cat ${i}_spacers/*_${j}.m8 | awk '$11 <= 0.05' > Spacers/${i}_${j}.m8
    done
done

