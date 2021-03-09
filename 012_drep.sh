#!/bin/bash

# conda activate drep

#sed 's|^|PLSDB/fastas/|;s|$|.fna|' PLSDB/acc_plasmid.txt > plsdb_where
#dRep dereplicate PLSDB/drep -l 1 -pa 0.9 -sa 0.95 --ignoreGenomeQuality --S_algorithm "fastANI" -cm "total" -nc 0.9 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -g plsdb_where
#rm plsdb_where

sed 's|^|Archaea/fastas/|;s|$|.fna|' Archaea/acc_plasmid.txt > archaea_where
dRep dereplicate Archaea/drep -l 1 -pa 0.9 -sa 0.95 --ignoreGenomeQuality --S_algorithm "fastANI" -cm "total" -nc 0.9 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -g archaea_where
rm archaea_where

