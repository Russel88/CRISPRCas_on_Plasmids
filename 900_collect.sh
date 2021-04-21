#!/bin/bash

mkdir -p Collect

cat Cas/Archaea_prediction.tab Cas/PLSDB_prediction.tab > Collect/Cas_plasmid.tab
cat Cas/Archaea_Hosts_prediction.tab Cas/Hosts_prediction.tab > Collect/Cas_host.tab

cat <(awk '{print $2,$1,$3,$4,$9+1,$20,$21}' CRISPRs/Archaea/crisprs_final_typed.tab) <(awk -F'\t' '{print $2,$1,$3,$4,$7,$13,$14}' CRISPRs/Archaea/RepeatMatch/arrays_typed.tab) <(awk '{print $2,$1,$3,$4,$9+1,$20,$21}' CRISPRs/PLSDB/crisprs_final_typed.tab) <(awk -F'\t' '{print $2,$1,$3,$4,$7,$13,$14}' CRISPRs/PLSDB/RepeatMatch/arrays_typed.tab) > Collect/CRISPR_plasmid.tab
cat <(awk '{print $2,$1,$3,$4,$9+1,$20,$21}' CRISPRs/Archaea_Hosts/crisprs_final_typed.tab) <(awk -F'\t' '{print $2,$1,$3,$4,$7,$13,$14}' CRISPRs/Archaea_Hosts/RepeatMatch/arrays_typed.tab) <(awk '{print $2,$1,$3,$4,$9+1,$20,$21}' CRISPRs/Hosts/crisprs_final_typed.tab) <(awk -F'\t' '{print $2,$1,$3,$4,$7,$13,$14}' CRISPRs/Hosts/RepeatMatch/arrays_typed.tab) > Collect/CRISPR_host.tab

cat CRISPRCas/Archaea_crispr_cas_distance_nearest.tab CRISPRCas/PLSDB_crispr_cas_distance_nearest.tab > Collect/CRISPRCas_plasmid.tab
cat CRISPRCas/Archaea_Hosts_crispr_cas_distance_nearest.tab CRISPRCas/Hosts_crispr_cas_distance_nearest.tab > Collect/CRISPRCas_host.tab

join -t$'\t' -j1 <(tail -n+2 PLSDB/plsdb_host_etc.tsv | sort -k1,1) <(cut -f1,5 PLSDB/plsdb_assembly_metadata.tsv | tail -n+2 | sort -k1,1) \
    | sort -k4,4 -u | awk -F'\t' 'BEGIN{OFS="\t"}{print $5,$2,$3,$4}' > Collect/cell_plsdb.tab

join -t$'\t' -j1 <(tail -n+2 Archaea/archaea_host_etc.tsv | sort -k1,1) <(cut -f1,4 Archaea/archaea_assembly_metadata.tsv | tail -n+2 | sort -k1,1) \
    | sort -k4,4 -u | awk -F'\t' 'BEGIN{OFS="\t"}{print $5,$2,$3,$4}' > Collect/cell_archaea.tab


cut -f2,5,9,25,43 PLSDB/plsdb.tsv > Collect/meta_plsdb.tab

cat <(cut -f1,2 PLSDB/gtdb/gtdbtk.bac120.summary.tsv | tail -n+2) <(cut -f1,2 Archaea/gtdb/gtdbtk.ar122.summary.tsv | tail -n+2) |  sed 's/.fna//;s/d__//;s/;[a-z]__/\t/g' > Collect/tax.tab

cp Archaea/gtdb/gtdbtk.ar122.classify.tree Collect/archaea.tree
cp PLSDB/gtdb/gtdbtk.bac120.classify.tree Collect/plsdb.tree

cp PLSDB/drep/data_tables/Cdb.csv Collect/PLSDB_cluster.csv
cp PLSDB/drep/data_tables/Wdb.csv Collect/PLSDB_drep.csv

cp Archaea/drep/data_tables/Cdb.csv Collect/Archaea_cluster.csv
cp Archaea/drep/data_tables/Wdb.csv Collect/Archaea_drep.csv

cp CRISPRs/Repeat_sim.csv Collect/

cat PLSDB/size.tab Archaea/size.tab > Collect/size_plasmid.tab
cat Hosts/size.tab Archaea_Hosts/size.tab > Collect/size_host.tab

cp Spacers/*.m8 Collect/

ln -s ${PWD}/PLSDB/network/scoring_similarity/bindash_sub.tab.gz Collect/bindash_sub.tab.gz

cat <(tail -q -n+2 PLSDB/mob/*) <(tail -q -n+2 Archaea/mob/*) | cut -f1,6,14 | sed 's/\.fna//' > Collect/mobility.tab
