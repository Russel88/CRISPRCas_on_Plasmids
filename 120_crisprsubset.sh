#!/bin/bash

# conda activate ccp2

subs1 () {

    # Good arrays
    awk '$14 == 4 || ($19 == "False" && $15 >= 70 && $16 < 50 && $18/sqrt($9) < 3)' CRISPRs/${1}/crisprs.tab > CRISPRs/${1}/crisprs_first_subset.tab

    # Quarantined arrays
    awk '$14 != 4 && $9 > 1 && ($19 == "True" || $15 < 70 || $16 >= 50 || $18/sqrt($9) >= 3)' CRISPRs/${1}/crisprs.tab > CRISPRs/${1}/crisprs_quarantine.tab

    # Get repeats from high quality arrays
    awk '$14 == 4' CRISPRs/${1}/crisprs_first_subset.tab | awk '{print $7}' | sort | uniq | awk '{sub("^",">"NR"\n",$1)}1' > CRISPRs/Repeats/${1}.fasta
}

subs2 () {

    # Rescue quarantined arrays by repeat matching
    awk '{print $1,$7}' CRISPRs/${1}/crisprs_quarantine.tab | sed 's|^|>|;s| |\n|' > CRISPRs/${1}/repeats_quarantine.fa

    blastn -db CRISPRs/Repeats/nr_${2} -query CRISPRs/${1}/repeats_quarantine.fa -task blastn-short -perc_identity 90 -qcov_hsp_perc 90 -dust no -soft_masking false -outfmt=6 \
        | cut -f1 | sort | uniq > CRISPRs/${1}/rescue_repeats.txt

    # Rescue quarantiend arrays by adjacency to cas operons
    ## Get predicted cas operons
    awk '$4 != "False" && $4 != "Partial"' Cas/${1}_prediction.tab > Cas/${1}_prediction_sub.tab

    ## Tables with cas operon positions
    tail -n+2 Cas/${1}_prediction_sub.tab | awk '{print $1,sub("@.*","",$1),$1,$2,$3}' | awk '{print $3,$1,$4,$5}' | sort -k1,1 > Cas/${1}_pos.tab

    ## Tables with array positions
    awk '{print $2,$1,$3,$4}' CRISPRs/${1}/crisprs_quarantine.tab | sort -k1,1 > CRISPRs/${1}/crisprs_pos.tab

    ## Merge arrays and operon tables and get array id when distance is within 1000 bp
    join -j 1 CRISPRs/${1}/crisprs_pos.tab Cas/${1}_pos.tab \
        | awk '{print $2,$3-$7,$6-$4}' \
        | awk '($2 > 0 && $2 <= 1000) || ($3 > 0 && $3 <= 1000) || ($2 < 0 && $3 < 0)' \
        | cut -d' ' -f1 | sort | uniq > CRISPRs/${1}/rescue_cas.txt

    ## Get rescued arrays
    grep -f <(cat CRISPRs/${1}/rescue_repeats.txt CRISPRs/${1}/rescue_cas.txt) CRISPRs/${1}/crisprs_quarantine.tab > CRISPRs/${1}/crisprs_rescued.tab

    # Make final CRISPR array table
    cat CRISPRs/${1}/crisprs_first_subset.tab CRISPRs/${1}/crisprs_rescued.tab > CRISPRs/${1}/crisprs_final.tab

    # Make masking file
    awk '{print $2,$3,$4}' CRISPRs/${1}/crisprs_final.tab | sed 's| |,|g' > CRISPRs/${1}/crisprs_mask.csv

}

subs1 "PLSDB"
subs1 "Hosts"
subs1 "Archaea"
subs1 "Archaea_Hosts"

# Make repeat database - PLSDB
cat CRISPRs/Repeats/CCF_DB.fasta CRISPRs/Repeats/Hosts.fasta CRISPRs/Repeats/PLSDB.fasta > CRISPRs/Repeats/all_plsdb.fasta
cd-hit-est -c 1 -i CRISPRs/Repeats/all_plsdb.fasta -o CRISPRs/Repeats/nr_plsdb.fasta
makeblastdb -dbtype nucl -in CRISPRs/Repeats/nr_plsdb.fasta -out CRISPRs/Repeats/nr_plsdb

# Make repeat database - Archaea
cat CRISPRs/Repeats/CCF_DB.fasta CRISPRs/Repeats/Archaea_Hosts.fasta CRISPRs/Repeats/Archaea.fasta > CRISPRs/Repeats/all_arch.fasta
cd-hit-est -c 1 -i CRISPRs/Repeats/all_arch.fasta -o CRISPRs/Repeats/nr_arch.fasta
makeblastdb -dbtype nucl -in CRISPRs/Repeats/nr_arch.fasta -out CRISPRs/Repeats/nr_arch

subs2 "PLSDB" "plsdb"
subs2 "Hosts" "plsdb"
subs2 "Archaea" "arch"
subs2 "Archaea_Hosts" "arch"

