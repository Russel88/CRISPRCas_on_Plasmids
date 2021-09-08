#!/bin/bash

# conda activate gtdbtk

mkdir -p Archaea/gcf
echo "Archaea"
for i in $(cut -f21 Archaea/archaea_assembly_metadata.tsv | tail -n+2 | sort | uniq | grep "ftp")
    do
    NAME=$(echo $i | sed 's|.*\/||;s|\(GC._[0-9]*\)\.[0-9]_.*|\1|')
    echo $NAME
    if [[ ! -f "Archaea/gcf/${NAME}.fna.gz" ]]
    then
        wget -q -nd -A'genomic.fna.gz' -R'from_genomic.fna.gz' -r $(echo $i | sed 's|GC._.*||')
        for k in ${NAME}*.gz
        do
            mv $k Archaea/gcf/${NAME}.fna.gz
        done
    fi
done

mkdir -p PLSDB/gcf
echo "PLSDB"
for i in $(cut -f22 PLSDB/plsdb_assembly_metadata.tsv | tail -n+2 | sort | uniq | grep "ftp")
    do
    NAME=$(echo $i | sed 's|.*\/||;s|\(GC._[0-9]*\)\.[0-9]_.*|\1|')
    echo $NAME
    if [[ ! -f "PLSDB/gcf/${NAME}.fna.gz" ]]
    then
        wget -q -nd -A'genomic.fna.gz' -R'from_genomic.fna.gz' -r $(echo $i | sed 's|GC._.*||')
        for k in ${NAME}*.gz
        do
            mv $k PLSDB/gcf/${NAME}.fna.gz
        done
    fi
done

nice -n10 gtdbtk classify_wf --genome_dir Archaea/gcf --out_dir Archaea/gtdb --cpus 10 --extension ".fna.gz"
nice -n10 gtdbtk classify_wf --genome_dir PLSDB/gcf --out_dir PLSDB/gtdb --cpus 20 --extension ".fna.gz"
