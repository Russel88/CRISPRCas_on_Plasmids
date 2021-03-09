#!/bin/bash

# conda activate crisprcasfinder

# git clone https://github.com/dcouvin/CRISPRCasFinder.git

# PLSDB
mkdir -p CRISPRs/PLSDB/CRISPRFinder

cat PLSDB/acc_plasmid.txt \
	| parallel --nice 10 -j40 'mkdir CRISPRs/PLSDB/CRISPRFinder/{} && \
		cd CRISPRs/PLSDB/CRISPRFinder/{} && \
			perl ../../../../CRISPRCasFinder/CRISPRCasFinder.pl -in ../../../../PLSDB/fastas/{}.fna -out 'Here' &> {}.log'

# Archaea
mkdir -p CRISPRs/Archaea/CRISPRFinder

cat Archaea/acc_plasmid.txt \
        | parallel --nice 10 -j40 'mkdir CRISPRs/Archaea/CRISPRFinder/{} && \
                cd CRISPRs/Archaea/CRISPRFinder/{} && \
                        perl ../../../../CRISPRCasFinder/CRISPRCasFinder.pl -in ../../../../Archaea/fastas/{}.fna -out 'Here' &> {}.log'

# PLSDB Hosts
mkdir -p CRISPRs/Hosts/CRISPRFinder

cat PLSDB/acc_host.txt \
        | parallel --nice 10 -j40 'mkdir CRISPRs/Hosts/CRISPRFinder/{} && \
                cd CRISPRs/Hosts/CRISPRFinder/{} && \
                        perl ../../../../CRISPRCasFinder/CRISPRCasFinder.pl -in ../../../../Hosts/fastas/{}.fna -out 'Here' &> {}.log'

# Archaea Hosts
mkdir -p CRISPRs/Archaea_Hosts/CRISPRFinder

cat Archaea/acc_host.txt \
        | parallel --nice 10 -j40 'mkdir CRISPRs/Archaea_Hosts/CRISPRFinder/{} && \
                cd CRISPRs/Archaea_Hosts/CRISPRFinder/{} && \
                        perl ../../../../CRISPRCasFinder/CRISPRCasFinder.pl -in ../../../../Archaea_Hosts/fastas/{}.fna -out 'Here' &> {}.log'

