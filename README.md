# CRISPR-Cas on Plasmids
This is the code repository for reproducing analyses in the paper:

[CRISPR-Cas systems are widespread accessory elements across bacterial and archaeal plasmids](https://doi.org/10.1093/nar/gkab859)

# Software
* See Conda/*.yml files for most software
* BLAST 2.8.1+
* CRISPRCasFinder 4.2.17
* R 3.6.3
    * ape 5.3
    * foreach 1.5.0
    * ggplot2 3.3.0
    * MASS 7.3-51.6
    * mixtools 1.2.0
    * eulerr 6.1.0
    * indicspecies 1.7.9

# Code
* 0xx: Downloading and preparation
    * 000_archaea.sh - Archaeal plasmids
    * 001_download.sh - Download PLSDB info and plasmid sequences
    * 002_host_information.py - Get assembly info for plasmids
    * 003_assembly_information.py - Download assembly reports
    * 004_host_table.py - Create table connecting plasmids and their chromosome hosts
    * 005_download_host.sh - Download host chromosome sequences
    * 006_prodigal_plasmids.sh - Find ORFs for plasmids
    * 007_prodigal_hosts.sh - Find ORFs for chromosomes
    * 008_mobtyper.sh - Mobility and Inc typing of plasmids
    * 009_genepos.sh - Get positions of high confidence ORFs
    * 010_dbDownload.sh - Get CRISPR repeat db
    * 011_tax.sh - Taxonomy of chromosomes
    * 012_drep.sh - Dereplicate plasmids
    * 013_size.sh - Sizes of sequences

* 1xx: CRISPRs
    * 101_CRISPRFinder.sh - Run CRISPRFinder
    * 102_extract_crispr.sh - Extract CRISPR tables
    * 103_crispy_stats.py - Calculate stats for CRISPRs
    * 104_orf_overlap.py - Determine overlap with ORFs (not very efficient...)
    * 105_crisprtab.sh - Combine all in one table
    * 120_crisprsubset.sh - Subset CRISPRs based on decision tree
    * 121_mask.sh - Mask the sequences
    * 122_repeatmatch.sh - BLAST repeats against masked sequences
    * 123_repeatmatch.py - Combine BLAST results into putative CRISPRs
    * 124_orf_array.py - Determine if new CRISPRs overlap with ORFs
    * 125_extract.sh - Extract new CRISPRs of good quality
    * 126_mask_again.sh - Re-mask the sequences for spacer-matching later
    * 127_repeat_sim.py - Calculate repeat similarity for typing CRISPRs distant to cas operons

* 2xx Cas Operons
    * 200_castyping.sh - Run CCTyper to find cas operons
    * 201_cas_tab.sh - Collect cas operons

* 3xx CRISPR-Cas
    * 300_crispr_cas.sh - Combine CRISPR and cas
    * 301_addtype.sh - Predict the subtype based on CRISPR repeat sequence

* 4xx Spacers
    * 400_spacers.sh - Get spacers from CRISPRFinder results
    * 401_spacers_new_arrays.sh - Get spacers from new arrays
    * 402_combine.sh - Combine all spacers
    * 403_spacer_match.sh - Download virus genomes and submit spacer matching scripts
    * 404/405/406/407_spacer_plsdb/hosts/archaea/archaea_hosts.qsh - Spacer matching script for computing cluster
    * 408_spacers.sh - Collect spacer results and filter by E-value
    * 409_spacers.sh - Get only ORF matching spacer matches

* 900_collect.sh - Collect files for analysis in R

* R/*.R - R scripts for statistics and figures
