# CRISPR-Cas on Plasmids
This is the code repository for reproducing analyses in the paper:

[CRISPR-Cas systems are widespread accessory elements across plasmids]()

# Code Structure
* 0xx: Downloading and preparation
* 1xx: CRISPRs
	* 10x: CRISPRFinder and thresholds
	* 12x: Repeat matching, SRUs and mini-arrays (run after Cas operon scripts)
* 2xx Cas Operons
* 3xx CRISPR-Cas
* 4xx Spacers

# Dir Structure
* PLSDB: PLSDB files, metadata, etc.
* Hosts
* CRISPRs: CRISPR array predictions
* Cas: Cas operon predictions
* CRISPRCas

# Software
* See conda_*.yml files for most software
* BLAST 2.8.1+
* CRISPRCasFinder 4.2.17
