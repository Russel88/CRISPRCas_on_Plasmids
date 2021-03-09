#!/bin/bash

# conda activate ccp2

mkdir -p PLSDB/network
# Get scripts from https://github.com/macman123/plasmid_network_analysis
git clone https://github.com/macman123/plasmid_network_analysis PLSDB/network

# Get paths to plasmids
ls PLSDB/fastas/ | sed 's/^/..\/..\/fastas\//' > PLSDB/network/scoring_similarity/plasmid_paths.txt

# Run
cd PLSDB/network/scoring_similarity/
chmod +x split_in_batches.sh
chmod +x pairwise_dists.sh
sed -i 's/bindash=\/path\/to\/bindash/bindash=bindash/' pairwise_dists.sh
./split_in_batches.sh
./pairwise_dists.sh
paste <(cut -f1,2 bindash_out.tsv | sed 's|.fna||g;s|......fastas.||g') <(cut -f5 bindash_out.tsv | awk -F'/' '{print 1-$1/$2}') > bindash.tab
awk '$3 < 1' bindash.tab > bindash_sub.tab
gzip bindash_out.tsv
gzip bindash.tab
gzip bindash_sub.tab
