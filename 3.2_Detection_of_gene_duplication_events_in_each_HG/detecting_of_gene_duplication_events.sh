#!/bin/bash

#
# Predict gene trees for each homologous groups
#

./scripts/predicting_phylogenetic_trees/predicting_phylogenetic_trees.py --input groups/groupsID.dump.seq.mci.I11_including_ALL_singletons --tmp tmp/ --output output/ --threads 2 --min 3 --extra n -p

#
# Process gene trees with Notung 2.6
#

./scripts/run_Notung_on_folder.py -i output/ -s speciesTree.nwk -o output_Notung/ -p .nt_ali_cleaned.phy_phyml_tree.txt

#
# Predict gene duplication events
#

./scripts/detect_duplication_events_-_for_all_species.py -i output_Notung/ -o recent_duplications.100pb.bootstrap.0.9.tsv -b 0.9 -p '.0'
./scripts/detect_duplication_events_-_for_NIATT_specific.py -i output_Notung/ -o recent_duplications_NIATT_only_100pb.bootstrap.0.9.tsv -b 0.9 -p '.0'
