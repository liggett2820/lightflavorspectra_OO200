#!/bin/bash
# makeEmbedFileLists_OO200.bash -- lightflavorspectra_OO200/submodule/MuDstProcessEmbedding
#
# Ported from lightflavorspectra_etof/submodule/MuDstProcessEmbedding/makeEmbedFileLists_19p6.bash,
# which is itself the README's own worked example of "how to make a filelist" for a new energy.
#
# REAL VALUES CONFIRMED (per Andrew, 2026-08-05+08-06), replacing all earlier placeholders:
#   - RequestID: 20232003 (single combined ID, not one per species)
#   - Production directory: production_OO_200GeV_2021
#   - RCF root: /star/data105/embedding (confirmed present; identical content also exists
#     under /sdcc/lustre02/star/data105/embedding -- same lustre storage under two mount
#     paths, use whichever is reachable from wherever this actually runs)
#
# NAMING CONVENTION IS DIFFERENT FROM THE 19.6 GeV REFERENCE SCRIPT: this production
# doesn't use simple `Species_*_RequestID` directories -- confirmed via `find
# -ipath "*production_OO_200GeV_2021*20232003*"` that each species has an "RFF" tag plus
# a large number of numbered sub-requests (e.g. PiplusRFF_20232003, PiplusRFF_100_20232003,
# PiplusRFF_101_20232003, ... PiplusRFF_202_20232003 -- 12 directories for pi+ alone, all
# under the SAME RequestID 20232003, presumably split by momentum/eta sub-range at request
# time). The `Species*_20232003` glob below (no underscore between the species tag and the
# wildcard) matches ALL of them at once -- both the bare/base directory and every numbered
# variant -- since `find` expands a multi-directory glob into separate search roots.
#
# PROTON HAS TWO TAG FAMILIES: the find output shows BOTH `Proton_*_20232003` (4 dirs)
# AND `ProtonRFF_*_20232003` (5 dirs) -- every other species only has the RFF-tagged
# family. Confirmed with Andrew (2026-08-06): only ProtonRFF* should be used, matching
# every other species -- the plain Proton_* set is deliberately excluded below.

echo Removing old 200GeV Embedding Filelists...
mkdir -p ./fileLists/200GeV_OO200/full/
echo Making new 200GeV Embedding Filelists...

EMBDIR=/star/data105/embedding/production_OO_200GeV_2021

echo Starting PionPlus...
find $EMBDIR/PiplusRFF*_20232003 -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pip_20232003.list

echo Starting PionMinus...
find $EMBDIR/PiminusRFF*_20232003 -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pim_20232003.list

echo Starting KaonPlus...
find $EMBDIR/KplusRFF*_20232003 -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_Kp_20232003.list

echo Starting KaonMinus...
find $EMBDIR/KminusRFF*_20232003 -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_Km_20232003.list

echo Starting ProtonPlus...
find $EMBDIR/ProtonRFF*_20232003 -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pro_20232003.list

echo Starting ProtonMinus...
find $EMBDIR/AntiProtonRFF*_20232003 -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pbar_20232003.list


echo Number of pip embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pip_20232003.list | wc -l

echo Number of pim embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pim_20232003.list | wc -l

echo Number of Kp embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_Kp_20232003.list | wc -l

echo Number of Km embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_Km_20232003.list | wc -l

echo Number of pro embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pro_20232003.list | wc -l

echo Number of pbar embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pbar_20232003.list | wc -l
