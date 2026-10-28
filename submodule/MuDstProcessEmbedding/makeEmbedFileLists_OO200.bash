#!/bin/bash
# makeEmbedFileLists_OO200.bash -- lightflavorspectra_OO200/submodule/MuDstProcessEmbedding
#
# Ported from lightflavorspectra_etof/submodule/MuDstProcessEmbedding/makeEmbedFileLists_19p6.bash,
# which is itself the README's own worked example of "how to make a filelist" for a new energy.
#
# Changed 2026-10-28: SL24y only. Set REQ to the P24iy embedding RequestID once it exists.
# Production directory: production_OO_200GeV_2021 under /star/data105/embedding (also
# mounted as /sdcc/lustre02/star/data105/embedding).
#
# The directory globs below follow the earlier P23ic request (20232003; in git history):
# each species was split into many numbered sub-directories with an "RFF" tag
# (PiplusRFF_20232003, PiplusRFF_100_20232003, ...), and protons also had an untagged
# Proton_* family that was excluded. CHECK the P24iy directory names with
#   find /star/data105/embedding/production_OO_200GeV_2021 -maxdepth 1 -name "*_$REQ"
# and adjust the species tags before using the lists.

REQ=SET_P24IY_REQUEST_ID
if [ "$REQ" = "SET_P24IY_REQUEST_ID" ]; then
  echo "makeEmbedFileLists_OO200.bash: set REQ to the P24iy embedding RequestID first"; exit 1
fi

echo Removing old 200GeV Embedding Filelists...
mkdir -p ./fileLists/200GeV_OO200/full/
echo Making new 200GeV Embedding Filelists...

EMBDIR=/star/data105/embedding/production_OO_200GeV_2021

echo Starting PionPlus...
find $EMBDIR/PiplusRFF*_${REQ} -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pip_${REQ}.list

echo Starting PionMinus...
find $EMBDIR/PiminusRFF*_${REQ} -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pim_${REQ}.list

echo Starting KaonPlus...
find $EMBDIR/KplusRFF*_${REQ} -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_Kp_${REQ}.list

echo Starting KaonMinus...
find $EMBDIR/KminusRFF*_${REQ} -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_Km_${REQ}.list

echo Starting ProtonPlus...
find $EMBDIR/ProtonRFF*_${REQ} -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pro_${REQ}.list

echo Starting ProtonMinus...
find $EMBDIR/AntiProtonRFF*_${REQ} -print | grep -i MuDst.root | awk '{print $0}' > ./fileLists/200GeV_OO200/full/200GeV_pbar_${REQ}.list


echo Number of pip embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pip_${REQ}.list | wc -l

echo Number of pim embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pim_${REQ}.list | wc -l

echo Number of Kp embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_Kp_${REQ}.list | wc -l

echo Number of Km embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_Km_${REQ}.list | wc -l

echo Number of pro embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pro_${REQ}.list | wc -l

echo Number of pbar embedding files:
cat ./fileLists/200GeV_OO200/full/200GeV_pbar_${REQ}.list | wc -l
