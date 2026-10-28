#!/bin/bash
# runBatchSpecificRequest_OO200GeV.bash -- lightflavorspectra_OO200/submodule/MuDstProcessEmbedding
#
# Ported from lightflavorspectra_etof/submodule/MuDstProcessEmbedding/runBatchSpecificRequest_19p6GeV.bash.
# Dispatches runSchedSubmission_OO200.bash once per (species, charge, RequestID) for the
# P24iy/SL24y O+O embedding. Add one line per RequestID if the request is split (e.g. a
# separate high-pT sample; see the 19.6 GeV example in the sibling etof repo).
#
# Changed 2026-10-28: SL24y only. Set REQ to the P24iy request's RequestID once it exists.
# (The earlier P23ic request, 20232003, is in git history.)

REQ=SET_P24IY_REQUEST_ID
LIB=SL24y
if [ "$REQ" = "SET_P24IY_REQUEST_ID" ]; then
  echo "runBatchSpecificRequest_OO200GeV.bash: set REQ to the P24iy embedding RequestID first"; exit 1
fi

./runSchedSubmission_OO200.bash 0  1 $REQ 200 $LIB   # pi+
./runSchedSubmission_OO200.bash 0 -1 $REQ 200 $LIB   # pi-

./runSchedSubmission_OO200.bash 1  1 $REQ 200 $LIB   # K+
./runSchedSubmission_OO200.bash 1 -1 $REQ 200 $LIB   # K-

./runSchedSubmission_OO200.bash 2  1 $REQ 200 $LIB   # proton
./runSchedSubmission_OO200.bash 2 -1 $REQ 200 $LIB   # antiproton

exit
