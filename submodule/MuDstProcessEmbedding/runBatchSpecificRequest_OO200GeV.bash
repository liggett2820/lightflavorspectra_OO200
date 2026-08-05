#!/bin/bash
# runBatchSpecificRequest_OO200GeV.bash -- lightflavorspectra_OO200/submodule/MuDstProcessEmbedding
#
# Ported from lightflavorspectra_etof/submodule/MuDstProcessEmbedding/runBatchSpecificRequest_19p6GeV.bash.
# Dispatches runSchedSubmission_OO200.bash once per (species, charge, RequestID) --
# add one line per RequestID your P23ic embedding request actually has (some requests
# split a species/charge across more than one RequestID, e.g. a normal-pT sample plus
# a separate high-pT sample -- see the original 19.6 GeV example for that shape).
#
# RequestID filled in as 20232003 (per Andrew, 2026-08-05) -- ASSUMED to be a single
# combined RequestID covering all 6 species/charges rather than one per species. If your
# P23ic request actually split species across separate RequestIDs, or split normal-pT vs
# a separate high-pT sample, update the relevant line(s) below and add extra lines the
# same way (see runBatchSpecificRequest_19p6GeV.bash in the sibling etof repo for that
# shape). library=SL23c below matches what you told me your P23ic embedding was produced
# with -- change it if that's wrong.

./runSchedSubmission_OO200.bash 0  1 20232003 200 SL23c   # pi+
./runSchedSubmission_OO200.bash 0 -1 20232003 200 SL23c   # pi-

./runSchedSubmission_OO200.bash 1  1 20232003 200 SL23c   # K+
./runSchedSubmission_OO200.bash 1 -1 20232003 200 SL23c   # K-

./runSchedSubmission_OO200.bash 2  1 20232003 200 SL23c   # proton
./runSchedSubmission_OO200.bash 2 -1 20232003 200 SL23c   # antiproton

exit
