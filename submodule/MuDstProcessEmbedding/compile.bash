#!/bin/bash
# compile.bash -- lightflavorspectra_OO200/submodule/MuDstProcessEmbedding
#
# Ported from lightflavorspectra_etof/submodule/MuDstProcessEmbedding/compile.bash.
#
# WHAT CHANGED FROM THE ORIGINAL, AND WHY:
#   - Dropped the MattMcFDEvent/MattMcFDTrack (feed-down-embedding) .so/.h copies --
#     this repo hasn't ported FeedDownEmbedding/HybridFeeddownMaker's RCF-side maker
#     yet (hybrid feed-down needs a SEPARATE embedding request with its own MuDst
#     sample, decay-chain-matched via a mother-particle geant ID -- not the plain
#     TPC-efficiency sample this pipeline targets). Add those back the same way
#     as the rest of this file once that request exists and gets ported.
#   - Dropped the StRefMultCorr_SL19b/SL21a/SL22dev/7p7GeV/9p2_11p5_17p3GeV .so copies
#     -- none of those match O+O 200 GeV.
#     RESOLVED 2026-09-01: checked STAR's public star-sw (StRoot/StRefMultCorr/Param.h)
#     against this repo's own xml/runPicoBinner_OO200_SDCC_template.xml run whitelist
#     (trgsetupname=production_OO_200GeV_2021). A real, official calibration already
#     exists there for O+O 200 GeV -- Param.h's mParamStr_ref6[0]/mParamStr_ref7[0]
#     ("Run 21 O+O 200 GeV, Trigger ID = 860001/860002/860011/860012"), selected via
#     StRefMultCorr("refmult6") or StRefMultCorr("totnMIP") for run numbers in
#     [22130029, 22144006]. All runs in this repo's picoBinner XML whitelist fall inside
#     that window (since Oct 2026 the whitelist is the 56 good trigger-860003 runs,
#     22133032-22136010; runlists/OO200_860003_runs.txt), so
#     it IS the right calibration IF this build were to use StRefMultCorr at all --
#     confirm your trigger IDs match 860001/860002/860011/860012 too if you ever do.
#     BUT: checked against two independent sources whether this embedding Maker SHOULD
#     use it, and both say no. (1) This repo's own data-side pipeline
#     (source/PicoBinner.cxx, macros/SetCutClass.C) deliberately never defines
#     _STREFMULTCORR_, using raw refMult centrality bins instead ("matching the
#     original's _OO_200_COL_ config"). (2) The real reference repo this was ported from
#     (lightflavorspectra_etof)'s actual top-level makefile_toggles.h has a genuine RCF
#     branch for _OO_200_COL_ (#ifndef _MAC_OSX_ / #ifdef _OO_200_COL_), and it ALSO never
#     defines _STREFMULTCORR_ there -- unlike every other collision energy in that same
#     file, which all pair _STREFMULTCORR_ with an SL/energy variant on RCF. So "no
#     StRefMultCorr for O+O 200 GeV, on RCF too" is a real, twice-confirmed, deliberate
#     choice, not an oversight -- using it here while the data side doesn't would bin
#     embedding efficiency into different centrality classes than the spectra it
#     corrects, a real physics bug. See submodule/MuDstProcessEmbedding/makefile_toggles_RCF.h
#     for where this is now encoded (as "leave _STREFMULTCORR_ undefined", not by
#     removing StRefMultCorr from doMuDst.C's _SL*_STREFMULTCORR_ gSystem->Load blocks,
#     which already no-op correctly when their guarding macro is never defined).
#   - submodule/ copies trimmed to what this repo actually has: ParticleInfo and
#     PicoDstReader_SL24y (instead of the original's SL19b/SL22b/SL22c/SL23a/SL23d
#     readers and a generic submodule/StRoot).
#   - Changed 2026-10-28: SL24y only. Stages submodule/PicoDstReader_SL24y/, the same
#     reader the PicoBinner/SDCC build uses (the production and the P24iy embedding are
#     both SL24y). The older P23ic-era reader and setup were removed; they are in git
#     history before this change.
#   - IMPORTANT: this maker runs under STAR's official offline software (root4star,
#     StChain, StMuDstMaker) at RCF/SDCC -- a fundamentally different runtime than the
#     rest of this repo, which ../../makefile_toggles.h deliberately hardcodes to a
#     single LOCAL-only build (`_MAC_OSX_` always defined, no _STREFMULTCORR_ toggle
#     at all -- see that file's own header comment). The `_MAC_OSX_` toggle there
#     means "not RCF" upstream, so building this maker against that toggle file as-is
#     will select the WRONG code branches in CutClass.h/HistogramUtilities.h (it skips
#     StPicoEvent*/StPicoTrack*-pointer overloads the RCF/StMuDst side may need, per
#     that file's own comment).
#     RESOLVED 2026-09-01: added ./makefile_toggles_RCF.h (this directory), the
#     RCF-flavored second toggles file this TODO asked for -- see its own header
#     comment for the full reasoning per toggle. The `cp ../../makefile_toggles.h ...`
#     lines below have been repointed at it.
#     headers/CutClass.h selects the SL24y headers for this build through
#     _PICO_READER_SL24y_ (defined in makefile_toggles_RCF.h).

echo --------------------------------------------------
echo NOTE: this must be run at RCF/SDCC inside the STAR software environment
echo   "(after starver SL24y) -- not in this cloud sandbox."

mkdir -p ./bin_temp/
cp ../../bin/HistogramUtilities_cxx.so ./bin_temp/
cp ../../bin/PhysMath_cxx.so ./bin_temp/
cp ../../bin/namespaces_cxx.so ./bin_temp/
echo NOTE: ------  It is okay to not have HistoUtil+PhysMath or namespaces but Not both
cp ../../bin/CutClass_cxx.so ./bin_temp/
cp ../../bin/SetCutClass_C.so ./bin_temp/
cp ../../bin/ParticleInfo_cxx.so ./bin_temp/
# No StRefMultCorr .so copy here -- _STREFMULTCORR_ is deliberately left undefined for
# this build (see header note above); doMuDst.C's _SL*_STREFMULTCORR_-gated
# gSystem->Load calls correctly no-op when their guarding macro is never defined.
cp ./makefile_toggles_RCF.h ./makefile_toggles.h
cp ../../bin/MattMcEvent_cxx.so ./bin_temp/
cp ../../bin/MattMcTrack_cxx.so ./bin_temp/
cp ../../bin/Helix_cxx.so ./bin_temp/

mkdir -p ./headers/
cp ../../headers/HistogramUtilities.h ./headers/
cp ../../headers/PhysMath.h ./headers/
cp ../../headers/CutClass.h ./headers/
cp ./makefile_toggles_RCF.h ./headers/makefile_toggles.h
cp ../../headers/MattMcEvent.h ./headers/
cp ../../headers/MattMcTrack.h ./headers/
cp ../../headers/Helix.h ./headers/

mkdir -p ./submodule/
cp -r ../../submodule/ParticleInfo ./submodule/
cp -r ../../submodule/PicoDstReader_SL24y ./submodule/
# NOTE: StMuDstMaker itself (StChain/StMuDstMaker/StMuTrack/...) comes from the STAR
# software environment via `starver`, not from this repo -- nothing to copy for that.

cons
