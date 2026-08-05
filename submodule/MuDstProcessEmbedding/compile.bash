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
#     TPC-efficiency P23ic sample this pipeline targets). Add those back the same way
#     as the rest of this file once that request exists and gets ported.
#   - Dropped the StRefMultCorr_SL19b/SL21a/SL22dev/7p7GeV/9p2_11p5_17p3GeV .so copies
#     -- none of those match O+O 200 GeV. TODO: you need to determine which
#     StRefMultCorr class (if any) applies to your run period/energy and copy the
#     matching .so here instead, following the same `cp ../../bin/StRefMultCorr_*.so
#     ./bin_temp/` pattern. This is real STAR-run-period info I don't have -- check
#     with your run's centrality/refmult documentation or ask in your working group.
#   - submodule/ copies trimmed to what this repo actually has (ParticleInfo,
#     PicoDstReader_SL24y) instead of the original's PicoDstReader_SL19b/SL22b/SL22c/
#     SL23a/SL23d + a generic submodule/StRoot (StRefMultCorr source) -- none of those
#     exist in this repo. If the StRefMultCorr TODO above needs its source (not just a
#     compiled .so), you'll need to vendor in whatever submodule/StRoot equivalent
#     provides it.
#   - IMPORTANT: this maker runs under STAR's official offline software (root4star,
#     StChain, StMuDstMaker) at RCF/SDCC -- a fundamentally different runtime than the
#     rest of this repo, which ../../makefile_toggles.h deliberately hardcodes to a
#     single LOCAL-only build (`_MAC_OSX_` always defined, no _STREFMULTCORR_ toggle
#     at all -- see that file's own header comment). The `_MAC_OSX_` toggle there
#     means "not RCF" upstream, so building this maker against that toggle file as-is
#     will select the WRONG code branches in CutClass.h/HistogramUtilities.h (it skips
#     StPicoEvent*/StPicoTrack*-pointer overloads the RCF/StMuDst side may need, per
#     that file's own comment). TODO: you likely need a second, RCF-flavored
#     makefile_toggles.h (no _MAC_OSX_, plus whatever _STREFMULTCORR_ variant you
#     determine above) used only when compiling this submodule at RCF -- don't just
#     reuse ../../makefile_toggles.h unmodified for this step.

echo --------------------------------------------------
echo NOTE: this must be run at RCF/SDCC inside the STAR software environment
echo   "(after starver <library>, e.g. starver SL23c) -- not in this cloud sandbox."

mkdir -p ./bin_temp/
cp ../../bin/HistogramUtilities_cxx.so ./bin_temp/
cp ../../bin/PhysMath_cxx.so ./bin_temp/
cp ../../bin/namespaces_cxx.so ./bin_temp/
echo NOTE: ------  It is okay to not have HistoUtil+PhysMath or namespaces but Not both
cp ../../bin/CutClass_cxx.so ./bin_temp/
cp ../../bin/SetCutClass_C.so ./bin_temp/
cp ../../bin/ParticleInfo_cxx.so ./bin_temp/
# TODO: cp ../../bin/StRefMultCorr_<YOUR_VARIANT>.so ./bin_temp/  -- see header note above
cp ../../makefile_toggles.h ./
cp ../../bin/MattMcEvent_cxx.so ./bin_temp/
cp ../../bin/MattMcTrack_cxx.so ./bin_temp/
cp ../../bin/Helix_cxx.so ./bin_temp/

mkdir -p ./headers/
cp ../../headers/HistogramUtilities.h ./headers/
cp ../../headers/PhysMath.h ./headers/
cp ../../headers/CutClass.h ./headers/
cp ../../makefile_toggles.h ./headers/
cp ../../headers/MattMcEvent.h ./headers/
cp ../../headers/MattMcTrack.h ./headers/
cp ../../headers/Helix.h ./headers/

mkdir -p ./submodule/
cp -r ../../submodule/ParticleInfo ./submodule/
cp -r ../../submodule/PicoDstReader_SL24y ./submodule/
# NOTE: StMuDstMaker itself (StChain/StMuDstMaker/StMuTrack/...) comes from the STAR
# software environment via `starver`, not from this repo -- nothing to copy for that.

cons
