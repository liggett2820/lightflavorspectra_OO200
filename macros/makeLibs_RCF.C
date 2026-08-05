// makeLibs_RCF.C -- lightflavorspectra_OO200
//
// Ported from lightflavorspectra_etof/macros/makeLibs_RCF.C. This is the RCF-side build
// macro macros/makefile's basic/col27/col19p6/etc. targets already call (`root -l -q -b
// macros/makeLibs_RCF.C`) -- it produces the bin/*.so files that
// submodule/MuDstProcessEmbedding/compile.bash then copies into its own bin_temp/ before
// compiling StMuAnalysisMaker with `cons` under STAR's official software. This file was
// simply missing from this fork until now (macros/makeLibs.C, the LOCAL/interactive build
// script actually used for the rest of this repo's macros, is a different thing -- see
// that file's own header comment).
//
// WHAT CHANGED FROM THE ORIGINAL, AND WHY:
//   - Trimmed to only the files that actually exist in this repo: dropped
//     PicoQA.cxx/runByRunQA.cxx/PicoDST_MattDST_Converter.cxx (never ported here) and
//     MattMcFDEvent.cxx/MattMcFDTrack.cxx (feed-down-embedding TTree schema -- not
//     ported yet, matches the same scoping decision documented in compile.bash's header:
//     hybrid feed-down needs a separate embedding request that hasn't been made/ported
//     yet).
//   - Kept source/Helix.cxx as a SEPARATE CompileMacro call (matching the original
//     exactly), even though this repo's own macros/makeLibs.C bundles Helix.cxx into
//     namespaces.cxx's single translation unit instead (see that file's comment) --
//     don't "fix" this to match makeLibs.C's convention. submodule/MuDstProcessEmbedding/
//     compile.bash (also ported from the original) expects a standalone Helix_cxx.so to
//     exist in bin/, and the original repo already runs both this file's separate Helix
//     compile AND namespaces.cxx's bundled one side by side in real production use --
//     don't assume that's a bug just because it looks redundant.
//   - StRefMultCorr: none of the original's #ifdef _SL19b_STREFMULTCORR_ / _SL21a_.../
//     _SL22dev_.../_7p7GeV_.../_9p2_11p5_17p3GeV_... variants apply to O+O 200 GeV --
//     all commented out below. TODO: same open question flagged in compile.bash's header
//     -- determine which StRefMultCorr class (if any) applies to your run period and add
//     a matching #ifdef block + gSystem->Load() line here, plus the matching
//     ./bin/StRefMultCorr_<variant>.so build rule in macros/makefile (see that file's
//     existing ./bin/StRefMultCorr_SL19b.so-style rules for the pattern to copy).
//
// Usage (run from the REPO ROOT, inside the STAR software environment at RCF, after
// `starver <library>`, e.g. `starver SL23c`):
//   root -l -q -b macros/makeLibs_RCF.C

#include "../makefile_toggles.h"

void makeLibs_RCF(TString opt=""){

  if (opt == "clean"){
    gSystem->Exec("rm -f bin/*");
    return;
  }

  if (gSystem->OpenDirectory("bin") == 0)
    gSystem->MakeDirectory("bin");

  gSystem->SetBuildDir("$PWD/bin/",true);
  gSystem->SetIncludePath(TString::Format("-I$ROOTSYS/include -I%s/headers",gSystem->pwd()));

  gSystem->Load("./bin/libStPicoDst.so");

  gSystem->CompileMacro("submodule/ParticleInfo/ParticleInfo/ParticleInfo.cxx","gk");

  // TODO: StRefMultCorr -- see header note above. Example shape (uncomment/adapt once
  // you know which variant applies and have added the matching makefile build rule):
  // #ifdef _SL23c_STREFMULTCORR_
  //   gSystem->Load("bin/StRefMultCorr_SL23c.so");
  // #endif

  gSystem->CompileMacro("source/namespaces.cxx","gk");
  gSystem->CompileMacro("source/CutClass.cxx","gk");
  gSystem->CompileMacro("source/CalibrationClass.cxx","gk");

  gSystem->CompileMacro("source/PicoBinner.cxx","gk");
  gSystem->CompileMacro("source/MattMcEvent.cxx","gk");
  gSystem->CompileMacro("source/MattMcTrack.cxx","gk");
  gSystem->CompileMacro("source/Helix.cxx","gk");

  gSystem->CompileMacro("macros/SetCutClass.C","gk");
}
