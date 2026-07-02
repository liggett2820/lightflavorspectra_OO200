// makeLibs.C -- lightflavorspectra_OO200
//
// Simplified build macro. Compiles PicoBinner and ZFitter plus their dependencies.
// As RawSpectraModifier/SpectraFitter get ported, add their CompileMacro lines here
// the same way the original repo's macros/makeLibs.C and macros/makeLibs_RCF.C do.
//
// NOTE ON ZFitter'S BUILD: the original repo built ZFitter via a separate real
// Makefile rather than ACLiC (see macros/makeLibs.C there -- its CompileMacro line for
// ZFitter.cxx is commented out). This repo instead compiles it with gSystem->
// CompileMacro like everything else, for consistency with how PicoBinner is already
// built here -- ACLiC handles ZFitter.cxx's size fine, just slower to build (it's a
// ~35,000 line single translation unit once its #include chain is counted). Expect
// this step to take several minutes the first time.
//
// Run this from the REPO ROOT (not from macros/), e.g.:
//   root -l -q -b macros/makeLibs.C
//
// Prerequisite: submodule/PicoDstReader_SL24y/libStPicoDst.so must already be built
// (see the README / chat instructions for the one-time `make` step in that directory)
// and copied into ./bin/.

#include "../makefile_toggles.h"

void makeLibs(TString opt = ""){

  if (opt == "clean"){
    gSystem->Exec("rm -f bin/*");
    return;
  }

  if (gSystem->OpenDirectory("bin") == 0)
    gSystem->MakeDirectory("bin");

  gSystem->SetBuildDir("$PWD/bin/", true);
  gSystem->SetIncludePath(TString::Format("-I$ROOTSYS/include -I%s/headers", gSystem->pwd()));

  // Built by submodule/PicoDstReader_SL24y's own Makefile, then copied here.
  gSystem->Load("./bin/libStPicoDst.so");

  // namespaces.cxx pulls in Helix.cxx, PhysMath.cxx, HistogramUtilities.cxx, and
  // PhysMathHistogramUtilities.cxx as a single translation unit -- this matches the
  // original repo's build exactly (those files are NOT compiled separately there
  // either; see source/namespaces.cxx).
  gSystem->CompileMacro("source/namespaces.cxx", "gk");

  gSystem->CompileMacro("submodule/ParticleInfo/ParticleInfo/ParticleInfo.cxx", "gk");

  gSystem->CompileMacro("source/CutClass.cxx", "gk");
  gSystem->CompileMacro("source/CalibrationClass.cxx", "gk");
  gSystem->CompileMacro("source/PicoBinner.cxx", "gk");
  gSystem->CompileMacro("macros/SetCutClass.C", "gk");

  // ZFitter dependencies -- must be compiled before ZFitter.cxx itself, since
  // headers/ZFitter.h includes GausMixModel.h and Histo2D.h.
  gSystem->CompileMacro("source/GausMixModel.cxx", "gk");
  gSystem->CompileMacro("source/Histo2D.cxx", "gk");
  gSystem->CompileMacro("source/ZFitter.cxx", "gk");

  // RawSpectraModifier dependencies -- must be compiled before RawSpectraModifier.cxx
  // itself, since headers/RawSpectraModifier.h includes SpectraStructure.h and
  // SplineFitter.h.
  gSystem->CompileMacro("source/SpectraStructure.cxx", "gk");
  gSystem->CompileMacro("source/SplineFitter.cxx", "gk");
  gSystem->CompileMacro("source/RawSpectraModifier.cxx", "gk");

  // SpectraFitter dependencies -- must be compiled before SpectraFitter.cxx itself,
  // since headers/SpectraFitter.h includes FitUtilities.h (in addition to CutClass.h,
  // PhysMath.h, ParticleInfo.h, all already compiled above).
  gSystem->CompileMacro("source/FitUtilities.cxx", "gk");
  gSystem->CompileMacro("source/SpectraFitter.cxx", "gk");
}
