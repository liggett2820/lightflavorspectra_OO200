// LoadZFitterLibs.C -- lightflavorspectra_OO200
//
// Same pattern as LoadPicoBinnerLibs.C: run this as its own macro FIRST, in the SAME
// root invocation as RunZFitter.C, so ZFitter's (and its dependencies') compiled
// symbols are already loaded and resolvable by the time RunZFitter.C is parsed, e.g.:
//
//   root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(0,0,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'

void LoadZFitterLibs(){
  gSystem->Load("../bin/namespaces_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/CutClass_cxx.so");
  gSystem->Load("../bin/GausMixModel_cxx.so");
  gSystem->Load("../bin/Histo2D_cxx.so");
  gSystem->Load("../bin/ZFitter_cxx.so");
  gSystem->Load("../bin/SetCutClass_C.so");
}
