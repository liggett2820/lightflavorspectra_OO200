// LoadEfficiencyFitterLibs.C -- lightflavorspectra_OO200
//
// Same pattern as LoadRawSpectraModifierLibs.C / LoadZFitterLibs.C: run this as its
// own macro FIRST, in the SAME root invocation as RunEfficiencyFitter.C, e.g.:
//
//   root -l -q -b LoadEfficiencyFitterLibs.C 'RunEfficiencyFitter.C("pionYield.root","kaonYield.root","efficiency_out.root")'

void LoadEfficiencyFitterLibs(){
  gSystem->Load("../bin/namespaces_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/SplineFitter_cxx.so");
  gSystem->Load("../bin/GraphFitter_cxx.so");
  gSystem->Load("../bin/EfficiencyFitter_cxx.so");
}
