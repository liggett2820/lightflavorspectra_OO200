// LoadRawSpectraModifierLibs.C -- lightflavorspectra_OO200
//
// Same pattern as LoadPicoBinnerLibs.C / LoadZFitterLibs.C: run this as its own macro
// FIRST, in the SAME root invocation as RunRawSpectraModifier.C, e.g.:
//
//   root -l -q -b LoadRawSpectraModifierLibs.C 'RunRawSpectraModifier.C("pionYield.root","kaonYield.root","protonYield.root","spectra.root","efficiency.root","spectra_modified.root")'

void LoadRawSpectraModifierLibs(){
  gSystem->Load("../bin/namespaces_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/CutClass_cxx.so");
  gSystem->Load("../bin/SpectraStructure_cxx.so");
  gSystem->Load("../bin/SplineFitter_cxx.so");
  gSystem->Load("../bin/RawSpectraModifier_cxx.so");
}
