// LoadSpectraFitterLibs.C -- lightflavorspectra_OO200
//
// Same pattern as LoadPicoBinnerLibs.C / LoadZFitterLibs.C / LoadRawSpectraModifierLibs.C:
// run this as its own macro FIRST, in the SAME root invocation as RunSpectraFitter.C, e.g.:
//
//   root -l -q -b LoadSpectraFitterLibs.C 'RunSpectraFitter.C("spectra_modified.root","fit_output.root")'

void LoadSpectraFitterLibs(){
  gSystem->Load("../bin/namespaces_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/CutClass_cxx.so");
  gSystem->Load("../bin/FitUtilities_cxx.so");
  gSystem->Load("../bin/SpectraFitter_cxx.so");
}
