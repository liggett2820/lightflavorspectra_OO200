// LoadDrawingLibs.C -- lightflavorspectra_OO200
//
// Same pattern as LoadPicoBinnerLibs.C/LoadZFitterLibs.C/LoadRawSpectraModifierLibs.C/
// LoadSpectraFitterLibs.C: run this as its own macro FIRST, in the SAME root invocation
// as DrawSpectraAndDnDy.C, so its (much smaller) dependency -- just ParticleInfo, used
// only to reconstruct the exact object-naming convention SpectraFitter::writeOutputs()
// used -- is already loaded and resolvable by the time DrawSpectraAndDnDy.C is parsed:
//
//   root -l -q -b LoadDrawingLibs.C 'DrawSpectraAndDnDy.C("fit_output.root")'
//
// Deliberately does NOT load CutClass/SpectraFitter/etc. -- this macro only reads
// already-written TGraphErrors/TH1D objects back out of SpectraFitter's output file,
// it doesn't need the fitting machinery itself.

void LoadDrawingLibs(){
  gSystem->Load("../bin/namespaces_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
}
