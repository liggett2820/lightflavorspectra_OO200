// LoadPicoBinnerLibs.C -- lightflavorspectra_OO200
//
// Separated out from RunPicoBinner.C on purpose. gSystem->Load() calls and calls to
// symbols from those libraries can't safely share one function body/translation unit
// with cling (it JIT-compiles and link-resolves a whole unit before any statement in
// it actually runs). Run this as its own macro FIRST, in the SAME root invocation as
// RunPicoBinner.C, e.g.:
//
//   root -l -q -b LoadPicoBinnerLibs.C 'RunPicoBinner.C("file.list","out.root",9,"PION","true","false","false","false")'
//
// Each file argument to `root` is processed as its own separate step within the same
// session, so by the time RunPicoBinner.C is parsed, these libraries are already
// loaded and their symbols already resolvable.

void LoadPicoBinnerLibs(){
  gSystem->Load("../bin/namespaces_cxx.so");
  gSystem->Load("../bin/ParticleInfo_cxx.so");
  gSystem->Load("../bin/CalibrationClass_cxx.so");
  gSystem->Load("../bin/CutClass_cxx.so");
  gSystem->Load("../bin/PicoBinner_cxx.so");
  gSystem->Load("../bin/SetCutClass_C.so");
}
