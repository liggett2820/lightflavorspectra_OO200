// LoadEventDisplayLibs.C -- lightflavorspectra_OO200
// Run FIRST, in the same root session as MakeEventDisplay_OO200.C, from the repo root:
//   root -l -b -q macros/LoadEventDisplayLibs.C 'macros/MakeEventDisplay_OO200.C+("filelist.list")'
// (Same split as LoadPicoBinnerLibs.C: the library must be loaded before the display
//  macro is compiled.)
void LoadEventDisplayLibs(){
  gSystem->Load("./bin/libStPicoDst_SL23c.so");
  gSystem->AddIncludePath("-I./submodule/PicoDstReader_SL23c");
}
