// LoadCutFlowLibs.C -- lightflavorspectra_OO200
// Loads the SL24y PicoDst reader library for macros/MakeEventCutFlow_OO200.C (the O+O
// production picoDsts are P24iy/SL24y). Run FIRST, in the same root session, from the repo root:
//   root -l -b -q macros/LoadCutFlowLibs.C 'macros/MakeEventCutFlow_OO200.C+("files.list","cutflow_part0.root")'
// bin/libStPicoDst.so comes from the README's "One-time build setup", step 1
// (submodule/PicoDstReader_SL24y). LoadEventDisplayLibs.C loads the SL23c reader instead,
// which is for the embedding (P23ic/SL23c) files, not the production picoDsts.
void LoadCutFlowLibs(){
  if(gSystem->Load("./bin/libStPicoDst.so") < 0)
    printf("LoadCutFlowLibs: could not load ./bin/libStPicoDst.so -- build it first (README, One-time build setup, step 1)\n");
  gSystem->AddIncludePath("-I./submodule/PicoDstReader_SL24y");
}
