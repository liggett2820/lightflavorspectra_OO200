#include "TStopwatch.h"
#include<string>
#include "./headers/makefile_toggles.h"

class StMaker;
class StChain;
class StMuDstMaker;


void doMuDst(const Char_t *  a_inputFileList, const Char_t* a_outputFileName, Char_t* a_partStr, int a_nEvents = 10000000, int a_motherID = -999){
  TStopwatch*   stopWatch = new TStopwatch();
  stopWatch->Start();
  
  gSystem->Load("./bin_temp/namespaces_cxx.so");
  gSystem->Load("./bin_temp/ParticleInfo_cxx.so");
  gSystem->Load("./bin_temp/namespaces_cxx.so");
  #ifdef _SL19b_STREFMULTCORR_  
    gSystem->Load("./bin_temp/StRefMultCorr_SL19b.so");
  #endif
  #ifdef _SL21a_STREFMULTCORR_  
    gSystem->Load("./bin_temp/StRefMultCorr_SL21a.so");
  #endif
  #ifdef _SL22dev_STREFMULTCORR_  
    gSystem->Load("./bin_temp/StRefMultCorr_SL22dev.so");
  #endif
  #ifdef _7p7GeV_STREFMULTCORR_
    gSystem->Load("./bin_temp/StRefMultCorr_7p7GeV.so");
  #endif
  #ifdef _9p2_11p5_17p3GeV_STREFMULTCORR_
    gSystem->Load("./bin_temp/StRefMultCorr_9p2_11p5_17p3GeV.so");
  #endif

  gSystem->Load("./bin_temp/CutClass_cxx.so");
  gSystem->Load("./bin_temp/SetCutClass_C.so");
  gSystem->Load("./bin_temp/Helix_cxx.so");
  gSystem->Load("./bin_temp/MattMcEvent_cxx.so");
  gSystem->Load("./bin_temp/MattMcFDEvent_cxx.so");
  gSystem->Load("./bin_temp/MattMcTrack_cxx.so");
  gSystem->Load("./bin_temp/MattMcFDTrack_cxx.so");

  gROOT->Macro("loadMuDst.C");
  gSystem->Load("StMuAnalysisMaker");
  
  CutClass *cuts = new CutClass();
  SetCutClass(cuts);
  StChain* chain = new StChain();
  //StMuDstMaker(int mode, int nameMode, const char* dirName="./", const char* fileName="", const char* filter=".", int maxfiles=10, const char* name="MuDst");
  StMuDstMaker* MuDstMaker = new StMuDstMaker(0, 0, "", a_inputFileList, "MuDst", 10000); //last num is maxfiles
  MuDstMaker->SetStatus("*", 1); //Turn all of the chains on (0 for all off)
  // MuDstMaker->SetStatus("MuEvent", 1);
  // MuDstMaker->SetStatus("PrimaryVertices", 1);
  // MuDstMaker->SetStatus("PrimaryTracks", 1);
  // MuDstMaker->SetStatus("GlobalTracks", 1);
  // MuDstMaker->SetStatus("CovGlobTrack", 1);
  // MuDstMaker->SetStatus("StMuMc*", 1); 
  // MuDstMaker->SetStatus("BTof*", 1);
   
  StMuAnalysisMaker* anaMaker = new StMuAnalysisMaker();
  anaMaker->setOutputName(a_outputFileName);
  cout << " SetOutputName : " << a_outputFileName << endl;
  anaMaker->setCutClass(cuts);
  #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
  cout << "Setting Mother Feed Down ID: " << a_motherID << endl;
  anaMaker->set_mother_id(a_motherID);
  #else
  anaMaker->setPart(a_partStr);
  //anaMaker->setUseCutClassMtM0Edges(false); // because we use fit functions, you shouldn't use the same structure if you have wide bins
  //anaMaker->setEta();
  //setMtM0 will set the mTm0 binning to be the same as in SetCutClass.C
  //You can give an integer argument and it will multiply the number of bins by that number (default of 1)
  //anaMaker->setMtM0(1); // will only work if using CutClass MtM0 edges
  #endif
  cout << "Finished with StMuAnalysisMaker initialization" << endl;

  int init_status = chain->Init();
  if(init_status) cout << "ERROR: Chain->Init() returned bad status" << endl;
  cout << "Finished chain->Init();" << endl;


  int total = MuDstMaker->chain()->GetEntries();
  cout << " Total Events to Process = " << total << endl;
  if(a_nEvents>total) a_nEvents = total;
  #ifdef _MuEmbDebug_
    if(a_nEvents > 100 && total > 100) a_nEvents = 100;
  #endif

  for (int eventIndex = 0; eventIndex < a_nEvents; eventIndex++){
    if (eventIndex % 100 == 0) cout << "Working on eventNumber " << eventIndex << endl;
    chain->Clear();
    int make_status = chain->Make(eventIndex);
    if (make_status != kStOK){
      continue; //this is a continue in the example code
    }
  }

  cout << "Finished Event Loop" << endl;
  chain->Finish();
  cout << "Completed Finish() with total number of events  " << a_nEvents << endl;
  delete chain;
  
  stopWatch->Stop();
  stopWatch->Print();
}
