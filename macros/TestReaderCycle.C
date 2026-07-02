// TestReaderCycle.C -- diagnostic only, not part of the analysis.
//
// Isolates one question: does repeatedly constructing a fresh StPicoDstReader per
// file, reading all its events, then Finish()+delete-ing it before moving to the next
// file, crash on its own -- independent of any PicoBinner histogram/analysis code?
//
// This only loads libStPicoDst.so (must already be built into ./bin/ by the
// submodule's own Makefile -- same prerequisite as makeLibs.C) and touches a couple
// of StPicoEvent/StPicoDst fields per event, just enough to force real memory reads.
//
// Run from the REPO ROOT:
//   root -l -q -b macros/TestReaderCycle.C'("/Volumes/IDISKK/filelist_OO_fastoffline.list")'

#include "../submodule/PicoDstReader_SL23a/StPicoDstReader.h"
#include "../submodule/PicoDstReader_SL23a/StPicoDst.h"
#include "../submodule/PicoDstReader_SL23a/StPicoEvent.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TSystem.h"

void TestReaderCycle(std::string filelist){

  gSystem->Load("./bin/libStPicoDst.so");

  std::vector<std::string> fileList;
  std::ifstream inputStream(filelist.c_str());
  if(!inputStream){ std::cout << "cannot open " << filelist << std::endl; return; }
  std::string line;
  while(std::getline(inputStream, line)){
    size_t pos = line.find_first_of(" ");
    if(pos != std::string::npos) line.erase(pos, line.length()-pos);
    if(line.find(".picoDst.root") == std::string::npos) continue;
    TFile* ftmp = TFile::Open(line.c_str());
    if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) fileList.push_back(line);
    if(ftmp) ftmp->Close();
  }
  std::cout << "Found " << fileList.size() << " valid files." << std::endl;

  Long64_t globalEventIndex = 0;
  int modNum = 1;

  // ONE reader for the whole run -- StPicoDst's picoArrays pointer is a process-wide
  // static (see StPicoDst.h), so constructing/destroying multiple StPicoDstReader
  // instances corrupts that shared state after a few cycles. SwitchToFile() resets the
  // internal TChain and points it at the next single file while reusing the SAME
  // TClonesArrays/branch addresses this one reader was built with.
  StPicoDstReader* reader = new StPicoDstReader(fileList[0].c_str());
  StPicoDst* dst = reader->picoDst();

  for(size_t fileIndex = 0; fileIndex < fileList.size(); fileIndex++){
    reader->SwitchToFile(fileList[fileIndex].c_str());
    Long64_t nEvents = (Long64_t) reader->chain()->GetEntries();

    for(Long64_t localIdx = 0; localIdx < nEvents; localIdx++, globalEventIndex++){
      int eventIndex = (int) globalEventIndex;
      if(eventIndex % modNum == 0){
        std::cout << "Working on event " << eventIndex << "  (file #" << fileIndex
                   << ": " << fileList[fileIndex] << ")" << std::endl;
        if(eventIndex == 10*modNum) modNum *= 10;
      }
      Bool_t ok = reader->readPicoEvent(localIdx);
      if(!ok){ std::cout << "  read failed at event " << eventIndex << std::endl; break; }
      StPicoEvent* ev = dst->event();
      if(!ev){ std::cout << "  null event at " << eventIndex << std::endl; continue; }
      // touch real memory so this isn't optimized away / doesn't just skip past bad reads
      volatile double vz = ev->primaryVertex().Z();
      volatile int nTracks = dst->numberOfTracks();
      (void)vz; (void)nTracks;
    }
  }

  reader->Finish();
  delete reader;

  std::cout << "TestReaderCycle finished cleanly. Total events read: " << globalEventIndex << std::endl;
}
