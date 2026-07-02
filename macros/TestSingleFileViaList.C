// TestSingleFileViaList.C -- diagnostic only, not part of the analysis.
//
// Reads ONE file through the exact same ".list"-triggered code path that
// StPicoDstReader::Init() uses for the original all-files-in-one-chain design (as
// opposed to the bare-single-file-path branch our per-file redesign uses). If this
// shows the same corruption as TestReaderCycle.C did on the same file, that proves the
// corruption is inherent to the file's bytes on disk, not an artifact of how we're
// constructing readers.
//
// Usage: first make a one-line list file, e.g.:
//   echo "/Volumes/IDISKK/OO_fastoffline_fxt-blue/st_physics_27035059_raw_1200011.picoDst.root" > /tmp/onefile.list
// Then, from the REPO ROOT:
//   root -l -q -b macros/TestSingleFileViaList.C'("/tmp/onefile.list")'

#include "../submodule/PicoDstReader_SL23a/StPicoDstReader.h"
#include "../submodule/PicoDstReader_SL23a/StPicoDst.h"
#include "../submodule/PicoDstReader_SL23a/StPicoEvent.h"
#include <iostream>
#include <string>
#include "TSystem.h"

void TestSingleFileViaList(std::string oneFileListPath){

  gSystem->Load("./bin/libStPicoDst.so");

  StPicoDstReader* reader = new StPicoDstReader(oneFileListPath.c_str());
  reader->Init();
  StPicoDst* dst = reader->picoDst();
  Long64_t n = reader->chain()->GetEntries();
  std::cout << "Entries: " << n << std::endl;

  int nGood = 0, nNull = 0, nFailed = 0;
  for(Long64_t i = 0; i < n; i++){
    Bool_t ok = reader->readPicoEvent(i);
    if(!ok){ std::cout << "read failed at " << i << std::endl; nFailed++; break; }
    StPicoEvent* ev = dst->event();
    if(!ev){ nNull++; continue; }
    nGood++;
    if(i % 100 == 0) std::cout << "ok at event " << i << std::endl;
  }

  std::cout << "TestSingleFileViaList finished. good=" << nGood << " null=" << nNull
             << " failed=" << nFailed << std::endl;

  reader->Finish();
  delete reader;
}
