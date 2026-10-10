// CheckDetectorHits_OO200.C -- lightflavorspectra_OO200
//
// Quick check of which detector hit collections are filled in the P24iy O+O picoDsts, split by
// whether the event fired trigger 860003 (minbias-bbc_etof). Written to answer "is the EPD read
// out for 860003?" after a bare-ROOT TTree::GetEntries("EpdHit_>0") check gave 0 for EPD, BTOF
// and eTOF alike (not trustworthy without the StPicoDst dictionaries). Uses the compiled reader.
//
// Usage (RCF, repo root, inside the SL7 container, SL24y). Argument: one picoDst path/URL or a .list:
//   root -l -b -q macros/LoadCutFlowLibs.C 'macros/CheckDetectorHits_OO200.C+("root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/production_OO_200GeV_2021/ReversedFullField/P24iy/2021/133/22133043/st_physics_22133043_raw_0000007.picoDst.root")'

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "TString.h"
#include "TChain.h"
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"

void CheckDetectorHits_OO200(const char* input, int maxFiles = 5){
  std::vector<std::string> files;
  TString in(input);
  if(in.EndsWith(".list")){
    std::ifstream fl(input); std::string line;
    while(std::getline(fl, line)){
      std::string f = line.substr(0, line.find_first_of(" \t"));
      if(f.find(".picoDst.root") != std::string::npos) files.push_back(f);
    }
  } else files.push_back(input);
  if(maxFiles > 0 && (int)files.size() > maxFiles) files.resize(maxFiles);

  const char* name[5] = {"Track", "BTofHit", "ETofHit", "EpdHit", "BbcHit"};
  // [0] all events, [1] 860003 events
  long long nEv[2] = {0, 0}, nWith[2][5], sum[2][5];
  for(int t = 0; t < 2; t++) for(int d = 0; d < 5; d++){ nWith[t][d] = 0; sum[t][d] = 0; }

  for(size_t f = 0; f < files.size(); f++){
    StPicoDstReader* reader = new StPicoDstReader(files[f].c_str());
    reader->Init();
    reader->SetStatus("*", 0);
    reader->SetStatus("Event", 1);
    for(int d = 0; d < 5; d++) reader->SetStatus(name[d], 1);
    Long64_t n = reader->chain()->GetEntries();
    for(Long64_t i = 0; i < n; i++){
      if(!reader->readPicoEvent(i)) continue;
      StPicoDst* dst = reader->picoDst();
      StPicoEvent* ev = dst->event();
      if(!ev) continue;
      unsigned int cnt[5] = {dst->numberOfTracks(), dst->numberOfBTofHits(), dst->numberOfETofHits(),
                             dst->numberOfEpdHits(), dst->numberOfBbcHits()};
      for(int t = 0; t < 2; t++){
        if(t == 1 && !ev->isTrigger(860003)) continue;
        nEv[t]++;
        for(int d = 0; d < 5; d++){ if(cnt[d] > 0) nWith[t][d]++; sum[t][d] += cnt[d]; }
      }
    }
    reader->Finish();
    delete reader;
  }

  std::cout << "\nFiles read: " << files.size() << std::endl;
  for(int t = 0; t < 2; t++){
    std::cout << (t == 0 ? "All events: " : "Trigger 860003: ") << nEv[t] << std::endl;
    for(int d = 0; d < 5; d++)
      std::cout << Form("  %-8s events with hits: %8lld (%5.1f%%)   mean hits/event: %8.2f", name[d], nWith[t][d],
                        nEv[t] ? 100.*nWith[t][d]/nEv[t] : 0., nEv[t] ? double(sum[t][d])/nEv[t] : 0.) << std::endl;
  }
}
