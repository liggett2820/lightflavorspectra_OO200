// macros/GetRefMultHisto.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS (2026-07-23): Andrew asked for a macro to plot refMult straight from
// raw fastoffline PicoDst files -- before PicoBinner has been run on this dataset, so
// there's no yieldHistos_*.root file yet for macros/PresentEventQA.C to read a "refMult"
// histogram out of (that macro's own refMult panel reads a histogram PicoBinner already
// booked and filled; this macro fills one from scratch instead).
//
// DIFFERENT DATASET (clarified 2026-07-24): this macro targets the 2026 O+O fastoffline
// sample, NOT the 2021 O+O data that the rest of this repo's pipeline
// (SetCutClass.C/PicoBinner.cxx/the yieldHistos_OOGeV_* files) has been analyzing. That
// means its trigger ID is its own, independent choice -- it is NOT supposed to track
// whatever SetCutClass.C configures for the 2021 dataset, and earlier revisions of this
// file that said otherwise were wrong. Default is trigger 6 for the 2026 sample; do not
// resync this to SetCutClass.C's trigger ID going forward.
//
// Deliberately minimal and dependency-free: no CutClass/SetCutClass, no pileup, vertex,
// or track cuts, no centrality binning -- just event->refMult() per event, straight from
// the picoDst, same scope as the reference standalone reader this was modeled on
// (OO_analysis/getrefmulthisto.cpp). This is meant as a fast raw look at the sample
// BEFORE any of that machinery is wired up, not a replacement for PicoBinner's own
// (cut-aware) refMult histogram.
//
// TRIGGER CUT (added 2026-07-23, corrected 2026-07-24): the one exception to "no cuts"
// above. Applies event->isTrigger(a_triggerId) (StPicoEvent.h line 68) directly against
// this 2026-dataset picoDst, default a_triggerId=6 (see DIFFERENT DATASET note above --
// this is unrelated to macros/SetCutClass.C's trigger ID, which is for the separate 2021
// dataset and its own PicoBinner.cxx pipeline). Pass a_triggerId<0 to disable (plot every
// event regardless of trigger, e.g. to compare against the cut version).
//
// FILE-BY-FILE READING (NOT copied from the OO_analysis reference macro, which builds
// ONE StPicoDstReader/TChain directly from the whole filelist): source/PicoBinner.cxx
// (this repo's own real analysis driver) documents a reproducible segfault from exactly
// that single-TChain-over-everything pattern -- it crashes at the exact event index
// corresponding to a file-boundary crossing, deep inside TBranchElement::GetEntry, on
// the real fastoffline filelist (see PicoBinner.cxx's "FILE-BY-FILE READING" comment,
// ~line 1161, for the full misdiagnosis-then-root-cause story: it was traced to how
// ROOT's TChain/TBranchElement machinery -- or the STAR-provided StPicoDstReader/
// StPicoDst dictionary code -- handles a file-transition boundary in this environment,
// NOT the SetStatus branch-enable pattern originally suspected). PicoBinner.cxx's fix,
// replicated here: parse the filelist by hand, and construct a FRESH StPicoDstReader for
// each individual file, fully Finish()+deleted before moving to the next -- this
// guarantees a file-boundary transition inside a single TChain instance never happens,
// sidestepping the bug entirely without changing which events get read.
//
// FILELIST FORMAT: one picoDst.root path per line (optionally followed by whitespace and
// a second token, e.g. an event count some external formatters append -- that second
// token is discarded here, matching PicoBinner.cxx's own parsing). Each candidate line is
// required to contain ".picoDst.root" and to open cleanly (TFile::Open succeeds, not a
// zombie, GetNkeys()>0) before it's added to the read list -- a corrupt or missing file
// in the list is silently skipped (reported in the "Found N valid picoDst files" count),
// not a hard error.
//
// Booked 1000 bins, 0-1000 -- same convention as PicoBinner.cxx's own refMult histogram
// (source/PicoBinner.cxx ~line 396), for consistency if you later want to overlay this
// against a PicoBinner-produced one. a_xMaxDisplay only clips the drawn range (like
// macros/PresentEventQA.C's refMult panel does at 140) -- it doesn't rebin or drop any
// underlying entries, and the full-range histogram is still what's written to
// a_outFile.
//
// USAGE (run from within macros/ or anywhere -- no ../bin/ compiled-library dependency,
// only the PicoDstReader_SL24y submodule's own headers/dictionary):
//   root -l -q 'GetRefMultHisto.C("/path/to/filelist_OO_fastoffline.list")'
//   root -l -q 'GetRefMultHisto.C("/path/to/filelist.list","./refmult_qa","refMult_fastoffline.root",200)'  // clip display to refMult<200
//   root -l -q 'GetRefMultHisto.C("/path/to/filelist.list",".","refMult_fastoffline.root",200,-1)'  // trigger cut OFF, for comparison

#include "../submodule/PicoDstReader_SL24y/StPicoDstReader.h"
#include "../submodule/PicoDstReader_SL24y/StPicoDst.h"
#include "../submodule/PicoDstReader_SL24y/StPicoEvent.h"

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

void GetRefMultHisto(string a_filelist,
                      string a_outDir = ".",
                      string a_outFile = "refMult_fastoffline.root",
                      double a_xMaxDisplay = 200,
                      int a_triggerId = 6){ // 2026 O+O fastoffline dataset's own trigger -- unrelated to SetCutClass.C's (2021 dataset). See DIFFERENT DATASET / TRIGGER CUT header comments. Pass <0 to disable

  gSystem->mkdir(a_outDir.c_str(), true);

  // ---- Parse the filelist ourselves -- see FILE-BY-FILE READING header comment above
  // for why this can't just be handed to a single StPicoDstReader/TChain. Logic mirrors
  // source/PicoBinner.cxx's own filelist parsing exactly (~line 1186-1204).
  cout << "Parsing filelist: " << a_filelist << endl;
  vector<string> picoFileList;
  {
    std::ifstream inputStream(a_filelist.c_str());
    if(!inputStream){
      cout << "GetRefMultHisto :: Error - cannot open filelist " << a_filelist << endl;
      return;
    }
    std::string line;
    while(std::getline(inputStream, line)){
      // Some external formatters emit "file NumEvents" -- keep only the first token.
      size_t pos = line.find_first_of(" ");
      if(pos != std::string::npos) line.erase(pos, line.length()-pos);
      if(line.find(".picoDst.root") == std::string::npos) continue;
      TFile* ftmp = TFile::Open(line.c_str());
      if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) picoFileList.push_back(line);
      if(ftmp) ftmp->Close();
    }
  }
  cout << "Found " << picoFileList.size() << " valid picoDst files in filelist." << endl;
  if(picoFileList.empty()){
    cout << "GetRefMultHisto :: no valid picoDst files found -- nothing to do." << endl;
    return;
  }

  TH1I* refMult = new TH1I("refMult","Multiplicity;refMult;Number of Events",1000,0,1000);

  int modNum = 1;
  Long64_t globalEventIndex = 0; // running count across all files, purely for the print cadence
  Long64_t nTriggerRejected = 0; // events skipped by the trigger cut -- see TRIGGER CUT header comment

  // ---- FILE-BY-FILE EVENT LOOP -- see header comment: never let a single
  // StPicoDstReader/TChain span more than one file.
  for(size_t fileIndex = 0; fileIndex < picoFileList.size(); fileIndex++){

    StPicoDstReader* picoReader = new StPicoDstReader(picoFileList[fileIndex].c_str());
    picoReader->Init();
    picoReader->SetStatus("*",0);
    picoReader->SetStatus("Event",1); // refMult lives on StPicoEvent -- no track branches needed
    StPicoDst* dst = picoReader->picoDst();
    Long64_t eventsInThisFile = (Long64_t) picoReader->chain()->GetEntries();
    cout << "File " << fileIndex+1 << "/" << picoFileList.size() << " (" << picoFileList[fileIndex]
         << "): " << eventsInThisFile << " event(s)" << endl;

    for(Long64_t localEventIndex = 0; localEventIndex < eventsInThisFile; localEventIndex++, globalEventIndex++){
      if(globalEventIndex % modNum == 0){
        cout << "Working on event " << globalEventIndex << endl;
        if(globalEventIndex == 10*modNum) modNum *= 10;
      }

      Bool_t readEvent = picoReader->readPicoEvent(localEventIndex);
      if(!readEvent){
        std::cerr << "GetRefMultHisto :: readPicoEvent failed at file " << fileIndex
                   << ", local event " << localEventIndex << " -- stopping this file early." << std::endl;
        break;
      }

      StPicoEvent* event = dst->event();
      if(!event){
        cerr << "GetRefMultHisto :: event " << globalEventIndex << " doesn't exist?!" << endl;
        continue;
      }

      // TRIGGER CUT -- see this file's header comment for why this is applied here
      // rather than relying on CutClass (whose configured trigger requirement never
      // actually gets checked in PicoBinner.cxx's real event loop).
      if(a_triggerId >= 0 && !event->isTrigger((unsigned int) a_triggerId)){
        nTriggerRejected++;
        continue;
      }

      refMult->Fill(event->refMult());
    }

    picoReader->Finish();
    delete picoReader;
  }
  cout << "Finished looping over all " << globalEventIndex << " events" << endl;
  if(a_triggerId >= 0){
    cout << "Trigger cut (isTrigger(" << a_triggerId << ")): rejected " << nTriggerRejected
         << " / " << globalEventIndex << " event(s); " << refMult->GetEntries()
         << " event(s) filled into refMult." << endl;
  }else{
    cout << "Trigger cut disabled (a_triggerId<0) -- all " << refMult->GetEntries()
         << " event(s) filled into refMult." << endl;
  }

  // ---- Write the raw histogram (full 0-1000 range, unclipped) ----
  string outFilePath = a_outDir + "/" + a_outFile;
  TFile* outFile = new TFile(outFilePath.c_str(),"RECREATE");
  outFile->cd();
  refMult->Write();
  outFile->Close();
  cout << "Wrote " << outFilePath << endl;

  // ---- Draw the plot ----
  gStyle->SetOptStat(1);
  TCanvas* c = new TCanvas("GetRefMultHisto","refMult",900,700);
  c->SetLogy();
  refMult->SetTitle(a_triggerId >= 0 ? Form("refMult, fastoffline, trigger %d",a_triggerId) : "refMult, fastoffline, no trigger cut");
  refMult->GetXaxis()->SetRangeUser(0, a_xMaxDisplay);
  refMult->Draw();

  string plotPath = a_outDir + "/refMult_fastoffline.png";
  c->SaveAs(plotPath.c_str());
  cout << "Wrote " << plotPath << endl;
}
