// macros/CheckZdcNeutronCapture_OO200.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS (2026-08-28): Andrew asked whether any neutrons were captured in the
// ZDC. StPicoEvent exposes ZdcSumAdcEast()/ZdcSumAdcWest() (summed east/west ZDC ADC per
// event) but nothing in this repo reads them yet -- this is a first, uncalibrated look:
// per-event ZDC ADC distributions plus a plain nonzero-ADC event count, straight from the
// raw picoDst, same "before PicoBinner has been run on this dataset" scope as
// GetRefMultHisto.C (which this macro is modeled on almost line-for-line).
//
// WHAT "NONZERO ADC" DOES AND DOESN'T TELL YOU: this is a naive proxy, not a calibrated
// neutron count. ZdcSumAdcEast()/West() are backed by UShort_t fields (StPicoEvent.h
// mZdcSumAdcEast/West, both unsigned) -- an unsigned type can't represent a
// pedestal-subtracted baseline that dips below zero, so it's NOT established here whether
// "0" means "genuinely zero signal" or "at/below whatever pedestal this production already
// baked in". No pedestal value or single-neutron peak position is asserted anywhere in
// this macro -- it just reports how many events have ADC>0 on each arm and plots the raw
// distributions so you can look at the peak structure yourself. Do not treat "N events
// nonzero" as a calibrated neutron multiplicity without checking that assumption.
//
// SAME DATASET / TRIGGER CUT AS GetRefMultHisto.C: this targets the 2026 O+O fastoffline
// sample, same trigger ID default (6) and same "pass <0 to disable" convention -- see that
// macro's header for why this trigger ID is independent of SetCutClass.C's (2021 dataset).
//
// FILE-BY-FILE READING (copied from GetRefMultHisto.C, not written fresh): see that
// macro's "FILE-BY-FILE READING" header comment for the full story -- a single
// StPicoDstReader/TChain spanning this filelist's file boundaries segfaults deep in
// TBranchElement::GetEntry (root-caused in source/PicoBinner.cxx). Fix replicated here
// unchanged: hand-parsed filelist, a fresh StPicoDstReader per file, Finish()+deleted
// before the next file -- never one TChain across file boundaries.
//
// BINNING: a_adcMax and a_adcBinWidth (both guesses, no calibration was available to base
// them on) set the histogram range and bin width independently -- default 0-4000 at 20
// ADC/bin. ROOT still records anything above a_adcMax in the overflow bin rather than
// silently dropping it, and this macro explicitly checks and reports overflow counts after
// the loop, so nothing is lost quietly -- if overflow is nonzero, rerun with a_adcMax
// raised. See a_adcBinWidth's own doc comment (next to the function signature below) for
// why you may want a much finer bin width than the default for peak-finding.
//
// USAGE (run from the REPO ROOT; the "+" IS REQUIRED -- see the note above the function
// below for why plain interpretation fails in this environment):
//
// 2026-08-28, third attempt: ACLiC compiles fine now (no more libc++ ABI symbol failures --
// confirms ACLiC is the right approach), but its own generated link command never mentions
// libStPicoDst.so, so StPicoDstReader/StPicoEvent symbols fail at hard link time. Fix:
// TSystem::AddLinkedLibs() -- the real API TSystem::CompileMacro() (what "+" calls) consults
// to build its link line -- called BEFORE the ACLiC compile, in the SAME root session, via a
// second -e (root runs multiple -e statements and the trailing macro in one session, in
// order). That link fix alone still wasn't enough, though: AddLinkedLibs finds the symbols
// at BUILD time via -L, but the resulting .so only records the bare name "libStPicoDst.so"
// as a RUNTIME dependency, and macOS's dyld doesn't remember -L paths -- it only searches a
// fixed set of standard locations plus DYLD_LIBRARY_PATH. CONFIRMED WORKING 2026-08-28 (real
// run over the full fastoffline filelist, 3017062 events, completed cleanly) -- all three
// pieces together:
//   export DYLD_LIBRARY_PATH="/Users/aliggett/repos/lightflavorspectra_OO200/bin:$DYLD_LIBRARY_PATH"
//   root -l -q -e 'gSystem->AddLinkedLibs("-L/Users/aliggett/repos/lightflavorspectra_OO200/bin -lStPicoDst")' 'macros/CheckZdcNeutronCapture_OO200.C+("/Volumes/IDISKK/filelist_OO_fastoffline.list")'
//   root -l -q -e 'gSystem->AddLinkedLibs("-L/Users/aliggett/repos/lightflavorspectra_OO200/bin -lStPicoDst")' 'macros/CheckZdcNeutronCapture_OO200.C+("/path/to/filelist.list",".","zdc_qa.root",4000,-1)'  // trigger cut OFF
//   root -l -q -e 'gSystem->AddLinkedLibs("-L/Users/aliggett/repos/lightflavorspectra_OO200/bin -lStPicoDst")' 'macros/CheckZdcNeutronCapture_OO200.C+("/Volumes/IDISKK/filelist_OO_fastoffline.list",".","zdc_qa_zoomed.root",300,6,2)'  // zoomed peak search: 0-300 ADC at 2 ADC/bin

#include "../submodule/PicoDstReader_SL24y/StPicoDstReader.h"
#include "../submodule/PicoDstReader_SL24y/StPicoDst.h"
#include "../submodule/PicoDstReader_SL24y/StPicoEvent.h"

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

// RUN THIS MACRO WITH ACLiC (the "+" suffix -- see USAGE above), not plain interpretation.
// 2026-08-28, two dead ends on plain "root -l -q 'file.C(args)'" interpretation before
// landing here:
//   1) gSystem->Load() as the function's first statement: unresolved-symbol linker errors
//      on every StPicoDstReader/StPicoEvent/StPicoDst method -- cling compiles the whole
//      function body as one unit when it parses the file, before that embedded Load()
//      statement has actually run.
//   2) gSystem->Load() moved to true top level (outside the function, so cling runs it
//      immediately as it reads the file): got past (1), but hit a much bigger wall --
//      cling's JIT failed to materialize a huge batch of std::string/std::vector symbols
//      (distinctive "B8ne190102"-tagged libc++ ABI) alongside the StPico ones, consistent
//      with cling (built via this Mac's conda `root`) and libStPicoDst.so (built with the
//      plain system compiler) using ABI-incompatible libc++ builds -- not something fixable
//      from inside this macro.
// ACLiC sidesteps both: it compiles this whole file for real via the system compiler (the
// same one that presumably built libStPicoDst.so) into an actual .so, rather than having
// cling JIT/interpret it -- no whole-function-as-one-unit timing issue, and no cling-vs-
// precompiled-library libc++ ABI mismatch. Its build command links with
// "-Wl,-undefined,dynamic_lookup" (lazy, runtime symbol resolution), which is exactly why
// gSystem->Load() belongs back INSIDE the function as its first statement for ACLiC: real
// compiled code executes it sequentially at runtime, before the StPico calls after it run
// -- unlike cling interpretation, where the whole function is JIT'd as one unit up front.
// Same prerequisite TestReaderCycle.C/TestSingleFileViaList.C document: libStPicoDst.so
// must already be built into ./bin/ by the submodule's own Makefile, and this must be run
// from the REPO ROOT (relative path below is relative to the cwd ROOT was started in).
void CheckZdcNeutronCapture_OO200(string a_filelist = "/Volumes/IDISKK/filelist_OO_fastoffline.list",
                                   string a_outDir = ".",
                                   string a_outFile = "zdc_neutron_capture_check.root",
                                   double a_adcMax = 4000,
                                   int a_triggerId = 6, // 2026 O+O fastoffline dataset's own trigger, same convention as GetRefMultHisto.C. Pass <0 to disable
                                   double a_adcBinWidth = 20.0){ // ADC per bin. 2026-08-28: the first run (a_adcMax=4000 at this same 20 ADC/bin default) showed a smooth featureless falling tail, no resolvable 1n/2n/3n peak structure -- but the nonzero-ADC stats from that run (east mean 69, west mean 115, both under ~1000-2100 max) suggest any real peak spacing could be well under 20 ADC, which would blur multiple peaks into exactly that kind of smooth continuum. Pass a much smaller value (e.g. 2) together with a small a_adcMax (e.g. 300) to zoom in and check whether peaks are actually there.

  gSystem->Load("./bin/libStPicoDst.so");
  gSystem->mkdir(a_outDir.c_str(), true);

  // ---- Parse the filelist ourselves -- see FILE-BY-FILE READING header comment above.
  cout << "Parsing filelist: " << a_filelist << endl;
  vector<string> picoFileList;
  {
    std::ifstream inputStream(a_filelist.c_str());
    if(!inputStream){
      cout << "CheckZdcNeutronCapture_OO200 :: Error - cannot open filelist " << a_filelist << endl;
      return;
    }
    std::string line;
    while(std::getline(inputStream, line)){
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
    cout << "CheckZdcNeutronCapture_OO200 :: no valid picoDst files found -- nothing to do." << endl;
    return;
  }

  int nAdcBins = (int)(a_adcMax / a_adcBinWidth); // see a_adcBinWidth's own doc comment above
  TH1D* hZdcEast = new TH1D("hZdcEast","ZDC East Sum ADC;Sum ADC (east);Events",nAdcBins,0,a_adcMax);
  TH1D* hZdcWest = new TH1D("hZdcWest","ZDC West Sum ADC;Sum ADC (west);Events",nAdcBins,0,a_adcMax);

  int modNum = 1;
  Long64_t globalEventIndex = 0;
  Long64_t nTriggerRejected = 0;
  Long64_t nEastNonzero = 0, nWestNonzero = 0, nEitherNonzero = 0, nBothNonzero = 0;
  double eastMin = -1, eastMax = -1, eastSumNonzero = 0; // min/max/sum over NONZERO east values only
  double westMin = -1, westMax = -1, westSumNonzero = 0; // same, west

  // ---- FILE-BY-FILE EVENT LOOP -- see header comment: never let a single
  // StPicoDstReader/TChain span more than one file.
  for(size_t fileIndex = 0; fileIndex < picoFileList.size(); fileIndex++){

    StPicoDstReader* picoReader = new StPicoDstReader(picoFileList[fileIndex].c_str());
    picoReader->Init();
    picoReader->SetStatus("*",0);
    picoReader->SetStatus("Event",1); // ZDC sums live on StPicoEvent -- no track branches needed
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
        std::cerr << "CheckZdcNeutronCapture_OO200 :: readPicoEvent failed at file " << fileIndex
                   << ", local event " << localEventIndex << " -- stopping this file early." << std::endl;
        break;
      }

      StPicoEvent* event = dst->event();
      if(!event){
        cerr << "CheckZdcNeutronCapture_OO200 :: event " << globalEventIndex << " doesn't exist?!" << endl;
        continue;
      }

      if(a_triggerId >= 0 && !event->isTrigger((unsigned int) a_triggerId)){
        nTriggerRejected++;
        continue;
      }

      double eastAdc = (double) event->ZdcSumAdcEast();
      double westAdc = (double) event->ZdcSumAdcWest();
      hZdcEast->Fill(eastAdc);
      hZdcWest->Fill(westAdc);

      bool eastHit = eastAdc > 0.0;
      bool westHit = westAdc > 0.0;
      if(eastHit){
        nEastNonzero++;
        eastSumNonzero += eastAdc;
        if(eastMin < 0 || eastAdc < eastMin) eastMin = eastAdc;
        if(eastAdc > eastMax) eastMax = eastAdc;
      }
      if(westHit){
        nWestNonzero++;
        westSumNonzero += westAdc;
        if(westMin < 0 || westAdc < westMin) westMin = westAdc;
        if(westAdc > westMax) westMax = westAdc;
      }
      if(eastHit || westHit) nEitherNonzero++;
      if(eastHit && westHit) nBothNonzero++;
    }

    picoReader->Finish();
    delete picoReader;
  }
  cout << "Finished looping over all " << globalEventIndex << " events" << endl;
  if(a_triggerId >= 0){
    cout << "Trigger cut (isTrigger(" << a_triggerId << ")): rejected " << nTriggerRejected
         << " / " << globalEventIndex << " event(s)." << endl;
  }else{
    cout << "Trigger cut disabled (a_triggerId<0)." << endl;
  }
  Long64_t nAnalyzed = globalEventIndex - nTriggerRejected;

  // ---- Overflow check -- see BINNING header note: report rather than silently clip. ----
  double eastOverflow = hZdcEast->GetBinContent(hZdcEast->GetNbinsX()+1);
  double westOverflow = hZdcWest->GetBinContent(hZdcWest->GetNbinsX()+1);

  cout << endl << "==================================================================" << endl;
  cout << "ZDC neutron-capture check (naive ADC>0 proxy -- see header caveat)" << endl;
  cout << "==================================================================" << endl;
  cout << "Events analyzed (after trigger cut): " << nAnalyzed << endl;
  cout << "East ZDC ADC>0:  " << nEastNonzero << " / " << nAnalyzed
       << " (" << Form("%.2f", nAnalyzed>0 ? 100.0*nEastNonzero/nAnalyzed : 0.0) << "%)" << endl;
  cout << "West ZDC ADC>0:  " << nWestNonzero << " / " << nAnalyzed
       << " (" << Form("%.2f", nAnalyzed>0 ? 100.0*nWestNonzero/nAnalyzed : 0.0) << "%)" << endl;
  cout << "Either side>0:   " << nEitherNonzero << " / " << nAnalyzed
       << " (" << Form("%.2f", nAnalyzed>0 ? 100.0*nEitherNonzero/nAnalyzed : 0.0) << "%)" << endl;
  cout << "Both sides>0:    " << nBothNonzero << " / " << nAnalyzed
       << " (" << Form("%.2f", nAnalyzed>0 ? 100.0*nBothNonzero/nAnalyzed : 0.0) << "%)" << endl;
  if(nEastNonzero > 0){
    cout << "East (nonzero only): min=" << eastMin << " max=" << eastMax
         << " mean=" << eastSumNonzero/nEastNonzero << endl;
  }
  if(nWestNonzero > 0){
    cout << "West (nonzero only): min=" << westMin << " max=" << westMax
         << " mean=" << westSumNonzero/nWestNonzero << endl;
  }
  if(eastOverflow > 0 || westOverflow > 0){
    cout << endl << "WARNING: " << eastOverflow << " east / " << westOverflow << " west event(s) fell above"
         << " the histogram's a_adcMax=" << a_adcMax << " range (in the overflow bin, not lost, but not"
         << " shown on the plot) -- rerun with a larger a_adcMax to see them." << endl;
  }
  cout << endl << "Reminder: ADC>0 is NOT a calibrated single-neutron threshold -- it's unknown from"
       << " this macro alone whether 0 means true zero signal or an already-applied production-level"
       << " pedestal floor. Look at the saved histograms' peak structure before treating these counts"
       << " as a neutron multiplicity." << endl;

  // ---- Write histograms ----
  string outFilePath = a_outDir + "/" + a_outFile;
  TFile* outFile = new TFile(outFilePath.c_str(),"RECREATE");
  outFile->cd();
  hZdcEast->Write();
  hZdcWest->Write();
  outFile->Close();
  cout << endl << "Wrote " << outFilePath << endl;

  // ---- Draw the plot ----
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("CheckZdcNeutronCapture","ZDC Sum ADC",900,700);
  c->SetLogy();
  hZdcEast->SetLineColor(kRed+1);
  hZdcWest->SetLineColor(kAzure+2);
  hZdcEast->SetLineWidth(2);
  hZdcWest->SetLineWidth(2);
  hZdcEast->SetTitle(a_triggerId >= 0 ? Form("ZDC Sum ADC, Fastoffline, Trigger %d",a_triggerId) : "ZDC Sum ADC, Fastoffline, No Trigger Cut");
  hZdcEast->Draw("hist");
  hZdcWest->Draw("hist same");

  TLegend* leg = new TLegend(0.55,0.75,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hZdcEast,"East","l");
  leg->AddEntry(hZdcWest,"West","l");
  leg->Draw();

  string plotPath = a_outDir + "/zdc_sum_adc_fastoffline.png";
  c->SaveAs(plotPath.c_str());
  cout << "Wrote " << plotPath << endl;
}
