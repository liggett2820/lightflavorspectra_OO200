// MakeRefMultCorrInputs_OO200.C -- lightflavorspectra_OO200
//
// Step 1 of building an StRefMultCorr-style centrality for O+O 200 GeV (Run 21), trigger 860003.
// Fills the event-level histograms that the corrections are fitted from; the fits are done in
// macros/FitRefMultCorr_OO200.C. Reads only the Event branch, so jobs are short.
//
// Event selection: trigger 860003, |Vz| <= 30 cm (StRefMultCorr range for O+O), Vr <= 1 cm.
// No pile-up cut is applied here except where noted, because the pile-up cut is one of the
// things being derived. The run list (56 good 860003 runs) is applied by the catalog query in
// xml/runRefMultCorrInputs_OO200_SDCC.xml.
//
// Histograms:
//   hRefMult_vs_Vz          refMult vs Vz (1 cm bins)                 -> Vz correction
//   hRefMult_vs_ZDCx        refMult vs ZDCx, |Vz| <= 10 cm            -> luminosity correction
//   hRefMult_vs_nBTOFMatch  refMult vs nBTOFMatch, |Vz| <= 30 cm      -> pile-up band
//   hRefMult_Vz_ZDCx        refMult x Vz x ZDCx, nBTOFMatch <= refMult + 100 (current cut)
//                           -> corrected refMult distribution for the NBD refit
//   pRefMult_vs_run, pNBTOFMatch_vs_run, pZDCx_vs_run               -> run-by-run stability
//   hVz, hZDCx, hEvents
//
// Usage (repo root, RCF, inside the SL7 container, SL24y):
//   root -l -b -q macros/LoadCutFlowLibs.C 'macros/MakeRefMultCorrInputs_OO200.C+("files.list","rmc_part0.root")'
// Output is additive (hadd).

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TChain.h"
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"

namespace {
  const unsigned int kTrigger = 860003;
  const double kVzMax   = 30.0;   // cm
  const double kVrMax   = 1.0;    // cm
  const double kVzLumi  = 10.0;   // cm, Vz window for the luminosity dependence
  const int    kTofLine = 100;    // current pile-up cut: nBTOFMatch <= refMult + 100
  // run index: (day - 130) * 100 + run number within the day (as in the ZDC macros)
  int runIndex(int runId){ int day = (runId % 1000000) / 1000; return (day - 130) * 100 + (runId % 1000); }
}

void MakeRefMultCorrInputs_OO200(const char* fileList, const char* outName = "rmc_OO200.root", int maxFiles = -1){
  std::vector<std::string> files; std::ifstream in(fileList); std::string line;
  while(std::getline(in, line)){
    std::string f = line.substr(0, line.find_first_of(" \t"));
    if(f.find(".picoDst.root") != std::string::npos) files.push_back(f);
  }
  if(files.empty()){ std::cout << "No .picoDst.root files in " << fileList << std::endl; return; }
  if(maxFiles > 0 && (int)files.size() > maxFiles) files.resize(maxFiles);

  TFile* out = TFile::Open(outName, "RECREATE");
  TH1D* hEvents = new TH1D("hEvents", "Events;;Events", 4, 0, 4);
  hEvents->GetXaxis()->SetBinLabel(1, "860003");
  hEvents->GetXaxis()->SetBinLabel(2, "|Vz| <= 30 cm");
  hEvents->GetXaxis()->SetBinLabel(3, "Vr <= 1 cm");
  hEvents->GetXaxis()->SetBinLabel(4, "nBTOFMatch <= refMult+100");

  TH1D* hVz   = new TH1D("hVz", "860003, Vr <= 1 cm;V_{z} (cm);Events", 600, -30, 30);
  TH1D* hZDCx = new TH1D("hZDCx", "860003, |V_{z}| <= 30 cm, Vr <= 1 cm;ZDCx (kHz);Events", 300, 0, 3);

  TH2F* hRefVz  = new TH2F("hRefMult_vs_Vz", "860003, Vr <= 1 cm;V_{z} (cm);refMult", 60, -30, 30, 300, 0, 300);
  TH2F* hRefZdc = new TH2F("hRefMult_vs_ZDCx", "860003, |V_{z}| <= 10 cm, Vr <= 1 cm;ZDCx (kHz);refMult", 60, 0, 3, 300, 0, 300);
  TH2F* hRefTof = new TH2F("hRefMult_vs_nBTOFMatch", "860003, |V_{z}| <= 30 cm, Vr <= 1 cm;nBTOFMatch;refMult", 300, 0, 300, 300, 0, 300);
  TH3F* h3      = new TH3F("hRefMult_Vz_ZDCx", "860003, |V_{z}| <= 30 cm, Vr <= 1 cm, nBTOFMatch <= refMult+100;refMult;V_{z} (cm);ZDCx (kHz)",
                           300, 0, 300, 60, -30, 30, 30, 0, 3);

  TProfile* pRefRun = new TProfile("pRefMult_vs_run",    "860003, |V_{z}| <= 30 cm, Vr <= 1 cm;run index = (day-130)*100 + run-in-day;#LTrefMult#GT", 700, 0, 700);
  TProfile* pTofRun = new TProfile("pNBTOFMatch_vs_run", "860003, |V_{z}| <= 30 cm, Vr <= 1 cm;run index;#LTnBTOFMatch#GT", 700, 0, 700);
  TProfile* pZdcRun = new TProfile("pZDCx_vs_run",       "860003, |V_{z}| <= 30 cm, Vr <= 1 cm;run index;#LTZDCx#GT (kHz)", 700, 0, 700);

  long long nPass = 0;
  for(size_t f = 0; f < files.size(); f++){
    StPicoDstReader* reader = new StPicoDstReader(files[f].c_str());
    reader->Init();
    reader->SetStatus("*", 0);
    reader->SetStatus("Event", 1);
    Long64_t n = reader->chain()->GetEntries();
    for(Long64_t i = 0; i < n; i++){
      if(!reader->readPicoEvent(i)) continue;
      StPicoEvent* ev = reader->picoDst()->event();
      if(!ev || !ev->isTrigger(kTrigger)) continue;
      hEvents->Fill(0.5);
      TVector3 pv = ev->primaryVertex();
      const double vz = pv.Z();
      if(TMath::Abs(vz) > kVzMax) continue;
      hEvents->Fill(1.5);
      if(TMath::Sqrt(pv.X()*pv.X() + pv.Y()*pv.Y()) > kVrMax) continue;
      hEvents->Fill(2.5);

      const int    rm  = ev->refMult();
      const int    tof = ev->nBTOFMatch();
      const double zdc = ev->ZDCx()/1000.0;   // Hz -> kHz
      const int    ri  = runIndex(ev->runId());

      hVz->Fill(vz); hZDCx->Fill(zdc);
      hRefVz->Fill(vz, rm);
      hRefTof->Fill(tof, rm);
      if(TMath::Abs(vz) <= kVzLumi) hRefZdc->Fill(zdc, rm);
      pRefRun->Fill(ri, rm); pTofRun->Fill(ri, tof); pZdcRun->Fill(ri, zdc);

      if(tof > rm + kTofLine) continue;
      hEvents->Fill(3.5);
      h3->Fill(rm, vz, zdc);
      nPass++;
    }
    reader->Finish();
    delete reader;
    if(f % 10 == 0) std::cout << "file " << f+1 << "/" << files.size() << "  events so far: " << nPass << std::endl;
  }

  out->Write();
  std::cout << "\n860003 events: " << hEvents->GetBinContent(1)
            << ", |Vz|<=30: " << hEvents->GetBinContent(2)
            << ", Vr<=1: " << hEvents->GetBinContent(3)
            << ", pile-up line: " << hEvents->GetBinContent(4) << std::endl;
  out->Close();
  std::cout << "Wrote " << outName << std::endl;
}
