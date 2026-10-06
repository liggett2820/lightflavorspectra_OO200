// MakeZdcSpectatorNeutrons_OO200.C -- lightflavorspectra_OO200
//
// ZDC spectator-neutron histograms for O+O 200 GeV (Run 21), straight from picoDsts, for the
// trigger-860003 events that pass the analysis event cuts (same cuts as
// macros/MakeEventCutFlow_OO200.C: |Vz| <= 2 cm, Vr <= 1 cm, nBTOFMatch <= refMult + 100,
// refMult <= 260).
//
// Stored, per side (East, West):
//   - unattenuated ZDC ADC (StPicoEvent::zdcUnAttenuatedEast/West) vs refMult       [TH2F]
//   - the same per centrality class (classes as in macros/SetCutClass.C:
//     refMult >= {44,37,28,17,5,0} for 0-5, 5-10, 10-20, 20-40, 40-80, 80-100%)    [TH1D x 6]
//   - attenuated ZDC ADC sum (ZdcSumAdcEast/West) per class, for comparison          [TH1D x 6]
//   - unattenuated ADC vs run, 80-100% events only, for the per-run 1n-peak check   [TH2F]
//   - East vs West unattenuated ADC for 80-100% events                              [TH2F]
//   - <ADC> vs refMult                                                               [TProfile]
// plus a counter of events with zero unattenuated ADC on a side (to catch productions where the
// unattenuated value was not filled).
//
// Neutron counting is done afterwards in macros/PlotZdcSpectatorNeutrons_OO200.C. The reference
// 1n peak positions (East ~125, West ~111 ADC) come from D. Wen's O+O UPC talk (Oct 6, 2026),
// fitted on UPC triggers; they are cross-checked here on our own peripheral events, not assumed.
//
// Usage (repo root, RCF, inside the SL7 container, SL24y):
//   root -l -b -q macros/LoadCutFlowLibs.C 'macros/MakeZdcSpectatorNeutrons_OO200.C+("files.list","zdc_part0.root")'
// The output is additive (hadd).

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TChain.h"
#include "TString.h"
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"

namespace {
  const unsigned int kTrigger = 860003;
  const double kVzMax   = 2.0;    // cm
  const double kVrMax   = 1.0;    // cm
  const int    kTofLine = 100;    // reject nBTOFMatch > refMult + 100
  const int    kRefMax  = 260;    // reject refMult > 260
  const int    kNCent   = 6;
  const int    kCentEdge[kNCent] = {44, 37, 28, 17, 5, 0};   // class i: refMult >= kCentEdge[i]
  const char*  kCentLab[kNCent]  = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
  const char*  kSide[2] = {"East", "West"};
  const int    kAdcBins = 1024;  const double kAdcMax = 4096.;   // 4 ADC per bin
  // run index for the per-run histograms: (day - 130) * 100 + (run number within the day)
  int runIndex(int runId){ int day = (runId % 1000000) / 1000; return (day - 130) * 100 + (runId % 1000); }
  int centClass(int refMult){ for(int i = 0; i < kNCent; i++) if(refMult >= kCentEdge[i]) return i; return -1; }
}

void MakeZdcSpectatorNeutrons_OO200(const char* fileList, const char* outName = "zdc_OO200.root", int maxFiles = -1){
  std::vector<std::string> files; std::ifstream in(fileList); std::string line;
  while(std::getline(in, line)){
    std::string f = line.substr(0, line.find_first_of(" \t"));
    if(f.find(".picoDst.root") != std::string::npos) files.push_back(f);
  }
  if(files.empty()){ std::cout << "No .picoDst.root files in " << fileList << std::endl; return; }
  if(maxFiles > 0 && (int)files.size() > maxFiles) files.resize(maxFiles);

  TFile* out = TFile::Open(outName, "RECREATE");
  TH1D* hEvents = new TH1D("hEvents", "Events;;Events", 4, 0, 4);
  hEvents->GetXaxis()->SetBinLabel(1, "860003, all event cuts");
  hEvents->GetXaxis()->SetBinLabel(2, "UA East == 0");
  hEvents->GetXaxis()->SetBinLabel(3, "UA West == 0");
  hEvents->GetXaxis()->SetBinLabel(4, "both UA == 0");
  TH1D* hCent = new TH1D("hCentClass", "Events per centrality class;;Events", kNCent, 0, kNCent);
  for(int i = 0; i < kNCent; i++) hCent->GetXaxis()->SetBinLabel(i+1, kCentLab[i]);

  TH2F* hUAvsRef[2]; TProfile* pUAvsRef[2]; TH2F* hUAvsRun[2];
  TH1D* hUA[2][kNCent]; TH1D* hSum[2][kNCent];
  for(int s = 0; s < 2; s++){
    hUAvsRef[s] = new TH2F(Form("hZdcUA%s_vs_refMult", kSide[s]), Form("ZDC %s unattenuated ADC vs refMult;refMult;ZDC %s unattenuated ADC", kSide[s], kSide[s]),
                           150, 0, 150, kAdcBins/2, 0, kAdcMax);
    pUAvsRef[s] = new TProfile(Form("pZdcUA%s_vs_refMult", kSide[s]), Form("<ZDC %s unattenuated ADC> vs refMult;refMult;<ADC>", kSide[s]), 150, 0, 150);
    hUAvsRun[s] = new TH2F(Form("hZdcUA%s_vs_run_peripheral", kSide[s]), Form("ZDC %s unattenuated ADC, 80-100%%;run index = (day-130)*100 + run-in-day;ADC", kSide[s]),
                           700, 0, 700, 250, 0, 1000);
    for(int c = 0; c < kNCent; c++){
      hUA[s][c]  = new TH1D(Form("hZdcUA%s_cent%d", kSide[s], c),  Form("ZDC %s unattenuated ADC, %s;ADC;Events", kSide[s], kCentLab[c]), kAdcBins, 0, kAdcMax);
      hSum[s][c] = new TH1D(Form("hZdcSum%s_cent%d", kSide[s], c), Form("ZDC %s ADC sum (attenuated), %s;ADC;Events", kSide[s], kCentLab[c]), kAdcBins, 0, kAdcMax);
    }
  }
  TH2F* hEW = new TH2F("hZdcUA_EastVsWest_peripheral", "ZDC unattenuated ADC, 80-100%;East ADC;West ADC", 250, 0, 1000, 250, 0, 1000);

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
      TVector3 pv = ev->primaryVertex();
      if(TMath::Abs(pv.Z()) > kVzMax) continue;
      if(TMath::Sqrt(pv.X()*pv.X() + pv.Y()*pv.Y()) > kVrMax) continue;
      const int rm = ev->refMult();
      if(ev->nBTOFMatch() > rm + kTofLine || rm > kRefMax) continue;
      const int c = centClass(rm);
      if(c < 0) continue;
      nPass++;
      hEvents->Fill(0.5); hCent->Fill(c + 0.5);

      const double ua[2]  = {ev->zdcUnAttenuatedEast(), ev->zdcUnAttenuatedWest()};
      const double sum[2] = {ev->ZdcSumAdcEast(), ev->ZdcSumAdcWest()};
      if(ua[0] <= 0) hEvents->Fill(1.5);
      if(ua[1] <= 0) hEvents->Fill(2.5);
      if(ua[0] <= 0 && ua[1] <= 0) hEvents->Fill(3.5);
      for(int s = 0; s < 2; s++){
        hUAvsRef[s]->Fill(rm, ua[s]); pUAvsRef[s]->Fill(rm, ua[s]);
        hUA[s][c]->Fill(ua[s]); hSum[s][c]->Fill(sum[s]);
        if(c == kNCent - 1) hUAvsRun[s]->Fill(runIndex(ev->runId()), ua[s]);
      }
      if(c == kNCent - 1) hEW->Fill(ua[0], ua[1]);
    }
    reader->Finish();
    delete reader;
    if(f % 10 == 0) std::cout << "file " << f+1 << "/" << files.size() << "  events passing so far: " << nPass << std::endl;
  }

  out->Write();
  std::cout << "\nEvents passing (860003 + event cuts): " << nPass << std::endl;
  for(int c = 0; c < kNCent; c++) std::cout << Form("  %-8s %12.0f", kCentLab[c], hCent->GetBinContent(c+1)) << std::endl;
  std::cout << Form("Unattenuated ADC == 0: East %.0f, West %.0f, both %.0f", hEvents->GetBinContent(2), hEvents->GetBinContent(3), hEvents->GetBinContent(4)) << std::endl;
  out->Close();
  std::cout << "Wrote " << outName << std::endl;
}
