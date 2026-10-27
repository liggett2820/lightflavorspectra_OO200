// MakeZdcSpectatorNeutrons_OO200.C -- lightflavorspectra_OO200
//
// ZDC spectator-neutron histograms for O+O 200 GeV (Run 21), straight from picoDsts, for the
// minbias-bbc_etof (860003) events that pass the SAME event selection and centrality as the spectra
// production (PicoBinner with macros/SetCutClass.C, Oct 26 2026):
//   |Vz| <= 30 cm, Vr <= 1 cm,
//   pile-up: the pass-1 refMult-nBTOFMatch pol4 band (as CutClass::setPileUpPolyBand), plus
//            nBTOFMatch <= refMult + 100 and refMult <= 260,
//   centrality: corrected EPD totnMIP (headers/EpdCentrality_OO200.h), classes
//            totnMIP_corr >= {364,310,234,130,36,0} for 0-5, 5-10, 10-20, 20-40, 40-80, 80-100%.
// Until Oct 26 2026 this macro used |Vz| <= 2 cm, the linear pile-up lines only and refMult classes
// {44,37,28,17,5,0}; the July ZDC plots were made that way.
//
// Stored, per side (East, West):
//   - unattenuated ZDC ADC (StPicoEvent::zdcUnAttenuatedEast/West) per centrality class  [TH1D x 6]
//   - the same in 10% sub-classes of 40-100% (totnMIP_corr >= {95,70,51,36,23,0}
//     for 40-50 ... 90-100%, edges from the same NBD-Glauber fit)                       [TH1D x 6]
//   - attenuated ZDC ADC sum (ZdcSumAdcEast/West) per class, for comparison              [TH1D x 6]
//   - unattenuated ADC vs corrected totnMIP, and vs refMult                              [TH2F, TProfile]
//   - unattenuated ADC vs run, 80-100% events only, for the per-run 1n-peak check        [TH2F]
//   - East vs West unattenuated ADC for 80-100% events                                   [TH2F]
// plus events per class and sub-class, and counters of events with zero unattenuated ADC on a side.
// Histogram names of the six classes are unchanged, so macros/PlotZdcSpectatorNeutrons_OO200.C reads
// the output as before.
//
// Usage (repo root, RCF, inside the SL7 container, SL24y):
//   root -l -b -q macros/LoadCutFlowLibs.C 'macros/MakeZdcSpectatorNeutrons_OO200.C+("files.list","zdc_part0.root")'
// The output is additive (hadd).

#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cstdio>
#include <unistd.h>
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
#include "StPicoEpdHit.h"
#include "../headers/EpdCentrality_OO200.h"

namespace {
  const unsigned int kTrigger = 860003;
  const double kVzMax   = 30.0;   // cm
  const double kVrMax   = 1.0;    // cm
  const int    kTofLine = 100;    // reject nBTOFMatch > refMult + 100
  const int    kRefMax  = 260;    // reject refMult > 260
  // pass-1 refMult-nBTOFMatch pile-up band, as in SetCutClass.C / MakeEpdCentInputs_OO200.C
  const double kBandUp[5] = {5.91052, 2.22934, -0.0221511, 0.000203259, -7.82252e-07};
  const double kBandLo[5] = {-4.67124, -0.259727, 0.0274938, -0.00035788, 1.67535e-06};
  const double kBandXm    = 85.;
  double bandEdge(const double* p, double t){
    const double tt = (t <= kBandXm) ? t : kBandXm;
    double v = p[0] + p[1]*tt + p[2]*tt*tt + p[3]*tt*tt*tt + p[4]*tt*tt*tt*tt;
    if(t > kBandXm) v += (p[1] + 2*p[2]*kBandXm + 3*p[3]*kBandXm*kBandXm + 4*p[4]*kBandXm*kBandXm*kBandXm) * (t - kBandXm);
    return v;
  }
  bool inBand(int tof, int rm){ return rm >= bandEdge(kBandLo, tof) && rm <= bandEdge(kBandUp, tof); }

  const int    kNCent = 6;
  const double kCentEdge[kNCent] = {364, 310, 234, 130, 36, 0};   // class i: totnMIP_corr >= kCentEdge[i]
  const char*  kCentLab[kNCent]  = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
  const int    kNSub = 6;                                         // 10% sub-classes of 40-100%
  const double kSubEdge[kNSub] = {95, 70, 51, 36, 23, 0};         // used only below 130 (40%)
  const char*  kSubLab[kNSub]  = {"40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};
  const char*  kSide[2] = {"East", "West"};
  const int    kAdcBins = 1024;  const double kAdcMax = 4096.;   // 4 ADC per bin
  int runIndex(int runId){ int day = (runId % 1000000) / 1000; return (day - 130) * 100 + (runId % 1000); }
  int centClass(double m){ for(int i = 0; i < kNCent; i++) if(m >= kCentEdge[i]) return i; return -1; }
  int subClass(double m){ if(m >= kCentEdge[3]) return -1; for(int i = 0; i < kNSub; i++) if(m >= kSubEdge[i]) return i; return -1; }
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
  TH1D* hCut = new TH1D("hCutFlow", "Event selection;;Events", 6, 0, 6);
  const char* cutLab[6] = {"860003", "|Vz| <= 30 cm", "Vr <= 1 cm", "pile-up band", "linear lines, refMult <= 260", "in a class"};
  for(int i = 0; i < 6; i++) hCut->GetXaxis()->SetBinLabel(i+1, cutLab[i]);
  TH1D* hCent = new TH1D("hCentClass", "Events per centrality class (EPD);;Events", kNCent, 0, kNCent);
  for(int i = 0; i < kNCent; i++) hCent->GetXaxis()->SetBinLabel(i+1, kCentLab[i]);
  TH1D* hSub = new TH1D("hSubClass", "Events per 10% sub-class of 40-100% (EPD);;Events", kNSub, 0, kNSub);
  for(int i = 0; i < kNSub; i++) hSub->GetXaxis()->SetBinLabel(i+1, kSubLab[i]);
  TH1D* hMip = new TH1D("hTotnMIPCorr", "Corrected EPD totnMIP;totnMIP_{corr};Events", 650, 0, 650);

  TH2F* hUAvsRef[2]; TProfile* pUAvsRef[2]; TH2F* hUAvsMip[2]; TProfile* pUAvsMip[2]; TH2F* hUAvsRun[2];
  TH1D* hUA[2][kNCent]; TH1D* hSum[2][kNCent]; TH1D* hUASub[2][kNSub];
  for(int s = 0; s < 2; s++){
    hUAvsRef[s] = new TH2F(Form("hZdcUA%s_vs_refMult", kSide[s]), Form("ZDC %s unattenuated ADC vs refMult;refMult;ZDC %s unattenuated ADC", kSide[s], kSide[s]),
                           150, 0, 150, kAdcBins/2, 0, kAdcMax);
    pUAvsRef[s] = new TProfile(Form("pZdcUA%s_vs_refMult", kSide[s]), Form("<ZDC %s unattenuated ADC> vs refMult;refMult;<ADC>", kSide[s]), 150, 0, 150);
    hUAvsMip[s] = new TH2F(Form("hZdcUA%s_vs_totnMIP", kSide[s]), Form("ZDC %s unattenuated ADC vs corrected totnMIP;totnMIP_{corr};ZDC %s unattenuated ADC", kSide[s], kSide[s]),
                           130, 0, 650, kAdcBins/2, 0, kAdcMax);
    pUAvsMip[s] = new TProfile(Form("pZdcUA%s_vs_totnMIP", kSide[s]), Form("<ZDC %s unattenuated ADC> vs corrected totnMIP;totnMIP_{corr};<ADC>", kSide[s]), 130, 0, 650);
    hUAvsRun[s] = new TH2F(Form("hZdcUA%s_vs_run_peripheral", kSide[s]), Form("ZDC %s unattenuated ADC, 80-100%%;run index = (day-130)*100 + run-in-day;ADC", kSide[s]),
                           700, 0, 700, 250, 0, 1000);
    for(int c = 0; c < kNCent; c++){
      hUA[s][c]  = new TH1D(Form("hZdcUA%s_cent%d", kSide[s], c),  Form("ZDC %s unattenuated ADC, %s;ADC;Events", kSide[s], kCentLab[c]), kAdcBins, 0, kAdcMax);
      hSum[s][c] = new TH1D(Form("hZdcSum%s_cent%d", kSide[s], c), Form("ZDC %s ADC sum (attenuated), %s;ADC;Events", kSide[s], kCentLab[c]), kAdcBins, 0, kAdcMax);
    }
    for(int c = 0; c < kNSub; c++)
      hUASub[s][c] = new TH1D(Form("hZdcUA%s_sub%d", kSide[s], c), Form("ZDC %s unattenuated ADC, %s;ADC;Events", kSide[s], kSubLab[c]), kAdcBins, 0, kAdcMax);
  }
  TH2F* hEW = new TH2F("hZdcUA_EastVsWest_peripheral", "ZDC unattenuated ADC, 80-100%;East ADC;West ADC", 250, 0, 1000, 250, 0, 1000);

  long long nPass = 0;
  for(size_t f = 0; f < files.size(); f++){
    StPicoDstReader* reader = new StPicoDstReader(files[f].c_str());
    reader->Init();
    reader->SetStatus("*", 0);
    reader->SetStatus("Event", 1);
    reader->SetStatus("EpdHit*", 1);
    Long64_t n = reader->chain()->GetEntries();
    for(Long64_t i = 0; i < n; i++){
      if(!reader->readPicoEvent(i)) continue;
      StPicoDst* dst = reader->picoDst();
      StPicoEvent* ev = dst->event();
      if(!ev || !ev->isTrigger(kTrigger)) continue;
      hCut->Fill(0.5);
      TVector3 pv = ev->primaryVertex();
      if(TMath::Abs(pv.Z()) > kVzMax) continue;
      hCut->Fill(1.5);
      if(TMath::Sqrt(pv.X()*pv.X() + pv.Y()*pv.Y()) > kVrMax) continue;
      hCut->Fill(2.5);
      const int rm  = ev->refMult();
      const int tof = ev->nBTOFMatch();
      if(!inBand(tof, rm)) continue;
      hCut->Fill(3.5);
      if(tof > rm + kTofLine || rm > kRefMax) continue;
      hCut->Fill(4.5);

      double raw = 0;
      for(int h = 0; h < (int) dst->numberOfEpdHits(); h++){
        StPicoEpdHit* hit = dst->epdHit(h);
        if(!hit) continue;
        raw += EpdCentOO200::truncMip(hit->nMIP());
      }
      const double mip = EpdCentOO200::correct(raw, pv.Z(), ev->ZDCx());
      const int c  = centClass(mip);
      if(c < 0) continue;
      const int sc = subClass(mip);
      hCut->Fill(5.5);
      nPass++;
      hEvents->Fill(0.5); hCent->Fill(c + 0.5); hMip->Fill(mip);
      if(sc >= 0) hSub->Fill(sc + 0.5);

      const double ua[2]  = {ev->zdcUnAttenuatedEast(), ev->zdcUnAttenuatedWest()};
      const double sum[2] = {ev->ZdcSumAdcEast(), ev->ZdcSumAdcWest()};
      if(ua[0] <= 0) hEvents->Fill(1.5);
      if(ua[1] <= 0) hEvents->Fill(2.5);
      if(ua[0] <= 0 && ua[1] <= 0) hEvents->Fill(3.5);
      for(int s = 0; s < 2; s++){
        hUAvsRef[s]->Fill(rm, ua[s]); pUAvsRef[s]->Fill(rm, ua[s]);
        hUAvsMip[s]->Fill(mip, ua[s]); pUAvsMip[s]->Fill(mip, ua[s]);
        hUA[s][c]->Fill(ua[s]); hSum[s][c]->Fill(sum[s]);
        if(sc >= 0) hUASub[s][sc]->Fill(ua[s]);
        if(c == kNCent - 1) hUAvsRun[s]->Fill(runIndex(ev->runId()), ua[s]);
      }
      if(c == kNCent - 1) hEW->Fill(ua[0], ua[1]);
    }
    reader->Finish();
    delete reader;
    if(f % 10 == 0) std::cout << "file " << f+1 << "/" << files.size() << "  events passing so far: " << nPass << std::endl;
  }

  out->Write();
  std::cout << "\nEvent selection:" << std::endl;
  for(int i = 0; i < 6; i++) std::cout << Form("  %-30s %12.0f", cutLab[i], hCut->GetBinContent(i+1)) << std::endl;
  std::cout << "Events per class (EPD):" << std::endl;
  for(int c = 0; c < kNCent; c++) std::cout << Form("  %-8s %12.0f  (%5.2f%%)", kCentLab[c], hCent->GetBinContent(c+1), nPass > 0 ? 100.0*hCent->GetBinContent(c+1)/nPass : 0.0) << std::endl;
  for(int c = 0; c < kNSub; c++) std::cout << Form("  %-8s %12.0f", kSubLab[c], hSub->GetBinContent(c+1)) << std::endl;
  std::cout << Form("Unattenuated ADC == 0: East %.0f, West %.0f, both %.0f", hEvents->GetBinContent(2), hEvents->GetBinContent(3), hEvents->GetBinContent(4)) << std::endl;
  out->Close();
  std::cout << "Wrote " << outName << std::endl;
  // The output file is already closed here. Batch use only (root -b -q): this ends the session.
  // ROOT aborts or hangs at teardown (double-registered PicoDst dictionaries), which kept SUMS
  // from copying the output back; same fix as MakeEpdCentInputs_OO200.C.
  std::cout.flush(); fflush(stdout); fflush(stderr);
  _exit(0);
}
