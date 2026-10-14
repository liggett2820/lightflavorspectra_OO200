// MakeEpdCentInputs_OO200.C -- lightflavorspectra_OO200
//
// EPD (totnMIP) counterpart of macros/MakeRefMultCorrInputs_OO200.C, for an EPD-based centrality
// in O+O 200 GeV (Run 21), trigger 860003. refMult cannot resolve centrality beyond ~75% (the 80%
// edge falls at refMult = 0 in the Oct 14 2026 NBD fit); the EPD (2.1 < |eta| < 5.1) can, and it
// avoids the autocorrelation of selecting on the mid-rapidity tracks the spectra are measured with.
//
// totnMIP follows the Run 21 O+O definition in StRefMultCorr("totnmip") (star-sw main, Apr 2026):
//   loop over all EPD hits (both sides, all rows): nMIP < 0.3 -> skip, nMIP > 6 -> add 6, else add nMIP.
// (StPicoEpdHit::nMIP() already returns 0 for tiles flagged bad.) The official definition's
// triggers are 860001/2/11/12; its cuts (80%..0-5%) are 37, 46, 55, 66, 77, 91, 105, 122, 141, 162,
// 186, 213, 242, 276, 315, 367. Here the same quantity is filled for 860003 so that our own
// NBD-Glauber fit can be compared with that table. Events that fired 860001 or 860002 are filled
// into a separate set (suffix _off) for a trigger comparison at low multiplicity only.
//
// Event selection: |Vz| <= 30 cm, Vr <= 1 cm; the refMult pile-up band from pass 1
// (MakeRefMultCorrInputs_OO200.C) is applied where noted. No totnMIP pile-up cut yet: deriving it
// (totnMIP vs nBTOFMatch band) is one of the things these histograms are for.
//
// Histograms (860003 unless _off):
//   hTotnMIP_vs_Vz            totnMIP vs Vz (1 cm bins)                -> Vz correction
//   hTotnMIP_vs_ZDCx          totnMIP vs ZDCx, |Vz| <= 10 cm           -> luminosity correction
//   hTotnMIP_vs_nBTOFMatch    totnMIP vs nBTOFMatch                    -> totnMIP pile-up band
//   hTotnMIP_vs_refMult       totnMIP vs refMult, refMult band applied -> refMult/EPD comparison
//   hTotnMIP_Vz_ZDCx          totnMIP x Vz x ZDCx, refMult band applied -> corrected totnMIP for the NBD fit
//   hTotnMIPEast_vs_West      East vs West sums, refMult band applied
//   hTotnMIPnoRow1_vs_Vz      totnMIP without the innermost row (row 1), refMult band applied
//                             -> check for spectator/fragment hits near beam rapidity
//   pRowMIP_East/West_vs_refMult  mean truncated nMIP per row (1-16) in refMult slices
//   hTotnMIP_vs_Vz_off, hTotnMIP_vs_refMult_off   same, events that fired 860001 or 860002
//   pTotnMIP_vs_run, pNEpdHits_vs_run, hEvents, hTrigOverlap
//
// Usage (repo root, RCF, inside the SL7 container, SL24y):
//   root -l -b -q macros/LoadCutFlowLibs.C 'macros/MakeEpdCentInputs_OO200.C+("files.list","epd_part0.root")'
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
#include "TProfile2D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TChain.h"
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoEpdHit.h"

namespace {
  const unsigned int kTrigger = 860003;
  const unsigned int kOff1 = 860001, kOff2 = 860002;   // official O+O minbias triggers (comparison only)
  const double kVzMax   = 30.0;   // cm
  const double kVrMax   = 1.0;    // cm
  const double kVzLumi  = 10.0;   // cm
  const double kMipMin  = 0.3;    // StRefMultCorr totnMIP threshold
  const double kMipMax  = 6.0;    // StRefMultCorr totnMIP truncation
  const int    kNMip    = 650;    // totnMIP histogram range (1 nMIP bins), as oo200_run21_totnMIPBins
  // pass-1 refMult pile-up band (same constants as MakeRefMultCorrInputs_OO200.C)
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
  int runIndex(int runId){ int day = (runId % 1000000) / 1000; return (day - 130) * 100 + (runId % 1000); }
  double truncMip(double m){ if(m < kMipMin) return 0.; return (m > kMipMax) ? kMipMax : m; }
}

void MakeEpdCentInputs_OO200(const char* fileList, const char* outName = "epd_OO200.root", int maxFiles = -1){
  std::vector<std::string> files; std::ifstream in(fileList); std::string line;
  while(std::getline(in, line)){
    std::string f = line.substr(0, line.find_first_of(" \t"));
    if(f.find(".picoDst.root") != std::string::npos) files.push_back(f);
  }
  if(files.empty()){ std::cout << "No .picoDst.root files in " << fileList << std::endl; return; }
  if(maxFiles > 0 && (int)files.size() > maxFiles) files.resize(maxFiles);

  TFile* out = TFile::Open(outName, "RECREATE");
  TH1D* hEvents = new TH1D("hEvents", "Events;;Events", 5, 0, 5);
  hEvents->GetXaxis()->SetBinLabel(1, "860003");
  hEvents->GetXaxis()->SetBinLabel(2, "|Vz| <= 30 cm");
  hEvents->GetXaxis()->SetBinLabel(3, "Vr <= 1 cm");
  hEvents->GetXaxis()->SetBinLabel(4, "refMult pile-up band (pass 1)");
  hEvents->GetXaxis()->SetBinLabel(5, "totnMIP = 0");
  TH2F* hTrig = new TH2F("hTrigOverlap", "|Vz| <= 30 cm, Vr <= 1 cm;fired 860003;fired 860001 or 860002", 2, 0, 2, 2, 0, 2);

  TH2F* hMipVz   = new TH2F("hTotnMIP_vs_Vz", "860003, Vr <= 1 cm;V_{z} (cm);totnMIP", 60, -30, 30, kNMip, 0, kNMip);
  TH2F* hMipZdc  = new TH2F("hTotnMIP_vs_ZDCx", "860003, |V_{z}| <= 10 cm;ZDCx (kHz);totnMIP", 60, 0, 3, kNMip, 0, kNMip);
  TH2F* hMipTof  = new TH2F("hTotnMIP_vs_nBTOFMatch", "860003, |V_{z}| <= 30 cm;nBTOFMatch;totnMIP", 300, 0, 300, kNMip, 0, kNMip);
  TH2F* hMipRef  = new TH2F("hTotnMIP_vs_refMult", "860003, refMult pile-up band;refMult;totnMIP", 300, 0, 300, kNMip, 0, kNMip);
  TH3F* h3       = new TH3F("hTotnMIP_Vz_ZDCx", "860003, refMult pile-up band;totnMIP;V_{z} (cm);ZDCx (kHz)",
                            kNMip, 0, kNMip, 60, -30, 30, 15, 0, 3);
  TH2F* hEW      = new TH2F("hTotnMIPEast_vs_West", "860003, refMult pile-up band;totnMIP West;totnMIP East", 400, 0, 400, 400, 0, 400);
  TH2F* hNoRow1  = new TH2F("hTotnMIPnoRow1_vs_Vz", "860003, refMult pile-up band, row 1 excluded;V_{z} (cm);totnMIP (rows 2-16)", 60, -30, 30, kNMip, 0, kNMip);
  TProfile2D* pRowE = new TProfile2D("pRowMIP_East_vs_refMult", "860003, East;refMult;EPD row;#LTtruncated nMIP per event#GT", 30, 0, 150, 16, 0.5, 16.5);
  TProfile2D* pRowW = new TProfile2D("pRowMIP_West_vs_refMult", "860003, West;refMult;EPD row;#LTtruncated nMIP per event#GT", 30, 0, 150, 16, 0.5, 16.5);

  TH2F* hMipVzOff  = new TH2F("hTotnMIP_vs_Vz_off", "860001 or 860002, Vr <= 1 cm;V_{z} (cm);totnMIP", 60, -30, 30, kNMip, 0, kNMip);
  TH2F* hMipRefOff = new TH2F("hTotnMIP_vs_refMult_off", "860001 or 860002, |V_{z}| <= 30 cm;refMult;totnMIP", 300, 0, 300, kNMip, 0, kNMip);

  TProfile* pMipRun = new TProfile("pTotnMIP_vs_run", "860003, |V_{z}| <= 30 cm;run index = (day-130)*100 + run-in-day;#LTtotnMIP#GT", 700, 0, 700);
  TProfile* pHitRun = new TProfile("pNEpdHits_vs_run", "860003, |V_{z}| <= 30 cm;run index;#LTnumber of EPD hits#GT", 700, 0, 700);

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
      if(!ev) continue;
      const bool isMine = ev->isTrigger(kTrigger);
      const bool isOff  = ev->isTrigger(kOff1) || ev->isTrigger(kOff2);
      if(!isMine && !isOff) continue;
      if(isMine) hEvents->Fill(0.5);
      TVector3 pv = ev->primaryVertex();
      const double vz = pv.Z();
      if(TMath::Abs(vz) > kVzMax) continue;
      if(isMine) hEvents->Fill(1.5);
      if(TMath::Sqrt(pv.X()*pv.X() + pv.Y()*pv.Y()) > kVrMax) continue;
      if(isMine) hEvents->Fill(2.5);
      hTrig->Fill(isMine ? 1.5 : 0.5, isOff ? 1.5 : 0.5);

      // totnMIP as in StRefMultCorr("totnmip"); also East/West, without row 1, and per row
      double tot = 0., totE = 0., totW = 0., totNoRow1 = 0.;
      double rowE[17] = {0}, rowW[17] = {0};
      const int nHits = dst->numberOfEpdHits();
      for(int h = 0; h < nHits; h++){
        StPicoEpdHit* hit = dst->epdHit(h);
        if(!hit) continue;
        const double m = truncMip(hit->nMIP());
        if(m <= 0.) continue;
        const int row = hit->row();
        tot += m;
        if(hit->side() < 0){ totE += m; if(row >= 1 && row <= 16) rowE[row] += m; }
        else               { totW += m; if(row >= 1 && row <= 16) rowW[row] += m; }
        if(row != 1) totNoRow1 += m;
      }

      const int    rm  = ev->refMult();
      const int    tof = ev->nBTOFMatch();
      const double zdc = ev->ZDCx()/1000.0;   // Hz -> kHz

      if(isOff){
        hMipVzOff->Fill(vz, tot);
        hMipRefOff->Fill(rm, tot);
      }
      if(!isMine) continue;

      hMipVz->Fill(vz, tot);
      hMipTof->Fill(tof, tot);
      if(TMath::Abs(vz) <= kVzLumi) hMipZdc->Fill(zdc, tot);
      const int ri = runIndex(ev->runId());
      pMipRun->Fill(ri, tot); pHitRun->Fill(ri, nHits);

      if(!inBand(tof, rm)) continue;
      hEvents->Fill(3.5);
      if(tot <= 0.) hEvents->Fill(4.5);
      hMipRef->Fill(rm, tot);
      h3->Fill(tot, vz, zdc);
      hEW->Fill(totW, totE);
      hNoRow1->Fill(vz, totNoRow1);
      for(int r = 1; r <= 16; r++){ pRowE->Fill(rm, r, rowE[r]); pRowW->Fill(rm, r, rowW[r]); }
      nPass++;
    }
    reader->Finish();
    delete reader;
    if(f % 10 == 0) std::cout << "file " << f+1 << "/" << files.size() << "  860003 events in band so far: " << nPass << std::endl;
  }

  out->Write();
  std::cout << "\n860003 events: " << hEvents->GetBinContent(1)
            << ", |Vz|<=30: " << hEvents->GetBinContent(2)
            << ", Vr<=1: " << hEvents->GetBinContent(3)
            << ", refMult band: " << hEvents->GetBinContent(4)
            << ", of which totnMIP = 0: " << hEvents->GetBinContent(5) << std::endl;
  out->Close();
  std::cout << "Wrote " << outName << std::endl;
}
