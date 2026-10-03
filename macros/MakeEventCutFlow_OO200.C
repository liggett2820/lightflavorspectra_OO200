// MakeEventCutFlow_OO200.C -- lightflavorspectra_OO200
//
// Event bookkeeping for the O+O 200 GeV (Run 21) analysis, straight from picoDsts:
//   1. Sequential cut flow: all events -> trigger 860003 -> |Vz| <= 2 cm -> Vr <= 1 cm
//      -> nBTOFMatch <= refMult + 100 -> refMult <= 260
//   2. Vz before and after the Vz cut (trigger-selected events)
//   3. TPC performance vs Vz: <refMult> and <nHitsFit> of primary tracks vs Vz, for
//      triggered events passing every cut except Vz (why the Vz window is tight)
//   4. Luminosity: ZDC coincidence rate (ZDCx) per event and per run, and <refMult> vs ZDCx
//   5. Every trigger ID present in the files, with event counts
//
// Cut values mirror macros/SetCutClass.C + source/CutClass.cxx + source/PicoBinner.cxx at
// commit 8296655 (trigger 860003, setZRange(-2,2), setRadialCut(1.0) about (0,0),
// setTofMatch(0), setUseT0(false), pile-up: nBTOFMatch > refMult + 100 (default
// m_pileUpLineCuts) or refMult > 260 (setPileUpCut(260))). If those change, change kCuts below.
//
// The output is additive: run on sub-lists in parallel and hadd the outputs, then plot with
// macros/PlotEventCutFlow_OO200.C.
//
// Usage (repo root, RCF, inside the SL7 container, SL24y):
//   root -l -b -q macros/LoadEventDisplayLibs.C 'macros/MakeEventCutFlow_OO200.C+("files.list","cutflow_part0.root")'
// Set doTracks = false to skip the track loop (much faster; drops the <nHitsFit> vs Vz profile).

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include "TChain.h"
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"

namespace {
  const unsigned int kTrigger = 860003;
  const double kVzMax   = 2.0;    // cm
  const double kVrMax   = 1.0;    // cm, about (0,0)
  const int    kTofLine = 100;    // reject nBTOFMatch > refMult + 100
  const int    kRefMax  = 260;    // reject refMult > 260
  const char*  kSteps[6] = {"All events", "Trigger 860003", "|V_{z}| #leq 2 cm",
                            "V_{r} #leq 1 cm", "nBTOFMatch #leq refMult+100", "refMult #leq 260"};
}

void MakeEventCutFlow_OO200(const char* fileList, const char* outName = "cutflow_OO200.root",
                            bool doTracks = true, int maxFiles = -1){
  std::vector<std::string> files; std::ifstream in(fileList); std::string line;
  while(std::getline(in,line)){
    std::string f = line.substr(0, line.find_first_of(" \t"));
    if(f.find(".picoDst.root") != std::string::npos) files.push_back(f);
  }
  if(files.empty()){ std::cout << "No .picoDst.root files in " << fileList << std::endl; return; }
  if(maxFiles > 0 && (int)files.size() > maxFiles) files.resize(maxFiles);

  TFile* out = TFile::Open(outName, "RECREATE");
  TH1D* hFlow = new TH1D("hCutFlow", "Sequential event cut flow;;Events", 6, 0, 6);
  for(int i=0;i<6;i++) hFlow->GetXaxis()->SetBinLabel(i+1, kSteps[i]);

  TH1D* hVzTrig  = new TH1D("hVz_trigger", "Trigger 860003, before the V_{z} cut;V_{z} (cm);Events", 800, -200, 200);
  TH1D* hVzFinal = new TH1D("hVz_final",   "After all event cuts;V_{z} (cm);Events", 800, -200, 200);
  TH1D* hVrTrig  = new TH1D("hVr_trigger", "Trigger 860003, |V_{z}| #leq 2 cm;V_{r} (cm);Events", 300, 0, 3);
  TH2D* hVxVy    = new TH2D("hVxVy_trigger", "Trigger 860003, |V_{z}| #leq 2 cm;V_{x} (cm);V_{y} (cm)", 300, -1.5, 1.5, 300, -1.5, 1.5);

  // TPC performance vs Vz: triggered events passing every cut except Vz
  TProfile* pRefVz  = new TProfile("pRefMult_vs_Vz",  "All cuts except V_{z};V_{z} (cm);#LTrefMult#GT", 160, -200, 200);
  TProfile* pHitsVz = new TProfile("pNHitsFit_vs_Vz", "Primary tracks |#eta| < 0.5, p_{T} > 0.2 GeV/c, all cuts except V_{z};V_{z} (cm);#LTnHitsFit#GT", 160, -200, 200);
  TProfile* pTrkVz  = new TProfile("pNGoodTracks_vs_Vz", "Tracks passing the analysis track cuts, |#eta| < 0.5;V_{z} (cm);#LTN_{tracks}#GT", 160, -200, 200);

  // Luminosity
  TH1D* hZdc        = new TH1D("hZDCx_final", "After all event cuts;ZDC coincidence rate (kHz);Events", 400, 0, 400);
  TProfile* pRefZdc = new TProfile("pRefMult_vs_ZDCx", "After all event cuts;ZDC coincidence rate (kHz);#LTrefMult#GT", 80, 0, 400);
  std::map<int, double> runEvents, runZdcSum;      // per run, after all cuts
  std::map<unsigned int, double> trigCounts;       // every trigger ID, all events

  for(size_t f = 0; f < files.size(); f++){
    StPicoDstReader* reader = new StPicoDstReader(files[f].c_str());
    reader->Init();
    reader->SetStatus("*", 0);
    reader->SetStatus("Event", 1);
    if(doTracks) reader->SetStatus("Track", 1);
    Long64_t n = reader->chain()->GetEntries();
    for(Long64_t i = 0; i < n; i++){
      if(!reader->readPicoEvent(i)) continue;
      StPicoDst* dst = reader->picoDst();
      StPicoEvent* ev = dst->event();
      if(!ev) continue;

      hFlow->Fill(0.5);
      std::vector<unsigned int> ids = ev->triggerIds();
      for(size_t t = 0; t < ids.size(); t++) trigCounts[ids[t]] += 1;
      if(!ev->isTrigger(kTrigger)) continue;
      hFlow->Fill(1.5);

      TVector3 pv = ev->primaryVertex();
      const double vz = pv.Z();
      const double vr = TMath::Sqrt(pv.X()*pv.X() + pv.Y()*pv.Y());
      const int rm = ev->refMult();
      const int tofm = ev->nBTOFMatch();
      const bool passVr = (vr <= kVrMax);
      const bool passPile = (tofm <= rm + kTofLine) && (rm <= kRefMax);
      hVzTrig->Fill(vz);

      // performance vs Vz: everything except the Vz cut
      if(passVr && passPile){
        pRefVz->Fill(vz, rm);
        if(doTracks){
          int nGood = 0;
          for(UInt_t t = 0; t < dst->numberOfTracks(); t++){
            StPicoTrack* tr = dst->track(t);
            if(!tr || !tr->isPrimary()) continue;
            TVector3 p = tr->pMom();
            if(TMath::Abs(p.Eta()) > 0.5 || p.Perp() < 0.2) continue;
            pHitsVz->Fill(vz, tr->nHitsFit());
            if(tr->nHitsFit() >= 10 && tr->nHitsDedx() >= 10 && tr->nHitsPoss() > 0 &&
               ((double)tr->nHits())/((double)tr->nHitsPoss()) >= 0.52 && tr->gDCA(pv).Mag() <= 1.0) nGood++;
          }
          pTrkVz->Fill(vz, nGood);
        }
      }

      if(TMath::Abs(vz) > kVzMax) continue;
      hFlow->Fill(2.5);
      hVrTrig->Fill(vr); hVxVy->Fill(pv.X(), pv.Y());
      if(!passVr) continue;
      hFlow->Fill(3.5);
      if(tofm > rm + kTofLine) continue;
      hFlow->Fill(4.5);
      if(rm > kRefMax) continue;
      hFlow->Fill(5.5);

      hVzFinal->Fill(vz);
      const double zdc = ev->ZDCx()/1000.0; // Hz -> kHz
      hZdc->Fill(zdc); pRefZdc->Fill(zdc, rm);
      runEvents[ev->runId()] += 1; runZdcSum[ev->runId()] += zdc;
    }
    reader->Finish();
    delete reader;
    if(f % 50 == 0) std::cout << "file " << f+1 << "/" << files.size() << "  events so far: " << hFlow->GetBinContent(1) << std::endl;
  }

  // per-run tree (additive under hadd)
  out->cd();
  TTree* tRun = new TTree("runs", "Per-run counts after all event cuts");
  int run; double nEv, zdcSum;
  tRun->Branch("run", &run, "run/I"); tRun->Branch("nEvents", &nEv, "nEvents/D"); tRun->Branch("zdcSum_kHz", &zdcSum, "zdcSum_kHz/D");
  for(std::map<int,double>::iterator it = runEvents.begin(); it != runEvents.end(); ++it){
    run = it->first; nEv = it->second; zdcSum = runZdcSum[run]; tRun->Fill();
  }
  TTree* tTrig = new TTree("triggers", "Every trigger ID in the files, event counts");
  unsigned int tid; double tN;
  tTrig->Branch("trigId", &tid, "trigId/i"); tTrig->Branch("nEvents", &tN, "nEvents/D");
  for(std::map<unsigned int,double>::iterator it = trigCounts.begin(); it != trigCounts.end(); ++it){
    tid = it->first; tN = it->second; tTrig->Fill();
  }
  out->Write();
  std::cout << "\nCut flow (" << files.size() << " files):" << std::endl;
  for(int i=1;i<=6;i++) std::cout << Form("  %-32s %14.0f", hFlow->GetXaxis()->GetBinLabel(i), hFlow->GetBinContent(i)) << std::endl;
  out->Close();
  std::cout << "Wrote " << outName << std::endl;
}
