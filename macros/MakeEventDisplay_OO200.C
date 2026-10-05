// MakeEventDisplay_OO200.C -- lightflavorspectra_OO200
//
// Draws one central O+O event (Run 21, sqrt(sNN) = 200 GeV) from picoDsts: primary-track
// helices in the transverse (x-y) view and a side (z-y) view, with the TPC outline.
// Tracks are reconstructed from each primary track's momentum, charge and the event's
// magnetic field (helix from the primary vertex out to the TPC outer radius), so this is a
// track-level display, not a hit-level one.
//
// Event choice: the analysis event cuts (trigger 860003, |Vz| <= 2 cm, Vr <= 1 cm, pile-up
// rejection refMult <= 260 and nBTOFMatch <= refMult + 100), 0-5% centrality (refMult >= 44),
// then the event with the most good tracks among the first `maxEvents` scanned.
// Track cuts: primary, nHitsFit >= 10, nHitsDedx >= 10, nHitsFit/nHitsPoss >= 0.52, gDCA <= 1 cm.
//
// Usage (from the repo root, RCF build):
//   root -l -b -q macros/LoadEventDisplayLibs.C 'macros/MakeEventDisplay_OO200.C+("filelist.list")'
// Writes OO200_EventDisplay_xy.png (square end view, for the outline slide),
// OO200_EventDisplay_xy_zy.png (both views) and prints the chosen run/event ID.

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPolyLine.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "TColor.h"
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"

namespace {
  const double kTpcRin  = 50.0;   // cm, TPC inner field cage (approx.)
  const double kTpcRout = 200.0;  // cm, TPC outer field cage (approx.)
  const double kTpcHalfZ = 210.0; // cm, TPC half-length (approx.)

  struct TrackPath { std::vector<double> x, y, z; int q; };

  bool goodEvent(StPicoEvent* ev){
    if(!ev->isTrigger(860003)) return false;
    TVector3 pv = ev->primaryVertex();
    if(TMath::Abs(pv.Z()) > 2.0) return false;
    if(TMath::Sqrt(pv.X()*pv.X()+pv.Y()*pv.Y()) > 1.0) return false;
    int rm = ev->refMult();
    if(rm > 260) return false;
    if((int)ev->nBTOFMatch() > rm + 100) return false;
    return rm >= 44; // 0-5%
  }

  bool goodTrack(StPicoTrack* t, const TVector3& pv){
    if(!t || !t->isPrimary()) return false;
    if(t->nHitsFit() < 10 || t->nHitsDedx() < 10) return false;
    if(t->nHitsPoss() <= 0 || ((double)t->nHits())/((double)t->nHitsPoss()) < 0.52) return false;
    if(t->gDCA(pv).Mag() > 1.0) return false;
    return true;
  }

  // Helix from the vertex, stepped in turning angle until it leaves the TPC.
  TrackPath helix(const TVector3& p, int q, const TVector3& v, double bTesla){
    TrackPath path; path.q = q;
    double pt = p.Perp();
    double R = (bTesla != 0) ? pt/(0.299792458*TMath::Abs(bTesla))*100.0 : 1e9; // cm
    double h = (q*bTesla > 0) ? -1.0 : 1.0; // sense of rotation
    double phi0 = TMath::ATan2(p.Y(), p.X());
    double dzda = R*p.Z()/pt;               // dz per radian of turning
    double x = v.X(), y = v.Y(), z = v.Z();
    for(double a = 0; a < TMath::Pi(); a += 0.002){
      x = v.X() + R*h*(TMath::Sin(phi0 + h*a) - TMath::Sin(phi0));
      y = v.Y() - R*h*(TMath::Cos(phi0 + h*a) - TMath::Cos(phi0));
      z = v.Z() + dzda*a;
      if(TMath::Sqrt(x*x+y*y) > kTpcRout || TMath::Abs(z) > kTpcHalfZ) break;
      path.x.push_back(x); path.y.push_back(y); path.z.push_back(z);
    }
    return path;
  }

  void drawFrame(bool sideView){
    static int nFrames = 0;
    TH2F* fr = new TH2F(Form("evdispFrame_%d",nFrames++),"",10,-230,230,10,-230,230);
    fr->SetStats(0);
    fr->GetXaxis()->SetAxisColor(0); fr->GetYaxis()->SetAxisColor(0);
    fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetLabelSize(0);
    fr->GetXaxis()->SetTickLength(0); fr->GetYaxis()->SetTickLength(0);
    fr->Draw("A");
    int rim = TColor::GetColor("#5B6F8C");
    if(!sideView){
      TEllipse* o = new TEllipse(0,0,kTpcRout,kTpcRout); o->SetFillStyle(0); o->SetLineColor(rim); o->SetLineWidth(2); o->Draw();
      TEllipse* i = new TEllipse(0,0,kTpcRin,kTpcRin);   i->SetFillStyle(0); i->SetLineColor(rim); i->SetLineWidth(2); i->Draw();
      for(int s=0;s<12;s++){ // sector boundaries
        double ph = (s*30.0+15.0)*TMath::DegToRad();
        TLine* l = new TLine(kTpcRin*TMath::Cos(ph),kTpcRin*TMath::Sin(ph),kTpcRout*TMath::Cos(ph),kTpcRout*TMath::Sin(ph));
        l->SetLineColor(rim); l->SetLineStyle(3); l->Draw();
      }
    } else {
      TBox* b1 = new TBox(-kTpcHalfZ,kTpcRin,kTpcHalfZ,kTpcRout);   b1->SetFillStyle(0); b1->SetLineColor(rim); b1->SetLineWidth(2); b1->Draw();
      TBox* b2 = new TBox(-kTpcHalfZ,-kTpcRout,kTpcHalfZ,-kTpcRin); b2->SetFillStyle(0); b2->SetLineColor(rim); b2->SetLineWidth(2); b2->Draw();
      TLine* m = new TLine(0,-kTpcRout,0,kTpcRout); m->SetLineColor(rim); m->SetLineStyle(3); m->Draw();
    }
  }
}

void MakeEventDisplay_OO200(const char* fileList, int maxEvents = 5000, const char* outPrefix = "OO200_EventDisplay"){
  gStyle->SetOptStat(0);
  // file-by-file reading (one reader per file), as in PicoBinner
  std::vector<std::string> files; std::ifstream in(fileList); std::string line;
  while(std::getline(in,line)){ if(line.find(".root") != std::string::npos) files.push_back(line); }
  if(files.empty()){ std::cout << "No .root files in " << fileList << std::endl; return; }

  int bestN = -1, bestRun = 0, bestEvt = 0, bestRefMult = 0; double bestB = 0;
  TVector3 bestPV; std::vector<TVector3> bestP; std::vector<int> bestQ;
  int scanned = 0;
  for(size_t f = 0; f < files.size() && scanned < maxEvents; f++){
    StPicoDstReader* reader = new StPicoDstReader(files[f].c_str());
    reader->Init();
    Long64_t n = reader->chain()->GetEntries();
    for(Long64_t i = 0; i < n && scanned < maxEvents; i++, scanned++){
      if(!reader->readPicoEvent(i)) continue;
      StPicoDst* dst = reader->picoDst();
      StPicoEvent* ev = dst->event();
      if(!ev || !goodEvent(ev)) continue;
      TVector3 pv = ev->primaryVertex();
      std::vector<TVector3> P; std::vector<int> Q;
      for(UInt_t t = 0; t < dst->numberOfTracks(); t++){
        StPicoTrack* tr = dst->track(t);
        if(!goodTrack(tr,pv)) continue;
        P.push_back(tr->pMom()); Q.push_back(tr->charge());
      }
      if((int)P.size() > bestN){
        bestN = P.size(); bestP = P; bestQ = Q; bestPV = pv;
        bestRun = ev->runId(); bestEvt = ev->eventId(); bestRefMult = ev->refMult();
        bestB = ev->bField()/10.0; // kG -> T
      }
    }
    delete reader;
  }
  if(bestN < 0){ std::cout << "No 0-5% event passing cuts in the first " << scanned << " events." << std::endl; return; }
  std::cout << "Chosen event: run " << bestRun << ", event " << bestEvt << ", refMult " << bestRefMult
            << ", good tracks " << bestN << ", B = " << bestB << " T (scanned " << scanned << " events)" << std::endl;

  std::vector<TrackPath> paths;
  for(size_t k = 0; k < bestP.size(); k++) paths.push_back(helix(bestP[k],bestQ[k],bestPV,bestB));
  int cPos = TColor::GetColor("#FFBF00"); // Aggie Gold
  int cNeg = TColor::GetColor("#7FB8FF");
  int bg   = TColor::GetColor("#022851"); // Aggie Blue

  // square end view
  TCanvas* c1 = new TCanvas("c1","",1200,1200);
  c1->SetFillColor(bg); c1->SetFrameFillColor(bg); c1->SetFrameLineColor(bg);
  c1->SetMargin(0.02,0.02,0.02,0.02);
  drawFrame(false);
  for(size_t k = 0; k < paths.size(); k++){
    if(paths[k].x.size() < 2) continue;
    TPolyLine* pl = new TPolyLine(paths[k].x.size(), &paths[k].x[0], &paths[k].y[0]);
    pl->SetLineColor(paths[k].q > 0 ? cPos : cNeg); pl->SetLineWidth(2); pl->Draw();
  }
  TLatex lab; lab.SetTextColor(0); lab.SetTextFont(42); lab.SetTextSize(0.028);
  lab.DrawLatexNDC(0.04,0.95,"STAR  O+O  #sqrt{s_{NN}} = 200 GeV");
  lab.DrawLatexNDC(0.04,0.91,Form("0-5%%, refMult = %d, %d primary tracks", bestRefMult, bestN));
  c1->SaveAs(Form("%s_xy.png",outPrefix));

  // both views
  TCanvas* c2 = new TCanvas("c2","",2000,1000);
  c2->SetFillColor(bg);
  c2->Divide(2,1,0.001,0.001);
  for(int v = 1; v <= 2; v++){
    c2->cd(v); gPad->SetFillColor(bg); gPad->SetFrameFillColor(bg); gPad->SetFrameLineColor(bg); gPad->SetMargin(0.02,0.02,0.02,0.02);
    drawFrame(v == 2);
    for(size_t k = 0; k < paths.size(); k++){
      if(paths[k].x.size() < 2) continue;
      TPolyLine* pl = (v == 1) ? new TPolyLine(paths[k].x.size(), &paths[k].x[0], &paths[k].y[0])
                               : new TPolyLine(paths[k].z.size(), &paths[k].z[0], &paths[k].y[0]);
      pl->SetLineColor(paths[k].q > 0 ? cPos : cNeg); pl->SetLineWidth(2); pl->Draw();
    }
  }
  c2->cd(1); lab.DrawLatexNDC(0.04,0.95,"STAR  O+O  #sqrt{s_{NN}} = 200 GeV");
  lab.DrawLatexNDC(0.04,0.91,Form("0-5%%, refMult = %d, %d primary tracks", bestRefMult, bestN));
  c2->cd(2); lab.DrawLatexNDC(0.04,0.95,"side view (z-y)");
  c2->SaveAs(Form("%s_xy_zy.png",outPrefix));
}
