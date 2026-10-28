// PlotZdcSmd_OO200.C -- lightflavorspectra_OO200
//
// First look at the ZDC shower-max detector (SMD) histograms from macros/MakeZdcSpectatorNeutrons_OO200.C
// (Oct 28 2026 version). Plain ROOT, no STAR libraries. Three pictures:
//   OO200_Smd_strips.png  mean ADC per strip, per side and plane, for 0n / 1n / 2n+ events
//                         (0n shows the pedestal level; 1n should show a bump where the neutron hits)
//   OO200_Smd_xy.png      shower centroid y vs x (strip units, raw ADC weights) for 1n and 2n+, per side
//   OO200_Smd_run.png     mean centroid x and y per run, 1n events (beam-centre stability)
// and prints the mean 0n ADC per strip (pedestal estimate) and the 1n centroid means.
//
// Usage: root -l -b -q 'macros/PlotZdcSmd_OO200.C("zdc_smdtest.root")'
// Written for ROOT 5 (C++98).

#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"

void PlotZdcSmd_OO200(const char* inName = "zdc_smdtest.root"){
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open(inName);
  if(!f || f->IsZombie()){ std::cout << "cannot open " << inName << std::endl; return; }
  const char* side[2]  = {"East", "West"};
  const char* plane[2] = {"Vertical", "Horizontal"};
  const char* gate[3]  = {"0n", "1n", "2n"};
  const int   col[3]   = {kGray + 2, kBlue + 1, kRed + 1};

  // 1. mean ADC per strip
  TCanvas* c1 = new TCanvas("c1", "", 1400, 900); c1->Divide(2, 2);
  for(int s = 0; s < 2; s++) for(int p = 0; p < 2; p++){
    c1->cd(1 + 2 * s + p);
    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    double ymax = 0; TProfile* pr[3];
    for(int g = 0; g < 3; g++){
      TH2F* h = (TH2F*) f->Get(Form("hSmd%s%s_%s", side[s], plane[p], gate[g]));
      pr[g] = h ? h->ProfileX(Form("pr_%d_%d_%d", s, p, g)) : 0;
      if(pr[g] && pr[g]->GetMaximum() > ymax) ymax = pr[g]->GetMaximum();
    }
    for(int g = 0; g < 3; g++){
      if(!pr[g]) continue;
      pr[g]->SetTitle(Form("ZDC-SMD %s %s;strip;mean ADC", side[s], plane[p]));
      pr[g]->SetLineColor(col[g]); pr[g]->SetMarkerColor(col[g]); pr[g]->SetMarkerStyle(20); pr[g]->SetLineWidth(2);
      pr[g]->SetMinimum(0); pr[g]->SetMaximum(1.2 * ymax);
      pr[g]->Draw(g == 0 ? "E1" : "E1 same");
      leg->AddEntry(pr[g], g == 2 ? "2n+" : gate[g], "lp");
      if(g == 0){
        std::cout << Form("0n mean ADC, %s %s:", side[s], plane[p]);
        for(int k = 1; k <= 8; k++) std::cout << Form(" %6.1f", pr[g]->GetBinContent(k));
        std::cout << std::endl;
      }
    }
    leg->Draw();
  }
  c1->SaveAs("OO200_Smd_strips.png");

  // 2. centroid maps
  TCanvas* c2 = new TCanvas("c2", "", 1400, 1200); c2->Divide(2, 2);
  for(int s = 0; s < 2; s++){
    TH2F* h1 = 0;
    for(int c = 0; c < 6; c++){
      TH2F* h = (TH2F*) f->Get(Form("hSmdXY%s_1n_cent%d", side[s], c));
      if(!h) continue;
      if(!h1){ h1 = (TH2F*) h->Clone(Form("xy1n_%d", s)); h1->SetTitle(Form("ZDC-SMD %s centroid, 1n, all classes;x (vertical strip index);y (horizontal strip index)", side[s])); }
      else h1->Add(h);
    }
    c2->cd(1 + s); gPad->SetRightMargin(0.13); if(h1) h1->Draw("colz");
    if(h1) std::cout << Form("1n centroid %s: <x> = %.2f (rms %.2f), <y> = %.2f (rms %.2f) strips, %.0f events",
                             side[s], h1->GetMean(1), h1->GetRMS(1), h1->GetMean(2), h1->GetRMS(2), h1->GetEntries()) << std::endl;
    TH2F* h2 = (TH2F*) f->Get(Form("hSmdXY%s_2n", side[s]));
    c2->cd(3 + s); gPad->SetRightMargin(0.13); if(h2) h2->Draw("colz");
  }
  c2->SaveAs("OO200_Smd_xy.png");

  // 3. centroid vs run
  TCanvas* c3 = new TCanvas("c3", "", 1400, 900); c3->Divide(2, 2);
  const char* xy[2] = {"X", "Y"};
  for(int s = 0; s < 2; s++) for(int a = 0; a < 2; a++){
    c3->cd(1 + 2 * s + a);
    TH2F* h = (TH2F*) f->Get(Form("hSmd%s%s_vs_run_1n", xy[a], side[s]));
    if(!h) continue;
    TProfile* pr = h->ProfileX(Form("prRun_%d_%d", s, a));
    pr->SetTitle(Form("ZDC-SMD %s mean centroid %s, 1n;run index = (day-130)*100 + run-in-day;<%s> (strip index)", side[s], a == 0 ? "x" : "y", a == 0 ? "x" : "y"));
    pr->SetMarkerStyle(20); pr->SetMarkerSize(0.7);
    pr->GetXaxis()->SetRangeUser(320, 620); pr->SetMinimum(-0.5); pr->SetMaximum(7.5);
    pr->Draw("E1");
  }
  c3->SaveAs("OO200_Smd_run.png");
}
