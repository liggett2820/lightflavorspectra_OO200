// PlotEpdVzCorr_OO200.C -- lightflavorspectra_OO200
//
// Redraws plot (a) of FitEpdCent_OO200.C, the EPD Vz correction, with the legend moved off the
// points and room for the y-axis title. Same input, same calculation (meanAbove, quantile and the
// |Vz| <= 2 cm reference copied from FitEpdCent_OO200.C); writes only OO200_EPD_VzCorr_v2.png and
// prints the reference and largest correction so they can be checked against
// OO200_EpdCent_params.txt (227.9, 0.35%). Nothing else is rewritten.
//
// Usage (container, SL24y, repo root):
//   root -l -b -q 'macros/PlotEpdVzCorr_OO200.C+("epd_OO200.root")'
// Written for ROOT 5 (C++98).

#include <cmath>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"

void PlotEpdVzCorr_OO200(const char* inName = "epd_OO200.root");

#ifndef __CINT__

namespace {
  const double kMipMean  = 100.;
  const double kVzCenter = 2.;
  double meanAbove(TH1* h, double xmin, double& err){
    double s = 0, sx = 0, sxx = 0;
    for(int i = 1; i <= h->GetNbinsX(); i++){
      if(h->GetBinLowEdge(i) < xmin) continue;
      double x = h->GetBinCenter(i), w = h->GetBinContent(i);
      s += w; sx += w*x; sxx += w*x*x;
    }
    if(s <= 0){ err = 0; return 0; }
    double m = sx/s; double var = sxx/s - m*m;
    err = (var > 0) ? std::sqrt(var/s) : 0;
    return m;
  }
  double quantile(TH1* h, double q){
    double tot = h->Integral(1, h->GetNbinsX());
    if(tot <= 0) return 0;
    double c = 0;
    for(int i = 1; i <= h->GetNbinsX(); i++){
      c += h->GetBinContent(i);
      if(c >= q*tot) return h->GetBinCenter(i);
    }
    return h->GetBinCenter(h->GetNbinsX());
  }
}

void PlotEpdVzCorr_OO200(const char* inName){
  gStyle->SetOptStat(0);
  TFile* in = TFile::Open(inName);
  if(!in || in->IsZombie()){ std::cout << "Cannot open " << inName << std::endl; return; }
  TH2F* hVz = (TH2F*) in->Get("hTotnMIP_vs_Vz");
  if(!hVz){ std::cout << "Missing hTotnMIP_vs_Vz in " << inName << std::endl; return; }

  const int nVz = hVz->GetNbinsX();
  std::vector<double> vzC(nVz), mean(nVz), meanErr(nVz), q99(nVz);
  double sRef = 0, wRef = 0, sQ = 0, nQ = 0;
  for(int i = 1; i <= nVz; i++){
    TH1D* p = hVz->ProjectionY(Form("pVz_%d", i), i, i);
    vzC[i-1] = hVz->GetXaxis()->GetBinCenter(i);
    double e; mean[i-1] = meanAbove(p, kMipMean, e); meanErr[i-1] = e;
    q99[i-1] = quantile(p, 0.99);
    if(TMath::Abs(vzC[i-1]) <= kVzCenter && e > 0){ sRef += mean[i-1]/(e*e); wRef += 1./(e*e); sQ += q99[i-1]; nQ++; }
    delete p;
  }
  const double refCenter = (wRef > 0) ? sRef/wRef : 0;
  const double q99Center = (nQ > 0) ? sQ/nQ : 0;
  TGraphErrors* gVz = new TGraphErrors();
  TGraph* gQ = new TGraph();
  double maxDev = 0, yLo = 1, yHi = 1;
  for(int i = 0; i < nVz; i++){
    if(mean[i] > 0){
      const double c = refCenter/mean[i], ce = c*meanErr[i]/mean[i];
      if(TMath::Abs(c - 1.) > maxDev) maxDev = TMath::Abs(c - 1.);
      int n = gVz->GetN(); gVz->SetPoint(n, vzC[i], c); gVz->SetPointError(n, 0, ce);
      if(c - ce < yLo) yLo = c - ce;
      if(c + ce > yHi) yHi = c + ce;
    }
    if(q99[i] > 0){
      const double c = q99Center/q99[i];
      int n = gQ->GetN(); gQ->SetPoint(n, vzC[i], c);
      if(c < yLo) yLo = c;
      if(c > yHi) yHi = c;
    }
  }
  TF1* fVz = new TF1("fVz", "pol6", -30, 30);
  gVz->Fit(fVz, "Q0");
  std::cout << Form("reference <totnMIP> (|Vz| <= %.0f cm, totnMIP >= %.0f) = %.2f, largest |c_vz - 1| = %.3f%%",
                    kVzCenter, kMipMean, refCenter, 100.*maxDev) << std::endl;

  TCanvas* c1 = new TCanvas("c1", "", 900, 600);
  c1->SetLeftMargin(0.15); c1->SetRightMargin(0.04); c1->SetBottomMargin(0.12);
  gVz->SetMarkerStyle(20);
  gVz->SetTitle("EPD V_{z} correction;V_{z} (cm);c_{vz} = #LTtotnMIP#GT_{center} / #LTtotnMIP#GT(V_{z})");
  // headroom above the points; the legend goes in the empty lower middle (|Vz| < 15 cm, below the 99th-percentile row)
  const double span = yHi - yLo;
  gVz->SetMinimum(yLo - 0.10*span); gVz->SetMaximum(yHi + 0.10*span);
  gVz->GetXaxis()->SetLimits(-33, 33);
  gVz->GetYaxis()->SetTitleOffset(1.6); gVz->GetYaxis()->SetDecimals();
  gVz->Draw("AP"); fVz->SetLineColor(kRed); fVz->Draw("same");
  gQ->SetMarkerStyle(24); gQ->SetMarkerColor(kBlue); gQ->Draw("P same");
  TLegend* l1 = new TLegend(0.36, 0.15, 0.74, 0.33);
  l1->SetBorderSize(0); l1->SetFillStyle(0);
  l1->AddEntry(gVz, Form("mean, totnMIP #geq %.0f", kMipMean), "p"); l1->AddEntry(fVz, "pol6 (guide only)", "l");
  l1->AddEntry(gQ, "99th percentile", "p"); l1->Draw();
  c1->Print("OO200_EPD_VzCorr_v2.png");
}

#endif
