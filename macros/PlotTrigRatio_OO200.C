// PlotTrigRatio_OO200.C -- lightflavorspectra_OO200
//
// Efficiency of trigger 860003 (minimum bias + eTOF) relative to the official O+O minimum-bias
// triggers 860001/860002, as a function of EPD totnMIP. In the inputs every 860003 event (|Vz| <= 30 cm,
// Vr <= 1 cm) also fired 860001 or 860002 (hTrigOverlap), so 860003 is a strict subset and
// N(860003)/N(860001|860002) per totnMIP bin is a binomial efficiency (Clopper-Pearson errors).
// The official triggers are used only as the reference here (comparison).
//
// Input: merged output of MakeEpdCentInputs_OO200.C (hTotnMIP_vs_Vz, hTotnMIP_vs_Vz_off; raw totnMIP,
// no pile-up band, |Vz| <= 30 cm, Vr <= 1 cm for both).
// Output: OO200_EPD_TrigRatio.png, OO200_EPD_TrigRatio.root (ratio graph, both spectra) and a table of the
// ratio in the official StRefMultCorr("totnmip") centrality classes on stdout.
//
// Usage (container, SL24y, repo root): root -l -b -q 'macros/PlotTrigRatio_OO200.C++("epd_OO200.root")'
// ROOT 5 / C++98.

#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

void PlotTrigRatio_OO200(const char* inName = "epd_OO200.root");

#ifndef __CINT__
namespace {
  const int    kNOff = 16;
  const double kOffPct[kNOff] = {80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5};
  const double kOffCut[kNOff] = {37, 46, 55, 66, 77, 91, 105, 122, 141, 162, 186, 213, 242, 276, 315, 367};
}

void PlotTrigRatio_OO200(const char* inName){
  gStyle->SetOptStat(0);
  TFile* in = TFile::Open(inName);
  if(!in || in->IsZombie()){ std::cout << "Cannot open " << inName << std::endl; return; }
  TH2F* h2 = (TH2F*) in->Get("hTotnMIP_vs_Vz");
  TH2F* h2o = (TH2F*) in->Get("hTotnMIP_vs_Vz_off");
  if(!h2 || !h2o){ std::cout << "Missing hTotnMIP_vs_Vz or hTotnMIP_vs_Vz_off" << std::endl; return; }
  TH1D* p3 = h2->ProjectionY("p3");
  TH1D* po = h2o->ProjectionY("po");

  // variable binning on the 1-unit totnMIP axis (edges are integers, so they align with the input bins)
  std::vector<double> e;
  for(int x = 0; x < 60; x += 3) e.push_back(x);
  for(int x = 60; x < 200; x += 10) e.push_back(x);
  for(int x = 200; x < 400; x += 25) e.push_back(x);
  for(int x = 400; x <= 650; x += 50) e.push_back(x);
  const int nb = e.size() - 1;
  TH1D* r3 = (TH1D*) p3->Rebin(nb, "hTotnMIP_860003", &e[0]);
  TH1D* ro = (TH1D*) po->Rebin(nb, "hTotnMIP_860001or2", &e[0]);
  // sanity: 860003 must not exceed the reference in any bin
  int bad = 0;
  for(int i = 1; i <= nb; i++) if(r3->GetBinContent(i) > ro->GetBinContent(i)) bad++;
  if(bad) std::cout << "WARNING: " << bad << " bins with N(860003) > N(860001|2); not a strict subset there" << std::endl;

  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  g->SetName("gRatio_860003_over_mb");
  g->Divide(r3, ro, "cp");

  // ---- plot: top spectra (per unit totnMIP, log), bottom ratio
  TCanvas* c = new TCanvas("c", "", 1000, 850);
  TPad* top = new TPad("top", "", 0, 0.42, 1, 1);  top->SetBottomMargin(0.02); top->SetLogy(); top->Draw();
  TPad* bot = new TPad("bot", "", 0, 0, 1, 0.42);  bot->SetTopMargin(0.03); bot->SetBottomMargin(0.22); bot->Draw();
  top->cd();
  TH1D* d3 = (TH1D*) r3->Clone("d3"); d3->Scale(1., "width");
  TH1D* dO = (TH1D*) ro->Clone("dO"); dO->Scale(1., "width");
  dO->SetLineColor(kGray+2); dO->SetLineWidth(2); d3->SetLineColor(kAzure+2); d3->SetLineWidth(2);
  dO->SetTitle(";;Events per unit totnMIP");
  dO->GetXaxis()->SetLabelSize(0); dO->GetYaxis()->SetTitleSize(0.055); dO->GetYaxis()->SetLabelSize(0.05); dO->GetYaxis()->SetTitleOffset(0.85);
  dO->GetXaxis()->SetRangeUser(0, 600);
  dO->Draw("hist"); d3->Draw("hist same");
  TLegend* l = new TLegend(0.55, 0.72, 0.88, 0.88); l->SetBorderSize(0);
  l->AddEntry(dO, "860001 or 860002 (min. bias)", "l"); l->AddEntry(d3, "860003 (min. bias + eTOF)", "l"); l->Draw();
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.05);
  tx.DrawLatex(0.12, 0.92, "O+O 200 GeV, Run 21, |V_{z}| #leq 30 cm, V_{r} #leq 1 cm");

  bot->cd();
  TH1D* frame = new TH1D("frame", ";EPD totnMIP;N(860003) / N(min. bias)", 10, 0, 600);
  frame->SetMinimum(0.5); frame->SetMaximum(1.05);
  frame->GetXaxis()->SetTitleSize(0.08); frame->GetXaxis()->SetLabelSize(0.07);
  frame->GetYaxis()->SetTitleSize(0.07); frame->GetYaxis()->SetLabelSize(0.065); frame->GetYaxis()->SetTitleOffset(0.6);
  frame->GetYaxis()->SetNdivisions(505);
  frame->Draw();
  TLine one; one.SetLineStyle(2); one.SetLineColor(kGray+1); one.DrawLine(0, 1, 600, 1);
  // official class edges every 10%
  TLatex lab; lab.SetTextSize(0.05); lab.SetTextColor(kGray+2); lab.SetTextAlign(21);
  TLine edge; edge.SetLineStyle(3); edge.SetLineColor(kGray+1);
  for(int k = 0; k < kNOff; k++){
    if(((int)kOffPct[k]) % 10 != 0) continue;
    edge.DrawLine(kOffCut[k], 0.5, kOffCut[k], 1.05);
    lab.DrawLatex(kOffCut[k], 0.52, Form("%.0f%%", kOffPct[k]));
  }
  g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->SetMarkerColor(kAzure+2); g->SetLineColor(kAzure+2);
  g->Draw("P same");
  c->Print("OO200_EPD_TrigRatio.png");

  // ---- table in the official classes
  std::cout << "\nN(860003)/N(860001|860002) in official totnMIP classes (raw totnMIP, |Vz|<=30, Vr<=1)\n";
  std::cout << "class        totnMIP range     N(860003)     N(min bias)   ratio\n";
  double lo = 0, hi = 0;
  for(int k = -1; k < kNOff; k++){
    // k = -1: below the 80% edge; otherwise [cut_k, cut_{k+1})
    lo = (k < 0) ? 0. : kOffCut[k];
    hi = (k + 1 < kNOff) ? kOffCut[k+1] : 1e9;
    int b1 = p3->FindBin(lo + 1e-6), b2 = (hi > 1e8) ? p3->GetNbinsX() + 1 : p3->FindBin(hi - 1e-6);
    double n3 = p3->Integral(b1, b2), no = po->Integral(b1, b2);
    TString cls = (k < 0) ? Form(">%.0f%%", kOffPct[0]) : Form("%.0f-%.0f%%", (k + 1 < kNOff) ? kOffPct[k+1] : 0., kOffPct[k]);
    std::cout << Form("%-12s [%5.0f,%5s)   %12.0f  %12.0f   %.4f\n", cls.Data(), lo, (hi > 1e8) ? "inf" : Form("%.0f", hi), n3, no, no > 0 ? n3/no : 0.);
  }
  TFile* out = TFile::Open("OO200_EPD_TrigRatio.root", "RECREATE");
  g->Write(); r3->Write(); ro->Write(); out->Close();
  std::cout << "Wrote OO200_EPD_TrigRatio.png and OO200_EPD_TrigRatio.root" << std::endl;
}
#endif // __CINT__
