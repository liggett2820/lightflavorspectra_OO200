// PlotZdcSpectatorNeutrons_OO200.C -- lightflavorspectra_OO200
//
// Reads the hadd-ed output of macros/MakeZdcSpectatorNeutrons_OO200.C and
//   1. fits the single-neutron (1n) peak of the unattenuated ZDC ADC, per side, in the 80-100%
//      class (pedestal tail + 1n Gaussian + 2n Gaussian at twice the 1n mean, sqrt(2) the width),
//      and compares it with the reference from D. Wen's O+O UPC talk (East ~125, West ~111 ADC);
//   2. repeats a 1n Gaussian fit run by run (80-100%) to check gain stability;
//   3. converts <ADC> per centrality class into an approximate <free neutrons> per side,
//      <N_n> ~ <ADC> / ADC_1n (approximate: the ZDC response is not linear at high signal);
//   4. draws the per-class ADC spectra with lines at n x ADC_1n.
// Outputs: OO200_ZdcNeutrons.txt, OO200_Zdc_1nFit.png, OO200_Zdc_1nVsRun.png, OO200_Zdc_SpectraByClass.png
//
// Usage: root -l -b -q 'macros/PlotZdcSpectatorNeutrons_OO200.C+("zdc_OO200_all.root")'
//   useReference = true skips the fit and uses the reference 1n positions for step 3.

#include <fstream>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"

namespace {
  const int kNCent = 6;
  const char* kCentLab[kNCent] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
  const char* kSide[2] = {"East", "West"};
  const double kRef1n[2] = {125., 111.};   // D. Wen, O+O UPC talk, Oct 6 2026 (UPC triggers)
  const int kColor[kNCent] = {kRed+1, kOrange+7, kGreen+2, kAzure+2, kBlue+3, kGray+2};
  // pedestal tail + 1n + 2n (2n mean = 2 x 1n mean, width = sqrt(2) x 1n width)
  double zdcModel(double* x, double* p){
    const double v = x[0];
    const double ped = p[0] * TMath::Exp(-v / p[1]);
    const double g1  = p[2] * TMath::Gaus(v, p[3], p[4]);
    const double g2  = p[5] * TMath::Gaus(v, 2.0 * p[3], TMath::Sqrt(2.0) * p[4]);
    return ped + g1 + g2;
  }
}

void PlotZdcSpectatorNeutrons_OO200(const char* inName = "zdc_OO200_all.root", bool useReference = false){
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open(inName);
  if(!f || f->IsZombie()){ std::cout << "Cannot open " << inName << std::endl; return; }
  std::ofstream txt("OO200_ZdcNeutrons.txt");

  TH1D* hEv = (TH1D*)f->Get("hEvents");
  const double nAll = hEv ? hEv->GetBinContent(1) : 0;
  std::cout << "Events (860003 + event cuts): " << nAll << std::endl;
  txt << "events_passing\t" << (long long)nAll << "\n";
  if(hEv && nAll > 0){
    std::cout << Form("Unattenuated ADC == 0: East %.2f%%, West %.2f%%, both %.2f%%",
                      100*hEv->GetBinContent(2)/nAll, 100*hEv->GetBinContent(3)/nAll, 100*hEv->GetBinContent(4)/nAll) << std::endl;
    txt << "frac_UA_zero_East\t" << hEv->GetBinContent(2)/nAll << "\nfrac_UA_zero_West\t" << hEv->GetBinContent(3)/nAll << "\n";
  }

  // ---- 1. 1n peak in 80-100%
  double adc1n[2] = {kRef1n[0], kRef1n[1]};
  double adc1nErr[2] = {0, 0};
  TCanvas* c1 = new TCanvas("c1", "", 1400, 600); c1->Divide(2, 1);
  for(int s = 0; s < 2; s++){
    c1->cd(s+1); gPad->SetLogy();
    TH1D* h = (TH1D*)f->Get(Form("hZdcUA%s_cent%d", kSide[s], kNCent-1));
    if(!h || h->GetEntries() < 100){ std::cout << "Too few 80-100% entries, " << kSide[s] << std::endl; continue; }
    h->SetLineColor(kBlack); h->GetXaxis()->SetRangeUser(0, 800); h->SetTitle(Form("ZDC %s, 80-100%%", kSide[s]));
    h->Draw("hist");
    TF1* fm = new TF1(Form("fZdc%s", kSide[s]), zdcModel, 40, 600, 6);
    fm->SetParNames("A_ped", "tau_ped", "A_1n", "mu_1n", "sigma_1n", "A_2n");
    const double peak = h->GetBinContent(h->FindBin(kRef1n[s]));
    fm->SetParameters(peak, 40., peak, kRef1n[s], 50., 0.3*peak);
    fm->SetParLimits(1, 5., 300.); fm->SetParLimits(3, 0.6*kRef1n[s], 1.5*kRef1n[s]); fm->SetParLimits(4, 15., 120.);
    fm->SetParLimits(0, 0., 1e12); fm->SetParLimits(2, 0., 1e12); fm->SetParLimits(5, 0., 1e12);
    int status = h->Fit(fm, "RQN0");
    fm->SetLineColor(kRed+1); fm->Draw("same");
    const double mu = fm->GetParameter(3), muE = fm->GetParError(3), sg = fm->GetParameter(4);
    const double chi2ndf = fm->GetNDF() > 0 ? fm->GetChisquare()/fm->GetNDF() : -1;
    std::cout << Form("%s 1n fit: mu = %.1f +- %.1f ADC, sigma = %.1f, chi2/ndf = %.1f, status %d (reference %.0f)",
                      kSide[s], mu, muE, sg, chi2ndf, status, kRef1n[s]) << std::endl;
    txt << "fit_1n_" << kSide[s] << "\tmu=" << mu << "\terr=" << muE << "\tsigma=" << sg << "\tchi2ndf=" << chi2ndf << "\tstatus=" << status << "\treference=" << kRef1n[s] << "\n";
    for(int k = 1; k <= 4; k++){ TLine* l = new TLine(k*mu, 1, k*mu, h->GetMaximum()); l->SetLineStyle(2); l->SetLineColor(kRed+1); l->Draw(); }
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
    lat.DrawLatex(0.45, 0.85, Form("1n fit: #mu = %.0f #pm %.0f ADC (ref. %.0f)", mu, muE, kRef1n[s]));
    if(!useReference && status == 0 && muE > 0 && muE < 0.1*mu){ adc1n[s] = mu; adc1nErr[s] = muE; }
    else std::cout << "  -> using the reference 1n position for " << kSide[s] << std::endl;
  }
  c1->SaveAs("OO200_Zdc_1nFit.png");

  // ---- 2. 1n peak vs run (80-100%)
  TCanvas* c2 = new TCanvas("c2", "", 1400, 500);
  TLegend* lg2 = new TLegend(0.75, 0.75, 0.9, 0.9); lg2->SetBorderSize(0);
  txt << "\nrun_index\tside\tmu_1n\terr\tentries\n";
  bool first = true;
  for(int s = 0; s < 2; s++){
    TH2F* h2 = (TH2F*)f->Get(Form("hZdcUA%s_vs_run_peripheral", kSide[s]));
    if(!h2) continue;
    TGraphErrors* g = new TGraphErrors(); int np = 0;
    for(int bx = 1; bx <= h2->GetNbinsX(); bx++){
      TH1D* p = h2->ProjectionY(Form("p%s_%d", kSide[s], bx), bx, bx);
      if(p->GetEntries() < 500){ delete p; continue; }
      TF1 gfit("gfit", "gaus", adc1n[s] - 60, adc1n[s] + 60);
      gfit.SetParameters(p->GetMaximum(), adc1n[s], 40);
      if(p->Fit(&gfit, "RQN0") == 0){
        const int ri = (int)h2->GetXaxis()->GetBinLowEdge(bx);
        g->SetPoint(np, ri, gfit.GetParameter(1)); g->SetPointError(np, 0, gfit.GetParError(1)); np++;
        txt << ri << "\t" << kSide[s] << "\t" << gfit.GetParameter(1) << "\t" << gfit.GetParError(1) << "\t" << (long long)p->GetEntries() << "\n";
      }
      delete p;
    }
    if(np == 0) continue;
    g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->SetMarkerColor(s == 0 ? kAzure+2 : kOrange+7); g->SetLineColor(g->GetMarkerColor());
    g->SetTitle(";run index = (day #minus 130) #times 100 + run in day;1n peak (ADC)");
    g->Draw(first ? "AP" : "P"); first = false;
    lg2->AddEntry(g, kSide[s], "p");
  }
  lg2->Draw(); c2->SaveAs("OO200_Zdc_1nVsRun.png");

  // ---- 3. <N_n> per class
  std::cout << "\n<free neutrons> per side ~ <ADC>/ADC_1n (approximate; ZDC nonlinear at high signal)" << std::endl;
  std::cout << Form("ADC_1n used: East %.1f, West %.1f", adc1n[0], adc1n[1]) << std::endl;
  txt << "\nADC_1n_used\tEast=" << adc1n[0] << "\tWest=" << adc1n[1] << "\n";
  txt << "class\tevents\tmeanADC_East\tNn_East\tmeanADC_West\tNn_West\tNn_East_ref\tNn_West_ref\n";
  for(int c = 0; c < kNCent; c++){
    TH1D* hE = (TH1D*)f->Get(Form("hZdcUAEast_cent%d", c));
    TH1D* hW = (TH1D*)f->Get(Form("hZdcUAWest_cent%d", c));
    if(!hE || !hW) continue;
    const double mE = hE->GetMean(), mW = hW->GetMean();
    std::cout << Form("  %-8s  N_ev = %10.0f   East <ADC> = %6.1f -> <N_n> = %5.2f   West <ADC> = %6.1f -> <N_n> = %5.2f",
                      kCentLab[c], hE->GetEntries(), mE, mE/adc1n[0], mW, mW/adc1n[1]) << std::endl;
    txt << kCentLab[c] << "\t" << (long long)hE->GetEntries() << "\t" << mE << "\t" << mE/adc1n[0] << "\t" << mW << "\t" << mW/adc1n[1]
        << "\t" << mE/kRef1n[0] << "\t" << mW/kRef1n[1] << "\n";
  }
  txt.close();

  // ---- 4. spectra by class
  TCanvas* c3 = new TCanvas("c3", "", 1400, 600); c3->Divide(2, 1);
  for(int s = 0; s < 2; s++){
    c3->cd(s+1); gPad->SetLogy();
    TLegend* lg = new TLegend(0.62, 0.55, 0.89, 0.89); lg->SetBorderSize(0);
    double ymax = 0;
    std::vector<TH1D*> hs;
    for(int c = 0; c < kNCent; c++){
      TH1D* h = (TH1D*)f->Get(Form("hZdcUA%s_cent%d", kSide[s], c));
      if(!h || h->Integral() <= 0){ hs.push_back(0); continue; }
      TH1D* hn = (TH1D*)h->Clone(Form("n%s%d", kSide[s], c)); hn->Rebin(4); hn->Scale(1.0/hn->Integral());
      ymax = TMath::Max(ymax, hn->GetMaximum()); hs.push_back(hn);
    }
    bool drawn = false;
    for(int c = 0; c < kNCent; c++){
      if(!hs[c]) continue;
      hs[c]->SetLineColor(kColor[c]); hs[c]->SetLineWidth(2);
      hs[c]->SetTitle(Form("ZDC %s unattenuated ADC by centrality;ADC;normalized", kSide[s]));
      hs[c]->SetMaximum(2*ymax); hs[c]->SetMinimum(1e-6); hs[c]->GetXaxis()->SetRangeUser(0, 3000);
      hs[c]->Draw(drawn ? "hist same" : "hist"); drawn = true;
      lg->AddEntry(hs[c], kCentLab[c], "l");
    }
    for(int k = 1; k <= 20; k++){ TLine* l = new TLine(k*adc1n[s], 1e-6, k*adc1n[s], 2*ymax); l->SetLineStyle(3); l->SetLineColor(kGray); l->Draw(); }
    lg->Draw();
  }
  c3->SaveAs("OO200_Zdc_SpectraByClass.png");
  std::cout << "Wrote OO200_ZdcNeutrons.txt and OO200_Zdc_{1nFit,1nVsRun,SpectraByClass}.png" << std::endl;
}
