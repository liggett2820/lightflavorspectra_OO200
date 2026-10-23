//____________________________________________________________________________________________________
// Is the EPD centrality, defined on |Vz| <= 30 cm, valid for the spectra's |Vz| <= 2 cm events?
// Applies one set of totnMIP class edges to the corrected totnMIP at |Vz| <= 30 cm (hTotnMIP) and at
// |Vz| <= 2 cm (hTotnMIPCorr_Vz2), both from macros/FitEpdCent_OO200.C, and compares the fraction of
// recorded 860003 events in each class. If the corrected distribution has the same shape at both vertex
// ranges, every ratio is 1 within statistics.
//
// Edges: lowest totnMIP_corr included in each class, central -> peripheral, 5% steps (comma separated).
// Default: fine-scan best point (24, 6, 0.20, 0) for 0-80% (StNbdFitMaker printout) and the coarse best
// point's CentCuts_OO200.C edges for 80-95% (fine-point cuts not run yet). Provisional; pass new ones.
// Usage (from the repo top, where hTotnMIPCorr_OO200.root is):
//   root -l -b -q 'macros/CompareVzCent_OO200.C+'
//   root -l -b -q 'macros/CompareVzCent_OO200.C+("hTotnMIPCorr_OO200.root","365,311,...")'
// Output: OO200_EPD_VzCent.png, OO200_EPD_VzCent.txt
//____________________________________________________________________________________________________
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <cstdio>
#include <iostream>

void CompareVzCent_OO200(const char* inName = "hTotnMIPCorr_OO200.root",
                         const char* edgeList = "364,310,269,234,203,175,151,130,111,95,82,70,60,52,44,36,28,21,14")
{
  gStyle->SetOptStat(0);
  TFile* in = TFile::Open(inName);
  if(!in || in->IsZombie()){ std::cout << "cannot open " << inName << std::endl; return; }
  TH1D* hAll = (TH1D*) in->Get("hTotnMIP");
  TH1D* hV2  = (TH1D*) in->Get("hTotnMIPCorr_Vz2");
  if(!hAll || !hV2){ std::cout << "hTotnMIP / hTotnMIPCorr_Vz2 missing" << std::endl; return; }
  hAll = (TH1D*) hAll->Clone("hAll"); hV2 = (TH1D*) hV2->Clone("hV2");
  hAll->SetDirectory(0); hV2->SetDirectory(0); in->Close();

  // parse edges
  TObjArray* tok = TString(edgeList).Tokenize(",");
  const int nE = tok->GetEntries();
  double* edge = new double[nE + 2];
  edge[0] = 1e9;                                   // upper edge of 0-5%
  for(int i = 0; i < nE; i++) edge[i+1] = ((TObjString*) tok->At(i))->GetString().Atof();
  edge[nE+1] = 0;                                  // lower edge of the last class
  const int nC = nE + 1;                           // classes
  delete tok;

  const double totA = hAll->Integral(0, hAll->GetNbinsX()+1);
  const double totB = hV2->Integral(0, hV2->GetNbinsX()+1);
  // event count for the statistical errors: the corrected histograms are filled from TH3 bins with weight =
  // number of events (split across neighbouring bins), so the integral is the event count; Sumw2-based
  // effective entries would badly undercount it
  const double nEffB = totB;

  FILE* out = fopen("OO200_EPD_VzCent.txt", "w");
  FILE* fps[2] = {stdout, out};
  double maxDev = 0, maxPull = 0;
  for(int k = 0; k < 2; k++){
    fprintf(fps[k], "Corrected totnMIP, 860003: |Vz|<=30 cm %.4g events, |Vz|<=2 cm %.4g events (%.1f%%)\n",
            totA, totB, 100.*totB/totA);
    fprintf(fps[k], "mean totnMIP: %.3f (30 cm), %.3f (2 cm)\n", hAll->GetMean(), hV2->GetMean());
    fprintf(fps[k], "class      totnMIP range     frac(30cm)  frac(2cm)   ratio   +-stat   pull\n");
  }
  for(int c = 0; c < nC; c++){
    const double lo = edge[c+1], hi = edge[c];
    double a = 0, b = 0;
    for(int ib = 0; ib <= hAll->GetNbinsX()+1; ib++){
      const double x = (ib == 0) ? -1 : (ib > hAll->GetNbinsX() ? 1e8 : hAll->GetXaxis()->GetBinLowEdge(ib));
      if(x >= lo && x < hi){ a += hAll->GetBinContent(ib); b += hV2->GetBinContent(ib); }
    }
    const double fa = a/totA, fb = b/totB;
    const double r  = (fa > 0) ? fb/fa : 0;
    // statistical error of fb (binomial with the effective |Vz|<=2 count); the 2 cm events are a subset
    // of the 30 cm ones and fa is ~7x better known, so its error is neglected
    const double efb = TMath::Sqrt(fb*(1-fb)/nEffB);
    const double er  = (fa > 0) ? efb/fa : 0;
    const double pull = (er > 0) ? (r - 1)/er : 0;
    if(TMath::Abs(r-1) > maxDev) maxDev = TMath::Abs(r-1);
    if(TMath::Abs(pull) > maxPull) maxPull = TMath::Abs(pull);
    TString rng = (c == 0) ? Form(">= %.0f", lo) : Form("[%.0f, %.0f)", lo, hi);
    for(int k = 0; k < 2; k++)
      fprintf(fps[k], "%3d-%3d%%   %-16s  %8.5f   %8.5f   %6.4f  %6.4f  %6.1f\n",
              5*c, 5*(c+1), rng.Data(), fa, fb, r, er, pull);
  }
  for(int k = 0; k < 2; k++)
    fprintf(fps[k], "largest |ratio-1| = %.4f, largest |pull| = %.1f\n", maxDev, maxPull);
  fclose(out);

  // plot: normalized shapes and their ratio, with the class edges
  TH1D* nA = (TH1D*) hAll->Clone("nA"); nA->Scale(1./totA);
  TH1D* nB = (TH1D*) hV2->Clone("nB");  nB->Scale(1./totB);
  TH1D* rr = (TH1D*) nB->Clone("rr");   rr->Divide(nA);
  for(int ib = 1; ib <= rr->GetNbinsX(); ib++){
    const double nb = hV2->GetBinContent(ib) * nEffB / totB;   // effective events in the bin
    rr->SetBinError(ib, (nb > 0) ? rr->GetBinContent(ib)/TMath::Sqrt(nb) : 0);
  }
  rr->Rebin(5); rr->Scale(1./5);
  TCanvas* cv = new TCanvas("cv", "", 900, 900);
  TPad* p1 = new TPad("p1", "", 0, 0.45, 1, 1); p1->SetBottomMargin(0.02); p1->SetLogy(); p1->Draw();
  TPad* p2 = new TPad("p2", "", 0, 0, 1, 0.45); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.2); p2->Draw();
  p1->cd();
  nA->SetTitle(""); nA->GetXaxis()->SetRangeUser(0, 650); nA->GetXaxis()->SetLabelSize(0);
  nA->GetYaxis()->SetTitle("fraction of events per unit totnMIP");
  nA->SetMinimum(1e-7); nA->SetLineColor(kBlack); nA->Draw("hist");
  nB->SetLineColor(kOrange+7); nB->Draw("hist same");
  TLegend* lg = new TLegend(0.55, 0.72, 0.88, 0.88); lg->SetBorderSize(0);
  lg->AddEntry(nA, "|V_{z}| #leq 30 cm (centrality fit)", "l"); lg->AddEntry(nB, "|V_{z}| #leq 2 cm (spectra)", "l"); lg->Draw();
  p2->cd();
  rr->GetXaxis()->SetRangeUser(0, 650); rr->SetMinimum(0.8); rr->SetMaximum(1.2);
  rr->SetMarkerStyle(20); rr->SetMarkerSize(0.5); rr->SetTitle("");
  rr->GetYaxis()->SetTitle("2 cm / 30 cm"); rr->GetXaxis()->SetTitle("totnMIP_{corr}");
  rr->GetXaxis()->SetLabelSize(0.06); rr->GetXaxis()->SetTitleSize(0.07);
  rr->GetYaxis()->SetLabelSize(0.06); rr->GetYaxis()->SetTitleSize(0.07); rr->GetYaxis()->SetTitleOffset(0.65);
  rr->Draw("p");
  TLine one; one.SetLineStyle(2); one.DrawLine(0, 1, 650, 1);
  TLine ed; ed.SetLineColor(kGray+1); ed.SetLineStyle(3);
  TLatex lab; lab.SetTextSize(0.045); lab.SetTextAlign(21);
  int nl = 0;
  for(int c = 1; c < nC; c++){
    if(c % 2 == 1) continue;                       // label every 10%
    ed.DrawLine(edge[c], 0.8, edge[c], 1.2);
    lab.DrawLatex(edge[c], (nl % 2 == 0) ? 0.82 : 0.85, Form("%d%%", 5*c));
    nl++;
  }
  cv->SaveAs("OO200_EPD_VzCent.png");
  delete[] edge;
}
