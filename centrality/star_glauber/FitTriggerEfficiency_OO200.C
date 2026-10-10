// FitTriggerEfficiency_OO200.C -- lightflavorspectra_OO200/centrality/star_glauber
//
// Step 4 of the StRefMultCorr-style centrality for O+O 200 GeV (Run 21), trigger 860003:
// the trigger/vertex efficiency at low multiplicity, eps(m) = data / (NBD x Glauber MC), from the
// Ratio file the NBD fit writes for the best point (hRefMult = data, hRefMultSim = MC normalized
// above the fit threshold, hRatio = MC/data). The event weight for centrality-integrated or
// peripheral quantities is w(m) = 1/eps(m) (StRefMultCorr::getWeight, StCentrality::GetReweighting).
//
// Two functional forms are fitted for m below the fit threshold:
//   (A) eps = 1 - exp(-p0 * m^p1)            the StCentrality form (lowrw/highrw use its +-2 sigma);
//                                            eps(0) = 0 by construction, so it is fitted for m >= 1
//   (B) eps = 1 - p0 * exp(-p1 * m^p2)       allows eps(0) > 0, which O+O needs (data/MC ~ 0.7 at 0)
// Their difference is a systematic on the weights. Parameters +-2 sigma are printed for (A) in the
// same convention as Init_*() in StCentrality.cxx (lowrw: p0 + 2s, p1 - 2s; highrw: p0 - 2s, p1 + 2s).
//
// Usage (container, SL24y, from centrality/star_glauber):
//   root -l -b -q 'FitTriggerEfficiency_OO200.C+("Ratio_best_mc15/Ratio_npp2.750_k2.500_x0.170_eff0.060.root", 15)'
// The second argument is the NBD fit threshold (refMult >= multCut); eps is only meaningful below it.
// Outputs: OO200_TrigEff_params.txt, OO200_TrigEff.png. Written for ROOT 5 (C++98).

#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"

void FitTriggerEfficiency_OO200(const char* ratioFile = "Ratio_best_mc15/Ratio_npp2.750_k2.500_x0.170_eff0.060.root",
                                int multCut = 15, int plotMax = 40){
  gStyle->SetOptStat(0);
  TFile* in = TFile::Open(ratioFile);
  if(!in || in->IsZombie()){ std::cout << "Cannot open " << ratioFile << std::endl; return; }
  TH1D* hData = (TH1D*) in->Get("hRefMult");
  TH1D* hSim  = (TH1D*) in->Get("hRefMultSim");
  if(!hData || !hSim){ std::cout << "hRefMult / hRefMultSim missing in " << ratioFile << std::endl; return; }

  // eps(m) = data/MC per unit bin, binomial-like error from the two histograms
  TGraphErrors* g = new TGraphErrors(); g->SetName("gTrigEff");
  for(int i = 1; i <= hData->GetNbinsX(); i++){
    double m = hData->GetXaxis()->GetBinLowEdge(i);     // integer multiplicity
    if(m > plotMax) break;
    double d = hData->GetBinContent(i), s = hSim->GetBinContent(i);
    double de = hData->GetBinError(i),  se = hSim->GetBinError(i);
    if(d <= 0 || s <= 0) continue;
    double e = d/s;
    double ee = e*TMath::Sqrt((de*de)/(d*d) + (se*se)/(s*s));
    int n = g->GetN(); g->SetPoint(n, m, e); g->SetPointError(n, 0, ee);
  }
  if(g->GetN() < 5){ std::cout << "Too few points" << std::endl; return; }

  // (A) StCentrality form, m >= 1
  TF1* fA = new TF1("fA", "1 - TMath::Exp(-[0]*TMath::Power(x,[1]))", 1, multCut);
  fA->SetParameters(1.0, 0.5);
  g->Fit(fA, "Q0R");
  // (B) form with eps(0) > 0
  TF1* fB = new TF1("fB", "1 - [0]*TMath::Exp(-[1]*TMath::Power(x,[2]))", 0, multCut);
  fB->SetParameters(0.3, 0.3, 1.0);
  fB->SetParLimits(0, 0., 1.); fB->SetParLimits(2, 0.1, 3.);
  g->Fit(fB, "Q0R");

  std::ofstream par("OO200_TrigEff_params.txt");
  par << "# Trigger/vertex efficiency eps(m) = data / (NBD x Glauber MC), O+O 200 GeV, trigger 860003\n";
  par << "# input: " << ratioFile << ", fit threshold refMult >= " << multCut << "\n";
  par << "# event weight w(m) = 1/eps(m) for m < " << multCut << ", w = 1 above\n\n";
  par << "[A] eps = 1 - exp(-p0 * m^p1), fitted for 1 <= m < " << multCut << "\n";
  par << Form("p0 = %.4f +- %.4f\np1 = %.4f +- %.4f\nchi2/ndf = %.1f/%d\n", fA->GetParameter(0), fA->GetParError(0),
              fA->GetParameter(1), fA->GetParError(1), fA->GetChisquare(), fA->GetNDF());
  par << Form("lowrw  (p0 + 2s, p1 - 2s): %.4f, %.4f\n", fA->GetParameter(0) + 2*fA->GetParError(0), fA->GetParameter(1) - 2*fA->GetParError(1));
  par << Form("highrw (p0 - 2s, p1 + 2s): %.4f, %.4f\n", fA->GetParameter(0) - 2*fA->GetParError(0), fA->GetParameter(1) + 2*fA->GetParError(1));
  par << "\n[B] eps = 1 - p0 * exp(-p1 * m^p2), fitted for 0 <= m < " << multCut << "\n";
  par << Form("p0 = %.4f +- %.4f\np1 = %.4f +- %.4f\np2 = %.4f +- %.4f\nchi2/ndf = %.1f/%d\n",
              fB->GetParameter(0), fB->GetParError(0), fB->GetParameter(1), fB->GetParError(1),
              fB->GetParameter(2), fB->GetParError(2), fB->GetChisquare(), fB->GetNDF());
  par << "\n# m   eps(data/MC)   err   fitA   fitB\n";
  for(int i = 0; i < g->GetN(); i++){
    double m = g->GetX()[i];
    par << Form("%3.0f  %.4f  %.4f  %.4f  %.4f\n", m, g->GetY()[i], g->GetEY()[i], m >= 1 ? fA->Eval(m) : 0., fB->Eval(m));
  }
  // fraction of events lost below the threshold: sum(MC) - sum(data) over m < multCut, relative to all MC
  double sData = 0, sSim = 0, sSimAll = hSim->Integral();
  for(int i = 1; i <= hData->GetNbinsX(); i++){
    if(hData->GetXaxis()->GetBinLowEdge(i) >= multCut) break;
    sData += hData->GetBinContent(i); sSim += hSim->GetBinContent(i);
  }
  par << Form("\nevents below threshold: data %.4g, MC %.4g -> overall efficiency of the recorded sample: %.4f\n",
              sData, sSim, sSimAll > 0 ? (sSimAll - (sSim - sData))/sSimAll : 0.);
  par.close();

  TCanvas* c = new TCanvas("c", "", 900, 600);
  g->SetMarkerStyle(20);
  g->SetTitle(Form("Trigger/vertex efficiency, O+O 200 GeV (860003);refMult;data / MC (NBD#times Glauber, normalized at refMult #geq %d)", multCut));
  g->GetYaxis()->SetRangeUser(0, 1.3);
  g->Draw("AP");
  fA->SetLineColor(kRed); fA->Draw("same");
  fB->SetLineColor(kBlue); fB->SetLineStyle(2); fB->Draw("same");
  TLine* l1 = new TLine(0, 1, plotMax, 1); l1->SetLineStyle(3); l1->Draw();
  TLine* l2 = new TLine(multCut, 0, multCut, 1.3); l2->SetLineStyle(2); l2->SetLineColor(kGray+2); l2->Draw();
  TLegend* leg = new TLegend(0.45, 0.15, 0.88, 0.35);
  leg->AddEntry(g, "data / MC", "p");
  leg->AddEntry(fA, Form("1 - exp(-%.2f m^{%.2f})", fA->GetParameter(0), fA->GetParameter(1)), "l");
  leg->AddEntry(fB, Form("1 - %.2f exp(-%.2f m^{%.2f})", fB->GetParameter(0), fB->GetParameter(1), fB->GetParameter(2)), "l");
  leg->Draw();
  c->Print("OO200_TrigEff.png");

  std::cout << "Wrote OO200_TrigEff_params.txt and OO200_TrigEff.png" << std::endl;
  std::cout << "eps(0) from data/MC: " << (g->GetN() ? g->GetY()[0] : 0) << ";  form B gives " << fB->Eval(0) << std::endl;
}
