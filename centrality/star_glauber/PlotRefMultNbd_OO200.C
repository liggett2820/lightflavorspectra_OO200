// PlotRefMultNbd_OO200.C -- lightflavorspectra_OO200, centrality/star_glauber
// Fit plot for the refMult NBD x Glauber best point: corrected refMult (data) and the model on a log
// scale, data/model below, the fit-range start and the class edges marked.
// Input: the best point's RatioChi2Files/Ratio_npp..._k..._x..._eff....root ("hRefMult" = data,
// "hRefMultSim" = model, already normalized to data above the fit-range start by StNbdFitMaker::Fit).
// Usage: root -l -b -q 'PlotRefMultNbd_OO200.C("RatioChi2Files_x/Ratio_npp2.550_k2.500_x0.230_eff0.040.root",15,"43,35,26,14","OO200_NBD_refMult_fit.png")'
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <iostream>

void PlotRefMultNbd_OO200(const char* ratioFile, double multCut = 15, const char* edges = "43,35,26,14",
                          const char* outName = "OO200_NBD_refMult_fit.png", double xmax = 100)
{
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  TFile* in = TFile::Open(ratioFile);
  if(in == 0 || in->IsZombie()){ std::cout << "cannot open " << ratioFile << std::endl; return; }
  TH1D* hD = (TH1D*) in->Get("hRefMult");
  TH1D* hM = (TH1D*) in->Get("hRefMultSim");
  if(hD == 0 || hM == 0){ std::cout << "hRefMult / hRefMultSim missing" << std::endl; return; }
  hD = (TH1D*) hD->Clone("hData"); hM = (TH1D*) hM->Clone("hModel");
  hD->SetDirectory(0); hM->SetDirectory(0);

  TCanvas* c = new TCanvas("c", "", 900, 900);
  TPad* up = new TPad("up", "", 0, 0.32, 1, 1); up->SetBottomMargin(0.02); up->SetLeftMargin(0.13); up->SetLogy(); up->Draw();
  TPad* dn = new TPad("dn", "", 0, 0, 1, 0.32); dn->SetTopMargin(0.02); dn->SetBottomMargin(0.28); dn->SetLeftMargin(0.13); dn->SetGridy(); dn->Draw();

  up->cd();
  hD->GetXaxis()->SetRangeUser(0, xmax);
  hD->SetMarkerStyle(20); hD->SetMarkerSize(0.6); hD->SetMarkerColor(1); hD->SetLineColor(1);
  hD->GetYaxis()->SetTitle("Events"); hD->GetYaxis()->SetTitleSize(0.05); hD->GetXaxis()->SetLabelSize(0);
  hD->SetMinimum(10); hD->SetMaximum(5*hD->GetMaximum());
  hD->Draw("P");
  hM->SetLineColor(632+1); hM->SetLineWidth(2); hM->Draw("hist same");
  up->Update();
  double yLo = TMath::Power(10, up->GetUymin()), yHi = TMath::Power(10, up->GetUymax());
  TLine* lc = new TLine(multCut, yLo, multCut, yHi); lc->SetLineStyle(2); lc->SetLineColor(600); lc->Draw();
  TObjArray* tok = TString(edges).Tokenize(",");
  TLine* lEdge = 0;
  for(int i = 0; i < tok->GetEntries(); i++){
    double e = ((TObjString*) tok->At(i))->GetString().Atof();
    TLine* l = new TLine(e, yLo, e, yHi); l->SetLineStyle(3); l->SetLineColor(920+2); l->SetLineWidth(2); l->Draw();
    lEdge = l;
  }
  TLegend* lg = new TLegend(0.55, 0.72, 0.88, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
  lg->AddEntry(hD, "corrected refMult, minbias-bbc_etof", "p");
  lg->AddEntry(hM, "Glauber #times NBD", "l");
  lg->AddEntry(lc, Form("fit range starts (refMult = %.0f)", multCut), "l");
  if(lEdge) lg->AddEntry(lEdge, Form("class edges (%s)", edges), "l");
  lg->Draw();

  dn->cd();
  TH1D* hR = (TH1D*) hD->Clone("hRatio"); hR->SetDirectory(0);
  hR->Divide(hM);
  hR->GetXaxis()->SetRangeUser(0, xmax);
  hR->SetMinimum(0); hR->SetMaximum(1.5);
  hR->GetYaxis()->SetTitle("data / model"); hR->GetXaxis()->SetTitle("refMult (corrected)");
  hR->GetXaxis()->SetLabelSize(0.09); hR->GetXaxis()->SetTitleSize(0.1); hR->GetXaxis()->SetTitleOffset(1.0);
  hR->GetYaxis()->SetLabelSize(0.08); hR->GetYaxis()->SetTitleSize(0.09); hR->GetYaxis()->SetTitleOffset(0.6); hR->GetYaxis()->SetNdivisions(505);
  hR->Draw("P");
  TLine* l1 = new TLine(0, 1, xmax, 1); l1->SetLineStyle(2); l1->Draw();
  c->SaveAs(outName);
}
