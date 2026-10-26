// PlotTpcBtofMatch_OO200.C -- lightflavorspectra_OO200
//
// TPC-BTOF consistency at mid-rapidity, one species at a time:
//   page 1: BTOF matching efficiency vs mT-m0, h+ and h-, all centralities
//           eff = nBTOFTracks / nTPCTracks_ForBTOF (PicoBinner: tracks with |nSigma_TPC| < cut,
//           numerator also has a good BTOF hit). Binomial errors.
//   page 2: h+  TPC vs BTOF/eff yield per centrality (top) and their ratio (bottom)
//   page 3: h-  same
//   page 4: raw h-/h+ ratio from TPC and from BTOF/eff, per centrality
// The ratio panel uses only rapidity bins where TPC, BTOF and the efficiency all exist, with the
// efficiency applied per rapidity bin (CommonRatio). The ratio BTOF/(eff x TPC) should be 1 where both detectors measure the same bins; a step
// there is the TPC-BTOF discontinuity.
//
// Inputs (counts, before any correction):
//   spectraFile: ZFitter output with DeDx_FitData/<Sp>/dEdxSpectra_<Sp>{Plus,Minus}_CentNN and
//                BTOF_FitData/<Sp>/InvBetaBTOFSpectra_<Sp>{Plus,Minus}_CentNN (TH2D: x = y, y = mT-m0,
//                content = fitted counts), e.g. spectra_OO200_Pion_Merged.root
//   yieldFile:   PicoBinner output with nTPCTracks{Plus,Minus}_ForBTOF_<Sp>_CentN and
//                nBTOFTracks{Plus,Minus}_<Sp>_CentN (TH2I: x = y, y = mT-m0), any subdirectory
// Mid-rapidity: rapidity bins lying inside |y| <= yMax. Yields are summed over the filled ones and
// divided by their total dy and by d(mT-m0), so the TPC and BTOF binnings can differ. Centrality
// classes are labelled from the six-class scheme (0-5 ... 80-100%).
//
// Usage (ROOT 5 or 6, no analysis libraries needed):
//   root -l -b -q 'PlotTpcBtofMatch_OO200.C+("Pion","spectra_OO200_Pion_Merged.root","yieldHistos_OO200_pion.root",0.25)'
// Output: <outPrefix>_<Sp>.pdf (4 pages), <outPrefix>_<Sp>_pN.png, <outPrefix>_<Sp>.txt (numbers)

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include <cstdio>
#include <iostream>

namespace TpcBtof {

const int kMaxCent = 6;
const char* kCentLabel[kMaxCent] = {"0-5%","5-10%","10-20%","20-40%","40-80%","80-100%"};
const int kCentColor[kMaxCent]   = {632+1, 800+7, 416+2, 860+2, 600+2, 880+1}; // red, orange, green, azure, blue, violet

double Mass(const TString& sp){
  if(sp == "Pion")   return 0.13957;
  if(sp == "Kaon")   return 0.493677;
  if(sp == "Proton") return 0.938272;
  return 0;
}

TObject* FindRecursive(TDirectory* d, const char* name){
  if(d == 0) return 0;
  TObject* o = d->Get(name);
  if(o) return o;
  TIter next(d->GetListOfKeys());
  TKey* k;
  while((k = (TKey*)next())){
    if(TString(k->GetClassName()).BeginsWith("TDirectory")){
      TObject* r = FindRecursive((TDirectory*)k->ReadObj(), name);
      if(r) return r;
    }
  }
  return 0;
}

// empty TH1D with the binning of an axis (fixed or variable bins)
TH1D* FromAxis(const TAxis* a, const char* name){
  TH1D* h = 0;
  if(a->GetXbins()->GetSize() > 0) h = new TH1D(name, "", a->GetNbins(), a->GetXbins()->GetArray());
  else                             h = new TH1D(name, "", a->GetNbins(), a->GetXmin(), a->GetXmax());
  h->SetDirectory(0);
  h->Sumw2();
  return h;
}

// rapidity bins fully inside |y| <= yMax
bool InWindow(const TAxis* ax, int b, double yMax){
  return ax->GetBinLowEdge(b) >= -yMax - 1e-6 && ax->GetBinUpEdge(b) <= yMax + 1e-6;
}

// fitted counts -> dN/(dy d(mT-m0)) at mid-rapidity; empty (unfitted) bins are skipped and dy
// counts only the filled rapidity bins of each mT-m0 bin
TH1D* MidRapDensity(const TH2* h, double yMax, const char* name){
  const TAxis* ax = h->GetXaxis();
  const TAxis* ay = h->GetYaxis();
  TH1D* out = FromAxis(ay, name);
  for(int j = 1; j <= ay->GetNbins(); j++){
    double sum = 0, err2 = 0, dy = 0;
    for(int i = 1; i <= ax->GetNbins(); i++){
      if(!InWindow(ax, i, yMax)) continue;
      double c = h->GetBinContent(i, j);
      if(c <= 0) continue;
      sum += c; err2 += h->GetBinError(i, j)*h->GetBinError(i, j); dy += ax->GetBinWidth(i);
    }
    if(sum <= 0 || dy <= 0) continue;
    double w = ay->GetBinWidth(j);
    out->SetBinContent(j, sum/(dy*w));
    out->SetBinError(j, TMath::Sqrt(err2)/(dy*w));
  }
  return out;
}

// matched / all at mid-rapidity, binomial errors
TH1D* MatchEff(const TH2* nTpc, const TH2* nTof, double yMax, const char* name){
  const TAxis* ax = nTpc->GetXaxis();
  const TAxis* ay = nTpc->GetYaxis();
  TH1D* out = FromAxis(ay, name);
  for(int j = 1; j <= ay->GetNbins(); j++){
    double n = 0, k = 0;
    for(int i = 1; i <= ax->GetNbins(); i++){
      if(!InWindow(ax, i, yMax)) continue;
      n += nTpc->GetBinContent(i, j); k += nTof->GetBinContent(i, j);
    }
    if(n < 10) continue;
    double e = k/n;
    out->SetBinContent(j, e);
    out->SetBinError(j, TMath::Sqrt(TMath::Max(e*(1 - e), 1.0/n)/n));
  }
  return out;
}

TGraphErrors* ToGraph(const TH1D* h, int color, int marker, double xMax = 1e9){
  TGraphErrors* g = new TGraphErrors();
  for(int j = 1; j <= h->GetNbinsX(); j++){
    if(h->GetBinContent(j) <= 0 || h->GetBinCenter(j) > xMax) continue;
    int n = g->GetN();
    g->SetPoint(n, h->GetBinCenter(j), h->GetBinContent(j));
    g->SetPointError(n, 0.5*h->GetBinWidth(j), h->GetBinError(j));
  }
  g->SetMarkerColor(color); g->SetLineColor(color); g->SetMarkerStyle(marker); g->SetMarkerSize(0.8);
  return g;
}

// numerator density at the denominator's bin centres (numerator bin containing each centre)
TGraphErrors* RatioAtCentres(const TH1D* num, const TH1D* den, int color, int marker){
  TGraphErrors* g = new TGraphErrors();
  for(int j = 1; j <= den->GetNbinsX(); j++){
    double d = den->GetBinContent(j);
    if(d <= 0) continue;
    int k = num->FindFixBin(den->GetBinCenter(j));
    if(k < 1 || k > num->GetNbinsX()) continue;
    double v = num->GetBinContent(k);
    if(v <= 0) continue;
    double r = v/d;
    double rel = TMath::Sqrt(TMath::Power(num->GetBinError(k)/v, 2) + TMath::Power(den->GetBinError(j)/d, 2));
    int n = g->GetN();
    g->SetPoint(n, den->GetBinCenter(j), r);
    g->SetPointError(n, 0.5*den->GetBinWidth(j), r*rel);
  }
  g->SetMarkerColor(color); g->SetLineColor(color); g->SetMarkerStyle(marker); g->SetMarkerSize(0.8);
  return g;
}

// BTOF yield divided by matching efficiency (same binning, bin by bin)
TH1D* Corrected(const TH1D* tof, const TH1D* eff, const char* name){
  TH1D* out = (TH1D*) tof->Clone(name);
  out->SetDirectory(0);
  out->Reset();
  for(int j = 1; j <= tof->GetNbinsX(); j++){
    double v = tof->GetBinContent(j);
    int k = eff->FindFixBin(tof->GetBinCenter(j));
    double e = (k >= 1 && k <= eff->GetNbinsX()) ? eff->GetBinContent(k) : 0;
    if(v <= 0 || e <= 0) continue;
    double rel = TMath::Sqrt(TMath::Power(tof->GetBinError(j)/v, 2) + TMath::Power(eff->GetBinError(k)/e, 2));
    out->SetBinContent(j, v/e);
    out->SetBinError(j, v/e*rel);
  }
  return out;
}

// BTOF/(eff x TPC) with matched rapidity coverage: for each TPC mT-m0 bin, only the rapidity bins
// where the TPC yield, the BTOF yield and the matching efficiency all exist are used, and the
// efficiency is applied per rapidity bin before summing. Without this, the two detectors can be
// averaged over different rapidity bins where the acceptance changes fast (BTOF threshold).
TGraphErrors* CommonRatio(const TH2* hT, const TH2* hB, const TH2* nT, const TH2* nB, double yMax, int color, int marker){
  TGraphErrors* g = new TGraphErrors();
  const TAxis* axT = hT->GetXaxis(); const TAxis* ayT = hT->GetYaxis();
  const TAxis* axB = hB->GetXaxis(); const TAxis* ayB = hB->GetYaxis();
  const TAxis* axN = nT->GetXaxis(); const TAxis* ayN = nT->GetYaxis();
  for(int j = 1; j <= ayT->GetNbins(); j++){
    double x = ayT->GetBinCenter(j);
    int k  = ayB->FindFixBin(x);
    int kn = ayN->FindFixBin(x);
    if(k < 1 || k > ayB->GetNbins() || kn < 1 || kn > ayN->GetNbins()) continue;
    double sT = 0, eT2 = 0, sB = 0, eB2 = 0;
    int nUsed = 0;
    for(int i = 1; i <= axT->GetNbins(); i++){
      if(!InWindow(axT, i, yMax)) continue;
      double yc = axT->GetBinCenter(i);
      int ib = axB->FindFixBin(yc), in = axN->FindFixBin(yc);
      if(ib < 1 || ib > axB->GetNbins() || in < 1 || in > axN->GetNbins()) continue;
      double t = hT->GetBinContent(i, j), b = hB->GetBinContent(ib, k);
      double nt = nT->GetBinContent(in, kn), nb = nB->GetBinContent(in, kn);
      if(t <= 0 || b <= 0 || nt < 10 || nb <= 0) continue;
      double e = nb/nt, ee2 = TMath::Max(e*(1 - e), 1.0/nt)/nt;
      double tD = t/ayT->GetBinWidth(j);
      double bD = b/ayB->GetBinWidth(k)/e;
      sT += tD; eT2 += TMath::Power(hT->GetBinError(i, j)/ayT->GetBinWidth(j), 2);
      sB += bD; eB2 += bD*bD*(TMath::Power(hB->GetBinError(ib, k)/b, 2) + ee2/(e*e));
      nUsed++;
    }
    if(nUsed == 0 || sT <= 0) continue;
    double r = sB/sT;
    double rel = TMath::Sqrt(eB2/(sB*sB) + eT2/(sT*sT));
    int n = g->GetN();
    g->SetPoint(n, x, r);
    g->SetPointError(n, 0.5*ayT->GetBinWidth(j), r*rel);
  }
  g->SetMarkerColor(color); g->SetLineColor(color); g->SetMarkerStyle(marker); g->SetMarkerSize(0.8);
  return g;
}

void WeightedMean(const TGraphErrors* g, double& mean, double& err, double& chi2ndf){
  double sw = 0, swx = 0;
  for(int i = 0; i < g->GetN(); i++){
    double e = g->GetEY()[i]; if(e <= 0) continue;
    sw += 1/(e*e); swx += g->GetY()[i]/(e*e);
  }
  mean = sw > 0 ? swx/sw : 0; err = sw > 0 ? 1/TMath::Sqrt(sw) : 0;
  double c2 = 0; int n = 0;
  for(int i = 0; i < g->GetN(); i++){
    double e = g->GetEY()[i]; if(e <= 0) continue;
    c2 += TMath::Power((g->GetY()[i] - 1)/e, 2); n++;
  }
  chi2ndf = n > 0 ? c2/n : 0;
}

} // namespace TpcBtof

void PlotTpcBtofMatch_OO200(const char* species = "Pion",
                            const char* spectraFile = "spectra_OO200_Pion_Merged.root",
                            const char* yieldFile = "yieldHistos_OO200_pion.root",
                            double yMax = 0.25,
                            const char* outPrefix = "TpcBtof_OO200")
{
  using namespace TpcBtof;
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  const TString sp(species);
  const char* pm[2] = {"Plus","Minus"};
  const char* sym[2] = {"+","-"};
  TString symName = sp == "Pion" ? "#pi" : (sp == "Kaon" ? "K" : "p");

  TFile* fs = TFile::Open(spectraFile);
  TFile* fy = TFile::Open(yieldFile);
  if(fs == 0 || fs->IsZombie()){ std::cout << "cannot open " << spectraFile << std::endl; return; }
  if(fy == 0 || fy->IsZombie()){ std::cout << "cannot open " << yieldFile << std::endl; return; }

  TH1D* tpc[2][kMaxCent]; TH1D* tof[2][kMaxCent]; TH1D* eff[2][kMaxCent]; TH1D* tofC[2][kMaxCent];
  TH2* h2T[2][kMaxCent]; TH2* h2B[2][kMaxCent]; TH2* n2T[2][kMaxCent]; TH2* n2B[2][kMaxCent];
  int nCent = 0;
  for(int c = 0; c < kMaxCent; c++){
    bool ok = true;
    for(int q = 0; q < 2; q++){
      tpc[q][c] = tof[q][c] = eff[q][c] = tofC[q][c] = 0;
      TH2* hT = (TH2*) fs->Get(Form("DeDx_FitData/%s/dEdxSpectra_%s%s_Cent%02d", species, species, pm[q], c));
      TH2* hB = (TH2*) fs->Get(Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%s%s_Cent%02d", species, species, pm[q], c));
      TH2* nT = (TH2*) FindRecursive(fy, Form("nTPCTracks%s_ForBTOF_%s_Cent%d", pm[q], species, c));
      TH2* nB = (TH2*) FindRecursive(fy, Form("nBTOFTracks%s_%s_Cent%d", pm[q], species, c));
      if(hT == 0 || hB == 0 || nT == 0 || nB == 0){
        if(c == 0) std::cout << "missing for cent 0, charge " << pm[q] << ": " << (hT ? "" : "dEdxSpectra ") << (hB ? "" : "InvBetaBTOFSpectra ")
                             << (nT ? "" : "nTPCTracks_ForBTOF ") << (nB ? "" : "nBTOFTracks") << std::endl;
        ok = false; break;
      }
      h2T[q][c] = hT; h2B[q][c] = hB; n2T[q][c] = nT; n2B[q][c] = nB;
      tpc[q][c]  = MidRapDensity(hT, yMax, Form("tpc_%s_%d", pm[q], c));
      tof[q][c]  = MidRapDensity(hB, yMax, Form("tof_%s_%d", pm[q], c));
      eff[q][c]  = MatchEff(nT, nB, yMax, Form("eff_%s_%d", pm[q], c));
      tofC[q][c] = Corrected(tof[q][c], eff[q][c], Form("tofC_%s_%d", pm[q], c));
      if(c == 0 && q == 0){
        TString w = "";
        for(int i = 1; i <= hT->GetNbinsX(); i++) if(InWindow(hT->GetXaxis(), i, yMax)) w += Form(" [%.2f,%.2f]", hT->GetXaxis()->GetBinLowEdge(i), hT->GetXaxis()->GetBinUpEdge(i));
        std::cout << "TPC rapidity bins used:" << w << std::endl;
        w = "";
        for(int i = 1; i <= hB->GetNbinsX(); i++) if(InWindow(hB->GetXaxis(), i, yMax)) w += Form(" [%.2f,%.2f]", hB->GetXaxis()->GetBinLowEdge(i), hB->GetXaxis()->GetBinUpEdge(i));
        std::cout << "BTOF rapidity bins used:" << w << std::endl;
      }
    }
    if(!ok) break;
    nCent = c + 1;
  }
  if(nCent == 0){ std::cout << "nothing found; check species name and files" << std::endl; return; }
  std::cout << species << ": " << nCent << " centrality classes" << std::endl;

  TString base = Form("%s_%s", outPrefix, species);
  FILE* txt = fopen(Form("%s.txt", base.Data()), "w");
  fprintf(txt, "# %s, |y| <= %.2f\n", species, yMax);

  // ---------- page 1: matching efficiency
  TCanvas* c1 = new TCanvas("c1", "", 1400, 650);
  c1->Divide(2, 1);
  for(int q = 0; q < 2; q++){
    c1->cd(q+1); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetGridy();
    TH1D* fr = new TH1D(Form("fr1_%d", q), Form(";m_{T}-m_{0} (GeV/c^{2});BTOF matching efficiency (%s^{%s}, |y| #leq %.2f)", symName.Data(), sym[q], yMax), 10, 0, 3);
    fr->SetDirectory(0); fr->SetMinimum(0); fr->SetMaximum(1.0); fr->Draw();
    TLegend* lg = new TLegend(0.62, 0.15, 0.88, 0.45); lg->SetBorderSize(0); lg->SetFillStyle(0);
    for(int c = 0; c < nCent; c++){
      TGraphErrors* g = ToGraph(eff[q][c], kCentColor[c], 20);
      g->Draw("P same"); lg->AddEntry(g, kCentLabel[c], "p");
      // plateau: weighted mean over mT-m0 0.5-1.5
      double sw = 0, swx = 0;
      for(int j = 1; j <= eff[q][c]->GetNbinsX(); j++){
        double x = eff[q][c]->GetBinCenter(j), e = eff[q][c]->GetBinError(j);
        if(x < 0.5 || x > 1.5 || e <= 0) continue;
        sw += 1/(e*e); swx += eff[q][c]->GetBinContent(j)/(e*e);
      }
      fprintf(txt, "matchEff %s%s %-8s plateau(mT-m0 0.5-1.5) = %.4f\n", species, pm[q], kCentLabel[c], sw > 0 ? swx/sw : -1);
    }
    lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.045); t.DrawLatex(0.15, 0.92, Form("O+O 200 GeV, %s^{%s}: nBTOF / nTPC", symName.Data(), sym[q]));
  }
  c1->Print(Form("%s.pdf(", base.Data()));
  c1->SaveAs(Form("%s_p1.png", base.Data()));

  // ---------- pages 2-3: TPC vs BTOF/eff per centrality
  for(int q = 0; q < 2; q++){
    TCanvas* c2 = new TCanvas(Form("c2_%d", q), "", 1500, 950);
    c2->Divide(3, 2, 0.002, 0.002);
    for(int c = 0; c < nCent; c++){
      c2->cd(c+1);
      TPad* up = new TPad(Form("up_%d_%d", q, c), "", 0, 0.38, 1, 1);
      TPad* dn = new TPad(Form("dn_%d_%d", q, c), "", 0, 0, 1, 0.38);
      up->SetBottomMargin(0.02); up->SetLeftMargin(0.15); up->SetLogy();
      dn->SetTopMargin(0.02); dn->SetBottomMargin(0.28); dn->SetLeftMargin(0.15); dn->SetGridy();
      up->Draw(); dn->Draw();

      // x range: up to where both detectors have points, plus some BTOF reach
      double xMax = 0;
      for(int j = tof[q][c]->GetNbinsX(); j >= 1; j--) if(tof[q][c]->GetBinContent(j) > 0){ xMax = tof[q][c]->GetBinLowEdge(j+1); break; }
      xMax = TMath::Min(TMath::Max(xMax, 0.5), 2.5);

      up->cd();
      TGraphErrors* gT = ToGraph(tpc[q][c], 1, 20);
      TGraphErrors* gB = ToGraph(tof[q][c], 920+2, 24);
      TGraphErrors* gC = ToGraph(tofC[q][c], 632+1, 25);
      double yMin = 1e30, yMx = 0;
      TGraphErrors* gs[3] = {gT, gB, gC};
      for(int k = 0; k < 3; k++) for(int i = 0; i < gs[k]->GetN(); i++){ double v = gs[k]->GetY()[i]; if(v > 0){ yMin = TMath::Min(yMin, v); yMx = TMath::Max(yMx, v); } }
      TH1D* fr = new TH1D(Form("fr2_%d_%d", q, c), ";;dN/dy d(m_{T}-m_{0}) (counts)", 10, 0, xMax);
      fr->SetDirectory(0); fr->SetMinimum(yMin*0.3); fr->SetMaximum(yMx*5);
      fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetTitleSize(0.06); fr->GetYaxis()->SetLabelSize(0.055); fr->GetYaxis()->SetTitleOffset(1.1);
      fr->Draw();
      gT->Draw("P same"); gB->Draw("P same"); gC->Draw("P same");
      TLatex t; t.SetNDC(); t.SetTextSize(0.07); t.DrawLatex(0.55, 0.82, Form("%s^{%s} %s", symName.Data(), sym[q], kCentLabel[c]));
      if(c == 0){
        TLegend* lg = new TLegend(0.5, 0.45, 0.88, 0.78); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.055);
        lg->AddEntry(gT, "TPC dE/dx", "p"); lg->AddEntry(gB, "BTOF 1/#beta (raw)", "p"); lg->AddEntry(gC, "BTOF / match. eff.", "p"); lg->Draw();
      }

      dn->cd();
      TGraphErrors* gR = CommonRatio(h2T[q][c], h2B[q][c], n2T[q][c], n2B[q][c], yMax, 632+1, 25);
      TH1D* fr2 = new TH1D(Form("fr3_%d_%d", q, c), ";m_{T}-m_{0} (GeV/c^{2});BTOF_{corr} / TPC", 10, 0, xMax);
      fr2->SetDirectory(0); fr2->SetMinimum(0.5); fr2->SetMaximum(1.5);
      fr2->GetXaxis()->SetTitleSize(0.1); fr2->GetXaxis()->SetLabelSize(0.09); fr2->GetYaxis()->SetTitleSize(0.085);
      fr2->GetYaxis()->SetLabelSize(0.08); fr2->GetYaxis()->SetTitleOffset(0.75); fr2->GetYaxis()->SetNdivisions(505);
      fr2->Draw();
      TLine* l1 = new TLine(0, 1, xMax, 1); l1->SetLineStyle(2); l1->Draw();
      gR->Draw("P same");

      double mean, err, chi2; WeightedMean(gR, mean, err, chi2);
      double xLo = gR->GetN() ? gR->GetX()[0] : 0, xHi = gR->GetN() ? gR->GetX()[gR->GetN()-1] : 0;
      fprintf(txt, "overlap %s%s %-8s n=%2d mT-m0 %.3f-%.3f  <BTOFcorr/TPC> = %.4f +- %.4f  chi2/n vs 1 = %.2f\n",
              species, pm[q], kCentLabel[c], gR->GetN(), xLo, xHi, mean, err, chi2);
      for(int i = 0; i < gR->GetN(); i++)
        fprintf(txt, "   %s%s %-8s mT-m0 %.3f  ratio %.4f +- %.4f\n", species, pm[q], kCentLabel[c], gR->GetX()[i], gR->GetY()[i], gR->GetEY()[i]);
    }
    c2->Print(Form("%s.pdf", base.Data()));
    c2->SaveAs(Form("%s_p%d.png", base.Data(), q+2));
  }

  // ---------- page 4: h-/h+ raw ratio
  TCanvas* c4 = new TCanvas("c4", "", 1500, 900);
  c4->Divide(3, 2, 0.002, 0.002);
  for(int c = 0; c < nCent; c++){
    c4->cd(c+1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetGridy();
    TGraphErrors* rT = RatioAtCentres(tpc[1][c], tpc[0][c], 1, 20);
    TGraphErrors* rB = RatioAtCentres(tofC[1][c], tofC[0][c], 632+1, 25);
    TH1D* fr = new TH1D(Form("fr4_%d", c), Form(";m_{T}-m_{0} (GeV/c^{2});%s^{-}/%s^{+} (raw)", symName.Data(), symName.Data()), 10, 0, 2.5);
    fr->SetDirectory(0);
    double lo = sp == "Proton" ? 0.4 : 0.7, hi = sp == "Proton" ? 1.2 : 1.3;
    fr->SetMinimum(lo); fr->SetMaximum(hi);
    fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->Draw();
    rT->Draw("P same"); rB->Draw("P same");
    TLatex t; t.SetNDC(); t.SetTextSize(0.055); t.DrawLatex(0.18, 0.85, kCentLabel[c]);
    if(c == 0){
      TLegend* lg = new TLegend(0.5, 0.72, 0.88, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
      lg->AddEntry(rT, "TPC", "p"); lg->AddEntry(rB, "BTOF / match. eff.", "p"); lg->Draw();
    }
    double mT_, eT_, c2T_, mB_, eB_, c2B_;
    WeightedMean(rT, mT_, eT_, c2T_); WeightedMean(rB, mB_, eB_, c2B_);
    fprintf(txt, "antiratio %s %-8s TPC <h-/h+> = %.4f +- %.4f (n=%d)   BTOF <h-/h+> = %.4f +- %.4f (n=%d)\n",
            species, kCentLabel[c], mT_, eT_, rT->GetN(), mB_, eB_, rB->GetN());
  }
  c4->Print(Form("%s.pdf)", base.Data()));
  c4->SaveAs(Form("%s_p4.png", base.Data()));
  fclose(txt);
  std::cout << "wrote " << base << ".pdf, _p1-4.png, .txt" << std::endl;
}
