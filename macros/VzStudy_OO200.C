// VzStudy_OO200.C -- lightflavorspectra_OO200
//
// Vz-window study: are the raw yields and the BTOF matching the same in every Vz slice?
// Reads the VzStudy/ directory that PicoBinner writes (since 2026-10-26):
//   vzStudy_Events<Sp>                       TH2I  Vz (12 x 5 cm, -30..30) x centrality index
//   vzStudy_nTPCTracks{Plus,Minus}_ForBTOF_<Sp>  TH3I  y x mT-m0 x Vz  (TPC tracks inside the nSigma cut)
//   vzStudy_nBTOFTracks{Plus,Minus}_<Sp>      TH3I  same tracks with a good BTOF hit
// Each Vz slice is compared with the central one, |Vz| < 5 cm (the two middle slices):
//   page 1: per-event TPC track yield at mid-rapidity vs x, slice / central, h+ and h-
//   page 2: BTOF matching efficiency at mid-rapidity vs x, slice / central
//   page 3: per-event TPC track yield vs y for xLo <= x <= xHi, slice / central (rapidity reach)
// These are raw counts (no efficiency correction), so they test the detector acceptance and
// matching across Vz, which is what a wider Vz window changes. Yields per event use all
// centrality classes in each slice (the class mix was checked to be the same for |Vz| <= 2 and <= 30).
//
// Usage (ROOT 5 or 6, no analysis libraries; the file is the hadd of the PicoBinner outputs):
//   root -l -b -q 'VzStudy_OO200.C+("VzStudy_Pion.root","Pion",0.5,"pT",0.2,1.0)'
//   args: file, species (Pion/Kaon/Proton), |y| window, x axis ("pT" or "mTm0"),
//         x range for page 3 and for the summary numbers (in the x variable chosen)
// Output: VzStudy_<Sp>.pdf (3 pages), VzStudy_<Sp>_pN.png, VzStudy_<Sp>.txt

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH3.h"
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

namespace VzStudy {

const int kNVz = 12;
bool   gUsePt = true;
double gMass  = 0;

double Mass(const TString& sp){
  if(sp == "Pion")   return 0.13957;
  if(sp == "Kaon")   return 0.493677;
  if(sp == "Proton") return 0.938272;
  return 0;
}
double XV(double mtm0){
  if(!gUsePt) return mtm0;
  double mt = mtm0 + gMass;
  return TMath::Sqrt(TMath::Max(mt*mt - gMass*gMass, 0.0));
}
const char* XTitle(){ return gUsePt ? "p_{T} (GeV/c)" : "m_{T}-m_{0} (GeV/c^{2})"; }

// colour per |Vz| slice (0-5 ... 25-30 cm); open markers for Vz < 0
int SliceColor(int iVz){
  const int col[6] = {1, 860+2, 416+2, 800+7, 632+1, 880+1};
  int a = iVz < kNVz/2 ? (kNVz/2 - 1 - iVz) : (iVz - kNVz/2);
  return col[a];
}
int SliceMarker(int iVz){ return iVz < kNVz/2 ? 24 : 20; }

TObject* FindInDir(TDirectory* d, const char* prefix){
  TIter next(d->GetListOfKeys());
  TKey* k;
  while((k = (TKey*)next())){
    if(TString(k->GetName()).BeginsWith(prefix)) return k->ReadObj();
  }
  return 0;
}

bool InWindow(const TAxis* ax, int b, double yMax){
  return ax->GetBinLowEdge(b) >= -yMax - 1e-6 && ax->GetBinUpEdge(b) <= yMax + 1e-6;
}

// counts summed over |y| <= yMax, one value per mT-m0 bin, for Vz bin k
void SumMidRap(const TH3* h, int k, double yMax, double* out, int nx){
  for(int j = 1; j <= nx; j++){
    double s = 0;
    for(int i = 1; i <= h->GetXaxis()->GetNbins(); i++){
      if(!InWindow(h->GetXaxis(), i, yMax)) continue;
      s += h->GetBinContent(i, j, k);
    }
    out[j-1] = s;
  }
}

// counts summed over xLo <= x <= xHi (x in the chosen variable), one value per y bin, for Vz bin k
void SumXRange(const TH3* h, int k, double xLo, double xHi, double* out, int ny){
  for(int i = 1; i <= ny; i++){
    double s = 0;
    for(int j = 1; j <= h->GetYaxis()->GetNbins(); j++){
      double x = XV(h->GetYaxis()->GetBinCenter(j));
      if(x < xLo || x > xHi) continue;
      s += h->GetBinContent(i, j, k);
    }
    out[i-1] = s;
  }
}

} // namespace VzStudy

void VzStudy_OO200(const char* fileName = "VzStudy_Pion.root",
                   const char* species = "Pion",
                   double yMax = 0.5,
                   const char* xVar = "pT",
                   double xLo = 0.2, double xHi = 1.0)
{
  using namespace VzStudy;
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  const TString sp(species);
  gUsePt = (TString(xVar) == "pT");
  gMass  = Mass(sp);
  const char* pm[2] = {"Plus","Minus"};
  const char* sym[2] = {"+","-"};
  TString symName = sp == "Pion" ? "#pi" : (sp == "Kaon" ? "K" : "p");

  TFile* f = TFile::Open(fileName);
  if(f == 0 || f->IsZombie()){ std::cout << "cannot open " << fileName << std::endl; return; }
  TDirectory* d = (TDirectory*) f->Get("VzStudy");
  if(d == 0){ std::cout << "no VzStudy directory in " << fileName << " (PicoBinner older than 2026-10-26?)" << std::endl; return; }
  TH2* hEv = (TH2*) FindInDir(d, "vzStudy_Events");
  TH3* hT[2]; TH3* hB[2];
  for(int q = 0; q < 2; q++){
    hT[q] = (TH3*) d->Get(Form("vzStudy_nTPCTracks%s_ForBTOF_%s", pm[q], species));
    hB[q] = (TH3*) d->Get(Form("vzStudy_nBTOFTracks%s_%s", pm[q], species));
    if(hT[q] == 0 || hB[q] == 0){ std::cout << "missing vzStudy histograms for " << species << pm[q] << std::endl; return; }
  }
  if(hEv == 0 || hEv->GetXaxis()->GetNbins() != kNVz || hT[0]->GetZaxis()->GetNbins() != kNVz){
    std::cout << "unexpected Vz binning" << std::endl; return;
  }

  // events per Vz slice (all classes)
  double nEv[kNVz];
  double nEvTot = 0;
  for(int k = 0; k < kNVz; k++){
    nEv[k] = 0;
    for(int c = 1; c <= hEv->GetYaxis()->GetNbins(); c++) nEv[k] += hEv->GetBinContent(k+1, c);
    nEvTot += nEv[k];
  }
  const int kC1 = kNVz/2 - 1, kC2 = kNVz/2; // the two |Vz| < 5 cm slices
  double nEvC = nEv[kC1] + nEv[kC2];
  if(nEvC <= 0){ std::cout << "no events at |Vz| < 5 cm" << std::endl; return; }

  TString base = Form("VzStudy_%s", species);
  FILE* txt = fopen(Form("%s.txt", base.Data()), "w");
  fprintf(txt, "# %s, |y| <= %.2f, x = %s, summary range %.2f-%.2f\n", species, yMax, xVar, xLo, xHi);
  fprintf(txt, "# events per Vz slice (all classes), total %.0f\n", nEvTot);
  for(int k = 0; k < kNVz; k++) fprintf(txt, "events Vz %+5.0f..%+5.0f  %.0f\n", hEv->GetXaxis()->GetBinLowEdge(k+1), hEv->GetXaxis()->GetBinUpEdge(k+1), nEv[k]);

  const int nx = hT[0]->GetYaxis()->GetNbins();
  const int ny = hT[0]->GetXaxis()->GetNbins();
  double* a = new double[nx]; double* b = new double[nx]; double* ca = new double[nx]; double* cb = new double[nx];
  double* ya = new double[ny]; double* yc = new double[ny];

  const char* pageTitle[2] = {"TPC tracks per event, slice / |V_{z}|<5 cm", "BTOF matching efficiency, slice / |V_{z}|<5 cm"};
  for(int page = 0; page < 2; page++){
    TCanvas* c = new TCanvas(Form("c%d", page), "", 1400, 650);
    c->Divide(2, 1);
    for(int q = 0; q < 2; q++){
      c->cd(q+1); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetGridy();
      double xMax = XV(hT[q]->GetYaxis()->GetXmax()); if(xMax > 2.5) xMax = 2.5;
      TH1D* fr = new TH1D(Form("fr%d_%d", page, q), Form(";%s;%s (%s^{%s}, |y| #leq %.2f)", XTitle(), page == 0 ? "per-event yield / |V_{z}|<5 cm" : "matching eff. / |V_{z}|<5 cm", symName.Data(), sym[q], yMax), 10, 0, xMax);
      fr->SetDirectory(0); fr->SetMinimum(page == 0 ? 0.6 : 0.8); fr->SetMaximum(page == 0 ? 1.4 : 1.2); fr->Draw();
      TLine* l1 = new TLine(0, 1, xMax, 1); l1->SetLineStyle(2); l1->Draw();
      TLegend* lg = new TLegend(0.14, 0.14, 0.62, 0.34); lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.03);
      // central reference
      SumMidRap(hT[q], kC1+1, yMax, a, nx); SumMidRap(hT[q], kC2+1, yMax, b, nx);
      for(int j = 0; j < nx; j++) ca[j] = a[j] + b[j];
      SumMidRap(hB[q], kC1+1, yMax, a, nx); SumMidRap(hB[q], kC2+1, yMax, b, nx);
      for(int j = 0; j < nx; j++) cb[j] = a[j] + b[j];
      for(int k = 0; k < kNVz; k++){
        if(k == kC1 || k == kC2 || nEv[k] <= 0) continue;
        SumMidRap(hT[q], k+1, yMax, a, nx);
        SumMidRap(hB[q], k+1, yMax, b, nx);
        TGraphErrors* g = new TGraphErrors();
        double sw = 0, swx = 0;
        for(int j = 0; j < nx; j++){
          if(a[j] < 20 || ca[j] < 20) continue;
          double r, e;
          if(page == 0){
            r = (a[j]/nEv[k]) / (ca[j]/nEvC);
            e = r*TMath::Sqrt(1/a[j] + 1/ca[j]);
          }else{
            if(b[j] <= 0 || cb[j] <= 0) continue;
            double e1 = b[j]/a[j], e2 = cb[j]/ca[j];
            r = e1/e2;
            e = r*TMath::Sqrt((1 - e1)/b[j] + (1 - e2)/cb[j]);
          }
          double lo = XV(hT[q]->GetYaxis()->GetBinLowEdge(j+1)), hi = XV(hT[q]->GetYaxis()->GetBinUpEdge(j+1));
          int n = g->GetN();
          g->SetPoint(n, 0.5*(lo + hi), r); g->SetPointError(n, 0.5*(hi - lo), e);
          double xc = 0.5*(lo + hi);
          if(xc >= xLo && xc <= xHi && e > 0){ sw += 1/(e*e); swx += r/(e*e); }
        }
        g->SetMarkerStyle(SliceMarker(k)); g->SetMarkerColor(SliceColor(k)); g->SetLineColor(SliceColor(k)); g->SetMarkerSize(0.7);
        g->Draw("P same");
        lg->AddEntry(g, Form("%+.0f..%+.0f", hEv->GetXaxis()->GetBinLowEdge(k+1), hEv->GetXaxis()->GetBinUpEdge(k+1)), "p");
        fprintf(txt, "%s %s%s Vz %+5.0f..%+5.0f  <slice/central> (x %.2f-%.2f) = %.4f +- %.4f\n",
                page == 0 ? "yield" : "matchEff", species, pm[q], hEv->GetXaxis()->GetBinLowEdge(k+1), hEv->GetXaxis()->GetBinUpEdge(k+1),
                xLo, xHi, sw > 0 ? swx/sw : -1, sw > 0 ? 1/TMath::Sqrt(sw) : -1);
      }
      lg->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.DrawLatex(0.13, 0.92, Form("O+O 200 GeV %s^{%s}: %s", symName.Data(), sym[q], pageTitle[page]));
    }
    c->Print(Form("%s.pdf%s", base.Data(), page == 0 ? "(" : ""));
    c->SaveAs(Form("%s_p%d.png", base.Data(), page+1));
  }

  // page 3: rapidity reach
  TCanvas* c3 = new TCanvas("c3", "", 1400, 650);
  c3->Divide(2, 1);
  for(int q = 0; q < 2; q++){
    c3->cd(q+1); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetGridy();
    double yLo = hT[q]->GetXaxis()->GetXmin(), yHi = hT[q]->GetXaxis()->GetXmax();
    TH1D* fr = new TH1D(Form("fr3_%d", q), Form(";y;per-event yield / |V_{z}|<5 cm (%s^{%s}, %.2f #leq %s #leq %.2f)", symName.Data(), sym[q], xLo, gUsePt ? "p_{T}" : "m_{T}-m_{0}", xHi), 10, yLo, yHi);
    fr->SetDirectory(0); fr->SetMinimum(0); fr->SetMaximum(2.0); fr->Draw();
    TLine* l1 = new TLine(yLo, 1, yHi, 1); l1->SetLineStyle(2); l1->Draw();
    SumXRange(hT[q], kC1+1, xLo, xHi, ya, ny);
    for(int i = 0; i < ny; i++) yc[i] = ya[i];
    SumXRange(hT[q], kC2+1, xLo, xHi, ya, ny);
    for(int i = 0; i < ny; i++) yc[i] += ya[i];
    TLegend* lg = new TLegend(0.3, 0.72, 0.88, 0.88); lg->SetNColumns(3); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.03);
    for(int k = 0; k < kNVz; k++){
      if(k == kC1 || k == kC2 || nEv[k] <= 0) continue;
      SumXRange(hT[q], k+1, xLo, xHi, ya, ny);
      TGraphErrors* g = new TGraphErrors();
      for(int i = 0; i < ny; i++){
        if(ya[i] < 20 || yc[i] < 20) continue;
        double r = (ya[i]/nEv[k]) / (yc[i]/nEvC);
        int n = g->GetN();
        g->SetPoint(n, hT[q]->GetXaxis()->GetBinCenter(i+1), r);
        g->SetPointError(n, 0.5*hT[q]->GetXaxis()->GetBinWidth(i+1), r*TMath::Sqrt(1/ya[i] + 1/yc[i]));
      }
      g->SetMarkerStyle(SliceMarker(k)); g->SetMarkerColor(SliceColor(k)); g->SetLineColor(SliceColor(k)); g->SetMarkerSize(0.7);
      g->Draw("P same");
      lg->AddEntry(g, Form("%+.0f..%+.0f", hEv->GetXaxis()->GetBinLowEdge(k+1), hEv->GetXaxis()->GetBinUpEdge(k+1)), "p");
    }
    lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.DrawLatex(0.13, 0.92, Form("O+O 200 GeV %s^{%s}: rapidity reach by V_{z} slice", symName.Data(), sym[q]));
  }
  c3->Print(Form("%s.pdf)", base.Data()));
  c3->SaveAs(Form("%s_p3.png", base.Data()));
  fclose(txt);
  std::cout << "wrote " << base << ".pdf, _p1-3.png, .txt" << std::endl;
}
