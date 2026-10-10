// FitRefMultCorr_OO200.C -- lightflavorspectra_OO200
//
// Step 2 of building an StRefMultCorr-style centrality for O+O 200 GeV (Run 21), trigger 860003.
// Reads the merged output of macros/MakeRefMultCorrInputs_OO200.C and derives:
//
//   (a) Vz correction      c_vz(Vz) = <refMult>(center) / <refMult>(Vz), from refMult >= kMultMin,
//                          in 1 cm bins over |Vz| <= 30 cm, plus a pol6 fit for reference.
//                          Cross-check: the 99th percentile of refMult per Vz bin.
//   (b) Luminosity corr.   <refMult>(refMult >= kMultMin) vs ZDCx for |Vz| <= 10 cm, fitted with
//                          p0 + p1*ZDCx; c_lumi(ZDCx) = p0 / (p0 + p1*ZDCx) (scaled to zero rate).
//   (c) Pile-up band       refMult vs nBTOFMatch: per nBTOFMatch slice, truncated mean and RMS of
//                          refMult; upper/lower edges mean +- kNSigma*RMS fitted with pol4 in
//                          nBTOFMatch. Coefficients printed in StRefMultCorr order (b0..b4 upper,
//                          c0..c4 lower: keep refMult inside [c(nBTOFMatch), b(nBTOFMatch)]).
//                          Beyond the last fitted slice the edges continue as straight lines
//                          (tangent there), because an unconstrained pol4 turns over.
//                          If the input has hRefMult_Vz_ZDCx_band (pass 2 of the inputs job, with the
//                          pass-1 band applied), (d) uses it; otherwise it uses hRefMult_Vz_ZDCx.
//   (d) Corrected refMult  refMult + U[0,1) times c_vz times c_lumi, rebuilt from the 3D histogram
//                          (refMult x Vz x ZDCx) by spreading each cell uniformly over its corrected
//                          interval. Written as hRefMult (unit bins) for the NBD rescan.
//
// Outputs:  OO200_RefMultCorr_params.txt           all parameters, human readable
//           hRefMultCorr_OO200.root                hRefMult (corrected, |Vz| <= 30 cm),
//                                                  hRefMultRaw (uncorrected, same events),
//                                                  hRefMultCorr_Vz2 (corrected, |Vz| <= 2 cm)
//           OO200_RMC_{VzCorr,LumiCorr,PileUp,RefMultCorr}.png
//
// Usage (container, SL24y, repo root):
//   root -l -b -q 'macros/FitRefMultCorr_OO200.C+("rmc_OO200_all.root")'
// Written for ROOT 5 (C++98): no auto / nullptr / range-for.

#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"

namespace {
  const double kMultMin   = 15.;   // refMult threshold for the Vz and luminosity means (fit range of the NBD fit)
  const double kVzCenter  = 2.;    // |Vz| <= 2 cm defines the reference
  const double kNSigma    = 5.;    // pile-up band half-width in RMS
  const int    kTofSlice  = 2;     // nBTOFMatch bins per slice
  const double kMinSliceEntries = 2000.;

  // mean of a 1D histogram for x >= xmin
  double meanAbove(TH1* h, double xmin, double& err){
    double s = 0, sx = 0, sxx = 0;
    for(int i = 1; i <= h->GetNbinsX(); i++){
      double x = h->GetBinLowEdge(i);           // refMult is an integer: use the bin low edge
      if(x < xmin) continue;
      double w = h->GetBinContent(i);
      s += w; sx += w*x; sxx += w*x*x;
    }
    if(s <= 0){ err = 0; return 0; }
    double m = sx/s; double var = sxx/s - m*m;
    err = (var > 0) ? std::sqrt(var/s) : 0;
    return m;
  }
  // quantile of a 1D histogram
  double quantile(TH1* h, double q){
    double tot = h->Integral(1, h->GetNbinsX());
    if(tot <= 0) return 0;
    double c = 0;
    for(int i = 1; i <= h->GetNbinsX(); i++){
      c += h->GetBinContent(i);
      if(c >= q*tot) return h->GetBinLowEdge(i);
    }
    return h->GetBinLowEdge(h->GetNbinsX());
  }
  // pile-up band edge: pol4 in nBTOFMatch up to the last fitted slice p[5], continued linearly
  // (tangent at p[5]) beyond it, so the edges cannot turn over where there are no data to fit
  double bandFunc(double* x, double* p){
    const double xm = p[5];
    const double xx = (x[0] <= xm) ? x[0] : xm;
    double v = p[0] + p[1]*xx + p[2]*xx*xx + p[3]*xx*xx*xx + p[4]*xx*xx*xx*xx;
    if(x[0] > xm) v += (p[1] + 2*p[2]*xm + 3*p[3]*xm*xm + 4*p[4]*xm*xm*xm) * (x[0] - xm);
    return v;
  }
  // add weight w spread uniformly over [lo, hi) into unit-width histogram hOut
  void spread(TH1D* hOut, double lo, double hi, double w){
    if(w <= 0 || hi <= lo) return;
    double width = hi - lo;
    int b1 = hOut->FindBin(lo), b2 = hOut->FindBin(hi - 1e-9);
    for(int b = b1; b <= b2; b++){
      double e1 = TMath::Max(lo, hOut->GetBinLowEdge(b));
      double e2 = TMath::Min(hi, hOut->GetBinLowEdge(b) + hOut->GetBinWidth(b));
      if(e2 > e1) hOut->AddBinContent(b, w*(e2 - e1)/width);
    }
  }
}

void FitRefMultCorr_OO200(const char* inName = "rmc_OO200_all.root"){
  gStyle->SetOptStat(0);
  TFile* in = TFile::Open(inName);
  if(!in || in->IsZombie()){ std::cout << "Cannot open " << inName << std::endl; return; }
  TH2F* hRefVz  = (TH2F*) in->Get("hRefMult_vs_Vz");
  TH2F* hRefZdc = (TH2F*) in->Get("hRefMult_vs_ZDCx");
  TH2F* hRefTof = (TH2F*) in->Get("hRefMult_vs_nBTOFMatch");
  TH3F* h3      = (TH3F*) in->Get("hRefMult_Vz_ZDCx");
  // pass 2 of MakeRefMultCorrInputs_OO200.C also stores the 3D histogram after the pass-1 pile-up band
  TH3F* h3band  = (TH3F*) in->Get("hRefMult_Vz_ZDCx_band");
  TH1D* hEv     = (TH1D*) in->Get("hEvents");
  if(!hRefVz || !hRefZdc || !hRefTof || !h3){ std::cout << "Missing input histograms in " << inName << std::endl; return; }

  std::ofstream par("OO200_RefMultCorr_params.txt");
  par << "# StRefMultCorr-style parameters, O+O 200 GeV Run 21, trigger 860003, refMult, |Vz| <= 30 cm\n";
  if(hEv) par << "# events: 860003 " << hEv->GetBinContent(1) << ", |Vz|<=30 " << hEv->GetBinContent(2)
              << ", Vr<=1 " << hEv->GetBinContent(3) << ", current pile-up line " << hEv->GetBinContent(4) << "\n";

  // ------------------------------------------------------------------ (a) Vz correction
  const int nVz = hRefVz->GetNbinsX();
  std::vector<double> vzC(nVz), mean(nVz), meanErr(nVz), q99(nVz);
  double sRef = 0, wRef = 0, sQ = 0, nQ = 0;
  for(int i = 1; i <= nVz; i++){
    TH1D* p = hRefVz->ProjectionY(Form("pVz_%d", i), i, i);
    vzC[i-1] = hRefVz->GetXaxis()->GetBinCenter(i);
    double e; mean[i-1] = meanAbove(p, kMultMin, e); meanErr[i-1] = e;
    q99[i-1] = quantile(p, 0.99);
    if(TMath::Abs(vzC[i-1]) <= kVzCenter && e > 0){ sRef += mean[i-1]/(e*e); wRef += 1./(e*e); sQ += q99[i-1]; nQ++; }
    delete p;
  }
  const double refCenter = (wRef > 0) ? sRef/wRef : 0;
  const double q99Center = (nQ > 0) ? sQ/nQ : 0;
  TGraphErrors* gVz = new TGraphErrors(); gVz->SetName("gVzCorr");
  TGraph* gQ = new TGraph(); gQ->SetName("gVzCorrQ99");
  std::vector<double> cvz(nVz, 1.);
  for(int i = 0; i < nVz; i++){
    if(mean[i] > 0){
      cvz[i] = refCenter/mean[i];
      int n = gVz->GetN(); gVz->SetPoint(n, vzC[i], cvz[i]); gVz->SetPointError(n, 0, cvz[i]*meanErr[i]/mean[i]);
    }
    if(q99[i] > 0){ int n = gQ->GetN(); gQ->SetPoint(n, vzC[i], q99Center/q99[i]); }
  }
  TF1* fVz = new TF1("fVz", "pol6", -30, 30);
  gVz->Fit(fVz, "Q0");
  par << "\n[Vz correction]  c_vz = <refMult>(|Vz|<=" << kVzCenter << ") / <refMult>(Vz), refMult >= " << kMultMin
      << ", reference <refMult> = " << refCenter << "\n";
  par << "pol6:";
  for(int k = 0; k < 7; k++) par << " " << fVz->GetParameter(k);
  par << "\n# Vz bin center (cm), c_vz (mean method), c_vz (99th-percentile method)\n";
  for(int i = 0; i < nVz; i++) par << Form("%6.1f  %.5f  %.5f\n", vzC[i], cvz[i], q99[i] > 0 ? q99Center/q99[i] : 0.);

  TCanvas* c1 = new TCanvas("c1", "", 900, 600);
  gVz->SetMarkerStyle(20); gVz->SetTitle("Vz correction;V_{z} (cm);c_{vz} = #LTrefMult#GT_{center} / #LTrefMult#GT(V_{z})");
  gVz->Draw("AP"); fVz->SetLineColor(kRed); fVz->Draw("same");
  gQ->SetMarkerStyle(24); gQ->SetMarkerColor(kBlue); gQ->Draw("P same");
  TLegend* l1 = new TLegend(0.15, 0.75, 0.6, 0.88);
  l1->AddEntry(gVz, Form("mean, refMult #geq %.0f", kMultMin), "p"); l1->AddEntry(fVz, "pol6 fit", "l");
  l1->AddEntry(gQ, "99th percentile", "p"); l1->Draw();
  c1->Print("OO200_RMC_VzCorr.png");

  // ------------------------------------------------------------------ (b) luminosity correction
  TGraphErrors* gL = new TGraphErrors(); gL->SetName("gLumi");
  for(int i = 1; i <= hRefZdc->GetNbinsX(); i++){
    TH1D* p = hRefZdc->ProjectionY(Form("pZdc_%d", i), i, i);
    double e; double m = meanAbove(p, kMultMin, e);
    double nAbove = p->Integral(p->FindBin(kMultMin), p->GetNbinsX());
    if(nAbove > 1000 && m > 0){ int n = gL->GetN(); gL->SetPoint(n, hRefZdc->GetXaxis()->GetBinCenter(i), m); gL->SetPointError(n, 0, e); }
    delete p;
  }
  TF1* fL = new TF1("fL", "pol1", 0, 3);
  double p0 = refCenter, p1 = 0;
  if(gL->GetN() >= 2){ gL->Fit(fL, "Q0"); p0 = fL->GetParameter(0); p1 = fL->GetParameter(1); }
  par << "\n[Luminosity correction]  <refMult>(refMult >= " << kMultMin << ", |Vz| <= 10 cm) = p0 + p1*ZDCx[kHz]\n";
  par << "p0 = " << p0 << " +- " << fL->GetParError(0) << "\np1 = " << p1 << " +- " << fL->GetParError(1) << " per kHz"
      << "  (= " << p1/1000. << " per Hz)\n";
  par << "c_lumi(ZDCx) = p0 / (p0 + p1*ZDCx);  ZDCx range in the data: see plot (lever arm is small at ~1 kHz)\n";
  TCanvas* c2 = new TCanvas("c2", "", 900, 600);
  gL->SetMarkerStyle(20); gL->SetTitle(Form("Luminosity dependence, |V_{z}| #leq 10 cm;ZDCx (kHz);#LTrefMult#GT (refMult #geq %.0f)", kMultMin));
  gL->Draw("AP"); fL->SetLineColor(kRed); fL->Draw("same");
  c2->Print("OO200_RMC_LumiCorr.png");

  // ------------------------------------------------------------------ (c) pile-up band
  TGraph* gUp = new TGraph(); TGraph* gLo = new TGraph(); TGraph* gMu = new TGraph();
  const int nTof = hRefTof->GetNbinsX();
  for(int i = 1; i <= nTof; i += kTofSlice){
    TH1D* p = hRefTof->ProjectionY(Form("pTof_%d", i), i, i + kTofSlice - 1);
    if(p->GetEntries() < kMinSliceEntries){ delete p; continue; }
    double m = p->GetMean(), s = p->GetRMS();
    for(int it = 0; it < 5; it++){            // truncated mean/RMS within +-3 sigma
      p->GetXaxis()->SetRangeUser(TMath::Max(0., m - 3*s), m + 3*s);
      m = p->GetMean(); s = p->GetRMS();
    }
    double x = 0.5*(hRefTof->GetXaxis()->GetBinLowEdge(i) + hRefTof->GetXaxis()->GetBinUpEdge(i + kTofSlice - 1));
    int n = gMu->GetN();
    gMu->SetPoint(n, x, m); gUp->SetPoint(n, x, m + kNSigma*s); gLo->SetPoint(n, x, m - kNSigma*s);
    delete p;
  }
  TF1* fUp = new TF1("fUp", "pol4", 0, 300); TF1* fLo = new TF1("fLo", "pol4", 0, 300);
  if(gUp->GetN() > 5){ gUp->Fit(fUp, "Q0"); gLo->Fit(fLo, "Q0"); }
  const double xmFit = gUp->GetN() ? gUp->GetX()[gUp->GetN()-1] : 300.;
  TF1* fUpE = new TF1("fUpE", bandFunc, 0, 300, 6); TF1* fLoE = new TF1("fLoE", bandFunc, 0, 300, 6);
  for(int k = 0; k < 5; k++){ fUpE->SetParameter(k, fUp->GetParameter(k)); fLoE->SetParameter(k, fLo->GetParameter(k)); }
  fUpE->SetParameter(5, xmFit); fLoE->SetParameter(5, xmFit);
  par << "\n[Pile-up band]  refMult vs nBTOFMatch, mean +- " << kNSigma << " x truncated RMS per " << kTofSlice
      << "-unit nBTOFMatch slice (slices with >= " << kMinSliceEntries << " events), fitted with pol4\n";
  par << "upper (b0..b4):"; for(int k = 0; k < 5; k++) par << " " << fUp->GetParameter(k); par << "\n";
  par << "lower (c0..c4):"; for(int k = 0; k < 5; k++) par << " " << fLo->GetParameter(k); par << "\n";
  par << "fit range in nBTOFMatch: 0 - " << (gUp->GetN() ? gUp->GetX()[gUp->GetN()-1] : 0) << " (beyond it: straight line, tangent at the last fitted point)\n";
  par << "band edges at the last fitted point: upper " << fUpE->Eval(xmFit) << ", lower " << fLoE->Eval(xmFit) << "\n";
  // fraction of events outside the band
  double inBand = 0, all = 0;
  for(int i = 1; i <= nTof; i++){
    double x = hRefTof->GetXaxis()->GetBinCenter(i);
    for(int j = 1; j <= hRefTof->GetNbinsY(); j++){
      double y = hRefTof->GetYaxis()->GetBinLowEdge(j), w = hRefTof->GetBinContent(i, j);
      all += w; if(y >= fLoE->Eval(x) && y <= fUpE->Eval(x)) inBand += w;
    }
  }
  par << "events inside the band: " << inBand << " of " << all << " (" << (all > 0 ? 100.*inBand/all : 0) << "%)\n";
  TCanvas* c3 = new TCanvas("c3", "", 900, 700); c3->SetLogz();
  hRefTof->SetTitle("refMult vs nBTOFMatch, |V_{z}| #leq 30 cm;nBTOFMatch;refMult");
  hRefTof->Draw("colz");
  fUpE->SetLineColor(kRed); fLoE->SetLineColor(kRed); fUpE->SetNpx(600); fLoE->SetNpx(600); fUpE->Draw("same"); fLoE->Draw("same");
  TF1* fOld = new TF1("fOld", "x - 100", 100, 300); fOld->SetLineColor(kBlack); fOld->SetLineStyle(2); fOld->Draw("same");
  TLegend* l3 = new TLegend(0.15, 0.78, 0.6, 0.88);
  l3->AddEntry(fUpE, Form("new band (#pm%.0f#sigma)", kNSigma), "l"); l3->AddEntry(fOld, "current cut: nBTOFMatch #leq refMult+100", "l"); l3->Draw();
  c3->Print("OO200_RMC_PileUp.png");

  // ------------------------------------------------------------------ (d) corrected refMult
  TH3F* h3use = h3band ? h3band : h3;
  par << "\n[Corrected refMult input] " << (h3band ? "hRefMult_Vz_ZDCx_band (pass-1 pile-up band applied)" : "hRefMult_Vz_ZDCx (no band; only nBTOFMatch <= refMult+100)") << "\n";
  TH1D* hCorr  = new TH1D("hRefMult", "Corrected refMult (Vz, luminosity, [0,1) smearing), |V_{z}| #leq 30 cm;refMult_{corr};Events", 300, 0, 300);
  TH1D* hRaw   = new TH1D("hRefMultRaw", "Uncorrected refMult, same events;refMult;Events", 300, 0, 300);
  TH1D* hCorr2 = new TH1D("hRefMultCorr_Vz2", "Corrected refMult, |V_{z}| #leq 2 cm;refMult_{corr};Events", 300, 0, 300);
  for(int iz = 1; iz <= h3use->GetNbinsZ(); iz++){
    double zdc = h3use->GetZaxis()->GetBinCenter(iz);
    double cl = (p0 + p1*zdc > 0) ? p0/(p0 + p1*zdc) : 1.;
    for(int iy = 1; iy <= h3use->GetNbinsY(); iy++){
      double vz = h3use->GetYaxis()->GetBinCenter(iy);
      int ivz = hRefVz->GetXaxis()->FindBin(vz);
      double c = cl * ((ivz >= 1 && ivz <= nVz) ? cvz[ivz-1] : 1.);
      for(int ix = 1; ix <= h3use->GetNbinsX(); ix++){
        double w = h3use->GetBinContent(ix, iy, iz);
        if(w <= 0) continue;
        double r = h3use->GetXaxis()->GetBinLowEdge(ix);    // integer refMult
        hRaw->Fill(r + 0.5, w);
        spread(hCorr, r*c, (r + 1.)*c, w);
        if(TMath::Abs(vz) <= kVzCenter) spread(hCorr2, r*c, (r + 1.)*c, w);
      }
    }
  }
  TFile* out = TFile::Open("hRefMultCorr_OO200.root", "RECREATE");
  hCorr->Write(); hRaw->Write(); hCorr2->Write(); gVz->Write(); gQ->Write(); gL->Write(); gUp->Write("gPileUpUpper"); gLo->Write("gPileUpLower"); fUpE->Write("fPileUpUpper"); fLoE->Write("fPileUpLower");
  fVz->Write(); fL->Write(); fUp->Write(); fLo->Write();
  out->Close();
  par << "\n[Corrected refMult]  hRefMult in hRefMultCorr_OO200.root: " << hCorr->Integral() << " events (|Vz| <= 30 cm), "
      << hCorr2->Integral() << " with |Vz| <= 2 cm\n";
  par.close();

  TCanvas* c4 = new TCanvas("c4", "", 900, 600); c4->SetLogy();
  TH1D* r1 = (TH1D*) hRaw->Clone("r1");  r1->Scale(1./r1->Integral());
  TH1D* c1n = (TH1D*) hCorr->Clone("c1n"); c1n->Scale(1./c1n->Integral());
  TH1D* c2n = (TH1D*) hCorr2->Clone("c2n"); if(c2n->Integral() > 0) c2n->Scale(1./c2n->Integral());
  r1->SetLineColor(kBlack); c1n->SetLineColor(kRed); c2n->SetLineColor(kBlue);
  r1->SetTitle("refMult before and after corrections;refMult;normalized");
  r1->GetXaxis()->SetRangeUser(0, 150); r1->Draw("hist"); c1n->Draw("hist same"); c2n->Draw("hist same");
  TLegend* l4 = new TLegend(0.55, 0.7, 0.88, 0.88);
  l4->AddEntry(r1, "raw, |V_{z}| #leq 30 cm", "l"); l4->AddEntry(c1n, "corrected, |V_{z}| #leq 30 cm", "l"); l4->AddEntry(c2n, "corrected, |V_{z}| #leq 2 cm", "l"); l4->Draw();
  c4->Print("OO200_RMC_RefMultCorr.png");

  std::cout << "Reference <refMult> (refMult >= " << kMultMin << ", |Vz| <= " << kVzCenter << "): " << refCenter << std::endl;
  std::cout << "Luminosity: p0 = " << p0 << ", p1 = " << p1 << " per kHz" << std::endl;
  std::cout << "Pile-up band keeps " << (all > 0 ? 100.*inBand/all : 0) << "% of events" << std::endl;
  std::cout << "Wrote OO200_RefMultCorr_params.txt, hRefMultCorr_OO200.root and OO200_RMC_*.png" << std::endl;
}
