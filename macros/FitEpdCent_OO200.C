// FitEpdCent_OO200.C -- lightflavorspectra_OO200
//
// Step 2 of an EPD (totnMIP) centrality for O+O 200 GeV (Run 21), trigger 860003: the EPD
// counterpart of macros/FitRefMultCorr_OO200.C. Reads the merged output of
// macros/MakeEpdCentInputs_OO200.C and derives (pass 1):
//
//   (a) Vz correction      c_vz(Vz) = <totnMIP>(|Vz| <= 2 cm) / <totnMIP>(Vz), from totnMIP >= kMipMean,
//                          in 1 cm bins over |Vz| <= 30 cm, plus a pol6 fit for reference.
//                          Cross-check: the 99th percentile of totnMIP per Vz bin.
//   (b) Luminosity corr.   <totnMIP>(totnMIP >= kMipMean) vs ZDCx for |Vz| <= 10 cm, p0 + p1*ZDCx;
//                          c_lumi(ZDCx) = p0 / (p0 + p1*ZDCx).
//   (c) Pile-up band       totnMIP vs nBTOFMatch: truncated mean +- kNSigma*RMS per slice, pol4 edges,
//                          continued linearly beyond the last fitted slice (same as the refMult band).
//                          Not applied to (d) here: the 3D input has no nBTOFMatch axis. Pass 2 of the
//                          inputs job would apply it, as was done for refMult.
//   (d) Corrected totnMIP  from hTotnMIP_Vz_ZDCx (refMult pile-up band already applied), each cell spread
//                          uniformly over its corrected interval. Written as hTotnMIP (650 unit bins) and,
//                          for a drop-in StNbdFitMaker input, as hRefMult (1000 unit bins, 0-1000).
//                          NOTE: StNbdFitMaker generates integer NBD multiplicities; totnMIP is a continuous
//                          sum of truncated nMIP with Landau-like fluctuations. Whether a plain NBD scan
//                          describes it is an open question; check how the official totnMIP table was made.
//   (e) Comparison with the official totnMIP cut table (triggers 860001/2/11/12, comparison only):
//                          fraction of recorded events above each official cut, for 860003 (raw and
//                          corrected) and for 860001/2. If the cut at X% is right and every event that
//                          central is recorded, frac(X) = X / eps, so X / frac estimates eps, the fraction
//                          of the inelastic cross section in the recorded sample (printed for 0-40%).
//   (f) QA                 <totnMIP> per run (runs > 3 sigma and > 2% from the mean are listed),
//                          totnMIP vs refMult, East vs West, totnMIP with and without row 1.
//
// Outputs:  OO200_EpdCent_params.txt   all parameters and tables, human readable
//           hTotnMIPCorr_OO200.root    hTotnMIP, hTotnMIPRaw, hTotnMIPCorr_Vz2, hTotnMIPOff, hRefMult (NBD input)
//           OO200_EPD_{VzCorr,LumiCorr,PileUp,TotnMIP,vsRefMult,EastWest,RunQA}.png
//
// Usage (container, SL24y, repo root):
//   root -l -b -q 'macros/FitEpdCent_OO200.C+("epd_OO200.root")'
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
#include "TProfile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"

namespace {
  const double kMipMean   = 100.;  // totnMIP threshold for the Vz and luminosity means (~ 50% centrality in the official table)
  const double kVzCenter  = 2.;    // |Vz| <= 2 cm defines the reference
  const double kNSigma    = 5.;    // pile-up band half-width in RMS
  const int    kTofSlice  = 2;     // nBTOFMatch bins per slice
  const double kMinSliceEntries = 2000.;
  const double kMinRunEntries   = 1000.;
  // official StRefMultCorr("totnmip") O+O cut table (star-sw main, Apr 2026): lower edge of each class
  const int    kNOff = 16;
  const double kOffPct[kNOff] = {80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5};
  const double kOffCut[kNOff] = {37, 46, 55, 66, 77, 91, 105, 122, 141, 162, 186, 213, 242, 276, 315, 367};

  // mean of a 1D histogram for x >= xmin (totnMIP is continuous: bin centers)
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
  // fraction of entries (in range, no under/overflow) with bin low edge >= cut
  double fracAbove(TH1* h, double cut){
    double tot = h->Integral(1, h->GetNbinsX());
    if(tot <= 0) return 0;
    int b = h->FindBin(cut + 1e-9);
    return h->Integral(b, h->GetNbinsX())/tot;
  }
  double bandFunc(double* x, double* p){
    const double xm = p[5];
    const double xx = (x[0] <= xm) ? x[0] : xm;
    double v = p[0] + p[1]*xx + p[2]*xx*xx + p[3]*xx*xx*xx + p[4]*xx*xx*xx*xx;
    if(x[0] > xm) v += (p[1] + 2*p[2]*xm + 3*p[3]*xm*xm + 4*p[4]*xm*xm*xm) * (x[0] - xm);
    return v;
  }
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

void FitEpdCent_OO200(const char* inName = "epd_OO200.root"){
  gStyle->SetOptStat(0);
  TFile* in = TFile::Open(inName);
  if(!in || in->IsZombie()){ std::cout << "Cannot open " << inName << std::endl; return; }
  TH2F* hVz   = (TH2F*) in->Get("hTotnMIP_vs_Vz");
  TH2F* hZdc  = (TH2F*) in->Get("hTotnMIP_vs_ZDCx");
  TH2F* hTof  = (TH2F*) in->Get("hTotnMIP_vs_nBTOFMatch");
  TH2F* hRef  = (TH2F*) in->Get("hTotnMIP_vs_refMult");
  TH3F* h3    = (TH3F*) in->Get("hTotnMIP_Vz_ZDCx");
  TH2F* hEW   = (TH2F*) in->Get("hTotnMIPEast_vs_West");
  TH2F* hNoR1 = (TH2F*) in->Get("hTotnMIPnoRow1_vs_Vz");
  TH2F* hVzOff = (TH2F*) in->Get("hTotnMIP_vs_Vz_off");
  TH2F* hTrig = (TH2F*) in->Get("hTrigOverlap");
  TProfile* pRun = (TProfile*) in->Get("pTotnMIP_vs_run");
  TH1D* hEv   = (TH1D*) in->Get("hEvents");
  if(!hVz || !hZdc || !hTof || !h3){ std::cout << "Missing input histograms in " << inName << std::endl; return; }

  std::ofstream par("OO200_EpdCent_params.txt");
  par << "# EPD totnMIP centrality inputs, pass 1, O+O 200 GeV Run 21, trigger 860003, |Vz| <= 30 cm, Vr <= 1 cm\n";
  if(hEv) par << "# events: 860003 " << hEv->GetBinContent(1) << ", |Vz|<=30 " << hEv->GetBinContent(2)
              << ", Vr<=1 " << hEv->GetBinContent(3) << ", refMult band " << hEv->GetBinContent(4)
              << ", of which totnMIP = 0: " << hEv->GetBinContent(5) << "\n";
  if(hTrig) par << "# trigger overlap (|Vz|<=30, Vr<=1): 860003 only " << hTrig->GetBinContent(2,1)
                << ", 860001/2 only " << hTrig->GetBinContent(1,2) << ", both " << hTrig->GetBinContent(2,2) << "\n";

  // ------------------------------------------------------------------ (a) Vz correction
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
  TGraphErrors* gVz = new TGraphErrors(); gVz->SetName("gVzCorr");
  TGraph* gQ = new TGraph(); gQ->SetName("gVzCorrQ99");
  std::vector<double> cvz(nVz, 1.);
  double maxDev = 0;
  for(int i = 0; i < nVz; i++){
    if(mean[i] > 0){
      cvz[i] = refCenter/mean[i];
      if(TMath::Abs(cvz[i] - 1.) > maxDev) maxDev = TMath::Abs(cvz[i] - 1.);
      int n = gVz->GetN(); gVz->SetPoint(n, vzC[i], cvz[i]); gVz->SetPointError(n, 0, cvz[i]*meanErr[i]/mean[i]);
    }
    if(q99[i] > 0){ int n = gQ->GetN(); gQ->SetPoint(n, vzC[i], q99Center/q99[i]); }
  }
  TF1* fVz = new TF1("fVz", "pol6", -30, 30);
  gVz->Fit(fVz, "Q0");
  par << "\n[Vz correction]  c_vz = <totnMIP>(|Vz|<=" << kVzCenter << ") / <totnMIP>(Vz), totnMIP >= " << kMipMean
      << ", reference <totnMIP> = " << refCenter << ", largest |c_vz - 1| = " << 100.*maxDev << "%\n";
  par << "pol6:";
  for(int k = 0; k < 7; k++) par << " " << fVz->GetParameter(k);
  par << "\n# Vz bin center (cm), c_vz (mean method), c_vz (99th-percentile method)\n";
  for(int i = 0; i < nVz; i++) par << Form("%6.1f  %.5f  %.5f\n", vzC[i], cvz[i], q99[i] > 0 ? q99Center/q99[i] : 0.);

  TCanvas* c1 = new TCanvas("c1", "", 900, 600);
  gVz->SetMarkerStyle(20); gVz->SetTitle("EPD Vz correction;V_{z} (cm);c_{vz} = #LTtotnMIP#GT_{center} / #LTtotnMIP#GT(V_{z})");
  gVz->Draw("AP"); fVz->SetLineColor(kRed); fVz->Draw("same");
  gQ->SetMarkerStyle(24); gQ->SetMarkerColor(kBlue); gQ->Draw("P same");
  TLegend* l1 = new TLegend(0.15, 0.75, 0.6, 0.88);
  l1->AddEntry(gVz, Form("mean, totnMIP #geq %.0f", kMipMean), "p"); l1->AddEntry(fVz, "pol6 fit", "l");
  l1->AddEntry(gQ, "99th percentile", "p"); l1->Draw();
  c1->Print("OO200_EPD_VzCorr.png");

  // ------------------------------------------------------------------ (b) luminosity correction
  TGraphErrors* gL = new TGraphErrors(); gL->SetName("gLumi");
  for(int i = 1; i <= hZdc->GetNbinsX(); i++){
    TH1D* p = hZdc->ProjectionY(Form("pZdc_%d", i), i, i);
    double e; double m = meanAbove(p, kMipMean, e);
    double nAbove = p->Integral(p->FindBin(kMipMean), p->GetNbinsX());
    if(nAbove > 1000 && m > 0){ int n = gL->GetN(); gL->SetPoint(n, hZdc->GetXaxis()->GetBinCenter(i), m); gL->SetPointError(n, 0, e); }
    delete p;
  }
  TF1* fL = new TF1("fL", "pol1", 0, 3);
  double p0 = refCenter, p1 = 0;
  if(gL->GetN() >= 2){ gL->Fit(fL, "Q0"); p0 = fL->GetParameter(0); p1 = fL->GetParameter(1); }
  par << "\n[Luminosity correction]  <totnMIP>(totnMIP >= " << kMipMean << ", |Vz| <= 10 cm) = p0 + p1*ZDCx[kHz]\n";
  par << "p0 = " << p0 << " +- " << fL->GetParError(0) << "\np1 = " << p1 << " +- " << fL->GetParError(1) << " per kHz"
      << "  (relative slope " << (p0 != 0 ? 100.*p1/p0 : 0) << "% per kHz)\n";
  TCanvas* c2 = new TCanvas("c2", "", 900, 600);
  gL->SetMarkerStyle(20); gL->SetTitle(Form("EPD luminosity dependence, |V_{z}| #leq 10 cm;ZDCx (kHz);#LTtotnMIP#GT (totnMIP #geq %.0f)", kMipMean));
  gL->Draw("AP"); fL->SetLineColor(kRed); fL->Draw("same");
  c2->Print("OO200_EPD_LumiCorr.png");

  // ------------------------------------------------------------------ (c) pile-up band, totnMIP vs nBTOFMatch
  TGraph* gUp = new TGraph(); TGraph* gLo = new TGraph(); TGraph* gMu = new TGraph();
  const int nTof = hTof->GetNbinsX();
  for(int i = 1; i <= nTof; i += kTofSlice){
    TH1D* p = hTof->ProjectionY(Form("pTof_%d", i), i, i + kTofSlice - 1);
    if(p->GetEntries() < kMinSliceEntries){ delete p; continue; }
    double m = p->GetMean(), s = p->GetRMS();
    for(int it = 0; it < 5; it++){
      p->GetXaxis()->SetRangeUser(TMath::Max(0., m - 3*s), m + 3*s);
      m = p->GetMean(); s = p->GetRMS();
    }
    double x = 0.5*(hTof->GetXaxis()->GetBinLowEdge(i) + hTof->GetXaxis()->GetBinUpEdge(i + kTofSlice - 1));
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
  par << "\n[Pile-up band]  totnMIP vs nBTOFMatch, mean +- " << kNSigma << " x truncated RMS per " << kTofSlice
      << "-unit nBTOFMatch slice (slices with >= " << kMinSliceEntries << " events), pol4\n";
  par << "upper (b0..b4):"; for(int k = 0; k < 5; k++) par << " " << fUp->GetParameter(k); par << "\n";
  par << "lower (c0..c4):"; for(int k = 0; k < 5; k++) par << " " << fLo->GetParameter(k); par << "\n";
  par << "fit range in nBTOFMatch: 0 - " << xmFit << " (beyond it: straight line, tangent at the last fitted point)\n";
  double inBand = 0, all = 0;
  for(int i = 1; i <= nTof; i++){
    double x = hTof->GetXaxis()->GetBinCenter(i);
    for(int j = 1; j <= hTof->GetNbinsY(); j++){
      double y = hTof->GetYaxis()->GetBinCenter(j), w = hTof->GetBinContent(i, j);
      all += w; if(y >= fLoE->Eval(x) && y <= fUpE->Eval(x)) inBand += w;
    }
  }
  par << "events inside the band: " << inBand << " of " << all << " (" << (all > 0 ? 100.*inBand/all : 0) << "%)\n";
  TCanvas* c3 = new TCanvas("c3", "", 900, 700); c3->SetLogz();
  hTof->SetTitle("totnMIP vs nBTOFMatch, |V_{z}| #leq 30 cm;nBTOFMatch;totnMIP");
  hTof->GetXaxis()->SetRangeUser(0, 150);
  hTof->Draw("colz");
  fUpE->SetLineColor(kRed); fLoE->SetLineColor(kRed); fUpE->SetNpx(600); fLoE->SetNpx(600); fUpE->Draw("same"); fLoE->Draw("same");
  TLegend* l3 = new TLegend(0.15, 0.82, 0.5, 0.88);
  l3->AddEntry(fUpE, Form("pass-1 band (#pm%.0f#sigma)", kNSigma), "l"); l3->Draw();
  c3->Print("OO200_EPD_PileUp.png");

  // ------------------------------------------------------------------ (d) corrected totnMIP
  const int nMipBins = h3->GetNbinsX();
  const double mipMax = h3->GetXaxis()->GetXmax();
  TH1D* hCorr  = new TH1D("hTotnMIP", "Corrected totnMIP (Vz, luminosity), refMult band, |V_{z}| #leq 30 cm;totnMIP_{corr};Events", nMipBins, 0, mipMax);
  TH1D* hRaw   = new TH1D("hTotnMIPRaw", "Uncorrected totnMIP, same events;totnMIP;Events", nMipBins, 0, mipMax);
  TH1D* hCorr2 = new TH1D("hTotnMIPCorr_Vz2", "Corrected totnMIP, |V_{z}| #leq 2 cm;totnMIP_{corr};Events", nMipBins, 0, mipMax);
  for(int iz = 1; iz <= h3->GetNbinsZ(); iz++){
    double zdc = h3->GetZaxis()->GetBinCenter(iz);
    double cl = (p0 + p1*zdc > 0) ? p0/(p0 + p1*zdc) : 1.;
    for(int iy = 1; iy <= h3->GetNbinsY(); iy++){
      double vz = h3->GetYaxis()->GetBinCenter(iy);
      int ivz = hVz->GetXaxis()->FindBin(vz);
      double c = cl * ((ivz >= 1 && ivz <= nVz) ? cvz[ivz-1] : 1.);
      for(int ix = 1; ix <= nMipBins; ix++){
        double w = h3->GetBinContent(ix, iy, iz);
        if(w <= 0) continue;
        double lo = h3->GetXaxis()->GetBinLowEdge(ix), hi = h3->GetXaxis()->GetBinUpEdge(ix);
        hRaw->Fill(0.5*(lo + hi), w);
        spread(hCorr, lo*c, hi*c, w);
        if(TMath::Abs(vz) <= kVzCenter) spread(hCorr2, lo*c, hi*c, w);
      }
    }
  }
  // drop-in input for StNbdFitMaker (same name/binning as hRefMult_OO200.root): see the NOTE at the top
  TH1D* hNbd = new TH1D("hRefMult", "Corrected totnMIP for the NBD scan;totnMIP_{corr};Events", 1000, 0, 1000);
  for(int i = 1; i <= hCorr->GetNbinsX(); i++) hNbd->SetBinContent(hNbd->FindBin(hCorr->GetBinCenter(i)), hCorr->GetBinContent(i));
  hNbd->SetEntries(hCorr->Integral());
  // official-trigger distribution (no refMult band, no corrections), comparison only
  TH1D* hOff = hVzOff ? hVzOff->ProjectionY("hTotnMIPOff") : 0;
  TH1D* hRawNoBand = hVz->ProjectionY("hTotnMIPRawNoBand");

  // ------------------------------------------------------------------ (e) official cut table comparison
  par << "\n[Official totnMIP cut table (860001/2/11/12), comparison only]\n";
  par << "# fraction of recorded events with totnMIP >= cut; eps = X/frac for X <= 40%\n";
  par << "#  class   cut   frac(860003 corr)  frac(860003 raw, no band)  frac(860001/2 raw)   eps(860003)  eps(860001/2)\n";
  for(int k = kNOff - 1; k >= 0; k--){
    double fC = fracAbove(hCorr, kOffCut[k]);
    double fR = fracAbove(hRawNoBand, kOffCut[k]);
    double fO = hOff ? fracAbove(hOff, kOffCut[k]) : 0;
    TString e1 = (kOffPct[k] <= 40 && fC > 0) ? Form("%.3f", 0.01*kOffPct[k]/fC) : "  -  ";
    TString e2 = (kOffPct[k] <= 40 && fO > 0) ? Form("%.3f", 0.01*kOffPct[k]/fO) : "  -  ";
    par << Form("  0-%2.0f%%  %4.0f      %.4f             %.4f                     %.4f            %s        %s\n",
                kOffPct[k], kOffCut[k], fC, fR, fO, e1.Data(), e2.Data());
  }
  par << "# reading: eps roughly constant over 0-40% means the shapes agree there; eps(860003) vs eps(860001/2)\n"
      << "# compares how much of the cross section each trigger records (assumes the official percentiles are right).\n";

  // ------------------------------------------------------------------ (f) run QA
  TGraphErrors* gRun = new TGraphErrors(); gRun->SetName("gRunQA");
  if(pRun){
    double s = 0, sw = 0;
    for(int i = 1; i <= pRun->GetNbinsX(); i++){
      if(pRun->GetBinEntries(i) < kMinRunEntries || pRun->GetBinError(i) <= 0) continue;
      double e = pRun->GetBinError(i);
      s += pRun->GetBinContent(i)/(e*e); sw += 1./(e*e);
    }
    const double runMean = (sw > 0) ? s/sw : 0;
    par << "\n[Run QA]  <totnMIP> per run (|Vz| <= 30 cm, no band), weighted mean over runs " << runMean << "\n";
    par << "# runs more than 3 sigma AND more than 2% from the mean:\n";
    int nFlag = 0;
    for(int i = 1; i <= pRun->GetNbinsX(); i++){
      if(pRun->GetBinEntries(i) < kMinRunEntries || pRun->GetBinError(i) <= 0) continue;
      int idx = (int) pRun->GetBinLowEdge(i);
      int runId = 22000000 + (130 + idx/100)*1000 + idx%100;
      double m = pRun->GetBinContent(i), e = pRun->GetBinError(i);
      int n = gRun->GetN(); gRun->SetPoint(n, idx, m); gRun->SetPointError(n, 0, e);
      double dev = (runMean > 0) ? (m - runMean)/runMean : 0;
      if(TMath::Abs(m - runMean) > 3*e && TMath::Abs(dev) > 0.02){
        par << Form("  run %d  <totnMIP> = %.2f +- %.2f  (%+.1f%%, %.0f events)\n", runId, m, e, 100.*dev, pRun->GetBinEntries(i));
        nFlag++;
      }
    }
    par << "# flagged runs: " << nFlag << " of " << gRun->GetN() << "\n";
    TCanvas* c7 = new TCanvas("c7", "", 1100, 500);
    gRun->SetMarkerStyle(20); gRun->SetMarkerSize(0.7);
    gRun->SetTitle("#LTtotnMIP#GT per run, 860003, |V_{z}| #leq 30 cm;run index = (day-130)*100 + run-in-day;#LTtotnMIP#GT");
    gRun->Draw("AP");
    if(gRun->GetN() > 0){
      TLine* lm = new TLine(gRun->GetXaxis()->GetXmin(), runMean, gRun->GetXaxis()->GetXmax(), runMean);
      lm->SetLineColor(kRed); lm->Draw();
    }
    c7->Print("OO200_EPD_RunQA.png");
  }

  // ------------------------------------------------------------------ plots: distributions and correlations
  TCanvas* c4 = new TCanvas("c4", "", 900, 600); c4->SetLogy();
  TH1D* r1 = (TH1D*) hRaw->Clone("r1");  r1->Scale(1./r1->Integral());
  TH1D* c1n = (TH1D*) hCorr->Clone("c1n"); c1n->Scale(1./c1n->Integral());
  r1->SetLineColor(kBlack); c1n->SetLineColor(kRed);
  r1->SetTitle("totnMIP, 860003;totnMIP;normalized");
  r1->GetXaxis()->SetRangeUser(0, 600); r1->Draw("hist"); c1n->Draw("hist same");
  TLegend* l4 = new TLegend(0.55, 0.7, 0.88, 0.88);
  l4->AddEntry(r1, "860003 raw", "l"); l4->AddEntry(c1n, "860003 corrected", "l");
  if(hOff && hOff->Integral() > 0){
    TH1D* o1 = (TH1D*) hOff->Clone("o1"); o1->Scale(1./o1->Integral()); o1->SetLineColor(kBlue); o1->SetLineStyle(2);
    o1->Draw("hist same"); l4->AddEntry(o1, "860001/2 raw (comparison)", "l");
  }
  l4->Draw();
  c4->Print("OO200_EPD_TotnMIP.png");

  if(hRef){
    TCanvas* c5 = new TCanvas("c5", "", 900, 700); c5->SetLogz();
    hRef->SetTitle("totnMIP vs refMult, 860003, refMult band;refMult;totnMIP");
    hRef->GetXaxis()->SetRangeUser(0, 80); hRef->GetYaxis()->SetRangeUser(0, 600);
    hRef->Draw("colz");
    TProfile* pr = hRef->ProfileX("pTotnMIP_vs_refMult_prof"); pr->SetMarkerStyle(20); pr->SetMarkerSize(0.6); pr->Draw("same");
    c5->Print("OO200_EPD_vsRefMult.png");
    // totnMIP for refMult = 0 and 1: what the 70-100% class looks like in the EPD
    TH1D* r0 = hRef->ProjectionY("hTotnMIP_refMult0", 1, 1);
    TH1D* r01 = hRef->ProjectionY("hTotnMIP_refMult01", 1, 2);
    par << "\n[refMult = 0 and 0-1 events in the EPD]\n";
    par << "refMult = 0:   " << r0->GetEntries() << " events, <totnMIP> = " << r0->GetMean() << ", fraction with totnMIP = 0: "
        << (r0->Integral() > 0 ? r0->GetBinContent(1)/r0->Integral() : 0) << "\n";
    par << "refMult = 0-1: " << r01->GetEntries() << " events, <totnMIP> = " << r01->GetMean() << "\n";
  }
  if(hEW){
    TCanvas* c6 = new TCanvas("c6", "", 800, 700); c6->SetLogz();
    hEW->GetXaxis()->SetRangeUser(0, 350); hEW->GetYaxis()->SetRangeUser(0, 350);
    hEW->Draw("colz");
    c6->Print("OO200_EPD_EastWest.png");
    par << "\n[East/West]  <West> = " << hEW->GetMean(1) << ", <East> = " << hEW->GetMean(2)
        << ", correlation = " << hEW->GetCorrelationFactor() << "\n";
  }
  if(hNoR1){
    TH1D* nr = hNoR1->ProjectionY("hTotnMIPnoRow1");
    par << "[Row 1]  <totnMIP without row 1> = " << nr->GetMean() << " vs <totnMIP> (same events) = " << hRaw->GetMean()
        << "  -> row 1 carries " << (hRaw->GetMean() > 0 ? 100.*(1. - nr->GetMean()/hRaw->GetMean()) : 0) << "% of the sum\n";
  }

  TFile* out = TFile::Open("hTotnMIPCorr_OO200.root", "RECREATE");
  hCorr->Write(); hRaw->Write(); hCorr2->Write(); hNbd->Write(); hRawNoBand->Write();
  if(hOff) hOff->Write();
  gVz->Write(); gQ->Write(); gL->Write(); gRun->Write();
  gUp->Write("gPileUpUpper"); gLo->Write("gPileUpLower"); fUpE->Write("fPileUpUpper"); fLoE->Write("fPileUpLower");
  fVz->Write(); fL->Write(); fUp->Write(); fLo->Write();
  out->Close();
  par << "\n[Corrected totnMIP]  hTotnMIP: " << hCorr->Integral() << " events (|Vz| <= 30 cm), "
      << hCorr2->Integral() << " with |Vz| <= 2 cm; <totnMIP> raw " << hRaw->GetMean() << ", corrected " << hCorr->GetMean() << "\n";
  par.close();

  std::cout << "Reference <totnMIP> (totnMIP >= " << kMipMean << ", |Vz| <= " << kVzCenter << "): " << refCenter
            << ", largest Vz correction " << 100.*maxDev << "%" << std::endl;
  std::cout << "Luminosity: p0 = " << p0 << ", p1 = " << p1 << " per kHz" << std::endl;
  std::cout << "totnMIP pile-up band keeps " << (all > 0 ? 100.*inBand/all : 0) << "% of events" << std::endl;
  std::cout << "Wrote OO200_EpdCent_params.txt, hTotnMIPCorr_OO200.root and OO200_EPD_*.png" << std::endl;
}
