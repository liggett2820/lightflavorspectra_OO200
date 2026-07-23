// CompareNHitsMaxByCharge.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: direct test of the leading hypothesis for the 2-3 diagonal
// "streaks" seen along each side of the eta_nHitsMax trapezoid (see
// PresentNHitsEta.C's middle panel, and ScanNHitsMaxVsEta.C for the already-
// understood main-ceiling / inner-outer-boundary / TPC-corner features). If those
// streaks are genuine fixed TPC structural crossings (e.g. field-cage support ribs,
// wheel-mount hardware) smeared into diagonal lines by track curvature, then positive
// and negative tracks -- which bend oppositely in STAR's field -- should cross the
// same structure at slightly different eta, shifting the streaks between charge
// signs. A purely eta-only geometric effect (no curvature dependence at all) would
// show identical streaks in both.
//
// PREREQUISITE: needs eta_nHitsMax_chargePlus/eta_nHitsMax_chargeMinus, added to
// PicoBinner.cxx 2026-07-23 -- rerun PicoBinner against a species file to populate
// these; they won't exist in any yield file produced before that change (run
// FindHistoInFile.C to check what's actually in a given file).
//
// Draws three things:
//   (1) the two 2D histograms (chargePlus, chargeMinus) side by side, same z-scale,
//       so you can visually compare the streak positions directly.
//   (2) the same continuous peak-vs-eta scan ScanNHitsMaxVsEta.C does, but as two
//       overlaid curves (blue=positive, red=negative) restricted to the eta range
//       where the streaks live (roughly |eta| in [0.9,1.6], past the TPC corner) --
//       this is the quantitative version of (1): if the ceiling trace itself is
//       offset in eta between the two charges, that's the curvature-shift signature.
//
// Usage:
//   root -l -q 'CompareNHitsMaxByCharge.C("/Users/aliggett/data/OO/yieldHistos_OO200_pion.root")'

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

// Same sliding-window peak scan as ScanNHitsMaxVsEta.C, factored out so it can run
// on both charge signs and return a TGraph for overlay comparison.
TGraph* scanPeakVsEta(TH2* a_hist, double a_etaLow, double a_etaHigh,
                       double a_etaStep, double a_etaHalfWidth, double a_minEntries){
  vector<double> etaVals, peakVals;
  for(double etaCenter = a_etaLow; etaCenter <= a_etaHigh; etaCenter += a_etaStep){
    int binLow  = a_hist->GetXaxis()->FindBin(etaCenter - a_etaHalfWidth);
    int binHigh = a_hist->GetXaxis()->FindBin(etaCenter + a_etaHalfWidth);
    TH1D* proj = a_hist->ProjectionY("chargeScanTmpProj", binLow, binHigh);
    double totalEntries = proj->Integral();
    if(totalEntries >= a_minEntries){
      int peakBin = proj->GetMaximumBin();
      double peakVal = proj->GetXaxis()->GetBinCenter(peakBin);
      etaVals.push_back(etaCenter);
      peakVals.push_back(peakVal);
    }
    delete proj;
  }
  if(etaVals.empty()) return nullptr;
  return new TGraph((int)etaVals.size(), &etaVals[0], &peakVals[0]);
}

void CompareNHitsMaxByCharge(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* hPlus  = (TH2*) inFile->Get("eta_nHitsMax_chargePlus");
  TH2* hMinus = (TH2*) inFile->Get("eta_nHitsMax_chargeMinus");
  if(!hPlus || !hMinus){
    cout << "!!! eta_nHitsMax_chargePlus/Minus not found in " << a_inputFile
         << " -- this file predates the 2026-07-23 PicoBinner.cxx change that adds "
         << "them. Rerun PicoBinner against this species first, or run "
         << "FindHistoInFile.C to double check what's actually in this file."
         << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outDir.c_str(), true);

  // (1) side-by-side 2D comparison, shared z-scale so colors mean the same thing in
  // both panels.
  double globalMax = TMath::Max(hPlus->GetMaximum(), hMinus->GetMaximum());
  hPlus->SetMaximum(globalMax);
  hMinus->SetMaximum(globalMax);
  hPlus->SetMinimum(1);
  hMinus->SetMinimum(1);

  TCanvas* c1 = new TCanvas("CompareNHitsMaxByCharge_2D","N_{hits}^{max} vs #eta by charge sign",1600,800);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogz(1);
  hPlus->SetTitle("Positive tracks;#eta;N_{hits}^{max}");
  hPlus->Draw("colz");
  c1->cd(2);
  gPad->SetLogz(1);
  hMinus->SetTitle("Negative tracks;#eta;N_{hits}^{max}");
  hMinus->Draw("colz");
  string outPath2D = a_outDir + "/CompareNHitsMaxByCharge_2D.png";
  c1->SaveAs(outPath2D.c_str());
  cout << "Wrote " << outPath2D << endl;

  // (2) overlaid peak-vs-eta scan, restricted to the streak region (roughly
  // |eta| in [0.9,1.6] on the positive side, mirrored on the negative side).
  double etaStep = 0.02, etaHalfWidth = 0.025, minEntries = 500;

  TGraph* peakPlusRight  = scanPeakVsEta(hPlus,  0.90, 1.65, etaStep, etaHalfWidth, minEntries);
  TGraph* peakMinusRight = scanPeakVsEta(hMinus, 0.90, 1.65, etaStep, etaHalfWidth, minEntries);
  TGraph* peakPlusLeft   = scanPeakVsEta(hPlus,  -1.65, -0.90, etaStep, etaHalfWidth, minEntries);
  TGraph* peakMinusLeft  = scanPeakVsEta(hMinus, -1.65, -0.90, etaStep, etaHalfWidth, minEntries);

  TCanvas* c2 = new TCanvas("CompareNHitsMaxByCharge_Scan","peak N_{hits}^{max} vs #eta, by charge sign",1600,800);
  c2->Divide(2,1);

  c2->cd(1);
  TH2F* frameLeft = new TH2F("frameLeft",";#eta;peak (mode) N_{hits}^{max}",10,-1.7,-0.85,10,0,80);
  frameLeft->Draw();
  if(peakPlusLeft){ peakPlusLeft->SetMarkerStyle(20); peakPlusLeft->SetMarkerColor(kBlue+2); peakPlusLeft->SetLineColor(kBlue+2); peakPlusLeft->Draw("P same"); }
  if(peakMinusLeft){ peakMinusLeft->SetMarkerStyle(20); peakMinusLeft->SetMarkerColor(kRed+1); peakMinusLeft->SetLineColor(kRed+1); peakMinusLeft->Draw("P same"); }
  TLegend* legLeft = new TLegend(0.15,0.15,0.45,0.3);
  legLeft->AddEntry(peakPlusLeft,"positive tracks","p");
  legLeft->AddEntry(peakMinusLeft,"negative tracks","p");
  legLeft->Draw();

  c2->cd(2);
  TH2F* frameRight = new TH2F("frameRight",";#eta;peak (mode) N_{hits}^{max}",10,0.85,1.7,10,0,80);
  frameRight->Draw();
  if(peakPlusRight){ peakPlusRight->SetMarkerStyle(20); peakPlusRight->SetMarkerColor(kBlue+2); peakPlusRight->SetLineColor(kBlue+2); peakPlusRight->Draw("P same"); }
  if(peakMinusRight){ peakMinusRight->SetMarkerStyle(20); peakMinusRight->SetMarkerColor(kRed+1); peakMinusRight->SetLineColor(kRed+1); peakMinusRight->Draw("P same"); }
  TLegend* legRight = new TLegend(0.55,0.7,0.85,0.85);
  legRight->AddEntry(peakPlusRight,"positive tracks","p");
  legRight->AddEntry(peakMinusRight,"negative tracks","p");
  legRight->Draw();

  string outPathScan = a_outDir + "/CompareNHitsMaxByCharge_Scan.png";
  c2->SaveAs(outPathScan.c_str());
  cout << "Wrote " << outPathScan << endl;

  cout << "\nIf positive/negative points visibly separate in eta (rather than lying "
       << "on top of each other) in the scan plot, that's the curvature-shift "
       << "signature -- the streaks come from crossing a fixed structural boundary "
       << "at a curvature-dependent eta. If they overlap, the streaks aren't "
       << "curvature-driven and the boundary hypothesis needs a different test "
       << "(e.g. pT-binned instead of just charge-sign)." << endl;
}
