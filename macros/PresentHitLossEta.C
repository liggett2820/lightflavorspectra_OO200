// PresentHitLossEta.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: tests whether the fitMaxRatio sag at large |eta| (see
// PresentNHitsEta.C's third panel) is a real, growing hit-finding inefficiency, or
// just a shrinking-denominator artifact of NHitsMax collapsing toward the TPC's
// acceptance edge. The ratio NHitsFit/NHitsMax will drop near the edge either way --
// a small, constant absolute hit loss looks huge once NHitsMax itself is small. This
// macro instead plots the ABSOLUTE hit loss, <NHitsMax> - <NHitsFit> vs eta: if that's
// flat, the ratio's decline is mostly geometric bookkeeping, not a growing physical
// effect (e.g. beampipe-material multiple scattering/energy loss knocking out hits).
// If it visibly grows toward large |eta|, that's real evidence worth chasing further
// (see also: splitting fitMaxRatio by pT, since multiple scattering scales as 1/p).
//
// Why profile subtraction is valid here (not just approximate): eta_nHitsFit and
// eta_nHitsMax are filled from the SAME track, with the SAME eta value, in the SAME
// loop iteration in PicoBinner.cxx (both gated on a_makeBasicHistos, back to back --
// see the "FILL TRACK HISTOGRAMS" block). So every entry in a given eta bin of one
// histogram has a corresponding entry from the identical track in that same eta bin
// of the other. By linearity of expectation, mean(NHitsMax) - mean(NHitsFit) per eta
// bin equals mean(NHitsMax - NHitsFit) per eta bin -- this is exact, not approximate,
// even though we never booked a directly-correlated (eta, NHitsMax-NHitsFit)
// histogram. That means this test can run right now against existing yield files --
// no PicoBinner.cxx change or rerun required.
//
// Reads (exact names from PicoBinner.cxx, written at top level like eta_pT):
//   eta_nHitsFit -- TH2 of eta vs N_hits^fit
//   eta_nHitsMax -- TH2 of eta vs N_hits^max
// Missing either one aborts with a clear message -- see FindHistoInFile.C to check
// what's actually in your file.
//
// Usage:
//   root -l -q 'PresentHitLossEta.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentHitLossEta.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
using namespace std;

void PresentHitLossEta(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* etaNHitsFit = (TH2*) inFile->Get("eta_nHitsFit");
  TH2* etaNHitsMax = (TH2*) inFile->Get("eta_nHitsMax");

  if(!etaNHitsFit || !etaNHitsMax){
    cout << "!!! Need both eta_nHitsFit and eta_nHitsMax in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced)."
         << endl;
    return;
  }

  // ProfileX gives <NHitsFit>(eta) and <NHitsMax>(eta) with proper errors-on-the-mean
  // per eta bin -- exactly the per-bin quantities we need to difference.
  TProfile* profFit = etaNHitsFit->ProfileX("prof_nHitsFit");
  TProfile* profMax = etaNHitsMax->ProfileX("prof_nHitsMax");

  int nBins = profFit->GetNbinsX();
  if(profMax->GetNbinsX() != nBins){
    cout << "!!! eta_nHitsFit and eta_nHitsMax have different eta binning -- "
         << "can't subtract bin-by-bin. Check PicoBinner.cxx booking." << endl;
    return;
  }

  TH1D* hitLoss = new TH1D("hitLossEta",
    ";#eta;#LT N_{hits}^{max} #GT - #LT N_{hits}^{fit} #GT",
    nBins, profFit->GetXaxis()->GetXmin(), profFit->GetXaxis()->GetXmax());

  for(int iBin = 1; iBin <= nBins; iBin++){
    // Skip empty bins -- an empty profile bin has content 0 and error 0, which would
    // otherwise show up as a spurious "zero hit loss" point instead of a gap.
    if(profFit->GetBinEntries(iBin) < 1 || profMax->GetBinEntries(iBin) < 1) continue;

    double diff    = profMax->GetBinContent(iBin) - profFit->GetBinContent(iBin);
    double errFit  = profFit->GetBinError(iBin);
    double errMax  = profMax->GetBinError(iBin);
    double errDiff = TMath::Sqrt(errFit*errFit + errMax*errMax);

    hitLoss->SetBinContent(iBin, diff);
    hitLoss->SetBinError(iBin, errDiff);
  }

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas("PresentHitLossEta","Absolute hit loss vs #eta",900,700);
  hitLoss->SetTitle("Absolute TPC Hit Loss vs #eta");
  hitLoss->SetMarkerStyle(20);
  hitLoss->SetMarkerColor(kBlue+1);
  hitLoss->SetLineColor(kBlue+1);
  hitLoss->GetYaxis()->SetTitleOffset(1.3);
  hitLoss->Draw("E1");

  string outPath = a_outDir + "/PresentHitLossEta.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
