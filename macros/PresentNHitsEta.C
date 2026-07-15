// PresentNHitsEta.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready QA plot of TPC hit counts vs pseudorapidity --
// N_hits^fit (hits actually used in the track fit), N_hits^max (the geometric ceiling
// on how many padrows that track's trajectory could possibly cross), and their ratio
// (hit-finding efficiency), all three side by side so you can tell acceptance and
// efficiency effects apart at a glance instead of needing two separate macros/plots.
//
// A track's nHitsFit should stay roughly flat and well above the isGoodTrack() cut
// across the accepted eta range; a sag near the eta edges (where fewer TPC padrows are
// crossed) or an unexpected dip in the middle would flag a tracking/acceptance problem
// worth chasing before trusting anything built on top of these tracks. Plotting
// nHitsMax alongside it tells you whether a low nHitsFit region is acceptance-limited
// (nHitsMax itself is small there -- geometry, not a cut problem) or cut-limited
// (nHitsMax is fine, something else is throwing hits away).
//
// The third panel (fitMaxRatio = N_hits^fit / N_hits^max, TRUE double division
// regardless of how CutClass's fitMaxRatio cut is currently enforced) isolates hit-
// finding EFFICIENCY from the pure geometric acceptance already covered by the first
// two panels. eta_nHitsMax shows a hard kinematic ceiling that falls off toward large
// |eta| -- that's expected TPC acceptance, not a problem. This panel answers a
// different question: for the hits that were geometrically available, what fraction
// did the fit actually use? If THIS ratio also sags at large |eta| (rather than
// staying flat until the ceiling itself forces it down), that points to a real
// hit-finding efficiency loss at those eta values -- e.g. multiple scattering or
// energy loss in whatever beampipe material the track crossed at that eta (STAR's
// beampipe isn't one uniform material along z), knocking real, geometrically-
// available hits out of the fit. That's a distinct failure mode from acceptance.
//
// Reads (exact names from PicoBinner.cxx -- see FindHistoInFile.C if your file uses
// different ones, e.g. because a_makeBasicHistos was false):
//   eta_nHitsFit    -- TH2 of eta vs N_hits^fit, written at top level (same as etaPhi/
//                      eta_pT), before any outFile->cd() into a subdirectory.
//   eta_nHitsMax    -- TH2 of eta vs N_hits^max, same top-level placement.
//   eta_fitMaxRatio -- TH2D of eta vs N_hits^fit/N_hits^max, same top-level placement.
// Any of the three missing (e.g. an older file that predates one being added) is
// tolerated -- this only aborts if ALL THREE are absent.
//
// Usage:
//   root -l -q 'PresentNHitsEta.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentNHitsEta.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentNHitsEta(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* etaNHitsFit    = (TH2*) inFile->Get("eta_nHitsFit");
  TH2* etaNHitsMax    = (TH2*) inFile->Get("eta_nHitsMax");
  TH2* etaFitMaxRatio = (TH2*) inFile->Get("eta_fitMaxRatio");

  if(!etaNHitsFit && !etaNHitsMax && !etaFitMaxRatio){
    cout << "!!! None of eta_nHitsFit, eta_nHitsMax, eta_fitMaxRatio found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced, or this "
         << "file predates these histograms being added to PicoBinner.cxx)."
         << endl;
    return;
  }
  if(!etaNHitsFit)    cout << "!!! eta_nHitsFit not found -- skipping that panel." << endl;
  if(!etaNHitsMax)    cout << "!!! eta_nHitsMax not found -- skipping that panel." << endl;
  if(!etaFitMaxRatio) cout << "!!! eta_fitMaxRatio not found -- skipping that panel "
                            << "(this file may predate it being added to PicoBinner.cxx)." << endl;

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentNHitsEta","N_{hits}^{fit}, N_{hits}^{max} & their ratio vs #eta",2400,700);
  c->Divide(3,1);

  // Extra right margin on all pads -- otherwise the z-axis (palette) power-of-10
  // labels get clipped at the edge of the pad.
  c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(etaNHitsFit){
    etaNHitsFit->SetTitle("TPC Hits Used in Fit vs #eta");
    etaNHitsFit->GetXaxis()->SetTitle("#eta");
    etaNHitsFit->GetYaxis()->SetTitle("N_{hits}^{fit}");
    etaNHitsFit->Draw("COLZ");
  }

  c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(etaNHitsMax){
    etaNHitsMax->SetTitle("TPC NHitsMax vs #eta");
    etaNHitsMax->GetXaxis()->SetTitle("#eta");
    etaNHitsMax->GetYaxis()->SetTitle("N_{hits}^{max}");
    etaNHitsMax->Draw("COLZ");
  }

  c->cd(3); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(etaFitMaxRatio){
    etaFitMaxRatio->SetTitle("Hit-Finding Efficiency vs #eta");
    etaFitMaxRatio->GetXaxis()->SetTitle("#eta");
    etaFitMaxRatio->GetYaxis()->SetTitle("N_{hits}^{fit}/N_{hits}^{max}");
    etaFitMaxRatio->Draw("COLZ");
  }

  string outPath = a_outDir + "/PresentNHitsEta.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
