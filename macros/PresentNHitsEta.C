// PresentNHitsEta.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready QA plot of TPC hit counts vs pseudorapidity --
// N_hits^fit (hits actually used in the track fit) side by side with N_hits^max (the
// geometric ceiling on how many padrows that track's trajectory could possibly cross).
// A track's nHitsFit should stay roughly flat and well above the isGoodTrack() cut
// across the accepted eta range; a sag near the eta edges (where fewer TPC padrows are
// crossed) or an unexpected dip in the middle would flag a tracking/acceptance problem
// worth chasing before trusting anything built on top of these tracks. Plotting
// nHitsMax alongside it tells you whether a low nHitsFit region is acceptance-limited
// (nHitsMax itself is small there -- geometry, not a cut problem) or cut-limited
// (nHitsMax is fine, something else is throwing hits away) -- see CutClass's
// fitMaxRatio cut, which compares exactly these two quantities per track.
//
// Reads (exact names from PicoBinner.cxx -- see FindHistoInFile.C if your file uses
// different ones, e.g. because a_makeBasicHistos was false):
//   eta_nHitsFit -- TH2 of eta vs N_hits^fit, written at top level (same as etaPhi/
//                   eta_pT), before any outFile->cd() into a subdirectory.
//   eta_nHitsMax -- TH2 of eta vs N_hits^max, same top-level placement.
// Either one missing (e.g. an older file that predates eta_nHitsMax being added) is
// tolerated -- this only aborts if BOTH are absent.
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

  TH2* etaNHitsFit = (TH2*) inFile->Get("eta_nHitsFit");
  TH2* etaNHitsMax = (TH2*) inFile->Get("eta_nHitsMax");

  if(!etaNHitsFit && !etaNHitsMax){
    cout << "!!! Neither eta_nHitsFit nor eta_nHitsMax found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced, or this "
         << "file predates these histograms being added to PicoBinner.cxx)."
         << endl;
    return;
  }
  if(!etaNHitsFit) cout << "!!! eta_nHitsFit not found -- plotting eta_nHitsMax only." << endl;
  if(!etaNHitsMax) cout << "!!! eta_nHitsMax not found -- plotting eta_nHitsFit only." << endl;

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentNHitsEta","N_{hits}^{fit} & N_{hits}^{max} vs #eta",1600,700);
  c->Divide(2,1);

  // Extra right margin on both pads -- otherwise the z-axis (palette) power-of-10
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

  string outPath = a_outDir + "/PresentNHitsEta.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
