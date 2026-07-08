// PresentNHitsEta.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready plot of TPC hits-used-in-fit (N_hits^fit) vs
// pseudorapidity -- a standard tracking QA check. A track's hit count should stay
// roughly flat and well above the isGoodTrack() N_hits^fit cut across the accepted eta
// range; a sag near the eta edges (where fewer TPC padrows are crossed) or an
// unexpected dip in the middle would flag a tracking/acceptance problem worth chasing
// before trusting anything built on top of these tracks.
//
// Reads (exact name from PicoBinner.cxx -- see FindHistoInFile.C if your file uses a
// different one, e.g. because a_makeBasicHistos was false):
//   eta_nHitsFit -- TH2 of eta vs N_hits^fit, written at top level (same as etaPhi/
//                   eta_pT), before any outFile->cd() into a subdirectory.
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

  TH2* etaNHits = (TH2*) inFile->Get("eta_nHitsFit");

  if(!etaNHits){
    cout << "!!! eta_nHitsFit not found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced, or this "
         << "file predates the eta_nHitsFit histogram being added to PicoBinner.cxx)."
         << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentNHitsEta","N_{hits}^{fit} vs #eta",900,700);

  // Extra right margin -- otherwise the z-axis (palette) power-of-10 labels get
  // clipped at the edge of the pad.
  c->cd(); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  etaNHits->SetTitle("TPC Hits Used in Fit vs #eta");
  etaNHits->GetXaxis()->SetTitle("#eta");
  etaNHits->GetYaxis()->SetTitle("N_{hits}^{fit}");
  etaNHits->Draw("COLZ");

  string outPath = a_outDir + "/PresentNHitsEta.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
