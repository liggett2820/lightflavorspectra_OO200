// PresentNHitsMaxEta.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: companion QA plot to PresentNHitsEta.C. N_hits^max is the geometric
// ceiling on how many TPC padrows a track *could* cross given its eta/pT (the
// denominator of CutClass's nHitsFit/nHitsMax fitMaxRatio cut) -- plotting it vs eta
// tells you whether a low-yield region (e.g. the large-|eta|, low-pT corner) is
// acceptance-limited (nHitsMax itself is small there, so no cut choice can recover it)
// or cut-limited (nHitsMax is plenty large, but nHitsFit/the fitMaxRatio cut is throwing
// tracks away that the detector otherwise could have reconstructed). Compare this plot's
// shape directly against PresentNHitsEta.C's eta_nHitsFit plot for the same file.
//
// Reads (exact name from PicoBinner.cxx -- see FindHistoInFile.C if your file uses a
// different one, e.g. because a_makeBasicHistos was false):
//   eta_nHitsMax -- TH2 of eta vs N_hits^max, written at top level (same as etaPhi/
//                   eta_pT/eta_nHitsFit), before any outFile->cd() into a subdirectory.
//
// Usage:
//   root -l -q 'PresentNHitsMaxEta.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentNHitsMaxEta.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentNHitsMaxEta(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* etaNHitsMax = (TH2*) inFile->Get("eta_nHitsMax");

  if(!etaNHitsMax){
    cout << "!!! eta_nHitsMax not found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced, or this "
         << "file predates the eta_nHitsMax histogram being added to PicoBinner.cxx)."
         << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentNHitsMaxEta","N_{hits}^{max} vs #eta",900,700);

  // Extra right margin -- otherwise the z-axis (palette) power-of-10 labels get
  // clipped at the edge of the pad.
  c->cd(); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  etaNHitsMax->SetTitle("TPC Geometric Max Hits vs #eta");
  etaNHitsMax->GetXaxis()->SetTitle("#eta");
  etaNHitsMax->GetYaxis()->SetTitle("N_{hits}^{max}");
  etaNHitsMax->Draw("COLZ");

  string outPath = a_outDir + "/PresentNHitsMaxEta.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
