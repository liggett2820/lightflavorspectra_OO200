// PresentFitMaxRatioEta.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: isolates hit-finding EFFICIENCY (N_hits^fit / N_hits^max) from the
// pure geometric acceptance already covered by PresentNHitsEta.C. eta_nHitsMax shows a
// hard kinematic ceiling that falls off toward large |eta| -- that's expected TPC
// acceptance, not a problem. This plot answers a different question: for the hits that
// were geometrically available, what fraction did the fit actually use? If THIS ratio
// also sags at large |eta| (rather than staying flat until the ceiling itself forces it
// down), that points to a real hit-finding efficiency loss at those eta values -- e.g.
// multiple scattering or energy loss in whatever beampipe material the track crossed at
// that eta (STAR's beampipe isn't one uniform material along z), knocking real,
// geometrically-available hits out of the fit. That's a distinct failure mode from
// acceptance, and it matters more than usual right now because CutClass's fitMaxRatio
// cut is currently enforced via integer division (see SetCutClass.C / PicoBinner.cxx
// comments on the 2026-07-11 revert), making it a de facto "reject unless the ratio is
// exactly 1.0" requirement -- so even a small efficiency dip at large |eta| gets
// punished much harder than the configured 0.5001 threshold would suggest.
//
// Note the histogram itself (eta_fitMaxRatio, booked in PicoBinner.cxx) is filled with
// the TRUE double division regardless of how the cut is currently enforced -- it needs
// to be graduated to be useful as a diagnostic, so don't expect it to just show 0s and
// 1s even while the cut itself is effectively binary.
//
// Reads (exact name from PicoBinner.cxx, written at top level like eta_pT/eta_nHitsFit):
//   eta_fitMaxRatio -- TH2D of eta vs N_hits^fit/N_hits^max.
// Missing (e.g. an older file that predates this histogram being added) just aborts
// with a clear message -- see FindHistoInFile.C to check what's actually in your file.
//
// Usage:
//   root -l -q 'PresentFitMaxRatioEta.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentFitMaxRatioEta.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentFitMaxRatioEta(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* etaFitMaxRatio = (TH2*) inFile->Get("eta_fitMaxRatio");
  if(!etaFitMaxRatio){
    cout << "!!! eta_fitMaxRatio not found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(this file may predate eta_fitMaxRatio being added to PicoBinner.cxx)."
         << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentFitMaxRatioEta","N_{hits}^{fit}/N_{hits}^{max} vs #eta",900,700);
  gPad->SetLogz();
  gPad->SetRightMargin(0.15);

  etaFitMaxRatio->SetTitle("Hit-Finding Efficiency vs #eta");
  etaFitMaxRatio->GetXaxis()->SetTitle("#eta");
  etaFitMaxRatio->GetYaxis()->SetTitle("N_{hits}^{fit}/N_{hits}^{max}");
  etaFitMaxRatio->Draw("COLZ");

  string outPath = a_outDir + "/PresentFitMaxRatioEta.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
