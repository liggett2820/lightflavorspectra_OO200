// PresentAcceptance.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready acceptance plots from a PicoBinner yield file --
// where in (y,pT) and (eta,pT) space the detector actually has coverage. Useful for
// explaining why some rapidity bins in the fitted spectra are empty (see the
// rapidity-bin-walker bug this repo's ZFitter had -- InspectRapBins.C/
// InspectRawYieldRapBins.C are the companion diagnostics for that).
//
// Reads (exact names from PicoBinner.cxx):
//   pion_y_pT  -- y vs pT assuming the pion mass hypothesis, any track passing basic
//                 cuts (not gated on TOF -- this is closest to "TPC-only acceptance")
//   eta_pT     -- eta vs pT, same idea but detector-frame (no mass hypothesis needed)
//
// Usage:
//   root -l -q 'PresentAcceptance.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentAcceptance.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentAcceptance(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* pionYPt      = (TH2*) inFile->Get("pion_y_pT");
  TH2* etaPt        = (TH2*) inFile->Get("eta_pT");

  if(!pionYPt && !etaPt){
    cout << "!!! Neither pion_y_pT nor eta_pT found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file." << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentAcceptance","Acceptance overview",1000,500);
  c->Divide(2,1);

  // y/eta are dimensionless (no unit); pT is in GeV/c (filled from PicoBinner's pT
  // variable, same convention as the dE/dx-vs-momentum plots in PresentPID.C).
  // Extra right margin -- otherwise the z-axis (palette) power-of-10 labels get
  // clipped/overlap the next pad in this Divide(2,1) grid.
  // Both (pseudo)rapidity axes are display-clipped to [-2,2] -- just clips the plotted
  // range, doesn't rebin or re-cut the underlying data.
  c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(pionYPt){ pionYPt->SetTitle("Acceptance: y (pion mass assumed) vs pT, TPC-only cuts"); pionYPt->GetYaxis()->SetTitle("p_{T} (GeV/c)"); pionYPt->GetXaxis()->SetRangeUser(-2, 2); pionYPt->Draw("COLZ"); }

  c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(etaPt){ etaPt->SetTitle("Acceptance: #eta vs pT"); etaPt->GetYaxis()->SetTitle("p_{T} (GeV/c)"); etaPt->GetXaxis()->SetRangeUser(-2, 2); etaPt->Draw("COLZ"); }

  string outPath = a_outDir + "/PresentAcceptance.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
