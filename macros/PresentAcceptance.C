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
//   pion_y_pT_tof, kaon_y_pT_tof, proton_y_pT_tof -- same as pion_y_pT but additionally
//                 requiring a BTOF or ETOF match, showing how acceptance shrinks once
//                 TOF is required
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
  TH2* protonYPtTof = (TH2*) inFile->Get("proton_y_pT_tof");

  if(!pionYPt && !etaPt){
    cout << "!!! Neither pion_y_pT nor eta_pT found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file." << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentAcceptance","Acceptance overview",1400,500);
  c->Divide(3,1);

  // y/eta are dimensionless (no unit); pT is in GeV/c (filled from PicoBinner's pT
  // variable, same convention as the dE/dx-vs-momentum plots in PresentPID.C).
  // Extra right margin -- otherwise the z-axis (palette) power-of-10 labels get
  // clipped/overlap the next pad in this Divide(3,1) grid.
  // Rapidity (y) axis capped at 0 on the high end for display -- keeps whatever low
  // edge the histogram was actually booked with (species-dependent, from CutClass),
  // just clips the plotted range rather than re-binning/re-cutting the data itself.
  // eta_pT is left alone -- its x-axis is detector eta, not rapidity.
  c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(pionYPt){ pionYPt->SetTitle("Acceptance: y (pion mass assumed) vs pT, TPC-only cuts"); pionYPt->GetYaxis()->SetTitle("p_{T} (GeV/c)"); pionYPt->GetXaxis()->SetRangeUser(pionYPt->GetXaxis()->GetXmin(), 0); pionYPt->Draw("COLZ"); }

  c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(etaPt){ etaPt->SetTitle("Acceptance: #eta vs pT"); etaPt->GetYaxis()->SetTitle("p_{T} (GeV/c)"); etaPt->Draw("COLZ"); }

  c->cd(3); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(protonYPtTof){ protonYPtTof->SetTitle("Acceptance: y (proton mass assumed) vs pT, TOF-matched"); protonYPtTof->GetYaxis()->SetTitle("p_{T} (GeV/c)"); protonYPtTof->GetXaxis()->SetRangeUser(protonYPtTof->GetXaxis()->GetXmin(), 0); protonYPtTof->Draw("COLZ"); }

  string outPath = a_outDir + "/PresentAcceptance.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
