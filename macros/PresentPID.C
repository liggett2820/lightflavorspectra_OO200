// PresentPID.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready PID plots from a PicoBinner yield file --
// dE/dx vs momentum (TPC) and mass vs momentum (BTOF), the two classic "does the
// detector separate species" plots. These are the histograms most people expect to
// see first when you say "here's the PID."
//
// Reads (exact names from PicoBinner.cxx -- see FindHistoInFile.C if your file uses
// different ones, e.g. because a_makeBasicHistos was false):
//   dEdxVsMom_Plus, dEdxVsMom_Minus     -- TPC dE/dx vs p_tot, log-binned in p
//   massVsMom_bTOF_Plus, massVsMom_bTOF_Minus -- BTOF mass vs p_tot
//
// Usage:
//   root -l -q 'PresentPID.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentPID.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include <iostream>
using namespace std;

// Draws `text` at data coordinates (x,y) on the CURRENT pad, bold black,
// with the given TLatex alignment code (22 = centered, 32 = right-aligned/
// vertically centered, i.e. text extends to the LEFT of (x,y)).
// (An earlier version of this drew a white pixel halo behind the text to
// guard against low contrast, but the fixed-pixel offset produced a
// smeared/blobby double-image on the thin "#pi" glyph -- especially bad on
// the log-scaled dE/dx pads -- and in practice these labels land mostly over
// sparse/white background anyway, so plain bold text reads cleanly without
// the extra machinery.)
static void DrawOutlinedLabel(double x, double y, const char* text, double textSize = 0.07, int align = 22){
  TLatex lat;
  lat.SetTextFont(62); // bold Helvetica
  lat.SetTextSize(textSize);
  lat.SetTextAlign(align);
  lat.SetTextColor(kBlack);
  lat.DrawLatex(x, y, text);
}

// Converts a data-space x on the CURRENT pad to the x that sits `pixels` to
// its left on screen -- works uniformly whether the pad's x-axis is linear
// or log-scaled, since it round-trips through actual screen pixels rather
// than subtracting a fixed data-space delta (which would be meaningless on
// a log axis).
static double ShiftXLeftPixels(double x, int pixels){
  return gPad->PixeltoX(gPad->XtoPixel(x) - pixels);
}

// Species labels for the dE/dx-vs-p pads (TPC, log-log). `a_isPlus` selects
// which pad's tracks these label: the Plus pad holds pi+/K+/p, the Minus pad
// holds pi-/K-/pbar. Proton uses the pbar bar-notation (no charge sign) on
// the Minus pad since "p-" isn't how antiprotons are written; pi/K use an
// explicit +/- superscript since both charge states are ordinary particles.
//
// Placement is physics-ordered, not guessed: at a FIXED low momentum, heavier
// species sit at HIGHER dE/dx (proton > kaon > pion) because dE/dx ~ 1/beta^2
// and a heavier particle is less relativistic at the same p. Equivalently, at
// a fixed dE/dx, the momentum where each species crosses that value scales
// with mass (p = beta*gamma*m at fixed beta*gamma), so pion crosses first
// (lowest p), then kaon, then proton (highest p). That ordering is what was
// wrong before (labels sat on the wrong bands / at a common height that only
// the heaviest species actually reached) -- fixing the ordering is the fix.
//
// The three (p, dE/dx) anchors below were computed from a Bethe-Bloch curve
// per species (PDG masses: pi=0.13957, K=0.49368, p=0.93827 GeV/c^2) with the
// overall normalization set so the minimum-ionizing plateau lands at ~2.7
// keV/cm, the typical STAR TPC truncated-mean MIP value. Each anchor sits
// right on its own band in the region where all three are still visibly
// separated from the merged MIP blob (below p ~ 0.3 GeV/c); the label itself
// is then right-aligned and nudged left by a fixed pixel offset so it reads
// as sitting NEXT TO the curve rather than printed on top of it. Actual
// calibration (gas gain, truncation %) varies run to run, so nudge these
// (x,y) pairs by eye after running against real data if a label drifts off
// its band -- what must NOT change is the relative order (pi lowest, K
// middle, p highest).
static void DrawDedxPIDLabels(bool a_isPlus){
  const int leftPx = 14; // gap between label and curve, in screen pixels
  const int protonLeftPx = 7; // proton label nudged closer to its curve (moved right per user feedback -- was sitting too far left at leftPx=14)
  const int align  = 32; // right-aligned, vertically centered
  DrawOutlinedLabel(ShiftXLeftPixels(0.082, leftPx), 7.0  * 1.15, a_isPlus ? "#pi^{+}" : "#pi^{-}", 0.07, align);
  DrawOutlinedLabel(ShiftXLeftPixels(0.149, leftPx), 18.0 * 1.15, a_isPlus ? "K^{+}"   : "K^{-}",   0.07, align);
  DrawOutlinedLabel(ShiftXLeftPixels(0.164, protonLeftPx), 42.0 * 1.15, a_isPlus ? "p"       : "#bar{p}", 0.07, align);
}

// Species labels for the BTOF mass-vs-p pads (linear x, linear y). Same
// charge convention as DrawDedxPIDLabels above. Unlike dE/dx, a BTOF mass
// measurement is calibration-independent by construction (m/q is
// reconstructed directly from beta and p), so these use the exact PDG
// masses rather than data-driven positions -- no run-to-run tuning needed.
// Placed at p=4 GeV/c, comfortably inside the 0-8 GeV/c axis range and past
// the low-p pileup where all species are still resolved.
static void DrawBtofPIDLabels(bool a_isPlus){
  DrawOutlinedLabel(4.0, 0.93827 + 0.13, a_isPlus ? "p"       : "#bar{p}");
  DrawOutlinedLabel(4.0, 0.49368 + 0.13, a_isPlus ? "K^{+}"   : "K^{-}");
  DrawOutlinedLabel(4.0, 0.13957 + 0.13, a_isPlus ? "#pi^{+}" : "#pi^{-}");
}

void PresentPID(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* dEdxPlus   = (TH2*) inFile->Get("dEdxVsMom_Plus");
  TH2* dEdxMinus  = (TH2*) inFile->Get("dEdxVsMom_Minus");
  TH2* massPlus   = (TH2*) inFile->Get("massVsMom_bTOF_Plus");
  TH2* massMinus  = (TH2*) inFile->Get("massVsMom_bTOF_Minus");

  if(!dEdxPlus && !massPlus){
    cout << "!!! Neither dEdxVsMom_Plus nor massVsMom_bTOF_Plus found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced)." << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentPID","PID Overview",1400,1000);
  c->Divide(2,2);

  // Units: StPicoTrack::dEdx() is stored in keV/cm (StPicoTrack.cxx::setDedx() scales
  // the native GeV/cm value by 1e6); mass comes from PhysMath::mass_fromBeta(p,beta)
  // with p in GeV/c, i.e. natural units (c=1), so mass is in GeV/c^2.
  // Extra right margin on every COLZ pad below -- otherwise the z-axis (palette)
  // power-of-10 labels get clipped/overlap the next pad in a Divide() grid.
  c->cd(1); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(dEdxPlus){ dEdxPlus->SetTitle("dE/dx vs p (Positive Tracks)"); dEdxPlus->GetYaxis()->SetTitle("dE/dx (keV/cm)"); dEdxPlus->Draw("COLZ"); DrawDedxPIDLabels(true); }

  c->cd(2); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(dEdxMinus){ dEdxMinus->SetTitle("dE/dx vs p (Negative Tracks)"); dEdxMinus->GetYaxis()->SetTitle("dE/dx (keV/cm)"); dEdxMinus->Draw("COLZ"); DrawDedxPIDLabels(false); }

  c->cd(3); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(massPlus){ massPlus->SetTitle("BTOF Mass vs p (Positive Tracks)"); massPlus->GetYaxis()->SetTitle("m/q (GeV/c^{2})"); massPlus->Draw("COLZ"); DrawBtofPIDLabels(true); }

  c->cd(4); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(massMinus){ massMinus->SetTitle("BTOF Mass vs p (Negative Tracks)"); massMinus->GetYaxis()->SetTitle("m/q (GeV/c^{2})"); massMinus->Draw("COLZ"); DrawBtofPIDLabels(false); }

  string outPath = a_outDir + "/PresentPID.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
