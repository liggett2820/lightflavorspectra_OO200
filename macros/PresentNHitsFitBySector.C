// PresentNHitsFitBySector.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: companion to PresentNHitsMaxBySector.C. That macro tests whether
// the eta~1.1-1.2 hit-loss bump traces back to the geometric ceiling (N_hits^max)
// varying by TPC sector. This one tests the other half: whether N_hits^fit itself --
// how many of the possible hits actually get fit -- varies by sector at a given eta,
// independent of the ceiling. A sector with a dead or noisy padrow region, for
// example, would show a normal nHitsMax ceiling (geometry is fine) but a depressed
// nHitsFit relative to other sectors at the same eta (something about that sector's
// hit-finding is worse) -- a pattern PresentNHitsMaxBySector.C alone can't distinguish
// from a purely geometric effect.
//
// Sectors are defined purely by phi (12 wedges of pi/6 each, phi in [-pi,pi)) -- NOT
// split further by west/east TPC end, since eta's sign already distinguishes that.
// Same convention as PresentNHitsMaxBySector.C/eta_nHitsMax_bySector.
//
// Reads (exact names from PicoBinner.cxx, added 2026-07-16):
//   eta_nHitsFit_sector0 .. eta_nHitsFit_sector11 -- TH2I of eta vs N_hits^fit, one per
//   phi sector, written at top level (same placement as eta_nHitsFit/eta_nHitsMax
//   themselves). These histograms do NOT exist in any yield file produced before this
//   macro's companion PicoBinner.cxx change (2026-07-16) -- rerun PicoBinner to
//   populate them, same as was needed for the sector0..11 eta_nHitsMax histograms.
//
// Usage:
//   root -l -q 'PresentNHitsFitBySector.C("/path/to/yieldHistos_OO200_pion.root")'
//   root -l -q 'PresentNHitsFitBySector.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentNHitsFitBySector(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* sectorHist[12] = {NULL};
  int  nFound = 0;
  int  globalMax = 0; // track the tallest bin across all sectors so every pad shares one z-scale

  for(int sectorIndex = 0; sectorIndex < 12; sectorIndex++){
    sectorHist[sectorIndex] = (TH2*) inFile->Get(Form("eta_nHitsFit_sector%d",sectorIndex));
    if(sectorHist[sectorIndex]){
      nFound++;
      int thisMax = (int) sectorHist[sectorIndex]->GetMaximum();
      if(thisMax > globalMax) globalMax = thisMax;
    }
  }

  if(nFound == 0){
    cout << "!!! None of eta_nHitsFit_sector0..11 found in " << a_inputFile << " -- "
         << "this file predates the 2026-07-16 PicoBinner.cxx change that adds the "
         << "per-sector nHitsFit histograms. Rerun PicoBinner against this species first, or "
         << "run FindHistoInFile.C to double check what's actually in this file."
         << endl;
    return;
  }
  if(nFound < 12){
    cout << "!!! Only " << nFound << "/12 eta_nHitsFit_sector* histograms found -- "
         << "missing sectors will be left blank below." << endl;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentNHitsFitBySector",
                            "N_{hits}^{fit} vs #eta, by TPC sector (30#circ each)",2400,1800);
  c->Divide(4,3); // 12 pads, one per sector, in reading order (sector 0 top-left)

  for(int sectorIndex = 0; sectorIndex < 12; sectorIndex++){
    c->cd(sectorIndex+1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(sectorHist[sectorIndex]){
      sectorHist[sectorIndex]->SetTitle(Form("Sector %d",sectorIndex));
      sectorHist[sectorIndex]->GetXaxis()->SetTitle("#eta");
      sectorHist[sectorIndex]->GetYaxis()->SetTitle("N_{hits}^{fit}");
      // Shared z-scale across all 12 pads (see globalMax above) so sector-to-sector
      // differences show up as color differences, not just each pad auto-scaling to
      // its own max and hiding a real deficit.
      sectorHist[sectorIndex]->SetMaximum(globalMax);
      sectorHist[sectorIndex]->SetMinimum(1); // keep log-z well-defined; empty bins stay white
      sectorHist[sectorIndex]->Draw("COLZ");
    }
  }

  string outPath = a_outDir + "/PresentNHitsFitBySector.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
