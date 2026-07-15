// PresentNHitsMaxBySector.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: requested at the 2026-07-15 advisor meeting, following up on the
// eta~1.1-1.2 hit-loss bump seen in PresentHitLossEta.C / the third panel of
// PresentNHitsEta.C. That bump's shape (rises then falls, rather than monotonic with
// |eta|) already ruled out the beampipe-material hypothesis (STAR's beryllium section
// covers this analysis's full eta range), leaving a TPC geometric/mechanical
// transition (e.g. barrel-to-endcap exit crossover, sector boundaries) as the leading
// candidate. This macro is the direct test of the "sector boundaries" half of that:
// it plots N_hits^max vs eta separately for each of the 12 azimuthal (30-degree) TPC
// sectors, so you can see by eye whether the bump (or any other feature) is uniform
// across all 12 sectors (-> a purely eta-dependent geometric effect, e.g. the
// barrel/endcap crossover) or concentrated in a subset of sectors (-> a sector-
// specific hardware/alignment issue instead).
//
// Sectors are defined purely by phi (12 wedges of pi/6 each, phi in [-pi,pi)) -- NOT
// split further by west/east TPC end, since eta's sign already distinguishes that.
// Each sector's plot spans the full eta range, same binning as the combined
// eta_nHitsMax histogram, so they're directly comparable to each other and to the
// all-sectors-combined plot in PresentNHitsEta.C.
//
// Reads (exact names from PicoBinner.cxx):
//   eta_nHitsMax_sector0 .. eta_nHitsMax_sector11 -- TH2I of eta vs N_hits^max, one per
//   phi sector, written at top level (same placement as eta_nHitsMax itself). These
//   histograms do NOT exist in any yield file produced before this macro's companion
//   PicoBinner.cxx change (2026-07-15) -- rerun PicoBinner to populate them.
//
// Usage:
//   root -l -q 'PresentNHitsMaxBySector.C("/path/to/yieldHistos_OO200_pion.root")'
//   root -l -q 'PresentNHitsMaxBySector.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentNHitsMaxBySector(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* sectorHist[12] = {NULL};
  int  nFound = 0;
  int  globalMax = 0; // track the tallest bin across all sectors so every pad shares one z-scale

  for(int sectorIndex = 0; sectorIndex < 12; sectorIndex++){
    sectorHist[sectorIndex] = (TH2*) inFile->Get(Form("eta_nHitsMax_sector%d",sectorIndex));
    if(sectorHist[sectorIndex]){
      nFound++;
      int thisMax = (int) sectorHist[sectorIndex]->GetMaximum();
      if(thisMax > globalMax) globalMax = thisMax;
    }
  }

  if(nFound == 0){
    cout << "!!! None of eta_nHitsMax_sector0..11 found in " << a_inputFile << " -- "
         << "this file predates the 2026-07-15 PicoBinner.cxx change that adds the "
         << "per-sector histograms. Rerun PicoBinner against this species first, or "
         << "run FindHistoInFile.C to double check what's actually in this file."
         << endl;
    return;
  }
  if(nFound < 12){
    cout << "!!! Only " << nFound << "/12 eta_nHitsMax_sector* histograms found -- "
         << "missing sectors will be left blank below." << endl;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentNHitsMaxBySector",
                            "N_{hits}^{max} vs #eta, by TPC sector (30#circ each)",2400,1800);
  c->Divide(4,3); // 12 pads, one per sector, in reading order (sector 0 top-left)

  for(int sectorIndex = 0; sectorIndex < 12; sectorIndex++){
    c->cd(sectorIndex+1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(sectorHist[sectorIndex]){
      sectorHist[sectorIndex]->SetTitle(Form("Sector %d",sectorIndex));
      sectorHist[sectorIndex]->GetXaxis()->SetTitle("#eta");
      sectorHist[sectorIndex]->GetYaxis()->SetTitle("N_{hits}^{max}");
      // Shared z-scale across all 12 pads (see globalMax above) so sector-to-sector
      // differences show up as color differences, not just each pad auto-scaling to
      // its own max and hiding a real deficit.
      sectorHist[sectorIndex]->SetMaximum(globalMax);
      sectorHist[sectorIndex]->SetMinimum(1); // keep log-z well-defined; empty bins stay white
      sectorHist[sectorIndex]->Draw("COLZ");
    }
  }

  string outPath = a_outDir + "/PresentNHitsMaxBySector.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
