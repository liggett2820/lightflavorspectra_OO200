// InspectRawYieldRapBins.C
//
// Companion to InspectRapBins.C, but for PicoBinner's raw yield output instead of
// ZFitter's fitted spectra output. PicoBinner writes the raw pre-fit yield as a TH3I
// per (species, charge, centrality) named zTPCPlus_<Species>_Cent<N> /
// zTPCMinus_<Species>_Cent<N> (source/PicoBinner.cxx, ~line 795), axes (rapidity,
// mT-m0, Z_TPC) -- this is what ZFitter actually reads and fits bin-by-bin. Prints
// which rapidity (x) bins have any nonzero raw content, summing over both the
// mT-m0 and Z_TPC axes, so it's directly comparable to InspectRapBins.C's output on
// the fitted dEdxSpectra_* histogram (same rapidity axis/binning, since CutClass's
// getRapMtM0BinStructure defines both).
//
// PicoBinner.cxx writes these under "yields/<Species>_space" (e.g. "yields/Proton_space"
// for a_partIndex's neutral-charge name, source/PicoBinner.cxx ~line 1930), not at the
// file's top level -- pass that as a_dirPath (default already set to the Proton case).
//
// Usage:
//   root -l -q 'InspectRawYieldRapBins.C("/Users/aliggett/data/OO/yieldHistos_OOGeV_proton_2026_07_01.root","zTPCPlus_Proton_Cent00")'

#include "TFile.h"
#include "TDirectory.h"
#include "TH3.h"
#include <iostream>
#include <string>
using namespace std;

void InspectRawYieldRapBins(string a_file, string a_histoName, string a_dirPath = "yields/Proton_space"){
  TFile* f = TFile::Open(a_file.c_str());
  if(!f || f->IsZombie()){ cout << "Could not open " << a_file << endl; return; }

  TH3I* h = (TH3I*) f->Get(a_histoName.c_str());
  if(!h){
    TDirectory* d = (TDirectory*) f->Get(a_dirPath.c_str());
    if(!d){
      cout << "Not found at top level, and directory not found either: " << a_dirPath << endl;
      cout << "Try `root -l " << a_file << "` then `.ls` / TBrowser to find the real path." << endl;
      return;
    }
    h = (TH3I*) d->Get(a_histoName.c_str());
    if(!h){
      cout << "Histogram not found: " << a_dirPath << "/" << a_histoName << endl;
      return;
    }
  }

  int nx = h->GetNbinsX(); // rapidity
  int ny = h->GetNbinsY(); // mT-m0
  int nz = h->GetNbinsZ(); // Z_TPC

  cout << a_file << "  " << a_histoName << endl;
  cout << "NBinsX (rapidity)=" << nx << "  (range " << h->GetXaxis()->GetXmin() << " to " << h->GetXaxis()->GetXmax() << ")" << endl;
  cout << "NBinsY (mT-m0)=" << ny << "  NBinsZ (Z_TPC)=" << nz << endl;

  int nPopulated = 0;
  double totalIntegral = h->Integral(); // Integral() only sums the *regular* (non-under/overflow) bins
  cout << "Total integral (regular bins only): " << totalIntegral << endl;

  for(int i = 1; i <= nx; i++){
    double s = 0;
    for(int j = 1; j <= ny; j++){
      for(int k = 1; k <= nz; k++){
        s += h->GetBinContent(i,j,k);
      }
    }
    if(s != 0){
      nPopulated++;
      cout << "  rapBin " << i << "  y=" << h->GetXaxis()->GetBinCenter(i)
           << "  [" << h->GetXaxis()->GetBinLowEdge(i) << "," << h->GetXaxis()->GetBinUpEdge(i) << "]"
           << "  sum=" << s << endl;
    }
  }
  cout << nPopulated << " of " << nx << " rapidity bins populated." << endl;

  // KEY CHECK: if the real rapidity distribution doesn't actually overlap with this
  // histogram's configured axis range [xmin,xmax], ROOT silently dumps those entries
  // into the rapidity (x) under/overflow bins instead of any "regular" bin -- they
  // never show up in the per-bin printout above or in Integral(), so a small/narrow
  // "populated" bin count above does NOT by itself mean the real distribution is
  // actually that narrow. Sum the x-underflow (bin index 0) and x-overflow (bin index
  // nx+1) slices across all mT-m0/Z_TPC bins to check.
  double underflowSum = 0, overflowSum = 0;
  for(int j = 0; j <= ny+1; j++){
    for(int k = 0; k <= nz+1; k++){
      underflowSum += h->GetBinContent(0,j,k);
      overflowSum  += h->GetBinContent(nx+1,j,k);
    }
  }
  cout << "Rapidity-axis UNDERFLOW sum (y < " << h->GetXaxis()->GetXmin() << "): " << underflowSum << endl;
  cout << "Rapidity-axis OVERFLOW sum  (y > " << h->GetXaxis()->GetXmax() << "): " << overflowSum << endl;
  if(overflowSum > 3*totalIntegral || underflowSum > 3*totalIntegral){
    cout << "  ==> WARNING: under/overflow dwarfs the regular-bin integral. The real rapidity"
         << " distribution likely extends well outside this histogram's configured axis"
         << " range [" << h->GetXaxis()->GetXmin() << "," << h->GetXaxis()->GetXmax()
         << "] -- check SetCutClass.C's setVariableRapMtM0BinningInfo() range for this species/detector." << endl;
  }
}
