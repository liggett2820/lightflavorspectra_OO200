// InspectRapBins.C
//
// Quick diagnostic: prints which rapidity bins of a dEdxSpectra_*/InvBetaBTOFSpectra_*
// 2D histogram (rapidity x mT-m0) actually carry nonzero yield, and the real rapidity
// value at each populated bin's center. Written as a proper macro (not a `root -e`
// one-liner) because opening these files via a bare `-e` multi-statement string was
// segfaulting in Cling on this machine, while the exact same TFile::Open/Get calls
// inside CompareSpectra_OO200_vs_etof.C work fine.
//
// Usage:
//   root -l -q 'InspectRapBins.C("/Users/aliggett/data/OO/etof_compare/spectra_200GeV_2024_09_12_000_Proton_TPC_RapBin_-1_-1_-1_-1_-1_-1.root","DeDx_FitData/ProtonPlus","dEdxSpectra_ProtonPlus_Cent00")'

#include "TFile.h"
#include "TDirectory.h"
#include "TH2D.h"
#include <iostream>
#include <string>
using namespace std;

void InspectRapBins(string a_file, string a_dirPath, string a_histoName){
  TFile* f = TFile::Open(a_file.c_str());
  if(!f || f->IsZombie()){ cout << "Could not open " << a_file << endl; return; }

  TDirectory* d = (TDirectory*) f->Get(a_dirPath.c_str());
  if(!d){ cout << "Directory not found: " << a_dirPath << endl; return; }

  TH2D* h = (TH2D*) d->Get(a_histoName.c_str());
  if(!h){ cout << "Histogram not found: " << a_dirPath << "/" << a_histoName << endl; return; }

  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();
  cout << a_file << "  " << a_dirPath << "/" << a_histoName << endl;
  cout << "NBinsX=" << nx << "  (rapidity range " << h->GetXaxis()->GetXmin() << " to " << h->GetXaxis()->GetXmax() << ")" << endl;
  cout << "NBinsY=" << ny << "  (mT-m0 range " << h->GetYaxis()->GetXmin() << " to " << h->GetYaxis()->GetXmax() << ")" << endl;

  int nPopulated = 0;
  for(int i = 1; i <= nx; i++){
    double s = 0;
    for(int j = 1; j <= ny; j++) s += h->GetBinContent(i,j);
    if(s != 0){
      nPopulated++;
      cout << "  rapBin " << i << "  y=" << h->GetXaxis()->GetBinCenter(i)
           << "  [" << h->GetXaxis()->GetBinLowEdge(i) << "," << h->GetXaxis()->GetBinUpEdge(i) << "]"
           << "  sum=" << s << endl;
    }
  }
  cout << nPopulated << " of " << nx << " rapidity bins populated." << endl;
}
