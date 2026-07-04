// DumpMTM0ZTPCProj2D.C -- lightflavorspectra_OO200
//
// Dumps the RAW (pre-fit) per-(rapidity bin, mT-m0 bin) dE/dx-amplitude counts that
// TPC fits actually feed on -- the exact same data drawn+fit in each
// dEdxFits_Cent0X/RapBin_XX_mTm0Bin_YY.png image -- out to a plain CSV. This is one
// half of a two-repo comparison: run this once per repo (each with its OWN compiled
// Histo2D library loaded -- etof and OO200 each compile a class named Histo2D
// separately, so you cannot load both repos' libraries in the same ROOT session), then
// diff the two resulting CSVs with DiffMTM0ZTPCDumps.C, which needs no library at all.
//
// The data comes from the Histo2D objects ZFitter writes to:
//   DeDx_FitData/<Species>/Cent/mTm0_ZTPC_Proj2D/mTm0_ZTPC_<Species><Charge>_Cent%02d_rapIndex%03d
// Histo2D is a custom (non-TH1/TH2) class -- one object per (species,charge,cent,
// rapidity bin), holding a TH1D dE/dx-amplitude distribution per mT-m0 bin
// (GetInternalYProjection(mtm0BinIndex)). There's no public getter for its bin count,
// so a_numMtm0Bins must be supplied (matches whatever SetCutClass set for this
// species/detector -- 150 for Proton TPC in both repos as of 2026-07).
//
// Usage (run from lightflavorspectra_OO200/macros/, loading LoadZFitterLibs.C first
// in the SAME invocation, same pattern as running RunZFitter.C itself):
//   root -l -q LoadZFitterLibs.C 'DumpMTM0ZTPCProj2D.C("/Users/aliggett/data/OO/etof_compare/spectra_OO200_Proton_TPC_RapBin_-1_-1_-1_-1_-1_-1.root","Proton","Plus",{0,1,2,3,4},150,"/Users/aliggett/data/OO/etof_compare/dump_oo200_ProtonPlus_TPC.csv")'

#include "../headers/Histo2D.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
using namespace std;

void DumpMTM0ZTPCProj2D(string a_specFile, string a_particleName, string a_charge,
                          vector<int> a_centralities, int a_numMtm0Bins, string a_outCsv){

  TFile* f = TFile::Open(a_specFile.c_str());
  if(!f || f->IsZombie()){ cout << "Could not open " << a_specFile << endl; return; }

  string dirPath = Form("DeDx_FitData/%s/Cent/mTm0_ZTPC_Proj2D", a_particleName.c_str());
  TDirectory* d = (TDirectory*) f->Get(dirPath.c_str());
  if(!d){ cout << "Directory not found: " << dirPath << " (in " << a_specFile << ")" << endl; return; }

  ofstream out(a_outCsv.c_str());
  out << "centrality,rapIndex,mtm0Bin,rapLow,rapHigh,integral\n";

  int totalObjectsFound = 0;
  for(int cent : a_centralities){
    int rapIndex = 0;
    int consecutiveMisses = 0;
    while(consecutiveMisses < 5){
      string histName = Form("mTm0_ZTPC_%s%s_Cent%02d_rapIndex%03d", a_particleName.c_str(), a_charge.c_str(), cent, rapIndex);
      TObject* obj = d->Get(histName.c_str());
      if(!obj){
        consecutiveMisses++;
        rapIndex++;
        continue;
      }
      consecutiveMisses = 0;
      totalObjectsFound++;

      Histo2D* h2d = (Histo2D*) obj;
      double rapLow = 0, rapHigh = 0;
      string title = h2d->GetTitle();
      size_t bracketPos = title.find("y=[");
      if(bracketPos != string::npos){
        sscanf(title.c_str() + bracketPos, "y=[%lf,%lf]", &rapLow, &rapHigh);
      }

      for(int mtm0Bin = 0; mtm0Bin < a_numMtm0Bins; mtm0Bin++){
        TH1D* proj = h2d->GetInternalYProjection(mtm0Bin);
        double integral = proj ? proj->Integral() : 0;
        if(integral != 0){
          out << cent << "," << rapIndex << "," << mtm0Bin << "," << rapLow << "," << rapHigh << "," << integral << "\n";
        }
      }

      rapIndex++;
    }
  }

  out.close();
  cout << "Found " << totalObjectsFound << " rapidity-bin Histo2D objects across "
       << a_centralities.size() << " centralities." << endl;
  cout << "Wrote " << a_outCsv << endl;
}
