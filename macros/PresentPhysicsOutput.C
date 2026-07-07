// PresentPhysicsOutput.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready view of the actual raw data ZFitter reads and
// fits bin-by-bin -- not as visually striking as the PID plots, but it's the honest
// "here is the raw yield my analysis starts from" plot, and it's exactly the histogram
// whose sparse/populated-bin structure this repo spent a long debugging session on
// (see the rapidity-bin-walker AND/OR fix in the ZFitter_*.cxx source files, and
// InspectRawYieldRapBins.C).
//
// Reads (exact Form() string from PicoBinner.cxx): zTPCPlus_<Species>_Cent<N> /
// zTPCMinus_<Species>_Cent<N> -- a TH3 with axes (rapidity, mT-m0, Z_TPC), <N> NOT
// zero-padded (Cent0, Cent1, ... not Cent00, Cent01).
//
// Projects down to (rapidity, mT-m0), summing over the Z_TPC axis -- this is what
// InspectRawYieldRapBins.C already inspects numerically; this macro just draws it.
//
// CAVEAT: ROOT's TH3::Project3D(option) axis-order convention (which letter becomes
// the new 2D histogram's x-axis vs y-axis) is easy to get backwards without testing
// live against ROOT -- if rapidity/mT-m0 come out swapped on the axes, flip "xy" to
// "yx" below.
//
// Usage:
//   root -l -q 'PresentPhysicsOutput.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentPhysicsOutput.C("/path/to/file.root","Proton",0)'
//   root -l -q 'PresentPhysicsOutput.C("/path/to/file.root","Proton",0,"/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentPhysicsOutput(string a_inputFile, string a_species = "Proton", int a_centIndex = 0, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  string plusName  = Form("zTPCPlus_%s_Cent%d",a_species.c_str(),a_centIndex);
  string minusName = Form("zTPCMinus_%s_Cent%d",a_species.c_str(),a_centIndex);

  TH3* zTPCPlus  = (TH3*) inFile->Get(plusName.c_str());
  TH3* zTPCMinus = (TH3*) inFile->Get(minusName.c_str());

  if(!zTPCPlus && !zTPCMinus){
    cout << "!!! Neither " << plusName << " nor " << minusName << " found in "
         << a_inputFile << " -- run FindHistoInFile.C to see what's actually in this "
         << "file (species name/casing and centrality index must match exactly)." << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentPhysicsOutput","Raw yield overview",1400,600);
  c->Divide(2,1);

  c->cd(1); gPad->SetLogz();
  if(zTPCPlus){
    TH2* projPlus = (TH2*) zTPCPlus->Project3D("xy");
    projPlus->SetTitle(Form("Raw TPC yield: %s+, Cent %d (y vs m_{T}-m_{0}, summed over Z_{TPC})",a_species.c_str(),a_centIndex));
    projPlus->Draw("COLZ");
  }

  c->cd(2); gPad->SetLogz();
  if(zTPCMinus){
    TH2* projMinus = (TH2*) zTPCMinus->Project3D("xy");
    projMinus->SetTitle(Form("Raw TPC yield: %s-, Cent %d (y vs m_{T}-m_{0}, summed over Z_{TPC})",a_species.c_str(),a_centIndex));
    projMinus->Draw("COLZ");
  }

  string outPath = a_outDir + "/PresentPhysicsOutput.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
