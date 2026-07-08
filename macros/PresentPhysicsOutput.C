// PresentPhysicsOutput.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready view of the actual raw data ZFitter reads and
// fits bin-by-bin -- not as visually striking as the PID plots, but it's the honest
// "here is the raw yield my analysis starts from" plot, and it's exactly the histogram
// whose sparse/populated-bin structure this repo spent a long debugging session on
// (see the rapidity-bin-walker AND/OR fix in the ZFitter_*.cxx source files, and
// InspectRawYieldRapBins.C).
//
// Reads (exact Form() string from PicoBinner.cxx): yields/<Species>_space/
// zTPCPlus_<Species>_Cent<N> / zTPCMinus_<Species>_Cent<N> -- a TH3 with axes
// (rapidity, mT-m0, Z_TPC), <N> NOT zero-padded (Cent0, Cent1, ... not Cent00, Cent01).
// Confirmed via FindHistoInFile.C against a real production yield file that these
// live inside the "yields/<Species>_space" subdirectory PicoBinner.cxx writes them
// into (outFile->mkdir("yields/" + species + "_space"); outFile->cd(...)) -- NOT at
// top level, same class of directory-prefix miss as PresentCalibrationQA.C's calib_*
// lookup.
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

  string speciesDir = "yields/" + a_species + "_space/";
  string plusName  = speciesDir + Form("zTPCPlus_%s_Cent%d",a_species.c_str(),a_centIndex);
  string minusName = speciesDir + Form("zTPCMinus_%s_Cent%d",a_species.c_str(),a_centIndex);

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

  // Extra right margin on both COLZ pads -- otherwise the z-axis (palette) power-of-10
  // labels get clipped/overlap the next pad in this Divide(2,1) grid.
  // Rapidity (y, the projected x-axis) capped at 0 on the high end for display -- keeps
  // whatever low edge the source histogram was booked with, just clips the plotted
  // range rather than re-binning/re-cutting the underlying data.
  c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(zTPCPlus){
    TH2* projPlus = (TH2*) zTPCPlus->Project3D("xy");
    projPlus->SetTitle(Form("Raw TPC yield: %s+, Cent %d (y vs m_{T}-m_{0}, summed over Z_{TPC})",a_species.c_str(),a_centIndex));
    projPlus->GetXaxis()->SetRangeUser(projPlus->GetXaxis()->GetXmin(), 0);
    projPlus->Draw("COLZ");
  }

  c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(zTPCMinus){
    TH2* projMinus = (TH2*) zTPCMinus->Project3D("xy");
    projMinus->SetTitle(Form("Raw TPC yield: %s-, Cent %d (y vs m_{T}-m_{0}, summed over Z_{TPC})",a_species.c_str(),a_centIndex));
    projMinus->GetXaxis()->SetRangeUser(projMinus->GetXaxis()->GetXmin(), 0);
    projMinus->Draw("COLZ");
  }

  string outPath = a_outDir + "/PresentPhysicsOutput.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
