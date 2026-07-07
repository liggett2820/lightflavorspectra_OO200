// PresentCalibrationQA.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready calibration QA plots from a PicoBinner yield
// file -- do the Bichsel-curve dE/dx predictions actually line up with the data, and
// does the BTOF-only-selected "isolated" dE/dx band land where the TPC dE/dx cut
// expects it to.
//
// Reads (exact names/Form() strings from PicoBinner.cxx -- only present if
// _PICOBINNER_VZ_CALIBRATION_ was defined in makefile_toggles.h when this file was
// produced, which it is by default in this repo):
//   calib_<Species>_VzIndex<NN> -- TH3 (eta, beta*gamma from p/m, ln(dE/dx)-prediction)
//                 for one V_Z calibration slice; NN defaults to 00 (the first V_Z bin,
//                 always present regardless of how many V_Z bins were configured)
//   dEdx_<Species>Plus_Isolated -- dE/dx vs p for tracks BTOF-identified as <Species>,
//                 with no other species contaminating that BTOF acceptance window
//
// CAVEAT: ROOT's TH3::Project3D(option) axis-order convention (which letter becomes
// the new 2D histogram's x-axis vs y-axis) is easy to get backwards without testing
// live against ROOT -- if beta*gamma/ln(dE/dx ratio) come out swapped on the axes,
// flip "zy" to "yz" below.
//
// Usage:
//   root -l -q 'PresentCalibrationQA.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentCalibrationQA.C("/path/to/file.root","Proton")'
//   root -l -q 'PresentCalibrationQA.C("/path/to/file.root","Proton","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentCalibrationQA(string a_inputFile, string a_species = "Proton", string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  string calibName = "calib_" + a_species + "_VzIndex00";
  string isolatedName = "dEdx_" + a_species + "Plus_Isolated";

  TH3* calib3D = (TH3*) inFile->Get(calibName.c_str());
  TH2* isolated = (TH2*) inFile->Get(isolatedName.c_str());

  if(!calib3D && !isolated){
    cout << "!!! Neither " << calibName << " nor " << isolatedName << " found in "
         << a_inputFile << " -- run FindHistoInFile.C to see what's actually in this "
         << "file (the calibration histos only exist if _PICOBINNER_VZ_CALIBRATION_ "
         << "was on, and the species name/casing must match exactly)." << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentCalibrationQA","Calibration QA overview",1400,600);
  c->Divide(2,1);

  c->cd(1); gPad->SetLogz();
  if(calib3D){
    // Project down to (beta*gamma, ln(dE/dx)-prediction), summing over eta -- this is
    // the actual calibration check: a flat band centered at 0 means the Bichsel
    // prediction matches the data at this V_Z slice.
    TH2* calibProj = (TH2*) calib3D->Project3D("zy");
    calibProj->SetTitle(Form("%s calibration (Vz bin 0): ln(dE/dx / prediction) vs #beta#gamma",a_species.c_str()));
    calibProj->Draw("COLZ");
  }

  c->cd(2); gPad->SetLogx(); gPad->SetLogz();
  if(isolated){
    isolated->SetTitle(Form("Isolated (BTOF-tagged) dE/dx for %s+",a_species.c_str()));
    isolated->Draw("COLZ");
  }

  string outPath = a_outDir + "/PresentCalibrationQA.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
