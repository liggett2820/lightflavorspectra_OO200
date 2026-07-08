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
//   Vz_Calibration/calib_<Species>_VzIndex<NN> -- TH3 (eta, beta*gamma from p/m,
//                 ln(dE/dx)-prediction) for one V_Z calibration slice, inside the
//                 "Vz_Calibration" subdirectory PicoBinner.cxx writes it into; NN
//                 defaults to 00 (the first V_Z bin, always present regardless of how
//                 many V_Z bins were configured)
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
#include "TLatex.h"
#include <iostream>
using namespace std;

void PresentCalibrationQA(string a_inputFile, string a_species = "Proton", string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  string calibName = "calib_" + a_species + "_VzIndex00";
  // calib_* histograms live inside a "Vz_Calibration" subdirectory in the file
  // (PicoBinner.cxx does outFile->mkdir("Vz_Calibration"); outFile->cd("Vz_Calibration");
  // before writing them), unlike the isolated dE/dx histograms below, which are written
  // at top level -- confirmed via FindHistoInFile.C against a real production yield
  // file (calib_Pion_VzIndex00 showed up at /Vz_Calibration/calib_Pion_VzIndex00, not
  // at top level). Get() needs the directory prefix, or it silently returns null.
  string calibPath = "Vz_Calibration/" + calibName;
  string isolatedName = "dEdx_" + a_species + "Plus_Isolated";

  TH3* calib3D = (TH3*) inFile->Get(calibPath.c_str());
  TH2* isolated = (TH2*) inFile->Get(isolatedName.c_str());

  if(!calib3D && !isolated){
    cout << "!!! Neither " << calibPath << " nor " << isolatedName << " found in "
         << a_inputFile << " -- run FindHistoInFile.C to see what's actually in this "
         << "file (the calibration histos only exist if _PICOBINNER_VZ_CALIBRATION_ "
         << "was on, and the species name/casing must match exactly)." << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentCalibrationQA","Calibration QA overview",1400,600);
  c->Divide(2,1);

  // Extra right margin on both COLZ pads -- otherwise the z-axis (palette) power-of-10
  // labels get clipped/overlap the next pad in this Divide(2,1) grid.
  c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(calib3D){
    // Project down to (beta*gamma, ln(dE/dx)-prediction), summing over eta -- this is
    // the actual calibration check: a flat band centered at 0 means the Bichsel
    // prediction matches the data at this V_Z slice.
    TH2* calibProj = (TH2*) calib3D->Project3D("zy");
    calibProj->SetTitle(Form("%s calibration (Vz bin 0): ln(dE/dx / prediction) vs #beta#gamma",a_species.c_str()));
    if(calibProj->GetEntries() == 0){
      // A found-but-empty histogram looks identical to a blank pad once drawn with
      // COLZ -- flag it explicitly rather than silently rendering nothing. Common
      // causes: this V_Z slice/species combo genuinely has no isolated tracks in this
      // file, or PicoBinner's calibration_matrix (source/PicoBinner.cxx) didn't book
      // this species when the file was made in single-species batch mode.
      cout << "!!! " << calibPath << " was found but has 0 entries -- panel 1 will be "
           << "blank. Check the file was produced with this species actually toggled "
           << "on in PicoBinner's calibration_matrix." << endl;
    }
    calibProj->Draw("COLZ");
  } else {
    // calib3D itself wasn't found at all (as opposed to found-but-empty above) --
    // this is the more likely reason panel 1 comes out blank: either
    // _PICOBINNER_VZ_CALIBRATION_ was off when this file was produced, the species
    // name/casing in calibName doesn't match what's actually in the file, or (for a
    // single-species production run) PicoBinner's calibration_matrix didn't book this
    // species at all. Run FindHistoInFile.C on this file to see what calib_* names
    // actually exist before assuming this macro is broken.
    cout << "!!! " << calibPath << " not found in " << a_inputFile
         << " -- panel 1 will be blank. Run FindHistoInFile.C to see what calib_* "
         << "histograms actually exist in this file." << endl;
    TLatex* missingLabel = new TLatex(0.5, 0.5, Form("%s not found", calibPath.c_str()));
    missingLabel->SetNDC();
    missingLabel->SetTextAlign(22);
    missingLabel->SetTextSize(0.045);
    missingLabel->Draw();
  }

  // log y only on panel 2 -- its y-axis is dE/dx (always positive, same as
  // PresentPID.C's dE/dx panels). Panel 1's y-axis is ln(dE/dx/prediction), which is
  // centered on 0 and goes negative, so a log scale isn't valid there.
  c->cd(2); gPad->SetLogx(); gPad->SetLogy(); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(isolated){
    isolated->SetTitle(Form("Isolated (BTOF-tagged) dE/dx for %s+",a_species.c_str()));
    isolated->Draw("COLZ");
  }

  string outPath = a_outDir + "/PresentCalibrationQA.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
