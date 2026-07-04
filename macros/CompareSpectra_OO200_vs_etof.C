// CompareSpectra_OO200_vs_etof.C
//
// Verifies that the lightflavorspectra_OO200 port reproduces the original
// lightflavorspectra_etof repo's results, bin-for-bin, on the same O+O 200 GeV
// Proton yields file. Compares the ZFitter-stage output spectra histograms written
// by each repo's RunZFitter(_Collider) macro:
//
//   TPC stage:  DeDx_FitData/<Species>/dEdxSpectra_<Species><Charge>_Cent%02d       (TH2D: rapBin x mTm0Bin)
//   BTOF stage: BTOF_FitData/<Species>/InvBetaBTOFSpectra_<Species><Charge>_Cent%02d (TH2D: rapBin x mTm0Bin)
//
// This is plain histogram reading -- no compiled ZFitter/CutClass libraries are
// needed, so this can be run directly with no LoadZFitterLibs.C step, from any
// directory:
//
//   root -l -q 'CompareSpectra_OO200_vs_etof.C("/Users/aliggett/data/OO/etof_compare/spectra_200GeV_2024_09_12_000_Proton_TPC_RapBin_-1_-1_-1_-1_-1_-1.root","/Users/aliggett/data/OO/etof_compare/spectra_200GeV_2024_09_12_000_Proton_BTOF_RapBin_-1_-1_-1_-1_-1_-1.root","/Users/aliggett/data/OO/etof_compare/spectra_OO200_Proton_TPC_RapBin_-1_-1_-1_-1_-1_-1.root","/Users/aliggett/data/OO/etof_compare/spectra_OO200_Proton_BTOF_RapBin_-1_-1_-1_-1_-1_-1.root")'
//
// (keep it all on one line -- ROOT's ACLiC arg parser chokes on line-continuation
// backslashes inside a quoted macro invocation)
//
// Check the actual filenames with `ls` first -- currentDate/rap-bin fields embedded
// in the names must match exactly what your two RunZFitter runs actually produced.
// Pass "" for either the TPC or BTOF pair to skip that comparison.
//
// For every (charge, centrality) bin pair it finds in both files, this numerically
// compares every populated 2D bin and reports the fraction outside a_relTolerance
// (default 5%), plus the worst single-bin deviation. All of a stage's rapidity-
// integrated ratio curves (one per charge x centrality) are then overlaid onto a
// SINGLE combined plot per stage -- two plots total, TPC and BTOF -- rather than one
// plot per bin.
//
// Ends by printing a single PASS/FAIL summary line across everything compared.

#include "TFile.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

TH2D* getHisto2D(TFile* f, string a_dirPath, string a_histoName){
  if(!f) return nullptr;
  TDirectory* dir = (TDirectory*) f->Get(a_dirPath.c_str());
  if(!dir){
    cout << "  [MISSING DIR]   " << a_dirPath << "  (in " << f->GetName() << ")" << endl;
    return nullptr;
  }
  TH2D* h = (TH2D*) dir->Get(a_histoName.c_str());
  if(!h){
    cout << "  [MISSING HISTO] " << a_dirPath << "/" << a_histoName << "  (in " << f->GetName() << ")" << endl;
  }
  return h;
}

// Compares one pair of same-named 2D histograms bin-by-bin and prints a numeric
// summary. On success, fills a_outRatio with the rapidity-integrated (summed over
// mT-m0) ratio histogram (etof/OO200) for the caller to collect and overlay later;
// leaves it null if there was nothing to compare. Returns true if every populated
// bin agreed within a_relTolerance.
bool compareOneHisto(TH2D* a_hEtof, TH2D* a_hOO200, string a_label, double a_relTolerance, TH1D*& a_outRatio){
  a_outRatio = nullptr;

  if(!a_hEtof || !a_hOO200){
    cout << "  SKIP (missing on one side): " << a_label << endl;
    return true; // don't fail the whole run over a bin/charge that wasn't fit
  }

  int nx = a_hEtof->GetNbinsX();
  int ny = a_hEtof->GetNbinsY();
  if(nx != a_hOO200->GetNbinsX() || ny != a_hOO200->GetNbinsY()){
    cout << "  [BINNING MISMATCH] " << a_label << ": etof(" << nx << "x" << ny
         << ") vs OO200(" << a_hOO200->GetNbinsX() << "x" << a_hOO200->GetNbinsY()
         << ") -- can't compare bin-by-bin, check SetCutClass binning matches" << endl;
    return false;
  }

  int nComparedBins = 0, nBadBins = 0, nOneSidedBins = 0;
  double worstDev = 0;
  int worstX = -1, worstY = -1;
  double sumRatio = 0, sumRatio2 = 0;

  for(int ix = 1; ix <= nx; ix++){
    for(int iy = 1; iy <= ny; iy++){
      double vEtof  = a_hEtof->GetBinContent(ix,iy);
      double vOO200 = a_hOO200->GetBinContent(ix,iy);
      if(vEtof == 0 && vOO200 == 0) continue; // both empty -- not a meaningful bin

      if(vEtof == 0 || vOO200 == 0){
        // one side has content, the other doesn't -- always a mismatch
        nOneSidedBins++;
        nComparedBins++;
        continue;
      }

      double ratio = vEtof / vOO200;
      nComparedBins++;
      sumRatio  += ratio;
      sumRatio2 += ratio*ratio;
      double dev = fabs(ratio - 1.0);
      if(dev > a_relTolerance) nBadBins++;
      if(dev > worstDev){ worstDev = dev; worstX = ix; worstY = iy; }
    }
  }

  nBadBins += nOneSidedBins;

  if(nComparedBins == 0){
    cout << "  " << a_label << ": no populated bins in either file -- nothing to compare" << endl;
    return true;
  }

  int nTwoSidedBins = nComparedBins - nOneSidedBins;
  double meanRatio = sumRatio / (nTwoSidedBins > 0 ? nTwoSidedBins : 1);
  double rmsRatio  = sqrt(max(0.0, sumRatio2/(nTwoSidedBins > 0 ? nTwoSidedBins : 1) - meanRatio*meanRatio));

  bool pass = (nBadBins == 0);

  cout << "  " << a_label << ": " << nComparedBins << " populated bins, "
       << nOneSidedBins << " one-sided-only, "
       << nBadBins << " outside " << (a_relTolerance*100) << "% tol, "
       << "mean ratio=" << meanRatio << ", RMS=" << rmsRatio
       << ", worst bin (" << worstX << "," << worstY << ") dev=" << worstDev
       << "  -> " << (pass ? "PASS" : "FAIL") << endl;

  // rapidity-integrated (summed over mT-m0) ratio, to be overlaid with the other
  // bins for this stage on one combined plot
  TH1D* pEtof  = a_hEtof->ProjectionX(Form("px_etof_%s",  a_label.c_str()));
  TH1D* pOO200 = a_hOO200->ProjectionX(Form("px_oo200_%s", a_label.c_str()));
  TH1D* pRatio = (TH1D*) pEtof->Clone(Form("ratio1D_%s", a_label.c_str()));
  pRatio->Divide(pOO200);
  pRatio->SetDirectory(0);
  a_outRatio = pRatio;

  return pass;
}

// Colors/markers for up to 10 (charge x centrality) combinations on one combined plot
int stylePalette[10] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kViolet+1, kAzure+2, kSpring+4};

void runStage(TFile* a_fEtof, TFile* a_fOO200, string a_stageLabel, string a_histoPrefix,
              string a_fitDataDirName, string a_particleName, vector<int> a_centralities,
              double a_relTolerance, string a_outDir, bool& a_allPass){
  if(!a_fEtof || !a_fOO200) return;

  cout << "=== " << a_stageLabel << " spectra (" << a_particleName << ") ===" << endl;
  string dirPath = Form("%s/%s", a_fitDataDirName.c_str(), a_particleName.c_str());
  vector<string> charges = {"Plus","Minus"};

  vector<TH1D*> ratios;
  vector<string> ratioLabels;

  for(auto& charge : charges){
    for(int cent : a_centralities){
      string histoName = Form("%s_%s%s_Cent%02d", a_histoPrefix.c_str(), a_particleName.c_str(), charge.c_str(), cent);
      TH2D* hEtof  = getHisto2D(a_fEtof,  dirPath, histoName);
      TH2D* hOO200 = getHisto2D(a_fOO200, dirPath, histoName);
      string label = Form("%s_%s_%s_Cent%02d", a_stageLabel.c_str(), a_particleName.c_str(), charge.c_str(), cent);

      TH1D* ratio = nullptr;
      bool pass = compareOneHisto(hEtof, hOO200, label, a_relTolerance, ratio);
      if(!pass) a_allPass = false;

      if(ratio){
        ratios.push_back(ratio);
        ratioLabels.push_back(Form("%s %s", charge.c_str(), Form("Cent%02d", cent)));
      }
    }
  }

  if(ratios.empty()){
    cout << "  (nothing to plot for " << a_stageLabel << ")" << endl;
    return;
  }

  gSystem->mkdir(a_outDir.c_str(), true);

  TCanvas* c = new TCanvas(Form("c_%s_combined", a_stageLabel.c_str()),
                            Form("%s ratio summary", a_stageLabel.c_str()), 900, 700);
  c->SetGridy();

  double yMin = 0.5, yMax = 1.5;
  ratios[0]->SetTitle(Form("%s: etof / OO200  (%s, all charges/centralities)", a_stageLabel.c_str(), a_particleName.c_str()));
  ratios[0]->GetYaxis()->SetTitle("etof / OO200");
  ratios[0]->GetYaxis()->SetRangeUser(yMin, yMax);

  TLegend* leg = new TLegend(0.72, 0.15, 0.97, 0.9);
  leg->SetTextSize(0.025);

  for(size_t i = 0; i < ratios.size(); i++){
    int color = stylePalette[i % 10];
    ratios[i]->SetLineColor(color);
    ratios[i]->SetMarkerColor(color);
    ratios[i]->SetMarkerStyle(20 + (i / 10)); // cycle marker shape if ever >10 curves
    ratios[i]->SetMarkerSize(1.6);
    ratios[i]->Draw(i == 0 ? "P" : "P SAME");
    leg->AddEntry(ratios[i], ratioLabels[i].c_str(), "p");
  }

  TLine* line = new TLine(ratios[0]->GetXaxis()->GetXmin(), 1, ratios[0]->GetXaxis()->GetXmax(), 1);
  line->SetLineStyle(2);
  line->Draw();
  leg->Draw();

  c->SaveAs(Form("%s/%s_ratio_summary.png", a_outDir.c_str(), a_stageLabel.c_str()));
}

void CompareSpectra_OO200_vs_etof(
    string a_etofFile_TPC   = "",
    string a_etofFile_BTOF  = "",
    string a_oo200File_TPC  = "",
    string a_oo200File_BTOF = "",
    string a_particleName   = "Proton",
    vector<int> a_centralities = {0,1,2,3,4},
    string a_outDir         = "./spectraComparison_OO200_vs_etof/",
    double a_relTolerance   = 0.05){

  gStyle->SetOptStat(0);
  bool allPass = true;

  TFile* fEtofTPC   = a_etofFile_TPC   != "" ? TFile::Open(a_etofFile_TPC.c_str())   : nullptr;
  TFile* fEtofBTOF  = a_etofFile_BTOF  != "" ? TFile::Open(a_etofFile_BTOF.c_str())  : nullptr;
  TFile* fOO200TPC  = a_oo200File_TPC  != "" ? TFile::Open(a_oo200File_TPC.c_str())  : nullptr;
  TFile* fOO200BTOF = a_oo200File_BTOF != "" ? TFile::Open(a_oo200File_BTOF.c_str()) : nullptr;

  if(a_etofFile_TPC  != "" && (!fEtofTPC   || fEtofTPC->IsZombie()))   { cout << "Could not open " << a_etofFile_TPC   << endl; return; }
  if(a_etofFile_BTOF != "" && (!fEtofBTOF  || fEtofBTOF->IsZombie()))  { cout << "Could not open " << a_etofFile_BTOF  << endl; return; }
  if(a_oo200File_TPC != "" && (!fOO200TPC  || fOO200TPC->IsZombie()))  { cout << "Could not open " << a_oo200File_TPC  << endl; return; }
  if(a_oo200File_BTOF!= "" && (!fOO200BTOF || fOO200BTOF->IsZombie())) { cout << "Could not open " << a_oo200File_BTOF << endl; return; }

  gSystem->mkdir(a_outDir.c_str(), true);

  runStage(fEtofTPC,  fOO200TPC,  "TPC",  "dEdxSpectra",         "DeDx_FitData", a_particleName, a_centralities, a_relTolerance, a_outDir, allPass);
  runStage(fEtofBTOF, fOO200BTOF, "BTOF", "InvBetaBTOFSpectra",  "BTOF_FitData", a_particleName, a_centralities, a_relTolerance, a_outDir, allPass);

  cout << endl << "=================================================" << endl;
  cout << (allPass ? "ALL COMPARISONS PASSED within tolerance." : "SOME COMPARISONS FAILED -- see FAIL/MISSING/MISMATCH lines above.") << endl;
  cout << "Plots written to: " << a_outDir << endl;
  cout << "  " << a_outDir << "/TPC_ratio_summary.png" << endl;
  cout << "  " << a_outDir << "/BTOF_ratio_summary.png" << endl;
  cout << "=================================================" << endl;
}
