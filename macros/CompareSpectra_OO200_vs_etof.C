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
// NOTE ON WHAT'S PLOTTED: both runs were done with combineImageDirectories=true and
// all a_rapBin_* args left at -1 ("AllRap" mode), meaning rapidity was NOT split into
// differential bins -- every populated 2D (rapidity x mT-m0) bin collapses onto a
// SINGLE rapidity slice. So there is no real rapidity axis to plot a curve against;
// the only real degree of freedom across these runs is centrality. This plots, per
// stage (TPC, BTOF), one point per (charge, centrality) -- 2 x 5 = 10 points -- at
// x=centrality bin, y=etof/OO200 (mean ratio across that bin's populated 2D cells),
// with the RMS-based spread as the error bar. Two plots total.
//
// Ends by printing a single PASS/FAIL summary line across everything compared.

#include "TFile.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
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
// summary. Fills a_outMeanRatio/a_outMeanRatioErr (standard error of the mean, across
// all populated 2D bins) for the caller to turn into one plotted point; returns false
// if there was nothing usable to compare (missing histos, binning mismatch, or no
// populated bins), true otherwise (including a real PASS/FAIL verdict via a_outPass).
bool compareOneHisto(TH2D* a_hEtof, TH2D* a_hOO200, string a_label, double a_relTolerance,
                     double& a_outMeanRatio, double& a_outMeanRatioErr, bool& a_outPass){
  a_outPass = true;

  if(!a_hEtof || !a_hOO200){
    cout << "  SKIP (missing on one side): " << a_label << endl;
    return false;
  }

  int nx = a_hEtof->GetNbinsX();
  int ny = a_hEtof->GetNbinsY();
  if(nx != a_hOO200->GetNbinsX() || ny != a_hOO200->GetNbinsY()){
    cout << "  [BINNING MISMATCH] " << a_label << ": etof(" << nx << "x" << ny
         << ") vs OO200(" << a_hOO200->GetNbinsX() << "x" << a_hOO200->GetNbinsY()
         << ") -- can't compare bin-by-bin, check SetCutClass binning matches" << endl;
    a_outPass = false;
    return false;
  }

  int nComparedBins = 0, nBadBins = 0, nOneSidedBins = 0;
  double worstDev = -1;
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
    return false;
  }

  int nTwoSidedBins = nComparedBins - nOneSidedBins;
  double meanRatio = sumRatio / (nTwoSidedBins > 0 ? nTwoSidedBins : 1);
  double rmsRatio  = sqrt(max(0.0, sumRatio2/(nTwoSidedBins > 0 ? nTwoSidedBins : 1) - meanRatio*meanRatio));
  double semRatio  = rmsRatio / sqrt((double) max(1, nTwoSidedBins));

  bool pass = (nBadBins == 0);
  a_outPass = pass;

  cout << "  " << a_label << ": " << nComparedBins << " populated bins, "
       << nOneSidedBins << " one-sided-only, "
       << nBadBins << " outside " << (a_relTolerance*100) << "% tol, "
       << "mean ratio=" << meanRatio << ", RMS=" << rmsRatio
       << ", worst bin (" << worstX << "," << worstY << ") dev=" << (worstDev < 0 ? 0 : worstDev)
       << "  -> " << (pass ? "PASS" : "FAIL") << endl;

  a_outMeanRatio    = meanRatio;
  a_outMeanRatioErr = semRatio;
  return true;
}

void runStage(TFile* a_fEtof, TFile* a_fOO200, string a_stageLabel, string a_histoPrefix,
              string a_fitDataDirName, string a_particleName, vector<int> a_centralities,
              double a_relTolerance, string a_outDir, bool& a_allPass){
  if(!a_fEtof || !a_fOO200) return;

  cout << "=== " << a_stageLabel << " spectra (" << a_particleName << ") ===" << endl;
  string dirPath = Form("%s/%s", a_fitDataDirName.c_str(), a_particleName.c_str());
  vector<string> charges = {"Plus","Minus"};
  int chargeColors[2] = {kRed, kBlue};
  int chargeMarkers[2] = {20, 21}; // filled circle, filled square
  double chargeXOffset[2] = {-0.08, 0.08}; // nudge apart so same-centrality points don't overlap

  vector<TGraphErrors*> graphs;

  for(int chargeIdx = 0; chargeIdx < (int) charges.size(); chargeIdx++){
    string charge = charges[chargeIdx];
    vector<double> xVals, yVals, yErrs;

    for(int cent : a_centralities){
      string histoName = Form("%s_%s%s_Cent%02d", a_histoPrefix.c_str(), a_particleName.c_str(), charge.c_str(), cent);
      TH2D* hEtof  = getHisto2D(a_fEtof,  dirPath, histoName);
      TH2D* hOO200 = getHisto2D(a_fOO200, dirPath, histoName);
      string label = Form("%s_%s_%s_Cent%02d", a_stageLabel.c_str(), a_particleName.c_str(), charge.c_str(), cent);

      double meanRatio = 0, meanRatioErr = 0;
      bool havePoint = false, pass = true;
      havePoint = compareOneHisto(hEtof, hOO200, label, a_relTolerance, meanRatio, meanRatioErr, pass);
      if(!pass) a_allPass = false;

      if(havePoint){
        xVals.push_back(cent + chargeXOffset[chargeIdx]);
        yVals.push_back(meanRatio);
        yErrs.push_back(meanRatioErr);
      }
    }

    if(!xVals.empty()){
      TGraphErrors* g = new TGraphErrors((int) xVals.size(), xVals.data(), yVals.data(), nullptr, yErrs.data());
      g->SetName(Form("g_%s_%s", a_stageLabel.c_str(), charge.c_str()));
      g->SetLineColor(chargeColors[chargeIdx]);
      g->SetMarkerColor(chargeColors[chargeIdx]);
      g->SetMarkerStyle(chargeMarkers[chargeIdx]);
      g->SetMarkerSize(1.8);
      graphs.push_back(g);
    }
  }

  if(graphs.empty()){
    cout << "  (nothing to plot for " << a_stageLabel << ")" << endl;
    return;
  }

  gSystem->mkdir(a_outDir.c_str(), true);

  TCanvas* c = new TCanvas(Form("c_%s_combined", a_stageLabel.c_str()),
                            Form("%s ratio summary", a_stageLabel.c_str()), 900, 700);
  c->SetGridy();

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(Form("%s: etof / OO200 (%s);Centrality bin;etof / OO200", a_stageLabel.c_str(), a_particleName.c_str()));

  TLegend* leg = new TLegend(0.72, 0.78, 0.97, 0.9);
  vector<string> chargeLabels = {"Proton+", "Proton-"};
  for(size_t i = 0; i < graphs.size(); i++){
    mg->Add(graphs[i], "P");
    leg->AddEntry(graphs[i], chargeLabels[i].c_str(), "p");
  }

  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0.5, 1.5);
  mg->GetXaxis()->SetLimits(-0.5, (double) *max_element(a_centralities.begin(), a_centralities.end()) + 0.5);
  gPad->Modified();
  gPad->Update();

  TLine* line = new TLine(mg->GetXaxis()->GetXmin(), 1, mg->GetXaxis()->GetXmax(), 1);
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
