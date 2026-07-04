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
// PLOT LAYOUT: one plot per (stage, populated rapidity bin), Proton+ only. Each plot
// has one point per centrality bin (etof/OO200 mean ratio, RMS-based error bar),
// restricted to that single rapidity row of the 2D (rapidity x mT-m0) histogram.
// "Populated rapidity bin" is auto-detected from the data (union, across all
// centralities and both files, of rapidity bins with any nonzero content) -- both
// runs so far were done with combineImageDirectories=true / all rapBin args at -1
// ("AllRap" mode), so there is currently only one such bin, but this generalizes
// cleanly to a real differential-rapidity run later.
//
// Ends by printing a single PASS/FAIL summary line across everything compared.

#include "TFile.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <vector>
#include <string>
#include <set>
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

// Compares bin content along the mT-m0 axis at a single fixed rapidity bin
// (a_ixFixed) between two same-binned 2D histograms, prints a numeric summary, and
// fills a_outMeanRatio/a_outMeanRatioErr (standard error of the mean across that
// row's populated cells). Returns false if there was nothing usable to compare.
bool compareRowStats(TH2D* a_hEtof, TH2D* a_hOO200, int a_ixFixed, string a_label, double a_relTolerance,
                     double& a_outMeanRatio, double& a_outMeanRatioErr, bool& a_outPass){
  a_outPass = true;

  if(!a_hEtof || !a_hOO200) return false;

  int ny = a_hEtof->GetNbinsY();
  int nComparedBins = 0, nBadBins = 0, nOneSidedBins = 0;
  double worstDev = -1;
  double sumRatio = 0, sumRatio2 = 0;

  for(int iy = 1; iy <= ny; iy++){
    double vEtof  = a_hEtof->GetBinContent(a_ixFixed,iy);
    double vOO200 = a_hOO200->GetBinContent(a_ixFixed,iy);
    if(vEtof == 0 && vOO200 == 0) continue; // both empty -- not a meaningful bin

    if(vEtof == 0 || vOO200 == 0){
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
    if(dev > worstDev) worstDev = dev;
  }

  nBadBins += nOneSidedBins;

  if(nComparedBins == 0){
    cout << "  " << a_label << ": no populated mT-m0 bins at this rapidity bin -- nothing to compare" << endl;
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
       << ", worst dev=" << (worstDev < 0 ? 0 : worstDev)
       << "  -> " << (pass ? "PASS" : "FAIL") << endl;

  a_outMeanRatio    = meanRatio;
  a_outMeanRatioErr = semRatio;
  return true;
}

void runStage(TFile* a_fEtof, TFile* a_fOO200, string a_stageLabel, string a_histoPrefix,
              string a_fitDataDirName, string a_particleName, string a_charge,
              vector<int> a_centralities, double a_relTolerance, string a_outDir, bool& a_allPass){
  if(!a_fEtof || !a_fOO200) return;

  cout << "=== " << a_stageLabel << " spectra (" << a_particleName << a_charge << ") ===" << endl;
  string dirPath = Form("%s/%s", a_fitDataDirName.c_str(), a_particleName.c_str());

  // Grab all centralities' histograms once, and auto-detect which rapidity (x) bins
  // actually have any content in either file, at any centrality.
  vector<TH2D*> etofHistos(a_centralities.size(), nullptr);
  vector<TH2D*> oo200Histos(a_centralities.size(), nullptr);
  set<int> populatedRapBins;
  int nx = -1;

  for(size_t i = 0; i < a_centralities.size(); i++){
    int cent = a_centralities[i];
    string histoName = Form("%s_%s%s_Cent%02d", a_histoPrefix.c_str(), a_particleName.c_str(), a_charge.c_str(), cent);
    TH2D* hEtof  = getHisto2D(a_fEtof,  dirPath, histoName);
    TH2D* hOO200 = getHisto2D(a_fOO200, dirPath, histoName);
    etofHistos[i]  = hEtof;
    oo200Histos[i] = hOO200;
    if(!hEtof || !hOO200) continue;

    if(nx == -1) nx = hEtof->GetNbinsX();
    if(hEtof->GetNbinsX() != hOO200->GetNbinsX() || hEtof->GetNbinsY() != hOO200->GetNbinsY()){
      cout << "  [BINNING MISMATCH] Cent" << Form("%02d",cent) << ": etof(" << hEtof->GetNbinsX() << "x" << hEtof->GetNbinsY()
           << ") vs OO200(" << hOO200->GetNbinsX() << "x" << hOO200->GetNbinsY() << ")" << endl;
      a_allPass = false;
      continue;
    }

    for(int ix = 1; ix <= hEtof->GetNbinsX(); ix++){
      for(int iy = 1; iy <= hEtof->GetNbinsY(); iy++){
        if(hEtof->GetBinContent(ix,iy) != 0 || hOO200->GetBinContent(ix,iy) != 0){
          populatedRapBins.insert(ix);
          break;
        }
      }
    }
  }

  if(populatedRapBins.empty()){
    cout << "  (no populated rapidity bins found -- nothing to plot for " << a_stageLabel << ")" << endl;
    return;
  }

  gSystem->mkdir(a_outDir.c_str(), true);

  for(int ixRap : populatedRapBins){
    vector<double> xVals, yVals, yErrs;

    for(size_t i = 0; i < a_centralities.size(); i++){
      int cent = a_centralities[i];
      string label = Form("%s_%s%s_RapBin%02d_Cent%02d", a_stageLabel.c_str(), a_particleName.c_str(), a_charge.c_str(), ixRap, cent);

      double meanRatio = 0, meanRatioErr = 0;
      bool pass = true;
      bool havePoint = compareRowStats(etofHistos[i], oo200Histos[i], ixRap, label, a_relTolerance, meanRatio, meanRatioErr, pass);
      if(!pass) a_allPass = false;

      if(havePoint){
        xVals.push_back(cent);
        yVals.push_back(meanRatio);
        yErrs.push_back(meanRatioErr);
      }
    }

    if(xVals.empty()) continue;

    // find a representative histogram to pull rapidity-bin edges from, for the title
    TH2D* refHisto = nullptr;
    for(auto* h : etofHistos) if(h){ refHisto = h; break; }
    double yLo = refHisto ? refHisto->GetXaxis()->GetBinLowEdge(ixRap) : 0;
    double yHi = refHisto ? refHisto->GetXaxis()->GetBinUpEdge(ixRap) : 0;

    TGraphErrors* g = new TGraphErrors((int) xVals.size(), xVals.data(), yVals.data(), nullptr, yErrs.data());
    g->SetName(Form("g_%s_RapBin%02d", a_stageLabel.c_str(), ixRap));
    g->SetTitle(Form("%s: etof / OO200, %s%s, RapBin %02d (y #in [%.2f, %.2f]);Centrality bin;etof / OO200",
                      a_stageLabel.c_str(), a_particleName.c_str(), a_charge.c_str(), ixRap, yLo, yHi));
    g->SetLineColor(kRed);
    g->SetMarkerColor(kRed);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.8);

    TCanvas* c = new TCanvas(Form("c_%s_RapBin%02d", a_stageLabel.c_str(), ixRap),
                              Form("%s RapBin %02d ratio summary", a_stageLabel.c_str(), ixRap), 900, 700);
    c->SetGridy();

    g->Draw("AP");
    g->GetYaxis()->SetRangeUser(0.5, 1.5);
    g->GetXaxis()->SetLimits(a_centralities.front() - 0.5, a_centralities.back() + 0.5);
    gPad->Modified();
    gPad->Update();

    TLine* line = new TLine(g->GetXaxis()->GetXmin(), 1, g->GetXaxis()->GetXmax(), 1);
    line->SetLineStyle(2);
    line->Draw();

    string outName = Form("%s/%s_%s%s_RapBin%02d_ratio_summary.png", a_outDir.c_str(),
                           a_stageLabel.c_str(), a_particleName.c_str(), a_charge.c_str(), ixRap);
    c->SaveAs(outName.c_str());
  }
}

void CompareSpectra_OO200_vs_etof(
    string a_etofFile_TPC   = "",
    string a_etofFile_BTOF  = "",
    string a_oo200File_TPC  = "",
    string a_oo200File_BTOF = "",
    string a_particleName   = "Proton",
    string a_charge         = "Plus",
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

  runStage(fEtofTPC,  fOO200TPC,  "TPC",  "dEdxSpectra",        "DeDx_FitData", a_particleName, a_charge, a_centralities, a_relTolerance, a_outDir, allPass);
  runStage(fEtofBTOF, fOO200BTOF, "BTOF", "InvBetaBTOFSpectra", "BTOF_FitData", a_particleName, a_charge, a_centralities, a_relTolerance, a_outDir, allPass);

  cout << endl << "=================================================" << endl;
  cout << (allPass ? "ALL COMPARISONS PASSED within tolerance." : "SOME COMPARISONS FAILED -- see FAIL/MISSING/MISMATCH lines above.") << endl;
  cout << "Plots written to: " << a_outDir << " (one PNG per stage x populated rapidity bin)" << endl;
  cout << "=================================================" << endl;
}
