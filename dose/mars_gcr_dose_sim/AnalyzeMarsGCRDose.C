// AnalyzeMarsGCRDose.C
//
// Post-processing for the MarsGCRDoseSim campaign (see main.cc). Reads the
// per-(Z,A,energy,layer) dose_Gy_per_primary / doseEq_Sv_per_primary ntuple
// written by RunAction for each physics-list variant, combines species and
// energy points into ONE dose-equivalent-vs-depth curve per variant, and
// derives a "days to NASA 600 mSv career limit vs depth" curve.
//
// UNITS: depth axis is LINEAR regolith depth in cm (not areal density,
// g/cm^2) -- converted from the Geant4 app's areal-density-defined layer
// geometry via the regolith bulk density (see DetectorParameters.hh /
// RunAction.cc). Total scanned depth is 100 g/cm^2 / 1.52 g/cm^3 = 65.79 cm.
//
// WEIGHTING -- read this before trusting the absolute curve shape:
//
//   Species weights: PDG Table 30.1 relative abundances, with documented
//   within-group split approximations for C/O, Ne/F, Mg/Na, Si/Al, Fe/Ni --
//   IDENTICAL table to GCRSpectrum.hh (kept in sync manually; if you change
//   one, change the other). See GCRSpectrum.hh's header comment for the
//   full sourcing/approximation writeup, not repeated here.
//
//   Energy weights: EVERY one of the 50 log-spaced energy grid points gets
//   equal weight, per species. Because the grid is log-spaced (not
//   linear-spaced), equal-weight-per-point approximates a flux spectrum
//   with dJ/dlogE ~ constant, i.e. dJ/dE ~ 1/E -- a real but *rough*
//   stand-in for the true GCR differential spectrum shape (which is peaked
//   around a few hundred MeV/u to a few GeV/u and falls off more steeply
//   above that, not a clean 1/E law at all energies). This under-weights
//   the peak-energy region and over-weights the extreme ends of the 10
//   MeV/u - 200 GeV/u range relative to a real GCR spectrum. Flagged
//   clearly rather than silently accepted -- if you want to test
//   sensitivity to this, EnergyWeight() below is the one place to edit,
//   and doing so does NOT require re-running the (expensive) Geant4
//   campaign, only re-running this macro.
//
// ABSOLUTE SCALE -- anchored to real measured data, not derived from the
//   (relative-only) weights above:
//   The species+energy-weighted curve from Geant4 has an arbitrary overall
//   normalization (it's dose PER PRIMARY, summed with relative weights that
//   aren't an absolute flux). To get real mSv/year, this macro rescales the
//   ENTIRE weighted curve by a single constant so that its value at the
//   shallowest scored depth (layer 0, nominally 0-1.32 cm, i.e.
//   effectively unshielded) matches the actual MEASURED Mars-surface GCR
//   dose-equivalent rate from MSL/Curiosity's RAD instrument:
//     0.64 +/- 0.12 mSv/day = 233.7 mSv/year
//     (Hassler et al. 2014, Science, "Mars' Surface Radiation Environment
//     Measured with the Mars Science Laboratory's Curiosity Rover"; value
//     cross-checked against Reitz et al. 2016, DLR/ICES, which reports the
//     same 0.64+/-0.12 mSv/day figure for Aug 7 2012 - Jun 1 2013.)
//   This means Geant4 supplies the SHAPE (how dose changes as regolith
//   shielding is added on top of the real Mars surface condition) and the
//   real RAD measurement supplies the ABSOLUTE SCALE at zero added
//   shielding -- deliberately avoiding inventing an absolute GCR flux
//   normalization from the simplified reference-ion/relative-abundance
//   approach alone.
//
// NASA CAREER LIMIT: 600 mSv effective dose, flat for all astronauts
// regardless of age/sex (NASA's 2021 policy revision, replacing the prior
// age/sex-dependent REID-based limits) -- same number already marked on
// Andrew's existing motivation slide.
//
// STATUS: this macro has not been run (depends on MarsGCRDoseSim's .root
// output, which doesn't exist until you build and run that application).
// Paren/brace/bracket balance was checked with a script before delivery,
// same as the other ROOT macros this session, but logical correctness
// (branch names, TTree access pattern) has not been exercised against a
// real file.
//
// Usage:
//   root -l -q 'AnalyzeMarsGCRDose.C("MarsGCRDose_QMD.root","MarsGCRDose_INCLXX.root")'

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPad.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>

using namespace std;

//============================================================================
// Species relative-abundance weights -- see file header. Keyed by Z (unique
// among our 8 reference ions, no need for the (Z,A) pair here).
//============================================================================
map<int, double> SpeciesRelAbundanceByZ(){
  return {
    {1,  550.0},  // H
    {2,  34.0},   // He
    {6,  1.20},   // C
    {8,  1.00},   // O
    {10, 0.27},   // Ne
    {12, 0.187},  // Mg
    {14, 0.1615}, // Si
    {26, 0.114},  // Fe
  };
}

// Equal weight per (log-spaced) energy grid point -- see file header caveat.
double EnergyWeight(double /*energyMeVPerNucleon*/){
  return 1.0;
}

// Real MSL/Curiosity RAD-measured Mars surface GCR dose-equivalent rate --
// see file header for citation.
const double kRADSurfaceDoseRate_mSvPerDay = 0.64;
const double kRADSurfaceDoseRate_SvPerYear = (kRADSurfaceDoseRate_mSvPerDay / 1000.0) * 365.25;

const double kNASACareerLimit_Sv = 0.600;
const int    kNLayers = 50;

// Total linear depth spanned by the 50-layer stack: 100 g/cm^2 (areal
// density, the geometry-defining parameter in DetectorParameters.hh) /
// 1.52 g/cm^3 (nominal Mars regolith bulk density) = 65.79 cm.
const double kTotalDepthCm = (kNLayers * 2.0) / 1.52;

//============================================================================
// LoadWeightedDepthCurve: read one physics-variant's ROOT file, produce the
// species+energy-weighted (but not yet RAD-anchored) doseEq_Sv_per_primary
// vs. layer-depth arrays.
//============================================================================
bool LoadWeightedDepthCurve(const string& a_rootFile,
                             vector<double>& a_depthMid_cm,
                             vector<double>& a_weightedDoseEq_perPrimary){
  TFile* f = TFile::Open(a_rootFile.c_str(), "READ");
  if(!f || f->IsZombie()){
    cout << "ERROR: could not open " << a_rootFile << endl;
    return false;
  }
  TTree* tree = (TTree*) f->Get("MarsGCRDose");
  if(!tree){
    cout << "ERROR: TTree \"MarsGCRDose\" not found in " << a_rootFile << endl;
    f->Close();
    return false;
  }

  Int_t    Z_in, A_in, physicsVariantCode_in, layerIndex_in, nPrimaries_in;
  Double_t energyMeVPerNucleon_in, layerDepthLo_in, layerDepthHi_in;
  Double_t dose_Gy_perPrim_in, doseEq_Sv_perPrim_in;
  tree->SetBranchAddress("Z", &Z_in);
  tree->SetBranchAddress("A", &A_in);
  tree->SetBranchAddress("energyMeVPerNucleon", &energyMeVPerNucleon_in);
  tree->SetBranchAddress("physicsVariantCode", &physicsVariantCode_in);
  tree->SetBranchAddress("layerIndex", &layerIndex_in);
  tree->SetBranchAddress("layerDepthLo_cm", &layerDepthLo_in);
  tree->SetBranchAddress("layerDepthHi_cm", &layerDepthHi_in);
  tree->SetBranchAddress("nPrimaries", &nPrimaries_in);
  tree->SetBranchAddress("dose_Gy_per_primary", &dose_Gy_perPrim_in);
  tree->SetBranchAddress("doseEq_Sv_per_primary", &doseEq_Sv_perPrim_in);

  map<int, double> relAbundance = SpeciesRelAbundanceByZ();

  vector<double> weightedSum(kNLayers, 0.0);
  vector<double> weightSum(kNLayers, 0.0);     // normalization (sum of weights actually seen)
  vector<double> depthMid(kNLayers, 0.0);
  vector<bool>   depthSet(kNLayers, false);

  Long64_t nEntries = tree->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    tree->GetEntry(i);
    if(layerIndex_in < 0 || layerIndex_in >= kNLayers) continue;

    auto it = relAbundance.find(Z_in);
    if(it == relAbundance.end()){
      cout << "WARNING: Z=" << Z_in << " not in the species weight table -- skipping entry." << endl;
      continue;
    }
    double w = it->second * EnergyWeight(energyMeVPerNucleon_in);

    weightedSum[layerIndex_in] += w * doseEq_Sv_perPrim_in;
    weightSum[layerIndex_in]   += w;

    if(!depthSet[layerIndex_in]){
      depthMid[layerIndex_in] = 0.5 * (layerDepthLo_in + layerDepthHi_in);
      depthSet[layerIndex_in] = true;
    }
  }
  f->Close();

  a_depthMid_cm.clear();
  a_weightedDoseEq_perPrimary.clear();
  for(int i = 0; i < kNLayers; i++){
    if(!depthSet[i] || weightSum[i] <= 0.0){
      cout << "WARNING: layer " << i << " has no entries in " << a_rootFile
           << " -- omitting from the curve (did every species/energy run complete?)" << endl;
      continue;
    }
    a_depthMid_cm.push_back(depthMid[i]);
    a_weightedDoseEq_perPrimary.push_back(weightedSum[i] / weightSum[i]);
  }
  return !a_depthMid_cm.empty();
}

void AnalyzeMarsGCRDose(string a_qmdFile = "MarsGCRDose_QMD.root",
                         string a_inclxxFile = "MarsGCRDose_INCLXX.root",
                         string a_outputDir = "."){
  vector<double> depthQMD, weightedQMD, depthINCLXX, weightedINCLXX;
  bool haveQMD    = LoadWeightedDepthCurve(a_qmdFile, depthQMD, weightedQMD);
  bool haveINCLXX = LoadWeightedDepthCurve(a_inclxxFile, depthINCLXX, weightedINCLXX);

  if(!haveQMD && !haveINCLXX){
    cout << "ERROR: neither input file produced a usable curve -- nothing to plot." << endl;
    return;
  }

  // RAD anchor: rescale so the shallowest depth point equals the real
  // measured Mars-surface GCR dose-equivalent rate (see file header).
  vector<double> doseRateQMD_SvPerYear, doseRateINCLXX_SvPerYear;
  if(haveQMD){
    double scale = kRADSurfaceDoseRate_SvPerYear / weightedQMD.front();
    for(double w : weightedQMD) doseRateQMD_SvPerYear.push_back(w * scale);
    cout << "QMD anchor scale factor = " << scale << endl;
  }
  if(haveINCLXX){
    double scale = kRADSurfaceDoseRate_SvPerYear / weightedINCLXX.front();
    for(double w : weightedINCLXX) doseRateINCLXX_SvPerYear.push_back(w * scale);
    cout << "INCLXX anchor scale factor = " << scale << endl;
  }

  // Derived: days to NASA 600 mSv career limit, at this dose rate.
  auto DaysToLimit = [&](double doseRate_SvPerYear){
    if(doseRate_SvPerYear <= 0.0) return -1.0;
    return kNASACareerLimit_Sv * 365.25 / doseRate_SvPerYear;
  };
  vector<double> daysQMD, daysINCLXX;
  for(double d : doseRateQMD_SvPerYear) daysQMD.push_back(DaysToLimit(d));
  for(double d : doseRateINCLXX_SvPerYear) daysINCLXX.push_back(DaysToLimit(d));

  //--------------------------------------------------------------------
  // Plot 1: dose-equivalent rate vs. regolith depth (mSv/year), the direct
  // Mars-regolith replacement for the borrowed Al-shielding slide figure.
  //--------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const char* colQMD    = "#2a78d6"; // dataviz palette slot 1, blue
  const char* colINCLXX = "#e34948"; // dataviz palette slot 8, red
  const char* SURFACE = "#fcfcfb";
  const char* TEXT_PRIMARY = "#0b0b0b";
  const char* TEXT_SECONDARY = "#52514e";
  const char* TEXT_MUTED = "#898781";
  const char* BASELINE = "#c3c2b7";

  TCanvas* c1 = new TCanvas("c_MarsGCRDose", "Mars GCR Dose vs Regolith Depth", 950, 750);
  c1->SetFillColor(TColor::GetColor(SURFACE));
  c1->SetFrameFillColor(TColor::GetColor(SURFACE));
  c1->SetLeftMargin(0.13); c1->SetRightMargin(0.06);
  c1->SetBottomMargin(0.14); c1->SetTopMargin(0.12); // only the title lives up here now, no subtitle

  double yMax = 0.0;
  for(double d : doseRateQMD_SvPerYear) yMax = max(yMax, d * 1000.0);
  for(double d : doseRateINCLXX_SvPerYear) yMax = max(yMax, d * 1000.0);
  yMax = max(yMax, 700.0); // headroom above the 600 mSv reference line

  TH1F* frame1 = c1->DrawFrame(0.0, 0.0, kTotalDepthCm, yMax * 1.15);
  frame1->GetXaxis()->SetTitle("Mars Regolith Shielding Depth (cm)");
  frame1->GetYaxis()->SetTitle("Dose Equivalent Rate (mSv/year)");
  frame1->GetXaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame1->GetYaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame1->GetXaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame1->GetYaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame1->GetXaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame1->GetYaxis()->SetAxisColor(TColor::GetColor(BASELINE));

  vector<TGraph*> keepAlive;
  TLegend* leg1 = new TLegend(0.55, 0.72, 0.93, 0.86);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.026);
  leg1->SetTextColor(TColor::GetColor(TEXT_SECONDARY));

  if(haveQMD){
    vector<double> y_mSv;
    for(double d : doseRateQMD_SvPerYear) y_mSv.push_back(d * 1000.0);
    TGraph* g = new TGraph((int) depthQMD.size(), &depthQMD[0], &y_mSv[0]);
    keepAlive.push_back(g);
    g->SetLineColor(TColor::GetColor(colQMD));
    g->SetLineWidth(3);
    g->Draw("L SAME");
    leg1->AddEntry(g, "Geant4 (QMD ion physics)", "l");
  }
  if(haveINCLXX){
    vector<double> y_mSv;
    for(double d : doseRateINCLXX_SvPerYear) y_mSv.push_back(d * 1000.0);
    TGraph* g = new TGraph((int) depthINCLXX.size(), &depthINCLXX[0], &y_mSv[0]);
    keepAlive.push_back(g);
    g->SetLineColor(TColor::GetColor(colINCLXX));
    g->SetLineWidth(3);
    g->Draw("L SAME");
    leg1->AddEntry(g, "Geant4 (INCLXX ion physics)", "l");
  }

  TLine* limitLine = new TLine(0.0, 600.0, kTotalDepthCm, 600.0);
  limitLine->SetLineColor(TColor::GetColor("#B8B8B8"));
  limitLine->SetLineStyle(2);
  limitLine->SetLineWidth(2);
  limitLine->Draw();
  TLatex limitLabel;
  limitLabel.SetTextColor(TColor::GetColor(TEXT_MUTED));
  limitLabel.SetTextSize(0.024);
  limitLabel.DrawLatex(0.02 * kTotalDepthCm, 610.0, "NASA career limit: 600 mSv total (not annual)");

  leg1->Draw();

  TLatex title1; title1.SetNDC(); title1.SetTextColor(TColor::GetColor(TEXT_PRIMARY));
  title1.SetTextFont(62); title1.SetTextSize(0.035); title1.SetTextAlign(21);
  title1.DrawLatex(0.5, 0.955, "Mars GCR Dose Equivalent vs. Regolith Shielding Depth");

  // No in-plot subtitle -- the method/caveat detail that used to live here
  // belongs in the slide's own figure caption instead.

  c1->RedrawAxis();
  string out1 = a_outputDir + "/MarsGCRDoseVsDepth.png";
  c1->SaveAs(out1.c_str());
  cout << "Wrote " << out1 << endl;

  //--------------------------------------------------------------------
  // Plot 2: derived -- days to NASA 600 mSv career limit vs. depth.
  //--------------------------------------------------------------------
  TCanvas* c2 = new TCanvas("c_MarsGCRDaysToLimit", "Days to NASA Career Limit", 950, 750);
  c2->SetLogy();
  c2->SetFillColor(TColor::GetColor(SURFACE));
  c2->SetFrameFillColor(TColor::GetColor(SURFACE));
  c2->SetLeftMargin(0.13); c2->SetRightMargin(0.06);
  c2->SetBottomMargin(0.14); c2->SetTopMargin(0.12); // only the title lives up here now, no subtitle

  double yMaxDays = 1.0;
  for(double d : daysQMD) if(d > 0) yMaxDays = max(yMaxDays, d);
  for(double d : daysINCLXX) if(d > 0) yMaxDays = max(yMaxDays, d);

  TH1F* frame2 = c2->DrawFrame(0.0, 10.0, kTotalDepthCm, yMaxDays * 3.0);
  frame2->GetXaxis()->SetTitle("Mars Regolith Shielding Depth (cm)");
  frame2->GetYaxis()->SetTitle("Days on Mars Surface to Reach 600 mSv Career Limit");
  frame2->GetXaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame2->GetYaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame2->GetXaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame2->GetYaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame2->GetXaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame2->GetYaxis()->SetAxisColor(TColor::GetColor(BASELINE));

  TLegend* leg2 = new TLegend(0.16, 0.72, 0.54, 0.86);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.026);
  leg2->SetTextColor(TColor::GetColor(TEXT_SECONDARY));

  if(haveQMD){
    TGraph* g = new TGraph((int) depthQMD.size(), &depthQMD[0], &daysQMD[0]);
    keepAlive.push_back(g);
    g->SetLineColor(TColor::GetColor(colQMD));
    g->SetLineWidth(3);
    g->Draw("L SAME");
    leg2->AddEntry(g, "Geant4 (QMD ion physics)", "l");
  }
  if(haveINCLXX){
    TGraph* g = new TGraph((int) depthINCLXX.size(), &depthINCLXX[0], &daysINCLXX[0]);
    keepAlive.push_back(g);
    g->SetLineColor(TColor::GetColor(colINCLXX));
    g->SetLineWidth(3);
    g->Draw("L SAME");
    leg2->AddEntry(g, "Geant4 (INCLXX ion physics)", "l");
  }
  leg2->Draw();

  TLatex title2; title2.SetNDC(); title2.SetTextColor(TColor::GetColor(TEXT_PRIMARY));
  title2.SetTextFont(62); title2.SetTextSize(0.035); title2.SetTextAlign(21);
  title2.DrawLatex(0.5, 0.955, "Days to NASA Career Radiation Limit vs. Regolith Depth");

  // No in-plot subtitle -- the method/caveat detail that used to live here
  // belongs in the slide's own figure caption instead.

  c2->RedrawAxis();
  string out2 = a_outputDir + "/MarsGCRDaysToLimitVsDepth.png";
  c2->SaveAs(out2.c_str());
  cout << "Wrote " << out2 << endl;
}
