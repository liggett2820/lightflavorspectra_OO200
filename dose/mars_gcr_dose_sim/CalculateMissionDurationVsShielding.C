// CalculateMissionDurationVsShielding.C
//
// Answers the actual mission-design question: given the NASA 600 mSv career
// radiation limit, how long can a crew stay on the Mars surface (and what's
// the resulting total round-trip mission length) as a function of ADDED
// REGOLITH SHIELDING DEPTH (cm) -- combining:
//
//   (1) the real transit dose debit -- MSL/RAD's actual MEASURED cruise
//       dose-equivalent rate, 1.8 mSv/day (Zeitlin et al. 2013, Science,
//       "Measurements of Energetic Particle Radiation in Transit to Mars on
//       the Mars Science Laboratory"), over the transit duration(s) below --
//       charged against the budget BEFORE any surface time is counted, and
//   (2) the Mars-surface dose-equivalent-rate-vs-regolith-depth curve from
//       this app's own Geant4 campaign (AnalyzeMarsGCRDose.C's
//       LoadWeightedDepthCurve pipeline: FTFP_BERT_HP + QMD/INCLXX ion
//       physics, RAD-surface-anchored, GCR only), re-used here unmodified.
//
// UNITS: regolith depth axis is linear cm (not areal density), same
// convention as AnalyzeMarsGCRDose.C -- see that file's header for the
// 100 g/cm^2 = 65.79 cm derivation.
//
// TRANSIT DURATION IS THE DOMINANT, UNRESOLVED FREE PARAMETER HERE -- READ
// THIS BEFORE TRUSTING ANY CURVE BELOW:
//   1.8 mSv/day is the ONE real, measured transit dose rate available (from
//   MSL's actual 253-day cruise, robotic, uncrewed cruise-stage shielding).
//   There is no measured dose rate for any other transit duration or
//   spacecraft design, so this macro applies that SAME 1.8 mSv/day rate to
//   every transit-duration scenario below -- i.e. it assumes dose rate is
//   set by GCR flux + a fixed (MSL-like) shielding mass, not by trajectory,
//   and ignores any solar-cycle-phase dependence of the ambient GCR flux.
//   Real crewed transit dose rates could differ from 1.8 mSv/day in either
//   direction (more consumables/water mass = more incidental shielding;
//   different trajectory/launch-window solar conditions = different flux).
//   Three ONE-WAY transit durations are scanned to bound this uncertainty --
//   chosen to straddle the break-even point (round-trip transit dose = the
//   ENTIRE 600 mSv budget happens at 600/1.8/2 = 166.7 days one-way, at this
//   fixed dose rate), not just clustered in the "commonly cited" 6-9 month
//   range, because every point in that range turns out to already exceed the
//   budget on transit alone (see KEY FINDING below):
//     130 days -- aggressive/fast-transit class (e.g. higher-thrust
//                 propulsion concepts discussed in NASA mission-architecture
//                 studies), leaves a real surface-dose budget
//     160 days -- just under the 166.7-day break-even point -- a thin but
//                 nonzero surface budget, where shielding depth matters most
//     253 days -- MSL's actual, REAL, measured cruise duration (the only
//                 scenario for which 1.8 mSv/day is not itself an
//                 extrapolation) -- ALREADY OVER BUDGET on transit alone
//   All three assume symmetric outbound/return legs (round-trip transit =
//   2x one-way) -- real missions need not be symmetric, but no asymmetric
//   trajectory data is used or invented here.
//
// KEY FINDING TO WATCH FOR: at the MSL-realistic 253-day-each-way transit
// duration, round-trip transit dose ALONE (2 x 253 x 1.8 mSv/day = 910.8 mSv)
// EXCEEDS the entire 600 mSv career limit before any surface time is counted.
// That means NO amount of Mars-surface regolith shielding can fix the budget
// in that scenario -- the binding constraint is transit shielding/duration,
// not surface shielding. This is a real, direct consequence of the real
// numbers, not a bug -- it's WHY the scenario list above deliberately
// includes faster transit times too, so the plot actually shows a regime
// where regolith depth matters, alongside the realistic regime where it
// currently doesn't.
//
// STATUS: depends on MarsGCRDoseSim's .root output (does not exist until
// that Geant4 campaign is built and run -- see main.cc for compute-cost
// warning). Not run. Paren/brace/bracket balance checked with a script
// before delivery, same as this session's other ROOT macros; logical
// correctness against a real file has not been exercised.
//
// Usage:
//   root -l -q 'CalculateMissionDurationVsShielding.C("MarsGCRDose_QMD.root","MarsGCRDose_INCLXX.root")'

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

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cmath>

using namespace std;

//============================================================================
// Same species relative-abundance weights as AnalyzeMarsGCRDose.C (kept in
// sync manually -- see that file / GCRSpectrum.hh for full sourcing).
//============================================================================
map<int, double> SpeciesRelAbundanceByZ(){
  return {
    {1,  550.0}, {2, 34.0}, {6, 1.20}, {8, 1.00},
    {10, 0.27}, {12, 0.187}, {14, 0.1615}, {26, 0.114},
  };
}
double EnergyWeight(double /*energyMeVPerNucleon*/){ return 1.0; }

const double kRADSurfaceDoseRate_mSvPerDay = 0.64; // Hassler et al. 2014, Science
const double kNASACareerLimit_mSv = 600.0;
const int    kNLayers = 50;
const double kTotalDepthCm = (kNLayers * 2.0) / 1.52; // see AnalyzeMarsGCRDose.C header

// Real MSL/RAD-measured transit dose rate -- see file header caveat on
// applying it across multiple transit-duration scenarios.
const double kTransitDoseRate_mSvPerDay = 1.8; // Zeitlin et al. 2013, Science

struct TransitScenario {
  string label;
  double oneWayDays;
};
vector<TransitScenario> TransitScenarios(){
  return {
    {"130 d one-way (fast-transit)",     130.0},
    {"160 d one-way (near break-even)",  160.0},
    {"253 d one-way (MSL-actual)",       253.0},
  };
}

const double kStandardSurfaceStayDays = 500.0; // conjunction-class orbital-mechanics floor

//============================================================================
// LoadWeightedDepthCurve -- identical logic/branch names to
// AnalyzeMarsGCRDose.C, duplicated here so this macro can run standalone.
//============================================================================
bool LoadWeightedDepthCurve(const string& a_rootFile,
                             vector<double>& a_depthMid_cm,
                             vector<double>& a_doseRate_mSvPerDay){
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
  vector<double> weightSum(kNLayers, 0.0);
  vector<double> depthMid(kNLayers, 0.0);
  vector<bool>   depthSet(kNLayers, false);

  Long64_t nEntries = tree->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    tree->GetEntry(i);
    if(layerIndex_in < 0 || layerIndex_in >= kNLayers) continue;

    auto it = relAbundance.find(Z_in);
    if(it == relAbundance.end()) continue;
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
  vector<double> weighted;
  for(int i = 0; i < kNLayers; i++){
    if(!depthSet[i] || weightSum[i] <= 0.0) continue;
    a_depthMid_cm.push_back(depthMid[i]);
    weighted.push_back(weightedSum[i] / weightSum[i]);
  }
  if(a_depthMid_cm.empty()) return false;

  // RAD anchor: rescale so depth=0 matches the real measured Mars-surface
  // rate, then convert Sv/primary(relative) -> mSv/day.
  double scale_SvPerYear = ((kRADSurfaceDoseRate_mSvPerDay / 1000.0) * 365.25) / weighted.front();
  a_doseRate_mSvPerDay.clear();
  for(double w : weighted){
    double SvPerYear = w * scale_SvPerYear;
    a_doseRate_mSvPerDay.push_back(SvPerYear * 1000.0 / 365.25);
  }
  return true;
}

void CalculateMissionDurationVsShielding(string a_qmdFile = "MarsGCRDose_QMD.root",
                                          string a_inclxxFile = "MarsGCRDose_INCLXX.root",
                                          string a_outputDir = "."){
  vector<double> depthQMD, doseRateQMD, depthINCLXX, doseRateINCLXX;
  bool haveQMD    = LoadWeightedDepthCurve(a_qmdFile, depthQMD, doseRateQMD);
  bool haveINCLXX = LoadWeightedDepthCurve(a_inclxxFile, depthINCLXX, doseRateINCLXX);

  if(!haveQMD && !haveINCLXX){
    cout << "ERROR: neither input file produced a usable curve -- nothing to plot." << endl;
    return;
  }
  // Prefer QMD as the primary curve for this plot (matches AnalyzeMarsGCRDose.C's
  // default-first pattern); fall back to INCLXX alone if QMD isn't available.
  vector<double>& depth = haveQMD ? depthQMD : depthINCLXX;
  vector<double>& doseRateSurface_mSvPerDay = haveQMD ? doseRateQMD : doseRateINCLXX;
  string primaryLabel = haveQMD ? "QMD" : "INCLXX";

  cout << fixed << setprecision(1);
  cout << "\n=== Mission duration vs. regolith shielding depth ===" << endl;
  cout << "Surface dose-rate curve: Geant4 (" << primaryLabel << " ion physics), "
       << depth.size() << " depth points, 0 - " << kTotalDepthCm << " cm" << endl;
  cout << "Transit dose rate applied to ALL scenarios: " << kTransitDoseRate_mSvPerDay
       << " mSv/day (real MSL/RAD measurement, 253-day cruise -- see file header caveat)" << endl;

  vector<TransitScenario> scenarios = TransitScenarios();
  vector<vector<double>> achievableSurfaceDays(scenarios.size());
  vector<double> remainingBudget_mSv(scenarios.size());

  for(size_t s = 0; s < scenarios.size(); s++){
    double roundTripTransitDays = 2.0 * scenarios[s].oneWayDays;
    double transitDose_mSv = roundTripTransitDays * kTransitDoseRate_mSvPerDay;
    remainingBudget_mSv[s] = kNASACareerLimit_mSv - transitDose_mSv;

    cout << "\n-- Scenario: " << scenarios[s].label << " --" << endl;
    cout << "  Round-trip transit: " << roundTripTransitDays << " days, dose = "
         << transitDose_mSv << " mSv" << endl;
    if(remainingBudget_mSv[s] <= 0.0){
      cout << "  >>> Transit dose ALONE exceeds the 600 mSv career limit by "
           << (-remainingBudget_mSv[s]) << " mSv."
           << " No surface time (at ANY regolith depth) fits in the budget for this scenario." << endl;
    } else {
      cout << "  Remaining surface-dose budget after transit: " << remainingBudget_mSv[s] << " mSv" << endl;
    }

    for(size_t i = 0; i < depth.size(); i++){
      double days = (remainingBudget_mSv[s] <= 0.0) ? 0.0
                    : remainingBudget_mSv[s] / doseRateSurface_mSvPerDay[i];
      achievableSurfaceDays[s].push_back(days);
    }

    // Report bare-surface (depth=0) and max-shielding (deepest scanned point)
    // achievable stay, plus whether/where the standard 500-day stay becomes reachable.
    if(!depth.empty()){
      cout << "  Achievable surface stay: " << achievableSurfaceDays[s].front()
           << " days at " << depth.front() << " cm regolith, "
           << achievableSurfaceDays[s].back() << " days at " << depth.back()
           << " cm regolith (deepest scanned)" << endl;
      bool crosses500 = false;
      for(size_t i = 0; i < depth.size(); i++){
        if(achievableSurfaceDays[s][i] >= kStandardSurfaceStayDays){
          cout << "  >>> Standard " << kStandardSurfaceStayDays
               << "-day stay becomes achievable at regolith depth >= " << depth[i]
               << " cm (this curve/anchor/weighting)." << endl;
          crosses500 = true;
          break;
        }
      }
      if(!crosses500 && remainingBudget_mSv[s] > 0.0){
        cout << "  Standard " << kStandardSurfaceStayDays
             << "-day stay is NOT reached even at the deepest scanned depth ("
             << kTotalDepthCm << " cm)." << endl;
      }
    }
  }

  //--------------------------------------------------------------------
  // Plot: achievable surface stay (days) vs. regolith shielding depth (cm),
  // one curve per transit-duration scenario.
  //--------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const char* SURFACE = "#fcfcfb";
  const char* TEXT_PRIMARY = "#0b0b0b";
  const char* TEXT_SECONDARY = "#52514e";
  const char* TEXT_MUTED = "#898781";
  const char* BASELINE = "#c3c2b7";
  const vector<string> seriesColors = {"#2a78d6", "#e34948", "#c8a13d"}; // blue, red, gold

  TCanvas* c = new TCanvas("c_MissionDuration", "Achievable Surface Stay vs Regolith Depth", 950, 750);
  c->SetFillColor(TColor::GetColor(SURFACE));
  c->SetFrameFillColor(TColor::GetColor(SURFACE));
  c->SetLeftMargin(0.13); c->SetRightMargin(0.06);
  c->SetBottomMargin(0.14); c->SetTopMargin(0.12); // only the title lives up here now, no subtitle

  double yMax = kStandardSurfaceStayDays;
  for(auto& v : achievableSurfaceDays) for(double d : v) yMax = max(yMax, d);
  yMax *= 1.2;

  TH1F* frame = c->DrawFrame(0.0, 0.0, kTotalDepthCm, yMax);
  frame->GetXaxis()->SetTitle("Added Mars Regolith Shielding Depth (cm)");
  frame->GetYaxis()->SetTitle("Achievable Surface Stay (days) Within Remaining Career-Dose Budget");
  frame->GetXaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetYaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetXaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetYaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetXaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetYaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetYaxis()->SetTitleSize(0.03);

  TLegend* leg = new TLegend(0.16, 0.72, 0.60, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.024);
  leg->SetTextColor(TColor::GetColor(TEXT_SECONDARY));

  vector<TGraph*> keepAlive;
  for(size_t s = 0; s < scenarios.size(); s++){
    TGraph* g = new TGraph((int) depth.size(), &depth[0], &achievableSurfaceDays[s][0]);
    keepAlive.push_back(g);
    g->SetLineColor(TColor::GetColor(seriesColors[s % seriesColors.size()].c_str()));
    g->SetLineWidth(3);
    if(remainingBudget_mSv[s] <= 0.0) g->SetLineStyle(3); // dotted: budget already blown by transit
    g->Draw("L SAME");
    leg->AddEntry(g, scenarios[s].label.c_str(), "l");
  }

  TLine* stdStayLine = new TLine(0.0, kStandardSurfaceStayDays, kTotalDepthCm, kStandardSurfaceStayDays);
  stdStayLine->SetLineColor(TColor::GetColor("#B8B8B8"));
  stdStayLine->SetLineStyle(2);
  stdStayLine->SetLineWidth(2);
  stdStayLine->Draw();
  TLatex stdStayLabel;
  stdStayLabel.SetTextColor(TColor::GetColor(TEXT_MUTED));
  stdStayLabel.SetTextSize(0.022);
  stdStayLabel.DrawLatex(0.02 * kTotalDepthCm, kStandardSurfaceStayDays + 0.02 * yMax,
                          "Standard conjunction-class stay: 500 days");

  leg->Draw();

  TLatex title; title.SetNDC(); title.SetTextColor(TColor::GetColor(TEXT_PRIMARY));
  title.SetTextFont(62); title.SetTextSize(0.032); title.SetTextAlign(21);
  title.DrawLatex(0.5, 0.955, "Achievable Mars Surface Stay vs. Regolith Depth");

  // No in-plot subtitle -- the method/caveat detail that used to live here
  // belongs in the slide's own figure caption instead.

  c->RedrawAxis();
  string out = a_outputDir + "/MissionDurationVsShieldingDepth.png";
  c->SaveAs(out.c_str());
  cout << "\nWrote " << out << endl;
}
