// CalculateTotalPenetrationDepth.C
//
// ============================================================================
// READ THIS BEFORE PUTTING ANY NUMBER FROM THIS MACRO ON A SLIDE
// ============================================================================
// This macro produces a TOY / APPROXIMATE estimate of the *total* (cumulative,
// multi-interaction) penetration depth of a secondary neutron in Mars
// regolith, by combining two measurements that are NOT the same physical
// process. Be able to explain this distinction if asked.
//
//   1. THE SOLID PART: the mean depth to a neutron's FIRST interaction in
//      regolith, from mars_neutron_attenuation's own Geant4 Monte Carlo
//      (FTFP_BERT_HP), which properly simulates neutron-nucleus hadronic
//      final states. Read from MarsNeutronAttenuation.root.
//
//   2. THE SHAKY PART: an assumed FRACTIONAL ENERGY LOSS PER INTERACTION,
//      derived from published nucleon "rapidity shift" (baryon stopping)
//      measurements in SYMMETRIC-ish HEAVY-ION COLLISIONS. These numbers
//      describe how much rapidity a LEADING BARYON loses over the course of
//      ONE FULL NUCLEUS-NUCLEUS COLLISION -- a different collision system
//      than a single secondary neutron scattering off an individual light
//      regolith nucleus (O, Si, Mg, Fe).
//
//   3. THE <nu>-CORRECTION (partial fix for item 2's granularity mismatch):
//      every measured Delta_y below is the CUMULATIVE shift after a
//      projectile nucleon crosses the ENTIRE target nucleus -- i.e. after
//      <nu> sequential nucleon-nucleon sub-collisions, not after just one.
//      Using a raw measured Delta_y as a per-Geant4-step energy loss (as an
//      earlier version of this macro did) conflates "one full nuclear
//      passage" with "one elementary hadron-nucleus collision." This version
//      instead estimates <nu> with the standard optical-limit Glauber
//      approximation for a head-on (b=0) passage through a uniform-density
//      sphere,
//        <nu> = sigma_NN,inel * T_Au(0),  T_Au(0) = 3*A_Au/(2*pi*R_Au^2),
//        R_Au = r0 * A_Au^(1/3), r0 = 1.2 fm, sigma_NN,inel = 32 mb (SPS/AGS
//        energies) --> <nu> ~ 6.2.
//      This is a simple, transparent, order-of-magnitude estimate -- NOT a
//      published Glauber Monte Carlo number -- applied uniformly to every
//      anchor point below even though they span different systems (Au+Au,
//      Pb+Pb, S+Au) of similar (heavy) target size. It fixes the granularity
//      mismatch; it does NOT fix collective/thermalized-medium effects (no
//      counterpart in dilute single-nucleus scattering), the copious-
//      particle-production-dominated SPS/RHIC energy regime (nothing like a
//      sub-GeV-to-few-GeV regolith neutron), or the fact that Delta_y is an
//      ensemble-averaged final-state observable, not a single-particle
//      energy-loss law.
//
//   4. ENERGY-DEPENDENT Delta_y, FROM REAL DATA AT MULTIPLE ENERGIES: rather
//      than holding Delta_y fixed at its ~200 GeV/nucleon value across the
//      entire 10 MeV-200 GeV chain, this version looks up Delta_y AT THE
//      NEUTRON'S CURRENT ENERGY at every step, interpolated from real
//      measurements from TWO independent sources, both converted to fixed-
//      target lab kinetic energy per nucleon (Delta_y is Lorentz-invariant
//      under boosts along the beam axis, so no ambiguity there; sqrt(s_NN)
//      values are converted via KE = s/(2 m_N) - 2 m_N):
//        AGS E917, Au+Au, most central (Barrette et al., nucl-ex/0003007):
//          6.0 GeV/nucleon:  Delta_y = 0.72 +/- 0.01
//          8.0 GeV/nucleon:  Delta_y = 0.78 +/- 0.01
//          10.8 GeV/nucleon: Delta_y = 0.88 +/- 0.01
//        Gao, Liu, Sun, Sun, Lacey (arXiv:1607.00611) -- three-source
//        Gaussian fits to net-proton rapidity distributions, AGS-RHIC, most
//        central Au+Au/Pb+Pb, converted from their quoted sqrt(s_NN):
//          1.19 GeV/nucleon (sqrt(s_NN)=2.4):  Delta_y = 0.295 +/- 0.042
//          3.25 GeV/nucleon (sqrt(s_NN)=3.1):  Delta_y = 0.612 +/- 0.036
//          5.03 GeV/nucleon (sqrt(s_NN)=3.6):  Delta_y = 0.674 +/- 0.036
//          7.09 GeV/nucleon (sqrt(s_NN)=4.1):  Delta_y = 0.752 +/- 0.066
//          9.90 GeV/nucleon (sqrt(s_NN)=4.7):  Delta_y = 0.966 +/- 0.066
//          11.45 GeV/nucleon (sqrt(s_NN)=5.0): Delta_y = 0.903 +/- 0.060
//          157.6 GeV/nucleon (sqrt(s_NN)=17.3): Delta_y = 1.393 +/- 0.069
//        SPS, S+Au (NA35, CDS 271384):
//          200 GeV/nucleon:  Delta_y = 1.98 +/- 0.05
//      NOTE the genuine ~20% discrepancy at ~158 GeV/nucleon: E917 quotes
//      1.71 for SPS Pb+Pb there, Gao et al.'s own fit to (presumably the
//      same underlying published) data gives 1.393 -- different rapidity-
//      distribution fitting methodology (single vs. three-source Gaussian),
//      not a typo. Rather than silently picking one, the anchor used below
//      is their unweighted average (1.55), which is itself an ad hoc choice
//      -- flagging, not resolving, the spread. Also note the AGS-region
//      points are not perfectly monotonic (Gao et al.'s 9.90 GeV point,
//      0.966, sits ABOVE both its 11.45 GeV neighbor, 0.903, and E917's
//      10.8 GeV point, 0.88) -- real data from two independent
//      analyses, left as-is rather than smoothed over.
//      REMAINING GAP: ~11.5 to 158 GeV/nucleon (over a decade) has NO data
//      point in this table. The actual NA49 SPS energy-scan measurements
//      (net-proton stopping at 20A, 30A, 40A, 80A GeV) exist in the
//      literature and would close this -- e.g. a compiled sqrt(s_NN)-vs-
//      Delta_y figure is shown in arXiv:2512.09630's Fig. 3 -- but the
//      numbers are published only as plotted points there, not as text/
//      table data extractable with the tools available in this pass. Only
//      the two endpoints bounding this gap are used here, log-linearly
//      interpolated. Below the lowest anchor (1.19 GeV/nucleon): ramped
//      linearly to the PHYSICAL constraint Delta_y(0)=0 (a nucleon at rest
//      has no rapidity to lose) -- a much shorter, more defensible
//      extrapolation than the earlier version's ramp from 6 GeV. Above 200
//      GeV/nucleon: flat-clamped (never actually reached -- the neutron
//      grid tops out at 200 GeV). Delta_y -> 0 as energy -> 0 being a
//      built-in physical constraint means the chain naturally stops
//      mattering at low energy without an arbitrary hard cutoff.
//
//      TREAT EVERY NUMBER OUT OF THIS MACRO AS AN ORDER-OF-MAGNITUDE TOY
//      ESTIMATE FOR A MOTIVATION/BACKGROUND SLIDE, NOT A VALIDATED PHYSICS
//      RESULT. If asked in the LFSUPC working group "where does this curve
//      come from", the honest answer is "an admittedly approximate analytic
//      chain model stitching Geant4 output to heavy-ion stopping data at a
//      handful of energies", not "Geant4" -- only the first-interaction
//      curve (dashed, open markers on the plot below) is a Geant4 result.
//
// A MORE RIGOROUS ALTERNATIVE, IF YOU WANT ONE LATER:
//   Track each primary through EVERY interaction (not just the first) inside
//   Geant4 itself, using its own validated neutron-nucleus physics at every
//   step, until the primary is absorbed or drops below a low-energy cutoff --
//   no borrowed heavy-ion numbers needed at all. mars_neutron_attenuation's
//   SteppingAction.cc already detects each interaction; it would just need to
//   stop killing the track at the first one and instead accumulate path
//   length across the full history. Ask if you want that built -- it directly
//   answers "total penetration depth" without the caveats above.
//
// Usage:
//   root -l -q 'CalculateTotalPenetrationDepth.C("MarsNeutronAttenuation.root")'
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace {

const double kNeutronMassMeV = 939.565;

// --- Mean number of sub-collisions per nuclear passage, <nu> --------------
// Optical-limit Glauber estimate for a head-on (b=0) passage of a projectile
// nucleon through a uniform-density gold nucleus. See header item 3 for the
// full derivation and its caveats -- this is a transparent order-of-
// magnitude estimate, not a published Glauber Monte Carlo result.
const double kR0Fm = 1.2;                    // fm, standard hard-sphere nuclear radius parameter
const double kAuMassNumber = 197.0;
const double kSigmaNNInelFm2 = 3.2;          // 32 mb = 3.2 fm^2, SPS-energy (sqrt(s_NN)~19-20 GeV) N-N inelastic xsec
const double kRAuFm = kR0Fm * cbrt(kAuMassNumber); // R_Au = r0 * A_Au^(1/3) ~ 6.984 fm
const double kThicknessAu0Fm2 = 3.0 * kAuMassNumber / (2.0 * M_PI * kRAuFm * kRAuFm); // T_Au(0), fm^-2
const double kMeanNCollisions = kSigmaNNInelFm2 * kThicknessAu0Fm2; // <nu> ~ 6.2

// --- Real Delta_y_total measurements at multiple beam energies, GeV/nucleon
// (kinetic energy, fixed-target lab frame) -- see header item 4 for sources,
// the E917/Gao-et-al. ~158 GeV discrepancy (averaged here), and the
// remaining unfilled ~11.5-158 GeV gap. Ascending order required
// (InterpolateDeltaYTotal assumes it).
const vector<double> kAnchorKE_GeV =
    {1.19, 3.25, 5.03, 6.0, 7.09, 8.0, 9.90, 10.8, 11.45, 158.0, 200.0};
const vector<double> kAnchorDeltaYTotal =
    {0.295, 0.612, 0.674, 0.72, 0.752, 0.78, 0.966, 0.88, 0.903, 1.55, 1.98};
// Representative relative uncertainty applied as a single multiplicative
// band around the whole interpolated curve (NOT a per-point propagation of
// each anchor's own +/-0.01 to +/-0.05 error bar -- that would need a full
// covariance treatment for a toy model that doesn't warrant one). +/-5% is a
// round number a bit larger than NA35's own ~2.5% relative uncertainty, as a
// crude stand-in for the interpolation's additional (unquantified) systematic.
const double kBandRelativeUncertainty = 0.05;

// Interpolated (or, below the lowest anchor / above the highest, ramped-to-
// zero / flat-clamped) TOTAL rapidity loss at a given kinetic energy -- see
// header item 4 for the physical justification of each regime.
double InterpolateDeltaYTotal(double KE_GeV){
  if(KE_GeV <= 0.0) return 0.0;
  if(KE_GeV <= kAnchorKE_GeV.front()){
    // No data below the lowest anchor -- linear ramp to the physical
    // constraint Delta_y(0)=0 (can't lose more rapidity than you have).
    return kAnchorDeltaYTotal.front() * (KE_GeV / kAnchorKE_GeV.front());
  }
  if(KE_GeV >= kAnchorKE_GeV.back()) return kAnchorDeltaYTotal.back(); // flat-clamp; grid never actually exceeds this

  size_t hi = (size_t)(upper_bound(kAnchorKE_GeV.begin(), kAnchorKE_GeV.end(), KE_GeV) - kAnchorKE_GeV.begin());
  size_t lo = hi - 1;
  double lx0 = log(kAnchorKE_GeV[lo]), lx1 = log(kAnchorKE_GeV[hi]);
  double t = (log(KE_GeV) - lx0) / (lx1 - lx0); // linear in ln(KE)
  return kAnchorDeltaYTotal[lo] + t * (kAnchorDeltaYTotal[hi] - kAnchorDeltaYTotal[lo]); // linear in Delta_y
}

// Log-log linear interpolation of the Geant4 first-interaction-depth table;
// flat-clamped outside the table's energy range (never extrapolated).
double InterpolateLogLog(const vector<double>& xs, const vector<double>& ys, double x){
  if(xs.empty()) return -1.0;
  if(x <= xs.front()) return ys.front();
  if(x >= xs.back())  return ys.back();
  size_t hi = (size_t)(upper_bound(xs.begin(), xs.end(), x) - xs.begin());
  size_t lo = hi - 1;
  double lx0 = log(xs[lo]), lx1 = log(xs[hi]);
  double ly0 = log(ys[lo]), ly1 = log(ys[hi]);
  double t = (log(x) - lx0) / (lx1 - lx0);
  return exp(ly0 + t * (ly1 - ly0));
}

// Chains interactions from E0_MeV (kinetic energy) down to the Geant4 data
// floor (xs.front(), ~10 MeV), re-evaluating Delta_y AT THE NEUTRON'S CURRENT
// ENERGY at every step (not just its starting energy) via
// InterpolateDeltaYTotal(), divided by <nu> for the per-sub-collision value.
// No separate low-energy cutoff is needed: Delta_y_total -> 0 as KE -> 0 is a
// built-in physical constraint, so steps taken near the floor are naturally
// tiny and the chain tapers off on its own. bandScale multiplies the whole
// interpolated Delta_y_total curve (1.0 central; 1+/-kBandRelativeUncertainty
// for the band -- see header item 4). Returns {totalDepth_cm, nInteractions}.
pair<double,int> TotalDepth(double E0_MeV, double bandScale,
                             const vector<double>& xs, const vector<double>& ys){
  double KE = E0_MeV;
  double depth = 0.0;
  int nInt = 0;
  const int kMaxIter = 200000;
  const double floorMeV = xs.empty() ? 0.0 : xs.front();

  // ">=" (not ">"): at E0 exactly equal to the Geant4 data floor -- which is
  // true for the lowest point on the energy grid, since that grid point IS
  // the floor -- a strict ">" skipped the loop body entirely and returned
  // depth=0. That fake zero is what was dragging the total-depth curve's
  // leftmost point down to the bottom of the (log-scale) plot instead of to
  // its real value, which is just the first-interaction depth (the chain
  // takes exactly one step at the floor energy, then the neutron's energy
  // drops below the floor and there is no more Geant4 data to chain from).
  while(KE >= floorMeV && nInt < kMaxIter){
    double lambda = InterpolateLogLog(xs, ys, KE);
    if(lambda <= 0) break;
    depth += lambda;
    nInt++;

    double deltaYTotal = InterpolateDeltaYTotal(KE / 1000.0) * bandScale; // MeV -> GeV
    double deltaYStep = deltaYTotal / kMeanNCollisions;
    double Etot = (KE + kNeutronMassMeV) * exp(-deltaYStep);
    KE = Etot - kNeutronMassMeV;
  }

  return make_pair(depth, nInt);
}

} // namespace

void CalculateTotalPenetrationDepth(string a_rootFile = "MarsNeutronAttenuation.root",
                                     string a_outputDir = "."){
  cout << "============================================================\n"
       << "CalculateTotalPenetrationDepth.C -- TOY MODEL. See file header\n"
       << "before trusting or presenting these numbers as physics.\n"
       << "Delta_y_total(KE) anchors [GeV/nucleon : Delta_y], AGS E917 +\n"
       << "Gao et al. (arXiv:1607.00611) + NA35:\n"
       << "  1.19:0.295  3.25:0.612  5.03:0.674  6.0:0.72  7.09:0.752\n"
       << "  8.0:0.78  9.90:0.966  10.8:0.88  11.45:0.903\n"
       << "  158:1.55(avg of 1.71,1.39 -- see header)  200:1.98\n"
       << "  (below 1.19 GeV/nucleon: linear ramp to Delta_y(0)=0, unmeasured;\n"
       << "   ~11.5-158 GeV/nucleon: no data, log-linear interpolation only)\n"
       << "Glauber-estimated <nu> (sub-collisions/nuclear passage) = "
       << kMeanNCollisions << "  [R_Au=" << kRAuFm << " fm, T_Au(0)="
       << kThicknessAu0Fm2 << " fm^-2, sigma_NN,inel=" << kSigmaNNInelFm2 << " fm^2]\n"
       << "Per-sub-collision Delta_y = Delta_y_total(current KE) / <nu>,\n"
       << "re-evaluated at EVERY step using the neutron's CURRENT energy.\n"
       << "Band: Delta_y_total(KE) scaled by 1 +/- " << kBandRelativeUncertainty
       << " (representative, not a rigorous per-point error propagation).\n"
       << "============================================================\n";

  TFile* f = TFile::Open(a_rootFile.c_str(), "READ");
  if(!f || f->IsZombie()){ cout << "ERROR: could not open " << a_rootFile << endl; return; }
  TTree* tree = (TTree*) f->Get("NeutronAttenuation");
  if(!tree){
    cout << "ERROR: TTree \"NeutronAttenuation\" not found in " << a_rootFile << endl;
    f->Close();
    return;
  }

  // stdErrOfMean_cm is the same per-energy-point statistical uncertainty
  // (standard error of the mean over that point's Geant4 primaries) that
  // PlotNeutronAttenuationLength.C already plots -- read here too so the
  // rigorous first-interaction curve carries real error bars instead of
  // implying (by omission) that it has none, next to the toy model's
  // shaded systematic band.
  Double_t energyMeV_in, meanDepth_in, stdErr_in;
  tree->SetBranchAddress("energyMeV", &energyMeV_in);
  tree->SetBranchAddress("meanFirstInteractionDepth_cm", &meanDepth_in);
  tree->SetBranchAddress("stdErrOfMean_cm", &stdErr_in);

  vector<double> xs, ys, yErrs;
  Long64_t nEntries = tree->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    tree->GetEntry(i);
    xs.push_back(energyMeV_in);
    ys.push_back(meanDepth_in);
    yErrs.push_back(stdErr_in);
  }
  f->Close();

  if(xs.empty()){ cout << "ERROR: no entries read -- nothing to compute." << endl; return; }

  // Sort ascending defensively (main.cc already loops the grid in ascending
  // order, but don't assume it for an arbitrary input file).
  vector<size_t> order(xs.size());
  for(size_t i = 0; i < order.size(); i++) order[i] = i;
  sort(order.begin(), order.end(), [&](size_t a, size_t b){ return xs[a] < xs[b]; });
  vector<double> sxs, sys, syErrs;
  for(size_t idx : order){ sxs.push_back(xs[idx]); sys.push_back(ys[idx]); syErrs.push_back(yErrs[idx]); }

  vector<double> energy, totalCentral, totalLow, totalHigh, firstOnly, firstErr;
  vector<int> nIntCentral;
  for(size_t i = 0; i < sxs.size(); i++){
    double e0 = sxs[i];
    pair<double,int> c  = TotalDepth(e0, 1.0, sxs, sys);
    pair<double,int> lo = TotalDepth(e0, 1.0 - kBandRelativeUncertainty, sxs, sys);
    pair<double,int> hi = TotalDepth(e0, 1.0 + kBandRelativeUncertainty, sxs, sys);
    energy.push_back(e0);
    totalCentral.push_back(c.first);
    totalLow.push_back(min(lo.first, hi.first));
    totalHigh.push_back(max(lo.first, hi.first));
    nIntCentral.push_back(c.second);
    // e0 is drawn from sxs itself, so this interpolation lands exactly on
    // a grid point (t=0 in InterpolateLogLog) and firstOnly[i] == sys[i]
    // exactly -- so syErrs[i] is the right error to pair with it, no
    // separate interpolation of the uncertainty needed.
    firstOnly.push_back(InterpolateLogLog(sxs, sys, e0));
    firstErr.push_back(syErrs[i]);
  }

  cout << "\nenergyMeV\tfirstInteraction_cm (+/- stdErr)\ttotalDepth_cm [band]\tnInteractions\n";
  for(size_t i = 0; i < energy.size(); i++){
    cout << energy[i] << "\t" << firstOnly[i] << " +/- " << firstErr[i] << "\t"
         << totalCentral[i] << " [" << totalLow[i] << ", " << totalHigh[i] << "]"
         << "\t" << nIntCentral[i] << "\n";
  }

  // --- Plot ---
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const char* SURFACE = "#fcfcfb";
  const char* TEXT_PRIMARY = "#0b0b0b";
  const char* TEXT_SECONDARY = "#52514e";
  const char* TEXT_MUTED = "#898781";
  const char* BASELINE = "#c3c2b7";
  const char* SERIES_FIRST = "#2a78d6"; // blue -- rigorous Geant4 first-interaction result
  const char* SERIES_TOTAL = "#c8553d"; // warm red-orange -- toy total-depth chain model

  TCanvas* c = new TCanvas("c_TotalDepth", "Total penetration depth (toy model)", 1000, 780);
  c->SetLogx();
  c->SetLogy();
  c->SetFillColor(TColor::GetColor(SURFACE));
  c->SetFrameFillColor(TColor::GetColor(SURFACE));
  c->SetLeftMargin(0.13); c->SetRightMargin(0.06);
  c->SetBottomMargin(0.16); c->SetTopMargin(0.12); // only the title lives up here now, no subtitle

  double xMin = energy.front() * 0.7, xMax = energy.back() * 1.4;
  double yMinData = 1e9, yMax = -1e9;
  for(size_t i = 0; i < energy.size(); i++){
    yMinData = min(yMinData, min(firstOnly[i], totalLow[i]));
    yMax = max(yMax, totalHigh[i]);
  }
  yMax = yMax * 1.5;
  // Reserve ~1.3 extra decades below the lowest data point (rather than
  // hugging it at yMinData*0.6) so there is a guaranteed-empty band across
  // the FULL x-range -- both curves' minimum is set by the Geant4
  // first-interaction depth, which never goes anywhere near this low -- to
  // anchor the in-plot legend and reference-line entries without having to
  // predict exactly where the curves fall at every x.
  double yMin = max(yMinData / 20.0, 0.1);

  TH1F* frame = c->DrawFrame(xMin, yMin, xMax, yMax);
  frame->GetXaxis()->SetTitle("Initial Neutron Kinetic Energy (MeV)");
  frame->GetYaxis()->SetTitle("Penetration Depth (cm)");
  frame->GetXaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetYaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetXaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetYaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetXaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetYaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetXaxis()->SetTitleOffset(1.35); // push the x title down, away from the tick-number row
  frame->GetXaxis()->SetLabelOffset(0.006);

  // Uncertainty band (from the +/-0.05 on Delta_y) for the toy total-depth curve.
  int n = (int) energy.size();
  TGraph* band = new TGraph(2 * n);
  for(int i = 0; i < n; i++) band->SetPoint(i, energy[i], totalHigh[i]);
  for(int i = 0; i < n; i++) band->SetPoint(n + i, energy[n - 1 - i], totalLow[n - 1 - i]);
  band->SetFillColorAlpha(TColor::GetColor(SERIES_TOTAL), 0.18);
  band->SetLineWidth(0);
  band->Draw("F SAME");

  TGraph* gTotal = new TGraph(n, &energy[0], &totalCentral[0]);
  gTotal->SetLineColor(TColor::GetColor(SERIES_TOTAL));
  gTotal->SetMarkerColor(TColor::GetColor(SERIES_TOTAL));
  gTotal->SetMarkerStyle(20); gTotal->SetMarkerSize(0.8); gTotal->SetLineWidth(2);
  gTotal->Draw("LP SAME");

  // TGraphErrors, not TGraph: firstErr is the real per-point statistical
  // uncertainty (standard error of the mean) from the Geant4 run, same
  // quantity PlotNeutronAttenuationLength.C already shows. No x errors.
  TGraphErrors* gFirst = new TGraphErrors(n, &energy[0], &firstOnly[0], nullptr, &firstErr[0]);
  gFirst->SetLineColor(TColor::GetColor(SERIES_FIRST));
  gFirst->SetMarkerColor(TColor::GetColor(SERIES_FIRST));
  gFirst->SetMarkerStyle(24); gFirst->SetMarkerSize(0.7);
  gFirst->SetLineWidth(2); gFirst->SetLineStyle(2);
  gFirst->Draw("LP SAME");

  // Marks the lowest real Delta_y data point (Gao et al., ~1.19 GeV/nucleon)
  // -- below this line, Delta_y_total is an unmeasured linear ramp to zero.
  double lowestAnchorMeV = kAnchorKE_GeV.front() * 1000.0;
  TLine* floorLine = new TLine(lowestAnchorMeV, yMin, lowestAnchorMeV, yMax);
  floorLine->SetLineColor(TColor::GetColor(TEXT_MUTED));
  floorLine->SetLineStyle(3);
  floorLine->Draw("SAME");

  // Marks the highest real Delta_y data point (SPS/NA35, 200 GeV/nucleon) --
  // above this line Delta_y_total is flat-clamped at that last measured
  // value. In practice this sits right at (or just past) the neutron energy
  // grid's own top edge, since the grid tops out at 200 GeV too -- see file
  // header item 4 -- so it will land close to the plot's last data points
  // rather than out in open space the way the floor line does.
  double highestAnchorMeV = kAnchorKE_GeV.back() * 1000.0;
  TLine* ceilingLine = new TLine(highestAnchorMeV, yMin, highestAnchorMeV, yMax);
  ceilingLine->SetLineColor(TColor::GetColor(TEXT_MUTED));
  ceilingLine->SetLineStyle(3);
  ceilingLine->Draw("SAME");

  // Bottom-right corner is empty in both curves across their full energy
  // range (first-interaction depth plateaus low there, total-depth model is
  // only large at HIGH y, and neither curve dips back down at high x), so
  // the legend no longer sits on top of either series the way the old
  // top-left placement did.
  // Text size kept small and margin tight -- TLegend does not wrap or clip
  // its entry labels, so at the default text size/margin these ran past the
  // right edge of the legend box and off the edge of the frame. Measured
  // against LiberationSans (metric-compatible with ROOT's default Helvetica,
  // font 42) at this canvas width, both labels below now end well short of
  // the frame's right border (NDC x=0.94).
  TLegend* leg = new TLegend(0.60, 0.185, 0.90, 0.345);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.018);
  leg->SetMargin(0.12);
  leg->AddEntry(gFirst, "First-interaction depth (Geant4)", "lp");
  leg->AddEntry(gTotal, "Total depth -- toy stopping-chain model", "lp");
  // Lowercase delta ("#delta", not "#Delta"): the rapidity-loss literature
  // this session has been citing (Videbaek, BRAHMS, Gao et al.) uniformly
  // writes this quantity as lowercase deltay -- capital Delta was wrong here.
  leg->AddEntry(floorLine, "Lowest measured #delta y (1.19 GeV/A)", "l");
  leg->AddEntry(ceilingLine, "Highest measured #delta y (200 GeV/A)", "l");
  leg->Draw();

  // No in-plot subtitle -- the method/caveat detail that used to live here
  // belongs in the slide's own figure caption instead, not baked into the
  // image (and shouldn't reference this macro, since the audience isn't
  // looking at the code).
  TLatex title; title.SetNDC(); title.SetTextColor(TColor::GetColor(TEXT_PRIMARY));
  title.SetTextFont(62); title.SetTextSize(0.030); title.SetTextAlign(21);
  title.DrawLatex(0.5, 0.955, "Total Neutron Penetration Depth in Mars Regolith (Toy Stopping-Chain Model)");

  c->RedrawAxis();
  string out = a_outputDir + "/TotalPenetrationDepthToyModel.png";
  c->SaveAs(out.c_str());
  cout << "\nWrote " << out << endl;
}
