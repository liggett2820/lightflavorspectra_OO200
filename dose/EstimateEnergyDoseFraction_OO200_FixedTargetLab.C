// EstimateEnergyDoseFraction_OO200_FixedTargetLab.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS (2026-08-31): follow-up to
// EstimateEnergyDoseFraction_OO200_LiteratureStopping.C, after a real physics
// gap surfaced in conversation: that macro's energy budget is built entirely
// in the NN CENTER-OF-MASS/COLLIDER frame (E_beam=sqrt(s_NN)/2=100 GeV given
// to EVERY nucleon, both nuclei symmetrically) -- correct and self-consistent
// for reproducing STAR/BRAHMS's own numbers, but NOT the frame in which
// "dose deposited in a target" (a GCR ion striking spacecraft shielding or
// tissue, which is at rest) is physically defined. Andrew confirmed he wants
// THIS version: a genuine fixed-target/lab-frame energy budget, with the
// stopping loss rebuilt from the frame-invariant rapidity shift rather than
// reusing the collider-frame 73 GeV figure unchanged.
//
// THE KEY PHYSICS POINT: energy is NOT frame-invariant, but the rapidity
// SHIFT a stopped participant undergoes (delta_y = 2.05, BRAHMS/Videbaek, see
// EstimateBaryonStoppingFromLiterature_AuAu200.C) IS frame-invariant --
// that's the actual portable physics. Because E = m_p*cosh(y) is a NONLINEAR
// function of y, applying the SAME delta_y at a different starting rapidity
// (the fixed-target lab frame's much higher beam rapidity, vs. the collider
// frame's) gives a DIFFERENT fractional energy loss, not just a different
// absolute one. Worked out below: ~87% in this frame, not 73%.
//
// EXACT KINEMATICS (fixed-target, one nucleon at rest, verified numerically
// -- see the python check this macro's numbers were validated against,
// referenced in the session that wrote this):
//   Projectile lab energy needed to reproduce sqrt(s_NN)=SQRT_S_NN_GEV
//   (matching the SAME STAR O+O collision energy the rest of this repo uses,
//   just reinterpreted in the fixed-target frame instead of the collider
//   frame) from s = 2*m_p*(m_p + E1_lab):
//     E1_lab = sqrt_s^2/(2*m_p) - m_p  =  21321.0 GeV  ~= 21.32 TeV
//   Beam rapidity in this frame (exact, y=0.5*ln((E+p)/(E-p))):
//     y_beam,lab = 10.7246   (vs. y_beam,CM = 5.3623 = ln(sqrt_s/m_p), the
//     value EstimateBaryonStoppingFromLiterature_AuAu200.C already uses)
//   A PROJECTILE participant loses the same delta_y=2.05 off ITS OWN beam
//   rapidity: y0_proj,lab = y_beam,lab - delta_y = 8.6746
//     -> T_proj,after = m_p*cosh(y0_proj,lab) - m_p = 2743.8 GeV
//     -> fractional loss = (T1_lab - T_proj,after)/T1_lab = 87.13%
//   A TARGET participant starts AT REST (y=0 in this frame) -- mirroring the
//   same symmetric stopping picture through the CM frame (a participant on
//   either side moves delta_y toward the other side's rapidity) gives, after
//   the algebra, the clean closed form:
//     y0_targ,lab = delta_y   (exactly -- the boost rapidity from lab to CM
//     equals y_beam,CM in this construction, which cancels the y_beam,CM
//     terms in the mirror formula; verified numerically, not asserted)
//     -> T_targ,after = m_p*(cosh(delta_y) - 1) = 2.7655 GeV
//   The target participant's kinetic-energy PICKUP (2.77 GeV, from being
//   struck while at rest) is what a real fixed-target dose calculation
//   should show as recoil energy -- and it is utterly negligible next to the
//   projectile's own ~21.3 TeV/18.6 TeV scale (0.013%), which is exactly
//   what "the target barely moves compared to the beam" should look like.
//   Kept in the accounting anyway (not dropped) since it costs nothing and
//   makes the energy budget exactly closed rather than approximately closed.
//
// ASYMMETRIC MODEL (the real structural difference from the collider-frame
// macro, not just different numbers): a fixed-target event is NOT symmetric
// between the two nuclei the way a collider event is. Only the PROJECTILE's
// 16 nucleons carry real kinetic energy going in (T1_lab each); the TARGET's
// 16 nucleons start at rest (T=0). Npart/Nspec themselves are unchanged from
// the existing Glauber table -- participant/spectator COUNTING is a
// transverse-geometry question (nuclear density profiles, impact parameter,
// inelastic cross section), not a kinematic one, so it does not depend on
// which frame the ENERGY accounting is done in. What changes is how that
// same Npart/Nspec is split and weighted: for a symmetric system (O+O, same
// species both sides), Npart splits evenly Npart_proj=Npart_targ=Npart/2,
// same for Nspec -- and only the PROJECTILE half of Nspec carries meaningful
// energy (target spectators are just sitting there, T=0, same as before the
// collision).
//
// ENERGY BUDGET PER CENTRALITY BIN (all kinetic energy, T = E - m_p, so
// "at rest" = 0 exactly rather than m_p -- rest mass is not what's being
// dosed and would otherwise swamp the target-side numbers, which are only a
// few GeV):
//   E_spectator[i]   = (Nspec[i]/2) * T1_lab
//     (target spectators contribute 0; only projectile spectators carry KE)
//   E_participant[i] = (Npart[i]/2) * (T_proj_after + T_targ_after)
//     (both sides' participants, after stopping)
//   E_produced[i]    = (Npart[i]/2) * (T1_lab - T_proj_after - T_targ_after)
//     (per-projectile-participant net energy converted to new particles,
//     i.e. the projectile's own loss minus what it handed directly to its
//     target partner as recoil KE -- derived from total energy conservation,
//     see the derivation this macro's header comment traces)
//   E_total[i] = 16 * T1_lab, CONSTANT across every bin (Npart[i]+Nspec[i]=32
//   always) -- same role as the old "3200 GeV" check, now ~341.1 TeV.
//
// WHAT DID NOT CHANGE FROM THE COLLIDER-FRAME MACRO:
//   - Npart/Nspec per bin: same Glauber-fit-calibrated numbers, same reader.
//   - delta_y = 2.05 +/- 0.17 (BRAHMS/Videbaek): frame-invariant, reused as-is.
//   - The cross-section-weighted min-bias average across the six centrality
//     bins (CENT_EDGES_FRAC weighting): this is a statement about how often
//     each centrality bin occurs, unrelated to which energy frame the
//     per-bin fractions are computed in, so the same weighting scheme
//     applies unchanged.
//   - The two literature caveats from EstimateBaryonStoppingFromLiterature_
//     AuAu200.C: delta_y=2.05 is an Au+Au (not O+O) measurement, and only
//     the 0-5% bin has a direct literature basis (the rest assume
//     centrality-independence). Both still apply here -- this macro changes
//     which FRAME the energy accounting is done in, not which underlying
//     physics measurement it's built from.
//
// NEW CAVEAT SPECIFIC TO THIS VERSION: this still uses the SAME sqrt(s_NN)=
// 200 GeV physics (Npart/Nspec, delta_y) as the collider-frame macro --
// it is a kinematic REINTERPRETATION of that same STAR-energy collision in a
// different frame, not an independent measurement at a different energy
// regime. It answers "what does a sqrt(s_NN)=200 GeV collision's energy
// split look like if you insist on viewing it in the frame where the target
// is at rest" -- which is what "GeV per collision" MEANS for a real
// fixed-target GCR-on-shielding event. It is explicitly NOT a claim about
// what fraction of a real astronaut's total mission dose comes from this
// energy (this is one narrow energy slice; see the flux-weighting caveat
// below), and NOT a claim that stopping itself works identically at 21 TeV
// projectile energy vs. published measurements (delta_y's own energy
// dependence beyond matching sqrt(s_NN) is not modeled -- it is taken as the
// same fixed number either way, exactly as the collider-frame macro already
// assumed centrality-independence of the same input).
//
// NOT ATTEMPTED HERE: folding in the actual GCR differential flux at the
// ~21.3 TeV/nucleon equivalent energy to get this energy slice's share of
// TOTAL accumulated dose across the real GCR spectrum. That's flux-weighting
// on top of this per-collision energy split, a separate calculation Andrew
// has not asked for yet (see prior discussion this macro followed).
//
// OUTPUT: console table (same layout as the collider-frame macro, TeV/GeV
// units called out explicitly since they now differ by side) and
// OO200_EnergyDoseFractionFixedTargetLab.png -- same stacked-bar style/
// palette as the collider-frame macro's PNG, retitled, for direct visual
// comparison between the two frames.
//
// USAGE:
//   root -l -b -q 'EstimateEnergyDoseFraction_OO200_FixedTargetLab.C()'
// or with an explicit Glauber ROOT file path:
//   root -l -b -q 'EstimateEnergyDoseFraction_OO200_FixedTargetLab.C("glauber_oo200_results_NBDfit.root")'
//
// STATUS: UNTESTED -- no ROOT in the sandbox that wrote this (same caveat as
// every other C++ macro in this session). The underlying kinematics were
// checked numerically with a standalone script before being written into
// this macro's constants/formulas, but the macro itself, as ROOT/cling C++,
// has not been run.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TStyle.h"
#include "TSystem.h"

using namespace std;

//============================================================================
// Palette -- identical to EstimateEnergyDoseFraction_OO200_LiteratureStopping.C
// for direct visual comparison between the two frames.
//============================================================================
const int COLOR_FIREBALL      = TColor::GetColor("#C0392B"); // brick red
const int COLOR_PARTICIPANT   = TColor::GetColor("#1E2761"); // navy
const int COLOR_SPECTATOR     = TColor::GetColor("#E8A33D"); // amber
const int COLOR_TEXT_ON_DARK  = kWhite;
const int COLOR_TEXT_ON_LIGHT = TColor::GetColor("#1A1A1A");
const int COLOR_CAVEAT        = TColor::GetColor("#7A7A7A");

const int N_CENT_BINS = 6;
const char* CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
const double CENT_EDGES_FRAC[N_CENT_BINS + 1] = {0.0, 0.05, 0.10, 0.20, 0.40, 0.80, 1.00};

// Cited inputs -- see header. Same sqrt(s_NN) and delta_y as the rest of the
// repo; this macro's whole point is REINTERPRETING them in a different frame,
// not changing them.
const double M_NUCLEON_GEV   = 0.938;  // same value used throughout dose/
const double SQRT_S_NN_GEV   = 200.0;  // same STAR O+O collision energy
const double DELTA_Y         = 2.05;   // BRAHMS/Videbaek, frame-invariant
const double DELTA_Y_ERR     = 0.17;

//============================================================================
// Exact fixed-target kinematics -- see header derivation. Computed here, not
// hardcoded, so it stays correct if SQRT_S_NN_GEV/M_NUCLEON_GEV are ever
// revised.
//============================================================================
double ProjectileLabEnergy(double sqrtS, double m){
  // s = 2*m*(m+E1_lab)  =>  E1_lab = s/(2m) - m
  return sqrtS*sqrtS/(2.0*m) - m;
}

double BeamRapidityLab(double E1_lab, double m){
  double p1_lab = sqrt(E1_lab*E1_lab - m*m);
  return 0.5*log((E1_lab+p1_lab)/(E1_lab-p1_lab));
}

double BeamRapidityCM(double sqrtS, double m){
  return acosh(sqrtS/(2.0*m)); // exact; matches ln(sqrtS/m) to high precision here
}

// Given delta_y, returns the projectile participant's and target participant's
// post-stopping KINETIC energy (T = E - m_p) in this lab frame. See header
// for the derivation (mirrored stopping through the CM frame for the target
// side; y0_targ_lab reduces to exactly delta_y in this construction).
void StoppedParticipantEnergies(double dy, double yBeamLab, double yBeamCM, double m,
                                 double& tProjAfter, double& tTargAfter){
  double y0ProjLab = yBeamLab - dy;
  tProjAfter = m*cosh(y0ProjLab) - m;

  double yBoost = yBeamLab - yBeamCM;
  double y0TargLab = -(yBeamCM - dy) + yBoost; // = dy exactly when yBoost==yBeamCM (verified numerically)
  tTargAfter = m*cosh(y0TargLab) - m;
}

//============================================================================
// Glauber centrality-bin reader -- identical to the collider-frame macro's
// (same TTree schema, "GlauberCentralityBins").
//============================================================================
struct GlauberBinRow {
  string cent;
  double meanNpart = 0.0, meanNspec = 0.0;
  bool foundNpart = false, foundNspec = false;
};

vector<GlauberBinRow> ReadGlauberCentralityTable(const string& a_path){
  vector<GlauberBinRow> out;
  TFile* f = TFile::Open(a_path.c_str(), "READ");
  if(!f || f->IsZombie()) return out;

  TTree* t = (TTree*) f->Get("GlauberCentralityBins");
  if(!t){ f->Close(); return out; }

  Char_t cent_in[16];
  Double_t meanNpart_in, meanNspec_in;
  t->SetBranchAddress("cent", cent_in);
  t->SetBranchAddress("meanNpart", &meanNpart_in);
  t->SetBranchAddress("meanNspec", &meanNspec_in);

  Long64_t nEntries = t->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    GlauberBinRow row;
    row.cent = string(cent_in);
    row.meanNpart = meanNpart_in;
    row.meanNspec = meanNspec_in;
    row.foundNpart = true;
    row.foundNspec = true;
    out.push_back(row);
  }
  f->Close();
  return out;
}

void DrawPercentLabel(double x, double y, double val, int color, double textSize = 0.028){
  if(val < 3.0) return;
  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextSize(textSize);
  lat.SetTextAlign(22);
  lat.SetTextColor(color);
  lat.DrawLatex(x, y, Form("%.0f%%", val));
}

//============================================================================
// Main macro
//============================================================================
void EstimateEnergyDoseFraction_OO200_FixedTargetLab(
    string a_glauberRoot = "glauber_oo200_results_NBDfit.root",
    string a_outputDir   = "."
){
  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outputDir.c_str(), kTRUE);

  //--------------------------------------------------------------------------
  // Exact fixed-target kinematics for this collision energy.
  //--------------------------------------------------------------------------
  double E1_lab = ProjectileLabEnergy(SQRT_S_NN_GEV, M_NUCLEON_GEV);
  double T1_lab = E1_lab - M_NUCLEON_GEV;
  double yBeamLab = BeamRapidityLab(E1_lab, M_NUCLEON_GEV);
  double yBeamCM  = BeamRapidityCM(SQRT_S_NN_GEV, M_NUCLEON_GEV);

  double tProjAfter, tTargAfter;
  StoppedParticipantEnergies(DELTA_Y, yBeamLab, yBeamCM, M_NUCLEON_GEV, tProjAfter, tTargAfter);

  cout << "==================================================================" << endl;
  cout << "EstimateEnergyDoseFraction_OO200_FixedTargetLab" << endl;
  cout << "O+O sqrt(s_NN)=" << SQRT_S_NN_GEV << " GeV, reinterpreted in the FIXED-TARGET LAB FRAME" << endl;
  cout << "(target at rest -- the frame in which dose is actually deposited/measured)" << endl;
  cout << "==================================================================" << endl;
  cout << "Projectile lab energy T1_lab = " << T1_lab << " GeV = " << T1_lab/1000.0 << " TeV" << endl;
  cout << "y_beam,lab = " << yBeamLab << "   y_beam,CM = " << yBeamCM
       << "   (collider-frame macro uses y_beam,CM=100 GeV/nucleon directly)" << endl;
  cout << "delta_y = " << DELTA_Y << " +/- " << DELTA_Y_ERR << " (BRAHMS/Videbaek, frame-invariant, reused as-is)" << endl;
  cout << "  -> projectile participant retains " << tProjAfter << " GeV = " << tProjAfter/1000.0
       << " TeV  (" << 100.0*(T1_lab-tProjAfter)/T1_lab << "% lost -- vs. 73% in the collider-frame macro)" << endl;
  cout << "  -> target participant gains " << tTargAfter << " GeV from rest (negligible next to the TeV scale above)" << endl;
  cout << "==================================================================" << endl;

  //--------------------------------------------------------------------------
  // Load Npart/Nspec per bin -- required, no fabrication if missing.
  //--------------------------------------------------------------------------
  vector<GlauberBinRow> rows = ReadGlauberCentralityTable(a_glauberRoot);
  if(rows.empty()){
    cerr << "ERROR: could not read any centrality-bin rows from '" << a_glauberRoot << "'."
         << " Run FitGlauberNBDToRefMult_OO200.C (or GenerateGlauberMC_OO200.C for the purely"
         << " geometric file) first, or pass the right path." << endl;
    return;
  }

  double meanNpart[N_CENT_BINS], meanNspec[N_CENT_BINS];
  for(int i = 0; i < N_CENT_BINS; i++){
    const GlauberBinRow* match = NULL;
    for(size_t j = 0; j < rows.size(); j++){
      if(rows[j].cent == CENT_LABELS[i]){ match = &rows[j]; break; }
    }
    if(!match || !match->foundNpart || !match->foundNspec){
      cerr << "ERROR: '" << a_glauberRoot << "' has no mean_Npart/mean_Nspec for cent bin '"
           << CENT_LABELS[i] << "' -- aborting rather than guessing." << endl;
      return;
    }
    meanNpart[i] = match->meanNpart;
    meanNspec[i] = match->meanNspec;
  }

  //--------------------------------------------------------------------------
  // Energy accounting per bin -- ASYMMETRIC projectile/target split, see
  // header. Npart/Nspec are split evenly between the two identical (O+O)
  // nuclei; only the projectile side carries meaningful kinetic energy.
  //--------------------------------------------------------------------------
  double eSpectator[N_CENT_BINS], eParticipant[N_CENT_BINS], eProduced[N_CENT_BINS], eTotal[N_CENT_BINS];
  double fireballFrac[N_CENT_BINS], participantFrac[N_CENT_BINS], spectatorFrac[N_CENT_BINS];

  printf("%-10s%10s%10s%16s%16s%16s%14s\n", "Cent", "Npart", "Nspec", "E_spec[TeV]", "E_part[TeV]", "E_prod[TeV]", "E_tot[TeV]");
  for(int i = 0; i < N_CENT_BINS; i++){
    double npartProj = meanNpart[i]/2.0, npartTarg = meanNpart[i]/2.0;
    double nspecProj = meanNspec[i]/2.0; // target spectators carry ~0, dropped from the energy sum

    eSpectator[i]   = nspecProj * T1_lab;
    eParticipant[i] = npartProj*tProjAfter + npartTarg*tTargAfter;
    eProduced[i]    = npartProj*T1_lab - npartProj*tProjAfter - npartTarg*tTargAfter;
    eTotal[i]       = eSpectator[i] + eParticipant[i] + eProduced[i];

    fireballFrac[i]    = eProduced[i]    / eTotal[i];
    participantFrac[i] = eParticipant[i] / eTotal[i];
    spectatorFrac[i]   = eSpectator[i]   / eTotal[i];

    printf("%-10s%10.2f%10.2f%16.3f%16.3f%16.3f%14.3f\n", CENT_LABELS[i], meanNpart[i], meanNspec[i],
           eSpectator[i]/1000.0, eParticipant[i]/1000.0, eProduced[i]/1000.0, eTotal[i]/1000.0);
  }

  cout << endl << "Consistency check: E_tot should equal 16*T1_lab = " << 16.0*T1_lab/1000.0
       << " TeV in EVERY bin (Npart+Nspec=32 always) -- if any row above deviates, something's off." << endl;

  cout << endl;
  printf("%-10s%16s%16s%16s\n", "Cent", "fireball %", "participant %", "spectator %");
  for(int i = 0; i < N_CENT_BINS; i++){
    printf("%-10s%16.1f%16.1f%16.1f\n", CENT_LABELS[i], 100.0*fireballFrac[i], 100.0*participantFrac[i], 100.0*spectatorFrac[i]);
  }

  //--------------------------------------------------------------------------
  // Cross-section-weighted min-bias average -- same weighting scheme as the
  // collider-frame macro (a statement about how often each centrality bin
  // occurs, independent of which energy frame the fractions were computed in).
  //--------------------------------------------------------------------------
  double binWidth[N_CENT_BINS], sumBinWidth = 0.0;
  for(int i = 0; i < N_CENT_BINS; i++){
    binWidth[i] = CENT_EDGES_FRAC[i + 1] - CENT_EDGES_FRAC[i];
    sumBinWidth += binWidth[i];
  }
  if(fabs(sumBinWidth - 1.0) > 1e-9){
    cout << "WARNING: CENT_EDGES_FRAC bin widths sum to " << sumBinWidth
         << ", not 1.0 -- check the edges array before trusting the min-bias average." << endl;
  }

  double minBiasFireballFrac = 0.0, minBiasParticipantFrac = 0.0, minBiasSpectatorFrac = 0.0;
  for(int i = 0; i < N_CENT_BINS; i++){
    minBiasFireballFrac    += binWidth[i] * fireballFrac[i];
    minBiasParticipantFrac += binWidth[i] * participantFrac[i];
    minBiasSpectatorFrac   += binWidth[i] * spectatorFrac[i];
  }

  cout << endl << "Cross-section-weighted min-bias (0-100%) average, weighted by each bin's"
       << " percentile width (5,5,10,20,40,20):" << endl;
  printf("%-10s%16.1f%16.1f%16.1f\n", "0-100%", 100.0*minBiasFireballFrac,
         100.0*minBiasParticipantFrac, 100.0*minBiasSpectatorFrac);
  cout << "(For comparison, the collider-frame macro's min-bias fireball fraction uses the SAME"
       << " weighting scheme but ~73% CM-frame stopping instead of ~87% lab-frame stopping --"
       << " expect this fireball_frac to read noticeably higher.)" << endl;

  //--------------------------------------------------------------------------
  // Cross-section-weighted CONTRIBUTION per bin (added 2026-09-01 per
  // Andrew's request) -- NOT the same thing as fireballFrac[i] above.
  // fireballFrac[i] answers "of THIS bin's own energy, what fraction is
  // fireball" (always sums to 100% within a bin, regardless of how much of
  // the total cross section that bin represents). This block instead scales
  // each bin's fraction by its own cross-section share (binWidth[i], same
  // CENT_EDGES_FRAC used for the min-bias average above), answering "how
  // many percentage-points of the OVERALL 0-100% total does this bin
  // contribute." Summed over all six bins, each column here reproduces the
  // min-bias numbers directly by construction:
  //   sum_i fireballFrac[i]*binWidth[i] = minBiasFireballFrac
  //--------------------------------------------------------------------------
  double fireballContrib[N_CENT_BINS], participantContrib[N_CENT_BINS], spectatorContrib[N_CENT_BINS];
  for(int i = 0; i < N_CENT_BINS; i++){
    fireballContrib[i]    = fireballFrac[i]    * binWidth[i];
    participantContrib[i] = participantFrac[i] * binWidth[i];
    spectatorContrib[i]   = spectatorFrac[i]   * binWidth[i];
  }
  cout << endl << "Cross-section-weighted CONTRIBUTION per bin (this bin's fraction x its own"
       << " percentile width -- these six rows sum to the min-bias row above):" << endl;
  printf("%-10s%16s%16s%16s\n", "Cent", "fireball pp", "participant pp", "spectator pp");
  double sumFbContrib = 0.0, sumPtContrib = 0.0, sumSpContrib = 0.0;
  for(int i = 0; i < N_CENT_BINS; i++){
    printf("%-10s%16.2f%16.2f%16.2f\n", CENT_LABELS[i], 100.0*fireballContrib[i],
           100.0*participantContrib[i], 100.0*spectatorContrib[i]);
    sumFbContrib += fireballContrib[i]; sumPtContrib += participantContrib[i]; sumSpContrib += spectatorContrib[i];
  }
  printf("%-10s%16.2f%16.2f%16.2f  (sum -- should match the min-bias row above)\n", "SUM",
         100.0*sumFbContrib, 100.0*sumPtContrib, 100.0*sumSpContrib);

  //--------------------------------------------------------------------------
  // Plot: stacked bar, same equal-width geometry as before -- what changed
  // is that each of the six centrality bars' HEIGHT is now that bin's
  // cross-section-weighted CONTRIBUTION (fireballContrib/participantContrib/
  // spectatorContrib above), not its own 100%-normalized internal split. So
  // 0-5% (5% of the cross section) now tops out around 5 on this axis while
  // 40-80% (40% of the cross section) tops out around 40 -- the bars' total
  // heights are what carry the "how much of the overall total does this bin
  // represent" information, while the color proportions WITHIN each bar
  // still show that bin's own internal fireball/participant/spectator split.
  // The Min-Bias bar is unchanged (it was already the full weighted total).
  //--------------------------------------------------------------------------
  const int N_PLOT_BINS = N_CENT_BINS + 1;
  const int MINBIAS_BIN = N_CENT_BINS + 1;
  TCanvas* c = new TCanvas("cEnergyDoseFractionLab", "Energy Dose Fraction (Fixed-Target Lab Frame)", 950, 700);
  gPad->SetBottomMargin(0.24);
  gPad->SetGridy();

  TH1D* hStackTop = new TH1D("hStackTopLab", "", N_PLOT_BINS, 0, N_PLOT_BINS);
  TH1D* hStackMid = new TH1D("hStackMidLab", "", N_PLOT_BINS, 0, N_PLOT_BINS);
  TH1D* hStackBot = new TH1D("hStackBotLab", "", N_PLOT_BINS, 0, N_PLOT_BINS);

  for(int i = 0; i < N_CENT_BINS; i++){
    double fb = 100.0 * fireballContrib[i];
    double pt = 100.0 * participantContrib[i];
    double sp = 100.0 * spectatorContrib[i];
    hStackTop->SetBinContent(i + 1, fb + pt + sp);
    hStackMid->SetBinContent(i + 1, fb + pt);
    hStackBot->SetBinContent(i + 1, fb);
    hStackTop->GetXaxis()->SetBinLabel(i + 1, CENT_LABELS[i]);
  }
  {
    double fb = 100.0 * minBiasFireballFrac;
    double pt = 100.0 * minBiasParticipantFrac;
    double sp = 100.0 * minBiasSpectatorFrac;
    hStackTop->SetBinContent(MINBIAS_BIN, fb + pt + sp);
    hStackMid->SetBinContent(MINBIAS_BIN, fb + pt);
    hStackBot->SetBinContent(MINBIAS_BIN, fb);
    hStackTop->GetXaxis()->SetBinLabel(MINBIAS_BIN, "0-100%");
  }

  double barWidth = 0.6, barOffset = (1.0 - barWidth) / 2.0;
  hStackTop->SetBarWidth(barWidth); hStackTop->SetBarOffset(barOffset);
  hStackMid->SetBarWidth(barWidth); hStackMid->SetBarOffset(barOffset);
  hStackBot->SetBarWidth(barWidth); hStackBot->SetBarOffset(barOffset);

  hStackTop->SetFillColor(COLOR_SPECTATOR);
  hStackMid->SetFillColor(COLOR_PARTICIPANT);
  hStackBot->SetFillColor(COLOR_FIREBALL);

  hStackTop->SetMaximum(100.0);
  hStackTop->SetMinimum(0.0);
  hStackTop->GetXaxis()->SetLabelSize(0.040);
  hStackTop->GetYaxis()->SetTitle("Contribution to Total Kinetic Energy [percentage points]");
  hStackTop->GetXaxis()->SetTitle("Centrality");
  hStackTop->GetXaxis()->CenterTitle(kTRUE);
  hStackTop->SetTitle("O+O #sqrt{s_{NN}}=200 GeV -- Cross-Section-Weighted Energy Dose Contribution (Fixed-Target Lab Frame)");

  hStackTop->Draw("bar");
  hStackMid->Draw("bar same");
  hStackBot->Draw("bar same");

  TLine* sepLine = new TLine((double) N_CENT_BINS, 0.0, (double) N_CENT_BINS, 100.0);
  sepLine->SetLineColor(TColor::GetColor("#B8B8B8"));
  sepLine->SetLineStyle(2);
  sepLine->SetLineWidth(1);
  sepLine->Draw();

  // Labels use the same CONTRIBUTION values as the bar heights above (not
  // fireballFrac[i] raw) so the printed numbers match what's actually drawn.
  for(int i = 0; i < N_CENT_BINS; i++){
    double x = hStackTop->GetXaxis()->GetBinCenter(i + 1);
    double fb = 100.0 * fireballContrib[i];
    double pt = 100.0 * participantContrib[i];
    double sp = 100.0 * spectatorContrib[i];
    DrawPercentLabel(x, fb / 2.0, fb, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt / 2.0, pt, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt + sp / 2.0, sp, COLOR_TEXT_ON_LIGHT);
  }
  {
    double x = hStackTop->GetXaxis()->GetBinCenter(MINBIAS_BIN);
    double fb = 100.0 * minBiasFireballFrac;
    double pt = 100.0 * minBiasParticipantFrac;
    double sp = 100.0 * minBiasSpectatorFrac;
    DrawPercentLabel(x, fb / 2.0, fb, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt / 2.0, pt, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt + sp / 2.0, sp, COLOR_TEXT_ON_LIGHT);
  }

  TLegend* leg = new TLegend(0.12, 0.08, 0.90, 0.16);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hStackBot, "Fireball-produced (net energy converted to new particles)", "f");
  leg->AddEntry(hStackMid, "Participant-retained (both sides, after stopping)", "f");
  leg->AddEntry(hStackTop, "Spectator (projectile side, full T1_lab, untouched)", "f");
  leg->Draw();

  TLatex caveat;
  caveat.SetNDC();
  caveat.SetTextFont(42);
  caveat.SetTextSize(0.020);
  caveat.SetTextColor(COLOR_CAVEAT);
  caveat.SetTextAlign(21);
  caveat.DrawLatex(0.5, 0.02,
    Form("Bar height = this bin's cross-section-weighted CONTRIBUTION to the 0-100%% total, not its own 100%%-normalized split. T1_lab=%.1f TeV, delta_y=%.2f (BRAHMS Au+Au) -- see macro header.", T1_lab/1000.0, DELTA_Y));

  string outPath = a_outputDir + "/OO200_EnergyDoseFractionFixedTargetLab.png";
  c->SaveAs(outPath.c_str());
  cout << endl << "Wrote " << outPath << endl;

  cout << endl << "Reminder: this is the SAME sqrt(s_NN)=200 GeV STAR O+O collision as the collider-frame"
       << " macro, reinterpreted in the frame where the target is at rest -- not a different energy"
       << " regime, and not yet weighted by how often a real GCR event at this energy actually occurs"
       << " (that needs the GCR differential flux at ~" << E1_lab/1000.0 << " TeV/nucleon folded in separately)." << endl;
}
