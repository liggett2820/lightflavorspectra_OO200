// MarsGCRPenetrationDepth.cpp
//
// PURPOSE: Continuous-slowing-down-approximation (CSDA) penetration depth of
// galactic cosmic ray (GCR) primary ions into Mars regolith, as a function of
// incident kinetic energy per nucleon. General physics calculation, not part
// of the lightflavorspectra_OO200 STAR analysis pipeline.
//
// NOTE ON VERIFICATION STATUS: the physics/integration code below (Bragg's
// rule, Sternheimer density effect, Bethe-Bloch, CSDA integration) is
// UNCHANGED from an earlier plain-C++ version of this file that was compiled
// with g++ and run directly in a sandbox with no ROOT installed. That run
// caught a real integration bug (see CSDARange_gcm2 below) via a numeric
// sanity check against NIST PSTAR values, and a follow-up exact z^2
// charge-scaling self-consistency check across all 8 species. This version
// swaps the output stage from a CSV file to a ROOT TFile/TTree, which means
// it now requires ROOT to compile -- so unlike that earlier verified run,
// THIS file has not been compiled or executed by me; only the output block
// changed, and that change was not exercised. A backup of the last
// g++-compiled, numerically-verified CSV-output version is kept alongside
// this file as MarsGCRPenetrationDepth_CSV_verified_reference.cpp.bak.
//
// SCOPE (as explicitly chosen by Andrew):
//   - Species: representative GCR reference ions H through Fe (H, He, C, O,
//     Ne, Mg, Si, Fe) -- the standard anchor-ion set used in NASA space
//     radiation literature (e.g. Badhwar-O'Neill / HZETRN reference nuclei),
//     not an exhaustive Z=1..26 scan.
//   - Target: Mars regolith/soil ONLY (i.e. depth already at the surface,
//     post-atmosphere) -- NOT atmospheric attenuation.
//   - Method: CSDA range from electromagnetic (Bethe-Bloch) stopping power.
//     EXPLICIT LIMITATION, chosen knowingly: this does NOT include nuclear
//     fragmentation/interaction losses. For light ions (H, He) at GCR energies
//     electromagnetic stopping dominates and this is a reasonable estimate.
//     For heavy ions (Ne and up, especially Fe) nuclear interactions remove a
//     significant fraction of primaries well before they range out
//     electromagnetically -- real transport codes (e.g. HZETRN) show heavy-ion
//     "depth-dose" is governed as much by the nuclear interaction mean free
//     path (attenuation of the primary beam) as by electromagnetic range. This
//     calculation gives the electromagnetic CSDA range only, which should be
//     read as an UPPER BOUND on how far a heavy-ion primary travels before
//     either stopping or fragmenting -- NOT a full transport-corrected result.
//
// MARS REGOLITH COMPOSITION -- SOURCED, NOT ASSUMED:
//   NASA NTRS technical report "Chemical, Mineralogical, and Physical
//   Properties of Martian Dust and Soil" (Gusev Crater "Panda" subclass
//   average, wt% oxides):
//     SiO2 46.52, Al2O3 10.46, FeO 12.18, Fe2O3 4.20, MgO 8.93, CaO 6.27,
//     Na2O 3.02, K2O 0.41, TiO2 0.87, SO3 4.90, Cl 0.61
//   (https://ntrs.nasa.gov/api/citations/20170005414/downloads/20170005414.pdf)
//   Sums to 98.37% (remainder is trace oxides/elements not itemized in the
//   source, e.g. MnO, P2O5) -- renormalized to 100% below, documented at the
//   point of renormalization rather than silently absorbed.
//
// MARS REGOLITH BULK DENSITY -- SOURCED, NOT ASSUMED, WITH A REAL RANGE:
//   From the JSC Mars-1 simulant paper (Allen et al., LPSC 1998,
//   https://www.lpi.usra.edu/meetings/LPSC98/pdf/1690.pdf), citing actual
//   Martian in-situ measurements:
//     Mars Pathfinder average bulk soil density: 1.52 g/cm^3  <- used as NOMINAL
//     Viking 1 drift material:                   1.2 +/- 0.2 g/cm^3
//     Viking XRFS fines:                          1.10 +/- 0.15 g/cm^3
//   This is loose regolith/soil, NOT competent bedrock (which would be denser,
//   ~2.7-3.3 g/cm^3 for basalt, and would give a correspondingly SHORTER
//   physical depth for the same areal density in g/cm^2 -- the areal-density
//   range in g/cm^2 is density-independent; only the g/cm^2 -> cm conversion
//   depends on which density is assumed). This macro reports depth at the
//   Pathfinder nominal density AND at the Viking-fines low end (1.10) and
//   Viking-drift high end (1.60, using the upper edge of Viking drift's
//   1.2+/-0.2) so the density-driven spread is visible rather than hidden
//   behind one silently-chosen number.
//
// PHYSICS:
//   Mean mass stopping power via the Bethe-Bloch formula (MeV cm^2/g):
//     -dE/dx = K z^2 (Z/A) (1/beta^2) [ 0.5 ln(2 me_c2 beta^2 gamma^2 Tmax / I^2) - beta^2 - delta/2 ]
//   K = 0.307075 MeV mol^-1 cm^2 (standard constant = 4 pi N_A r_e^2 me c^2).
//   Target (Z/A) and mean excitation energy I are combined from the elemental
//   composition via Bragg's additivity rule (standard practice for compounds/
//   mixtures without directly measured I).
//   Elemental I values (eV) are the standard ICRU-37 tabulated values.
//   Density-effect correction delta(beta*gamma) uses the Sternheimer generic-
//   material parameterization (Sternheimer, Berger & Seltzer, At. Data Nucl.
//   Data Tables 30, 261 (1984)) -- the standard fallback recipe used when no
//   material-specific fitted Sternheimer parameters exist (true here: nobody
//   has published fitted Sternheimer coefficients for "Mars regolith" as a
//   named material). Without SOME density-effect treatment, dE/dx would be
//   spuriously overestimated at high energy (unphysical unbounded relativistic
//   rise), which would UNDERESTIMATE range/depth at the GeV/nucleon end of the
//   GCR spectrum -- this correction exists specifically to avoid that bias.
//   No shell correction, no Barkas/Bloch correction -- both are small at the
//   MeV/nucleon-and-up energies tabulated here and are omitted for simplicity;
//   flagged as a (minor) known limitation, not silently ignored.
//
// CSDA range: R(T0) = Integral from E_CUTOFF to T0 of dE' / (dE/dx)(E'),
// in g/cm^2, then divided by density to get cm. E_CUTOFF = 1 MeV/nucleon --
// below this the Bethe-Bloch formula (no shell correction) is not reliable,
// but the contribution to total range from below 1 MeV/nucleon is negligible
// for the GCR energies tabulated (>= 10 MeV/nucleon) since dE/dx is very large
// (short range) in that low-energy region anyway.
//
// GCR REFERENCE IONS: H, He, C, O, Ne, Mg, Si, Fe -- most abundant isotope
// mass numbers used (A = 1, 4, 12, 16, 20, 24, 28, 56).
//
// OUTPUT: MarsGCRPenetrationDepth.root -- one TTree named "MarsGCRPenetrationDepth"
// with one entry per (species, energy grid point), branches:
//   species/C (char[8]), Z/I, A/D, KE_MeV_per_nucleon/D, KE_total_MeV/D,
//   range_g_cm2/D, depth_cm_nominal_rho1p52/D, depth_cm_lowdensity_rho1p10/D,
//   depth_cm_highdensity_rho1p60/D
// Read back with, e.g.:
//   TFile* f = TFile::Open("MarsGCRPenetrationDepth.root");
//   TTree* t = (TTree*) f->Get("MarsGCRPenetrationDepth");
//
// STATUS: the physics/integration code is carried over unchanged from a
// version that WAS compiled with g++ and run in-sandbox -- see the
// verification note above. This ROOT-output version has NOT been compiled
// or run by me (no ROOT in this sandbox); review before trusting its output.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <cstring>

#include "TFile.h"
#include "TTree.h"

using namespace std;

//============================================================================
// Physical constants
//============================================================================
const double K_BETHE      = 0.307075;   // MeV mol^-1 cm^2
const double ME_C2_MEV    = 0.510999;   // electron rest energy, MeV
const double U_MEV        = 931.494;    // atomic mass unit, MeV/c^2
const double E_CUTOFF_MEV_PER_U = 1.0;  // low-energy integration cutoff, see header

//============================================================================
// Elemental data: Z, standard atomic weight (g/mol), ICRU-37 mean excitation
// energy I (eV). Only elements appearing in the Mars regolith composition
// below are included.
//============================================================================
struct ElementData { int Z; double atomicWeight; double I_eV; };
map<string, ElementData> ElementTable(){
  map<string, ElementData> t;
  t["O"]  = {8,  15.999, 95.0};
  t["Si"] = {14, 28.085, 173.0};
  t["Al"] = {13, 26.982, 166.0};
  t["Fe"] = {26, 55.845, 286.0};
  t["Mg"] = {12, 24.305, 156.0};
  t["Ca"] = {20, 40.078, 191.0};
  t["Na"] = {11, 22.990, 149.0};
  t["K"]  = {19, 39.098, 190.0};
  t["Ti"] = {22, 47.867, 233.0};
  t["S"]  = {16, 32.060, 180.0};
  t["Cl"] = {17, 35.450, 174.0};
  return t;
}

//============================================================================
// Oxide composition of Mars regolith (Gusev Crater "Panda" average, see
// header for source). Each entry: element symbol, atom count of that element
// per formula unit, oxygen count per formula unit, wt% of this oxide in the
// bulk composition (as reported by the source -- NOT yet renormalized).
//============================================================================
struct OxideEntry { string element; int elementCount; int oxygenCount; double wtPercentOxide; };
vector<OxideEntry> MarsRegolithOxides(){
  return {
    {"Si", 1, 2, 46.52},  // SiO2
    {"Al", 2, 3, 10.46},  // Al2O3
    {"Fe", 1, 1, 12.18},  // FeO
    {"Fe", 2, 3, 4.20},   // Fe2O3
    {"Mg", 1, 1, 8.93},   // MgO
    {"Ca", 1, 1, 6.27},   // CaO
    {"Na", 2, 1, 3.02},   // Na2O
    {"K",  2, 1, 0.41},   // K2O
    {"Ti", 1, 2, 0.87},   // TiO2
    {"S",  1, 3, 4.90},   // SO3 (S itself is not a target constituent used in
                          // Bragg's rule below except via this oxide's O and,
                          // notionally, S -- S is added to the element table
                          // via its own entry, matching how FeO/Fe2O3 both
                          // route into "Fe")
  };
}
// Cl reported as elemental wt%, not an oxide, in APXS-style data (see source).
const double CL_ELEMENTAL_WTPCT = 0.61;

//============================================================================
// DecomposeToElementalMassFractions: oxide wt% -> renormalized elemental mass
// fractions (fractions sum to 1.0 over the full composition, computed from
// atomic weights and stoichiometry rather than hand-entered, to avoid
// transcription error).
//============================================================================
map<string, double> DecomposeToElementalMassFractions(){
  map<string, ElementData> elems = ElementTable();
  map<string, double> elementWtPct; // running wt%, not yet renormalized

  for(const OxideEntry& ox : MarsRegolithOxides()){
    double elementMass = ox.elementCount * elems[ox.element].atomicWeight;
    double oxygenMass  = ox.oxygenCount  * elems["O"].atomicWeight;
    double molarMass   = elementMass + oxygenMass;
    double elementFracInOxide = elementMass / molarMass;
    double oxygenFracInOxide  = oxygenMass  / molarMass;
    elementWtPct[ox.element] += ox.wtPercentOxide * elementFracInOxide;
    elementWtPct["O"]        += ox.wtPercentOxide * oxygenFracInOxide;
  }
  elementWtPct["Cl"] += CL_ELEMENTAL_WTPCT;

  double total = 0.0;
  for(auto& kv : elementWtPct) total += kv.second;
  cout << "Mars regolith elemental decomposition: raw sum = " << total
       << " wt% (source oxides summed to 98.37%, renormalizing to 100%)." << endl;

  map<string, double> massFrac;
  for(auto& kv : elementWtPct) massFrac[kv.first] = kv.second / total; // renormalize to sum=1
  return massFrac;
}

//============================================================================
// Bragg's additivity rule: composite (Z/A) and mean excitation energy I for
// the regolith mixture, from elemental mass fractions.
//============================================================================
struct TargetMixture { double ZA; double I_eV; };
TargetMixture ComputeBraggMixture(){
  map<string, ElementData> elems = ElementTable();
  map<string, double> w = DecomposeToElementalMassFractions();

  double ZAmix = 0.0;
  double weightedLnI = 0.0;
  for(auto& kv : w){
    const string& el = kv.first;
    double wi = kv.second;
    double ZAi = elems[el].Z / elems[el].atomicWeight;
    ZAmix += wi * ZAi;
    weightedLnI += wi * ZAi * log(elems[el].I_eV);
  }
  double Imix = exp(weightedLnI / ZAmix);

  cout << fixed << setprecision(5);
  cout << "Composite target (Z/A) = " << ZAmix << " mol/g" << endl;
  cout << "Composite target I     = " << Imix << " eV" << endl;
  cout.unsetf(ios::fixed);

  return {ZAmix, Imix};
}

//============================================================================
// Sternheimer generic-material density-effect parameters (see header for
// citation). rho in g/cm^3, I in eV.
//============================================================================
struct SternheimerParams { double Cbar, X0, X1, a, m; };
SternheimerParams ComputeSternheimerGeneric(double ZAmix, double I_eV, double rho_gcm3){
  double hbarOmegaP_eV = 28.816 * sqrt(rho_gcm3 * ZAmix); // plasma energy, eV
  double Cbar = 2.0 * log(I_eV / hbarOmegaP_eV) + 1.0;
  double X0, X1;
  const double m = 3.0;
  if(I_eV < 100.0){
    X1 = 2.0;
    X0 = (Cbar < 3.681) ? 0.2 : (0.326 * Cbar - 1.0);
  } else {
    X1 = 3.0;
    X0 = (Cbar < 5.215) ? 0.2 : (0.326 * Cbar - 1.5);
  }
  double a = (Cbar - 4.606 * X0) / pow(X1 - X0, m);
  return {Cbar, X0, X1, a, m};
}

double DensityEffectDelta(double betaGamma, const SternheimerParams& s){
  double X = log10(betaGamma);
  if(X < s.X0) return 0.0; // delta0 = 0, non-conductor assumption
  if(X < s.X1) return 4.606 * X - s.Cbar + s.a * pow(s.X1 - X, s.m);
  return 4.606 * X - s.Cbar;
}

//============================================================================
// Bethe-Bloch mass stopping power, MeV cm^2/g, for an ion of charge z and
// mass A_ion (amu) at kinetic energy T (MeV per nucleon) in the target
// mixture (ZAmix, I_eV, sternheimer params for that mixture's density).
//============================================================================
double MassStoppingPower(double T_MeV_per_u, int z, double A_ion_u,
                          double ZAmix, double I_eV, const SternheimerParams& s){
  double M_ion_MeV = A_ion_u * U_MEV;           // total rest mass of the ion
  double totalKE_MeV = T_MeV_per_u * A_ion_u;    // total kinetic energy
  double gamma = 1.0 + totalKE_MeV / M_ion_MeV;
  double beta2 = 1.0 - 1.0 / (gamma * gamma);
  if(beta2 <= 0.0) return 1e12; // guard, shouldn't trigger above E_CUTOFF
  double beta = sqrt(beta2);
  double betaGamma = beta * gamma;

  double meOverM = ME_C2_MEV / M_ion_MeV;
  double Tmax = (2.0 * ME_C2_MEV * beta2 * gamma * gamma)
              / (1.0 + 2.0 * gamma * meOverM + meOverM * meOverM);

  double I_MeV = I_eV * 1.0e-6;
  double delta = DensityEffectDelta(betaGamma, s);

  double logTerm = 0.5 * log(2.0 * ME_C2_MEV * beta2 * gamma * gamma * Tmax / (I_MeV * I_MeV));
  double bracket = logTerm - beta2 - delta / 2.0;

  double dEdx = K_BETHE * z * z * ZAmix * (1.0 / beta2) * bracket;
  return dEdx; // MeV cm^2/g
}

//============================================================================
// CSDA range via straightforward fine-step numerical integration in
// log-energy (dE/dx varies smoothly, no adaptive stepping needed).
//============================================================================
double CSDARange_gcm2(double T0_MeV_per_u, int z, double A_ion_u,
                       double ZAmix, double I_eV, const SternheimerParams& s){
  if(T0_MeV_per_u <= E_CUTOFF_MEV_PER_U) return 0.0;
  const int N_STEPS = 20000;
  double logLo = log(E_CUTOFF_MEV_PER_U);
  double logHi = log(T0_MeV_per_u);
  double dLog = (logHi - logLo) / N_STEPS;

  double range = 0.0;
  for(int i = 0; i < N_STEPS; i++){
    double logE_a = logLo + i * dLog;
    double logE_b = logLo + (i + 1) * dLog;
    double Ea = exp(logE_a), Eb = exp(logE_b);
    double dEdx_a = MassStoppingPower(Ea, z, A_ion_u, ZAmix, I_eV, s);
    double dEdx_b = MassStoppingPower(Eb, z, A_ion_u, ZAmix, I_eV, s);
    // trapezoid in E (not log E): dR = dE / (dE/dx). Standard trapezoidal rule:
    // integral over [Ea,Eb] of f(E) dE ~= (Eb-Ea)/2 * (f(Ea)+f(Eb)), with
    // f(E) = 1/(dE/dx)(E). (An earlier version of this line had an extra
    // duplicated *0.5 here, silently halving every range/depth value in the
    // output -- caught by sanity-checking the 100 MeV proton range against
    // the well-known NIST PSTAR value in water (~7.7 g/cm^2): this material's
    // lower Z/A and higher I should give a LONGER range than water, not a
    // shorter one, and the halved output was doing the opposite. Fixed.)
    range += 0.5 * ((Eb - Ea) / dEdx_a + (Eb - Ea) / dEdx_b);
  }
  return range; // g/cm^2
}

//============================================================================
// GCR reference ions
//============================================================================
struct GCRIon { string name; int Z; double A; };
vector<GCRIon> GCRReferenceIons(){
  return {
    {"H",  1, 1.0},
    {"He", 2, 4.0},
    {"C",  6, 12.0},
    {"O",  8, 16.0},
    {"Ne", 10, 20.0},
    {"Mg", 12, 24.0},
    {"Si", 14, 28.0},
    {"Fe", 26, 56.0},
  };
}

//============================================================================
// Main
//============================================================================
int main(){
  const double RHO_NOMINAL = 1.52; // g/cm^3, Mars Pathfinder average bulk soil density
  const double RHO_LOW     = 1.10; // g/cm^3, Viking XRFS fines
  const double RHO_HIGH    = 1.60; // g/cm^3, upper edge of Viking 1 drift material (1.2+0.2*2 rounded)

  cout << "=== MarsGCRPenetrationDepth ===" << endl;
  TargetMixture mix = ComputeBraggMixture();
  SternheimerParams sSternheimer = ComputeSternheimerGeneric(mix.ZA, mix.I_eV, RHO_NOMINAL);
  cout << "Sternheimer generic params: Cbar=" << sSternheimer.Cbar
       << " X0=" << sSternheimer.X0 << " X1=" << sSternheimer.X1
       << " a=" << sSternheimer.a << " m=" << sSternheimer.m << endl;

  // Energy grid: 10 MeV/nucleon to 200 GeV/nucleon, log-spaced -- spans the
  // bulk of the GCR energy spectrum relevant to surface dose, extended up to
  // 200 GeV/nucleon per request (was 10 GeV/nucleon). Point count raised from
  // 40 to 50 to keep roughly the same log-spacing density per decade as
  // before (was ~13.3 points/decade over 3 decades; now ~11.6/decade over
  // 4.3 decades) rather than thinning out the curve over the wider range.
  vector<double> energyGrid_MeV_per_u;
  const int N_ENERGY_POINTS = 50;
  double eLo = 10.0, eHi = 200000.0;
  for(int i = 0; i < N_ENERGY_POINTS; i++){
    double frac = (double) i / (N_ENERGY_POINTS - 1);
    energyGrid_MeV_per_u.push_back(eLo * pow(eHi / eLo, frac));
  }

  TFile* outFile = new TFile("MarsGCRPenetrationDepth.root", "RECREATE");
  TTree* tree = new TTree("MarsGCRPenetrationDepth",
                           "Mars regolith CSDA GCR penetration depth vs. energy per nucleon");

  Char_t species_out[8];
  Int_t    Z_out;
  Double_t A_out;
  Double_t KE_MeV_per_nucleon_out;
  Double_t KE_total_MeV_out;
  Double_t range_g_cm2_out;
  Double_t depth_cm_nominal_rho1p52_out;
  Double_t depth_cm_lowdensity_rho1p10_out;
  Double_t depth_cm_highdensity_rho1p60_out;

  tree->Branch("species",                       species_out,                     "species/C");
  tree->Branch("Z",                              &Z_out,                          "Z/I");
  tree->Branch("A",                              &A_out,                          "A/D");
  tree->Branch("KE_MeV_per_nucleon",             &KE_MeV_per_nucleon_out,         "KE_MeV_per_nucleon/D");
  tree->Branch("KE_total_MeV",                   &KE_total_MeV_out,               "KE_total_MeV/D");
  tree->Branch("range_g_cm2",                    &range_g_cm2_out,                "range_g_cm2/D");
  tree->Branch("depth_cm_nominal_rho1p52",       &depth_cm_nominal_rho1p52_out,   "depth_cm_nominal_rho1p52/D");
  tree->Branch("depth_cm_lowdensity_rho1p10",    &depth_cm_lowdensity_rho1p10_out, "depth_cm_lowdensity_rho1p10/D");
  tree->Branch("depth_cm_highdensity_rho1p60",   &depth_cm_highdensity_rho1p60_out, "depth_cm_highdensity_rho1p60/D");

  for(const GCRIon& ion : GCRReferenceIons()){
    cout << "--- " << ion.name << " (Z=" << ion.Z << ", A=" << ion.A << ") ---" << endl;
    strncpy(species_out, ion.name.c_str(), sizeof(species_out) - 1);
    species_out[sizeof(species_out) - 1] = '\0';
    Z_out = ion.Z;
    A_out = ion.A;

    for(double T : energyGrid_MeV_per_u){
      double range_gcm2 = CSDARange_gcm2(T, ion.Z, ion.A, mix.ZA, mix.I_eV, sSternheimer);
      KE_MeV_per_nucleon_out = T;
      KE_total_MeV_out = T * ion.A;
      range_g_cm2_out = range_gcm2;
      depth_cm_nominal_rho1p52_out = range_gcm2 / RHO_NOMINAL;
      depth_cm_lowdensity_rho1p10_out = range_gcm2 / RHO_LOW;
      depth_cm_highdensity_rho1p60_out = range_gcm2 / RHO_HIGH;
      tree->Fill();
    }
    // print a few representative rows to stdout for a quick sanity read
    for(double Tcheck : {100.0, 1000.0}){
      double range_gcm2 = CSDARange_gcm2(Tcheck, ion.Z, ion.A, mix.ZA, mix.I_eV, sSternheimer);
      cout << "  " << Tcheck << " MeV/u: range=" << range_gcm2 << " g/cm^2, depth(nominal rho)="
           << range_gcm2 / RHO_NOMINAL << " cm" << endl;
    }
  }

  outFile->cd();
  tree->Write();
  outFile->Close();
  cout << "Wrote MarsGCRPenetrationDepth.root (TTree \"MarsGCRPenetrationDepth\", "
       << GCRReferenceIons().size() * energyGrid_MeV_per_u.size() << " entries)" << endl;
  return 0;
}
