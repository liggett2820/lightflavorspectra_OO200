// ============================================================================
// GCR_InteractionLength_MarsRegolith.C
//
// Calculates and plots nuclear interaction (reaction) mean free paths for
// representative Galactic Cosmic Ray (GCR) primary ions in Mars regolith,
// and uses them to show how quickly primaries undergo at least one nuclear
// interaction with depth -- i.e., how much of the GCR-induced radiation
// field at depth should be expected to come from secondary particles /
// nuclear fragmentation products rather than surviving primaries.
//
// PHYSICS MODEL
// --------------
// Nuclear reaction cross section: Bradt-Peters geometric overlap formula,
//
//     sigma_R(Ap,At) = pi * r0^2 * ( Ap^(1/3) + At^(1/3) - b0 )^2
//
// the standard "black disk / geometric overlap" approximation for heavy-ion
// total reaction cross sections. It is energy-independent and is a good
// approximation above roughly ~200 MeV/nucleon, which covers the bulk of
// the GCR flux that dominates dose (the GCR spectrum peaks in dose-relevant
// flux around a few hundred MeV/n to a few GeV/n).
// Reference: Bradt & Peters, Phys. Rev. 77, 54 (1950); geometric
// cross-section systematics are also summarized in standard nuclear/particle
// physics texts (e.g. Povh, Rith, Scholz, Zetsche, "Particle and Nuclear
// Physics").
//
// r0 and b0 below are representative literature values for this class of
// parameterization. Published fits differ (Silberberg-Tsao; Tripathi et al.
// 1996/1999, used in NASA's HZETRN; Sihver et al. add energy-dependent
// corrections) -- if your write-up needs to match a specific published
// cross-section model, replace sigmaReaction_mb() below; the rest of the
// macro (mean free path, survival fraction) is unaffected by that choice.
//
// TARGET COMPOSITION (Mars regolith)
// -----------------------------------
// Representative bulk Mars regolith composition, built from typical
// oxide-weight-percent soil averages broadly consistent with multi-lander
// APXS bulk-soil measurements (Viking/Pathfinder/MER-era Mars soil
// geochemistry summaries, e.g. Taylor & McLennan (2009)-type reviews):
//
//   SiO2 45%   FeO  17%   Al2O3 9%   CaO 7%
//   MgO   8%   SO3   6%   Na2O  2%   K2O 0.5%
//   TiO2  1%   Cl  0.5%   other 4%  (other lumped at <A> ~ 23)
//
// Elemental mass fractions are derived from these oxide fractions via
// standard molar-mass ratios (shown next to each entry below). ADJUST
// regolithComposition[] if you have a specific measured composition (e.g.
// for a specific landing site) you need to match.
//
// OUTPUT
// ------
//  1) Console table: reaction cross section / interaction mean free path
//     (g/cm^2 and cm) for each GCR reference ion in Mars regolith, plus
//     the fraction of each ion species that has undergone >=1 nuclear
//     interaction by several representative depths.
//  2) Canvas 1: interaction length (cm) per ion species, with a shaded
//     band showing a typical near-surface regolith shielding-depth range.
//  3) Canvas 2: survival fraction (fraction of primaries that have NOT
//     yet undergone a nuclear interaction) vs. depth, for each ion --
//     (1 - survival) is the fraction that HAS produced secondary /
//     fragmentation products by that depth.
//
// USAGE
//   root -l GCR_InteractionLength_MarsRegolith.C
//   (written for a standard ROOT6/Cling install; if you're on an old
//   ROOT5/CINT setup the range-based for-loops below will need rewriting
//   as index loops)
// ============================================================================

#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
const double kAvogadro = 6.02214076e23;      // mol^-1

// Bradt-Peters geometric cross-section parameters (representative values --
// see header comment; adjust to match a specific published parameterization
// if needed).
const double kR0_fm = 1.20;                  // fm
const double kB0_fm = 1.00;                  // fm  (nuclear "overlap"/transparency term)
const double kFm2_to_mb = 10.0;              // 1 fm^2 = 10 mb

// Assumed bulk Mars regolith density range (loose, near-surface regolith).
// Broadly consistent with Mars lander/rover bulk-density estimates (e.g.
// Golombek et al.); adjust for compacted vs. loose regolith.
const double kRegolithDensity_gcm3      = 1.60;  // "nominal" loose regolith
const double kRegolithDensity_gcm3_low  = 1.30;  // loose/dusty
const double kRegolithDensity_gcm3_high = 1.90;  // more compacted

// ---------------------------------------------------------------------------
// GCR reference primary ions (standard reduced set used in space radiation
// transport studies -- protons and He dominate by number; C/O/Mg/Si/Fe are
// the next-most dose-relevant "HZE" (High-Z, High-Energy) contributors).
// ---------------------------------------------------------------------------
struct Ion {
  TString name;
  double  Z;
  double  A;
  double  relAbundance; // approximate relative abundance BY NUMBER in the
                         // GCR flux (raw, unnormalized -- normalized to sum
                         // to 1 over just this 7-ion subset at run time).
};

// Approximate GCR relative abundances by number, representative of the
// integral GCR spectrum above roughly a few hundred MeV/nucleon (consistent
// with standard cosmic-ray composition summaries, e.g. Simpson (1983) and
// similar reviews used in space-radiation references): protons and He
// dominate by number (~87% / ~12%), with all heavier nuclei combined
// contributing only ~1%. Within the heavies, C and O are the most abundant,
// and Fe shows the well-known nucleosynthetic "iron peak" -- comparable to
// or exceeding lighter heavies like Si despite its much larger mass. These
// are representative/order-of-magnitude values -- replace with a specific
// GCR model (e.g. Badhwar-O'Neill) if you need to match a particular
// reference or solar-modulation epoch.
std::vector<Ion> gcrIons = {
  {"H",  1,  1, 87.0},
  {"He", 2,  4, 12.0},
  {"C",  6, 12,  0.16},
  {"O",  8, 16,  0.15},
  {"Mg",12, 24,  0.04},
  {"Si",14, 28,  0.03},
  {"Fe",26, 56,  0.02},
};

// ---------------------------------------------------------------------------
// Mars regolith elemental composition, derived from oxide wt% via molar-mass
// ratios (see header comment for the oxide table this comes from):
//   Si : from SiO2 (45%) * (28.09/60.08)     = 21.0 %
//   Fe : from FeO  (17%) * (55.85/71.85)     = 13.2 %
//   Al : from Al2O3(9%)  * (2*26.98/101.96)  =  4.8 %
//   Ca : from CaO  (7%)  * (40.08/56.08)     =  5.0 %
//   Mg : from MgO  (8%)  * (24.31/40.30)     =  4.8 %
//   S  : from SO3  (6%)  * (32.07/80.06)     =  2.4 %
//   Na : from Na2O (2%)  * (2*22.99/61.98)   =  1.5 %
//   K  : from K2O  (0.5%)* (2*39.10/94.20)   =  0.4 %
//   Ti : from TiO2 (1%)  * (47.87/79.87)     =  0.6 %
//   Cl : direct estimate                      =  0.5 %
//   O  : remainder from all the above oxides =  41.8 %
//   Other (trace elements, <A> ~ 23)          =  4.0 %
//  (fractions sum to ~100%, rounding)
// ---------------------------------------------------------------------------
struct TargetElement {
  TString name;
  double  A;          // mass number (approx, dominant isotope)
  double  massFrac;   // mass fraction (0-1) in bulk regolith
};

std::vector<TargetElement> regolithComposition = {
  {"O",    16, 0.418},
  {"Si",   28, 0.210},
  {"Fe",   56, 0.132},
  {"Ca",   40, 0.050},
  {"Al",   27, 0.048},
  {"Mg",   24, 0.048},
  {"S",    32, 0.024},
  {"Na",   23, 0.015},
  {"Ti",   48, 0.006},
  {"Cl",   35, 0.005},
  {"K",    39, 0.004},
  {"Other",23, 0.040},
};

// ---------------------------------------------------------------------------
// Bradt-Peters geometric reaction cross section [mb]
// ---------------------------------------------------------------------------
double sigmaReaction_mb(double Ap, double At){
  double r = kR0_fm * ( std::pow(Ap,1.0/3.0) + std::pow(At,1.0/3.0) - kB0_fm );
  if (r < 0) r = 0; // guard against unphysical (very light-on-light) inputs
  double sigma_fm2 = TMath::Pi() * r * r;
  return sigma_fm2 * kFm2_to_mb; // mb
}

// ---------------------------------------------------------------------------
// Nuclear interaction mean free path in Mars regolith, in g/cm^2, for a
// projectile of mass number Ap, summed over all target elements weighted
// by their mass fraction:
//
//   1/lambda [cm^2/g] = Sum_i ( w_i * N_A * sigma_i(Ap,At_i) / At_i )
//
// sigma_i is converted from mb to cm^2 (1 mb = 1e-27 cm^2).
// ---------------------------------------------------------------------------
double interactionLength_gcm2(double Ap){
  const double mb_to_cm2 = 1.0e-27;
  double invLambda = 0.0; // cm^2/g
  for (const auto& elem : regolithComposition){
    double sigma_mb  = sigmaReaction_mb(Ap, elem.A);
    double sigma_cm2 = sigma_mb * mb_to_cm2;
    invLambda += elem.massFrac * kAvogadro * sigma_cm2 / elem.A;
  }
  return 1.0/invLambda; // g/cm^2
}

// ---------------------------------------------------------------------------
// Main macro
// ---------------------------------------------------------------------------
void GCR_InteractionLength_MarsRegolith(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  std::cout << "============================================================\n";
  std::cout << " GCR nuclear interaction lengths in Mars regolith\n";
  std::cout << " (Bradt-Peters geometric cross section; nominal regolith\n";
  std::cout << "  density = " << kRegolithDensity_gcm3 << " g/cm^3)\n";
  std::cout << "============================================================\n";
  std::cout << std::left  << std::setw(6)  << "Ion"
            << std::right << std::setw(10) << "A"
            << std::setw(16) << "lambda[g/cm2]"
            << std::setw(14) << "lambda[cm]"
            << std::setw(16) << "lambda_low[cm]"
            << std::setw(16) << "lambda_high[cm]" << "\n";

  std::vector<double> lambda_cm;
  std::vector<double> lambda_gcm2_all;

  for (const auto& ion : gcrIons){
    double lam_gcm2    = interactionLength_gcm2(ion.A);
    double lam_cm      = lam_gcm2 / kRegolithDensity_gcm3;
    double lam_cm_low  = lam_gcm2 / kRegolithDensity_gcm3_high; // denser -> shorter physical length
    double lam_cm_high = lam_gcm2 / kRegolithDensity_gcm3_low;  // looser -> longer physical length

    lambda_cm.push_back(lam_cm);
    lambda_gcm2_all.push_back(lam_gcm2);

    std::cout << std::left  << std::setw(6)  << ion.name
              << std::right << std::setw(10) << ion.A
              << std::setw(16) << std::fixed << std::setprecision(2) << lam_gcm2
              << std::setw(14) << std::fixed << std::setprecision(2) << lam_cm
              << std::setw(16) << std::fixed << std::setprecision(2) << lam_cm_low
              << std::setw(16) << std::fixed << std::setprecision(2) << lam_cm_high
              << "\n";
  }

  // Concrete numbers to quote directly in a written justification: fraction
  // of each ion species that has undergone >=1 nuclear interaction by a
  // given near-surface regolith depth.
  std::vector<double> refDepths_cm = {10, 30, 50, 100};
  std::cout << "\nFraction of primaries that have undergone >=1 nuclear\n"
            << "interaction by a given regolith depth (nominal density):\n";
  std::cout << std::left << std::setw(6) << "Ion";
  for (double d : refDepths_cm)
    std::cout << std::right << std::setw(12) << (TString::Format("%.0f cm",d)).Data();
  std::cout << "\n";
  for (size_t i = 0; i < gcrIons.size(); ++i){
    std::cout << std::left << std::setw(6) << gcrIons[i].name;
    for (double d : refDepths_cm){
      double frac_interacted = 1.0 - std::exp(-d/lambda_cm[i]);
      std::cout << std::right << std::setw(12) << std::fixed << std::setprecision(3) << frac_interacted;
    }
    std::cout << "\n";
  }
  std::cout << "============================================================\n\n";

  // -------------------------------------------------------------------
  // GCR-abundance-weighted "system" interaction length.
  //
  // This combines the per-species mean free paths the same way the
  // regolith's own mixed-target mean free path was combined above: the
  // effective interaction RATE per unit areal density is the abundance-
  // weighted average of each species' own rate,
  //
  //   1/lambda_system [cm^2/g] = Sum_i ( f_i * (1/lambda_i) )
  //
  // where f_i = relAbundance_i / Sum(relAbundance), i.e. each ion's share
  // of the (7-species) GCR flux BY NUMBER.
  //
  // IMPORTANT CAVEAT: because protons so overwhelmingly dominate the GCR
  // flux by number (~87%), this number-weighted average is pulled strongly
  // toward the (long) proton interaction length, even though the rare
  // heavy ions interact far more readily per particle and are the ones
  // driving disproportionate secondary/fragment production. A "system"
  // interaction length weighted by number is therefore NOT the same
  // statement as "how much of the radiation field is secondaries" -- for
  // that, look at the per-species survival-fraction curves (Canvas 2) and
  // the "fraction interacted by depth" table above, especially for Fe/Si/Mg.
  // -------------------------------------------------------------------
  double sumAbundance = 0.0;
  for (const auto& ion : gcrIons) sumAbundance += ion.relAbundance;

  double invLambdaSys_cm2g = 0.0;
  for (size_t i = 0; i < gcrIons.size(); ++i){
    double fi = gcrIons[i].relAbundance / sumAbundance;
    invLambdaSys_cm2g += fi / lambda_gcm2_all[i];
  }
  double lambdaSys_gcm2 = 1.0 / invLambdaSys_cm2g;
  double lambdaSys_cm   = lambdaSys_gcm2 / kRegolithDensity_gcm3;

  std::cout << "GCR number-abundance-weighted SYSTEM interaction length:\n";
  std::cout << "  lambda_system = " << std::fixed << std::setprecision(2)
            << lambdaSys_gcm2 << " g/cm^2  ("
            << lambdaSys_cm << " cm at nominal density "
            << kRegolithDensity_gcm3 << " g/cm^3)\n";
  std::cout << "  (dominated by protons, ~"
            << std::setprecision(1) << (100.0*gcrIons[0].relAbundance/sumAbundance)
            << "% of this flux by number -- see caveat in the code comments\n"
            << "   above the calculation; this is NOT a dose- or fragmentation-\n"
            << "   weighted average.)\n";
  std::cout << "============================================================\n\n";

  // -------------------------------------------------------------------
  // Canvas 1: interaction length (cm) per species, with a shaded band
  // for a typical near-surface regolith shielding-depth range.
  // -------------------------------------------------------------------
  int nIon = (int)gcrIons.size();
  TGraph* grLen = new TGraph(nIon);
  for (int i = 0; i < nIon; ++i) grLen->SetPoint(i, i, lambda_cm[i]);
  grLen->SetTitle("Nuclear interaction length of GCR ions in Mars regolith;;Interaction length (cm)");
  grLen->SetMarkerStyle(21);
  grLen->SetMarkerSize(1.6);
  grLen->SetMarkerColor(kAzure+2);
  grLen->SetLineColor(kAzure+2);
  grLen->SetLineWidth(2);

  TCanvas* c1 = new TCanvas("c1","GCR interaction length in Mars regolith",800,600);
  grLen->Draw("APL");
  grLen->GetXaxis()->Set(nIon, -0.5, nIon-0.5);
  for (int i = 0; i < nIon; ++i)
    grLen->GetXaxis()->SetBinLabel(grLen->GetXaxis()->FindBin(i), gcrIons[i].name);
  grLen->GetXaxis()->LabelsOption("h");

  // Illustrative "typical near-surface shielding depth of interest" band --
  // adjust kShieldDepthMin/Max_cm to your specific scenario (e.g. a habitat
  // regolith berm thickness).
  double kShieldDepthMin_cm = 30.0;
  double kShieldDepthMax_cm = 100.0;
  TBox* band = new TBox(-0.5, kShieldDepthMin_cm, nIon-0.5, kShieldDepthMax_cm);
  band->SetFillColorAlpha(kGray+1, 0.35);
  band->SetLineColor(kGray+1);
  band->Draw("same");

  TLatex* bandLabel = new TLatex(nIon-0.6, kShieldDepthMax_cm*1.03, "typical shielding depth range");
  bandLabel->SetTextSize(0.030);
  bandLabel->SetTextAlign(31);
  bandLabel->Draw();

  TLine* sysLine = new TLine(-0.5, lambdaSys_cm, nIon-0.5, lambdaSys_cm);
  sysLine->SetLineColor(kRed+1);
  sysLine->SetLineStyle(2);
  sysLine->SetLineWidth(2);
  sysLine->Draw();
  TLatex* sysLabel = new TLatex(-0.4, lambdaSys_cm*1.06, "number-weighted system average");
  sysLabel->SetTextSize(0.028);
  sysLabel->SetTextColor(kRed+1);
  sysLabel->Draw();

  c1->SetLogy();
  c1->SetGridy();
  c1->Update();

  // -------------------------------------------------------------------
  // Canvas 2: survival fraction (NOT yet interacted) vs depth, per ion.
  // (1 - this) is the fraction that HAS produced secondaries/fragments
  // by that depth -- the direct evidence for secondary/fragmentation-
  // driven radiation increase.
  // -------------------------------------------------------------------
  const int nDepth = 200;
  double maxDepth_cm = 150.0;

  TCanvas* c2 = new TCanvas("c2","Survival fraction vs depth",800,600);
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle("Fraction of GCR primaries NOT yet nuclear-interacted;Depth in Mars regolith (cm);Surviving primary fraction");

  TLegend* leg = new TLegend(0.65,0.55,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  int colors[] = {kBlack, kRed+1, kOrange+1, kGreen+2, kCyan+2, kAzure+2, kMagenta+2};

  for (int i = 0; i < nIon; ++i){
    TGraph* g = new TGraph(nDepth);
    for (int j = 0; j < nDepth; ++j){
      double depth = maxDepth_cm * (j+1) / nDepth;
      double surv  = std::exp(-depth/lambda_cm[i]);
      g->SetPoint(j, depth, surv);
    }
    g->SetLineColor(colors[i % 7]);
    g->SetLineWidth(2);
    mg->Add(g, "L");
    leg->AddEntry(g, gcrIons[i].name, "l");
  }

  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0,1);
  leg->Draw();

  TLine* half = new TLine(0,0.5,maxDepth_cm,0.5);
  half->SetLineStyle(2);
  half->SetLineColor(kGray+2);
  half->Draw();

  c2->SetGridx();
  c2->SetGridy();
  c2->Update();

  c1->SaveAs("GCR_InteractionLength_MarsRegolith.png");
  c2->SaveAs("GCR_SurvivalFraction_vs_Depth_MarsRegolith.png");

  std::cout << "Saved plots:\n"
            << "  GCR_InteractionLength_MarsRegolith.png\n"
            << "  GCR_SurvivalFraction_vs_Depth_MarsRegolith.png\n";
}
