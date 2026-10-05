// ExtractNetProtonStoppingFraction_OO200.C -- lightflavorspectra_OO200
//
// PURPOSE: data-driven complement to dose/EstimateNpartNspecFractions_OO200.py's
// GEOMETRIC (Glauber Monte Carlo) participant/spectator nucleon split. That script
// answers "how many nucleons overlapped," but a Glauber MC has no fragmentation/
// de-excitation physics for the spectator remnant, and doesn't tell us how BARYON
// NUMBER (i.e. dose-carrying baryons) is actually distributed in rapidity after the
// collision. This macro instead reads this analysis's own measured net-proton
// dN/dy(y) -- p minus pbar, per Andrew's explicit choice -- and uses it as a
// DATA-DRIVEN proxy for the "stopped" (participant/fireball-associated) vs.
// "un-stopped" (spectator-associated, i.e. still moving near beam rapidity)
// baryon split. This is directly the p/pbar analogue of GCR fragmentation
// physics for spacecraft shielding: how much of the incident baryon number gets
// redistributed to mid-rapidity (dose deposited locally) vs. continues forward
// near the beam/GCR-primary rapidity (dose deposited downstream).
//
// WHY NET-PROTON (p - pbar) AND NOT RAW PROTON: at 200 GeV, pbar/p pair production
// at midrapidity is not negligible, so raw proton dN/dy conflates the ORIGINAL
// nucleons (from the two O-16 beams) with newly pair-produced protons. Net-proton
// cancels the pair-produced component (to the extent p and pbar are produced in
// equal numbers, standard assumption) and isolates transported/stopped baryon
// number. This gives a hard validation constraint: integrated over ALL rapidity,
// net-proton yield must equal exactly 16 (8 protons/nucleus x 2 nuclei), independent
// of centrality, since baryon (here: proton) number is conserved in each event. See
// the sum-rule check below.
//
// WHY A DOUBLE/TRIPLE-GAUSSIAN EXTRAPOLATION (Andrew's explicit choice): this
// analysis's measured rapidity acceptance is only |y| < 1.55 (rapidityMin_CC/
// rapidityMax_CC, headers/globalDefinitions.h), while beam rapidity at 200 GeV is
// y_beam = ln(sqrt(s_NN)/m_p) = ln(200/0.938) ~ 5.36 -- i.e. the actual
// spectator-associated peak (expected near/approaching y_beam at these energies,
// per RHIC systematics of incomplete baryon stopping at top energy) is far outside
// the measured window. A functional form is needed to extrapolate the measured
// shoulder out to y_beam. Two symmetric forms are fit and compared here:
//   - "double" Gaussian:  two mirror-symmetric Gaussians centered at +-y0 (shared
//     amplitude/width by y -> -y symmetry of a symmetric A+A system). No separate
//     y=0 component -- if there IS a resolvable central excess, this form will fit
//     it poorly and lose to the triple-Gaussian below.
//   - "triple" Gaussian: the double-Gaussian side pair PLUS a third Gaussian
//     centered at y=0 for the participant/fireball-"stopped" net-baryon yield.
// The macro fits BOTH per centrality bin and picks the preferred model via a
// standard nested-model F-test on the fit chi2/ndf (see selectPreferredModel()).
//
// IMPORTANT CAVEAT -- READ BEFORE TRUSTING THE EXTRAPOLATED FRACTIONS: the fit is
// constrained by data over |y|<1.55, a small fraction of the 0-to-y_beam=5.36
// range being integrated over. The side-Gaussian centroid y0 and width, which
// mostly control the extrapolated spectator-associated yield, are only weakly
// constrained by a shoulder (if any) that may not even be fully resolved within
// this acceptance. Parameter limits below keep the fit in a physically sensible
// regime (0 <= y0 <= y_beam), but the resulting participant/spectator SPLIT should
// be read as a first-pass, largely-extrapolated estimate, with the sum-rule check
// (integral == 16) as the main quantitative sanity check available until wider
// rapidity acceptance (e.g. published STAR forward data, or a transport-model
// cross-check) is in hand. This mirrors the same caveat already carried in
// EstimateNpartNspecFractions_OO200.py for its centrality proxy: a genuine
// first-pass estimate, not yet a calibrated final result.
//
// A SECOND CAVEAT ON THE SUM RULE: strictly, only total baryon number (protons +
// neutrons) is exactly conserved event-by-event; net-proton number specifically can
// in principle be shifted by isospin/charge-exchange processes (p <-> n via charged
// pion exchange). For an isospin-symmetric initial state (O-16: 8p+8n each beam)
// this effect is expected to average out to leading order across an ensemble, and
// net-proton sum-rule checks of this kind are standard practice in the heavy-ion
// literature -- but it is not an exact identity at the percent level, so treat
// deviations of order a few percent from 16 as expected, not necessarily a fit
// failure.
//
// WHY THIS IS A ROOT C++ MACRO, NOT PYTHON: this repo's real dN/dy(y) output lives
// in ROOT TH1D histograms inside a ROOT-format fit_output.root (SpectraFitter's
// writeOutputs(), see below for exact paths/naming). This cloud sandbox has no ROOT
// installed, no network route to install it (apt-get is blocked here), and no
// pip-installable ROOT-independent reader either -- `pip install uproot` fails with
// "No matching distribution found" here, same as the `pip install root/pyroot`
// attempts already documented in EstimateNpartNspecFractions_OO200.py's commit
// history. So, unlike that script, this one could NOT be reimplemented and
// numerically validated in pure Python in this session. It is written directly as
// a ROOT macro, following this repo's existing macro conventions (see
// macros/RunSpectraFitter.C), and should be run and checked by Andrew wherever
// ROOT + a real fit_output.root are available (RCF, or a local ROOT install) --
// same pattern already used for this session's earlier RCF-targeted C++ edits.
//
// INPUT: a_fitOutputFile must be the ROOT file written by SpectraFitter::
// writeOutputs() (i.e. the a_outputFile argument of macros/RunSpectraFitter.C,
// commonly "fit_output.root"). Required histograms, read per centrality bin
// centIndex in [0, a_numCent):
//   ParticleYields/ProtonPlus/dNdy_ProtonPlus_Cent%02d_Nominal   (proton dN/dy(y))
//   ParticleYields/ProtonMinus/dNdy_ProtonMinus_Cent%02d_Nominal (antiproton dN/dy(y))
// (naming confirmed against SpectraFitter.cxx's fitSingleRapiditySpectra()/
// writeOutputs(), and ParticleInfo::GetParticleName(2,+-1) = "ProtonPlus"/
// "ProtonMinus" for partIndex=2="Proton"; bin content AND bin error are both
// explicitly SetBinContent/SetBinError upstream, so TH1::Add's automatic
// quadrature error propagation for the net-proton subtraction below is correct).
//
// OUTPUT (written under a_outputDir):
//   NetProtonStopping_Summary.csv       -- one row per centrality bin: chi2/ndf for
//     both models, which model was selected, integrated net-proton yield (fit) vs.
//     the 16 sum-rule target, and the participant/fireball- vs. spectator-associated
//     fractional split.
//   NetProtonStopping_Cent%02d_Fit.png  -- per-centrality diagnostic plot: measured
//     net-proton dN/dy(y) points, selected fit curve over the measured range (solid)
//     and extrapolated out to +-y_beam (dashed), with individual Gaussian components
//     overlaid when the triple model is selected.
//   NetProtonStopping_FractionSummary.png -- fireball/participant vs.
//     spectator-associated fraction vs. centrality, same visual style as
//     OO200_Npart_Nspec_byCentrality_WoodsSaxon.png for direct side-by-side
//     comparison with the geometric Glauber result.
//
// USAGE (ROOT, at RCF or anywhere with ROOT + a real fit_output.root):
//   root -l -b -q 'ExtractNetProtonStoppingFraction_OO200.C("fit_output.root")'
// or with explicit output dir / centrality count:
//   root -l -b -q 'ExtractNetProtonStoppingFraction_OO200.C("fit_output.root","./netproton_stopping_output",6)'
//
// NEXT STEPS once this has been run against real data:
//   1. Cross-check the sum-rule integral (should be ~16) as the primary sanity
//      check on the extrapolation before trusting the fractional split.
//   2. Compare the spectator-associated fraction here against the GEOMETRIC
//      Nspec/(2*A) fraction from EstimateNpartNspecFractions_OO200.py's Woods-Saxon
//      Glauber result, bin-by-bin -- they answer related but distinct questions
//      (geometric nucleon count vs. data-driven baryon-number transport) and should
//      be broadly consistent in trend (more peripheral -> larger spectator share)
//      even if not numerically identical.
//   3. Once real forward-rapidity data (STAR EPD/ZDC-adjacent, or a published
//      wider-acceptance net-proton measurement) is available, use it to directly
//      constrain y0/sigma of the side Gaussians instead of leaving them almost
//      entirely extrapolation-driven.
//   4. Fold in the "fireball-produced particles" (thermal secondary pi/K/p, NOT the
//      transported nucleons this macro targets) as the third dose source named on
//      the WG slide -- that is a separate multiplicity estimate (e.g. from measured
//      mid-rapidity pi/K yields, or GlauberClass::ProduceParticles()'s two-component
//      model), not something net-proton dN/dy(y) alone can give.
//   5. Turn the resulting participant/spectator-associated BARYON fractions into
//      actual fractional DOSE (this macro stops at baryon-number bookkeeping) --
//      needs a dose-per-baryon or LET/fluence-to-dose conversion, same open item
//      already flagged in EstimateNpartNspecFractions_OO200.py.

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TColor.h"
#include "TMatrixDSym.h"

using namespace std;

//============================================================================
// Physics constants
//============================================================================
const double M_PROTON       = 0.938272; // GeV, PDG proton mass
const double SQRT_S_NN      = 200.0;    // GeV, this analysis's collision energy
const int    NET_PROTON_SUM_RULE = 16;  // O+O: 8 protons/nucleus x 2 nuclei; see
                                         // header caveat on exactness

const double Y_ACCEPT_MIN = -1.55; // matches globalDefinitions.h rapidityMin_CC
const double Y_ACCEPT_MAX =  1.55; // matches globalDefinitions.h rapidityMax_CC

const int N_CENT_BINS_DEFAULT = 6;
const char* CENT_LABELS[6] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};

// "Midnight Executive" palette, matching this repo's existing WG-deck figures
// (dose/OO200_Npart_Nspec_byCentrality_WoodsSaxon.png) for direct visual comparison.
const int COLOR_PARTICIPANT = TColor::GetColor("#1E2761"); // navy
const int COLOR_SPECTATOR   = TColor::GetColor("#E8A33D"); // amber
const int COLOR_DATA        = kBlack;
const int COLOR_FIT_TOTAL   = TColor::GetColor("#C0392B"); // brick red, distinct from both fractions

//============================================================================
// Fit functional forms
//============================================================================

// "Double" model: two mirror-symmetric Gaussians at +-y0, shared amplitude/width.
// par[0] = A     (shared amplitude of each side Gaussian)
// par[1] = y0    (>=0, centroid of the +y side Gaussian; -y0 for the other)
// par[2] = sigma (shared width)
double DoubleGaussSymmetric(double* x, double* par){
  double y = x[0];
  double A = par[0], y0 = par[1], sigma = par[2];
  if(sigma <= 0) return 0;
  double gPlus  = A * TMath::Gaus(y, +y0, sigma, kFALSE);
  double gMinus = A * TMath::Gaus(y, -y0, sigma, kFALSE);
  return gPlus + gMinus;
}

// "Triple" model: DoubleGaussSymmetric's side pair PLUS a central Gaussian at y=0
// for the participant/fireball-"stopped" net-proton yield.
// par[0] = Ac     (central-Gaussian amplitude, y=0)
// par[1] = sigmaC (central-Gaussian width)
// par[2] = As     (shared amplitude of each side Gaussian)
// par[3] = y0     (>=0, side-Gaussian centroid)
// par[4] = sigmaS (shared side-Gaussian width)
double TripleGaussSymmetric(double* x, double* par){
  double y = x[0];
  double Ac = par[0], sigmaC = par[1];
  double As = par[2], y0 = par[3], sigmaS = par[4];
  if(sigmaC <= 0 || sigmaS <= 0) return 0;
  double central = Ac * TMath::Gaus(y, 0.0, sigmaC, kFALSE);
  double gPlus   = As * TMath::Gaus(y, +y0, sigmaS, kFALSE);
  double gMinus  = As * TMath::Gaus(y, -y0, sigmaS, kFALSE);
  return central + gPlus + gMinus;
}

// Individual-component helpers, used only for drawing the decomposed curves and
// for splitting the integrated yield into participant- vs spectator-associated
// pieces after a fit converges (par[] passed in verbatim from the fit result).
double CentralComponentOnly(double* x, double* par){
  // par[0]=Ac, par[1]=sigmaC (same indexing as the first two TripleGaussSymmetric pars)
  double y = x[0];
  double Ac = par[0], sigmaC = par[1];
  if(sigmaC <= 0) return 0;
  return Ac * TMath::Gaus(y, 0.0, sigmaC, kFALSE);
}
double SideComponentOnly(double* x, double* par){
  // par[0]=As, par[1]=y0, par[2]=sigmaS
  double y = x[0];
  double As = par[0], y0 = par[1], sigmaS = par[2];
  if(sigmaS <= 0) return 0;
  return As * TMath::Gaus(y, +y0, sigmaS, kFALSE) + As * TMath::Gaus(y, -y0, sigmaS, kFALSE);
}

//============================================================================
// Per-centrality-bin result bundle
//============================================================================
struct StoppingFractionResult {
  int    centIndex;
  bool   ok;              // false if histograms were missing/empty -- skip in summary
  bool   usedTripleModel;
  double chi2Double, ndfDouble;
  double chi2Triple, ndfTriple;
  double pValueFTest;
  double integralTotal, integralTotalErr;      // fit integral over [-y_beam,+y_beam]
  double fracParticipant, fracParticipantErr;   // central-component share of integralTotal
  double fracSpectator,  fracSpectatorErr;      // side-component share of integralTotal
};

//============================================================================
// selectPreferredModel: standard nested-model F-test between the double (3-param)
// and triple (5-param) fits. Returns true if the triple model is preferred.
//============================================================================
bool selectPreferredModel(double chi2Double, int ndfDouble,
                           double chi2Triple, int ndfTriple,
                           double centralAmplitude, double& pValueOut){
  pValueOut = 1.0;
  int deltaNdf = ndfDouble - ndfTriple; // extra parameters used by the triple model (nominally 2)
  if(deltaNdf <= 0 || ndfTriple <= 0 || chi2Triple <= 0) return false;
  if(centralAmplitude <= 0) return false; // unphysical/non-converged central component
  double fStat = ((chi2Double - chi2Triple) / deltaNdf) / (chi2Triple / ndfTriple);
  if(fStat <= 0) return false;
  pValueOut = 1.0 - TMath::FDistI(fStat, deltaNdf, ndfTriple);
  return (pValueOut < 0.05);
}

//============================================================================
// Main macro
//============================================================================
void ExtractNetProtonStoppingFraction_OO200(string a_fitOutputFile,
                                             string a_outputDir = "./netproton_stopping_output",
                                             int    a_numCent = N_CENT_BINS_DEFAULT){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gSystem->mkdir(a_outputDir.c_str(), kTRUE);

  double yBeam = TMath::Log(SQRT_S_NN / M_PROTON);
  cout << "==================================================================" << endl;
  cout << " ExtractNetProtonStoppingFraction_OO200" << endl;
  cout << " sqrt(s_NN)       = " << SQRT_S_NN << " GeV" << endl;
  cout << " y_beam           = " << yBeam << endl;
  cout << " Measured acceptance: [" << Y_ACCEPT_MIN << ", " << Y_ACCEPT_MAX << "]" << endl;
  cout << " Net-proton sum rule target = " << NET_PROTON_SUM_RULE << endl;
  cout << "==================================================================" << endl;

  TFile* inFile = TFile::Open(a_fitOutputFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){
    cerr << "ERROR: could not open input file '" << a_fitOutputFile << "'" << endl;
    return;
  }

  vector<StoppingFractionResult> results;
  ofstream csv((a_outputDir + "/NetProtonStopping_Summary.csv").c_str());
  csv << "centIndex,centLabel,model,chi2ndfDouble,chi2ndfTriple,pValueFTest,"
      << "integralTotal,integralTotalErr,fracParticipant,fracParticipantErr,"
      << "fracSpectator,fracSpectatorErr\n";

  for(int centIndex = 0; centIndex < a_numCent; centIndex++){

    StoppingFractionResult res;
    res.centIndex = centIndex;
    res.ok = false;
    res.usedTripleModel = false;

    TString hNameProton    = Form("ParticleYields/ProtonPlus/dNdy_ProtonPlus_Cent%02d_Nominal", centIndex);
    TString hNameAntiproton = Form("ParticleYields/ProtonMinus/dNdy_ProtonMinus_Cent%02d_Nominal", centIndex);

    TH1D* hProton     = (TH1D*)inFile->Get(hNameProton);
    TH1D* hAntiproton = (TH1D*)inFile->Get(hNameAntiproton);

    if(!hProton || !hAntiproton){
      cerr << "WARNING: centIndex " << centIndex << " (" << CENT_LABELS[centIndex]
           << ") -- missing '" << (!hProton ? hNameProton : hNameAntiproton)
           << "', skipping this bin." << endl;
      results.push_back(res);
      continue;
    }

    TH1D* hNet = (TH1D*)hProton->Clone(Form("dNdy_NetProton_Cent%02d", centIndex));
    hNet->SetDirectory(0);
    hNet->Add(hAntiproton, -1.0); // errors combine in quadrature automatically;
                                  // valid because bin errors are explicitly set
                                  // upstream in SpectraFitter::fitSingleRapiditySpectra
    hNet->SetTitle(Form("Net-Proton dN/dy -- %s;Rapidity y;dN/dy (Net Proton)", CENT_LABELS[centIndex]));

    if(hNet->Integral() <= 0){
      cerr << "WARNING: centIndex " << centIndex << " (" << CENT_LABELS[centIndex]
           << ") -- net-proton histogram is empty/non-positive, skipping." << endl;
      results.push_back(res);
      continue;
    }

    double hMax = hNet->GetMaximum();
    double centralBinContent = hNet->GetBinContent(hNet->FindBin(0.0));

    // ---- Double-Gaussian fit ----
    TF1* fDouble = new TF1(Form("fDouble_Cent%02d", centIndex), DoubleGaussSymmetric,
                            Y_ACCEPT_MIN, Y_ACCEPT_MAX, 3);
    fDouble->SetParameters(0.5 * hMax, 1.2, 1.0);
    fDouble->SetParLimits(0, 0.0, 10.0 * hMax);       // A >= 0
    fDouble->SetParLimits(1, 0.0, yBeam);              // 0 <= y0 <= y_beam
    fDouble->SetParLimits(2, 0.05, 2.0 * yBeam);       // sigma > 0, physically bounded
    TFitResultPtr fitDoubleResult = hNet->Fit(fDouble, "RSQ0"); // R=range, S=save, Q=quiet, 0=no auto-draw

    // ---- Triple-Gaussian fit ----
    TF1* fTriple = new TF1(Form("fTriple_Cent%02d", centIndex), TripleGaussSymmetric,
                            Y_ACCEPT_MIN, Y_ACCEPT_MAX, 5);
    fTriple->SetParameters(TMath::Max(centralBinContent, 0.1 * hMax), 0.8,
                            0.3 * hMax, 1.3, 1.0);
    fTriple->SetParLimits(0, 0.0, 10.0 * hMax);        // Ac >= 0
    fTriple->SetParLimits(1, 0.05, 2.0 * yBeam);       // sigmaC > 0
    fTriple->SetParLimits(2, 0.0, 10.0 * hMax);        // As >= 0
    fTriple->SetParLimits(3, 0.0, yBeam);              // 0 <= y0 <= y_beam
    fTriple->SetParLimits(4, 0.05, 2.0 * yBeam);       // sigmaS > 0
    TFitResultPtr fitTripleResult = hNet->Fit(fTriple, "RSQ0");

    res.chi2Double = fDouble->GetChisquare();
    res.ndfDouble  = fDouble->GetNDF();
    res.chi2Triple = fTriple->GetChisquare();
    res.ndfTriple  = fTriple->GetNDF();

    double pValue = 1.0;
    bool useTriple = selectPreferredModel(res.chi2Double, (int)res.ndfDouble,
                                           res.chi2Triple, (int)res.ndfTriple,
                                           fTriple->GetParameter(0), pValue);
    res.pValueFTest = pValue;
    res.usedTripleModel = useTriple;

    TF1* fSelected = useTriple ? fTriple : fDouble;

    // ---- Integrate the SELECTED model over the full rapidity range to y_beam ----
    // (Gaussian tails beyond y_beam are physically meaningless -- baryon number
    // cannot exceed beam rapidity -- so the integral is bounded there rather than
    // taken to +-infinity.)
    //
    // IMPORTANT: TF1::IntegralError() with no explicit params/covariance argument
    // silently reads the LAST fit's global fitter state, which by this point in
    // the loop is whichever of fDouble/fTriple was fit most recently -- NOT
    // necessarily the model actually selected. To avoid that stale-state bug,
    // pass each model's own TFitResultPtr parameters/covariance explicitly, so
    // the error is always computed from the correct fit regardless of call order.
    double integralDouble    = fDouble->Integral(-yBeam, yBeam);
    double integralDoubleErr = fDouble->IntegralError(-yBeam, yBeam,
                                  fitDoubleResult->GetParams(),
                                  fitDoubleResult->GetCovarianceMatrix().GetMatrixArray());
    double integralTriple    = fTriple->Integral(-yBeam, yBeam);
    double integralTripleErr = fTriple->IntegralError(-yBeam, yBeam,
                                  fitTripleResult->GetParams(),
                                  fitTripleResult->GetCovarianceMatrix().GetMatrixArray());

    double integralTotal = useTriple ? integralTriple : integralDouble;
    double integralErr   = useTriple ? integralTripleErr : integralDoubleErr;

    res.integralTotal = integralTotal;
    res.integralTotalErr = integralErr;

    // ---- Split into participant/fireball- vs spectator-associated fractions ----
    double integralCentral = 0, integralSide = 0;
    double integralCentralErr = 0, integralSideErr = 0;

    if(useTriple){
      TF1 fCentralOnly("fCentralOnly", CentralComponentOnly, -yBeam, yBeam, 2);
      fCentralOnly.SetParameters(fTriple->GetParameter(0), fTriple->GetParameter(1));
      TF1 fSideOnly("fSideOnly", SideComponentOnly, -yBeam, yBeam, 3);
      fSideOnly.SetParameters(fTriple->GetParameter(2), fTriple->GetParameter(3), fTriple->GetParameter(4));

      integralCentral = fCentralOnly.Integral(-yBeam, yBeam);
      integralSide     = fSideOnly.Integral(-yBeam, yBeam);
      // Approximate component errors by propagating the fractional error on the
      // total fit integral, weighted by each component's share -- the full
      // covariance-aware split would require passing fTriple's covariance matrix
      // into a combined error propagation, which ROOT's IntegralError doesn't do
      // per-component out of the box. This is a reasonable first-pass estimate;
      // treat these component errors as approximate, not exact.
      double fracErrApprox = (integralTotal > 0) ? (integralErr / integralTotal) : 0;
      integralCentralErr = integralCentral * fracErrApprox;
      integralSideErr     = integralSide * fracErrApprox;
    } else {
      // Double-Gaussian model: by construction there is no separate y=0 component
      // -- ALL of the fit yield is formally "side"/spectator-associated-shaped.
      // This does not necessarily mean zero participant-associated yield exists;
      // it means this simpler model did not resolve one as statistically
      // preferred over the triple model (see header caveat). Flag this explicitly
      // rather than silently reporting 0% participant fraction as if it were a
      // confident physics result.
      integralCentral = 0;
      integralSide = integralTotal;
      integralCentralErr = 0;
      integralSideErr = integralErr;
    }

    res.fracParticipant = (integralTotal != 0) ? (integralCentral / integralTotal) : 0;
    res.fracSpectator    = (integralTotal != 0) ? (integralSide / integralTotal) : 0;
    res.fracParticipantErr = (integralTotal != 0) ? (integralCentralErr / integralTotal) : 0;
    res.fracSpectatorErr    = (integralTotal != 0) ? (integralSideErr / integralTotal) : 0;
    res.ok = true;
    results.push_back(res);

    double chi2ndfDouble = (res.ndfDouble > 0) ? res.chi2Double / res.ndfDouble : -1;
    double chi2ndfTriple = (res.ndfTriple > 0) ? res.chi2Triple / res.ndfTriple : -1;

    cout << "------------------------------------------------------------------" << endl;
    cout << " Cent " << CENT_LABELS[centIndex] << " (index " << centIndex << ")" << endl;
    cout << "   chi2/ndf: double=" << chi2ndfDouble << "  triple=" << chi2ndfTriple
         << "  F-test p=" << pValue << "  -> using " << (useTriple ? "TRIPLE" : "DOUBLE") << " model" << endl;
    cout << "   Integrated net-proton yield [-y_beam,y_beam] = " << integralTotal
         << " +- " << integralErr << "   (sum-rule target = " << NET_PROTON_SUM_RULE
         << ", deviation = " << (100.0 * (integralTotal - NET_PROTON_SUM_RULE) / NET_PROTON_SUM_RULE)
         << "%)" << endl;
    cout << "   Participant/fireball-associated fraction = " << 100.0 * res.fracParticipant
         << "% +- " << 100.0 * res.fracParticipantErr << "%" << endl;
    cout << "   Spectator-associated fraction             = " << 100.0 * res.fracSpectator
         << "% +- " << 100.0 * res.fracSpectatorErr << "%" << endl;

    csv << centIndex << "," << CENT_LABELS[centIndex] << ","
        << (useTriple ? "triple" : "double") << ","
        << chi2ndfDouble << "," << chi2ndfTriple << "," << pValue << ","
        << integralTotal << "," << integralErr << ","
        << res.fracParticipant << "," << res.fracParticipantErr << ","
        << res.fracSpectator << "," << res.fracSpectatorErr << "\n";

    // ---- Diagnostic per-centrality plot ----
    TCanvas* c = new TCanvas(Form("cFit_Cent%02d", centIndex), "", 900, 650);
    hNet->SetMarkerStyle(20);
    hNet->SetMarkerColor(COLOR_DATA);
    hNet->SetLineColor(COLOR_DATA);
    hNet->GetXaxis()->SetRangeUser(-yBeam * 1.05, yBeam * 1.05);
    hNet->Draw("PE");

    // Solid curve over the measured acceptance, dashed extrapolation beyond it.
    TF1* fMeasuredRange = (TF1*)fSelected->Clone(Form("fMeasured_Cent%02d", centIndex));
    fMeasuredRange->SetRange(Y_ACCEPT_MIN, Y_ACCEPT_MAX);
    fMeasuredRange->SetLineColor(COLOR_FIT_TOTAL);
    fMeasuredRange->SetLineWidth(3);
    fMeasuredRange->SetLineStyle(1);
    fMeasuredRange->Draw("SAME");

    TF1* fExtrapLow = (TF1*)fSelected->Clone(Form("fExtrapLow_Cent%02d", centIndex));
    fExtrapLow->SetRange(-yBeam, Y_ACCEPT_MIN);
    fExtrapLow->SetLineColor(COLOR_FIT_TOTAL);
    fExtrapLow->SetLineWidth(2);
    fExtrapLow->SetLineStyle(2);
    fExtrapLow->Draw("SAME");

    TF1* fExtrapHigh = (TF1*)fSelected->Clone(Form("fExtrapHigh_Cent%02d", centIndex));
    fExtrapHigh->SetRange(Y_ACCEPT_MAX, yBeam);
    fExtrapHigh->SetLineColor(COLOR_FIT_TOTAL);
    fExtrapHigh->SetLineWidth(2);
    fExtrapHigh->SetLineStyle(2);
    fExtrapHigh->Draw("SAME");

    TLegend* leg = new TLegend(0.15, 0.70, 0.55, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hNet, "Net-proton (p - p#bar{})", "PE");
    leg->AddEntry(fMeasuredRange, useTriple ? "Fit (triple-Gauss, measured range)" : "Fit (double-Gauss, measured range)", "L");
    leg->AddEntry(fExtrapHigh, "Extrapolation to y_{beam}", "L");
    leg->Draw();

    TLatex tl;
    tl.SetNDC();
    tl.SetTextSize(0.035);
    tl.DrawLatex(0.15, 0.92, Form("O+O #sqrt{s_{NN}}=200 GeV, %s -- Participant %.0f%%, Spectator %.0f%%",
                                   CENT_LABELS[centIndex], 100.0 * res.fracParticipant, 100.0 * res.fracSpectator));

    c->SaveAs(Form("%s/NetProtonStopping_Cent%02d_Fit.png", a_outputDir.c_str(), centIndex));
    delete c;
  }
  csv.close();
  cout << "==================================================================" << endl;
  cout << "Wrote " << a_outputDir << "/NetProtonStopping_Summary.csv" << endl;

  // ---- Summary plot: fireball/participant vs spectator-associated fraction vs centrality ----
  vector<double> xVals, fracPart, fracPartErr, fracSpec, fracSpecErr;
  for(size_t i = 0; i < results.size(); i++){
    if(!results[i].ok) continue;
    xVals.push_back((double)results[i].centIndex);
    fracPart.push_back(100.0 * results[i].fracParticipant);
    fracPartErr.push_back(100.0 * results[i].fracParticipantErr);
    fracSpec.push_back(100.0 * results[i].fracSpectator);
    fracSpecErr.push_back(100.0 * results[i].fracSpectatorErr);
  }

  if(xVals.size() > 0){
    TCanvas* cSum = new TCanvas("cSummary", "", 900, 650);
    TGraphErrors* gPart = new TGraphErrors(xVals.size(), &xVals[0], &fracPart[0], 0, &fracPartErr[0]);
    TGraphErrors* gSpec = new TGraphErrors(xVals.size(), &xVals[0], &fracSpec[0], 0, &fracSpecErr[0]);
    gPart->SetMarkerStyle(21); gPart->SetMarkerColor(COLOR_PARTICIPANT); gPart->SetLineColor(COLOR_PARTICIPANT); gPart->SetLineWidth(2);
    gSpec->SetMarkerStyle(21); gSpec->SetMarkerColor(COLOR_SPECTATOR);   gSpec->SetLineColor(COLOR_SPECTATOR);   gSpec->SetLineWidth(2);

    gPart->SetTitle(";Centrality Bin Index (0=Most Central);Fraction of Net-Proton Yield [%]");
    gPart->GetYaxis()->SetRangeUser(0, 100);
    gPart->GetXaxis()->SetLimits(-0.5, (double)a_numCent - 0.5);
    gPart->Draw("APL");
    gSpec->Draw("PL SAME");

    TLegend* legSum = new TLegend(0.15, 0.78, 0.55, 0.90);
    legSum->SetBorderSize(0);
    legSum->SetFillStyle(0);
    legSum->AddEntry(gPart, "Participant/fireball-associated (net-proton, data-driven)", "PL");
    legSum->AddEntry(gSpec, "Spectator-associated (net-proton, data-driven)", "PL");
    legSum->Draw();

    cSum->SaveAs(Form("%s/NetProtonStopping_FractionSummary.png", a_outputDir.c_str()));
    delete cSum;
    cout << "Wrote " << a_outputDir << "/NetProtonStopping_FractionSummary.png" << endl;
  } else {
    cout << "No centrality bins produced a valid fit -- summary plot skipped." << endl;
  }

  inFile->Close();
  cout << "==================================================================" << endl;
  cout << "Done. Cross-check the sum-rule column in the CSV (target = "
       << NET_PROTON_SUM_RULE << ") before trusting the fractional split." << endl;
}
