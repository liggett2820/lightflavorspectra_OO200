// GenerateGlauberMC_OO200.C -- lightflavorspectra_OO200
//
// PURPOSE: C++/ROOT replacement for EstimateNpartNspecFractions_OO200.py, per
// Andrew's request to keep the Glauber-refMult pipeline in C++/ROOT rather
// than mixing in Python (same reasoning already applied to
// PlotDoseFractionComparison_OO200.C -- see that file's header: the Python
// version was only used once because the cloud sandbox that first wrote this
// analysis had no ROOT install). This macro reproduces the EXACT SAME physics
// as the Python script line-for-line (same source citations against
// GlauberClass.cxx, same constants), just as a native ROOT macro:
//   - GlauberClass::GetNuclearRadius()        -> R(A) = 1.25 * A^(1/3) fm
//                                                 (used ONLY for b_max)
//   - GlauberClass::GetWoodsSaxonParameters() -> O-16 3-param-Fermi params
//     {p0, w=-0.051, c=2.608 fm, a=0.513 fm, B20=B40=0} (H. De Vries et al.,
//     At. Data Nucl. Data Tables 36 (1987))
//   - GlauberClass::WoodsSaxonFunc()          -> rho(r)/p0 =
//       (1 + w*(r/c)^2 if r<c else 1+w) / (1 + exp((r-c)/a))
//     B20=B40=0 for O-16 -> spherically symmetric, no theta dependence
//     (RotateNucleons() no-ops), so this is a purely radial 3-param Fermi
//     shape, NOT a plain 2-parameter Woods-Saxon (the w term is a real
//     core-density correction that matters for a light nucleus like O-16).
//   - GlauberClass::SetImpactParameter()      -> b sampled from PDF ~ b on
//     [0, b_max] (matches TF1("2*pi*x"))
//   - GlauberClass::GenerateNucleonPositions()-> b_max = 1.3*(R_A+R_B); nucleus
//     B shifted by +b along x only
//   - GlauberClass::CountNpartNcoll()         -> two nucleons collide if their
//     TRANSVERSE (x,y only, z ignored) separation < r_eff, where
//     r_eff = sqrt(sigma_NN_mb / (10*pi)) fm (standard "black disk" test)
//   - GlauberModel.cxx                        -> only Npart>0 ("successful"/
//     reacting) events are kept -- this defines sigma_reaction
//
// SAMPLING METHOD: r is drawn via TH1D::GetRandom() off a histogram filled
// with rho(r)*r^2 (the radial density weighted by the spherical volume
// element -- same combination GlauberClass::WoodsSaxon2D_PDF_Func integrates
// against dtheta, and the same quantity the Python version tabulates into an
// inverse CDF by hand). Direction is drawn isotropically (cos(theta) uniform
// on [-1,1], phi uniform on [0,2pi]), exact given B20=B40=0 -- only the
// resulting TRANSVERSE (x,y) position is kept, since CountNpartNcoll()/the
// collision test above ignores z entirely.
//
// SCOPE NOTE (same as the Python version): this does NOT use a negative-
// binomial multiplicity model tied to real refMult -- that is exactly what
// FitGlauberNBDToRefMult_OO200.C does downstream, consuming this macro's
// per-cell CSV output. This macro is the purely-geometric Glauber MC input to
// that fit, nothing more.
//
// OUTPUT -- CHANGED from JSON+CSV to a single native ROOT file, per Andrew's
// explicit instruction to get Python and JSON entirely out of this pipeline
// (EstimateNpartNspecFractions_OO200.py -- the original source of this
// schema -- and PlotDoseFractionComparison_OO200.py are retired for the same
// reason; see their replacement macros' headers). Both former JSON/CSV
// outputs are now TTrees in one glauber_oo200_results_woodssaxon.root:
//   TTree "GlauberCentralityBins" -- one entry per centrality bin. Branches:
//     cent/C, nEvents/D, meanB_fm/D, meanNpart/D, meanNspec/D, meanNcoll/D,
//     NpartRMS/D, NspecRMS/D, fracNpartOfTotal/D, fracNspecOfTotal/D,
//     refMultEdgeHi/D (always -999 here -- this macro is purely geometric,
//     no refMult involved; see FitGlauberNBDToRefMult_OO200.C for the
//     variant that fills this field for real). THIS SCHEMA IS SHARED with
//     FitGlauberNBDToRefMult_OO200.C's own GlauberCentralityBins output (that
//     macro fills meanB_fm/NpartRMS/NspecRMS with -999 instead, since it
//     doesn't compute those) -- every downstream macro that reads a
//     GlauberCentralityBins tree can point at EITHER file interchangeably,
//     same drop-in intent the JSON schema used to serve.
//   TTree "GlauberCells" -- one entry per distinct (Npart,Ncoll) pair among
//     successful events, with its event count. Branches: Npart/I, Nspec/I,
//     Ncoll/I, weight/L. Replaces GlauberNpartNcollWeights_OO200.csv for
//     FitGlauberNBDToRefMult_OO200.C's NBD fit (same collapsed-cell idea,
//     see that macro's header for why summing weight*NBD_PMF over these
//     cells is exactly equivalent to summing over every raw MC event).
//   TTree "GlauberConfig" -- one entry, the run configuration/diagnostics
//     (model, A_nucleons, sigma_NN_mb, N_thrown, N_success, WS_w, WS_c_fm,
//     WS_a_fm, R_nuclearRadius_fm_bmax_only, b_max_fm, r_eff_fm,
//     sigma_reaction_mb). Informational only -- no downstream macro reads
//     this back (same role the JSON "config" sub-object played: printed and
//     archived, not consumed programmatically elsewhere).
// Every downstream macro that reads this -- EstimateFireballFractionFromRefMult_OO200.C,
// PlotDoseFractionComparison_OO200.C, FitGlauberNBDToRefMult_OO200.C -- has
// been updated to open this .root file and TTree::GetEntry() instead of
// hand-scanning JSON text; see each of those files' own header for its
// updated reader.
//
// PERFORMANCE: default a_nEvents=2,000,000, each doing a 16x16 transverse
// pairwise collision test (256 distance evaluations) plus 32 TH1::GetRandom()
// calls -- call this with the trailing '+' (ACLiC) or it will be slow
// interpreted:
//   root -l -b -q 'GenerateGlauberMC_OO200.C+()'
//
// USAGE:
//   root -l -b -q 'GenerateGlauberMC_OO200.C+()'
// or with explicit event count / output dir:
//   root -l -b -q 'GenerateGlauberMC_OO200.C+(2000000, ".")'
//
// STATUS: UNTESTED -- no ROOT in the sandbox that wrote this (same caveat as
// every other C++ macro in this session). The physics/constants are a direct
// line-for-line port of the already-reasoned-through Python version (see that
// script's own header for the original derivation and citations) -- the new
// risk here is purely mechanical (ROOT API usage, sampling correctness), not a
// re-derivation of the physics. Andrew HAS run an earlier (JSON+CSV-writing)
// version of this macro against real ROOT before -- glauber_oo200_results_woodssaxon.json
// on disk reflects that real run. The OUTPUT-WRITING block was just changed
// from JSON+CSV to a single .root file (3 TTrees, see header) and has NOT
// been exercised against real ROOT -- re-run before trusting
// glauber_oo200_results_woodssaxon.root or its downstream consumers.
// Compare its printed sigma_reaction and per-bin <Npart>/<Nspec> against the
// physics-numbers already validated in the old JSON as the first sanity check
// -- they should agree exactly (same event loop, only the output format
// changed), then delete the old .json/.csv once the new .root is confirmed good.

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStopwatch.h"

using namespace std;

//============================================================================
// Configuration -- matches EstimateNpartNspecFractions_OO200.py exactly.
//============================================================================
const int    A_NUCLEONS  = 16;   // oxygen-16, both beams
const int    N_CENT_BINS = 6;
const char*  CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
const double CENT_EDGES[N_CENT_BINS + 1] = {0.0, 0.05, 0.10, 0.20, 0.40, 0.80, 1.00};

// O-16 3-parameter Fermi charge-density params, verbatim from
// GlauberClass::GetWoodsSaxonParameters(16) -- {p0, w, c, a, B20, B40, gamma2, gamma4}
const double WS_W = -0.051; // central-density correction (negative = core depletion)
const double WS_C = 2.608;  // half-density radius, fm
const double WS_A = 0.513;  // diffuseness/skin depth, fm
const double R_GRID_MAX = 18.0; // fm -- matches the TF2 range GlauberClass itself uses
const int    R_GRID_NBINS = 200000; // fine-binned PDF histogram, matches Python's grid density

//============================================================================
// rho(r)/p0 -- GlauberClass::WoodsSaxonFunc() with B20=B40=0 (O-16 is
// spherically symmetric), i.e. no theta dependence.
//============================================================================
double WoodsSaxonShape(double r){
  double numerator = (r < WS_C) ? (1.0 + WS_W * (r / WS_C) * (r / WS_C)) : (1.0 + WS_W);
  double denominator = 1.0 + TMath::Exp((r - WS_C) / WS_A);
  return numerator / denominator;
}

//============================================================================
// Sample one nucleon's TRANSVERSE (x,y) position: r ~ rho(r)*r^2 via the
// histogram's built-in inverse-CDF sampling, direction isotropic. z is not
// computed at all -- CountNpartNcoll()'s transverse-only collision test never
// uses it (see header).
//============================================================================
inline void SampleNucleonXY(TH1D* a_hR, double& a_x, double& a_y){
  double r = a_hR->GetRandom();
  double cosTheta = gRandom->Uniform(-1.0, 1.0);
  double sinTheta = TMath::Sqrt(1.0 - cosTheta * cosTheta);
  double phi = gRandom->Uniform(0.0, 2.0 * TMath::Pi());
  double rhoPerp = r * sinTheta;
  a_x = rhoPerp * TMath::Cos(phi);
  a_y = rhoPerp * TMath::Sin(phi);
}

//============================================================================
// Main macro
//============================================================================
void GenerateGlauberMC_OO200(
    Long64_t a_nEvents     = 2000000,
    string   a_outputDir   = ".",
    double   a_sigmaNN_mb  = 42.0,   // NN inelastic cross section, matches this repo's OO200 usage
    UInt_t   a_seed        = 12345
){
  gSystem->mkdir(a_outputDir.c_str(), kTRUE);
  gRandom->SetSeed(a_seed);

  double R_nucRad = 1.25 * TMath::Power((double) A_NUCLEONS, 1.0 / 3.0); // GetNuclearRadius(16), b_max only
  double bMax = 1.3 * (R_nucRad + R_nucRad);
  double rEff = TMath::Sqrt(a_sigmaNN_mb / (10.0 * TMath::Pi()));
  double rEff2 = rEff * rEff;

  cout << "==================================================================" << endl;
  cout << "GenerateGlauberMC_OO200 -- Woods-Saxon (3-param Fermi) O-16 Glauber MC" << endl;
  cout << "  w=" << WS_W << "  c=" << WS_C << " fm  a=" << WS_A << " fm" << endl;
  cout << "  GetNuclearRadius(16) [b_max only] = " << R_nucRad << " fm" << endl;
  cout << "  b_max = " << bMax << " fm    r_eff (sigma_NN=" << a_sigmaNN_mb << " mb) = " << rEff << " fm" << endl;
  cout << "  Thrown events = " << a_nEvents << endl;
  cout << "==================================================================" << endl;

  //--------------------------------------------------------------------------
  // Build the r ~ rho(r)*r^2 sampling histogram (fine-binned inverse-CDF via
  // TH1::GetRandom()), and print the p0 normalization as a sanity check
  // (matches ComputeWoodsSaxonNormalizationParameter -- not needed for
  // sampling shape itself, GetRandom() only cares about relative bin
  // content, but useful for cross-checking against the Python run's printout).
  //--------------------------------------------------------------------------
  TH1D* hR = new TH1D("hR_WoodsSaxonPDF", "", R_GRID_NBINS, 0.0, R_GRID_MAX);
  double binWidth = R_GRID_MAX / R_GRID_NBINS;
  double volIntegralUnnorm = 0.0;
  for(int b = 1; b <= R_GRID_NBINS; b++){
    double r = hR->GetXaxis()->GetBinCenter(b);
    double w = WoodsSaxonShape(r) * r * r;
    hR->SetBinContent(b, w);
    volIntegralUnnorm += w * binWidth;
  }
  volIntegralUnnorm *= 4.0 * TMath::Pi();
  double p0Norm = (volIntegralUnnorm > 0) ? (A_NUCLEONS / volIntegralUnnorm) : 0.0;
  cout << "Computed p0 normalization (sanity check) = " << p0Norm << endl;

  //--------------------------------------------------------------------------
  // Event loop
  //--------------------------------------------------------------------------
  vector<double> bSuccess; bSuccess.reserve(a_nEvents / 2);
  vector<int> npartSuccess, ncollSuccess;
  npartSuccess.reserve(a_nEvents / 2);
  ncollSuccess.reserve(a_nEvents / 2);

  double xA[64], yA[64], xB[64], yB[64]; // sized generously beyond A_NUCLEONS=16
  bool collidedA[64], collidedB[64];

  Long64_t nSuccess = 0;
  TStopwatch sw; sw.Start();

  for(Long64_t ev = 0; ev < a_nEvents; ev++){
    double u = gRandom->Uniform(0.0, 1.0);
    double b = bMax * TMath::Sqrt(u); // inverse-CDF sample for PDF ~ b

    for(int i = 0; i < A_NUCLEONS; i++) SampleNucleonXY(hR, xA[i], yA[i]);
    for(int i = 0; i < A_NUCLEONS; i++){ SampleNucleonXY(hR, xB[i], yB[i]); xB[i] += b; }

    for(int i = 0; i < A_NUCLEONS; i++){ collidedA[i] = false; }
    for(int j = 0; j < A_NUCLEONS; j++){ collidedB[j] = false; }

    int nColl = 0;
    for(int i = 0; i < A_NUCLEONS; i++){
      for(int j = 0; j < A_NUCLEONS; j++){
        double dx = xA[i] - xB[j];
        double dy = yA[i] - yB[j];
        if(dx * dx + dy * dy < rEff2){
          nColl++;
          collidedA[i] = true;
          collidedB[j] = true;
        }
      }
    }

    int nPart = 0;
    for(int i = 0; i < A_NUCLEONS; i++) if(collidedA[i]) nPart++;
    for(int j = 0; j < A_NUCLEONS; j++) if(collidedB[j]) nPart++;

    if(nPart > 0){
      bSuccess.push_back(b);
      npartSuccess.push_back(nPart);
      ncollSuccess.push_back(nColl);
      nSuccess++;
    }

    if(ev > 0 && ev % 500000 == 0){
      cout << "  ... " << ev << " / " << a_nEvents << " thrown (" << nSuccess << " successful so far), "
           << sw.RealTime() << " s elapsed" << endl;
      sw.Continue();
    }
  }
  sw.Stop();
  cout << "Event loop done in " << sw.RealTime() << " s (" << sw.CpuTime() << " s CPU)." << endl;

  double successRate = (double) nSuccess / (double) a_nEvents;
  double sigmaReactionFm2 = successRate * TMath::Pi() * bMax * bMax;
  double sigmaReactionMb = sigmaReactionFm2 * 10.0;

  cout << "Successful (reacting) events = " << nSuccess << " (" << 100.0 * successRate << "%)" << endl;
  cout << "sigma_reaction (geometric)   = " << sigmaReactionMb << " mb" << endl;

  //--------------------------------------------------------------------------
  // Rank successful events by impact parameter (ascending = most central
  // first) and bin at the standard centrality percentile edges.
  //--------------------------------------------------------------------------
  vector<size_t> order(nSuccess);
  iota(order.begin(), order.end(), 0);
  sort(order.begin(), order.end(), [&](size_t a, size_t b_idx){ return bSuccess[a] < bSuccess[b_idx]; });

  cout << endl;
  printf("%-10s%12s%10s%10s%10s%10s%12s%12s\n", "Cent bin", "N_events", "<b> fm", "<Npart>", "<Nspec>", "<Ncoll>", "Npart RMS", "Nspec RMS");

  struct BinResult {
    string cent; long long nEvents; double meanB, meanNpart, meanNspec, meanNcoll, rmsNpart, rmsNspec;
    double fracNpart, fracNspec;
  };
  vector<BinResult> results;

  for(int i = 0; i < N_CENT_BINS; i++){
    long long lo = (long long) llround(CENT_EDGES[i] * nSuccess);
    long long hi = (long long) llround(CENT_EDGES[i + 1] * nSuccess);
    long long n = hi - lo;
    if(n <= 0){
      cerr << "WARNING: centrality bin " << CENT_LABELS[i] << " is empty -- check N_EVENTS/statistics." << endl;
      continue;
    }

    double sumB = 0, sumNpart = 0, sumNspec = 0, sumNcoll = 0;
    vector<double> npartVals, nspecVals;
    npartVals.reserve(n); nspecVals.reserve(n);
    for(long long k = lo; k < hi; k++){
      size_t idx = order[k];
      double npart = npartSuccess[idx];
      double nspec = 2.0 * A_NUCLEONS - npart;
      sumB += bSuccess[idx];
      sumNpart += npart;
      sumNspec += nspec;
      sumNcoll += ncollSuccess[idx];
      npartVals.push_back(npart);
      nspecVals.push_back(nspec);
    }
    double meanB = sumB / n, meanNpart = sumNpart / n, meanNspec = sumNspec / n, meanNcoll = sumNcoll / n;

    double sqDevNpart = 0, sqDevNspec = 0;
    for(long long k = 0; k < n; k++){
      sqDevNpart += (npartVals[k] - meanNpart) * (npartVals[k] - meanNpart);
      sqDevNspec += (nspecVals[k] - meanNspec) * (nspecVals[k] - meanNspec);
    }
    double rmsNpart = (n > 1) ? TMath::Sqrt(sqDevNpart / (n - 1)) : 0.0;
    double rmsNspec = (n > 1) ? TMath::Sqrt(sqDevNspec / (n - 1)) : 0.0;

    printf("%-10s%12lld%10.3f%10.2f%10.2f%10.2f%12.2f%12.2f\n",
           CENT_LABELS[i], n, meanB, meanNpart, meanNspec, meanNcoll, rmsNpart, rmsNspec);

    BinResult r;
    r.cent = CENT_LABELS[i]; r.nEvents = n; r.meanB = meanB;
    r.meanNpart = meanNpart; r.meanNspec = meanNspec; r.meanNcoll = meanNcoll;
    r.rmsNpart = rmsNpart; r.rmsNspec = rmsNspec;
    r.fracNpart = meanNpart / (2.0 * A_NUCLEONS);
    r.fracNspec = meanNspec / (2.0 * A_NUCLEONS);
    results.push_back(r);
  }

  cout << endl << "Total nucleons in system (2*A) = " << 2 * A_NUCLEONS << endl << endl;
  printf("%-10s%12s%12s\n", "Cent bin", "Npart frac", "Nspec frac");
  for(size_t i = 0; i < results.size(); i++){
    printf("%-10s%12.3f%12.3f\n", results[i].cent.c_str(), results[i].fracNpart, results[i].fracNspec);
  }

  //--------------------------------------------------------------------------
  // Write glauber_oo200_results_woodssaxon.root -- ONE ROOT file, three
  // TTrees, replacing the former JSON config+centrality_bins and the CSV
  // cell table (see header for the schema and why -- Python/JSON are
  // retired from this pipeline per Andrew's instruction).
  //--------------------------------------------------------------------------
  string rootPath = a_outputDir + "/glauber_oo200_results_woodssaxon.root";
  TFile* outFile = new TFile(rootPath.c_str(), "RECREATE");
  outFile->cd();

  // ---- GlauberCentralityBins: one entry per centrality bin ----
  TTree* tBins = new TTree("GlauberCentralityBins",
                            "Woods-Saxon Glauber MC, b-percentile-ranked <Npart>/<Nspec>/<Ncoll> per bin");
  Char_t   cent_out[16];
  Double_t nEvents_out, meanB_out, meanNpart_out, meanNspec_out, meanNcoll_out;
  Double_t npartRMS_out, nspecRMS_out, fracNpart_out, fracNspec_out, refMultEdgeHi_out;
  tBins->Branch("cent",             cent_out,           "cent/C");
  tBins->Branch("nEvents",          &nEvents_out,       "nEvents/D");
  tBins->Branch("meanB_fm",         &meanB_out,         "meanB_fm/D");
  tBins->Branch("meanNpart",        &meanNpart_out,     "meanNpart/D");
  tBins->Branch("meanNspec",        &meanNspec_out,     "meanNspec/D");
  tBins->Branch("meanNcoll",        &meanNcoll_out,     "meanNcoll/D");
  tBins->Branch("NpartRMS",         &npartRMS_out,      "NpartRMS/D");
  tBins->Branch("NspecRMS",         &nspecRMS_out,      "NspecRMS/D");
  tBins->Branch("fracNpartOfTotal", &fracNpart_out,     "fracNpartOfTotal/D");
  tBins->Branch("fracNspecOfTotal", &fracNspec_out,     "fracNspecOfTotal/D");
  tBins->Branch("refMultEdgeHi",    &refMultEdgeHi_out, "refMultEdgeHi/D");

  for(size_t i = 0; i < results.size(); i++){
    const BinResult& r = results[i];
    strncpy(cent_out, r.cent.c_str(), sizeof(cent_out) - 1);
    cent_out[sizeof(cent_out) - 1] = '\0';
    nEvents_out    = (double) r.nEvents;
    meanB_out      = r.meanB;
    meanNpart_out  = r.meanNpart;
    meanNspec_out  = r.meanNspec;
    meanNcoll_out  = r.meanNcoll;
    npartRMS_out   = r.rmsNpart;
    nspecRMS_out   = r.rmsNspec;
    fracNpart_out  = r.fracNpart;
    fracNspec_out  = r.fracNspec;
    refMultEdgeHi_out = -999.0; // not applicable -- purely geometric, no refMult here
    tBins->Fill();
  }
  tBins->Write();

  // ---- GlauberCells: every distinct (Npart,Ncoll) pair among ALL
  // successful events (not just one centrality bin) with its event count --
  // replaces GlauberNpartNcollWeights_OO200.csv for
  // FitGlauberNBDToRefMult_OO200.C's NBD fit. See header for why this
  // collapsed form is exactly what that fit needs. ----
  map<pair<int,int>, long long> cellCounts;
  for(Long64_t i = 0; i < nSuccess; i++){
    cellCounts[make_pair(npartSuccess[i], ncollSuccess[i])]++;
  }
  TTree* tCells = new TTree("GlauberCells",
                             "Distinct (Npart,Ncoll) cells among successful events, with event-count weight");
  Int_t   npart_cell_out, nspec_cell_out, ncoll_cell_out;
  Long64_t weight_cell_out;
  tCells->Branch("Npart",  &npart_cell_out, "Npart/I");
  tCells->Branch("Nspec",  &nspec_cell_out, "Nspec/I");
  tCells->Branch("Ncoll",  &ncoll_cell_out, "Ncoll/I");
  tCells->Branch("weight", &weight_cell_out, "weight/L");

  long long totalWeightCheck = 0;
  for(map<pair<int,int>, long long>::const_iterator it = cellCounts.begin(); it != cellCounts.end(); ++it){
    npart_cell_out = it->first.first;
    ncoll_cell_out = it->first.second;
    nspec_cell_out = 2 * A_NUCLEONS - npart_cell_out;
    weight_cell_out = it->second;
    totalWeightCheck += weight_cell_out;
    tCells->Fill();
  }
  tCells->Write();
  cout << endl << cellCounts.size() << " distinct (Npart,Ncoll) cells, " << totalWeightCheck
       << " successful events total -- should equal N_success=" << nSuccess << " above." << endl;

  // ---- GlauberConfig: one entry, run configuration/diagnostics. Not read
  // back by any downstream macro -- archival, same role the JSON "config"
  // sub-object played. ----
  TTree* tConfig = new TTree("GlauberConfig", "Run configuration/diagnostics for this Glauber MC");
  Char_t   model_out[64];
  Int_t    aNucleons_out;
  Double_t sigmaNN_out;
  Long64_t nThrown_out, nSuccess_out;
  Double_t wsW_out, wsC_out, wsA_out, rNucRad_out, bMax_out, rEff_out, sigmaReaction_out;
  strncpy(model_out, "WoodsSaxon (3-param Fermi, O-16)", sizeof(model_out) - 1);
  model_out[sizeof(model_out) - 1] = '\0';
  tConfig->Branch("model",          model_out,           "model/C");
  tConfig->Branch("A_nucleons",     &aNucleons_out,      "A_nucleons/I");
  tConfig->Branch("sigma_NN_mb",    &sigmaNN_out,        "sigma_NN_mb/D");
  tConfig->Branch("N_thrown",       &nThrown_out,        "N_thrown/L");
  tConfig->Branch("N_success",      &nSuccess_out,       "N_success/L");
  tConfig->Branch("WS_w",           &wsW_out,            "WS_w/D");
  tConfig->Branch("WS_c_fm",        &wsC_out,            "WS_c_fm/D");
  tConfig->Branch("WS_a_fm",        &wsA_out,            "WS_a_fm/D");
  tConfig->Branch("R_nuclearRadius_fm_bmax_only", &rNucRad_out, "R_nuclearRadius_fm_bmax_only/D");
  tConfig->Branch("b_max_fm",       &bMax_out,           "b_max_fm/D");
  tConfig->Branch("r_eff_fm",       &rEff_out,           "r_eff_fm/D");
  tConfig->Branch("sigma_reaction_mb", &sigmaReaction_out, "sigma_reaction_mb/D");

  aNucleons_out = A_NUCLEONS; sigmaNN_out = a_sigmaNN_mb;
  nThrown_out = a_nEvents; nSuccess_out = nSuccess;
  wsW_out = WS_W; wsC_out = WS_C; wsA_out = WS_A;
  rNucRad_out = R_nucRad; bMax_out = bMax; rEff_out = rEff; sigmaReaction_out = sigmaReactionMb;
  tConfig->Fill();
  tConfig->Write();

  outFile->Close();
  cout << "Wrote " << rootPath << " (TTrees \"GlauberCentralityBins\" [" << results.size()
       << " entries], \"GlauberCells\" [" << cellCounts.size() << " entries], \"GlauberConfig\" [1 entry])" << endl;

  cout << endl << "Reminder: this is a purely geometric Glauber MC (no refMult/detector"
       << " multiplicity model folded in yet) -- see FitGlauberNBDToRefMult_OO200.C for that step."
       << endl;
}
