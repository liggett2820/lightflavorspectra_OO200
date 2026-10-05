// FitGlauberNBDToRefMult_OO200.C -- lightflavorspectra_OO200
//
// PURPOSE: the real NBD-Glauber centrality-determination fit Andrew asked for
// ("I need it based in a glauber model of my ref mult plot from one of my
// yield histos files") -- the standard STAR/PHOBOS/ALICE technique, distinct
// from and more work than the two things already delivered this session:
//   - dose/EstimateFireballFractionFromRefMult_OO200.C combines a REAL refMult
//     histogram with Npart/Nspec numbers from a GENERIC Glauber MC run (b
//     ranked, percentile-cut, not tied to this detector's actual measured
//     multiplicity shape at all).
//   - dose/EstimateNpartNspecFractions_OO200.py's Npart/Nspec-by-centrality
//     numbers are likewise purely geometric: rank successful Glauber MC events
//     by IMPACT PARAMETER and cut at percentile boundaries -- no refMult
//     involved.
// THIS macro instead folds a two-component negative-binomial (NBD) ancestor
// multiplicity model through the Glauber MC and FITS its free parameters so
// the predicted multiplicity distribution's SHAPE matches the real measured
// refMult histogram. Centrality edges and <Npart>/<Ncoll>/<Nspec> per bin are
// then derived from that fit, not from geometry alone -- this is what makes it
// "based on a Glauber model of the refMult plot" rather than just "next to" it.
//
// PHYSICS MODEL (two-component ancestor/NBD, the standard method -- see e.g.
// the approach used by TGlauberMC's own NBD-fit utilities and standard STAR/
// PHOBOS centrality-determination papers):
//   For a Glauber MC event with Npart participants and Ncoll binary collisions,
//   the mean of its predicted multiplicity is
//     nbar(Npart,Ncoll) = npp * [ (1-x)*Npart/2 + x*Ncoll ]
//   i.e. a mix of "soft" (proportional to participant pairs, Npart/2) and
//   "hard" (proportional to binary collisions, Ncoll) particle production,
//   weighted by x in [0,1]. The predicted per-event multiplicity is then drawn
//   from a Negative Binomial Distribution with that mean and a SHARED shape
//   parameter k (same k for every event -- k controls the extra fluctuation
//   beyond Poisson, a single "one ancestor's" multiplicity-fluctuation width).
//   npp, k, x are the free parameters fit here. This is a SHAPE fit: an
//   overall scale factor (folding in refMult's limited acceptance/efficiency,
//   which npp alone cannot fully absorb since it also sets the shape via the
//   Npart/Ncoll mix) is profiled out analytically each iteration rather than
//   fit as a 4th Minuit parameter -- see ProfiledScale() below.
//
// WHY A WEIGHTED-CELL SUM INSTEAD OF RE-THROWING MC EVENTS PER FIT ITERATION:
// Npart and Ncoll are both small non-negative integers, so tabulating
// (Npart, Ncoll, weight=event count) over the Glauber MC's full successful-
// event sample collapses millions of raw events into a few thousand distinct
// cells (done in EstimateNpartNspecFractions_OO200.py's newly-added CSV
// export -- see that script's added block and GlauberNpartNcollWeights_OO200.csv).
// Summing weight*NBD_PMF(m; nbar(Npart,Ncoll), k) over these cells is exactly
// mathematically equivalent to summing over every individual MC event (each
// event's predicted multiplicity distribution IS an NBD with that event's own
// nbar), but touches thousands of terms instead of millions per chi2
// evaluation -- essential for this to run in a Minuit loop in finite time.
//
// FIT RANGE -- WHY NOT THE WHOLE refMult HISTOGRAM: real refMult distributions
// are truncated/distorted at low multiplicity by trigger and vertex-finding
// inefficiency (an event with very few tracks is less likely to satisfy the
// trigger or reconstruct a good vertex at all) -- fitting that region would
// bias npp/k/x toward compensating for detector inefficiency rather than
// genuine multiplicity-fluctuation physics. This macro therefore fits only
// refMult in [a_refMultFitLo, a_refMultFitHi) by default, skipping the low end
// entirely; a_refMultFitLo defaults to 15 (adjust based on where YOUR
// histogram visibly turns over/flattens -- not re-derived from source here).
// a_refMultFitHi defaults to -1 = auto-detect: the highest refMult bin with at
// least a_minBinCountForFit (default 20) events, to avoid fitting into the
// sparse/noisy statistical tail. PRINT AND LOOK AT YOUR HISTOGRAM before
// trusting these defaults for a real O+O-sized (small) system -- there is no
// prior O+O refMult distribution in this repo to calibrate them against.
//
// ONCE FIT: centrality edges are set the way it's actually meant to be done --
// NOT by reading percentiles off the raw (trigger-truncated) DATA histogram,
// but off the FITTED model extrapolated down to refMult=0 (the NBD_PMF
// naturally extends there), which estimates the TRUE total reaction cross
// section including the low-multiplicity events the trigger/vertex cuts lose.
// Cutting percentiles on the raw truncated data would systematically misplace
// the low-refMult edges (see header note above on trigger inefficiency) -- see
// DerivedCentralityBin() and the main macro body below for exactly how the
// extrapolated distribution is used for both the refMult edges AND the
// Npart/Ncoll/Nspec means per bin.
//
// INPUTS -- CHANGED from CSV+JSON to a single ROOT file, per Andrew's
// explicit instruction to get Python and JSON out of this pipeline entirely
// (see GenerateGlauberMC_OO200.C's header for the full schema and the retired
// EstimateNpartNspecFractions_OO200.py/PlotDoseFractionComparison_OO200.py):
//   a_glauberRoot -- glauber_oo200_results_woodssaxon.root, from
//                    GenerateGlauberMC_OO200.C. TWO trees read from this ONE
//                    file (previously two separate files, a CSV and a JSON):
//                      "GlauberCells" -- Npart/I, Nspec/I, Ncoll/I, weight/L,
//                        replaces GlauberNpartNcollWeights_OO200.csv for the
//                        NBD fit itself.
//                      "GlauberCentralityBins" -- replaces the old-vs-new
//                        comparison table's JSON source.
//   a_yieldFile  -- a PicoBinner output ROOT file with a top-level "refMult"
//                   TH1I histogram (same source/convention as
//                   EstimateFireballFractionFromRefMult_OO200.C -- see that
//                   macro's header for exactly how/where refMult is booked).
//                   Andrew confirmed: /Users/aliggett/data/OO/yieldHistos_OO200_pion.root
//
// OUTPUTS (under a_outputDir):
//   OO200_NBDGlauberFit_RefMult.png   -- data vs. best-fit prediction overlay
//     (log-y), with the derived centrality-bin edges marked as vertical lines.
//   glauber_oo200_results_NBDfit.root -- CHANGED from JSON to ROOT, same
//     reasoning as the input side. TTree "GlauberCentralityBins" -- SAME
//     schema as GenerateGlauberMC_OO200.C's own GlauberCentralityBins tree
//     (cent, meanNpart, meanNspec, meanNcoll, fracNpartOfTotal,
//     fracNspecOfTotal, refMultEdgeHi, nEvents; meanB_fm/NpartRMS/NspecRMS
//     filled as -999 sentinels here since this fit doesn't compute them) --
//     any downstream macro that reads a GlauberCentralityBins tree can point
//     at either this file or the geometric one interchangeably, same drop-in
//     intent the JSON schema used to serve. Plus TTree "GlauberConfig" (one
//     entry: npp/k/x + Minuit errors, chi2, ndof, refMult fit range) --
//     informational only, not read back downstream.
//   Console: fit parameters + Minuit errors, chi2/ndf, the derived refMult
//     edges next to SetCutClass.C's currently-assumed {44,37,28,17,5,0} as a
//     direct sanity check, and old-geometric-vs-new-fit-derived <Npart>/<Nspec>
//     per bin (read back from a_glauberRoot's GlauberCentralityBins tree if
//     present -- soft-skipped if that file isn't found, same graceful-
//     degradation policy as PlotDoseFractionComparison_OO200.C's optional
//     Glauber overlay).
//
// USAGE (ROOT, from dose/ -- run GenerateGlauberMC_OO200.C FIRST if
// glauber_oo200_results_woodssaxon.root doesn't exist yet):
//   root -l -b -q 'GenerateGlauberMC_OO200.C+()'
//   root -l -b -q 'FitGlauberNBDToRefMult_OO200.C("glauber_oo200_results_woodssaxon.root","/Users/aliggett/data/OO/yieldHistos_OO200_pion.root")'
// PERFORMANCE: compile it (note the trailing '+') if Migrad is slow --
//   root -l -b -q 'FitGlauberNBDToRefMult_OO200.C+("glauber_oo200_results_woodssaxon.root","/Users/aliggett/data/OO/yieldHistos_OO200_pion.root")'
// Each chi2 evaluation is O(nCells * nBinsFit); a_minCellWeight (default 5)
// drops low-weight cells to bound this -- the printed "retained weight
// fraction" should stay above ~99% or the threshold is cutting real signal.
//
// STATUS: UNTESTED, same caveat as every ROOT macro in this session -- no ROOT
// in the sandbox that wrote this. This one carries MORE fit-stability risk
// than anything else delivered so far (Minuit convergence depends on
// reasonable starting values/limits and on a_refMultFitLo/Hi actually avoiding
// the trigger-biased region of YOUR real histogram -- neither was checked
// against real data here). Run it, look at the printed chi2/ndf and whether
// any fitted parameter sits on its limit (printed explicitly), and look at the
// overlay plot before trusting the derived Npart/Nspec numbers for anything
// downstream. This is a genuine first-pass calibration, not a final centrality
// definition -- exactly the caveat this repo already carries for the purely-
// geometric Glauber numbers it replaces. Andrew HAS run an earlier
// (CSV+JSON-based) version of this fit against real data -- the real
// glauber_oo200_results_NBDfit.json on disk reflects that run's actual fitted
// npp/k/x and per-bin Npart/Nspec. The I/O was just switched to ROOT and has
// NOT been re-exercised -- re-run and confirm the fit converges to the same
// npp/k/x (same chi2/ndf) before trusting glauber_oo200_results_NBDfit.root.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;

//============================================================================
// Palette -- same hex values used across this repo's dose/ figures.
//============================================================================
const int COLOR_DATA = kBlack;
const int COLOR_FIT  = TColor::GetColor("#C0392B"); // brick red, matches "fireball" elsewhere
const int COLOR_EDGE = TColor::GetColor("#1E2761"); // navy, matches "participant" elsewhere

//============================================================================
// A weighted Glauber MC cell: all successful events with this exact
// (Npart, Ncoll) pair, collapsed to one entry with a count.
//============================================================================
struct GlauberCell {
  int    Npart, Nspec, Ncoll;
  double weight;
  double nbar; // filled in per fit-parameter evaluation, see BuildPrediction()
};

//============================================================================
// Globals used by the Minuit chi2 functor (ROOT::Math::Functor needs a plain
// callable; keeping the fit inputs at file scope is the simplest correct way
// to do this in a single macro, matching the level of complexity already used
// elsewhere in this repo for hand-rolled readers/fitters).
//============================================================================
vector<GlauberCell> g_cells;
int    g_fitLo = 0, g_fitHi = 0; // integer refMult range used in the chi2 sum
vector<double> g_dataContent;    // data histogram content, bins [g_fitLo, g_fitHi)
vector<double> g_dataError;      // sqrt(content), floor 1 where content>0

//============================================================================
// NBD_PMF -- P(m; nbar, k), evaluated via log-gamma to avoid overflow for
// large k/m. m is treated as a real number via the Gamma-function
// generalization of the binomial coefficient (standard trick; m is always an
// integer refMult value in practice here).
//============================================================================
double NBD_PMF(double m, double nbar, double k){
  if(nbar <= 0.0 || k <= 0.0) return 0.0;
  double logP = TMath::LnGamma(m + k) - TMath::LnGamma(k) - TMath::LnGamma(m + 1.0)
              + k * TMath::Log(k / (k + nbar))
              + m * TMath::Log(nbar / (k + nbar));
  return TMath::Exp(logP);
}

//============================================================================
// BuildPrediction -- fills nbar for every cell (given npp,k,x) and returns
// the (unnormalized) predicted content in integer bins [mLo, mHi).
//============================================================================
vector<double> BuildPrediction(vector<GlauberCell>& cells, double npp, double k, double x, int mLo, int mHi){
  vector<double> pred(mHi - mLo, 0.0);
  for(size_t c = 0; c < cells.size(); c++){
    double nbar = npp * ((1.0 - x) * cells[c].Npart / 2.0 + x * cells[c].Ncoll);
    cells[c].nbar = nbar;
    if(nbar <= 0.0) continue;
    for(int m = mLo; m < mHi; m++){
      pred[m - mLo] += cells[c].weight * NBD_PMF((double) m, nbar, k);
    }
  }
  return pred;
}

//============================================================================
// ProfiledScale -- closed-form linear scale factor minimizing
// sum[(data_i - scale*pred_i)^2 / err_i^2] over the fit-range bins, so the
// overall normalization (which folds in refMult's acceptance/efficiency,
// separate from the shape parameters) doesn't need its own Minuit parameter.
//============================================================================
double ProfiledScale(const vector<double>& pred){
  double num = 0.0, den = 0.0;
  for(size_t i = 0; i < pred.size(); i++){
    if(g_dataContent[i] <= 0.0) continue;
    double e2 = g_dataError[i] * g_dataError[i];
    num += g_dataContent[i] * pred[i] / e2;
    den += pred[i] * pred[i] / e2;
  }
  return (den > 0.0) ? (num / den) : 0.0;
}

//============================================================================
// Chi2Function -- the Minuit objective. par = {npp, k, x}.
//============================================================================
double Chi2Function(const double* par){
  double npp = par[0], k = par[1], x = par[2];
  vector<double> pred = BuildPrediction(g_cells, npp, k, x, g_fitLo, g_fitHi);
  double scale = ProfiledScale(pred);
  double chi2 = 0.0;
  for(size_t i = 0; i < pred.size(); i++){
    if(g_dataContent[i] <= 0.0) continue;
    double e2 = g_dataError[i] * g_dataError[i];
    double d = g_dataContent[i] - scale * pred[i];
    chi2 += d * d / e2;
  }
  return chi2;
}

//============================================================================
// ReadGlauberCells -- reads the "GlauberCells" TTree (Npart/I, Nspec/I,
// Ncoll/I, weight/L), written by GenerateGlauberMC_OO200.C. Replaces the
// former hand-rolled CSV reader (GlauberNpartNcollWeights_OO200.csv, from the
// now-retired EstimateNpartNspecFractions_OO200.py). Returns an EMPTY vector
// (not fabricated) on failure -- caller must check and abort.
//============================================================================
vector<GlauberCell> ReadGlauberCells(const string& a_path, double a_minCellWeight, double& a_totalWeightOut, double& a_keptWeightOut){
  vector<GlauberCell> cells;
  a_totalWeightOut = 0.0;
  a_keptWeightOut = 0.0;

  TFile* f = TFile::Open(a_path.c_str(), "READ");
  if(!f || f->IsZombie()){
    cerr << "ERROR: could not open '" << a_path << "'." << endl;
    return cells;
  }
  TTree* t = (TTree*) f->Get("GlauberCells");
  if(!t){
    cerr << "ERROR: no \"GlauberCells\" TTree in '" << a_path << "'"
         << " (expected from GenerateGlauberMC_OO200.C -- run that macro first)." << endl;
    f->Close();
    return cells;
  }

  Int_t npart_in, nspec_in, ncoll_in;
  Long64_t weight_in;
  t->SetBranchAddress("Npart", &npart_in);
  t->SetBranchAddress("Nspec", &nspec_in);
  t->SetBranchAddress("Ncoll", &ncoll_in);
  t->SetBranchAddress("weight", &weight_in);

  Long64_t nEntries = t->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    double w = (double) weight_in;
    a_totalWeightOut += w;
    if(w < a_minCellWeight) continue; // drop rare low-weight cells, see header perf note
    a_keptWeightOut += w;
    GlauberCell cell;
    cell.Npart  = npart_in;
    cell.Nspec  = nspec_in;
    cell.Ncoll  = ncoll_in;
    cell.weight = w;
    cell.nbar   = 0.0;
    cells.push_back(cell);
  }
  f->Close();
  return cells;
}

//============================================================================
// Read mean_Npart/mean_Nspec per "cent" label out of the existing geometric
// Glauber ROOT file's "GlauberCentralityBins" tree, for the console
// cross-check table. Replaces the former hand-rolled JSON line-scanning
// reader. Soft-fails (empty vector) if the file/tree is missing -- this is
// the one input in this macro allowed to be absent (see header: the old-vs-
// new table is a secondary comparison, not the fit itself).
//============================================================================
struct OldGlauberRow { double meanNpart = 0, meanNspec = 0; bool found = false; };
vector<pair<string, OldGlauberRow>> ReadOldGeometricGlauber(const string& a_path){
  vector<pair<string, OldGlauberRow>> out;
  TFile* f = TFile::Open(a_path.c_str(), "READ");
  if(!f || f->IsZombie()){
    cout << "NOTE: '" << a_path << "' not found -- skipping old-vs-new comparison table." << endl;
    return out;
  }
  TTree* t = (TTree*) f->Get("GlauberCentralityBins");
  if(!t){
    cout << "NOTE: no \"GlauberCentralityBins\" tree in '" << a_path
         << "' -- skipping old-vs-new comparison table." << endl;
    f->Close();
    return out;
  }

  Char_t cent_in[16];
  Double_t meanNpart_in, meanNspec_in;
  t->SetBranchAddress("cent", cent_in);
  t->SetBranchAddress("meanNpart", &meanNpart_in);
  t->SetBranchAddress("meanNspec", &meanNspec_in);

  Long64_t nEntries = t->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    OldGlauberRow row;
    row.meanNpart = meanNpart_in;
    row.meanNspec = meanNspec_in;
    row.found = true;
    out.push_back(make_pair(string(cent_in), row));
  }
  f->Close();
  return out;
}

//============================================================================
// Main macro
//============================================================================
void FitGlauberNBDToRefMult_OO200(
    string a_glauberRoot       = "glauber_oo200_results_woodssaxon.root", // GlauberCells + GlauberCentralityBins, both read from this ONE file -- see header
    string a_yieldFile         = "/Users/aliggett/data/OO/yieldHistos_OO200_pion.root",
    string a_outputDir         = ".",
    int    a_refMultFitLo      = 15,   // skip the trigger/vertex-inefficiency region -- see header
    int    a_refMultFitHi      = -1,   // -1 = auto-detect from the data histogram, see header
    int    a_minBinCountForFit = 20,   // used only for the auto-detect above
    double a_minCellWeight     = 5.0   // drop rarer (Npart,Ncoll) cells, see header perf note
){
  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outputDir.c_str(), kTRUE);

  //--------------------------------------------------------------------------
  // Load Glauber MC cells
  //--------------------------------------------------------------------------
  double totalWeight = 0.0, keptWeight = 0.0;
  g_cells = ReadGlauberCells(a_glauberRoot, a_minCellWeight, totalWeight, keptWeight);
  if(g_cells.empty()){
    cerr << "ERROR: no usable Glauber cells read from '" << a_glauberRoot << "' -- aborting."
         << " (Run GenerateGlauberMC_OO200.C first -- it writes this file's \"GlauberCells\" tree.)" << endl;
    return;
  }
  double keptFrac = (totalWeight > 0) ? (100.0 * keptWeight / totalWeight) : 0.0;
  cout << "Read " << g_cells.size() << " Glauber cells from '" << a_glauberRoot << "' ("
       << keptFrac << "% of total MC weight retained after the weight>=" << a_minCellWeight
       << " cut -- should be close to 100%)." << endl;
  if(keptFrac < 99.0){
    cout << "WARNING: retained weight fraction is below 99% -- a_minCellWeight may be dropping"
         << " real signal, not just noise. Consider lowering it." << endl;
  }

  //--------------------------------------------------------------------------
  // Load real refMult histogram
  //--------------------------------------------------------------------------
  TFile* inFile = TFile::Open(a_yieldFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){
    cerr << "ERROR: could not open '" << a_yieldFile << "'." << endl;
    return;
  }
  TH1* hRefMult = (TH1*) inFile->Get("refMult");
  if(!hRefMult){
    cerr << "ERROR: no top-level 'refMult' histogram in '" << a_yieldFile << "'"
         << " (see EstimateFireballFractionFromRefMult_OO200.C's header for how it's booked)." << endl;
    inFile->Close();
    return;
  }
  hRefMult->SetDirectory(0);
  inFile->Close();

  //--------------------------------------------------------------------------
  // Auto-detect the upper fit edge if requested
  //--------------------------------------------------------------------------
  if(a_refMultFitHi < 0){
    int hiBin = 1;
    for(int b = 1; b <= hRefMult->GetNbinsX(); b++){
      if(hRefMult->GetBinContent(b) >= a_minBinCountForFit) hiBin = b;
    }
    a_refMultFitHi = (int) hRefMult->GetXaxis()->GetBinUpEdge(hiBin);
    cout << "Auto-detected refMult fit upper edge = " << a_refMultFitHi
         << " (highest bin with >= " << a_minBinCountForFit << " events)." << endl;
  }
  if(a_refMultFitHi <= a_refMultFitLo){
    cerr << "ERROR: refMult fit range [" << a_refMultFitLo << "," << a_refMultFitHi
         << ") is empty or inverted -- check a_refMultFitLo/a_refMultFitHi against your histogram." << endl;
    return;
  }
  cout << "Fitting refMult in [" << a_refMultFitLo << ", " << a_refMultFitHi << ")." << endl;

  g_fitLo = a_refMultFitLo;
  g_fitHi = a_refMultFitHi;
  int nBinsFit = g_fitHi - g_fitLo;
  g_dataContent.assign(nBinsFit, 0.0);
  g_dataError.assign(nBinsFit, 0.0);
  for(int m = g_fitLo; m < g_fitHi; m++){
    int b = hRefMult->GetXaxis()->FindBin((double) m);
    double c = hRefMult->GetBinContent(b);
    g_dataContent[m - g_fitLo] = c;
    g_dataError[m - g_fitLo] = (c > 0.0) ? TMath::Sqrt(c) : 1.0;
  }

  //--------------------------------------------------------------------------
  // Data-driven starting guess for npp: match the data's mean multiplicity in
  // the top decile of the fit range to the Glauber cells' mean (Npart/2,Ncoll)
  // in that same rough region, so Migrad starts in a physically sensible spot
  // instead of an arbitrary hardcoded constant.
  //--------------------------------------------------------------------------
  double dataHighEnd = g_fitLo + 0.9 * (g_fitHi - g_fitLo);
  double sumM = 0.0, sumW = 0.0;
  for(int m = (int) dataHighEnd; m < g_fitHi; m++){
    double c = g_dataContent[m - g_fitLo];
    sumM += c * m; sumW += c;
  }
  double dataMeanHighEnd = (sumW > 0) ? (sumM / sumW) : (double) g_fitHi;

  double maxNpart = 0, maxNcoll = 0;
  for(size_t c = 0; c < g_cells.size(); c++){
    maxNpart = max(maxNpart, (double) g_cells[c].Npart);
    maxNcoll = max(maxNcoll, (double) g_cells[c].Ncoll);
  }
  double softHardAtMax = 0.5 * (0.5 * maxNpart) + 0.5 * maxNcoll; // x=0.5 seed mix
  double nppSeed = (softHardAtMax > 0) ? (dataMeanHighEnd / softHardAtMax) : 1.0;
  cout << "Data-driven npp starting guess = " << nppSeed
       << " (from mean refMult=" << dataMeanHighEnd << " in the top ~10% of the fit range"
       << " vs. max-Npart/Ncoll Glauber cell)." << endl;

  //--------------------------------------------------------------------------
  // Minuit fit: par = {npp, k, x}
  //--------------------------------------------------------------------------
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  if(!min){
    cout << "NOTE: Minuit2 not available, falling back to legacy Minuit." << endl;
    min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
  }
  if(!min){
    cerr << "ERROR: could not create any Minuit minimizer -- is ROOT's minimizer plugin available?" << endl;
    return;
  }
  min->SetMaxFunctionCalls(100000);
  min->SetMaxIterations(10000);
  min->SetTolerance(0.01);
  min->SetPrintLevel(1);

  ROOT::Math::Functor f(&Chi2Function, 3);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "npp", nppSeed, nppSeed * 0.05, 1e-6, 1000.0);
  min->SetLimitedVariable(1, "k",   1.5,     0.1,            1e-3, 200.0);
  min->SetLimitedVariable(2, "x",   0.1,     0.02,           0.0,  1.0);

  bool fitOk = min->Minimize();
  if(!fitOk){
    cerr << "WARNING: Minuit did not report full convergence -- inspect the printout above"
         << " and treat the result below as provisional." << endl;
  }

  const double* par = min->X();
  const double* err = min->Errors();
  double nppFit = par[0], kFit = par[1], xFit = par[2];
  double chi2min = Chi2Function(par);
  int ndof = nBinsFit - 3 - 1; // 3 shape params + 1 profiled scale
  double chi2ndf = (ndof > 0) ? (chi2min / ndof) : -1.0;

  cout << "==================================================================" << endl;
  cout << "NBD-Glauber fit result:" << endl;
  cout << "  npp = " << nppFit << " +/- " << err[0] << endl;
  cout << "  k   = " << kFit   << " +/- " << err[1] << endl;
  cout << "  x   = " << xFit   << " +/- " << err[2] << endl;
  cout << "  chi2/ndf = " << chi2min << " / " << ndof << " = " << chi2ndf << endl;
  if(xFit < 1e-3 || xFit > 0.999) cout << "  WARNING: x is at/near its [0,1] limit." << endl;
  if(kFit > 199.0) cout << "  WARNING: k is at/near its upper limit -- consider raising it." << endl;
  cout << "==================================================================" << endl;

  double scaleAtBest = ProfiledScale(BuildPrediction(g_cells, nppFit, kFit, xFit, g_fitLo, g_fitHi));

  //--------------------------------------------------------------------------
  // Derive centrality edges + <Npart>/<Ncoll>/<Nspec> per bin from the fitted
  // model EXTRAPOLATED down to refMult=0 -- see header for why the raw
  // (trigger-truncated) data histogram is NOT used for this percentile cut.
  //--------------------------------------------------------------------------
  const int N_CENT_BINS = 6;
  const char* CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
  const double CENT_EDGES_FRAC[N_CENT_BINS + 1] = {0.0, 0.05, 0.10, 0.20, 0.40, 0.80, 1.00};

  // Sort cells by descending nbar (most-central first), matching the
  // b-percentile convention used elsewhere in this repo, but ranked by the
  // FITTED multiplicity proxy instead of impact parameter.
  BuildPrediction(g_cells, nppFit, kFit, xFit, 0, 1); // fills cell.nbar as a side effect
  vector<GlauberCell> sortedCells = g_cells;
  sort(sortedCells.begin(), sortedCells.end(), [](const GlauberCell& a, const GlauberCell& b){ return a.nbar > b.nbar; });

  double totalKeptWeight = 0.0;
  for(size_t c = 0; c < sortedCells.size(); c++) totalKeptWeight += sortedCells[c].weight;

  cout << endl << "Derived centrality bins (NBD-Glauber-fit-calibrated):" << endl;
  printf("%-10s%14s%14s%14s%12s\n", "Cent", "<Npart>", "<Nspec>", "<Ncoll>", "N_MCevt");

  vector<pair<string, OldGlauberRow>> oldGlauber = ReadOldGeometricGlauber(a_glauberRoot);

  vector<double> binNpartMean(N_CENT_BINS), binNspecMean(N_CENT_BINS), binNcollMean(N_CENT_BINS);
  vector<double> binNevt(N_CENT_BINS);
  vector<double> binNbarLo(N_CENT_BINS), binNbarHi(N_CENT_BINS);

  size_t idx = 0;
  double cumWeight = 0.0;
  for(int i = 0; i < N_CENT_BINS; i++){
    double targetHi = CENT_EDGES_FRAC[i + 1] * totalKeptWeight;
    double sumNpart = 0, sumNspec = 0, sumNcoll = 0, sumW = 0;
    double nbarLo = -1, nbarHi = -1;
    while(idx < sortedCells.size() && cumWeight < targetHi){
      sumNpart += sortedCells[idx].weight * sortedCells[idx].Npart;
      sumNspec += sortedCells[idx].weight * sortedCells[idx].Nspec;
      sumNcoll += sortedCells[idx].weight * sortedCells[idx].Ncoll;
      sumW += sortedCells[idx].weight;
      if(nbarHi < 0) nbarHi = sortedCells[idx].nbar;
      nbarLo = sortedCells[idx].nbar;
      cumWeight += sortedCells[idx].weight;
      idx++;
    }
    binNpartMean[i] = (sumW > 0) ? sumNpart / sumW : 0.0;
    binNspecMean[i] = (sumW > 0) ? sumNspec / sumW : 0.0;
    binNcollMean[i] = (sumW > 0) ? sumNcoll / sumW : 0.0;
    binNevt[i] = sumW;
    binNbarLo[i] = nbarLo; binNbarHi[i] = nbarHi;
    printf("%-10s%14.2f%14.2f%14.2f%12.0f\n", CENT_LABELS[i], binNpartMean[i], binNspecMean[i], binNcollMean[i], sumW);
  }

  cout << endl << "Derived refMult edges (fitted-nbar boundary of each bin, approximate --"
       << " read the actual refMult value on the overlay plot at each dashed line):" << endl;
  for(int i = 0; i < N_CENT_BINS; i++){
    cout << "  " << CENT_LABELS[i] << ": nbar in [" << binNbarLo[i] << ", " << binNbarHi[i] << ")" << endl;
  }
  cout << endl << "Compare against SetCutClass.C's currently-assumed refMult cuts: {44, 37, 28, 17, 5, 0}"
       << " (see EstimateFireballFractionFromRefMult_OO200.C's header for that mapping)." << endl;

  if(!oldGlauber.empty()){
    cout << endl << "Old geometric (b-percentile) vs. new NBD-fit-derived <Npart>/<Nspec>:" << endl;
    printf("%-10s%16s%16s%16s%16s\n", "Cent", "Npart(old)", "Npart(new)", "Nspec(old)", "Nspec(new)");
    for(int i = 0; i < N_CENT_BINS; i++){
      double oldNpart = 0, oldNspec = 0; bool found = false;
      for(size_t j = 0; j < oldGlauber.size(); j++){
        if(oldGlauber[j].first == CENT_LABELS[i]){ oldNpart = oldGlauber[j].second.meanNpart; oldNspec = oldGlauber[j].second.meanNspec; found = true; break; }
      }
      if(found) printf("%-10s%16.2f%16.2f%16.2f%16.2f\n", CENT_LABELS[i], oldNpart, binNpartMean[i], oldNspec, binNspecMean[i]);
      else      printf("%-10s%16s%16.2f%16s%16.2f\n", CENT_LABELS[i], "n/a", binNpartMean[i], "n/a", binNspecMean[i]);
    }
  }

  //--------------------------------------------------------------------------
  // Overlay plot: data vs. best-fit prediction, full histogram range, log-y.
  //--------------------------------------------------------------------------
  int plotLo = 0;
  int plotHi = hRefMult->GetXaxis()->GetXmax();
  vector<double> fullPred = BuildPrediction(g_cells, nppFit, kFit, xFit, plotLo, plotHi);

  TH1D* hPred = new TH1D("hPred", "", plotHi - plotLo, plotLo, plotHi);
  for(int m = plotLo; m < plotHi; m++) hPred->SetBinContent(m - plotLo + 1, scaleAtBest * fullPred[m - plotLo]);
  hPred->SetLineColor(COLOR_FIT);
  hPred->SetLineWidth(2);

  hRefMult->SetMarkerStyle(20);
  hRefMult->SetMarkerColor(COLOR_DATA);
  hRefMult->SetLineColor(COLOR_DATA);
  hRefMult->GetXaxis()->SetTitle("refMult");
  hRefMult->GetYaxis()->SetTitle("Events");
  hRefMult->SetTitle("O+O #sqrt{s_{NN}}=200 GeV -- refMult vs. NBD-Glauber Fit");

  TCanvas* c = new TCanvas("cNBDGlauberFit", "NBD-Glauber Fit to refMult", 900, 650);
  c->SetLogy();
  // Zoom to where there's actually data -- the booked histogram runs to 1000
  // but a small O+O system never gets near that; without this the fit region
  // and edge lines are squeezed into an unreadable sliver at the left edge.
  double xZoomHi = 1.6 * g_fitHi;
  hRefMult->GetXaxis()->SetRangeUser(0.0, xZoomHi);
  hRefMult->Draw("PE");
  hPred->Draw("HIST SAME");

  // mark the derived edges: reconstruct refMult values from the fitted
  // prediction's own cumulative (descending) distribution, matching the same
  // percentile boundaries used for the Npart/Nspec table above.
  vector<double> cumFromHigh(plotHi - plotLo, 0.0);
  double running = 0.0;
  for(int m = plotHi - 1; m >= plotLo; m--){
    running += scaleAtBest * fullPred[m - plotLo];
    cumFromHigh[m - plotLo] = running;
  }
  double totalPredArea = running;
  double ymax = hRefMult->GetMaximum();
  vector<int> edgeRefMult(N_CENT_BINS, plotLo);
  for(int i = 0; i < N_CENT_BINS; i++){
    // cumFromHigh[m] = fraction of events with refMult >= m (maximal at
    // m=plotLo, decreasing as m increases -- see the loop that filled it
    // above). Bin i covers the most-central CENT_EDGES_FRAC[i+1] fraction of
    // events (e.g. i=0 -> top 5%), so its lower/peripheral-side edge is the
    // refMult where the "at-least-this-central" cumulative first drops to
    // CENT_EDGES_FRAC[i+1] of the total -- NOT its complement. (A first pass
    // at this used (1.0 - targetFrac) here, which finds the edge for the
    // complementary fraction instead -- caught after a real run: it put every
    // derived edge down near refMult 0-30, not tracking SetCutClass.C's
    // {44,37,28,17,5,0} scale at all, and collapsed the last bin to whatever
    // plotLo initialized it to since the complementary threshold for the
    // 100% bin is 0 and is essentially never satisfied by an m>plotLo.)
    double targetFrac = CENT_EDGES_FRAC[i + 1];
    for(int m = plotLo; m < plotHi; m++){
      if(cumFromHigh[m - plotLo] <= targetFrac * totalPredArea){ edgeRefMult[i] = m; break; }
    }
    TLine* ln = new TLine(edgeRefMult[i], 0.5, edgeRefMult[i], ymax);
    ln->SetLineColor(COLOR_EDGE);
    ln->SetLineStyle(2);
    ln->Draw();
  }

  TLegend* leg = new TLegend(0.55, 0.65, 0.88, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hRefMult, "Data (measured refMult)", "pe");
  leg->AddEntry(hPred, "NBD-Glauber fit (scaled)", "l");
  leg->Draw();

  TLatex fitLabel;
  fitLabel.SetNDC();
  fitLabel.SetTextSize(0.028);
  fitLabel.DrawLatex(0.55, 0.60, Form("npp=%.3f  k=%.2f  x=%.2f", nppFit, kFit, xFit));
  fitLabel.DrawLatex(0.55, 0.56, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2min, ndof, chi2ndf));

  string outPng = a_outputDir + "/OO200_NBDGlauberFit_RefMult.png";
  c->SaveAs(outPng.c_str());
  cout << "Wrote " << outPng << endl;

  cout << endl << "Derived refMult edges from the fit (compare to SetCutClass.C's {44,37,28,17,5,0}):" << endl;
  for(int i = 0; i < N_CENT_BINS; i++) cout << "  " << CENT_LABELS[i] << ": refMult >= " << edgeRefMult[i] << endl;

  //--------------------------------------------------------------------------
  // Write glauber_oo200_results_NBDfit.root -- CHANGED from JSON to ROOT (see
  // header). TTree "GlauberCentralityBins" uses the SAME schema as
  // GenerateGlauberMC_OO200.C's own output tree (meanB_fm/NpartRMS/NspecRMS
  // filled as -999 sentinels here, since this fit doesn't compute them; every
  // downstream macro reading a GlauberCentralityBins tree can point at either
  // file interchangeably). TTree "GlauberConfig" carries the fit parameters,
  // matching the former JSON "config" sub-object -- informational only, not
  // read back downstream.
  //--------------------------------------------------------------------------
  string rootPath = a_outputDir + "/glauber_oo200_results_NBDfit.root";
  TFile* fitOutFile = new TFile(rootPath.c_str(), "RECREATE");
  fitOutFile->cd();

  TTree* tBins = new TTree("GlauberCentralityBins",
                            "NBD-Glauber fit to real refMult -- <Npart>/<Nspec>/<Ncoll> per bin, "
                            "derived from the fitted model extrapolated to refMult=0");
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

  for(int i = 0; i < N_CENT_BINS; i++){
    strncpy(cent_out, CENT_LABELS[i], sizeof(cent_out) - 1);
    cent_out[sizeof(cent_out) - 1] = '\0';
    nEvents_out       = binNevt[i];
    meanB_out         = -999.0; // not computed by this fit -- see header
    meanNpart_out     = binNpartMean[i];
    meanNspec_out     = binNspecMean[i];
    meanNcoll_out     = binNcollMean[i];
    npartRMS_out      = -999.0; // not computed by this fit -- see header
    nspecRMS_out      = -999.0;
    fracNpart_out     = binNpartMean[i] / 32.0; // 2*A_NUCLEONS, matches GenerateGlauberMC_OO200.C's convention
    fracNspec_out     = binNspecMean[i] / 32.0;
    refMultEdgeHi_out = (double) edgeRefMult[i];
    tBins->Fill();
  }
  tBins->Write();

  TTree* tConfig = new TTree("GlauberConfig", "NBD-Glauber fit parameters and configuration");
  Char_t   model_out[64];
  Char_t   glauberRoot_out[256], yieldFile_out[256];
  Int_t    refMultFitLo_out, refMultFitHi_out, ndof_out;
  Double_t npp_out, nppErr_out, k_out, kErr_out, x_out, xErr_out, chi2_out, chi2ndf_out;
  strncpy(model_out, "NBD-Glauber fit to real refMult (two-component ancestor model)", sizeof(model_out) - 1);
  model_out[sizeof(model_out) - 1] = '\0';
  strncpy(glauberRoot_out, a_glauberRoot.c_str(), sizeof(glauberRoot_out) - 1);
  glauberRoot_out[sizeof(glauberRoot_out) - 1] = '\0';
  strncpy(yieldFile_out, a_yieldFile.c_str(), sizeof(yieldFile_out) - 1);
  yieldFile_out[sizeof(yieldFile_out) - 1] = '\0';
  tConfig->Branch("model",          model_out,        "model/C");
  tConfig->Branch("glauber_root",   glauberRoot_out,  "glauber_root/C");
  tConfig->Branch("yield_file",     yieldFile_out,    "yield_file/C");
  tConfig->Branch("refMultFitLo",   &refMultFitLo_out, "refMultFitLo/I");
  tConfig->Branch("refMultFitHi",   &refMultFitHi_out, "refMultFitHi/I");
  tConfig->Branch("npp",            &npp_out,          "npp/D");
  tConfig->Branch("npp_err",        &nppErr_out,       "npp_err/D");
  tConfig->Branch("k",              &k_out,            "k/D");
  tConfig->Branch("k_err",          &kErr_out,         "k_err/D");
  tConfig->Branch("x",              &x_out,            "x/D");
  tConfig->Branch("x_err",          &xErr_out,         "x_err/D");
  tConfig->Branch("chi2",           &chi2_out,         "chi2/D");
  tConfig->Branch("ndof",           &ndof_out,         "ndof/I");
  tConfig->Branch("chi2_per_ndof",  &chi2ndf_out,      "chi2_per_ndof/D");

  refMultFitLo_out = g_fitLo; refMultFitHi_out = g_fitHi;
  npp_out = nppFit; nppErr_out = err[0];
  k_out = kFit; kErr_out = err[1];
  x_out = xFit; xErr_out = err[2];
  chi2_out = chi2min; ndof_out = ndof; chi2ndf_out = chi2ndf;
  tConfig->Fill();
  tConfig->Write();

  fitOutFile->Close();
  cout << "Wrote " << rootPath << " (TTrees \"GlauberCentralityBins\" [" << N_CENT_BINS
       << " entries], \"GlauberConfig\" [1 entry])" << endl;

  cout << endl << "Reminder: this is a first-pass NBD-Glauber calibration -- check the chi2/ndf,"
       << " whether any parameter sits on its limit, and the overlay plot before trusting the"
       << " derived Npart/Nspec/refMult-edge numbers for anything downstream." << endl;
}
