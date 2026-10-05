// EstimateFireballFractionFromRefMult_OO200.C -- lightflavorspectra_OO200
//
// PURPOSE: a FIRST-PASS particle-COUNT fraction (explicitly NOT yet a dose
// fraction -- see caveat below) for fireball-produced particles vs. participant
// vs. spectator nucleons, per centrality bin, using ONLY data that already exists
// today -- no fit_output.root, no SpectraFitter run, no resolution of the
// StRefMultCorr/EPD open question needed:
//   fireball proxy      = <refMult> (mean measured charged-track multiplicity per
//                          event), read directly from PicoBinner's own output.
//   participant proxy   = <Npart>  (Woods-Saxon Glauber MC, from
//                          GenerateGlauberMC_OO200.C)
//   spectator proxy     = <Nspec>  (same Glauber result)
//
// WHY THIS IS A PARTICLE-COUNT FRACTION, NOT A DOSE FRACTION: it treats one
// reconstructed charged track and one transported/spectator nucleon as
// contributing EQUALLY -- "one unit" each. Real dose depends on species and
// energy (a proton or antiproton deposits more than a pion; annihilation adds
// more for antiprotons; a Glauber Npart nucleon and a refMult track aren't even
// the same KIND of "particle"). Converting this to an actual dose fraction needs
// the per-species, per-energy treatment already built in
// ExtractDoseFractionSummary_OO200.C / ExtractDoubleDifferentialCrossSections_OO200.C
// -- this macro is a deliberately cruder, but immediately computable, complement
// to those, not a replacement. Andrew confirmed this scope explicitly ("don't
// care about specific species yet") before this was written.
//
// A SECOND CAVEAT -- ACCEPTANCE MISMATCH: refMult counts charged tracks in
// whatever limited detector acceptance/pT threshold STAR's refMult definition
// uses (TPC pseudorapidity + track-quality cuts -- the exact cuts weren't
// re-derived from source for this macro). Npart/Nspec, by contrast, are full-4pi
// NUCLEON counts from the Glauber MC. So this fraction combines a
// limited-acceptance PARTICLE count with a full-phase-space NUCLEON count --
// another reason to read this as a rough first pass, not a precision result.
//
// REFMULT HISTOGRAM -- SOURCED DIRECTLY FROM PicoBinner.cxx, NOT GUESSED:
//   Booked:  source/PicoBinner.cxx:396  -- new TH1I("refMult","Multiplicity; refMult ;
//            Number of Events",1000,0,1000)
//   Filled:  source/PicoBinner.cxx:1452 -- refMult->Fill(refMultVar), where
//            refMultVar = event->refMult() (line 1433) -- RAW refMult, since
//            _STREFMULTCORR_ is not defined in this build (matches
//            macros/SetCutClass.C's documented choice).
//   Written: source/PicoBinner.cxx:1977 -- HistogramUtilities::ConditionalWrite(refMult),
//            called right after outFile->cd() with NO intervening mkdir/cd, so it
//            lives at the TOP LEVEL ("/") of the output file, not in a subdirectory.
//            Cross-confirmed by macros/PresentEventQA.C:44's own retrieval:
//            (TH1*) inFile->Get("refMult").
// So: TH1* refMult = (TH1*) inFile->Get("refMult");  -- exactly what this macro does.
//
// CENTRALITY-ON-REFMULT BIN EDGES AND PARTICIPANT/SPECTATOR NUMBERS -- READ AT
// RUNTIME FROM a_glauberRoot, NOT HARDCODED (changed from this macro's first
// version, which hardcoded six copy-pasted Npart/Nspec pairs and SetCutClass.C's
// assumed edges -- that version's own header already flagged the staleness
// risk: "re-derive/update these six pairs of numbers if that file is
// regenerated with different statistics." It since was regenerated, by a
// proper NBD-Glauber fit to this analysis's own real refMult data (see
// FitGlauberNBDToRefMult_OO200.C) rather than a purely geometric b-percentile
// MC -- reading the ROOT file live instead of re-copying six numbers by hand
// removes that staleness risk for good.):
//   a_glauberRoot defaults to glauber_oo200_results_NBDfit.root, which carries
//   BOTH the refMult edge ("refMult_edge_hi", despite the name this is the
//   LOWER/minimum refMult required to belong to this bin or a more central
//   one -- same convention as SetCutClass.C's centCutsArray, see that file's
//   CutClass::centralityIndex() logic) AND mean_Npart/mean_Nspec per bin, so
//   both halves of this macro's output are now self-consistent (same fitted
//   centrality definition), instead of mixing SetCutClass.C's assumed edges
//   for the refMult side with a separately-computed geometric Glauber MC for
//   the nucleon side. If a_glauberRoot doesn't carry refMult_edge_hi (e.g. you
//   point this at the older purely-geometric glauber_oo200_results_woodssaxon.root
//   instead), this macro falls back to SetCutClass.C's assumed edges
//   ({44,37,28,17,5,0}, still real/documented, just not from this file) with a
//   printed warning -- but mean_Npart/mean_Nspec are REQUIRED from the ROOT
//   file; if the file or a bin is missing, this macro aborts loudly rather
//   than falling back to stale numbers (same no-fabrication policy as every
//   other macro in this pipeline).
//
// KNOWN CAVEAT ON THE DEFAULT FILE'S 0-5%/5-10% BINS: the NBD-Glauber fit that
// produced glauber_oo200_results_NBDfit.root visibly overshoots the real
// refMult distribution's tail (refMult gtr-sim 50, growing through 90+) -- almost
// entirely inside the 0-5% bin (refMult >= 43) and the top of 5-10%. The
// mean_Npart/mean_Nspec for those two bins come from ranking Glauber MC cells
// by fitted nbar and cutting by cross-section-weight percentile, which is
// mostly set by the well-fit bulk region rather than reading the tail shape
// directly -- so it's not as exposed as the refMult edge itself, but treat the
// 0-5%/5-10% numbers here with that in mind, more than the well-fit 10-80%
// bins. See this session's fit diagnostics for the full discussion.
//
// INPUT: a_yieldFile -- any PicoBinner output file with a top-level "refMult"
// histogram, e.g. the real file this was written against:
//   /Users/aliggett/data/OO/yieldHistos_OO200_pion.root
// (any of the pion/kaon/proton yield files should carry an identical refMult
// histogram, since it's an event-level quantity independent of a_partIndex --
// worth spot-checking that they actually agree if you have more than one).
//
// OUTPUT: a_outputFile (default FireballFractionFromRefMult_Summary.root), containing
// one TTree named "FireballFractionFromRefMult" with one entry per centrality bin and
// branches: centIndex/I, centLabel/C, refMultLoEdge/D, refMultHiEdge/D, meanRefMult/D,
// nEventsInBin/D, meanNpart/D, meanNspec/D, fireball_frac/D, participant_frac/D,
// spectator_frac/D. Read it back with e.g.
//   TFile* f = TFile::Open("FireballFractionFromRefMult_Summary.root");
//   TTree* t = (TTree*) f->Get("FireballFractionFromRefMult");
//   t->Scan("centLabel:fireball_frac:participant_frac:spectator_frac");
// Also prints, per bin, what fraction of TOTAL events fell in that refMult range --
// compare against this analysis's intended {5,10,20,40,80,100}% design as a sanity
// check that PicoBinner's actual data lands where the cuts were meant to put it.
//
// USAGE:
//   root -l -b -q 'EstimateFireballFractionFromRefMult_OO200.C("yieldHistos_OO200_pion.root")'
// (writes FireballFractionFromRefMult_Summary.root in the current directory,
// reading centrality edges + Npart/Nspec from glauber_oo200_results_NBDfit.root
// in the current directory by default -- pass a third string argument to point
// at a different Glauber ROOT file, or a second argument to change the output path)
//
// UNTESTED -- same caveat as every ROOT macro in this directory: no ROOT in this
// sandbox. The Glauber-table reader was just switched from hand-rolled JSON
// line-scanning to a TTree read (Python/JSON retired from this pipeline, per
// Andrew's instruction) -- this file's own TTree-writing output logic is
// unchanged and was already real-data-tested, but the READER side of this
// specific change has not been re-exercised. Re-run after regenerating
// glauber_oo200_results_NBDfit.root and confirm the printed <Npart>/<Nspec>
// per bin match what the old JSON-based run reported.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

const int N_CENT_BINS = 6;
const char* CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};

// Fallback-only refMult edges, used ONLY if a_glauberRoot doesn't carry
// "refMult_edge_hi" for a bin (e.g. pointed at the purely-geometric
// glauber_oo200_results_woodssaxon.root instead of the NBD-fit one) --
// SetCutClass.C's centCutsArray, see the header for the full derivation.
// The corresponding upper edges aren't needed as a separate table: refMultHi[i]
// is always just refMultLo[i-1] (or -1/no-limit for the 0-5% bin), computed
// below from whichever refMultLo[] values end up in use.
const double DEFAULT_REFMULT_LO[N_CENT_BINS] = {44, 37, 28, 17, 5, 0};

//============================================================================
// One row of the Glauber centrality-bin table: cent label, mean_Npart,
// mean_Nspec, and refMultEdgeHi (sentinel -999 for the purely-geometric
// woodssaxon file, real for the NBD-fit file) -- matches the shared
// "GlauberCentralityBins" TTree schema both glauber_oo200_results_NBDfit.root
// and glauber_oo200_results_woodssaxon.root write (see
// GenerateGlauberMC_OO200.C's header for the full schema). Replaces the
// former hand-rolled JSON line-scanning reader now that Python/JSON are
// retired from this pipeline.
//============================================================================
struct GlauberBinRow {
  string cent;
  double meanNpart = 0.0, meanNspec = 0.0, refMultEdgeHi = -999.0;
  bool foundNpart = false, foundNspec = false, foundEdge = false;
};

vector<GlauberBinRow> ReadGlauberCentralityTable(const string& a_path){
  vector<GlauberBinRow> out;
  TFile* f = TFile::Open(a_path.c_str(), "READ");
  if(!f || f->IsZombie()) return out; // caller checks emptiness and decides how to react

  TTree* t = (TTree*) f->Get("GlauberCentralityBins");
  if(!t){ f->Close(); return out; }

  Char_t cent_in[16];
  Double_t meanNpart_in, meanNspec_in, refMultEdgeHi_in;
  t->SetBranchAddress("cent", cent_in);
  t->SetBranchAddress("meanNpart", &meanNpart_in);
  t->SetBranchAddress("meanNspec", &meanNspec_in);
  t->SetBranchAddress("refMultEdgeHi", &refMultEdgeHi_in);

  Long64_t nEntries = t->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    GlauberBinRow row;
    row.cent = string(cent_in);
    row.meanNpart = meanNpart_in;
    row.meanNspec = meanNspec_in;
    row.foundNpart = true;
    row.foundNspec = true;
    // refMultEdgeHi/D is always present as a branch, but carries the -999
    // sentinel for a purely-geometric file (see GenerateGlauberMC_OO200.C) --
    // treat that sentinel the same way the old "field absent" case was
    // treated: foundEdge=false, triggering this macro's documented fallback
    // to SetCutClass.C's assumed edges.
    row.refMultEdgeHi = refMultEdgeHi_in;
    row.foundEdge = (refMultEdgeHi_in > -998.0);
    out.push_back(row);
  }
  f->Close();
  return out;
}

void EstimateFireballFractionFromRefMult_OO200(string a_yieldFile,
                                                 string a_outputFile  = "FireballFractionFromRefMult_Summary.root",
                                                 string a_glauberRoot = "glauber_oo200_results_NBDfit.root"){

  //--------------------------------------------------------------------------
  // Load the Glauber centrality-bin table (Npart/Nspec required; refMult
  // edges optional, see header) and match it onto CENT_LABELS order. No
  // fabrication: abort if the file or any bin's Npart/Nspec is missing.
  //--------------------------------------------------------------------------
  vector<GlauberBinRow> glauberRows = ReadGlauberCentralityTable(a_glauberRoot);
  if(glauberRows.empty()){
    cerr << "ERROR: could not read any centrality-bin rows from '" << a_glauberRoot << "'."
         << " Run FitGlauberNBDToRefMult_OO200.C (or GenerateGlauberMC_OO200.C for the purely"
         << " geometric file) first, or pass the right path as the 3rd argument." << endl;
    return;
  }

  double meanNpart[N_CENT_BINS], meanNspec[N_CENT_BINS];
  double refMultLo[N_CENT_BINS], refMultHi[N_CENT_BINS];
  bool usedFallbackEdges = false;
  for(int i = 0; i < N_CENT_BINS; i++){
    const GlauberBinRow* match = NULL;
    for(size_t j = 0; j < glauberRows.size(); j++){
      if(glauberRows[j].cent == CENT_LABELS[i]){ match = &glauberRows[j]; break; }
    }
    if(!match || !match->foundNpart || !match->foundNspec){
      cerr << "ERROR: '" << a_glauberRoot << "' has no mean_Npart/mean_Nspec for cent bin '"
           << CENT_LABELS[i] << "' -- aborting rather than guessing." << endl;
      return;
    }
    meanNpart[i] = match->meanNpart;
    meanNspec[i] = match->meanNspec;

    if(match->foundEdge){
      refMultLo[i] = match->refMultEdgeHi;
    } else {
      refMultLo[i] = DEFAULT_REFMULT_LO[i];
      usedFallbackEdges = true;
    }
  }
  for(int i = 0; i < N_CENT_BINS; i++){
    refMultHi[i] = (i == 0) ? -1.0 : refMultLo[i - 1];
  }
  if(usedFallbackEdges){
    cout << "NOTE: '" << a_glauberRoot << "' has no refMult_edge_hi field for at least one bin --"
         << " falling back to SetCutClass.C's assumed edges {44,37,28,17,5,0} for those." << endl;
  }
  cout << "Loaded Glauber centrality table from '" << a_glauberRoot << "':" << endl;
  for(int i = 0; i < N_CENT_BINS; i++){
    cout << "  " << CENT_LABELS[i] << ": refMult >= " << refMultLo[i]
         << "   <Npart>=" << meanNpart[i] << "   <Nspec>=" << meanNspec[i] << endl;
  }

  TFile* inFile = TFile::Open(a_yieldFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){
    cerr << "ERROR: could not open input file '" << a_yieldFile << "'" << endl;
    return;
  }

  TH1* hRefMult = (TH1*) inFile->Get("refMult");
  if(!hRefMult){
    cerr << "ERROR: no top-level 'refMult' histogram found in '" << a_yieldFile << "'." << endl;
    cerr << "       (Confirm this is a PicoBinner output file -- see header for how it's booked.)" << endl;
    inFile->Close();
    return;
  }

  double totalEvents = hRefMult->Integral(1, hRefMult->GetNbinsX());
  if(totalEvents <= 0){
    cerr << "ERROR: refMult histogram is empty." << endl;
    inFile->Close();
    return;
  }

  cout << "==================================================================" << endl;
  cout << " EstimateFireballFractionFromRefMult_OO200" << endl;
  cout << " Input: " << a_yieldFile << "   Total events in refMult histo: " << totalEvents << endl;
  cout << "==================================================================" << endl;

  TFile* outFile = new TFile(a_outputFile.c_str(), "RECREATE");
  outFile->cd();

  TTree* tree = new TTree("FireballFractionFromRefMult",
                           "Fireball/participant/spectator PARTICLE-COUNT fractions "
                           "(refMult + Woods-Saxon Glauber, not yet dose -- see macro header) "
                           "per centrality bin");

  Int_t  centIndex_out;
  Char_t centLabel_out[16];
  Double_t refMultLoEdge_out, refMultHiEdge_out, meanRefMult_out, nEventsInBin_out;
  Double_t meanNpart_out, meanNspec_out;
  Double_t fireball_frac_out, participant_frac_out, spectator_frac_out;

  tree->Branch("centIndex",        &centIndex_out,        "centIndex/I");
  tree->Branch("centLabel",        centLabel_out,          "centLabel/C");
  tree->Branch("refMultLoEdge",    &refMultLoEdge_out,    "refMultLoEdge/D");
  tree->Branch("refMultHiEdge",    &refMultHiEdge_out,    "refMultHiEdge/D");
  tree->Branch("meanRefMult",      &meanRefMult_out,      "meanRefMult/D");
  tree->Branch("nEventsInBin",     &nEventsInBin_out,     "nEventsInBin/D");
  tree->Branch("meanNpart",        &meanNpart_out,        "meanNpart/D");
  tree->Branch("meanNspec",        &meanNspec_out,        "meanNspec/D");
  tree->Branch("fireball_frac",    &fireball_frac_out,    "fireball_frac/D");
  tree->Branch("participant_frac", &participant_frac_out, "participant_frac/D");
  tree->Branch("spectator_frac",   &spectator_frac_out,   "spectator_frac/D");

  for(int centIndex = 0; centIndex < N_CENT_BINS; centIndex++){

    double lo = refMultLo[centIndex];
    double hi = refMultHi[centIndex];

    int loBin = hRefMult->GetXaxis()->FindBin(lo);
    int hiBin = (hi < 0) ? hRefMult->GetNbinsX() : (hRefMult->GetXaxis()->FindBin(hi) - 1);
    if(hiBin < loBin){
      cerr << "  WARNING: centIndex " << centIndex << " -- empty bin range, skipping." << endl;
      continue;
    }

    double nEventsInBin = hRefMult->Integral(loBin, hiBin);
    double sumWeighted = 0.0;
    for(int b = loBin; b <= hiBin; b++){
      sumWeighted += hRefMult->GetBinContent(b) * hRefMult->GetXaxis()->GetBinCenter(b);
    }
    double meanRefMult = (nEventsInBin > 0) ? (sumWeighted / nEventsInBin) : 0.0;

    double meanNpartBin = meanNpart[centIndex];
    double meanNspecBin = meanNspec[centIndex];

    double total = meanRefMult + meanNpartBin + meanNspecBin;
    double fireball_frac    = (total > 0) ? meanRefMult / total : 0.0;
    double participant_frac = (total > 0) ? meanNpartBin / total : 0.0;
    double spectator_frac    = (total > 0) ? meanNspecBin / total : 0.0;

    double eventFracOfTotal = 100.0 * nEventsInBin / totalEvents;

    cout << "------------------------------------------------------------------" << endl;
    cout << " Cent " << CENT_LABELS[centIndex] << " (refMult in ["
         << lo << "," << (hi < 0 ? "inf" : Form("%.0f", hi)) << ")): "
         << nEventsInBin << " events (" << eventFracOfTotal << "% of total, design target "
         << (centIndex == 0 ? 5 : centIndex == 1 ? 5 : centIndex == 2 ? 10 : centIndex == 3 ? 20 : centIndex == 4 ? 40 : 20)
         << "%)" << endl;
    cout << "   <refMult> = " << meanRefMult << "   fireball=" << 100.0*fireball_frac
         << "%  participant=" << 100.0*participant_frac << "%  spectator=" << 100.0*spectator_frac << "%" << endl;

    centIndex_out = centIndex;
    strncpy(centLabel_out, CENT_LABELS[centIndex], sizeof(centLabel_out) - 1);
    centLabel_out[sizeof(centLabel_out) - 1] = '\0';
    refMultLoEdge_out    = lo;
    refMultHiEdge_out    = hi;
    meanRefMult_out       = meanRefMult;
    nEventsInBin_out      = nEventsInBin;
    meanNpart_out          = meanNpartBin;
    meanNspec_out          = meanNspecBin;
    fireball_frac_out     = fireball_frac;
    participant_frac_out  = participant_frac;
    spectator_frac_out     = spectator_frac;
    tree->Fill();
  }

  outFile->cd();
  tree->Write();
  outFile->Close();
  inFile->Close();
  cout << "==================================================================" << endl;
  cout << "Wrote " << a_outputFile << " (TTree \"FireballFractionFromRefMult\", "
       << N_CENT_BINS << " entries)" << endl;
  cout << "Check the event-fraction-vs-design-target line above first -- if PicoBinner's" << endl;
  cout << "actual refMult distribution doesn't land close to {5,5,10,20,40,20}% per bin," << endl;
  cout << "something upstream (cut config, luminosity mix) doesn't match what SetCutClass.C assumed." << endl;
}
