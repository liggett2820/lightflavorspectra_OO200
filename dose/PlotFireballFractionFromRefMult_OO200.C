// PlotFireballFractionFromRefMult_OO200.C -- lightflavorspectra_OO200
//
// PURPOSE: plots the fireball/participant/spectator PARTICLE-COUNT fraction
// (refMult-based proxy, per EstimateFireballFractionFromRefMult_OO200.C) vs.
// centrality, using the TTree that macro already wrote to
// FireballFractionFromRefMult_Summary.root -- this is REAL data that has
// already been run once against dose/FireballFractionFromRefMult_Summary.root
// in this repo (confirmed present on disk).
//
// WHY THIS EXISTS INSTEAD OF PlotDoseFractionComparison_OO200.C: that macro
// needs DoseFractionSummary.csv, which needs ExtractDoseFractionSummary_OO200.C,
// which needs fit_output.root -- the output of macros/RunSpectraFitter.C run
// against macros/RunRawSpectraModifier.C's spectra_modified.root. NEITHER of
// those exists anywhere in this repo yet (checked: no fit_output.root,
// no spectra_modified.root, no NetProtonStopping_Summary.csv, no
// DoseFractionSummary.csv on disk as of this macro being written). That chain
// needs a full corrected-spectra STAR analysis pass (embedding etc.) that
// hasn't been run for O+O 200 yet. This macro does NOT fabricate that missing
// result -- it plots the different, cruder, but already-real
// FireballFractionFromRefMult_Summary.root output instead, clearly labeled as
// a particle-count proxy, not a true dose fraction (same distinction the
// source macro's own header draws -- see its "WHY THIS IS A PARTICLE-COUNT
// FRACTION, NOT A DOSE FRACTION" section). Once fit_output.root exists, re-run
// ExtractNetProtonStoppingFraction_OO200.C -> ExtractDoseFractionSummary_OO200.C
// -> PlotDoseFractionComparison_OO200.C for the real dose-weighted result and
// treat this plot as the interim/first-pass check.
//
// Reads: a_summaryFile (default FireballFractionFromRefMult_Summary.root),
// TTree "FireballFractionFromRefMult", branches centIndex/I, centLabel/C,
// refMultLoEdge/D, refMultHiEdge/D, meanRefMult/D, nEventsInBin/D,
// meanNpart/D, meanNspec/D, fireball_frac/D, participant_frac/D,
// spectator_frac/D -- exact schema from
// EstimateFireballFractionFromRefMult_OO200.C's TTree::Branch calls. Prints a
// loud error and returns with NO plot if the tree can't be read -- no
// synthetic-placeholder fallback, matching PlotDoseFractionComparison_OO200.C's
// policy (see that file's header for the earlier incident this guards
// against).
//
// Output: OO200_FireballFractionFromRefMult.png under a_outputDir -- single
// stacked-bar panel, fireball/participant/spectator fraction by centrality,
// same color convention as PlotDoseFractionComparison_OO200.C (brick red
// fireball #C0392B, navy participant #1E2761, amber spectator #E8A33D) for
// visual continuity across this repo's figures. A caveat line is drawn
// directly on the canvas so this can't be mistaken for the dose-weighted
// result if it's viewed on its own.
//
// USAGE (ROOT, from dose/ where FireballFractionFromRefMult_Summary.root
// already lives):
//   root -l -b -q 'PlotFireballFractionFromRefMult_OO200.C()'
// or with explicit paths:
//   root -l -b -q 'PlotFireballFractionFromRefMult_OO200.C("FireballFractionFromRefMult_Summary.root","./plots")'
//
// STATUS: UNTESTED -- no ROOT in the sandbox that wrote this (same caveat as
// every other C++ macro in this session). Run it and report back what prints
// before trusting the figure.

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;

//============================================================================
// Palette -- identical hex values to PlotDoseFractionComparison_OO200.C, see
// that file's header for the CVD-validation provenance.
//============================================================================
const int COLOR_FIREBALL      = TColor::GetColor("#C0392B"); // brick red
const int COLOR_PARTICIPANT   = TColor::GetColor("#1E2761"); // navy
const int COLOR_SPECTATOR     = TColor::GetColor("#E8A33D"); // amber
const int COLOR_TEXT_ON_DARK  = kWhite;
const int COLOR_TEXT_ON_LIGHT = TColor::GetColor("#1A1A1A");
const int COLOR_CAVEAT        = TColor::GetColor("#7A7A7A"); // neutral gray

struct FracRow {
  int    centIndex;
  string centLabel;
  double meanRefMult, nEventsInBin, meanNpart, meanNspec;
  double fireballFrac, participantFrac, spectatorFrac;
};

//============================================================================
// ReadFractionTree -- returns an EMPTY vector (not a fabricated one) on any
// failure. Caller must check and abort -- no synthetic-placeholder fallback,
// see header.
//============================================================================
vector<FracRow> ReadFractionTree(const string& a_path){
  vector<FracRow> rows;

  TFile* f = TFile::Open(a_path.c_str(), "READ");
  if(!f || f->IsZombie()){
    cerr << "ERROR: could not open '" << a_path << "'." << endl;
    return rows;
  }

  TTree* t = (TTree*) f->Get("FireballFractionFromRefMult");
  if(!t){
    cerr << "ERROR: no TTree named 'FireballFractionFromRefMult' found in '" << a_path << "'." << endl;
    f->Close();
    return rows;
  }

  Int_t    centIndex_in;
  Char_t   centLabel_in[16];
  Double_t refMultLoEdge_in, refMultHiEdge_in, meanRefMult_in, nEventsInBin_in;
  Double_t meanNpart_in, meanNspec_in;
  Double_t fireball_frac_in, participant_frac_in, spectator_frac_in;

  t->SetBranchAddress("centIndex",        &centIndex_in);
  t->SetBranchAddress("centLabel",        centLabel_in);
  t->SetBranchAddress("refMultLoEdge",    &refMultLoEdge_in);
  t->SetBranchAddress("refMultHiEdge",    &refMultHiEdge_in);
  t->SetBranchAddress("meanRefMult",      &meanRefMult_in);
  t->SetBranchAddress("nEventsInBin",     &nEventsInBin_in);
  t->SetBranchAddress("meanNpart",        &meanNpart_in);
  t->SetBranchAddress("meanNspec",        &meanNspec_in);
  t->SetBranchAddress("fireball_frac",    &fireball_frac_in);
  t->SetBranchAddress("participant_frac", &participant_frac_in);
  t->SetBranchAddress("spectator_frac",   &spectator_frac_in);

  Long64_t nEntries = t->GetEntries();
  if(nEntries <= 0){
    cerr << "ERROR: TTree 'FireballFractionFromRefMult' in '" << a_path << "' has no entries." << endl;
    f->Close();
    return rows;
  }

  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    FracRow row;
    row.centIndex       = centIndex_in;
    row.centLabel        = string(centLabel_in);
    row.meanRefMult       = meanRefMult_in;
    row.nEventsInBin      = nEventsInBin_in;
    row.meanNpart          = meanNpart_in;
    row.meanNspec          = meanNspec_in;
    row.fireballFrac      = fireball_frac_in;
    row.participantFrac   = participant_frac_in;
    row.spectatorFrac      = spectator_frac_in;
    rows.push_back(row);
  }
  f->Close();

  sort(rows.begin(), rows.end(), [](const FracRow& a, const FracRow& b){ return a.centIndex < b.centIndex; });
  return rows;
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
void PlotFireballFractionFromRefMult_OO200(
    string a_summaryFile = "FireballFractionFromRefMult_Summary.root",
    string a_outputDir   = "."
){
  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outputDir.c_str(), kTRUE);

  vector<FracRow> rows = ReadFractionTree(a_summaryFile);
  if(rows.empty()){
    cerr << "ERROR: no usable rows read from '" << a_summaryFile << "' -- aborting."
         << " (This macro does NOT fabricate placeholder data. Run"
         << " EstimateFireballFractionFromRefMult_OO200.C against a real PicoBinner"
         << " yield file first.)" << endl;
    return;
  }
  cout << "Read " << rows.size() << " centrality bins from '" << a_summaryFile << "'." << endl;

  int n = (int) rows.size();

  TCanvas* c = new TCanvas("cFireballFractionFromRefMult", "Fireball Fraction from refMult", 950, 700);
  gPad->SetBottomMargin(0.24); // room for legend + caveat line below the axis
  gPad->SetGridy();

  TH1D* hStackTop = new TH1D("hStackTop", "", n, 0, n); // fireball+participant+spectator
  TH1D* hStackMid = new TH1D("hStackMid", "", n, 0, n); // fireball+participant
  TH1D* hStackBot = new TH1D("hStackBot", "", n, 0, n); // fireball only

  for(int i = 0; i < n; i++){
    double fb = 100.0 * rows[i].fireballFrac;
    double pt = 100.0 * rows[i].participantFrac;
    double sp = 100.0 * rows[i].spectatorFrac;
    hStackTop->SetBinContent(i+1, fb + pt + sp);
    hStackMid->SetBinContent(i+1, fb + pt);
    hStackBot->SetBinContent(i+1, fb);
    hStackTop->GetXaxis()->SetBinLabel(i+1, rows[i].centLabel.c_str());
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
  hStackTop->GetXaxis()->SetLabelSize(0.045);
  hStackTop->GetYaxis()->SetTitle("Fraction of (#LTrefMult#GT + #LTN_{part}#GT + #LTN_{spec}#GT) [%]");
  hStackTop->GetXaxis()->SetTitle("Centrality");
  hStackTop->SetTitle("O+O #sqrt{s_{NN}}=200 GeV -- Particle-Count Fraction by Source (refMult Proxy)");

  hStackTop->Draw("bar");
  hStackMid->Draw("bar same");
  hStackBot->Draw("bar same");

  for(int i = 0; i < n; i++){
    double x = hStackTop->GetXaxis()->GetBinCenter(i+1);
    double fb = 100.0 * rows[i].fireballFrac;
    double pt = 100.0 * rows[i].participantFrac;
    double sp = 100.0 * rows[i].spectatorFrac;
    DrawPercentLabel(x, fb / 2.0, fb, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt / 2.0, pt, COLOR_TEXT_ON_DARK);
    DrawPercentLabel(x, fb + pt + sp / 2.0, sp, COLOR_TEXT_ON_LIGHT);
  }

  TLegend* leg = new TLegend(0.12, 0.06, 0.90, 0.16);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(1);
  leg->AddEntry(hStackBot, "Fireball-produced (proxy: #LTrefMult#GT)", "f");
  leg->AddEntry(hStackMid, "Participant-associated (proxy: #LTN_{part}#GT, Woods-Saxon Glauber)", "f");
  leg->AddEntry(hStackTop, "Spectator-associated (proxy: #LTN_{spec}#GT, Woods-Saxon Glauber)", "f");
  leg->Draw();

  TLatex caveat;
  caveat.SetNDC();
  caveat.SetTextFont(42);
  caveat.SetTextSize(0.024);
  caveat.SetTextColor(COLOR_CAVEAT);
  caveat.SetTextAlign(21);
  caveat.DrawLatex(0.5, 0.015,
    "Particle-count proxy (one track / one nucleon = one unit) -- NOT a dose fraction. "
    "Real dose split pending fit_output.root (see ExtractDoseFractionSummary_OO200.C).");

  string outPath = a_outputDir + "/OO200_FireballFractionFromRefMult.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;

  // ---- Text table to stdout ----
  cout << endl;
  printf("%-10s%14s%14s%14s%12s%14s%15s%13s\n",
         "Cent", "<refMult>", "<N_part>", "<N_spec>", "N_evt", "Fireball %", "Participant %", "Spectator %");
  for(int i = 0; i < n; i++){
    printf("%-10s%14.2f%14.2f%14.2f%12.0f%12.1f%15.1f%13.1f\n",
           rows[i].centLabel.c_str(), rows[i].meanRefMult, rows[i].meanNpart, rows[i].meanNspec,
           rows[i].nEventsInBin, 100.0*rows[i].fireballFrac, 100.0*rows[i].participantFrac,
           100.0*rows[i].spectatorFrac);
  }
  cout << endl << "Reminder: this is the refMult-based particle-count proxy, not a dose fraction." << endl;
  cout << "See this macro's header and EstimateFireballFractionFromRefMult_OO200.C's header for scope." << endl;
}
