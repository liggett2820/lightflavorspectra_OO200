// PlotDoseFractionComparison_OO200.C -- lightflavorspectra_OO200
//
// ROOT/C++ port of PlotDoseFractionComparison_OO200.py, per Andrew's explicit
// request to keep the dose-fraction pipeline in C++/ROOT going forward rather
// than mixing in Python for the plotting-only step (the Python version was used
// once because this cloud sandbox that wrote it has no ROOT install, and this
// particular step's inputs -- a CSV and a JSON, not a ROOT file -- happened not
// to require one). PlotDoseFractionComparison_OO200.py is now RETIRED --
// this macro fully supersedes it, and the Glauber input below was switched
// from JSON to a native ROOT TTree at the same time (see next paragraph),
// finishing the Python/JSON removal from this whole pipeline.
//
// DELIBERATE DIFFERENCE FROM THE PYTHON VERSION -- READ THIS FIRST:
// The Python script auto-generated a SYNTHETIC PLACEHOLDER CSV and plotted it,
// unprompted, whenever the real DoseFractionSummary.csv wasn't found. That
// silent fallback produced a rendered, real-looking chart from fabricated
// numbers, which is exactly what got caught and had to be deleted earlier in
// this session (see NetProtonStopping/DoseFraction chat history: "so these
// pseudo plots are incorrect then"). This C++ version does NOT do that -- if
// the real CSV is missing, it prints a loud error and returns without
// producing any plot at all. Do not add a synthetic-placeholder fallback back
// into this macro.
//
// Reads:
//   DoseFractionSummary.csv       -- from ExtractDoseFractionSummary_OO200.C.
//     Columns (exact order, see that macro): centIndex,centLabel,sigma_bin_mb,
//     fireball_dNdy,fireball_dNdy_err,participant_dNdy,participant_dNdy_err,
//     spectator_dNdy,spectator_dNdy_err,total_dNdy,fireball_frac,
//     fireball_frac_err,participant_frac,participant_frac_err,spectator_frac,
//     spectator_frac_err,fireball_sigma_mb,participant_sigma_mb,
//     spectator_sigma_mb. Only centIndex/centLabel/the three *_frac columns are
//     used here (matches the Python version's scope -- no error bars drawn on
//     the bars there either, just value labels).
//   glauber_oo200_results_woodssaxon.root -- from GenerateGlauberMC_OO200.C's
//     Woods-Saxon Glauber MC (CHANGED from JSON to ROOT -- see that macro's
//     header for the full schema and why). Only the "GlauberCentralityBins"
//     tree's cent/fracNspecOfTotal branches are used. This file is OPTIONAL:
//     if missing, the bottom comparison panel is skipped with an on-canvas
//     message (same graceful degradation as the Python version) -- this is
//     the one place a missing input does NOT abort the whole macro, since
//     it's a secondary overlay, not the primary result.
//   TTree reader below (ReadGlauberSpecFractions) -- replaces the former
//   hand-rolled line-scanning JSON reader now that Python/JSON are retired
//   from this pipeline.
//
// Output: OO200_DoseFractionComparison.png under a_outputDir, 2-panel:
//   Top:    manually-stacked bar (fireball/participant/spectator dose fraction
//           by centrality). ROOT's THStack has a lazy-axis-creation quirk that
//           makes custom bin labels fiddly, so this draws three full-height-to-
//           short-height histograms in back-to-front order instead (draw the
//           tallest cumulative height first in the TOP segment's color, then
//           progressively shorter cumulative heights on top in each lower
//           segment's color) -- a standard ROOT manual-stack idiom, and it
//           keeps full control of the x-axis bin labels.
//   Bottom: paired bars, data-driven (net-proton) vs geometric (Glauber)
//           spectator-associated fraction, via TH1::SetBarWidth/SetBarOffset.
//
// Colors: identical hex values to the Python version and to
// OO200_Npart_Nspec_byCentrality_WoodsSaxon.png (navy participant, amber
// spectator) for visual continuity across this repo's figures; brick red
// fireball was CVD-validated against both via this repo's dataviz skill in the
// original Python version's development -- same palette, so that validation
// still applies here.
//
// USAGE (ROOT, wherever a real DoseFractionSummary.csv exists):
//   root -l -b -q 'PlotDoseFractionComparison_OO200.C()'
// or with explicit paths:
//   root -l -b -q 'PlotDoseFractionComparison_OO200.C("DoseFractionSummary.csv","glauber_oo200_results_woodssaxon.root","./plots")'
//
// STATUS: UNTESTED (no ROOT in the sandbox that wrote this -- same caveat as
// every other C++ macro in this session). ReadGlauberSpecFractions() was just
// switched from JSON line-scanning to a TTree read and has not been
// exercised against a real glauber_oo200_results_woodssaxon.root -- confirm
// the bottom panel's geometric-spectator bars still populate before trusting
// them. Run it and report back what prints before trusting the figure.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

using namespace std;

//============================================================================
// Palette -- identical hex values to the Python version, see header
//============================================================================
const int COLOR_FIREBALL     = TColor::GetColor("#C0392B"); // brick red
const int COLOR_PARTICIPANT  = TColor::GetColor("#1E2761"); // navy
const int COLOR_SPECTATOR    = TColor::GetColor("#E8A33D"); // amber
const int COLOR_GEOM_SPEC    = TColor::GetColor("#B8B8B8"); // neutral gray
const int COLOR_TEXT_ON_DARK  = kWhite;
const int COLOR_TEXT_ON_LIGHT = TColor::GetColor("#1A1A1A");

//============================================================================
// Row structs
//============================================================================
struct DoseFractionRow {
  int    centIndex;
  string centLabel;
  double fireballFrac, fireballFracErr;
  double participantFrac, participantFracErr;
  double spectatorFrac, spectatorFracErr;
};

struct GlauberSpecRow {
  string centLabel;
  double fracNspecOfTotal;
};

//============================================================================
// ReadDoseFractionCSV -- hand-rolled, matches ExtractDoseFractionSummary_OO200.C's
// exact 19-column schema. Returns an EMPTY vector (not a fabricated one) on any
// failure -- caller must check and abort, per this macro's no-synthetic-fallback
// policy (see header).
//============================================================================
vector<DoseFractionRow> ReadDoseFractionCSV(const string& a_path){
  vector<DoseFractionRow> rows;
  ifstream f(a_path.c_str());
  if(!f.is_open()){
    cerr << "ERROR: could not open '" << a_path << "'" << endl;
    return rows;
  }

  string line;
  if(!getline(f, line)){ // header line, discarded (column order is hardcoded below, matching the known schema)
    cerr << "ERROR: '" << a_path << "' is empty (no header line)." << endl;
    return rows;
  }

  int lineNum = 1;
  while(getline(f, line)){
    lineNum++;
    if(line.empty()) continue;
    vector<string> tok;
    stringstream ss(line);
    string field;
    while(getline(ss, field, ',')) tok.push_back(field);
    if(tok.size() != 19){
      cerr << "WARNING: '" << a_path << "' line " << lineNum << " has " << tok.size()
           << " fields, expected 19 -- skipping this row." << endl;
      continue;
    }
    DoseFractionRow row;
    row.centIndex          = atoi(tok[0].c_str());
    row.centLabel           = tok[1];
    row.fireballFrac        = atof(tok[10].c_str());
    row.fireballFracErr     = atof(tok[11].c_str());
    row.participantFrac     = atof(tok[12].c_str());
    row.participantFracErr  = atof(tok[13].c_str());
    row.spectatorFrac       = atof(tok[14].c_str());
    row.spectatorFracErr    = atof(tok[15].c_str());
    rows.push_back(row);
  }
  f.close();

  // sort by centIndex, matching the Python version's explicit sort
  for(size_t i = 0; i < rows.size(); i++)
    for(size_t j = i+1; j < rows.size(); j++)
      if(rows[j].centIndex < rows[i].centIndex) std::swap(rows[i], rows[j]);

  return rows;
}

//============================================================================
// ReadGlauberSpecFractions -- reads the "GlauberCentralityBins" TTree's
// cent/fracNspecOfTotal branches (see GenerateGlauberMC_OO200.C's header for
// the full shared schema). Replaces the former hand-rolled line-scanning
// JSON reader now that Python/JSON are retired from this pipeline. Returns
// an empty vector if the file/tree is missing (soft-fail -- the bottom panel
// is optional, unlike the CSV above).
//============================================================================
vector<GlauberSpecRow> ReadGlauberSpecFractions(const string& a_path){
  vector<GlauberSpecRow> out;
  TFile* f = TFile::Open(a_path.c_str(), "READ");
  if(!f || f->IsZombie()){
    cout << "WARNING: Glauber ROOT file '" << a_path << "' not found -- bottom comparison panel will be skipped." << endl;
    return out;
  }
  TTree* t = (TTree*) f->Get("GlauberCentralityBins");
  if(!t){
    cout << "WARNING: no \"GlauberCentralityBins\" tree in '" << a_path
         << "' -- bottom comparison panel will be skipped." << endl;
    f->Close();
    return out;
  }

  Char_t cent_in[16];
  Double_t fracNspec_in;
  t->SetBranchAddress("cent", cent_in);
  t->SetBranchAddress("fracNspecOfTotal", &fracNspec_in);

  Long64_t nEntries = t->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    GlauberSpecRow row;
    row.centLabel = string(cent_in);
    row.fracNspecOfTotal = fracNspec_in;
    out.push_back(row);
  }
  f->Close();

  if(out.empty()){
    cout << "WARNING: no entries in '" << a_path
         << "'s \"GlauberCentralityBins\" tree -- bottom comparison panel will be skipped." << endl;
  }
  return out;
}

//============================================================================
// Draws a value-label ("NN%") centered at (x, y) if val is large enough to be
// legible (>=3, matching the Python version's slivers-too-thin-to-label skip).
//============================================================================
void DrawPercentLabel(double x, double y, double val, int color, double textSize = 0.028){
  if(val < 3.0) return;
  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextSize(textSize);
  lat.SetTextAlign(22); // centered
  lat.SetTextColor(color);
  lat.DrawLatex(x, y, Form("%.0f%%", val));
}

//============================================================================
// Main macro
//============================================================================
void PlotDoseFractionComparison_OO200(
    string a_doseFractionCsv = "DoseFractionSummary.csv",
    string a_glauberRoot     = "glauber_oo200_results_woodssaxon.root",
    string a_outputDir       = "."
){
  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outputDir.c_str(), kTRUE);

  vector<DoseFractionRow> rows = ReadDoseFractionCSV(a_doseFractionCsv);
  if(rows.empty()){
    cerr << "ERROR: no usable rows read from '" << a_doseFractionCsv << "' -- aborting."
         << " (This macro does NOT fabricate placeholder data -- see header. Run"
         << " ExtractDoseFractionSummary_OO200.C against a real fit_output.root first.)" << endl;
    return;
  }
  cout << "Read " << rows.size() << " centrality bins from '" << a_doseFractionCsv << "'." << endl;

  vector<GlauberSpecRow> glauber = ReadGlauberSpecFractions(a_glauberRoot);

  int n = (int) rows.size();

  TCanvas* c = new TCanvas("cDoseFractionComparison", "Dose Fraction Comparison", 950, 1050);
  c->Divide(1, 2);

  //--------------------------------------------------------------------------
  // Top panel: manually-stacked bar (fireball bottom, participant middle,
  // spectator top) -- see header for why this uses back-to-front overdraw
  // instead of THStack.
  //--------------------------------------------------------------------------
  c->cd(1);
  gPad->SetBottomMargin(0.22); // room for the legend drawn below the axis
  gPad->SetGridy();

  TH1D* hStackTop = new TH1D("hStackTop", "", n, 0, n);      // full height: fireball+participant+spectator
  TH1D* hStackMid = new TH1D("hStackMid", "", n, 0, n);      // fireball+participant
  TH1D* hStackBot = new TH1D("hStackBot", "", n, 0, n);      // fireball only

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
  hStackTop->GetYaxis()->SetTitle("Fraction of Net-Proton-Normalized Dose-Carrier Yield [%]");
  hStackTop->GetXaxis()->SetTitle("Centrality");
  hStackTop->SetTitle("O+O #sqrt{s_{NN}}=200 GeV -- Dose-Fraction Split by Source");

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

  TLegend* leg1 = new TLegend(0.12, 0.02, 0.90, 0.14);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetNColumns(1);
  leg1->AddEntry(hStackBot, "Fireball-produced (#pi, K, #bar{p})", "f");
  leg1->AddEntry(hStackMid, "Participant-associated (stopped net-p)", "f");
  leg1->AddEntry(hStackTop, "Spectator-associated (forward net-p)", "f");
  leg1->Draw();

  //--------------------------------------------------------------------------
  // Bottom panel: paired bars, data-driven vs geometric spectator-associated
  // fraction. Skipped gracefully (with an on-canvas message) if the Glauber
  // ROOT file wasn't found -- see header, this is the one soft-fail input.
  //--------------------------------------------------------------------------
  c->cd(2);

  if(!glauber.empty()){
    gPad->SetGridy();

    TH1D* hDataSpec = new TH1D("hDataSpec", "", n, 0, n);
    TH1D* hGeomSpec = new TH1D("hGeomSpec", "", n, 0, n);
    vector<bool> haveGeom(n, false);

    for(int i = 0; i < n; i++){
      hDataSpec->SetBinContent(i+1, 100.0 * rows[i].spectatorFrac);
      hDataSpec->GetXaxis()->SetBinLabel(i+1, rows[i].centLabel.c_str());
      double geomVal = 0.0;
      for(size_t j = 0; j < glauber.size(); j++){
        if(glauber[j].centLabel == rows[i].centLabel){
          geomVal = 100.0 * glauber[j].fracNspecOfTotal;
          haveGeom[i] = true;
          break;
        }
      }
      hGeomSpec->SetBinContent(i+1, geomVal);
    }

    double pairWidth = 0.35;
    hDataSpec->SetBarWidth(pairWidth); hDataSpec->SetBarOffset(0.5 - pairWidth);
    hGeomSpec->SetBarWidth(pairWidth); hGeomSpec->SetBarOffset(0.5);

    hDataSpec->SetFillColor(COLOR_SPECTATOR);
    hGeomSpec->SetFillColor(COLOR_GEOM_SPEC);

    hDataSpec->SetMaximum(100.0);
    hDataSpec->SetMinimum(0.0);
    hDataSpec->GetXaxis()->SetLabelSize(0.045);
    hDataSpec->GetYaxis()->SetTitle("Spectator-Associated Fraction [%]");
    hDataSpec->GetXaxis()->SetTitle("Centrality");
    hDataSpec->SetTitle("Spectator Fraction: Data-Driven (Net-Proton) vs. Geometric (Glauber)");

    hDataSpec->Draw("bar");
    hGeomSpec->Draw("bar same");

    for(int i = 0; i < n; i++){
      double x = hDataSpec->GetXaxis()->GetBinCenter(i+1);
      double dataVal = 100.0 * rows[i].spectatorFrac;
      TLatex latData;
      latData.SetTextFont(42); latData.SetTextSize(0.03); latData.SetTextAlign(22);
      latData.SetTextColor(COLOR_TEXT_ON_LIGHT);
      latData.DrawLatex(x - pairWidth/2.0 - 0.05, dataVal + 3.0, Form("%.0f%%", dataVal));
      if(haveGeom[i]){
        double geomVal = hGeomSpec->GetBinContent(i+1);
        TLatex latGeom;
        latGeom.SetTextFont(42); latGeom.SetTextSize(0.03); latGeom.SetTextAlign(22);
        latGeom.SetTextColor(COLOR_TEXT_ON_LIGHT);
        latGeom.DrawLatex(x + pairWidth/2.0 + 0.05, geomVal + 3.0, Form("%.0f%%", geomVal));
      }
    }

    TLegend* leg2 = new TLegend(0.12, 0.83, 0.88, 0.92);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(hDataSpec, "Data-driven (net-proton, spectator-associated)", "f");
    leg2->AddEntry(hGeomSpec, "Geometric (Woods-Saxon Glauber, N_{spec}/2A)", "f");
    leg2->Draw();
  } else {
    TLatex noGlauber;
    noGlauber.SetTextAlign(22);
    noGlauber.SetTextSize(0.04);
    noGlauber.DrawLatexNDC(0.5, 0.5, "Glauber comparison ROOT file not found -- panel skipped.");
  }

  string outPath = a_outputDir + "/OO200_DoseFractionComparison.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;

  // ---- Text table to stdout, same relief-for-contrast-WARN purpose as the Python version ----
  cout << endl;
  printf("%-10s%12s%15s%13s\n", "Cent", "Fireball %", "Participant %", "Spectator %");
  for(int i = 0; i < n; i++){
    printf("%-10s%12.1f%15.1f%13.1f\n", rows[i].centLabel.c_str(),
           100.0*rows[i].fireballFrac, 100.0*rows[i].participantFrac, 100.0*rows[i].spectatorFrac);
  }
}
