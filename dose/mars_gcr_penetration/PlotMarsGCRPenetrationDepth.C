// PlotMarsGCRPenetrationDepth.C
//
// Reads MarsGCRPenetrationDepth.root (TTree "MarsGCRPenetrationDepth", written
// by MarsGCRPenetrationDepth.cpp) and draws the same log-log, multi-species
// penetration-depth-vs-energy chart that plot_mars_gcr_depth.py previously
// produced -- ROOT/C++ replacement per Andrew's request to use C++ instead of
// Python for future deliverables in this project.
//
// STATUS: like every ROOT macro in this session, this has NOT been compiled
// or run by me -- this sandbox has no ROOT installed. The physics/data it
// plots comes from MarsGCRPenetrationDepth.root, itself produced by code
// whose CSV-output predecessor WAS compiled and numerically verified (see
// MarsGCRPenetrationDepth.cpp header) but whose ROOT-output stage was not.
//
// Usage:
//   root -l -q 'PlotMarsGCRPenetrationDepth.C("MarsGCRPenetrationDepth.root")'
//
// Output: MarsGCRPenetrationDepth.png in the current directory.

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLine.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

//============================================================================
// dataviz skill default categorical palette, slots 1-8, fixed order --
// same colors used in the earlier matplotlib version, for continuity.
//============================================================================
struct SpeciesStyle { string name; string hex; };
vector<SpeciesStyle> SpeciesOrder(){
  return {
    {"H",  "#2a78d6"}, // slot 1 blue
    {"He", "#eb6834"}, // slot 2 orange
    {"C",  "#1baf7a"}, // slot 3 aqua
    {"O",  "#eda100"}, // slot 4 yellow
    {"Ne", "#e87ba4"}, // slot 5 magenta
    {"Mg", "#008300"}, // slot 6 green
    {"Si", "#4a3aa7"}, // slot 7 violet
    {"Fe", "#e34948"}, // slot 8 red
  };
}

const char* TEXT_PRIMARY   = "#0b0b0b";
const char* TEXT_SECONDARY = "#52514e";
const char* TEXT_MUTED     = "#898781";
const char* GRIDLINE       = "#e1e0d9";
const char* BASELINE       = "#c3c2b7";
const char* SURFACE        = "#fcfcfb";

//============================================================================
// ReadSpeciesSeries: pull (KE_MeV_per_nucleon, depth_cm_nominal_rho1p52) pairs
// for one species out of the TTree, in energy order (the tree is already
// written species-major / energy-ascending by MarsGCRPenetrationDepth.cpp,
// but we don't rely on that -- we filter by species name explicitly).
//============================================================================
bool ReadSpeciesSeries(TTree* a_tree, const string& a_species,
                        vector<double>& a_E, vector<double>& a_depth){
  Char_t species_in[8];
  Double_t KE_MeV_per_nucleon_in, depth_cm_nominal_rho1p52_in;
  a_tree->SetBranchAddress("species", species_in);
  a_tree->SetBranchAddress("KE_MeV_per_nucleon", &KE_MeV_per_nucleon_in);
  a_tree->SetBranchAddress("depth_cm_nominal_rho1p52", &depth_cm_nominal_rho1p52_in);

  a_E.clear();
  a_depth.clear();
  Long64_t nEntries = a_tree->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    a_tree->GetEntry(i);
    if(a_species != string(species_in)) continue;
    a_E.push_back(KE_MeV_per_nucleon_in);
    a_depth.push_back(depth_cm_nominal_rho1p52_in);
  }
  return !a_E.empty();
}

void PlotMarsGCRPenetrationDepth(string a_rootFile = "MarsGCRPenetrationDepth.root",
                                  string a_outputDir = "."){
  TFile* f = TFile::Open(a_rootFile.c_str(), "READ");
  if(!f || f->IsZombie()){
    cout << "ERROR: could not open " << a_rootFile << endl;
    return;
  }
  TTree* tree = (TTree*) f->Get("MarsGCRPenetrationDepth");
  if(!tree){
    cout << "ERROR: TTree \"MarsGCRPenetrationDepth\" not found in " << a_rootFile << endl;
    f->Close();
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* c = new TCanvas("c_MarsGCR", "Mars GCR Penetration Depth", 950, 750);
  c->SetLogx();
  c->SetLogy();
  c->SetFillColor(TColor::GetColor(SURFACE));
  c->SetFrameFillColor(TColor::GetColor(SURFACE));
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.14); // extra room for end-of-line species labels
  c->SetBottomMargin(0.16);
  c->SetTopMargin(0.12); // only the title lives up here now, no subtitle

  vector<SpeciesStyle> speciesOrder = SpeciesOrder();

  // First pass: read every series and track the global depth range, so the
  // frame can be sized before any graph is drawn.
  map<string, vector<double>> E_by_species, depth_by_species;
  double globalMinDepth = 1e300, globalMaxDepth = -1e300;
  double globalMinE = 1e300, globalMaxE = -1e300;
  for(const SpeciesStyle& sp : speciesOrder){
    vector<double> E, depth;
    if(!ReadSpeciesSeries(tree, sp.name, E, depth)){
      cout << "WARNING: no entries found for species \"" << sp.name << "\" -- skipping." << endl;
      continue;
    }
    E_by_species[sp.name] = E;
    depth_by_species[sp.name] = depth;
    for(double e : E){ globalMinE = min(globalMinE, e); globalMaxE = max(globalMaxE, e); }
    for(double d : depth){ globalMinDepth = min(globalMinDepth, d); globalMaxDepth = max(globalMaxDepth, d); }
  }
  if(E_by_species.empty()){
    cout << "ERROR: no species series read from tree -- nothing to plot." << endl;
    f->Close();
    return;
  }

  // Frame with extra headroom for end-of-line labels and the right-margin
  // extension used for the label text itself (x range extended ~3x beyond
  // the data so labels drawn just past the last point aren't clipped).
  double xMax = globalMaxE * 3.0;
  TH1F* frame = c->DrawFrame(globalMinE * 0.9, globalMinDepth * 0.5, xMax, globalMaxDepth * 2.0);
  frame->GetXaxis()->SetTitle("Kinetic Energy (MeV / Nucleon)");
  frame->GetYaxis()->SetTitle("CSDA Penetration Depth Into Mars Regolith (cm)");
  frame->GetXaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetYaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetXaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetYaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetXaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetYaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetXaxis()->SetTitleSize(0.038);
  frame->GetYaxis()->SetTitleSize(0.038);

  gPad->SetGridx(false);
  gPad->SetGridy(false);

  TLegend* leg = new TLegend(0.14, 0.68, 0.40, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  leg->SetTextColor(TColor::GetColor(TEXT_SECONDARY));
  leg->SetNColumns(2);
  leg->SetHeader("GCR reference ion");

  vector<TGraph*> keepAlive; // ROOT graphs must outlive the Draw calls

  for(const SpeciesStyle& sp : speciesOrder){
    if(E_by_species.find(sp.name) == E_by_species.end()) continue;
    const vector<double>& E = E_by_species[sp.name];
    const vector<double>& depth = depth_by_species[sp.name];

    TGraph* g = new TGraph((int) E.size(), &E[0], &depth[0]);
    keepAlive.push_back(g);
    int col = TColor::GetColor(sp.hex.c_str());
    g->SetLineColor(col);
    g->SetLineWidth(3);
    g->Draw("L SAME");
    leg->AddEntry(g, sp.name.c_str(), "l");

    // direct end-of-line label, same rationale as the matplotlib version:
    // relief for the palette's contrast-WARN slots (aqua/yellow/magenta),
    // and generally easier to read than a legend alone on an 8-line log-log chart.
    TLatex* lab = new TLatex(E.back() * 1.15, depth.back(), sp.name.c_str());
    lab->SetTextColor(col);
    lab->SetTextFont(62); // bold
    lab->SetTextSize(0.03);
    lab->SetTextAlign(12); // left, vcenter
    lab->Draw();
  }

  leg->Draw();

  TLatex title;
  title.SetNDC();
  title.SetTextColor(TColor::GetColor(TEXT_PRIMARY));
  title.SetTextFont(62);
  title.SetTextSize(0.035);
  title.SetTextAlign(21); // center, bottom
  title.DrawLatex(0.5, 0.955, "GCR Reference-Ion Penetration Depth vs. Energy");

  // No in-plot subtitle -- the method detail that used to live here belongs
  // in the slide's own figure caption instead. The bottom CAVEAT text below
  // is left in place -- it's a physics-validity warning about the plotted
  // curves themselves (electromagnetic-only, no nuclear fragmentation), not
  // a methods subtitle, so it doesn't fit the same "caption, not caveat on
  // the image" reasoning. Flag if you'd like that removed too.

  TLatex caveat;
  caveat.SetNDC();
  caveat.SetTextColor(TColor::GetColor(TEXT_MUTED));
  caveat.SetTextFont(42);
  caveat.SetTextSize(0.020);
  caveat.SetTextAlign(21);
  caveat.DrawLatex(0.5, 0.02,
    "CAVEAT: electromagnetic stopping only -- omits nuclear fragmentation, which materially shortens the");
  caveat.DrawLatex(0.5, 0.005,
    "effective range of heavier ions (Ne and up) at GCR energies. Reads as an upper bound, not a full transport result.");

  c->RedrawAxis();

  string outPath = a_outputDir + "/MarsGCRPenetrationDepth.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;

  f->Close();
}
