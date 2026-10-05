// PlotNeutronAttenuationLength.C
//
// Reads MarsNeutronAttenuation.root (TTree "NeutronAttenuation", written by
// MarsNeutronAttenuation, see main.cc) and plots the neutron attenuation
// length (mean depth to first interaction) vs. incident neutron energy in
// Mars regolith, log-log, with statistical error bars (standard error of
// the mean from RunAction.cc).
//
// Also re-checks and reports the transmittedFraction column at each energy
// (the "was the slab thick enough" diagnostic -- see DetectorParameters.hh
// and RunAction.cc) so a bad point isn't silently plotted as if it were
// trustworthy.
//
// STATUS: not run -- depends on MarsNeutronAttenuation's .root output.
// Paren/brace/bracket balance checked with a script before delivery, same
// as this session's other ROOT macros; logical correctness against a real
// file has not been exercised.
//
// Usage:
//   root -l -q 'PlotNeutronAttenuationLength.C("MarsNeutronAttenuation.root")'

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TColor.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

void PlotNeutronAttenuationLength(string a_rootFile = "MarsNeutronAttenuation.root",
                                   string a_outputDir = "."){
  TFile* f = TFile::Open(a_rootFile.c_str(), "READ");
  if(!f || f->IsZombie()){
    cout << "ERROR: could not open " << a_rootFile << endl;
    return;
  }
  TTree* tree = (TTree*) f->Get("NeutronAttenuation");
  if(!tree){
    cout << "ERROR: TTree \"NeutronAttenuation\" not found in " << a_rootFile << endl;
    f->Close();
    return;
  }

  Double_t energyMeV_in, meanDepth_in, stdErr_in, transmittedFrac_in;
  Int_t nEvents_in, nInteracted_in, nTransmitted_in;
  tree->SetBranchAddress("energyMeV", &energyMeV_in);
  tree->SetBranchAddress("nEvents", &nEvents_in);
  tree->SetBranchAddress("nInteracted", &nInteracted_in);
  tree->SetBranchAddress("nTransmitted", &nTransmitted_in);
  tree->SetBranchAddress("meanFirstInteractionDepth_cm", &meanDepth_in);
  tree->SetBranchAddress("stdErrOfMean_cm", &stdErr_in);
  tree->SetBranchAddress("transmittedFraction", &transmittedFrac_in);

  vector<double> energy, depth, depthErr;
  int nFlagged = 0;
  Long64_t nEntries = tree->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    tree->GetEntry(i);
    energy.push_back(energyMeV_in);
    depth.push_back(meanDepth_in);
    depthErr.push_back(stdErr_in);
    if(transmittedFrac_in > 0.02){
      nFlagged++;
      cout << "WARNING: E=" << energyMeV_in << " MeV had "
           << (transmittedFrac_in * 100.0) << "% of primaries transit the "
           << "whole slab without interacting -- this point's depth is "
           << "likely a truncated underestimate (see DetectorParameters.hh)." << endl;
    }
  }
  f->Close();

  if(energy.empty()){
    cout << "ERROR: no entries read from " << a_rootFile << " -- nothing to plot." << endl;
    return;
  }
  if(nFlagged > 0){
    cout << nFlagged << " of " << energy.size() << " energy points were flagged above -- "
         << "consider re-running with a thicker slab (MarsNeutron::kSlabArealThickness) "
         << "if this matters for the energies you care about." << endl;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const char* SURFACE = "#fcfcfb";
  const char* TEXT_PRIMARY = "#0b0b0b";
  const char* TEXT_SECONDARY = "#52514e";
  const char* TEXT_MUTED = "#898781";
  const char* BASELINE = "#c3c2b7";
  const char* SERIES = "#2a78d6"; // dataviz palette slot 1, blue

  TCanvas* c = new TCanvas("c_NeutronAttenuation", "Neutron attenuation length vs energy", 950, 750);
  c->SetLogx();
  c->SetLogy();
  c->SetFillColor(TColor::GetColor(SURFACE));
  c->SetFrameFillColor(TColor::GetColor(SURFACE));
  c->SetLeftMargin(0.13); c->SetRightMargin(0.06);
  c->SetBottomMargin(0.14); c->SetTopMargin(0.12); // only the title lives up here now, no subtitle

  double xMin = energy.front() * 0.7, xMax = energy.back() * 1.4;
  double yMin = 1e9, yMax = -1e9;
  for(size_t i = 0; i < depth.size(); i++){
    yMin = min(yMin, depth[i] - depthErr[i]);
    yMax = max(yMax, depth[i] + depthErr[i]);
  }
  yMin = max(yMin * 0.7, 0.1);
  yMax = yMax * 1.4;

  TH1F* frame = c->DrawFrame(xMin, yMin, xMax, yMax);
  frame->GetXaxis()->SetTitle("Neutron Kinetic Energy (MeV)");
  frame->GetYaxis()->SetTitle("Mean Depth to First Interaction (cm) -- Attenuation Length");
  frame->GetXaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetYaxis()->SetTitleColor(TColor::GetColor(TEXT_SECONDARY));
  frame->GetXaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetYaxis()->SetLabelColor(TColor::GetColor(TEXT_MUTED));
  frame->GetXaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetYaxis()->SetAxisColor(TColor::GetColor(BASELINE));
  frame->GetYaxis()->SetTitleSize(0.032);

  TGraphErrors* g = new TGraphErrors((int) energy.size(), &energy[0], &depth[0], nullptr, &depthErr[0]);
  g->SetLineColor(TColor::GetColor(SERIES));
  g->SetMarkerColor(TColor::GetColor(SERIES));
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.9);
  g->SetLineWidth(2);
  g->Draw("LP SAME");

  TLatex title; title.SetNDC(); title.SetTextColor(TColor::GetColor(TEXT_PRIMARY));
  title.SetTextFont(62); title.SetTextSize(0.035); title.SetTextAlign(21);
  title.DrawLatex(0.5, 0.955, "Neutron Attenuation Length in Mars Regolith vs. Energy");

  // No in-plot subtitle -- the method/caveat detail that used to live here
  // belongs in the slide's own figure caption instead.

  c->RedrawAxis();
  string out = a_outputDir + "/NeutronAttenuationLengthVsEnergy.png";
  c->SaveAs(out.c_str());
  cout << "Wrote " << out << endl;
}
