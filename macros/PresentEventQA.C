// PresentEventQA.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready event-level sanity plots from a PicoBinner yield
// file -- vertex distributions, multiplicity, and events-per-centrality-bin. These are
// the "the event sample looks reasonable" plots that should come before anyone trusts
// the physics downstream.
//
// Reads (exact names from PicoBinner.cxx):
//   dca       -- per-track DCA (distance of closest approach) to the primary vertex
//   xyHisto   -- event vertex X vs Y
//   refMult   -- raw reference multiplicity per event
//   centEvents[<Species>] -- events per centrality bin. If PicoBinner was run for a
//                 single species (a_partIndex >= 0, the convention this repo's SDCC
//                 jobs use -- see xml/runPicoBinner_OO200_SDCC_template.xml), this
//                 histogram is named "centEvents<Species>" (e.g. "centEventsProton"),
//                 not plain "centEvents" -- pass a_speciesSuffix accordingly.
//
// Usage:
//   root -l -q 'PresentEventQA.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root","Proton")'
//   root -l -q 'PresentEventQA.C("/path/to/file.root","")'    // if centEvents has no species suffix
//   root -l -q 'PresentEventQA.C("/path/to/file.root","Proton","/path/to/outDir")'

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include <iostream>
using namespace std;

void PresentEventQA(string a_inputFile, string a_speciesSuffix = "Proton", string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH1* dca       = (TH1*) inFile->Get("dca");
  TH2* xyHisto   = (TH2*) inFile->Get("xyHisto");
  TH1* refMult   = (TH1*) inFile->Get("refMult");
  // Oct 2026 productions: corrected EPD totnMIP (the centrality variable since then), name has the species suffix
  TH1* totnMIP   = (TH1*) inFile->Get(("totnMIPCorr"+a_speciesSuffix).c_str());
  TH1* centEvents = (TH1*) inFile->Get(("centEvents"+a_speciesSuffix).c_str());
  if(!centEvents) centEvents = (TH1*) inFile->Get("centEvents"); // try the no-suffix name too

  if(!dca && !xyHisto && !refMult && !centEvents){
    cout << "!!! None of dca/xyHisto/refMult/centEvents(" << a_speciesSuffix
         << ") found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file." << endl;
    return;
  }

  gStyle->SetOptStat(1);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentEventQA","Event QA Overview",1400,1000);
  c->Divide(2,2);

  // dca is track-level (booked 0-3.5, unitless in PicoBinner.cxx, but the DCA cut
  // itself is quoted in cm elsewhere in this repo -- e.g. SetCutClass.C's "DCA
  // tightened to 1cm" -- so label the axis accordingly here.
  c->cd(1);
  if(dca){
    dca->SetTitle("Track DCA");
    dca->GetXaxis()->SetTitle("DCA [cm]");
    // Display range capped at 3cm -- just clips the plotted range (booked out to 3.5),
    // doesn't rebin or drop any underlying entries.
    dca->GetXaxis()->SetRangeUser(0, 3);
    dca->Draw();
  }

  // Extra right margin -- otherwise the z-axis (palette) power-of-10 labels get
  // clipped/overlap the next pad in this Divide(2,2) grid.
  c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(xyHisto){
    xyHisto->SetTitle("Event Vertex X vs Y");
    xyHisto->Draw("COLZ");
    // The default stats box (gStyle->SetOptStat(1) above) sits top-right, right on top
    // of the COLZ palette -- move it to the top-left corner instead, out of the way.
    gPad->Update();
    TPaveStats* xyStats = (TPaveStats*) xyHisto->FindObject("stats");
    if(xyStats){
      xyStats->SetX1NDC(0.15); xyStats->SetX2NDC(0.40);
      xyStats->SetY1NDC(0.75); xyStats->SetY2NDC(0.90);
      gPad->Modified();
    }
  }

  // Centrality variable with the class edges. Files with totnMIPCorr<Species> (Oct 2026 on) use the EPD
  // classes of SetCutClass.C, totnMIP_corr >= {364,310,234,130,36} at 5/10/20/40/80% (the 80-100% edge, 0,
  // sits on the axis). Older files: raw refMult with {44,37,28,17,5}. Note totnMIPCorr is filled for every
  // event passing the vertex cuts, before the pile-up and centrality selections.
  c->cd(3); gPad->SetLogy();
  TH1* centVar = totnMIP ? totnMIP : refMult;
  if(centVar){
    double centCuts[5]     = {44, 37, 28, 17, 5};
    if(totnMIP){ centCuts[0] = 364; centCuts[1] = 310; centCuts[2] = 234; centCuts[3] = 130; centCuts[4] = 36; }
    int    centPercents[5] = {5, 10, 20, 40, 80};
    centVar->SetTitle(totnMIP ? "Corrected EPD totnMIP" : "refMult");
    centVar->GetXaxis()->SetRangeUser(0, totnMIP ? 650 : 140);
    centVar->Draw();
    gPad->Update();
    double yMin = TMath::Power(10, gPad->GetUymin());
    double yMax = TMath::Power(10, gPad->GetUymax());
    for(int iii = 0; iii < 5; iii++){
      TLine* cutLine = new TLine(centCuts[iii], yMin, centCuts[iii], yMax);
      cutLine->SetLineColor(kRed);
      cutLine->SetLineStyle(2);
      cutLine->SetLineWidth(2);
      cutLine->Draw("SAME");

      TLatex* cutLabel = new TLatex(centCuts[iii], yMax, Form("%d%%",centPercents[iii]));
      cutLabel->SetTextColor(kRed);
      cutLabel->SetTextSize(0.035);
      cutLabel->SetTextAlign(21); // centered horizontally, bottom-aligned
      cutLabel->Draw("SAME");
    }
  }

  c->cd(4);
  if(centEvents){
    centEvents->SetTitle("Events per Centrality Bin");
    // Force the y-axis to start at 0 -- ROOT's default autoscale starts near the
    // lowest bin content, which zero-suppresses the axis and exaggerates bin-to-bin
    // differences that are actually small relative to the total.
    centEvents->SetMinimum(0);
    if(centEvents->GetNbinsX() == 6){
      const char* lab[6] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
      for(int iii = 0; iii < 6; iii++) centEvents->GetXaxis()->SetBinLabel(iii+1, lab[iii]);
      centEvents->GetXaxis()->SetTitle("");
    }
    centEvents->Draw();
  }

  string outPath = a_outDir + "/PresentEventQA.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
