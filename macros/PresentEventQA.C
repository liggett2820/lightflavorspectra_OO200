// PresentEventQA.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready event-level sanity plots from a PicoBinner yield
// file -- vertex distributions, multiplicity, and events-per-centrality-bin. These are
// the "the event sample looks reasonable" plots that should come before anyone trusts
// the physics downstream.
//
// Reads (exact names from PicoBinner.cxx):
//   zVertex   -- event Z vertex position
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

  TH1* zVertex   = (TH1*) inFile->Get("zVertex");
  TH2* xyHisto   = (TH2*) inFile->Get("xyHisto");
  TH1* refMult   = (TH1*) inFile->Get("refMult");
  TH1* centEvents = (TH1*) inFile->Get(("centEvents"+a_speciesSuffix).c_str());
  if(!centEvents) centEvents = (TH1*) inFile->Get("centEvents"); // try the no-suffix name too

  if(!zVertex && !xyHisto && !refMult && !centEvents){
    cout << "!!! None of zVertex/xyHisto/refMult/centEvents(" << a_speciesSuffix
         << ") found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file." << endl;
    return;
  }

  gStyle->SetOptStat(1);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentEventQA","Event QA overview",1400,1000);
  c->Divide(2,2);

  // zVertex's X axis already carries units (cm) from PicoBinner.cxx; its Y axis (a
  // plain bin count) isn't set there, so label it here for consistency with
  // xyHisto/refMult/centEvents, which already have it baked into their titles.
  c->cd(1);
  if(zVertex){
    zVertex->SetTitle("Event Z vertex");
    zVertex->GetYaxis()->SetTitle("Number of Events");
    // Display range capped at [-4,4]cm -- just clips the plotted range (SetCutClass.C's
    // setZRange(-2,2) event-level cut is well inside this), doesn't rebin or drop any
    // underlying entries.
    zVertex->GetXaxis()->SetRangeUser(-4, 4);
    zVertex->Draw();
  }

  // Extra right margin -- otherwise the z-axis (palette) power-of-10 labels get
  // clipped/overlap the next pad in this Divide(2,2) grid.
  c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
  if(xyHisto){
    xyHisto->SetTitle("Event vertex X vs Y");
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

  // Display range capped at 140 -- refMult is booked 0-1000 in PicoBinner.cxx; this
  // just clips the plotted range, doesn't rebin or drop any underlying entries.
  c->cd(3); gPad->SetLogy();
  if(refMult){
    refMult->SetTitle("refMult");
    refMult->GetXaxis()->SetRangeUser(0, 140);
    refMult->Draw();

    // Overlay the official 5-bin centrality cuts from SetCutClass.C (2026-07-03:
    // {44,37,28,17,5} at the 5/10/20/40/80% boundaries -- these are the values that
    // actually produced this file's centIndex, NOT re-derived here).
    double centCuts[5]     = {44, 37, 28, 17, 5};
    int    centPercents[5] = {5, 10, 20, 40, 80};
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
    centEvents->SetTitle("Events per centrality bin");
    // Force the y-axis to start at 0 -- ROOT's default autoscale starts near the
    // lowest bin content, which zero-suppresses the axis and exaggerates bin-to-bin
    // differences that are actually small relative to the total.
    centEvents->SetMinimum(0);
    centEvents->Draw();
  }

  string outPath = a_outDir + "/PresentEventQA.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
