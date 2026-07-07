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

  c->cd(1);
  if(zVertex){ zVertex->SetTitle("Event Z vertex"); zVertex->Draw(); }

  c->cd(2); gPad->SetLogz();
  if(xyHisto){ xyHisto->SetTitle("Event vertex X vs Y"); xyHisto->Draw("COLZ"); }

  c->cd(3); gPad->SetLogy();
  if(refMult){ refMult->SetTitle("refMult"); refMult->Draw(); }

  c->cd(4);
  if(centEvents){ centEvents->SetTitle("Events per centrality bin"); centEvents->Draw(); }

  string outPath = a_outDir + "/PresentEventQA.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
