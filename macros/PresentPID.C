// PresentPID.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: presentation-ready PID plots from a PicoBinner yield file --
// dE/dx vs momentum (TPC) and mass vs momentum (BTOF), the two classic "does the
// detector separate species" plots. These are the histograms most people expect to
// see first when you say "here's the PID."
//
// Reads (exact names from PicoBinner.cxx -- see FindHistoInFile.C if your file uses
// different ones, e.g. because a_makeBasicHistos was false):
//   dEdxVsMom_Plus, dEdxVsMom_Minus     -- TPC dE/dx vs p_tot, log-binned in p
//   massVsMom_bTOF_Plus, massVsMom_bTOF_Minus -- BTOF mass vs p_tot
//
// Usage:
//   root -l -q 'PresentPID.C("/path/to/yieldHistos_OOGeV_proton_2026_07_01.root")'
//   root -l -q 'PresentPID.C("/path/to/file.root","/path/to/outDir")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void PresentPID(string a_inputFile, string a_outDir = "."){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* dEdxPlus   = (TH2*) inFile->Get("dEdxVsMom_Plus");
  TH2* dEdxMinus  = (TH2*) inFile->Get("dEdxVsMom_Minus");
  TH2* massPlus   = (TH2*) inFile->Get("massVsMom_bTOF_Plus");
  TH2* massMinus  = (TH2*) inFile->Get("massVsMom_bTOF_Minus");

  if(!dEdxPlus && !massPlus){
    cout << "!!! Neither dEdxVsMom_Plus nor massVsMom_bTOF_Plus found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced)." << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("PresentPID","PID overview",1400,1000);
  c->Divide(2,2);

  c->cd(1); gPad->SetLogx(); gPad->SetLogz();
  if(dEdxPlus){ dEdxPlus->SetTitle("dE/dx vs p (positive tracks)"); dEdxPlus->Draw("COLZ"); }

  c->cd(2); gPad->SetLogx(); gPad->SetLogz();
  if(dEdxMinus){ dEdxMinus->SetTitle("dE/dx vs p (negative tracks)"); dEdxMinus->Draw("COLZ"); }

  c->cd(3); gPad->SetLogz();
  if(massPlus){ massPlus->SetTitle("BTOF mass vs p (positive tracks)"); massPlus->Draw("COLZ"); }

  c->cd(4); gPad->SetLogz();
  if(massMinus){ massMinus->SetTitle("BTOF mass vs p (negative tracks)"); massMinus->Draw("COLZ"); }

  string outPath = a_outDir + "/PresentPID.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
