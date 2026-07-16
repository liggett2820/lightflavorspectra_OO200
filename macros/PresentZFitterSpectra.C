// PresentZFitterSpectra.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: the original request (before the pion TPC crash/Bichsel-curve
// detour took over) was a macro to plot the pion spectra straight out of ZFitter's
// output -- raw yield vs mT-m0, overlaying all rapidity bins, combining TPC and BTOF
// on the same plot. This is that macro, generalized to any of the three species
// (pion/kaon/proton) this analysis fits, since the plotting logic doesn't depend on
// species beyond the name string.
//
// WHAT IT READS (names/paths confirmed directly from ZFitter::makeSpectra() in
// source/ZFitter.cxx, not guessed):
//   TPC:  DeDx_FitData/<Species>/dEdxSpectra_<Species>_Plus_Cent%02d
//         DeDx_FitData/<Species>/dEdxSpectra_<Species>_Minus_Cent%02d
//   BTOF: BTOF_FitData/<Species>/InvBetaBTOFSpectra_<Species>_Plus_Cent%02d
//         BTOF_FitData/<Species>/InvBetaBTOFSpectra_<Species>_Minus_Cent%02d
//   ETOF: ETOF_FitData/<Species>/InvBetaETOFSpectra_<Species>_Plus_Cent%02d
//         ETOF_FitData/<Species>/InvBetaETOFSpectra_<Species>_Minus_Cent%02d
//         (only drawn if present -- this analysis has m_hasETOF=true, but no ETOF
//         ZFitter run has necessarily been done yet, so this is optional at draw time)
// All three are TH2D with X-axis = rapidity bin (title "y_{...}"), Y-axis = mT-m0
// [GeV/c^2] (title "m_{T}-m_{...}"), bin content = fitted yield/amplitude -- confirmed
// from the SetBinContent(rapBin, mTm0Bin, value) call and SetTitle() strings at the
// point these histograms are created in ZFitter.cxx. "Overlaying all rapidity bins"
// means, for a fixed species/charge/centrality: take a Y-axis projection (ProjectionY)
// at each non-empty X bin (rapidity bin) to get a 1D yield-vs-(mT-m0) spectrum, then
// draw all of those on one pad.
//
// NAMING CAVEAT (found and fixed 2026-07-15, see source/RawSpectraModifier.cxx): the
// Plus/Minus underscore in the names above matches what ZFitter.cxx actually writes.
// RawSpectraModifier.cxx used to look for these without the underscore (e.g.
// "dEdxSpectra_PionPlus_Cent00") and silently got NULL back -- if you ever see NULL
// histograms reading a ZFitter output file with a name that looks *almost* right,
// double check the underscore placement against ZFitter.cxx's own SetName() calls,
// not against a reader macro's assumptions.
//
// LEGEND CHOICE: a fit typically has dozens of rapidity bins (the pion TPC crash log
// referenced RapBin 35, i.e. at least 36 bins) -- one legend entry per bin would be
// unreadable. Instead, each rapidity bin's curve is colored along a continuous
// rainbow gradient (dark = most negative/backward rapidity in this bin range, bright
// = most positive/forward), with a compact color-key drawn in the corner instead of a
// per-bin legend, plus a two-entry marker-style legend (filled circle = TPC, open
// circle = BTOF, open square = ETOF) since marker style (not color) distinguishes
// detector.
//
// Usage:
//   root -l -q 'PresentZFitterSpectra.C("/path/to/zfitter_output_file.root",0,0)'
//     -> pion (partIndex 0), centrality bin 0, both charges, output to "."
//   root -l -q 'PresentZFitterSpectra.C("/path/to/file.root",1,2,"/path/to/outDir")'
//     -> kaon (partIndex 1), centrality bin 2, output to /path/to/outDir
//
// a_partIndex: 0=Pion, 1=Kaon, 2=Proton (matches ParticleInfo's speciesIndex ordering)
// a_centIndex: centrality bin, matches ZFitter's nCentToRun indexing (0 = most central
//              in this analysis's 6-bin scheme, per README's centrality section)

#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPad.h"
#include <iostream>
#include <vector>
using namespace std;

//_______________________________________________________________________________
// One pad's worth of work: given the Plus (or Minus) TPC/BTOF/ETOF TH2Ds for one
// species+centrality, project out every non-empty rapidity bin and draw them all
// overlaid, TPC/BTOF/ETOF distinguished by marker style, rapidity bin distinguished
// by color along a rainbow gradient.
void drawOneChargeSign(TH2* a_tpcHist, TH2* a_btofHist, TH2* a_etofHist,
                        string a_speciesName, string a_chargeLabel, int a_centIndex){

  gPad->SetLogy();
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);

  if(!a_tpcHist && !a_btofHist && !a_etofHist){
    TLatex missingLabel;
    missingLabel.SetTextSize(0.05);
    missingLabel.DrawLatexNDC(0.15,0.5,Form("No TPC/BTOF/ETOF spectra found for %s %s, Cent%02d",
                                             a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex));
    return;
  }

  // Use whichever histogram is available to set up the frame/axis ranges and to
  // determine how many rapidity bins exist (TPC, BTOF, and ETOF share the same
  // rapidity binning by construction -- all cloned from the same structHisto-family
  // template in ZFitter.cxx -- so any of them is a valid reference).
  TH2* referenceHist = a_tpcHist ? a_tpcHist : (a_btofHist ? a_btofHist : a_etofHist);
  int nRapBins = referenceHist->GetNbinsX();

  bool drewFrame = false;
  double globalMaxY = 0;
  double globalMinPositiveY = 1e30;

  // First pass: figure out the y-axis range across every projection we're about to
  // draw, so the frame can be sized once up front (ROOT's default auto-ranging on
  // the first Draw() call doesn't know about the later ones).
  vector<TH1D*> tpcProjections, btofProjections, etofProjections;
  for(int rapBin = 1; rapBin <= nRapBins; rapBin++){
    TH1D* tpcProj  = a_tpcHist  ? a_tpcHist->ProjectionY(Form("tpcProj_%s_%s_Cent%02d_rap%03d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,rapBin),rapBin,rapBin) : NULL;
    TH1D* btofProj = a_btofHist ? a_btofHist->ProjectionY(Form("btofProj_%s_%s_Cent%02d_rap%03d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,rapBin),rapBin,rapBin) : NULL;
    TH1D* etofProj = a_etofHist ? a_etofHist->ProjectionY(Form("etofProj_%s_%s_Cent%02d_rap%03d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex,rapBin),rapBin,rapBin) : NULL;

    if(tpcProj  && tpcProj->Integral()  <= 0){ delete tpcProj;  tpcProj  = NULL; }
    if(btofProj && btofProj->Integral() <= 0){ delete btofProj; btofProj = NULL; }
    if(etofProj && etofProj->Integral() <= 0){ delete etofProj; etofProj = NULL; }

    tpcProjections.push_back(tpcProj);
    btofProjections.push_back(btofProj);
    etofProjections.push_back(etofProj);

    for(TH1D* proj : {tpcProj,btofProj,etofProj}){
      if(!proj) continue;
      for(int b = 1; b <= proj->GetNbinsX(); b++){
        double content = proj->GetBinContent(b);
        if(content <= 0) continue;
        if(content > globalMaxY) globalMaxY = content;
        if(content < globalMinPositiveY) globalMinPositiveY = content;
      }
    }
  }

  if(globalMaxY <= 0){
    TLatex emptyLabel;
    emptyLabel.SetTextSize(0.05);
    emptyLabel.DrawLatexNDC(0.15,0.5,Form("All rapidity bins empty for %s %s, Cent%02d",
                                          a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex));
    return;
  }

  int nColors = 200;
  gStyle->SetPalette(kRainBow);

  for(int rapBin = 1; rapBin <= nRapBins; rapBin++){
    double colorFrac = (nRapBins > 1) ? (double)(rapBin-1)/(double)(nRapBins-1) : 0.5;
    int colorIndex = TColor::GetColorPalette((int)(colorFrac*(nColors-1)));

    TH1D* tpcProj  = tpcProjections[rapBin-1];
    TH1D* btofProj = btofProjections[rapBin-1];
    TH1D* etofProj = etofProjections[rapBin-1];

    if(tpcProj){
      tpcProj->SetMarkerStyle(20); // filled circle = TPC
      tpcProj->SetMarkerColor(colorIndex);
      tpcProj->SetLineColor(colorIndex);
      tpcProj->SetMarkerSize(0.7);
      tpcProj->GetYaxis()->SetRangeUser(globalMinPositiveY*0.3,globalMaxY*3);
      tpcProj->SetTitle(Form("%s %s, Cent%02d -- raw ZFitter yield vs m_{T}-m_{0}, all rapidity bins;m_{T}-m_{0} [GeV/c^{2}];Yield",
                              a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex));
      tpcProj->Draw(drewFrame ? "P SAME" : "P");
      drewFrame = true;
    }
    if(btofProj){
      btofProj->SetMarkerStyle(24); // open circle = BTOF
      btofProj->SetMarkerColor(colorIndex);
      btofProj->SetLineColor(colorIndex);
      btofProj->SetMarkerSize(0.7);
      if(!drewFrame){
        btofProj->GetYaxis()->SetRangeUser(globalMinPositiveY*0.3,globalMaxY*3);
        btofProj->SetTitle(Form("%s %s, Cent%02d -- raw ZFitter yield vs m_{T}-m_{0}, all rapidity bins;m_{T}-m_{0} [GeV/c^{2}];Yield",
                                 a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex));
      }
      btofProj->Draw(drewFrame ? "P SAME" : "P");
      drewFrame = true;
    }
    if(etofProj){
      etofProj->SetMarkerStyle(25); // open square = ETOF
      etofProj->SetMarkerColor(colorIndex);
      etofProj->SetLineColor(colorIndex);
      etofProj->SetMarkerSize(0.7);
      if(!drewFrame){
        etofProj->GetYaxis()->SetRangeUser(globalMinPositiveY*0.3,globalMaxY*3);
        etofProj->SetTitle(Form("%s %s, Cent%02d -- raw ZFitter yield vs m_{T}-m_{0}, all rapidity bins;m_{T}-m_{0} [GeV/c^{2}];Yield",
                                 a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex));
      }
      etofProj->Draw(drewFrame ? "P SAME" : "P");
      drewFrame = true;
    }
  }

  // Marker-style legend (TPC vs BTOF vs ETOF) -- NOT one entry per rapidity bin, see
  // header comment. Built from dummy histograms so the legend swatches show the
  // correct marker style independent of which bins actually had data.
  TLegend* styleLegend = new TLegend(0.70,0.75,0.93,0.91);
  styleLegend->SetBorderSize(0);
  styleLegend->SetFillStyle(0);
  styleLegend->SetTextSize(0.035);
  TH1D* tpcDummy = new TH1D(Form("tpcDummy_%s_%s_%d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex),"",1,0,1);
  TH1D* btofDummy = new TH1D(Form("btofDummy_%s_%s_%d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex),"",1,0,1);
  TH1D* etofDummy = new TH1D(Form("etofDummy_%s_%s_%d",a_speciesName.c_str(),a_chargeLabel.c_str(),a_centIndex),"",1,0,1);
  tpcDummy->SetMarkerStyle(20); tpcDummy->SetMarkerColor(kBlack);
  btofDummy->SetMarkerStyle(24); btofDummy->SetMarkerColor(kBlack);
  etofDummy->SetMarkerStyle(25); etofDummy->SetMarkerColor(kBlack);
  if(a_tpcHist)  styleLegend->AddEntry(tpcDummy,"TPC (dE/dx)","P");
  if(a_btofHist) styleLegend->AddEntry(btofDummy,"BTOF (1/#beta)","P");
  if(a_etofHist) styleLegend->AddEntry(etofDummy,"ETOF (1/#beta)","P");
  styleLegend->Draw();

  // Compact color-key text instead of a per-bin legend (see header comment).
  TLatex colorKey;
  colorKey.SetNDC();
  colorKey.SetTextSize(0.03);
  colorKey.DrawLatex(0.15,0.85,"color: rapidity bin index 1 (dark) #rightarrow last bin (bright)");
}

//_______________________________________________________________________________
void PresentZFitterSpectra(string a_zFitterOutputFile, int a_partIndex = 0, int a_centIndex = 0, string a_outDir = "."){

  const string speciesNames[3] = {"Pion","Kaon","Proton"}; // matches ParticleInfo's speciesIndex ordering for the 3 species this analysis fits
  if(a_partIndex < 0 || a_partIndex > 2){
    cout << "!!! a_partIndex must be 0 (Pion), 1 (Kaon), or 2 (Proton). Got " << a_partIndex << "." << endl;
    return;
  }
  string speciesName = speciesNames[a_partIndex];

  TFile* inFile = new TFile(a_zFitterOutputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_zFitterOutputFile << endl;
    return;
  }

  TH2* tpcPlus  = (TH2*) inFile->Get(Form("DeDx_FitData/%s/dEdxSpectra_%s_Plus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex));
  TH2* tpcMinus = (TH2*) inFile->Get(Form("DeDx_FitData/%s/dEdxSpectra_%s_Minus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex));
  TH2* btofPlus  = (TH2*) inFile->Get(Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%s_Plus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex));
  TH2* btofMinus = (TH2*) inFile->Get(Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%s_Minus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex));
  TH2* etofPlus  = (TH2*) inFile->Get(Form("ETOF_FitData/%s/InvBetaETOFSpectra_%s_Plus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex));
  TH2* etofMinus = (TH2*) inFile->Get(Form("ETOF_FitData/%s/InvBetaETOFSpectra_%s_Minus_Cent%02d",speciesName.c_str(),speciesName.c_str(),a_centIndex));

  if(!tpcPlus && !tpcMinus && !btofPlus && !btofMinus && !etofPlus && !etofMinus){
    cout << "!!! No TPC, BTOF, or ETOF spectra found for " << speciesName << ", Cent" << Form("%02d",a_centIndex)
         << " in " << a_zFitterOutputFile << ". Either that centrality/species wasn't fit yet in this ZFitter "
         << "run, or this file predates a completed run for it (check the ZFitter log for errors)." << endl;
    return;
  }
  if(!btofPlus && !btofMinus) cout << "Note: no BTOF spectra found -- this looks like a TPC-only ZFitter run so far. Drawing TPC only." << endl;
  if(!etofPlus && !etofMinus) cout << "Note: no ETOF spectra found -- this analysis has ETOF available (m_hasETOF), but no ETOF ZFitter run has completed yet, or none exists in this file." << endl;

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(Form("PresentZFitterSpectra_%s_Cent%02d",speciesName.c_str(),a_centIndex),
                            Form("%s ZFitter spectra, Cent%02d",speciesName.c_str(),a_centIndex),2000,900);
  c->Divide(2,1);

  c->cd(1);
  drawOneChargeSign(tpcPlus,btofPlus,etofPlus,speciesName,"Plus",a_centIndex);

  c->cd(2);
  drawOneChargeSign(tpcMinus,btofMinus,etofMinus,speciesName,"Minus",a_centIndex);

  string outPath = a_outDir + Form("/PresentZFitterSpectra_%s_Cent%02d.png",speciesName.c_str(),a_centIndex);
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
