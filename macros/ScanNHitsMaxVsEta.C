// ScanNHitsMaxVsEta.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: follow-up to ProjectNHitsMaxSlices.C, which showed that the six
// hand-picked eta slices aren't uniform -- eta=0/0.5/0.8 each have a single sharp
// spike in NHitsMax (a near-deterministic geometric ceiling), but eta=1.1 is instead
// a broad, smeared plateau (near-peak fraction collapses from ~0.40 to ~0.12), before
// re-sharpening by eta=1.4 at a much lower ceiling value. Per Andrew (2026-07-23):
// the drop in "spikiness" around |eta|~0.9 is the track's path crossing the physical
// CORNER of the TPC (where the barrel meets the endcap plane -- exactly clearing it
// or not becomes sensitive to small differences in dip angle/vertex-z rather than
// eta alone, smearing the ceiling into a spread of values instead of one), and the
// ceiling value around NHitsMax~40 that shows up past that point is the INNER/OUTER
// TPC SECTOR BOUNDARY (a track that exits before completing the outer sector's
// padrows is capped there instead of at the full inner+outer total near 72).
//
// This macro scans that behavior CONTINUOUSLY across eta (rather than six fixed
// slices) using a sliding window, and plots two profiles stacked on one canvas:
//   (top)    peak (mode) NHitsMax value vs eta -- should sit flat at the full
//            inner+outer ceiling (~72) in the barrel, and show where/how sharply it
//            steps down to the inner-only-outer-partial regime (~40) once the track
//            starts exiting before finishing the outer sector.
//   (bottom) "spikiness" = fraction of that window's entries within +-1 hit of the
//            peak -- should stay high (single dominant ceiling value) except right
//            around the TPC-corner crossing, where it should dip sharply (many
//            comparably-populated ceiling values instead of one).
// Reference lines are drawn at |eta|=0.9 (TPC corner, per Andrew) and NHitsMax=40
// (inner/outer boundary, per Andrew) so the physical explanation is visible directly
// on the plot rather than requiring a caption.
//
// Reads (same histogram as ProjectNHitsMaxSlices.C / PresentHitLossEta.C):
//   eta_nHitsMax -- TH2I of eta vs N_hits^max, written at top level by PicoBinner.cxx
//
// Usage:
//   root -l -q 'ScanNHitsMaxVsEta.C("/Users/aliggett/data/OO/yieldHistos_OO200_pion.root")'
//   root -l -q 'ScanNHitsMaxVsEta.C("/path/to/file.root","/path/to/outDir",0.02,0.025)'
//
// a_etaStep: spacing between scan points. a_etaHalfWidth: half-width of the sliding
// window used at each scan point (same convention/default as ProjectNHitsMaxSlices.C
// so the two are directly comparable). Points whose window has fewer than
// a_minEntries total entries are skipped (avoids noisy/meaningless peak-finding in
// near-empty acceptance-edge bins).

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

void ScanNHitsMaxVsEta(string a_inputFile,
                        string a_outDir = ".",
                        double a_etaStep = 0.02,
                        double a_etaHalfWidth = 0.025,
                        double a_minEntries = 1000){

  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* etaNHitsMax = (TH2*) inFile->Get("eta_nHitsMax");
  if(!etaNHitsMax){
    cout << "!!! eta_nHitsMax not found in " << a_inputFile
         << " -- run FindHistoInFile.C to see what's actually in this file "
         << "(a_makeBasicHistos may have been false when it was produced)."
         << endl;
    return;
  }

  double etaMin = etaNHitsMax->GetXaxis()->GetXmin();
  double etaMax = etaNHitsMax->GetXaxis()->GetXmax();

  vector<double> etaVals, peakVals, spikinessVals;

  for(double etaCenter = etaMin + a_etaHalfWidth; etaCenter <= etaMax - a_etaHalfWidth; etaCenter += a_etaStep){
    int binLow  = etaNHitsMax->GetXaxis()->FindBin(etaCenter - a_etaHalfWidth);
    int binHigh = etaNHitsMax->GetXaxis()->FindBin(etaCenter + a_etaHalfWidth);

    TH1D* proj = etaNHitsMax->ProjectionY("scanTmpProj", binLow, binHigh);
    double totalEntries = proj->Integral();
    if(totalEntries < a_minEntries){
      delete proj;
      continue;
    }

    int peakBin = proj->GetMaximumBin();
    double peakVal = proj->GetXaxis()->GetBinCenter(peakBin);
    double peakEntries = proj->GetBinContent(peakBin);
    double nearPeakEntries = peakEntries
      + proj->GetBinContent(peakBin - 1)
      + proj->GetBinContent(peakBin + 1);
    double nearPeakFrac = nearPeakEntries / totalEntries;

    etaVals.push_back(etaCenter);
    peakVals.push_back(peakVal);
    spikinessVals.push_back(nearPeakFrac);

    delete proj;
  }

  cout << "Scanned " << etaVals.size() << " eta points (step=" << a_etaStep
       << ", half-width=" << a_etaHalfWidth << ", min entries=" << a_minEntries << ")" << endl;

  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outDir.c_str(), true);

  TGraph* peakGraph = new TGraph((int)etaVals.size(), &etaVals[0], &peakVals[0]);
  peakGraph->SetTitle(";#eta;peak (mode) N_{hits}^{max}");
  peakGraph->SetMarkerStyle(20);
  peakGraph->SetMarkerSize(0.5);
  peakGraph->SetMarkerColor(kBlue+2);
  peakGraph->SetLineColor(kBlue+2);

  TGraph* spikeGraph = new TGraph((int)etaVals.size(), &etaVals[0], &spikinessVals[0]);
  spikeGraph->SetTitle(";#eta;fraction within #pm1 hit of peak (\"spikiness\")");
  spikeGraph->SetMarkerStyle(20);
  spikeGraph->SetMarkerSize(0.5);
  spikeGraph->SetMarkerColor(kRed+1);
  spikeGraph->SetLineColor(kRed+1);

  TCanvas* c = new TCanvas("ScanNHitsMaxVsEta","NHitsMax ceiling + spikiness vs eta",900,900);
  c->Divide(1,2);

  c->cd(1);
  peakGraph->Draw("AP");
  peakGraph->GetYaxis()->SetRangeUser(0,80);
  // Reference line at NHitsMax=40 -- per Andrew, the inner/outer TPC sector boundary:
  // a track that exits before completing the outer sector's padrows is capped near
  // here instead of at the full inner+outer ceiling (~72).
  TLine* innerOuterLine = new TLine(etaMin, 40, etaMax, 40);
  innerOuterLine->SetLineColor(kGray+2);
  innerOuterLine->SetLineStyle(2);
  innerOuterLine->Draw("same");
  TLatex latTop;
  latTop.SetTextSize(0.035);
  latTop.SetTextColor(kGray+2);
  latTop.DrawLatex(etaMin+0.05, 41.5, "inner/outer TPC sector boundary");
  // Reference lines at |eta|=0.9 -- per Andrew, where the track's path crosses the
  // physical corner of the TPC (barrel/endcap edge).
  TLine* cornerLineNeg = new TLine(-0.9, 0, -0.9, 80);
  TLine* cornerLinePos = new TLine(0.9, 0, 0.9, 80);
  for(TLine* l : {cornerLineNeg, cornerLinePos}){
    l->SetLineColor(kGreen+2);
    l->SetLineStyle(2);
    l->Draw("same");
  }
  latTop.SetTextColor(kGreen+2);
  latTop.DrawLatex(0.92, 5, "TPC corner");
  latTop.DrawLatex(-0.9-0.30, 5, "TPC corner");

  c->cd(2);
  spikeGraph->Draw("AP");
  spikeGraph->GetYaxis()->SetRangeUser(0,0.6);
  TLine* cornerLineNeg2 = new TLine(-0.9, 0, -0.9, 0.6);
  TLine* cornerLinePos2 = new TLine(0.9, 0, 0.9, 0.6);
  for(TLine* l : {cornerLineNeg2, cornerLinePos2}){
    l->SetLineColor(kGreen+2);
    l->SetLineStyle(2);
    l->Draw("same");
  }
  TLatex latBot;
  latBot.SetTextSize(0.035);
  latBot.SetTextColor(kGreen+2);
  latBot.DrawLatex(0.92, 0.55, "TPC corner");
  latBot.DrawLatex(-0.9-0.30, 0.55, "TPC corner");

  string outPath = a_outDir + "/ScanNHitsMaxVsEta.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
