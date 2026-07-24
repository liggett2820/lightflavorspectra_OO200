// ProjectNHitsMaxSlices.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: direct follow-up to the "why are there red streaks along the sides
// of the NHitsMax-vs-eta trapezoid" question (see PresentNHitsEta.C /
// PresentNHitsMaxBySector.C). The working explanation is that track->nHitsMax() is a
// per-track quantity computed upstream by STAR's own tracking software (read straight
// off the PicoDST, see submodule/PicoDstReader_SL24y/StPicoTrack.h -- nothing in this
// repo computes it), and it's nearly deterministic in eta: almost all tracks at a
// given eta cross very close to the same number of padrows (set by fixed TPC
// geometry), so the 2D histogram's statistical weight piles up on that geometric
// ceiling, saturating the color scale there and leaving a much sparser "fan" of
// off-nominal tracks (vertex-z offset, sector-boundary clipping, dead padrows) below
// it. This macro tests that directly rather than just eyeballing the 2D plot: take a
// few narrow eta slices and look at the 1D NHitsMax distribution within each one. If
// the explanation is right, each slice should show a tall, narrow spike at the
// ceiling value with a long, low tail below -- not a broad, flat-ish distribution.
//
// Reads (same histogram PresentHitLossEta.C / PresentNHitsEta.C already use):
//   eta_nHitsMax -- TH2I of eta vs N_hits^max, written at top level by PicoBinner.cxx
//
// For each requested eta slice this prints, and annotates on the plot: the peak
// (mode) NHitsMax value, and what fraction of that slice's total entries sit at or
// within 1 hit of the peak -- a direct, quantitative "how spike-like is this" number,
// not just a visual impression.
//
// Usage:
//   root -l -q 'ProjectNHitsMaxSlices.C("/path/to/yieldHistos_OO200_pion.root")'
//   root -l -q 'ProjectNHitsMaxSlices.C("/path/to/file.root","/path/to/outDir")'
//
// Default slices cover both the flat-plateau region (eta near 0) and the region where
// the eta~1.1-1.2 hit-loss bump sits, so you can compare "normal" vs "bump" behavior
// directly. Pass your own via a_etaCenters/a_etaHalfWidth if you want different ones.

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

void ProjectNHitsMaxSlices(string a_inputFile,
                            string a_outDir = ".",
                            vector<double> a_etaCenters = {0.0, 0.5, 0.8, 1.1, 1.4, 1.6},
                            double a_etaHalfWidth = 0.025){

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

  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outDir.c_str(), true);

  int nSlices = (int) a_etaCenters.size();
  int nCols = (nSlices <= 3) ? nSlices : (int)TMath::Ceil(nSlices/2.0);
  int nRows = (nSlices <= 3) ? 1 : 2;

  TCanvas* c = new TCanvas("ProjectNHitsMaxSlices","NHitsMax 1D slices vs eta",
                            450*nCols, 450*nRows);
  c->Divide(nCols, nRows);

  cout << "eta slice results (half-width " << a_etaHalfWidth << "):" << endl;

  for(int i = 0; i < nSlices; i++){
    double etaCenter = a_etaCenters[i];
    double etaLow  = etaCenter - a_etaHalfWidth;
    double etaHigh = etaCenter + a_etaHalfWidth;

    int binLow  = etaNHitsMax->GetXaxis()->FindBin(etaLow);
    int binHigh = etaNHitsMax->GetXaxis()->FindBin(etaHigh);

    TH1D* proj = etaNHitsMax->ProjectionY(Form("proj_eta_%.3f", etaCenter), binLow, binHigh);
    proj->SetTitle(Form("#eta #in [%.3f, %.3f];N_{hits}^{max};entries", etaLow, etaHigh));
    proj->SetDirectory(0);

    double totalEntries = proj->Integral();
    int peakBin = proj->GetMaximumBin();
    double peakVal = proj->GetXaxis()->GetBinCenter(peakBin);
    double peakEntries = proj->GetBinContent(peakBin);

    // entries within +-1 hit of the peak, as a fraction of the slice's total -- the
    // "how spike-like is this" number.
    double nearPeakEntries = peakEntries
      + proj->GetBinContent(peakBin - 1)
      + proj->GetBinContent(peakBin + 1);
    double nearPeakFrac = (totalEntries > 0) ? nearPeakEntries / totalEntries : 0.0;

    cout << "  eta=" << etaCenter << ": peak at NHitsMax=" << peakVal
         << ", total entries=" << totalEntries
         << ", fraction within +-1 hit of peak=" << nearPeakFrac << endl;

    c->cd(i+1);
    gPad->SetLogy(1);
    proj->Draw("hist");

    TLine* peakLine = new TLine(peakVal, 0, peakVal, peakEntries);
    peakLine->SetLineColor(kRed);
    peakLine->SetLineStyle(2);
    peakLine->SetLineWidth(2);
    peakLine->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.DrawLatex(0.15, 0.85, Form("peak=%.0f", peakVal));
    lat.DrawLatex(0.15, 0.79, Form("#pm1 hit frac=%.2f", nearPeakFrac));
  }

  string outPath = a_outDir + "/ProjectNHitsMaxSlices.png";
  c->SaveAs(outPath.c_str());
  cout << "Wrote " << outPath << endl;
}
