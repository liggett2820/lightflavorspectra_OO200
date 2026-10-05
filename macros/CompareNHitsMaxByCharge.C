// CompareNHitsMaxByCharge.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: direct test of the leading hypothesis for the 2-3 diagonal
// "streaks" seen along each side of the eta_nHitsMax trapezoid (see
// PresentNHitsEta.C's middle panel, and ScanNHitsMaxVsEta.C for the already-
// understood main-ceiling / inner-outer-boundary / TPC-corner features). If those
// streaks are genuine fixed TPC structural crossings (e.g. field-cage support ribs,
// wheel-mount hardware) smeared into diagonal lines by track curvature, then positive
// and negative tracks -- which bend oppositely in STAR's field -- should cross the
// same structure at slightly different eta, shifting the streaks between charge
// signs. A purely eta-only geometric effect (no curvature dependence at all) would
// show identical streaks in both.
//
// PREREQUISITE: needs eta_nHitsMax_chargePlus/eta_nHitsMax_chargeMinus, added to
// PicoBinner.cxx 2026-07-23 -- rerun PicoBinner against a species file to populate
// these; they won't exist in any yield file produced before that change (run
// FindHistoInFile.C to check what's actually in a given file).
//
// Draws three things:
//   (1) the two 2D histograms (chargePlus, chargeMinus) side by side, same z-scale,
//       so you can visually compare the streak positions directly.
//   (2) the same continuous peak-vs-eta scan ScanNHitsMaxVsEta.C does, but as two
//       overlaid curves (blue=positive, red=negative) restricted to the eta range
//       where the streaks live (roughly |eta| in [0.9,1.6], past the TPC corner) --
//       this is the quantitative version of (1): if the ceiling trace itself is
//       offset in eta between the two charges, that's the curvature-shift signature.
//
// Usage:
//   root -l -q 'CompareNHitsMaxByCharge.C("/Users/aliggett/data/OO/yieldHistos_OO200_pion.root")'
//   root -l -q 'CompareNHitsMaxByCharge.C("/Users/aliggett/data/OO/yieldHistos_OO200_pion.root",".",100)'  // looser min-entries-per-window
//
// REVISION 1 (2026-07-29): first real run (default a_minEntries=500) showed the 2D
// side-by-side comparison overlapping (no visible charge-dependent streak shift), but
// the quantitative peak-vs-eta scan came back nearly all negative-track points in the
// |eta| in [0.9,1.65] streak region -- positive tracks only cleared the 500-entries-
// per-window bar at 2-3 points total, nowhere near enough to actually compare trends
// against the negative-track curve. Andrew asked to lower the threshold so positive
// tracks get a fair shot at the scan too. Exposed the previously-hardcoded
// minEntries=500 as a new a_minEntries parameter (default kept at 500, so omitting it
// reproduces the original run exactly) instead of picking a new hardcoded value --
// whether 500 was simply too strict, or there's a genuine positive/negative statistics
// asymmetry in this eta band worth its own follow-up, is easier to tell by trying a
// couple of values than by guessing one number up front.
//
// REVISION 2 (2026-07-29, same day): a_minEntries=100 and a_minEntries=20 produced
// BYTE-IDENTICAL scan output -- same point count, same positions, confirmed both
// visually and by diffing the PNGs. That rules out "500 was just a bit too strict":
// if there were a continuum of positive-track windows sitting somewhere between 20 and
// 500 entries, dropping the bar 25x (500->20) would have picked up at least some of
// them. It didn't, which means positive-track windows in the |eta| in [0.9,1.65]
// streak region are bimodal -- either comfortably above the threshold already, or
// essentially empty -- not marginal. Andrew asked to see the raw numbers directly
// instead of continuing to guess threshold values. Added scanTotalEntriesVsEta() (the
// same sliding window as scanPeakVsEta(), but returning the untouched, un-thresholded
// entry COUNT instead of a threshold-gated peak position) and a third output plot,
// CompareNHitsMaxByCharge_TotalEntries.png, overlaying raw total-entries-per-window vs
// eta for both charges on one log-y panel spanning the full streak region -- so the
// actual shape of the positive/negative statistics gap (gradual falloff vs a hard
// cliff vs a real acceptance gap) is visible directly, rather than inferred from which
// scan points happened to clear an arbitrary bar.

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

// Same sliding-window peak scan as ScanNHitsMaxVsEta.C, factored out so it can run
// on both charge signs and return a TGraph for overlay comparison.
TGraph* scanPeakVsEta(TH2* a_hist, double a_etaLow, double a_etaHigh,
                       double a_etaStep, double a_etaHalfWidth, double a_minEntries){
  vector<double> etaVals, peakVals;
  for(double etaCenter = a_etaLow; etaCenter <= a_etaHigh; etaCenter += a_etaStep){
    int binLow  = a_hist->GetXaxis()->FindBin(etaCenter - a_etaHalfWidth);
    int binHigh = a_hist->GetXaxis()->FindBin(etaCenter + a_etaHalfWidth);
    TH1D* proj = a_hist->ProjectionY("chargeScanTmpProj", binLow, binHigh);
    double totalEntries = proj->Integral();
    if(totalEntries >= a_minEntries){
      int peakBin = proj->GetMaximumBin();
      double peakVal = proj->GetXaxis()->GetBinCenter(peakBin);
      etaVals.push_back(etaCenter);
      peakVals.push_back(peakVal);
    }
    delete proj;
  }
  if(etaVals.empty()) return nullptr;
  return new TGraph((int)etaVals.size(), &etaVals[0], &peakVals[0]);
}

//_______________________________________________________________________________
// NEW in REVISION 2 (2026-07-29) -- see header comment above. Same sliding window as
// scanPeakVsEta(), but returns the RAW total-entries count for every window across the
// full range, no threshold and no peak-finding -- this is the direct look at whatever
// caused positive tracks to nearly vanish from the thresholded scan. A window with
// literally zero entries is floored to 0.3 (not 0, and not skipped) so it still plots
// as a visible point near the bottom of a log-y axis instead of silently disappearing
// or breaking the log scale -- the "count" axis label makes clear this is a display
// floor, not a real value.
TGraph* scanTotalEntriesVsEta(TH2* a_hist, double a_etaLow, double a_etaHigh,
                               double a_etaStep, double a_etaHalfWidth){
  vector<double> etaVals, countVals;
  for(double etaCenter = a_etaLow; etaCenter <= a_etaHigh; etaCenter += a_etaStep){
    int binLow  = a_hist->GetXaxis()->FindBin(etaCenter - a_etaHalfWidth);
    int binHigh = a_hist->GetXaxis()->FindBin(etaCenter + a_etaHalfWidth);
    TH1D* proj = a_hist->ProjectionY("chargeTotalEntriesTmpProj", binLow, binHigh);
    double totalEntries = proj->Integral();
    etaVals.push_back(etaCenter);
    countVals.push_back((totalEntries > 0) ? totalEntries : 0.3);
    delete proj;
  }
  if(etaVals.empty()) return nullptr;
  return new TGraph((int)etaVals.size(), &etaVals[0], &countVals[0]);
}

void CompareNHitsMaxByCharge(string a_inputFile, string a_outDir = ".", double a_minEntries = 500){
  TFile* inFile = new TFile(a_inputFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_inputFile << endl;
    return;
  }

  TH2* hPlus  = (TH2*) inFile->Get("eta_nHitsMax_chargePlus");
  TH2* hMinus = (TH2*) inFile->Get("eta_nHitsMax_chargeMinus");
  if(!hPlus || !hMinus){
    cout << "!!! eta_nHitsMax_chargePlus/Minus not found in " << a_inputFile
         << " -- this file predates the 2026-07-23 PicoBinner.cxx change that adds "
         << "them. Rerun PicoBinner against this species first, or run "
         << "FindHistoInFile.C to double check what's actually in this file."
         << endl;
    return;
  }

  gStyle->SetOptStat(0);
  gSystem->mkdir(a_outDir.c_str(), true);

  // (1) side-by-side 2D comparison, shared z-scale so colors mean the same thing in
  // both panels.
  double globalMax = TMath::Max(hPlus->GetMaximum(), hMinus->GetMaximum());
  hPlus->SetMaximum(globalMax);
  hMinus->SetMaximum(globalMax);
  hPlus->SetMinimum(1);
  hMinus->SetMinimum(1);

  TCanvas* c1 = new TCanvas("CompareNHitsMaxByCharge_2D","N_{hits}^{max} vs #eta by Charge Sign",1600,800);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogz(1);
  hPlus->SetTitle("Positive Tracks;#eta;N_{hits}^{max}");
  hPlus->Draw("colz");
  c1->cd(2);
  gPad->SetLogz(1);
  hMinus->SetTitle("Negative Tracks;#eta;N_{hits}^{max}");
  hMinus->Draw("colz");
  string outPath2D = a_outDir + "/CompareNHitsMaxByCharge_2D.png";
  c1->SaveAs(outPath2D.c_str());
  cout << "Wrote " << outPath2D << endl;

  // (2) overlaid peak-vs-eta scan, restricted to the streak region (roughly
  // |eta| in [0.9,1.6] on the positive side, mirrored on the negative side).
  double etaStep = 0.02, etaHalfWidth = 0.025, minEntries = a_minEntries;

  TGraph* peakPlusRight  = scanPeakVsEta(hPlus,  0.90, 1.65, etaStep, etaHalfWidth, minEntries);
  TGraph* peakMinusRight = scanPeakVsEta(hMinus, 0.90, 1.65, etaStep, etaHalfWidth, minEntries);
  TGraph* peakPlusLeft   = scanPeakVsEta(hPlus,  -1.65, -0.90, etaStep, etaHalfWidth, minEntries);
  TGraph* peakMinusLeft  = scanPeakVsEta(hMinus, -1.65, -0.90, etaStep, etaHalfWidth, minEntries);

  // NOTE (2026-07-29, overlap-diagnostic): the Scan plot appeared to show almost no
  // positive-track points even at very loose thresholds (down to 20), while the
  // separately-added TotalEntries plot proved there's no real statistics deficit for
  // positive tracks anywhere (millions of entries per window, full eta range). Since
  // both scanPeakVsEta() and scanTotalEntriesVsEta() use identical ProjectionY+
  // Integral logic over the same eta windows, that ruled out a statistics/threshold
  // bug. The remaining candidate: peakPlusRight/Left and peakMinusRight/Left are
  // nearly IDENTICAL (same peak NHitsMax value at each eta point, since the ceiling is
  // set by fixed TPC geometry, not charge), so the red (minus) markers -- drawn second,
  // "P same", same solid style/size as blue -- sit exactly on top of and fully hide the
  // blue (plus) markers underneath. That would look exactly like "almost all red, a
  // few blue" even though every point is really two coincident markers. Print a direct
  // point-by-point comparison to confirm this before trusting the visual read.
  auto dumpPeakOverlap = [](const char* label, TGraph* gPlus, TGraph* gMinus){
    if(!gPlus || !gMinus){
      cout << "CompareNHitsMaxByCharge: " << label << " -- plus or minus graph is NULL ("
           << "plus=" << (gPlus ? gPlus->GetN() : -1) << " pts, "
           << "minus=" << (gMinus ? gMinus->GetN() : -1) << " pts)" << endl;
      return;
    }
    cout << "CompareNHitsMaxByCharge: " << label << " -- plus has " << gPlus->GetN()
         << " pts, minus has " << gMinus->GetN() << " pts" << endl;
    int nCommon = TMath::Min(gPlus->GetN(), gMinus->GetN());
    int nCoincident = 0;
    for(int i = 0; i < nCommon; i++){
      double xp = gPlus->GetX()[i],  yp = gPlus->GetY()[i];
      double xm = gMinus->GetX()[i], ym = gMinus->GetY()[i];
      bool coincident = (TMath::Abs(xp - xm) < 1e-6) && (TMath::Abs(yp - ym) < 1e-6);
      if(coincident) nCoincident++;
      if(i < 3 || i == nCommon/2 || i >= nCommon - 3){
        cout << "    eta~" << xp << ": plus peak=" << yp << ", minus peak=" << ym
             << (coincident ? "  [COINCIDENT]" : "  [DIFFERS]") << endl;
      }
    }
    cout << "CompareNHitsMaxByCharge: " << label << " -- " << nCoincident << "/" << nCommon
         << " points exactly coincident (plus and minus peak at the same NHitsMax value)" << endl;
  };
  dumpPeakOverlap("peak scan, right (eta in [0.90,1.65])", peakPlusRight, peakMinusRight);
  dumpPeakOverlap("peak scan, left (eta in [-1.65,-0.90])", peakPlusLeft, peakMinusLeft);

  TCanvas* c2 = new TCanvas("CompareNHitsMaxByCharge_Scan","Peak N_{hits}^{max} vs #eta, by Charge Sign",1600,800);
  c2->Divide(2,1);

  // Marker styles deliberately differ (solid blue circle vs open red circle, red drawn
  // larger) so that a positive-track point sitting exactly under a negative-track point
  // remains visible -- with matching solid markers of the same size, an exact overlap
  // makes the blue point completely invisible, which is what was creating the illusion
  // of "almost no positive tracks" in earlier renders of this plot.
  c2->cd(1);
  TH2F* frameLeft = new TH2F("frameLeft",";#eta;Peak (Mode) N_{hits}^{max}",10,-1.7,-0.85,10,0,80);
  frameLeft->Draw();
  if(peakPlusLeft){ peakPlusLeft->SetMarkerStyle(20); peakPlusLeft->SetMarkerSize(1.0); peakPlusLeft->SetMarkerColor(kBlue+2); peakPlusLeft->SetLineColor(kBlue+2); peakPlusLeft->Draw("P same"); }
  if(peakMinusLeft){ peakMinusLeft->SetMarkerStyle(24); peakMinusLeft->SetMarkerSize(1.4); peakMinusLeft->SetMarkerColor(kRed+1); peakMinusLeft->SetLineColor(kRed+1); peakMinusLeft->Draw("P same"); }
  TLegend* legLeft = new TLegend(0.15,0.15,0.45,0.3);
  legLeft->AddEntry(peakPlusLeft,"positive tracks","p");
  legLeft->AddEntry(peakMinusLeft,"negative tracks","p");
  legLeft->Draw();

  c2->cd(2);
  TH2F* frameRight = new TH2F("frameRight",";#eta;Peak (Mode) N_{hits}^{max}",10,0.85,1.7,10,0,80);
  frameRight->Draw();
  if(peakPlusRight){ peakPlusRight->SetMarkerStyle(20); peakPlusRight->SetMarkerSize(1.0); peakPlusRight->SetMarkerColor(kBlue+2); peakPlusRight->SetLineColor(kBlue+2); peakPlusRight->Draw("P same"); }
  if(peakMinusRight){ peakMinusRight->SetMarkerStyle(24); peakMinusRight->SetMarkerSize(1.4); peakMinusRight->SetMarkerColor(kRed+1); peakMinusRight->SetLineColor(kRed+1); peakMinusRight->Draw("P same"); }
  TLegend* legRight = new TLegend(0.55,0.7,0.85,0.85);
  legRight->AddEntry(peakPlusRight,"positive tracks","p");
  legRight->AddEntry(peakMinusRight,"negative tracks","p");
  legRight->Draw();

  string outPathScan = a_outDir + "/CompareNHitsMaxByCharge_Scan.png";
  c2->SaveAs(outPathScan.c_str());
  cout << "Wrote " << outPathScan << endl;

  // (3) NEW in REVISION 2 (2026-07-29) -- raw total-entries-per-window vs eta, no
  // threshold, both charges overlaid on one log-y panel spanning the full streak
  // region. See header comment for why: the thresholded scan above couldn't
  // distinguish "positive tracks are genuinely sparse here" from "500 was just too
  // strict a bar," and lowering the bar 25x (500->20) changed nothing -- this plot
  // shows the actual entry counts directly instead of inferring them from which
  // points survived a cut.
  TGraph* totalPlus  = scanTotalEntriesVsEta(hPlus,  -1.65, 1.65, etaStep, etaHalfWidth);
  TGraph* totalMinus = scanTotalEntriesVsEta(hMinus, -1.65, 1.65, etaStep, etaHalfWidth);

  // NOTE (2026-07-29, same-day fix #1): the first version of this panel used a raw
  // `new TH2F(...); Draw();` as the axis frame under SetLogy() -- rendered a totally
  // empty plot. Switched to pad->DrawFrame(...) (the log-y TGraph-overlay idiom already
  // proven elsewhere in this repo, see PresentRawSpectraModifierOutput.C's
  // drawSpectraIntoPad()) -- but that ALSO came back empty, byte-identical to the
  // TH2F version. That ruled out the frame-drawing method as the actual cause.
  //
  // NOTE (2026-07-29, same-day fix #2): the real bug was the hardcoded y-axis ceiling
  // of 1e5. This histogram's total-entries-per-window is summed over the ENTIRE
  // NHitsMax range (unlike the threshold checks in scanPeakVsEta(), which only ever
  // confirmed a window had "at least" 500/100/20 entries, never an upper bound) -- for
  // a dataset this size, midrapidity windows plausibly hold far more than 1e5 entries
  // each, silently clipping every single point off the top of the frame with no error
  // printed (ROOT just doesn't draw points outside the frame's declared range).
  // Replaced the guessed fixed ceiling with an auto-computed one, same pattern already
  // used in ScanNHitsMaxVsEta.C/drawSpectraIntoPad() (derive the range from the actual
  // data instead of assuming a number) -- scans both graphs' real Y values for the max
  // and gives it 2x headroom, so this can't silently under-range again regardless of
  // how large this or a future input file's statistics are.
  // TEMP DIAGNOSTIC (2026-07-29) -- ground-truth dump of the graph objects themselves,
  // before any drawing/frame code touches them, since several rendering-side fixes in
  // a row have all produced byte-identical output. This will confirm or rule out the
  // data side directly instead of continuing to guess about pad/frame mechanics.
  cout << "CompareNHitsMaxByCharge: totalPlus is " << (totalPlus ? "non-null" : "NULL")
       << ", N points = " << (totalPlus ? totalPlus->GetN() : -1) << endl;
  cout << "CompareNHitsMaxByCharge: totalMinus is " << (totalMinus ? "non-null" : "NULL")
       << ", N points = " << (totalMinus ? totalMinus->GetN() : -1) << endl;
  if(totalPlus && totalPlus->GetN() > 0){
    cout << "  totalPlus first point: (" << totalPlus->GetX()[0] << ", " << totalPlus->GetY()[0] << ")" << endl;
    int midP = totalPlus->GetN()/2;
    cout << "  totalPlus middle point: (" << totalPlus->GetX()[midP] << ", " << totalPlus->GetY()[midP] << ")" << endl;
    int lastP = totalPlus->GetN()-1;
    cout << "  totalPlus last point: (" << totalPlus->GetX()[lastP] << ", " << totalPlus->GetY()[lastP] << ")" << endl;
  }
  if(totalMinus && totalMinus->GetN() > 0){
    cout << "  totalMinus first point: (" << totalMinus->GetX()[0] << ", " << totalMinus->GetY()[0] << ")" << endl;
    int midM = totalMinus->GetN()/2;
    cout << "  totalMinus middle point: (" << totalMinus->GetX()[midM] << ", " << totalMinus->GetY()[midM] << ")" << endl;
    int lastM = totalMinus->GetN()-1;
    cout << "  totalMinus last point: (" << totalMinus->GetX()[lastM] << ", " << totalMinus->GetY()[lastM] << ")" << endl;
  }

  double maxTotalEntries = 1.0;
  if(totalPlus)  for(int p = 0; p < totalPlus->GetN();  p++) maxTotalEntries = TMath::Max(maxTotalEntries, totalPlus->GetY()[p]);
  if(totalMinus) for(int p = 0; p < totalMinus->GetN(); p++) maxTotalEntries = TMath::Max(maxTotalEntries, totalMinus->GetY()[p]);
  cout << "CompareNHitsMaxByCharge: max total-entries-per-window found = " << maxTotalEntries << endl;

  // NOTE (2026-07-29, same-day fix #3): auto-ranging (fix #2) confirmed the real max
  // is ~1.58e7 -- so the auto-computed frame (0.2 to ~3.16e7) DID draw real points, but
  // almost the entire visible range (0.2 to 1e5) turned out to be empty/near-empty
  // territory, squeezing every actually-informative point into a thin sliver at the
  // very top -- which is why it still visually read as "blank" at a glance. Andrew
  // asked to zoom into the populated region directly: fixed range 1e5-1e8 (comfortably
  // bracketing the observed ~1.58e7 max with room to spare for other species/files
  // with somewhat higher statistics) instead of the auto-computed one. This trades the
  // auto-range's future-proofing (fix #2) for readability now that the actual scale of
  // this dataset is known -- if a future input's max total entries exceeds ~1e8, this
  // range will need revisiting (the console print above still reports the true max
  // every run, so that's easy to notice).
  TCanvas* c3 = new TCanvas("CompareNHitsMaxByCharge_TotalEntries","Raw Entries per Window vs #eta, by Charge Sign",1200,800);
  c3->SetLogy();
  TH1F* frameTotal = c3->DrawFrame(-1.7, 1e5, 1.7, 1e8,
                        ";#eta;Entries per Window (Log Scale)");
  // NOTE (2026-07-29, same-day fix #4): fix #3's DrawFrame(...,1e5,...,1e8,...) call
  // alone did NOT stick -- the rendered plot still came back ranged 1 to 1e5 (confirmed
  // by cropping and directly reading the axis tick labels), even though the source and
  // console output both showed fix #3's code was what actually ran. DrawFrame() returns
  // a TH1F whose OWN display min/max (a separate pair of knobs from the axis range
  // passed into the constructor call) can still fall back to auto-computed values under
  // SetLogy() in some ROOT versions. Setting them explicitly and redundantly here, plus
  // forcing the pad to recompute its range before drawing anything else, so there's no
  // remaining path for it to silently fall back to an auto/default range again.
  frameTotal->SetMinimum(1e5);
  frameTotal->SetMaximum(1e8);
  frameTotal->GetYaxis()->SetRangeUser(1e5, 1e8);
  c3->Modified();
  c3->Update();
  cout << "CompareNHitsMaxByCharge: TotalEntries frame range is now ["
       << frameTotal->GetMinimum() << ", " << frameTotal->GetMaximum() << "]" << endl;
  if(totalPlus){  totalPlus->SetMarkerStyle(20);  totalPlus->SetMarkerColor(kBlue+2); totalPlus->SetLineColor(kBlue+2);  totalPlus->Draw("PL same"); }
  if(totalMinus){ totalMinus->SetMarkerStyle(20); totalMinus->SetMarkerColor(kRed+1); totalMinus->SetLineColor(kRed+1); totalMinus->Draw("PL same"); }
  // NOTE: the old a_minEntries=500 reference line was dropped here -- with the
  // fix #3 range now starting at 1e5, y=500 sits far below the visible window and
  // would never actually appear on the plot.
  TLegend* legTotal = new TLegend(0.15,0.78,0.4,0.88);
  legTotal->AddEntry(totalPlus,"positive tracks","p");
  legTotal->AddEntry(totalMinus,"negative tracks","p");
  legTotal->Draw();

  string outPathTotal = a_outDir + "/CompareNHitsMaxByCharge_TotalEntries.png";
  c3->SaveAs(outPathTotal.c_str());
  cout << "Wrote " << outPathTotal << endl;

  cout << "\nIf positive/negative points visibly separate in eta (rather than lying "
       << "on top of each other) in the scan plot, that's the curvature-shift "
       << "signature -- the streaks come from crossing a fixed structural boundary "
       << "at a curvature-dependent eta. If they overlap, the streaks aren't "
       << "curvature-driven and the boundary hypothesis needs a different test "
       << "(e.g. pT-binned instead of just charge-sign)." << endl;
}
