// macros/FitZdcSingleNeutronPeak_OO200.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS (2026-08-28): follow-up to macros/CheckZdcNeutronCapture_OO200.C. Andrew
// asked how to convert raw ZDC Sum ADC into an actual neutron multiplicity. There is no
// official STAR calibration constant available for this dataset (confirmed: "I don't know
// of a calibration constant for this dataset"), so this macro derives one empirically from
// the data itself: it locates the single-neutron ("1n") peak in the already-produced
// zoomed ZDC ADC spectrum, fits it, and uses that fitted centroid (above an estimated
// pedestal) as an empirical ADC-per-neutron scale, ADC_1n. That scale is then used to
// relabel the Sum ADC axis of the already-produced histograms into a neutron-multiplicity
// axis. This is explicitly a SELF-DERIVED, DATA-DRIVEN estimate, not an official
// STAR/collaboration calibration -- every printed/plotted result says so.
//
// INPUT FILES (reads two files already produced by CheckZdcNeutronCapture_OO200.C -- this
// macro does NOT touch the picoDst filelist again, deliberately, per Andrew's request):
//   1) a_zoomedFile (default "zdc_qa_zoomed.root"): hZdcEast/hZdcWest booked 0-300 ADC,
//      2 ADC/bin (from a_adcMax=300, a_adcBinWidth=2 run) -- fine enough binning to
//      resolve the 1n peak, used ONLY for the peak fit.
//   2) a_fullRangeFile (default "zdc_neutron_capture_check.root"): hZdcEast/hZdcWest
//      booked 0-4000 ADC, 20 ADC/bin (the original default-binning run) -- coarser, but
//      covers the FULL Sum ADC range actually populated in the data (the first full run's
//      plot showed real entries out to ~2000-2100 ADC on both arms). This is the histogram
//      actually relabeled into the neutron-multiplicity axis and written out, specifically
//      because the zoomed 0-300 ADC histogram alone would silently truncate every event
//      above ~300/ADC_1n neutrons -- see CAVEATS below.
//
// METHODOLOGY:
//   Pedestal estimate: the zoomed spectrum shows a large pedestal spike sitting at/near
//   ADC=0 (visible in both the first full-range run and the zoomed rerun). Rather than
//   force a Gaussian onto what looked like a spike clipped against the ADC=0 edge (a
//   genuine electronics-pedestal Gaussian can't be trusted to fit well against a hard
//   boundary), this macro just takes the bin of maximum content within
//   [0, a_pedestalSearchMax] as the pedestal estimate (ADC_ped) for each side. Printed to
//   console either way, so it's checkable.
//
//   1n peak (REVISED 2026-08-28 after the first version's fit visibly failed -- see below):
//   the first version of this macro fit a single Gaussian over a fixed, eyeballed window
//   ([8,40] ADC). That failed badly on real data: East came back mean=8.7, sigma=10.85,
//   chi2/ndf=190.7/13 -- the fitted mean sat right at the window's lower edge and sigma was
//   almost as large as the mean itself, both signs the fit was pulled onto the still-falling
//   tail of the huge pedestal spike (confirmed visually: the linear-scale QA plot showed the
//   pedestal bin near 2x10^4 vs. the real 1n bump only ~2x10^3, so a wide fit window starting
//   right after the pedestal still contains far more pedestal-tail counts than bump counts).
//   Fixed by finding the bump PROGRAMMATICALLY instead of trusting an eyeballed range: scan
//   for the bin of maximum content in (a_bumpScanMin, a_bumpScanMax] (default 6-45 ADC, i.e.
//   starting right after the pedestal search region and ending before the ADC~50-70 dip/
//   second-hump region visible in the already-produced log-scale plot), then fit a Gaussian
//   in a NARROW, symmetric window of +/-a_fitHalfWidth (default 12 ADC) centered on that
//   located bin. This keeps the fit from reaching back into the pedestal tail or forward into
//   the next (2n+) hump. CHECK THE OUTPUT PNG (now log-y, fit overlaid on the histogram, with
//   the scan/fit window marked) BEFORE TRUSTING THE PRINTED ADC_1n VALUE -- a chi2/ndf
//   warning is also printed if the fit quality still looks bad. If the located bump or fit
//   still look wrong in the PNG, narrow a_fitHalfWidth or adjust a_bumpScanMin/a_bumpScanMax
//   to bracket just the real bump as it appears there.
//
//   ADC_1n (effective, pedestal-subtracted single-neutron scale) = (fitted 1n peak mean)
//   - ADC_ped, computed separately for East and West (the two arms are NOT assumed to
//   share a scale -- the earlier full-range plot showed East/West have different ADC
//   ranges/means, consistent with the fxt-blue trigger's forward-track asymmetry Andrew
//   pointed out, not necessarily a shared gain).
//
//   Multiplicity relabeling: for a Sum ADC histogram with UNIFORM binning, n = (ADC -
//   ADC_ped)/ADC_1n is an affine (linear) map of the x-axis, so a bin-for-bin content copy
//   into a new histogram with rescaled axis limits reproduces the identical per-event
//   assignment as would applying that formula event-by-event -- no need to re-read
//   per-event data. Applied here to a_fullRangeFile's hZdcEast/hZdcWest (NOT the zoomed
//   ones -- see INPUT FILES above) so the multiplicity histogram covers the full range of
//   Sum ADC actually seen in the data, at the coarser 20-ADC/bin resolution (i.e. bin
//   widths in neutron units are 20/ADC_1n, roughly 0.5-1 neutron/bin depending on side --
//   fine for distinguishing low multiplicities, blurring at higher n where physically the
//   peaks overlap anyway from photostatistics).
//
// CAVEATS (read before using ADC_1n or the multiplicity histograms downstream):
//   - This is an empirical, self-derived calibration from this exact dataset's own ADC
//     spectrum, NOT an official STAR/collaboration ZDC calibration constant. Treat any
//     neutron count derived from it as approximate.
//   - No true pedestal-width/resolution term is propagated -- ADC_ped is a single bin
//     center, not a fitted mean+sigma, because the pedestal looked clipped/asymmetric
//     against the ADC=0 edge rather than a clean Gaussian (see METHODOLOGY above). If the
//     printed pedestal estimate looks wrong (e.g. sitting well above where the PNG shows
//     the spike), that's a sign to inspect the zoomed histogram's low-ADC bins directly.
//   - The 1n peak fit window is now programmatically centered on the tallest bin found in
//     the scan region (see REVISED METHODOLOGY), which fixed the first version's failure --
//     but "tallest bin in range" can still land somewhere wrong if the real spectrum shape
//     is unusual. Always check the fit-overlay PNG (and its printed chi2/ndf) before
//     trusting ADC_1n.
//   - The n=0 bin of the resulting multiplicity histogram is dominated by whatever landed
//     in/near the pedestal itself (i.e. "no resolvable neutron," not necessarily "exactly
//     zero neutrons" -- ZDC threshold effects near pedestal are not modeled here).
//
// USAGE: this macro only touches TFile/TH1/TF1/TCanvas/TLegend -- it never includes or
// links against libStPicoDst.so, unlike CheckZdcNeutronCapture_OO200.C. Reasoning (not
// yet confirmed in this environment): none of the earlier ACLiC/DYLD_LIBRARY_PATH/
// AddLinkedLibs build gymnastics documented in that macro were about ROOT/cling itself --
// they were specifically about StPicoDstReader/libStPicoDst.so symbol resolution. Since
// this macro never references those symbols, plain cling interpretation should work
// without any of that. Try the plain form first; if it unexpectedly hits a linker/symbol
// error, fall back to the ACLiC form (second line below) the same way
// CheckZdcNeutronCapture_OO200.C needed to.
//
//   root -l -q 'macros/FitZdcSingleNeutronPeak_OO200.C("zdc_qa_zoomed.root","zdc_neutron_capture_check.root")'
//
//   // fallback ONLY if the plain form errors out -- not expected to be needed:
//   root -l -q -e 'gSystem->AddLinkedLibs("-L./bin -lStPicoDst")' 'macros/FitZdcSingleNeutronPeak_OO200.C+("zdc_qa_zoomed.root","zdc_neutron_capture_check.root")'
//
//   // to narrow the fit window (radius in ADC around the auto-located bump) if it still
//   // looks wrong in zdc_1n_peak_fit_qa.png, e.g. a tighter +/-8 ADC window:
//   root -l -q 'macros/FitZdcSingleNeutronPeak_OO200.C("zdc_qa_zoomed.root","zdc_neutron_capture_check.root",".", "zdc_neutron_multiplicity.root",6.0,6.0,45.0,8.0)'

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

// Finds the bin of maximum content in (0, xMax] on hist, returns its bin CENTER.
// Used for the pedestal estimate -- see METHODOLOGY header comment for why this is a
// simple max-bin lookup rather than a Gaussian fit.
double FindPedestalEstimate(TH1* hist, double xMax){
  int binMax = hist->GetXaxis()->FindBin(xMax);
  int iBest = 1;
  double best = -1;
  for(int i = 1; i <= binMax; i++){
    double c = hist->GetBinContent(i);
    if(c > best){ best = c; iBest = i; }
  }
  return hist->GetBinCenter(iBest);
}

// Finds the bin of maximum content in (scanMin, scanMax] on hist, returns its bin CENTER.
// Used to locate the real 1n bump programmatically rather than trusting an eyeballed fit
// window -- see the REVISED METHODOLOGY header comment for why the first version's fixed,
// wide window failed.
double FindBumpCenter(TH1* hist, double scanMin, double scanMax){
  int binLo = hist->GetXaxis()->FindBin(scanMin) + 1; // start just past scanMin's own bin
  int binHi = hist->GetXaxis()->FindBin(scanMax);
  int iBest = binLo;
  double best = -1;
  for(int i = binLo; i <= binHi; i++){
    double c = hist->GetBinContent(i);
    if(c > best){ best = c; iBest = i; }
  }
  return hist->GetBinCenter(iBest);
}

// Fits a Gaussian to hist in [fitMin,fitMax] (a narrow window centered on a programmatically
// located bump -- see FindBumpCenter), draws hist+fit, returns the fitted mean. Prints
// chi2/ndf and warns if it looks bad (loose threshold -- this is a QA aid, not a hard cut;
// always look at the PNG too).
double FitSingleNeutronPeak(TH1* hist, double fitMin, double fitMax, const char* fname){
  TF1* f = new TF1(fname, "gaus", fitMin, fitMax);
  f->SetParameter(1, 0.5*(fitMin+fitMax)); // seed the mean at the window center (the located bump)
  int fitStatus = hist->Fit(f, "RQS", "", fitMin, fitMax);
  double mean = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double chi2 = f->GetChisquare();
  int ndf = f->GetNDF();
  cout << "  Gaussian fit [" << fitMin << "," << fitMax << "]: mean=" << mean
       << "  sigma=" << sigma << "  chi2/ndf=" << chi2 << "/" << ndf
       << "  status=" << fitStatus << endl;
  if(fitStatus != 0){
    cout << "FitZdcSingleNeutronPeak_OO200 :: WARNING - fit '" << fname
         << "' did not converge cleanly (status " << fitStatus
         << "). Printed ADC_1n for this side is unreliable -- inspect the PNG." << endl;
  }
  if(ndf > 0 && chi2/ndf > 5.0){
    cout << "FitZdcSingleNeutronPeak_OO200 :: WARNING - fit '" << fname
         << "' has a large chi2/ndf (" << chi2/ndf << "). The Gaussian may be a poor "
         << "description of this window (e.g. still catching a falling shoulder rather than "
         << "an isolated bump) -- inspect the PNG before trusting ADC_1n for this side."
         << endl;
  }
  return mean;
}

// Bin-for-bin copy of srcHist into a new histogram whose x-axis has been affinely
// rescaled n = (ADC - pedestal)/adc1n. Valid because srcHist has uniform binning -- see
// METHODOLOGY header comment. srcHist must NOT be modified by the caller afterward (we
// read its bin contents/errors directly).
TH1D* MakeMultiplicityHisto(TH1* srcHist, double pedestal, double adc1n, const char* newName, const char* newTitle){
  int nb = srcHist->GetNbinsX();
  double oldXmin = srcHist->GetXaxis()->GetXmin();
  double oldXmax = srcHist->GetXaxis()->GetXmax();
  double newXmin = (oldXmin - pedestal) / adc1n;
  double newXmax = (oldXmax - pedestal) / adc1n;
  TH1D* h = new TH1D(newName, newTitle, nb, newXmin, newXmax);
  for(int i = 0; i <= nb+1; i++){ // include under/overflow so Integral()/entries stay consistent
    h->SetBinContent(i, srcHist->GetBinContent(i));
    h->SetBinError(i, srcHist->GetBinError(i));
  }
  h->SetEntries(srcHist->GetEntries());
  return h;
}

void FitZdcSingleNeutronPeak_OO200(string a_zoomedFile = "zdc_qa_zoomed.root",
                                    string a_fullRangeFile = "zdc_neutron_capture_check.root",
                                    string a_outDir = ".",
                                    string a_outFile = "zdc_neutron_multiplicity.root",
                                    double a_pedestalSearchMax = 6.0,  // ADC; look for the pedestal bin in (0, this]
                                    double a_bumpScanMin = 6.0,        // ADC; start of the window scanned for the real 1n bump (after the pedestal)
                                    double a_bumpScanMax = 45.0,       // ADC; end of that scan window (before the ADC~50-70 dip seen in the log-scale plot)
                                    double a_fitHalfWidth = 12.0){     // ADC; Gaussian fit window is [bumpCenter-this, bumpCenter+this] -- see REVISED METHODOLOGY

  gSystem->mkdir(a_outDir.c_str(), true);

  // ---- Open the fine-binned (zoomed) file, used only for the peak fit ----
  TFile* fZoom = TFile::Open(a_zoomedFile.c_str());
  if(!fZoom || fZoom->IsZombie()){
    cout << "FitZdcSingleNeutronPeak_OO200 :: Error - cannot open " << a_zoomedFile << endl;
    return;
  }
  TH1D* hZoomEast = (TH1D*) fZoom->Get("hZdcEast");
  TH1D* hZoomWest = (TH1D*) fZoom->Get("hZdcWest");
  if(!hZoomEast || !hZoomWest){
    cout << "FitZdcSingleNeutronPeak_OO200 :: Error - hZdcEast/hZdcWest not found in "
         << a_zoomedFile << endl;
    return;
  }
  hZoomEast->SetDirectory(0);
  hZoomWest->SetDirectory(0);

  // ---- Open the full-range (coarse) file, this is what actually gets relabeled ----
  TFile* fFull = TFile::Open(a_fullRangeFile.c_str());
  if(!fFull || fFull->IsZombie()){
    cout << "FitZdcSingleNeutronPeak_OO200 :: Error - cannot open " << a_fullRangeFile << endl;
    return;
  }
  TH1D* hFullEast = (TH1D*) fFull->Get("hZdcEast");
  TH1D* hFullWest = (TH1D*) fFull->Get("hZdcWest");
  if(!hFullEast || !hFullWest){
    cout << "FitZdcSingleNeutronPeak_OO200 :: Error - hZdcEast/hZdcWest not found in "
         << a_fullRangeFile << endl;
    return;
  }
  hFullEast->SetDirectory(0);
  hFullWest->SetDirectory(0);

  cout << "==== Pedestal estimate (max bin in (0," << a_pedestalSearchMax << "] ADC) ====" << endl;
  double pedEast = FindPedestalEstimate(hZoomEast, a_pedestalSearchMax);
  double pedWest = FindPedestalEstimate(hZoomWest, a_pedestalSearchMax);
  cout << "  East pedestal estimate: " << pedEast << " ADC" << endl;
  cout << "  West pedestal estimate: " << pedWest << " ADC" << endl;

  cout << "==== Locating the real 1n bump (max bin in (" << a_bumpScanMin << ","
       << a_bumpScanMax << "] ADC) ====" << endl;
  double bumpEast = FindBumpCenter(hZoomEast, a_bumpScanMin, a_bumpScanMax);
  double bumpWest = FindBumpCenter(hZoomWest, a_bumpScanMin, a_bumpScanMax);
  cout << "  East bump located at: " << bumpEast << " ADC" << endl;
  cout << "  West bump located at: " << bumpWest << " ADC" << endl;

  double fitMinEast = std::max(a_bumpScanMin, bumpEast - a_fitHalfWidth);
  double fitMaxEast = std::min(a_bumpScanMax, bumpEast + a_fitHalfWidth);
  double fitMinWest = std::max(a_bumpScanMin, bumpWest - a_fitHalfWidth);
  double fitMaxWest = std::min(a_bumpScanMax, bumpWest + a_fitHalfWidth);

  cout << "==== Single-neutron (1n) peak fit, zoomed histograms (narrow window around located bump) ====" << endl;
  cout << "East:" << endl;
  double meanEast = FitSingleNeutronPeak(hZoomEast, fitMinEast, fitMaxEast, "f1nEast");
  cout << "West:" << endl;
  double meanWest = FitSingleNeutronPeak(hZoomWest, fitMinWest, fitMaxWest, "f1nWest");

  double adc1nEast = meanEast - pedEast;
  double adc1nWest = meanWest - pedWest;
  cout << "==== Empirical ADC_1n (pedestal-subtracted single-neutron scale) ====" << endl;
  cout << "  East ADC_1n = " << adc1nEast << "  (NOT an official STAR calibration constant -- self-derived, see header comment)" << endl;
  cout << "  West ADC_1n = " << adc1nWest << "  (NOT an official STAR calibration constant -- self-derived, see header comment)" << endl;
  if(adc1nEast <= 0 || adc1nWest <= 0){
    cout << "FitZdcSingleNeutronPeak_OO200 :: WARNING - a nonpositive ADC_1n means the peak fit "
         << "landed at or below the pedestal estimate. The located bump/fit window is almost "
         << "certainly wrong for this side -- check the PNG and rerun with adjusted "
         << "a_bumpScanMin/a_bumpScanMax/a_fitHalfWidth before trusting anything downstream."
         << endl;
  }

  // ---- Build multiplicity histograms from the FULL-RANGE file (not the zoomed one) ----
  TH1D* hMultEast = MakeMultiplicityHisto(hFullEast, pedEast, adc1nEast,
    "hNeutronMultEast", "ZDC East -- empirical neutron multiplicity (self-derived calibration, not official);N_{n} (east, empirical);Events");
  TH1D* hMultWest = MakeMultiplicityHisto(hFullWest, pedWest, adc1nWest,
    "hNeutronMultWest", "ZDC West -- empirical neutron multiplicity (self-derived calibration, not official);N_{n} (west, empirical);Events");

  // ---- Write everything out ----
  string outFilePath = a_outDir + "/" + a_outFile;
  TFile* outFile = new TFile(outFilePath.c_str(), "RECREATE");
  outFile->cd();
  hMultEast->Write();
  hMultWest->Write();
  outFile->Close();
  cout << "Wrote " << outFilePath << " (hNeutronMultEast, hNeutronMultWest)" << endl;

  // ---- QA plot: fit overlay on the zoomed histograms, LOG-Y (required -- the pedestal bin
  // is ~2x10^4 vs. the real bump's ~2x10^3, so a linear axis flattens the bump to invisible
  // and hides exactly the kind of bad fit the first version of this macro produced) ----
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas* cFit = new TCanvas("cFitZdc1n", "ZDC 1n peak fit QA", 1200, 600);
  cFit->Divide(2,1);
  cFit->cd(1)->SetLogy();
  hZoomEast->SetLineColor(kRed+1);
  hZoomEast->SetTitle(Form("East;Sum ADC;Events   (ped=%.1f, bump@%.1f, fit[%.1f,%.1f], ADC_{1n}=%.1f)",
                           pedEast, bumpEast, fitMinEast, fitMaxEast, adc1nEast));
  hZoomEast->Draw();
  cFit->cd(2)->SetLogy();
  hZoomWest->SetLineColor(kBlue+1);
  hZoomWest->SetTitle(Form("West;Sum ADC;Events   (ped=%.1f, bump@%.1f, fit[%.1f,%.1f], ADC_{1n}=%.1f)",
                           pedWest, bumpWest, fitMinWest, fitMaxWest, adc1nWest));
  hZoomWest->Draw();
  string fitPngPath = a_outDir + "/zdc_1n_peak_fit_qa.png";
  cFit->SaveAs(fitPngPath.c_str());
  cout << "Wrote " << fitPngPath << " -- CHECK THIS before trusting ADC_1n above (see header comment)." << endl;

  // ---- Result plot: the derived multiplicity histograms ----
  TCanvas* cMult = new TCanvas("cNeutronMult", "ZDC empirical neutron multiplicity", 900, 700);
  cMult->SetLogy();
  hMultEast->SetLineColor(kRed+1);
  hMultWest->SetLineColor(kBlue+1);
  hMultEast->SetTitle("ZDC empirical neutron multiplicity (self-derived, NOT an official calibration);N_{n} (empirical);Events");
  double xMaxDisplay = 40; // just for the drawn range, same convention as GetRefMultHisto.C's a_xMaxDisplay
  hMultEast->GetXaxis()->SetRangeUser(0, xMaxDisplay);
  hMultEast->Draw("hist");
  hMultWest->Draw("hist same");
  TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
  leg->AddEntry(hMultEast, "East", "l");
  leg->AddEntry(hMultWest, "West", "l");
  leg->Draw();
  string multPngPath = a_outDir + "/zdc_neutron_multiplicity.png";
  cMult->SaveAs(multPngPath.c_str());
  cout << "Wrote " << multPngPath << endl;
}
