// FitProtonDeuteronRatio_OO200.C
//
// *** VERSION 3 -- fixes a divergence bug found by actually running v2. ***
// v2 used ONE hardcoded location/scale guess pair applied blindly to every
// histogram (Plus proton, Minus proton, Plus/Minus deuteron seeds all shared
// a_protonLocGuess=0.0 / a_protonScaleGuess=0.3 / etc). Running it against the
// real file showed the proton calibration fit diverging to physically
// meaningless values (location=9962, scale=803660) even though the histogram's
// own axis only spans [-2,2] and its actual peak (confirmed directly:
// hz->GetXaxis()->GetBinCenter(hz->GetMaximumBin())) sits at Z=-0.917 --
// nowhere near the runaway fit result. A blind, shared initial guess with no
// parameter limit on location let MIGRAD wander off with nothing to anchor it.
//
// v3 fixes this at the root: every fit now derives its OWN location guess
// (peak bin of that specific histogram) and scale guess (that histogram's RMS)
// directly from the data, constrains the location parameter to stay within
// that histogram's own real axis range (so a repeat of the 9962 runaway is not
// possible), and uses the histogram's own range as the fit range by default
// instead of a hardcoded constant. It also adds an explicit minimum-entries
// gate (default 10) so a case like zTPC_DeuteronMinus_Isolated (confirmed by
// the user to have exactly 1 entry) is refused outright rather than fed into a
// 4-parameter fit that cannot possibly be meaningful with that little data.
//
// CONFIRMED PHYSICAL MODEL (source-traced, PicoBinner.cxx, + verified live against
// the real file by the user directly in ROOT)
// ---------------------------------------------------------------------------------
// This file is a DEUTERON-DRIVEN PicoBinner run (yields/Deuteron_space). Per
// PicoBinner.cxx:46-54 and the booking code at :803-823 plus the fill logic at
// :1887-1898 (isSpecificParticleByBTOF/ETOF, CutClass.cxx:1533-1543):
//
//   zTPC_Deuteron<Charge>_AllCent   -- the REAL, INCLUSIVE TPC dE/dx-vs-kinematics
//                                      histogram in the deuteron mass hypothesis
//                                      (y_d, mT-m_d, Z_TPC,d). Contains the actual
//                                      mixture of real deuterons PLUS proton (and
//                                      other) contamination bleeding into that band.
//                                      This is the histogram to FIT for a real yield.
//                                      CONFIRMED: 15883 entries (Plus), 8913 (Minus).
//
//   zTPC_Proton<Charge>_Isolated    -- a PURE, TOF-tagged (BTOF/ETOF 1/beta cut,
//                                      independent of TPC dE/dx) proton calibration
//                                      sample, with its TPC value filled under the
//                                      SAME deuteron kinematic hypothesis. This is a
//                                      SHAPE TEMPLATE for what proton contamination
//                                      looks like in the deuteron band -- NOT an
//                                      independent proton yield measurement.
//                                      CONFIRMED: 1239 entries (Plus), 65 (Minus).
//                                      Confirmed peak (Plus): Z=-0.917.
//                                      CONFIRMED no "zTPC_Proton<Charge>_AllCent"
//                                      exists (AllCent is driven-species-only).
//
//   zTPC_Deuteron<Charge>_Isolated  -- same TOF-tagging mechanism applied to the
//                                      driven species itself: a PURE, TOF-tagged
//                                      deuteron sample. Used here only as a SHAPE
//                                      SEED for the deuteron component (not fixed).
//                                      CONFIRMED: 109 entries (Plus), 1 (!) entry
//                                      (Minus) -- the Minus/antideuteron calibration
//                                      sample is statistically meaningless; expect
//                                      it to be refused by the min-entries gate.
//
// CORRECTED FIT STRATEGY (unchanged from v2, only the numerics changed in v3)
// -----------------------------------------------------------------------------
// For each charge sign:
//   1. Fit zTPC_Proton<Charge>_Isolated ALONE with a single skew-normal (auto-seeded
//      from that histogram's own peak/RMS). This fixes the proton contamination
//      peak's shape (xi, omega, alpha).
//   2. Fit zTPC_Deuteron<Charge>_Isolated ALONE with a single skew-normal (auto-
//      seeded), used only to SEED (not constrain) the deuteron component in step 3.
//      If this fit is refused (insufficient entries) or fails, step 3 falls back to
//      auto-seeding the deuteron component directly from the AllCent histogram
//      itself, since deuterons should be visible in their own inclusive band.
//   3. Fit zTPC_Deuteron<Charge>_AllCent (the real inclusive data) with a two-
//      component skew-normal sum: the proton component's xi/omega/alpha are FIXED
//      to the step-1 calibration values (only its amplitude floats); the deuteron
//      component floats freely (auto-seeded).
//
// WHAT THE RESULTING "PROTON YIELD" ACTUALLY MEANS
// ---------------------------------------------------
// The proton number out of step 3 is "protons detected within the deuteron TPC
// dE/dx band, evaluated at deuteron kinematic acceptance (y_d, mT-m_d binning)" --
// NOT a general/independent proton yield in proton-native kinematics. Do not use it
// as a standalone proton measurement outside this specific ratio.
//
// KNOWN LIMITATIONS
// --------------------
// - The main fit (step 3) treats the calibration shape (xi, omega, alpha from step 1)
//   as EXACTLY fixed, with no uncertainty propagated from that calibration into the
//   step-3 proton yield's reported error. protonInBandYieldErr is a lower bound, not
//   the full uncertainty. A rigorous treatment would re-fit step 3 with the shape
//   parameters shifted by +/-1 sigma from step 1 and take the spread as a systematic
//   -- not implemented here.
// - Auto-seeding from peak-bin/RMS is a heuristic. For a low-statistics or multi-
//   peaked histogram it can still seed badly -- the printed "auto-guess: loc=...
//   scale=..." lines let you sanity-check every seed before trusting the fit result,
//   and the printed fit status/parameter values after each fit let you catch a repeat
//   of the v2 divergence (values wildly outside the histogram's own printed range).
// - This macro remains UNTESTED end-to-end as of this version (no ROOT available in
//   the sandbox that wrote it). Run it and report back exactly what prints.
// - Whether the etof-era ZFitter (source not available in either repo mirror) used
//   this exact fixed-shape-template strategy is NOT confirmed -- this is a
//   defensible, standard PID-fitting technique, but a reconstruction, not a
//   verified match to whatever the original analysis did.
// - Look at the diagnostic PNG canvases before trusting any yield: the proton
//   component (shape fixed, amplitude free) should visibly track a real bump in
//   the AllCent data, not just be forced to near-zero amplitude by the fit.
//
// BTOF is attempted with the analogous logic (Student-t components, matching
// ZFitter's m_useStudentTDistributionsForTOF=true convention) IF matching
// zBTOF_Deuteron<Charge>_AllCent / zBTOF_Proton<Charge>_Isolated / zBTOF_Deuteron<Charge>_Isolated
// histograms exist under any of several attempted naming patterns. If not found
// (as was the case in the v2 run -- all three MISSING for both charges), BTOF is
// skipped for that charge with an explicit message.
//
// USAGE
// -----
//   root -l -b -q 'FitProtonDeuteronRatio_OO200.C("/Users/aliggett/data/OO/yieldHistos_OOGeV_deuteron_2025_08_28.root")'

#include "TFile.h"
#include "TH3.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TDirectoryFile.h"
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Fit component shapes
// ---------------------------------------------------------------------------

Double_t SkewNormalPDF(Double_t x, Double_t xi, Double_t omega, Double_t alpha){
  Double_t z = (x - xi) / omega;
  Double_t phi = TMath::Gaus(z, 0.0, 1.0, kTRUE);
  Double_t PHI = 0.5 * (1.0 + TMath::Erf(alpha * z / TMath::Sqrt2()));
  return (2.0 / omega) * phi * PHI;
}

Double_t SkewNormalComponentOnly(Double_t *x, Double_t *par){
  // par[0..3] = Amp, xi, omega, alpha
  return par[0] * SkewNormalPDF(x[0], par[1], par[2], par[3]);
}

Double_t TwoSpeciesSkewNormal(Double_t *x, Double_t *par){
  // par[0..3] = proton: Amp, xi, omega, alpha (xi/omega/alpha FIXED by caller)
  // par[4..7] = deuteron: Amp, xi, omega, alpha (all free)
  return par[0] * SkewNormalPDF(x[0], par[1], par[2], par[3])
       + par[4] * SkewNormalPDF(x[0], par[5], par[6], par[7]);
}

Double_t StudentTPDF(Double_t x, Double_t mu, Double_t sigma, Double_t ndf){
  Double_t T = (x - mu) / sigma;
  return TMath::Student(T, ndf) / sigma;
}

Double_t StudentTComponentOnly(Double_t *x, Double_t *par){
  return par[0] * StudentTPDF(x[0], par[1], par[2], par[3]);
}

Double_t TwoSpeciesStudentT(Double_t *x, Double_t *par){
  return par[0] * StudentTPDF(x[0], par[1], par[2], par[3])
       + par[4] * StudentTPDF(x[0], par[5], par[6], par[7]);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const Double_t AUTO = -999.0; // sentinel meaning "derive this from the histogram itself"
const Double_t MIN_ENTRIES_FOR_FIT = 10.0; // below this, a 4-parameter fit is not meaningful

TH3* FindTH3(TDirectoryFile *a_dir, const std::vector<std::string> &a_candidateNames,
             std::string &a_foundName){
  for(size_t iii = 0; iii < a_candidateNames.size(); iii++){
    TObject *obj = a_dir->Get(a_candidateNames[iii].c_str());
    if(obj){
      TH3 *h3 = dynamic_cast<TH3*>(obj);
      if(h3){ a_foundName = a_candidateNames[iii]; return h3; }
    }
  }
  a_foundName = "";
  return NULL;
}

TH1D* ProjectOntoZ(TH3 *a_hist, const char* a_projName){
  std::cout << "  [" << a_hist->GetName() << "] entries=" << a_hist->GetEntries()
            << " axis titles: x=\"" << a_hist->GetXaxis()->GetTitle() << "\" y=\""
            << a_hist->GetYaxis()->GetTitle() << "\" z=\""
            << a_hist->GetZaxis()->GetTitle() << "\"" << std::endl;
  TH1D *proj = (TH1D*) a_hist->Project3D("z");
  proj->SetName(a_projName);
  return proj;
}

struct CalibFitResult {
  Bool_t converged;
  Bool_t attempted;
  Double_t amp, ampErr;
  Double_t xi, xiErr;
  Double_t omega, omegaErr;
  Double_t alpha, alphaErr;
};

// Single-component fit -- used for both the proton contamination-shape
// calibration and the deuteron shape seed. Auto-derives its own initial
// guesses from the histogram (peak bin, RMS) unless overridden, and refuses
// to fit at all below MIN_ENTRIES_FOR_FIT.
CalibFitResult FitSingleSpecies(TH1D *a_proj, const char* a_label, Bool_t a_isTPC,
                                 Double_t a_locGuessOverride, Double_t a_scaleGuessOverride,
                                 TDirectoryFile *a_outDir){
  CalibFitResult result;
  result.converged = kFALSE; result.attempted = kFALSE;
  result.amp = result.ampErr = result.xi = result.xiErr = 0.0;
  result.omega = result.omegaErr = result.alpha = result.alphaErr = 0.0;

  if(!a_proj){
    std::cout << "  [" << a_label << "] missing projection -- cannot calibrate." << std::endl;
    return result;
  }
  if(a_proj->GetEntries() < MIN_ENTRIES_FOR_FIT){
    std::cout << "  [" << a_label << "] INSUFFICIENT STATISTICS (entries=" << a_proj->GetEntries()
              << " < " << MIN_ENTRIES_FOR_FIT << ") for a 4-parameter shape fit -- refusing to fit." << std::endl;
    return result;
  }
  result.attempted = kTRUE;

  Double_t xmin = a_proj->GetXaxis()->GetXmin();
  Double_t xmax = a_proj->GetXaxis()->GetXmax();
  Double_t locGuess = (a_locGuessOverride != AUTO) ? a_locGuessOverride
                     : a_proj->GetXaxis()->GetBinCenter(a_proj->GetMaximumBin());
  Double_t scaleGuess = (a_scaleGuessOverride != AUTO) ? a_scaleGuessOverride
                       : TMath::Max(a_proj->GetRMS(), (xmax - xmin) * 0.02);

  std::cout << "  [" << a_label << "] range=[" << xmin << "," << xmax << "]"
            << "  auto-guess: loc=" << locGuess << " scale=" << scaleGuess << std::endl;

  TF1 *fit = a_isTPC
    ? new TF1(Form("calibFit_%s", a_label), SkewNormalComponentOnly, xmin, xmax, 4)
    : new TF1(Form("calibFit_%s", a_label), StudentTComponentOnly, xmin, xmax, 4);

  if(a_isTPC){
    fit->SetParameters(a_proj->GetMaximum(), locGuess, scaleGuess, 0.0);
    fit->SetParNames("Amp","xi","omega","alpha");
    fit->SetParLimits(3, -50.0, 50.0);
  } else {
    fit->SetParameters(a_proj->GetMaximum(), locGuess, scaleGuess, 5.0);
    fit->SetParNames("Amp","mu","sigma","ndf");
    fit->SetParLimits(3, 1.0, 200.0);
  }
  fit->SetParLimits(0, 0.0, 1e9);
  fit->SetParLimits(1, xmin, xmax);                 // <-- prevents the v2 runaway (loc=9962)
  fit->SetParLimits(2, (xmax - xmin) * 1e-4, (xmax - xmin) * 2.0);

  TFitResultPtr fitResult = a_proj->Fit(fit, "S Q R");
  result.converged = (fitResult.Get() != NULL) && (fitResult->Status() == 0);
  result.amp = fit->GetParameter(0); result.ampErr = fit->GetParError(0);
  result.xi = fit->GetParameter(1); result.xiErr = fit->GetParError(1);
  result.omega = fit->GetParameter(2); result.omegaErr = fit->GetParError(2);
  result.alpha = fit->GetParameter(3); result.alphaErr = fit->GetParError(3);

  std::cout << "  [" << a_label << "] calib fit status=" << (result.converged ? "OK" : "FAILED/NOT-CONVERGED")
            << "  loc=" << result.xi << "+/-" << result.xiErr
            << "  scale=" << result.omega << "+/-" << result.omegaErr
            << "  shape=" << result.alpha << "+/-" << result.alphaErr << std::endl;

  TCanvas *canv = new TCanvas(Form("canv_%s", a_label), a_label, 800, 600);
  a_proj->SetMarkerStyle(20); a_proj->SetMarkerSize(0.7);
  a_proj->SetTitle(Form("%s (Calibration/Seed Fit);Z-Value;Counts / Bin", a_label));
  a_proj->Draw("PE");
  fit->SetLineColor(kRed+1); fit->SetLineWidth(2); fit->Draw("SAME");
  a_outDir->cd(); canv->Write(); canv->SaveAs(Form("%s.png", canv->GetName()));

  return result;
}

struct MainFitResult {
  Bool_t attempted;
  Bool_t converged;
  Double_t protonInBandYield, protonInBandYieldErr;
  Double_t deuteronYield, deuteronYieldErr;
  Double_t ratio_p_over_d, ratio_p_over_d_err;
};

MainFitResult FitInclusiveBand(TH1D *a_proj, const char* a_label, Bool_t a_isTPC,
                                const CalibFitResult &a_protonShape,
                                Double_t a_deuteronLocOverride, Double_t a_deuteronScaleOverride,
                                TDirectoryFile *a_outDir){
  MainFitResult result;
  result.attempted = kFALSE; result.converged = kFALSE;
  result.protonInBandYield = result.protonInBandYieldErr = 0.0;
  result.deuteronYield = result.deuteronYieldErr = 0.0;
  result.ratio_p_over_d = result.ratio_p_over_d_err = 0.0;

  if(!a_proj || a_proj->GetEntries() < MIN_ENTRIES_FOR_FIT){
    std::cout << "  [" << a_label << "] missing or insufficient-statistics inclusive projection -- skipping." << std::endl;
    return result;
  }
  if(!a_protonShape.converged){
    std::cout << "  [" << a_label << "] SKIPPED: proton calibration shape did not converge --"
              << " refusing to fix an unreliable shape into the main fit." << std::endl;
    return result;
  }
  result.attempted = kTRUE;

  Double_t xmin = a_proj->GetXaxis()->GetXmin();
  Double_t xmax = a_proj->GetXaxis()->GetXmax();
  Double_t binWidth = a_proj->GetXaxis()->GetBinWidth(1);

  Double_t deuteronLocGuess = (a_deuteronLocOverride != AUTO) ? a_deuteronLocOverride
                             : a_proj->GetXaxis()->GetBinCenter(a_proj->GetMaximumBin());
  Double_t deuteronScaleGuess = (a_deuteronScaleOverride != AUTO) ? a_deuteronScaleOverride
                               : TMath::Max(a_proj->GetRMS(), (xmax - xmin) * 0.02);

  std::cout << "  [" << a_label << "] range=[" << xmin << "," << xmax << "]"
            << "  deuteron auto-guess: loc=" << deuteronLocGuess << " scale=" << deuteronScaleGuess
            << "  (proton shape fixed: xi=" << a_protonShape.xi << " omega=" << a_protonShape.omega
            << " alpha=" << a_protonShape.alpha << ")" << std::endl;

  TF1 *fit = a_isTPC
    ? new TF1(Form("mainFit_%s", a_label), TwoSpeciesSkewNormal, xmin, xmax, 8)
    : new TF1(Form("mainFit_%s", a_label), TwoSpeciesStudentT, xmin, xmax, 8);

  // proton component: shape fixed from calibration, amplitude free
  fit->SetParameter(0, a_proj->GetMaximum() * 0.3);
  fit->SetParLimits(0, 0.0, 1e9);
  fit->FixParameter(1, a_protonShape.xi);
  fit->FixParameter(2, a_protonShape.omega);
  fit->FixParameter(3, a_protonShape.alpha);

  // deuteron component: free, auto-seeded (or overridden), constrained to the real range
  fit->SetParameter(4, a_proj->GetMaximum());
  fit->SetParLimits(4, 0.0, 1e9);
  fit->SetParameter(5, deuteronLocGuess);
  fit->SetParLimits(5, xmin, xmax);
  fit->SetParameter(6, deuteronScaleGuess);
  fit->SetParLimits(6, (xmax - xmin) * 1e-4, (xmax - xmin) * 2.0);
  fit->SetParameter(7, a_isTPC ? 0.0 : 5.0);
  if(a_isTPC) fit->SetParLimits(7, -50.0, 50.0);
  else fit->SetParLimits(7, 1.0, 200.0);

  fit->SetParNames("protonAmp","protonShape_xi_FIXED","protonShape_omega_FIXED","protonShape_alpha_FIXED",
                    "deuteronAmp","deuteronLoc","deuteronScale","deuteronShape");

  TFitResultPtr fitResult = a_proj->Fit(fit, "S Q R");
  result.converged = (fitResult.Get() != NULL) && (fitResult->Status() == 0);

  Double_t protonAmp = fit->GetParameter(0), protonAmpErr = fit->GetParError(0);
  Double_t deuteronAmp = fit->GetParameter(4), deuteronAmpErr = fit->GetParError(4);

  result.protonInBandYield = protonAmp / binWidth;
  result.protonInBandYieldErr = protonAmpErr / binWidth; // NOTE: excludes calibration-shape uncertainty, see header
  result.deuteronYield = deuteronAmp / binWidth;
  result.deuteronYieldErr = deuteronAmpErr / binWidth;

  if(result.deuteronYield > 0){
    result.ratio_p_over_d = result.protonInBandYield / result.deuteronYield;
    Double_t relP = (result.protonInBandYield > 0) ? result.protonInBandYieldErr / result.protonInBandYield : 0.0;
    Double_t relD = result.deuteronYieldErr / result.deuteronYield;
    result.ratio_p_over_d_err = result.ratio_p_over_d * TMath::Sqrt(relP*relP + relD*relD);
  }

  std::cout << "  [" << a_label << "] main fit status=" << (result.converged ? "OK" : "FAILED/NOT-CONVERGED")
            << "  proton-in-band-yield=" << result.protonInBandYield << "+/-" << result.protonInBandYieldErr
            << "  deuteron-yield=" << result.deuteronYield << "+/-" << result.deuteronYieldErr
            << "  p/d=" << result.ratio_p_over_d << "+/-" << result.ratio_p_over_d_err << std::endl;

  TCanvas *canv = new TCanvas(Form("canv_%s", a_label), a_label, 900, 700);
  a_proj->SetMarkerStyle(20); a_proj->SetMarkerSize(0.7);
  a_proj->SetTitle(Form("%s;Z-Value;Counts / Bin", a_label));
  a_proj->Draw("PE");
  fit->SetLineColor(kBlack); fit->SetLineWidth(2); fit->Draw("SAME");

  TF1 *protonComp = a_isTPC
    ? new TF1(Form("protonComp_%s", a_label), SkewNormalComponentOnly, xmin, xmax, 4)
    : new TF1(Form("protonComp_%s", a_label), StudentTComponentOnly, xmin, xmax, 4);
  protonComp->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3));
  protonComp->SetLineColor(kAzure+2); protonComp->SetLineStyle(2); protonComp->Draw("SAME");

  TF1 *deuteronComp = a_isTPC
    ? new TF1(Form("deuteronComp_%s", a_label), SkewNormalComponentOnly, xmin, xmax, 4)
    : new TF1(Form("deuteronComp_%s", a_label), StudentTComponentOnly, xmin, xmax, 4);
  deuteronComp->SetParameters(fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(6), fit->GetParameter(7));
  deuteronComp->SetLineColor(kRed+1); deuteronComp->SetLineStyle(2); deuteronComp->Draw("SAME");

  TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg->AddEntry(a_proj, "Data (inclusive)", "PE");
  leg->AddEntry(fit, "Total fit", "L");
  leg->AddEntry(protonComp, "Proton (shape fixed)", "L");
  leg->AddEntry(deuteronComp, "Deuteron (free)", "L");
  leg->Draw();

  a_outDir->cd(); canv->Write(); canv->SaveAs(Form("%s.png", canv->GetName()));

  return result;
}

// ---------------------------------------------------------------------------
// Main macro
// ---------------------------------------------------------------------------
void FitProtonDeuteronRatio_OO200(
    std::string a_yieldFile,
    std::string a_outputFile = "ProtonDeuteronRatio_AllCent.root"
){
  std::cout << "=== FitProtonDeuteronRatio_OO200 (v3 -- auto-seeded, range-constrained) ===" << std::endl;
  std::cout << "Input file: " << a_yieldFile << std::endl;
  std::cout << "NOTE: this macro is UNTESTED end-to-end (no ROOT available in the sandbox that wrote it)." << std::endl;
  std::cout << "NOTE: all-centrality-combined ONLY (per-Cent bins in this file use a different" << std::endl;
  std::cout << "      centrality scheme than the current OO200 analysis -- confirmed by inspection)." << std::endl;
  std::cout << "NOTE: 'proton yield' below means protons detected WITHIN the deuteron TPC/BTOF" << std::endl;
  std::cout << "      band at deuteron kinematic acceptance -- see header comment before using it." << std::endl;

  TFile *inFile = TFile::Open(a_yieldFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){ std::cout << "ERROR: could not open " << a_yieldFile << std::endl; return; }

  TDirectoryFile *yieldsDir = (TDirectoryFile*) inFile->Get("yields");
  if(!yieldsDir){ std::cout << "ERROR: no 'yields' directory found." << std::endl; return; }
  TDirectoryFile *deuteronSpaceDir = (TDirectoryFile*) yieldsDir->Get("Deuteron_space");
  if(!deuteronSpaceDir){ std::cout << "ERROR: no 'yields/Deuteron_space' directory found." << std::endl; return; }

  TFile *outFile = new TFile(a_outputFile.c_str(), "RECREATE");

  TTree *tree = new TTree("ProtonDeuteronRatio", "All-centrality proton-in-deuteron-band / deuteron yield ratio");
  Char_t chargeBuf[16], detBuf[16], noteBuf[256];
  Int_t protonCalibAttempted_out, protonCalibConverged_out, deuteronCalibAttempted_out, deuteronCalibConverged_out;
  Int_t mainFitAttempted_out, mainFitConverged_out;
  Double_t protonCalibXi_out, protonCalibOmega_out, protonCalibAlpha_out;
  Double_t deuteronCalibXi_out, deuteronCalibOmega_out, deuteronCalibAlpha_out;
  Double_t protonInBandYield_out, protonInBandYieldErr_out;
  Double_t deuteronYield_out, deuteronYieldErr_out;
  Double_t ratio_out, ratioErr_out;
  tree->Branch("chargeSign", chargeBuf, "chargeSign/C");
  tree->Branch("detector", detBuf, "detector/C");
  tree->Branch("note", noteBuf, "note/C");
  tree->Branch("protonCalibAttempted", &protonCalibAttempted_out, "protonCalibAttempted/I");
  tree->Branch("protonCalibConverged", &protonCalibConverged_out, "protonCalibConverged/I");
  tree->Branch("protonCalibXi", &protonCalibXi_out, "protonCalibXi/D");
  tree->Branch("protonCalibOmega", &protonCalibOmega_out, "protonCalibOmega/D");
  tree->Branch("protonCalibAlpha", &protonCalibAlpha_out, "protonCalibAlpha/D");
  tree->Branch("deuteronCalibAttempted", &deuteronCalibAttempted_out, "deuteronCalibAttempted/I");
  tree->Branch("deuteronCalibConverged", &deuteronCalibConverged_out, "deuteronCalibConverged/I");
  tree->Branch("deuteronCalibXi", &deuteronCalibXi_out, "deuteronCalibXi/D");
  tree->Branch("deuteronCalibOmega", &deuteronCalibOmega_out, "deuteronCalibOmega/D");
  tree->Branch("deuteronCalibAlpha", &deuteronCalibAlpha_out, "deuteronCalibAlpha/D");
  tree->Branch("mainFitAttempted", &mainFitAttempted_out, "mainFitAttempted/I");
  tree->Branch("mainFitConverged", &mainFitConverged_out, "mainFitConverged/I");
  tree->Branch("protonInBandYield", &protonInBandYield_out, "protonInBandYield/D");
  tree->Branch("protonInBandYieldErr", &protonInBandYieldErr_out, "protonInBandYieldErr/D");
  tree->Branch("deuteronYield", &deuteronYield_out, "deuteronYield/D");
  tree->Branch("deuteronYieldErr", &deuteronYieldErr_out, "deuteronYieldErr/D");
  tree->Branch("ratio_p_over_d", &ratio_out, "ratio_p_over_d/D");
  tree->Branch("ratio_p_over_d_err", &ratioErr_out, "ratio_p_over_d_err/D");

  const char* charges[2] = {"Plus", "Minus"};
  const char* noteText = "proton yield = in-deuteron-band contamination at deuteron kinematics, "
                          "NOT a native proton yield; error excludes calibration-shape uncertainty";

  for(int iCharge = 0; iCharge < 2; iCharge++){
    std::string chargeStr = charges[iCharge];

    // ---------------- TPC ----------------
    std::cout << "--- TPC, charge=" << chargeStr << " ---" << std::endl;
    std::string foundName;
    CalibFitResult protonCalib, deuteronCalib;
    MainFitResult mainResult;
    protonCalib.converged = kFALSE; protonCalib.attempted = kFALSE;
    deuteronCalib.converged = kFALSE; deuteronCalib.attempted = kFALSE;
    mainResult.attempted = kFALSE; mainResult.converged = kFALSE;

    std::vector<std::string> protonIsoCandidates;
    protonIsoCandidates.push_back(Form("zTPC_Proton%s_Isolated", chargeStr.c_str()));
    TH3 *hProtonIso = FindTH3(deuteronSpaceDir, protonIsoCandidates, foundName);
    if(hProtonIso){
      std::cout << "Found: " << foundName << std::endl;
      TH1D *proj = ProjectOntoZ(hProtonIso, Form("proj_ProtonCalib_TPC_%s", chargeStr.c_str()));
      protonCalib = FitSingleSpecies(proj, Form("TPC_ProtonCalib_%s", chargeStr.c_str()), kTRUE,
                                      AUTO, AUTO, outFile);
    } else {
      std::cout << "WARNING: no zTPC_Proton" << chargeStr << "_Isolated found -- cannot calibrate proton shape." << std::endl;
    }

    std::vector<std::string> deuteronIsoCandidates;
    deuteronIsoCandidates.push_back(Form("zTPC_Deuteron%s_Isolated", chargeStr.c_str()));
    TH3 *hDeuteronIso = FindTH3(deuteronSpaceDir, deuteronIsoCandidates, foundName);
    Double_t deuteronSeedLoc = AUTO, deuteronSeedScale = AUTO;
    if(hDeuteronIso){
      std::cout << "Found: " << foundName << std::endl;
      TH1D *proj = ProjectOntoZ(hDeuteronIso, Form("proj_DeuteronSeed_TPC_%s", chargeStr.c_str()));
      deuteronCalib = FitSingleSpecies(proj, Form("TPC_DeuteronSeed_%s", chargeStr.c_str()), kTRUE,
                                        AUTO, AUTO, outFile);
      if(deuteronCalib.converged){ deuteronSeedLoc = deuteronCalib.xi; deuteronSeedScale = deuteronCalib.omega; }
    } else {
      std::cout << "WARNING: no zTPC_Deuteron" << chargeStr << "_Isolated found -- main fit will auto-seed deuteron component from the AllCent histogram itself." << std::endl;
    }

    std::vector<std::string> allCentCandidates;
    allCentCandidates.push_back(Form("zTPC_Deuteron%s_AllCent", chargeStr.c_str()));
    TH3 *hAllCent = FindTH3(deuteronSpaceDir, allCentCandidates, foundName);
    if(hAllCent){
      std::cout << "Found: " << foundName << std::endl;
      TH1D *proj = ProjectOntoZ(hAllCent, Form("proj_Inclusive_TPC_%s", chargeStr.c_str()));
      mainResult = FitInclusiveBand(proj, Form("TPC_Inclusive_%s", chargeStr.c_str()), kTRUE,
                                     protonCalib, deuteronSeedLoc, deuteronSeedScale, outFile);
    } else {
      std::cout << "WARNING: no zTPC_Deuteron" << chargeStr << "_AllCent found -- cannot run inclusive fit." << std::endl;
    }

    strncpy(chargeBuf, chargeStr.c_str(), 15); chargeBuf[15]=0;
    strncpy(detBuf, "TPC", 15); detBuf[15]=0;
    strncpy(noteBuf, noteText, 255); noteBuf[255]=0;
    protonCalibAttempted_out = protonCalib.attempted ? 1 : 0;
    protonCalibConverged_out = protonCalib.converged ? 1 : 0;
    protonCalibXi_out = protonCalib.xi; protonCalibOmega_out = protonCalib.omega; protonCalibAlpha_out = protonCalib.alpha;
    deuteronCalibAttempted_out = deuteronCalib.attempted ? 1 : 0;
    deuteronCalibConverged_out = deuteronCalib.converged ? 1 : 0;
    deuteronCalibXi_out = deuteronCalib.xi; deuteronCalibOmega_out = deuteronCalib.omega; deuteronCalibAlpha_out = deuteronCalib.alpha;
    mainFitAttempted_out = mainResult.attempted ? 1 : 0;
    mainFitConverged_out = mainResult.converged ? 1 : 0;
    protonInBandYield_out = mainResult.protonInBandYield; protonInBandYieldErr_out = mainResult.protonInBandYieldErr;
    deuteronYield_out = mainResult.deuteronYield; deuteronYieldErr_out = mainResult.deuteronYieldErr;
    ratio_out = mainResult.ratio_p_over_d; ratioErr_out = mainResult.ratio_p_over_d_err;
    tree->Fill();

    // ---------------- BTOF (opportunistic) ----------------
    std::cout << "--- BTOF, charge=" << chargeStr << " ---" << std::endl;
    CalibFitResult protonCalibBTOF, deuteronCalibBTOF;
    MainFitResult mainResultBTOF;
    protonCalibBTOF.converged = kFALSE; protonCalibBTOF.attempted = kFALSE;
    deuteronCalibBTOF.converged = kFALSE; deuteronCalibBTOF.attempted = kFALSE;
    mainResultBTOF.attempted = kFALSE; mainResultBTOF.converged = kFALSE;

    std::vector<std::string> protonIsoBTOFCandidates;
    protonIsoBTOFCandidates.push_back(Form("zBTOF_Proton%s_Isolated", chargeStr.c_str()));
    protonIsoBTOFCandidates.push_back(Form("zBTOF%s_Proton_Isolated", chargeStr.c_str()));
    TH3 *hProtonIsoBTOF = FindTH3(deuteronSpaceDir, protonIsoBTOFCandidates, foundName);

    std::vector<std::string> deuteronIsoBTOFCandidates;
    deuteronIsoBTOFCandidates.push_back(Form("zBTOF_Deuteron%s_Isolated", chargeStr.c_str()));
    deuteronIsoBTOFCandidates.push_back(Form("zBTOF%s_Deuteron_Isolated", chargeStr.c_str()));
    std::string foundNameD;
    TH3 *hDeuteronIsoBTOF = FindTH3(deuteronSpaceDir, deuteronIsoBTOFCandidates, foundNameD);

    std::vector<std::string> allCentBTOFCandidates;
    allCentBTOFCandidates.push_back(Form("zBTOF_Deuteron%s_AllCent", chargeStr.c_str()));
    allCentBTOFCandidates.push_back(Form("zBTOF%s_Deuteron_AllCent", chargeStr.c_str()));
    std::string foundNameA;
    TH3 *hAllCentBTOF = FindTH3(deuteronSpaceDir, allCentBTOFCandidates, foundNameA);

    if(!hProtonIsoBTOF || !hDeuteronIsoBTOF || !hAllCentBTOF){
      std::cout << "INFO: BTOF triplet incomplete for charge=" << chargeStr
                << " (protonIso=" << (hProtonIsoBTOF?"found":"MISSING")
                << ", deuteronIso=" << (hDeuteronIsoBTOF?"found":"MISSING")
                << ", allCent=" << (hAllCentBTOF?"found":"MISSING")
                << ") -- BTOF fit skipped, not fabricated." << std::endl;
    } else {
      std::cout << "Found: " << foundName << ", " << foundNameD << ", " << foundNameA << std::endl;
      TH1D *projP = ProjectOntoZ(hProtonIsoBTOF, Form("proj_ProtonCalib_BTOF_%s", chargeStr.c_str()));
      protonCalibBTOF = FitSingleSpecies(projP, Form("BTOF_ProtonCalib_%s", chargeStr.c_str()), kFALSE,
                                          AUTO, AUTO, outFile);
      TH1D *projD = ProjectOntoZ(hDeuteronIsoBTOF, Form("proj_DeuteronSeed_BTOF_%s", chargeStr.c_str()));
      deuteronCalibBTOF = FitSingleSpecies(projD, Form("BTOF_DeuteronSeed_%s", chargeStr.c_str()), kFALSE,
                                            AUTO, AUTO, outFile);
      Double_t seedLocB = AUTO, seedScaleB = AUTO;
      if(deuteronCalibBTOF.converged){ seedLocB = deuteronCalibBTOF.xi; seedScaleB = deuteronCalibBTOF.omega; }
      TH1D *projA = ProjectOntoZ(hAllCentBTOF, Form("proj_Inclusive_BTOF_%s", chargeStr.c_str()));
      mainResultBTOF = FitInclusiveBand(projA, Form("BTOF_Inclusive_%s", chargeStr.c_str()), kFALSE,
                                         protonCalibBTOF, seedLocB, seedScaleB, outFile);
    }

    strncpy(chargeBuf, chargeStr.c_str(), 15); chargeBuf[15]=0;
    strncpy(detBuf, "BTOF", 15); detBuf[15]=0;
    strncpy(noteBuf, noteText, 255); noteBuf[255]=0;
    protonCalibAttempted_out = protonCalibBTOF.attempted ? 1 : 0;
    protonCalibConverged_out = protonCalibBTOF.converged ? 1 : 0;
    protonCalibXi_out = protonCalibBTOF.xi; protonCalibOmega_out = protonCalibBTOF.omega; protonCalibAlpha_out = protonCalibBTOF.alpha;
    deuteronCalibAttempted_out = deuteronCalibBTOF.attempted ? 1 : 0;
    deuteronCalibConverged_out = deuteronCalibBTOF.converged ? 1 : 0;
    deuteronCalibXi_out = deuteronCalibBTOF.xi; deuteronCalibOmega_out = deuteronCalibBTOF.omega; deuteronCalibAlpha_out = deuteronCalibBTOF.alpha;
    mainFitAttempted_out = mainResultBTOF.attempted ? 1 : 0;
    mainFitConverged_out = mainResultBTOF.converged ? 1 : 0;
    protonInBandYield_out = mainResultBTOF.protonInBandYield; protonInBandYieldErr_out = mainResultBTOF.protonInBandYieldErr;
    deuteronYield_out = mainResultBTOF.deuteronYield; deuteronYieldErr_out = mainResultBTOF.deuteronYieldErr;
    ratio_out = mainResultBTOF.ratio_p_over_d; ratioErr_out = mainResultBTOF.ratio_p_over_d_err;
    tree->Fill();
  }

  outFile->cd();
  tree->Write();
  outFile->Close();
  inFile->Close();

  std::cout << "=== Done. Results written to " << a_outputFile << " ===" << std::endl;
  std::cout << "Inspect with: root -l " << a_outputFile
            << " then: ((TTree*)_file0->Get(\"ProtonDeuteronRatio\"))->Scan(\"chargeSign:detector:mainFitAttempted:mainFitConverged:protonInBandYield:deuteronYield:ratio_p_over_d:ratio_p_over_d_err\")"
            << std::endl;
}
