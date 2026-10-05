// ExtractDoubleDifferentialCrossSections_OO200.C -- lightflavorspectra_OO200
//
// PURPOSE: converts this analysis's per-rapidity-bin invariant-yield spectra fits
// (BlastWave for protons/antiprotons, Tsallis for pions/kaons -- see
// macros/RunSpectraFitter.C's functName_Nom assignment) into lab-frame double-
// differential production cross sections d^2(sigma)/dE dOmega(E, theta), the
// quantity a radiation-transport code (HZETRN, FLUKA, GEANT4) actually consumes as
// input, rather than the (pT or mT, y) form heavy-ion spectra are usually quoted in.
//
// THIS IS DELIBERATELY SCOPED TO PIONS, KAONS, PROTONS, AND ANTIPROTONS ONLY -- the
// species with a measured pT/mT shape at each rapidity bin to convert. It does NOT
// attempt a double-differential for the net-proton participant/spectator split
// (that's ExtractDoseFractionSummary_OO200.C's job, at the bulk/integrated level
// only) -- there is no measured pT shape at the forward rapidities the spectator
// component lives at, so a genuine (E,theta) differential for that component would
// need a physically distinct model (e.g. Fermi-motion pT smearing of a beam-
// rapidity nucleon, using a local Fermi-gas estimate off the same Woods-Saxon
// density already used in EstimateNpartNspecFractions_OO200.py:
// p_F(r) = hbar*c*(3*pi^2*rho(r))^(1/3) -- not yet implemented anywhere in this
// repo).
//
// THE (mT-m0, y) -> (E, theta) TRANSFORM:
//   Let x = mT - m0 (this analysis's fit variable), mT = x + m0, pT = sqrt(mT^2-m0^2).
//   At fixed rapidity y: pz = mT*sinh(y), E = mT*cosh(y), p = sqrt(pT^2+pz^2).
//   Lab polar angle from the beam axis: theta = atan2(pT, pz).
//   The fit functions here are already in the standard STAR "invariant yield"
//   convention: f(x,y) = E*d^3N/dp^3 = (1/(2*pi*mT)) * d^2N/(dmT dy) -- confirmed
//   directly from source for BlastWave (blastWaveModelFit_NotInvariant returns
//   2*pi*(x+par[5])*<invariant integrand>, i.e. NotInvariant = 2*pi*mT*invariant,
//   and dNdy_total_fit = Integral(NotInvariant) reproduces dN/dy exactly -- see
//   SpectraFitter.cxx's fitSingleRapiditySpectra). The other five fit shapes
//   (BoseEinstein/PtExpo/Tsallis/MtM0Expo/Boltzmann/Thermal) are fit and integrated
//   by the SAME shared code path with the same NotInvariant=2*pi*mT*invariant
//   relationship assumed generically (Par0IsDnDy is the only per-function branch);
//   this was verified directly only for BlastWave's own source this session -- if
//   PhysMath::getTsallisFunction() (used for pions/kaons, per RunSpectraFitter.C's
//   functName_Nom) turns out to use a different x convention, this whole macro's
//   pion/kaon output needs revisiting. Worth a quick sanity check (integrate this
//   macro's pion/kaon "NotInvariant"-equivalent construction against the existing
//   dNdy_PionPlus_Cent%02d_Nominal histograms -- they should match) before trusting
//   the pion/kaon numbers.
//   Standard change-of-variables result (d^3p = mT dmT dphi E dy, integrated over
//   phi assuming a phi-isotropic/azimuthally-averaged measurement):
//     d^2N/(dE dOmega) = p * E * f(x,y)
//   and multiplying by this centrality bin's share of the reaction cross section
//   (same sigma_bin = CENT_BIN_WIDTH_FRAC[centIndex] * sigma_reaction_mb convention
//   as ExtractDoseFractionSummary_OO200.C):
//     d^2(sigma)/(dE dOmega) = sigma_bin_mb * p * E * f(x,y)      [mb / (GeV * sr)]
//
// WHY THIS DOESN'T PRODUCE A REGULARLY-BINNED TH2D(E,theta): each rapidity bin has
// its OWN independently-fit invariant-yield shape (this analysis fits mT-m0 spectra
// separately in each of the 31 rapidity bins, |y|<1.55 width 0.1 -- there is no
// continuous y-dependence to interpolate between bins with, only 31 discrete
// shapes). Sweeping x at FIXED y traces a CURVE through (E,theta) space (theta
// depends on x too, through pz = mT sinh(y)), not a horizontal or vertical line --
// so the natural output of this transform is 31 curves per (species, charge,
// centrality), each exactly computed, NOT a resampled/interpolated uniform grid.
// Inventing a rebinning/interpolation scheme to force a uniform TH2D would risk
// silently distorting the normalization in a macro nobody has run yet to catch the
// mistake -- so this macro deliberately stops at the exact, ungridded (E, theta,
// value) points via TGraph2D, and leaves any regridding needed for a specific
// transport code's input format as a follow-up once someone can validate it against
// real numbers. See "NEXT STEPS" below.
//
// INPUT: a_fitOutputFile, same SpectraFitter::writeOutputs() ROOT file as the other
// two scripts in this directory. Reads, per (species, charge, centIndex, rapIndex):
//   ParticleYields/<Name>/dNdy_<Name>_Cent%02d_Nominal   -- for the rapidity bin's
//     actual y value (GetXaxis()->GetBinCenter(rapIndex+1); avoids needing the live
//     SpectraFitter object's getRapidity() accessor, matching the approach already
//     used in ExtractNetProtonStoppingFraction_OO200.C)
//   SpectraFits/<Name>/<Name>_<FunctionName>_Cent%02d_yIndex%02d
//     -- the invariant-yield fit TF1 itself. FunctionName = "Tsallis" for pions and
//     kaons, "BlastWave" for protons/antiprotons (macros/RunSpectraFitter.C's
//     functName_Nom). Object naming confirmed directly against
//     SpectraFitter::fitSingleRapiditySpectra() source this session:
//     fitFunct->SetName(Form("%s_%s_Cent%02d_yIndex%02d", ParticleName, FunctionName,
//     centIndex, rapIndex)) (SpectraFitter.cxx ~line 4813), and this exact TF1 is
//     what gets Clone()'d (no rename) into m_nominalFits_Isolated_Plus/Minus, which
//     is what writeOutputs() writes under SpectraFits/<Name>/.
//   Species/name list: PionPlus, PionMinus, KaonPlus, KaonMinus, ProtonPlus,
//     ProtonMinus (particleSpeciesName[] confirmed in submodule/ParticleInfo/
//     ParticleInfo/ParticleInfo.cxx: index 0=Pion, 1=Kaon, 2=Proton).
//
// OUTPUT: a_outputFile (default DoubleDifferentialCrossSections_OO200.root),
// organized as <ParticleName>/Cent%02d/ containing:
//   RapBin%02d_EThetaSigma   -- TGraph2D(E [GeV], theta [deg], d2sigma/dEdOmega
//                                [mb/(GeV*sr)]) for that one rapidity bin, exact
//                                (no interpolation).
//   Combined_EThetaSigma     -- all 31 rapidity bins' points concatenated into one
//                                TGraph2D, for a single overview plot.
//   cCombined (canvas, .png alongside the ROOT file) -- theta (x) vs E (y), points
//                                colored by log10(d2sigma/dEdOmega), drawn PCOL.
//
// PARTICLE MASSES: standard PDG values hardcoded below (pi+-: 0.13957 GeV, K+-:
// 0.493677 GeV, p/pbar: 0.938272 GeV) rather than pulled from
// ParticleInfo::GetParticleMass() to keep this macro standalone (no ParticleInfo/
// SpectraFitter class dependency, just a raw TFile read) -- cross-check against
// ParticleInfo::GetParticleMass() if this repo uses non-PDG values anywhere.
//
// USAGE:
//   root -l -b -q 'ExtractDoubleDifferentialCrossSections_OO200.C("fit_output.root")'
//
// UNTESTED -- same caveat as the other two scripts here: no ROOT in this sandbox.
// The physics transform above was derived and checked by hand (d^3p = p^2 dp dOmega
// -> dE dOmega via E dE = p dp is a standard, straightforward identity), and every
// object name was traced directly to its construction in SpectraFitter.cxx this
// session rather than guessed -- but this has not been run against real data.
//
// NEXT STEPS:
//   1. Sanity check: integrate this macro's own p*E*f(x,y) reconstruction of
//      NotInvariant-equivalent yield back over x at fixed y and confirm it
//      reproduces dNdy_<Name>_Cent%02d_Nominal's bin content at that y -- validates
//      both the Jacobian algebra AND the Tsallis x-convention assumption above in
//      one check.
//   2. If a specific transport code needs a uniformly-gridded TH2D(E,theta) rather
//      than these exact scattered points, grid/interpolate deliberately (e.g.
//      Delaunay triangulation via TGraph2D's own Interpolate(), or a physically
//      motivated interpolation IN (y, mT) SPACE across adjacent rapidity bins
//      before transforming to (E,theta) -- not by interpolating directly in
//      (E,theta), which would smear across the curved rapidity-bin loci) -- once
//      there's a real spectrum to validate the interpolation against.
//   3. Same Fermi-motion extension noted above, to eventually give the
//      spectator-associated component an actual (E,theta) differential instead of
//      the bulk-only number ExtractDoseFractionSummary_OO200.C reports.

#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"

using namespace std;

const int N_CENT_BINS = 6;
const char* CENT_LABELS[N_CENT_BINS] = {"0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%"};
const double CENT_BIN_WIDTH_FRAC[N_CENT_BINS] = {0.05, 0.05, 0.10, 0.20, 0.40, 0.20}; // matches the other two scripts

const int N_RAP_BINS = 31; // matches globalDefinitions.h nRapidityBins_CC

struct SpeciesInfo {
  const char* name;         // ParticleInfo::GetParticleName(partIndex,charge), e.g. "PionPlus"
  const char* functionName; // macros/RunSpectraFitter.C's functName_Nom for this species
  double mass;              // GeV, PDG
};

// Only the species with a measured pT/mT shape -- see header for why net-proton's
// participant/spectator pieces are excluded here.
const int N_SPECIES = 6;
const SpeciesInfo SPECIES_LIST[N_SPECIES] = {
  {"PionPlus",    "Tsallis",   0.13957},
  {"PionMinus",   "Tsallis",   0.13957},
  {"KaonPlus",    "Tsallis",   0.493677},
  {"KaonMinus",   "Tsallis",   0.493677},
  {"ProtonPlus",  "BlastWave", 0.938272},
  {"ProtonMinus", "BlastWave", 0.938272},
};

void ExtractDoubleDifferentialCrossSections_OO200(string a_fitOutputFile,
                                                    double a_sigmaReactionMb = 1173.0, // Woods-Saxon Glauber default -- see the other two scripts
                                                    double a_xMax = 3.0,               // mT-m0 upper sweep limit, GeV -- see header on trusting extrapolation beyond this
                                                    int    a_nXPoints = 150,
                                                    string a_outputFile = "DoubleDifferentialCrossSections_OO200.root",
                                                    string a_outputImageDir = "./doublediff_images"){

  gSystem->mkdir(a_outputImageDir.c_str(), kTRUE);

  TFile* inFile = TFile::Open(a_fitOutputFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){
    cerr << "ERROR: could not open input file '" << a_fitOutputFile << "'" << endl;
    return;
  }

  TFile* outFile = new TFile(a_outputFile.c_str(), "RECREATE");

  for(int s = 0; s < N_SPECIES; s++){
    const SpeciesInfo& sp = SPECIES_LIST[s];
    cout << "==================================================================" << endl;
    cout << " Species: " << sp.name << "  (fit: " << sp.functionName << ", mass " << sp.mass << " GeV)" << endl;

    outFile->mkdir(sp.name);

    for(int centIndex = 0; centIndex < N_CENT_BINS; centIndex++){

      TH1D* hDNdy = (TH1D*)inFile->Get(Form("ParticleYields/%s/dNdy_%s_Cent%02d_Nominal", sp.name, sp.name, centIndex));
      if(!hDNdy){
        cerr << "  WARNING: centIndex " << centIndex << " -- missing dNdy histogram, skipping this bin." << endl;
        continue;
      }

      double sigma_bin_mb = CENT_BIN_WIDTH_FRAC[centIndex] * a_sigmaReactionMb;

      outFile->mkdir(Form("%s/Cent%02d", sp.name, centIndex));
      outFile->cd(Form("%s/Cent%02d", sp.name, centIndex));

      vector<double> combinedE, combinedTheta, combinedSigma;
      int nBinsWithFits = 0;

      for(int rapIndex = 0; rapIndex < N_RAP_BINS; rapIndex++){

        // dN/dy histogram bin content check -- skip empty/unfit bins rather than
        // evaluating a fit that may not exist or may not have converged there.
        if(hDNdy->GetBinContent(rapIndex + 1) <= 0 && hDNdy->GetBinError(rapIndex + 1) <= 0) continue;

        double y = hDNdy->GetXaxis()->GetBinCenter(rapIndex + 1);

        TF1* fitFunct = (TF1*)inFile->Get(Form("SpectraFits/%s/%s_%s_Cent%02d_yIndex%02d",
                                                 sp.name, sp.name, sp.functionName, centIndex, rapIndex));
        if(!fitFunct){
          continue; // this rapidity bin's fit wasn't written/didn't converge -- silently skip, not an error
        }

        double xMaxThisFit = TMath::Min(a_xMax, fitFunct->GetXmax());

        vector<double> eVals, thetaVals, sigmaVals;
        eVals.reserve(a_nXPoints);
        thetaVals.reserve(a_nXPoints);
        sigmaVals.reserve(a_nXPoints);

        for(int xi = 0; xi < a_nXPoints; xi++){
          double x = xMaxThisFit * (double)xi / (double)(a_nXPoints - 1); // 0 .. xMaxThisFit
          double mT = x + sp.mass;
          double pT2 = mT*mT - sp.mass*sp.mass;
          if(pT2 < 0) continue; // guards x=0 rounding
          double pT = sqrt(pT2);
          double pz = mT * TMath::SinH(y);
          double E  = mT * TMath::CosH(y);
          double p  = sqrt(pT*pT + pz*pz);
          double thetaRad = atan2(pT, pz); // 0 = forward beam direction, pi = backward
          double thetaDeg = thetaRad * TMath::RadToDeg();

          double invariantYield = fitFunct->Eval(x); // f(x,y) = E d^3N/dp^3, per-event
          double d2sigma_dEdOmega = sigma_bin_mb * p * E * invariantYield; // mb / (GeV * sr)

          if(!TMath::Finite(d2sigma_dEdOmega)) continue;

          eVals.push_back(E);
          thetaVals.push_back(thetaDeg);
          sigmaVals.push_back(d2sigma_dEdOmega);
        }

        if(eVals.empty()) continue;

        TGraph2D* gBin = new TGraph2D(Form("RapBin%02d_EThetaSigma", rapIndex),
                                       Form("%s Cent %s y=%.2f;E [GeV];#theta [deg];d^{2}#sigma/dEd#Omega [mb/(GeV sr)]",
                                            sp.name, CENT_LABELS[centIndex], y),
                                       eVals.size(), &eVals[0], &thetaVals[0], &sigmaVals[0]);
        gBin->Write();
        nBinsWithFits++;

        for(size_t i = 0; i < eVals.size(); i++){
          combinedE.push_back(eVals[i]);
          combinedTheta.push_back(thetaVals[i]);
          combinedSigma.push_back(sigmaVals[i]);
        }
      } // rapidity loop

      if(combinedE.empty()){
        cerr << "  WARNING: centIndex " << centIndex << " -- no rapidity bins had usable fits, nothing written." << endl;
        continue;
      }

      TGraph2D* gCombined = new TGraph2D("Combined_EThetaSigma",
                                          Form("%s Cent %s -- All Rapidity Bins;E [GeV];#theta [deg];d^{2}#sigma/dEd#Omega [mb/(GeV sr)]",
                                               sp.name, CENT_LABELS[centIndex]),
                                          combinedE.size(), &combinedE[0], &combinedTheta[0], &combinedSigma[0]);
      gCombined->Write();

      cout << "  Cent " << CENT_LABELS[centIndex] << ": " << nBinsWithFits << "/" << N_RAP_BINS
           << " rapidity bins had usable fits, " << combinedE.size() << " total (E,theta) points." << endl;

      // ---- Diagnostic overview plot: theta vs E, colored by log10(d2sigma/dEdOmega) ----
      TCanvas* c = new TCanvas(Form("cCombined_%s_Cent%02d", sp.name, centIndex), "", 900, 650);
      gCombined->SetMarkerStyle(20);
      gCombined->SetMarkerSize(0.6);
      gCombined->Draw("PCOL"); // colored markers by Z value; if this draw option isn't
                                 // available in your ROOT version, "P0" (plain markers,
                                 // no color-by-Z) is the safe fallback.
      c->SetLogz();
      c->SaveAs(Form("%s/DoubleDiff_%s_Cent%02d.png", a_outputImageDir.c_str(), sp.name, centIndex));
      delete c;

    } // centrality loop
  } // species loop

  outFile->Close();
  inFile->Close();
  cout << "==================================================================" << endl;
  cout << "Wrote " << a_outputFile << " and diagnostic images to " << a_outputImageDir << endl;
  cout << "Run the sanity check in NEXT STEPS #1 before trusting these numbers." << endl;
}
