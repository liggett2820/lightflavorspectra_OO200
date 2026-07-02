// ConvertSpectraToCrossSection.C -- lightflavorspectra_OO200
//
// Converts SpectraFitter's per-event invariant-yield spectra (TotalSpectra/<Particle>/
// totalSpectra_<Particle>_Cent_<CC>_yIndex<YY>, TGraphErrors, y-axis =
// (1/N_evt)(1/2*pi*m_T) d^2N/dm_T dy -- see macros/RunSpectraFitter.C and
// source/SpectraFitter.cxx) into an actual cross section, in both the standard
// conventions:
//
//   Invariant:     E d^3(sigma)/dp^3        = (1/2*pi*m_T) d^2(sigma)/dm_T dy   [mb GeV^-2]
//   Non-invariant: d^2(sigma)/dp_T dy                                          [mb GeV^-1]
//
// (These two differ only by a factor of 2*pi*p_T -- the "invariant yield per event"
// stored by SpectraFitter is already invariant under whether you differentiate wrt
// p_T or m_T, since dm_T/dp_T = p_T/m_T cancels the difference; see e.g. any PDG
// kinematics reference. p_T at each point is reconstructed from the graph's x-axis
// (m_T-m_0) plus the particle's PDG mass via m_T = x + m_0, p_T = sqrt(m_T^2 - m_0^2).)
//
// THE CROSS SECTION NORMALIZATION (per centrality bin): the stored quantity is
// "per event *in that centrality bin*", so converting a single bin to a cross section
// means multiplying by that bin's SHARE of the total reaction cross section, i.e.
//   sigma_bin = (percentile width of the bin) * sigma_reaction_total
// e.g. the 0-5% (most central) bin gets 5% of sigma_reaction_total, the 10-20% bin
// gets 10%, etc. The percentile widths used here are copied verbatim from
// source/CutClass.cxx's setStRefMultCorrCentralities() (see its own header comment,
// "Centrality binning:", around line 400) -- NOT re-derived or assumed -- so they are
// guaranteed to match this repo's own centIndex convention (centIndex 0 = 0-5% most
// central ... centIndex 8 = 70-80% most peripheral for the 9-bin scheme this repo
// defaults to; see CutClass::centralityIndex(), "centIndex = 8 - getCentralityBin9()").
//
// sigma_reaction_total (mb, with its uncertainty) is NOT computed here -- it comes from
// a separate Glauber Monte Carlo calculation (this repo has no Glauber code of its
// own). Use ../../glaubermchardness/macros/BootstrapReactionCrossSection.C to produce
// it (see that macro's header comment for why GlauberModel.cxx's own printed cross
// section isn't reusable directly, and how this one re-derives + bootstraps it) --
// point a_sigmaReactionFile at that macro's output file.
//
// USAGE (run from within macros/, loading LoadDrawingLibs.C first -- same dependency,
// just ParticleInfo, as DrawSpectraAndDnDy.C):
//   root -l -q -b LoadDrawingLibs.C 'ConvertSpectraToCrossSection.C("fit_output.root","../../glaubermchardness/data/sigmaReaction_OO200.root")'
//   root -l -q -b LoadDrawingLibs.C 'ConvertSpectraToCrossSection.C("fit_output.root","sigmaReaction.root","crossSections.root",1)'
//
// a_centScheme: 0 (default) = the 5-bin scheme (0-5,5-10,10-20,20-40,40-80%) SetCutClass.C
// currently configures; 1 = the old 9-bin scheme (0-5,5-10,10-20,20-30,...,70-80%); 2 =
// the old 16-bin scheme (all 5%-wide). Must match whatever a_numCent RunSpectraFitter.C
// was actually run with -- passing the wrong one will silently apply the wrong bin widths.

#include "../makefile_toggles.h"
// Needed so cling knows this type when the macro is interpreted -- gSystem->Load() in
// LoadDrawingLibs.C loads the compiled implementation, but doesn't by itself guarantee
// the class declaration is known to the interpreter (same reasoning as RunPicoBinner.C).
#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TParameter.h"
#include "TMath.h"
#include "TSystem.h"
#include <vector>
#include <string>

// Fixed array bounds copied from headers/SpectraFitter.h -- see DrawSpectraAndDnDy.C's
// header comment for the same reasoning; kept independent here since this macro doesn't
// load SpectraFitter.h either.
static const int g_MaxCent    = 16;
static const int g_MaxRapBins = 35;

void ConvertSpectraToCrossSection(string  a_fitOutputFile,                    // SpectraFitter's writeOutputs() output, e.g. "fit_output.root"
                                   string  a_sigmaReactionFile,                // BootstrapReactionCrossSection.C's output (see header comment)
                                   string  a_outputFile     = "crossSections.root",
                                   int     a_centScheme     = 0){              // 0=5-bin (current default, see header comment), 1=old 9-bin, 2=old 16-bin

  // ---- Centrality percentile-width table (fractional, e.g. 0.05 = 5%) ----
  // centEdges5 matches SetCutClass.C's current setCentralities(5,...) call exactly
  // (0-5,5-10,10-20,20-40,40-80%). centEdges9/centEdges16 are copied verbatim from
  // source/CutClass.cxx's setStRefMultCorrCentralities() m_centPercents (see that
  // function's header comment for the full bin-index-to-% table) and kept here for
  // anyone still working with older 9-bin/16-bin output. centIndex i covers
  // [edges[i], edges[i+1])% most-central.
  std::vector<double> centEdges5  = {0,5,10,20,40,80};
  std::vector<double> centEdges9  = {0,5,10,20,30,40,50,60,70,80};
  std::vector<double> centEdges16 = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80};
  std::vector<double>& centEdges = (a_centScheme == 2) ? centEdges16 : (a_centScheme == 1) ? centEdges9 : centEdges5;
  int nCentBinsInScheme = (int) centEdges.size() - 1;

  // ---- Load sigma_reaction (mb) + its uncertainty ----
  TFile* sigmaFile = TFile::Open(a_sigmaReactionFile.c_str(), "READ");
  if(!sigmaFile || sigmaFile->IsZombie()){
    cout << "ERROR: could not open a_sigmaReactionFile " << a_sigmaReactionFile << endl;
    return;
  }
  TParameter<Double_t>* sigmaParam    = (TParameter<Double_t>*) sigmaFile->Get("sigmaReaction_mb");
  TParameter<Double_t>* sigmaErrParam = (TParameter<Double_t>*) sigmaFile->Get("sigmaReaction_mb_err");
  if(!sigmaParam){
    cout << "ERROR: \"sigmaReaction_mb\" not found in " << a_sigmaReactionFile
         << " -- is this really BootstrapReactionCrossSection.C's output file?" << endl;
    return;
  }
  double sigmaReaction_mb     = sigmaParam->GetVal();
  double sigmaReaction_mb_err = sigmaErrParam ? sigmaErrParam->GetVal() : 0.0;
  double sigmaReaction_relErr = (sigmaReaction_mb > 0) ? sigmaReaction_mb_err / sigmaReaction_mb : 0.0;
  cout << "Using sigma_reaction = " << sigmaReaction_mb << " +/- " << sigmaReaction_mb_err
       << " mb (" << 100.0*sigmaReaction_relErr << "% relative)" << endl;
  sigmaFile->Close();

  // ---- Open input spectra / output file ----
  TFile* inFile = TFile::Open(a_fitOutputFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){
    cout << "ERROR: could not open a_fitOutputFile " << a_fitOutputFile << endl;
    return;
  }
  TDirectory* specTopDir = inFile->GetDirectory("TotalSpectra");
  if(!specTopDir){
    cout << "ERROR: no \"TotalSpectra\" directory in " << a_fitOutputFile
         << " -- is this really SpectraFitter's writeOutputs() output?" << endl;
    return;
  }

  TFile* outFile = new TFile(a_outputFile.c_str(), "RECREATE");
  outFile->mkdir("CrossSection_Invariant");
  outFile->mkdir("CrossSection_NonInvariant");

  ParticleInfo* particleInfo = new ParticleInfo(); // default ctor -- matches SpectraFitter's own m_partInfo, so names match exactly

  int nConverted = 0;

  for(int partIndex = 0; partIndex < 3; partIndex++){ // pion/kaon/proton -- only species RunSpectraFitter.C ever fits
    double mass = particleInfo->GetParticleMass(partIndex); // GeV

    for(int chargeSign = 0; chargeSign < 2; chargeSign++){
      int charge = (chargeSign == 0) ? 1 : -1;
      TString particleName = particleInfo->GetParticleName(partIndex, charge); // e.g. "PionPlus"

      TDirectory* particleDir = specTopDir->GetDirectory(particleName.Data());
      if(!particleDir) continue;

      outFile->cd("CrossSection_Invariant");
      gDirectory->mkdir(particleName.Data());
      outFile->cd("CrossSection_NonInvariant");
      gDirectory->mkdir(particleName.Data());

      for(int centIndex = 0; centIndex < g_MaxCent; centIndex++){

        if(centIndex >= nCentBinsInScheme) break; // beyond the chosen scheme's bin count

        double binWidthFrac = (centEdges[centIndex+1] - centEdges[centIndex]) / 100.0; // e.g. 0.05
        double sigmaBin_mb     = binWidthFrac * sigmaReaction_mb;
        double sigmaBin_mb_err = binWidthFrac * sigmaReaction_mb_err; // width is exact/fixed, so relative error carries through unchanged

        for(int yIndex = 0; yIndex < g_MaxRapBins; yIndex++){
          TString objName = Form("totalSpectra_%s_Cent_%02d_yIndex%02d", particleName.Data(), centIndex, yIndex);
          TGraphErrors* g = (TGraphErrors*) particleDir->Get(objName.Data());
          if(!g || g->GetN() == 0) continue;

          TGraphErrors* gInv    = (TGraphErrors*) g->Clone(Form("crossSection_Invariant_%s_Cent_%02d_yIndex%02d",    particleName.Data(), centIndex, yIndex));
          TGraphErrors* gNonInv = (TGraphErrors*) g->Clone(Form("crossSection_NonInvariant_%s_Cent_%02d_yIndex%02d", particleName.Data(), centIndex, yIndex));
          gInv->SetTitle(Form("Invariant Cross Section %s Cent%02d yIndex%02d;m_{T}-m_{0} (GeV/c^{2});E d^{3}#sigma/dp^{3} (mb GeV^{-2})", particleName.Data(), centIndex, yIndex));
          gNonInv->SetTitle(Form("d^{2}#sigma/dp_{T}dy %s Cent%02d yIndex%02d;p_{T} (GeV/c);d^{2}#sigma/dp_{T}dy (mb GeV^{-1})", particleName.Data(), centIndex, yIndex));

          for(int p = 0; p < g->GetN(); p++){
            double mTm0  = g->GetX()[p];
            double mT    = mTm0 + mass;
            double pT    = TMath::Sqrt(TMath::Max(0.0, mT*mT - mass*mass));
            double yVal  = g->GetY()[p];
            double yErr  = g->GetErrorY(p);

            // Invariant: E d3sigma/dp3 = (invariant yield per event) * sigma_bin
            double invY    = yVal * sigmaBin_mb;
            double invYErr = TMath::Sqrt(TMath::Power(yErr*sigmaBin_mb, 2) + TMath::Power(yVal*sigmaBin_mb_err, 2));
            gInv->SetPoint(p, mTm0, invY);
            gInv->SetPointError(p, g->GetErrorX(p), invYErr);

            // Non-invariant: d2sigma/dpT dy = 2*pi*pT * (invariant cross section)
            double jac        = 2.0 * TMath::Pi() * pT;
            double nonInvY     = invY * jac;
            double nonInvYErr  = invYErr * jac;
            gNonInv->SetPoint(p, pT, nonInvY);
            gNonInv->SetPointError(p, g->GetErrorX(p), nonInvYErr); // x-error carried through unchanged (dpT/dmT ~ 1 near midrapidity; not re-derived here)
          }

          outFile->cd(Form("CrossSection_Invariant/%s", particleName.Data()));
          gInv->Write();
          outFile->cd(Form("CrossSection_NonInvariant/%s", particleName.Data()));
          gNonInv->Write();
          nConverted++;
        } // yIndex
      } // centIndex
    } // chargeSign
  } // partIndex

  outFile->Close();
  inFile->Close();

  cout << "-------------  End of ConvertSpectraToCrossSection Macro  -------------" << endl;
  cout << "Converted " << nConverted << " spectra to cross sections (both conventions) -> " << a_outputFile << endl;
  if(nConverted == 0) cout << "WARNING: nothing converted -- check a_fitOutputFile actually contains TotalSpectra/ objects." << endl;
}
