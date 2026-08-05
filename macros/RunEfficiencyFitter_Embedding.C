// RunEfficiencyFitter_Embedding.C -- lightflavorspectra_OO200
//
// Stage 4 of the embedding pipeline (StMuAnalysisMaker at RCF -> hadd -> RunEmbeddingBinner.C
// locally -> THIS macro -> RawSpectraModifier). Fits TPC efficiency and energy loss from
// the P23ic/SL23c embedding sample (already reprocessed by RunEmbeddingBinner.C) using
// EfficiencyFitter's already-ported fitEfficiency(doTPCEff=true)/fitEnergyLoss() -- these
// are NOT new code, they're the same class RunEfficiencyFitter.C (the existing BTOF-only
// driver in this repo) already uses for the data-driven BTOF path; this macro just calls
// the embedding-data path of the same class instead.
//
// This is a NEW file (no direct _Collider.C equivalent to port from -- the original repo's
// RunEfficiencyFitter_<energy>GeV.C macros do TPC eff/energy loss/BTOF eff/feed-down all
// together in one enormous per-energy macro; this repo keeps BTOF-only (RunEfficiencyFitter.C)
// and embedding-based (this file) separate, matching how this repo's RunRawSpectraModifier.C
// port already documents that a_efficiencyFile is meant to hold BOTH -- see that macro's
// PREREQUISITES comment). Written following the same conventions/fixes already established
// in RunEfficiencyFitter.C: setUseIntegralChiSqr(false) (real, confirmed-necessary fix --
// see that macro's own header comment for the segfault it prevents), setNoFitNumBins(),
// nCentBins=6 (matches SetCutClass.C).
//
// OUTPUT: writes a SEPARATE file from RunEfficiencyFitter.C's BTOF-only output (fitter->
// write() presumably recreates its output file, so writing directly into the same file
// would clobber the existing BTOF efficiency objects rather than adding to them -- not
// verified against EfficiencyFitter::write()'s exact TFile-open mode, so don't assume
// UPDATE is safe there). RawSpectraModifier's loadAndApplyTPCEffAndEnergyLossAndBTOFEffFile()
// only takes ONE a_efficiencyFile, so once both this and the BTOF-only output exist,
// `hadd efficiency_OO200_combined.root efficiency_OO200_BTOF.root efficiency_OO200_Embedding.root`
// combines them (their object paths don't overlap -- BTOFEfficiencyFits/... vs
// EfficiencyFits/.../EnergyLossFits/...) -- pass the hadd'd file to RunRawSpectraModifier.C.
//
// PREREQUISITES: run RunEmbeddingBinner.C first for every species/charge you want fit here.

#include "../makefile_toggles.h"

void RunEfficiencyFitter_Embedding(string a_outFileName = "efficiency_OO200_Embedding.root"){

  // TODO: fill in the reprocessed EmbeddingBinner output for each species/charge you have
  // (RunEmbeddingBinner.C's output, e.g. "embed_OO200_pip_Reprocessed.root"). Leave a
  // species/charge as "" to skip it -- loadEmbeddingFiles/fitEfficiency/fitEnergyLoss are
  // only called below if the corresponding path is non-empty.
  string embFile_pip  = "";
  string embFile_pim  = "";
  string embFile_Kp   = "";
  string embFile_Km   = "";
  string embFile_pro  = "";
  string embFile_pbar = "";

  int nCentBins = 6; // matches SetCutClass.C / RunRawSpectraModifier.C / RunEfficiencyFitter.C
  bool convert16to9CentBins = false;
  int  mtm0Rebin = 1;

  cout << "RunEfficiencyFitter_Embedding starting..." << endl;

  EfficiencyFitter* fitter = new EfficiencyFitter();
  fitter->setImageDir("./EfficiencyFitterImages_Embedding/");
  // Same fix as RunEfficiencyFitter.C (BTOF-only driver) -- EfficiencyFitter's own
  // constructor defaults m_useIntegralChiSqr to true, which segfaults on a
  // non-converging fit (fitFunctSys stays NULL, then gets Clone()'d with no check).
  // Every original per-energy macro upstream explicitly disables this.
  fitter->setUseIntegralChiSqr(false);
  fitter->setNoFitNumBins(25, 120); // TPC, BTOF -- BTOF unused here but matches upstream convention

  const int PION   = 0;
  const int KAON   = 1;
  const int PROTON = 2;

  struct SpeciesJob{ int partIndex; int charge; string file; const char* label; };
  vector<SpeciesJob> jobs = {
    {PION,   1, embFile_pip,  "pi+"},
    {PION,  -1, embFile_pim,  "pi-"},
    {KAON,   1, embFile_Kp,   "K+"},
    {KAON,  -1, embFile_Km,   "K-"},
    {PROTON, 1, embFile_pro,  "proton"},
    {PROTON,-1, embFile_pbar, "antiproton"},
  };

  //##################   LOAD EMBEDDING + FIT TPC EFFICIENCY + ENERGY LOSS   ##################
  for(unsigned int jobIndex = 0; jobIndex < jobs.size(); jobIndex++){
    if(jobs[jobIndex].file == ""){
      cout << "Skipping " << jobs[jobIndex].label << " -- no embedding file given yet (see TODO above)." << endl;
      continue;
    }
    cout << "\nLoading embedding for " << jobs[jobIndex].label << "..." << endl;
    // a_inFileNameTwo="" -- RunEmbeddingBinner.C already combines every RequestID for this
    // species/charge into one file, so there's nothing left to add a second file for here.
    fitter->loadEmbeddingFiles(jobs[jobIndex].file, "", jobs[jobIndex].partIndex, jobs[jobIndex].charge,
                                mtm0Rebin, convert16to9CentBins);

    cout << "Fitting AllCent TPC efficiency for " << jobs[jobIndex].label << "..." << endl;
    // fitEfficiency(partIndex, charge, centIndex, doTPCEff, doDataDrivenBTOF, doLegendreResiduleFits)
    // doTPCEff=true, doDataDrivenBTOF=false -- this is the embedding-data path, not the
    // BTOF-only data-driven path RunEfficiencyFitter.C uses.
    fitter->fitEfficiency(jobs[jobIndex].partIndex, jobs[jobIndex].charge, -1, true, false, false);
    // RawSpectraModifier has no "AllCent" fallback for TPC efficiency either (same issue
    // RunEfficiencyFitter.C's header comment flags for BTOF) -- copy the AllCent fit into
    // every Cent%02d slot so loadAndApplyTPCEffAndEnergyLossAndBTOFEffFile() finds what it
    // expects.
    fitter->copyAllCentTPCFitsToSpecCent(jobs[jobIndex].partIndex, jobs[jobIndex].charge);

    cout << "Fitting AllCent energy loss for " << jobs[jobIndex].label << "..." << endl;
    fitter->fitEnergyLoss(jobs[jobIndex].partIndex, jobs[jobIndex].charge, -1);
    // NOTE: fitEnergyLoss's AllCent-vs-Cent%02d naming wasn't traced as deeply as TPC
    // efficiency's here -- if RawSpectraModifier warns it can't find a Cent%02d-named
    // energy-loss object when you run RunRawSpectraModifier.C, check EfficiencyFitter.cxx's
    // fitEnergyLoss()/write() for whether it needs the same copyAllCent...ToSpecCent()
    // treatment and this macro will need a matching call added.
  }

  cout << "\nWriting output file: " << a_outFileName << endl;
  fitter->write(a_outFileName);

  cout << "RunEfficiencyFitter_Embedding done." << endl;
}
