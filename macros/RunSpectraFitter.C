// RunSpectraFitter.C -- lightflavorspectra_OO200
//
// Ported from the original repo's macros/RunSpectraFitter_Collider.C, the
// ColliderCenter-eventConfig run macro (matching this repo's O+O 200 GeV setup -- same
// reasoning as RunZFitter.C/RunRawSpectraModifier.C: SpectraFitter itself has no
// FixedTarget/ColliderCenter branching, so this only affects which macro's
// CONFIGURATION VALUES are the right starting point).
//
// WHAT CHANGED FROM THE ORIGINAL, AND WHY IT'S SAFE:
//
// The original is a 7-energy-scan macro (27/19.6/14.6/13.7/11.5/9.2/7.7 GeV -- note
// 13.7 GeV, not 17.3 GeV, and its baseDir[3] is an unpopulated "" placeholder: this
// energy slot was never actually filled in by the original author). It is ALSO, unlike
// RawSpectraModifier's run macro, mostly DEAD CODE: the vast majority of the 872-line
// file is wrapped in `/* ... */` block comments, including several entire alternative
// fitting call sequences (Legendre-constrained simultaneous fits, the
// RapidityCentralityLooper variant, all of fitPionSpectra/fitKaonSpectra/
// fitProtonSpectra/fitProtonSpectraFXT/fitProtonSpectraIndivBlastWave/
// fitProtonSpectraThermal/fitProtonDnDy/fitParticlesSimultaneouslyByRapidity/fitDnDy).
//
// This was verified two ways, since a naive grep of `fitter->` calls is NOT sufficient
// here (some calls with no per-line `//` marker turned out to be inside a `/* */` block
// opened many lines earlier): (1) a script that tracks true C++ comment/string state
// character-by-character through the whole file, and (2) direct reading of every line
// range the script flagged as live, to confirm by eye. The conclusion: the ENTIRE live
// fitting call sequence in the Collider-mode macro is:
//
//   CutClass setup, then fitter->readCuts/setBlastWaveSeeds/setOldSpectraLoadingConvention
//     /readSpectraFile x3 (pion/kaon/proton) /setImageDir
//   -> a double loop over particle/centrality/charge/rapidity-bin calling
//      fitter->fitSingleRapiditySpectra(...) TWICE per bin: once "Nominal" (Tsallis for
//      pion and kaon, BlastWave for proton) and once "Systematic" (BoseEinstein for
//      pion, MtM0Expo for kaon, Thermal for proton)
//   -> fitter->writeOutputs(...)
//
// Two things worth flagging about that live loop, preserved exactly below:
//   - A `doConstraint`-gated per-bin parameter-constraint lookup is threaded through the
//     loop (reads pre-fit TF1 constraint functions from a file keyed by energy). In the
//     original, doConstraint is HARDCODED FALSE and the per-energy constraint file names
//     are Matthew-specific absolute paths that are never opened -- this branch is
//     provably always a no-op (confirmed: constraintFile stays NULL, paramFunctions all
//     stay NULL, numLoadedConstraintFunctions stays 0). It's dropped here rather than
//     carried forward as inert always-false plumbing, matching the same "genuinely dead,
//     don't carry it" precedent used for RawSpectraModifier's 27-GeV trim table and this
//     file's own huge unused historical TF1-constraint-object block (~120 lines building
//     TF1s like p_T_7p7_Plus/pro_T_thermal_19p6/etc. that are ONLY read by the dead
//     doConstrain loop above -- also not carried forward).
//   - The `if(partIndex==0){...}if(partIndex==1){...}else if(partIndex==2){...}` chain
//     (missing "else" before the partIndex==1 check) is preserved verbatim below even
//     though it looks like a typo, per the physics-fidelity mandate -- it is harmless
//     here since partIndex only ever ranges over 0/1/2 and each branch unconditionally
//     overwrites functName_Nom/functName_Sys, so the missing "else" changes nothing
//     observable, but the exact structure is kept anyway rather than "fixed."
//
// specificParticleIndex/specificCentralityIndex/specificCharge (hardcoded to
// proton-only/all-centralities/both-charges in this particular historical run of the
// macro) are exposed as function arguments here instead, defaulting to "all" (-1/-1/0)
// since there's no reason a fresh O+O 200 GeV analysis should start out proton-restricted.
//
// PREREQUISITE: a_spectraFile is RawSpectraModifier's writeRapiditySpectra(...) output
// (see RunRawSpectraModifier.C).
//
// USAGE (run from within macros/, loading LoadSpectraFitterLibs.C first):
//   root -l -q -b LoadSpectraFitterLibs.C 'RunSpectraFitter.C("spectra_modified.root","fit_output.root")'

#include "../makefile_toggles.h"
#include "../headers/SpectraFitter.h"

void RunSpectraFitter(string a_spectraFile,   // RawSpectraModifier's writeRapiditySpectra(...) output
                       string a_outputFile,   // where the fit results get written
                       string a_imageDir = "./spectrafit_images",  // directory for diagnostic fit-quality plots
                       int    a_numCent = 6, // matches CutClass's 6-bin scheme (0-5,5-10,10-20,20-40,40-80,80-100) set in SetCutClass.C -- bumped from 5 on 2026-07-09 when the 80-100% bin was added (see README's "Centrality binning" section)
                       bool   a_corrected = true,
                       bool   a_doInclusiveFits = false, // false does PRIMORDIAL fits (original repo's own comment)
                       TString a_TPC_TOF = "TPC_BTOF",
                       bool   a_oldSpectraLoadingConvention = false,
                       double a_blastWaveTempSeed  = 0.1,  // temp_kin seed, verbatim from the original (all 7 energies used the same seed)
                       double a_blastWaveBetaSSeed = 0.7,  // beta_S seed, verbatim from the original
                       int    a_specificParticleIndex   = -1, // -1 = pion+kaon+proton (0/1/2 to restrict to one)
                       int    a_specificCentralityIndex = -1, // -1 = all centrality bins
                       int    a_specificCharge          =  0){ // 0 = both +1 and -1

  CutClass *cuts = new CutClass();
  SetCutClass(cuts);

  SpectraFitter* fitter = new SpectraFitter();
  fitter->readCuts(cuts);
  fitter->setBlastWaveSeeds(a_blastWaveTempSeed, a_blastWaveBetaSSeed);
  fitter->setOldSpectraLoadingConvention(a_oldSpectraLoadingConvention);

  //########## DO NOT CHANGE THE ORDER OF THESE FUNCTIONS!!! ###############
  fitter->readSpectraFile(a_spectraFile, 0, a_TPC_TOF, a_numCent, a_corrected, a_doInclusiveFits); // pion
  fitter->readSpectraFile(a_spectraFile, 1, a_TPC_TOF, a_numCent, a_corrected, a_doInclusiveFits); // kaon
  fitter->readSpectraFile(a_spectraFile, 2, a_TPC_TOF, a_numCent, a_corrected, a_doInclusiveFits); // proton

  //must come after reading to know which dir to make
  fitter->setImageDir(a_imageDir);

  // Fit ranges: fixed verbatim -- identical across all 7 original energies.
  double lowMtM0FitRange_Nominal[3][2]    = {{0.0,0.0},{0.0,0.0},{0.0,0.0}};
  double lowMtM0FitRange_Systematic[3][2] = {{0.0,0.0},{0.0,0.0},{0.0,0.0}};
  double highMtM0FitRange_Nominal[3][2]    = {{5.0,5.0},{5.0,5.0},{5.0,5.0}};
  double highMtM0FitRange_Systematic[3][2] = {{5.0,5.0},{5.0,5.0},{5.0,5.0}};
  // "HIGH FITTING" delta-range block -- the live one of three alternatives in the
  // original (the "LOW FITTING" and "FULL FITTING" blocks were both dead/commented).
  double lowMtM0DeltaFitRange_Nominal[3][2]     = {{-1,-1},{-1,-1},{-1,-1}};
  double lowMtM0DeltaFitRange_Systematic[3][2]  = {{-1,-1},{-1,-1},{-1,-1}};
  double highMtM0DeltaFitRange_Nominal[3][2]    = {{0.6,0.6},{0.6,0.6},{0.6,0.6}};
  double highMtM0DeltaFitRange_Systematic[3][2] = {{0.6,0.6},{0.6,0.6},{0.6,0.6}};

  //STEP: core per-particle/centrality/charge/rapidity-bin fitting loop. This is the
  // ONLY live fitting call sequence in the original Collider-mode macro (see the header
  // comment above) -- fitSingleRapiditySpectra_RapidityCentralityLooper,
  // fitSpectra_Simultaneous_SingleCharge_Legendre_Constraints, fitPionSpectra,
  // fitKaonSpectra, fitProtonSpectra, fitProtonSpectraFXT, fitProtonSpectraIndivBlastWave,
  // fitProtonSpectraThermal, fitProtonDnDy, fitParticlesSimultaneouslyByRapidity, and
  // fitDnDy were all confirmed dead (block-commented) in the source macro and are not
  // called here either.
  for(int partIndex = 0; partIndex < 3; partIndex++){
    if(a_specificParticleIndex >= 0 && a_specificParticleIndex != partIndex) continue;
    for(int centIndex = 0; centIndex < a_numCent; centIndex++){
      if(a_specificCentralityIndex >= 0 && a_specificCentralityIndex != centIndex) continue;
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        int charge = 1;
        if(pmIndex == 1) charge = -1;
        if(a_specificCharge != 0 && charge != a_specificCharge) continue;
        for(unsigned int yIndex = 0; yIndex < fitter->getNRapidityBins(partIndex,0); yIndex++){

          bool   constrained[8]       = {false,false,false,false,false,false,false,false};
          double constrained_value[8] = {0,0,0,0,0,0,0,0};
          bool   constrained_sys[8]       = {false,false,false,false,false,false,false,false};
          double constrained_value_sys[8] = {0,0,0,0,0,0,0,0};

          string functName_Nom;
          string functName_Sys;
          // Preserved verbatim, including the missing "else" before the partIndex==1
          // check -- harmless here since partIndex only ranges over 0/1/2 and each
          // branch unconditionally sets both strings (see header comment above).
          if(partIndex == 0){
            functName_Nom = "Tsallis";
            functName_Sys = "BoseEinstein";
          }if(partIndex == 1){
            functName_Nom = "Tsallis";
            functName_Sys = "MtM0Expo";
          }else if(partIndex == 2){
            functName_Nom = "BlastWave";
            functName_Sys = "Thermal";
          }

          fitter->fitSingleRapiditySpectra(partIndex, charge, centIndex, yIndex, functName_Nom, true,  0, 10, 5, NULL, NULL, NULL, constrained,     constrained_value,     lowMtM0DeltaFitRange_Nominal[partIndex][pmIndex],    highMtM0DeltaFitRange_Nominal[partIndex][pmIndex]);
          fitter->fitSingleRapiditySpectra(partIndex, charge, centIndex, yIndex, functName_Sys, false, 0, 10, 5, NULL, NULL, NULL, constrained_sys, constrained_value_sys, lowMtM0DeltaFitRange_Systematic[partIndex][pmIndex], highMtM0DeltaFitRange_Systematic[partIndex][pmIndex]);

        }//rapidity loop
      }//pmIndex
    }//centIndex
  }//particle loop

  fitter->writeOutputs(a_outputFile, "RECREATE");

  cout << "-------------  End of SpectraFitter Macro  -------------" << endl;
  cout << "   Please search for ERROR and WARNING to make sure the code completed successfully " << endl;
  return;
}
