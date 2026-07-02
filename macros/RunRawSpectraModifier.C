// RunRawSpectraModifier.C -- lightflavorspectra_OO200
//
// Ported from the original repo's macros/RunRawSpectraModifier_Collider.C, the
// ColliderCenter-eventConfig run macro (matching this repo's O+O 200 GeV setup, same
// reasoning as RunZFitter.C/headers/ZFitter.h -- RawSpectraModifier itself has no
// FixedTarget/ColliderCenter branching at all, so this choice only affects which
// macro's CONFIGURATION VALUES are the right starting point to carry over).
//
// WHAT CHANGED FROM THE ORIGINAL, AND WHY IT'S SAFE:
//
// Like RunZFitter_Collider.C, this is a 7-ENERGY SCAN macro (7.7/9.2/11.5/14.6/17.3/
// 19.6/27 GeV) with per-energy arrays of hardcoded absolute file paths and several
// per-energy config values. Going through every value used by the live code path
// (i.e. excluding a large `/* ... */`-commented-out block of hand-tuned per-rapidity-
// bin trim data for the 27 GeV analysis specifically, which was already fully dead in
// the original and is not carried forward here at all):
//
//   - Values IDENTICAL across all 7 energies were kept verbatim as fixed local
//     variables: nCentBins=9, converted16CentBinsTo9CentBins=false, pTLossMode=true,
//     useTPCEffSystematic/useBTOFEffSystematic/useEnergyLossSystematic=true,
//     useTPCEffDataPoints/useBTOFEffDataPoints=true (original comment: "should be
//     false in final analysis" -- carried forward as-is, your call when to flip it),
//     useCent16ForCent9=false, overallEfficiency=1.0/overallEfficiencyErr=0.0,
//     offsets_tpc_accept=all zero, maxFracError=10.0, removeByPreviousYieldDiffFactTwo
//     =false, mTm0_tolerence=0.01, startMom_tpc={0,0,0} (always -- a no-op "start from
//     the beginning" cut), startMom_eTOF={0.6,0.5,1.0}, stopMom_eTOF={3,3,3},
//     lowEtaCutTPC/highEtaCutTPC/lowEtaCutBTOF/highEtaCutBTOF={-5,-5,-5}/{5,5,5}
//     (i.e. no eta cut), doEnergyLoss/doBTOFEff/doTPCEff=true, and the entire
//     startMtM0_matrix_bTOF table (a single hardcoded matrix used identically for
//     every energy in the original, not indexed by energy at all).
//   - Values that DO vary by energy in the original (startMom_bTOF, stopMom_bTOF,
//     lowRapCut/highRapCut, doHybridFeedDown, and of course all the file paths) have
//     no "right" answer for a brand-new O+O 200 GeV dataset that wasn't part of the
//     original scan, so they are exposed as explicit function arguments here instead
//     of silently reusing one energy's tuned numbers. Defaults below are the most
//     common values across the original 7 entries (documented inline per argument),
//     but you should sanity-check them against this dataset's actual rapidity
//     acceptance/momentum ranges before trusting the output.
//   - trimData/pointByPointData: 5 of the original's 7 energies used an EMPTY point-
//     by-point cleanup list (only 7.7 GeV had one, itself just a handful of hand-
//     picked bad points found by eye in that specific dataset) -- so defaulting to
//     empty here (matching the majority behavior) rather than guessing is the
//     faithful choice; pass your own once you've identified problem points in this
//     dataset's spectra, matching the original class's trimSpectra() interface.
//   - mod->setCuts(...)/mod->setStarVer(...) were commented out (unused) in the
//     original Collider macro too -- not resurrected here, so this macro has no
//     CutClass dependency at all (matching RawSpectraModifier.h's constructor, which
//     needs neither to run the live call sequence below).
//
// PREREQUISITES: the a_spectraFile argument is ZFitter's makeSpectra(...) output (see
// RunZFitter.C); the three a_picoBinnerFile arguments are PicoBinner's per-particle
// yield-histogram outputs (see RunPicoBinner.C); a_efficiencyFile is an embedding-
// derived TPC/BTOF-efficiency + energy-loss correction file for THIS dataset -- there
// is no such file shipped with this repo, you'll need to produce/supply your own.
//
// USAGE (run from within macros/, loading LoadRawSpectraModifierLibs.C first):
//   root -l -q -b LoadRawSpectraModifierLibs.C 'RunRawSpectraModifier.C("pionYield.root","kaonYield.root","protonYield.root","spectra.root","efficiency.root","spectra_modified.root")'

#include "../makefile_toggles.h"
#include "../headers/RawSpectraModifier.h"

void RunRawSpectraModifier(string a_pionYieldFile, string a_kaonYieldFile, string a_protonYieldFile,
                            string a_spectraFile,       // ZFitter's makeSpectra(...) output
                            string a_efficiencyFile,    // embedding-derived TPC/BTOF eff + energy-loss corrections for THIS dataset
                            string a_outputFile,        // where the modified/corrected spectra get written
                            bool   a_rawOnly = false,   // true skips all corrections -- just writes the raw rapidity spectra
                            bool   a_doHybridFeedDown = false, // true requires a HybridFeedDownMaker output -- see loadAndApplyHybridFeedDownCorrections
                            vector<double> a_startMom_bTOF = {0.55, 0.4, 0.8},  // pi/K/p -- most common value across the original's 7 energies
                            vector<double> a_stopMom_bTOF  = {4.0, 4.0, 4.0},   // pi/K/p -- most common value across the original's 7 energies
                            vector<double> a_lowRapCut  = {-1.05, -0.95, -0.95}, // pi/K/p -- widest (least restrictive) value seen in the original
                            vector<double> a_highRapCut = { 1.05,  0.95,  0.95}, // pi/K/p -- widest (least restrictive) value seen in the original
                            vector< vector<double> > a_trimData = {},
                            vector< vector<double> > a_pointByPointData = {}){

  if(a_spectraFile == a_outputFile){
    cout << " Don't overwrite your ZFitter spectra file -- a_outputFile must differ from a_spectraFile. Returning." << endl;
    return;
  }

  // ################  FIXED CONFIGURATION (identical across all 7 original entries)  ################
  int  nCentBins = 5; // matches CutClass's 5-bin scheme (0-5,5-10,10-20,20-40,40-80) set in SetCutClass.C
  bool converted16CentBinsTo9CentBins = false;
  bool pTLossMode = true;

  bool useTPCEffSystematic     = true;
  bool useBTOFEffSystematic    = true;
  bool useEnergyLossSystematic = true;
  bool useTPCEffDataPoints     = true; // original repo's own comment: "should be false in final analysis"
  bool useBTOFEffDataPoints    = true;

  bool useCent16ForCent9 = false;
  double overallEfficiency    = 1.0;
  double overallEfficiencyErr = 0.0;
  int offsets_tpc_accept[9] = {0,0,0,0,0,0,0,0,0};
  double maxFracError = 10.0;
  bool removeByPreviousYieldDiffFactTwo = false;
  double mTm0_tolerence = 0.01;

  bool doEnergyLoss = true;
  bool doTPCEff     = true;
  bool doBTOFEff    = true;

  double startMom_tpc[3]  = {0.0, 0.0, 0.0};   // no-op: start from the beginning
  double startMom_eTOF[3] = {0.6, 0.5, 1.0};
  double stopMom_eTOF[3]  = {3.0, 3.0, 3.0};
  double lowEtaCutTPC[3]   = {-5,-5,-5};
  double highEtaCutTPC[3]  = { 5, 5, 5};
  double lowEtaCutBTOF[3]  = {-5,-5,-5};
  double highEtaCutBTOF[3] = { 5, 5, 5};
  double minMtM0[3] = {0.0, 0.0, 0.0};
  double maxMtM0[3] = {5.0, 5.0, 5.0};

  // startMtM0_matrix_bTOF: single hardcoded matrix used identically for every energy
  // in the original macro (not indexed by energy there either). {part, rapidity,
  // bTofStartMtM0}. Carried over verbatim.
  vector<vector<double>> startMtM0_matrix_bTOF = {
    {0, 0.0, 0.15}, {0, 0.1, 0.15}, {0, 0.2, 0.15}, {0, 0.3, 0.15}, {0, 0.4, 0.15},
    {0, 0.5, 0.15}, {0, 0.6, 0.15}, {0, 0.7, 0.15}, {0, 0.8, 0.15}, {0, 0.9, 0.15},
    {0, 1.0, 0.20}, {0, 1.1, 0.31}, {0, 1.2, 1.0},  {0, 1.3, 1.0},  {0, 1.4, 1.0},
    {0, 1.5, 1.0},  {0, 1.6, 1.0},

    {1, 0.0, 0.13}, {1, 0.1, 0.115},{1, 0.2, 0.10}, {1, 0.3, 0.19}, {1, 0.4, 0.08},
    {1, 0.5, 0.07}, {1, 0.6, 0.07}, {1, 0.7, 0.07}, {1, 0.8, 0.06}, {1, 0.9, 0.06},
    {1, 1.0, 0.06}, {1, 1.1, 0.06},

    {2, 0.0, 0.1},  {2, 0.1, 0.1},  {2, 0.2, 0.1},  {2, 0.3, 0.1},  {2, 0.4, 0.1},
    {2, 0.5, 0.1},  {2, 0.6, 0.15}, {2, 0.7, 0.18}, {2, 0.8, 0.25}, {2, 0.9, 0.37},
    {2, 1.0, 0.65}, {2, 1.1, 1.0},  {2, 1.2, 1.0}
  };

  RawSpectraModifier* mod = new RawSpectraModifier();
  mod->setChrisCorrectionMode(false);
  mod->setPtEnergyLossMode(pTLossMode);

  //########## DO NOT CHANGE THE ORDER OF THESE FUNCTIONS!!! ###############
  mod->loadDataFile(a_pionYieldFile,   0, nCentBins, converted16CentBinsTo9CentBins);
  mod->loadDataFile(a_kaonYieldFile,   1, nCentBins, converted16CentBinsTo9CentBins);
  mod->loadDataFile(a_protonYieldFile, 2, nCentBins, converted16CentBinsTo9CentBins);

  mod->loadSpectraFile(a_spectraFile.c_str());
  mod->convertSpectraToInvariant(overallEfficiency, overallEfficiencyErr);

  //STEP 1: Make Raw Spectra TGraphs
  mod->makeRapiditySpectraTGraphs();

  //STEP 2: Apply the feeddown/knockout/TPC efficiencies (IN ORDER!) -- feed-down
  // before TPC efficiency, energy loss before TPC efficiency, hybrid feed-down last
  // (it doesn't go through StarSim, unlike the others).
  if(!a_rawOnly){
    mod->loadAndApplyTPCEffAndEnergyLossAndBTOFEffFile(a_efficiencyFile, useTPCEffSystematic, useBTOFEffSystematic, useEnergyLossSystematic,
             offsets_tpc_accept, useCent16ForCent9, useTPCEffDataPoints, useBTOFEffDataPoints, doEnergyLoss, doTPCEff, doBTOFEff);
    // a_momStop_tpc is deliberately handed a_startMom_bTOF's values -- the TPC
    // spectrum's upper momentum bound IS the BTOF spectrum's starting momentum (the
    // TPC-to-BTOF handoff point), matching the original macro's own call exactly.
    mod->cleanSpectra(startMom_tpc, &a_startMom_bTOF[0], startMom_eTOF, &a_startMom_bTOF[0], &a_stopMom_bTOF[0], stopMom_eTOF,
                    lowEtaCutTPC, highEtaCutTPC, lowEtaCutBTOF, highEtaCutBTOF,
                    minMtM0, maxMtM0, maxFracError, &a_lowRapCut[0], &a_highRapCut[0], removeByPreviousYieldDiffFactTwo, startMtM0_matrix_bTOF);
    mod->trimSpectra(a_trimData, a_pointByPointData, mTm0_tolerence);
    mod->populateInclusiveSpectra(); // corrected spectra without feed-down corrections
    if(a_doHybridFeedDown) mod->loadAndApplyHybridFeedDownCorrections(a_efficiencyFile, false); // bool is use-point-eval instead of integral
  }

  mod->writeMemberData(a_outputFile.c_str(), "RECREATE");
  mod->writeRapiditySpectra(a_outputFile.c_str(), "UPDATE");

  cout << "-------------  End of RawSpectraModifier Macro  -------------" << endl;
  cout << "   Please search for ERROR and WARNING to make sure the code completed successfully " << endl;
  return;
}
