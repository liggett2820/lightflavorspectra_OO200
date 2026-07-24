// RunEfficiencyFitter.C -- lightflavorspectra_OO200
//
// Ported from the original repo's macros/RunEfficiencyFitter_Collider.C, TRIMMED TO
// BTOF-EFFICIENCY-ONLY (the original's "a_fitMode==4 / B" path). See
// claude/OO200_analysis_log.md for the full story of why this exists: RawSpectraModifier
// ::loadAndApplyTPCEffAndEnergyLossAndBTOFEffFile() needs an externally-produced
// efficiency file containing BTOFEfficiencyGraphs/btofFromDataEfficiencyGraph_... and
// BTOFEfficiencyFits/btofEfficiencyFit_DataDriven_... objects, and this "Efficiency
// Fitter" tool (EfficiencyFitter class, ported verbatim from lightflavorspectra_etof)
// is what produces them. Its BTOF-only mode is entirely PicoBinner-yield-file-driven
// (the SAME yieldHistos_OO200_*.root files RunRawSpectraModifier.C's a_pionYieldFile/
// a_kaonYieldFile args and RunSpectraFitter.C already use) -- it does NOT need an
// embedding/simulation sample, unlike the TPC-efficiency and energy-loss pieces of the
// original tool, which is why only the BTOF path is ported here.
//
// WHAT'S DELIBERATELY NOT PORTED (out of scope for now, not silently dropped):
//   - TPC efficiency fitting (fitMode T/TBF/TBFE/etc.) -- needs an embedding sample,
//     which doesn't exist for this dataset yet.
//   - Energy-loss fitting (fitMode E) -- same reason.
//   - Feed-down fitting (fitMode F) -- same reason.
//   - Knockout background (fitMode K) -- a SEPARATE fit mode from plain BTOF
//     efficiency (uses the ProtonDCA histograms loadPicoBinnerDataFile() also reads
//     for proton, but never calls the knockout-specific fit function). Not wired
//     here regardless of whether a proton file is supplied -- if/when knockout
//     background fitting is needed, port that piece the same way as the others.
// If/when an O+O 200 embedding sample exists, those pieces can be ported the same way
// this one was: read the *_Collider.C block for that fit mode, confirm its actual repo-
// specific dependencies, port verbatim, wire into makeLibs.C.
//
// WHAT CHANGED FROM THE ORIGINAL, AND WHY:
//   - a_numCentBins defaults to 6, NOT the original's numCentBins_picoBinner=9 --
//     matches THIS repo's SetCutClass.C (cuts->setCentralities(6,...), see
//     macros/SetCutClass.C) and RunRawSpectraModifier.C's own nCentBins=6. Confirmed
//     via source/PicoBinner.cxx:178 (`int numCentBins = a_cutClass->getNCentBins();`)
//     that PicoBinner's raw output bin count comes from CutClass, not a fixed
//     constant -- so this is NOT a value that should be blindly copied from the
//     etof scan, and isn't.
//   - a_convert16to9CentBins defaults to false -- PicoBinner's output here is
//     natively 6 bins (via CutClass), never 16, so no down-conversion applies. (The
//     original's convert16to9_picoBinner flag exists for OTHER energies in that scan
//     that ran with a 16-bin CutClass config and wanted 9 bins out -- doesn't apply
//     to this repo's fixed 6-bin config at all.)
//   - a_fitAllCentBTOF / a_useAllCentBTOFFits are exposed as explicit arguments
//     rather than hardcoded, because this is a genuine statistics/physics call that
//     depends on how much data lands in each of this (small-system O+O) dataset's
//     centrality bins -- see the long comment above their declaration below.
//   - Proton is OPTIONAL, gated on a_protonPicoBinnerFile != "" (2026-07-20: the
//     current proton PicoBinner yield file, yieldHistos_OOGeV_proton_2026_07_01.root,
//     is confirmed corrupt at Cent5 -- see project log -- so the default is "",
//     skipping proton entirely, same convention RunRawSpectraModifier.C already uses
//     for an optional species file. Pass a real proton file once one exists and
//     pion/kaon/proton all get fit together in the same run/output file.)
//   - No per-energy file-path arrays: just plain file-path arguments (pion, kaon,
//     optionally proton).
//   - FIX (2026-07-21): the port originally never called fitter->setUseIntegralChiSqr(),
//     leaving EfficiencyFitter's own constructor default (m_useIntegralChiSqr=true) in
//     effect -- this silently enables a genuinely buggy code path (confirmed by reading
//     source/EfficiencyFitter.cxx directly): when fitEfficiency()'s GraphFitter-based
//     integral-chi-square fit fails to converge (GraphFitter::minimize() status!=0),
//     fitFunctSys stays NULL and a "WARNING: BAD FIT...." line prints, but the very next
//     block (a_doDataDrivenBTOF && allCentMode) unconditionally calls
//     fitFunctSys->Clone() with no NULL check -- a straight segfault on any non-
//     converging fit. Confirmed EVERY original per-energy macro in the etof repo
//     (RunEfficiencyFitter_Collider.C and all the single-energy variants) explicitly
//     sets this to false (their shared "useIntegralGraphFitter = NONONONONONONONONONO"
//     convention, NONONONONONONONONONO being a locally-declared name for `false`) --
//     i.e. this bug is real and universally avoided upstream, not something specific to
//     this port. The line below fixes that gap. (A previous comment on a_doLegendreBTOFFits
//     below incorrectly implied IT was the relevant "buggy, avoid" toggle -- it isn't;
//     that parameter controls an unrelated, separate Legendre-polynomial residual refit
//     applied after the main systematic fit, see fitEfficiency()'s a_doLegendreResiduleFits
//     doc in headers/EfficiencyFitter.h. Corrected below.)
//   - FIX (2026-07-21, part 2): a_mtm0RebinPion_Plus/etc above (the mT-m0 rebin factors)
//     turned out to have NO effect on AllCent-mode BTOF fits -- confirmed empirically
//     (chi^2/ndf's DENOMINATOR, i.e. the bin count, was identical between a rebin=1 run
//     and a rebin=2 run) and then traced to the actual cause in EfficiencyFitter.cxx:
//     loadPicoBinnerDataFile()'s AllCent path rebins the already-a_mtm0Rebin-rebinned
//     BTOF histogram a SECOND time, unconditionally, via HistogramUtilities::
//     rebinEqualByCDF(..., m_btof_eff_no_fit_num_bins, ...) -- silently overriding
//     whatever the first rebin produced. m_btof_eff_no_fit_num_bins defaults to 80 in
//     EfficiencyFitter's own constructor and is only changed via setNoFitNumBins(),
//     which this port never called. Confirmed every original per-energy macro in the
//     etof repo calls fitter->setNoFitNumBins(25, 120) -- i.e. upstream always uses 120
//     BTOF bins, never the 80 this port was silently stuck at. New a_btofNoFitNumBins
//     argument (default 120, matching upstream) below wires this up. The TPC count (25)
//     is hardcoded rather than exposed, since TPC efficiency isn't fit by this BTOF-only
//     port at all -- 25 is passed only to match upstream convention in case TPC fitting
//     is added here later.
//
// Usage (run from macros/, load LoadEfficiencyFitterLibs.C FIRST in the SAME
// invocation):
//   root -l -q -b LoadEfficiencyFitterLibs.C 'RunEfficiencyFitter.C("yieldHistos_OO200_pion.root","yieldHistos_OO200_kaon.root","efficiency_OO200_BTOF.root")'
// ...or, once a good proton yield file exists, pass it as the 4th argument to fit all
// three species together in one output file:
//   root -l -q -b LoadEfficiencyFitterLibs.C 'RunEfficiencyFitter.C("yieldHistos_OO200_pion.root","yieldHistos_OO200_kaon.root","efficiency_OO200_BTOF.root","yieldHistos_OO200_proton.root")'

#include "../makefile_toggles.h"

void RunEfficiencyFitter(string a_pionPicoBinnerFile,
                          string a_kaonPicoBinnerFile,
                          string a_outFileName,
                          string a_protonPicoBinnerFile    = "",     // optional -- see header note above. Leave "" to skip proton (current default,
                                                                      // since the current proton yield file is corrupt); pass a good file to fit
                                                                      // all 3 species together.
                          string a_imageDir               = "./EfficiencyFitterImages/",
                          int    a_numCentBins             = 6,      // see header note above -- matches SetCutClass.C / RunRawSpectraModifier.C
                          bool   a_convert16to9CentBins    = false,  // see header note above
                          int    a_mtm0RebinPion_Plus      = 1,      // rebin factor for the mT-m0 axis used by the Data-Driven BTOF efficiency
                          int    a_mtm0RebinPion_Minus     = 1,      // fit -- original repo used 1 (no rebin) for every energy it ran, so that's
                          int    a_mtm0RebinKaon_Plus      = 1,      // the default here too. Bump these up if the O+O 200 statistics in a given
                          int    a_mtm0RebinKaon_Minus     = 1,      // mT-m0 bin turn out too sparse to fit cleanly once you look at the output.
                          int    a_mtm0RebinProton_Plus    = 1,      // same, for proton (only used if a_protonPicoBinnerFile != "").
                          int    a_mtm0RebinProton_Minus   = 1,
                          // ####  PHYSICS/STATISTICS CHOICE -- READ BEFORE TRUSTING THE OUTPUT  ####
                          // The BTOF efficiency can be fit two ways:
                          //   a_fitAllCentBTOF=true (default, matches the original tool's own default):
                          //     fit ONE combined-centrality efficiency curve per species/charge/rapidity-bin
                          //     using ALL centralities' data pooled together (most statistics, most stable
                          //     fit, but assumes the BTOF efficiency doesn't depend on centrality).
                          //   a_fitAllCentBTOF=false:
                          //     fit each of the a_numCentBins centralities independently (captures real
                          //     centrality dependence if it exists, but each fit has ~1/6th the statistics
                          //     of the AllCent version -- may be noisy for peripheral bins in a
                          //     small system like O+O).
                          // Either way, RawSpectraModifier needs genuine Cent%02d-named objects to exist
                          // (it has no "AllCent" fallback -- confirmed via its loadAndApplyTPCEffAnd
                          // EnergyLossAndBTOFEffFile() read pattern). So when a_fitAllCentBTOF=true,
                          // a_useAllCentBTOFFits MUST also be true (default here) -- that copies the one
                          // AllCent fit into every Cent%02d slot via copyAllCentBTOFFitsToSpecCent(), so
                          // RawSpectraModifier finds what it expects, just with all centralities sharing
                          // an identical fit. If you want genuinely per-centrality fits instead, set
                          // a_fitAllCentBTOF=false (a_useAllCentBTOFFits is then ignored -- the per-cent
                          // loop runs directly).
                          bool   a_fitAllCentBTOF          = true,
                          bool   a_useAllCentBTOFFits      = true,
                          bool   a_doLegendreBTOFFits      = false, // original per-energy-macro default (a separate, optional Legendre-polynomial
                                                                     // residual refit AFTER the main systematic fit -- NOT the useIntegralGraphFitter
                                                                     // toggle, see the 2026-07-21 FIX note in the header comment above and the
                                                                     // fitter->setUseIntegralChiSqr(false) call below for that actual fix)
                          int    a_btofNoFitNumBins       = 120){    // FIX (2026-07-21, part 2, see header comment above): the ACTUAL knob for
                                                                      // AllCent-mode BTOF fit bin count -- a_mtm0RebinPion_Plus/etc above do NOT
                                                                      // control this, despite their name/doc comment. Default 120 matches every
                                                                      // original per-energy macro upstream.

  bool doProton = (a_protonPicoBinnerFile != "");

  cout << "RunEfficiencyFitter (BTOF-only) starting..." << endl;
  cout << "  a_pionPicoBinnerFile:      " << a_pionPicoBinnerFile << endl;
  cout << "  a_kaonPicoBinnerFile:      " << a_kaonPicoBinnerFile << endl;
  cout << "  a_protonPicoBinnerFile:    " << (doProton ? a_protonPicoBinnerFile : "(none -- proton skipped)") << endl;
  cout << "  a_outFileName:             " << a_outFileName << endl;
  cout << "  a_numCentBins:             " << a_numCentBins << endl;
  cout << "  a_convert16to9CentBins:    " << (a_convert16to9CentBins ? "true" : "false") << endl;
  cout << "  a_fitAllCentBTOF:          " << (a_fitAllCentBTOF ? "true" : "false") << endl;
  cout << "  a_useAllCentBTOFFits:      " << (a_useAllCentBTOFFits ? "true" : "false") << endl;
  cout << "  a_doLegendreBTOFFits:      " << (a_doLegendreBTOFFits ? "true" : "false") << endl;
  cout << "  a_btofNoFitNumBins:        " << a_btofNoFitNumBins << endl;

  if(a_fitAllCentBTOF && !a_useAllCentBTOFFits){
    cout << "WARNING: a_fitAllCentBTOF=true with a_useAllCentBTOFFits=false will NOT produce any"
         << " Cent%02d-named BTOF efficiency objects -- RawSpectraModifier will find nothing to read"
         << " and will fall back to skipping the BTOF correction (see its fracApplied_BTOFEff logic)."
         << " You almost certainly want a_useAllCentBTOFFits=true here, or a_fitAllCentBTOF=false"
         << " instead. Proceeding anyway since you asked for it explicitly." << endl;
  }

  gSystem->mkdir(a_imageDir.c_str(), kTRUE);

  EfficiencyFitter* fitter = new EfficiencyFitter();
  fitter->setImageDir(a_imageDir);
  // FIX (2026-07-21, see header comment above for the full story): explicitly disable
  // the buggy integral-chi-square GraphFitter path -- EfficiencyFitter's own
  // constructor defaults m_useIntegralChiSqr to true, and every original per-energy
  // macro in the source etof repo overrides it to false. Without this call, a
  // non-converging fit (real for this dataset -- hit on KaonPlus, AllCent, y=-0.9)
  // segfaults inside fitEfficiency() on a NULL fitFunctSys->Clone().
  fitter->setUseIntegralChiSqr(false);
  // FIX (2026-07-21, part 2, see header comment above for the full story): wire up the
  // BTOF "no-fit" bin count -- this is the actual knob controlling AllCent-mode BTOF
  // fit bin count, NOT a_mtm0RebinPion_Plus/etc above (confirmed those have zero effect
  // on this path). TPC count (25) hardcoded to match upstream convention -- irrelevant
  // to this port since TPC efficiency isn't fit here at all.
  fitter->setNoFitNumBins(25, a_btofNoFitNumBins);

  // partIndex convention matches the rest of this repo: 0=Pion, 1=Kaon, 2=Proton
  // (proton loaded/fit only if a_protonPicoBinnerFile != "" -- see header note above).
  const int PION   = 0;
  const int KAON   = 1;
  const int PROTON = 2;

  //##################      LOAD PICOBINNER DATA (BOTH CHARGES, PION+KAON[+PROTON])   ##################
  cout << "\nLoading PicoBinner data files for Data-Driven BTOF efficiency..." << endl;
  fitter->loadPicoBinnerDataFile(a_pionPicoBinnerFile, PION,  1, a_numCentBins, a_mtm0RebinPion_Plus,  a_convert16to9CentBins);
  fitter->loadPicoBinnerDataFile(a_pionPicoBinnerFile, PION, -1, a_numCentBins, a_mtm0RebinPion_Minus, a_convert16to9CentBins);
  fitter->loadPicoBinnerDataFile(a_kaonPicoBinnerFile, KAON,  1, a_numCentBins, a_mtm0RebinKaon_Plus,  a_convert16to9CentBins);
  fitter->loadPicoBinnerDataFile(a_kaonPicoBinnerFile, KAON, -1, a_numCentBins, a_mtm0RebinKaon_Minus, a_convert16to9CentBins);
  if(doProton){
    // NOTE: for partIndex==2 (proton) with charge>0, EfficiencyFitter::loadPicoBinnerDataFile()
    // also attempts to load ProtonDCA/rap_mTm0_dca_Proton{Plus,Minus}_Cent_%02d histograms
    // (used for knockout-background fitting, which this macro does NOT fit -- see header note
    // above). That load is guarded internally (skips cleanly if those histograms aren't
    // present), so it's harmless here even though we never use the result.
    fitter->loadPicoBinnerDataFile(a_protonPicoBinnerFile, PROTON,  1, a_numCentBins, a_mtm0RebinProton_Plus,  a_convert16to9CentBins);
    fitter->loadPicoBinnerDataFile(a_protonPicoBinnerFile, PROTON, -1, a_numCentBins, a_mtm0RebinProton_Minus, a_convert16to9CentBins);
  }

  //##################      FIT BTOF EFFICIENCY   #######################
  // fitEfficiency(partIndex, charge, centIndex, doTPCEff, doDataDrivenBTOF, doLegendreResiduleFits)
  // -- doTPCEff is always false here (BTOF-only port); doDataDrivenBTOF is always true.
  if(a_fitAllCentBTOF){
    cout << "\nFitting AllCent (combined-centrality) BTOF efficiency..." << endl;
    fitter->fitEfficiency(PION, 1, -1, false, true, a_doLegendreBTOFFits);
    fitter->fitEfficiency(PION,-1, -1, false, true, a_doLegendreBTOFFits);
    fitter->fitEfficiency(KAON, 1, -1, false, true, a_doLegendreBTOFFits);
    fitter->fitEfficiency(KAON,-1, -1, false, true, a_doLegendreBTOFFits);
    if(doProton){
      fitter->fitEfficiency(PROTON, 1, -1, false, true, a_doLegendreBTOFFits);
      fitter->fitEfficiency(PROTON,-1, -1, false, true, a_doLegendreBTOFFits);
    }

    if(a_useAllCentBTOFFits){
      cout << "Copying AllCent BTOF fits into each Cent%02d slot (a_useAllCentBTOFFits=true)..." << endl;
      fitter->copyAllCentBTOFFitsToSpecCent(PION, 1);
      fitter->copyAllCentBTOFFitsToSpecCent(PION,-1);
      fitter->copyAllCentBTOFFitsToSpecCent(KAON, 1);
      fitter->copyAllCentBTOFFitsToSpecCent(KAON,-1);
      if(doProton){
        fitter->copyAllCentBTOFFitsToSpecCent(PROTON, 1);
        fitter->copyAllCentBTOFFitsToSpecCent(PROTON,-1);
      }
    }
  }else{
    cout << "\nFitting per-centrality BTOF efficiency (a_fitAllCentBTOF=false)..." << endl;
    for(int centIndex = 0; centIndex < a_numCentBins; centIndex++){
      cout << "  Cent" << Form("%02d",centIndex) << "..." << endl;
      fitter->fitEfficiency(PION, 1, centIndex, false, true, a_doLegendreBTOFFits);
      fitter->fitEfficiency(PION,-1, centIndex, false, true, a_doLegendreBTOFFits);
      fitter->fitEfficiency(KAON, 1, centIndex, false, true, a_doLegendreBTOFFits);
      fitter->fitEfficiency(KAON,-1, centIndex, false, true, a_doLegendreBTOFFits);
      if(doProton){
        fitter->fitEfficiency(PROTON, 1, centIndex, false, true, a_doLegendreBTOFFits);
        fitter->fitEfficiency(PROTON,-1, centIndex, false, true, a_doLegendreBTOFFits);
      }
    }
  }

  cout << "\nWriting output file: " << a_outFileName << endl;
  fitter->write(a_outFileName);

  cout << "RunEfficiencyFitter (BTOF-only) done." << endl;
}
