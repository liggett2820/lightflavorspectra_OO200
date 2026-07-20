#ifndef _ZFITTER_
#define _ZFITTER_

// =====================================================================================
// ZFitter.h -- ported from lightflavorspectra_etof for LOCAL, 200 GeV O+O, pi/K/p ONLY use.
//
// WHAT CHANGED FROM THE ORIGINAL (read this before trusting any output):
//
//   The original repo's RunZFitter_Collider.C is the run macro whose eventConfig
//   ("Collider"/ColliderCenter) matches this repo's O+O 200 GeV configuration (as
//   opposed to RunZFitter_Collider_Liggett.C, which is a personal variant for a
//   completely different analysis: 3.2/3.3/6.2/9.7 GeV FixedTarget deuteron/triton/
//   helion/alpha light-nuclei fits). Cross-referencing exactly which ZFitter methods
//   RunZFitter_Collider.C actually calls against where those methods are DEFINED
//   (grep across all fourteen source/ZFitter_*.cxx files ZFitter.cxx #includes) showed
//   that a large fraction of the original class is legacy/superseded/out-of-scope code
//   that compiles but is never exercised by this analysis path. Declarations for that
//   dead code have been removed from this header (and the corresponding .cxx files are
//   simply not being ported), specifically:
//
//   1. FixedTarget-only TPC fitting: fitTPCFXT_SimulCent, fitTPC_FXT_SimulCent_Part2,
//      setFXTTPCMeans (source/ZFitter_TPC_FXT_SimulCent.cxx, entirely FXT-specific --
//      dead for ColliderCenter, same reasoning as PicoBinner's FixedTarget branch
//      removal).
//   2. Legacy/superseded duplicate implementations of the SAME physics that the actual
//      run macro doesn't call, because a newer version replaced it (each pair below:
//      dead first, live/kept second):
//        fitTPCPionKaon, fitTPCPionKaon_SimulCent, fitTPCPionKaon_SimulCent_Part2
//          -- superseded by fitTPCPionKaon_SimulCent_RapidityLooper (kept)
//        fitTPCProton, fitTPCProton_SimulCent
//          -- superseded by fitTPCProton_SimulCent_RapidityLooper (kept)
//        fitBTOFColliderCenter (non-SimulCent)
//          -- superseded by fitBTOF_SimulCent (kept)
//      (RunZFitter_Collider.C's calls to the old versions are literally commented out
//      in favor of the RapidityLooper/_SimulCent calls -- this isn't a guess.)
//   3. ETOF-specific FITTING methods -- fitInvBetaVsMom_ETOF, fitETOF, fitETOFCent
//      (source/ZFitter_eTOFFuncts.cxx) -- are called nowhere in RunZFitter_Collider.C
//      (commented out there too), even though O+O 200 DOES have ETOF and its data
//      still gets loaded/stored (that data-loading code lives in the kept
//      ZFitter.cxx/loadDataHistograms and is NOT removed). Only the actual ETOF fit
//      DRIVERS are unused for this exact analysis path, so ZFitter_eTOFFuncts.cxx is
//      not being ported at this stage.
//
//   IMPORTANT CORRECTION vs. an earlier draft of this file: deuteron/triton/helion/
//   alpha ("NonMixing") fitting methods (fitBTOF_NonMixing, fitBTOF_NonMixing_Liggett,
//   fitBTOFCent_NonMixing, fitBTOF_withBackground, fitBTOFCent_withBackground,
//   fitTPC_NonMixing, fitTPC_NonMixing_Liggett, fitTPCCent_NonMixing,
//   fitTPCCent_NonMixing_Liggett), the Liggett draw variants (drawBichselOnDeDxPlot_
//   Liggett, drawInvBetaPredictionPlot_Liggett), and setFXTTPCMeans were INITIALLY
//   assumed to live in the same dedicated files as the other FXT/Liggett-only code
//   and dropped accordingly -- that assumption was wrong. Grepping the actual ported
//   .cxx files turned up real definitions for all of these INSIDE
//   ZFitter_TPCFuncts.cxx, ZFitter_bTOFFuncts.cxx, and ZFitter.cxx itself -- i.e. the
//   same files that also define the genuinely-live fitTPC/fitTPCCent/fitSoloTPC/
//   fitBTOF/fitBTOFCent methods this repo needs. Those files are organized by
//   detector/topic, not by liveness, so their dead and live methods are physically
//   interleaved throughout. Since those files are being ported anyway (needed for the
//   live methods) and hand-excising individual dead function bodies out of
//   multi-thousand-line files carries real risk of a mismatched brace silently
//   corrupting a neighboring live function, these declarations are being kept here to
//   match the implementations that are, in fact, present and compiled. They are not
//   called anywhere in this repo's run macro (this repo's PicoBinner only fills
//   pi/K/p yield histograms, so these light-nuclei/background/Liggett fitters have no
//   data to operate on even if invoked) -- present but inert, exactly like the
//   internal 9-species contamination-modeling machinery inside the live methods.
//
//   NOTHING about pi/K/p physics/fit logic itself has been changed or simplified --
//   only whole methods/files confirmed dead for this exact run-macro call chain have
//   been dropped. Internal contamination-modeling machinery that pi/K/p fits actually
//   depend on (references to electron/kaon/other species inside the KEPT methods) is
//   left completely untouched, matching the same caution applied throughout
//   PicoBinner's port.
// =====================================================================================

#include <vector>
#include <string>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TString.h"
#include "TPaveStats.h"
#include "TGraphSmooth.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompLU.h"

#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"
#include "../headers/CutClass.h"
#include "../headers/HistogramUtilities.h"
#include "../headers/PhysMath.h"
#include "../makefile_toggles.h"
#include "../headers/GausMixModel.h"
#include "../headers/Histo2D.h"

// _MAC_OSX_ is always defined in this repo (see makefile_toggles.h note) so this stays.
#ifdef _MAC_OSX_
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Fit/ParameterSettings.h"
#endif

using namespace std;

//NOTE: define _ZFITTER_DEBUG_ in makefiletoggles.h to debug issues

const int g_ZFitter_MaxNumRapBins = 35;
const int g_ZFitter_MaxNumMtM0Bins = 300;

// 2026-07: everything the per-bin diagnostic plotting code (in
// ZFitter_PionKaonTPC_SimulCent_SingleLoop.cxx / ZFitter_ProtonTPC_SimulCent_SingleLoop.cxx)
// needs to redraw a single bin's fit-overlay canvas + pull plot, WITHOUT the minimizer,
// the histogram-loading machinery, or any other live-fit state still being around --
// this is what gets written to / read back from the diagnostic snapshot TTree (see
// setSaveDiagnosticSnapshots() below and source/ZFitter_DiagnosticSnapshot.cxx).
//
// "species" here always means 3 slots in a fixed order: slot 0 = pion, slot 1 = kaon,
// slot 2 = electron (driverType 0, PionKaonTPC's skew-pion + 2-Gaussian fit) or proton
// (driverType 1, ProtonTPC's 3-skew-normal fit). Gaussian species just don't use their
// gamma slot (left at 0) -- which render function to call, and which combined-function
// builder it uses internally, is decided entirely by driverType, not by anything stored
// per-species here.
//
// The text-box contents (topTextLines / bottomTextLines) are stored PRE-FORMATTED,
// generated with the exact same Form(...) calls the live plotting code already uses,
// captured at fill-time while the minimizer/paramNames/limits are still valid -- rather
// than trying to re-derive fixed-flags / "Ignore" skips / AT-LIMIT percentages later
// from raw numbers, which would risk silently drifting from what a live run actually
// shows. The render functions just AddText() these lines back in, unchanged.
struct ZFitterDiagnosticSnapshot{
  int driverType; // 0 = PionKaonTPC (pion skew + kaon/electron Gaussian), 1 = ProtonTPC (3 skew)
  string particleName; // m_partInfo->GetParticleName(m_currentPartIndex, charge) at fill time
  int charge;
  int centIndex;
  int rapBin;
  int mtm0Bin;
  int fittingRound;
  string roundTag; // ProtonTPC only (appended to the output filename there); empty for PionKaonTPC

  double lowRange;
  double highRange;
  double currentLowFitRange;
  double currentHighFitRange;
  double sigmaGuess;
  bool   drawInitialSeeds;
  bool   saveNoLogImages;

  double N[3];
  double N_err[3];
  double mu[3];
  double mu_err[3];
  double sigma[3];
  double sigma_err[3];
  double gamma[3]; // unused (0) for Gaussian species
  double gamma_err[3];

  // Seed/initial values used to draw the optional "_init" overlay curves (only drawn
  // when drawInitialSeeds is true). For PionKaonTPC's electron slot this is already the
  // doElectronPionRatio-resolved amplitude (that branch is resolved once at fill time,
  // not re-derived here).
  double seedN[3];
  double seedMu[3];
  double seedSigma[3];
  double seedGamma[3];
  bool   doFit[3]; // doFitPion/doFitKaon/doFitElectron(or Proton) -- ProtonTPC zeroes a seed curve's amplitude when false

  int color[3]; // m_partInfo->GetParticleColor(<global particle index>) at fill time -- slot order doesn't match global particle index for PionKaonTPC's electron slot, so this can't be re-derived from driverType alone

  vector<string> topTextLines;    // fitParsTxt: final fit values + (fixed) tags
  vector<string> bottomTextLines; // fitTxt: initial/seed values + limits (+ AT LIMIT for PionKaonTPC)

  TH1D* hist; //the raw dE/dx histogram that was fit for this bin/centrality -- owned by the caller (either the live driver's m_currentHistosToFit_ByCent, or a clone read back from the snapshot tree); render functions only read it, never delete it
};

class ZFitter{
public:

  //########## DEF IN ZFITTER.CXX ############
  ZFitter();

  //import the output of PicoBinner (has histograms with rap_mTm0_dEdx binning)
  //this is the file the has the dEdx momentom plots and stuff; the one where you passed true as the fourth argument of runPicoBinner
  void loadBasicDataHistograms(string a_filename,  bool a_hasETOF=true);

  //This reads output of PicoBinner, used if the histograms are in different root files (like when processing on RCF)
  //all particles must have the same number of centralities
  void loadDataHistograms(string a_filename, int a_partIndex, int a_numCent, int a_rebinX=1, int a_rebinY=1,  int a_rebinZTPC=1, int a_rebinZBTOF=1, bool a_convert16to9CentBins = false, bool a_combine5and10CentBins = false, bool a_combine40and50_60and70_centBins=false, int a_rebinZETOF = 1, vector< vector<double> > a_dynamicRebin_Plus = {{}},  vector< vector<double> > a_dynamicRebin_Minus = {{}});

  void loadCalibrationHistogram(string a_filename, int a_partIndex);

  // must be done after loadDataHistograms. Creates a file with TH2Ds instead of TH2Is with the rebinning done once. This saves time.
  void createQuickStartFile(string a_filename);

  void loadQuickStartFile(string a_filename, int a_numCent, int a_partIndex = -1, int a_detectorIndex = -1);  //can be used instead of loadDataHistograms if you created a quick file before
  void loadTextFileInformation(string a_textFileList);


  TH1D* grabSpecificHistogram(int a_partIndex, int a_charge, int a_centIndex, int a_detectorIndex, int a_rapBin, int a_mTm0Bin, int a_rebinCDFNumBins = -1, bool a_useProminenceAllCent = false, double a_lowBound = -999, double a_highBound = -999, bool a_divideByNumEvents = false, bool a_propagateNumEvtError = true);
  void  createProjectionFile( string a_outFileName, int a_partIndex, int a_detectorIndex, int a_rebinCDFNumBins = -1);

  void makeFitDataHistos(int a_partIndex = -1, int a_detectorIndex = -1);
  void makeAllCentHistos();
  void setParticleInfo(string a_starVer);
  void setCuts(CutClass* a_cuts){m_cuts = a_cuts;};
  void setParticleToggles(bool *partToggles);
  void setParticleIndex(int a_index){m_currentPartIndex = a_index;};
  void setTesting(bool a_testing = false, int a_div = 3, int a_divMtM0 = 3){m_testing = a_testing; m_divTesting = a_div; m_divTestingMtM0 = a_divMtM0;};
  void setTestingRapBins(vector<int> a_rapBinsToTest){m_rapBinsToTest = a_rapBinsToTest;};
  void setMimicInfo(int* a_partIndicies, double* momentumStart, double* momentumStop);
  void setFullMimicToggle(bool* a_fullParticleMimic);
  void setSeparatePM_bToF_allCent(bool a_separatePM_bToF_allCent = false){m_separatePM_bToF_allCent = a_separatePM_bToF_allCent;};
  void setNdofStudentTMinimum(double a_ndofMin = 1.0){m_NdofStudentTMinimum = a_ndofMin;};

  void usePionMinusFitAsSeedForPionPlus(bool a_toggle = true){m_usePionMinusFitAsSeedForPionPlus = a_toggle;};
  void useMirroredRapidityAsSeed(bool a_toggle = true){m_useMirroredRapidityAsSeed = a_toggle;};
  void useStudentTDistributionsForTPC(bool a_toggle = true){m_useStudentTDistributionsForTPC = a_toggle;};
  void useStudentTDistributionsForTOF(bool a_toggle = true){m_useStudentTDistributionsForTOF = a_toggle;};
  void useMatrixInitializationForTPC(bool a_toggle = true){m_useMatrixInitializationForTPC = a_toggle;};
  void useRootMinimizer(bool a_toggle = true){m_useRootMinimizer = a_toggle;};
  void useGausMixModel(bool a_toggle = true){m_useGausMixModel = a_toggle;};
  void allowCentIntegratedFitsToVary(bool a_toggle = true){m_allowCentIntegratedFitsToVary = a_toggle;};
  void allowSpecCentFitsToVary(bool a_toggle = true){m_allowSpecCentFitsToVary = a_toggle;};
  void allowSpecCentFitsToVary_BTOF(bool a_toggle = true){m_allowSpecCentFitsToVary_BTOF = a_toggle;};
  void allowDeDxParameterExtrapolation(bool a_toggle = true){m_allowDeDxParameterExtrapolation = a_toggle;};
  void useRLBiasedSingleParticleZoom(bool a_toggle = true){m_useRLBiasedSingleParticleZoom = a_toggle;};
  void calculateBTOFIsolatedParameterDifferencesToTotalFits(bool a_toggle = true){m_calculateBTOFIsolatedParameterDifferencesToTotalFits = a_toggle;};
  void usePearsonIVinSingleFits(bool a_toggle = true){m_usePearsonIVinSingle = a_toggle;};
  void useIdenticalRangeOfAllCent(bool a_toggle = true){m_useIdenticalRangeOfAllCent = a_toggle;};
  void massSquaredMode(bool a_toggle = true){m_massSquaredMode = a_toggle;};
  void storeFitsNotStatusZero(bool a_toggle = true){m_storeFitsNotStatusZero = a_toggle;};

  void setSeedFromOppositeChargeDirection(int* a_seedDirections){for(int iii = 0; iii < 9; iii++) m_seedFromOppositeChargeDirection[iii] = a_seedDirections[iii]; };
  void setWeightedAverageMean(bool* a_weightedAverage_mean){for(int iii=0; iii< 9; iii++) m_weightedAverage_mean[iii] = a_weightedAverage_mean[iii]; };
  void setWeightedAverageStdDev(bool* a_weightedAverage_stdDev){for(int iii=0; iii< 9; iii++) m_weightedAverage_stdDev[iii] = a_weightedAverage_stdDev[iii]; };
  void drawInitialSeedsToFits(bool a_toggle = true){m_drawInitialSeedsToFits = a_toggle;};
  //where do you want the pictures of the fits to be generated
  // a_detectorIndex (0=TPC,1=BTOF, <0 = unknown/legacy caller): 2026-07, only creates
  // the detector-relevant image subdirectories instead of the full menu for both
  // detectors -- see the .cxx for the verification behind which folders are actually
  // live. <0 preserves the old fully-inclusive behavior for any caller this wasn't
  // verified against.
  void setImageDir(string a_prePath, string a_imgDirName, int a_detectorIndex = -1);
  void setTextFileDir(string a_textFileDir);

  //LEGACY used to constrain electron parameters to single value
  void setElectronParams(double a_mu, double a_sig, double a_mean, double a_stdDev);

  void setElectronExtrapolationParams(bool a_useElecExtrap, double a_lowMom, double a_highMom){
        m_useElectronExtrapolation = a_useElecExtrap; m_electronExtrapolation_LowMom = a_lowMom; m_electronExtrapolation_HighMom = a_highMom;};
  void setParticleIndividualFitMomHighBounds(vector<double> a_mom_dEdx_plus, vector<double> a_mom_dEdx_minus, vector<double> a_mom_invBeta_plus, vector<double> a_mom_invBeta_minus);
  void setParticleIndividualFitMomHighBoundsETOF(vector<double> a_mom_invBeta_plus, vector<double> a_mom_invBeta_minus);
  void setParticleMaxMomForMeanAndSigmaParToVaryDeDx(vector<double> a_mom_dEdx_plus, vector<double> a_mom_dEdx_minus);
  void setPearsonZTPCMaxMomStop(vector<double> a_mom_plus, vector<double> a_mom_minus);
  void setParticleFitMatrix(vector< vector< bool > > a_plusMatrix, vector< vector< bool > > a_minusMatrix);
  void setParticleFitMatrix(bool a_plusMatrix[9][9], bool a_minusMatrix[9][9]);
  void setPartilceMomStopMatrix(vector< vector< double > > a_plusMatrix_TPC, vector< vector< double > > a_minusMatrix_TPC, vector< vector< double > > a_plusMatrix_BTOF, vector< vector< double > > a_minusMatrix_BTOF);
  void setPartilceMomStopMatrixETOF(vector< vector< double > > a_plusMatrix_ETOF, vector< vector< double > > a_minusMatrix_ETOF);
  void setPartilceExtrapolationMatrix(vector< vector< double > > a_plusMatrix_TPC, vector< vector< double > > a_minusMatrix_TPC, vector< vector< double > > a_plusMatrix_BTOF, vector< vector< double > > a_minusMatrix_BTOF);
  void setPartilceExtrapolationMatrixETOF(vector< vector< double > > a_plusMatrix_ETOF, vector< vector< double > > a_minusMatrix_ETOF);
  void setBTOFExtrapolationRangeMatrix(int a_mainPartIndex, int a_subPartIndex, int a_charge, double a_lowMom, double a_highMom);
  void setTPCExtrapolationRangeMatrix(int a_mainPartIndex, int a_subPartIndex, int a_charge, double a_lowMom, double a_highMom);
  // NonMixing/withBackground/Liggett/FXT-mean methods: implementations are physically
  // colocated inside ZFitter_TPCFuncts.cxx / ZFitter_bTOFFuncts.cxx / ZFitter.cxx
  // (files kept for their live methods) -- see header note above. Not called anywhere
  // in this repo's run macro.
  void fitTPC_NonMixing(int a_particleIndex, double a_momentum_startUsingBkgd = 0.0, bool a_fitMinus = false, bool a_constrainMinusToPlus = true);
  void fitTPC_NonMixing_Liggett(int a_particleIndex, double a_momentum_startUsingBkgd = 0.0, bool a_fitMinus = false, bool a_constrainMinusToPlus = true);

  void fitTPCCent_NonMixing(int a_centIndex, int a_particleIndex, double a_momentum_startUsingBkgd = 0.0, bool a_fitMinus = false);
  void fitTPCCent_NonMixing_Liggett(int a_centIndex, int a_particleIndex, double a_momentum_startUsingBkgd = 0.0, bool a_fitMinus = false);

  void setFXTTPCMeans(double *a_FXTMeanOverwriteMinMom, double *a_FXTMeanOverwriteMaxMom);

  void fitBTOF_withBackground(int a_partIndex, vector<double> a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT = false, bool a_useExpBackground = true, bool a_fitBackgroundOnRight = false, double a_lowMassOverRide = -999, double a_highMassOverRide = -999);
  void fitBTOFCent_withBackground(int a_centIndex, int a_partIndex, vector<double> a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT = false, bool a_useExpBackground = true, bool a_fitBackgroundOnRight = false, double a_lowMassOverRide = -999, double a_highMassOverRide = -999);

  void fitBTOF_NonMixing(int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT = true, double a_lowMassOverRide = -999, double a_highMassOverRide = -999);
  void fitBTOF_NonMixing_Liggett(int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT = true, double a_lowMassOverRide = -999, double a_highMassOverRide = -999);

  void fitBTOFCent_NonMixing(int a_centIndex, int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus = false, bool a_constrainMinusToPlus = true, bool a_useStudentT = true , double a_lowMassOverRide = -999, double a_highMassOverRide = -999);

  void setCanvasParameters(int a_windowWidth,int a_windowHeight, int a_canvasWidth, int a_canvasHeight);
  void setSaveNoLogImages(bool a_toggle = true){m_saveNoLogImages = a_toggle;};
  // Perf toggle added 2026-07: set false to skip generating the per-bin diagnostic
  // fit-overlay PNGs (TPC Pion/Kaon, TPC Proton, and BTOF fitters) -- this is by far the
  // most expensive part of a ZFitter run (canvas render + PNG encode for every
  // rapidity/mTm0/centrality/step bin combination, vastly more than the fit itself), and
  // turning it off does not change the fitted spectra output, only skips the QA plots.
  void setSaveDiagnosticImages(bool a_toggle = true){m_saveDiagnosticImages = a_toggle;};
  // Perf toggle added 2026-07: independent of setSaveDiagnosticImages() above -- either,
  // both, or neither can be on. When on, the per-bin diagnostic plotting code writes a
  // lightweight "snapshot" (the raw fit histogram + everything needed to redraw the same
  // plot -- see ZFitterDiagnosticSnapshot above) into a_snapshotFilePath's TTree instead
  // of (or alongside) rendering a PNG live: no TCanvas, no TF1 construction, no PNG
  // encode -- just a histogram clone + a TTree::Fill(), so this is cheap even when
  // setSaveDiagnosticImages(false) is skipping the expensive part of a run. Run
  // macros/MakeDiagnosticImagesFromSnapshot.C afterward to turn a snapshot file into the
  // same PNGs a live setSaveDiagnosticImages(true) run would have produced, without
  // re-running the fit. See source/ZFitter_DiagnosticSnapshot.cxx.
  void setSaveDiagnosticSnapshots(bool a_toggle, string a_snapshotFilePath){
    m_saveDiagnosticSnapshots = a_toggle; m_diagnosticSnapshotFilePath = a_snapshotFilePath;};
  // Flush and close the snapshot TTree/TFile if setSaveDiagnosticSnapshots(true,...) ever
  // opened one. Call once, after all fitXXX(...) calls are done (see RunZFitter.C) --
  // safe/no-op to call even if snapshots were never turned on.
  void closeDiagnosticSnapshotFile();
  // Shared rendering code for both live plotting and macros/MakeDiagnosticImagesFromSnapshot.C
  // -- defined in source/ZFitter_DiagnosticSnapshot.cxx. Produces byte-identical PNGs to
  // the pre-2026-07 inline plotting code these replaced (see that file's header comment).
  void renderPionKaonTPCDiagnosticPlot(const ZFitterDiagnosticSnapshot& a_snapshot);
  void renderProtonTPCDiagnosticPlot(const ZFitterDiagnosticSnapshot& a_snapshot);
  // Lazily opens m_diagnosticSnapshotFile/m_diagnosticSnapshotTree on first use (only
  // called when m_saveDiagnosticSnapshots is true) and books all branches. Internal
  // glue used from within the two live fit drivers -- not meant to be called directly
  // from a run macro.
  void ensureDiagnosticSnapshotTreeOpen();
  // Clones a_snapshot.hist with a unique name, copies every other field into the bound
  // "current row" scratch members, and Fill()s the tree. Internal glue, same caveat as
  // ensureDiagnosticSnapshotTreeOpen() above.
  void fillDiagnosticSnapshotTreeRow(const ZFitterDiagnosticSnapshot& a_snapshot);
  void setConvertInvBetaToMassSquared(bool a_toggle = true){m_convertInvBetaToMassSquared = a_toggle;}; // for mapping invBeta BTOF into mass^2
  void setUseColliderStopTable(bool a_toggle = true){m_useColliderStopTable = a_toggle;}; // this will stop fitting tpc when it isn't needed anymore

  void setParticleFitMatrixElementPlus(int a_particleOfInterest, int a_secondaryParticleIndex, bool a_toggle){m_partToFitMatrix_Plus[a_particleOfInterest][a_secondaryParticleIndex] = a_toggle;};
  void setParticleFitMatrixElementMinus(int a_particleOfInterest, int a_secondaryParticleIndex, bool a_toggle){m_partToFitMatrix_Minus[a_particleOfInterest][a_secondaryParticleIndex] = a_toggle;};
  void setParticleDeDxMomStopMatrixElementPlus(int a_particleOfInterest, int a_secondaryParticleIndex, double a_momToStopAt){m_partDeDxMomStopMatrix_Plus[a_particleOfInterest][a_secondaryParticleIndex] = a_momToStopAt;};
  void setParticleDeDxMomStopMatrixElementMinus(int a_particleOfInterest, int a_secondaryParticleIndex, double a_momToStopAt){m_partDeDxMomStopMatrix_Minus[a_particleOfInterest][a_secondaryParticleIndex] = a_momToStopAt;};
  void setParticleDeDxMomStartExtrapolationMatrixElementPlus(int a_particleOfInterest, int a_secondaryParticleIndex, double a_momToStopAt){m_partDeDx_MomStartExtrapolation_Matrix_Plus[a_particleOfInterest][a_secondaryParticleIndex] = a_momToStopAt;};
  void setParticleDeDxMomStartExtrapolationMatrixElementMinus(int a_particleOfInterest, int a_secondaryParticleIndex, double a_momToStopAt){m_partDeDx_MomStartExtrapolation_Matrix_Minus[a_particleOfInterest][a_secondaryParticleIndex] = a_momToStopAt;};
  void setParticleBTOFMomStartExtrapolationMatrixElementPlus(int a_particleOfInterest, int a_secondaryParticleIndex, double a_momToStopAt){m_partBTOF_MomStartExtrapolation_Matrix_Plus[a_particleOfInterest][a_secondaryParticleIndex] = a_momToStopAt;};
  void setParticleBTOFMomStartExtrapolationMatrixElementMinus(int a_particleOfInterest, int a_secondaryParticleIndex, double a_momToStopAt){m_partBTOF_MomStartExtrapolation_Matrix_Minus[a_particleOfInterest][a_secondaryParticleIndex] = a_momToStopAt;};

  void setBTOFExtrapolationFunctionVersion(int a_extrap_fn_vers){m_extrap_fn_vers = a_extrap_fn_vers;};

  //THIS CREATES THE OUTPUT FILE (ALSO CONTAINS FIT INFO)
  void makeSpectra(string a_outFileName);

  void setDoKaonDecay(bool a_doKaonDecay){m_doKaonDecay = a_doKaonDecay;};


  //########## DEF IN ZFITTER_TPCFuncts.CXX ############

  //this fits the isolated ZTPC histograms in each particles space
  void fitSoloTPC(int a_ZTPC_Rebin = 1, bool a_saveImages = true);

  //runs the fitting procedure for DeDx
  //fills the empty EtaPtBins from the solo dE/dx fitter (low Pt)
  void fitTPC(bool a_saveImages = true, bool a_useProminenceAllCent = true);

  //fit as function of centrality using results from fitZTPC
  void fitTPCCent(int a_centIndex = -1, bool a_saveImages = true);

  //########## DEF IN ZFITTER_PionKaonTPC_SimulCent_SingleLoop.CXX ############
  // (the only live TPC pi/K fit driver -- see header comment at top of file)
  void fitTPCPionKaon_SimulCent_RapidityLooper(int a_partIndexSpace, int a_charge);

  //########## DEF IN ZFITTER_ProtonTPC_SimulCent_SingleLoop.CXX ############
  // (the only live TPC proton fit driver -- see header comment at top of file)
  void fitTPCProton_SimulCent_RapidityLooper(int a_charge);


  //########## DEF IN ZFITTER_bTOFFuncts.CXX ############

  void fitBTOF(bool a_saveImages = true);
  void fitBTOFCent(int a_centIndex = -1, bool a_saveImages = true);

  //########## DEF IN ZFITTER_BTOF_ColliderCenter_SimulCent.CXX ############
  // (the only live BTOF fit driver -- see header comment at top of file)
  void fitBTOF_SimulCent(int a_partIndexSpace, int a_charge);


  //########## DEF IN ZFITTER_dEdxFuncts.CXX ############

  //this disregards variations in sigma for different rapidities
  //will overwrite electron fits in fitSoloZTPC
  void fitZTPC_Electrons();

  //this fits the dE/dx of single particle with bTof info then zooms in on the full particle dE/dx at low momentum
  void fitDeDx_Solo(int a_rebinMom);

  //LEGACY,BUT WORKS... this fits the dEdx_LogLog Spectra ! RIGHT NOW ONLY FITS PLUS !
  void fitDeDxVsMom();

  //this fits dEdx as a function of betaGamma
  void fitBetaGammaDeDx(bool a_isBTOF = true, double a_betaGammaFitStart = 0.4, bool a_saveImages = false);
  //dEdx fits from Beta Gamma in eta slices from BTOF
  void fitEtaBetaGammaDeDx(string a_storageFileName,bool a_isBTOF = true, double a_betaGammaFitStart = 0.4, bool a_saveImages = false); //set to false to fit etof data (doesn't overwrite btof data)

  void loadEtaBetaGammaDeDxFile(string a_storageFileName, int a_lowIndexToLoad, int a_highIndexToLoad);// if you did the function above, then you don't need to do it again, so just load it. Must be called after loadBasicHistos

  //This creates the data to estimate electron contributions under other bands
  //Step1: this fits a dE/dx projection from a_lowMomRange to a_highMomRange

  void electronStudy(double a_lowMomRange = 0.29, // low momentum value, where electrons first poke out of pion band
                     double a_highMomRange = 0.33, // high mom val before electrons go into kaon band (used for dEdx Projection)
                     double beginDeDxSoloPion = 2.8, // fits the pion first, this is the low dEdx range, recommended slightly after pion peak
                     double endDeDxSoloPion = 3.6, // end dEdx fit range for solo pion, should be right before electron shelf
                     double dEdxStop = 4.75 // this should be slightly before the dip between electrons and kaons
                     );

  void fitElectronDeDxByRapidity(int a_partIndex);


  void findPartShiftsAndRatiosWithRap(vector< double > a_momToFitByParticle);

  void fitIsolatedElectronsByRapidity(int a_partIndex);

  void drawBichselOnDeDxPlot();
  void drawInvBetaPredictionPlot();
  // Liggett-variant draws: implementation colocated in ZFitter.cxx itself (kept for
  // the live methods above) -- see header note near the top of this file. Not called
  // by this repo's (non-Liggett) run macro.
  void drawBichselOnDeDxPlot_Liggett();
  void drawInvBetaPredictionPlot_Liggett();


  void findMassSquaredWidths();

  //draws current status of fits (soloZTPC)
  void drawPartdEdxFits(int a_partIndex);

  void draw_ZVar_mTm0_data_with_fit_data(int a_partIndex, int a_detIndex, int a_centIndex, int a_rap_index, int a_charge);


#ifdef _MAC_OSX_
  //used for the root minimizer
  double chiSquared_TPC(const double* a_params);
  double chiSquared_TOF(const double* a_params);

  double simultaneous_centrality_ZTPC_chisqr(const double* a_param);
  double simultaneous_centrality_ZTOF_chisqr(const double* a_param);
  double simultaneous_centrality_ZTPC_chisqr_E_K_Ratio_Linked(const double* a_param);
  double simultaneous_centrality_ZTPC_chisqr_E_Pi_Ratio_Linked(const double* a_param);
  double simultaneous_centrality_ZTPC_chisqr_tof_ratio_constrain(const double* a_param);
  double simultaneous_centrality_ZTPC_OneSkewNormal_chisqr(const double* a_param);
  double simultaneous_centrality_ZTPC_TwoSkewNormal_chisqr(const double* a_param);
  double simultaneous_centrality_ZTPC_skewNormalAndTwoGaus_chisqr(const double* a_param);
  double simultaneous_centrality_ZTPC_skewNormalAndTwoGaus_sigmaLinked_chisqr(const double* a_param);
  double simultaneous_centrality_ZTPC_skewNormalAndTwoGaus_Proton_chisqr(const double* a_param);
  double simultaneous_centrality_ZTPC_ThreeSkewNormal_Proton_chisqr(const double* a_param);

#endif

private:
  //MEMBER VARIABLES
  CutClass* m_cuts;
  ParticleInfo* m_partInfo;
  int m_numCentralities;
  int m_numParticles;
  int m_currentPartIndex; // which rapidiy space are we in?
  int m_currentCentIndex;
  int m_currentRapBin;
  int m_currentMtM0Bin;
  double m_currentRapidity;
  double m_currentMtM0;
  double m_currentIntegralofHistogram;
  int m_currentPartCharge;
  int m_currentPlusMinusIndex;
  double m_currentLowFitRange;
  double m_currentHighFitRange;
  int m_currentLowFitBin[16];
  int m_currentHighFitBin[16];


  double m_currentMomOfBin;
  int    m_currentNumPartToFit;
  int    m_currentFreeParameters;
  bool   m_inAllCentMode;
  bool   m_separatePM_bToF_allCent;
  bool   m_dEdxMomFittingCompleted;

  double m_canvasParameters[4]; // window width , window height, canvas width, canvas height
  string m_textFileDir;
  int m_extrap_fn_vers;

  #ifdef _MAC_OSX_
    ROOT::Math::Minimizer* m_minimizer;
  #endif

  bool m_useGausMixModel;
  GausMixModel* m_gausMixModel;

  TH1D* m_currentHistoToFit;
  TH1D* m_currentHistoToFitPlus;
  TH1D* m_currentHistoToFitMinus;
  TH1D* m_currentHistosToFit_ByCent[16];
  TH1D* m_currentHistosToFit_ByCent_InvBeta[16]; // when doing BTOF Fitting, some of the code stores both m2 and invBeta histograms
  double m_current_TOF_Ratios_ByCent[16][3]; //cent then K/pi K/p pi/p

  bool   m_initialParticleToggles[9];
  bool   m_partToFitMatrix_Plus[9][9];
  bool   m_partToFitMatrix_Minus[9][9];
  int    m_particleLeftRightSeedDirection[9]; // -1 is left, 1 is right
  double m_partDeDxMomStopMatrix_Plus[9][9];
  double m_partDeDxMomStopMatrix_Minus[9][9];
  double m_partBTOFMomStopMatrix_Plus[9][9];
  double m_partBTOFMomStopMatrix_Minus[9][9];
  double m_partETOFMomStopMatrix_Plus[9][9];
  double m_partETOFMomStopMatrix_Minus[9][9];
  double m_partDeDx_MomStartExtrapolation_Matrix_Plus[9][9];
  double m_partDeDx_MomStartExtrapolation_Matrix_Minus[9][9];
  double m_partBTOF_MomStartExtrapolation_Matrix_Plus[9][9];
  double m_partBTOF_MomStartExtrapolation_Matrix_Minus[9][9];
  double m_partETOF_MomStartExtrapolation_Matrix_Plus[9][9];
  double m_partETOF_MomStartExtrapolation_Matrix_Minus[9][9];

  //ranges to use for parabola fits to fit parameters
  double m_TPC_Extrapolation_Range_Matrix[9][9][2][2]; // partIndex partIndex2 (start/stop mom) (plus/minus)
  double m_BTOF_Extrapolation_Range_Matrix[9][9][2][2]; // partIndex partIndex2 (start/stop mom) (plus/minus)

  bool m_partHasBeenExtrapolatedCharge[9][2]; //normal extrapolation
  bool m_partHasBeenExtrapolated_Solo_Charge[9][2]; // diff between full fit vs solo (BTOF Isolated) fit extrapolation

  bool m_particleLoaded[9];


  //########     TOGGLES ##################
  bool m_usePionMinusFitAsSeedForPionPlus;
  bool m_useMirroredRapidityAsSeed;  // for ColliderCenter, ex. use y=0.5 data as seed to y= -0.5
  bool m_useStudentTDistributionsForTPC;
  bool m_useStudentTDistributionsForTOF;
  bool m_useMatrixInitializationForTPC;
  bool m_useElectronExtrapolation; // these three are for constraining the electron amplitues under the kaons
  double m_electronExtrapolation_LowMom;
  double m_electronExtrapolation_HighMom;
  bool m_electronExtrapolationComplete;
  bool m_electronStudyComplete;
  bool m_bichselShiftsComplete;
  bool m_usePearsonIVinSingle;
  bool m_useIdenticalRangeOfAllCent;
  bool m_saveNoLogImages;
  // Perf toggle added 2026-07: gates the per-bin diagnostic TCanvas/TF1-overlay/PNG
  // generation in the TPC (Pion/Kaon and Proton) and BTOF fitters -- NOT any actual
  // fitting/minimizer/parameter-storage logic, so turning this off does not change the
  // fitted spectra output at all, only skips writing out the per-bin QA plots. Defaults
  // to true so existing behavior (every bin's fit gets a PNG) is unchanged unless a
  // caller opts in via setSaveDiagnosticImages(false). See RunZFitter.C.
  bool m_saveDiagnosticImages;

  // 2026-07: independent snapshot-instead-of/alongside-images toggle -- see
  // setSaveDiagnosticSnapshots() above for the full explanation. m_diagnosticSnapshotFile/
  // Tree are opened lazily (ensureDiagnosticSnapshotTreeOpen()) on first fill, not here.
  bool m_saveDiagnosticSnapshots;
  string m_diagnosticSnapshotFilePath;
  TFile* m_diagnosticSnapshotFile;
  TTree* m_diagnosticSnapshotTree;
  // "current row" scratch space that the tree's branches are bound to -- repopulated and
  // Fill()'d once per bin/centrality when snapshots are on. TH1D* branch needs its own
  // pointer slot separate from the struct (see .cxx) since the struct itself isn't a
  // branch-address target.
  ZFitterDiagnosticSnapshot m_diagnosticSnapshotRow;
  TH1D* m_diagnosticSnapshotHistPtr;
  // bool[3] doesn't get a clean ROOT TTree leaflist type -- this stable-address int[3]
  // shadow of m_diagnosticSnapshotRow.doFit is what the "doFit" branch is actually bound
  // to (see ZFitter_DiagnosticSnapshot.cxx).
  int m_diagnosticSnapshotDoFitInt[3];

  bool m_useRootMinimizer; // not yet complete
  bool m_useGeneralETOFFitDataConstraint;
  bool m_allowCentIntegratedFitsToVary;
  bool m_allowSpecCentFitsToVary;
  bool m_allowSpecCentFitsToVary_BTOF;
  bool m_allowDeDxParameterExtrapolation;
  bool m_useRLBiasedSingleParticleZoom;
  bool m_calculateBTOFIsolatedParameterDifferencesToTotalFits; // correct BTOF isolated dEdx fits
  bool m_drawInitialSeedsToFits;
  bool m_massSquaredMode;
  bool m_doIntegralOfBinForFits;
  bool m_useFitInfoForOppositeCharge;
  bool m_storeFitsNotStatusZero;
  bool m_convertInvBetaToMassSquared;
  bool m_useColliderStopTable;

  int m_seedFromOppositeChargeDirection[9];
  vector<bool> m_initialParticleToggles_Plus;
  vector<bool> m_initialParticleToggles_Minus;
  vector< vector<double> > m_colliderStopTable;

  bool m_weightedAverage_mean[9];
  bool m_weightedAverage_stdDev[9];


  double m_ZTPCshift;
  TF1* m_betaGammaBTOF_dEdx_FitFunction[2]; //for plus minus
  vector< TF1* > m_eta_betaGammaBTOF_dEdx_FitFunction[2]; //for plus minus
  TF1* m_betaGammaETOF_dEdx_FitFunction[2]; //for plus minus
  vector< TF1* > m_eta_betaGammaETOF_dEdx_FitFunction[2]; //for plus minus

  TF1* m_bichselCurves[9];
  TF1* m_invBetaCurves[9];
  double m_isolatedZTPC_starts[9];
  double m_ZTPCStops[9];
  double m_dEdxSeparatedStops[9][2]; // If the current momentum is less than these values, the program will zoom into that peak and fit with one logNormal
  double m_dEdxMaximumMomForTotalFit[9][2]; //this is the last momentum in which a particle can really be seen before disappearing under another particle peak
                                             // after this point, the fits will rely on means and stddev from bTof isolation fits or not fit at all
  double m_bTOFSeparatedStops[9][2]; // If the current momentum is less than these values, the program will zoom into that peak and fit with one gaus
  double m_eTOFSeparatedStops[9][2]; // If the current momentum is less than these values, the program will zoom into that peak and fit with one gaus
  double m_PearsonZTPCMaxMomStop[9][2]; // maximum momentums for the pearsonIV fits (include skewness and kurtosis for first bins)
  double m_lowMassCut[9];
  double m_highMassCut[9];
  vector<double> m_currentParameters;
  vector<int> m_currentLiveParticleIndicies;
  vector<string> m_currentParamNames;
  vector<bool> m_paramToggles; // see if the variable is fixed or not
  double m_NdofStudentTMinimum;
  vector<double> m_previousGoodBTOFParams; // storage for better initialization

  // 0: don't fit (no parameters to fit)
  // 1: use solo dEdx fits (only amplitude free) or you are in centrality fits and only amplitudes need to be fit
  // 2: fit the whole log normal/ gaus (all 3 params free)
  // 3: in a region of separation (zoom in and fit with one distribution, fix all parameters)
  // 4: get bTOF amplitudes of particle and particle referenced in m_particleMimicIndicies
  //    contrain the ratio of amplitudes in fit, the mean varies with the mimic particle mean, the stddev varies with the mimic particle stddev
  //    no free parameters (all parameters "fixed")
  // 5: for electrons, look at mimic particle's deviation from bichsel curve and add it to your fixed mean
  //    the mean varies with the mimic particle mean
  //    only amplitude as free parameter
  // 6: constrains the mean and width but let the amplitude be free. the mean is a shift of the mimic (amplitude free, mean and width fixed)
  int m_currentParticleTags[9];
  vector<string> m_particleTagDescriptions; //descriptions for printout
  int m_particleMimicIndicies[9]; // if tag = 4, then the particle will mimic the one referenced here
  bool m_fullParticleMimic[9];
  double m_particleMimicLowMom[9];
  double m_particleMimicHighMom[9];
  double m_currentConstrainedAmpRaios[9];
  double m_currentMeansForTagFive[9];

  bool m_doKaonDecay;
  vector< TLine* > m_TLinePtrs; // prediciton lines to avoid memeory leaks (easy to delete)

  TH3D* m_rap_mTm0_dEdx_Plus_allCent[9]; //by particle rapidity space
  TH3D* m_rap_mTm0_invBetaBTOF_Plus_allCent[9]; //by particle rapidity space
  TH3D* m_rap_mTm0_invBetaETOF_Plus_allCent[9]; //by particle rapidity space

  TH3D* m_rap_mTm0_dEdx_Minus_allCent[9]; //by particle rapidity space
  TH3D* m_rap_mTm0_invBetaBTOF_Minus_allCent[9]; //by particle rapidity space
  TH3D* m_rap_mTm0_invBetaETOF_Minus_allCent[9]; //by particle rapidity space

  TH3D* m_rap_mTm0_dEdx_allCent[9]; //by particle rapidity space
  TH3D* m_rap_mTm0_invBetaBTOF_allCent[9]; //by particle rapidity space
  TH3D* m_rap_mTm0_invBetaETOF_allCent[9]; //by particle rapidity space

  //dE/dx Mom Fit Data (NO ETA PT BINNING)
  //particle, amp/mean/sigma
  TH1D* m_mom_dEdxFitData[9][3];
  TH1D* m_mom_dEdxFitData_BTOFIsolated_Plus[9][3];
  TH1D* m_mom_dEdxFitData_BTOFIsolated_Minus[9][3];

  TH1D* m_mom_invBetaBTOFFitData[9][3];
  TH1D* m_mom_invBetaETOFFitData[9][3];

  //isolated electron histograms
  TH2I* m_electron_rap_dEdx[9][2];
  TH1D* m_electronDeDxParamsByRapidity[9][3][2];
  TF1*  m_electronDeDxParamFitsByRapidity[9][3][2];

  //For BetaGamma dEdx Fitting //structure: N, mu, sig
  TH1D* m_betaGammaBTOF_DeDx_Data_Plus[3];
  TH1D* m_betaGammaBTOF_DeDx_Data_Minus[3];
  TH1D* m_betaGammaETOF_DeDx_Data_Plus[3];
  TH1D* m_betaGammaETOF_DeDx_Data_Minus[3];

  TGraph* m_betaGammaBTOF_DeDx_Plus_SmoothFitData[3]; //{N,mu,sig}
  TGraph* m_betaGammaBTOF_DeDx_Minus_SmoothFitData[3]; //{N,mu,sig}
  TGraph* m_betaGammaETOF_DeDx_Plus_SmoothFitData[3]; //{N,mu,sig}
  TGraph* m_betaGammaETOF_DeDx_Minus_SmoothFitData[3]; //{N,mu,sig}

  //For Eta BetaGamma dEdx Fitting //structure: N, mu, sig
  TH2D* m_eta_betaGammaBTOF_DeDx_Data_Plus[3];
  TH2D* m_eta_betaGammaBTOF_DeDx_Data_Minus[3];
  TH2D* m_eta_betaGammaETOF_DeDx_Data_Plus[3];
  TH2D* m_eta_betaGammaETOF_DeDx_Data_Minus[3];

  vector< vector< TGraph* > > m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData;//{N,mu,sig} //eta
  vector< vector< TGraph* > > m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData;//{N,mu,sig} //eta
  vector< vector< TGraph* > > m_eta_betaGammaETOF_DeDx_Plus_SmoothFitData; //{N,mu,sig}  //eta
  vector< vector< TGraph* > > m_eta_betaGammaETOF_DeDx_Minus_SmoothFitData; //{N,mu,sig} //eta

  //#############          CENTRALITY INTEGRATED FIT DATA       ##########################
  TH2D* m_Fit_Data_ZTPC[9][9][5][2]; // partIndex for rapidity space, particle fit, amp mean sig ndf, charge 0=plus 1-minus
  TH2D* m_Fit_Data_Solo_ZTPC[9][9][4][2];// partIndex for rapidity space, particle fit, amp mean sig ndf, charge 0=plus 1-minus
  TH2D* m_Fit_Data_Solo_ZTPCShifts[9][9][2];// partIndex for rapidity space, particle fit, charge 0=plus 1-minus
  TH2D* m_Fit_Data_ZbTOF[9][9][5][2]; // particle rap space, particle , parameter, charge 0=plus 1-minus
  TH2D* m_Fit_Data_ZeTOF[9][9][5][2]; // particle rap space, particle , parameter, charge 0=plus 1-minus

  TH2D* m_Fit_Data_Cent_ZTPC[9][9][16][5][2]; // partIndex for rapidity space, particle fit, centrality amp mean sig ndf, charge 0=plus 1-minus
  TH2D* m_Fit_Data_Cent_ZbTOF[9][9][16][5][2]; // particle rap space, particle , centrality, parameter, charge 0=plus 1-minus
  TH2D* m_Fit_Data_Cent_ZeTOF[9][9][16][5][2]; // particle rap space, particle , centrality, parameter, charge 0=plus 1-minus

  TH2D* m_ZTPC_LowRange[9][2]; // particle charge 0=plus 1=minus
  TH2D* m_ZTPC_HighRange[9][2]; // particle charge 0=plus 1=minus
  TH2D* m_ZbTOF_LowRange[9][2]; // particle charge 0=plus 1=minus
  TH2D* m_ZbTOF_HighRange[9][2]; // particle charge 0=plus 1=minus
  TH2D* m_ZbTOF_partBinning_LowRange[9][9][2]; // particle binning, particle charge 0=plus 1=minus
  TH2D* m_ZbTOF_partBinning_HighRange[9][9][2]; // particle binning, particle charge 0=plus 1=minus

  //if using non-mixing fit functions, you will have background parameters
  TH2D* m_Background_FitData_ZTPC[9][4][2]; // particle of interest,  parameter, plus minus
  TH2D* m_Background_FitData_ZbTOF[9][4][2]; // particle of interest,  parameter, plus minus
  TH2D* m_Background_FitData_ZeTOF[9][4][2]; // particle of interest,  parameter, plus minus

  //if using fit functions with background, you will have background parameters
  TH2D* m_Background_FitData_partBinning_ZbTOF[9][9][5][2]; // binning particle, particle of interest,  parameter, plus minus

  TH2D* m_ChiSqr_ZTPC[9][2];
  TH2D* m_ChiSqr_ZbTOF[9][2];
  TH2D* m_ChiSqr_ZeTOF[9][2];
  TH2D* m_FitStatus_ZTPC[9][2];
  TH2D* m_FitStatus_ZbTOF[9][2];
  TH2D* m_FitStatus_ZeTOF[9][2];


  //centrality specific fit data

  TH2D* m_ChiSqr_Cent_ZTPC[9][16][2]; //particle centrality charge
  TH2D* m_ChiSqr_Cent_ZbTOF[9][16][2];
  TH2D* m_ChiSqr_Cent_ZeTOF[9][16][2];
  TH2D* m_FitStatus_Cent_ZTPC[9][16][2];
  TH2D* m_FitStatus_Cent_ZbTOF[9][16][2];
  TH2D* m_FitStatus_Cent_ZeTOF[9][16][2];

  Histo2D* m_mTm0_ZTPC_by_Rapidity[9][16][g_ZFitter_MaxNumRapBins][2];
  Histo2D* m_mTm0_ZBTOF_by_Rapidity[9][16][g_ZFitter_MaxNumRapBins][2];
  Histo2D* m_mTm0_ZETOF_by_Rapidity[9][16][g_ZFitter_MaxNumRapBins][2];



  TH2D* m_Spectra_Cent_ZTPC[9][16][2];
  TH2D* m_Spectra_Cent_ZbTOF[9][16][2];
  TH2D* m_Spectra_Cent_ZeTOF[9][16][2];



  //dEdx amplitude storage for dEdx constraints
  TH2D* m_dEdxYieldsCentPlus[9][9][16]; // particle space ; particle of interest ; centrality
  TH2D* m_dEdxYieldsCentMinus[9][9][16]; // particle space ;  particle of interest ; centrality

  //TOF amplitude storage for dEdx constraints
  TH2D* m_InvBetaBTOFYieldsCentPlus[9][9][16]; // particle space ; particle of interest ; centrality
  TH2D* m_InvBetaBTOFYieldsCentMinus[9][9][16]; // particle space ;  particle of interest ; centrality


  //filled by findPartShiftsAndRatiosWithRap();
  TH1D* m_bichselShifts_Plus[9][9]; //[particle space][particle]
  TH1D* m_bichselShifts_Minus[9][9]; //[particle space][particle]
  TH1D* m_tpcResolution_Plus[9][9]; //[particle space][particle]
  TH1D* m_tpcResolution_Minus[9][9]; //[particle space][particle]
  TH1D* m_quickYieldNums_Plus[9][9]; //[particle space][particle]
  TH1D* m_quickYieldNums_Minus[9][9]; //[particle space][particle]

  double m_mass2Widths[9][3]; //particle mom range
  double m_mass2WidthsMomentum[3];

  bool m_testing;  // will choose a subset of things to fit so it takes shorter time
  bool m_hasETOF;  // some datasets don't have an etof detector
  int  m_divTesting; // number of rapidity bins to do when testing
  int m_divTestingMtM0; // doe every m_divTestingMtM0 mTm0 bin
  vector< int > m_rapBinsToTest;
  string m_imagePreDir;
  string m_imgDirName;
  double m_basicElectronPionRatio;
  double m_basicKaonPionRatio;
  double m_electronParams[6];
  double m_dEdxParams_Electrons[2][4];// {positron, electron} {mu,sig, mean,stdDev}

  TF1* m_logNormFunct;
  TF1* m_gausFunct;
  TF1* m_student_T_Funct;

  TH1D* m_rebinZTPCHistoStorage[1000];
  double m_rebinZTPCHistoIntegralStorage[1000];
  TH1D* m_rebinZTPCHistoStorage_ByCent[16][1000];
  TH1D* m_electronRatioFitHistos_ByCent[16][3];



  //######### LOADED DATA FROM PICOBINNER ##########

  TH1I* m_centEvents[9];
  TH2I* m_mom_dEdx_Plus;
  TH2I* m_mom_dEdx_Minus;
  TH2I* m_betaGammaBTOF_dEdx_Plus;
  TH2I* m_betaGammaBTOF_dEdx_Minus;
  TH2I* m_betaGammaETOF_dEdx_Plus;
  TH2I* m_betaGammaETOF_dEdx_Minus;
  TH3I* m_eta_betaGammaBtof_dEdx_Plus;
  TH3I* m_eta_betaGammaBtof_dEdx_Minus;
  TH3I* m_eta_betaGammaEtof_dEdx_Plus;
  TH3I* m_eta_betaGammaEtof_dEdx_Minus;
  TH2I* m_mom_dEdx_BTOFIsolated_Plus[9]; // by particle
  TH2I* m_mom_dEdx_BTOFIsolated_Minus[9]; // by particle
  TH2I* m_mom_invBetaBTOF_Plus;
  TH2I* m_mom_invBetaBTOF_Minus;
  TH2I* m_mom_invBetaETOF_Plus;
  TH2I* m_mom_invBetaETOF_Minus;

  TH2I* m_electronReduct;

  TH3D* m_rap_mTm0_dEdx_allCent_PionElectron[9]; //by particle rapidity space



  TH3D*  m_rap_mTm0_dEdx_AllCentRead_Plus[9];
  TH3D*  m_rap_mTm0_dEdx_AllCentRead_Minus[9];

  TH3D*  m_rap_mTm0_dEdx_Plus_allCentIsolated[9][9]; //by particle rapidity space, isolated particle
  TH3D*  m_rap_mTm0_dEdx_Minus_allCentIsolated[9][9]; //by particle rapidity space, isolated particle
  TH3D*  m_rap_mTm0_dEdx_Plus[9][16]; // by particle rapidity space, by centrality bin
  TH3D*  m_rap_mTm0_invBetaBTOF_Plus[9][16];// by particle rapidity space, by centrality bin
  TH3D*  m_rap_mTm0_invBetaETOF_Plus[9][16];// by particle rapidity space, by centrality bin
  TH3D*  m_rap_mTm0_dEdx_Minus[9][16];// by particle rapidity space, by centrality bin
  TH3D*  m_rap_mTm0_invBetaBTOF_Minus[9][16];// by particle rapidity space, by centrality bin
  TH3D*  m_rap_mTm0_invBetaETOF_Minus[9][16];// by particle rapidity space, by centrality bin



  //########      PRIVATE FUNCTIONS    #############
  //########  THAT DON'T NEED TO BE CALLED IN THE MACRO #########

  void initZTPCParams();
  void initZTPCParams_student();
  void initBTOFParams();
  void initBTOFParams_student();
  void initETOFParams();
  void initETOFParams_student();
  void initZTPCFunction(TF1* a_funct);
  void initZTPCFunction_student(TF1* a_funct);
  void initBTOFFunction(TF1* a_funct);
  void initBTOFFunction_student(TF1* a_funct);
  void initETOFFunction(TF1* a_funct);
  void initETOFFunction_student(TF1* a_funct);
  void initdEdxMomParams();
  void initDEDXFunction(TF1* a_funct);
  void initInvBetaMomParams();
  void initInvBetaFunction(TF1* a_funct);
  #ifdef _MAC_OSX_
  void initMinimizer_TPC();
  void initMinimizer_TOF();
  #endif

  vector<TF1*> drawIndividualFits_LogNorm();
  vector<TF1*> drawIndividualFits_Gaus();
  vector<TF1*> drawIndividualFits_Student();
  void DrawPredictionLines_TOF(TCanvas *canvas);
  void DrawPredictionLines_dEdx(TCanvas *canvas);
  void printCurrentParams();
  //this will get you the position of the counts in m_currentParmters
  int getCurrParamIndex(int a_partIndex);
  double gaus(double a_mean, double a_sig, double a_x);
  double gausAmp(const double* a_params, double a_x);
  double logNorm(double a_mu, double a_sig, double a_x);
  double logNormAmp(double* a_params, double a_x);
  double studentT(double a_mean, double a_sig, double a_ndf, double a_x);
  double studentTAmp(double* a_params, double a_x);
  double dEdxBetaGammaOperator(Double_t *x,  Double_t *par);

  double predictZTPC(int a_partIndexOfInterest, int a_partIndexZTPCSpace, double a_mom);
  double predictZTOF(int a_partIndexOfInterest, int a_partIndexZTPCSpace, double a_mom);
  double predictZTPC_fromEtaBetaGammaFits(int a_partIndexOfInterest, int a_partIndexZTPCSpace, int a_charge, double a_rapidity, double a_mTm0);
  void getRangeDeDx(double* a_returnData,int a_partIndex); // returns low then high as a_returnData[0] and [1]
  void getRangeDeDx(double &a_lowDeDx, double &a_highDeDx); // given tags and multiple particles
  void getRangeZTPC(double &a_lowZTPC, double &a_highZTPC);
  void getRangeZTPC_Simple(double &a_lowZTPC, double &a_highZTPC, double a_nSigma = 4.0); // no particle tags or matrix toggle stuff
  void getRangeZTOF(double &a_lowZTOF, double &a_highZTOF, double a_lowMassSquaredOverRide = -999, double a_highMassSquaredOverRide = -999);
  void getRangeZTOF_Simple(double &a_lowZTOF, double &a_highZTOF, double a_lowMassOverRide = -999, double a_highMassOverRide = -999); // no particle tags/ matrix toggle stuff
  bool hasTOFRatioToConstrainDeDx(int a_partIndex); // sees if there are actual numbers in the tof fits at the current bin// used for dedx tag 4
  TLine* PredictionLine(TLine* line, Double_t prediction, Double_t height);

  void getRangeZTOF_mass2(double &a_lowZTOF, double &a_highZTOF, int a_partIndex);

  //for the minimizer
  double ZTPCFunct(double a_dEdx);
  double ZTPCFunct_student(double a_dEdx);
  double ZTOFFunct(double a_invBeta);
  double ZTOFFunct_student(double a_ZTOF);
  double dEdxBetaGammaFunct(double a_betaGamma);
  double ZTPC_prediction_betagamma_fit(int a_partIndex,double a_y, double a_mTm0);
  double chiSqr(const double* a_params); //not used right now
  double ZTPCOperator(double* a_x,const double* a_params);
  double ZTOFOperator(double* a_x,const double* a_params);
  double ZTOFOperator_student(double* a_x, const double* a_params);
  double dEdxFunct(double a_dEdx);
  double dEdxOperator(double* a_x,double* a_params);
  double invBetaFunct(double a_invBeta);
  double invBetaOperator(double* a_x,double* a_params);

  bool areDeDxParamsRidiculous();
  bool areBTOFParamsRidiculous();
  bool doThisTPCBinByColliderStopTable(int a_partIndex,double a_rapidity,double a_momentum);


  //when paramTag = 3, zoom into dEdx when the particle is separated from other contamination and fit alone
  //returns N,mu,sig
  vector<double> fitIsolation_dEdx(int a_particleIndex);
  vector<double> fitIsolation_ZTPC_Pearson(); // just for particle of interest
  vector<double> fitIsolation_ZTPC_withBackground(vector<double> a_fixedParameters);
  vector<double> fitIsolation_BTOF(int a_particleIndex);
  vector<double> fitIsolation_BTOF_withBackground(vector<double> a_fixedParameters, bool a_useStudentT);
  vector<double> fitIsolation_ETOF(int a_particleIndex);
  vector<double> fit_BTOF_withBackground(int a_particleIndex,vector<double> a_fixedParameters, vector<double> a_momentum_startUsingBkgd, bool a_useStudentT, bool a_useExpBackground = true, bool a_fitBackgroundOnRight = false);


  void makeElectronExtrapolationFunction();


  // Some functions for the single loop pion kaon TPC SimulCent
  bool calculateAverageElectronWidth();
  bool calculateElectronRatios(); // returns status (true = good)
  void calculateMeanShifts(double &a_delta_Pi,double &a_delta_K,double &a_delta_E);
  bool calculatePionGammaLowPtExtrapolation();
  bool calculateElectronMeanExtrapolation();
  void calculateProtonMeanShifts(double &a_delta_Pi,double &a_delta_K,double &a_delta_P);
  bool calculateProtonSkewExtrapolation(double &a_skew);
  void fitTPCPionKaon_SimulCent_ByRapidity();
  void fitTPCProton_SimulCent_ByRapidity();

  //Used Inside FitTPC to fix the means and widths
  void extrapolateMeanAndWidthParams_TPC(int a_partIndex);
  void extrapolateMeanAndWidthParams_TPC_SecondOrder();
  void extrapolateMeanAndWidthParams_BTOF(int a_partIndex);
  void extrapolateMeanAndWidthParams_BTOF_m2();
  void extrapolateMeanAndWidthParams_FromSoloTPC_TPC(int a_partIndex);
  bool extrapolate_BTOF_Fit_Parameters(int a_particleIndexToExtrapolate, double a_lowMomFitRange, double a_highMomFitRange, int a_functionVersion, int a_centIndex = -1, bool a_meanIsInvBeta = false, int a_specificParameterIndex = -999, int a_convertInvBetaToM2 = false);
  //using the dEdxVsMom plot, it will genereate a ZTPC histo for the current bin
  TH1D* getZTPC_Example();

  // For specifically the Histo2D class
  void make_TPC_Histo2Ds(int a_detIndex = -1); // -1 for all detectors, 0 = TPC, 1 = BTOF, 2 = ETOF


};

#endif
