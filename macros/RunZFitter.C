// RunZFitter.C -- lightflavorspectra_OO200
//
// Ported from the original repo's macros/RunZFitter_Collider.C, which is the run
// macro whose eventConfig ("Collider"/ColliderCenter) matches this repo's O+O 200 GeV
// configuration (RunZFitter_Collider_Liggett.C is a different, FXT light-nuclei
// analysis and was not used as a basis here -- see headers/ZFitter.h for the fuller
// version of this reasoning).
//
// WHAT CHANGED FROM THE ORIGINAL, AND WHY IT'S SAFE:
//
// RunZFitter_Collider.C is actually a 7-ENERGY SCAN macro (27/19.6/14.6/17.3/11.5/9.2/
// 7.7 GeV) selected by an a_energyIndex argument, with a per-energy array of hardcoded
// absolute file paths (Matthew's /home/matthewharasty/ResearchData/... directories)
// and two per-energy booleans (hasETOF, convertInvBetaToMassSquared). Checking every
// other configuration value used in that macro -- every fit-matrix, momentum-range,
// rebinning array, mimic-particle setting, and every fitter->setXXX(...) toggle --
// confirmed NONE of them actually vary by a_energyIndex; they are single hardcoded
// values applied identically across all 7 energies. That means this port only had two
// real jobs: (1) replace the 7-way path lookup table with plain arguments, since this
// repo only ever runs one configuration (O+O, 200 GeV, local laptop, not SDCC), and
// (2) pick the correct fixed values for hasETOF and convertInvBetaToMassSquared for
// this configuration -- hasETOF = true (PicoBinner.cxx defines _HAS_ETOF_ for this
// repo, unlike the original's 27 GeV entry, which was the one energy without ETOF) and
// convertInvBetaToMassSquared = true (identical across all 7 original entries, so
// unambiguous). Every fit-matrix/range/toggle value below is copied verbatim from the
// original macro -- none of it was re-tuned or guessed for 200 GeV, since nothing in
// the original macro indicated these values were energy-dependent in the first place.
//
// Also dropped, and why:
//   - createProblemProjections / createAllProjectionFiles branches: debug-only
//     utilities (both toggled off by default in the original too) that call
//     grabSpecificHistogram/createProjectionFile with a hardcoded list of Matthew's
//     specific problem-bin coordinates from the BES-II analysis -- not meaningful for
//     this dataset. The underlying ZFitter methods are still ported and available if
//     you want to use them by hand later.
//   - create2DProjectionFile branch: calls a ZFitter::create2DProjectionFile method
//     whose defining file (source/ZFitter_2DFits.cxx) was never in scope for this port
//     (see headers/ZFitter.h) -- it isn't declared on this repo's ZFitter class at all.
//   - The commented-out `doingDeuteronTriton` branch (deuteron/triton NonMixing fits):
//     was ALREADY fully commented out in the original macro (hardcoded to the `else`
//     path via a comment-delimiter trick) since doingDeuteronTriton is always false for
//     this analysis. Only the live `else` content (BTOF then TPC fitting loops, pi/K/p
//     only) is kept here; the dead deuteron/triton calls (which also reference the
//     NonMixing methods this repo doesn't drive from any macro) were dropped rather
//     than carried forward as inert comments.
//   - Commented-out calls to methods this port never uses either way (electronStudy,
//     findPartShiftsAndRatiosWithRap, fitElectronDeDxByRapidity,
//     setElectronExtrapolationParams, fitDeDx_Solo, fitDeDxVsMom,
//     fitInvBetaVsMom_ETOF) were left out of this file entirely rather than carried
//     forward as dead comments, since they were already inert in the original too.
//
// USAGE (run from within macros/, loading LoadZFitterLibs.C first in the same
// invocation, matching the LoadPicoBinnerLibs.C / RunPicoBinner.C pattern):
//
//   root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(0,0,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'
//
// a_partIndex:    0=Pion, 1=Kaon, 2=Proton
// a_detectorIndex: 0=TPC, 1=BTOF  (ETOF fit drivers were not ported -- see headers/ZFitter.h)

#include "../makefile_toggles.h"
// Needed so cling knows these types when this macro is interpreted in a fresh ROOT
// session -- gSystem->Load() (done in LoadZFitterLibs.C, run first) loads the compiled
// implementations, but doesn't by itself guarantee the class declarations are known to
// the interpreter. SetCutClass() itself is declared inside the SetCutClass_C dictionary
// loaded by LoadZFitterLibs.C -- do NOT #include "SetCutClass.C" here too, or you get
// the same ambiguous-duplicate-declaration problem documented in RunPicoBinner.C.
#include "../headers/CutClass.h"
#include "../headers/ZFitter.h"

void RunZFitter(int a_partIndex, int a_detectorIndex,
                string a_pionYieldFile, string a_kaonYieldFile, string a_protonYieldFile,
                string a_outputDir, string a_etaBetaGammaDeDxFile,
                string a_quickStartFile = "",
                bool   a_useQuickLoad = true,
                bool   a_createQuickStart = false,
                bool   a_loadEtaBetaGamma = true,
                bool   a_loadTextData = false,
                string a_textFileList = "",
                int a_rapBin_1=-1,int a_rapBin_2=-1,int a_rapBin_3=-1,
                int a_rapBin_4=-1,int a_rapBin_5=-1,int a_rapBin_6=-1){

  // ################  FIXED CONFIGURATION FOR THIS REPO (O+O, 200 GeV, ColliderCenter)  ################
  bool convertInvBetaToMassSquared = true; // identical across all 7 entries in the original macro
  bool hasETOF                     = true; // this repo's PicoBinner defines _HAS_ETOF_ (unlike the
                                            // original's one ETOF-less entry, 27 GeV)

  string partNames[3] = {"Pion","Kaon","Proton"};
  string detNames[4]  = {"TPC","BTOF","ETOF","NA"};
  bool combineImageDirectories = true;
  string imageDirName;
  if(combineImageDirectories){
    imageDirName = Form("yieldFits_OO200_%s_%s_AllRap",partNames[a_partIndex].c_str(),detNames[a_detectorIndex].c_str());
  }else{
    imageDirName = Form("yieldFits_OO200_%s_%s_RapBin_%d_%d_%d_%d_%d_%d",partNames[a_partIndex].c_str(),detNames[a_detectorIndex].c_str(),
            a_rapBin_1,a_rapBin_2,a_rapBin_3,a_rapBin_4,a_rapBin_5,a_rapBin_6);
  }
  string spectraFileName = Form("spectra_OO200_%s_%s_RapBin_%d_%d_%d_%d_%d_%d.root",
             partNames[a_partIndex].c_str(),detNames[a_detectorIndex].c_str(),a_rapBin_1,a_rapBin_2,a_rapBin_3,a_rapBin_4,a_rapBin_5,a_rapBin_6);
  string textDirName = a_outputDir + "yieldFits_OO200_text";

  cout << "FILES BEING USED:" << endl;
  cout << "  pion yield file:   " << a_pionYieldFile << endl;
  cout << "  kaon yield file:   " << a_kaonYieldFile << endl;
  cout << "  proton yield file: " << a_protonYieldFile << endl;
  cout << "  eta-betaGamma dEdx file: " << a_etaBetaGammaDeDxFile << endl;
  cout << "  output dir:        " << a_outputDir << endl;
  cout << "  image dir name:    " << imageDirName << endl;
  cout << "  spectra file name: " << spectraFileName << endl;

  int  nCentToRun = 5; // matches CutClass's 5-bin scheme (0-5,5-10,10-20,20-40,40-80) set in SetCutClass.C
  bool onlyDoOneCent = false;

  vector<int> rapBinsToTest = {};
  if(a_rapBin_1 != -1) rapBinsToTest.push_back(a_rapBin_1);
  if(a_rapBin_2 != -1) rapBinsToTest.push_back(a_rapBin_2);
  if(a_rapBin_3 != -1) rapBinsToTest.push_back(a_rapBin_3);
  if(a_rapBin_4 != -1) rapBinsToTest.push_back(a_rapBin_4);
  if(a_rapBin_5 != -1) rapBinsToTest.push_back(a_rapBin_5);
  if(a_rapBin_6 != -1) rapBinsToTest.push_back(a_rapBin_6);

  bool testing             = false;
  int  divNumForTesting    = 1;
  int  mTm0DivNumForTesting = 1;

  bool fitPionTPC    = false;
  bool fitPionBTOF   = false;
  bool fitKaonTPC    = false;
  bool fitKaonBTOF   = false;
  bool fitProtonTPC  = false;
  bool fitProtonBTOF = false;
  if(a_partIndex == 0){
    if(a_detectorIndex == 0) fitPionTPC = true;
    if(a_detectorIndex == 1) fitPionBTOF = true;
  }else if(a_partIndex == 1){
    if(a_detectorIndex == 0) fitKaonTPC = true;
    if(a_detectorIndex == 1) fitKaonBTOF = true;
  }else if(a_partIndex == 2){
    if(a_detectorIndex == 0) fitProtonTPC = true;
    if(a_detectorIndex == 1) fitProtonBTOF = true;
  }

  int  selectiveCharged = 0; // 0 for both
  bool massSquaredMode  = true; // different than convertInvBetaToMassSquared above

  bool doKaonDecay       = false;
  bool combine5and10Cent = false; // legacy
  bool convert16to9      = false; // must edit CutClass after if you flip this

  int rapidityRebin = 1; // desyncs from CutClass if changed from 1 -- be careful
  int mTm0Rebin_byParticle[9] = {1,1,1,0,1,1,5,5,0};

  // particleMimicIndicies/fullParticleMimic/mom start-stop are part of the internal
  // 9-species contamination-modeling machinery this repo preserves in full -- see
  // headers/ZFitter.h. Values copied verbatim.
  int particleMimicIndicies[9] = {0,0,0,0,0,0,0,0,0}; // everything mimics a pion by default
  bool fullParticleMimic[9] = {true,true,true,true,true,true,true,true,true};
  double particleMimicMomStart[9] = {0 ,0.0 ,0.0,0 ,0 ,0 ,0 ,0 ,0};
  double particleMimicMomStop[9]  = {0 ,0.0 ,0.0 ,0 ,0 ,0 ,0 ,0 ,0}; // zeros mean don't mimic

  bool useRootMinimizer              = false;
  bool allowCentIntegratedFitsToVary = true;
  bool allowSpecCent_DeDxToVary      = false;
  bool allowSpecCent_ToVary_BTOF     = false;

  // partToggles: which of the 9 species are tracked at all internally. electron/
  // deuteron/triton stay true (contamination-modeling support for the pi/K/p bands);
  // helion/alpha/muon stay false, matching the original macro exactly.
  bool partToggles[9] = {true,true,true,   //pion,kaon,proton
                         true,true,true,   //electron, deuteron, triton
                         false,false,false}; //helion, alpha, muon

  bool usePionMinusAsSeedForPionPlus = false;
  bool useStudentTForTOF             = true;
  bool useStudentTForTPC             = false; // only implemented for solo de/dx right now
  bool separatePM_bToF_AllCent       = true;
  bool drawInitialSeedsToFits        = true;

  // stop momentum to zoom in and fit with single bump (per particle-space, plus/minus)
  vector<double> a_mom_dEdx_plus      = {0.3, 0.4, 0.8, 0.0, 1.4, 1.4, 0.5, 0.5, 0.1};
  vector<double> a_mom_dEdx_minus     = {0.3, 0.4, 0.8, 0.0, 1.4, 1.4, 0.5, 0.5, 0.1};
  vector<double> a_mom_invBeta_plus   = {0.65, 0.65, 1.2, 0.0, 1.3, 1.3, 0.0, 0.0, 0.0};
  vector<double> a_mom_invBeta_minus  = {0.65, 0.65, 1.2, 0.0, 1.3, 1.3, 0.0, 0.0, 0.0};
  vector<double> mom_Pearson_plus     = {0.2, 0.4, 0.4, 0.0, 0.9, 0.9, 0.5, 0.5, 0.1};
  vector<double> mom_Pearson_minus    = {0.2, 0.4, 0.4, 0.0, 0.9, 0.9, 0.5, 0.5, 0.1};
  vector<double> a_mom_dEdx_maxVary_plus  = {0.85, 0.85, 1.4, 0.0, 1.9, 2.0, 2.0, 2.0, 0.0};
  vector<double> a_mom_dEdx_maxVary_minus = {0.85, 0.85, 1.4, 0.0, 1.9, 2.0, 2.0, 2.0, 0.0};

  vector< vector< double > > dynamicRebinInfo_Plus = {
    {0,1,2,1.3,1.765},
    {0,1,4,1.73,4.5},
    {1,1,2,0.9,1.47},
    {1,1,4,1.47,4.5},
    {2,1,2,0.6,0.95},
    {2,1,3,0.92,1.38},
    {2,1,4,1.35,4.5}
  };
  vector< vector< double > > dynamicRebinInfo_Minus = {
    {0,1,2,1.3,1.765},
    {0,1,4,1.73,4.5},
    {1,1,2,0.9,1.47},
    {1,1,4,1.47,4.5},
    {2,0,4,0,5},
    {2,0,2,0.5,5},
    {2,0,2,1,5},
    {2,1,4,0.0,5},
    {2,1,2,0.6,5},
    {2,1,2,0.89,5},
    {2,1,2,1.32,5}
  };

  vector<vector< bool > > plusFitMatrix = {
    {true ,true ,true,true ,false,false,false,false,false},
    {true ,true ,true,true ,false,false,false,false,false},
    {true ,true ,true ,false,false,false,false,false,false},
    {false,false,false,false,false,false,false,false,false},
    {false,false,true ,false,true ,false,false,false,false},
    {false,false,false,false,true ,true,false,false,false},
    {false,false,false,false,false,false,true,false,false},
    {false,false,false,false,false,false,false,true,false},
    {false,false,false,false,false,false,false,false,false}
  };
  vector<vector< bool > > minusFitMatrix = {
    {true ,true ,true,true,false,false,false,false,false},
    {true ,true ,true,true,false,false,false,false,false},
    {true ,true ,true ,false,false,false,false,false,false},
    {false,false,false,false,false,false,false,false,false},
    {false,false,true ,false,true ,false,false,false,false},
    {false,false,false,false,true ,true ,false,false,false},
    {false,false,false,false,false,false,true ,false,false},
    {false,false,false,false,false,false,false,true ,false},
    {false,false,false,false,false,false,false,false,false}
  };

  vector<vector< double > > dEdxStopPlusMatrix = {
    {0.9,0.9  ,0  ,0.45  ,0  ,0  ,0  ,0  ,0  },
    {0.9  ,0.9  ,0  ,0.45  ,0  ,0  ,0  ,0  ,0  },
    {3.5  ,3.5  ,3.5  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0  ,0  ,5.0  ,0  ,5.0  ,0  ,0  ,0  ,0  },
    {0  ,0  ,0  ,0  ,4.5  ,4.5  ,0  ,0  ,0  },
    {0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  }
  };
  vector<vector< double > > dEdxStopMinusMatrix = {
    {0.9  ,0.9  ,0  ,0.45  ,0  ,0  ,0  ,0  ,0  },
    {0.9  ,0.9  ,0  ,0.45  ,0  ,0  ,0  ,0  ,0  },
    {3.5  ,3.5  ,3.5  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0    ,0    ,0  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0    ,0    ,5.0  ,0  ,5.0  ,0  ,0  ,0  ,0  },
    {0    ,0    ,0  ,0  ,4.5  ,4.5  ,0  ,0  ,0  },
    {0    ,0    ,0  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0    ,0    ,0  ,0  ,0  ,0  ,0  ,0  ,0  },
    {0    ,0    ,0  ,0  ,0  ,0  ,0  ,0  ,0  }
  };

  vector<vector< double > > BTOFStopPlusMatrix = {
    {3.3,3.3,5,5,5,5,5,5,5},
    {3.3,3.3,5,5,5,5,5,5,5},
    {5,2.0,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5}
  };
  vector<vector< double > > BTOFStopMinusMatrix = {
    {3.3,3.3,5,5,5,5,5,5,5},
    {3.3,3.3,5,5,5,5,5,5,5},
    {5,2.0,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5}
  };

  // used for extrapolation and BTOF isolated solo fit corrections, even if not used for total fits
  vector<vector< double > > TPCExtrapolationPlusMatrix = {
    {0.7,0.7,5,5,5,5,5,5,5},
    {0.7,0.7,5,5,5,5,5,5,5},
    {1.05,1.05,1.05,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5}
  };
  vector<vector< double > > TPCExtrapolationMinusMatrix = {
    {0.7,0.7,5,5,5,5,5,5,5},
    {0.7,0.7,5,5,5,5,5,5,5},
    {1.05,1.05,1.05,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5}
  };

  vector<vector< double > > BTOFExtrapolationPlusMatrix = {
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {2.9, 5 ,2.9,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,4.4,5,4.4,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5}
  };
  vector<vector< double > > BTOFExtrapolationMinusMatrix = {
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {2.9,5,2.9,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,4.4,5,4.4,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5},
    {5,5,5,5,5,5,5,5,5}
  };

  CutClass* cuts = new CutClass();
  SetCutClass(cuts);

  ZFitter* fitter = new ZFitter();
  fitter->setCuts(cuts);

  fitter->setTesting(testing,divNumForTesting,mTm0DivNumForTesting);
  fitter->setTestingRapBins(rapBinsToTest);
  fitter->setParticleInfo(cuts->getStarver());
  fitter->setParticleToggles(partToggles);
  fitter->setMimicInfo(particleMimicIndicies, particleMimicMomStart, particleMimicMomStop);
  fitter->setFullMimicToggle(fullParticleMimic);
  fitter->setDoKaonDecay(doKaonDecay);
  fitter->setSeparatePM_bToF_allCent(separatePM_bToF_AllCent);
  fitter->useRootMinimizer(useRootMinimizer);
  fitter->allowCentIntegratedFitsToVary(allowCentIntegratedFitsToVary);
  fitter->allowSpecCentFitsToVary_BTOF(allowSpecCent_DeDxToVary);
  fitter->allowSpecCentFitsToVary_BTOF(allowSpecCent_ToVary_BTOF);
  fitter->useRLBiasedSingleParticleZoom(false);
  fitter->allowDeDxParameterExtrapolation(false);
  fitter->calculateBTOFIsolatedParameterDifferencesToTotalFits(false);
  fitter->usePionMinusFitAsSeedForPionPlus(usePionMinusAsSeedForPionPlus);
  fitter->useMirroredRapidityAsSeed(true);
  fitter->useStudentTDistributionsForTOF(useStudentTForTOF);
  fitter->useStudentTDistributionsForTPC(useStudentTForTPC);
  fitter->useMatrixInitializationForTPC(false); // Buggy - don't use, matches original
  fitter->drawInitialSeedsToFits(drawInitialSeedsToFits);
  fitter->usePearsonIVinSingleFits(true);
  fitter->useIdenticalRangeOfAllCent(true);
  fitter->setNdofStudentTMinimum(1.0);
  fitter->massSquaredMode(massSquaredMode);
  fitter->storeFitsNotStatusZero(true);
  fitter->setSaveNoLogImages(false);
  fitter->setConvertInvBetaToMassSquared(convertInvBetaToMassSquared);

  fitter->loadBasicDataHistograms(a_pionYieldFile, hasETOF); // PION file must have simple histos in it

  fitter->setParticleIndividualFitMomHighBounds(a_mom_dEdx_plus, a_mom_dEdx_minus, a_mom_invBeta_plus, a_mom_invBeta_minus);
  fitter->setParticleMaxMomForMeanAndSigmaParToVaryDeDx(a_mom_dEdx_maxVary_plus, a_mom_dEdx_maxVary_minus);
  fitter->setParticleFitMatrix(plusFitMatrix, minusFitMatrix);
  fitter->setPartilceMomStopMatrix(dEdxStopPlusMatrix, dEdxStopMinusMatrix, BTOFStopPlusMatrix, BTOFStopMinusMatrix);
  fitter->setPartilceExtrapolationMatrix(TPCExtrapolationPlusMatrix, TPCExtrapolationMinusMatrix, BTOFExtrapolationPlusMatrix, BTOFExtrapolationMinusMatrix);
  fitter->setPearsonZTPCMaxMomStop(mom_Pearson_plus, mom_Pearson_minus);

  //############  TPC EXTRAPOLATION RANGES ###############
  //PION PLUS
  fitter->setTPCExtrapolationRangeMatrix(0,0,1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(0,1,1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(0,3,1,0.3,0.36);
  //KAON PLUS
  fitter->setTPCExtrapolationRangeMatrix(1,0,1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(1,1,1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(1,3,1,0.3,0.36);
  //PROTON PLUS
  fitter->setTPCExtrapolationRangeMatrix(2,0,1,1.4,2.0);
  fitter->setTPCExtrapolationRangeMatrix(2,1,1,1.0,1.8);
  fitter->setTPCExtrapolationRangeMatrix(2,2,1,1.7,2.4);
  //PION MINUS
  fitter->setTPCExtrapolationRangeMatrix(0,0,-1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(0,1,-1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(1,3,1,0.3,0.36);
  //KAON MINUS
  fitter->setTPCExtrapolationRangeMatrix(1,0,-1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(1,1,-1,0.5,0.57);
  fitter->setTPCExtrapolationRangeMatrix(1,3,1,0.3,0.36);
  //PROTON MINUS
  fitter->setTPCExtrapolationRangeMatrix(2,0,-1,1.4,2.0);
  fitter->setTPCExtrapolationRangeMatrix(2,1,-1,1.0,1.8);
  fitter->setTPCExtrapolationRangeMatrix(2,2,-1,1.7,2.4);

  //############  BTOF EXTRAPOLATION RANGES ###############
  //PION PLUS
  fitter->setBTOFExtrapolationRangeMatrix(0,0,1,1.2,1.8);
  fitter->setBTOFExtrapolationRangeMatrix(0,1,1,1.2,1.8);
  //KAON PLUS
  fitter->setBTOFExtrapolationRangeMatrix(1,0,1,1.2,1.8);
  fitter->setBTOFExtrapolationRangeMatrix(1,1,1,1.2,1.8);
  //PROTON PLUS
  fitter->setBTOFExtrapolationRangeMatrix(2,0,1,1.4,2.0);
  fitter->setBTOFExtrapolationRangeMatrix(2,1,1,1.0,1.8);
  fitter->setBTOFExtrapolationRangeMatrix(2,2,1,1.7,2.4);
  //PION MINUS
  fitter->setBTOFExtrapolationRangeMatrix(0,0,-1,1.2,1.8);
  fitter->setBTOFExtrapolationRangeMatrix(0,1,-1,1.2,1.8);
  //KAON MINUS
  fitter->setBTOFExtrapolationRangeMatrix(1,0,-1,1.2,1.8);
  fitter->setBTOFExtrapolationRangeMatrix(1,1,-1,1.2,1.8);
  //PROTON MINUS
  fitter->setBTOFExtrapolationRangeMatrix(2,0,-1,1.4,2.0);
  fitter->setBTOFExtrapolationRangeMatrix(2,1,-1,1.0,1.8);
  fitter->setBTOFExtrapolationRangeMatrix(2,2,-1,1.7,2.4);

  bool loaded[9];
  for(int iii = 0; iii < 9; iii++) loaded[iii] = false;

  if(!a_useQuickLoad){
    if(a_pionYieldFile != ""){
      fitter->loadDataHistograms(a_pionYieldFile,0,nCentToRun,rapidityRebin,mTm0Rebin_byParticle[0],1,1,convert16to9,combine5and10Cent,false,1,dynamicRebinInfo_Plus,dynamicRebinInfo_Minus);
      loaded[0] = true;
    }
    if(a_kaonYieldFile != ""){
      fitter->loadDataHistograms(a_kaonYieldFile,1,nCentToRun,rapidityRebin,mTm0Rebin_byParticle[1],1,1,convert16to9,combine5and10Cent,false,1,dynamicRebinInfo_Plus,dynamicRebinInfo_Minus);
      loaded[1] = true;
    }
    if(a_protonYieldFile != ""){
      fitter->loadDataHistograms(a_protonYieldFile,2,nCentToRun,rapidityRebin,mTm0Rebin_byParticle[2],1,1,convert16to9,combine5and10Cent,false,1,dynamicRebinInfo_Plus,dynamicRebinInfo_Minus);
      loaded[2] = true;
    }
    fitter->makeAllCentHistos();

    if(a_createQuickStart){
      fitter->createQuickStartFile(a_quickStartFile);
      return;
    }
  }else{
    if(a_pionYieldFile != "")   loaded[0] = true;
    if(a_kaonYieldFile != "")   loaded[1] = true;
    if(a_protonYieldFile != "") loaded[2] = true;
    fitter->loadQuickStartFile(a_quickStartFile,9,a_partIndex,a_detectorIndex);
  }

  fitter->setImageDir(a_outputDir, imageDirName);
  fitter->setTextFileDir(textDirName);
  fitter->makeFitDataHistos(a_partIndex,a_detectorIndex);

  if(a_loadTextData){
    fitter->loadTextFileInformation(a_textFileList);
  }

  fitter->drawBichselOnDeDxPlot();
  fitter->drawInvBetaPredictionPlot();

  if(!a_loadEtaBetaGamma) fitter->fitEtaBetaGammaDeDx(a_etaBetaGammaDeDxFile,true,0.4,false);
  else fitter->loadEtaBetaGammaDeDxFile(a_etaBetaGammaDeDxFile,2,22);

  if(!a_loadTextData){
    // ORDER MATTERS! DON'T CHANGE! (BTOF fits feed seeds/constraints the TPC fits use)
    for(int partIndex = 0; partIndex < 9; partIndex++){
      if(!loaded[partIndex]) continue;
      fitter->setParticleIndex(partIndex);
      if(partIndex == 0 && fitPionBTOF){
        if(selectiveCharged >= 0) fitter->fitBTOF_SimulCent(0, 1);
        if(selectiveCharged <= 0) fitter->fitBTOF_SimulCent(0, -1);
      }
      if(partIndex == 1 && fitKaonBTOF){
        if(selectiveCharged >= 0) fitter->fitBTOF_SimulCent(1, 1);
        if(selectiveCharged <= 0) fitter->fitBTOF_SimulCent(1, -1);
      }
      if(partIndex == 2 && fitProtonBTOF){
        if(selectiveCharged >= 0) fitter->fitBTOF_SimulCent(2, 1);
        if(selectiveCharged <= 0) fitter->fitBTOF_SimulCent(2, -1);
      }
    }
  }
  cout << "FINISHED BTOF FITTING EACH CENTRALITY" << endl;

  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(!loaded[partIndex]) continue;
    fitter->setParticleIndex(partIndex);
    if(partIndex == 0 && fitPionTPC){
      if(selectiveCharged >= 0) fitter->fitTPCPionKaon_SimulCent_RapidityLooper(partIndex, 1);
      if(selectiveCharged <= 0) fitter->fitTPCPionKaon_SimulCent_RapidityLooper(partIndex, -1);
    }
    if(partIndex == 1 && fitKaonTPC){
      if(selectiveCharged >= 0) fitter->fitTPCPionKaon_SimulCent_RapidityLooper(partIndex, 1);
      if(selectiveCharged <= 0) fitter->fitTPCPionKaon_SimulCent_RapidityLooper(partIndex, -1);
    }
    if(partIndex == 2 && fitProtonTPC){
      if(selectiveCharged >= 0) fitter->fitTPCProton_SimulCent_RapidityLooper(1);
      if(selectiveCharged <= 0) fitter->fitTPCProton_SimulCent_RapidityLooper(-1);
    }
  }

  cout << "Attempting to make the spectra output file..." << endl;
  fitter->makeSpectra(a_outputDir + spectraFileName);
  cout << "Done With ZFitter..." << endl;
  return;
}
