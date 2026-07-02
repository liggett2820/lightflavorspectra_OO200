#ifndef _ZFITTER_CXX_
#define _ZFITTER_CXX_
#include "../headers/ZFitter.h"

//THIS PUTS THE OTHER FILES RIGHT HERE DONT PUT THEM IN THE MAKEFILE
//
// Ported from the original repo's 13-file #include list (matching the "one giant
// translation unit" build pattern documented in headers/ZFitter.h). Six files are
// intentionally NOT ported/included here -- see the header comment at the top of
// headers/ZFitter.h for the full reasoning and a careful cross-file dependency audit
// (several functions had to be extracted out of these excluded files and copied
// verbatim into the files below, because "live" and "dead" code turned out to be
// interleaved within the same physical files rather than cleanly separated by file):
//   ZFitter_eTOFFuncts.cxx            -- ETOF fit drivers, never called by this repo's run macro
//   ZFitter_PionKaonTPC.cxx           -- superseded by fitTPCPionKaon_SimulCent_RapidityLooper
//   ZFitter_PionKaonTPC_SimulCent.cxx -- superseded by fitTPCPionKaon_SimulCent_RapidityLooper
//                                        (its 2 still-live chi-squared functors were extracted
//                                        into ZFitter_PionKaonTPC_SimulCent_SingleLoop.cxx)
//   ZFitter_ProtonTPC.cxx             -- superseded by fitTPCProton_SimulCent_RapidityLooper
//   ZFitter_BTOF_ColliderCenter.cxx   -- superseded by fitBTOF_SimulCent
//                                        (its still-live extrapolate_BTOF_Fit_Parameters was
//                                        extracted into ZFitter_BTOF_ColliderCenter_SimulCent.cxx)
//   ZFitter_ProtonTPC_SimulCent.cxx   -- superseded by fitTPCProton_SimulCent_RapidityLooper
//   ZFitter_TPC_FXT_SimulCent.cxx     -- FixedTarget-only, not applicable to O+O 200 (ColliderCenter)
#include "../source/ZFitter_TPCFuncts.cxx"
#include "../source/ZFitter_bTOFFuncts.cxx"
#include "../source/ZFitter_dEdxFuncts.cxx"
#include "../source/ZFitter_PionKaonTPC_SimulCent_SingleLoop.cxx"
#include "../source/ZFitter_BTOF_ColliderCenter_SimulCent.cxx"
#include "../source/ZFitter_ProtonTPC_SimulCent_SingleLoop.cxx"


//=========================================================================================================================================================================
//=========================================================================================================================================================================

ZFitter::ZFitter(){
  m_partInfo = new ParticleInfo();
  m_particleTagDescriptions.push_back("Particle Not Fit (0 parameters free, 0 fixed)");                             //Particle Tag 0
  m_particleTagDescriptions.push_back("Scale Parameter to Fit. Mean and StdDev Fixed (1 free, 2 fixed)"); //Particle Tag 1
  m_particleTagDescriptions.push_back("Scale, Mean, and StdDev to Fit (3 free, 0 fixed)");                         //Particle Tag 2
  m_particleTagDescriptions.push_back("All Parameters Fixed (0 free, 3 fixed)");                                    //Particle Tag 3
  m_particleTagDescriptions.push_back("Scale constrained by BTOF ratio, Mean and StdDev mimic another particle's parameters (0 free, 3 fixed)"); //Particle Tag 4
  m_particleTagDescriptions.push_back("Scale to be Fit, Mean mimics another particle's mean, StdDev is Fixed (1 free, 2 fixed)"); //Particle Tag 5
  m_particleTagDescriptions.push_back("Scale to be Fit, Mean and StdDev mimic another particle's parameters (1 free, 2 fixed)"); //Particle Tag 6
  m_particleTagDescriptions.push_back("Scale constrained by BTOF ratio, Mean and StdDev free (2 free, 1 fixed)"); //Particle Tag 7
  m_particleTagDescriptions.push_back("Scale constrained by BTOF ratio, Mean and StdDev fixed (0 free, 3 fixed)"); //Particle Tag 8
  m_particleTagDescriptions.push_back("Scale to be Fit, Mean is fixed and StdDev is free (2 free, 1 fixed)"); //Particle Tag 9
  m_particleTagDescriptions.push_back("Scale to be Fit, Mean is fixed and StdDev mimics (1 free, 2 fixed)"); //Particle Tag 10
  m_particleTagDescriptions.push_back("Scale to be Fit, Mean is free and StdDev mimics (2 free, 1 fixed)"); //Particle Tag 11
  m_particleTagDescriptions.push_back("Scale constrained by BTOF ratio, Mean is fixed and StdDev mimics (0 free, 3 fixed)"); //Particle Tag 12
  m_particleTagDescriptions.push_back("Scale constrained by BTOF ratio, Mean and StdDev are fixed(0 free, 3 fixed)"); //Particle Tag 13

  m_isolatedZTPC_starts[0]=0.0;
  m_isolatedZTPC_starts[1]=0.0;
  m_isolatedZTPC_starts[2]=0.0;
  m_isolatedZTPC_starts[3]=0.0;
  m_isolatedZTPC_starts[4]=0.0;
  m_isolatedZTPC_starts[5]=0.0;
  m_isolatedZTPC_starts[6]=0.0;
  m_isolatedZTPC_starts[7]=0.0;
  m_isolatedZTPC_starts[8]=0.1;


  m_initialParticleToggles[0] = true;
  m_initialParticleToggles[1] = true;
  m_initialParticleToggles[2] = true;
  m_initialParticleToggles[3] = true;
  m_initialParticleToggles[4] = true;
  m_initialParticleToggles[5] = true;
  m_initialParticleToggles[6] = true;
  m_initialParticleToggles[7] = true;
  m_initialParticleToggles[8] = true;

  // -1 left, 0 middle, 1 right
  m_particleLeftRightSeedDirection[0] = -1;
  m_particleLeftRightSeedDirection[1] = 1;
  m_particleLeftRightSeedDirection[2] = 1;
  m_particleLeftRightSeedDirection[3] = 0;
  m_particleLeftRightSeedDirection[4] = 1;
  m_particleLeftRightSeedDirection[5] = 1;
  m_particleLeftRightSeedDirection[6] = 1;
  m_particleLeftRightSeedDirection[7] = 1;
  m_particleLeftRightSeedDirection[8] = 1;

  cout << " Initial Last Momentums to fit dEdx with only one gaussian: " << endl;
  for(unsigned int iii = 0; iii < 9; iii++) cout << "\t" << m_partInfo->GetParticleName(iii) << "  isolated fit stops at p= " << m_dEdxSeparatedStops[0] << endl;
  cout << endl;

  //######## SETUP PARTICLE FIT MATRIX ###################
  //NOTE: EDIT WITH setParticleFitMatrixElementPlus or Minus in Macro ONLY
  // which particles will be fit with particle of interest
  for(int iii = 0; iii < 9; iii++){
    m_particleLoaded[iii] = false;
    m_currentParticleTags[iii] = 0;
    m_partHasBeenExtrapolatedCharge[iii][0] = false;
    m_partHasBeenExtrapolatedCharge[iii][1] = false;
    m_partHasBeenExtrapolated_Solo_Charge[iii][0] = false;
    m_partHasBeenExtrapolated_Solo_Charge[iii][1] = false;
    m_dEdxSeparatedStops[iii][0] = -1;
    m_dEdxSeparatedStops[iii][1] = -1;
    m_bTOFSeparatedStops[iii][0] = -1;
    m_bTOFSeparatedStops[iii][1] = -1;
    m_eTOFSeparatedStops[iii][0] = -1;
    m_eTOFSeparatedStops[iii][1] = -1;
    m_PearsonZTPCMaxMomStop[iii][0] = -1;
    m_PearsonZTPCMaxMomStop[iii][1] = -1;

    m_dEdxMaximumMomForTotalFit[iii][0] = 999;
    m_dEdxMaximumMomForTotalFit[iii][1] = 999;
    m_electron_rap_dEdx[iii][0] = NULL;
    for(int jjj = 0; jjj < 3; jjj++){
      m_electronDeDxParamsByRapidity[iii][jjj][0] = NULL;
      m_electronDeDxParamsByRapidity[iii][jjj][1] = NULL;
      m_electronDeDxParamFitsByRapidity[iii][jjj][0] = NULL;
      m_electronDeDxParamFitsByRapidity[iii][jjj][1] = NULL;
    }

    for(int jjj = 0; jjj < 9; jjj++){
      m_partToFitMatrix_Plus[iii][jjj] = false;
      m_partToFitMatrix_Minus[iii][jjj] = false;
      m_partDeDxMomStopMatrix_Plus[iii][jjj] = 0.0;
      m_partDeDxMomStopMatrix_Minus[iii][jjj] = 0.0;
      m_partBTOFMomStopMatrix_Plus[iii][jjj] = 0.0;
      m_partBTOFMomStopMatrix_Minus[iii][jjj] = 0.0;
      m_partETOFMomStopMatrix_Plus[iii][jjj] = 0.0;
      m_partETOFMomStopMatrix_Minus[iii][jjj] = 0.0;
      m_partDeDx_MomStartExtrapolation_Matrix_Plus[iii][jjj] = 100.0; // initialize as super high momentum so default is no extrapolation
      m_partDeDx_MomStartExtrapolation_Matrix_Minus[iii][jjj] = 100.0;
      m_partBTOF_MomStartExtrapolation_Matrix_Plus[iii][jjj] = 100.0;
      m_partBTOF_MomStartExtrapolation_Matrix_Minus[iii][jjj] = 100.0;
      m_partETOF_MomStartExtrapolation_Matrix_Plus[iii][jjj] = 100.0;
      m_partETOF_MomStartExtrapolation_Matrix_Minus[iii][jjj] = 100.0;

      m_BTOF_Extrapolation_Range_Matrix[iii][jjj][0][0] = 999;
      m_BTOF_Extrapolation_Range_Matrix[iii][jjj][0][1] = 999;
      m_BTOF_Extrapolation_Range_Matrix[iii][jjj][1][0] = 999;
      m_BTOF_Extrapolation_Range_Matrix[iii][jjj][1][1] = 999;

      if(iii == jjj){
        m_partToFitMatrix_Plus[iii][jjj] = true;
        m_partToFitMatrix_Minus[iii][jjj] = true;
      }
    }
  }



  cout << "####    CURRENT POSITIVE PARTICLE FIT TOGGLES  ######" << endl;
  for(int iii = 0; iii < 9; iii++){
    cout << " While fitting " << m_partInfo->GetParticleName(iii) << " : " ;
    for(int jjj = 0; jjj < 9; jjj++){
      if(m_partToFitMatrix_Plus[iii][jjj] ) cout << m_partInfo->GetParticleName(jjj) << "  ";
    }
    cout << " will be fit" <<endl;
  }

  cout << "####    CURRENT NEGATIVE PARTICLE FIT TOGGLES  ######" << endl;
  for(int iii = 0; iii < 9; iii++){
    cout << " While fitting " << m_partInfo->GetParticleName(iii) << " : " ;
    for(int jjj = 0; jjj < 9; jjj++){
      if(m_partToFitMatrix_Minus[iii][jjj] ) cout << m_partInfo->GetParticleName(jjj) << "  ";
    }
    cout << " will be fit" <<endl;
  }




  m_lowMassCut[0]=0.0;
  m_lowMassCut[1]=0.3561;
  m_lowMassCut[2]=0.7135;
  m_lowMassCut[3]=0.0;
  m_lowMassCut[4]=1.6;
  m_lowMassCut[5]=2.7;
  m_lowMassCut[6]=2.7;
  m_lowMassCut[7]=3.6;
  m_lowMassCut[8]=0.97;

  m_highMassCut[0]=0.3561;
  m_highMassCut[1]=0.7135;
  m_highMassCut[2]=1.264;
  m_highMassCut[3]=0.065;
  m_highMassCut[4]=1.9;
  m_highMassCut[5]=2.9;
  m_highMassCut[6]=2.9;
  m_highMassCut[7]=3.8;
  m_highMassCut[8]=0.121;

  m_canvasParameters[0] = 1300;
  m_canvasParameters[1] = 800;
  m_canvasParameters[2] = 1200;
  m_canvasParameters[3] = 700;

  m_currentCentIndex = -1;
  m_currentRapidity = 0.0;
  m_currentMtM0 = 0.0;
  m_inAllCentMode = false;
  m_textFileDir = "YouDidntSetATextFileDir";
  m_logNormFunct = new TF1("m_logNormFunct",
    "[0]*0.3989422804 * TMath::Exp(-0.5*( TMath::Log(x) - [1] )*( TMath::Log(x) - [1] )/([2]*[2])) /([2] * x)",
    0,200);
  m_logNormFunct->SetNpx(500);
  m_logNormFunct->SetLineColor(kGreen);

  m_gausFunct = new TF1("m_gausFunct",
    "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",
    -200,200);
  m_gausFunct->SetNpx(500);
  m_gausFunct->SetLineColor(kGreen);

  m_student_T_Funct = new TF1("m_student_T_Funct",
    "[0]*TMath::Student((x - [1])/[2], [3])/[2]",
    -10,10);

  TVirtualFitter::SetMaxIterations( 1000000 );



  m_testing = false;
  m_divTesting = 1;
  m_divTestingMtM0 = 1;
  m_dEdxMomFittingCompleted = false;
  m_electronExtrapolationComplete = false;
  m_usePionMinusFitAsSeedForPionPlus = false;
  m_useMirroredRapidityAsSeed = false;
  m_useStudentTDistributionsForTPC = false;
  m_useStudentTDistributionsForTOF = true;
  m_useRootMinimizer = false;
  m_useGausMixModel = false;
  m_allowCentIntegratedFitsToVary = false;
  m_gausMixModel = new GausMixModel();
  m_useGeneralETOFFitDataConstraint = false;
  m_drawInitialSeedsToFits = false;
  m_bichselShiftsComplete = false;
  m_electronStudyComplete = false;
  m_allowSpecCentFitsToVary = false;
  m_allowSpecCentFitsToVary_BTOF = false;
  m_useRLBiasedSingleParticleZoom = false;
  m_usePearsonIVinSingle = false;
  m_calculateBTOFIsolatedParameterDifferencesToTotalFits = false;
  m_useIdenticalRangeOfAllCent = true;
  m_NdofStudentTMinimum = 1.0;
  m_massSquaredMode = false;
  m_useFitInfoForOppositeCharge = false;
  m_doIntegralOfBinForFits = false;
  m_storeFitsNotStatusZero = true;
  m_saveNoLogImages = false;
  m_convertInvBetaToMassSquared = false;
  m_useColliderStopTable = false;

  m_seedFromOppositeChargeDirection[0] = -1;
  for(int iii = 1; iii < 9; iii++) m_seedFromOppositeChargeDirection[iii] = 1;

  for( int partIndex = 0; partIndex < 9; partIndex++){
    m_weightedAverage_mean[partIndex] = false;
    m_weightedAverage_stdDev[partIndex] = false;
  }

  m_extrap_fn_vers = 0;
  m_basicElectronPionRatio = 0.0;
  m_basicKaonPionRatio = 0.0;
  m_currentIntegralofHistogram = 0.0;

  m_betaGammaBTOF_dEdx_FitFunction[0] = NULL;
  m_betaGammaBTOF_dEdx_FitFunction[1] = NULL;

  for(int partIndex = 0; partIndex<9; partIndex++){
    m_particleMimicIndicies[partIndex] = 0;
    m_particleMimicLowMom[partIndex] = 0;
    m_particleMimicHighMom[partIndex] = 0;
    m_fullParticleMimic[partIndex] = true;
    m_currentConstrainedAmpRaios[partIndex] = 0.0;
  }

  m_colliderStopTable = {
    //particle rapidity bTofStartMom easyDeDxMom hardDeDxMom doStep0and1
    {0, 1, 0.0, 0.24, 0.40, 0.60, 1},
    {0, 1, 0.1, 0.24, 0.40, 0.60, 1},
    {0, 1, 0.2, 0.24, 0.40, 0.60, 1},
    {0, 1, 0.3, 0.24, 0.40, 0.60, 1},
    {0, 1, 0.4, 0.24, 0.40, 0.60, 1},
    {0, 1, 0.5, 0.24, 0.41, 0.60, 1},
    {0, 1, 0.6, 0.24, 0.43, 0.60, 1},      // last rap for left side electron fit
    {0, 1, 0.7, 0.25, 0.43, 0.60, 1},
    {0, 1, 0.8, 0.26, 0.43, 0.60, 1},
    {0, 1, 0.9, 0.34, 0.44, 0.60, 1},
    {0, 1, 1.0, 0.42, 0.44, 0.60, 1},
    {0, 1, 1.1, 0.65, 0.44, 0.70, 1},
    {0, 1, 1.2, 3.00, 0.45, 0.75, 1},
    {0, 1, 1.3, 3.00, 0.45, 0.75, 1},
    {0, 1, 1.4, 3.00, 0.45, 0.75, 1},

    {0, -1, 0.0, 0.24, 0.40, 0.60, 1},
    {0, -1, 0.1, 0.24, 0.40, 0.60, 1},
    {0, -1, 0.2, 0.24, 0.40, 0.60, 1},
    {0, -1, 0.3, 0.24, 0.40, 0.60, 1},
    {0, -1, 0.4, 0.24, 0.40, 0.60, 1},
    {0, -1, 0.5, 0.24, 0.41, 0.60, 1},
    {0, -1, 0.6, 0.24, 0.43, 0.60, 1},      // last rap for left side electron fit
    {0, -1, 0.7, 0.25, 0.43, 0.60, 1},
    {0, -1, 0.8, 0.26, 0.43, 0.60, 1},
    {0, -1, 0.9, 0.34, 0.44, 0.60, 1},
    {0, -1, 1.0, 0.42, 0.44, 0.60, 1},
    {0, -1, 1.1, 0.65, 0.44, 0.70, 1},
    {0, -1, 1.2, 3.00, 0.45, 0.75, 1},
    {0, -1, 1.3, 3.00, 0.45, 0.75, 1},
    {0, -1, 1.4, 3.00, 0.45, 0.75, 1},


    //particle rapidity bTofStartMom easyDeDxMom hardDeDxMom doStep0and1

    // easy dEdx can fit electrons without constraints
    {1, 1, 0.0, 0.27, 0.42, 0.5, 0 },
    {1, 1, 0.1, 0.27, 0.41, 0.5, 0 },
    {1, 1, 0.2, 0.28, 0.41, 0.5, 0 },
    {1, 1, 0.3, 0.29, 0.41, 0.5, 0 },
    {1, 1, 0.4, 0.30, 0.41, 0.5, 0 },
    {1, 1, 0.5, 0.36, 0.42, 0.5, 0 },
    {1, 1, 0.6, 0.45, 0.43, 0.5, 1 },
    {1, 1, 0.7, 0.54, 0.44, 0.5, 1 },       // this is more difficult to do ratio method
    {1, 1, 0.8, 0.68, 0.44, 0.5, 1 },
    {1, 1, 0.9, 0.86, 0.44, 0.5, 1 },
    {1, 1, 1.0, 1.10, 0.44, 0.5, 1 },
    {1, 1, 1.1, 2.10, 0.44, 0.5, 1 },

    {1, -1, 0.0, 0.27, 0.42, 0.5, 0 },
    {1, -1, 0.1, 0.27, 0.41, 0.5, 0 },
    {1, -1, 0.2, 0.28, 0.41, 0.5, 0 },
    {1, -1, 0.3, 0.29, 0.41, 0.5, 0 },
    {1, -1, 0.4, 0.30, 0.41, 0.5, 0 },
    {1, -1, 0.5, 0.36, 0.42, 0.5, 0 },
    {1, -1, 0.6, 0.45, 0.43, 0.5, 1 },
    {1, -1, 0.7, 0.54, 0.44, 0.5, 1 },       // this is more difficult to do ratio method
    {1, -1, 0.8, 0.68, 0.44, 0.5, 1 },
    {1, -1, 0.9, 0.86, 0.44, 0.5, 1 },
    {1, -1, 1.0, 1.10, 0.44, 0.5, 1 },
    {1, -1, 1.1, 2.10, 0.44, 0.5, 1 },

    // easy dEdx can be fit with single bump
    {2, 1,  0.0, 0.35, 0.65, 1.00,0.45},
    {2, 1,  0.1, 0.35, 0.65, 1.00,0.45},
    {2, 1,  0.2, 0.36, 0.67, 1.00,0.45},
    {2, 1,  0.3, 0.41, 0.70, 0.92,0.5},
    {2, 1,  0.4, 0.50, 0.70, 0.92,0.6},
    {2, 1,  0.5, 0.65, 0.62, 0.92,0.7},
    {2, 1,  0.6, 0.81, 0.40, 0.92,0.9},
    {2, 1,  0.7, 1.00, 0.40, 1.03,1.1},
    {2, 1,  0.8, 1.27, 0.40, 1.03,1.1},
    {2, 1,  0.9, 1.62, 0.40, 1.03,1.1},
    {2, 1,  1.0, 2.20, 0.40, 0.50,1.1},
    {2, 1,  1.1, 4.00, 0.40, 0.50,1.1},

    {2, -1, 0.0, 0.35, 0.65, 1.00,0.45},
    {2, -1, 0.1, 0.35, 0.65, 1.00,0.45},
    {2, -1, 0.2, 0.36, 0.67, 1.00,0.45},
    {2, -1, 0.3, 0.41, 0.70, 0.92,0.5},
    {2, -1, 0.4, 0.50, 0.70, 0.92,0.6},
    {2, -1, 0.5, 0.65, 0.62, 0.92,0.7},
    {2, -1, 0.6, 0.81, 0.40, 0.92,0.9},
    {2, -1, 0.7, 1.00, 0.40, 1.03,1.1},
    {2, -1, 0.8, 1.27, 0.40, 1.03,1.1},
    {2, -1, 0.9, 1.62, 0.40, 1.03,1.1},
    {2, -1, 1.0, 2.20, 0.40, 0.50,1.1},
    {2, -1, 1.1, 4.00, 0.40, 0.50,1.1}

  };


  for(int iii = 0; iii < 9; iii++){
    for(int jjj = 0; jjj < 5; jjj++){
      m_Background_FitData_ZTPC[iii][jjj][0] = NULL;
      m_Background_FitData_ZTPC[iii][jjj][1] = NULL;
      m_Background_FitData_ZbTOF[iii][jjj][0] = NULL;
      m_Background_FitData_ZbTOF[iii][jjj][1] = NULL;
      m_Background_FitData_ZeTOF[iii][jjj][0] = NULL;
      m_Background_FitData_ZeTOF[iii][jjj][1] = NULL;
    }
  }

  for(int iii = 0; iii < 3; iii++) for(int jjj = 0; jjj < 16; jjj++) m_electronRatioFitHistos_ByCent[jjj][iii] = NULL;

  for(int iii = 0; iii < 1000; iii++){
    m_rebinZTPCHistoStorage[iii] = NULL;
    m_rebinZTPCHistoIntegralStorage[iii] = 1.0;
    for(int centIndex = 0; centIndex < 16; centIndex++){
      m_rebinZTPCHistoStorage_ByCent[centIndex][iii] = NULL;
    }
  }

  for(int partBaseIndex = 0; partBaseIndex < 9; partBaseIndex++){
    for(int partIndex = 0; partIndex < 9; partIndex++){
      for(int centIndex = 0; centIndex < 16; centIndex++){
        m_dEdxYieldsCentPlus[partBaseIndex][partIndex][centIndex] = NULL;
        m_dEdxYieldsCentMinus[partBaseIndex][partIndex][centIndex] = NULL;
        m_InvBetaBTOFYieldsCentPlus[partBaseIndex][partIndex][centIndex] = NULL;
        m_InvBetaBTOFYieldsCentMinus[partBaseIndex][partIndex][centIndex] = NULL;
      }
    }
  }


  for(int partIndex = 0; partIndex < 9; partIndex++){

    m_centEvents[partIndex] = NULL;
    m_rap_mTm0_dEdx_allCent_PionElectron[partIndex] = NULL; //by particle rapidity space
    m_rap_mTm0_dEdx_AllCentRead_Plus[partIndex]  = NULL;
    m_rap_mTm0_dEdx_AllCentRead_Minus[partIndex] = NULL;
    m_mom_dEdx_BTOFIsolated_Plus[partIndex]  = NULL; // by particle
    m_mom_dEdx_BTOFIsolated_Minus[partIndex] = NULL; // by particle


    for(int iii = 0; iii < 3; iii++){
      m_mom_dEdxFitData[partIndex][iii] = NULL;
      m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][iii] = NULL;
      m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][iii] = NULL;
      m_mom_invBetaBTOFFitData[partIndex][iii] = NULL;
      m_mom_invBetaETOFFitData[partIndex][iii] = NULL;
    }


    m_ZTPC_LowRange[partIndex][0] = NULL;
    m_ZTPC_HighRange[partIndex][0] = NULL;
    m_ZbTOF_LowRange[partIndex][0] = NULL;
    m_ZbTOF_HighRange[partIndex][0] = NULL;

    m_ZTPC_LowRange[partIndex][1] = NULL;
    m_ZTPC_HighRange[partIndex][1] = NULL;
    m_ZbTOF_LowRange[partIndex][1] = NULL;
    m_ZbTOF_HighRange[partIndex][1] = NULL;


    for(int partIndex2=0; partIndex2<9; partIndex2++){
      m_ZbTOF_partBinning_LowRange[partIndex][partIndex2][0]  = NULL;
      m_ZbTOF_partBinning_HighRange[partIndex][partIndex2][0] = NULL;
      m_ZbTOF_partBinning_LowRange[partIndex][partIndex2][1]  = NULL;
      m_ZbTOF_partBinning_HighRange[partIndex][partIndex2][1] = NULL;
    }

    m_rap_mTm0_dEdx_Plus_allCent[partIndex] = NULL;
    m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex] = NULL;
    m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex] = NULL;
    m_rap_mTm0_dEdx_Minus_allCent[partIndex] = NULL;
    m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex] = NULL;
    m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex] = NULL;
    m_rap_mTm0_dEdx_allCent[partIndex] = NULL;
    m_rap_mTm0_invBetaBTOF_allCent[partIndex] = NULL;
    m_rap_mTm0_invBetaETOF_allCent[partIndex] = NULL;


    for(int centIndex = 0; centIndex < 16; centIndex++){
      m_currentHistosToFit_ByCent[centIndex] = NULL;

      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        m_ChiSqr_Cent_ZTPC[partIndex][centIndex][pmIndex] = NULL;
        m_ChiSqr_Cent_ZbTOF[partIndex][centIndex][pmIndex] = NULL;
        m_ChiSqr_Cent_ZeTOF[partIndex][centIndex][pmIndex] = NULL;
        m_FitStatus_Cent_ZTPC[partIndex][centIndex][pmIndex] = NULL;
        m_FitStatus_Cent_ZbTOF[partIndex][centIndex][pmIndex] = NULL;
        m_FitStatus_Cent_ZeTOF[partIndex][centIndex][pmIndex] = NULL;
        m_Spectra_Cent_ZTPC[partIndex][centIndex][pmIndex] = NULL;
        m_Spectra_Cent_ZbTOF[partIndex][centIndex][pmIndex] = NULL;
        m_Spectra_Cent_ZeTOF[partIndex][centIndex][pmIndex] = NULL;
        for(int rap_index = 0; rap_index < g_ZFitter_MaxNumRapBins; rap_index++){
          m_mTm0_ZTPC_by_Rapidity[partIndex][centIndex][rap_index][pmIndex] = NULL;
          m_mTm0_ZBTOF_by_Rapidity[partIndex][centIndex][rap_index][pmIndex] = NULL;
          m_mTm0_ZETOF_by_Rapidity[partIndex][centIndex][rap_index][pmIndex] = NULL;
        }
      }

      m_rap_mTm0_dEdx_Plus[partIndex][centIndex] = NULL;
      m_rap_mTm0_invBetaBTOF_Plus[partIndex][centIndex] = NULL;
      m_rap_mTm0_invBetaETOF_Plus[partIndex][centIndex] = NULL;
      m_rap_mTm0_dEdx_Minus[partIndex][centIndex] = NULL;
      m_rap_mTm0_invBetaBTOF_Minus[partIndex][centIndex] = NULL;
      m_rap_mTm0_invBetaETOF_Minus[partIndex][centIndex] = NULL;



    }

    for(int iii = 0; iii < 5; iii++){
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        m_Background_FitData_ZTPC[partIndex][iii][pmIndex] = NULL;
        m_Background_FitData_ZbTOF[partIndex][iii][pmIndex] = NULL;
        m_Background_FitData_ZeTOF[partIndex][iii][pmIndex] = NULL;
        for(int partIndex2=0; partIndex2<9; partIndex2++){
          m_Background_FitData_partBinning_ZbTOF[partIndex][partIndex2][iii][pmIndex] = NULL;
        }
      }
    }
    for(int pmIndex = 0; pmIndex < 2; pmIndex++){
      m_ChiSqr_ZTPC[partIndex][pmIndex] = NULL;
      m_ChiSqr_ZbTOF[partIndex][pmIndex] = NULL;
      m_ChiSqr_ZeTOF[partIndex][pmIndex] = NULL;
      m_FitStatus_ZTPC[partIndex][pmIndex] = NULL;
      m_FitStatus_ZbTOF[partIndex][pmIndex] = NULL;
      m_FitStatus_ZeTOF[partIndex][pmIndex] = NULL;
    }

    for(int iii = 0; iii < 9; iii++){
      m_bichselShifts_Plus[partIndex][iii] = NULL;
      m_bichselShifts_Minus[partIndex][iii] = NULL;
      m_tpcResolution_Plus[partIndex][iii] = NULL;
      m_tpcResolution_Minus[partIndex][iii] = NULL;
      m_quickYieldNums_Plus[partIndex][iii] = NULL;
      m_quickYieldNums_Minus[partIndex][iii] = NULL;

      m_rap_mTm0_dEdx_Plus_allCentIsolated[partIndex][iii] = NULL;
      m_rap_mTm0_dEdx_Minus_allCentIsolated[partIndex][iii] = NULL;

      for(int pmIndex = 0; pmIndex < 2; pmIndex++){

         // partIndex for rapidity space, particle fit, amp mean sig ndf, charge 0=plus 1-minus
        m_Fit_Data_Solo_ZTPCShifts[partIndex][iii][pmIndex] = NULL;// partIndex for rapidity space, particle fit, charge 0=plus 1-minus

        for(int jjj = 0; jjj < 5; jjj++){
          m_Fit_Data_Solo_ZTPC[partIndex][iii][jjj][pmIndex] = NULL;// partIndex for rapidity space, particle fit, amp mean sig ndf, charge 0=plus 1-minus
          m_Fit_Data_ZbTOF[partIndex][iii][jjj][pmIndex] = NULL; // particle rap space, particle , parameter, charge 0=plus 1-minus
          m_Fit_Data_ZeTOF[partIndex][iii][jjj][pmIndex] = NULL; // particle rap space, particle , parameter, charge 0=plus 1-minus
          m_Fit_Data_ZTPC[partIndex][iii][jjj][pmIndex] = NULL;
          for(int cent_index = 0; cent_index < 9; cent_index++){
            m_Fit_Data_Cent_ZTPC[partIndex][iii][cent_index][jjj][pmIndex] = NULL;
            m_Fit_Data_Cent_ZbTOF[partIndex][iii][cent_index][jjj][pmIndex] = NULL;
            m_Fit_Data_Cent_ZeTOF[partIndex][iii][cent_index][jjj][pmIndex] = NULL;

          }

        }
        m_Fit_Data_ZTPC[partIndex][iii][4][pmIndex] = NULL;

      }

    }


  }//end partIndex loop

  for(int iii = 0; iii < 3; iii++){
    m_betaGammaBTOF_DeDx_Data_Plus[iii] = NULL;
    m_betaGammaBTOF_DeDx_Data_Minus[iii] = NULL;
    m_betaGammaETOF_DeDx_Data_Plus[iii] = NULL;
    m_betaGammaETOF_DeDx_Data_Minus[iii] = NULL;
    m_betaGammaBTOF_DeDx_Plus_SmoothFitData[iii] = NULL;
    m_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii] = NULL;
    m_betaGammaETOF_DeDx_Plus_SmoothFitData[iii] = NULL;
    m_betaGammaETOF_DeDx_Minus_SmoothFitData[iii] = NULL;
    m_eta_betaGammaBTOF_DeDx_Data_Plus[iii] = NULL;
    m_eta_betaGammaBTOF_DeDx_Data_Minus[iii] = NULL;
    m_eta_betaGammaETOF_DeDx_Data_Plus[iii] = NULL;
    m_eta_betaGammaETOF_DeDx_Data_Minus[iii] = NULL;
  }



  m_betaGammaBTOF_dEdx_FitFunction[0] = NULL;
  m_betaGammaETOF_dEdx_FitFunction[0] = NULL;
  m_betaGammaBTOF_dEdx_FitFunction[1] = NULL;
  m_betaGammaETOF_dEdx_FitFunction[1] = NULL;



  //For BetaGamma dEdx Fitting //structure: N, mu, sig

  m_mom_dEdx_Plus = NULL;
  m_mom_dEdx_Minus = NULL;
  m_betaGammaBTOF_dEdx_Plus = NULL;
  m_betaGammaBTOF_dEdx_Minus = NULL;
  m_betaGammaETOF_dEdx_Plus = NULL;
  m_betaGammaETOF_dEdx_Minus = NULL;
  m_eta_betaGammaBtof_dEdx_Plus = NULL;
  m_eta_betaGammaBtof_dEdx_Minus = NULL;
  m_eta_betaGammaEtof_dEdx_Plus = NULL;
  m_eta_betaGammaEtof_dEdx_Minus = NULL;
  m_mom_invBetaBTOF_Plus = NULL;
  m_mom_invBetaBTOF_Minus = NULL;
  m_mom_invBetaETOF_Plus = NULL;
  m_mom_invBetaETOF_Minus = NULL;
  m_electronReduct = NULL;

  m_doKaonDecay=false;


  for(int centIndex = 0; centIndex < 16; centIndex++){
    m_electronRatioFitHistos_ByCent[centIndex][0] = NULL;
    m_electronRatioFitHistos_ByCent[centIndex][1] = NULL;
    m_electronRatioFitHistos_ByCent[centIndex][2] = NULL;
  }


  #ifdef _ZFITTER_DEBUG_
  cout << " Done with ZFitter Constructor" << endl;
  #endif
}

//=========================================================================================================================================================================
//=========================================================================================================================================================================


void ZFitter::setParticleIndividualFitMomHighBounds(vector<double> a_mom_dEdx_plus, vector<double> a_mom_dEdx_minus, vector<double> a_mom_invBeta_plus, vector<double> a_mom_invBeta_minus){
  for(unsigned int iii = 0; iii < a_mom_dEdx_plus.size(); iii++){
    m_dEdxSeparatedStops[iii][0] = a_mom_dEdx_plus[iii];
  }
  for(unsigned int iii = 0; iii < a_mom_dEdx_minus.size(); iii++){
    m_dEdxSeparatedStops[iii][1] = a_mom_dEdx_minus[iii];
  }
  for(unsigned int iii = 0; iii < a_mom_invBeta_plus.size(); iii++){
    m_bTOFSeparatedStops[iii][0] = a_mom_invBeta_plus[iii];
  }
  for(unsigned int iii = 0; iii < a_mom_invBeta_minus.size(); iii++){
    m_bTOFSeparatedStops[iii][1] = a_mom_invBeta_minus[iii];
  }
  return;
}

void ZFitter::setParticleIndividualFitMomHighBoundsETOF(vector<double> a_mom_invBeta_plus, vector<double> a_mom_invBeta_minus){
  for(unsigned int iii = 0; iii < a_mom_invBeta_plus.size(); iii++){
    m_eTOFSeparatedStops[iii][0] = a_mom_invBeta_plus[iii];
  }
  for(unsigned int iii = 0; iii < a_mom_invBeta_minus.size(); iii++){
    m_eTOFSeparatedStops[iii][1] = a_mom_invBeta_minus[iii];
  }
  return;
}



void ZFitter::setParticleMaxMomForMeanAndSigmaParToVaryDeDx(vector<double> a_mom_dEdx_plus, vector<double> a_mom_dEdx_minus){
  for(unsigned int iii = 0; iii < a_mom_dEdx_plus.size(); iii++){
    m_dEdxMaximumMomForTotalFit[iii][0] = a_mom_dEdx_plus[iii];
  }
  for(unsigned int iii = 0; iii < a_mom_dEdx_minus.size(); iii++){
    m_dEdxMaximumMomForTotalFit[iii][1] = a_mom_dEdx_minus[iii];
  }
}

void ZFitter::setPearsonZTPCMaxMomStop(vector<double> a_mom_plus, vector<double> a_mom_minus){
  for(unsigned int iii = 0; iii < a_mom_plus.size(); iii++){
    m_PearsonZTPCMaxMomStop[iii][0] = a_mom_plus[iii];
    m_PearsonZTPCMaxMomStop[iii][1] = a_mom_minus[iii];
  }
}

void ZFitter::setParticleFitMatrix(vector< vector< bool > > a_plusMatrix, vector< vector< bool > > a_minusMatrix){
  for(unsigned int iii = 0; iii < 9; iii++){
    for(unsigned int jjj = 0; jjj < 9; jjj++){
       m_partToFitMatrix_Plus[iii][jjj] = a_plusMatrix[iii][jjj];
       m_partToFitMatrix_Minus[iii][jjj] = a_minusMatrix[iii][jjj];
    }
  }
}

void ZFitter::setParticleFitMatrix(bool a_plusMatrix[9][9], bool a_minusMatrix[9][9]){
  for(unsigned int iii = 0; iii < 9; iii++){
    for(unsigned int jjj = 0; jjj < 9; jjj++){
       m_partToFitMatrix_Plus[iii][jjj] = a_plusMatrix[iii][jjj];
       m_partToFitMatrix_Minus[iii][jjj] = a_minusMatrix[iii][jjj];
    }
  }
}

void ZFitter::setPartilceMomStopMatrix(vector< vector< double > > a_plusMatrix_TPC, vector< vector< double > > a_minusMatrix_TPC, vector< vector< double > > a_plusMatrix_BTOF, vector< vector< double > > a_minusMatrix_BTOF){
  for(unsigned int iii = 0; iii < 9; iii++){
    for(unsigned int jjj = 0; jjj < 9; jjj++){
       m_partDeDxMomStopMatrix_Plus[iii][jjj] = a_plusMatrix_TPC[iii][jjj];
       m_partDeDxMomStopMatrix_Minus[iii][jjj] = a_minusMatrix_TPC[iii][jjj];
       m_partBTOFMomStopMatrix_Plus[iii][jjj] = a_plusMatrix_BTOF[iii][jjj];
       m_partBTOFMomStopMatrix_Minus[iii][jjj] = a_minusMatrix_BTOF[iii][jjj];
    }
  }
}

void ZFitter::setPartilceMomStopMatrixETOF(vector< vector< double > > a_plusMatrix_ETOF, vector< vector< double > > a_minusMatrix_ETOF){
  for(unsigned int iii = 0; iii < 9; iii++){
    for(unsigned int jjj = 0; jjj < 9; jjj++){
       m_partETOFMomStopMatrix_Plus[iii][jjj] = a_plusMatrix_ETOF[iii][jjj];
       m_partETOFMomStopMatrix_Minus[iii][jjj] = a_minusMatrix_ETOF[iii][jjj];
    }
  }
}

void ZFitter::setPartilceExtrapolationMatrix(vector< vector< double > > a_plusMatrix_TPC, vector< vector< double > > a_minusMatrix_TPC, vector< vector< double > > a_plusMatrix_BTOF, vector< vector< double > > a_minusMatrix_BTOF){
  for(unsigned int iii = 0; iii < 9; iii++){
    for(unsigned int jjj = 0; jjj < 9; jjj++){
       m_partDeDx_MomStartExtrapolation_Matrix_Plus[iii][jjj] = a_plusMatrix_TPC[iii][jjj];
       m_partDeDx_MomStartExtrapolation_Matrix_Minus[iii][jjj] = a_minusMatrix_TPC[iii][jjj];
       m_partBTOF_MomStartExtrapolation_Matrix_Plus[iii][jjj] = a_plusMatrix_BTOF[iii][jjj];
       m_partBTOF_MomStartExtrapolation_Matrix_Minus[iii][jjj] = a_minusMatrix_BTOF[iii][jjj];
     }
   }
}

void ZFitter::setPartilceExtrapolationMatrixETOF(vector< vector< double > > a_plusMatrix_ETOF, vector< vector< double > > a_minusMatrix_ETOF){
  for(unsigned int iii = 0; iii < 9; iii++){
    for(unsigned int jjj = 0; jjj < 9; jjj++){
       m_partETOF_MomStartExtrapolation_Matrix_Plus[iii][jjj] = a_plusMatrix_ETOF[iii][jjj];
       m_partETOF_MomStartExtrapolation_Matrix_Minus[iii][jjj] = a_minusMatrix_ETOF[iii][jjj];
     }
   }
}


void ZFitter::setBTOFExtrapolationRangeMatrix(int a_mainPartIndex, int a_subPartIndex, int a_charge, double a_lowMom, double a_highMom){
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;
  m_BTOF_Extrapolation_Range_Matrix[a_mainPartIndex][a_subPartIndex][0][pmIndex] = a_lowMom;
  m_BTOF_Extrapolation_Range_Matrix[a_mainPartIndex][a_subPartIndex][1][pmIndex] = a_highMom;

}

void ZFitter::setTPCExtrapolationRangeMatrix(int a_mainPartIndex, int a_subPartIndex, int a_charge, double a_lowMom, double a_highMom){
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;
  m_TPC_Extrapolation_Range_Matrix[a_mainPartIndex][a_subPartIndex][0][pmIndex] = a_lowMom;
  m_TPC_Extrapolation_Range_Matrix[a_mainPartIndex][a_subPartIndex][1][pmIndex] = a_highMom;
}

void ZFitter::setCanvasParameters(int a_windowWidth,int a_windowHeight, int a_canvasWidth, int a_canvasHeight){
  m_canvasParameters[0] = a_windowWidth;
  m_canvasParameters[1] = a_windowHeight;
  m_canvasParameters[2] = a_canvasWidth;
  m_canvasParameters[3] = a_canvasHeight;
}



void ZFitter::setParticleToggles(bool *partToggles){
  int nParticles = 0;
  for(int iii=0; iii<9; iii++){
    m_initialParticleToggles[iii] = partToggles[iii];
    if(m_initialParticleToggles[iii]) nParticles = iii + 1;
  }
  m_numParticles = nParticles;
  if(m_numParticles < 4) m_numParticles = 4;

  return;
}

//=========================================================================================================================================================================
//=========================================================================================================================================================================


void ZFitter::loadBasicDataHistograms(string a_filename,  bool a_hasETOF){
  cout << "Start of loadBasicDataHistograms(" << a_filename << "," << a_hasETOF << ")" << endl;
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "Warning in ZFitter::loadDataHistograms.. couldn't find file" << endl;
    return;
  }
  m_hasETOF = a_hasETOF;

  TH1I* centEvents = (TH1I*) inFile->Get("centEvents");
  if(centEvents){ //if every particle was done at once in a single root file
    m_centEvents[0] = centEvents;
    m_centEvents[1] = centEvents;
    m_centEvents[2] = centEvents;
    m_centEvents[3] = centEvents;
    m_centEvents[4] = centEvents;
    m_centEvents[5] = centEvents;
    m_centEvents[6] = centEvents;
    m_centEvents[7] = centEvents;
    m_centEvents[8] = centEvents;
  }

  m_mom_dEdx_Plus  = (TH2I*) inFile->Get("dEdxVsMom_Plus");
  m_mom_dEdx_Minus = (TH2I*) inFile->Get("dEdxVsMom_Minus");
  m_electronReduct = (TH2I*) inFile->Get("electronSubtract_dEdx");

  m_betaGammaBTOF_dEdx_Plus   = (TH2I*) inFile->Get("betaGammaBtof_dEdx_Plus");
  m_betaGammaBTOF_dEdx_Minus  =  (TH2I*) inFile->Get("betaGammaBtof_dEdx_Minus");
  if(m_hasETOF) m_betaGammaETOF_dEdx_Plus  = (TH2I*) inFile->Get("betaGammaEtof_dEdx_Plus");
  if(m_hasETOF) m_betaGammaETOF_dEdx_Minus = (TH2I*) inFile->Get("betaGammaEtof_dEdx_Minus");

  m_eta_betaGammaBtof_dEdx_Plus = (TH3I*) inFile->Get("eta_betaGammaBtof_dEdx_Plus");
  m_eta_betaGammaBtof_dEdx_Minus = (TH3I*) inFile->Get("eta_betaGammaBtof_dEdx_Minus");
  if(m_hasETOF) m_eta_betaGammaEtof_dEdx_Plus  = (TH3I*) inFile->Get("eta_betaGammaEtof_dEdx_Plus");
  if(m_hasETOF) m_eta_betaGammaEtof_dEdx_Minus = (TH3I*) inFile->Get("eta_betaGammaEtof_dEdx_Minus");

  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex ==3 || partIndex == 8) continue;
    m_mom_dEdx_BTOFIsolated_Plus[partIndex] = (TH2I*) inFile->Get(Form("dEdx_%s_Isolated",
     m_partInfo->GetParticleName(partIndex,1).Data()) );
    m_mom_dEdx_BTOFIsolated_Minus[partIndex] = (TH2I*) inFile->Get(Form("dEdx_%s_Isolated",
     m_partInfo->GetParticleName(partIndex,-1).Data()));
  }

  m_mom_invBetaBTOF_Plus   = (TH2I*) inFile->Get("invBeta_bTOF_plus");
  m_mom_invBetaBTOF_Minus  = (TH2I*) inFile->Get("invBeta_bTOF_minus");
  m_mom_invBetaETOF_Plus   = (TH2I*) inFile->Get("invBeta_eTOF_plus");
  m_mom_invBetaETOF_Minus  = (TH2I*) inFile->Get("invBeta_eTOF_minus");
  #ifdef _ZFITTER_DEBUG_
    cout << "End of loadBasicDataHistograms(" << a_filename << "," << a_hasETOF << ")" << endl;
  #endif
}


//=========================================================================================================================================================================
//=========================================================================================================================================================================


//import the output of EtaPtBinner
void ZFitter::loadDataHistograms(string a_filename, int a_partIndex, int a_numCent, int a_rebinX, int a_rebinY, int a_rebinZTPC, int a_rebinZBTOF,  bool a_convert16to9CentBins, bool a_combine5and10CentBins, bool a_combine40and50_60and70_centBins,int a_rebinZETOF, vector< vector<double> > a_dynamicRebin_Plus,  vector< vector<double> > a_dynamicRebin_Minus){
  cout << "Start of loadDataHistograms(" << a_filename << "," << a_partIndex << "," << a_numCent << "," << a_rebinX << "," << a_rebinY << "," << a_rebinZTPC<< "," << a_rebinZBTOF << ")" << endl;
  m_numCentralities = a_numCent;
  m_particleLoaded[a_partIndex] = true;


  //################  FOR DYNAMIC REBIN / REMOVE UNNECESSARY INFORMATIONFROM VECTOR  ####################
  vector< vector< double > > rebinInfo_TPC_Plus   = a_dynamicRebin_Plus;
  vector< vector< double > > rebinInfo_bTOF_Plus  = a_dynamicRebin_Plus;
  vector< vector< double > > rebinInfo_eTOF_Plus  = a_dynamicRebin_Plus;
  vector< vector< double > > rebinInfo_TPC_Minus  = a_dynamicRebin_Minus;
  vector< vector< double > > rebinInfo_bTOF_Minus = a_dynamicRebin_Minus;
  vector< vector< double > > rebinInfo_eTOF_Minus = a_dynamicRebin_Minus;

  //  ######## Remove unwanted particles
  for(int rebinIndex = rebinInfo_TPC_Plus.size() - 1; rebinIndex >= 0 && rebinIndex < rebinInfo_TPC_Plus.size(); rebinIndex--){
    int attemptPartIndex = (int) rebinInfo_TPC_Plus[rebinIndex][0];
    if(attemptPartIndex != a_partIndex){
      rebinInfo_TPC_Plus.erase(rebinInfo_TPC_Plus.begin() + rebinIndex);
      rebinInfo_bTOF_Plus.erase(rebinInfo_bTOF_Plus.begin() + rebinIndex);
      rebinInfo_eTOF_Plus.erase(rebinInfo_eTOF_Plus.begin() + rebinIndex);
    }
  }
  for(int rebinIndex = rebinInfo_TPC_Minus.size() - 1; rebinIndex >= 0 && rebinIndex < rebinInfo_TPC_Minus.size(); rebinIndex--){
    int attemptPartIndex = (int) rebinInfo_TPC_Minus[rebinIndex][0];
    if(attemptPartIndex != a_partIndex){
      rebinInfo_TPC_Minus.erase(rebinInfo_TPC_Minus.begin() + rebinIndex);
      rebinInfo_bTOF_Minus.erase(rebinInfo_bTOF_Minus.begin() + rebinIndex);
      rebinInfo_eTOF_Minus.erase(rebinInfo_eTOF_Minus.begin() + rebinIndex);
    }
  }

  //  ######## Remove unwanted detectors
  for(int rebinIndex = rebinInfo_TPC_Plus.size() - 1; rebinIndex >= 0; rebinIndex--){
    int detIndex = (int) rebinInfo_TPC_Plus[rebinIndex][1];
    if(detIndex != 0){
      rebinInfo_TPC_Plus.erase(rebinInfo_TPC_Plus.begin() + rebinIndex);
    }else{
      rebinInfo_TPC_Plus[rebinIndex].erase(rebinInfo_TPC_Plus[rebinIndex].begin()); //partIndex
      rebinInfo_TPC_Plus[rebinIndex].erase(rebinInfo_TPC_Plus[rebinIndex].begin()); //detIndex
    }
  }
  for(int rebinIndex = rebinInfo_bTOF_Plus.size() - 1; rebinIndex >= 0; rebinIndex--){
    int detIndex = (int) rebinInfo_bTOF_Plus[rebinIndex][1];
    if(detIndex != 1){
      rebinInfo_bTOF_Plus.erase(rebinInfo_bTOF_Plus.begin() + rebinIndex);
    }else{
      rebinInfo_bTOF_Plus[rebinIndex].erase(rebinInfo_bTOF_Plus[rebinIndex].begin()); //partIndex
      rebinInfo_bTOF_Plus[rebinIndex].erase(rebinInfo_bTOF_Plus[rebinIndex].begin()); //detIndex
    }
  }
  for(int rebinIndex = rebinInfo_eTOF_Plus.size() - 1; rebinIndex >= 0; rebinIndex--){
    int detIndex = (int) rebinInfo_eTOF_Plus[rebinIndex][1];
    if(detIndex != 2){
      rebinInfo_eTOF_Plus.erase(rebinInfo_eTOF_Plus.begin() + rebinIndex);
    }else{
      rebinInfo_eTOF_Plus[rebinIndex].erase(rebinInfo_eTOF_Plus[rebinIndex].begin()); //partIndex
      rebinInfo_eTOF_Plus[rebinIndex].erase(rebinInfo_eTOF_Plus[rebinIndex].begin()); //detIndex
    }
  }

  for(int rebinIndex = rebinInfo_TPC_Minus.size() - 1; rebinIndex >= 0; rebinIndex--){
    int detIndex = (int) rebinInfo_TPC_Minus[rebinIndex][1];
    if(detIndex != 0){
      rebinInfo_TPC_Minus.erase(rebinInfo_TPC_Minus.begin() + rebinIndex);
    }else{
      rebinInfo_TPC_Minus[rebinIndex].erase(rebinInfo_TPC_Minus[rebinIndex].begin()); //partIndex
      rebinInfo_TPC_Minus[rebinIndex].erase(rebinInfo_TPC_Minus[rebinIndex].begin()); //detIndex
    }
  }

  for(int rebinIndex = rebinInfo_bTOF_Minus.size() - 1; rebinIndex >= 0; rebinIndex--){
    int detIndex = (int) rebinInfo_bTOF_Minus[rebinIndex][1];
    if(detIndex != 1){
      rebinInfo_bTOF_Minus.erase(rebinInfo_bTOF_Minus.begin() + rebinIndex);
    }else{
      rebinInfo_bTOF_Minus[rebinIndex].erase(rebinInfo_bTOF_Minus[rebinIndex].begin()); //partIndex
      rebinInfo_bTOF_Minus[rebinIndex].erase(rebinInfo_bTOF_Minus[rebinIndex].begin()); //detIndex
    }
  }

  for(int rebinIndex = rebinInfo_eTOF_Minus.size() - 1; rebinIndex >= 0; rebinIndex--){
    int detIndex = (int) rebinInfo_eTOF_Minus[rebinIndex][1];
    if(detIndex != 2){
      rebinInfo_eTOF_Minus.erase(rebinInfo_eTOF_Minus.begin() + rebinIndex);
    }else{
      rebinInfo_eTOF_Minus[rebinIndex].erase(rebinInfo_eTOF_Minus[rebinIndex].begin()); //partIndex
      rebinInfo_eTOF_Minus[rebinIndex].erase(rebinInfo_eTOF_Minus[rebinIndex].begin()); //detIndex
    }
  }





  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "Warning in ZFitter::loadDataHistograms.. couldn't find file " << a_filename.c_str() << endl;
    return;
  }


  TH1I* centEvents = (TH1I*) inFile->Get(Form("centEvents%s", m_partInfo->GetParticleName(a_partIndex,0).Data()));
  if(centEvents){ //if every particle was done at once in a single root file
    m_centEvents[a_partIndex] = centEvents;
  }else{
    cout << "May Be An Error: Can't Load " << Form("centEvents%s", m_partInfo->GetParticleName(a_partIndex,0).Data()) << endl;
  }

  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex ==3 || partIndex == 8){
      m_electron_rap_dEdx[partIndex][0] = NULL;
      m_electron_rap_dEdx[partIndex][1] = NULL;
      continue;
    }
    m_electron_rap_dEdx[partIndex][0] = (TH2I*) inFile->Get(Form("positrons_rap_dEdx_%s",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_electron_rap_dEdx[partIndex][1] = (TH2I*) inFile->Get(Form("electron_rap_dEdx_%s",m_partInfo->GetParticleName(partIndex,0).Data()));
  }


  inFile->cd();
  string firstDir = "yields";
  string secondDir = Form("%s_space",m_partInfo->GetParticleName(a_partIndex,0).Data());


  #ifdef _ZFITTER_DEBUG_
    cout << "Grabbing " <<  Form("/%s/%s/zTPC_AllCent", firstDir.c_str(),secondDir.c_str()) << "  from " << inFile->GetName() << endl;
  #endif
  TH3I* tempAllCentPlus   = (TH3I*) inFile->Get( Form("/%s/%s/zTPC_%s_AllCent", firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,1).Data()) );
  TH3I* tempAllCentMinus  =  (TH3I*) inFile->Get( Form("/%s/%s/zTPC_%s_AllCent", firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,-1).Data()) );
  #ifdef _ZFITTER_DEBUG_
    cout << "tempAllCentPlus: " << tempAllCentPlus << endl;
    cout << "tempAllCentMinus: " << tempAllCentMinus << endl;
  #endif
  TH3D* tempAllCentPlus_D  = HistogramUtilities::convertInt_to_DoubleHisto3D(tempAllCentPlus);
  TH3D* tempAllCentMinus_D = HistogramUtilities::convertInt_to_DoubleHisto3D(tempAllCentMinus);

  m_rap_mTm0_dEdx_AllCentRead_Plus[a_partIndex]  = HistogramUtilities::reBinVar3D(tempAllCentPlus_D, a_rebinX, a_rebinY, a_rebinZTPC);
  m_rap_mTm0_dEdx_AllCentRead_Minus[a_partIndex] = HistogramUtilities::reBinVar3D(tempAllCentMinus_D, a_rebinX, a_rebinY, a_rebinZTPC);
  if(rebinInfo_TPC_Plus.size() > 0){
    TH3D* tempPlus  = m_rap_mTm0_dEdx_AllCentRead_Plus[a_partIndex];
    m_rap_mTm0_dEdx_AllCentRead_Plus[a_partIndex] = HistogramUtilities::reBinBar3D_specificY(tempPlus, rebinInfo_TPC_Plus);
    #ifndef _SAVE_POINTERS_
      delete tempPlus;
    #endif
  }
  if(rebinInfo_TPC_Minus.size() > 0){
    TH3D* tempMinus = m_rap_mTm0_dEdx_AllCentRead_Minus[a_partIndex];
    m_rap_mTm0_dEdx_AllCentRead_Minus[a_partIndex] = HistogramUtilities::reBinBar3D_specificY(tempMinus, rebinInfo_TPC_Minus);
    #ifndef _SAVE_POINTERS_
      delete tempMinus;
    #endif
  }

  HistogramUtilities::makeSqrtNErrors3D(m_rap_mTm0_dEdx_AllCentRead_Plus[a_partIndex]);
  HistogramUtilities::makeSqrtNErrors3D(m_rap_mTm0_dEdx_AllCentRead_Minus[a_partIndex]);
  #ifndef _SAVE_POINTERS_
    delete tempAllCentPlus;
    delete tempAllCentMinus;
    delete tempAllCentPlus_D;
    delete tempAllCentMinus_D;
  #endif
  if(a_partIndex < 2) m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex] = (TH3D*) inFile->Get( Form("/%s/%s/zTPC_PionElectron_Isolated", firstDir.c_str(),secondDir.c_str()) );
  else                m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex] = NULL;
  if(a_partIndex < 2 && m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex]){
    m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex] = HistogramUtilities::reBinVar3D( m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex], a_rebinX, a_rebinY, a_rebinZTPC);
    if(rebinInfo_TPC_Plus.size() > 0){
      TH3D* tempPlus  = m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex];
      m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex] = HistogramUtilities::reBinBar3D_specificY(tempPlus, rebinInfo_TPC_Plus);
      #ifndef _SAVE_POINTERS_
        delete tempPlus;
      #endif
    }
    HistogramUtilities::makeSqrtNErrors3D(m_rap_mTm0_dEdx_allCent_PionElectron[a_partIndex]);
  }

  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    #ifdef _ZFITTER_DEBUG_
      cout << "  Reading Follow Histograms for Centrality " << centIndex << endl;
      cout << "Grabbing " <<Form("/%s/%s/zTPCPlus_%s_Cent%d",
        firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,0).Data(),centIndex) << endl;
    #endif
    vector< TH3I* > histosToDelete;
    vector< TH3D* > histosToDelete_double;

    TH3I* DHP = (TH3I*) inFile->Get( Form("/%s/%s/zTPCPlus_%s_Cent%d",
      firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,0).Data(),centIndex) );
    TH3I* DHM = (TH3I*) inFile->Get( Form("/%s/%s/zTPCMinus_%s_Cent%d",
      firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,0).Data(),centIndex) );
    TH3I* THP = (TH3I*) inFile->Get( Form("/%s/%s/zBTOFPlus_%s_Cent%d",
      firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,0).Data(),centIndex) );
    TH3I* THM = (TH3I*) inFile->Get( Form("/%s/%s/zBTOFMinus_%s_Cent%d",
      firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,0).Data(),centIndex) );
    TH3I* EHP = (TH3I*) inFile->Get( Form("/%s/%s/zETOFPlus_%s_Cent%d",
      firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,0).Data(),centIndex) );
    TH3I* EHM = (TH3I*) inFile->Get( Form("/%s/%s/zETOFMinus_%s_Cent%d",
      firstDir.c_str(),secondDir.c_str(),m_partInfo->GetParticleName(a_partIndex,0).Data(),centIndex) );

    histosToDelete.push_back(DHP);
    histosToDelete.push_back(DHM);
    histosToDelete.push_back(THP);
    histosToDelete.push_back(THM);
    histosToDelete.push_back(EHP);
    histosToDelete.push_back(EHM);

    #ifdef _ZFITTER_DEBUG_
      cout << "DHP: " << DHP << endl;
      cout << "DHM: " << DHM << endl;
      cout << "THP: " << THP << endl;
      cout << "THM: " << THM << endl;
    #endif
    TH3D* DHP_D = HistogramUtilities::convertInt_to_DoubleHisto3D(DHP);
    TH3D* DHM_D = HistogramUtilities::convertInt_to_DoubleHisto3D(DHM);
    TH3D* THP_D = HistogramUtilities::convertInt_to_DoubleHisto3D(THP);
    TH3D* THM_D = HistogramUtilities::convertInt_to_DoubleHisto3D(THM);
    TH3D* EHP_D = NULL;
    TH3D* EHM_D = NULL;
    if(m_hasETOF){
      #ifdef _ZFITTER_DEBUG_
        cout << "EHP: " << EHP << endl;
        cout << "EHM: " << EHM << endl;
      #endif
      EHP_D = HistogramUtilities::convertInt_to_DoubleHisto3D(EHP);
      EHM_D = HistogramUtilities::convertInt_to_DoubleHisto3D(EHM);
    }

    histosToDelete_double.push_back(DHP_D);
    histosToDelete_double.push_back(DHM_D);
    histosToDelete_double.push_back(THP_D);
    histosToDelete_double.push_back(THM_D);
    histosToDelete_double.push_back(EHP_D);
    histosToDelete_double.push_back(EHM_D);

    DHP_D = HistogramUtilities::reBinVar3D(DHP_D, a_rebinX, a_rebinY, a_rebinZTPC);
    DHM_D = HistogramUtilities::reBinVar3D(DHM_D, a_rebinX, a_rebinY, a_rebinZTPC);
    THP_D = HistogramUtilities::reBinVar3D(THP_D, a_rebinX, a_rebinY, a_rebinZBTOF);
    THM_D = HistogramUtilities::reBinVar3D(THM_D, a_rebinX, a_rebinY, a_rebinZBTOF);
    if(rebinInfo_TPC_Plus.size() > 0){
      histosToDelete_double.push_back(DHP_D);
      TH3D* originalHisto = DHP_D;
      DHP_D = HistogramUtilities::reBinBar3D_specificY(DHP_D, rebinInfo_TPC_Plus);
      #ifndef _SAVE_POINTERS_
        delete originalHisto;
      #endif
    }
    if(rebinInfo_TPC_Minus.size() > 0){
      histosToDelete_double.push_back(DHM_D);
      TH3D* originalHisto = DHM_D;
      DHM_D = HistogramUtilities::reBinBar3D_specificY(DHM_D, rebinInfo_TPC_Minus);
      #ifndef _SAVE_POINTERS_
        delete originalHisto;
      #endif
    }
    if(rebinInfo_bTOF_Plus.size() > 0){
      histosToDelete_double.push_back(THP_D);
      TH3D* originalHisto = THP_D;
      THP_D = HistogramUtilities::reBinBar3D_specificY(THP_D, rebinInfo_bTOF_Plus);
      #ifndef _SAVE_POINTERS_
        delete originalHisto;
      #endif
    }
    if(rebinInfo_bTOF_Minus.size() > 0){
      histosToDelete_double.push_back(THM_D);
      TH3D* originalHisto = THM_D;
      THM_D = HistogramUtilities::reBinBar3D_specificY(THM_D, rebinInfo_bTOF_Minus);
      #ifndef _SAVE_POINTERS_
        delete originalHisto;
      #endif
    }
    HistogramUtilities::makeSqrtNErrors3D(DHP_D);
    HistogramUtilities::makeSqrtNErrors3D(DHM_D);
    HistogramUtilities::makeSqrtNErrors3D(THP_D);
    HistogramUtilities::makeSqrtNErrors3D(THM_D);
    if(m_hasETOF){
      EHP_D = HistogramUtilities::reBinVar3D(EHP_D, a_rebinX, a_rebinY, a_rebinZETOF);
      EHM_D = HistogramUtilities::reBinVar3D(EHM_D, a_rebinX, a_rebinY, a_rebinZETOF);
      if(rebinInfo_eTOF_Plus.size() > 0){
        histosToDelete_double.push_back(EHP_D);
        TH3D* originalHisto = EHP_D;
        EHP_D = HistogramUtilities::reBinBar3D_specificY(EHP_D, rebinInfo_eTOF_Plus);
        #ifndef _SAVE_POINTERS_
          delete originalHisto;
        #endif
      }
      if(rebinInfo_eTOF_Minus.size() > 0){
        histosToDelete_double.push_back(EHM_D);
        TH3D* originalHisto = EHM_D;
        EHM_D = HistogramUtilities::reBinBar3D_specificY(EHM_D, rebinInfo_eTOF_Minus);
        #ifndef _SAVE_POINTERS_
          delete originalHisto;
        #endif
      }
      HistogramUtilities::makeSqrtNErrors3D(EHP_D);
      HistogramUtilities::makeSqrtNErrors3D(EHM_D);
    }

    #ifndef _SAVE_POINTERS_
      for(int iii = 0; iii < histosToDelete.size(); iii++){
        if(histosToDelete[iii]){delete histosToDelete[iii];}
        if(histosToDelete_double[iii]){delete histosToDelete_double[iii];}
      }
    #endif

    m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex]          = DHP_D;
    m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex]         = DHM_D;
    m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex]   = THP_D;
    m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex]  = THM_D;
    if(m_hasETOF){
      m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex]  = EHP_D;
      m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex] = EHM_D;
      #ifdef _ZFITTER_DEBUG_
        cout << "ETOF HAS BEEN READ?  ptr: " << EHP << endl;
      #endif
    }else{
      m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex]  = NULL;
      m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex] = NULL;
    }

  }

  for(int partIndex = 0; partIndex < 9; partIndex++){
    TH3I* isoTemp = (TH3I*) inFile->Get( Form("/%s/%s/zTPC_%s_Isolated",
     firstDir.c_str(),secondDir.c_str(), m_partInfo->GetParticleName(partIndex,1).Data()) );
    #ifdef _ZFITTER_DEBUG_
      cout << "Getting Isolated Particle Histograms for partIndex: " << partIndex << endl;
      cout << "isoTemp plus: " << isoTemp << endl;
    #endif
    TH3D* isoTemp_D = HistogramUtilities::convertInt_to_DoubleHisto3D(isoTemp);
    TH3D* isoTemp_D_Rebin = NULL;
    if(isoTemp_D){
      isoTemp_D_Rebin = HistogramUtilities::reBinVar3D(isoTemp_D, a_rebinX, a_rebinY, a_rebinZTPC);
      if(rebinInfo_TPC_Plus.size() > 0){
        TH3D* tempPlus  = isoTemp_D_Rebin;
        isoTemp_D_Rebin = HistogramUtilities::reBinBar3D_specificY(tempPlus, rebinInfo_TPC_Plus);
        #ifndef _SAVE_POINTERS_
          delete tempPlus;
        #endif
      }
      HistogramUtilities::makeSqrtNErrors3D(isoTemp_D_Rebin);
    }
    m_rap_mTm0_dEdx_Plus_allCentIsolated[a_partIndex][partIndex] = isoTemp_D_Rebin;
    #ifndef _SAVE_POINTERS_
      delete isoTemp;
      delete isoTemp_D;
    #endif


    isoTemp = (TH3I*) inFile->Get( Form("/%s/%s/zTPC_%s_Isolated",
     firstDir.c_str(),secondDir.c_str(), m_partInfo->GetParticleName(partIndex,-1).Data()) );
    //cout << "isoTemp minus: " << isoTemp << endl;
    isoTemp_D = HistogramUtilities::convertInt_to_DoubleHisto3D(isoTemp);
    isoTemp_D_Rebin = NULL;
    if(isoTemp_D){
      isoTemp_D_Rebin = HistogramUtilities::reBinVar3D(isoTemp_D, a_rebinX, a_rebinY, a_rebinZTPC);
      if(rebinInfo_TPC_Minus.size() > 0){
        TH3D* tempMinus  = isoTemp_D_Rebin;
        isoTemp_D_Rebin = HistogramUtilities::reBinBar3D_specificY(tempMinus, rebinInfo_TPC_Minus);
        #ifndef _SAVE_POINTERS_
          delete tempMinus;
        #endif
      }
      HistogramUtilities::makeSqrtNErrors3D(isoTemp_D_Rebin);
    }
    m_rap_mTm0_dEdx_Minus_allCentIsolated[a_partIndex][partIndex] = isoTemp_D_Rebin;
    #ifndef _SAVE_POINTERS_
      delete isoTemp;
      delete isoTemp_D;
    #endif
  }

  if(a_convert16to9CentBins){
    cout << "Converting 16 Cent Bins to 9" << endl;
    if(m_numCentralities != 16){
      cout << "you cant use the a_convert16to9CentBins unless you actually have 16 bins..." << endl;
      return;
    }

    for(int centIndex = 2; centIndex < m_numCentralities; centIndex+=2){
      int newCentIndex = 1+centIndex/2;
      m_rap_mTm0_dEdx_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex],m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex+1]);
      m_rap_mTm0_dEdx_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex],m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex+1]);
      m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex],m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex+1]);
      m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex],m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex+1]);
      if(m_hasETOF){
        m_rap_mTm0_invBetaETOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex],m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex+1]);
        m_rap_mTm0_invBetaETOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex],m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex+1]);
      }
    }
    #ifndef _SAVE_POINTERS_
      for(int centIndex = 9; centIndex < m_numCentralities; centIndex+=1){
        delete m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex];
        delete m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex];
        delete m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex];
        delete m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex];
        if(m_hasETOF){
          delete m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex];
          delete m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex];
        }
      }
    #endif

    m_numCentralities = 9;

  }

  if(a_combine5and10CentBins){
    cout << "Combining 0-5%% and 5-10%% Centrialities" << endl;

    for(int centIndex = 1; centIndex < m_numCentralities; centIndex++){
      int newCentIndex = centIndex - 1;
      if(centIndex == 1){
        cout << "adding centIndex " << centIndex-1 << " " << m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex-1]->GetTitle() << " and centIndex " << centIndex << " " << m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex]->GetTitle() << " into centIndex " << newCentIndex <<  endl;
        m_rap_mTm0_dEdx_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex-1],m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex]);
        cout << "adding centIndex " << centIndex-1 << " " << m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex-1]->GetTitle() << " and centIndex " << centIndex << " " << m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex]->GetTitle() << " into centIndex " << newCentIndex <<  endl;
        m_rap_mTm0_dEdx_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex-1],m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex]);
        cout << "adding centIndex " << centIndex-1 << " " << m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex-1]->GetTitle() << " and centIndex " << centIndex << " " << m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex]->GetTitle() << " into centIndex " << newCentIndex <<  endl;
        m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex-1],m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex]);
        cout << "adding centIndex " << centIndex-1 << " " << m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex-1]->GetTitle() << " and centIndex " << centIndex << " " << m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex]->GetTitle() << " into centIndex " << newCentIndex <<  endl;
        m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex-1],m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex]);
        if(m_hasETOF){
          cout << "adding centIndex " << centIndex-1 << " " << m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex-1]->GetTitle() << " and centIndex " << centIndex << " " << m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex]->GetTitle() << " into centIndex " << newCentIndex <<  endl;
          m_rap_mTm0_invBetaETOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex-1],m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex]);
          cout << "adding centIndex " << centIndex-1 << " " << m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex-1]->GetTitle() << " and centIndex " << centIndex << " " << m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex]->GetTitle() << " into centIndex " << newCentIndex <<  endl;
          m_rap_mTm0_invBetaETOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex-1],m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex]);
        }
      }else{
        cout << "pushing " << m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex]->GetTitle() << "from centIndex " << centIndex << " to " << newCentIndex << endl;
        m_rap_mTm0_dEdx_Plus[a_partIndex][newCentIndex] = m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex];
        cout << "pushing " << m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex]->GetTitle() << "from centIndex " << centIndex << " to " << newCentIndex << endl;
        m_rap_mTm0_dEdx_Minus[a_partIndex][newCentIndex] = m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex];
        cout << "pushing " << m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex]->GetTitle() << "from centIndex " << centIndex << " to " << newCentIndex << endl;
        m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][newCentIndex] = m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex];
        cout << "pushing " << m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex]->GetTitle() << "from centIndex " << centIndex << " to " << newCentIndex << endl;
        m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][newCentIndex] = m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex];
        if(m_hasETOF){
          cout << "pushing " << m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex]->GetTitle() << "from centIndex " << centIndex << " to " << newCentIndex << endl;
          m_rap_mTm0_invBetaETOF_Plus[a_partIndex][newCentIndex] = m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex];
          cout << "pushing " << m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex]->GetTitle() << "from centIndex " << centIndex << " to " << newCentIndex << endl;
          m_rap_mTm0_invBetaETOF_Minus[a_partIndex][newCentIndex] = m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex];
        }
      }
    }

    /*for(int centIndex = m_numCentralities-1; centIndex < m_numCentralities; centIndex+=1){
      delete m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex];
      delete m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex];
      delete m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex];
      delete m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex];
      if(m_hasETOF){
        delete m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex];
        delete m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex];
      }
    }*/

    m_numCentralities = m_numCentralities - 1;
  }

  if(a_combine40and50_60and70_centBins){
    cout << "Converting 40-50%% & 50-60%% and 60-70%% & 70-80%% Cent Bins" << endl;
    cout << "Initial number of centrality bins: " << m_numCentralities << endl;
    if(m_numCentralities != 9 && (m_numCentralities != 8 && !a_combine5and10CentBins)){
      cout << "you cant use the a_combine40and50_60and70_centBins unless you have 8 or 9 bins..." << endl;
      return;
    }

    for(int centIndex = m_numCentralities - 4; centIndex < m_numCentralities; centIndex+=2){
      int newCentIndex = (m_numCentralities-4)+(centIndex-(m_numCentralities-4))/2;
      cout << "for original centIndex " << centIndex << " mapping to new centIndex " << newCentIndex << endl;
      cout << m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex]->GetTitle() << " and " << m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex+1]->GetTitle() << endl;
      m_rap_mTm0_dEdx_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex],m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex+1]);
      cout << m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex]->GetTitle() << " and " << m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex+1]->GetTitle() << endl;
      m_rap_mTm0_dEdx_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex],m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex+1]);
      cout << m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex]->GetTitle() << " and " << m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex+1]->GetTitle() << endl;
      m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex],m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex+1]);
      cout << m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex]->GetTitle() << " and " << m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex+1]->GetTitle() << endl;
      m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex],m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex+1]);
      if(m_hasETOF){
        cout << m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex]->GetTitle() << " and " << m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex+1]->GetTitle() << endl;
        m_rap_mTm0_invBetaETOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex],m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex+1]);
        cout << m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex]->GetTitle() << " and " << m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex+1]->GetTitle() << endl;
        m_rap_mTm0_invBetaETOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex],m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex+1]);
      }
    }

    cout << "finished adding" << endl;

    for(int centIndex = m_numCentralities-2; centIndex < m_numCentralities; centIndex+=1){
      cout << "removing old centrality arrays for old centIndex: " << centIndex << endl;
      #ifndef _SAVE_POINTERS_
        delete m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex];
        delete m_rap_mTm0_dEdx_Minus[a_partIndex][centIndex];
        delete m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][centIndex];
        delete m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][centIndex];
        if(m_hasETOF){
          delete m_rap_mTm0_invBetaETOF_Plus[a_partIndex][centIndex];
          delete m_rap_mTm0_invBetaETOF_Minus[a_partIndex][centIndex];
        }
      #endif
    }

    m_numCentralities = m_numCentralities - 2;
    cout << "Final number of centrality bins: " << m_numCentralities << endl;

  }




  cout << "End of loadDataHistograms(...)" << endl;

}



//=========================================================================================================================================================================


void ZFitter::loadCalibrationHistogram(string a_filename, int a_partIndex){

  cout << "Start of loadCalibrationHistogram(" << a_filename << "," << a_partIndex << ")" << endl;

  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "Warning in ZFitter::loadCalibrationHistogram.. couldn't find file " << a_filename.c_str() << endl;
    return;
  }

  TH1I* centEvents = (TH1I*) inFile->Get(Form("centEvents%s", m_partInfo->GetParticleName(a_partIndex,0).Data()));
  if(centEvents){ //if every particle was done at once in a single root file
    m_centEvents[a_partIndex] = centEvents;
  }else{
    cout << "May Be An Error: Can't Load " << Form("centEvents%s", m_partInfo->GetParticleName(a_partIndex,0).Data()) << endl;
  }


}





//=========================================================================================================================================================================

void ZFitter::makeAllCentHistos(){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to makeAllCentHistos()" << endl;
    cout << "m_numParticles: " << m_numParticles << endl;
  #endif
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(!m_particleLoaded[partIndex]) continue;
    #ifdef _ZFITTER_DEBUG_
      cout << "     Making All Cent for Particle: " << m_partInfo->GetParticleName(partIndex,0) << endl;
    #endif
    m_rap_mTm0_dEdx_Plus_allCent[partIndex]        = (TH3D*) m_rap_mTm0_dEdx_Plus[partIndex][0]->Clone();
    m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex] = (TH3D*) m_rap_mTm0_invBetaBTOF_Plus[partIndex][0]->Clone();
    if(m_hasETOF) m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex] = (TH3D*) m_rap_mTm0_invBetaETOF_Plus[partIndex][0]->Clone();
    m_rap_mTm0_dEdx_Minus_allCent[partIndex]        = (TH3D*) m_rap_mTm0_dEdx_Minus[partIndex][0]->Clone();
    m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex] = (TH3D*) m_rap_mTm0_invBetaBTOF_Minus[partIndex][0]->Clone();
    if(m_hasETOF) m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex] = (TH3D*) m_rap_mTm0_invBetaETOF_Minus[partIndex][0]->Clone();
    for(int centIndex = 1; centIndex < m_numCentralities; centIndex++){
      m_rap_mTm0_dEdx_Plus_allCent[partIndex]->Add(m_rap_mTm0_dEdx_Plus[partIndex][centIndex],1.0);
      m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex]->Add(m_rap_mTm0_invBetaBTOF_Plus[partIndex][centIndex],1.0);
      if(m_hasETOF) m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex]->Add(m_rap_mTm0_invBetaETOF_Plus[partIndex][centIndex],1.0);
      m_rap_mTm0_dEdx_Minus_allCent[partIndex]->Add(m_rap_mTm0_dEdx_Minus[partIndex][centIndex],1.0);
      m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex]->Add(m_rap_mTm0_invBetaBTOF_Minus[partIndex][centIndex],1.0);
      if(m_hasETOF) m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex]->Add(m_rap_mTm0_invBetaETOF_Minus[partIndex][centIndex],1.0);
    }
    m_rap_mTm0_dEdx_allCent[partIndex] = (TH3D*) m_rap_mTm0_dEdx_Plus_allCent[partIndex]->Clone();
    m_rap_mTm0_dEdx_allCent[partIndex]->Add(m_rap_mTm0_dEdx_Minus_allCent[partIndex]);
    m_rap_mTm0_invBetaBTOF_allCent[partIndex] = (TH3D*) m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex]->Clone();
    m_rap_mTm0_invBetaBTOF_allCent[partIndex]->Add(m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex]);
    if(m_hasETOF){
      m_rap_mTm0_invBetaETOF_allCent[partIndex] = (TH3D*) m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex]->Clone();
      m_rap_mTm0_invBetaETOF_allCent[partIndex]->Add(m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex]);
    }
  }
  return;
}


//=========================================================================================================================================================================

void ZFitter::createQuickStartFile(string a_filename){

  cout << "Start of createQuickStartFile(" << a_filename << " ) " << endl;
  TFile* outFile = new TFile(a_filename.c_str(),"RECREATE");
  if(!outFile){
    cout << "Warning in ZFitter::createQuickStartFile.. couldn't make file " << a_filename.c_str() << endl;
    return;
  }

  for(unsigned int partIndex = 0; partIndex < 9; partIndex++){

    if(!m_particleLoaded[partIndex]) continue;
    cout << "partIndex: " << partIndex << endl;
    if(m_centEvents[partIndex]){
      m_centEvents[partIndex]->Write();
    }

    if(m_electron_rap_dEdx[partIndex][0]){
      m_electron_rap_dEdx[partIndex][0]->SetName(Form("ELEC_%d",partIndex));
      m_electron_rap_dEdx[partIndex][0]->Write();
    };
    if(m_electron_rap_dEdx[partIndex][1]){
      m_electron_rap_dEdx[partIndex][1]->SetName(Form("POSIT_%d",partIndex));
      m_electron_rap_dEdx[partIndex][1]->Write();
    };

    if(partIndex < 2 && m_rap_mTm0_dEdx_allCent_PionElectron[partIndex]){
      m_rap_mTm0_dEdx_allCent_PionElectron[partIndex]->SetName(Form("pionElectronAllCent_%02d",partIndex));
      m_rap_mTm0_dEdx_allCent_PionElectron[partIndex]->Write();
    }

    m_rap_mTm0_dEdx_AllCentRead_Plus[partIndex]->SetName(Form("ZTPC_AllCent_Read_%02d_Plus",partIndex));
    m_rap_mTm0_dEdx_AllCentRead_Minus[partIndex]->SetName(Form("ZTPC_AllCent_Read_%02d_Minus",partIndex));
    m_rap_mTm0_dEdx_AllCentRead_Plus[partIndex]->Write();
    m_rap_mTm0_dEdx_AllCentRead_Minus[partIndex]->Write();


    if(m_rap_mTm0_dEdx_allCent[partIndex]){
      m_rap_mTm0_dEdx_allCent[partIndex]->SetName(Form("ZTPC_AllCent_%02d",partIndex));
      m_rap_mTm0_dEdx_allCent[partIndex]->Write();
    }
    if(m_rap_mTm0_dEdx_Plus_allCent[partIndex]){
      m_rap_mTm0_dEdx_Plus_allCent[partIndex]->SetName(Form("ZTPC_AllCent_Plus_%02d",partIndex));
      m_rap_mTm0_dEdx_Plus_allCent[partIndex]->Write();
    }
    if(m_rap_mTm0_dEdx_Minus_allCent[partIndex]){
      m_rap_mTm0_dEdx_Minus_allCent[partIndex]->SetName(Form("ZTPC_AllCent_Minus_%02d",partIndex));
      m_rap_mTm0_dEdx_Minus_allCent[partIndex]->Write();
    }
    if(m_rap_mTm0_invBetaBTOF_allCent[partIndex]){
      m_rap_mTm0_invBetaBTOF_allCent[partIndex]->SetName(Form("ZBTOF_AllCent_%02d",partIndex));
      m_rap_mTm0_invBetaBTOF_allCent[partIndex]->Write();
    }
    if(m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex]){
      m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex]->SetName(Form("ZBTOF_AllCent_Plus_%02d",partIndex));
      m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex]->Write();
    }
    if(m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex]){
      m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex]->SetName(Form("ZBTOF_AllCent_Minus_%02d",partIndex));
      m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex]->Write();
    }
    if(m_hasETOF){
      if(m_rap_mTm0_invBetaETOF_allCent[partIndex]){
        m_rap_mTm0_invBetaETOF_allCent[partIndex]->SetName(Form("ZETOF_AllCent_%02d",partIndex));
        m_rap_mTm0_invBetaETOF_allCent[partIndex]->Write();
      }
      if(m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex]){
        m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex]->SetName(Form("ZETOF_AllCent_Plus_%02d",partIndex));
        m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex]->Write();
      }
      if(m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex]){
        m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex]->SetName(Form("ZETOF_AllCent_Minus_%02d",partIndex));
        m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex]->Write();
      }
    }



    for(unsigned int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      cout << "cent: " << centIndex << endl;
      if(m_rap_mTm0_dEdx_Plus[partIndex][centIndex]){
        m_rap_mTm0_dEdx_Plus[partIndex][centIndex]->SetName(Form("ZTPC_Plus_%02d_%02d",partIndex,centIndex));
        m_rap_mTm0_dEdx_Plus[partIndex][centIndex]->Write();
      }
      if(m_rap_mTm0_dEdx_Minus[partIndex][centIndex]){
        m_rap_mTm0_dEdx_Minus[partIndex][centIndex]->SetName(Form("ZTPC_Minus_%02d_%02d",partIndex,centIndex));
        m_rap_mTm0_dEdx_Minus[partIndex][centIndex]->Write();
      }
      if(m_rap_mTm0_invBetaBTOF_Plus[partIndex][centIndex]){
        m_rap_mTm0_invBetaBTOF_Plus[partIndex][centIndex]->SetName(Form("ZBTOF_Plus_%02d_%02d",partIndex,centIndex));
        m_rap_mTm0_invBetaBTOF_Plus[partIndex][centIndex]->Write();
      }
      if(m_rap_mTm0_invBetaBTOF_Minus[partIndex][centIndex]){
        m_rap_mTm0_invBetaBTOF_Minus[partIndex][centIndex]->SetName(Form("ZBTOF_Minus_%02d_%02d",partIndex,centIndex));
        m_rap_mTm0_invBetaBTOF_Minus[partIndex][centIndex]->Write();
      }
      if(m_hasETOF){
        if(m_rap_mTm0_invBetaETOF_Plus[partIndex][centIndex]){
          m_rap_mTm0_invBetaETOF_Plus[partIndex][centIndex]->SetName(Form("ZETOF_Plus_%02d_%02d",partIndex,centIndex));
          m_rap_mTm0_invBetaETOF_Plus[partIndex][centIndex]->Write();
        }
        if(m_rap_mTm0_invBetaETOF_Minus[partIndex][centIndex]){
          m_rap_mTm0_invBetaETOF_Minus[partIndex][centIndex]->SetName(Form("ZETOF_Minus_%02d_%02d",partIndex,centIndex));
          m_rap_mTm0_invBetaETOF_Minus[partIndex][centIndex]->Write();
        }
      }
    }// cent index

    for(unsigned int partIndex2 = 0; partIndex2 < 9; partIndex2++){
      cout << "partIndex2: " << partIndex2 << endl;
      if(!m_rap_mTm0_dEdx_Plus_allCentIsolated[partIndex][partIndex2]) continue;
      m_rap_mTm0_dEdx_Plus_allCentIsolated[partIndex][partIndex2]->SetName(Form("ZTPC_Iso_Plus%02d_%02d",partIndex,partIndex2));
      cout << "set plus name" << endl;
      m_rap_mTm0_dEdx_Minus_allCentIsolated[partIndex][partIndex2]->SetName(Form("ZTPC_Iso_Minus%02d_%02d",partIndex,partIndex2));
      cout << "set minus name" << endl;
      m_rap_mTm0_dEdx_Plus_allCentIsolated[partIndex][partIndex2]->Write();
      cout << "wrote plus" << endl;
      m_rap_mTm0_dEdx_Minus_allCentIsolated[partIndex][partIndex2]->Write();
      cout << "end of partIndex2: " << partIndex2 << " loop" << endl;
    }// end partIndex2
      cout << "end of partIndex: " << partIndex << " loop" << endl;

  }//end partIndex loop


  cout << "End of createQuickStartFile(...)" << endl;






}


//=========================================================================================================================================================================


void ZFitter::loadQuickStartFile(string a_filename, int a_numCent, int a_partIndex, int a_detectorIndex){
  cout << "Start of loadQuickStartFile(" << a_filename << " ) " << endl;

  m_numCentralities = a_numCent;
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "Warning in ZFitter::loadQuickStartFile.. couldn't find file " << a_filename.c_str() << endl;
    return;
  }

  for(unsigned int partIndex = 0; partIndex < 9; partIndex++){
    if(a_partIndex >= 0 && a_partIndex != partIndex) continue;
    m_centEvents[partIndex] = (TH1I*) inFile->Get(Form("centEvents%s", m_partInfo->GetParticleName(partIndex,0).Data()));

    if(a_detectorIndex < 0 || a_detectorIndex == 0){
      if(partIndex < 2){
        m_rap_mTm0_dEdx_allCent_PionElectron[partIndex] = (TH3D*) inFile->Get(Form("pionElectronAllCent_%02d",partIndex));
      }else{
        m_rap_mTm0_dEdx_allCent_PionElectron[partIndex] = NULL;
      }
      m_electron_rap_dEdx[partIndex][0]            = (TH2I*) inFile->Get(Form("ELEC_%d",partIndex));
      m_electron_rap_dEdx[partIndex][1]            = (TH2I*) inFile->Get(Form("POSIT_%d",partIndex));
      m_rap_mTm0_dEdx_AllCentRead_Plus[partIndex]  = (TH3D*) inFile->Get(Form("ZTPC_AllCent_Read_%02d_Plus",partIndex));
      m_rap_mTm0_dEdx_AllCentRead_Minus[partIndex] = (TH3D*) inFile->Get(Form("ZTPC_AllCent_Read_%02d_Minus",partIndex));
      m_rap_mTm0_dEdx_allCent[partIndex]           = (TH3D*) inFile->Get(Form("ZTPC_AllCent_%02d",partIndex));
      m_rap_mTm0_dEdx_Plus_allCent[partIndex]      = (TH3D*) inFile->Get(Form("ZTPC_AllCent_Plus_%02d",partIndex));
      m_rap_mTm0_dEdx_Minus_allCent[partIndex]     = (TH3D*) inFile->Get(Form("ZTPC_AllCent_Minus_%02d",partIndex));
      for(unsigned int partIndex2 = 0; partIndex2 < 9; partIndex2++){
        m_rap_mTm0_dEdx_Plus_allCentIsolated[partIndex][partIndex2]  = (TH3D*) inFile->Get(Form("ZTPC_Iso_Plus%02d_%02d",partIndex,partIndex2));
        m_rap_mTm0_dEdx_Minus_allCentIsolated[partIndex][partIndex2] = (TH3D*) inFile->Get(Form("ZTPC_Iso_Minus%02d_%02d",partIndex,partIndex2));
      }// end partIndex2
      if(m_rap_mTm0_dEdx_allCent[partIndex]){
        m_particleLoaded[partIndex] = true;
      }
    }

    if(a_detectorIndex < 0 || a_detectorIndex == 1){ // LOAD BTOF
      m_rap_mTm0_invBetaBTOF_allCent[partIndex]       = (TH3D*) inFile->Get(Form("ZBTOF_AllCent_%02d",partIndex));
      m_rap_mTm0_invBetaBTOF_Plus_allCent[partIndex]  = (TH3D*) inFile->Get(Form("ZBTOF_AllCent_Plus_%02d",partIndex));
      m_rap_mTm0_invBetaBTOF_Minus_allCent[partIndex] = (TH3D*) inFile->Get(Form("ZBTOF_AllCent_Minus_%02d",partIndex));
      if(m_rap_mTm0_invBetaBTOF_allCent[partIndex]){
        m_particleLoaded[partIndex] = true;
      }
    }
    if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
      m_rap_mTm0_invBetaETOF_allCent[partIndex]       = (TH3D*) inFile->Get(Form("ZETOF_AllCent_%02d",partIndex));
      m_rap_mTm0_invBetaETOF_Plus_allCent[partIndex]  = (TH3D*) inFile->Get(Form("ZETOF_AllCent_Plus_%02d",partIndex));
      m_rap_mTm0_invBetaETOF_Minus_allCent[partIndex] = (TH3D*) inFile->Get(Form("ZETOF_AllCent_Minus_%02d",partIndex));
      if(m_rap_mTm0_invBetaETOF_allCent[partIndex]){
        m_particleLoaded[partIndex] = true;
      }
    }

    for(unsigned int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      if(a_detectorIndex < 0 || a_detectorIndex == 0){
        m_rap_mTm0_dEdx_Plus[partIndex][centIndex]         = (TH3D*) inFile->Get(Form("ZTPC_Plus_%02d_%02d",partIndex,centIndex));
        m_rap_mTm0_dEdx_Minus[partIndex][centIndex]        = (TH3D*) inFile->Get(Form("ZTPC_Minus_%02d_%02d",partIndex,centIndex));
      }
      if(a_detectorIndex < 0 || a_detectorIndex == 1){
        m_rap_mTm0_invBetaBTOF_Plus[partIndex][centIndex]  = (TH3D*) inFile->Get(Form("ZBTOF_Plus_%02d_%02d",partIndex,centIndex));
        m_rap_mTm0_invBetaBTOF_Minus[partIndex][centIndex] = (TH3D*) inFile->Get(Form("ZBTOF_Minus_%02d_%02d",partIndex,centIndex));
      }
      if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
        m_rap_mTm0_invBetaETOF_Plus[partIndex][centIndex]  = (TH3D*) inFile->Get(Form("ZETOF_Plus_%02d_%02d",partIndex,centIndex));
        m_rap_mTm0_invBetaETOF_Minus[partIndex][centIndex] = (TH3D*) inFile->Get(Form("ZETOF_Minus_%02d_%02d",partIndex,centIndex));
      }


    }// cent index

  }//end partIndex loop




  cout << "End of loadQuickStartFile(...)" << endl;




}


void ZFitter::loadTextFileInformation(string a_textFileList){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::loadTextFileInformation(" << a_textFileList << ")" << endl;
  #endif

  ifstream filelist(a_textFileList.c_str(),ios::in);
  string filename = "";
  while(filelist.good()){
    filelist >> filename;
    #ifdef _ZFITTER_DEBUG_
      cout << "Reading File " << filename << endl;
    #endif
    ifstream file(filename.c_str(),ios::in);
    int partIndex,plusMinusIndex,detIndex,centIndex,rapBin,mTm0Bin;
    double value,error;
    while(file.good()){
      file >> partIndex;
      file >> plusMinusIndex;
      file >> detIndex;
      file >> centIndex;
      file >> rapBin;
      file >> mTm0Bin;
      file >> value;
      file >> error;

      #ifdef _ZFITTER_DEBUG_
        cout << " Read in : " << partIndex << " " <<  plusMinusIndex
             << " " <<  detIndex << " " <<  centIndex << " " <<  rapBin << " " <<  mTm0Bin << " " <<  value
             << " " <<  error << endl;
      #endif

      if(detIndex == 0){
        m_Spectra_Cent_ZTPC[partIndex][centIndex][plusMinusIndex]->SetBinContent(rapBin,mTm0Bin,value);
        m_Spectra_Cent_ZTPC[partIndex][centIndex][plusMinusIndex]->SetBinError(rapBin,mTm0Bin,error);
      }else if(detIndex == 1){
        m_Spectra_Cent_ZbTOF[partIndex][centIndex][plusMinusIndex]->SetBinContent(rapBin,mTm0Bin,value);
        m_Spectra_Cent_ZbTOF[partIndex][centIndex][plusMinusIndex]->SetBinError(rapBin,mTm0Bin,error);
      }else if(detIndex == 2){
        m_Spectra_Cent_ZeTOF[partIndex][centIndex][plusMinusIndex]->SetBinContent(rapBin,mTm0Bin,value);
        m_Spectra_Cent_ZeTOF[partIndex][centIndex][plusMinusIndex]->SetBinError(rapBin,mTm0Bin,error);
      }
    }







    file.close();
  }



  filelist.close();
}


//=========================================================================================================================================================================

TH1D* ZFitter::grabSpecificHistogram(int a_partIndex, int a_charge, int a_centIndex, int a_detectorIndex, int a_rapBin, int a_mTm0Bin, int a_rebinCDFNumBins, bool a_useProminenceAllCent, double a_lowBound, double a_highBound, bool a_divideByNumEvents, bool a_propagateNumEvtError){

  TH3D* histo3DPM = NULL;
  if(a_charge < 0){
    if(a_centIndex < 0){
      if(a_detectorIndex == 0){
        if(a_useProminenceAllCent){
          histo3DPM = m_rap_mTm0_dEdx_AllCentRead_Minus[a_partIndex];
        }else{
          histo3DPM = m_rap_mTm0_dEdx_Minus_allCent[a_partIndex];
        }
      }else if(a_detectorIndex == 1) histo3DPM = m_rap_mTm0_invBetaBTOF_Minus_allCent[a_partIndex];
      else if(a_detectorIndex == 2)  histo3DPM = m_rap_mTm0_invBetaETOF_Minus_allCent[a_partIndex];
    }else{
      // specific centrality
      if(a_detectorIndex == 0)      histo3DPM = m_rap_mTm0_dEdx_Minus[a_partIndex][a_centIndex];
      else if(a_detectorIndex == 1) histo3DPM = m_rap_mTm0_invBetaBTOF_Minus[a_partIndex][a_centIndex];
      else if(a_detectorIndex == 2) histo3DPM = m_rap_mTm0_invBetaETOF_Minus[a_partIndex][a_centIndex];
    }
  }else{
    if(a_centIndex < 0){
      if(a_detectorIndex == 0){
        if(a_useProminenceAllCent){
          histo3DPM = m_rap_mTm0_dEdx_AllCentRead_Plus[a_partIndex];
        }else{
          histo3DPM = m_rap_mTm0_dEdx_Plus_allCent[a_partIndex];
        }
      }else if(a_detectorIndex == 1) histo3DPM = m_rap_mTm0_invBetaBTOF_Plus_allCent[a_partIndex];
      else if(a_detectorIndex == 2)  histo3DPM = m_rap_mTm0_invBetaETOF_Plus_allCent[a_partIndex];
    }else{
      // specific centrality
      if(a_detectorIndex == 0)      histo3DPM = m_rap_mTm0_dEdx_Plus[a_partIndex][a_centIndex];
      else if(a_detectorIndex == 1) histo3DPM = m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][a_centIndex];
      else if(a_detectorIndex == 2) histo3DPM = m_rap_mTm0_invBetaETOF_Plus[a_partIndex][a_centIndex];
    }
  }

  if(!histo3DPM){
    cout << " The histogram doesn't exist" << endl;
    return NULL;
  }

  string detNames[3] = {"TPC","BTOF","ETOF"};

  TH1D* histo = histo3DPM->ProjectionZ(Form("%s_%s_Cent_%02d_RapBin_%02d_mTm0Bin_%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), detNames[a_detectorIndex].c_str(), a_centIndex, a_rapBin , a_mTm0Bin),
         a_rapBin, a_rapBin, a_mTm0Bin, a_mTm0Bin);
  HistogramUtilities::makeSqrtNErrors(histo);


  if(a_lowBound > -998.0){
    TH1D* origHistToFit = histo;
    histo = HistogramUtilities::truncateHistogram(histo,a_lowBound,a_highBound);
    #ifndef _SAVE_POINTERS_
      delete origHistToFit;
    #endif
  }

  if(a_rebinCDFNumBins > 0 && !HistogramUtilities::isEmpty(histo)){
    TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(histo, a_rebinCDFNumBins, false, true, true); // isDifferrential doLogCDF
    #ifndef _SAVE_POINTERS_
      delete histo;
    #endif
    histo = rebinHisto;
  }



  if(a_centIndex >= 0 && a_divideByNumEvents && m_centEvents[a_partIndex]){
    double N_Evt     = m_centEvents[a_partIndex]->GetBinContent(a_centIndex+1);
    double N_Evt_Err = sqrt(N_Evt);
    for(int binX = 1; binX <= histo->GetNbinsX(); binX++){
      double content   = histo->GetBinContent(binX);
      double error     = histo->GetBinError(binX);
      double newContent= content/N_Evt;
      double newErr    = newContent * PhysMath::addInQuadriture(error/content,N_Evt_Err/N_Evt);
      if(std::isfinite(newContent) && std::isfinite(newErr)){
        histo->SetBinContent(binX,newContent);
        if(a_propagateNumEvtError){
          histo->SetBinError(binX,newErr);
        }else{
          histo->SetBinError(binX,error/N_Evt);
        }
      }else{
        histo->SetBinContent(binX,0.0);
        histo->SetBinError(binX,0.0);
      }
    }
  }

  histo->Scale(1.0,"width");
  return histo;


}




void ZFitter::createProjectionFile(string a_outFileName, int a_partIndex, int a_detectorIndex, int a_rebinCDFNumBins){

  TH1D* histo_database[16][g_ZFitter_MaxNumRapBins][g_ZFitter_MaxNumMtM0Bins][2];
  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int rapIndex = 0; rapIndex < g_ZFitter_MaxNumRapBins; rapIndex++){
      for(int mTm0Index = 0; mTm0Index < g_ZFitter_MaxNumMtM0Bins; mTm0Index++){
        histo_database[centIndex][rapIndex][mTm0Index][0] = NULL;
        histo_database[centIndex][rapIndex][mTm0Index][1] = NULL;
      }
    }
  }

  int numRapBins = m_rap_mTm0_dEdx_Plus[a_partIndex][0]->GetNbinsX();
  int numMtM0Bins = m_rap_mTm0_dEdx_Plus[a_partIndex][0]->GetNbinsY();
  if(a_detectorIndex == 1){
    numRapBins = m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][0]->GetNbinsX();
    numMtM0Bins = m_rap_mTm0_invBetaBTOF_Plus[a_partIndex][0]->GetNbinsY();
  }else if(a_detectorIndex == 2){
    numRapBins = m_rap_mTm0_invBetaETOF_Plus[a_partIndex][0]->GetNbinsX();
    numMtM0Bins = m_rap_mTm0_invBetaETOF_Plus[a_partIndex][0]->GetNbinsY();
  }

  for(int centIndex = 0; centIndex < 16; centIndex++){
    if(!m_rap_mTm0_dEdx_Plus[a_partIndex][centIndex]) continue;
    for(int rapIndex = 0; rapIndex < numRapBins; rapIndex++){
      for(int mTm0Index = 0; mTm0Index < numMtM0Bins; mTm0Index++){
        histo_database[centIndex][rapIndex][mTm0Index][0] = grabSpecificHistogram(a_partIndex,  1, centIndex, a_detectorIndex, rapIndex+1, mTm0Index+1, a_rebinCDFNumBins);
        histo_database[centIndex][rapIndex][mTm0Index][1] = grabSpecificHistogram(a_partIndex, -1, centIndex, a_detectorIndex, rapIndex+1, mTm0Index+1, a_rebinCDFNumBins);
      }
    }
  }

  TFile* outFile = new TFile(a_outFileName.c_str(),"RECREATE");

  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int rapIndex = 0; rapIndex < g_ZFitter_MaxNumRapBins; rapIndex++){
      for(int mTm0Index = 0; mTm0Index < g_ZFitter_MaxNumMtM0Bins; mTm0Index++){
        HistogramUtilities::ConditionalWrite(histo_database[centIndex][rapIndex][mTm0Index][0]);
        HistogramUtilities::ConditionalWrite(histo_database[centIndex][rapIndex][mTm0Index][1]);
      }
    }
  }





}

//=========================================================================================================================================================================
// NOTE: ZFitter::create2DProjectionFile was intentionally dropped when porting this
// file from lightflavorspectra_etof (see README.md, "What's deliberately NOT in this
// repo" -- it's dead/SDCC-only for this configuration). Its declaration was already
// removed from headers/ZFitter.h; this leftover out-of-line definition (which no
// longer matched any header declaration and broke the ACLiC build) was removed here
// to match.
//=========================================================================================================================================================================


void ZFitter::setParticleInfo(string a_starVer){
  #ifdef _ZFITTER_DEBUG_
    cout << "ZFitter::setParticleInfo(" << a_starVer << ")" << endl;
  #endif

  m_partInfo = new ParticleInfo(a_starVer.c_str(),true,1);
  for(int partIndex = 0; partIndex < 9; partIndex++){
    m_bichselCurves[partIndex] = m_partInfo->GetBichselFunction(partIndex);
    m_invBetaCurves[partIndex] = m_partInfo->GetInverseBetaFunction(partIndex);
  }

  for(int partIndex = 0; partIndex < 9; partIndex++){

    m_invBetaCurves[partIndex] = new TF1(Form("invBeta_part%02d",partIndex),"sqrt([0]*[0] + x*x)/x", 0, 8);
    m_invBetaCurves[partIndex]->SetParameter(0,m_partInfo->GetParticleMass(partIndex));
    m_invBetaCurves[partIndex]->SetNpx(500);

  }
  #ifdef _ZFITTER_DEBUG_
    cout << "End of ZFitter::setParticleInfo(" << a_starVer << ")" << endl;
  #endif
}

//=========================================================================================================================================================================
void ZFitter::setImageDir(string a_prePath, string a_imgDirName){

  m_imagePreDir = a_prePath;
  m_imgDirName = a_imgDirName;
  cout << "m_imagePreDir: " << m_imagePreDir << endl;
  cout << "m_imgDirName: " << m_imgDirName << endl;
  gSystem->mkdir(Form("%s/%s",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  gSystem->mkdir(Form("%s/%s/eta_betaGamma_dEdx",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  gSystem->mkdir(Form("%s/%s/mom_invBetaBTOFFits",m_imagePreDir.c_str(),m_imgDirName.c_str()));
  gSystem->mkdir(Form("%s/%s/mom_invBetaETOFFits",m_imagePreDir.c_str(),m_imgDirName.c_str()));
  gSystem->mkdir(Form("%s/%s/dEdxMomFits",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  gSystem->mkdir(Form("%s/%s/electronStudy",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  gSystem->mkdir(Form("%s/%s/BichselShifts",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  gSystem->mkdir(Form("%s/%s/ParticleWidths",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  /*for(int partIndex = 0; partIndex < 9; partIndex++){
    if(!m_particleLoaded[partIndex]) continue;
    gSystem->mkdir(Form("%s/%s/%s",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    gSystem->mkdir(Form("%s/%s/%s/dEdxFits",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    gSystem->mkdir(Form("%s/%s/%s/invBetaBTOFFits",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    if(m_hasETOF) gSystem->mkdir(Form("%s/%s/%s/invBetaETOFFits",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      gSystem->mkdir(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data(),centIndex));
      if(m_hasETOF) gSystem->mkdir(Form("%s/%s/%s/invBetaETOFFits_Cent%02d",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data(),centIndex));
      gSystem->mkdir(Form("%s/%s/%s/dEdxFits_Cent%02d",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data(),centIndex));
    }
    //gSystem->mkdir(Form("%s/%s/%s/ZTPC_Zoomed",m_imagePreDir.c_str(), m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    //gSystem->mkdir(Form("%s/%s/%s/BTOF_Zoomed",m_imagePreDir.c_str(), m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    //if(m_hasETOF) gSystem->mkdir(Form("%s/%s/%s/ETOF_Zoomed",m_imagePreDir.c_str(), m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    //gSystem->mkdir(Form("%s/%s/%s/eta_betaGamma_dEdx",m_imagePreDir.c_str(), m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));
    gSystem->mkdir(Form("%s/%s/%s/ZTPC_BTOFIsolated/",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));

    for(int partIndex2 = 0; partIndex2 < 9; partIndex2++){
      gSystem->mkdir(Form("%s/%s/%s/ZTPC_BTOFIsolated/%s",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data(),m_partInfo->GetParticleName(partIndex2,0).Data()));
    }
    gSystem->mkdir(Form("%s/%s/%s/ZTPC_BTOFIsolated/PionElectron",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,0).Data()));

  }*/

  for(int chargePM = -1; chargePM < 2; chargePM+=2){
    for(int partIndex = 0; partIndex < 9; partIndex++){
      if(!m_particleLoaded[partIndex]) continue;
      gSystem->mkdir(Form("%s/%s/%s",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      gSystem->mkdir(Form("%s/%s/%s/dEdxFits",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      gSystem->mkdir(Form("%s/%s/%s/dEdxFits_Overview",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      gSystem->mkdir(Form("%s/%s/%s/dEdxParamExtrapolations",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      gSystem->mkdir(Form("%s/%s/%s/invBetaBTOFFits",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      gSystem->mkdir(Form("%s/%s/%s/invBetaBTOFFits_Overview",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      gSystem->mkdir(Form("%s/%s/%s/BTOFParamExtrapolations",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      gSystem->mkdir(Form("%s/%s/%s/TPCParamExtrapolations",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      if(m_hasETOF){
        gSystem->mkdir(Form("%s/%s/%s/invBetaETOFFits",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
        gSystem->mkdir(Form("%s/%s/%s/invBetaETOFFits_Overview",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      }
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        gSystem->mkdir(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
        if(m_hasETOF) gSystem->mkdir(Form("%s/%s/%s/invBetaETOFFits_Cent%02d",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
        gSystem->mkdir(Form("%s/%s/%s/dEdxFits_Cent%02d",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
      }
      gSystem->mkdir(Form("%s/%s/%s/ZTPC_BTOFIsolated/",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));
      for(int partIndex2 = 0; partIndex2 < 9; partIndex2++){
        gSystem->mkdir(Form("%s/%s/%s/ZTPC_BTOFIsolated/%s",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data(),m_partInfo->GetParticleName(partIndex2,chargePM).Data()));

      }
      gSystem->mkdir(Form("%s/%s/%s/ZTPC_BTOFIsolated/PionElectron",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(partIndex,chargePM).Data()));

    }
  }
}

void ZFitter::setTextFileDir(string a_textFileDir){
  m_textFileDir = a_textFileDir;
  gSystem->mkdir(a_textFileDir.c_str());
};

//=========================================================================================================================================================================

void ZFitter::setMimicInfo(int* a_partIndicies, double* a_momentumStart, double* a_momentumStop){
  cout << "MimicInfo: " << endl;
  for(int partIndex = 0; partIndex < 9 ; partIndex++){
    m_particleMimicIndicies[partIndex] = a_partIndicies[partIndex];
    m_particleMimicLowMom[partIndex] = a_momentumStart[partIndex];
    m_particleMimicHighMom[partIndex] = a_momentumStop[partIndex];
    cout << "   " << m_particleMimicIndicies[partIndex] << " " << m_particleMimicLowMom[partIndex] << " " << m_particleMimicHighMom[partIndex] << endl;
  }
}

//=========================================================================================================================================================================

void ZFitter::setFullMimicToggle(bool* a_fullParticleMimic){
  cout << "Full Mimic Toggle: " << endl;
  for(int partIndex = 0; partIndex < 9 ; partIndex++){
    m_fullParticleMimic[partIndex] = a_fullParticleMimic[partIndex];
    cout << "  " << m_fullParticleMimic[partIndex] << endl;
  }
}

//=========================================================================================================================================================================
void ZFitter::setElectronParams(double a_mu, double a_sig, double a_mean, double a_stdDev){
  m_electronParams[0] = a_mu;
  m_electronParams[1] = a_sig;
  m_electronParams[2] = a_mean;
  m_electronParams[3] = a_stdDev;
  m_electronParams[4] = TMath::Exp(a_mu+a_sig*a_sig/2.0);
  m_electronParams[5] = (TMath::Exp(a_sig*a_sig) - 1.0)
  *TMath::Exp(2.0*a_mu + a_sig*a_sig);
  cout << "Electron Paramters: ";
  for(int iii = 0; iii < 6; iii++) cout << m_electronParams[iii] << "  ";
  cout << endl;
}

//=========================================================================================================================================================================

//void ZFitter::setMaxForTotalFitMomenta(double* a_momentaPlus, double* a_momentaMinus){
//  for(int iii = 0; iii < 9; iii++){
//    m_dEdxMaximumMomForTotalFit[iii][0] = a_momentaPlus[iii];
//    m_dEdxMaximumMomForTotalFit[iii][1] = a_momentaMinus[iii];
//  }
//
//}

//=========================================================================================================================================================================


void ZFitter::makeFitDataHistos(int a_partIndex, int a_detectorIndex){
  for(int partBaseIndex = 0; partBaseIndex < m_numParticles; partBaseIndex++){
      cout << "Make fir data histos for partbaseindex " << partBaseIndex << " loaded? " << m_particleLoaded[partBaseIndex] << endl;
    if(!m_particleLoaded[partBaseIndex]) continue;
    if(a_partIndex >= 0 && a_partIndex != partBaseIndex) continue;
    //if(partBaseIndex ==3 || partBaseIndex == 8) continue;
    //TH2D* structHisto = HistogramUtilities::projection2D(m_rap_mTm0_dEdx_Plus[partBaseIndex][0], 1, 1);
    //HistogramUtilities::zeroBins2D(structHisto);
    //TH2D* structHistoInvBeta = HistogramUtilities::projection2D(m_rap_mTm0_invBetaBTOF_Plus[partBaseIndex][0], 1, 1);
    //HistogramUtilities::zeroBins2D(structHistoInvBeta);
    //TH2D* structHistoInvBetaETOF = HistogramUtilities::projection2D(m_rap_mTm0_invBetaETOF_Plus[partBaseIndex][0], 1, 1);
    //HistogramUtilities::zeroBins2D(structHistoInvBetaETOF);


    string partSymbol = m_partInfo->GetParticleSymbol(partBaseIndex).Data();
    TH2D* structHisto[2] = {nullptr,nullptr};
    if(a_detectorIndex < 0 || a_detectorIndex == 0){
      structHisto[0]  = HistogramUtilities::projection2D(m_rap_mTm0_dEdx_Plus[partBaseIndex][0], 1, 1);
      structHisto[1] = HistogramUtilities::projection2D(m_rap_mTm0_dEdx_Minus[partBaseIndex][0], 1, 1);
      HistogramUtilities::zeroBins2D(structHisto[0]);
      HistogramUtilities::zeroBins2D(structHisto[1]);
      structHisto[0]->SetTitle(Form("structHistoPlus; y_{%s}; m_{T}-m_{%s} [GeV/c^{2}]",partSymbol.c_str(),partSymbol.c_str()));
      structHisto[1]->SetTitle(Form("structHistoMinus; y_{%s}; m_{T}-m_{%s} [GeV/c^{2}]",partSymbol.c_str(),partSymbol.c_str()));
    }

    TH2D* structHistoInvBeta[2] = {nullptr,nullptr};
    if(a_detectorIndex < 0 || a_detectorIndex == 1){
      structHistoInvBeta[0] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaBTOF_Plus[partBaseIndex][0], 1, 1);
      structHistoInvBeta[1] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaBTOF_Minus[partBaseIndex][0], 1, 1);
      HistogramUtilities::zeroBins2D(structHistoInvBeta[0]);
      HistogramUtilities::zeroBins2D(structHistoInvBeta[1]);
      structHistoInvBeta[0]->SetTitle(Form("structHistoInvBetaPlus; y_{%s}; m_{T}-m_{%s} [GeV/c^{2}]",partSymbol.c_str(),partSymbol.c_str()));
      structHistoInvBeta[1]->SetTitle(Form("structHistoInvBetaMinus; y_{%s}; m_{T}-m_{%s} [GeV/c^{2}]",partSymbol.c_str(),partSymbol.c_str()));
    }

    TH2D* structHistoInvBetaETOF[2] = {nullptr,nullptr};
    if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
      structHistoInvBetaETOF[0] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaETOF_Plus[partBaseIndex][0], 1, 1);
      structHistoInvBetaETOF[1] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaETOF_Minus[partBaseIndex][0], 1, 1);
      HistogramUtilities::zeroBins2D(structHistoInvBetaETOF[0]);
      HistogramUtilities::zeroBins2D(structHistoInvBetaETOF[1]);
      structHistoInvBetaETOF[0]->SetTitle(Form("structHistoInvBetaETOFPlus; y_{%s}; m_{T}-m_{%s} [GeV/c^{2}]",partSymbol.c_str(),partSymbol.c_str()));
      structHistoInvBetaETOF[1]->SetTitle(Form("structHistoInvBetaETOFMinus; y_{%s}; m_{T}-m_{%s} [GeV/c^{2}]",partSymbol.c_str(),partSymbol.c_str()));
    }


    if(a_detectorIndex < 0 || a_detectorIndex == 0){
      m_ZTPC_LowRange[partBaseIndex][0] = (TH2D*) structHisto[0]->Clone();
      m_ZTPC_LowRange[partBaseIndex][0]->SetName("RangeLowZTPC_Plus");
      m_ZTPC_LowRange[partBaseIndex][0]->SetTitle("Low Fit Range for ZTPC Plus");
      m_ZTPC_HighRange[partBaseIndex][0] = (TH2D*) structHisto[0]->Clone();
      m_ZTPC_HighRange[partBaseIndex][0]->SetName("RangeHighZTPC_Plus");
      m_ZTPC_HighRange[partBaseIndex][0]->SetTitle("High Fit Range for ZTPC Plus");
      m_ZTPC_LowRange[partBaseIndex][1] = (TH2D*) structHisto[1]->Clone();
      m_ZTPC_LowRange[partBaseIndex][1]->SetName("RangeLowZTPC_Minus");
      m_ZTPC_LowRange[partBaseIndex][1]->SetTitle("Low Fit Range for ZTPC Minus");
      m_ZTPC_HighRange[partBaseIndex][1] = (TH2D*) structHisto[1]->Clone();
      m_ZTPC_HighRange[partBaseIndex][1]->SetName("RangeHighZTPC_Minus");
      m_ZTPC_HighRange[partBaseIndex][1]->SetTitle("High Fit Range for ZTPC Minus");
      m_ChiSqr_ZTPC[partBaseIndex][0] = (TH2D*) structHisto[0]->Clone();
      m_ChiSqr_ZTPC[partBaseIndex][0]->SetName("dEdxChiSqrPlus_AllCent");
      m_ChiSqr_ZTPC[partBaseIndex][0]->SetTitle("dE/dx Plus Chi Squared All Centrialities");
      m_ChiSqr_ZTPC[partBaseIndex][1] = (TH2D*) structHisto[1]->Clone();
      m_ChiSqr_ZTPC[partBaseIndex][1]->SetName("dEdxChiSqrMinus_AllCent");
      m_ChiSqr_ZTPC[partBaseIndex][1]->SetTitle("dE/dx Minus Chi Squared All Centrialities");
      m_FitStatus_ZTPC[partBaseIndex][0] = (TH2D*) structHisto[0]->Clone();
      m_FitStatus_ZTPC[partBaseIndex][0]->SetName("dEdxFitStatusPlus_AllCent");
      m_FitStatus_ZTPC[partBaseIndex][0]->SetTitle("dE/dx Plus Fit Status All Centrialities");
      m_FitStatus_ZTPC[partBaseIndex][1] = (TH2D*) structHisto[1]->Clone();
      m_FitStatus_ZTPC[partBaseIndex][1]->SetName("dEdxFitStatusMinus_AllCent");
      m_FitStatus_ZTPC[partBaseIndex][1]->SetTitle("dE/dx Minus Fit Status All Centrialities");
    }
    if(a_detectorIndex < 0 || a_detectorIndex == 1){
      m_ZbTOF_LowRange[partBaseIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
      m_ZbTOF_LowRange[partBaseIndex][0]->SetName("RangeLowBTOF_Plus");
      m_ZbTOF_LowRange[partBaseIndex][0]->SetTitle("Low Fit Range for bTOF Plus");
      m_ZbTOF_HighRange[partBaseIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
      m_ZbTOF_HighRange[partBaseIndex][0]->SetName("RangeHighBTOF_Plus");
      m_ZbTOF_HighRange[partBaseIndex][0]->SetTitle("High Fit Range for bTOF Plus");
      m_ZbTOF_LowRange[partBaseIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
      m_ZbTOF_LowRange[partBaseIndex][1]->SetName("RangeLowBTOF_Minus");
      m_ZbTOF_LowRange[partBaseIndex][1]->SetTitle("Low Fit Range for bTOF Minus");
      m_ZbTOF_HighRange[partBaseIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
      m_ZbTOF_HighRange[partBaseIndex][1]->SetName("RangeHighBTOF_Minus");
      m_ZbTOF_HighRange[partBaseIndex][1]->SetTitle("High Fit Range for bTOF Minus");
      m_ChiSqr_ZbTOF[partBaseIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
      m_ChiSqr_ZbTOF[partBaseIndex][0]->SetName("invBetaChiSqr_AllCent_Plus");
      m_ChiSqr_ZbTOF[partBaseIndex][0]->SetTitle("invBeta BTOF Plus Chi Squared All Centrialities");
      m_ChiSqr_ZbTOF[partBaseIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
      m_ChiSqr_ZbTOF[partBaseIndex][1]->SetName("invBetaChiSqr_AllCent_Minus");
      m_ChiSqr_ZbTOF[partBaseIndex][1]->SetTitle("invBeta BTOF Minus Chi Squared All Centrialities");
      m_FitStatus_ZbTOF[partBaseIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
      m_FitStatus_ZbTOF[partBaseIndex][0]->SetName("InvBetaBTOFFitStatusPlus_AllCent");
      m_FitStatus_ZbTOF[partBaseIndex][0]->SetTitle("InvBetaBTOF Plus Fit Status All Centrialities");
      m_FitStatus_ZbTOF[partBaseIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
      m_FitStatus_ZbTOF[partBaseIndex][1]->SetName("InvBetaBTOFFitStatusMinus_AllCent");
      m_FitStatus_ZbTOF[partBaseIndex][1]->SetTitle("InvBetaBTOF Minus Fit Status All Centrialities");
    }



    if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
      m_ChiSqr_ZeTOF[partBaseIndex][0] = (TH2D*) structHistoInvBetaETOF[0]->Clone();
      m_ChiSqr_ZeTOF[partBaseIndex][0]->SetName("invBetaETOFChiSqrPlus_AllCent");
      m_ChiSqr_ZeTOF[partBaseIndex][0]->SetTitle("invBeta ETOF Plus Chi Squared All Centrialities");
      m_ChiSqr_ZeTOF[partBaseIndex][1] = (TH2D*) structHistoInvBetaETOF[1]->Clone();
      m_ChiSqr_ZeTOF[partBaseIndex][1]->SetName("invBetaETOFChiSqrMinus_AllCent");
      m_ChiSqr_ZeTOF[partBaseIndex][1]->SetTitle("InvBetaETOF Minus Chi Squared All Centrialities");
      m_FitStatus_ZeTOF[partBaseIndex][0] = (TH2D*) structHistoInvBetaETOF[0]->Clone();
      m_FitStatus_ZeTOF[partBaseIndex][0]->SetName("InvBetaETOFFitStatusPlus_AllCent");
      m_FitStatus_ZeTOF[partBaseIndex][0]->SetTitle("InvBetaETOF Plus Fit Status All Centrialities");
      m_FitStatus_ZeTOF[partBaseIndex][1] = (TH2D*) structHistoInvBetaETOF[1]->Clone();
      m_FitStatus_ZeTOF[partBaseIndex][1]->SetName("InvBetaETOFFitStatusMinus_AllCent");
      m_FitStatus_ZeTOF[partBaseIndex][1]->SetTitle("InvBetaETOF Minus Fit Status All Centrialities");
    }


    for(int partIndex = 0; partIndex < m_numParticles; partIndex ++){
      cout << "Making Fit Data Histograms for " << m_partInfo->GetParticleName(partIndex,0) << " in " << m_partInfo->GetParticleName(partBaseIndex,0) << " Space" << endl;

      for(int pmIndex=0; pmIndex < 2; pmIndex++){
        int chargePM = 1;
        if(pmIndex == 1){
          chargePM = -1;
          if(!m_partToFitMatrix_Minus[partBaseIndex][partIndex]) continue;
        }else{
          if(!m_partToFitMatrix_Plus[partBaseIndex][partIndex]) continue;
        }

        if(a_detectorIndex < 0 || a_detectorIndex == 0){
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][0][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][0][pmIndex]->SetName(Form("dEdx_Amp_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][0][pmIndex]->SetTitle(Form("dEdx Amp for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][1][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][1][pmIndex]->SetName(Form("dEdx_Mean_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][1][pmIndex]->SetTitle(Form("dEdx Mean for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][2][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][2][pmIndex]->SetName(Form("dEdx_StdDev_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][2][pmIndex]->SetTitle(Form("dEdx StdDev for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][3][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][3][pmIndex]->SetName(Form("dEdx_NDF_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][3][pmIndex]->SetTitle(Form("dEdx NDF for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][4][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][4][pmIndex]->SetName(Form("dEdx_p5_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZTPC[partBaseIndex][partIndex][4][pmIndex]->SetTitle(Form("dEdx pearsonIV param 5 (nu) for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));

          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][0][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][0][pmIndex]->SetName(Form("dEdx_Solo_Amp_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][0][pmIndex]->SetTitle(Form("dEdx_Solo Amp for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][1][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][1][pmIndex]->SetName(Form("dEdx_Solo_Mean_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][1][pmIndex]->SetTitle(Form("dEdx_Solo Mean for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][2][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][2][pmIndex]->SetName(Form("dEdx_Solo_StdDev_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][2][pmIndex]->SetTitle(Form("dEdx StdDev for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][3][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][3][pmIndex]->SetName(Form("dEdx_Solo_NDF_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][3][pmIndex]->SetTitle(Form("dEdx NDF for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));

          m_Fit_Data_Solo_ZTPCShifts[partBaseIndex][partIndex][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
          m_Fit_Data_Solo_ZTPCShifts[partBaseIndex][partIndex][pmIndex]->SetName(Form("dEdx_Solo_ZTPCShifts_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_Solo_ZTPCShifts[partBaseIndex][partIndex][pmIndex]->SetTitle(Form("dEdx Shift from Prediction for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
        }

        if(a_detectorIndex < 0 || a_detectorIndex == 1){
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][0][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][0][pmIndex]->SetName(Form("invBeta_BTOF_Amp_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][0][pmIndex]->SetTitle(Form("InvBeta BTOF Amp for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][1][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][1][pmIndex]->SetName(Form("invBeta_BTOF_Mean_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][1][pmIndex]->SetTitle(Form("InvBeta BTOF Mean for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][2][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][2][pmIndex]->SetName(Form("invBeta_BTOF_StdDev_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][2][pmIndex]->SetTitle(Form("InvBeta BTOF StdDev for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][3][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][3][pmIndex]->SetName(Form("invBeta_BTOF_NDF_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][3][pmIndex]->SetTitle(Form("InvBeta BTOF NDF for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][4][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][4][pmIndex]->SetName(Form("invBeta_BTOF_p5_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZbTOF[partBaseIndex][partIndex][4][pmIndex]->SetTitle(Form("InvBeta BTOF p5 for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
        }



        if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][0][pmIndex] = (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][0][pmIndex]->SetName(Form("invBeta_ETOF_Amp_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][0][pmIndex]->SetTitle(Form("InvBeta ETOF Amp for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][1][pmIndex] = (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][1][pmIndex]->SetName(Form("invBeta_ETOF_Mean_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][1][pmIndex]->SetTitle(Form("InvBeta ETOF Mean for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][2][pmIndex]= (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][2][pmIndex]->SetName(Form("invBeta_ETOF_StdDev_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][2][pmIndex]->SetTitle(Form("InvBeta ETOF StdDev for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][3][pmIndex]= (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][3][pmIndex]->SetName(Form("invBeta_ETOF_NDF_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][3][pmIndex]->SetTitle(Form("InvBeta ETOF NDF for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][4][pmIndex]= (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][4][pmIndex]->SetName(Form("invBeta_ETOF_p5_%s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
          m_Fit_Data_ZeTOF[partBaseIndex][partIndex][4][pmIndex]->SetTitle(Form("InvBeta ETOF p5 for %s",m_partInfo->GetParticleName(partIndex,chargePM).Data()));
        }

        for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
          if(a_detectorIndex < 0 || a_detectorIndex == 0){
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][0][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][0][pmIndex]->SetName(Form("dEdx_Amp_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][0][pmIndex]->SetTitle(Form("dEdx Amp for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][1][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][1][pmIndex]->SetName(Form("dEdx_Mean_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][1][pmIndex]->SetTitle(Form("dEdx Mean for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][2][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][2][pmIndex]->SetName(Form("dEdx_StdDev_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][2][pmIndex]->SetTitle(Form("dEdx StdDev for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][3][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][3][pmIndex]->SetName(Form("dEdx_NDF_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][3][pmIndex]->SetTitle(Form("dEdx NDF for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][4][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][4][pmIndex]->SetName(Form("dEdx_p5_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][4][pmIndex]->SetTitle(Form("dEdx pearsonIV param 5 (nu) for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
          }
          if(a_detectorIndex < 0 || a_detectorIndex == 1){
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][0][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][0][pmIndex]->SetName(Form("invBeta_BTOF_Amp_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][0][pmIndex]->SetTitle(Form("InvBeta BTOF Amp for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][1][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][1][pmIndex]->SetName(Form("invBeta_BTOF_Mean_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][1][pmIndex]->SetTitle(Form("InvBeta BTOF Mean for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][2][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][2][pmIndex]->SetName(Form("invBeta_BTOF_StdDev_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][2][pmIndex]->SetTitle(Form("InvBeta BTOF StdDev for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][3][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][3][pmIndex]->SetName(Form("invBeta_BTOF_NDF_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][3][pmIndex]->SetTitle(Form("InvBeta BTOF NDF for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][4][pmIndex]= (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][4][pmIndex]->SetName(Form("invBeta_BTOF_p5_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZbTOF[partBaseIndex][partIndex][centIndex][4][pmIndex]->SetTitle(Form("InvBeta BTOF p5 for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
          }
          if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][0][pmIndex] = (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][0][pmIndex]->SetName(Form("invBeta_ETOF_Amp_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][0][pmIndex]->SetTitle(Form("InvBeta ETOF Amp for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][1][pmIndex] = (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][1][pmIndex]->SetName(Form("invBeta_ETOF_Mean_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][1][pmIndex]->SetTitle(Form("InvBeta ETOF Mean for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][2][pmIndex]= (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][2][pmIndex]->SetName(Form("invBeta_ETOF_StdDev_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][2][pmIndex]->SetTitle(Form("InvBeta ETOF StdDev for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][3][pmIndex]= (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][3][pmIndex]->SetName(Form("invBeta_ETOF_NDF_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][3][pmIndex]->SetTitle(Form("InvBeta ETOF NDF for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][4][pmIndex]= (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][4][pmIndex]->SetName(Form("invBeta_ETOF_p5_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
            m_Fit_Data_Cent_ZeTOF[partBaseIndex][partIndex][centIndex][4][pmIndex]->SetTitle(Form("InvBeta ETOF p5 for %s Cent %02d",m_partInfo->GetParticleName(partIndex,chargePM).Data(),centIndex));
          }
        } // end cent loop

      } // end plus minus loop
      cout << "Done Making Fit Data Histograms " << endl;



      cout << "Making chi sqrd / fit status / spectra histograms : " << structHisto << "  " << structHistoInvBeta << "  " << structHistoInvBetaETOF << endl;
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

        if(a_detectorIndex < 0 || a_detectorIndex == 0){
          TH2D* centPlusHistoDeDx = (TH2D*) structHisto[0]->Clone();
          centPlusHistoDeDx->SetName(Form("dEdx_Amp_%s_Cent_%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex));
          centPlusHistoDeDx->SetTitle(Form("dE/dx  Amp for %s Centrality %02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex));
          m_dEdxYieldsCentPlus[partBaseIndex][partIndex][centIndex] = centPlusHistoDeDx;
          TH2D* centMinusHistoDeDx= (TH2D*) structHisto[1]->Clone();
          centMinusHistoDeDx->SetName(Form("dEdx_Amp_%s_Cent_%02d",m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex));
          centMinusHistoDeDx->SetTitle(Form("dE/dx  Amp for %s Centrality %02d",m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex));
          m_dEdxYieldsCentMinus[partBaseIndex][partIndex][centIndex] = centMinusHistoDeDx;
        }
        if(a_detectorIndex < 0 || a_detectorIndex == 1){
          TH2D* centPlusHisto = (TH2D*) structHistoInvBeta[0]->Clone();
          centPlusHisto->SetName(Form("invBeta_BTOF_Amp_%s_Cent_%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex));
          centPlusHisto->SetTitle(Form("InvBeta BTOF  Amp for %s Centrality %02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex));
          m_InvBetaBTOFYieldsCentPlus[partBaseIndex][partIndex][centIndex] = centPlusHisto;
          TH2D* centMinusHisto = (TH2D*) structHistoInvBeta[1]->Clone();
          centMinusHisto->SetName(Form("invBeta_BTOF_Amp_%s_Cent_%02d",m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex));
          centMinusHisto->SetTitle(Form("InvBeta BTOF  Amp for %s Centrality %02d",m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex));
          m_InvBetaBTOFYieldsCentMinus[partBaseIndex][partIndex][centIndex] = centMinusHisto;
        }




        for(int pmIndex = 0; pmIndex < 2; pmIndex++){
          int chargePM = 1;
          if(pmIndex == 1) chargePM = -1;
          //####   Chi Squared ######
          if(a_detectorIndex < 0 || a_detectorIndex == 0){
            m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][pmIndex]->SetName(Form("ZTPC_ChiSqr_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("ZTPC %s Chi Squared Cent Bin %d", m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][pmIndex]->SetName(Form("ZTPC_FitStatus_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("ZTPC %s Fit Status Cent Bin %d", m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHisto[pmIndex]->Clone();
            m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][pmIndex]->SetName(Form("dEdxSpectra_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(),centIndex));
            m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("dE/dx Amp for %s of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(),centIndex));
          }
          if(a_detectorIndex < 0 || a_detectorIndex == 1){
            m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex]->SetName(Form("ZbTOF_ChiSqr_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("ZbTOF %s Chi Squared Cent Bin %d", m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex]->SetName(Form("ZbTOF_FitStatus_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("ZbTOF %s Fit Status Cent Bin %d", m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHistoInvBeta[pmIndex]->Clone();
            m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex]->SetName(Form("InvBetaBTOFSpectra_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(),centIndex));
            m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("Inverse Beta BTOF Amp for %s  of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(),centIndex));
          }
          if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
            m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex]->SetName(Form("ZbTOF_ChiSqr_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("ZbTOF %s Chi Squared Cent Bin %d", m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex]->SetName(Form("ZeTOF_FitStatus_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("ZeTOF %s Fit Status Cent Bin %d", m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(), centIndex));
            m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex] = (TH2D*) structHistoInvBetaETOF[pmIndex]->Clone();
            m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex]->SetName(Form("InvBetaETOFSpectra_%s_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(),centIndex));
            m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][pmIndex]->SetTitle(Form("Inverse Beta ETOF Amp for %s of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,chargePM).Data(),centIndex));
          }
        }




      } // end cent loop

    } // partIndex

    #ifndef _SAVE_POINTERS_
    if(structHisto[0]) delete structHisto[0];
    if(structHisto[1]) delete structHisto[1];
    if(structHistoInvBeta[0]) delete structHistoInvBeta[0];
    if(structHistoInvBeta[1]) delete structHistoInvBeta[1];
    if(m_hasETOF){
      if(structHistoInvBetaETOF[0]) delete structHistoInvBetaETOF[0];
      if(structHistoInvBetaETOF[1]) delete structHistoInvBetaETOF[1];
    }
    #endif
  } // partBaseIndex

  cout << " Making Beta Gamma Histograms " << endl;
  //Setup Fit Data Structures for BetaGamma dEdx Fitting
  if(!m_betaGammaBTOF_dEdx_Plus){
    cout << "  ERROR: Unable to locate betaGammaBtof_dEdx_Plus...    Your basic histogram file must be wrong." << endl;
  }else{
    TH1D* structHisto_betagamma_BTOF = m_betaGammaBTOF_dEdx_Plus->ProjectionX("structBetaGammaBTOFNoEta",1,1);
    HistogramUtilities::zeroBins(structHisto_betagamma_BTOF);
    TH1D* structHisto_betagamma_ETOF = NULL;
    if(m_hasETOF){
      structHisto_betagamma_ETOF = m_betaGammaETOF_dEdx_Plus->ProjectionX("structBetaGammaETOFNoEta",1,1);
      HistogramUtilities::zeroBins(structHisto_betagamma_ETOF);
    }

    m_betaGammaBTOF_DeDx_Data_Plus[0] = (TH1D*) structHisto_betagamma_BTOF->Clone();
    m_betaGammaBTOF_DeDx_Data_Plus[0]->SetName("betaGammaDeDx_BTOF_Plus_N");
    m_betaGammaBTOF_DeDx_Data_Plus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Plus; #beta#gamma; N");
    m_betaGammaBTOF_DeDx_Data_Plus[1] = (TH1D*) structHisto_betagamma_BTOF->Clone();
    m_betaGammaBTOF_DeDx_Data_Plus[1]->SetName("betaGammaDeDx_BTOF_Plus_Mu");
    m_betaGammaBTOF_DeDx_Data_Plus[1]->SetTitle("#mu for #beta#gamma dE/dx Plus; #beta#gamma; #mu");
    m_betaGammaBTOF_DeDx_Data_Plus[2] = (TH1D*) structHisto_betagamma_BTOF->Clone();
    m_betaGammaBTOF_DeDx_Data_Plus[2]->SetName("betaGammaDeDx_BTOF_Plus_Sig");
    m_betaGammaBTOF_DeDx_Data_Plus[2]->SetTitle("#sigma for #beta#gamma dE/dx Plus; #beta#gamma; #sigma");

    m_betaGammaBTOF_DeDx_Data_Minus[0] = (TH1D*) structHisto_betagamma_BTOF->Clone();
    m_betaGammaBTOF_DeDx_Data_Minus[0]->SetName("betaGammaDeDx_BTOF_Minus_N");
    m_betaGammaBTOF_DeDx_Data_Minus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Minus; #beta#gamma; N");
    m_betaGammaBTOF_DeDx_Data_Minus[1] = (TH1D*) structHisto_betagamma_BTOF->Clone();
    m_betaGammaBTOF_DeDx_Data_Minus[1]->SetName("betaGammaDeDx_BTOF_Minus_Mu");
    m_betaGammaBTOF_DeDx_Data_Minus[1]->SetTitle("#mu for #beta#gamma dE/dx Minus; #beta#gamma; #mu");
    m_betaGammaBTOF_DeDx_Data_Minus[2] = (TH1D*) structHisto_betagamma_BTOF->Clone();
    m_betaGammaBTOF_DeDx_Data_Minus[2]->SetName("betaGammaDeDx_BTOF_Minus_Sig");
    m_betaGammaBTOF_DeDx_Data_Minus[2]->SetTitle("#sigma for #beta#gamma dE/dx Minus; #beta#gamma; #sigma");

    if(m_hasETOF){
      m_betaGammaETOF_DeDx_Data_Plus[0] = (TH1D*) structHisto_betagamma_ETOF->Clone();
      m_betaGammaETOF_DeDx_Data_Plus[0]->SetName("betaGammaDeDx_ETOF_Plus_N");
      m_betaGammaETOF_DeDx_Data_Plus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Plus; #beta#gamma; N");
      m_betaGammaETOF_DeDx_Data_Plus[1] = (TH1D*) structHisto_betagamma_ETOF->Clone();
      m_betaGammaETOF_DeDx_Data_Plus[1]->SetName("betaGammaDeDx_ETOF_Plus_Mu");
      m_betaGammaETOF_DeDx_Data_Plus[1]->SetTitle("#mu for #beta#gamma dE/dx Plus; #beta#gamma; #mu");
      m_betaGammaETOF_DeDx_Data_Plus[2] = (TH1D*) structHisto_betagamma_ETOF->Clone();
      m_betaGammaETOF_DeDx_Data_Plus[2]->SetName("betaGammaDeDx_ETOF_Plus_Sig");
      m_betaGammaETOF_DeDx_Data_Plus[2]->SetTitle("#sigma for #beta#gamma dE/dx Plus; #beta#gamma; #sigma");

      m_betaGammaETOF_DeDx_Data_Minus[0] = (TH1D*) structHisto_betagamma_ETOF->Clone();
      m_betaGammaETOF_DeDx_Data_Minus[0]->SetName("betaGammaDeDx_ETOF_Minus_N");
      m_betaGammaETOF_DeDx_Data_Minus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Minus; #beta#gamma; N");
      m_betaGammaETOF_DeDx_Data_Minus[1] = (TH1D*) structHisto_betagamma_ETOF->Clone();
      m_betaGammaETOF_DeDx_Data_Minus[1]->SetName("betaGammaDeDx_ETOF_Minus_Mu");
      m_betaGammaETOF_DeDx_Data_Minus[1]->SetTitle("#mu for #beta#gamma dE/dx Minus; #beta#gamma; #mu");
      m_betaGammaETOF_DeDx_Data_Minus[2] = (TH1D*) structHisto_betagamma_ETOF->Clone();
      m_betaGammaETOF_DeDx_Data_Minus[2]->SetName("betaGammaDeDx_ETOF_Minus_Sig");
      m_betaGammaETOF_DeDx_Data_Minus[2]->SetTitle("#sigma for #beta#gamma dE/dx Minus; #beta#gamma; #sigma");
    }
  }

  cout << " making eta beta gamma stuff" << endl;
  if(m_eta_betaGammaBtof_dEdx_Plus){
    TH2I* structHisto_eta_betagamma_BTOF_int = HistogramUtilities::projection2D(m_eta_betaGammaBtof_dEdx_Plus, 1, 1);
    TH2D* structHisto_eta_betagamma_BTOF     = HistogramUtilities::convertInt_to_DoubleHisto2D(structHisto_eta_betagamma_BTOF_int);
    #ifndef _SAVE_POINTERS_
    delete structHisto_eta_betagamma_BTOF_int;
    #endif
    HistogramUtilities::zeroBins2D(structHisto_eta_betagamma_BTOF);

    TH2D* structHisto_eta_betagamma_ETOF = NULL;
    if(m_hasETOF){
       TH2I* structHisto_eta_betagamma_ETOF_int = HistogramUtilities::projection2D(m_eta_betaGammaEtof_dEdx_Plus, 1, 1);
       structHisto_eta_betagamma_ETOF = HistogramUtilities::convertInt_to_DoubleHisto2D(structHisto_eta_betagamma_ETOF_int);
       #ifndef _SAVE_POINTERS_
       delete structHisto_eta_betagamma_ETOF_int;
       #endif
       HistogramUtilities::zeroBins2D(structHisto_eta_betagamma_ETOF);

      m_eta_betaGammaETOF_DeDx_Data_Plus[0] = (TH2D*) structHisto_eta_betagamma_ETOF->Clone();
      m_eta_betaGammaETOF_DeDx_Data_Plus[0]->SetName("etaBetaGammaDeDx_ETOF_Plus_N");
      m_eta_betaGammaETOF_DeDx_Data_Plus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Plus;#eta;#beta#gamma");
      m_eta_betaGammaETOF_DeDx_Data_Plus[1] = (TH2D*) structHisto_eta_betagamma_ETOF->Clone();
      m_eta_betaGammaETOF_DeDx_Data_Plus[1]->SetName("etaBetaGammaDeDx_ETOF_Plus_Mu");
      m_eta_betaGammaETOF_DeDx_Data_Plus[1]->SetTitle("#mu for #beta#gamma dE/dx Plus;#eta;#beta#gamma");
      m_eta_betaGammaETOF_DeDx_Data_Plus[2] = (TH2D*) structHisto_eta_betagamma_ETOF->Clone();
      m_eta_betaGammaETOF_DeDx_Data_Plus[2]->SetName("etaBetaGammaDeDx_ETOF_Plus_Sig");
      m_eta_betaGammaETOF_DeDx_Data_Plus[2]->SetTitle("#sigma for #beta#gamma dE/dx Plus;#eta;#beta#gamma");

      m_eta_betaGammaETOF_DeDx_Data_Minus[0] = (TH2D*) structHisto_eta_betagamma_ETOF->Clone();
      m_eta_betaGammaETOF_DeDx_Data_Minus[0]->SetName("etaBetaGammaDeDx_ETOF_Minus_N");
      m_eta_betaGammaETOF_DeDx_Data_Minus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Minus;#eta;#beta#gamma");
      m_eta_betaGammaETOF_DeDx_Data_Minus[1] = (TH2D*) structHisto_eta_betagamma_ETOF->Clone();
      m_eta_betaGammaETOF_DeDx_Data_Minus[1]->SetName("etaBetaGammaDeDx_ETOF_Minus_Mu");
      m_eta_betaGammaETOF_DeDx_Data_Minus[1]->SetTitle("#mu for #beta#gamma dE/dx Minus;#eta;#beta#gamma");
      m_eta_betaGammaETOF_DeDx_Data_Minus[2] = (TH2D*) structHisto_eta_betagamma_ETOF->Clone();
      m_eta_betaGammaETOF_DeDx_Data_Minus[2]->SetName("etaBetaGammaDeDx_ETOF_Minus_Sig");
      m_eta_betaGammaETOF_DeDx_Data_Minus[2]->SetTitle("#sigma for #beta#gamma dE/dx Minus;#eta;#beta#gamma");

    }

    m_eta_betaGammaBTOF_DeDx_Data_Plus[0] = (TH2D*) structHisto_eta_betagamma_BTOF->Clone();
    m_eta_betaGammaBTOF_DeDx_Data_Plus[0]->SetName("etaBetaGammaDeDx_BTOF_Plus_N");
    m_eta_betaGammaBTOF_DeDx_Data_Plus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Plus;#eta;#beta#gamma");
    m_eta_betaGammaBTOF_DeDx_Data_Plus[1] = (TH2D*) structHisto_eta_betagamma_BTOF->Clone();
    m_eta_betaGammaBTOF_DeDx_Data_Plus[1]->SetName("etaBetaGammaDeDx_BTOF_Plus_Mu");
    m_eta_betaGammaBTOF_DeDx_Data_Plus[1]->SetTitle("#mu for #beta#gamma dE/dx Plus;#eta;#beta#gamma");
    m_eta_betaGammaBTOF_DeDx_Data_Plus[2] = (TH2D*) structHisto_eta_betagamma_BTOF->Clone();
    m_eta_betaGammaBTOF_DeDx_Data_Plus[2]->SetName("etaBetaGammaDeDx_BTOF_Plus_Sig");
    m_eta_betaGammaBTOF_DeDx_Data_Plus[2]->SetTitle("#sigma for #beta#gamma dE/dx Plus;#eta;#beta#gamma");

    m_eta_betaGammaBTOF_DeDx_Data_Minus[0] = (TH2D*) structHisto_eta_betagamma_BTOF->Clone();
    m_eta_betaGammaBTOF_DeDx_Data_Minus[0]->SetName("etaBetaGammaDeDx_BTOF_Minus_N");
    m_eta_betaGammaBTOF_DeDx_Data_Minus[0]->SetTitle("Number of Counts for #beta#gamma dE/dx Minus;#eta;#beta#gamma");
    m_eta_betaGammaBTOF_DeDx_Data_Minus[1] = (TH2D*) structHisto_eta_betagamma_BTOF->Clone();
    m_eta_betaGammaBTOF_DeDx_Data_Minus[1]->SetName("etaBetaGammaDeDx_BTOF_Minus_Mu");
    m_eta_betaGammaBTOF_DeDx_Data_Minus[1]->SetTitle("#mu for #beta#gamma dE/dx Minus;#eta;#beta#gamma");
    m_eta_betaGammaBTOF_DeDx_Data_Minus[2] = (TH2D*) structHisto_eta_betagamma_BTOF->Clone();
    m_eta_betaGammaBTOF_DeDx_Data_Minus[2]->SetName("etaBetaGammaDeDx_BTOF_Minus_Sig");
    m_eta_betaGammaBTOF_DeDx_Data_Minus[2]->SetTitle("#sigma for #beta#gamma dE/dx Minus;#eta;#beta#gamma");
    cout  << "some structs: " << structHisto_eta_betagamma_ETOF << "  " << structHisto_eta_betagamma_BTOF << endl;

    #ifndef _SAVE_POINTERS_
      if(m_hasETOF) delete structHisto_eta_betagamma_ETOF;
      delete structHisto_eta_betagamma_BTOF;
    #endif

  }

  make_TPC_Histo2Ds(0); // 0 for TPC, 1 for BTOF, 2 for ETOF
  make_TPC_Histo2Ds(1);
  make_TPC_Histo2Ds(2);
}

//=========================================================================================================================================================================
/*
void ZFitter::setupSpectra(int a_partIndex, int a_detectorIndex){
  for(int partBaseIndex = 0; partBaseIndex < m_numParticles; partBaseIndex++){
    if(a_partIndex >= 0 && a_partIndex != partBaseIndex) continue;
    if(!m_particleLoaded[partBaseIndex]) continue;
    #ifdef _ZFITTER_DEBUG_
      cout << "Setting up Spectra for " << m_partInfo->GetParticleName(partBaseIndex,0).Data() << endl;
    #endif
    TH2D* structHisto[2] = {nullptr,nullptr};
    if(a_detectorIndex < 0 || a_detectorIndex == 0){
      structHisto[0] = HistogramUtilities::projection2D(m_rap_mTm0_dEdx_Plus[partBaseIndex][0], 1, 1);
      structHisto[1] = HistogramUtilities::projection2D(m_rap_mTm0_dEdx_Minus[partBaseIndex][0], 1, 1);
      structHisto[0]->SetTitle(Form(";y_{%s};m_{T}-m_{%s} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partBaseIndex).Data(),m_partInfo->GetParticleSymbol(partBaseIndex).Data()));
      structHisto[1]->SetTitle(Form(";y_{%s};m_{T}-m_{%s} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partBaseIndex).Data(),m_partInfo->GetParticleSymbol(partBaseIndex).Data()));
      HistogramUtilities::zeroBins2D(structHisto[0]);
      HistogramUtilities::zeroBins2D(structHisto[1]);
    }

    TH2D* structHistoInvBeta[2] = {nullptr,nullptr};
    if(a_detectorIndex < 0 || a_detectorIndex == 1){
      structHistoInvBeta[0] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaBTOF_Plus[partBaseIndex][0], 1, 1);
      structHistoInvBeta[1] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaBTOF_Minus[partBaseIndex][0], 1, 1);
      structHistoInvBeta[0]->SetTitle(Form(";y_{%s};m_{T}-m_{%s} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partBaseIndex).Data(),m_partInfo->GetParticleSymbol(partBaseIndex).Data()));
      structHistoInvBeta[1]->SetTitle(Form(";y_{%s};m_{T}-m_{%s} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partBaseIndex).Data(),m_partInfo->GetParticleSymbol(partBaseIndex).Data()));
      HistogramUtilities::zeroBins2D(structHistoInvBeta[0]);
      HistogramUtilities::zeroBins2D(structHistoInvBeta[1]);
    }

    TH2D* structHistoInvBetaETOF[2] = {nullptr,nullptr};;
    if(m_hasETOF && (a_detectorIndex < 0 || a_detectorIndex == 2)){
      structHistoInvBetaETOF[0] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaETOF_Plus[partBaseIndex][0], 1, 1);
      structHistoInvBetaETOF[1] = HistogramUtilities::projection2D(m_rap_mTm0_invBetaETOF_Minus[partBaseIndex][0], 1, 1);
      structHistoInvBetaETOF[0]->SetTitle(Form(";y_{%s};m_{T}-m_{%s} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partBaseIndex).Data(),m_partInfo->GetParticleSymbol(partBaseIndex).Data()));
      structHistoInvBetaETOF[1]->SetTitle(Form(";y_{%s};m_{T}-m_{%s} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partBaseIndex).Data(),m_partInfo->GetParticleSymbol(partBaseIndex).Data()));
      HistogramUtilities::zeroBins2D(structHistoInvBetaETOF[0]);
      HistogramUtilities::zeroBins2D(structHistoInvBetaETOF[1]);
    }

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      //###### Plus ######
      // Chi Squared for Plus
      m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][0] = (TH2D*) structHisto[0]->Clone();
      m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][0]->SetName(Form("dEdxChiSqr_Plus_Cent%02d",centIndex));
      m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][0]->SetTitle(Form("dE/dx Chi Squared Plus of Centrality Bin %02d",centIndex));
      m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
      m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][0]->SetName(Form("InvBetaChiSqr_Plus_Cent%02d",centIndex));
      m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][0]->SetTitle(Form("Inverse Beta Chi Squared Plus of Centrality Bin %02d",centIndex));
      if(m_hasETOF){
        m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
        m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][0]->SetName(Form("InvBetaETOFChiSqr_Plus_Cent%02d",centIndex));
        m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][0]->SetTitle(Form("Inverse Beta ETOF Chi Squared Plus of Centrality Bin %02d",centIndex));
      }

      // Fit Status for Plus
      m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][0] = (TH2D*) structHisto[0]->Clone();
      m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][0]->SetName(Form("dEdxFitStatus_Plus_Cent%02d",centIndex));
      m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][0]->SetTitle(Form("dE/dx Fit Status Plus of Centrality Bin %02d",centIndex));
      m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
      m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][0]->SetName(Form("InvBetaBTOFFitStatus_Plus_Cent%02d",centIndex));
      m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][0]->SetTitle(Form("Inverse Beta Fit Status  Plus of Centrality Bin %02d",centIndex));
      if(m_hasETOF){
        m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
        m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][0]->SetName(Form("InvBetaETOFBTOFFitStatus_Plus_Cent%02d",centIndex));
        m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][0]->SetTitle(Form("Inverse Beta ETOF Fit Status  Plus of Centrality Bin %02d",centIndex));
      }

      // Spectra for Plus
      m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][0] = (TH2D*) structHisto[0]->Clone();
      m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][0]->SetName(Form("dEdxSpectra_%s_Plus_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][0]->SetTitle(Form("dE/dx Amp for %s Plus of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][0] = (TH2D*) structHistoInvBeta[0]->Clone();
      m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][0]->SetName(Form("InvBetaBTOFSpectra_%s_Plus_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][0]->SetTitle(Form("Inverse Beta BTOF Amp for %s Plus of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      if(m_hasETOF){
        m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][0] = (TH2D*) structHistoInvBetaETOF[0]->Clone();
        m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][0]->SetName(Form("InvBetaETOFSpectra_%s_Plus_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
        m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][0]->SetTitle(Form("Inverse Beta ETOF Amp for %s Plus of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      }


      //###### MINUS ######

      // Chi Squared for Minus
      m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][1] = (TH2D*) structHisto[1]->Clone();
      m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][1]->SetName(Form("dEdxChiSqr_Minus_Cent%02d",centIndex));
      m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][1]->SetTitle(Form("dE/dx Chi Squared Minus of Centrality Bin %02d",centIndex));
      m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
      m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][1]->SetName(Form("InvBetaChiSqr_Minus_Cent%02d",centIndex));
      m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][1]->SetTitle(Form("Inverse Beta Chi Squared Minus of Centrality Bin %02d",centIndex));
      if(m_hasETOF){
        m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
        m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][1]->SetName(Form("InvBetaETOFChiSqr_Minus_Cent%02d",centIndex));
        m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][1]->SetTitle(Form("Inverse Beta ETOF Chi Squared Minus of Centrality Bin %02d",centIndex));
      }

      // Fit Status for Minus
      m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][1] = (TH2D*) structHisto[1]->Clone();
      m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][1]->SetName(Form("dEdxFitStatus_Minus_Cent%02d",centIndex));
      m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][1]->SetTitle(Form("dE/dx Fit Status Minus of Centrality Bin %02d",centIndex));
      m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
      m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][1]->SetName(Form("InvBetaBTOFFitStatus_Minus_Cent%02d",centIndex));
      m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][1]->SetTitle(Form("Inverse Beta Fit Status  Minus of Centrality Bin %02d",centIndex));
      if(m_hasETOF){
        m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
        m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][1]->SetName(Form("InvBetaETOFBTOFFitStatus_Minus_Cent%02d",centIndex));
        m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][1]->SetTitle(Form("Inverse Beta ETOF Fit Status  Minus of Centrality Bin %02d",centIndex));
      }

      //Spectra for Minus
      m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][1] = (TH2D*) structHisto[1]->Clone();
      m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][1]->SetName(Form("dEdxSpectra_%s_Minus_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][1]->SetTitle(Form("dE/dx Amp for %s Minus of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][1] = (TH2D*) structHistoInvBeta[1]->Clone();
      m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][1]->SetName(Form("InvBetaBTOFSpectra_%s_Minus_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][1]->SetTitle(Form("Inverse Beta BTOF Amp for %s Minus of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      if(m_hasETOF){
        m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][1] = (TH2D*) structHistoInvBetaETOF[1]->Clone();
        m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][1]->SetName(Form("InvBetaETOFSpectra_%s_Minus_Cent%02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
        m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][1]->SetTitle(Form("Inverse Beta ETOF Amp for %s Minus of Centrality Bin %02d",m_partInfo->GetParticleName(partBaseIndex,0).Data(),centIndex));
      }

    }
  }//part loop
  cout << "Finished Setup Histograms" << endl;
}
*/


//=========================================================================================================================================================================

bool ZFitter::areDeDxParamsRidiculous(){
  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
     if(iii%3 == 0 && m_currentParameters[iii] < 0.0) return true;
     if(iii%3 == 2 && m_currentParameters[iii] <= 0.0) return true;
  }
  return false;
}

//=========================================================================================================================================================================

bool ZFitter::areBTOFParamsRidiculous(){
  int divNum = 3;
  if(m_useStudentTDistributionsForTOF){
    divNum = 4;
  }
  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
     if(iii%divNum == 0 && m_currentParameters[iii] < 0.0) return true;
     if(iii%divNum == 2  && m_currentParameters[iii] <= 0.0) return true;
     if(iii%divNum == 3  && m_currentParameters[iii] <= 0.0) return true;
  }
  return false;
}

//=========================================================================================================================================================================

bool ZFitter::doThisTPCBinByColliderStopTable(int a_partIndex,double a_rapidity,double a_momentum){
  int index = -1;
  for(unsigned int iii = 0; iii < m_colliderStopTable.size(); iii++){
    if(m_colliderStopTable[iii][0] == a_partIndex && fabs(a_rapidity - m_colliderStopTable[iii][0]) < 0.01 ){
      index = iii;
      break;
    }
  }
  if(a_momentum < m_colliderStopTable[index][2]) return true;
  double dist = m_colliderStopTable[index][3] - m_colliderStopTable[index][2];
  if(dist > 0.1){
    if(a_momentum < m_colliderStopTable[index][2] + 0.5*dist) return true;
    else return false;
  }else{
    if(a_momentum < m_colliderStopTable[index][4]) return true;
    else return false;
  }

  return false;
}


//=========================================================================================================================================================================


double ZFitter::gausAmp(const double* a_params, double a_y){
  return a_params[0] * 0.3989422804 * TMath::Exp(
   -0.5*( a_y - a_params[1] )*( a_y - a_params[1] )/(a_params[2]*a_params[2])
   ) /(a_params[2]);
}

double ZFitter::gaus(double a_mean, double a_sig, double a_x){
  return 0.3989422804 * TMath::Exp(-0.5*( a_x - a_mean )*( a_x - a_mean )/(a_sig*a_sig)) /a_sig;
}

//=========================================================================================================================================================================

double ZFitter::logNorm(double a_mu, double a_sig, double a_x){
  return 0.3989422804 * TMath::Exp(
   -0.5*( TMath::Log(a_x) - a_mu )*( TMath::Log(a_x) - a_mu )/(a_sig*a_sig)
   ) /(a_sig * a_x);
}

double ZFitter::logNormAmp(double* a_params, double a_x){
  return a_params[0]*0.3989422804 * TMath::Exp(
   -0.5*( TMath::Log(a_x) - a_params[1] )*( TMath::Log(a_x) - a_params[1] )/(a_params[2]*a_params[2])
   ) /(a_params[2] * a_x);
}

//=========================================================================================================================================================================

double ZFitter::studentT(double a_mean, double a_sig, double a_ndf, double a_x){
  return TMath::Student((a_x - a_mean)/a_sig, a_ndf)/a_sig;
}

double ZFitter::studentTAmp(double* a_params, double a_x){
  return a_params[0]*studentT(a_params[1],a_params[2],a_params[3], a_x);

}
//============================================================================================================
//double ZFitter::predictZTPC(int m_currentPartIndexOfInterest, int m_currentPartIndexZTPCSpace, double a_mom){
//  //double momentum = PhysMath::pTot_from_y_mTm0(a_rap,a_mTm0, m_partInfo->GetParticleMass(m_currentPartIndexZTPCSpace));
//  return log( m_bichselCurves[m_currentPartIndexOfInterest]->Eval(a_mom) / m_bichselCurves[m_currentPartIndexZTPCSpace]->Eval(a_mom));
//}


//double ZFitter::predictZTOF(int m_currentPartIndexOfInterest, int m_currentPartIndexZTPCSpace, double a_mom){
//  //double a_mom = PhysMath::pTot_from_y_mTm0(a_rap,a_mTm0, m_partInfo->GetParticleMass(m_currentPartIndexZTPCSpace));
//  return log( m_bichselCurves[m_currentPartIndexOfInterest]->Eval(a_mom) / m_bichselCurves[m_currentPartIndexZTPCSpace]->Eval(a_mom));
//}




//=========================================================================================================================================================================


//uses internal state parameters to return ampl for given dEdx
double ZFitter::ZTPCFunct(double a_ZTPC){
  double val = 0;
  for(int iii = 0; iii < m_currentNumPartToFit; iii++){
    val+=m_currentIntegralofHistogram * gausAmp(&m_currentParameters[3*iii], a_ZTPC);
    //val+=logNormAmp(&m_currentParameters[3*iii], a_ZTPC);
  }
  if(std::isnan(val)){
   #ifdef _ZFITTER_DEBUG_
   cout << "Eval ZTPCFunct at ZTPC=" << a_ZTPC << " Is Nan" << endl;
   printCurrentParams();
   #endif

   return 0.0;
 }
 return val;
}


//uses internal state parameters to return ampl for given dEdx
double ZFitter::ZTPCFunct_student(double a_ZTPC){
  double val = 0;
  for(int iii = 0; iii < m_currentNumPartToFit; iii++){
    val += m_currentIntegralofHistogram * studentTAmp(&m_currentParameters[4*iii], a_ZTPC);
  }
  if(std::isnan(val)){
    #ifdef _ZFITTER_DEBUG_
     cout << "Z_TPC = " << a_ZTPC << " Is Nan  with total integral as " << m_currentIntegralofHistogram <<  endl;
     cout << "     with   ";
     printCurrentParams();
    #endif
   return 0.0;

 }
 return val;
}


//=========================================================================================================================================================================

//uses internal state parameters to return ampl for given dEdx
double ZFitter::ZTOFFunct(double a_ZTOF){
  double val = 0;
  for(int iii = 0; iii < m_currentNumPartToFit; iii++){
    val += m_currentIntegralofHistogram * gausAmp(&m_currentParameters[3*iii], a_ZTOF);
  }
  if(std::isnan(val)){
    #ifdef _ZFITTER_DEBUG_
      cout << "invBeta = " << a_ZTOF << " Is Nan" << endl;
      cout << "     with   ";
      printCurrentParams();
   #endif
   return 0.0;
 }
 return val;
}

//=========================================================================================================================================================================

//uses internal state parameters to return ampl for given dEdx
double ZFitter::ZTOFFunct_student(double a_ZTOF){
  double val = 0;
  for(int iii = 0; iii < m_currentNumPartToFit; iii++){
    val += m_currentIntegralofHistogram * studentTAmp(&m_currentParameters[4*iii], a_ZTOF);
  }
  if(std::isnan(val)){
    #ifdef _ZFITTER_DEBUG_
     cout << "invBeta = " << a_ZTOF << " Is Nan  with total integral as " << m_currentIntegralofHistogram <<  endl;
     cout << "     with   ";
     printCurrentParams();
    #endif
   return 0.0;

 }
 return val;
}


//=========================================================================================================================================================================


double ZFitter::ZTPCOperator(double* a_x,const double* a_params){



  int paramIndexCounter = 0;
  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
    if(m_paramToggles[iii]){
      m_currentParameters[iii] = a_params[paramIndexCounter];
      //cout << "setting param: " << iii << " to " << a_params[paramIndexCounter] << endl;
      paramIndexCounter++;
      //if(iii%3==0 && m_currentParameters[iii] < 0.0){
        //m_currentParameters[iii] = 0.0;
        //cout << "Negative Amplitude, setting it to zero" << endl;
      //}

    }
  }

  //if the particle has a tof constraint, change the m_current parameters
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_currentParticleTags[partIndex] == 4){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      m_currentParameters[partLoc] = m_currentParameters[partMimicLoc] * m_currentConstrainedAmpRaios[partIndex];
      // mean = ZTPC + MimicShiftFromPrediction
      m_currentParameters[partLoc+1] = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin))
      + m_currentParameters[partMimicLoc+1]  -  log(m_bichselCurves[m_particleMimicIndicies[partIndex]]->Eval(m_currentMomOfBin));
                   //+ log(m_currentParameters[partMimicLoc+1] / m_bichselCurves[m_particleMimicIndicies[partIndex]]->Eval(m_currentMomOfBin));
      m_currentParameters[partLoc+2] = m_currentParameters[partMimicLoc+2];
      // stdDev = stdDevMimic * dEdx / dEdxMimic
      //m_currentParameters[partLoc+2] = m_currentParameters[partMimicLoc+2] * (m_currentParameters[partLoc+1] + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)))
      //                    / (m_currentParameters[partMimicLoc+1] + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)));
      //cout << "Tag4: " <<  m_currentParameters[partLoc] << " " <<  m_currentParameters[partLoc+1] << " " <<  m_currentParameters[partLoc+2] << endl;
    }else if(m_currentParticleTags[partIndex] == 5){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      // mean = ZTPC + MimicShiftFromPrediction
      m_currentParameters[partLoc+1] = m_currentMeansForTagFive[partIndex] + m_currentParameters[partMimicLoc+1]
      -  log(m_bichselCurves[m_particleMimicIndicies[partIndex]]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
    }else if(m_currentParticleTags[partIndex] == 6){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);

      // mean = ZTPC + MimicShiftFromPrediction
      m_currentParameters[partLoc+1] = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin))
      + m_currentParameters[partMimicLoc+1]  -  log(m_bichselCurves[m_particleMimicIndicies[partIndex]]->Eval(m_currentMomOfBin));
                   //+ log(m_currentParameters[partMimicLoc+1] / m_bichselCurves[m_particleMimicIndicies[partIndex]]->Eval(m_currentMomOfBin));
      m_currentParameters[partLoc+2] = m_currentParameters[partMimicLoc+2];
      // stdDev = stdDevMimic * dEdx / dEdxMimic
      //m_currentParameters[partLoc+2] = m_currentParameters[partMimicLoc+2] * (m_currentParameters[partLoc+1] + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)))
      //                    / (m_currentParameters[partMimicLoc+1] + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)));
    }else if(m_currentParticleTags[partIndex] == 7){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      m_currentParameters[partLoc] = m_currentParameters[partMimicLoc] * m_currentConstrainedAmpRaios[partIndex];
    }else if(m_currentParticleTags[partIndex] == 8){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      m_currentParameters[partLoc] = m_currentParameters[partMimicLoc] * m_currentConstrainedAmpRaios[partIndex];
    }else if(m_currentParticleTags[partIndex] == 10){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      m_currentParameters[partLoc+2] = m_currentParameters[partMimicLoc+2];
    }else if(m_currentParticleTags[partIndex] == 11){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      m_currentParameters[partLoc+2] = m_currentParameters[partMimicLoc+2];
    }else if(m_currentParticleTags[partIndex] == 12){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      m_currentParameters[partLoc] = m_currentParameters[partMimicLoc] * m_currentConstrainedAmpRaios[partIndex];
      m_currentParameters[partLoc+2] = m_currentParameters[partMimicLoc+2];
    }else if(m_currentParticleTags[partIndex] == 13){
      int partLoc = getCurrParamIndex(partIndex);
      int partMimicLoc = getCurrParamIndex(m_particleMimicIndicies[partIndex]);
      m_currentParameters[partLoc] = m_currentParameters[partMimicLoc] * m_currentConstrainedAmpRaios[partIndex];
    }//end particle tag if statements
  }//end partIndex loop


  //cout << "with dEdxFunt = " << dEdxFunct(a_x[0]) << endl;
  return ZTPCFunct(a_x[0]);
}

//=========================================================================================================================================================================

double ZFitter::ZTOFOperator(double* a_x,const double* a_params){
  int paramIndexCounter = 0;
  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
    if(m_paramToggles[iii]){
      m_currentParameters[iii] = a_params[paramIndexCounter];
      //cout << "setting param: " << iii << " to " << a_params[paramIndexCounter] << endl;
      paramIndexCounter++;
    }
  }
  return ZTOFFunct(a_x[0]);
}

//=========================================================================================================================================================================

double ZFitter::ZTOFOperator_student(double* a_x,const double* a_params){
  int paramIndexCounter = 0;
  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
    if(m_paramToggles[iii]){
      m_currentParameters[iii] = a_params[paramIndexCounter];
      //cout << "setting param: " << iii << " to " << a_params[paramIndexCounter] << endl;
      paramIndexCounter++;
    }
  }
  return ZTOFFunct_student(a_x[0]);
}
//=========================================================================================================================================================================
//=========================================================================================================================================================================

#ifdef _MAC_OSX_
double ZFitter::chiSquared_TPC(const double* a_params){
  double chiSqrd = 0.0;
  for(int binX = 1; binX <= m_currentHistoToFit->GetNbinsX(); binX++){
    double binContent = m_currentHistoToFit->GetBinContent(binX);
    if(binContent <= 0.0) continue;
    double binError = m_currentHistoToFit->GetBinError(binX);
    double binCenter = m_currentHistoToFit->GetBinCenter(binX);
    double prediction = ZTPCOperator(&binCenter,a_params);
    chiSqrd += (binContent - prediction)*(binContent - prediction)/(binError*binError);
  }
  return chiSqrd;
}


double ZFitter::chiSquared_TOF(const double* a_params){
  double chiSqrd = 0.0;
  for(int binX = 1; binX <= m_currentHistoToFit->GetNbinsX(); binX++){
    double binContent = m_currentHistoToFit->GetBinContent(binX);
    if(binContent <= 0.0) continue;
    double binError = m_currentHistoToFit->GetBinError(binX);
    double binCenter = m_currentHistoToFit->GetBinCenter(binX);
    double prediction = ZTOFOperator_student(&binCenter,a_params);
    chiSqrd += (binContent - prediction)*(binContent - prediction)/(binError*binError);
  }
  return chiSqrd;
}

#endif

//=========================================================================================================================================================================
//=========================================================================================================================================================================







double ZFitter::dEdxBetaGammaOperator(double *x,  double *par){

  double betaGamma = x[0];
  double m1 = par[0];
  double m2 = par[1];
  double m3 = par[2];
  double m4 = par[3];
  double m5 = par[4];
  double m6 = par[5];
  double m7 = par[6];

  double term1 = TMath::Power(1.0 + TMath::Power(1.0/betaGamma,2),m2);
  double term2 = TMath::Power(TMath::Abs(TMath::Log(m3*TMath::Power(betaGamma,2))),m4);
  term2 += m5 * TMath::Power(1.0+TMath::Power(1.0/betaGamma,2),m6);

  return m1 * term1 * term2 - m7;

}

//=========================================================================================================================================================================
//=========================================================================================================================================================================

double ZFitter::ZTPC_prediction_betagamma_fit(int a_partIndex,double a_y, double a_mTm0){
  //double eta =
  return 1.0;

}

//=========================================================================================================================================================================
//=========================================================================================================================================================================
vector<TF1*> ZFitter::drawIndividualFits_Gaus(){
  int graphIndex = 0;
  vector<TF1*> returnVect;
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_currentParticleTags[partIndex] != 0){
     TF1* indFunct = new TF1(Form("indFunct%02d",partIndex),
      "0.3989422804 * [0] * exp(  -0.5*( x - [1] )*( x - [1] )/([2]*[2])   ) /[2]",
      m_currentParameters[graphIndex*3+1] - 6.0*m_currentParameters[graphIndex*3+2],
      m_currentParameters[graphIndex*3+1] + 6.0*m_currentParameters[graphIndex*3+2]);
     indFunct->SetParameter(0,m_currentParameters[graphIndex*3] * m_currentIntegralofHistogram);
     indFunct->SetParameter(1,m_currentParameters[graphIndex*3+1]);
     indFunct->SetParameter(2,m_currentParameters[graphIndex*3+2]);

     indFunct->SetNpx(300);
     indFunct->SetLineColor(m_partInfo->GetParticleColor(partIndex));
     indFunct->Draw("same");
     graphIndex++;
     returnVect.push_back(indFunct);
   }
 }
 return returnVect;
}
//=========================================================================================================================================================================
//=========================================================================================================================================================================

vector<TF1*> ZFitter::drawIndividualFits_Student(){
  int graphIndex = 0;
  vector<TF1*> returnVect;
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_currentParticleTags[partIndex] != 0){
     TF1* indFunct = (TF1*) m_student_T_Funct->Clone();
     indFunct->SetName(Form("indFunct%02d",partIndex));
     indFunct->SetRange(-1.0,10.0);
     indFunct->SetParameter(0,m_currentParameters[graphIndex*4]* m_currentIntegralofHistogram);
     indFunct->SetParameter(1,m_currentParameters[graphIndex*4+1]);
     indFunct->SetParameter(2,m_currentParameters[graphIndex*4+2]);
     indFunct->SetParameter(3,m_currentParameters[graphIndex*4+3]);
     indFunct->SetNpx(300);
     indFunct->SetLineColor(m_partInfo->GetParticleColor(partIndex));
     indFunct->Draw("same");
     graphIndex++;
     returnVect.push_back(indFunct);
   }
 }
 return returnVect;
}

//=========================================================================================================================================================================
//=========================================================================================================================================================================

vector<TF1*> ZFitter::drawIndividualFits_LogNorm(){
  int graphIndex = 0;
  vector<TF1*> retVect;
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_currentParticleTags[partIndex] != 0){

     TF1* indFunct = new TF1("indFunct",
      "[0]*0.3989422804 * TMath::Exp(-0.5*( TMath::Log(x) - [1] )*( TMath::Log(x) - [1] )/([2]*[2]) ) /([2] * x)"
      ,0,200);
     indFunct->SetParameter(0,m_currentParameters[graphIndex*3] * m_currentIntegralofHistogram);
     indFunct->SetParameter(1,m_currentParameters[graphIndex*3+1]);
     indFunct->SetParameter(2,m_currentParameters[graphIndex*3+2]);

     indFunct->SetNpx(500);
     //     if(m_currentParticleTags[partIndex] == 1) indFunct->SetLineColor(kGreen);
     //     if(m_currentParticleTags[partIndex] == 2) indFunct->SetLineColor(kBlue);
     indFunct->SetLineColor(m_partInfo->GetParticleColor(partIndex));
     indFunct->Draw("same");
     graphIndex++;
     retVect.push_back(indFunct);
   }
 }
 return retVect;
}

//=========================================================================================================================================================================
//=========================================================================================================================================================================
void ZFitter::printCurrentParams(){
  cout << "  Current Parameters: ";
  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
    cout << m_currentParameters[iii] << "  ";
  }
  cout << endl;
}

//=========================================================================================================================================================================
//=========================================================================================================================================================================

int ZFitter::getCurrParamIndex(int a_partIndex){
  int currentCount = 0;
  for(int iii = 0; iii < m_numParticles; iii++){
    if(m_currentParticleTags[iii] != 0){
      if(iii == a_partIndex) return currentCount;
      currentCount += 3;
    }
  }

  cout << "Particle " << a_partIndex << "not Fit ###!@#3##" << endl;
  return -1;

}


//=========================================================================================================================================================================
//=========================================================================================================================================================================


double ZFitter::predictZTPC_fromEtaBetaGammaFits(int a_partIndexOfInterest, int a_partIndexZTPCSpace, int a_charge, double a_rapidity, double a_mTm0){
  double eta = PhysMath::eta_from_y_mTm0(a_rapidity, a_mTm0,m_partInfo->GetParticleMass(a_partIndexZTPCSpace));
  double momentum = PhysMath::pTot_from_y_mTm0(a_rapidity,a_mTm0, m_partInfo->GetParticleMass(a_partIndexZTPCSpace));
  //double betaGamma_ZTPCSpace = PhysMath::betagamma(momentum, m_partInfo->GetParticleMass(a_partIndexZTPCSpace));
  double betaGamma_OfInterest = PhysMath::betagamma(momentum, m_partInfo->GetParticleMass(a_partIndexOfInterest));
  cout << " eta: " << eta << " mom: " << momentum << " betaGamma_OfInterest: " << betaGamma_OfInterest << endl;
  TF1* fitToUse = NULL;
  if(a_charge > 0.0){
    int etaBin = m_eta_betaGammaBtof_dEdx_Plus->GetXaxis()->FindBin(eta);
    if(etaBin > 0 && etaBin <= m_eta_betaGammaBtof_dEdx_Plus->GetNbinsX()){
      fitToUse = m_eta_betaGammaBTOF_dEdx_FitFunction[0][etaBin-1];
    }
    cout  << " initial try: " << etaBin << " with " << fitToUse << endl;
    if(!fitToUse){
      //find the closest one
      int closestIndex = 300;
      int closestDiff = 300;
      for(int index = 0; index < m_eta_betaGammaBTOF_dEdx_FitFunction[0].size(); index++){
        if(m_eta_betaGammaBTOF_dEdx_FitFunction[0][index] && abs(etaBin - index - 1) < closestDiff){
          closestIndex = index;
          closestDiff = abs(etaBin - index - 1);
        }
      }
      cout  << " closest: " << closestIndex << " with " << closestDiff << endl;

      if(closestIndex != 300) fitToUse = m_eta_betaGammaBTOF_dEdx_FitFunction[0][closestIndex];
    }

  }else{
    int etaBin = m_eta_betaGammaBtof_dEdx_Minus->GetXaxis()->FindBin(eta);
    if(etaBin > 0 && etaBin <= m_eta_betaGammaBtof_dEdx_Plus->GetNbinsX()){
      fitToUse = m_eta_betaGammaBTOF_dEdx_FitFunction[1][etaBin-1];
    }
    if(!fitToUse){
      //find the closest one
      int closestIndex = 300;
      int closestDiff = 300;
      for(int index = 0; index < m_eta_betaGammaBTOF_dEdx_FitFunction[1].size(); index++){
        if(m_eta_betaGammaBTOF_dEdx_FitFunction[1][index] && abs(etaBin - index - 1) < closestDiff){
          closestIndex = index;
          closestDiff = abs(etaBin - index - 1);
        }
      }
      if(closestIndex != 300) fitToUse = m_eta_betaGammaBTOF_dEdx_FitFunction[1][closestIndex];
    }
  }
  cout << " fitFunct: " << fitToUse << endl;
  if(fitToUse){
    double initialPrediction = fitToUse->Eval(betaGamma_OfInterest) - log(m_bichselCurves[a_partIndexZTPCSpace]->Eval(momentum));
    if(a_partIndexOfInterest == 0) initialPrediction -= fabs(eta)*0.07;
    if(a_partIndexOfInterest == 1) initialPrediction += fabs(eta)*0.07;
    //dont do anything for proton
    return initialPrediction;
  }else{
    cout << "WARNING: the eta betagamma dedx prediction failed for eta: " << eta << " mom: " << momentum << endl;
    return log(m_bichselCurves[a_partIndexOfInterest]->Eval(momentum)) - log(m_bichselCurves[a_partIndexZTPCSpace]->Eval(momentum));
  }
}




//=========================================================================================================================================================================
//=========================================================================================================================================================================

void ZFitter::drawPartdEdxFits(int a_partIndex){
  if(!m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][0][0]) return;
  TCanvas* resultCanvas = new TCanvas();
  resultCanvas->Divide(2,2);

  resultCanvas->cd(1);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][0][0]->SetContour(99);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][0][0]->Draw("colz");
  resultCanvas->cd(2);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][1][0]->SetContour(99);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][1][0]->Draw("colz");
  resultCanvas->cd(3);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][2][0]->SetContour(99);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][2][0]->Draw("colz");
  resultCanvas->cd(4);
  m_ChiSqr_ZTPC[m_currentPartIndex][0]->SetContour(99);
  m_ChiSqr_ZTPC[m_currentPartIndex][0]->Draw("colz");



  if(!m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][0][1]) return;
  TCanvas* resultCanvas_2 = new TCanvas();
  resultCanvas_2->Divide(2,2);

  resultCanvas_2->cd(1);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][0][1]->SetContour(99);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][0][1]->Draw("colz");
  resultCanvas_2->cd(2);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][1][1]->SetContour(99);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][1][1]->Draw("colz");
  resultCanvas_2->cd(3);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][2][1]->SetContour(99);
  m_Fit_Data_ZTPC[m_currentPartIndex][a_partIndex][2][1]->Draw("colz");
  resultCanvas_2->cd(4);
  m_ChiSqr_ZTPC[m_currentPartIndex][1]->SetContour(99);
  m_ChiSqr_ZTPC[m_currentPartIndex][1]->Draw("colz");



}

//=========================================================================================================================================================================


void ZFitter::draw_ZVar_mTm0_data_with_fit_data(int a_partIndex, int a_detIndex, int a_centIndex, int a_rap_index, int a_charge){

  int pm_index = 0;
  if(a_charge == -1) pm_index = 1;

  string det_labels[3] = {"Z_{TPC}", "Z_{BTOF}", "Z_{ETOF}"};

  TCanvas* canv = new TCanvas();
  canv->SetWindowSize(1600,800);
  canv->SetCanvasSize(1600,800);

  canv->Divide(2,1);

  //Get the Histo2D data to draw first
  Histo2D* dataHist = nullptr;
  if(a_detIndex == 0) dataHist = m_mTm0_ZTPC_by_Rapidity[a_partIndex][a_centIndex][a_rap_index][pm_index];
  if(a_detIndex == 1) dataHist = m_mTm0_ZBTOF_by_Rapidity[a_partIndex][a_centIndex][a_rap_index][pm_index];
  if(a_detIndex == 2) dataHist = m_mTm0_ZETOF_by_Rapidity[a_partIndex][a_centIndex][a_rap_index][pm_index];
  if(!dataHist) return;

  //Generate TGraphErrors for the fit data of each particle with the mean as the y value and the stdDev as the error
  TGraphErrors* fitGraphs[9] = {nullptr};
  TGraphErrors* ampGraphs[9] = {nullptr};
  TGraphErrors* ratioGraphs[9] = {nullptr};
  double highestXValue = -999;
  double lowestXValue = 999;
  double highestYValue = -999;
  double lowestYValue = 999;
  bool allEmpty = true;
  for(int partIndex = 0; partIndex < 9; partIndex++){
    TH2D* fitDataAmp = nullptr;
    TH2D* fitDataMean = nullptr;
    TH2D* fitDataStdDev = nullptr;
    if(a_detIndex == 0){
      fitDataAmp    = m_Fit_Data_Cent_ZTPC[a_partIndex][partIndex][a_centIndex][0][pm_index];
      fitDataMean   = m_Fit_Data_Cent_ZTPC[a_partIndex][partIndex][a_centIndex][1][pm_index];
      fitDataStdDev = m_Fit_Data_Cent_ZTPC[a_partIndex][partIndex][a_centIndex][2][pm_index];
    }
    if(a_detIndex == 1){
      fitDataAmp    = m_Fit_Data_Cent_ZbTOF[a_partIndex][partIndex][a_centIndex][0][pm_index];
      fitDataMean   = m_Fit_Data_Cent_ZbTOF[a_partIndex][partIndex][a_centIndex][1][pm_index];
      fitDataStdDev = m_Fit_Data_Cent_ZbTOF[a_partIndex][partIndex][a_centIndex][2][pm_index];
    }
    if(a_detIndex == 2){
      fitDataAmp    = m_Fit_Data_Cent_ZbTOF[a_partIndex][partIndex][a_centIndex][0][pm_index];
      fitDataMean   = m_Fit_Data_Cent_ZeTOF[a_partIndex][partIndex][a_centIndex][1][pm_index];
      fitDataStdDev = m_Fit_Data_Cent_ZeTOF[a_partIndex][partIndex][a_centIndex][2][pm_index];
    }

    if(!fitDataMean || !fitDataStdDev) continue;
    fitGraphs[partIndex] = new TGraphErrors();
    ampGraphs[partIndex] = new TGraphErrors();
    for(int mTm0_binIndex = 0; mTm0_binIndex < fitDataMean->GetNbinsY(); mTm0_binIndex++){
      double amp     = fitDataAmp->GetBinContent(a_rap_index+1,mTm0_binIndex+1);
      double amp_err = fitDataAmp->GetBinError(a_rap_index+1,mTm0_binIndex+1);
      double mean    = fitDataMean->GetBinContent(a_rap_index+1,mTm0_binIndex+1);
      double stdDev  = fitDataStdDev->GetBinContent(a_rap_index+1,mTm0_binIndex+1);
      if(stdDev <= 0.0 || (stdDev == 1.0 && mean == 0.0)) continue;
      int numPoints = fitGraphs[partIndex]->GetN();
      fitGraphs[partIndex]->SetPoint(numPoints,fitDataMean->GetYaxis()->GetBinCenter(mTm0_binIndex+1),mean);
      fitGraphs[partIndex]->SetPointError(numPoints,0.0,stdDev);
      ampGraphs[partIndex]->SetPoint(numPoints,fitDataMean->GetYaxis()->GetBinCenter(mTm0_binIndex+1),amp);
      ampGraphs[partIndex]->SetPointError(numPoints,fitDataMean->GetYaxis()->GetBinCenter(mTm0_binIndex+1),amp_err);
      if(fitDataMean->GetYaxis()->GetBinCenter(mTm0_binIndex+1) > highestXValue) highestXValue = fitDataMean->GetYaxis()->GetBinCenter(mTm0_binIndex+1);
      if(fitDataMean->GetYaxis()->GetBinCenter(mTm0_binIndex+1) < lowestXValue) lowestXValue = fitDataMean->GetYaxis()->GetBinCenter(mTm0_binIndex+1);
      if(mean+stdDev > highestYValue) highestYValue = mean+stdDev;
      if(mean-stdDev < lowestYValue) lowestYValue = mean-stdDev;
    }
    fitGraphs[partIndex]->SetMarkerStyle(20);
    fitGraphs[partIndex]->SetMarkerSize(0.25);
    fitGraphs[partIndex]->SetMarkerColor(m_partInfo->GetParticleColor(partIndex));
    fitGraphs[partIndex]->SetLineColor(m_partInfo->GetParticleColor(partIndex));
    ampGraphs[partIndex]->SetMarkerStyle(20);
    ampGraphs[partIndex]->SetMarkerSize(0.25);
    ampGraphs[partIndex]->SetMarkerColor(m_partInfo->GetParticleColor(partIndex));
    ampGraphs[partIndex]->SetLineColor(m_partInfo->GetParticleColor(partIndex));
    if(fitGraphs[partIndex]->GetN() > 0) allEmpty = false;
  }

  if(allEmpty){
    for(int partIndex = 0; partIndex < 9; partIndex++){
      if(fitGraphs[partIndex]) delete fitGraphs[partIndex];
    }
    return;
  }

  //construct the ratios in comparison to the particle of interest
  double highestRatioValue = -999999;
  double lowestRatioValue = 999999;
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(a_partIndex == partIndex) continue;
    if(!ampGraphs[a_partIndex] || !ampGraphs[partIndex]) continue;
    ratioGraphs[partIndex] = new TGraphErrors();
    for(int pointIndex = 0; pointIndex < ampGraphs[partIndex]->GetN(); pointIndex++){
      double x_val = 0.0;
      double y_val = 0.0;
      double y_err = ampGraphs[partIndex]->GetErrorY(pointIndex);
      ampGraphs[partIndex]->GetPoint(pointIndex,x_val,y_val);

      int mainPointIndex = HistogramUtilities::getPointClosestToXVal(ampGraphs[a_partIndex], x_val);
      // SDCC/ROOT5 BUILD FIX (2026-07-02, not in the original): TGraphErrors::GetPointY()
      // doesn't exist in ROOT 5.34.38 (added upstream well after this repo's original Mac
      // build's ROOT 6) -- replaced with the classic GetPoint(index,x,y) form already used
      // a few lines above in this same function, which has always been part of TGraph's
      // API. Same value, just via the older accessor (the x output here is unused).
      double main_x_val_unused = 0.0;
      double main_y_val = 0.0;
      ampGraphs[a_partIndex]->GetPoint(mainPointIndex, main_x_val_unused, main_y_val);
      double main_y_err = ampGraphs[a_partIndex]->GetErrorY(mainPointIndex);

      double ratio = y_val/main_y_val;
      double ratio_err = ratio*sqrt(pow(y_err/y_val,2)+pow(main_y_err/main_y_val,2));
      if(!std::isfinite(ratio)){
        ratio = 1.0;
        ratio_err = 10.0;
      }
      if(!std::isfinite(ratio_err)){
        ratio_err = 10.0;
      }
      ratioGraphs[partIndex]->SetPoint(pointIndex,x_val,ratio);
      ratioGraphs[partIndex]->SetPointError(pointIndex,0.0,ratio_err);

      if(ratio > highestRatioValue) highestRatioValue = ratio;
      if(ratio < lowestRatioValue) lowestRatioValue = ratio;
    }
    ratioGraphs[partIndex]->SetMarkerStyle(20);
    ratioGraphs[partIndex]->SetMarkerSize(0.25);
    ratioGraphs[partIndex]->SetMarkerColor(m_partInfo->GetParticleColor(partIndex));
    ratioGraphs[partIndex]->SetLineColor(m_partInfo->GetParticleColor(partIndex));
  }



  double x_width = highestXValue - lowestXValue;
  lowestXValue -= 0.1*x_width;
  highestXValue += 0.1*x_width;
  double y_width = highestYValue - lowestYValue;
  //lowestYValue -= 0.1*y_width;
  //highestYValue += 0.1*y_width;
  lowestYValue = 0.8*lowestYValue;
  highestYValue = 1.2*highestYValue;

  canv->cd(1);
  dataHist->Draw("colz",true,"m_{T}-m_{0} [GeV/c^{2}]",det_labels[a_detIndex]);
  dataHist->SetXDrawRange(lowestXValue,highestXValue);
  dataHist->SetYDrawRange(lowestYValue,highestYValue);
  dataHist->SetZDrawRangeForLogScale();
  gPad->SetLogz();
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(fitGraphs[partIndex]) fitGraphs[partIndex]->Draw("PL");
  }

  canv->cd(2);
  bool firstToDraw = true;
  TH1D* graph_frame = nullptr;
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(!ratioGraphs[partIndex] || ratioGraphs[partIndex]->GetN() == 0) continue;
    if(firstToDraw){
      graph_frame = HistogramUtilities::generateGraphFrame(ratioGraphs[partIndex],lowestXValue,highestXValue,false,true,"m_{T}-m_{0} [GeV/c^{2}]",Form("Amp. Ratio to %s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()),false);
      firstToDraw = false;
      graph_frame->SetStats(kFALSE);
      graph_frame->GetYaxis()->SetRangeUser(0.8*lowestRatioValue,1.2*highestRatioValue);
      graph_frame->Draw();
    }
    ratioGraphs[partIndex]->Draw("P");
  }
  gPad->SetLogy();


  if(a_detIndex == 0){
    canv->SaveAs(Form("%s/%s/%s/dEdxFits_Overview/Cent%02d_RapBin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex, a_rap_index+1));
  }
  if(a_detIndex == 1){
    canv->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Overview/Cent%02d_RapBin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex, a_rap_index+1));
  }
  if(a_detIndex == 2){
    canv->SaveAs(Form("%s/%s/%s/invBetaETOFFits_Overview/Cent%02d_RapBin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex, a_rap_index+1));
  }
  canv->Clear();
  if(graph_frame) delete graph_frame;
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(ampGraphs[partIndex]) delete ampGraphs[partIndex];
    if(ratioGraphs[partIndex]) delete ratioGraphs[partIndex];
    if(fitGraphs[partIndex]) delete fitGraphs[partIndex];
  }


}



//=========================================================================================================================================================================
//=========================================================================================================================================================================
void ZFitter::getRangeDeDx(double* a_returnData,int a_partIndex){

  double lowDeDx = 0;
  double highDeDx = 0;

  if(a_partIndex == 0){
    lowDeDx = 0.0;
    double massHigh = 0.5*(m_partInfo->GetParticleMass(0) + m_partInfo->GetParticleMass(1));
    TF1* funct1 = m_partInfo->GetBichselFunction(massHigh,1.0);
    highDeDx = funct1->Eval(m_currentMomOfBin);
    #ifndef _SAVE_POINTERS_
      delete funct1;
    #endif
  }else if(a_partIndex == 1){
    double massLow = 0.5*(m_partInfo->GetParticleMass(0) + m_partInfo->GetParticleMass(1));
    double massHigh = 0.5*(m_partInfo->GetParticleMass(1) + m_partInfo->GetParticleMass(2));
    TF1* funct1 = m_partInfo->GetBichselFunction(massLow,1.0);
    TF1* funct2 = m_partInfo->GetBichselFunction(massHigh,1.0);
    lowDeDx = funct1->Eval(m_currentMomOfBin);
    highDeDx = funct2->Eval(m_currentMomOfBin);
    #ifndef _SAVE_POINTERS_
      delete funct1;
      delete funct2;
    #endif
  }else if(a_partIndex == 2){
    double massLow = 0.5*(m_partInfo->GetParticleMass(1) + m_partInfo->GetParticleMass(2));
    double massHigh = 0.5*(m_partInfo->GetParticleMass(2) + m_partInfo->GetParticleMass(4));
    TF1* funct1 = m_partInfo->GetBichselFunction(massLow,1.0);
    TF1* funct2 = m_partInfo->GetBichselFunction(massHigh,1.0);
    lowDeDx = funct1->Eval(m_currentMomOfBin);
    highDeDx = funct2->Eval(m_currentMomOfBin);
    #ifndef _SAVE_POINTERS_
      delete funct1;
      delete funct2;
    #endif
  }else if(a_partIndex > 5){
    cout << "update ZFitter::getRangeDeDx" << endl;
  }else{
    double massLow = 0.5*(m_partInfo->GetParticleMass(a_partIndex - 1) + m_partInfo->GetParticleMass(a_partIndex));
    double massHigh = 0.5*(m_partInfo->GetParticleMass(a_partIndex) + m_partInfo->GetParticleMass(a_partIndex + 1));
    TF1* funct1 = m_partInfo->GetBichselFunction(massLow,1.0);
    TF1* funct2 = m_partInfo->GetBichselFunction(massHigh,1.0);
    lowDeDx = funct1->Eval(m_currentMomOfBin);
    highDeDx = funct2->Eval(m_currentMomOfBin);
    #ifndef _SAVE_POINTERS_
      delete funct1;
      delete funct2;
    #endif
  }

  a_returnData[0] = lowDeDx;
  a_returnData[1] = highDeDx;

  return;
}


//=========================================================================================================================================================================
//=========================================================================================================================================================================

void ZFitter::getRangeDeDx(double &a_lowDeDx, double &a_highDeDx){
  double lowDeDx = 500;
  double highDeDx = -200;

  double allowedSigma = 2.5;
  double sigmaPercentage = 0.07;

  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(m_currentParticleTags[partIndex] == 0) continue;
    double center = m_bichselCurves[partIndex]->Eval(m_currentMomOfBin);
    double sigma = sigmaPercentage * center;
    double lower = center - allowedSigma * sigma;
    double higher = center + allowedSigma * sigma;
    if(lower < lowDeDx) lowDeDx = lower;
    if(higher > highDeDx) highDeDx = higher;
  }

  a_lowDeDx = lowDeDx;
  a_highDeDx = highDeDx;

  double firstLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::firstBinWithContent(m_currentHistoToFit) );
  double lastLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::lastBinWithContent(m_currentHistoToFit) + 1 );

  #ifdef _ZFITTER_DEBUG_
    cout << "##### getRangeDeDx:  " << "  LowCut: " << a_lowDeDx << " HighCut: " << a_highDeDx << " FirstContent: " << firstLocWithContent << " LastContent: " << lastLocWithContent << endl;
  #endif

  if(firstLocWithContent > a_lowDeDx) a_lowDeDx = firstLocWithContent;
  if(lastLocWithContent < a_highDeDx) a_highDeDx = lastLocWithContent;

  return;
}

//=========================================================================================================================================================================
//=========================================================================================================================================================================

void ZFitter::getRangeZTPC(double &a_lowZTPC, double &a_highZTPC){
  double centerForPOI = m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  if(centerForPOI <= 0){
    #ifdef _ZFITTER_DEBUG_
      cout << "getRangeZTPC : Bichsel Curve Doesn't work here: setting range to -0.4 to 0.4" << endl;
    #endif
    a_lowZTPC = -0.4;
    a_highZTPC = 0.4;
    return;
  }
  double lowDeDx = centerForPOI;
  double highDeDx = centerForPOI;

  double sigmaPercentage = 0.07;
  double allowedSigma = 2.5;

  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex == 3 || partIndex == 8) continue;
    /*if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
        || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex])){
      double center = m_bichselCurves[partIndex]->Eval(m_currentMomOfBin);
      if(center <= 0.0) continue;
      double sigma = sigmaPercentage * center;
      double lower = center - allowedSigma * sigma;
      double higher = center + allowedSigma * sigma;
      if(lower < highDeDx) highDeDx = lower;
      if(higher > lowDeDx) lowDeDx = higher;
    }else{ */
    if((m_currentPartCharge == 1 && m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
        || (m_currentPartCharge == -1 && m_partToFitMatrix_Minus[m_currentPartIndex][partIndex])){
      if(m_currentParticleTags[partIndex] == 0) continue;
      double center = m_bichselCurves[partIndex]->Eval(m_currentMomOfBin);
      if(center <= 0.0) continue;
      double sigma = sigmaPercentage * center;
      double lower = center - allowedSigma * sigma;
      double higher = center + allowedSigma * sigma;
      if(lower < lowDeDx) lowDeDx = lower;
      if(higher > highDeDx) highDeDx = higher;
    }


  }

  a_lowZTPC = log(lowDeDx/centerForPOI);
  a_highZTPC = log(highDeDx/centerForPOI);

  if(m_bichselShiftsComplete){
    a_lowZTPC += m_bichselShifts_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
    a_highZTPC += m_bichselShifts_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
  }

  double firstLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::firstBinWithContent(m_currentHistoToFit) );
  double lastLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::lastBinWithContent(m_currentHistoToFit) + 1 );



  if(firstLocWithContent > a_lowZTPC || std::isnan(a_lowZTPC)) a_lowZTPC = firstLocWithContent;
  if(lastLocWithContent < a_highZTPC || std::isnan(a_highZTPC)) a_highZTPC = lastLocWithContent;

  #ifdef _ZFITTER_DEBUG_
    cout << "##### getRangeZTPC: Center dEdx: " << centerForPOI << "  LowCut: " << a_lowZTPC << " HighCut: " << a_highZTPC << " FirstContent: " << firstLocWithContent << " LastContent: " << lastLocWithContent << endl;
  #endif

  return;
}

//=========================================================================================================================================================================
//=========================================================================================================================================================================



void ZFitter::getRangeZTPC_Simple(double &a_lowZTPC, double &a_highZTPC, double a_nSigma){

  a_lowZTPC = -a_nSigma*log(1.07);
  a_highZTPC = -a_lowZTPC;
  return;

}

//=========================================================================================================================================================================
//=========================================================================================================================================================================

void ZFitter::getRangeZTOF(double &a_lowZTOF, double &a_highZTOF, double a_lowMassSquaredOverRide, double a_highMassSquaredOverRide){
  double centerForPOI = PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), m_currentMomOfBin);
  if(m_massSquaredMode) centerForPOI = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0);
  double lowInvBeta = centerForPOI;
  double highInvBeta = centerForPOI;


  double allowedSigma = 2.5;
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex == 3 || partIndex == 8) continue;
    if(m_currentParticleTags[partIndex] == 0) continue;
    if((m_currentPartCharge == 1 && m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
        || (m_currentPartCharge == -1 && m_partToFitMatrix_Minus[m_currentPartIndex][partIndex])){
      double lowerEdge = PhysMath::invBeta(m_cuts->getLowMassCut(partIndex), m_currentMomOfBin);
      double higherEdge = PhysMath::invBeta(m_cuts->getHighMassCut(partIndex), m_currentMomOfBin);
      if(m_massSquaredMode){
        lowerEdge = m_cuts->getLowMassSquaredCut(partIndex);
        higherEdge = m_cuts->getHighMassSquaredCut(partIndex);
      }
      if(a_lowMassSquaredOverRide > -998){
        lowerEdge = PhysMath::invBeta(sqrt(a_lowMassSquaredOverRide), m_currentMomOfBin);
        higherEdge = PhysMath::invBeta(sqrt(a_highMassSquaredOverRide), m_currentMomOfBin);
        if(m_massSquaredMode){
          lowerEdge = a_lowMassSquaredOverRide;
          higherEdge = a_highMassSquaredOverRide;
        }
      }

      if(partIndex == 0 && !m_massSquaredMode) lowerEdge = 0.95;
      if(lowerEdge < lowInvBeta) lowInvBeta = lowerEdge;
      if(higherEdge > highInvBeta) highInvBeta = higherEdge;
    }
  }



  if(m_massSquaredMode){
    a_lowZTOF = lowInvBeta;
    a_highZTOF = highInvBeta;
  }else{
    a_lowZTOF =  lowInvBeta - centerForPOI;
    a_highZTOF = highInvBeta - centerForPOI;
  }

  double firstLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::firstBinWithContent(m_currentHistoToFit) );
  double lastLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::lastBinWithContent(m_currentHistoToFit) + 1 );

  #ifdef _ZFITTER_DEBUG_
    cout << "##### getRangeZTOF: Center InvBeta: " << centerForPOI << "  LowCut: " << a_lowZTOF << " HighCut: " << a_highZTOF << " FirstContent: " << firstLocWithContent << " LastContent: " << lastLocWithContent << endl;
  #endif


  if(firstLocWithContent > a_lowZTOF) a_lowZTOF = firstLocWithContent;
  if(lastLocWithContent < a_highZTOF) a_highZTOF = lastLocWithContent;
  return;
}

//============================================================================================================

void ZFitter::getRangeZTOF_Simple(double &a_lowZTOF, double &a_highZTOF, double a_lowMassOverRide, double a_highMassOverRide){
  double centerForPOI = PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), m_currentMomOfBin);

  double lower = PhysMath::invBeta(m_cuts->getLowMassCut(m_currentPartIndex), m_currentMomOfBin);
  double higher = PhysMath::invBeta(m_cuts->getHighMassCut(m_currentPartIndex), m_currentMomOfBin);

  if(a_lowMassOverRide > 0.0){
    lower = PhysMath::invBeta(a_lowMassOverRide, m_currentMomOfBin);
  }

  if(a_highMassOverRide > 0.0){
    higher = PhysMath::invBeta(a_highMassOverRide, m_currentMomOfBin);
  }

  a_lowZTOF =  lower - centerForPOI;
  a_highZTOF = higher - centerForPOI;


  double firstLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::firstBinWithContent(m_currentHistoToFit) );
  double lastLocWithContent = m_currentHistoToFit->GetBinLowEdge( HistogramUtilities::lastBinWithContent(m_currentHistoToFit) + 1 );

  #ifdef _ZFITTER_DEBUG_
    cout << "##### getRangeZTOF: Center InvBeta: " << centerForPOI << "  LowCut: " << a_lowZTOF << " HighCut: " << a_highZTOF << " FirstContent: " << firstLocWithContent << " LastContent: " << lastLocWithContent << endl;
  #endif

  if(firstLocWithContent > a_lowZTOF) a_lowZTOF = firstLocWithContent;
  if(lastLocWithContent < a_highZTOF) a_highZTOF = lastLocWithContent;

  return;
}

//============================================================================================================

bool ZFitter::hasTOFRatioToConstrainDeDx(int a_partIndex){
  double minPercentError = 0.5;
  TH2D* tofSpectraMimic = NULL;
  TH2D* tofSpectra = NULL;
  if(m_inAllCentMode){
    tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][m_particleMimicIndicies[a_partIndex]][0][m_currentPlusMinusIndex];
    tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][a_partIndex][0][m_currentPlusMinusIndex];
  }else{
    if(m_currentPartCharge == 1){
      tofSpectraMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][m_particleMimicIndicies[a_partIndex]][m_currentCentIndex];
      tofSpectra = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][a_partIndex][m_currentCentIndex];
    }else{
      tofSpectraMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][m_particleMimicIndicies[a_partIndex]][m_currentCentIndex];
      tofSpectra = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][a_partIndex][m_currentCentIndex];
    }
  }
  double tofAmpMimic = tofSpectraMimic->GetBinContent(tofSpectraMimic->GetXaxis()->FindBin(m_currentRapidity), tofSpectraMimic->GetYaxis()->FindBin(m_currentMtM0));
  double tofAmp =  tofSpectra->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
  double ratioError = sqrt(1.0/tofAmpMimic + 1.0/tofAmp);
  if(tofAmpMimic == 0 || tofAmp == 0 || std::isnan(ratioError) || std::isinf(ratioError) || ratioError > minPercentError){
    return false;
  }else{
    return true;
  }
}


//============================================================================================================

//######      MAKE SPECTRA    ###########
void ZFitter::makeSpectra(string a_outFileName){
  cout << "#####  WRITING SPECTRA ######" << endl;
  //write everything out

  TFile* outFile = new TFile(a_outFileName.c_str(),"RECREATE");

  outFile->cd();
  outFile->mkdir("DeDx_FitData");
  outFile->cd("DeDx_FitData");
  //if(m_dEdxChiSqrHisto) m_dEdxChiSqrHisto->Write();

  for(int partBaseIndex = 0; partBaseIndex< m_numParticles; partBaseIndex++){
    if(!m_particleLoaded[partBaseIndex]) continue;
    cout << "Writing dEdx Data for " << m_partInfo->GetParticleName(partBaseIndex,0) << endl;
    outFile->cd();
    outFile->mkdir(Form("DeDx_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
    outFile->mkdir(Form("DeDx_FitData/%s/Cent",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
    outFile->mkdir(Form("DeDx_FitData/%s/Cent/mTm0_ZTPC_Proj2D",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
    //outFile->cd(Form("DeDx_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      if(m_partToFitMatrix_Plus[partBaseIndex][partIndex]){
        if(m_Fit_Data_ZTPC[partBaseIndex][partIndex][0][0]){

          outFile->cd(Form("DeDx_FitData/%s/Cent/mTm0_ZTPC_Proj2D",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
          for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
            for(int rap_index = 0; rap_index < g_ZFitter_MaxNumRapBins; rap_index++){
              HistogramUtilities::ConditionalWrite(m_mTm0_ZTPC_by_Rapidity[partIndex][centIndex][rap_index][0]);
            }
          }

          for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][0][0]);
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][1][0]);
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][2][0]);
          }

          outFile->cd(Form("DeDx_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));

          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][0][0]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][1][0]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][2][0]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][3][0]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][4][0]);

          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][0][0]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][1][0]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][2][0]);

          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPCShifts[partBaseIndex][partIndex][0]);

          if(m_useStudentTDistributionsForTPC){
            HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][3][0]);
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][3][0]);
          }
          if(m_bichselShiftsComplete) HistogramUtilities::ConditionalWrite(m_quickYieldNums_Plus[partBaseIndex][partIndex]);
          if(m_bichselShiftsComplete) HistogramUtilities::ConditionalWrite(m_bichselShifts_Plus[partBaseIndex][partIndex]);

        }

      }
      if(m_partToFitMatrix_Minus[partBaseIndex][partIndex]){
        if(m_Fit_Data_ZTPC[partBaseIndex][partIndex][0][1]){

          outFile->cd(Form("DeDx_FitData/%s/Cent/mTm0_ZTPC_Proj2D",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
          for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
            for(int rap_index = 0; rap_index < g_ZFitter_MaxNumRapBins; rap_index++){
              HistogramUtilities::ConditionalWrite(m_mTm0_ZTPC_by_Rapidity[partIndex][centIndex][rap_index][1]);
            }
          }

          outFile->cd(Form("DeDx_FitData/%s/Cent",m_partInfo->GetParticleName(partBaseIndex,0).Data()));

          for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][0][1]);
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][1][1]);
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Cent_ZTPC[partBaseIndex][partIndex][centIndex][2][1]);
          }

          outFile->cd(Form("DeDx_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));

          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][0][1]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][1][1]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][2][1]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][3][1]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][4][1]);

          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][0][1]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][1][1]);
          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][2][1]);

          HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPCShifts[partBaseIndex][partIndex][1]);

          if(m_useStudentTDistributionsForTPC){
            HistogramUtilities::ConditionalWrite(m_Fit_Data_ZTPC[partBaseIndex][partIndex][3][1]);
            HistogramUtilities::ConditionalWrite(m_Fit_Data_Solo_ZTPC[partBaseIndex][partIndex][3][1]);
          }
          if(m_bichselShiftsComplete) HistogramUtilities::ConditionalWrite(m_quickYieldNums_Minus[partBaseIndex][partIndex]);
          if(m_bichselShiftsComplete) HistogramUtilities::ConditionalWrite(m_bichselShifts_Minus[partBaseIndex][partIndex]);
        }


      }



    }

    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][0][0]);
    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][1][0]);
    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][2][0]);
    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][3][0]);
    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][0][1]);
    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][1][1]);
    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][2][1]);
    HistogramUtilities::ConditionalWrite(m_Background_FitData_ZTPC[m_currentPartIndex][3][1]);




    HistogramUtilities::ConditionalWrite(m_ChiSqr_ZTPC[partBaseIndex][0]);
    HistogramUtilities::ConditionalWrite(m_ChiSqr_ZTPC[partBaseIndex][1]);
    HistogramUtilities::ConditionalWrite(m_FitStatus_ZTPC[partBaseIndex][0]);
    HistogramUtilities::ConditionalWrite(m_FitStatus_ZTPC[partBaseIndex][1]);
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      //if(m_dEdxChiSqrMinus[partBaseIndex][centIndex]) m_dEdxChiSqrMinus[partBaseIndex][centIndex]->Write();
      //if(m_dEdxChiSqrPlus[partBaseIndex][centIndex]) m_dEdxChiSqrPlus[partBaseIndex][centIndex]->Write();
      //if(m_dEdxFitStatusMinus[partBaseIndex][centIndex]) m_dEdxFitStatusMinus[partBaseIndex][centIndex]->Write();
      //if(m_dEdxFitStatusPlus[partBaseIndex][centIndex]) m_dEdxFitStatusPlus[partBaseIndex][centIndex]->Write();

      HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZTPC[partBaseIndex][centIndex][1]);
      HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][1]);
      HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZTPC[partBaseIndex][centIndex][1]);
      HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][1]);
      if(m_hasETOF){
        HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][0]);
        HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][0]);
        HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][1]);
        HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][1]);
      }
    }

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        HistogramUtilities::ConditionalWrite(m_dEdxYieldsCentPlus[partBaseIndex][partIndex][centIndex]);
      }
    }

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        HistogramUtilities::ConditionalWrite(m_dEdxYieldsCentMinus[partBaseIndex][partIndex][centIndex]);
      }
    }

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      HistogramUtilities::discardNanHistoData(m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][0]);
    }
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      HistogramUtilities::discardNanHistoData(m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][1]);
      HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZTPC[partBaseIndex][centIndex][1]);
    }

  }


  if(m_dEdxMomFittingCompleted){
    outFile->cd();
    outFile->mkdir("DeDx_mom_FitData");
    outFile->cd("DeDx_mom_FitData");
    for(int partIndex = 0; partIndex<m_numParticles; partIndex++){
      if(partIndex ==3 || partIndex== 8) continue;

        HistogramUtilities::ConditionalWrite(m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][0]);
        HistogramUtilities::ConditionalWrite(m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]);
        HistogramUtilities::ConditionalWrite(m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]);
        HistogramUtilities::ConditionalWrite(m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][0]);
        HistogramUtilities::ConditionalWrite(m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]);
        HistogramUtilities::ConditionalWrite(m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]);

    }
  }


  cout << "Writing dEdx BetaGamma Fit Data" << endl;
  outFile->cd();
  outFile->mkdir("betagamma_FitData");
  outFile->cd("betagamma_FitData");
  for(int iii = 0; iii < 3; iii++){
    HistogramUtilities::ConditionalWrite(m_betaGammaBTOF_DeDx_Data_Plus[iii]);
    HistogramUtilities::ConditionalWrite(m_betaGammaBTOF_DeDx_Data_Minus[iii]);
    if(m_hasETOF) HistogramUtilities::ConditionalWrite(m_betaGammaETOF_DeDx_Data_Plus[iii]);
    if(m_hasETOF) HistogramUtilities::ConditionalWrite(m_betaGammaETOF_DeDx_Data_Minus[iii]);
  }
  for(int iii = 0; iii < 3; iii++){
    HistogramUtilities::ConditionalWrite(m_eta_betaGammaBTOF_DeDx_Data_Plus[iii]);
    HistogramUtilities::ConditionalWrite(m_eta_betaGammaBTOF_DeDx_Data_Minus[iii]);
    if(m_hasETOF) HistogramUtilities::ConditionalWrite(m_eta_betaGammaETOF_DeDx_Data_Plus[iii]);
    if(m_hasETOF) HistogramUtilities::ConditionalWrite(m_eta_betaGammaETOF_DeDx_Data_Minus[iii]);
  }

  //Write functional fits
  HistogramUtilities::ConditionalWrite(m_betaGammaBTOF_dEdx_FitFunction[0]);
  HistogramUtilities::ConditionalWrite(m_betaGammaBTOF_dEdx_FitFunction[1]);
  for(int etaIndex = 0; etaIndex < (int) m_eta_betaGammaBTOF_dEdx_FitFunction[0].size(); etaIndex++){
      HistogramUtilities::ConditionalWrite(m_eta_betaGammaBTOF_dEdx_FitFunction[0][etaIndex]);
  }
  for(int etaIndex = 0; etaIndex < (int) m_eta_betaGammaBTOF_dEdx_FitFunction[1].size(); etaIndex++){
      HistogramUtilities::ConditionalWrite(m_eta_betaGammaBTOF_dEdx_FitFunction[1][etaIndex]);
  }
  for(int etaIndex = 0; etaIndex < (int)  m_eta_betaGammaETOF_dEdx_FitFunction[0].size(); etaIndex++){
      HistogramUtilities::ConditionalWrite(m_eta_betaGammaETOF_dEdx_FitFunction[0][etaIndex]);
  }
  for(int etaIndex = 0; etaIndex < (int) m_eta_betaGammaETOF_dEdx_FitFunction[1].size(); etaIndex++){
      HistogramUtilities::ConditionalWrite(m_eta_betaGammaETOF_dEdx_FitFunction[1][etaIndex]);
  }

  for(int iii = 0; iii < (int) m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData.size(); iii++){
   //cout << "size: " << m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[iii].size() << endl;
   for(int jjj = 0; jjj < (int) m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[iii].size(); jjj++){
      //cout << "Writing " << iii << " " << jjj << endl;
      //cout << "  ptr:" << m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[iii][jjj] << endl;
      //cout << "  name:" << m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[iii][jjj]->GetName() << endl;
    HistogramUtilities::ConditionalWrite(m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[iii][jjj]);
      //if(m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii][jjj]) m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii][jjj]->Write();
  }
}
for(int iii = 0; iii < (int) m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData.size(); iii++){
   //cout << "size: " << m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii].size() << endl;
 for(int jjj = 0; jjj < (int) m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii].size(); jjj++){
      //cout << "Writing " << iii << " " << jjj << endl;
      //cout << "  ptr:" << m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii][jjj] << endl;
      //cout << "  name:" << m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii][jjj]->GetName() << endl;
    HistogramUtilities::ConditionalWrite(m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[iii][jjj]);
  }
}

if(m_hasETOF){

  for(int iii = 0; iii < (int) m_eta_betaGammaETOF_DeDx_Plus_SmoothFitData.size(); iii++){
    for(int jjj = 0; jjj < (int) m_eta_betaGammaETOF_DeDx_Plus_SmoothFitData[iii].size(); jjj++){
      HistogramUtilities::ConditionalWrite(m_eta_betaGammaETOF_DeDx_Plus_SmoothFitData[iii][jjj]);
      HistogramUtilities::ConditionalWrite(m_eta_betaGammaETOF_DeDx_Minus_SmoothFitData[iii][jjj]);
    }
  }
}

outFile->cd();
outFile->mkdir("BTOF_FitData");
outFile->cd("BTOF_FitData");
for(int partBaseIndex = 0; partBaseIndex< m_numParticles; partBaseIndex++){
  if(!m_particleLoaded[partBaseIndex]) continue;
  cout << "Writing BTOF Fit Data for " << m_partInfo->GetParticleName(partBaseIndex,0) << endl;
  outFile->cd();
  outFile->mkdir(Form("BTOF_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
  outFile->cd(Form("BTOF_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));

  cout << "    writing background fit data" << endl;


  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][0][0]);
  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][1][0]);
  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][2][0]);
  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][3][0]);
  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][0][1]);
  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][1][1]);
  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][2][1]);
  HistogramUtilities::ConditionalWrite(m_Background_FitData_ZbTOF[partBaseIndex][3][1]);


  cout << "    writing normal fit data" << endl;
  for(int partIndex = 0; partIndex<m_numParticles; partIndex++){
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][0][0]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][1][0]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][2][0]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][3][0]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][4][0]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][0][1]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][1][1]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][2][1]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][3][1]);
    HistogramUtilities::ConditionalWrite(m_Fit_Data_ZbTOF[partBaseIndex][partIndex][4][1]);

  }

  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      HistogramUtilities::ConditionalWrite(m_InvBetaBTOFYieldsCentPlus[partBaseIndex][partIndex][centIndex]);
    }
  }

  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      HistogramUtilities::ConditionalWrite(m_InvBetaBTOFYieldsCentMinus[partBaseIndex][partIndex][centIndex]);
    }
  }

  HistogramUtilities::ConditionalWrite(m_ChiSqr_ZbTOF[partBaseIndex][0]);
  HistogramUtilities::ConditionalWrite(m_ChiSqr_ZbTOF[partBaseIndex][1]);
  HistogramUtilities::ConditionalWrite(m_FitStatus_ZbTOF[partBaseIndex][0]);
  HistogramUtilities::ConditionalWrite(m_FitStatus_ZbTOF[partBaseIndex][1]);
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][0]);
    HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZbTOF[partBaseIndex][centIndex][1]);
    HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][0]);
    HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZbTOF[partBaseIndex][centIndex][1]);
    HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][0]);
    HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][1]);
  }
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    if(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][0]) HistogramUtilities::discardNanHistoData(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][0]);
    HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][0]);
  }

  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    if(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][1]) HistogramUtilities::discardNanHistoData(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][1]);
    HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZbTOF[partBaseIndex][centIndex][1]);
  }

}

if(m_hasETOF){
  outFile->cd();
  outFile->mkdir("ETOF_FitData");
  outFile->cd("ETOF_FitData");
  for(int partBaseIndex = 0; partBaseIndex< m_numParticles; partBaseIndex++){
    if(!m_particleLoaded[partBaseIndex]) continue;
    cout << "Writing ETOF Fit Data for " << m_partInfo->GetParticleName(partBaseIndex,0) << endl;
    outFile->cd();
    outFile->mkdir(Form("ETOF_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
    outFile->cd(Form("ETOF_FitData/%s",m_partInfo->GetParticleName(partBaseIndex,0).Data()));
    for(int partIndex = 0; partIndex<m_numParticles; partIndex++){

      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][0][0]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][1][0]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][2][0]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][3][0]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][4][0]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][0][1]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][1][1]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][2][1]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][3][1]);
      HistogramUtilities::ConditionalWrite(m_Fit_Data_ZeTOF[partBaseIndex][partIndex][4][1]);
    }
    HistogramUtilities::ConditionalWrite(m_ChiSqr_ZeTOF[partBaseIndex][0]);
    HistogramUtilities::ConditionalWrite(m_ChiSqr_ZeTOF[partBaseIndex][1]);
    HistogramUtilities::ConditionalWrite(m_FitStatus_ZeTOF[partBaseIndex][0]);
    HistogramUtilities::ConditionalWrite(m_FitStatus_ZeTOF[partBaseIndex][1]);
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_ChiSqr_Cent_ZeTOF[partBaseIndex][centIndex][1]);
      HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_FitStatus_Cent_ZeTOF[partBaseIndex][centIndex][1]);
      if(m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][0]) HistogramUtilities::discardNanHistoData(m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][0]);
      HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][0]);
      if(m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][1]) HistogramUtilities::discardNanHistoData(m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][1]);
      HistogramUtilities::ConditionalWrite(m_Spectra_Cent_ZeTOF[partBaseIndex][centIndex][1]);
    }

  }
}
cout << "DONE Writing, Closing File" << endl;
outFile->Close();
cout << "FILE Closed" << endl;
return;
}



TH1D* ZFitter::getZTPC_Example(){
  double horLow = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin);
  double horHigh = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1);
  double vertLow = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin);
  double vertHigh = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1);
  double horCenter = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
  double vertCenter = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
  double possibleMom[4];

  possibleMom[0] =vertLow * TMath::CosH( horLow ) ;
  possibleMom[1] =vertHigh * TMath::CosH( horLow );
  possibleMom[2] =vertLow * TMath::CosH( horHigh) ;
  possibleMom[3] =vertHigh * TMath::CosH( horHigh);


  cout << "HorBin = " << m_currentRapBin << "  VertBin = " << m_currentMtM0Bin << endl;
  cout << "   momentum = " << possibleMom[0] << "  "<< possibleMom[1] << "  "<< possibleMom[2] << "  "<< possibleMom[3] << endl;

  double momLow = 1E7;
  double momHigh = -999;
  for(int iii = 0; iii < 4; iii++){
    if(possibleMom[iii] > momHigh) momHigh = possibleMom[iii];
    if(possibleMom[iii] < momLow) momLow = possibleMom[iii];
  }
  cout << "   momentum Low = " << momLow << endl;
  cout << "   momentum High= " << momHigh << endl;


  int firstBin = m_mom_dEdx_Plus->GetXaxis()->FindBin(momLow);
  int lastBin = m_mom_dEdx_Plus->GetXaxis()->FindBin(momHigh);
  cout << "   firstBin= " << firstBin << "    lastBin=" << lastBin << endl;
  //TH1D* proj = m_dEdx_LogLog->ProjectionY("dEdx_proj_loglog",1,1);
  //HistogramUtilities::zeroBins(proj);

  TH1D* shifted = new TH1D("dEdx_shifted_toZTPC","dEdx_shifted_toZTPC",500,-3.0,3.0);


  for(int momBin = firstBin; momBin <= lastBin; momBin++){
    double momentum = m_mom_dEdx_Plus->GetXaxis()->GetBinCenter(momBin);
    for(int dEdxBin = 1; dEdxBin <= m_mom_dEdx_Plus->GetNbinsY(); dEdxBin++){
      double dEdx = m_mom_dEdx_Plus->GetYaxis()->GetBinCenter(dEdxBin);

      double newContent = m_mom_dEdx_Plus->GetBinContent(momBin,dEdxBin);
      shifted->Fill(dEdx,newContent);
      //shifted->AddBinContent(shifted->GetXaxis()->FindBin(ZTPC), m_mom_dEdx_Plus->GetBinContent(momBin,dEdxBin));
    }
  }
  return shifted;
}


void ZFitter::fitZTPC_Electrons(){
 TFile f("ztpc_examples.root","RECREATE");
 int numHorBins = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetNbinsX();
 int numVertBins = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetNbinsY();
 for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
  m_currentRapBin = HorBinIndex;
  for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
    m_currentMtM0Bin = VertBinIndex;
    double horCenter = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
    double vertCenter = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
    double momentum = PhysMath::pTot_from_y_mTm0(horCenter, vertCenter, m_partInfo->GetParticleMass(m_currentPartIndex));
    TH1D* example = getZTPC_Example();
    cout << "    momentum= " << momentum << endl;
    example->SetName(Form("ztpc_example_Eta_%1.2f_pT_%1.3f",horCenter,vertCenter));
    example->SetTitle(Form("ZTPC Example Momentum = %f",momentum));
    f.cd();
    example->Write();

  }
}
return;
}


//--------------------------------------------------------------------
TLine*  ZFitter::PredictionLine(TLine* line, Double_t prediction, Double_t height){
  line->SetX1(prediction);
  line->SetX2(prediction);
  line->SetY1(0);
  line->SetY2(height);
  return line;
}



void ZFitter::DrawPredictionLines_dEdx(TCanvas *canvas){
  canvas->cd();

  //ZTPC predictions for line drawing
  double pionPredict = log(m_bichselCurves[0]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double kaonPredict = log(m_bichselCurves[1]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double protonPredict = log(m_bichselCurves[2]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double electronPredict = log(m_bichselCurves[3]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double deuteronPredict = log(m_bichselCurves[4]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double tritonPredict = log(m_bichselCurves[5]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double helionPredict = log(m_bichselCurves[6]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double alphaPredict = log(m_bichselCurves[7]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  double lineHeight = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(0.0));

  TLine *pionLine = new TLine();
  TLine *kaonLine = new TLine();
  TLine *protonLine = new TLine();
  TLine *electronLine = new TLine();
  TLine *deuteronLine = new TLine();
  TLine *tritonLine = new TLine();
  TLine *helionLine = new TLine();
  TLine *alphaLine = new TLine();

  m_TLinePtrs.push_back(pionLine);
  m_TLinePtrs.push_back(kaonLine);
  m_TLinePtrs.push_back(protonLine);
  m_TLinePtrs.push_back(electronLine);
  m_TLinePtrs.push_back(deuteronLine);
  m_TLinePtrs.push_back(tritonLine);
  m_TLinePtrs.push_back(helionLine);
  m_TLinePtrs.push_back(alphaLine);

  pionLine->SetLineColor(m_partInfo->GetParticleColor(PION));
  kaonLine->SetLineColor(m_partInfo->GetParticleColor(KAON));
  protonLine->SetLineColor(m_partInfo->GetParticleColor(PROTON));
  electronLine->SetLineColor(m_partInfo->GetParticleColor(ELECTRON));
  deuteronLine->SetLineColor(m_partInfo->GetParticleColor(DEUTERON));
  tritonLine->SetLineColor(m_partInfo->GetParticleColor(TRITON));
  helionLine->SetLineColor(m_partInfo->GetParticleColor(HELION));
  alphaLine->SetLineColor(m_partInfo->GetParticleColor(ALPHA));

  pionLine->SetLineWidth(2);
  kaonLine->SetLineWidth(2);
  protonLine->SetLineWidth(2);
  electronLine->SetLineWidth(2);
  deuteronLine->SetLineWidth(2);
  tritonLine->SetLineWidth(2);
  helionLine->SetLineWidth(2);
  alphaLine->SetLineWidth(2);

  PredictionLine(pionLine,pionPredict,lineHeight);
  PredictionLine(kaonLine,kaonPredict,lineHeight);
  PredictionLine(protonLine,protonPredict,lineHeight);
  PredictionLine(electronLine,electronPredict,lineHeight);
  PredictionLine(deuteronLine,deuteronPredict,lineHeight);
  PredictionLine(tritonLine,tritonPredict,lineHeight);
  PredictionLine(helionLine,helionPredict,lineHeight);
  PredictionLine(alphaLine,alphaPredict,lineHeight);


  if(pionPredict > m_currentHistoToFit->GetBinLowEdge(1)     && pionPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))     pionLine->Draw("SAME");
  if(kaonPredict > m_currentHistoToFit->GetBinLowEdge(1)     && kaonPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))     kaonLine->Draw("SAME");
  if(protonPredict > m_currentHistoToFit->GetBinLowEdge(1)   && protonPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))   protonLine->Draw("SAME");
  if(electronPredict > m_currentHistoToFit->GetBinLowEdge(1) && electronPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1)) electronLine->Draw("SAME");
  if(deuteronPredict > m_currentHistoToFit->GetBinLowEdge(1) && deuteronPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1)) deuteronLine->Draw("SAME");
  if(tritonPredict > m_currentHistoToFit->GetBinLowEdge(1)   && tritonPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))   tritonLine->Draw("SAME");
  if(helionPredict > m_currentHistoToFit->GetBinLowEdge(1)   && helionPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))   helionLine->Draw("SAME");
  if(alphaPredict > m_currentHistoToFit->GetBinLowEdge(1)    && alphaPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))    alphaLine->Draw("SAME");
}


void ZFitter::DrawPredictionLines_TOF(TCanvas *canvas){
  canvas->cd();

  //ZTPC predictions for line drawing
  double pionPredict = m_invBetaCurves[0]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double kaonPredict =  m_invBetaCurves[1]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double protonPredict =  m_invBetaCurves[2]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double electronPredict =  m_invBetaCurves[3]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double deuteronPredict =  m_invBetaCurves[4]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double tritonPredict =  m_invBetaCurves[5]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double helionPredict =  m_invBetaCurves[6]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double alphaPredict =  m_invBetaCurves[7]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
  double lineHeight = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(0.0));


  if(m_massSquaredMode){
    pionPredict = pow(m_partInfo->GetParticleMass(0),2.0);
    kaonPredict =  pow(m_partInfo->GetParticleMass(1),2.0);
    protonPredict =  pow(m_partInfo->GetParticleMass(2),2.0);
    electronPredict =  pow(m_partInfo->GetParticleMass(3),2.0);
    deuteronPredict =  pow(m_partInfo->GetParticleMass(4),2.0);
    tritonPredict =  pow(m_partInfo->GetParticleMass(5),2.0);
    helionPredict =  pow(m_partInfo->GetParticleMass(6),2.0);
    alphaPredict =  pow(m_partInfo->GetParticleMass(7),2.0);
    lineHeight = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0)));
  }

  TLine *pionLine = new TLine();
  TLine *kaonLine = new TLine();
  TLine *protonLine = new TLine();
  TLine *electronLine = new TLine();
  TLine *deuteronLine = new TLine();
  TLine *tritonLine = new TLine();
  TLine *helionLine = new TLine();
  TLine *alphaLine = new TLine();

  m_TLinePtrs.push_back(pionLine);
  m_TLinePtrs.push_back(kaonLine);
  m_TLinePtrs.push_back(protonLine);
  m_TLinePtrs.push_back(electronLine);
  m_TLinePtrs.push_back(deuteronLine);
  m_TLinePtrs.push_back(tritonLine);
  m_TLinePtrs.push_back(helionLine);
  m_TLinePtrs.push_back(alphaLine);

  pionLine->SetLineColor(m_partInfo->GetParticleColor(PION));
  kaonLine->SetLineColor(m_partInfo->GetParticleColor(KAON));
  protonLine->SetLineColor(m_partInfo->GetParticleColor(PROTON));
  electronLine->SetLineColor(m_partInfo->GetParticleColor(ELECTRON));
  deuteronLine->SetLineColor(m_partInfo->GetParticleColor(DEUTERON));
  tritonLine->SetLineColor(m_partInfo->GetParticleColor(TRITON));
  helionLine->SetLineColor(m_partInfo->GetParticleColor(HELION));
  alphaLine->SetLineColor(m_partInfo->GetParticleColor(ALPHA));

  pionLine->SetLineWidth(2);
  kaonLine->SetLineWidth(2);
  protonLine->SetLineWidth(2);
  electronLine->SetLineWidth(2);
  deuteronLine->SetLineWidth(2);
  tritonLine->SetLineWidth(2);
  helionLine->SetLineWidth(2);
  alphaLine->SetLineWidth(2);

  //GetThe right height
  PredictionLine(pionLine,pionPredict,lineHeight);
  PredictionLine(kaonLine,kaonPredict,lineHeight);
  PredictionLine(protonLine,protonPredict,lineHeight);
  PredictionLine(electronLine,electronPredict,lineHeight);
  PredictionLine(deuteronLine,deuteronPredict,lineHeight);
  PredictionLine(tritonLine,tritonPredict,lineHeight);
  PredictionLine(helionLine,helionPredict,lineHeight);
  PredictionLine(alphaLine,alphaPredict,lineHeight);


  if(pionPredict > m_currentHistoToFit->GetBinLowEdge(1)     && pionPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))     pionLine->Draw("SAME");
  if(kaonPredict > m_currentHistoToFit->GetBinLowEdge(1)     && kaonPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))     kaonLine->Draw("SAME");
  if(protonPredict > m_currentHistoToFit->GetBinLowEdge(1)   && protonPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))   protonLine->Draw("SAME");
  if(electronPredict > m_currentHistoToFit->GetBinLowEdge(1) && electronPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1)) electronLine->Draw("SAME");
  if(deuteronPredict > m_currentHistoToFit->GetBinLowEdge(1) && deuteronPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1)) deuteronLine->Draw("SAME");
  if(tritonPredict > m_currentHistoToFit->GetBinLowEdge(1)   && tritonPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))   tritonLine->Draw("SAME");
  if(helionPredict > m_currentHistoToFit->GetBinLowEdge(1)   && helionPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))   helionLine->Draw("SAME");
  if(alphaPredict > m_currentHistoToFit->GetBinLowEdge(1)    && alphaPredict < m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1))    alphaLine->Draw("SAME");

  //cout << "Drawing Prediction Lines with line height: " << lineHeight << endl;
}





void ZFitter::findPartShiftsAndRatiosWithRap(vector< double > a_momToFitByParticle){
  cout << "Finding Particle Shifts and Ratios by Fitting dE/dx by momentum" << endl;
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  //make the data histograms
  for(int partBaseIndex = 0; partBaseIndex < 9; partBaseIndex++){
    if(!m_particleLoaded[partBaseIndex]) continue;
    TH1D* structHisto = m_rap_mTm0_dEdx_allCent[partBaseIndex]->ProjectionX("structHistoBichShifts",1,1,1,1);
    HistogramUtilities::zeroBins(structHisto);

    for(int partIndex = 0; partIndex < 9; partIndex++){
      if(m_partToFitMatrix_Plus[partBaseIndex][partIndex]){
        m_quickYieldNums_Plus[partBaseIndex][partIndex] = (TH1D*) structHisto->Clone();
        m_quickYieldNums_Plus[partBaseIndex][partIndex]->SetName(Form("quickYields_%s_for_%s",m_partInfo->GetParticleName(partBaseIndex,1).Data(), m_partInfo->GetParticleName(partIndex,1).Data()));
        m_quickYieldNums_Plus[partBaseIndex][partIndex]->SetName(Form("Quick Yield for %s in %s Space",m_partInfo->GetParticleSymbol(partIndex,1).Data(), m_partInfo->GetParticleSymbol(partBaseIndex,1).Data()));
        m_bichselShifts_Plus[partBaseIndex][partIndex] = (TH1D*) structHisto->Clone();
        m_bichselShifts_Plus[partBaseIndex][partIndex]->SetName(Form("bichselShifts_%s_for_%s",m_partInfo->GetParticleName(partBaseIndex,1).Data(), m_partInfo->GetParticleName(partIndex,1).Data()));
        m_bichselShifts_Plus[partBaseIndex][partIndex]->SetName(Form("Deviation in Bichsel Curve for %s in %s Space",m_partInfo->GetParticleSymbol(partIndex,1).Data(), m_partInfo->GetParticleSymbol(partBaseIndex,1).Data()));
        m_tpcResolution_Plus[partBaseIndex][partIndex] = (TH1D*) structHisto->Clone();
        m_tpcResolution_Plus[partBaseIndex][partIndex]->SetName(Form("tpcResolution_%s_for_%s",m_partInfo->GetParticleName(partBaseIndex,1).Data(), m_partInfo->GetParticleName(partIndex,1).Data()));
        m_tpcResolution_Plus[partBaseIndex][partIndex]->SetName(Form("TPC Resolution for %s in %s Space",m_partInfo->GetParticleSymbol(partIndex,1).Data(), m_partInfo->GetParticleSymbol(partBaseIndex,1).Data()));

      }
      if(m_partToFitMatrix_Minus[partBaseIndex][partIndex]){
        m_bichselShifts_Minus[partBaseIndex][partIndex] = (TH1D*) structHisto->Clone();
        m_bichselShifts_Minus[partBaseIndex][partIndex]->SetName(Form("bichselShifts_%s_for_%s",m_partInfo->GetParticleName(partBaseIndex,-1).Data(), m_partInfo->GetParticleName(partIndex,-1).Data()));
        m_bichselShifts_Minus[partBaseIndex][partIndex]->SetName(Form("Deviation in Bichsel Curve for %s in %s Space",m_partInfo->GetParticleSymbol(partIndex,-1).Data(), m_partInfo->GetParticleSymbol(partBaseIndex,-1).Data()));
        m_quickYieldNums_Minus[partBaseIndex][partIndex] = (TH1D*) structHisto->Clone();
        m_quickYieldNums_Minus[partBaseIndex][partIndex]->SetName(Form("quickYields_%s_for_%s",m_partInfo->GetParticleName(partBaseIndex,-1).Data(), m_partInfo->GetParticleName(partIndex,-1).Data()));
        m_quickYieldNums_Minus[partBaseIndex][partIndex]->SetName(Form("Quick Yield for %s in %s Space",m_partInfo->GetParticleSymbol(partIndex,-1).Data(), m_partInfo->GetParticleSymbol(partBaseIndex,-1).Data()));
        m_tpcResolution_Minus[partBaseIndex][partIndex] = (TH1D*) structHisto->Clone();
        m_tpcResolution_Minus[partBaseIndex][partIndex]->SetName(Form("tpcResolution_%s_for_%s",m_partInfo->GetParticleName(partBaseIndex,-1).Data(), m_partInfo->GetParticleName(partIndex,-1).Data()));
        m_tpcResolution_Minus[partBaseIndex][partIndex]->SetName(Form("TPC Resolution for %s in %s Space",m_partInfo->GetParticleSymbol(partIndex,-1).Data(), m_partInfo->GetParticleSymbol(partBaseIndex,-1).Data()));
      }
    }
  }

  //make the data histograms
  for(int partBaseIndex = 0; partBaseIndex < 9; partBaseIndex++){
    if(!m_particleLoaded[partBaseIndex]) continue;

    for(int partIndex = 0; partIndex < 9; partIndex++){
      if(partIndex == 3 || partIndex == 8) continue;
      //if(partBaseIndex != partIndex) continue;
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        for(int iii = 0; iii < 9; iii++) m_currentParticleTags[iii] = 0;
        m_currentParticleTags[partIndex] = 1;

        m_currentPartCharge = -1;
        if(pmIndex == 0) m_currentPartCharge = 1;
        if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[partBaseIndex][partIndex])
           || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[partBaseIndex][partIndex])) continue;

        TH3D* histogram3D = NULL;
        if( m_currentPartCharge == 1 ){
          //histogram3D = m_rap_mTm0_dEdx_AllCentRead_Plus[partBaseIndex];
          histogram3D = m_rap_mTm0_dEdx_Plus_allCent[partBaseIndex];
        }else{
          //histogram3D = m_rap_mTm0_dEdx_AllCentRead_Minus[partBaseIndex];
          histogram3D = m_rap_mTm0_dEdx_Minus_allCent[partBaseIndex];
        }

        for(int rapBin = 1; rapBin <= histogram3D->GetNbinsX(); rapBin++){
          m_currentRapBin = rapBin;
          m_currentRapidity = histogram3D->GetXaxis()->GetBinCenter(m_currentRapBin);
          m_currentMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity, a_momToFitByParticle[partIndex], m_partInfo->GetParticleMass(partBaseIndex));
          m_currentMtM0Bin = histogram3D->GetYaxis()->FindBin(m_currentMtM0);

          if(m_currentMtM0Bin < 1 || m_currentMtM0Bin > histogram3D->GetNbinsY()){
            cout << " Momentum is outside mTm0 bin range..." << endl;
            continue;
          }

          cout << " Working on " << m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()
           << " in " << m_partInfo->GetParticleName(partBaseIndex,m_currentPartCharge).Data() << " space  : RapBin : " << m_currentRapBin << endl;

          int binToGet = histogram3D->GetYaxis()->FindBin(m_currentMtM0);
          m_currentHistoToFit = histogram3D->ProjectionZ(Form("proj_%02d_%02d_%d_%02d",partBaseIndex,partIndex,pmIndex,m_currentRapBin),m_currentRapBin,m_currentRapBin,m_currentMtM0Bin,m_currentMtM0Bin);



          //there is a wierd effect where the a dEdx bin near the electron peak gets orders of magnitude more counts then goes back
          // this finds that weird bin and zeros it
          for(unsigned int projBin = 1; projBin <= m_currentHistoToFit->GetNbinsX(); projBin++ ){
            //check if bin before and after have content
            if(m_currentHistoToFit->GetBinContent(projBin-1) <= 0.0 ||
               m_currentHistoToFit->GetBinContent(projBin+1) <= 0.0) continue;
            double average = (m_currentHistoToFit->GetBinContent(projBin-1) + m_currentHistoToFit->GetBinContent(projBin+1) ) * 0.5;
            if(m_currentHistoToFit->GetBinContent(projBin)/average > 1E2){
              m_currentHistoToFit->SetBinContent(projBin,0.0);
              m_currentHistoToFit->SetBinError(projBin,0.0);
            }
          }





          double sigAllow = 1.0;
          double prediction = m_bichselCurves[partIndex]->Eval(a_momToFitByParticle[partIndex]);
          double predictionBase = m_bichselCurves[partBaseIndex]->Eval(a_momToFitByParticle[partIndex]);
          double lowZTPC = log((1.0 - sigAllow*0.07)*prediction/predictionBase);
          double highZTPC = log((1.0 + (sigAllow + 1.5)*0.07)*prediction/predictionBase);


          double meanInitial = log(prediction/predictionBase);
          int meanBin = m_currentHistoToFit->FindBin(meanInitial);
          cout << " mean initial : "  << meanInitial << "     mean Bin: " << meanBin << endl;
          //correct the mean to the highest bin in the area
          /*for(unsigned int binDiff = 0; meanBin + binDiff + 2 <= m_currentHistoToFit->GetNbinsX() && meanBin - binDiff > 2; binDiff++){
            //check if bin is the max of the average around it
            double average = 0.25*(m_currentHistoToFit->GetBinContent(meanBin + binDiff-1) + m_currentHistoToFit->GetBinContent(meanBin + binDiff+1)
                              + m_currentHistoToFit->GetBinContent(meanBin + binDiff-2) + m_currentHistoToFit->GetBinContent(meanBin + binDiff+2)) ;
            if(m_currentHistoToFit->GetBinContent(meanBin + binDiff) > average){
              //meanBin = meanBin + binDiff;
              meanInitial = m_currentHistoToFit->GetBinCenter(meanBin + binDiff);
              break;
            }else{ // do the left hand search
              average = 0.25*(m_currentHistoToFit->GetBinContent(meanBin - binDiff-1) + m_currentHistoToFit->GetBinContent(meanBin - binDiff+1)
                                + m_currentHistoToFit->GetBinContent(meanBin - binDiff-2) + m_currentHistoToFit->GetBinContent(meanBin - binDiff+2)) ;
              if(m_currentHistoToFit->GetBinContent(meanBin - binDiff) > average){
                //meanBin = meanBin - binDiff;
                meanInitial = m_currentHistoToFit->GetBinCenter(meanBin - binDiff);
                break;
              }
            }
          }*/

          double maxRightBin = meanBin;
          double maxRightBinContent = m_currentHistoToFit->GetBinContent(maxRightBin);
          double maxLeftBin = meanBin;
          double maxLeftBinContent = m_currentHistoToFit->GetBinContent(maxLeftBin);

          for(int binDiff = 0; meanBin + binDiff <= m_currentHistoToFit->GetNbinsX() ; binDiff++){
            if(m_currentHistoToFit->GetBinContent(meanBin + binDiff) > maxRightBinContent){
              maxRightBinContent = m_currentHistoToFit->GetBinContent(meanBin + binDiff);
              maxRightBin = meanBin + binDiff;
            }else{
              break;
            }
          }

          for(int binDiff = 0; meanBin - binDiff > 0 ; binDiff--){
            if(m_currentHistoToFit->GetBinContent(meanBin - binDiff) > maxRightBinContent){
              maxRightBinContent = m_currentHistoToFit->GetBinContent(meanBin - binDiff);
              maxRightBin = meanBin - binDiff;
            }else{
              break;
            }
          }

          if(maxRightBinContent > maxLeftBinContent){
            meanInitial = m_currentHistoToFit->GetBinCenter(maxRightBin);
          }else{
            meanInitial = m_currentHistoToFit->GetBinCenter(maxLeftBin);
          }



          cout << " mean after : " << meanInitial << "     mean Bin after: " << meanBin << endl;


          double widthInitial = log(1.07);
          double countsInitial = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(meanInitial)) / m_currentHistoToFit->GetBinWidth(m_currentHistoToFit->FindBin(meanInitial)) /gaus(meanInitial, widthInitial, meanInitial);
          double x_low, x_mid, x_high = 0;
          if(m_particleLeftRightSeedDirection[partIndex] == -1){
            x_low = meanInitial - widthInitial;
            x_mid = meanInitial - 0.5*widthInitial;
            x_high = meanInitial;
          }else if(m_particleLeftRightSeedDirection[partIndex] == 1){
            x_low = meanInitial + widthInitial;
            x_mid = meanInitial + 0.5*widthInitial;
            x_high = meanInitial;
          }else{
            x_low = meanInitial - widthInitial;
            x_mid = meanInitial;
            x_high = meanInitial + widthInitial;
          }

          //double x_low = meanInitial - widthInitial;
          int binX_Low = m_currentHistoToFit->FindBin(x_low);
          int binX_Mid = m_currentHistoToFit->FindBin(x_mid);
          int binX_High =m_currentHistoToFit->FindBin(x_high);
          double y_low  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, binX_Low  , 2) / m_currentHistoToFit->GetBinWidth( binX_Low ) ;
          double y_mid  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, binX_Mid  , 2) / m_currentHistoToFit->GetBinWidth( binX_Mid ) ;
          double y_high  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, binX_High , 2) / m_currentHistoToFit->GetBinWidth(binX_High ) ;


          //double y_low = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(x_low)) / m_currentHistoToFit->GetBinWidth(m_currentHistoToFit->FindBin(x_low));
          ////double x_mid = meanInitial;
          //double y_mid = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(x_mid)) / m_currentHistoToFit->GetBinWidth(m_currentHistoToFit->FindBin(x_mid));
          ////double x_high = meanInitial + widthInitial;
          //double y_high = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(x_high)) / m_currentHistoToFit->GetBinWidth(m_currentHistoToFit->FindBin(x_high));

          double N_init, mu_init, sig_init = 0;
          PhysMath::get3PointGausParams(N_init, mu_init, sig_init, x_low , y_low , x_mid , y_mid , x_high, y_high);

          cout << "Before : N: " << countsInitial << "  mu: " << meanInitial << "   widthInitial: " << widthInitial << endl;
          cout << "After : N: " << N_init << "  mu: " << mu_init << "   widthInitial: " << sig_init << endl;

          if(sig_init > 0.2 || fabs(log10(N_init/countsInitial)) > 2){
            cout << " Won't be changing the seeds..." << endl;
            cout << " logPar = " << fabs(log10(N_init/countsInitial)) << " and it is greater than 2" << endl;
          }else{
            meanInitial = mu_init;
            widthInitial = sig_init;
            countsInitial = N_init;
          }

          cout << " ptr: " << m_currentHistoToFit << endl;

          if(countsInitial < 1){
            cout << "Not Enough Counts: " << countsInitial << endl;
            #ifndef _SAVE_POINTERS_
              delete m_currentHistoToFit;
            #endif
            continue;
          }
          cout << " here 1 " << endl;

          if(std::isnan(meanInitial) || std::isnan(widthInitial) || std::isnan(countsInitial)){
            cout << "Bad Initialization: " << endl;
            #ifndef _SAVE_POINTERS_
              delete m_currentHistoToFit;
            #endif
            continue;
          }

          HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
          m_currentHistoToFit->Scale(1.0,"width");

          cout << " here 2 " << endl;

          /*double gaussData[3];
          HistogramUtilities::getGausData(gaussData, m_currentHistoToFit, true, true, lowZTPC, highZTPC);
          countsInitial = gaussData[0];
          meanInitial = gaussData[1];
          widthInitial = gaussData[2];
          lowZTPC = log(1.0 - sigAllow*0.07) + meanInitial;
          highZTPC = log(1.0 + (sigAllow + 1.5)*0.07) + meanInitial;
          */
          lowZTPC = meanInitial - sigAllow*widthInitial;
          highZTPC = meanInitial + sigAllow*widthInitial;
          TF1* gausFunct = (TF1*) m_gausFunct->Clone();
          string gausName = Form("gaus_%02d,%02d,%02d,%02d",partBaseIndex,partIndex,pmIndex,rapBin);
          gausFunct->SetName(gausName.c_str());
          gausFunct->SetParameter(0, countsInitial);
          gausFunct->SetParameter(1, meanInitial);
          gausFunct->SetParameter(2, widthInitial);
          gausFunct->SetParLimits(0, 0.1*countsInitial,10.0*countsInitial);
          gausFunct->SetParLimits(1, meanInitial - 0.3, meanInitial + 0.5);
          gausFunct->SetParLimits(2, 0.3 * widthInitial,5.0 * widthInitial);
          gausFunct->SetParNames("N","#mu","#sigma");
          gausFunct->SetRange(lowZTPC,highZTPC);
          cout << " here 3 " << endl;
          //m_currentHistoToFit->GetXaxis()->SetRangeUser(lowZTPC,highZTPC);
          TF1* initFunct = (TF1*) gausFunct->Clone();
          initFunct->SetLineColor(kYellow);
          cout << " fitting for shifts" << endl;
          int status = m_currentHistoToFit->Fit(gausName.c_str(),"","",lowZTPC,highZTPC);
          TF1* fitResult = m_currentHistoToFit->GetFunction(gausName.c_str());
          if(!fitResult){
            cout << "There is no fit result, continuing" << endl;
            continue;
          }
          double shift = fitResult->GetParameter(1) - log(prediction/predictionBase);
          cout << "  pushing back results" << endl;
          if(m_currentPartCharge == 1){
            m_quickYieldNums_Plus[partBaseIndex][partIndex]->SetBinContent(m_currentRapBin, fitResult->GetParameter(0));
            m_quickYieldNums_Plus[partBaseIndex][partIndex]->SetBinError(m_currentRapBin, fitResult->GetParError(0));
            m_bichselShifts_Plus[partBaseIndex][partIndex]->SetBinContent(m_currentRapBin, shift);
            m_bichselShifts_Plus[partBaseIndex][partIndex]->SetBinError(m_currentRapBin, fitResult->GetParError(1));
            m_tpcResolution_Plus[partBaseIndex][partIndex]->SetBinContent(m_currentRapBin, fitResult->GetParameter(2));
            m_tpcResolution_Plus[partBaseIndex][partIndex]->SetBinError(m_currentRapBin, fitResult->GetParError(2));
          }else{
            m_quickYieldNums_Minus[partBaseIndex][partIndex]->SetBinContent(m_currentRapBin, fitResult->GetParameter(0));
            m_quickYieldNums_Minus[partBaseIndex][partIndex]->SetBinError(m_currentRapBin, fitResult->GetParError(0));
            m_bichselShifts_Minus[partBaseIndex][partIndex]->SetBinContent(m_currentRapBin, shift);
            m_bichselShifts_Minus[partBaseIndex][partIndex]->SetBinError(m_currentRapBin, fitResult->GetParError(1));
            m_tpcResolution_Minus[partBaseIndex][partIndex]->SetBinContent(m_currentRapBin, fitResult->GetParameter(2));
            m_tpcResolution_Minus[partBaseIndex][partIndex]->SetBinError(m_currentRapBin, fitResult->GetParError(2));
          }


          fittingCanvas->cd();
          m_currentHistoToFit->Draw("E");
          gPad->SetRightMargin(0.4);
          gPad->SetLogy();
          gPad->Update();
          //fitResult->SetNpx(500);
          //fitResult->Draw("same");
          gSystem->ProcessEvents();
          initFunct->Draw("SAME");
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/%s_RapBin_%02d_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex,m_currentPartCharge).Data(),m_currentRapBin,m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));

          gPad->SetLogy(false);
          gPad->Update();
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/NoLog_%s_RapBin_%02d_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex,m_currentPartCharge).Data(),m_currentRapBin,m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));
          cout << "  deleing stuff: " << initFunct << "  " << m_currentHistoToFit << endl;
          #ifndef _SAVE_POINTERS_
            delete initFunct;
            delete m_currentHistoToFit;
          #endif
        }// end rap bin

        cout << "  drawing stuff" << endl;
        if(m_currentPartCharge > 0){  //positive particles
          fittingCanvas->cd();
          gPad->SetLogy(false);
          m_bichselShifts_Plus[partBaseIndex][partIndex]->Draw("E");
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/BichselShift_%s_space_AllRap_for_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex, m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));

          fittingCanvas->cd();
          m_quickYieldNums_Plus[partBaseIndex][partIndex]->Draw("E");
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/QuickYield_%s_space_AllRap_for_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex, m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));

          fittingCanvas->cd();
          m_tpcResolution_Plus[partBaseIndex][partIndex]->Draw("E");
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/TPCResolution_%s_space_AllRap_for_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex, m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));


        }else{  //negative particles

          fittingCanvas->cd();
          gPad->SetLogy(false);
          m_bichselShifts_Minus[partBaseIndex][partIndex]->Draw("E");
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/BichselShift_%s_space_AllRap_for_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex, m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));

          fittingCanvas->cd();
          m_quickYieldNums_Minus[partBaseIndex][partIndex]->Draw("E");
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/QuickYield_%s_space_AllRap_for_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex, m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));

          fittingCanvas->cd();
          m_tpcResolution_Minus[partBaseIndex][partIndex]->Draw("E");
          fittingCanvas->SaveAs(Form("%s/%s/BichselShifts/TPCResolution_%s_space_AllRap_for_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partBaseIndex, m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data()));
        }



      }// end plus minus
    }// end inner part
  }// end outer part

  cout << " getting kaon ratio" << endl;

  m_basicKaonPionRatio = m_quickYieldNums_Plus[0][1]->GetBinContent(m_quickYieldNums_Plus[0][1]->FindBin(m_cuts->getYCM()))
                         /m_quickYieldNums_Plus[0][0]->GetBinContent(m_quickYieldNums_Plus[0][0]->FindBin(m_cuts->getYCM()));


  m_bichselShiftsComplete = true;
}


void ZFitter::findMassSquaredWidths(){

  if(!m_massSquaredMode){
    cout << "Call to ZFitter::findMassSquaredWidths, but you are not in mass squared mode... returning without doing anything" << endl;
    return;
  }

  cout << "starting ZFitter::findMassSquaredWidths" << endl;

  int yCM_bin = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetXaxis()->FindBin(m_cuts->getYCM());
  cout << "yCM_bin: " << yCM_bin << endl;
  int lowP_bin = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->FindBin(0.45);
  int midP_bin = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->FindBin(0.75);
  int highP_bin1 = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->FindBin(1.2);
  int highP_bin2 = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->FindBin(1.4);

  m_mass2WidthsMomentum[0] = PhysMath::pTot_from_y_mTm0(m_cuts->getYCM(), m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->GetBinCenter(lowP_bin), m_partInfo->GetParticleMass(0));
  m_mass2WidthsMomentum[1] = PhysMath::pTot_from_y_mTm0(m_cuts->getYCM(), m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->GetBinCenter(midP_bin), m_partInfo->GetParticleMass(0));
  m_mass2WidthsMomentum[2] = PhysMath::pTot_from_y_mTm0(m_cuts->getYCM(), (m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->GetBinLowEdge(highP_bin1)+m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->GetYaxis()->GetBinLowEdge(highP_bin2 + 1))/2.0, m_partInfo->GetParticleMass(0));

  cout << "lowP: " << m_mass2WidthsMomentum[0] << endl;
  cout << "midP: " << m_mass2WidthsMomentum[1] << endl;
  cout << "highP: " << m_mass2WidthsMomentum[2] << endl;

  TH1D *lowP_hist = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->ProjectionZ("BTOF_mass2_width_lowHist",yCM_bin, yCM_bin, lowP_bin, lowP_bin);
  TH1D *midP_hist = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->ProjectionZ("BTOF_mass2_width_midHist",yCM_bin, yCM_bin, midP_bin, midP_bin);
  TH1D *highP_hist = m_rap_mTm0_invBetaBTOF_Plus_allCent[0]->ProjectionZ("BTOF_mass2_width_highHist",yCM_bin, yCM_bin, highP_bin1, highP_bin2);

  lowP_hist->SetMarkerStyle(20);
  midP_hist->SetMarkerStyle(20);
  highP_hist->SetMarkerStyle(20);

  TCanvas *widthFittingCanvas = new TCanvas();


  for(int partIndex=0; partIndex<9; partIndex++){
    TF1 *partGauss = new TF1("partGauss","[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",pow(m_partInfo->GetParticleMass(partIndex)-0.1,2.0),pow(m_partInfo->GetParticleMass(partIndex)+0.1,2.0));
    cout << "function limits [" << pow(m_partInfo->GetParticleMass(partIndex)-0.05,2.0) << "," << pow(m_partInfo->GetParticleMass(partIndex)+0.05,2.0) <<"]" << endl;
    double amp = 1;
    double mean = pow(m_partInfo->GetParticleMass(partIndex),2.0);
    double sig = 0.05;

    double x1 = lowP_hist->GetBinCenter(lowP_hist->FindBin(mean - 0.25*sig));
    double x2 = lowP_hist->GetBinCenter(lowP_hist->FindBin(mean));
    double x3 = lowP_hist->GetBinCenter(lowP_hist->FindBin(mean + 0.25*sig));
    if(x1 == x2) x1 = lowP_hist->GetBinCenter(lowP_hist->FindBin(mean) - 1);
    if(x3 == x2) x3 = lowP_hist->GetBinCenter(lowP_hist->FindBin(mean) + 1);
    cout << "x1: " << x1 << endl;
    cout << "x2: " << x2 << endl;
    cout << "x3: " << x3 << endl;

    double y1 = HistogramUtilities::getBinContentAvg(lowP_hist, lowP_hist->FindBin(x1), 2);
    cout << "y1: " << y1 << endl;
    double y2 = HistogramUtilities::getBinContentAvg(lowP_hist, lowP_hist->FindBin(x2), 2);
    cout << "y2: " << y2 << endl;
    double y3 = HistogramUtilities::getBinContentAvg(lowP_hist, lowP_hist->FindBin(x3), 2);
    cout << "y3: " << y3 << endl;

    PhysMath::get3PointGausParams(amp, mean, sig,x1,y1,x2,y2,x3,y3);

    cout << "ampSeed: " << amp << endl;
    cout << "mean: " << mean << endl;
    cout << "sigmaSeed: " << sig << endl;

    partGauss->SetParameters(amp,mean,sig);
    partGauss->SetParLimits(1,mean-0.05,mean+0.05);
    partGauss->SetParLimits(2,0.1*sig,2.0*sig);

    lowP_hist->Fit(partGauss,"E");
    widthFittingCanvas->cd();
    widthFittingCanvas->SetLogy();
    lowP_hist->Draw("AP");
    partGauss->Draw("same");
    widthFittingCanvas->SaveAs(Form("%s/%s/ParticleWidths/%s_width_lowP.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partIndex, 1).Data()));
    widthFittingCanvas->Clear();
    m_mass2Widths[partIndex][0] = partGauss->GetParameter(2);
    midP_hist->Fit(partGauss,"E");
    widthFittingCanvas->SetLogy();
    midP_hist->Draw("AP");
    partGauss->Draw("same");
    widthFittingCanvas->SaveAs(Form("%s/%s/ParticleWidths/%s_width_midP.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partIndex, 1).Data()));
    widthFittingCanvas->Clear();
    m_mass2Widths[partIndex][1] = partGauss->GetParameter(2);
    highP_hist->Fit(partGauss,"E");
    widthFittingCanvas->SetLogy();
    highP_hist->Draw("AP");
    partGauss->Draw("same");
    widthFittingCanvas->SaveAs(Form("%s/%s/ParticleWidths/%s_width_highP.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
              m_partInfo->GetParticleName(partIndex, 1).Data()));
    widthFittingCanvas->Clear();
    m_mass2Widths[partIndex][2] = partGauss->GetParameter(2);
    #ifndef _SAVE_POINTERS_
      delete partGauss;
    #endif

  }

  return;

}

void ZFitter::getRangeZTOF_mass2(double &a_lowZTOF, double &a_highZTOF, int a_partIndex){

  if(!m_massSquaredMode) return;

  cout << "getting mass2 range for partIndex " << a_partIndex << " with a total momentum of " << m_currentMomOfBin << endl;

  //double mu = ((m_mass2Widths[a_partIndex][2]-m_mass2Widths[a_partIndex][0])/(m_mass2Widths[a_partIndex][1]-m_mass2Widths[a_partIndex][0])*(pow(m_mass2WidthsMomentum[1],2.0)-pow(m_mass2WidthsMomentum[0],2.0))+pow(m_mass2WidthsMomentum[0],2.0)-pow(m_mass2WidthsMomentum[2],2.0))/(2.0*((m_mass2Widths[a_partIndex][2]-m_mass2Widths[a_partIndex][0])/(m_mass2Widths[a_partIndex][1]-m_mass2Widths[a_partIndex][0])*(m_mass2WidthsMomentum[1]-m_mass2WidthsMomentum[0])-m_mass2WidthsMomentum[2] + m_mass2WidthsMomentum[0]));

  //cout << "mu: " << mu << endl;

  //double A = (m_mass2Widths[a_partIndex][1]-m_mass2Widths[a_partIndex][0])/(pow(m_mass2WidthsMomentum[1]-mu,2.0)-pow(m_mass2WidthsMomentum[0]-mu,2.0));

  //cout << "A: " << A << endl;
  //double Const = m_mass2Widths[a_partIndex][0] - A*pow(m_mass2WidthsMomentum[0]-mu,2.0);

  //cout << "Const: " << Const << endl;

//  double quadratic = A*pow(m_currentMomOfBin - mu,2.0) + Const;
  double quadratic = m_mass2Widths[a_partIndex][0]*(m_currentMomOfBin-m_mass2WidthsMomentum[1])*(m_currentMomOfBin-m_mass2WidthsMomentum[2])/((m_mass2WidthsMomentum[0]-m_mass2WidthsMomentum[1])*(m_mass2WidthsMomentum[0]-m_mass2WidthsMomentum[2]))
    +m_mass2Widths[a_partIndex][1]*(m_currentMomOfBin-m_mass2WidthsMomentum[0])*(m_currentMomOfBin-m_mass2WidthsMomentum[2])/((m_mass2WidthsMomentum[1]-m_mass2WidthsMomentum[0])*(m_mass2WidthsMomentum[1]-m_mass2WidthsMomentum[2]))
    +m_mass2Widths[a_partIndex][2]*(m_currentMomOfBin-m_mass2WidthsMomentum[0])*(m_currentMomOfBin-m_mass2WidthsMomentum[1])/((m_mass2WidthsMomentum[2]-m_mass2WidthsMomentum[0])*(m_mass2WidthsMomentum[2]-m_mass2WidthsMomentum[1]));
  cout << "quadratic: " << quadratic << endl;

  double quadraticDown = -999;
  if(a_partIndex > 0) quadraticDown = m_mass2Widths[a_partIndex-1][0]*(m_currentMomOfBin-m_mass2WidthsMomentum[1])*(m_currentMomOfBin-m_mass2WidthsMomentum[2])/((m_mass2WidthsMomentum[0]-m_mass2WidthsMomentum[1])*(m_mass2WidthsMomentum[0]-m_mass2WidthsMomentum[2]))
    +m_mass2Widths[a_partIndex-1][1]*(m_currentMomOfBin-m_mass2WidthsMomentum[0])*(m_currentMomOfBin-m_mass2WidthsMomentum[2])/((m_mass2WidthsMomentum[1]-m_mass2WidthsMomentum[0])*(m_mass2WidthsMomentum[1]-m_mass2WidthsMomentum[2]))
    +m_mass2Widths[a_partIndex-1][2]*(m_currentMomOfBin-m_mass2WidthsMomentum[0])*(m_currentMomOfBin-m_mass2WidthsMomentum[1])/((m_mass2WidthsMomentum[2]-m_mass2WidthsMomentum[0])*(m_mass2WidthsMomentum[2]-m_mass2WidthsMomentum[1]));
  double quadraticUp = m_mass2Widths[a_partIndex+1][0]*(m_currentMomOfBin-m_mass2WidthsMomentum[1])*(m_currentMomOfBin-m_mass2WidthsMomentum[2])/((m_mass2WidthsMomentum[0]-m_mass2WidthsMomentum[1])*(m_mass2WidthsMomentum[0]-m_mass2WidthsMomentum[2]))
    +m_mass2Widths[a_partIndex+1][1]*(m_currentMomOfBin-m_mass2WidthsMomentum[0])*(m_currentMomOfBin-m_mass2WidthsMomentum[2])/((m_mass2WidthsMomentum[1]-m_mass2WidthsMomentum[0])*(m_mass2WidthsMomentum[1]-m_mass2WidthsMomentum[2]))
    +m_mass2Widths[a_partIndex+1][2]*(m_currentMomOfBin-m_mass2WidthsMomentum[0])*(m_currentMomOfBin-m_mass2WidthsMomentum[1])/((m_mass2WidthsMomentum[2]-m_mass2WidthsMomentum[0])*(m_mass2WidthsMomentum[2]-m_mass2WidthsMomentum[1]));

  a_lowZTOF = pow(m_partInfo->GetParticleMass(a_partIndex),2.0) - 3*quadratic;
  a_highZTOF = pow(m_partInfo->GetParticleMass(a_partIndex),2.0) + 5*quadratic;

  if(a_partIndex == 0){
    if(a_highZTOF > pow(m_partInfo->GetParticleMass(1),2.0) - 3*quadraticUp) a_highZTOF = pow(m_partInfo->GetParticleMass(1),2.0) - 3*quadraticUp;
  }else if(a_partIndex == 1){
    if(a_highZTOF > pow(m_partInfo->GetParticleMass(2),2.0) - 3*quadraticUp) a_highZTOF = pow(m_partInfo->GetParticleMass(2),2.0) - 3*quadraticUp;
    if(a_lowZTOF < pow(m_partInfo->GetParticleMass(0),2.0) + 3*quadraticDown) a_lowZTOF = pow(m_partInfo->GetParticleMass(0),2.0) + 3*quadraticDown;
  }else if(a_partIndex == 2){
    if(a_highZTOF > pow(m_partInfo->GetParticleMass(4),2.0) - 3*quadraticUp) a_highZTOF = pow(m_partInfo->GetParticleMass(4),2.0) - 3*quadraticUp;
    if(a_lowZTOF < pow(m_partInfo->GetParticleMass(1),2.0) + 3*quadraticDown) a_lowZTOF = pow(m_partInfo->GetParticleMass(1),2.0) + 3*quadraticDown;
  }else if(a_partIndex == 4){
    if(a_highZTOF > pow(m_partInfo->GetParticleMass(5),2.0) - 3*quadraticUp) a_highZTOF = pow(m_partInfo->GetParticleMass(5),2.0) - 3*quadraticUp;
    if(a_lowZTOF < pow(m_partInfo->GetParticleMass(2),2.0) + 3*quadraticDown) a_lowZTOF = pow(m_partInfo->GetParticleMass(2),2.0) + 3*quadraticDown;
  }else if(a_partIndex == 5){
    if(a_highZTOF > pow(m_partInfo->GetParticleMass(6),2.0) - 3*quadraticUp) a_highZTOF = pow(m_partInfo->GetParticleMass(6),2.0) - 3*quadraticUp;
    if(a_lowZTOF < pow(m_partInfo->GetParticleMass(5),2.0) + 3*quadraticDown) a_lowZTOF = pow(m_partInfo->GetParticleMass(5),2.0) + 3*quadraticDown;
  }


  cout << "a_lowZTOF: " << a_lowZTOF << endl;
  cout << "a_highZTOF: " << a_highZTOF << endl;

  return;

}




void ZFitter::makeElectronExtrapolationFunction(){
  #ifdef _ZFITTER_DEBUG_
    cout << " Call to ZFitter::makeElectronExtrapolationFunction()" << endl;
  #endif
  TCanvas* fittingCanvas = new TCanvas();

  TF1* elecExtrap = new TF1("elecExtrap","[0]*exp(-[1]*x+[2])",0.01,5.0);


  TH1D* proj = NULL;
  if(m_currentPartCharge == 1){
    proj = m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->ProjectionY(Form("proj_plus_%02d_%02d",m_currentPartIndex,m_currentRapBin),
      m_currentRapBin,m_currentRapBin);
  }else if(m_currentPartCharge == -1){
    proj = m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->ProjectionY(Form("proj_plus_%02d_%02d",m_currentPartIndex,m_currentRapBin),
      m_currentRapBin,m_currentRapBin);
  }

  if(proj->Integral() <= 0.0){
    cout << " electron amp proj for extrapolation is empty, cant work" << endl;
    return;
  }

  proj->Scale(1.0,"width");
  for(int binX = 1; binX < proj->GetNbinsX(); binX++){
    double center = proj->GetBinCenter(binX);
    double content = proj->GetBinContent(binX);
    double error =  proj->GetBinError(binX);
    double newContent = content/center;
    double newError = error/center;
  }

  proj->GetYaxis()->SetTitle("1/(mt-m0) dN/dmtm0");


  double rapidity = m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->GetXaxis()->GetBinCenter(m_currentRapBin);

  int bin1 = proj->FindBin(PhysMath::mTm0_from_y_p(rapidity, m_electronExtrapolation_LowMom, m_partInfo->GetParticleMass(m_currentPartIndex)));
  while(proj->GetBinContent(bin1) <= 0 ) bin1++;
  double x1 = proj->GetBinCenter(bin1);
  double y1 = proj->GetBinContent(bin1);
  int bin2 = proj->FindBin(PhysMath::mTm0_from_y_p(rapidity, m_electronExtrapolation_HighMom, m_partInfo->GetParticleMass(m_currentPartIndex))) - 1;
  double x2 = proj->GetBinCenter(bin2);
  double y2 = proj->GetBinContent(bin2);
  if(y2 <= 0.0){
    cout << "high mom elec bin is empty... ERROR" << endl;
  }
  double b = log(y2/y1)/(x1-x2);
  double a = y1*exp(b*x1)/b;

  elecExtrap->SetParameter(0,a);
  elecExtrap->SetParameter(1,b);
  elecExtrap->SetParameter(2,0.0);
  elecExtrap->SetParLimits(0,0.2*a,5.0*a);
  elecExtrap->SetParLimits(1,0.1*b,10*b);
  elecExtrap->SetParLimits(2,-5,5);

  cout << "bin1: " << bin1 << "  x1:" << x1 << " y1:" << y1 << endl;
  cout << "bin2: " << bin2 << "  x2:" << x2 << " y2:" << y2 << endl;
  cout << "amp: " << a << "slope" << b << endl;

  int status = proj->Fit("elecExtrap","E","",x1,x2);
  if(status != 0){
    cout << "BAD ELECTRON EXTRAPOLATION" << endl;
    //break the code
    //cout << m_dEdxFitDataMinus[100][1][1]->GetName() << endl;
    //return;
  }
  TF1* fitResult = proj->GetFunction("elecExtrap");
  if(fitResult){
    fitResult->SetNpx(500);
    fitResult->SetRange(x1,10.0);
  }


  double biggestBinCont = 0.0;
  double smallestBinCont = 1E60;
  for(int binX = 1; binX < proj->GetNbinsX(); binX++){
    double content = proj->GetBinContent(binX);
    if(content > biggestBinCont) biggestBinCont = content;
    if(content < smallestBinCont && content > 0.0) smallestBinCont = content;
  }


  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  proj->GetXaxis()->SetRangeUser(0.0,x2+0.1);
  proj->GetYaxis()->SetRangeUser(smallestBinCont,biggestBinCont);
  gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  proj->SetStats(true);

  proj->Draw("E");
  if(fitResult) fitResult->Draw("same");

  fittingCanvas->SaveAs(Form("%s/%s/electronStudy/electron_extrapolation_%s_Rap_%02d.png",
   m_imagePreDir.c_str(),      m_imgDirName.c_str(),
   m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),     m_currentRapBin));

  if(fitResult){
    if(m_currentPartCharge == 1){
      for(int vertBin = bin2; vertBin <= m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->GetNbinsY(); vertBin++){
        double vertCenter = m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->GetYaxis()->GetBinCenter(vertBin);
        double vertWidth = m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->GetYaxis()->GetBinWidth(vertBin);
        double prediction = fitResult->Eval(vertCenter)*vertWidth*vertCenter;
        m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->SetBinContent(m_currentRapBin,vertBin, prediction);
      }
    }else if(m_currentPartCharge == -1){
      for(int vertBin = bin2; vertBin <= m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->GetNbinsY(); vertBin++){
        double vertCenter = m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->GetYaxis()->GetBinCenter(vertBin);
        double vertWidth =  m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->GetYaxis()->GetBinWidth(vertBin);
        double prediction = fitResult->Eval(vertCenter)*vertWidth*vertCenter;
        m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->SetBinContent(m_currentRapBin,vertBin, prediction);
      }
    }
    m_electronExtrapolationComplete = true;


  }
  #ifndef _SAVE_POINTERS_
    delete proj;
  #endif

  #ifdef _ZFITTER_DEBUG_
    cout << " End Call to ZFitter::makeElectronExtrapolationFunction()" << endl;
  #endif


  return;
}






void ZFitter::drawBichselOnDeDxPlot(){
  TCanvas* fittingCanvas = new TCanvas();
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  fittingCanvas->Divide(2,1);
  fittingCanvas->cd(1);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_dEdx_Plus->Draw("colz");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex == 8) continue;
    m_bichselCurves[partIndex]->Draw("same");
  }

  fittingCanvas->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_dEdx_Minus->Draw("colz");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex == 8) continue;
    m_bichselCurves[partIndex]->Draw("same");
  }

   fittingCanvas->SaveAs(Form("%s/%s/bichsel_plot.png", m_imagePreDir.c_str(), m_imgDirName.c_str()));

}

void ZFitter::drawBichselOnDeDxPlot_Liggett(){
  TCanvas* fittingCanvas = new TCanvas();
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  //fittingCanvas->Divide(2,1);
  fittingCanvas->cd();
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_dEdx_Plus->Draw("colz");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex == 8) continue;
    m_bichselCurves[partIndex]->Draw("same");
  }

  /*fittingCanvas->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_dEdx_Minus->Draw("colz");
  for(int partIndex = 4; partIndex < 9; partIndex++){
    if(partIndex == 8) continue;
    m_bichselCurves[partIndex]->Draw("same");
  }
*/
   fittingCanvas->SaveAs(Form("%s/%s/bichsel_plot.png", m_imagePreDir.c_str(), m_imgDirName.c_str()));

}

void ZFitter::drawInvBetaPredictionPlot(){
  TCanvas* fittingCanvas = new TCanvas();
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  fittingCanvas->Divide(2,1);
  fittingCanvas->cd(1);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_invBetaBTOF_Plus->Draw("colz");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    //if(partIndex == 8) continue;
    m_invBetaCurves[partIndex]->Draw("same");
  }

  fittingCanvas->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_invBetaBTOF_Minus->Draw("colz");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    //if(partIndex == 8) continue;
    m_invBetaCurves[partIndex]->Draw("same");
  }

   fittingCanvas->SaveAs(Form("%s/%s/invBetaCalib_plot.png", m_imagePreDir.c_str(), m_imgDirName.c_str()));

}

void ZFitter::drawInvBetaPredictionPlot_Liggett(){
  TCanvas* fittingCanvas = new TCanvas();
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  //fittingCanvas->Divide(2,1);
  fittingCanvas->cd();
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_invBetaBTOF_Plus->Draw("colz");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex == 8) continue;
    cout << "Print line for particle: " << partIndex << endl;
    m_invBetaCurves[partIndex]->Draw("same");
  }

  /*fittingCanvas->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetLogy();
  m_mom_invBetaBTOF_Minus->Draw("colz");
  for(int partIndex = 4; partIndex < 9; partIndex++){
    //if(partIndex == 8) continue;
    m_invBetaCurves[partIndex]->Draw("same");
    }*/

   fittingCanvas->SaveAs(Form("%s/%s/invBetaCalib_plot.png", m_imagePreDir.c_str(), m_imgDirName.c_str()));

}



void ZFitter::make_TPC_Histo2Ds(int a_detIndex){
  // -1 for all detectors, 0 = TPC, 1 = BTOF, 2 = ETOF
  for(int detIndex = 0; detIndex < 3; detIndex++){
    if(a_detIndex != -1 && a_detIndex != detIndex) continue;
    for(int partIndex = 0; partIndex < 9; partIndex++){
      if(partIndex == 3 && partIndex == 8) continue;
      for(int chargeIndex = 0; chargeIndex < 2; chargeIndex++){
        int charge = 1;
        if(chargeIndex == 1) charge = -1;

        int num_rap_bins = -1;
        vector<double> rap_bin_edges;
        vector<double> mtm0_bin_edges;
        string det_var_name = "";
        bool valid_read = false;
        if(detIndex == 0){
          det_var_name = "ZTPC";
          if(chargeIndex == 0){
            if(!m_rap_mTm0_dEdx_Plus[partIndex][0]) continue;
            num_rap_bins = m_rap_mTm0_dEdx_Plus[partIndex][0]->GetNbinsX();
            rap_bin_edges = HistogramUtilities::getBinEdges_X(m_rap_mTm0_dEdx_Plus[partIndex][0]);
            mtm0_bin_edges = HistogramUtilities::getBinEdges_Y(m_rap_mTm0_dEdx_Plus[partIndex][0]);
          }else{
            if(!m_rap_mTm0_dEdx_Minus[partIndex][0]) continue;
            num_rap_bins = m_rap_mTm0_dEdx_Minus[partIndex][0]->GetNbinsX();
            rap_bin_edges = HistogramUtilities::getBinEdges_X(m_rap_mTm0_dEdx_Minus[partIndex][0]);
            mtm0_bin_edges = HistogramUtilities::getBinEdges_Y(m_rap_mTm0_dEdx_Minus[partIndex][0]);
          }
        }else if(detIndex == 1){
          det_var_name = "ZBTOF";
          if(chargeIndex == 0){
            if(!m_rap_mTm0_invBetaBTOF_Plus[partIndex][0]) continue;
            num_rap_bins = m_rap_mTm0_invBetaBTOF_Plus[partIndex][0]->GetNbinsX();
            rap_bin_edges = HistogramUtilities::getBinEdges_X(m_rap_mTm0_invBetaBTOF_Plus[partIndex][0]);
            mtm0_bin_edges = HistogramUtilities::getBinEdges_Y(m_rap_mTm0_invBetaBTOF_Plus[partIndex][0]);
          }else{
            if(!m_rap_mTm0_invBetaBTOF_Minus[partIndex][0]) continue;
            num_rap_bins = m_rap_mTm0_invBetaBTOF_Minus[partIndex][0]->GetNbinsX();
            rap_bin_edges = HistogramUtilities::getBinEdges_X(m_rap_mTm0_invBetaBTOF_Minus[partIndex][0]);
            mtm0_bin_edges = HistogramUtilities::getBinEdges_Y(m_rap_mTm0_invBetaBTOF_Minus[partIndex][0]);
          }
        }else if(detIndex == 2){
          det_var_name = "ZETOF";
          if(chargeIndex == 0){
            if(!m_rap_mTm0_invBetaETOF_Plus[partIndex][0]) continue;
            num_rap_bins = m_rap_mTm0_invBetaETOF_Plus[partIndex][0]->GetNbinsX();
            rap_bin_edges = HistogramUtilities::getBinEdges_X(m_rap_mTm0_invBetaETOF_Plus[partIndex][0]);
            mtm0_bin_edges = HistogramUtilities::getBinEdges_Y(m_rap_mTm0_invBetaETOF_Plus[partIndex][0]);
          }else{
            if(!m_rap_mTm0_invBetaETOF_Minus[partIndex][0]) continue;
            num_rap_bins = m_rap_mTm0_invBetaETOF_Minus[partIndex][0]->GetNbinsX();
            rap_bin_edges = HistogramUtilities::getBinEdges_X(m_rap_mTm0_invBetaETOF_Minus[partIndex][0]);
            mtm0_bin_edges = HistogramUtilities::getBinEdges_Y(m_rap_mTm0_invBetaETOF_Minus[partIndex][0]);
          }
        }
        for(int cent_index = 0; cent_index < 9; cent_index++){
          for(int rap_index = 0; rap_index < num_rap_bins; rap_index++){
            //Histo2D* hist = new Histo2D("name", "title", ((int)rap_bin_edges.size())-1, &rap_bin_edges[0]);
            string hist_name = Form("mTm0_%s_%s_Cent%02d_rapIndex%03d", det_var_name.c_str(), m_partInfo->GetParticleName(partIndex,charge).Data(), cent_index, rap_index);
            string hist_title = Form("%s %s Cent %02d y=[%1.2f,%1.2f]", det_var_name.c_str(), m_partInfo->GetParticleName(partIndex,charge).Data(), cent_index, rap_bin_edges[rap_index], rap_bin_edges[rap_index+1]);
            Histo2D* hist = new Histo2D(hist_name, hist_title, ((int)mtm0_bin_edges.size())-1, &mtm0_bin_edges[0]); // this is a weak construction (adds TH1D projections as we go)
            if(detIndex == 0){
              m_mTm0_ZTPC_by_Rapidity[partIndex][cent_index][rap_index][chargeIndex] = hist;
            }
            if(detIndex == 1){
              m_mTm0_ZBTOF_by_Rapidity[partIndex][cent_index][rap_index][chargeIndex] = hist;
            }
            if(detIndex == 2){
              m_mTm0_ZETOF_by_Rapidity[partIndex][cent_index][rap_index][chargeIndex] = hist;
            }
          } // rap_index
        } // cent_index
      } // chargeIndex
    } // partIndex
  } // detIndex
}







#endif
