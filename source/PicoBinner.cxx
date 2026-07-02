// =====================================================================================
// PicoBinner.cxx -- ported from lightflavorspectra_etof for LOCAL, 200 GeV O+O ONLY use.
//
// WHAT CHANGED FROM THE ORIGINAL (read this before trusting any output):
//
//   1. PicoDstReader version: the original supported SL22b/SL22c/SL23a/SL23d/SL19b
//      readers via nested #ifdef. _OO_200_COL_ always selected SL23a, so that's the
//      only reader included here -- all other reader includes were deleted.
//   2. MattDST mode (_MATTDST_MODE_) was never defined for this analysis and has been
//      deleted; this file always reads PicoDsts.
//   3. FixedTarget branches deleted throughout. PicoBinner itself enforces
//      eventConfig=="ColliderCenter" for O+O (it exits otherwise), so every
//      eventConfig=="FixedTarget" branch in the original was dead code for this
//      analysis and has been removed.
//   4. _STREFMULTCORR_ was NOT defined for _OO_200_COL_ in the original
//      makefile_toggles.h -- centrality here always uses raw refMult, never
//      StRefMultCorr. The StRefMultCorr setup/branches have been removed; only the
//      raw-refMult path survives.
//   5. _JIA_CHEN_19p6_Modification_, _PICOBINNER_VZ_STUDY_, and
//      _PICOBINNER_2D_DEDX_BTOF_ were not defined and have been deleted outright.
//   6. _USE_CALIB_CLASS_ was not defined; only the _PICOBINNER_USE_BICHSEL_ path
//      (which WAS defined) for computing zTPC survives.
//   7. _PICO_BINNER_DEBUG_ verbose cout statements (never defined/compiled) removed
//      for readability.
//   8. The original called picoReader->Finish() and `delete picoReader` once after
//      the event loop, then had a SECOND Finish() call guarded by
//      `#ifndef _PICO_READER_FINISH_BYPASS_` at the very end of the function. Since
//      _PICO_READER_FINISH_BYPASS_ WAS defined, that second call never ran (and would
//      have been a use-after-free on the already-deleted picoReader if it had run).
//      That dead/dangerous second call has been deleted; only the first, correct
//      Finish()+delete survives.
//
// WHAT DID NOT CHANGE (deliberately preserved exactly, even though this build only
// analyzes pi/K/p):
//
//   - nParticles is still 9 and every "isolated"/contamination-modeling histogram
//     (ZTPCHistoIsolatedPlus/Minus, the zTPC_byPart_hist[9][2] calibration array, the
//     inner partIndex2 loops, etc.) still spans all 9 species. This machinery models
//     dE/dx and TOF band overlap from bystander species (electrons, muons, deuterons,
//     etc.) into the pion/kaon/proton bands. It was not verified safe to remove, so it
//     was left untouched to protect "must match exactly" fidelity for pi/K/p. Only the
//     DRIVER (which species get a full/isolated per-run analysis, via a_partIndex) is
//     ever pointed at pion(0)/kaon(1)/proton(2) -- exactly how the original was already
//     used in practice.
//   - All bin-edge/histogram-content/fit-input numeric logic is unchanged.
// =====================================================================================

/// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
/// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TGraph.h"

#include "../makefile_toggles.h"

// PicoDst headers -- SL23a only (the reader version _OO_200_COL_ always selected)
#include "../submodule/PicoDstReader_SL23a/StPicoDstReader.h"
#include "../submodule/PicoDstReader_SL23a/StPicoDst.h"
#include "../submodule/PicoDstReader_SL23a/StPicoEvent.h"
#include "../submodule/PicoDstReader_SL23a/StPicoTrack.h"
#include "../submodule/PicoDstReader_SL23a/StPicoBTofHit.h"
#include "../submodule/PicoDstReader_SL23a/StPicoBTowHit.h"
#include "../submodule/PicoDstReader_SL23a/StPicoEmcTrigger.h"
#include "../submodule/PicoDstReader_SL23a/StPicoBTofPidTraits.h"
#include "../submodule/PicoDstReader_SL23a/StPicoETofPidTraits.h"
#include "../submodule/PicoDstReader_SL23a/StPicoTrackCovMatrix.h"
#include "../submodule/PicoDstReader_SL23a/StPicoEpdHit.h"
#include "../submodule/PicoDstReader_SL23a/StPicoBbcHit.h"
#include "../submodule/PicoDstReader_SL23a/StPicoETofHit.h"

#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"
#include "../headers/HistogramUtilities.h"
#include "../headers/PhysMath.h"
#include "../headers/CutClass.h"
#include "../headers/CalibrationClass.h"
using namespace std;

double tofPathLength(TVector3 beginPoint, TVector3 endPoint, double curvature){
  double xdif = endPoint.X() - beginPoint.X();
  double ydif = endPoint.Y() - beginPoint.Y();
  double C = sqrt(xdif*xdif + ydif*ydif);
  double s_perp = C;
  if (curvature){
    double R = 1/curvature;
    if(C/(2*R) > 1 || C/(2*R)<-1) return -999.;
    s_perp = 2*R * asin(C/(2*R));
  }
  double s_z = fabs(endPoint.Z() - beginPoint.Z());
  return sqrt(s_perp*s_perp + s_z*s_z);
}

double lowestMass(int partBaseIndex, int a_charge, CutClass* a_cutClass, ParticleInfo *a_partInfo){
  double currentLowMass = 100.0;
  bool particleFitMatrix[9][9];
  if(a_charge>0) a_cutClass->getParticleFitMatrix(1,particleFitMatrix);
  else if(a_charge<0) a_cutClass->getParticleFitMatrix(-1,particleFitMatrix);
  for(int partIndex=0; partIndex<9; partIndex++){
    if(particleFitMatrix[partBaseIndex][partIndex]){
      if(a_partInfo->GetParticleMass(partIndex)<currentLowMass) currentLowMass = a_partInfo->GetParticleMass(partIndex);
    }
  }
  return currentLowMass;
}

double highestMass(int partBaseIndex, int a_charge, CutClass* a_cutClass, ParticleInfo *a_partInfo){
  double currentHighMass = 0.0;
  bool particleFitMatrix[9][9];
  if(a_charge>0) a_cutClass->getParticleFitMatrix(1,particleFitMatrix);
  else if (a_charge<0) a_cutClass->getParticleFitMatrix(-1,particleFitMatrix);
  for(int partIndex=0; partIndex<9; partIndex++){
    if(particleFitMatrix[partBaseIndex][partIndex]){
      if(a_partInfo->GetParticleMass(partIndex)>currentHighMass) currentHighMass = a_partInfo->GetParticleMass(partIndex);
    }
  }
  return currentHighMass;
}

//###########################################################################################################################################################
//###########################################################################################################################################################
void PicoBinner(string    a_filelist,
                string    a_outFileName,
                int       a_numParticles,   // pass 9 -- see header comment on why this must not be reduced
                CutClass* a_cutClass,
                int       a_partIndex       = -1,   // 0=pion, 1=kaon, 2=proton for this build
                bool      a_makeBasicHistos = true,
                bool      a_correctT0       = false,
                bool      a_ignoreNonParticleOfInterestByBTOF = false,
                bool      a_doMassSquared   = false,
                CalibrationClass* a_calib = NULL){

  const bool doProtonDCA   = true;
  const bool useHuiPileup  = a_cutClass->getUseHuiDT27GevPileup();
  const bool doMassSquaredVariableBinning = true; // if do mass squared make it variable binning
  if(useHuiPileup) for(int iii = 0; iii < 50; iii++) cout << "\n USING HUI PILEUP!!!!!! CHECK THE TOGGLES AT THE TOP OF PICO BINNER!!!!!!!!!! \n" <<endl;

  TString eventConfig  = a_cutClass->getEventConfig();
  Float_t yCM          = a_cutClass->getYCM();
  bool useLogMTBinning = a_cutClass->getLogMTBinning();

  if(eventConfig != "ColliderCenter"){
    cout << "PicoBinner :: Error - this build is O+O 200 GeV only and requires eventConfig==\"ColliderCenter\" from CutClass." << endl;
    exit(EXIT_FAILURE);
  }
  if(yCM != 0.0){
    cout << "PicoBinner :: Error - working in Collider mode but yCM!=0" << endl;
    exit(EXIT_FAILURE);
  }

  string filelistPath   = a_filelist;
  string outputFileName = a_outFileName;
  TString starver       = a_cutClass->getStarver();

  int nParticles  = a_numParticles;
  int numCentBins = a_cutClass->getNCentBins();

  //######## LOAD PARTICLE INFO for Bichsel Curves and InvBeta Curves
  ParticleInfo* particleInfo= new ParticleInfo(starver,true,1); //do not change the last two arguments unless you know what you are doing
  TF1* bichselCurves[9];
  TF1* invBetaCurves[9];
  for(int partIndex = 0; partIndex < 9; partIndex++){
    bichselCurves[partIndex] = particleInfo->GetBichselFunction(partIndex);
    invBetaCurves[partIndex] = particleInfo->GetInverseBetaFunction(partIndex);
  }

  //###################################################
  //###         SET UP SIMPLE HISTOGRAMS            ###
  //###################################################

  TH1I* zHisto      = NULL;
  TH2I* xyHisto     = NULL;
  TH1I* refMult     = NULL;
  TH1I* refMultHigh = NULL;
  TH1I* refMultLow  = NULL;
  TH2I* refMult_refMultHigh    = NULL;
  TH2I* refMult_refMultLow     = NULL;
  TH2I* refMultHigh_refMultLow = NULL;
  TH2I* tofMatchVsRefMult      = NULL;
  TH2I* tofMultVsRefMult       = NULL;

  TH1I* dcaHisto = NULL;
  TH2I* etaPhi   = NULL;
  TH2I* eta_pT   = NULL;
  TH2I* dEdx_Plus  = NULL;
  TH2I* dEdx_Minus = NULL;
  TH2I* dEdx_Plus_withoutBGCut  = NULL;
  TH2I* dEdx_Minus_withoutBGCut = NULL;
  TH2I* pion_y_pT               = NULL;
  TH2I* pion_y_pT_tof           = NULL;
  TH2I* kaon_y_pT_tof           = NULL;
  TH2I* proton_y_pT_tof         = NULL;
  TH2I* pion_y_pT_tof_2         = NULL;
  TH2I* kaon_y_pT_tof_2         = NULL;
  TH2I* proton_y_pT_tof_2       = NULL;
  TH2I* pion_y_pT_tof_3         = NULL;
  TH2I* kaon_y_pT_tof_3         = NULL;
  TH2I* proton_y_pT_tof_3       = NULL;

  TH2I* zTPC_byPart_hist[9][2];
  for(int iii = 0; iii < 9; iii++) for(int jjj = 0; jjj < 2; jjj++) zTPC_byPart_hist[iii][jjj] = NULL;
  TH2I* invBetaBTOF_Plus    = NULL;
  TH2I* invBetaBTOF_Minus   = NULL;
  TH2I* betaGammaBtof_dEdx_Plus  = NULL;
  TH2I* betaGammaBtof_dEdx_Minus = NULL;
  TH3I* eta_betaGammaBtof_dEdx_Plus    = NULL;
  TH3I* eta_betaGammaBtof_dEdx_Minus   = NULL;
  TH2I* invBetaGammaSqrBtof_dEdx_Plus  = NULL;
  TH2I* invBetaGammaSqrBtof_dEdx_Minus = NULL;
  TH2I* modDeDx_invBetaGammaSqrBtof_Plus   = NULL;
  TH2I* modDeDx_invBetaGammaSqrBtof_Minus  = NULL;
  TH2I* invBetaGammaPion_invBetaGammaBtof_Plus = NULL;
  TH2I* invBetaGammaPion_invBetaGammaBtof_Minus = NULL;
  TH2I* massVsMom_bTOF_Plus  = NULL;
  TH2I* massVsMom_bTOF_Minus = NULL;
  TH2I* btof_localY_localZ   = NULL;
  TH2I* pionElectronDeDx     = NULL;

  TH2I* invBetaETOFVsdEdx = NULL;
  TH2I* invBetaETOF_Plus = NULL;
  TH2I* invBetaETOF_Minus = NULL;
  TH2I* deltaXY_ETOF_Plus = NULL;
  TH2I* deltaXY_ETOF_Minus = NULL;
  TH2I* betaGammaEtof_dEdx_Plus = NULL;
  TH2I* betaGammaEtof_dEdx_Minus = NULL;
  TH2I* invBetaGammaSqrEtof_dEdx_Plus = NULL;
  TH2I* invBetaGammaSqrEtof_dEdx_Minus = NULL;
  TH2I* modDeDx_invBetaGammaSqrEtof_Plus = NULL;
  TH2I* modDeDx_invBetaGammaSqrEtof_Minus = NULL;
  TH2I* invBetaGammaPion_invBetaGammaEtof_Plus = NULL;
  TH2I* invBetaGammaPion_invBetaGammaEtof_Minus = NULL;

  TH3I* eta_betaGammaEtof_dEdx_Plus = NULL;
  TH3I* eta_betaGammaEtof_dEdx_Minus = NULL;
  TH2I* massVsMom_eTOF_Plus = NULL;
  TH2I* massVsMom_eTOF_Minus = NULL;
  TH2I* hit_vs_digi_mult_eTOF = NULL;
  TH1I* etof_goodEventFlagFraction = NULL;
  TH2I* etofGoodEventFlag = NULL;
  TH2I* eTOFNSigmaProtonPlusByETOFDet = NULL;

  string centEventsSuffix = "";
  if( a_partIndex >= 0 ){
    centEventsSuffix = particleInfo->GetParticleName(a_partIndex,0);
  }

  TH1I* centEvents = new TH1I(Form("centEvents%s",centEventsSuffix.c_str()),"Number of Events for Each Centrality Bin; Centrality Bin; Number of Events",
                        numCentBins, -0.5, ((double) numCentBins) - 0.5);
  std::vector< TH1I *> centVarHisto(numCentBins,(TH1I *)NULL);
  for (int centIndex = 0; centIndex < numCentBins; centIndex++){
    centVarHisto.at(centIndex) =
          new TH1I(Form("centralityVariable_Cent%d%s",centIndex,centEventsSuffix.c_str()),
            Form("CentralityVariable Distribution in Centrality Bin %d",centIndex),
            1000,0,1000);
  }

  #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
    std::vector<TH2I*> electron_rap_dEdx(nParticles,(TH2I *)NULL);
    std::vector<TH2I*> positron_rap_dEdx(nParticles,(TH2I *)NULL);
    std::vector<TH2I*> dEdxIsolated_Plus(nParticles,(TH2I *)NULL);
    std::vector<TH2I*> dEdxIsolated_Minus(nParticles,(TH2I *)NULL);
    std::vector<TH2I*> dEdxIsolated_Plus_Prominence(nParticles,(TH2I *)NULL);
    std::vector<TH2I*> dEdxIsolated_Minus_Prominence(nParticles,(TH2I *)NULL);
  #endif

  #ifdef _PICOBINNER_VZ_CALIBRATION_
    bool calibration_matrix[9][9] = { // given a particle of interest, which calibration histograms are necessary
      {true ,true ,true ,false,true ,false,false,false,false},
      {true ,true ,true ,false,true ,false,false,false,false},
      {true ,true ,true ,false,true ,false,false,false,false},
      {false,false,false,false,false,false,false,false,false},
      {false,false,true ,false,true ,true ,false,false,false},
      {false,false,true ,false,true ,true ,true ,false,false},
      {false,false,false,false,true ,true ,true ,true ,false},
      {false,false,false,false,true ,true ,true ,true ,false},
      {false,false,false,false,false,false,false,false,false}
    };
    const int numBetaGammaBins_Calib = 100;
    double betaGammaLimits[9][2] = {
      {0.5,5},
      {0.44,1.42},
      {0.1,1.18},
      {0,0},
      {0.1,1.6},
      {0.07,1.43},
      {0.07,1.43},
      {0.07,1.43},
      {0.0,0.0}
    };

    int vZCalib_ChargeToggles[9] = {-1,0,1,0,1,1,1,1,0};

    int numVzBinsForCalib = 9; // always ColliderCenter for this build

    TH1I* vZ_calib_struct = new TH1I("vZ_calib_struct","Bin Structure For V_{Z} dE/dx Calibration; V_{Z} [cm]; (N/A)",
                                     numVzBinsForCalib,a_cutClass->getLowZVertexCut(),a_cutClass->getHighZVertexCut());
    vector< vector< TH3I*> >  dEdx_calibration_byPart_byVz;
    vector< vector< TH3I*> >  dEdx_calibration_byPart_byVz_fromTOF;

    for(int iii = 0; iii < 9; iii++){
      vector< TH3I* > temp;
      vector< TH3I* > temp2;
      dEdx_calibration_byPart_byVz.push_back(temp);
      dEdx_calibration_byPart_byVz_fromTOF.push_back(temp2);
      for(int jjj = 0; jjj < numVzBinsForCalib; jjj++){
        dEdx_calibration_byPart_byVz[iii].push_back(NULL);
        dEdx_calibration_byPart_byVz_fromTOF[iii].push_back(NULL);
      }
    }

    int    numEtaBins_calib   = 15;
    double lowEtaRange_calib  = -1.35;
    double highEtaRange_calib =  1.35;

    for(int partIndex = 0; partIndex < 9; partIndex++){
      if( a_partIndex < 0 || calibration_matrix[a_partIndex][partIndex]){
        for(int vZ_index = 0; vZ_index < numVzBinsForCalib; vZ_index++){
          dEdx_calibration_byPart_byVz[partIndex][vZ_index] = new TH3I(Form("calib_%s_VzIndex%02d",particleInfo->GetParticleName(partIndex).Data(),vZ_index),
          Form("Calibration Histogram for %s V_{Z}=[%4.2f,%4.2f]; #eta; #beta#gamma (from p/m); ln(dE/dx) - prediction",
               particleInfo->GetParticleSymbol(partIndex).Data(),vZ_calib_struct->GetBinLowEdge(vZ_index+1),vZ_calib_struct->GetBinLowEdge(vZ_index+2)),
               numEtaBins_calib,lowEtaRange_calib,highEtaRange_calib,
               numBetaGammaBins_Calib,betaGammaLimits[partIndex][0],betaGammaLimits[partIndex][1],
               50,-0.35,0.35);
          dEdx_calibration_byPart_byVz_fromTOF[partIndex][vZ_index] = new TH3I(Form("calib_%s_VzIndex%02d_fromTOF",particleInfo->GetParticleName(partIndex).Data(),vZ_index),
          Form("Calibration Histogram for %s V_{Z}=[%4.2f,%4.2f] from TOF; #eta; #beta#gamma (from TOF); ln(dE/dx) - prediction",
               particleInfo->GetParticleSymbol(partIndex).Data(),vZ_calib_struct->GetBinLowEdge(vZ_index+1),vZ_calib_struct->GetBinLowEdge(vZ_index+2)),
               numEtaBins_calib,lowEtaRange_calib,highEtaRange_calib,
               40,betaGammaLimits[partIndex][1] - 0.15,14,
               50,-0.35,0.35);
        }
      }
    }
  #endif

  if(a_makeBasicHistos){
    //## EVENT LEVEL ### (ColliderCenter only)
    xyHisto = new TH2I("xyHisto","X vs Y Vertex Location;Event Vertex X [cm]; Event Vertex Y [cm]",1000,-5,5,1000,-5,5);
    zHisto                 = new TH1I("zVertex","Z Vertex Location; Z Vertex [cm]",4000,-400,400);
    refMult                = new TH1I("refMult","Multiplicity; refMult ; Number of Events",1000,0,1000);
    refMultHigh            = new TH1I("refMultHigh","Multiplicity; Number of Primary Tracks #eta=[-1,0]; Number of Events",1000,0,1000);
    refMultLow             = new TH1I("refMultLow","Multiplicity; Number of Primary Tracks #eta=[0,1]; Number of Events",1000,0,1000);;
    refMult_refMultHigh    = new TH2I("refMult_refMultHigh","Multiplicity;refMult; Number of Primary Tracks #eta=[0,1];",1000,0,1000,1000,0,1000);
    refMult_refMultLow     = new TH2I("refMult_refMultLow","Multiplicity;refMult; Number of Primary Tracks #eta=[0,1];",1000,0,1000,1000,0,1000);
    refMultHigh_refMultLow = new TH2I("refMult2High_refMultHigh","Multiplicity;Number of Primary Tracks #eta=[-1,0]; Number of Primary Tracks #eta=[0,1];",1000,0,1000,1000,0,1000);

    tofMatchVsRefMult  = new TH2I("tofMatchVsRefMult","TOF Matches Against RefMult; refMult; bTOF Matched Tracks",1000,0,1000,1000,0,1000);
    tofMultVsRefMult   = new TH2I("tofMultVsRefMult","TOF Multiplicity Against RefMult; refMult; bTOF Multiplicity",1000,0,1000,1000,0,1000);

    //## TRACK LEVEL ####
    dcaHisto = new TH1I("dca", "Event DCA; DCA; Counts", 1000, 0, 3.5);
    etaPhi   = new TH2I("etaPhi","Phi vs Eta; #eta; #phi",600, -2.0, 2.0,600,-3.1415,3.1415);
    eta_pT   = new TH2I("eta_pT",";#eta;p_{T}",350,-3.5,2.0,250,0,6.0);

    dEdx_Plus  = HistogramUtilities::make2DHistLogBinsInt("dEdxVsMom_Plus","dE/dx vs Primary Track Momentum;p_{tot} (GeV/c); dE/dx",
                 true, 1000,0.01, 20,true,1000,1.0,200);
    dEdx_Minus = HistogramUtilities::make2DHistLogBinsInt("dEdxVsMom_Minus","dE/dx vs Primary Track Momentum;p_{tot} (GeV/c); dE/dx",
                 true, 1000,0.01, 20,true,1000,1.0,200);
    dEdx_Plus_withoutBGCut  = HistogramUtilities::make2DHistLogBinsInt("dEdxVsMom_Plus_withoutBGCut","dE/dx vs Primary Track Momentum;p_{tot} (GeV/c); dE/dx",
                 true, 1000,0.01, 20,true,1000,1.0,200);
    dEdx_Minus_withoutBGCut = HistogramUtilities::make2DHistLogBinsInt("dEdxVsMom_Minus_withoutBGCut","dE/dx vs Primary Track Momentum;p_{tot} (GeV/c); dE/dx",
                 true, 1000,0.01, 20,true,1000,1.0,200);
    pion_y_pT        = new TH2I("pion_y_pT","Pion Acceptance; y; pT",250,-3.5,2.0,250,0,6);

    // ColliderCenter TOF acceptance binning (FixedTarget variant deleted -- dead for this build)
    pion_y_pT_tof          = new TH2I("pion_y_pT_tof","Pion Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    kaon_y_pT_tof          = new TH2I("kaon_y_pT_tof","Kaon Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    proton_y_pT_tof        = new TH2I("proton_y_pT_tof","Proton Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    pion_y_pT_tof_2        = new TH2I("pion_y_pT_tof_2","Pion Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    kaon_y_pT_tof_2        = new TH2I("kaon_y_pT_tof_2","Kaon Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    proton_y_pT_tof_2      = new TH2I("proton_y_pT_tof_2","Proton Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    pion_y_pT_tof_3        = new TH2I("pion_y_pT_tof_3","Pion Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    kaon_y_pT_tof_3        = new TH2I("kaon_y_pT_tof_3","Kaon Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);
    proton_y_pT_tof_3      = new TH2I("proton_y_pT_tof_3","Proton Acceptance in TOF; y; pT",500,-2.5,0,500,0,2.5);

    for(int partIndex = 0; partIndex < 9; partIndex++){
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        int charge = 1;
        if(pmIndex == 1) charge = -1;
        zTPC_byPart_hist[partIndex][pmIndex] = new TH2I(Form("zTPC_mom_%s",particleInfo->GetParticleName(partIndex,charge).Data()),
                                         Form("Calibration Z_TPC %s; momentum [GeV/c]; Z_{TPC}",particleInfo->GetParticleSymbol(partIndex,charge).Data()),200,0,2.5,200,-2,2);
      }
    }

    invBetaBTOF_Plus  = HistogramUtilities::make2DHistLogBinsInt("invBeta_bTOF_plus","bTOF #beta^{-1} vs Primary Track Momentum;p_{tot} (GeV/c); #beta^{-1}",
                       true,800,0.025,10,true, 800,0.7,5.0);
    invBetaBTOF_Minus = HistogramUtilities::make2DHistLogBinsInt("invBeta_bTOF_minus","bTOF #beta^{-1} vs Primary Track Momentum;p_{tot} (GeV/c); #beta^{-1}",
                       true,800,0.025,10,true, 800,0.7,5.0);
    betaGammaBtof_dEdx_Plus  = HistogramUtilities::make2DHistLogBinsInt("betaGammaBtof_dEdx_Plus", "Positive Primary Tracks from BTOF; #beta#gamma; dE/dx",
                                          true,800,0.05,100,true,800,0.8,200);
    betaGammaBtof_dEdx_Minus = HistogramUtilities::make2DHistLogBinsInt("betaGammaBtof_dEdx_Minus", "Negative Primary Tracks from BTOF; #beta#gamma; dE/dx",
                                          true,800,0.05,100,true,800,0.8,200);

    vector<double> invBetaGammmaSqr_binEdges = HistogramUtilities::makeEvenBinEdges(200, -0.15,0.003);
    vector<double> logExtendEdges            = HistogramUtilities::makeLogBinEdges(600,0.004, 120);
    for(unsigned int iii = 0; iii < logExtendEdges.size(); iii++) invBetaGammmaSqr_binEdges.push_back(logExtendEdges[iii]);
    vector<double> logDeDxEdges = HistogramUtilities::makeLogBinEdges(600,1, 200);
    vector<double> modDeDxEdges = HistogramUtilities::makeEvenBinEdges(600,-7, 16);
    invBetaGammaSqrBtof_dEdx_Plus = new TH2I("invBetaGammaSqrBtof_dEdx_Plus", "Positive Primary Tracks from BTOF; 1/(#beta#gamma)^{2}; dE/dx",
                                          invBetaGammmaSqr_binEdges.size()-1,&invBetaGammmaSqr_binEdges[0],logDeDxEdges.size()-1,&logDeDxEdges[0]);
    invBetaGammaSqrBtof_dEdx_Minus = (TH2I*) invBetaGammaSqrBtof_dEdx_Plus->Clone();
    invBetaGammaSqrBtof_dEdx_Minus->SetName("invBetaGammaSqrBtof_dEdx_Minus");
    invBetaGammaSqrBtof_dEdx_Minus->SetTitle("Negative Primary Tracks from BTOF; 1/(#beta#gamma)^{2}; dE/dx [keV/cm]");

    modDeDx_invBetaGammaSqrBtof_Plus = new TH2I("modDeDx_invBetaGammaSqrBtof_Plus","Positive Primary Tracks from BTOF;1/(#beta#gamma)^{2}; n#sigma_{dE/dx}",invBetaGammmaSqr_binEdges.size()-1,&invBetaGammmaSqr_binEdges[0],modDeDxEdges.size()-1,&modDeDxEdges[0]);
    modDeDx_invBetaGammaSqrBtof_Minus = new TH2I("modDeDx_invBetaGammaSqrBtof_Minus","Negative Primary Tracks from BTOF;1/(#beta#gamma)^{2}; n#sigma_{dE/dx}",invBetaGammmaSqr_binEdges.size()-1,&invBetaGammmaSqr_binEdges[0],modDeDxEdges.size()-1,&modDeDxEdges[0]);

    invBetaGammaPion_invBetaGammaBtof_Plus = new TH2I("invBetaGammaPion_invBetaGammaBtof_Plus","Pion BTOF; 1/#beta#gamma (from momentum) ;1/#beta#gamma (from BTOF)",
                                                      300,0.0,3,300,-0.4,5);
    invBetaGammaPion_invBetaGammaBtof_Minus = new TH2I("invBetaGammaPion_invBetaGammaBtof_Minus","Pion BTOF; 1/#beta#gamma (from momentum) ;1/#beta#gamma (from BTOF)",
                                                      300,0.0,3,300,-0.4,5);

    eta_betaGammaBtof_dEdx_Plus = HistogramUtilities::make3DHistLogBinsInt(
            "eta_betaGammaBtof_dEdx_Plus", "Positive Primary Tracks from BTOF; #eta;#beta#gamma;dE/dx",
            false,a_cutClass->getEtaPtBinStructure(2)->GetNbinsX() ,
            a_cutClass->getEtaPtBinStructure(2)->GetBinLowEdge(1),
            HistogramUtilities::getHighEdge(a_cutClass->getEtaPtBinStructure(2)),
            true,300,0.1,10,true,300,1.0,40);
    eta_betaGammaBtof_dEdx_Minus = HistogramUtilities::make3DHistLogBinsInt(
            "eta_betaGammaBtof_dEdx_Minus", "Negative Primary Tracks from BTOF; #eta;#beta#gamma;dE/dx",
            false,a_cutClass->getEtaPtBinStructure(2)->GetNbinsX() ,
            a_cutClass->getEtaPtBinStructure(2)->GetBinLowEdge(1),
            HistogramUtilities::getHighEdge(a_cutClass->getEtaPtBinStructure(2)),
            true,300,0.1,10,true,300,1.0,40);

    btof_localY_localZ   = new TH2I("btof_localY_localZ","bTOF Local Pad Distance; local y; local z",60, -3,3,60,-3,3);
    massVsMom_bTOF_Plus  = new TH2I("massVsMom_bTOF_Plus","Mass for Positive Particles BTOF ; p_{tot} (GeV/c); m/q", 800,0,8, 1000,0,4.2);
    massVsMom_bTOF_Minus = new TH2I("massVsMom_bTOF_Minus","Mass for Negative Particles BTOF ; p_{tot} (GeV/c); m/q", 800,0,8, 1000,0,4.2);
    pionElectronDeDx = HistogramUtilities::make2DHistLogBinsInt("pionElectronDeDx","dEdx for Pions and Electrons; P_{total} (GeV/c), dE/dx",true,1000,0.1,15,true,1000,0.8,220);

    //#####  ETOF SPECIFIC ###### (_HAS_ETOF_ IS defined for _OO_200_COL_)
    invBetaETOFVsdEdx = new TH2I("invBetaETOFVsdEdx","eTOF Inverse Beta vs dE/dx;dE/dx;#beta^{-1}",2000,0,80,2000,0.7,5);
    invBetaETOF_Plus = HistogramUtilities::make2DHistLogBinsInt("invBeta_eTOF_plus","eTOF #beta^{-1} vs Primary Track Momentum;p_{tot} (GeV/c); #beta^{-1}",
                       true,800,0.025,10,true, 800,0.7,5.0);
    invBetaETOF_Minus = HistogramUtilities::make2DHistLogBinsInt("invBeta_eTOF_minus","eTOF #beta^{-1} vs Primary Track Momentum;p_{tot} (GeV/c); #beta^{-1}",
                       true,800,0.025,10,true, 800,0.7,5.0);

    deltaXY_ETOF_Plus = new TH2I("deltaXY_ETOF_Plus"," DeltaXY for Positive ETOF Tracks; Delta X; Delta Y", 150, -15,15,150,-15,15);
    deltaXY_ETOF_Minus = new TH2I("deltaXY_ETOF_Minus"," DeltaXY for Negative ETOF Tracks; Delta X; Delta Y", 150, -15,15,150,-15,15);
    betaGammaEtof_dEdx_Plus = HistogramUtilities::make2DHistLogBinsInt("betaGammaEtof_dEdx_Plus", "Positive Primary Tracks from ETOF;#beta#gamma;dE/dx",
                                          true,800,0.05,100,true,600,0.8,200);
    betaGammaEtof_dEdx_Minus = HistogramUtilities::make2DHistLogBinsInt("betaGammaEtof_dEdx_Minus", "Negative Primary Tracks from ETOF;#beta#gamma;dE/dx",
                                          true,800,0.05,100,true,600,0.8,200);

    invBetaGammaSqrEtof_dEdx_Plus = (TH2I*) invBetaGammaSqrBtof_dEdx_Plus->Clone();
    invBetaGammaSqrEtof_dEdx_Plus->SetName("invBetaGammaSqrEtof_dEdx_Plus");
    invBetaGammaSqrEtof_dEdx_Plus->SetTitle("Positive Primary Tracks from ETOF; 1/(#beta#gamma)^{2}; dE/dx [keV/cm]");

    invBetaGammaSqrEtof_dEdx_Minus = (TH2I*) invBetaGammaSqrBtof_dEdx_Plus->Clone();
    invBetaGammaSqrEtof_dEdx_Minus->SetName("invBetaGammaSqrEtof_dEdx_Minus");
    invBetaGammaSqrEtof_dEdx_Minus->SetTitle("Negative Primary Tracks from ETOF; 1/(#beta#gamma)^{2}; dE/dx [keV/cm]");
    modDeDx_invBetaGammaSqrEtof_Plus = new TH2I("modDeDx_invBetaGammaSqrEtof_Plus","Positive Primary Tracks from ETOF;1/(#beta#gamma)^{2}; n#sigma_{dE/dx}",invBetaGammmaSqr_binEdges.size()-1,&invBetaGammmaSqr_binEdges[0],modDeDxEdges.size()-1,&modDeDxEdges[0]);
    modDeDx_invBetaGammaSqrEtof_Minus = new TH2I("modDeDx_invBetaGammaSqrEtof_Minus","Negative Primary Tracks from ETOF;1/(#beta#gamma)^{2}; n#sigma_{dE/dx}",invBetaGammmaSqr_binEdges.size()-1,&invBetaGammmaSqr_binEdges[0],modDeDxEdges.size()-1,&modDeDxEdges[0]);
    invBetaGammaPion_invBetaGammaEtof_Plus = new TH2I("invBetaGammaPion_invBetaGammaEtof_Plus","Pion ETOF; 1/#beta#gamma (from momentum) ;1/#beta#gamma (from BTOF)",
                                                    300,0.0,3,300,-0.4,5);
    invBetaGammaPion_invBetaGammaEtof_Minus = new TH2I("invBetaGammaPion_invBetaGammaEtof_Minus","Pion ETOF; 1/#beta#gamma (from momentum) ;1/#beta#gamma (from BTOF)",
                                                    300,0.0,3,300,-0.4,5);

    eta_betaGammaEtof_dEdx_Plus = HistogramUtilities::make3DHistLogBinsInt(
            "eta_betaGammaEtof_dEdx_Plus", "Positive Primary Tracks from ETOF; #eta;#beta#gamma;dE/dx",
            false,a_cutClass->getEtaPtBinStructure(3)->GetNbinsX() ,
            a_cutClass->getEtaPtBinStructure(3)->GetBinLowEdge(1),
            HistogramUtilities::getHighEdge(a_cutClass->getEtaPtBinStructure(3)),
            true,300,0.01,10,true,300,1.0,120);
    eta_betaGammaEtof_dEdx_Minus = HistogramUtilities::make3DHistLogBinsInt(
            "eta_betaGammaEtof_dEdx_Minus", "Negative Primary Tracks from ETOF; #eta;#beta#gamma;dE/dx",
            false,a_cutClass->getEtaPtBinStructure(3)->GetNbinsX() ,
            a_cutClass->getEtaPtBinStructure(3)->GetBinLowEdge(1),
            HistogramUtilities::getHighEdge(a_cutClass->getEtaPtBinStructure(3)),
            true,300,0.01,10,true,300,1.0,120);
    massVsMom_eTOF_Plus = new TH2I("massVsMom_eTOF_Plus","Mass for Positive Particles ETOF ; p_{tot} (GeV/c); m/q", 800,0,8, 1000,0,4.2);
    massVsMom_eTOF_Minus = new TH2I("massVsMom_eTOF_Minus","Mass for Negative Particles ETOF ; p_{tot} (GeV/c); m/q", 800,0,8, 1000,0,4.2);
    hit_vs_digi_mult_eTOF = new TH2I("hit_vs_digi_mult_eTOF","Multiplicity for eTOF; digiMult; hitMult",1201,-0.5,1200.5,601,-0.5,600.5);
    etof_goodEventFlagFraction = new TH1I("etof_goodEventFlagFraction","Fraction of Good ETOF Detectors (108 Max); Good/Total Fraction; Number of Events",108,0.0,1.0);
    etofGoodEventFlag = new TH2I("etofGoodEventFlag","Good eTOF Detectors;3*module+counter;sector",9,0,9,12,0,12);
    eTOFNSigmaProtonPlusByETOFDet = new TH2I("eTOFNSigmaProtonPlusByETOFDet","eTOF Calibration Quality;eTOF Detector ID;n_{#sigma p eTOF}",108,0,108,200, -5.0, 5.0);

    #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
      for(int partIndex = 0; partIndex < nParticles; partIndex++){
        const char *particleSymbol = particleInfo->GetParticleSymbol(partIndex,0).Data();
        const char *particleSymbolPlus = particleInfo->GetParticleSymbol(partIndex,1).Data();
        const char *particleSymbolMinus = particleInfo->GetParticleSymbol(partIndex,-1).Data();
        int nRapidityBins_btof  = a_cutClass->getRapMtM0BinStructure(partIndex,2)->GetNbinsX();
        double rapidityMin_btof = HistogramUtilities::getLowEdge(a_cutClass->getRapMtM0BinStructure(partIndex,2));
        double rapidityMax_btof = HistogramUtilities::getHighEdge(a_cutClass->getRapMtM0BinStructure(partIndex,2));
        const char *particleName = particleInfo->GetParticleName(partIndex,0).Data();

        if(a_partIndex < 0 || (a_partIndex >=0 && partIndex == a_partIndex)){
          electron_rap_dEdx.at(partIndex) = new TH2I(Form("electron_rap_dEdx_%s",particleName),
            Form("Electrons from BTOF in %s Space ;y_%s;dE/dx [kEV/cm]",particleName,particleSymbolMinus),nRapidityBins_btof,rapidityMin_btof,rapidityMax_btof,80,3.4,5.0);
          positron_rap_dEdx.at(partIndex) = new TH2I(Form("positrons_rap_dEdx_%s",particleName),
            Form("Positrons from BTOF in %s Space ;y_%s;dE/dx [kEV/cm]",particleName,particleSymbolPlus),nRapidityBins_btof,rapidityMin_btof,rapidityMax_btof,80,3.4,5.0);
        }

        dEdxIsolated_Plus.at(partIndex) = (TH2I*) dEdx_Plus->Clone();
        dEdxIsolated_Plus.at(partIndex)->
        SetName(Form("dEdx_%s_Isolated",particleInfo->GetParticleName(partIndex,1).Data()));
        dEdxIsolated_Plus.at(partIndex)->
        SetTitle(Form("dE/dx of %s ;p_{tot} (GeV/c); dE/dx (keV/cm)", particleSymbolPlus));

        dEdxIsolated_Minus.at(partIndex) = (TH2I*) dEdx_Minus->Clone();
        dEdxIsolated_Minus.at(partIndex)->
        SetName(Form("dEdx_%s_Isolated",particleInfo->GetParticleName(partIndex,-1).Data()));
        dEdxIsolated_Minus.at(partIndex)->
        SetTitle(Form("dE/dx of %s ;p_{tot} (GeV/c); dE/dx (keV/cm)", particleSymbolMinus));

        dEdxIsolated_Plus_Prominence.at(partIndex) = (TH2I*) dEdx_Plus->Clone();
        dEdxIsolated_Plus_Prominence.at(partIndex)->
        SetName(Form("dEdx_%s_Isolated_Prominence",particleInfo->GetParticleName(partIndex,1).Data()));
        dEdxIsolated_Plus_Prominence.at(partIndex)->
        SetTitle(Form("dE/dx of %s without other BTOF particles; p_{tot} (GeV/c); dE/dx (keV/cm)", particleSymbolPlus));

        dEdxIsolated_Minus_Prominence.at(partIndex) = (TH2I*) dEdx_Minus->Clone();
        dEdxIsolated_Minus_Prominence.at(partIndex)->
        SetName(Form("dEdx_%s_Isolated_Prominence",particleInfo->GetParticleName(partIndex,-1).Data()));
        dEdxIsolated_Minus_Prominence.at(partIndex)->
        SetTitle(Form("dE/dx of %s without other BTOF particles; p_{tot} (GeV/c); dE/dx (keV/cm)", particleSymbolMinus));
      }
    #endif

  }//end simple histo toggle

  //###################################################
  //###         SET UP YIELD  HISTOGRAMS            ###
  //###################################################

  std::vector<std::vector<TH3I *> > ZTPCHistoIsolatedPlus(nParticles,vector<TH3I *>(nParticles,(TH3I *)NULL));
  std::vector<std::vector<TH3I *> > ZTPCHistoIsolatedMinus(nParticles,vector<TH3I *>(nParticles,(TH3I *)NULL));
  //pion + electron histogram (hard to separate by TOF) just of pions and kaons
  std::vector< TH3I* > ZTPCHistoPionElectron(2,(TH3I *)NULL);
  std::vector< TH3I * > ZTPCHistoAllCentPlus(nParticles,(TH3I *)NULL);
  std::vector< TH3I * > ZTPCHistoAllCentMinus(nParticles,(TH3I *)NULL);

  //The Yield Histograms for the difference centrality bins and particle species
  std::vector<std::vector<TH3I *> > ZTPCHistoPlus(nParticles,vector<TH3I *>(numCentBins,(TH3I *)NULL));
  std::vector<std::vector<TH3I *> > ZBTOFHistoPlus(nParticles,vector<TH3I *>(numCentBins,(TH3I *)NULL));
  std::vector<std::vector<TH3I *> > ZTPCHistoMinus(nParticles,vector<TH3I *>(numCentBins,(TH3I *)NULL));
  std::vector<std::vector<TH3I *> > ZBTOFHistoMinus(nParticles,vector<TH3I *>(numCentBins,(TH3I *)NULL));

  //Number of tracks seen in the TPC and TOF - Will be used for TOF Matching Efficiency
  std::vector<std::vector<TH2I *> > TPCTracksPlus_ForBTOF(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  std::vector<std::vector<TH2I *> > BTOFTracksPlus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  std::vector<std::vector<TH2I *> > TPCTracksMinus_ForBTOF(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  std::vector<std::vector<TH2I *> > BTOFTracksMinus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));

  #ifdef _PICOBINNER_ETA_PT_EFF_
    std::vector<std::vector<TH2I *> > TPCTracksPlus_ForBTOF_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I *> > BTOFTracksPlus_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I *> > TPCTracksMinus_ForBTOF_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I *> > BTOFTracksMinus_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  #endif

  std::vector< TH3I* > DCA_by_cent_PionMinus(numCentBins, (TH3I *) NULL);
  std::vector< TH3I* > DCA_by_cent_ProtonPlus(numCentBins, (TH3I *) NULL);
  std::vector< TH3I* > DCA_by_cent_ProtonMinus(numCentBins, (TH3I *) NULL);
  std::vector< TH3I* > y_mTm0_DCA_Plus(nParticles, (TH3I *) NULL);
  std::vector< TH3I* > y_mTm0_DCA_Minus(nParticles, (TH3I *) NULL);
  std::vector< TH3I* > y_mTm0_DCA_woTOF_Plus(nParticles, (TH3I *) NULL);
  std::vector< TH3I* > y_mTm0_DCA_woTOF_Minus(nParticles, (TH3I *) NULL);

  std::vector<std::vector<TH3I *> > ZETOFHistoPlus(nParticles,vector<TH3I *>(numCentBins,(TH3I *)NULL));
  std::vector<std::vector<TH3I *> > ZETOFHistoMinus(nParticles,vector<TH3I *>(numCentBins,(TH3I *)NULL));
  std::vector<std::vector<TH2I *> > TPCTracksPlus_ForETOF(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  std::vector<std::vector<TH2I *> > ETOFTracksPlus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  std::vector<std::vector<TH2I *> > TPCTracksMinus_ForETOF(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  std::vector<std::vector<TH2I *> > ETOFTracksMinus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  #ifdef _PICOBINNER_ETA_PT_EFF_
    std::vector<std::vector<TH2I *> > TPCTracksPlus_ForETOF_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I *> > ETOFTracksPlus_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I *> > TPCTracksMinus_ForETOF_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I *> > ETOFTracksMinus_EtaPt(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  #endif

  #ifdef _PICOBINNER_SIMPLE_SPECTRA_
    std::vector<std::vector<TH2I*> > simple_spectra_TPC_Plus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I*> > simple_spectra_TPC_Minus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I*> > simple_spectra_bTOF_Plus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I*> > simple_spectra_bTOF_Minus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I*> > simple_spectra_eTOF_Plus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
    std::vector<std::vector<TH2I*> > simple_spectra_eTOF_Minus(nParticles,vector<TH2I *>(numCentBins,(TH2I *)NULL));
  #endif

  //Loop Over the Particle Species and centrality bins and define the histograms
  //NOTE: this still spans all 9 species -- see header comment, contamination modeling
  for (int partIndex=0; partIndex<nParticles; partIndex++){
    if(a_partIndex >= 0 && partIndex != a_partIndex) continue;
    cout << "Creating ZTPC and ZTOF for " << particleInfo->GetParticleName(partIndex,0).Data() << endl;
    //Skip Electrons and Muons
    if (partIndex == 3 || partIndex == 8) continue;

    int numMtM0Bins =  a_cutClass->getRapMtM0BinStructure(partIndex,0)->GetNbinsX();
    double MtM0min= HistogramUtilities::getLowEdge(a_cutClass->getRapMtM0BinStructure(partIndex,0));
    double MtM0max= HistogramUtilities::getHighEdge(a_cutClass->getRapMtM0BinStructure(partIndex,0));

    int nRapidityBins_tpc   = a_cutClass->getRapMtM0BinStructure(partIndex,1)->GetNbinsX();
    double rapidityMin_tpc  = HistogramUtilities::getLowEdge(a_cutClass->getRapMtM0BinStructure(partIndex,1));
    double rapidityMax_tpc  = HistogramUtilities::getHighEdge(a_cutClass->getRapMtM0BinStructure(partIndex,1));
    int nRapidityBins_btof  = a_cutClass->getRapMtM0BinStructure(partIndex,2)->GetNbinsX();
    double rapidityMin_btof = HistogramUtilities::getLowEdge(a_cutClass->getRapMtM0BinStructure(partIndex,2));
    double rapidityMax_btof = HistogramUtilities::getHighEdge(a_cutClass->getRapMtM0BinStructure(partIndex,2));

    double ZTPC_HW = 2.0;
    double ZBTOF_HW = 0.5;
    double ZBTOF_LowEdge = -ZBTOF_HW;
    double ZBTOF_HighEdge = ZBTOF_HW;
    if(a_doMassSquared){
      ZBTOF_LowEdge  = 0.5*pow(lowestMass(partIndex,1,a_cutClass,particleInfo),2.0);
      ZBTOF_HighEdge = 2.0*pow(highestMass(partIndex,1,a_cutClass,particleInfo),2.0);
      if(ZBTOF_LowEdge < 0.0) ZBTOF_LowEdge = 0.0;
    }

    double ZETOF_HW = 0.3;
    int nRapidityBins_etof  = a_cutClass->getRapMtM0BinStructure(partIndex,3)->GetNbinsX();
    double rapidityMin_etof = HistogramUtilities::getLowEdge(a_cutClass->getRapMtM0BinStructure(partIndex,3));
    double rapidityMax_etof = HistogramUtilities::getHighEdge(a_cutClass->getRapMtM0BinStructure(partIndex,3));

    //particleSymbols
    const char* particleSymbol      = particleInfo->GetParticleSymbol(partIndex,0).Data();
    const char* particleSymbolPlus  = particleInfo->GetParticleSymbol(partIndex,1).Data();
    const char* particleSymbolMinus = particleInfo->GetParticleSymbol(partIndex,-1).Data();
    string transverseAxisLabel = Form("m_{T}-m_{%s} (GeV/c^{2})",particleSymbol);

    int numZTPCBins = 600;
    int numZBTOFBins = 1000;
    int numZETOFBins = 1000;
    //get the bin edges from cut class for this particle
    vector< double > mtm0Edges = HistogramUtilities::getBinEdges(a_cutClass->getRapMtM0BinStructure(partIndex,0));
    vector< double > tpcRapEdges = HistogramUtilities::getBinEdges(a_cutClass->getRapMtM0BinStructure(partIndex,1));
    vector< double > btofRapEdges = HistogramUtilities::getBinEdges(a_cutClass->getRapMtM0BinStructure(partIndex,2));
    vector< double > etofRapEdges = HistogramUtilities::getBinEdges(a_cutClass->getRapMtM0BinStructure(partIndex,3));
    vector< double > btofMtM0Edges;
    vector< double > etofMtM0Edges;
    if(a_cutClass->getRapMtM0BinStructure(partIndex,4)) btofMtM0Edges = HistogramUtilities::getBinEdges(a_cutClass->getRapMtM0BinStructure(partIndex,4));
    else btofMtM0Edges = mtm0Edges;
    if(a_cutClass->getRapMtM0BinStructure(partIndex,5)) etofMtM0Edges = HistogramUtilities::getBinEdges(a_cutClass->getRapMtM0BinStructure(partIndex,5));
    else etofMtM0Edges = mtm0Edges;

    vector< double > zTPCEdges = HistogramUtilities::makeEvenBinEdges(numZTPCBins,-ZTPC_HW,ZTPC_HW);
    vector< double > zBTOFEdges;

    if(doMassSquaredVariableBinning && a_doMassSquared){
      double momentumSamplePoint = 0.3; // GeV/c
      double timingResolution = 35E-11 * 3.0E8; // 35 ps * c
      double lowestMassSquared = ZBTOF_LowEdge;
      if(lowestMassSquared < 0.003) lowestMassSquared =  -3.0*momentumSamplePoint*momentumSamplePoint/4.0;
      lowestMassSquared = -0.3;
      double highestMassSquared = ZBTOF_HighEdge;
      int numVariableBTOFBins = 1000;
      double topSampleVal = timingResolution*momentumSamplePoint*sqrt(highestMassSquared + momentumSamplePoint*momentumSamplePoint);
      double sampleWidth = topSampleVal/((double)numVariableBTOFBins);
      for(int edgeIndex = 0; edgeIndex <= numVariableBTOFBins; edgeIndex++){
        double edgeValue = pow(( sampleWidth*((double) edgeIndex) )/(timingResolution * momentumSamplePoint),2.0) - momentumSamplePoint*momentumSamplePoint;
        zBTOFEdges.push_back(edgeValue);
      }
      numZBTOFBins = numVariableBTOFBins;
    }else{
      zBTOFEdges = HistogramUtilities::makeEvenBinEdges(numZBTOFBins,ZBTOF_LowEdge,ZBTOF_HighEdge);
    }

    vector< double > zETOFEdges;
    if(doMassSquaredVariableBinning && a_doMassSquared){
      zETOFEdges = zBTOFEdges;
      numZETOFBins = numZBTOFBins;
    }else{
      zETOFEdges = HistogramUtilities::makeEvenBinEdges(numZETOFBins,-ZETOF_HW,ZETOF_HW);
    }

    //Define Isolated Histogram (contamination modeling -- spans all 9 species, see header note)
    for (int partIndex2=0; partIndex2 < nParticles; partIndex2++){
      const char *particleSymbol2 = particleInfo->GetParticleSymbol(partIndex2,1).Data();
      ZTPCHistoIsolatedPlus.at(partIndex).at(partIndex2) =
        new TH3I("","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size()-1,&mtm0Edges[0],numZTPCBins,&zTPCEdges[0]);
      ZTPCHistoIsolatedPlus.at(partIndex).at(partIndex2)->
        SetName(Form("zTPC_%s_Isolated",particleInfo->GetParticleName(partIndex2,1).Data()));
      ZTPCHistoIsolatedPlus.at(partIndex).at(partIndex2)->
        SetTitle(Form("Z_{TPC} of %s ;y_{%s}; %s; Z_{TPC,%s}",
        particleSymbol2,particleSymbol,transverseAxisLabel.c_str(),particleSymbol));
    }
    for (int partIndex2=0; partIndex2 < nParticles; partIndex2++){
      const char *particleSymbol2 = particleInfo->GetParticleSymbol(partIndex2,-1).Data();
      ZTPCHistoIsolatedMinus.at(partIndex).at(partIndex2) =
        new TH3I("","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size() - 1,&mtm0Edges[0],numZTPCBins,&zTPCEdges[0]);
      ZTPCHistoIsolatedMinus.at(partIndex).at(partIndex2)->
        SetName(Form("zTPC_%s_Isolated",particleInfo->GetParticleName(partIndex2,-1).Data()));
      ZTPCHistoIsolatedMinus.at(partIndex).at(partIndex2)->
        SetTitle(Form("Z_{TPC} of %s ;y_{%s}; %s; Z_{TPC,%s}",
        particleSymbol2,particleSymbol,transverseAxisLabel.c_str(),particleSymbol));
    }

    if(partIndex == 0 || partIndex == 1){
      ZTPCHistoPionElectron.at(partIndex) = new TH3I("","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size()-1,&mtm0Edges[0],numZTPCBins,&zTPCEdges[0]);
      ZTPCHistoPionElectron.at(partIndex)->SetName("zTPC_PionElectron_Isolated");
      ZTPCHistoPionElectron.at(partIndex)->SetTitle(Form("Z_{TPC} of Pion + Electrons ;y_{%s}; %s; Z_{TPC,%s}",
        particleSymbol,transverseAxisLabel.c_str(),particleSymbol));
    }

    ZTPCHistoAllCentPlus.at(partIndex) =
      new TH3I("","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size()-1,&mtm0Edges[0],numZTPCBins,&zTPCEdges[0]);
    ZTPCHistoAllCentPlus.at(partIndex)->
      SetName(Form("zTPC_%s_AllCent",particleInfo->GetParticleName(partIndex,1).Data()));
    ZTPCHistoAllCentMinus.at(partIndex) =
      new TH3I("","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size()-1,&mtm0Edges[0],numZTPCBins,&zTPCEdges[0]);
    ZTPCHistoAllCentMinus.at(partIndex)->
      SetName(Form("zTPC_%s_AllCent",particleInfo->GetParticleName(partIndex,-1).Data()));
    ZTPCHistoAllCentPlus.at(partIndex)->
      SetTitle(Form("Z_{TPC} of %s All Centralities;y_{%s};%s;Z_{TPC,%s}",
      particleSymbol,particleSymbol,transverseAxisLabel.c_str(),particleSymbol));
    ZTPCHistoAllCentMinus.at(partIndex)->
      SetTitle(Form("Z_{TPC} of %s All Centralities;y_{%s};%s;Z_{TPC,%s}",
      particleSymbol,particleSymbol,transverseAxisLabel.c_str(),particleSymbol));

    //FILL DATA BY CENTRALITY
    for(int centIndex = 0; centIndex < numCentBins; centIndex++){
      //Define Yield Histograms
      ZTPCHistoPlus.at(partIndex).at(centIndex) =
        new TH3I("","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size()-1,&mtm0Edges[0], numZTPCBins,&zTPCEdges[0]);
      ZBTOFHistoPlus.at(partIndex).at(centIndex) =
        new TH3I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0], numZBTOFBins,&zBTOFEdges[0]);

      ZTPCHistoMinus.at(partIndex).at(centIndex) =
        new TH3I("","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size()-1,&mtm0Edges[0], numZTPCBins,&zTPCEdges[0]);
      ZBTOFHistoMinus.at(partIndex).at(centIndex) =
        new TH3I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0], numZBTOFBins,&zBTOFEdges[0]);

      TPCTracksPlus_ForBTOF.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);
      BTOFTracksPlus.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);

      TPCTracksMinus_ForBTOF.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);
      BTOFTracksMinus.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);

      #ifdef _PICOBINNER_ETA_PT_EFF_
        TPCTracksPlus_ForBTOF_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);
        BTOFTracksPlus_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);

        TPCTracksMinus_ForBTOF_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);
        BTOFTracksMinus_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);
      #endif

      #ifdef _PICOBINNER_SIMPLE_SPECTRA_
        simple_spectra_TPC_Plus.at(partIndex).at(centIndex)   = new TH2I("tp","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size() - 1,&mtm0Edges[0]);
        simple_spectra_TPC_Minus.at(partIndex).at(centIndex)  = new TH2I("tm","",nRapidityBins_tpc,&tpcRapEdges[0],mtm0Edges.size() - 1,&mtm0Edges[0]);
        simple_spectra_bTOF_Plus.at(partIndex).at(centIndex)  = new TH2I("bp","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);
        simple_spectra_bTOF_Minus.at(partIndex).at(centIndex) = new TH2I("bm","",nRapidityBins_btof,&btofRapEdges[0],btofMtM0Edges.size() - 1,&btofMtM0Edges[0]);
        simple_spectra_TPC_Plus.at(partIndex).at(centIndex)->
          SetName(Form("SpectraTPC_%s_Cent%d",particleInfo->GetParticleName(partIndex,1).Data(),centIndex));
        simple_spectra_TPC_Plus.at(partIndex).at(centIndex)->
          SetTitle(Form("TPC Spectra for %s Centrality Bin: %d;y_{%s};%s",
            particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str()));
        simple_spectra_TPC_Minus.at(partIndex).at(centIndex)->
          SetName(Form("SpectraTPC_%s_Cent%d",particleInfo->GetParticleName(partIndex,-1).Data(),centIndex));
        simple_spectra_TPC_Minus.at(partIndex).at(centIndex)->
           SetTitle(Form("TPC Spectra for %s Centrality Bin: %d;y_{%s};%s",
            particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));
        simple_spectra_bTOF_Plus.at(partIndex).at(centIndex)->
          SetName(Form("SpectraBTOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,1).Data(),centIndex));
        simple_spectra_bTOF_Plus.at(partIndex).at(centIndex)->
          SetTitle(Form("bTOF Spectra for %s Centrality Bin: %d;y_{%s};%s",
            particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str()));
        simple_spectra_bTOF_Minus.at(partIndex).at(centIndex)->
          SetName(Form("SpectraBTOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,-1).Data(),centIndex));
        simple_spectra_bTOF_Minus.at(partIndex).at(centIndex)->
          SetTitle(Form("bTOF Spectra for %s Centrality Bin: %d;y_{%s};%s",
            particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));
      #endif

      ZTPCHistoPlus.at(partIndex).at(centIndex)->
        SetName(Form("zTPCPlus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      ZTPCHistoPlus.at(partIndex).at(centIndex)->
        SetTitle(Form("Z_{TPC} of %s Centrality Bin: %d;y_{%s};%s;Z_{TPC,%s}",
        particleSymbolPlus,centIndex,particleSymbolPlus,
        transverseAxisLabel.c_str(),particleSymbolPlus));

      ZTPCHistoMinus.at(partIndex).at(centIndex)->
        SetName(Form("zTPCMinus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      ZTPCHistoMinus.at(partIndex).at(centIndex)->
        SetTitle(Form("Z_{TPC} of %s Centrality Bin: %d;y_{%s};%s;Z_{TPC,%s}",
        particleSymbolMinus,centIndex,particleSymbolMinus,
        transverseAxisLabel.c_str(),particleSymbolMinus));
      ZBTOFHistoPlus.at(partIndex).at(centIndex)->
        SetName(Form("zBTOFPlus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      ZBTOFHistoMinus.at(partIndex).at(centIndex)->
        SetName(Form("zBTOFMinus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      if(!a_doMassSquared){
        ZBTOFHistoPlus.at(partIndex).at(centIndex)->
          SetTitle(Form("Z_{BTOF} of %s Centrality Bin: %d;y_{%s};%s;Z_{BTOF,%s}",
          particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str(),particleSymbolPlus));
        ZBTOFHistoMinus.at(partIndex).at(centIndex)->
          SetTitle(Form("Z_{BTOF} of %s Centrality Bin: %d;y_{%s};%s;Z_{BTOF,%s}",
          particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str(),particleSymbolMinus));
      }else{
        ZBTOFHistoPlus.at(partIndex).at(centIndex)->
          SetTitle(Form("m^{2} of %s Centrality Bin: %d;y_{%s};%s;m^{2}",
          particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str()));
        ZBTOFHistoMinus.at(partIndex).at(centIndex)->
          SetTitle(Form("m^{2} of %s Centrality Bin: %d;y_{%s};%s;m^{2}",
          particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));
      }

      //  2D Histograms
      TPCTracksPlus_ForBTOF.at(partIndex).at(centIndex)->
        SetName(Form("nTPCTracksPlus_ForBTOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      TPCTracksPlus_ForBTOF.at(partIndex).at(centIndex)->
        SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;y_{%s};%s",
        particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str()));

      BTOFTracksPlus.at(partIndex).at(centIndex)->
        SetName(Form("nBTOFTracksPlus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      BTOFTracksPlus.at(partIndex).at(centIndex)->
        SetTitle(Form("n TOF Tracks for %s Centrality Bin: %d;y_{%s};%s",
        particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str()));

      TPCTracksMinus_ForBTOF.at(partIndex).at(centIndex)->
        SetName(Form("nTPCTracksMinus_ForBTOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      TPCTracksMinus_ForBTOF.at(partIndex).at(centIndex)->
        SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;y_{%s};%s",
        particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));

      BTOFTracksMinus.at(partIndex).at(centIndex)->
        SetName(Form("nBTOFTracksMinus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      BTOFTracksMinus.at(partIndex).at(centIndex)->
        SetTitle(Form("n TOF Tracks for %s Centrality Bin: %d;y_{%s};%s",
        particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));

      #ifdef _PICOBINNER_ETA_PT_EFF_
        TPCTracksPlus_ForBTOF_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nTPCTracksPlus_ForBTOF_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        TPCTracksPlus_ForBTOF_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;#eta;p_{T}",
          particleSymbolPlus,centIndex));

        BTOFTracksPlus_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nBTOFTracksPlus_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        BTOFTracksPlus_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n TOF Tracks for %s Centrality Bin: %d;#eta;p_{T}",
          particleSymbolPlus,centIndex));

        TPCTracksMinus_ForBTOF_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nTPCTracksMinus_ForBTOF_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        TPCTracksMinus_ForBTOF_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;#eta;p_{T}",
          particleSymbolMinus,centIndex));

        BTOFTracksMinus_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nBTOFTracksMinus_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        BTOFTracksMinus_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n TOF Tracks for %s Centrality Bin: %d;#eta;p_{T}",
          particleSymbolMinus,centIndex));
      #endif

      ZETOFHistoPlus.at(partIndex).at(centIndex) =
        new TH3I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0], numZETOFBins,&zETOFEdges[0]);
      ZETOFHistoMinus.at(partIndex).at(centIndex) =
        new TH3I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0], numZETOFBins,&zETOFEdges[0]);
      TPCTracksPlus_ForETOF.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
      ETOFTracksPlus.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
      TPCTracksMinus_ForETOF.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
      ETOFTracksMinus.at(partIndex).at(centIndex) =
        new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
      #ifdef _PICOBINNER_ETA_PT_EFF_
        TPCTracksPlus_ForETOF_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
        ETOFTracksPlus_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
        TPCTracksMinus_ForETOF_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
        ETOFTracksMinus_EtaPt.at(partIndex).at(centIndex) =
          new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
      #endif

      ZETOFHistoPlus.at(partIndex).at(centIndex)->
        SetName(Form("zETOFPlus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      ZETOFHistoPlus.at(partIndex).at(centIndex)->
        SetTitle(Form("Z_{ETOF} of %s Centrality Bin: %d;y_{%s};%s;Z_{ETOF,%s}",
          particleSymbolPlus,centIndex,particleSymbolPlus,
          transverseAxisLabel.c_str(),particleSymbolPlus));
      ZETOFHistoMinus.at(partIndex).at(centIndex)->
        SetName(Form("zETOFMinus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      ZETOFHistoMinus.at(partIndex).at(centIndex)->
        SetTitle(Form("Z_{ETOF} of %s Centrality Bin: %d;y_{%s};%s;Z_{ETOF,%s}",
          particleSymbolMinus,centIndex,particleSymbolMinus,
          transverseAxisLabel.c_str(),particleSymbolMinus));
      TPCTracksPlus_ForETOF.at(partIndex).at(centIndex)->
        SetName(Form("nTPCTracksPlus_ForETOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      TPCTracksPlus_ForETOF.at(partIndex).at(centIndex)->
        SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;y_{%s};%s",
          particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str()));
      ETOFTracksPlus.at(partIndex).at(centIndex)->
        SetName(Form("nETOFTracksPlus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      ETOFTracksPlus.at(partIndex).at(centIndex)->
        SetTitle(Form("n ETOF Tracks for %s Centrality Bin: %d;y_{%s};%s",
          particleSymbolPlus,centIndex,particleSymbolPlus,transverseAxisLabel.c_str()));
      TPCTracksMinus_ForETOF.at(partIndex).at(centIndex)->
        SetName(Form("nTPCTracksMinus_ForETOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      TPCTracksMinus_ForETOF.at(partIndex).at(centIndex)->
        SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;y_{%s};%s",
          particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));
      ETOFTracksMinus.at(partIndex).at(centIndex)->
        SetName(Form("nETOFTracksMinus_%s_Cent%d",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
      ETOFTracksMinus.at(partIndex).at(centIndex)->
        SetTitle(Form("n ETOF Tracks for %s Centrality Bin: %d;y_{%s};%s",
          particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));

      #ifdef _PICOBINNER_ETA_PT_EFF_
        TPCTracksPlus_ForETOF_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nTPCTracksPlus_ForETOF_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        TPCTracksPlus_ForETOF_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;#eta;p_{T}",
            particleSymbolPlus,centIndex));
        ETOFTracksPlus_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nETOFTracksPlus_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        ETOFTracksPlus_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n ETOF Tracks for %s Centrality Bin: %d;#eta;p_{T}",
            particleSymbolPlus,centIndex));
        TPCTracksMinus_ForETOF_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nTPCTracksMinus_ForETOF_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        TPCTracksMinus_ForETOF_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n TPC Tracks for %s Centrality Bin: %d;#eta;p_{T}",
            particleSymbolMinus,centIndex));
        ETOFTracksMinus_EtaPt.at(partIndex).at(centIndex)->
          SetName(Form("nETOFTracksMinus_%s_Cent%d_EtaPt",particleInfo->GetParticleName(partIndex,0).Data(),centIndex));
        ETOFTracksMinus_EtaPt.at(partIndex).at(centIndex)->
          SetTitle(Form("n ETOF Tracks for %s Centrality Bin: %d;#eta;p_{T}",
            particleSymbolMinus,centIndex));
      #endif

      #ifdef _PICOBINNER_SIMPLE_SPECTRA_
        simple_spectra_eTOF_Plus.at(partIndex).at(centIndex)  = new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
        simple_spectra_eTOF_Minus.at(partIndex).at(centIndex) = new TH2I("","",nRapidityBins_etof,&etofRapEdges[0],etofMtM0Edges.size() - 1,&etofMtM0Edges[0]);
        simple_spectra_eTOF_Plus.at(partIndex).at(centIndex)->
          SetName(Form("SpectraETOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,1).Data(),centIndex));
        simple_spectra_eTOF_Plus.at(partIndex).at(centIndex)->
          SetTitle(Form("ETOF Spectra for %s Centrality Bin: %d;y_{%s};%s",
            particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));
        simple_spectra_eTOF_Minus.at(partIndex).at(centIndex)->
          SetName(Form("SpectraETOF_%s_Cent%d",particleInfo->GetParticleName(partIndex,-1).Data(),centIndex));
        simple_spectra_eTOF_Minus.at(partIndex).at(centIndex)->
          SetTitle(Form("ETOF Spectra for %s Centrality Bin: %d;y_{%s};%s",
            particleSymbolMinus,centIndex,particleSymbolMinus,transverseAxisLabel.c_str()));
      #endif
    }//End Loop Over Centrality Bins
  }//End Loop Over particle species

  int nRapidityBins_tpc_Proton  = a_cutClass->getRapMtM0BinStructure(2,1)->GetNbinsX();
  double rapidityMin_tpc_Proton = HistogramUtilities::getLowEdge(a_cutClass->getRapMtM0BinStructure(2,1));
  double rapidityMax_tpc_Proton = HistogramUtilities::getHighEdge(a_cutClass->getRapMtM0BinStructure(2,1));

  if(doProtonDCA){
   for(int centIndex = 0; centIndex < numCentBins; centIndex++){
      DCA_by_cent_ProtonPlus[centIndex] = new TH3I(
        Form("rap_mTm0_dca_ProtonPlus_Cent_%02d",centIndex),
        "DCA for Protons; y_{p} ; m_{T}-m_{p} ; DCA [cm]",
        nRapidityBins_tpc_Proton, rapidityMin_tpc_Proton, rapidityMax_tpc_Proton,
        40,0.0,0.08,
        60,0.0,3.0);

      DCA_by_cent_ProtonMinus[centIndex] = new TH3I(
        Form("rap_mTm0_dca_ProtonMinus_Cent_%02d",centIndex),
        "DCA for Anti-Protons; y_{p} ; m_{T}-m_{p} ; DCA [cm]",
        nRapidityBins_tpc_Proton, rapidityMin_tpc_Proton, rapidityMax_tpc_Proton,
        40,0.0,0.08,
        60,0.0,3.0);

      DCA_by_cent_PionMinus[centIndex] = new TH3I(
        Form("rap_mTm0_dca_PionMinus_Cent_%02d",centIndex),
        "DCA for Anti-Protons; y_{p} ; m_{T}-m_{p} ; DCA [cm]",
        nRapidityBins_tpc_Proton, rapidityMin_tpc_Proton, rapidityMax_tpc_Proton,
        40,0.0,0.08,
        60,0.0,3.0);
    }
  }

  for(int partIndex = 0; partIndex < nParticles; partIndex++){
    if(a_partIndex >= 0 && partIndex != a_partIndex) continue;
    if(partIndex == 3 || partIndex == 8){
      continue;
    }

    int nRapidityBins_tpc_Part  = a_cutClass->getRapMtM0BinStructure(partIndex,1)->GetNbinsX();
    double rapidityMin_tpc_Part = HistogramUtilities::getLowEdge(a_cutClass->getRapMtM0BinStructure(partIndex,1));
    double rapidityMax_tpc_Part = HistogramUtilities::getHighEdge(a_cutClass->getRapMtM0BinStructure(partIndex,1));

    y_mTm0_DCA_Plus[partIndex] = new TH3I(
        Form("rap_mTm0_dca_%s",particleInfo->GetParticleName(partIndex,1).Data()),
        Form("DCA for %s; y ; m_{T}-m_{0} ; DCA [cm]",particleInfo->GetParticleSymbol(partIndex,1).Data()),
        nRapidityBins_tpc_Part, rapidityMin_tpc_Part, rapidityMax_tpc_Part,
        75,0.0,1.5,
        60,0.0,3.0);

    y_mTm0_DCA_Minus[partIndex] = new TH3I(
        Form("rap_mTm0_dca_%s",particleInfo->GetParticleName(partIndex,-1).Data()),
        Form("DCA for %s; y ; m_{T}-m_{0} ; DCA [cm]",particleInfo->GetParticleSymbol(partIndex,-1).Data()),
        nRapidityBins_tpc_Part, rapidityMin_tpc_Part, rapidityMax_tpc_Part,
        75,0.0,1.5,
        60,0.0,3.0);

    y_mTm0_DCA_woTOF_Plus[partIndex] = new TH3I(
        Form("rap_mTm0_dca_woTOF_%s",particleInfo->GetParticleName(partIndex,1).Data()),
        Form("DCA for %s; y ; m_{T}-m_{0} ; DCA [cm]",particleInfo->GetParticleSymbol(partIndex,1).Data()),
        nRapidityBins_tpc_Part, rapidityMin_tpc_Part, rapidityMax_tpc_Part,
        75,0.0,1.5,
        60,0.0,3.0);

    y_mTm0_DCA_woTOF_Minus[partIndex] = new TH3I(
        Form("rap_mTm0_dca_woTOF_%s",particleInfo->GetParticleName(partIndex,-1).Data()),
        Form("DCA for %s; y ; m_{T}-m_{0} ; DCA [cm]",particleInfo->GetParticleSymbol(partIndex,-1).Data()),
        nRapidityBins_tpc_Part, rapidityMin_tpc_Part, rapidityMax_tpc_Part,
        75,0.0,1.5,
        60,0.0,3.0);
  }

  //###################################################
  //#########       BIN THE DATA                #######
  //###################################################

  // ----------------------------------------------------------------------------------
  // FILE-BY-FILE READING (changed from the original single-TChain-over-everything
  // approach):
  //
  // The original single-reader/single-TChain design (one StPicoDstReader built directly
  // from a_filelist, chaining all files together) segfaults reproducibly on the full
  // real fastoffline filelist, at the exact event index corresponding to a file-boundary
  // crossing, deep inside TBranchElement::GetEntry (called from
  // StPicoDstReader::readPicoEvent). This was initially misdiagnosed as being caused by
  // the SetStatus("*",0)+selective-re-enable pattern the original repo used (removed
  // above at one point) -- but after confirming that fix was genuinely compiled in and
  // loaded, the identical crash still occurred at the identical event, proving the
  // SetStatus pattern was not the actual cause. The real bug lives in how ROOT's
  // TChain/TBranchElement machinery (or the STAR-provided, unmodified
  // StPicoDstReader/StPicoDst dictionary code) handles a file-transition boundary in
  // this environment.
  //
  // Workaround: never let a single StPicoDstReader/TChain span more than one file. We
  // parse the filelist ourselves (replicating StPicoDstReader::Init's own per-line
  // parsing/corruption-check logic) and construct a fresh StPicoDstReader for each
  // individual file, fully finishing and deleting it before moving to the next. This
  // guarantees a file-boundary transition inside a single TChain instance never happens,
  // which sidesteps the bug entirely. This changes nothing about which events are read
  // or how they are analyzed -- it only changes how the input files are opened.
  // ----------------------------------------------------------------------------------
  cout << "Parsing filelist: " << a_filelist << endl;
  vector<string> picoFileList;
  {
    std::ifstream inputStream(a_filelist.c_str());
    if(!inputStream){
      cout << "PicoBinner :: Error - cannot open filelist " << a_filelist << endl;
      exit(EXIT_FAILURE);
    }
    std::string line;
    while(std::getline(inputStream, line)){
      // our external formatters may pass "file NumEvents" -- take only the first part
      size_t pos = line.find_first_of(" ");
      if(pos != std::string::npos) line.erase(pos, line.length()-pos);
      if(line.find(".picoDst.root") == std::string::npos) continue;
      TFile* ftmp = TFile::Open(line.c_str());
      if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) picoFileList.push_back(line);
      if(ftmp) ftmp->Close();
    }
  }
  cout << "Found " << picoFileList.size() << " valid picoDst files in filelist." << endl;

  //#######  TRACK INFO    ################
  bool preciseBTOF = false;
  bool   bTOF_identified_part[9];
  double nSigmaTPC_part[9];
  double zTPC_part[9];
  double nSigmaBTOF_part[9];
  double zBTOF_part[9];
  double modNSigBtof = -1;
  double invBeta = -1;
  double momentum = -1;
  double eta = -1;
  double dEdx = -1;
  double rapidity_part[9];
  double mTm0_part[9];
  double pT = -1;
  double pZ = -1;
  double mass = -1;
  double massSquared_signed = -999;
  double phi = -5;
  double betaGammaBtof = -1;
  double invBetaGammaBtof_signed = -999;
  double invBetaGammaSqrBtof_signed = -999;
  bool preciseETOF = false;
  double modNSigEtof = -1;
  double invBetaETOF = -1;
  double massETOF = -1;
  double massSquared_signed_ETOF = -1;
  double betaGammaEtof = -1;
  bool eTOF_identified_part[9];
  double nSigmaETOF_part[9];
  double zETOF_part[9];
  int eTOF_detID = -1;
  double invBetaGammaEtof_signed    = -999;
  double invBetaGammaSqrEtof_signed = -999;

  int modNum = 1;
  int currentRunNumber = -1;
  Long64_t globalEventIndex = 0; // running count across all files, purely for the "Working on event N" print/mod cadence
  //############################################################################################
  //##################      EVENT LOOP  (outer: one file at a time; inner: events in that file)#
  //############################################################################################
  for(size_t fileIndex = 0; fileIndex < picoFileList.size(); fileIndex++){

    StPicoDstReader* picoReader = new StPicoDstReader(picoFileList[fileIndex].c_str());
    picoReader->Init();
    StPicoDst* dst = picoReader->picoDst();
    Long64_t eventsInThisFile = (Long64_t) picoReader->chain()->GetEntries();

    for(Long64_t localEventIndex = 0; localEventIndex < eventsInThisFile; localEventIndex++, globalEventIndex++){
    int eventIndex = (int) globalEventIndex; // preserves original print/mod semantics across the whole run
    if( eventIndex % modNum == 0 ){
      cout << "Working on event " << eventIndex << endl;
      if (eventIndex == 10*modNum) modNum *= 10;
    }

    Bool_t readEvent = picoReader->readPicoEvent(localEventIndex);

    if( !readEvent ) {
      std::cerr << "Something went wrong. Nothing to analyze... read event was bad" << std::endl;
      std::cerr << "PicoReader Ptr: " << picoReader << endl;
      break;
    }

    StPicoEvent* event = dst->event();
    if(!event){
      cerr << "Event " << eventIndex << " doesn't exist?!?!" << endl;
      continue;
    }

    //#################     CORRECT T0  ############################
    if(a_correctT0){
      int nCan = 0;
      double tSum = 0.0;
      double t0[5000];
      memset(t0, 0.0, sizeof(t0));

      int nTracks_T0 = dst->numberOfTracks();
      for(int trackIndex=0; trackIndex<nTracks_T0; trackIndex++){
        StPicoTrack* track = dst->track(trackIndex);
        if(!track->isPrimary()) continue;
        Int_t bTofPidIndex = track->bTofPidTraitsIndex();
        if(bTofPidIndex == -1) continue;
        StPicoBTofPidTraits *traits = dst->btofPidTraits(bTofPidIndex);
        double mom = track->pPtot();
        double ch  = track->charge();
        double nSigE  = -999.;
        double nSigPi = -999.;
        double nSigK  = -999.;
        double nSigP  = -999.;
        nSigE  = track->nSigmaElectron();
        nSigPi = track->nSigmaPion();
        nSigK  = track->nSigmaKaon();
        nSigP  = track->nSigmaProton();

        double tof = traits->btof();
        if(tof == -999) continue;
        double trackL = tofPathLength(event->primaryVertex(), traits->btofHitPos(), track->helix(event->bField()).curvature());
        if(trackL == -999.) continue;

        double tofPi = trackL*sqrt(0.13957*0.13957+mom*mom)/(mom*(29.98));
        double tofP  = trackL*sqrt(0.9382723128*0.9382723128+mom*mom)/(mom*(29.98));

        if( mom>0.2 && mom<1.0 && fabs(nSigPi)<2.0 && nSigE<-3.0 && nSigK<-3.0 ){
          tSum    += tof - tofPi;
          t0[nCan] = tof - tofPi;
          nCan++;
        }
        else if( ch>0 && mom<2.0 && fabs(nSigP)<2.0 ){
          tSum    += tof - tofP;
          t0[nCan] = tof - tofP;
          nCan++;
        }
      }//end track loop

      int nTzero = nCan;
      while(nTzero>2){
        int    iMax  = 0;
        double dtMax = 0.0;
        for(int i = 0; i < nCan; i++){
          if(t0[i] == -9999.) continue;
          double tdiff = t0[i] - (tSum-t0[i])/(nTzero-1);
          if(fabs(tdiff)>dtMax){
            iMax = i;
            dtMax = fabs(tdiff);
          }
        }//end loop over candidates
        if(dtMax > sqrt(1.0+1.0/(nTzero-1))*0.5){
          tSum -= t0[iMax];
          t0[iMax] = -9999.;
          nTzero--;
        }else{
          break;
        }
      }

      double tStart = nTzero>0 ? tSum/nTzero : -9999.;
      event->setNTofT0(nTzero);
      for(int tr=0; tr<nTracks_T0; tr++){
        StPicoTrack* track = dst->track(tr);
        Int_t bTofPidIndex = track->bTofPidTraitsIndex();
        if(bTofPidIndex == -1) continue;
        StPicoBTofPidTraits *traits = dst->btofPidTraits(bTofPidIndex);
        double mom    = track->pPtot();
        double ch     = track->charge();
        double tof    = traits->btof();
        double beta   = traits->btofBeta();
        double trackL = tofPathLength(event->primaryVertex(), traits->btofHitPos(), track->helix(event->bField()).curvature());
        if(trackL == -999.){
          traits->setTOF(-999.);
          traits->setBeta(-999.);
        }else{
          traits->setTOF(tof - tStart);
          double betaCorr = trackL/((tof - tStart)*(29.98));
          traits->setBeta(betaCorr);
        }
      }//end track loop
    }//end if mCalcT0

    if(!a_cutClass->isGoodEvent(event->primaryVertex().X(), event->primaryVertex().Y(), event->primaryVertex().Z(), event->nBTOFMatch(), event->nTofT0() )){
     continue;
    }
    TVector3 eventVertex = event->primaryVertex();

    // NOTE: _STREFMULTCORR_ was never defined for _OO_200_COL_ in the original -- the
    // StRefMultCorr setup block that used to live here has been deleted; centrality
    // always uses raw refMult for this build (see the raw-refMult branch below).

    int numETOFGoodEventFlags = 0;
    for(int sector_index = 13; sector_index < 25; sector_index++){
      for(int mod_index = 0; mod_index < 9; mod_index++){
        if(event->eTofGoodEventFlag(sector_index,(mod_index/3)+1,(mod_index%3)+1)){
          numETOFGoodEventFlags++;
          etofGoodEventFlag->Fill(sector_index-13,mod_index);
        }
      }
    }
    double fracEtofGood = ((double) numETOFGoodEventFlags)/108.0;
    etof_goodEventFlagFraction->Fill(fracEtofGood);

    int centralityIndex = -1;
    int centVariable = 0;
    int refMultVar = 0;
    int refMultHighVar = 0;
    int refMultLowVar = 0;

    for (int trackIndex = 0; trackIndex < (int)dst->numberOfTracks(); trackIndex++){
      StPicoTrack* track = dst->track(trackIndex);
      if(!a_cutClass->isGoodTrack(track->nHitsFit(), track->nHits()/track->nHitsPoss(),track->nHitsDedx(), track->gDCA(eventVertex).Mag(), phi)) continue;
      double eta = track->pMom().Eta();
      if(eta >= -1.0 && eta <= 0.0) refMultHighVar++;
      if(eta >= 0.0 && eta <= 1.0) refMultLowVar++;
    }

    // ColliderCenter always, raw refMult always (see note above)
    refMultVar = event->refMult();
    if(a_cutClass->getVertexConfig() == "Center")     centVariable = refMultVar;
    if(a_cutClass->getVertexConfig() == "OffsetLow")  centVariable = refMultLowVar;
    if(a_cutClass->getVertexConfig() == "OffsetHigh") centVariable = refMultHighVar;
    centralityIndex = a_cutClass->centralityIndex(centVariable);

    //pileup cut
    if((useHuiPileup && a_cutClass->isPileUp(centVariable, event->btofTrayMultiplicity()))
      || (!useHuiPileup && a_cutClass->isPileUp(centVariable, event->nBTOFMatch())) ){
     continue;
    }
    if(!useHuiPileup && a_cutClass->isLinPileUp(centVariable) ){
      continue;
    }
    #ifdef _PICOBINNER_VZ_CALIBRATION_
      int vZ_calib_index = vZ_calib_struct->FindBin(event->primaryVertex().Z()) - 1;
    #endif

    if(a_makeBasicHistos){
      refMult->Fill(refMultVar);
      refMultHigh->Fill(refMultHighVar);
      refMultLow->Fill(refMultLowVar);
      refMult_refMultHigh->Fill(refMultVar,refMultHighVar);
      refMult_refMultLow->Fill(refMultVar,refMultLowVar);
      refMultHigh_refMultLow->Fill(refMultHighVar,refMultLowVar);
      hit_vs_digi_mult_eTOF->Fill(event->etofDigiMultiplicity(),event->etofHitMultiplicity());
      tofMatchVsRefMult->Fill(centVariable, event->nBTOFMatch());
      tofMultVsRefMult->Fill(centVariable, event->btofTrayMultiplicity());
    }
    if(centralityIndex == -1 || centralityIndex >= numCentBins){
      continue;
    }

    centEvents->Fill(centralityIndex);
    centVarHisto[centralityIndex]->Fill(centVariable);
    if(a_makeBasicHistos) zHisto->Fill(event->primaryVertex().Z());
    if(a_makeBasicHistos) xyHisto->Fill(event->primaryVertex().X(),event->primaryVertex().Y());

    //###########################################################################
    //############################      TRACK LOOP    ###########################
    //###########################################################################
    for(int trackIndex = 0; trackIndex < (int) dst->numberOfTracks(); trackIndex++){
      StPicoTrack* track = dst->track(trackIndex);
      momentum = (double) track->pPtot();
      pT       = track->pPt();
      pZ       = track->pMom().Z();
      dEdx     = track->dEdx();
      //FILL DCA for Proton and AntiProton
      if(doProtonDCA && (a_partIndex < 0 || a_partIndex == 2) && a_cutClass->isGoodTrack(track->nHitsFit(), track->nHits()/track->nHitsPoss(),track->nHitsDedx(), track->gDCA(eventVertex).Mag(), phi)){
        double rapidity = PhysMath::rapidity(particleInfo->GetParticleMass(2), pT, pZ);
        double mTm0     = PhysMath::mTm0(particleInfo->GetParticleMass(2), pT);
        Int_t bTofPidIndex = track->bTofPidTraitsIndex();
        bool hasBTOFTraits = (bTofPidIndex >=0);
        invBeta = -999;
        if(hasBTOFTraits){
          StPicoBTofPidTraits *traits = dst->btofPidTraits(bTofPidIndex);
          if(a_cutClass->isGoodBTof(traits->btofYLocal(),traits->btofZLocal())){
            invBeta =1.0/traits->btofBeta();
            if(!a_cutClass->isGoodTrack_GivenBTOFInvBetaGammaSqr(dEdx,PhysMath::invBetaGammaSqr_signed_fromInvBeta(invBeta))){
              invBeta = -999;
            }
          }
        }
        if(fabs(track->nSigmaProton()) < 2.0 && (!hasBTOFTraits || (hasBTOFTraits && invBeta > 0.0 && a_cutClass->isSpecificParticleByBTOF(invBeta,momentum,2)))){
          if(track->charge() > 0.0){
            DCA_by_cent_ProtonPlus[centralityIndex]->Fill(rapidity,mTm0,track->gDCA(eventVertex).Mag());
          }else{
            DCA_by_cent_ProtonMinus[centralityIndex]->Fill(rapidity,mTm0,track->gDCA(eventVertex).Mag());
          }
        }
        if(fabs(track->nSigmaPion()) < 2.0){
          DCA_by_cent_PionMinus[centralityIndex]->Fill(rapidity,mTm0,track->gDCA(eventVertex).Mag());
        }
      }//Proton DCA

      if(a_cutClass->isGoodTrack(track->nHitsFit(), track->nHits()/track->nHitsPoss(),track->nHitsDedx(), track->gDCA(eventVertex).Mag(), phi)){
        for(int partBaseIndex = 0; partBaseIndex < nParticles; partBaseIndex++){
          if(a_partIndex >= 0 && partBaseIndex != a_partIndex) continue;
          if(partBaseIndex == 3 || partBaseIndex == 8){
            continue;
          }

          bool validParticle = false;
          if(partBaseIndex == 0 && fabs(track->nSigmaPion()) < 2.0) validParticle = true;
          if(partBaseIndex == 1 && fabs(track->nSigmaKaon()) < 2.0) validParticle = true;
          if(partBaseIndex == 2 && fabs(track->nSigmaProton()) < 2.0) validParticle = true;
          if(partBaseIndex > 3  && fabs(log(dEdx/bichselCurves[partBaseIndex]->Eval(momentum))/log(1.07)) < 2.0) validParticle = true;
          Int_t bTofPidIndex = track->bTofPidTraitsIndex();
          bool hasBTOFTraits = (bTofPidIndex >=0);
          invBeta = -999;
          if(hasBTOFTraits){
            StPicoBTofPidTraits *traits = dst->btofPidTraits(bTofPidIndex);
            if(a_cutClass->isGoodBTof(traits->btofYLocal(),traits->btofZLocal())){
              invBeta =1.0/traits->btofBeta();
              if(!a_cutClass->isGoodTrack_GivenBTOFInvBetaGammaSqr(dEdx,PhysMath::invBetaGammaSqr_signed_fromInvBeta(invBeta))){
                invBeta = -999;
              }
            }
          }

          bool validParticleWithTOF = validParticle && (!hasBTOFTraits || (hasBTOFTraits && invBeta > 0.0 && a_cutClass->isSpecificParticleByBTOF(invBeta,momentum,partBaseIndex)));

          if(validParticle){
            double rapidity = PhysMath::rapidity(particleInfo->GetParticleMass(partBaseIndex), pT, pZ);
            double mTm0     = PhysMath::mTm0(particleInfo->GetParticleMass(partBaseIndex), pT);
            if(track->charge() > 0.0){
              y_mTm0_DCA_woTOF_Plus[partBaseIndex]->Fill(rapidity,mTm0,track->gDCA(eventVertex).Mag());
              if(validParticleWithTOF) y_mTm0_DCA_Plus[partBaseIndex]->Fill(rapidity,mTm0,track->gDCA(eventVertex).Mag());
            }else{
              y_mTm0_DCA_woTOF_Minus[partBaseIndex]->Fill(rapidity,mTm0,track->gDCA(eventVertex).Mag());
              if(validParticleWithTOF) y_mTm0_DCA_Minus[partBaseIndex]->Fill(rapidity,mTm0,track->gDCA(eventVertex).Mag());
            }
          }
        }
      }

      if(track->charge() > 0.0) dEdx_Plus_withoutBGCut->Fill(momentum,dEdx);
      else                      dEdx_Minus_withoutBGCut->Fill(momentum,dEdx);

      if(!a_cutClass->isGoodTrack(track->nHitsFit(), track->nHits()/track->nHitsPoss(),track->nHitsDedx(), track->gDCA(eventVertex).Mag(), phi)) continue;
      if(track->charge() == 0.0) continue;

      // RESET TRACK VARIABLES
      for(int partIndex= 0; partIndex < nParticles; partIndex++){
         bTOF_identified_part[partIndex] = false;
         eTOF_identified_part[partIndex] = false;
      }
      preciseBTOF = false;
      modNSigBtof = -999;
      invBeta = -1;
      mass = -1;
      massSquared_signed = -1;
      betaGammaBtof = -1;
      invBetaGammaBtof_signed = -999;
      invBetaGammaSqrBtof_signed = -999;
      preciseETOF = false;
      invBetaETOF = -1;
      massETOF = -1;
      modNSigEtof = -999;
      massSquared_signed_ETOF = -1;
      invBetaGammaBtof_signed    = -999;
      invBetaGammaSqrBtof_signed = -999;
      betaGammaEtof = -1;
      eTOF_detID = -1;

      eta = track->pMom().Eta();
      phi = TMath::ATan2(track->pMom().Y(), track->pMom().X());

      //########## Read BTOF Beta ###############
      bool passedBetaGammaCutsBTOF = true;
      bool passedElectronCutBTOF = false;
      Int_t bTofPidIndex = track->bTofPidTraitsIndex();
      if (bTofPidIndex!=-1){
        StPicoBTofPidTraits *traits = dst->btofPidTraits(bTofPidIndex);
        invBeta = 0.0;
        if(a_cutClass->isGoodBTof(traits->btofYLocal(),traits->btofZLocal())){
          invBeta = 1.0/traits->btofBeta();
          invBetaGammaSqrBtof_signed = PhysMath::invBetaGammaSqr_signed_fromInvBeta(invBeta);
          if(!a_cutClass->isGoodTrack_GivenBTOFInvBetaGammaSqr(dEdx,invBetaGammaSqrBtof_signed)){
            passedBetaGammaCutsBTOF = false;
            continue;
          }
          mass = PhysMath::mass_fromBeta(momentum, traits->btofBeta());
          massSquared_signed =  PhysMath::massSquared_signed_fromBeta(momentum,traits->btofBeta());
          betaGammaBtof = PhysMath::beta_to_betagamma(traits->btofBeta());
          invBetaGammaBtof_signed = PhysMath::signed_sqrt(invBetaGammaSqrBtof_signed);
          if(a_cutClass->isElectron_GivenBTOFInvBetaGammaSqr(dEdx,invBetaGammaSqrBtof_signed)){
            passedElectronCutBTOF = true;
          }
          preciseBTOF = true;
          btof_localY_localZ->Fill(traits->btofYLocal(),traits->btofZLocal());

          for(int partIndex= 0; partIndex < nParticles; partIndex++){
            if(partIndex == 3 || partIndex == 8) continue;
            if(a_cutClass->isSpecificParticleByBTOF(invBeta,momentum,partIndex)){
              bTOF_identified_part[partIndex] = true;
              #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
                if(track->charge() > 0.0) dEdxIsolated_Plus[partIndex]->Fill(momentum,dEdx);
                else dEdxIsolated_Minus[partIndex]->Fill(momentum,dEdx);
              #endif
            }
          }
          modNSigBtof = a_cutClass->getModNSigma_bTOF(dEdx,invBetaGammaSqrBtof_signed);
        }
      }

      //########## Read ETOF Beta ##############
      bool passedBetaGammaCutsETOF = true;
      bool passedElectronCutETOF = false;
      int eTofPidIndex = track->eTofPidTraitsIndex();
      if(eTofPidIndex!=-1){
        StPicoETofPidTraits* traits  = dst->etofPidTraits(eTofPidIndex);
        StPicoETofHit*       etofHit = dst->etofHit(traits->hitIndex());
        invBetaETOF = 0.0;
        if(a_cutClass->isGoodETof(traits,etofHit)){
          invBetaETOF                = 1.0/traits->beta();
          invBetaGammaSqrEtof_signed = PhysMath::invBetaGammaSqr_signed_fromInvBeta(invBetaETOF);
          if(!a_cutClass->isGoodTrack_GivenETOFInvBetaGammaSqr(dEdx,invBetaGammaSqrEtof_signed)){
            passedBetaGammaCutsETOF = false;
            continue;
          }
          massETOF                = PhysMath::mass_fromBeta(momentum, traits->beta());
          massSquared_signed_ETOF = PhysMath::massSquared_signed_fromBeta(momentum,traits->beta());
          betaGammaEtof           = PhysMath::beta_to_betagamma(traits->beta());
          invBetaGammaEtof_signed = PhysMath::signed_sqrt(invBetaGammaSqrEtof_signed);
          if(a_cutClass->isElectron_GivenETOFInvBetaGammaSqr(dEdx,invBetaGammaSqrEtof_signed)){
            passedElectronCutETOF = true;
          }
          if(track->charge() > 0.0){
            deltaXY_ETOF_Plus->Fill(traits->deltaX(),traits->deltaY());
          }else{
            deltaXY_ETOF_Minus->Fill(traits->deltaX(),traits->deltaY());
          }
          preciseETOF = true;
          for(int partIndex= 0; partIndex < nParticles; partIndex++){
            if(a_cutClass->isSpecificParticleByETOF(invBetaETOF,momentum,partIndex)){
              eTOF_identified_part[partIndex] = true;
              #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
                if(track->charge() > 0.0) dEdxIsolated_Plus[partIndex]->Fill(momentum,dEdx);
                else dEdxIsolated_Minus[partIndex]->Fill(momentum,dEdx);
              #endif
            }
          }
          modNSigEtof = a_cutClass->getModNSigma_eTOF(dEdx,invBetaGammaSqrEtof_signed);
          eTOF_detID = etofHit->geomId();
        }
      }

      nSigmaTPC_part[0] = track->nSigmaPion();
      nSigmaTPC_part[1] = track->nSigmaKaon();
      nSigmaTPC_part[2] = track->nSigmaProton();
      nSigmaTPC_part[3] = track->nSigmaElectron();
      for(int partIndex = 0; partIndex < nParticles; partIndex++){
         zBTOF_part[partIndex]        = -999;
         nSigmaBTOF_part[partIndex]   = -999;
         rapidity_part[partIndex]  = PhysMath::rapidity(particleInfo->GetParticleMass(partIndex), pT, pZ);
         mTm0_part[partIndex]      = PhysMath::mTm0(particleInfo->GetParticleMass(partIndex), pT);
         if(partIndex >3 ){
           nSigmaTPC_part[partIndex] = log(dEdx/bichselCurves[partIndex]->Eval(momentum))/log(1.07);
         }
         // _USE_CALIB_CLASS_ was not defined; only the Bichsel-curve path survives.
         zTPC_part[partIndex] = log(dEdx/bichselCurves[partIndex]->Eval(momentum));

         if(preciseBTOF){
            zBTOF_part[partIndex] = invBeta - PhysMath::invBeta(a_cutClass->getMassBTOF(partIndex),momentum);
            nSigmaBTOF_part[partIndex]   = zBTOF_part[partIndex]/0.02;
            if(a_doMassSquared) zBTOF_part[partIndex] = massSquared_signed;
         }
         if(preciseETOF){
            zETOF_part[partIndex] = invBetaETOF - PhysMath::invBeta(a_cutClass->getMassBTOF(partIndex),momentum);
            nSigmaETOF_part[partIndex]   = zETOF_part[partIndex]/0.015;
            if(a_doMassSquared) zETOF_part[partIndex] = massSquared_signed_ETOF;
         }else{
            zETOF_part[partIndex] = -999;
            nSigmaETOF_part[partIndex]   = -999;
         }
      }

      #ifdef _PICOBINNER_VZ_CALIBRATION_
        for(int partIndex = 0; partIndex < 9; partIndex++){
          if(dEdx_calibration_byPart_byVz[partIndex][vZ_calib_index]){
            if(vZCalib_ChargeToggles[partIndex] == 0 || (track->charge() > 0 && vZCalib_ChargeToggles[partIndex] > 0) || (track->charge() < 0 && vZCalib_ChargeToggles[partIndex] < 0)){
              dEdx_calibration_byPart_byVz[partIndex][vZ_calib_index]->Fill(eta,momentum/particleInfo->GetParticleMass(partIndex),log(dEdx/a_cutClass->getDeDxPrediciton_byPart(momentum,partIndex)));

              if(bTOF_identified_part[partIndex]){
                dEdx_calibration_byPart_byVz_fromTOF[partIndex][vZ_calib_index]->Fill(eta,betaGammaBtof,log(dEdx/a_cutClass->getDeDxPrediciton_byPart(momentum,partIndex)));
              }
              if(eTOF_identified_part[partIndex]){
                dEdx_calibration_byPart_byVz_fromTOF[partIndex][vZ_calib_index]->Fill(eta,betaGammaEtof,log(dEdx/a_cutClass->getDeDxPrediciton_byPart(momentum,partIndex)));
              }
            }
          }
        }
      #endif

      //## FILL TRACK HISTOGRAMS
      if(a_makeBasicHistos){
        etaPhi->Fill(eta,phi);
        eta_pT->Fill(eta,pT);
        dcaHisto->Fill(track->gDCA(eventVertex).Mag());
        pion_y_pT->Fill(rapidity_part[0],pT);

        // ColliderCenter binning only (FixedTarget fabs() variant deleted -- dead for this build)
        if(bTOF_identified_part[0] || eTOF_identified_part[0]) pion_y_pT_tof->Fill(rapidity_part[0],pT);
        if(bTOF_identified_part[1] || eTOF_identified_part[1]) kaon_y_pT_tof->Fill(rapidity_part[1],pT);
        if(bTOF_identified_part[2] || eTOF_identified_part[2]) proton_y_pT_tof->Fill(rapidity_part[2],pT);
        if(invBeta > 0 || invBetaETOF > 0) pion_y_pT_tof_2->Fill(rapidity_part[0],pT);
        if(invBeta > 0 || invBetaETOF > 0) kaon_y_pT_tof_2->Fill(rapidity_part[1],pT);
        if(invBeta > 0 || invBetaETOF > 0) proton_y_pT_tof_2->Fill(rapidity_part[2],pT);
        if(invBeta >= 0.0 || invBetaETOF >= 0.0) pion_y_pT_tof_3->Fill(rapidity_part[0],pT);
        if(invBeta >= 0.0 || invBetaETOF >= 0.0) kaon_y_pT_tof_3->Fill(rapidity_part[1],pT);
        if(invBeta >= 0.0 || invBetaETOF >= 0.0) proton_y_pT_tof_3->Fill(rapidity_part[2],pT);

        if(track->charge() > 0.0){
          for(int partIndex = 0 ; partIndex < 9; partIndex++) zTPC_byPart_hist[partIndex][0]->Fill(momentum, zTPC_part[partIndex]);
          dEdx_Plus->Fill(momentum,dEdx);
          if(preciseBTOF){
            invBetaBTOF_Plus->Fill(momentum, invBeta);
            massVsMom_bTOF_Plus->Fill(momentum, mass);
            betaGammaBtof_dEdx_Plus->Fill(betaGammaBtof,dEdx);
            invBetaGammaSqrBtof_dEdx_Plus->Fill(invBetaGammaSqrBtof_signed,dEdx);
            eta_betaGammaBtof_dEdx_Plus->Fill(eta,betaGammaBtof,dEdx);
            modDeDx_invBetaGammaSqrBtof_Plus->Fill(invBetaGammaSqrBtof_signed,modNSigBtof);
            invBetaGammaPion_invBetaGammaBtof_Plus->Fill(particleInfo->GetParticleMass(0)/momentum,invBetaGammaBtof_signed);
          }
          if(preciseETOF){
            invBetaETOF_Plus->Fill(momentum, invBetaETOF);
            massVsMom_eTOF_Plus->Fill(momentum, massETOF);
            betaGammaEtof_dEdx_Plus->Fill(betaGammaEtof,dEdx);
            invBetaGammaSqrEtof_dEdx_Plus->Fill(invBetaGammaSqrEtof_signed,dEdx);
            eta_betaGammaEtof_dEdx_Plus->Fill(eta,betaGammaEtof,dEdx);
            modDeDx_invBetaGammaSqrEtof_Plus->Fill(invBetaGammaSqrEtof_signed,modNSigEtof);
            invBetaGammaPion_invBetaGammaEtof_Plus->Fill(particleInfo->GetParticleMass(0)/momentum,invBetaGammaEtof_signed);
            eTOFNSigmaProtonPlusByETOFDet->Fill(eTOF_detID,nSigmaETOF_part[2]);
          }
        }else{
          for(int partIndex = 0; partIndex < 9; partIndex++) zTPC_byPart_hist[partIndex][1]->Fill(momentum,zTPC_part[partIndex]);
          dEdx_Minus->Fill(momentum,dEdx);
          if(preciseBTOF){
            invBetaBTOF_Minus->Fill(momentum, invBeta);
            massVsMom_bTOF_Minus->Fill(momentum, mass);
            betaGammaBtof_dEdx_Minus->Fill(betaGammaBtof,dEdx);
            invBetaGammaSqrBtof_dEdx_Minus->Fill(invBetaGammaSqrBtof_signed,dEdx);
            eta_betaGammaBtof_dEdx_Minus->Fill(eta,betaGammaBtof,dEdx);
            modDeDx_invBetaGammaSqrBtof_Minus->Fill(invBetaGammaSqrBtof_signed,modNSigBtof);
            invBetaGammaPion_invBetaGammaBtof_Minus->Fill(particleInfo->GetParticleMass(0)/momentum,invBetaGammaBtof_signed);
          }
          if(preciseETOF){
            invBetaETOF_Minus->Fill(momentum, invBetaETOF);
            massVsMom_eTOF_Minus->Fill(momentum, massETOF);
            betaGammaEtof_dEdx_Minus->Fill(betaGammaEtof,dEdx);
            invBetaGammaSqrEtof_dEdx_Minus->Fill(invBetaGammaSqrEtof_signed,dEdx);
            eta_betaGammaEtof_dEdx_Minus->Fill(eta,betaGammaEtof,dEdx);
            modDeDx_invBetaGammaSqrEtof_Minus->Fill(invBetaGammaSqrEtof_signed,modNSigEtof);
            invBetaGammaPion_invBetaGammaEtof_Minus->Fill(particleInfo->GetParticleMass(0)/momentum,invBetaGammaEtof_signed);
          }
        }

        #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
          for(int partIndex= 0; partIndex < nParticles; partIndex++){
            if( (!preciseBTOF) || (preciseBTOF && bTOF_identified_part[partIndex])){
              if(track->charge() > 0.0)  dEdxIsolated_Plus_Prominence[partIndex]->Fill(momentum,dEdx);
              else                       dEdxIsolated_Minus_Prominence[partIndex]->Fill(momentum,dEdx);
            }
          }
        #endif

      }//end basic histos if

      //####################################################################################################
      //######## FILL Particle Specific HISTOGRAMS #########
      //####################################################################################################
      for(int partBaseIndex = 0; partBaseIndex < nParticles; partBaseIndex++){
        if(a_partIndex >= 0 && partBaseIndex != a_partIndex) continue;
        if(partBaseIndex == 3 || partBaseIndex == 8){
          continue;
        }

        if(track->charge() > 0.0){
          if( (!a_ignoreNonParticleOfInterestByBTOF) || (a_ignoreNonParticleOfInterestByBTOF && ((!preciseBTOF) || (preciseBTOF && bTOF_identified_part[partBaseIndex])))){
            ZTPCHistoAllCentPlus[partBaseIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
          }
          ZTPCHistoPlus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
          if(preciseBTOF && !passedElectronCutBTOF){
            ZBTOFHistoPlus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zBTOF_part[partBaseIndex]);
            #ifdef _PICOBINNER_SIMPLE_SPECTRA_
              if(bTOF_identified_part[partBaseIndex]){
                simple_spectra_bTOF_Plus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
              }
            #endif
          }
          #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
            if(preciseBTOF && passedElectronCutBTOF && momentum < 1.0){
              positron_rap_dEdx[partBaseIndex]->Fill(rapidity_part[partBaseIndex],dEdx);
            }
          #endif
          if(preciseETOF && !passedElectronCutETOF){
            ZETOFHistoPlus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zETOF_part[partBaseIndex]);
            #ifdef _PICOBINNER_SIMPLE_SPECTRA_
              if(eTOF_identified_part[partBaseIndex]){
                simple_spectra_eTOF_Plus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
              }
            #endif
          }
        }else{
          if( (!a_ignoreNonParticleOfInterestByBTOF) || (a_ignoreNonParticleOfInterestByBTOF && ((!preciseBTOF) || (preciseBTOF && bTOF_identified_part[partBaseIndex])))){
            ZTPCHistoAllCentMinus[partBaseIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
          }
          ZTPCHistoMinus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
          if(preciseBTOF && !passedElectronCutBTOF){
            ZBTOFHistoMinus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zBTOF_part[partBaseIndex]);
            #ifdef _PICOBINNER_SIMPLE_SPECTRA_
              if(bTOF_identified_part[partBaseIndex]){
                simple_spectra_bTOF_Minus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
              }
            #endif
          }
          #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
            if(preciseBTOF && passedElectronCutBTOF && momentum < 1.0){
              electron_rap_dEdx[partBaseIndex]->Fill(rapidity_part[partBaseIndex],dEdx);
            }
          #endif
          if(preciseETOF && !passedElectronCutETOF){
            ZETOFHistoMinus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zETOF_part[partBaseIndex]);
            #ifdef _PICOBINNER_SIMPLE_SPECTRA_
              if(eTOF_identified_part[partBaseIndex]){
                simple_spectra_eTOF_Minus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
              }
            #endif
          }
        }
        #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
          if(preciseBTOF){
            for(int partIndex= 0; partIndex < nParticles; partIndex++){
              if(a_cutClass->isSpecificParticleByBTOF(invBeta,momentum,partIndex)){
                if(track->charge() > 0.0) ZTPCHistoIsolatedPlus[partBaseIndex][partIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
                else                      ZTPCHistoIsolatedMinus[partBaseIndex][partIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
              }
            }
            if(partBaseIndex < 2 && ( a_cutClass->isSpecificParticleByBTOF(invBeta,momentum,0) || a_cutClass->isSpecificParticleByBTOF(invBeta,momentum,3) )){
              ZTPCHistoPionElectron[partBaseIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
            }
          }
          if(preciseETOF){
            for(int partIndex= 0; partIndex < nParticles; partIndex++){
              if(a_cutClass->isSpecificParticleByETOF(invBetaETOF,momentum,partIndex)){
                if(track->charge() > 0.0) ZTPCHistoIsolatedPlus[partBaseIndex][partIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
                else                      ZTPCHistoIsolatedMinus[partBaseIndex][partIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
              }
            }
            if(partBaseIndex < 2 && ( a_cutClass->isSpecificParticleByETOF(invBetaETOF,momentum,0) || a_cutClass->isSpecificParticleByETOF(invBetaETOF,momentum,3))){
              ZTPCHistoPionElectron[partBaseIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex],zTPC_part[partBaseIndex]);
            }
          }
        #endif
        //FILL TOF EFF HISTOGRAMS
        if(fabs(nSigmaTPC_part[partBaseIndex]) < a_cutClass->getNSigmaCut(partBaseIndex)){
          #ifdef _PICOBINNER_SIMPLE_SPECTRA_
            if(track->charge() > 0.0) simple_spectra_TPC_Plus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            else                      simple_spectra_TPC_Minus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
          #endif

          if(track->charge()>0){
            TPCTracksPlus_ForBTOF[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            if(preciseBTOF) BTOFTracksPlus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            #ifdef _PICOBINNER_ETA_PT_EFF_
              TPCTracksPlus_ForBTOF_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
              if(preciseBTOF) BTOFTracksPlus_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
            #endif
            TPCTracksPlus_ForETOF[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            if(preciseETOF) ETOFTracksPlus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            #ifdef _PICOBINNER_ETA_PT_EFF_
              TPCTracksPlus_ForETOF_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
              if(preciseETOF) ETOFTracksPlus_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
            #endif
          }else{
            TPCTracksMinus_ForBTOF[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            if(preciseBTOF) BTOFTracksMinus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            #ifdef _PICOBINNER_ETA_PT_EFF_
              TPCTracksMinus_ForBTOF_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
              if(preciseBTOF) BTOFTracksMinus_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
            #endif
            if(preciseETOF) ETOFTracksMinus[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            TPCTracksMinus_ForETOF[partBaseIndex][centralityIndex]->Fill(rapidity_part[partBaseIndex],mTm0_part[partBaseIndex]);
            #ifdef _PICOBINNER_ETA_PT_EFF_
              TPCTracksMinus_ForETOF_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
              if(preciseETOF) ETOFTracksMinus_EtaPt[partBaseIndex][centralityIndex]->Fill(eta,pT);
            #endif
          }
        }
      }//end particle base loop
    }//End Track Loop

    } //end of per-file event loop

    picoReader->Finish();
    delete picoReader;
  } //end of file loop
  cout << "Finished looping over all " << globalEventIndex << " events" << endl;
  // NOTE: the original had a SECOND picoReader->Finish() call at the very end of the
  // function, guarded by `#ifndef _PICO_READER_FINISH_BYPASS_`. That toggle WAS
  // defined, so the second call never ran (it would have been a use-after-free on the
  // already-deleted picoReader above). That dead/dangerous code has been deleted. The
  // per-file Finish()+delete above (inside the file loop) is the only reader
  // teardown that happens now.

  //##################     WRITE HISTOGRAMS   ###################
  TFile* outFile = new TFile(outputFileName.c_str(),"RECREATE");
  outFile->cd();

  HistogramUtilities::ConditionalWrite(a_cutClass->getEventRejectionHistogram());
  HistogramUtilities::ConditionalWrite(a_cutClass->getTrackRejectionHistogram());
  HistogramUtilities::ConditionalWrite(a_cutClass->getCutInformationHistogram());
  HistogramUtilities::ConditionalWrite(centEvents);
  for(int iii = 0; iii < (int) centVarHisto.size(); iii++){
    HistogramUtilities::ConditionalWrite(centVarHisto[iii]);
  }

  if(a_makeBasicHistos){
    HistogramUtilities::ConditionalWrite(zHisto);
    HistogramUtilities::ConditionalWrite(xyHisto);
    HistogramUtilities::ConditionalWrite(refMult);
    HistogramUtilities::ConditionalWrite(refMultHigh);
    HistogramUtilities::ConditionalWrite(refMultLow);
    HistogramUtilities::ConditionalWrite(refMult_refMultHigh);
    HistogramUtilities::ConditionalWrite(refMult_refMultLow);
    HistogramUtilities::ConditionalWrite(refMultHigh_refMultLow);
    HistogramUtilities::ConditionalWrite(tofMatchVsRefMult);
    HistogramUtilities::ConditionalWrite(tofMultVsRefMult);
    HistogramUtilities::ConditionalWrite(dcaHisto);
    HistogramUtilities::ConditionalWrite(etaPhi);
    HistogramUtilities::ConditionalWrite(eta_pT);
    HistogramUtilities::ConditionalWrite(pion_y_pT);
    HistogramUtilities::ConditionalWrite(pion_y_pT_tof);
    HistogramUtilities::ConditionalWrite(kaon_y_pT_tof);
    HistogramUtilities::ConditionalWrite(proton_y_pT_tof);
    HistogramUtilities::ConditionalWrite(pion_y_pT_tof_2);
    HistogramUtilities::ConditionalWrite(kaon_y_pT_tof_2);
    HistogramUtilities::ConditionalWrite(proton_y_pT_tof_2);
    HistogramUtilities::ConditionalWrite(pion_y_pT_tof_3);
    HistogramUtilities::ConditionalWrite(kaon_y_pT_tof_3);
    HistogramUtilities::ConditionalWrite(proton_y_pT_tof_3);

    for(int iii = 0; iii < 9; iii++){
      HistogramUtilities::ConditionalWrite(zTPC_byPart_hist[iii][0]);
      HistogramUtilities::ConditionalWrite(zTPC_byPart_hist[iii][1]);
    }
    HistogramUtilities::ConditionalWrite(dEdx_Plus);
    HistogramUtilities::ConditionalWrite(dEdx_Minus);
    HistogramUtilities::ConditionalWrite(dEdx_Plus_withoutBGCut);
    HistogramUtilities::ConditionalWrite(dEdx_Minus_withoutBGCut);
    HistogramUtilities::ConditionalWrite(invBetaBTOF_Plus);
    HistogramUtilities::ConditionalWrite(invBetaBTOF_Minus);
    HistogramUtilities::ConditionalWrite(betaGammaBtof_dEdx_Plus);
    HistogramUtilities::ConditionalWrite(betaGammaBtof_dEdx_Minus);
    HistogramUtilities::ConditionalWrite(invBetaGammaSqrBtof_dEdx_Plus);
    HistogramUtilities::ConditionalWrite(invBetaGammaSqrBtof_dEdx_Minus);
    HistogramUtilities::ConditionalWrite(modDeDx_invBetaGammaSqrBtof_Plus);
    HistogramUtilities::ConditionalWrite(modDeDx_invBetaGammaSqrBtof_Minus);
    HistogramUtilities::ConditionalWrite(invBetaGammaPion_invBetaGammaBtof_Plus);
    HistogramUtilities::ConditionalWrite(invBetaGammaPion_invBetaGammaBtof_Minus);
    HistogramUtilities::ConditionalWrite(massVsMom_bTOF_Plus);
    HistogramUtilities::ConditionalWrite(massVsMom_bTOF_Minus);
    HistogramUtilities::ConditionalWrite(eta_betaGammaBtof_dEdx_Plus);
    HistogramUtilities::ConditionalWrite(eta_betaGammaBtof_dEdx_Minus);
    HistogramUtilities::ConditionalWrite(btof_localY_localZ);

    HistogramUtilities::ConditionalWrite(etof_goodEventFlagFraction);
    HistogramUtilities::ConditionalWrite(etofGoodEventFlag);
    HistogramUtilities::ConditionalWrite(eTOFNSigmaProtonPlusByETOFDet);
    HistogramUtilities::ConditionalWrite(hit_vs_digi_mult_eTOF);
    HistogramUtilities::ConditionalWrite(invBetaETOF_Plus);
    HistogramUtilities::ConditionalWrite(invBetaETOF_Minus);
    HistogramUtilities::ConditionalWrite(betaGammaEtof_dEdx_Plus);
    HistogramUtilities::ConditionalWrite(betaGammaEtof_dEdx_Minus);
    HistogramUtilities::ConditionalWrite(invBetaGammaSqrEtof_dEdx_Plus);
    HistogramUtilities::ConditionalWrite(invBetaGammaSqrEtof_dEdx_Minus);
    HistogramUtilities::ConditionalWrite(modDeDx_invBetaGammaSqrEtof_Plus);
    HistogramUtilities::ConditionalWrite(modDeDx_invBetaGammaSqrEtof_Minus);
    HistogramUtilities::ConditionalWrite(invBetaGammaPion_invBetaGammaEtof_Plus);
    HistogramUtilities::ConditionalWrite(invBetaGammaPion_invBetaGammaEtof_Minus);
    HistogramUtilities::ConditionalWrite(massVsMom_eTOF_Plus);
    HistogramUtilities::ConditionalWrite(massVsMom_eTOF_Minus);
    HistogramUtilities::ConditionalWrite(deltaXY_ETOF_Plus);
    HistogramUtilities::ConditionalWrite(deltaXY_ETOF_Minus);
    HistogramUtilities::ConditionalWrite(eta_betaGammaEtof_dEdx_Plus);
    HistogramUtilities::ConditionalWrite(eta_betaGammaEtof_dEdx_Minus);
    HistogramUtilities::ConditionalWrite(invBetaETOFVsdEdx);

    #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
      for(int partIndex = 0; partIndex < nParticles; partIndex++){
        if(partIndex == 3 || partIndex == 8) continue;
        HistogramUtilities::ConditionalWrite(electron_rap_dEdx[partIndex]);
        HistogramUtilities::ConditionalWrite(positron_rap_dEdx[partIndex]);
        HistogramUtilities::ConditionalWrite(dEdxIsolated_Plus[partIndex]);
        HistogramUtilities::ConditionalWrite(dEdxIsolated_Minus[partIndex]);
        if(a_ignoreNonParticleOfInterestByBTOF){
          HistogramUtilities::ConditionalWrite(dEdxIsolated_Plus_Prominence[partIndex]);
          HistogramUtilities::ConditionalWrite(dEdxIsolated_Minus_Prominence[partIndex]);
        }
      }
    #endif
  }

  TF1* calibFunct = a_cutClass->getInvBetaGammaSqrDeDxCalibFunct();
  TH1D* param_histo = new TH1D("calibration_InvBetaGammaSqr_dEdx_params","Params For Calibration Function",6,-0.5,5.5);
  param_histo->GetXaxis()->SetBinLabel(1,"param 0");
  param_histo->GetXaxis()->SetBinLabel(2,"param 1");
  param_histo->GetXaxis()->SetBinLabel(3,"param 2");
  param_histo->GetXaxis()->SetBinLabel(4,"param 3");
  param_histo->GetXaxis()->SetBinLabel(5,"param 4");
  param_histo->GetXaxis()->SetBinLabel(6,"Should be 1.0");
  param_histo->SetBinContent(1,calibFunct->GetParameter(0));
  param_histo->SetBinContent(2,calibFunct->GetParameter(1));
  param_histo->SetBinContent(3,calibFunct->GetParameter(2));
  param_histo->SetBinContent(4,calibFunct->GetParameter(3));
  param_histo->SetBinContent(5,calibFunct->GetParameter(4));
  param_histo->SetBinContent(6,1.0);
  HistogramUtilities::ConditionalWrite(param_histo);

  #ifdef _PICOBINNER_VZ_CALIBRATION_
    outFile->mkdir("Vz_Calibration");
    outFile->cd("Vz_Calibration");
    HistogramUtilities::ConditionalWrite(vZ_calib_struct);
    for(int partIndex = 0; partIndex < 9; partIndex++){
      for(int vZ_index = 0; vZ_index < numVzBinsForCalib; vZ_index++){
        HistogramUtilities::ConditionalWrite(dEdx_calibration_byPart_byVz[partIndex][vZ_index]);
        HistogramUtilities::ConditionalWrite(dEdx_calibration_byPart_byVz_fromTOF[partIndex][vZ_index]);
      }
    }
  #endif

  if(doProtonDCA && (a_partIndex < 0 || a_partIndex == 2 )){
    outFile->mkdir("ProtonDCA");
    outFile->cd("ProtonDCA");
    for(int centIndex = 0; centIndex < numCentBins; centIndex++){
      HistogramUtilities::ConditionalWrite(DCA_by_cent_ProtonPlus[centIndex]);
      HistogramUtilities::ConditionalWrite(DCA_by_cent_ProtonMinus[centIndex]);
      HistogramUtilities::ConditionalWrite(DCA_by_cent_PionMinus[centIndex]);
    }
  }

  //Write Y and mT-m0  in it's own directory
  outFile->mkdir("yields");
  outFile->cd("yields");
  for(int partIndex = 0; partIndex < nParticles; partIndex++){
    if(a_partIndex >= 0 && partIndex != a_partIndex) continue;
    if(partIndex ==3 || partIndex ==8) continue;
    outFile->mkdir("yields/" + particleInfo->GetParticleName(partIndex,0)+"_space");
    outFile->cd("yields/" + particleInfo->GetParticleName(partIndex,0)+"_space");
    #ifdef _PICOBINNER_BTOF_ISOLATED_HISTOS_
      for(int partIndex2 = 0; partIndex2 < nParticles; partIndex2++){
        HistogramUtilities::ConditionalWrite(ZTPCHistoIsolatedPlus[partIndex][partIndex2]);
        HistogramUtilities::ConditionalWrite(ZTPCHistoIsolatedMinus[partIndex][partIndex2]);
      }
      if(partIndex == 0 || partIndex == 1){ // only pions and kaons matter
        HistogramUtilities::ConditionalWrite(ZTPCHistoPionElectron[partIndex]);
      }
    #endif
    HistogramUtilities::ConditionalWrite(ZTPCHistoAllCentMinus[partIndex]);
    HistogramUtilities::ConditionalWrite(ZTPCHistoAllCentPlus[partIndex]);

    for(int centIndex = 0; centIndex < numCentBins; centIndex++){
      #ifdef _PICOBINNER_SIMPLE_SPECTRA_
        HistogramUtilities::ConditionalWrite(simple_spectra_TPC_Plus[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(simple_spectra_TPC_Minus[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(simple_spectra_bTOF_Plus[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(simple_spectra_bTOF_Minus[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(simple_spectra_eTOF_Plus[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(simple_spectra_eTOF_Minus[partIndex][centIndex]);
      #endif
      HistogramUtilities::ConditionalWrite(ZTPCHistoPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(ZTPCHistoMinus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(ZBTOFHistoPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(ZBTOFHistoMinus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(TPCTracksPlus_ForBTOF[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(BTOFTracksPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(TPCTracksMinus_ForBTOF[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(BTOFTracksMinus[partIndex][centIndex]);

      #ifdef _PICOBINNER_ETA_PT_EFF_
        HistogramUtilities::ConditionalWrite(TPCTracksPlus_ForBTOF_EtaPt[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(BTOFTracksPlus_EtaPt[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(TPCTracksMinus_ForBTOF_EtaPt[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(BTOFTracksMinus_EtaPt[partIndex][centIndex]);
      #endif

      HistogramUtilities::ConditionalWrite(ZETOFHistoPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(ZETOFHistoMinus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(TPCTracksPlus_ForETOF[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(ETOFTracksPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(TPCTracksMinus_ForETOF[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(ETOFTracksMinus[partIndex][centIndex]);

      #ifdef _PICOBINNER_ETA_PT_EFF_
        HistogramUtilities::ConditionalWrite(TPCTracksPlus_ForETOF_EtaPt[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(ETOFTracksPlus_EtaPt[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(TPCTracksMinus_ForETOF_EtaPt[partIndex][centIndex]);
        HistogramUtilities::ConditionalWrite(ETOFTracksMinus_EtaPt[partIndex][centIndex]);
      #endif
    }

    for(int partIndex = 0; partIndex < nParticles; partIndex++){
      HistogramUtilities::ConditionalWrite(y_mTm0_DCA_Plus[partIndex]);
      HistogramUtilities::ConditionalWrite(y_mTm0_DCA_Minus[partIndex]);
      HistogramUtilities::ConditionalWrite(y_mTm0_DCA_woTOF_Plus[partIndex]);
      HistogramUtilities::ConditionalWrite(y_mTm0_DCA_woTOF_Minus[partIndex]);
    }
  }

  outFile->Close();
  cout << "Wrote Histograms" << endl;
  cout << "Finished PicoBinner Function Call -- GoodBye" << endl;
  return;
}//end of function
