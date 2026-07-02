#ifndef __CUT_CLASS__
#define __CUT_CLASS__

#include <vector>
#include <string>
#include <iostream>

//#include "TMath.h"
//#include "TObject.h"


#include "../makefile_toggles.h"
#include "HistogramUtilities.h"
#include "PhysMath.h"

#ifndef _MAC_OSX_
  #ifdef _PICO_READER_
    #include "../submodule/PicoDstReader/StPicoEvent.h"
    #include "../submodule/PicoDstReader/StPicoTrack.h"
    #include "../submodule/PicoDstReader/StPicoBTofHit.h"
    #include "../submodule/PicoDstReader/StPicoBTofPidTraits.h"
    #include "../submodule/PicoDstReader/StPicoETofPidTraits.h"
    #include "../submodule/PicoDstReader/StPicoETofHit.h"
  #endif

  #ifdef _PICO_READER_SL22b_
    #include "../submodule/PicoDstReader_SL22b/StPicoEvent.h"
    #include "../submodule/PicoDstReader_SL22b/StPicoTrack.h"
    #include "../submodule/PicoDstReader_SL22b/StPicoBTofHit.h"
    #include "../submodule/PicoDstReader_SL22b/StPicoBTofPidTraits.h"
    #include "../submodule/PicoDstReader_SL22b/StPicoETofPidTraits.h"
    #include "../submodule/PicoDstReader_SL22b/StPicoETofHit.h"
  #endif

  #ifdef _PICO_READER_SL22c_
    #include "../submodule/PicoDstReader_SL22c/StPicoEvent.h"
    #include "../submodule/PicoDstReader_SL22c/StPicoTrack.h"
    #include "../submodule/PicoDstReader_SL22c/StPicoBTofHit.h"
    #include "../submodule/PicoDstReader_SL22c/StPicoBTofPidTraits.h"
    #include "../submodule/PicoDstReader_SL22c/StPicoETofPidTraits.h"
    #include "../submodule/PicoDstReader_SL22c/StPicoETofHit.h"
  #endif

  #ifdef _PICO_READER_SL23a_
    #include "../submodule/PicoDstReader_SL23a/StPicoEvent.h"
    #include "../submodule/PicoDstReader_SL23a/StPicoTrack.h"
    #include "../submodule/PicoDstReader_SL23a/StPicoBTofHit.h"
    #include "../submodule/PicoDstReader_SL23a/StPicoBTofPidTraits.h"
    #include "../submodule/PicoDstReader_SL23a/StPicoETofPidTraits.h"
    #include "../submodule/PicoDstReader_SL23a/StPicoETofHit.h"
  #endif

  #ifdef _PICO_READER_SL23d_
    #include "../submodule/PicoDstReader_SL23d/StPicoEvent.h"
    #include "../submodule/PicoDstReader_SL23d/StPicoTrack.h"
    #include "../submodule/PicoDstReader_SL23d/StPicoBTofHit.h"
    #include "../submodule/PicoDstReader_SL23d/StPicoBTofPidTraits.h"
    #include "../submodule/PicoDstReader_SL23d/StPicoETofPidTraits.h"
    #include "../submodule/PicoDstReader_SL23d/StPicoETofHit.h"
  #endif

  #ifdef _PICO_READER_SL19b_
    #include "../submodule/PicoDstReader_SL19b/StPicoEvent.h"
    #include "../submodule/PicoDstReader_SL19b/StPicoTrack.h"
    #include "../submodule/PicoDstReader_SL19b/StPicoBTofHit.h"
    #include "../submodule/PicoDstReader_SL19b/StPicoBTofPidTraits.h"
  #endif

  #ifdef _SL19b_STREFMULTCORR_
    #include "../submodule/StRoot/StRefMultCorr_SL19b/StRefMultCorr.h"
    #include "../submodule/StRoot/StRefMultCorr_SL19b/CentralityMaker.h"
  #endif

  #ifdef _SL21a_STREFMULTCORR_
    #include "../submodule/StRoot/StRefMultCorr_SL21a/StRefMultCorr.h"
    #include "../submodule/StRoot/StRefMultCorr_SL21a/CentralityMaker.h"
  #endif

  #ifdef _SL22dev_STREFMULTCORR_
    #include "../submodule/StRoot/StRefMultCorr_SL22dev/StRefMultCorr.h"
    #include "../submodule/StRoot/StRefMultCorr_SL22dev/CentralityMaker.h"
  #endif

  #ifdef _7p7GeV_STREFMULTCORR_
    #include "../submodule/StRoot/StRefMultCorr_7p7GeV/StRefMultCorr.h"
    #include "../submodule/StRoot/StRefMultCorr_7p7GeV/CentralityMaker.h"
  #endif

  #ifdef _9p2_11p5_17p3GeV_STREFMULTCORR_
   #include "../submodule/StRoot/StRefMultCorr_9p2_11p5_17p3GeV/StRefMultCorr.h"
   #include "../submodule/StRoot/StRefMultCorr_9p2_11p5_17p3GeV/CentralityMaker.h"
  #endif



#endif

// ---------------------------------------------------------------------------
// isGoodETof() below needs StPicoETofPidTraits/StPicoETofHit's real definitions
// regardless of _MAC_OSX_. In the original repo these types were only ever
// included inside the "#ifndef _MAC_OSX_" block above, which meant isGoodETof
// could never compile on a local/Mac build (see source/CutClass.cxx, where its
// definition -- and only its definition, nothing else it depends on -- was
// nested inside that same "#ifndef _MAC_OSX_ // no StPicoDst stuff" block).
// That combination (_HAS_ETOF_ + local Mac build) was never actually exercised
// in the original repo: its local/"Not on RCF" toggle branch has no _OO_200_COL_
// case at all, so _HAS_ETOF_ was never turned on together with a Mac build there.
// Since this repo needs real eTOF PID locally, these two includes and the
// isGoodETof() declaration are pulled out from behind that gate -- gated only on
// _HAS_ETOF_, matching how PicoBinner.cxx itself gates its ETOF code.
#ifdef _HAS_ETOF_
  #include "../submodule/PicoDstReader_SL23a/StPicoETofPidTraits.h"
  #include "../submodule/PicoDstReader_SL23a/StPicoETofHit.h"
#endif

#include "TH1D.h"


class CutClass{
public:
  CutClass();
  //~CutClass(){};

  //SETTERS
  void setTriggers(int a_numTriggers, unsigned int* a_triggers);
  void setHLTTriggers(int a_numHLTTriggers, unsigned int* a_HLTTriggers);
  #ifdef _HAS_ETOF_
    void setEToFTriggers(int a_numEToFTriggers, unsigned int* a_EToFTriggers);
  #endif

  // assumes 0->_%, so start at 5,10,... then lower bounds for those 350,300,___
  void setCentralities(int a_numCent, double* a_cent_edges, int* a_cent_percents); //not needed if using StRefMultCorr
  void setTrigToggle(bool a_toggle){m_cutOnTriggers = a_toggle;};

  void setBeamLocation(double a_x, double a_y){m_beamXCenter = a_x; m_beamYCenter = a_y;};
  //  void setBeamLocation(double a_x, double a_y);
  void setZRange(double a_zLow, double a_zHigh);
  void setDCA(double a_dca){m_DCA = a_dca;};
  void setPhiCut(bool a_useCut, double a_lowCut = 0.0, double a_highCut=0.0){
    m_usePhiCut = a_useCut;
    m_lowPhiCut = a_lowCut;
    m_highPhiCut = a_highCut;
  };
  void setEtaCut(bool a_useCut, double a_lowCut = 0.0, double a_highCut = 0.0){
    m_useEtaCut = a_useCut;
    m_lowEtaCut = a_lowCut;
    m_highEtaCut = a_highCut;
  };

  void setPtCut(bool a_useCut, double a_lowCut = 0.0){
    m_usePtCut = a_useCut;
    m_lowPtCut = a_lowCut;
  };

  void setPileUpCut(double a_cut=1e150){m_pileUpCut = a_cut;}
  void setPileUpLinearCuts(double a_highSlope, double a_highIntercept, double a_lowSlope, double a_lowIntercept){
       m_pileUpLineCuts[0] = a_highSlope;
       m_pileUpLineCuts[1] = a_highIntercept;
       m_pileUpLineCuts[2] = a_lowSlope;
       m_pileUpLineCuts[3] = a_lowIntercept;
  };
  void setPileUpLinearCutsByLowPoints(double a_x1, double a_y1, double a_x2, double a_y2,double a_x3, double a_y3,double a_x4, double a_y4){
       //m_pileUpLineCuts[0] = a_highSlope;
       //m_pileUpLineCuts[1] = a_highIntercept;
       m_pileUpLineCuts[2] = (a_y2 - a_y1)/(a_x2 - a_x1);
       m_pileUpLineCuts[3] = a_y2 - m_pileUpLineCuts[2]*a_x2;
       m_pileUpLineCuts[4] = (a_y4 - a_y3)/(a_x4 - a_x3);
       m_pileUpLineCuts[5] = a_y4 - m_pileUpLineCuts[4]*a_x4;
       cout << "Lower Pile-up Cuts:  bTofMatch < " << m_pileUpLineCuts[2] << " * refMult + " << m_pileUpLineCuts[3] << endl;
       cout << "Lower Pile-up Cuts:  bTofMatch < " << m_pileUpLineCuts[4] << " * refMult + " << m_pileUpLineCuts[5] << endl;
       if(m_pileUpLineCuts[2] < 0.0 || m_pileUpLineCuts[4] < 0.0) cout << "WARNING: Your pileup line has a negative slope, that isn't normal" << endl;


  };


  void useHuiDT27GeVPileup(bool a_toggle = true){m_useHui_dt_27GeVPileup = a_toggle;};

  void setRadialCut(double a_rad){m_radialCut = a_rad;};
  void setTofMatch(int a_tofMatches){ m_tofMatchCut = a_tofMatches; };

  void setUseT0(bool a_useT0){ m_useT0 = a_useT0; };
  void setTofT0Cut(int a_nTofT0Cut){ m_nTofT0Cut = a_nTofT0Cut; };
  void setNSigmaCuts(vector<double> a_nSigmaCuts){m_nSigmaCuts = a_nSigmaCuts;};

  void setNHitsFit(double a_hits){m_NHitsFit = a_hits;};
  void setNHitsDeDx(double a_hits){m_NHitsDeDx = a_hits;};
  void setFitMaxRatio(double a_ratio){m_fitMaxRatio = a_ratio;};
  void setBTOF(double a_localY, double a_localZ){m_btofLocY = a_localY; m_btofLocZ = a_localZ;};
  void setETOF(double a_meanDeltaX, double a_meanDeltaY, double a_widthDeltaX, double a_widthDeltaY){
               m_etofDetlaXY_Center_X = a_meanDeltaX;       m_etofDeltaXY_Center_Y = a_meanDeltaY;
               m_etofDeltaX_Ellipse_Cut = a_widthDeltaX;    m_etofDeltaY_Ellipse_Cut = a_widthDeltaY;};
  void setETOFMaxClusterSize(int a_cut){  m_etof_max_cluster_size = a_cut;};
  void setETOFMatchFlagCut(int a_cut){    m_etof_match_flag_cut = a_cut;};
  void setYCM(Float_t a_yCM);
  void setEventConfig(string a_eventConfig){m_eventConfig = a_eventConfig;};
  void setVertexConfig(string a_vertexConfig){m_vertexConfig = a_vertexConfig;}; // OffsetLow OffsetHigh Center
  void setStarver(string a_starver){m_starver = a_starver;};
  void setSqrts_NN(Double_t a_sqrts_NN){ m_sqrts_NN = a_sqrts_NN; };
  void setBeamEnergy(Double_t a_beamEnergy){ m_beamEnergy = a_beamEnergy; };

  void setEtaPtBinningInfo( int a_numPtBins, double a_pTMin, double a_pTMax,
                        int a_numEtaBins_TPC, double a_EtaMin_TPC, double a_EtaMax_TPC,
                        int a_numEtaBins_bTOF, double a_EtaMin_bTOF, double a_EtaMax_bTOF,
                        int a_numEtaBIns_eTOF, double a_EtaMin_eTOF, double a_EtaMax_eTOF);

  void setLogMTBinning(bool a_toggle){m_useLogMTBinning = a_toggle;};
  //must be by particle index in ParticleInfo
  void setRapMtM0BinningInfo(vector<int> a_numMtM0Bins, vector<double> a_mTm0Mins, vector<double> a_mTm0Maxs,
                              vector<int> a_numRapBins_tpc, vector<double> a_RapMins_tpc, vector<double> a_RapMaxs_tpc,
                              vector<int> a_numRapBins_btof, vector<double> a_RapMins_btof, vector<double> a_RapMaxs_btof,
                              vector<int> a_numRapBins_etof, vector<double> a_RapMins_etof, vector<double> a_RapMaxs_etof);

  void setVariableRapMtM0BinningInfo(int a_partIndex,
                                     int a_numRapBins_tpc, double a_lowRap_tpc, double a_highRap_tpc,
                                     int a_numMtM0Bins_tpc, vector<double> a_mTm0BinEdges_tpc,
                                     int a_numRapBins_btof, double a_lowRap_btof, double a_highRap_btof,
                                     int a_numMtM0Bins_btof, vector<double> a_mTm0BinEdges_btof,
                                     int a_numRapBins_etof, double a_lowRap_etof, double a_highRap_etof,
                                     int a_numMtM0Bins_etof, vector<double> a_mTm0BinEdges_etof);

  void setVariableRapMtM0BinningInfo(int a_partIndex,
                                     int a_numRapBins_tpc, double a_lowRap_tpc, double a_highRap_tpc,
                                     int a_numMtM0Bins_tpc, double* a_mTm0BinEdges_tpc,
                                     int a_numRapBins_btof, double a_lowRap_btof, double a_highRap_btof,
                                     int a_numMtM0Bins_btof, double* a_mTm0BinEdges_btof,
                                     int a_numRapBins_etof, double a_lowRap_etof, double a_highRap_etof,
                                     int a_numMtM0Bins_etof, double* a_mTm0BinEdges_etof);
  //void setRapidityBinningInfo(int a_numRapidityBins, double a_RapMin, double a_RapMax ){
  //              m_numRapBins = a_numRapidityBins; m_rapMin = a_RapMin; m_rapMax = a_RapMax; };

  void setMassCuts(vector<double> a_lowMassCut, vector<double> a_highMassCut){
              m_lowMassCut = a_lowMassCut; m_highMassCut = a_highMassCut;};
  void setParticleMass_bTOF(int a_partIndex,double a_mass){m_particleMasses_bTOF[a_partIndex] = a_mass;};
  void setParticleMass_eTOF(int a_partIndex,double a_mass){m_particleMasses_eTOF[a_partIndex] = a_mass;};
  void setInvBetaWidth_bTOF(int a_partIndex,double a_invBetaWidth){m_invBetaCut_bTOF[a_partIndex] = a_invBetaWidth;};
  void setInvBetaWidth_eTOF(int a_partIndex,double a_invBetaWidth){m_invBetaCut_eTOF[a_partIndex] = a_invBetaWidth;};
  void setParticleMasses_bTOF(vector<double> a_masses){ m_particleMasses_bTOF = a_masses;};
  void setParticleMasses_eTOF(vector<double> a_masses){ m_particleMasses_eTOF = a_masses;};
  void setInvBetaCuts_bTOF(vector<double> a_invBetaWidth){ m_invBetaCut_bTOF = a_invBetaWidth; };
  void setInvBetaCuts_eTOF(vector<double> a_invBetaWidth){ m_invBetaCut_eTOF = a_invBetaWidth; };
  void setModNSigmaCuts_bTOF(bool a_toggle,double a_lowStart, double a_lowCut, double a_highCut, double a_lowInvBGCut){
                   m_useModNSigma_InvBetaGammaSqrBTOFCuts = a_toggle;
                   m_lowInvBetaGammaSqr_ToUseModNSigmaBTOFCut = a_lowStart;
                   m_lowModNSigma_BTOF_Cut = a_lowCut;
                   m_highModNSigma_BTOF_Cut = a_highCut;
                   m_lowInvBetaGammaSqrBTOFCut = a_lowInvBGCut;};

  void setModNSigmaCuts_eTOF(bool a_toggle,double a_lowStart, double a_lowCut, double a_highCut, double a_lowInvBGCut){
                   m_useModNSigma_InvBetaGammaSqrETOFCuts = a_toggle;
                   m_lowInvBetaGammaSqr_ToUseModNSigmaETOFCut = a_lowStart;
                   m_lowModNSigma_ETOF_Cut = a_lowCut;
                   m_highModNSigma_ETOF_Cut = a_highCut;
                   m_lowInvBetaGammaSqrETOFCut = a_lowInvBGCut;};
  //these only cut TOF, not dE/dx
  void setLowInvBetaGammaSqrBTOFCut(double a_value){m_lowInvBetaGammaSqrBTOFCut = a_value;};
  void setLowInvBetaGammaSqrETOFCut(double a_value){m_lowInvBetaGammaSqrETOFCut = a_value;};


  void setParticleFitMatrix(bool a_particleFitMatrix[9][9], int a_charge);

  void setBTOFBetaGammaCuts(bool a_useCuts, vector<double> a_lowLineXPts, vector<double> a_lowLineYPts,
                            vector<double> a_highLineXPts, vector<double> a_highLineYPts);
  void setBTOFInvBetaGammaSqrCuts(bool a_useCuts, vector<double> a_lowLineXPts, vector<double> a_lowLineYPts,
                            vector<double> a_highLineXPts, vector<double> a_highLineYPts);
  void setBTOFElectronCut(bool a_toggle){m_useBTOFElectronCut = a_toggle;};
  void setETOFInvBetaGammaSqrCuts(bool a_useCuts, vector<double> a_lowLineXPts, vector<double> a_lowLineYPts,
                            vector<double> a_highLineXPts, vector<double> a_highLineYPts);
  void setETOFElectronCut(bool a_toggle){m_useETOFElectronCut = a_toggle;};


  void setMomStops_TPC( vector<double> a_momStops_TPC_Plus, vector<double> a_momStops_TPC_Minus ){
                        m_momStops_TPC_Plus = a_momStops_TPC_Plus;
                        m_momStops_TPC_Minus = a_momStops_TPC_Minus; };
  void setMomStops_bToF( vector<double> a_momStops_bToF_Plus, vector<double> a_momStops_bToF_Minus ){
                        m_momStops_bToF_Plus = a_momStops_bToF_Plus;
                        m_momStops_bToF_Minus = a_momStops_bToF_Minus; };
  void setMomStartStops_eToF(vector<double> a_momStarts_eToF_Plus, vector<double> a_momStarts_eToF_Minus,
                             vector<double> a_momStops_eToF_Plus, vector<double> a_momStops_eToF_Minus ){
                        m_momStarts_eToF_Plus = a_momStarts_eToF_Plus;
                        m_momStarts_eToF_Minus = a_momStarts_eToF_Minus;
                        m_momStops_eToF_Plus = a_momStops_eToF_Plus;
                        m_momStops_eToF_Minus = a_momStops_eToF_Minus; };


  #ifndef _MAC_OSX_
    bool isHLTEvent(StPicoEvent* a_event);
    bool isGoodEvent(StPicoEvent* a_event);
    #ifdef _HAS_ETOF_
      bool isEToFEvent(StPicoEvent* a_event);
    #endif
    bool isGoodTrack(StPicoTrack* a_track, TVector3* a_EventVertex, bool a_ignoreDCA = false);
    bool isGoodBTof(StPicoBTofPidTraits* a_traits);
  #endif

  // Pulled out from the _MAC_OSX_ gate above -- see the comment near the top of
  // this file, next to the StPicoETofPidTraits/StPicoETofHit includes.
  #ifdef _HAS_ETOF_
    bool isGoodETof(StPicoETofPidTraits* a_traits, StPicoETofHit* a_hit = NULL);
  #endif



  bool isGoodEvent(double a_x, double a_y, double a_z, int a_nBTOFMatch, int a_nTofT0 = 999);
  bool isGoodTrack(bool a_isPrimary,int a_nHitsFit, int a_nHitsMax, int a_nHitsDedx, double a_gDCA, double a_phi);
  bool isGoodTrack(int a_nHitsFit, double a_nHitsOverNPoss, int a_nHitsDedx, double a_gDCA, double a_phi);


  void makeCutInformationHistogram();
  TH1D* getCutInformationHistogram(){return m_cutInformation;};

  // call this to determine that you have good dEdx data after you determine it has goodBTOF
  //Uses btof to define betagamma region for precise tracks (Cuts out double pion peaks)
  bool isGoodTrack_GivenBTOFBetaGamma(double a_dEdx, double a_betaGamma);
  bool isGoodTrack_GivenBTOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr);
  bool isElectron_GivenBTOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr);// returns false if toggle is off
  bool isGoodTrack_GivenETOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr);
  bool isElectron_GivenETOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr);// returns false if toggle is off

  bool isGoodBTof(double a_localY, double a_localZ);//, double a_btofBeta, double a_dEdx);
  bool isSpecificParticleByTOF(double a_tofMass, int a_partIndex);
  bool isSpecificParticleByBTOF(double a_invBeta,double a_p, int a_partIndex);
  bool isSpecificParticleByETOF(double a_invBeta,double a_p, int a_partIndex);

  #ifdef _STREFMULTCORR_
    void   setStRefMultCorrToggle(bool a_toggle, bool a_use9bins = true){m_useStRefMultCorr = a_toggle; m_refMultCorr9Bins = a_use9bins;};
    void   setStRefMultCorrCentralities();
    void   setStRefMultCorrRunId(int a_runID){ m_refMultCorr->init(a_runID);};
    void   setStRefMultCorrEvent(StPicoEvent* a_event);
    void   setStRefMultCorrEvent(int a_refMult, double a_Vz, double a_ZDCx);
    bool   isBadStRefMultCorrRunId(int a_runID);
    double getCorrectedRefMult(){return m_refMultCorr->getRefMultCorr();};
    double getCorrectedWeight(){return m_refMultCorr->getWeight();};
  #endif

  bool isPileUp(int a_centVariable, int a_tofmatch);
  //doesn't work with refmult
  bool isLinPileUp(int a_centVariable){if(!m_useStRefMultCorr) return (a_centVariable > m_pileUpCut); else return false;};
  //bool isGoodDeDxBtof(double a_dEdx, double a_betaGama);

  // If using StRefMultCorr, 0 is most central (purpously reversed from standard)
  int    centralityIndex(int a_mult); //returns -1 if not found //a_mult is neglected if StRefMultCorr is defined
  int    getNCentBins(){return m_numCentralities;};
  vector<unsigned int> getTriggers(){return m_triggers;};
  bool   getLogMTBinning(){return m_useLogMTBinning;};
  bool   getVariableMtM0BinEdges(){return m_variableMtM0BinEdges;};
  TH1D*  getEtaPtBinStructure(int a_tag); //0=pt, eta-> 1=tpc, 2=btof, 3=etof
  TH1D*  getRapMtM0BinStructure(int a_particleIndex, int a_tag); //0=mtm0, rap-> 1=tpc, 2=btof, 3=etof 4 = mTm0_btof, 5 = mTm0_etof
  TH1D*  getEventRejectionHistogram(){return m_rejectedEvents;};
  TH1D*  getTrackRejectionHistogram(){return m_rejectedTracks;};
  int    getNRapidityBins(int a_particleIndex, int a_tag);
  double getLowRapidity(int a_particleIndex, int a_tag);
  double getHighRapidity(int a_particleIndex, int a_tag);
  double getRapidityRangeHigh(int a_yIndex, int a_partIndex, int a_tag);
  double getRapidityRangeLow(int a_yIndex, int a_partIndex, int a_tag);
  double getRapidityRangeCenter(int a_yIndex, int a_partIndex, int a_tag);
  int    getN_mTm0Bins(int a_particleIndex);
  double getLow_mTm0(int a_particleIndex);
  double getHigh_mTm0(int a_particleIndex);
  double getLowCentPercent(int a_centIndex){return m_centPercents[a_centIndex];};
  double getHighCentPercent(int a_centIndex){return m_centPercents[a_centIndex+1];};
  double getCentCutLowEdge(int a_centIndex){return m_centEdges.at(a_centIndex+1);};
  double getCentCutHighEdge(int a_centIndex){return m_centEdges.at(a_centIndex);};

  double getLinPileUpCut(){return m_pileUpCut;};
  bool   getUseHuiDT27GevPileup(){return m_useHui_dt_27GeVPileup;};
  float  getYCM(){return m_yCM;};

  string getEventConfig(){return m_eventConfig;};
  string getVertexConfig(){return m_vertexConfig;};
  string getStarver(){return m_starver;};
  //string getPartInfoStarver(){return m_partInfo->GetStarVer();};
  double getSqrts_NN(){ return m_sqrts_NN; };
  double getBeamEnergy(){ return m_beamEnergy; };

  double getDCACut(){ return m_DCA; };
  double getNHitsFitCut(){ return m_NHitsFit; };
  double getNHitsDeDxCut(){ return m_NHitsDeDx; };
  double getFitMaxRatioCut(){ return m_fitMaxRatio; };
  double getNSigmaCut(int a_partIndex){return m_nSigmaCuts[a_partIndex];};


  void getParticleFitMatrix(int a_charge, bool a_particleFitMatrix[9][9]);

  double getLowMassCut(int a_partIndex){return m_lowMassCut.at(a_partIndex);};
  double getHighMassCut(int a_partIndex){return m_highMassCut.at(a_partIndex);};
  double getLowMassSquaredCut(int a_partIndex);
  double getHighMassSquaredCut(int a_partIndex);
  double getMassBTOF(int a_partIndex){return m_particleMasses_bTOF[a_partIndex];};
  double getMassETOF(int a_partIndex){return m_particleMasses_eTOF[a_partIndex];};
  double getMomStops_TPC(int a_partIndex, int a_charge);
  double getMomStops_bToF(int a_partIndex, int a_charge);
  double getMomStart_eToF(int a_partIndex, int a_charge);
  double getMomStop_eToF(int a_partIndex, int a_charge);

  std::pair<double,double> getZVertexCuts();
  double getHighZVertexCut(){return m_ZHighCut;};
  double getLowZVertexCut(){return m_ZLowCut;};
  double getRadialVertexCut(){ return m_radialCut; };
  std::pair<double,double> getBeamSpotLocation();

  double getModNSigma_bTOF(double a_dEdx, double a_invBetaGammaSqr);
  double getModNSigma_eTOF(double a_dEdx, double a_invBetaGammaSqr);
  double getDeDxPrediction(double a_p, double a_m);
  double getDeDxPrediciton_byPart(double a_p,int a_partIndex);
  TF1*   getInvBetaGammaSqrDeDxCalibFunct(){return m_invBetaGammaSqr_dEdx_Calib_funct;};
  //ClassDef(CutClass, 1);

  //Rapidity and Mt-M0 Binning Info
  vector< TH1D* > m_mTm0_structs; // by particle index
  vector< TH1D* > m_tpc_rap_structs;
  vector< TH1D* > m_btof_rap_structs;
  vector< TH1D* > m_etof_rap_structs;
  vector< TH1D* > m_btof_mTm0_structs;
  vector< TH1D* > m_etof_mTm0_structs;




private:

  string m_eventConfig; // FixedTarget or ColliderCenter
  string m_vertexConfig;
  float  m_yCM;
  double m_sqrts_NN;
  double m_beamEnergy;
  string m_starver;
  double m_particleMasses[9];

  //####  BINNING INFORMATION ######
  //ETA and Pt
  TH1D* m_pT_struct;
  TH1D* m_eta_tpc_struct;
  TH1D* m_eta_btof_struct;
  TH1D* m_eta_etof_struct;

  TH1D* m_rejectedEvents; //histogram telling which cut killed your events
  TH1D* m_rejectedTracks; //histogram telling which cut killed your tracks

  TH1D* m_cutInformation; // stores cuts in a histogram so you can have a record of them

  TF1*  m_invBetaGammaSqr_dEdx_Calib_funct;

  int m_numCentralities;
  vector<double> m_centPercents;
  vector<double> m_centEdges;


  #ifdef _STREFMULTCORR_
    StRefMultCorr* m_refMultCorr;
    bool m_refMultCorr9Bins;
  #endif

  // TOGGLES
  bool m_useLogMTBinning;
  bool m_variableMtM0BinEdges;
  bool m_useStRefMultCorr;
  bool m_cutOnTriggers;
  bool m_usePhiCut;
  bool m_useEtaCut;
  bool m_usePtCut;
  bool m_useT0;
  bool m_useHui_dt_27GeVPileup; // specific toggle to copy Hui's analysis






  // EVENT INFO/CUTS
  double m_beamXCenter;
  double m_beamYCenter;
  double m_ZLowCut;     // Vz >= cut
  double m_ZHighCut;    // Vz <= cut
  double m_radialCut;   // Vr <= cut
  int    m_tofMatchCut; // n >= cut
  int    m_nTofT0Cut;   // n >= cut
  double m_pileUpCut;
  double m_pileUpLineCuts[6];


  // TRACK INFO/CUTS
  double m_DCA;
  int m_NHitsFit;   // n >= cut
  int m_NHitsDeDx;  // n >= cut
  double m_fitMaxRatio;

  double m_btofLocY;
  double m_btofLocZ;
  double m_etofDeltaCut; // defines a radius for a circle centered at 0,0 for delta X vs delta Y
  double m_etofDetlaXY_Center_X;
  double m_etofDeltaXY_Center_Y;
  double m_etofDeltaX_Ellipse_Cut;
  double m_etofDeltaY_Ellipse_Cut;
  int    m_etof_max_cluster_size; // n <= cut
  int    m_etof_match_flag_cut; // -1 for > 0 |  1 for single-single  | 2 for single-multiple

  double m_lowPhiCut; //run 18 had itpc not enabled so cut phi
  double m_highPhiCut;
  double m_lowEtaCut;
  double m_highEtaCut;
  double m_lowPtCut;



  vector< double > m_nSigmaCuts;
  vector<unsigned int> m_triggers;
  vector<unsigned int> m_HLTTriggers;
  #ifdef _HAS_ETOF_
    vector<unsigned int> m_EToFTriggers;
  #endif


  vector<double> m_lowMassCut;
  vector<double> m_highMassCut;
  vector<double> m_particleMasses_bTOF;
  vector<double> m_particleMasses_eTOF;
  vector<double> m_invBetaCut_bTOF;
  vector<double> m_invBetaCut_eTOF;


  bool m_particleFitMatrix_Plus[9][9];
  bool m_particleFitMatrix_Minus[9][9];

  bool m_useBTOFBetaGammaCuts;
  vector<double> m_BTOFBetaGammaLowLineSlopes;
  vector<double> m_BTOFBetaGammaLowLineIntercepts;
  vector<double> m_BTOFBetaGammaHighLineXPoints;
  vector<double> m_BTOFBetaGammaHighLineYPoints;

  bool m_useBTOFInvBetaGammaSqrCuts;
  vector<double> m_BTOFInvBetaGammaSqrLowLineSlopes;
  vector<double> m_BTOFInvBetaGammaSqrLowLineIntercepts;
  vector<double> m_BTOFInvBetaGammaSqrHighLineXPoints;
  vector<double> m_BTOFInvBetaGammaSqrHighLineYPoints;

  bool m_useETOFInvBetaGammaSqrCuts;
  vector<double> m_ETOFInvBetaGammaSqrLowLineSlopes;
  vector<double> m_ETOFInvBetaGammaSqrLowLineIntercepts;
  vector<double> m_ETOFInvBetaGammaSqrHighLineXPoints;
  vector<double> m_ETOFInvBetaGammaSqrHighLineYPoints;


  bool   m_useModNSigma_InvBetaGammaSqrBTOFCuts;
  double m_lowInvBetaGammaSqr_ToUseModNSigmaBTOFCut;
  double m_lowModNSigma_BTOF_Cut;
  double m_highModNSigma_BTOF_Cut;
  double m_lowInvBetaGammaSqrBTOFCut;

  bool   m_useBTOFElectronCut;
  double m_electron_lowInvBetaGammaSqrBTOFCut;
  double m_electron_highInvBetaGammaSqrBTOFCut;
  double m_electron_lowDeDxBTOFCut;
  double m_electron_highDeDxBTOFCut;

  bool   m_useModNSigma_InvBetaGammaSqrETOFCuts;
  double m_lowInvBetaGammaSqr_ToUseModNSigmaETOFCut;
  double m_lowModNSigma_ETOF_Cut;
  double m_highModNSigma_ETOF_Cut;
  double m_lowInvBetaGammaSqrETOFCut;

  bool   m_useETOFElectronCut;
  double m_electron_lowInvBetaGammaSqrETOFCut;
  double m_electron_highInvBetaGammaSqrETOFCut;
  double m_electron_lowDeDxETOFCut;
  double m_electron_highDeDxETOFCut;

  vector<double> m_momStops_TPC_Minus;
  vector<double> m_momStops_TPC_Plus;
  vector<double> m_momStops_bToF_Minus;
  vector<double> m_momStops_bToF_Plus;

  vector<double> m_momStarts_eToF_Minus;
  vector<double> m_momStarts_eToF_Plus;
  vector<double> m_momStops_eToF_Minus;
  vector<double> m_momStops_eToF_Plus;



};

#endif
