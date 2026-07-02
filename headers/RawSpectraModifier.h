#ifndef _RAWSPECTRAMODIFIER_
#define _RAWSPECTRAMODIFIER_

#include <iostream>
#include <string>
#include <vector>
//#include "../source/CutClass.cxx"
#include "../headers/PhysMath.h"
#include "../headers/HistogramUtilities.h"
#include "../headers/CutClass.h"
#include "../headers/SpectraStructure.h"
#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"
#include "../headers/SplineFitter.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualFitter.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLine.h"

using namespace std;

static const int g_RawSpecMod_MaxNumRapBins = 120;


class RawSpectraModifier{
public:

  RawSpectraModifier(); 

  //output of picoBinner
  void loadDataFile(string a_filename, int a_numPart, int a_numCent, bool a_converted16to9 = false, bool a_combine5and10CentBins = false, bool combine40and50_60and70_Cent = false, bool a_useETOF = false, bool a_populateJankSpectra = false, double a_tpcEffJank = 1.0);
  //output of ZFitter
  void loadSpectraFile(string a_filename);

  void loadAndApplyTPCEffAndEnergyLossAndBTOFEffFile(string a_filename, bool a_useTPCEffSystematic = false, bool a_useBTOFEffSystematic = false, bool a_useEnergyLossSystematic = false, int* m_offsets = NULL, bool a_useCent16ForCent9 = false, 
             bool a_useTPCEffDataPoints = false,bool a_useBTOFEffDataPoints = false, bool a_doEnergyLoss = true, bool a_doTPCEff = true, bool a_doBTOFEff = true, bool a_useEnergyLossErrorAsXErrorBars = false);

  void applyTPCEffConf(TGraphErrors* spectraGraph, TGraphAsymmErrors* effGr, TF1* effFunct);

  void loadAndApplyFeedDownCorrections(string a_filename, bool a_useSystematic = false, int a_backwardsCentNum = 0, int* m_offsets = NULL); //leave last to zero

  // Use the HybridFeedDownMaker output to correct for feeddown (this is a subtraction instead of a feed down fraciton)
  void loadAndApplyHybridFeedDownCorrections(string a_filename, bool a_useSinglePointInsteadOfIntegral = false); // you boolean only for testing
  
  void loadAndApplyKnockoutCorrections(string a_filename, int a_backwardsCentNum = 0, int m_offset = -999); //leave last to zero

  //void setCuts(CutClass* a_cuts){m_cuts = a_cuts;};
  void setStarVer(string a_starVer);
  void setChrisCorrectionMode(bool a_toggle = true){m_chrisCorrectionFileMode = a_toggle;};
  void setPtEnergyLossMode(bool a_toggle = true){m_pTLossMode = a_toggle;};
  //void setEqualCountsRebinForTPC(int a_numBins){m_equalCountsRebinForTPC = a_numBins;};
  void setDynamicRebinForJankSpectra(vector< vector<double> > a_dynamicRebin_Plus,  vector< vector<double> > a_dynamicRebin_Minus);  // partIndex detIndex numRebin lowmTm0Range highmTm0Range

  int convertRapidityIndex(int a_yIndexSource, int a_sourceDetectorIndex, int a_targetDetectorIndex, int a_partIndex);

  // this approriatly scales the counts to be invariant yields instead of raw counts
  //Only works on the TH2Ds
  void convertSpectraToInvariant(double a_efficiency = 1.0, //like if you do a global phi cut
                                 double a_effErr = 0.0);

  //NOTE: All this bToF stuff is legacy, use the Efficiency Fitter instead
  void makeBTofEff(string a_imageDir,double* a_momMatch,double* a_momStop, bool a_combinePlusMinus = true, bool a_usePlusEffForMinus = false);
  //Only works on the TH2Ds
  //fit tags are 0=normal, 1=combine plus and minus, 2=use plus for minus
  void fitBTofEff(string a_imageDir,double* a_momMatch,double* a_momStop, bool a_reduceFitRange = false, vector< int > a_fitTags = {0,0,0,0,0,0,0,0,0});
  //Only works on the TH2Ds
  void applyBTOFEff(bool a_useFit = true);

  //Only works on the TH2Ds
  //void scaleBTOF(double* a_momMatch, bool m_constrainScaling = true);

  void makeRapiditySpectraTGraphs();

  //not enough tpc info to really do this
  void fitETofEff(string a_imageDir,double* a_momMatch,double* a_momStop);
  

  void matchBTOFtoTPC(double* a_momMatch, bool a_matchBinByBin, double* a_lowRapToFit, double* a_highRapToFit);

  //This only works on the TGraphs
  //momentum cuts for matching
  void cleanSpectra(double* a_momStart_tpc,
                    double* a_momStart_btof, 
                    double* a_momStart_etof,
                    double* a_momStop_tpc,
                    double* a_momStop_btof, 
                    double* a_momStop_etof,
                    double* a_lowEtaCutTPC,
                    double* a_highEtaCutTPC,
                    double* a_lowEtaCutBTOF,
                    double* a_highEtaCutBTOF,
                    double* a_minMtM0,
                    double* a_maxMtM0,
                    double a_maxErrorFraction,
                    double* a_minRapidity,
                    double* a_maxRapidity,
                    bool a_removeByPreviousRatio,
                    vector<vector<double>> a_btof_mtm0_start_matrix);
                    

  //This is in addition to cleanSpectra, but more specific
  //TRIM DATA: use min and max mTm0 to delete beginning and end points for a specific spectrum
  //  first layer of vector is how many rap spectra you would like to modify
  //  second layer of vector is data on how to modify with this layout
  //        part centIndex charge rapidityValue minMtM0 maxMtM0
  //POINT BY POINT DATA: use  mTm0 and tolerence to delete a point within a specific spectrum
  //  first layer of vector is how many points you would like to modify
  //  second layer of vector is data on how to modify with this layout
  //        part centIndex charge rapidityValue mTm0Val

  void trimSpectra(vector< vector <double> > a_trimData, 
                   vector< vector<double> > a_pointByPointData, 
                   double a_mTm0Tolerence = 0.005); // tolerence will delete a_ponitByPoint with this precision

  void cleanSpectraFXT();


  //LEGACY
  // run everything before this function
  // this uses 2D delayney triangle interpolation to get 
  // rapidity mT-m0 from eta pT spectra
  void makeRapiditySpectraFromEtaSpectra(string a_filename, string a_option = "UPDATE"); // this will write the spectra to a file, you can choose "recreate" or "update" as options
  
  //LEGACY //if in rapidity space already, use this function
  //void makeRapiditySpectra(string a_filename, string a_option = "UPDATE");

  //run this before the final step of hybrid-feeddown correction to copy the current corrected spectra over to the inclusive spectra
  void populateInclusiveSpectra();


  //writes spectra to a file with raw or inclusive options
  void writeRapiditySpectra(string a_filename, string a_option = "UPDATE"); // inclusive is to be done before hybrid feeddown correction

  //writes modified copy to new root file
  void writeMemberData(string a_filename, string a_option = "UPDATE");



private:

  bool m_particleLoaded[9];
  bool m_loadedETOF;
  bool m_combined16To9;
  int m_numCentralities;
  int m_numParticles; 
  int m_numEtaBins_TPC[9];
  int m_numPtBins_TPC[9];
  int m_numEtaBins_BTOF[9];
  int m_numPtBins_BTOF[9];
  int m_numEtaBins_ETOF[9];
  int m_numPtBins_ETOF[9];

  bool m_BTOF_Eff_Calc;
  bool m_ETOF_Eff_Calc;
  int m_tofIndexOffset[9]; // yIndex TOF + offset = yIndex TPC
  int m_etofIndexOffset[9]; // yIndex ETOF + offset

  bool m_chrisCorrectionFileMode;
  //Chris's correction files have energy loss corrections as pT^reco - pT^emb in pT space
  //Our correction is mT^emb - mT^reco in mTm0 space

  bool m_pTLossMode; // Matt's pT^emb - pT^reco. this does not need to be true if chris correction mode
  int m_equalCountsRebinForTPC;
  string m_starVer;
  ParticleInfo* m_partInfo;
  //CutClass* m_cuts;
  vector< vector<double> > m_dynamicRebin_Plus; // for jank spectra
  vector< vector<double> > m_dynamicRebin_Minus; // for jank spectra


  TH1D* m_picoBinner_rapidity_structure[9][3]; // tpc btof etof
  TH1D* m_embedding_rapidity_structure[9];        
  TH1D* m_hybridFeedDown_rapidity_structure[9];   





  TH2D* m_BTOFEffPlusRatio[9][16];
  TH2D* m_BTOFEffMinusRatio[9][16];
  TH2D* m_BTOFEffPlus[9][16];
  TH2D* m_BTOFEffMinus[9][16];
  TH2D* m_ETOFEffPlus[9][16];
  TH2D* m_ETOFEffMinus[9][16];

  //the btof spectra are scaled to the dedx spectra at a given momentum, this stores that scale
  TH1D*  m_btofScaleHistoPlus[9][16]; //part cent
  TH1D*  m_btofScaleHistoMinus[9][16]; //part cent
  vector< double > m_matchMom_History; // by particle , this is copied from the last time someone uses scaleBTOF(double* a_matchMom)
  double m_momStart_tpc[9];
  double m_momStart_btof[9];
  double m_momStart_etof[9];
  double m_momStop_tpc[9];
  double m_momStop_btof[9];
  double m_momStop_etof[9];


  //order: part cent rap
  TGraphErrors* m_spectraPlusDeDx[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors* m_spectraMinusDeDx[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors* m_spectraPlusInvBeta[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors* m_spectraMinusInvBeta[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors* m_spectraPlusInvBetaETOF[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors* m_spectraMinusInvBetaETOF[9][16][g_RawSpecMod_MaxNumRapBins];

  //order: part cent rap
  TGraphErrors*  m_inclusiveSpectraPlusDeDx[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_inclusiveSpectraMinusDeDx[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_inclusiveSpectraPlusInvBeta[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_inclusiveSpectraMinusInvBeta[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_inclusiveSpectraPlusInvBetaETOF[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_inclusiveSpectraMinusInvBetaETOF[9][16][g_RawSpecMod_MaxNumRapBins];

  //order: part cent rap
  TGraphErrors*  m_rawSpectraPlusDeDx[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_rawSpectraMinusDeDx[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_rawSpectraPlusInvBeta[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_rawSpectraMinusInvBeta[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_rawSpectraPlusInvBetaETOF[9][16][g_RawSpecMod_MaxNumRapBins];
  TGraphErrors*  m_rawSpectraMinusInvBetaETOF[9][16][g_RawSpecMod_MaxNumRapBins];
  
  //order: part cent rap (+/-)    The spectra from feeddown, spectra*feeddownfraction
  TGraphErrors* m_feeddownSpectraDeDx[9][16][g_RawSpecMod_MaxNumRapBins][2];
  TGraphErrors* m_feeddownSpectraInvBeta[9][16][g_RawSpecMod_MaxNumRapBins][2];
  TGraphErrors* m_feeddownSpectraInvBetaETOF[9][16][g_RawSpecMod_MaxNumRapBins][2];

  SpectraStructure* m_spectraStruct[9][16][g_RawSpecMod_MaxNumRapBins][2]; // particle centrality rapidity charge



  //#################     LOADED DATA   #####################
  TH1D* m_centEvents[9];
  TH3D* m_rap_mTm0_dEdx_allCent[9];
  TH3D* m_rap_mTm0_invBetaBTOF_allCent[9];
  TH3D* m_rap_mTm0_invBetaETOF_allCent[9];

  //for Tof Eff// These are used for getting binning
  //order: part cent
  TH2D* m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[9][16];
  TH2D* m_rap_mTm0_NumTPCTracks_ForETOF_Plus[9][16];
  TH2D* m_rap_mTm0_NumBTOFTracks_Plus[9][16];
  TH2D* m_rap_mTm0_NumETOFTracks_Plus[9][16];
  TH2D* m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[9][16];
  TH2D* m_rap_mTm0_NumTPCTracks_ForETOF_Minus[9][16];
  TH2D* m_rap_mTm0_NumBTOFTracks_Minus[9][16];
  TH2D* m_rap_mTm0_NumETOFTracks_Minus[9][16];

  //member variables to fill with use of class functions
  TH2D*  m_dEdxSpectraPlus[9][16]; //part cent
  TH2D*  m_InvBetaBTOFSpectraPlus[9][16]; //part cent
  TH2D*  m_InvBetaETOFSpectraPlus[9][16]; //part cent

  TH2D*   m_dEdxSpectraMinus[9][16]; //part cent
  TH2D*   m_InvBetaBTOFSpectraMinus[9][16]; //part cent
  TH2D*   m_InvBetaETOFSpectraMinus[9][16]; //part cent

  TH2D* m_TPCEfficiency_Nominal_Plus[9][16]; //part cent
  TH2D* m_TPCEfficiency_Systematic_Plus[9][16]; //part cent
  TH2D* m_EnergyLoss_Nominal_Plus[9]; //part
  TH2D* m_EnergyLoss_Systematic_Plus[9]; //part

  TH2D* m_TPCEfficiency_Nominal_Minus[9][16]; //part cent
  TH2D* m_TPCEfficiency_Systematic_Minus[9][16]; //part cent
  TH2D* m_EnergyLoss_Nominal_Minus[9]; //part
  TH2D* m_EnergyLoss_Systematic_Minus[9]; //part




  double m_btofLowMtM0Ignore;
  double m_btofHighMtM0Ignore;
  TH1D* m_currentBTOF_Eff_ToFit;
  double btofEffOperator(double* a_x, const double* a_params);



  SpectraStructure* findAssociatedSpectraStructure(int a_partIndex, int a_centIndex, int a_charge, double a_rapidity, double a_rapidityTolerence = 0.0025); // the center should be within tolerence


};




#endif
