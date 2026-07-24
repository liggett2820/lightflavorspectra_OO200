
#ifndef _Efficiency_Fitter_
#define _Efficiency_Fitter_

#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TVirtualFitter.h"
#include "TLine.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TVirtualFitter.h"


#include "../makefile_toggles.h"

#include "../headers/PhysMath.h"
#include "../headers/HistogramUtilities.h"
#include "../headers/PhysMathHistogramUtilities.h"
//#include "../headers/CutClass.h"
#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"
#include "../headers/SplineFitter.h"
#include "../headers/GraphFitter.h"

static const int g_EffFit_MaxNumRapBins = 120;
static const int g_EffFit_MaxNumMtM0Bins = 150;
static const int g_EffFit_MaxNumVzBins = 28;

class EfficiencyFitter{
public:
	  EfficiencyFitter();

    //######## OPTIONS ##############
    void setFixedTargetMode(bool a_toggle = true){m_fixed_target_mode = a_toggle;};
    void testingOneCent(bool a_toggle = true){m_testingOneCent = a_toggle;};
    void testingMidRapidity(bool a_toggle = true,double a_ycm = 0.0){m_testingMidRapidity = a_toggle; m_ycm = a_ycm;};
    void setImageDir(string a_imageDir){m_imageDir = a_imageDir;  gSystem->mkdir(a_imageDir.c_str());};
    void setLegendreOrderForTPCEff(int a_partIndex, int a_charge, int a_order);
    void setLegendreOrderForBTOFEff(int a_partIndex, int a_charge, int a_order);
    void setPtLossMode(bool a_toggle = true){m_pTLossMode = a_toggle;};
    void setUseEnergyLossTProfiles(bool a_toggle = true){m_useEnergyLossTProfiles = a_toggle;};
    void setUseIntegralChiSqr(bool a_toggle = true){m_useIntegralChiSqr = a_toggle;};
    void setTpcEff_Systematic_SeedTuple(vector<vector<double>> a_tuple){m_tpcEff_Systematic_SeedTuple = a_tuple;}; // {partIndex,rapidity,p0,p1,p2}
    void setBtofEff_Systematic_SeedTuple(vector<vector<double>> a_tuple){m_btofEff_Systematic_SeedTuple = a_tuple;}; // {partIndex,rapidity,p0,p1,p2}
    void setTpcEff_Legendre_SeedTuples(vector<vector<double>> a_legendre_start_tuple,vector<vector<double>> a_line_start_tuple){m_TPC_Eff_lowMtM0Legendre_Tuple = a_legendre_start_tuple;m_TPC_Eff_highMtM0ToLine_Tuple = a_line_start_tuple;}; // {partIndex,rapidity,p0,p1,p2}
    void setRapBinsToTest(vector<double> a_rapidityValues){m_rapiditiesToTest = a_rapidityValues;};
    void setNoFitNumBins(int a_tpc_num_bins = 80, int a_btof_num_bins = 120){ m_tpc_eff_no_fit_num_bins = a_tpc_num_bins; m_btof_eff_no_fit_num_bins = a_btof_num_bins;};

    //#############   LOADERS   ##############
    //Load the embedding files from the submodule/MuDstProcessEmbedding   TPC Efficiency, Energy Loss, Btof Efficiency // if a_inFileNameTwo is not "", then it adds both plus and minus together
    void loadEmbeddingFiles(string a_inFileName, string a_inFileNameTwo, int a_partIndex, int a_charge, int a_mtm0Rebin = 1, bool a_convert16to9CentBins = false, bool a_savePtLossGausImages = false, bool a_fitEnergyLossByCent = false, int a_energyLossRebin = 1, bool a_doEnergyLossFits = true, vector<double> a_refMultCuts = {});

    void loadEmbeddingFiles_ForVzRecombination(int a_partIndex, int a_charge, vector<string> a_inFileNames,vector<vector<double>> a_VzRanges, string a_picobinner_file_for_Vz_data_histo, vector<double> a_refMultCuts = {});

    //void remap_tpc_eff_byRefMult(string a_inFileName, int a_partIndex, int a_charge, vector<double> a_refMultCuts); // a_refMultCuts should be in order like {293,256,238,...} 0-5% will be 256->293

    //Load the output of PicoBinner for bTOF Efficiency (data_driven)
    void loadPicoBinnerDataFile(string  a_inFileName, int a_partIndex, int a_charge, int a_numPicoBinnerCentBins = 16, int a_mtm0Rebin = 1, bool a_convert16to9CentBins = false, int a_knockout_dca_rebin = 1, int a_knockout_mTm0_rebin = 1);

    void loadHybridFeedDown(string a_filename, int a_partIndex, int a_charge);


    //###############  DATA MANIPULATORS #####################
    void copyAllCentTPCFitsToSpecCent(int a_partIndex, int a_charge);
    void copyAllCentBTOFFitsToSpecCent(int a_partIndex, int a_charge);

    //###############  FITTERS #####################
    void fitEfficiency(int a_partIndex, int a_charge, int a_centIndex = -1, bool a_doTPCEff = false, bool a_doDataDrivenBTOF = true, bool a_doLegendreResiduleFits = true); // for embedding data or data Driven / if you don't know which, do DataDriven
    void fitEnergyLoss(int a_partIndex, int a_charge, int a_centIndex = -1,
          vector<double> a_lowRangeExcludeRegionPt = {}, vector<double> a_highRangeExcludeRegionPt = {},
          vector<double> a_lowRangeExcludeRegionEta = {}, vector<double> a_highRangeExcludeRegionEta = {}, bool a_fixFitFunction_y0_toZero = true ); // centIndex<0 means fit allcent
    void fitFeedDown(int a_partIndex, int a_charge, int a_centIndex = -1); // negative centIndex for all cent fit
    void fitHybridFeedDown(int a_partIndex, int a_charge, int a_centIndex = -1); // negative centIndex for all cent fit
    void fitProtonKnockout(int a_centIndex = -1, double a_maxDCA = 3.0);
    //##############   WRITE THE DATA ##################
    void write(string a_outFileName);



    //public chi squared functions (don't use in macros)
    double legendre_ChiSqr(const double* a_param);
    double getChiSquared_DCA(const double* a_params);
    double protonDCAOperator(double* a_x, double* a_params);
    double TGraphAsymmErrors_ChiSqr(const double* a_param);
    double spline_operator_3(double* a_x, double* a_params);
    double spline_operator_5(double* a_x, double* a_params);
    double spline_operator_5_two_point(double* a_x, double* a_params); // the second is a difference in x locs





    //##############   THINGS YOU PROBABLY SHOULDNT USE ##################

    //LEGACY, USE HYBRID FEED DOWN IF YOU CAN. Load output of urqmd simulation from source/embedding/ProcessFeeddownEmbedding.cxx
    void loadFeedDownFile(string  a_inFileName, int a_partIndex, int a_charge, bool a_convert16to9CentBins);
    void compareTurnOnFunctions(); //Fits TPCxAcceptance with multiple sigmoid functions and determines the best chi^2
    //void fitTPCEff(int a_partIndex, int a_charge, int a_centIndex = -1);//, int* a_numFitPoints = NULL); // numFitPointsShould be an array that can override the number of fit points used for each rapidity index   (LEGACY)

    //this function excludes contamination regions: LEGACY (DONT USE, this is just a study! use fitEfficiency instead)
    void fitBTOFEfficiency(bool a_useEmbedding = true, bool a_defineContaminationRegionToSkip = false, int a_specificFitFunctionIndex = -1);

    //~EfficiencyFitter();
    //void fitTPCEff_Chris(int a_partIndex);







private:
  unsigned long m_uniqueIdentifier;
  //int m_numCentBins;
  //vector< bool > m_particleLoaded;
  bool m_particleLoaded_embedding[9][2];
  bool m_particleLoaded_btofEff[9][2];
  bool m_particleLoaded_feeddown[9][2];

  int m_polynomialOrderForTPCEff[9][2];
  int m_polynomialOrderForBTOFEff[9][2];

  string m_imageDir;

  ParticleInfo* m_partInfo;
  //CutClass* m_cuts;

  double m_ycm;
  bool m_fixed_target_mode;
  bool m_testingOneCent;
  bool m_testingMidRapidity;
  bool m_makeChrisConfidenceGraphs;
  bool m_pTLossMode;
  bool m_useEnergyLossTProfiles;
  bool m_useIntegralChiSqr;
  bool m_useChebyshevSpacing;
  vector<double> m_rapiditiesToTest;
  int  m_tpc_eff_no_fit_num_bins;
  int  m_btof_eff_no_fit_num_bins;


  TF1*               m_currentFitFunction;
  TGraphAsymmErrors* m_currentGraphToFit;
  TH1D* m_currentProtonDCA;
  TH1D* m_currentAntiProtonDCA;
  int m_currentDCALowerFitRangeBin;
  int m_currentDCAUpperFitRangeBin;
  int m_currentRapIndex;
  int m_currentCentIndex;

  int       m_spline_numKnots;
  double    m_splineParams[100];
  TSpline3* m_spline;
  TSpline5* m_spline5;


  vector< vector<double>> m_tpcEff_Systematic_SeedTuple; // {partIndex,rapidity,p0,p1,p2}
  vector< vector<double>> m_btofEff_Systematic_SeedTuple; // {partIndex,rapidity,p0,p1,p2}

  TH1D* m_knockout_mTm0_structure;
  TH1D* m_picoBinner_rapidity_structure[9][3];    // tpc btof etof
  TH1D* m_embedding_rapidity_structure[9];        // only one
  TH1D* m_embedding_pseudorapidity_structure[9];  // only one
  TH1D* m_hybridFeedDown_rapidity_structure[9];   // only one
  //TH1D* m_embedding_vZ_structure;

  //##################   Embedding Data for TPC Efficiency #####################
  //#### DATA TO FIT
  TGraphAsymmErrors* m_TPC_Efficiency_AllCent[9][g_EffFit_MaxNumRapBins][2];    //part rapidity charge {+1,-1}
  TGraphAsymmErrors* m_TPC_Efficiency_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TGraphAsymmErrors* m_TPC_Efficiency_AllCent_EtaPt[9][g_EffFit_MaxNumRapBins][2];    //part pseudorapidity charge {+1,-1}
  TGraphAsymmErrors* m_TPC_Efficiency_ByCent_EtaPt[9][16][g_EffFit_MaxNumRapBins][2]; //part cent pseudorapidity charge {+1,-1}
  //TGraph2DErrors* m_TPC_Efficiency_ByCent[9][16][2]; //part cent charge {+1,-1}
  //#### CHI^2 HISTOGRAMS ####
  TH2D* m_TPC_Efficiency_ChiSqrd_Compare[9][2];
  TH2D* m_TPC_Efficiency_ChiSqrd[9][2]; // part charge {+1,-1} TH2D is x=rapidity y= centrality 17th bin is AllCent
  TH2D* m_TPC_Efficiency_ChiSqrd_Sys[9][2]; // part charge {+1,-1}
  //#### NOMINAL
  TF1* m_TPC_Eff_Fits_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_TPC_Eff_Fits_AllCent[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TF1* m_TPC_Eff_Fits_HighPt_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_TPC_Eff_Fits_HighPt_AllCent[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TGraphErrors* m_TPC_Eff_FitErrGraphs_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TGraphErrors* m_TPC_Eff_FitErrGraphs_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_TPC_Eff_CovMat_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_TPC_Eff_CovMat_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_TPC_Eff_CovMat_HighPt_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_TPC_Eff_CovMat_HighPt_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part rap charge
  //#### SYSTEMATIC
  TF1* m_TPC_Eff_Fits_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_TPC_Eff_Fits_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TGraphErrors* m_TPC_Eff_FitErrGraphs_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TGraphErrors* m_TPC_Eff_FitErrGraphs_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_TPC_Eff_CovMat_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_TPC_Eff_CovMat_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part cent rap

  TF1* m_TPC_Eff_highMtM0ToLine[9]; // the mtm0 cutoff for switching to line fit

  vector< vector<double>> m_TPC_Eff_lowMtM0Legendre_Tuple; // the mtm0 cutoff for starting Legendre fit // {partIndex,rapidity,mTm0}
  vector< vector<double>> m_TPC_Eff_highMtM0ToLine_Tuple; // the mtm0 cutoff for switching to line fit // {partIndex,rapidity,mTm0}

  //##################   Embedding Data for BTOF Efficiency #####################
  //#### DATA TO FIT
  TGraphAsymmErrors* m_BTOF_Efficiency_AllCent[9][g_EffFit_MaxNumRapBins][2];    //part rapidity charge {+1,-1}
  TGraphAsymmErrors* m_BTOF_Efficiency_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  //#### CHI^2 HISTOGRAMS ####
  TH2D* m_BTOF_Efficiency_ChiSqrd_Compare[9][2];
  TH2D* m_BTOF_Efficiency_ChiSqrd[9][2]; // part charge {+1,-1} TH2D is x=rapidity y= centrality 17th bin is AllCent
  TH2D* m_BTOF_Efficiency_ChiSqrd_Sys[9][2]; // part charge {+1,-1}
  //#### NOMINAL
  TF1* m_BTOF_Eff_Fits_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_BTOF_Eff_Fits_AllCent[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TMatrixDSym* m_BTOF_Eff_CovMat_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_BTOF_Eff_CovMat_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  //#### SYSTEMATIC
  TF1* m_BTOF_Eff_Fits_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_BTOF_Eff_Fits_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TMatrixDSym* m_BTOF_Eff_CovMat_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_BTOF_Eff_CovMat_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part rap

  //##################   Embedding Data for Energy Loss #####################
  //### DATA TO FIT
  TGraphAsymmErrors* m_EnergyLoss_AllCent[9][g_EffFit_MaxNumRapBins][2];    //part rapidity charge {+1,-1}
  TGraphAsymmErrors* m_EnergyLoss_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TProfile*          m_EnergyLossProfile_AllCent[9][g_EffFit_MaxNumRapBins][2];    //part rapidity charge {+1,-1}
  TProfile*          m_EnergyLossProfile_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  //#### CHI^2 HISTOGRAMS ####
  TH2D* m_EnergyLoss_ChiSqrd[9][2]; // part charge {+1,-1} TH2D is x=rapidity y= centrality 17th bin is AllCent
  TH2D* m_EnergyLoss_ChiSqrd_Sys[9][2]; // part charge {+1,-1}
  //#### NOMINAL
  TF1* m_EnergyLoss_Fits_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_EnergyLoss_Fits_AllCent[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TGraphErrors* m_EnergyLoss_FitErrGraphs_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TGraphErrors* m_EnergyLoss_FitErrGraphs_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_EnergyLoss_CovMat_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_EnergyLoss_CovMat_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  //#### SYSTEMATIC
  TF1* m_EnergyLoss_Fits_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_EnergyLoss_Fits_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TGraphErrors* m_EnergyLoss_FitErrGraphs_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TGraphErrors* m_EnergyLoss_FitErrGraphs_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_EnergyLoss_CovMat_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym*  m_EnergyLoss_CovMat_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part cent rap

  TH3I* m_mTm0ReponseMatrixHisto[9][2]; // Centrality Integrated
  TH3I* m_mTm0ReponseMatrixHisto_BTOF[9][2]; // Centrality Integrated
  TH3D* m_mTm0ReponseMatrixHisto_PDF[9][2]; // Centrality Integrated
  TH3D* m_mTm0ReponseMatrixHisto_BTOF_PDF[9][2]; // Centrality Integrated




  //##################   FeedDown Data / Pure URQMD #####################
  TGraphAsymmErrors* m_FeedDown_primaryOverTotal_AllCent[9][g_EffFit_MaxNumRapBins][2]; //part rapidity charge {+1,-1}
  TGraphAsymmErrors* m_FeedDown_primaryOverTotal_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TGraphAsymmErrors* m_FeedDown_feeddownOverTotal_AllCent[9][g_EffFit_MaxNumRapBins][2]; //part rapidity charge {+1,-1}
  TGraphAsymmErrors* m_FeedDown_feeddownOverTotal_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  //#### CHI^2 HISTOGRAMS ####
  TH2D* m_FeedDown_ChiSqrd[9][2]; // part charge {+1,-1} TH2D is x=rapidity y= centrality 17th bin is AllCent
  TH2D* m_FeedDown_ChiSqrd_Sys[9][2]; // part charge {+1,-1}
  //#### NOMINAL
  TF1* m_FeedDown_Fits_AllCent[9][g_EffFit_MaxNumRapBins][2]; //part rapidity charge {+1,-1}
  TF1* m_FeedDown_Fits_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TMatrixDSym* m_FeedDown_CovMat_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_FeedDown_CovMat_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  //#### SYSTEMATIC
  TF1* m_FeedDown_Fits_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; //part rapidity charge {+1,-1}
  TF1* m_FeedDown_Fits_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TMatrixDSym* m_FeedDown_CovMat_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_FeedDown_CovMat_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part cent rap



  //##################   FeedDown Data / Hybrid Method #####################
  //TGraphAsymmErrors* m_FeedDown_Hybrid_AllCent[9][g_EffFit_MaxNumRapBins][2]; //part rapidity charge {+1,-1}
  TH1D* m_FeedDown_Hybrid_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TGraphErrors* m_FeedDown_Hybrid_ByCent_graphs[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}

  //#### CHI^2 HISTOGRAMS ####
  //TH2D* m_FeedDown_Hybrid_ChiSqrd[9][2]; // part charge {+1,-1} TH2D is x=rapidity y= centrality 17th bin is AllCent
  //TH2D* m_FeedDown_Hybrid_ChiSqrd_Sys[9][2]; // part charge {+1,-1}
  //#### NOMINAL
  //TF1* m_FeedDown_Hybrid_Fits_AllCent[9][g_EffFit_MaxNumRapBins][2]; //part rapidity charge {+1,-1}
  TF1* m_FeedDown_Hybrid_Fits_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TF1* m_FeedDown_Hybrid_Fits_ByCent_lowPt[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  TMatrixDSym* m_FeedDown_Hybrid_CovMat_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_FeedDown_Hybrid_CovMat_ByCent_lowPt[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap

  //TMatrixDSym* m_FeedDown_Hybrid_CovMat_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  //#### SYSTEMATIC
  //TF1* m_FeedDown_Hybrid_Fits_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; //part rapidity charge {+1,-1}
  //TF1* m_FeedDown_Hybrid_Fits_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  //TMatrixDSym* m_FeedDown_Hybrid_CovMat_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  //TMatrixDSym* m_FeedDown_Hybrid_CovMat_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part cent rap


  //##################   Data for BTOF Efficiency DATA DRIVEN #####################
  //#### DATA TO FIT
  TGraphAsymmErrors* m_BTOF_fromData_Efficiency_AllCent[9][g_EffFit_MaxNumRapBins][2];    //part rapidity charge {+1,-1}
  TGraphAsymmErrors* m_BTOF_fromData_Efficiency_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; //part cent rapidity charge {+1,-1}
  //#### CHI^2 HISTOGRAMS ####
  TH2D* m_BTOF_fromData_Efficiency_ChiSqrd_Compare[9][2];
  TH2D* m_BTOF_fromData_Efficiency_ChiSqrd[9][2]; // part charge {+1,-1} TH2D is x=rapidity y= centrality 17th bin is AllCent
  TH2D* m_BTOF_fromData_Efficiency_ChiSqrd_Sys[9][2]; // part charge {+1,-1}
  //#### NOMINAL
  TF1* m_BTOF_fromData_Eff_Fits_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_BTOF_fromData_Eff_Fits_AllCent[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TF1* m_BTOF_fromData_Eff_Fits_HighPt_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_BTOF_fromData_Eff_Fits_HighPt_AllCent[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TMatrixDSym* m_BTOF_fromData_Eff_CovMat_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_BTOF_fromData_Eff_CovMat_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_BTOF_fromData_Eff_CovMat_HighPt_ByCent[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_BTOF_fromData_Eff_CovMat_HighPt_AllCent[9][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TF1* m_BTOF_highMtM0ToLine[9];
  //#### SYSTEMATIC
  TF1* m_BTOF_fromData_Eff_Fits_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // particle centrality rapdiity
  TF1* m_BTOF_fromData_Eff_Fits_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // particle rapdiity
  TMatrixDSym* m_BTOF_fromData_Eff_CovMat_ByCent_Sys[9][16][g_EffFit_MaxNumRapBins][2]; // part cent rap
  TMatrixDSym* m_BTOF_fromData_Eff_CovMat_AllCent_Sys[9][g_EffFit_MaxNumRapBins][2]; // part cent rap


  //##################   Data for Proton DCA Knockout #####################
  //#### DATA TO FIT    (knockout is fraction as funct of mTm0    |    dca is dca fit for individual y & mTm0 fit)
  TGraphAsymmErrors* m_proton_knockout_AllCent[g_EffFit_MaxNumRapBins];    //rapidity
  TGraphAsymmErrors* m_proton_knockout_ByCent[16][g_EffFit_MaxNumRapBins]; //cent rapidity
  TH1D* m_proton_dca_AllCent[g_EffFit_MaxNumRapBins][g_EffFit_MaxNumMtM0Bins][2];    //rapidity charge {+1,-1}
  TH1D* m_proton_dca_ByCent[16][g_EffFit_MaxNumRapBins][g_EffFit_MaxNumMtM0Bins][2]; //cent rapidity charge {+1,-1}
  //#### CHI^2 HISTOGRAMS ####
  TH2D* m_proton_knockout_ChiSqrd; // TH2D is x=rapidity y= centrality 17th bin is AllCent
  TH2D* m_proton_dca_ChiSqrd[16]; //  TH2D is x=rapidity y= mTm0
  //#### NOMINAL
  TF1* m_proton_knockout_Fits_ByCent[16][g_EffFit_MaxNumRapBins]; // particle centrality rapdiity
  TF1* m_proton_knockout_Fits_AllCent[g_EffFit_MaxNumRapBins]; // particle rapdiity
  TF1* m_proton_dca_Fits_ByCent[16][g_EffFit_MaxNumRapBins][g_EffFit_MaxNumMtM0Bins]; // particle centrality rapdiity
  TF1* m_proton_dca_Fits_AllCent[g_EffFit_MaxNumRapBins][g_EffFit_MaxNumMtM0Bins]; // particle rapdiity
  TMatrixDSym* m_proton_knockout_CovMat_ByCent[16][g_EffFit_MaxNumRapBins]; // part cent rap
  TH2D* m_proton_knockout_histogram_ByCent[16]; //  TH2D is x=rapidity y= mTm0


  //##################   BTOF Efficiency Data  #####################
  TH1D* m_centEvents[9];

  //BTOF Efficiencies calculated by Data Driven Methods (nSigma cut and ratio)
  TH2D* m_BTOFEff_Data_AllCent_Plus[9];
  TH2D* m_BTOFEff_Data_AllCent_Minus[9];
  TH2D* m_BTOFEff_Data_AllCent_Combined[9];
  TH2D* m_BTOFEff_Data_Plus[9][16];
  TH2D* m_BTOFEff_Data_Minus[9][16];
  TH2D* m_BTOFEff_Data_Combined[9][16];
  TH2D* m_tpcTrack_For_BTOFEff_Data_Plus[9][16];
  TH2D* m_bTOFTrack_For_BTOFEff_Data_Plus[9][16];
  TH2D* m_tpcTrack_For_BTOFEff_Data_Minus[9][16];
  TH2D* m_bTOFTrack_For_BTOFEff_Data_Minus[9][16];

  TH1D* m_BTOFEff_Data_FitParams_Plus[9][16][8][6];//part cent fit param
  TH1D* m_BTOFEff_Data_FitParams_Minus[9][16][8][6];
  TH1D* m_BTOFEff_Data_FitParams_Combined[9][16][8][6];
  TH1D* m_BTOFEff_Data_FitParams_AllCent_Plus[9][8][6];
  TH1D* m_BTOFEff_Data_FitParams_AllCent_Minus[9][8][6];
  TH1D* m_BTOFEff_Data_FitParams_AllCent_Combined[9][8][6];

  TH2D*   m_BTOFEff_Data_ChiSqrd_AllCent_Plus[9]; // part rap
  TH2D*   m_BTOFEff_Data_ChiSqrd_AllCent_Minus[9]; // part rap
  TH2D*   m_BTOFEff_Data_ChiSqrd_AllCent_Combined[9]; // part rap

  //BTOF Efficiencies calculated by Embedding Samples
  TH2D* m_BTOFEff_Emb_AllCent_Plus[9];
  TH2D* m_BTOFEff_Emb_AllCent_Minus[9];
  TH2D* m_BTOFEff_Emb_AllCent_Combined[9];
  TH2D* m_BTOFEff_Emb_Plus[9][16];
  TH2D* m_BTOFEff_Emb_Minus[9][16];
  TH2D* m_BTOFEff_Emb_Combined[9][16];
  TH2D* m_tpcTrack_For_BTOFEff_Emb_Plus[9][16];
  TH2D* m_bTOFTrack_For_BTOFEff_Emb_Plus[9][16];
  TH2D* m_tpcTrack_For_BTOFEff_Emb_Minus[9][16];
  TH2D* m_bTOFTrack_For_BTOFEff_Emb_Minus[9][16];

  TH1D* m_BTOFEff_Emb_FitParams_Plus[9][16][8][6];//part cent fit param
  TH1D* m_BTOFEff_Emb_FitParams_Minus[9][16][8][6];
  TH1D* m_BTOFEff_Emb_FitParams_Combined[9][16][8][6];
  TH1D* m_BTOFEff_Emb_FitParams_AllCent_Plus[9][8][6];
  TH1D* m_BTOFEff_Emb_FitParams_AllCent_Minus[9][8][6];
  TH1D* m_BTOFEff_Emb_FitParams_AllCent_Combined[9][8][6];

  TH2D*   m_BTOFEff_Emb_ChiSqrd_AllCent_Plus[9]; // part rap
  TH2D*   m_BTOFEff_Emb_ChiSqrd_AllCent_Minus[9]; // part rap
  TH2D*   m_BTOFEff_Emb_ChiSqrd_AllCent_Combined[9]; // part rap


  vector< SplineFitter* > m_splinesToDelete;


  double btofEffOperator(double* a_x, const double* a_params);
  //void   makeTSpline3(string a_name, int a_numPoints, vector<double> a_x_knots, double a_low_deriv, double a_high_deriv);
  double m_btofHighMtM0Ignore;
  double m_btofLowMtM0Ignore;

  unsigned long getUniqueID(){return m_uniqueIdentifier++;};
  TF1* fitKnockoutTrendHistogram(int a_paramIndex,TH1D* a_histo);

  TF1* m_functToFitMinimizer;
  TGraphAsymmErrors* m_graphToFitMinimizer;




};






#endif
