#ifndef _SPECTRA_FITTER_
#define _SPECTRA_FITTER_

#include <string>
#include <vector>
#include <iostream>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TVirtualFitter.h"
#include <TCanvas.h>
#include <TString.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLine.h>


#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"
#include "../headers/FitUtilities.h"
#include "../headers/CutClass.h"
#include "../headers/PhysMath.h"
using namespace std;

static const int g_SpecFit_MaxNumRapBins = 35;



class SpectraFitter{
public:

  //#######################     SETUP   ######################################
  SpectraFitter();
  //~SpectraFitter();
  
  void readSpectraFile(string a_filename, int a_partIndex, TString TPC_TOF, int a_numCent = -1, bool a_corrected = true, bool a_inclusive = false); // inclusive has priority over corrected if true
  void readCuts(CutClass* a_cuts){m_cuts = a_cuts;};
  void setImageDir(string a_imageDir);
  //void setUseInclusiveSpectra(bool a_toggle = true){m_useInclusiveSpectra = a_toggle;};
  void setTemperatureSeed(double a_temp){m_temperatureSeed = a_temp;};
  void setLowTemperatureSeed(double a_temp){m_lowTemperatureSeed = a_temp;};
  void setBlastWaveSeeds(double a_temp, double a_betaS){m_blastWave_temperature_seed = a_temp; m_blastWave_betaS_seed = a_betaS;};
  void setDoulbeBoseForPion(bool a_toggle = true){ m_useDoubleBoseForPions = a_toggle; };
  void setUseProtonBetaSForAntiprotons(bool a_toggle = true){m_constrainAntiProtonBetaS_withParamsFromProtons = a_toggle;};
  void setBlastWaveNFix(double a_nValue = 0.5){
    m_blastWave_nFix = a_nValue;
    if(m_blastWave_nFix > 0) m_blastWave->FixParameter(3,m_blastWave_nFix);
  };
  void setVC(vector<double> a_VC){m_VC = a_VC;};
  void setOldSpectraLoadingConvention(bool a_toggle = true){m_oldSpectraLoadingConvention = a_toggle;};
  
  //void setLegendreTemperatureOrder(int a_partIndex, int a_charge, int a_order){ m_legendreTemperatureOrder[a_partIndex][a_charge > 0 ? 0 : 1] = a_order;}; // if the constraint is legendre order dependent
  //void setLegendreParTwoOrder(int a_partIndex, int a_charge, int a_order){ m_legendreParTwoOrder[a_partIndex][a_charge > 0 ? 0 : 1] = a_order;};


  // This function fits a single spectra with a fit function from the list. You can choose to initialize it / fix parameters, or let it initialize the parameters for you.
  // Options for fit function names:     BoseEinstein     PtExpo    Tsallis     MtM0Expo   Boltzmann   BlastWave    Thermal
  // Note: This function will change m_isTemperatureConstrained->False  m_temperatureSeed-> -1.0 if you do not supply parameters
  //       Set low limit to high limit to fix a parameter
  //       Right now this uses EX0, and doesn't use an integral chi^2
  void fitSingleRapiditySpectra(int a_partIndex, 
                                int a_charge, 
                                int a_centIndex,
                                int a_rapidityIndex, 
                                string a_functionName, 
                                bool a_isNominal = true, 
                                double a_lowMtM0FitRange = 0.0, 
                                double a_highMtM0FitRange = 10.0, 
                                int a_numParam = -1, 
                                double* a_params = NULL, 
                                double* a_param_lower_limits = NULL, 
                                double* a_param_upper_limits = NULL, 
                                bool* a_fixed_list = NULL, 
                                double* a_fixed_params = NULL,
                                double a_deltaLowMtM0FitRange = -1,
                                double a_deltaHighMtM0FitRange = -1
                                );
  void fitSingleRapiditySpectra_RapidityCentralityLooper(int a_partIndex, int a_charge, string a_functionName,int a_centIndex = -1, bool a_isNominal = true, double a_lowMtM0FitRange = 0.0, double a_highMtM0FitRange = 10.0, double a_deltaLowMtM0FitRange = -1, double a_deltaHighMtM0FitRange = -1); //CentIndex = -1 for all centralities
  

  void fitPionSpectra(bool a_constrainTemperature = true, int a_numTrialsForDnDyBootstrap = 2500, int a_centIndex = 0, int a_charge = 1); // bootstrap doesn't matter unless you specified
  void fitPionSpectraSimultaneous(bool saveFinalFits = false, bool seedAllWith_yCM = true, int a_numTrialsForDnDyBootstrap = 2500);

  //Fit all pion spectra simultaneously. dN/dy free and  T/Cosh(y-y_cm) temperature dependence
  void fitPionSpectraSimultaneous_BoseEinstein(bool a_isNominal = true, bool a_enforceDnDySymmetry = true, bool a_saveFinalFits = false);
  void fitPionSpectraSimultaneous_doubleBoseEinstein(bool saveFinalFits = false, bool seedAllWith_yCM = true, int a_numTrialsForDnDyBootstrap = 2500, bool fixPlustoMinus = false);
  void fitPionSpectraSimultaneous_doubleBoseEinstein_CoulombCorrected(bool saveFinalFits = false, bool seedAllWith_yCM = true);
  void fitPionSpectraSimultaneous_doubleBoseEinstein_CoulombCorrected_FreeTemps(bool saveFinalFits = false, bool seedAllWith_yCM = true);
  void fitPionSpectraWithDoubleBoseEinstein(double a_momOfHighTempConstraint);

  void fitKaonSpectra(bool a_constrainTemperature = true, int a_numTrialsForDnDyBootstrap = 2500);
  void fitKaonSpectraSimultaneous(bool saveFinalFits = false, bool seedAllWith_yCM = true, int a_numTrialsForDnDyBootstrap = 2500);
  void fitKaonSpectraSimultaneousSymmetric(bool saveFinalFits = false, bool seedAllWith_yCM = true, int a_numTrialsForDnDyBootstrap = 2500);

  //void fitProtonSpectra(int a_lowRapIndex, int a_highRapIndex, bool a_constrainTemperature = true);

  void fitProtonSpectra(int a_centIndex, int a_charge, bool a_constrainBetaS = true, int a_numTrialsForDnDyBootstrap = 2500, double a_lowRapBound_BetaS_Constraint = -5.0, double a_highRapBound_BetaS_Constraint = 5.0);
  void fitSpectra_Simultaneous_SingleCharge_Legendre_Constraints(int a_partIndex, int a_centIndex, int a_charge, bool a_isNominal, int a_lowRapIndex = -1, int a_highRapIndex = -1, double a_lowMtM0FitRange = 0.0, double a_highMtM0FitRange = 10.0); // This function requires that you run the fitSingleRapiditySpectra_RapidityCentralityLooper for protons first, hardcoded blastwave and thermal funct
  void fitProtonSpectraFXT(bool a_constrainBetaS = true, int a_numTrialsForDnDyBootstrap = 2500, double a_lowRapBound_BetaS_Constraint = -5.0, double a_highRapBound_BetaS_Constraint = 5.0);
  void fitProtonSpectraFXT_betaConstr(bool a_constrainBetaS = true, int a_numTrialsForDnDyBootstrap = 2500, double a_lowRapBound_BetaS_Constraint = -5.0, double a_highRapBound_BetaS_Constraint = 5.0);
  void fitProtonSpectraFXT_betaConstr_newIntegration(bool a_constrainBetaS = true, double a_lowRapBound_BetaS_Constraint = -5.0, double a_highRapBound_BetaS_Constraint = 5.0);
  void fitProtonSpectraFXT_betaConstr_newIntegration_forceSymmetry(bool a_constrainBetaS = true, double a_lowRapBound_BetaS_Constraint = -5.0, double a_highRapBound_BetaS_Constraint = 5.0);

  // Simultanteous Proton Fitter
  //void fitProtonSpectra(int a_numTrialsForDnDyBootstrap = 2500);
  // By Rapidity Bin, plus minus pi K p are fit together with a blastwave
  // MUST HAVE SAME RAPIDITY BINNING FOR NOW
  void fitParticlesSimultaneouslyByRapidity();

  double getBlastWaveNFix(){return m_blastWave_nFix;};
  int    getNRapidityBins(int a_partIndex, int a_detectorIndex);
  double getRapidity(int a_partIndex, int a_detectorIndex, int a_rapIndex);
  
  void fitProtonSpectraIndivBlastWave(bool a_constrainTemperature = true);
  void fitProtonSpectraThermal(bool a_constrainTemperature = true);

  void fitProtonDnDy();
  void fitDnDy(bool a_fitMirror = true);
  void compute_dNdEta();

  void writeOutputs(string a_filename, string a_option = "RECREATE");

  //####  CHI SQUARED FUNCITONS ####
  //Used by Root Minimizer, No need to call in a macro, must be public to work
  double simultaneous_blastWave_ChiSqr(const double* a_param); // all rapidities of protons
  double simultaneous_blastWave_betaConstr_ChiSqr(const double* a_param); // all rapidities of protons
  double simultaneous_blastWave_betaConstr_forceSymmetry_ChiSqr(const double* a_param); // all rapidities of protons
  double simultaneous_blastWave_Legendre_BetaS_And_Temp_Plus_ChiSqr(const double* a_param); // all rapidities of protons     NOTE: m_blastWave must have n & m set before calls made
  double simultaneous_blastWave_Legendre_BetaS_And_Temp_Minus_ChiSqr(const double* a_param); // all rapidities of protons     NOTE: m_blastWave must have n & m set before calls made
  double simultaneous_tsallis_Legendre_qm1_And_Temp_Plus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_tsallis_Legendre_qm1_And_Temp_Minus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_boseEinstein_Legendre_Temp_Plus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_boseEinstein_Legendre_Temp_Minus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_mTm0Exponential_Legendre_Temp_Plus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_mTm0Exponential_Legendre_Temp_Minus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_thermal_Legendre_Temp_Plus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_thermal_Legendre_Temp_Minus_ChiSqr(const double* a_param); // all rapidities of protons     (params set before calls)
  double simultaneous_thermal_ChiSqr(const double* a_param); // all rapidities of protons
  double simultaneous_thermal_FXT_ChiSqr(const double* a_param); // all rapidities of protons
  double simultaneousParticle_blastWave_ChiSqr(const double* a_param); // one rapidity all particles
  double simultaneous_doubleThermal_ChiSqr(const double* a_param); // all rapidities of pions
  double simultaneous_doubleBoseEinstein_ChiSqr(const double* a_param); // all rapidities of pions
  double simultaneous_doubleBoseEinstein_CoulombCorrected_ChiSqr(const double* a_param); // all rapidities of pions, both plus and minus
  double simultaneous_doubleBoseEinstein_CoulombCorrected_FreeTemps_ChiSqr(const double* a_param); // all rapidities of pions, both plus and minus
  double simultaneous_doubleBoseEinstein_CoulombCorrected_FreeTemps_noMem_ChiSqr(const double* a_param); // all rapidities of pions, both plus and minus
  double simultaneous_mTm0Exponential_BE_ChiSqr(const double* a_param); // all rapidities of pions
  double simultaneous_BoseEinstein_ChiSqr(const double* a_param); // all rapidities of pions
  double simultaneous_mTm0Exponential_ChiSqr(const double* a_param); // all rapidities of kaons
  double simultaneous_mTm0Exponential_symm_ChiSqr(const double* a_param); // all rapidities of kaons
  double simultaneous_thermal_K_ChiSqr(const double* a_param); // all rapidities of kaons
  double simultaneous_pTExponential_ChiSqr(const double* a_param); // all rapidities of kaons


  double boseEinstein_VC_fit_Plus_noMem(double *x, double *par);
  double boseEinstein_VC_fit_Minus_noMem(double *x, double *par);
  double doubleBoseEinstein_VC_fit_Ch(double *x, double *par);
  double doubleBoseEinstein_VC_fit_Ch_HADES(double *x, double *par);


  double integrate_mTm0_invariant_with_2pi_mT(double *x, double *par); //For dN/dy integral error


  //TESTERS
  void test_doubleBE_CC();
  void test_blastwave_function();


private:

  //TF1* m_fitFunctions[9][3]; // particle // nominal;systematic;...

  double m_nPart[7]    = {334.432,287.019,224.183,159.5835,101.4045,59.1425,20.184};
  double m_nPartErr[7] = {6.533,10.048,11.374192,8.6947497,6.3299037,5.707869,2.0499263};

  double m_Tpro_amp[7] = {0.252297,0.23596,0.20449,0.188534,0.163917,0.13915,0.127722};
  double m_Tpro_sig[7] = {1.57028,1.44414,1.37874,1.43276,1.37272,1.35113,1.46965};

  bool m_useDoubleBoseForPions;

  string        m_imageDir;
  bool          m_particleLoaded[9];
  ParticleInfo* m_partInfo;
  CutClass*     m_cuts;
  TH1D*         m_rapidity_structure[9][3]; // tpc btof etof

  int    m_currentPartIndex;
  int    m_currentCentIndex;
  int    m_currentRapIndex;
  int    m_currentMidRapIndex;
  double m_currentRapValue;
  double m_currentParticleMass;
  TF1*   m_currentFitForIntegration;
  double m_currentLowFitRange;
  double m_currentHighFitRange;
  int    m_currentCharge;
  int    m_pionFittingStage;
  int    m_kaonFittingStage;
  int    m_protonFittingStage;

  bool   m_isCurrentFitNominal;
  bool   m_oldSpectraLoadingConvention;
  double m_temperatureSeed;
  double m_lowTemperatureSeed;
  double m_blastWave_temperature_seed;
  double m_blastWave_betaS_seed;
  int    m_legendreTemperatureOrder[9][2];
  int    m_legendreParTwoOrder[9][2];// could be betaS or q-1 for tsallis
  TGraphErrors* m_currentSpectra;
  TGraphErrors* m_currentSpectra_NotInvariant;
  TGraphAsymmErrors* m_currentSpectraAsymm;
  TH1D* m_currentTemperatureHisto;
  TH1D* m_currentDnDyHisto;
  bool  m_isTemperatureConstrained;
  bool  m_isLowTemperatureConstrained;
  bool  m_constrainAntiProtonBetaS_withParamsFromProtons;
  bool  m_enforceDnDySymmetry;
  bool  m_useBootstrapperForIntegralError;

  bool m_isTPC;
  bool m_isBTOF;
  bool m_isETOF;

  // STORAGE FOR TPC BTOF AND ETOF SPECTRA
  TGraphErrors* m_spectra_TPC_Plus[9][16][g_SpecFit_MaxNumRapBins]; // part cent yindex
  TGraphErrors* m_spectra_TPC_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TGraphErrors* m_spectra_BTOF_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TGraphErrors* m_spectra_BTOF_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TGraphErrors* m_spectra_ETOF_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TGraphErrors* m_spectra_ETOF_Minus[9][16][g_SpecFit_MaxNumRapBins];

  //TGraphErrors* m_spectra_inclusive_TPC_Plus[9][16][g_SpecFit_MaxNumRapBins]; // part cent yindex
  //TGraphErrors* m_spectra_inclusive_TPC_Minus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors* m_spectra_inclusive_BTOF_Plus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors* m_spectra_inclusive_BTOF_Minus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors* m_spectra_inclusive_ETOF_Plus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors* m_spectra_inclusive_ETOF_Minus[9][16][g_SpecFit_MaxNumRapBins];

  // dNdy Storage Histograms
  TH1D* m_dNdy_Plus_Nominal[9][16]; // part cent
  TH1D* m_dNdyHigh_Plus_Nominal[9][16]; // part cent
  TH1D* m_dNdyLow_Plus_Nominal[9][16]; // part cent
  TH1D* m_dNdy_Plus_Systematic[9][16]; // part cent
  TH1D* m_dNdyHigh_Plus_Systematic[9][16]; // part cent
  TH1D* m_dNdyLow_Plus_Systematic[9][16]; // part cent

  TH1D* m_dNdy_Minus_Nominal[9][16];
  TH1D* m_dNdyHigh_Minus_Nominal[9][16];
  TH1D* m_dNdyLow_Minus_Nominal[9][16];
  TH1D* m_dNdy_Minus_Systematic[9][16];
  TH1D* m_dNdyHigh_Minus_Systematic[9][16];
  TH1D* m_dNdyLow_Minus_Systematic[9][16];

  TH1D* m_dNdy_Fiducial_Measured[9][16][2]; // part cent +/-
  TH1D* m_dNdy_Fiducial_And_Extrap_Nominal[9][16][2]; // part cent +/-  // This is fiducial measured + extrapFromFit
  TH1D* m_dNdy_Fiducial_And_Extrap_Systematic[9][16][2]; // part cent +/-  // This is fiducial measured + extrapFromFit

  //Fit Integrals for dN/dy
  TH1D* m_dNdy_LowMtM0_Extrap_Nominal[9][16][2]; // part cent +/-
  TH1D* m_dNdy_LowMtM0_Extrap_Systematic[9][16][2]; // part cent +/-
  TH1D* m_dNdy_MidMtM0_Extrap_Nominal[9][16][2]; // part cent +/-
  TH1D* m_dNdy_MidMtM0_Extrap_Systematic[9][16][2]; // part cent +/-
  TH1D* m_dNdy_HighMtM0_Extrap_Nominal[9][16][2]; // part cent +/-
  TH1D* m_dNdy_HighMtM0_Extrap_Systematic[9][16][2]; // part cent +/-
  TH1D* m_dNdy_Fiducial_Nominal[9][16][2]; // part cent +/-
  TH1D* m_dNdy_Fiducial_Systematic[9][16][2]; // part cent +/-



  //dNdy Storage Histograms //this is for the bootstrap method, it iterates many times and computes an average dN/dy for all the integrals it did
  TH1D* m_dNdy_Plus_Nominal_Avg[9][16]; // part cent
  TH1D* m_dNdy_Plus_Systematic_Avg[9][16];
  TH1D* m_dNdy_Minus_Nominal_Avg[9][16];
  TH1D* m_dNdy_Minus_Systematic_Avg[9][16];

  // dNdy Storage Histograms
  TGraphErrors* m_dNdySymmetric_Plus_Nominal[9][16]; // part cent
  TGraphErrors* m_dNdySymmetric_Plus_Systematic[9][16];
  TGraphErrors* m_dNdySymmetric_Minus_Nominal[9][16];
  TGraphErrors* m_dNdySymmetric_Minus_Systematic[9][16];
  TGraphErrors* m_dNdyGhostSymmetric_Plus_Nominal[9][16]; // part cent
  TGraphErrors* m_dNdyGhostSymmetric_Plus_Systematic[9][16];
  TGraphErrors* m_dNdyGhostSymmetric_Minus_Nominal[9][16];
  TGraphErrors* m_dNdyGhostSymmetric_Minus_Systematic[9][16];

  // Temperature Histograms (protons not included)
  TH1D* m_Temp_Plus_Nominal[9][16]; // part cent
  TH1D* m_TempLow_Plus_Nominal[9][16]; // part cent
  TH1D* m_Temp_Plus_Systematic[9][16];
  TH1D* m_TempLow_Plus_Systematic[9][16]; 
  TH1D* m_Temp_Minus_Nominal[9][16];
  TH1D* m_TempLow_Minus_Nominal[9][16];
  TH1D* m_Temp_Minus_Systematic[9][16];
  TH1D* m_TempLow_Minus_Systematic[9][16];

  //Simultaneous Particle BlastWave Fitting (all 6 +/- particles T and Beta_S by Rapidity) 
  TH1D* m_dNdy_Plus_SimPartBW[9][16]; // part cent
  TH1D* m_dNdy_Minus_SimPartBW[9][16]; // part cent
  TH1D* m_BetaS_SimPartBW[16]; //cent
  TH1D* m_Temp_SimPartBW[16]; //cent

  //Temperature Histogram for Protons
  TH1D* m_Proton_Temp_Plus_Nom;
  TH1D* m_Proton_Temp_Plus_Sys;
  TH1D* m_Proton_Temp_Minus_Nom;
  TH1D* m_Proton_Temp_Minus_Sys;

  TH1D* m_BetaS_Plus_Nominal[16]; // cent
  TH1D* m_BetaS_Minus_Nominal[16]; // cent

  //4pi Yields by Centrality
  TGraphErrors* m_FourPiYields_Nom_Plus[16];
  TGraphErrors* m_FourPiYields_Nom_Minus[16];
  TGraphErrors* m_FourPiYields_Sys_Plus[16];
  TGraphErrors* m_FourPiYields_Sys_Minus[16];

  //Storage of fit paramters for isolated rapidity fits

  TH1D* m_Isolated_Nominal_Fit_Params_Plus[9][16][8];
  TH1D* m_Isolated_Nominal_Fit_Params_Minus[9][16][8];
  TH1D* m_Isolated_Systematic_Fit_Params_Plus[9][16][8];
  TH1D* m_Isolated_Systematic_Fit_Params_Minus[9][16][8];

  //Storage of fit paramters for simultaneous (final) rapidity fits
  TH1D* m_Simultaneous_Nominal_Fit_Params_Plus[9][16][8];
  TH1D* m_Simultaneous_Nominal_Fit_Params_Minus[9][16][8];
  TH1D* m_Simultaneous_Systematic_Fit_Params_Plus[9][16][8];
  TH1D* m_Simultaneous_Systematic_Fit_Params_Minus[9][16][8];

  TGraphErrors*  m_totalSpectra_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TGraphErrors*  m_totalSpectra_NotInvariant_Plus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors*  m_totalSpectraInclusive_Plus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors*  m_totalSpectraInclusive_NotInvariant_Plus[9][16][g_SpecFit_MaxNumRapBins];
  double m_totalSpectra_Plus_RapVals[9][16][g_SpecFit_MaxNumRapBins]; // because etof might have spectra without dEdx
  TF1*   m_nominalFits_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_nominalFits_NotInvariant_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_NotInvariant_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_SimPartBWFits_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_nominalFits_Isolated_Plus[9][16][g_SpecFit_MaxNumRapBins];               //These isolated fits are from a single rapidity fit only (not simultaneous rapidity)
  TF1*   m_nominalFits_Isolated_NotInvariant_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_nominalFitsCovMat_Isolated_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_nominalFitsCorrMat_Isolated_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_Isolated_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_Isolated_NotInvariant_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_systematicFitsCovMat_Isolated_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_systematicFitsCorrMat_Isolated_Plus[9][16][g_SpecFit_MaxNumRapBins];

  TGraphErrors* m_totalSpectra_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TGraphErrors* m_totalSpectra_NotInvariant_Minus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors* m_totalSpectraInclusive_Minus[9][16][g_SpecFit_MaxNumRapBins];
  //TGraphErrors* m_totalSpectraInclusive_NotInvariant_Minus[9][16][g_SpecFit_MaxNumRapBins];
  double m_totalSpectra_Minus_RapVals[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_nominalFits_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_nominalFits_NotInvariant_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_NotInvariant_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_SimPartBWFits_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_nominalFits_Isolated_Minus[9][16][g_SpecFit_MaxNumRapBins];           //These isolated fits are from a single rapidity fit only (not simultaneous rapidity)
  TF1*   m_nominalFits_Isolated_NotInvariant_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_nominalFitsCovMat_Isolated_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_nominalFitsCorrMat_Isolated_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_Isolated_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*   m_systematicFits_Isolated_NotInvariant_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_systematicFitsCovMat_Isolated_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TH2D*  m_systematicFitsCorrMat_Isolated_Minus[9][16][g_SpecFit_MaxNumRapBins];

  // Physics should be symmetric, so you can combine the corrected spectra and fit them as symmetric
  TGraphErrors*  m_totalSpectra_Symmetric_Plus[9][16][g_SpecFit_MaxNumRapBins];
  double         m_totalSpectra_Symmetric_Plus_RapVals[9][16][g_SpecFit_MaxNumRapBins]; // because etof might have spectra without dEdx
  TF1*           m_nominalFits_Symmetric_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*           m_systematicFits_Symmetric_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*           m_SimPartBWFits_Symmetric_Plus[9][16][g_SpecFit_MaxNumRapBins];

  TGraphErrors*  m_totalSpectra_Symmetric_Minus[9][16][g_SpecFit_MaxNumRapBins];
  double         m_totalSpectra_Symmetric_Minus_RapVals[9][16][g_SpecFit_MaxNumRapBins];
  TF1*           m_nominalFits_Symmetric_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*           m_systematicFits_Symmetric_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TF1*           m_SimPartBWFits_Symmetric_Minus[9][16][g_SpecFit_MaxNumRapBins];


  //Bootstrapping Data (history of integrals because of sampling parameters for dN/dy);
  TH1D* m_bootstrapData_Nom_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TH1D* m_bootstrapData_Nom_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TH1D* m_bootstrapData_Sys_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TH1D* m_bootstrapData_Sys_Minus[9][16][g_SpecFit_MaxNumRapBins];
  TH1D* m_bootstrapData_SimPartBW_Plus[9][16][g_SpecFit_MaxNumRapBins];
  TH1D* m_bootstrapData_SimPartBW_Minus[9][16][g_SpecFit_MaxNumRapBins];

  //Fitting Functions
  //Pions
  TF1* m_boseEinstein; //Bose Einstein:  [0]=dN/dy [1]=kT  [2]=mass (fixed) 
  TF1* m_boseEinstein_CoulombCorrected; //Bose Einstein:  [0]=Amp [1]=kT  [2]=mass (fixed) [3]=VC
  TF1* m_pTExponential; //pT Exponential:  [0]=dN/dy [1]=kT 
  TF1* m_tsallis; // TSallis Function
  TF1* m_doubleBoseEinstein;
  TF1* m_doubleBoseEinstein_CoulombCorrected;
  TF1* m_doubleThermal; //double Thermal: [0]=dN/dy_low [1]=kT_low [2]=mass (fixed) [3]=dN/dy_high [4]=kT_high
  TF1* m_thermal_boseEinstein; //thermal Bose Einstein: [0]=dN/dy_low [1]=kT_low [2]=mass (fixed) [3]=A_high [4]=kT_high
  TF1* m_mTm0Exponential_BE; //mTm0 Exponential:  [0]=dN/dy_low [1]=kT_low [2]=mass (fixed) [3]=dN/dy_high [4]=kT_high

  vector<double> m_VC; //Coulomb potential for fixing with Coulomb corrected Bose-Einstein:  [cent]

  //Kaons
  TF1* m_mTm0Exponential; //mTm0 Exponential:  [0]=dN/dy [1]=kT [2]=mass (fixed) 
  TF1* m_boltzmann; //Boltzmann: [0]=dN/dy [1]=kT [2]=mass (fixed) 

  //Protons
  TF1* m_mTExponential; //mT Exponential: [0]=dN/dy [1]=kT [2]=mass (fixed)

  TF1* m_blastWave;/// BEWARE, DO NOT CLONE
  TF1* m_blastWaveIntegrand;
  double m_blastWave_nFix;
  TF1* m_thermal; //Thermal: [0]=dN/dy [1]=kT [2]=mass (fixed)


  //####  NOT INVARIANT FUNCTIONS (Scaled by 2PiMt)
  TF1* m_boseEinstein_NotInvariant;    //Bose Einstein:  [0]=dN/dy [1]=kT  [2]=mass (fixed) 
  TF1* m_boseEinstein_NotInvariant_WithLegendreTemp;
  TF1* m_pTExponential_NotInvariant;   //pT Exponential:  [0]=dN/dy [1]=kT
  //TF1* m_pTExponential_NotInvariant_WithLegendreTempAndBetaS; 
  TF1* m_tsallis_NotInvariant;
  TF1* m_tsallis_NotInvariant_WithLegendreTempAndQm1; 
  TF1* m_mTm0Exponential_NotInvariant; //mTm0 Exponential:  [0]=dN/dy [1]=kT [2]=mass (fixed) 
  TF1* m_mTm0Exponential_NotInvariant_WithLegendreTemp; 
  TF1* m_boltzmann_NotInvariant;       //Boltzmann: [0]=dN/dy [1]=kT [2]=mass (fixed) 
  TF1* m_mTExponential_NotInvariant;   //mT Exponential: [0]=dN/dy [1]=kT [2]=mass (fixed)
  TF1* m_blastWave_NotInvariant;       // BEWARE, DO NOT CLONE
  TF1* m_blastWave_NotInvariant_WithLegendreTempAndBetaS; //beware don't clone
  TF1* m_thermal_NotInvariant;         //Thermal: [0]=dN/dy [1]=kT [2]=mass (fixed)
  TF1* m_thermal_NotInvariant_WithLegendreTemp; //beware don't clone



  TF1* m_simultaneousBlastWave;
  
  ROOT::Math::Minimizer* m_minimizer;
  ROOT::Math::Functor* m_chiSqrdFunctor;

  ROOT::Math::Minimizer* m_minimizer_sys;
  ROOT::Math::Functor* m_chiSqrdFunctor_sys;

  //Temperature Constraint Functions
  TF1* m_constFunct;
  TF1* m_gausFunct;

  TF1* m_mirroredGausFunct;
  TF1* m_mirroredGausPlusGausFunct;


  
  //scales by 2*pi*mT assuming a_graph is in mT-m0 space
  void scaleTGraphByInvFactor(TGraphErrors* a_graph, double a_mass);
  double blastWaveModelFit(double *x, double *par);
  double blastWaveModelFit_NotInvariant(double *x, double *par);
  double blastWaveModelFit_NotInvariant_WithLegendreTempAndBetaS(double *x, double *par);
  double boseEinstein_VC_fit(double *x, double *par);
  double boseEinstein_VC_fit_Plus(double *x, double *par);
  double boseEinstein_VC_fit_Minus(double *x, double *par);
  //double boseEinstein_VC_fit_Plus_noMem(double *x, double *par);
  //double boseEinstein_VC_fit_Minus_noMem(double *x, double *par);
  double doubleBoseEinstein_VC_fit(double *x, double *par);
  //double doubleBoseEinstein_VC_fit_Ch(double *x, double *par);



  void init_boseEinstein(TF1* a_alternateFunct = NULL);
  void init_boseEinstein_CoulombCorrected();
  void init_doubleBoseEinstein(double a_dNdyHigh = -999, double a_tempHigh = -999);
  void init_doubleBoseEinstein_CoulombCorrected(double a_ampHigh = -999, double a_tempHigh = -999, double a_VC = -999);
  void init_doubleThermal();
  void init_doubleThermal_v2(double a_dNdyHigh = -999, double a_tempHigh = -999);
  void init_Thermal_boseEinstein();
  void init_pTExponential(TF1* a_alternateFunct = NULL);
  void init_tsallis(TF1* a_alternateFunct = NULL, int a_regionIndex = 0); // region 0=mid 1=low 2=high
  void init_mTExponential(TF1* a_alternateFunct = NULL);
  void init_mTm0Exponential(TF1* a_alternateFunct = NULL);
  void init_mTm0Exponential_BE();
  void init_mTm0Exponential_BE_v2(double a_dNdyHigh = -999, double a_tempHigh = -999);
  void init_boltzmann(TF1* a_alternateFunct = NULL);
  void init_thermal(TF1* a_alternateFunct = NULL);
  void init_blastWave(TF1* a_alternateFunct = NULL);


  void init_simul_blastWave();
  void init_simul_particle_blastWave();
  
  int m_currentDnDyParameters;
  vector< int > m_liveRapiditySpectra;
  vector< int > m_liveRapSpec_to_ParamIndex; //if enforcing dNdySymmetry, some will use same parameters
  vector< int > m_ParamIndex_to_LiveRapSpec;
  vector< double > m_liveRapidityValueSpectra;

  vector< int > m_liveRapidityPionSpectra;
  vector< double > m_liveRapidityValuePionSpectra;
  vector< int > m_liveRapidityKaonSpectra;
  vector< double > m_liveRapidityValueKaonSpectra;
  vector< int > m_liveRapidityProtonSpectra;
  vector< double > m_liveRapidityValueProtonSpectra;

  vector< vector< int > > m_liveRapidityPionSpectra_charge;
  vector< vector< double > > m_liveRapidityValuePionSpectra_charge;

  double m_Vc[20][7];

  void init_constFunctForTemp();
  void init_gausFunctForTemp();
  void init_doubleGausFunctForTemp();
  void init_doubleGausPlusGausFunctForDnDy();

  void makeSymmetricDnDy(TH1D* a_histo, TGraphErrors* a_mainGraph, TGraphErrors* a_ghostGraph);






};
#endif