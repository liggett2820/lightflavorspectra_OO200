#ifndef _CALIBRATION_CLASS_
#define _CALIBRATION_CLASS_

#include <iostream>
#include <string>
#include "TH1D.h"
#include "TF1.h"
#include "TProfile.h"
#include "PhysMath.h"
#include "HistogramUtilities.h"
#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"



static const int g_CalibClass_MaxNumVzBins  = 30;
static const int g_CalibClass_MaxNumEtaBins = 30;
static const int g_CalibClass_MaxNumBetaGammaBins = 200;
using namespace std;


//############################
// This class loads the output of PicoBinner (specifically some calibration histograms),
// fits the dE/dx trends with a functional form, and makes a calibration file that can be
// supplied to PicoBinner in a second pass for better dE/dx predictions and widths 
//##############################



class CalibrationClass{

public:
  CalibrationClass();
  
  //######   SETUP    #########
  void loadDeDxCalibration(string a_filename, bool a_byParamHistos = true); // from DeDxCalibration  
  //void loadBTOFCalibration(string a_filename); // Not Yet Implemented
  //void loadETOFCalibration(string a_filename); // Not Yet Implemented
  //simple profile and parabola fit for low pt
  void loadDeDxCalibrationProfile(string a_filename);

  void setParticleInfo(ParticleInfo* a_info);
  void setDoSigmaAdjustment(bool a_toggle = true);
  void setUsePionCalibForKaon(bool a_toggle = true);

  


  //######  OPERATIONS  ##########
  void   setParticleOfInterest(int a_partIndex);
  double getCalibratedDeDx(double a_Vz, double a_eta, double a_momentum); //the calibration curves
  double getCalibratedNSigma(double a_Vz, double a_eta, double a_momentum, double a_dEdx);
  double getCalibratedZTPC(double a_Vz, double a_eta, double a_momentum, double a_dEdx);

  //uses loadDeDxCalibrationProfile information (TProfile + Parabola)
  double getDeDxFromProfileParab(double a_momentum);

protected:

  ParticleInfo* m_partInfo;
  bool          m_doSigmaAdjustment;
  bool          m_usePionCalibForKaon;
  int           m_particleOfInterest;
  double        m_masses[9];
  bool          m_particle_loaded[9];
  // LOADED INFORATION
  TH1I*              m_vZ_calib_struct;
  TH1I*              m_eta_calib_struct;
  TF1*               m_mean_calibrationFunctions[9][g_CalibClass_MaxNumVzBins][g_CalibClass_MaxNumEtaBins];
  TF1*               m_meanPlusSigma_calibrationFunctions[9][g_CalibClass_MaxNumVzBins][g_CalibClass_MaxNumEtaBins];
  TF1*               m_sigma_calibrationFunctions[9][g_CalibClass_MaxNumVzBins][g_CalibClass_MaxNumEtaBins];
  int                m_minMaxVzEtaIndex[9][4];

  TProfile*   m_dEdx_invBetaGammaSqr_profile;
  TF1*        m_dEdx_invBetaGammaSqr_parabola;

};










#endif
