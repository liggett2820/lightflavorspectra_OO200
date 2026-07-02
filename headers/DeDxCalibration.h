#ifndef _DEDX_CALIBRATION_ 
#define _DEDX_CALIBRATION_

//#include "TF1.h"
//#include "TGraphErrors.h"
//#include "TGraphAsymmErrors.h"
//#include <TFitResult.h>
//#include <TFitResultPtr.h>
//#include "Math/Functor.h"
//#include "Math/Minimizer.h"
//#include "Math/Factory.h"
//#include "TVirtualFitter.h"

#include "TFile.h"
#include "TF1.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TFitResult.h"
#include "../headers/HistogramUtilities.h"
#include "../headers/PhysMath.h"
#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"


static const int g_DeDxCalib_MaxNumVzBins  = 30;
static const int g_DeDxCalib_MaxNumEtaBins = 30;
static const int g_DeDxCalib_MaxNumBetaGammaBins = 200;
using namespace std;


//############################
// This class loads the output of PicoBinner (specifically some calibration histograms),
// fits the dE/dx trends with a functional form, and makes a calibration file that can be
// supplied to PicoBinner in a second pass for better dE/dx predictions and widths 
//##############################


class DeDxCalibration{

public:
  DeDxCalibration();
  
  //######   SETUP    #########
  void loadPicoBinnerOutput(string a_filename, int a_partIndex, int a_betaGamma_rebin = 1);
  void loadSaveFile(string a_filename); // loads previously made mean/width graphs (not necessary, for development)
  void setInitialCalibrationCurve(TF1* a_funct); // override function. loadPicoBinnerOutput automatically loads a calibration function.
  void setVzStructureHistogram(TH1I* a_histo);   // override function. loadPicoBinnerOutput automatically loads a structure histogram.
  void setCdfRebinNumBins(int a_numBins);        //This is for ln(dE/dx) space, depending on statistics
  void setImageDir(string a_directory);          //setting this will toggle the writing of images
  void setTestingToggle(bool a_toggle = true);
  

  //######    CALIBRATION OPERATIONS     #########
  void extractMeansAndWidths(int a_partIndex, bool a_saveImages = false);
  void calibrateMeans(int a_partIndex, bool a_doOneSigma = false); // if true, does a mean+sigma fit
  void calibrateWidths(int a_partIndex);
  

  //######    WRITE OUTPUT     #########
  void makeCalibrationFile(string a_outputFileName);

protected:

  TH1D*              m_currentHistoToFit;
  TGraphAsymmErrors* m_currentGraphToFit;
  string             m_imageDir;
  ParticleInfo*      m_partInfo;
  int                m_numBins_cdf_rebin;
  bool               m_testing;
  double             m_maxBetaGammaForTPCFitting[9];

  //###########      LOADED INFORATION      #################
  TF1*               m_initialCalibCurve; // should be loaded from CutClass // this is an 1/(betagamma)^2 funciton
  TF1*               m_initialCalibCurve_betaGamma; // this is an betagamma funciton
  TH1I*              m_vZ_calib_struct;
  TH1I*              m_eta_calib_struct;
  TH1I*              m_betaGamma_calib_struct;
  TH1I*              m_betaGamma_calib_struct_fromTOF;
  TH1D*              m_lnDeDx_calibrationHistograms[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins][g_DeDxCalib_MaxNumBetaGammaBins];
  TH1D*              m_lnDeDx_calibrationHistograms_fromTOF[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins][g_DeDxCalib_MaxNumBetaGammaBins];


  //###########    CALIBRATION INFROMATION     ###############
  TGraphAsymmErrors* m_mean_calibrationGraphs[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins];
  TGraphAsymmErrors* m_meanPlusSigma_calibrationGraphs[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins];
  TGraphAsymmErrors* m_sigma_calibrationGraphs[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins];
  TF1*               m_mean_calibrationFunctions[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins];
  TF1*               m_meanPlusSigma_calibrationFunctions[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins];
  TF1*               m_sigma_calibrationFunctions[9][g_DeDxCalib_MaxNumVzBins][g_DeDxCalib_MaxNumEtaBins];


};




#endif
