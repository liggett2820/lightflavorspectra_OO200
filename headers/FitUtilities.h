#ifndef _FIT_UTLITIES_
#define _FIT_UTLITIES_

#include "TF1.h"
#include "TDecompChol.h"
#include "TRandom3.h"
#include "PhysMath.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <time.h>

using namespace std;

//defined in cxx
//TF1* spectraFitToScale_FitUtil;

namespace FitUtilities{


  //double scaleSpectraByInvFactor(double *x, double *par);

  //double scalingOperator(double a_x, double a_mass, TF1* a_function);

  //this assumes it is a spectra, does 2pi*mT scaling then integrates
  //cant use with C++ class functions as TF1s
  void bootstrapIntegral(double &a_integral, double &a_integralError ,TF1* a_fitFunction, 
                         TFitResult* a_fitResult,double a_mass,
                         double a_lowIntegralBound = -999, double a_highIntegralBound = -999);

  //returns list of integrals sampled
  vector<double> bootstrapIntegralSingleParameter(double &a_integral,double &a_avgIntegral, double &a_integralError,
              TF1* a_fitFunction, TFitResult* a_fitResult, double a_mass, int a_parameterToVary, int a_numTrials = 1000, double a_paramError = -1,
              double a_lowIntegralBound = -999, double a_highIntegralBound = -999, bool a_extrapolateToInf = true);

  void bootstrapIntegralMultiParameter(double &a_integral,double &a_avgIntegral, double &a_integralError,
              TF1* a_fitFunction, TFitResult* a_fitResult, double a_mass, int a_numTrials = 1000,
              double a_lowIntegralBound = -999, double a_highIntegralBound = -999);


  //#####  TESTERS ########

  void testBootstrapper();

};



#endif


