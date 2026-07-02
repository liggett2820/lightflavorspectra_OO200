#ifndef _Gaus_Mixture_Model_
#define _Gaus_Mixture_Model_

// Ported verbatim from lightflavorspectra_etof/headers/GausMixModel.h -- a standalone
// numerical utility class (Gaussian mixture model fitting for dE/dx-type peak
// decomposition) with no dependence on particle species/energy/eventConfig, so nothing
// needed to change for the 200 GeV O+O pi/K/p restriction.

#include <stdlib.h>
//#include <time.h>
#include <vector>
#include <string>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TPad.h"
#include "TRandom3.h"

#include "../headers/HistogramUtilities.h"
//#include "/scratch_altair/mharasty/repos/picodstanalysis/src/submodules/ParticleInfo/ParticleInfo/ParticleInfo.h"

using namespace std;

// Design Notes:
// KMEANS Not Yet Suitible for Variable bin width histograms
// doGMM must be a differential histogram (e.g. histo->Scale(1.0,"width"))
class GausMixModel{
public:
  GausMixModel();
  //GausMixModel();

  void importDataAsHistogram(TH1D* a_histo);
  void importDataAs2DHistogram(TH2D* a_histo);
  //void importDataAsArray(int a_numData, double* a_data);
  void setNumberOfGaussians(int a_numGaus);
  //void initializeParameters(int a_numParams, double* a_parameters);

  //This will perform a K-Means Algorithm to initialize the means of the distribution
  //This will return false if the KMeans Fails
  bool doKMeans(double* a_meansToFill);

  //K-Means for 2D...
  bool doKMeans2D(double* a_meansToFill);

  //This will call KMeans to initialize the means, then carry out a gausian mixture
  //model algorithm to find all the parameters of the gaussian
  //the result will be stored in a_parametersToFill
  //Note, this returns parameters as    Fraction Mean Variance, Fraction Mean Variance, ...
  //
  // SDCC/ROOT5-CINT BUILD FIX (2026-07-02, not in the original): a_initialMeans/
  // a_initialAmps/a_initialSigmas' default values (all NULL) were removed here. This
  // overload and the one below it (doGMM with a_avgSigma instead of a_initialSigmas)
  // originally both had every parameter after the first defaulted, and both share the
  // exact same first-3-parameter types (double*,double*,double*) -- when CINT's
  // rootcint builds this class's dictionary, it generates interactive-call stubs for
  // every "N arguments given, rest defaulted" arity, and for N=1,2,3 the two overloads'
  // stubs become genuinely ambiguous (this only surfaces under ROOT5/CINT, not ROOT6/
  // Cling, which is why it never showed up building on a Mac). This overload's only
  // call site (inside GausMixModel.cxx, from the other overload's body) always passes
  // all 7 arguments explicitly, so requiring the first 4 removes the ambiguous stub
  // arities without changing any real call site's behavior.
  void doGMM(double* a_parametersToFill,
             double* a_initialMeans,
             double* a_initialAmps, //each are [0,1] fractions
             double* a_initialSigmas,
             int a_maxTrials = 300,
             double a_tolerance = 1e-7,
             double a_minBinContent = -1 // this will zero bins below the threshold if positive
            );

  //This is the same above, but avgSigma as seed instead
  void doGMM(double* a_parametersToFill,
             double* a_initialMeans = NULL,
             double* a_initialAmps = NULL, //each are [0,1] fractions
             int a_maxTrials = 300,
             double a_tolerance = 1e-7,
             double a_avgSigma=-1,
             double a_minBinContent = -1 // this will zero bins below the threshold if positive
            );


  // ampli
  void doGMM2D(double* a_parametersToFill, // [amp,muX,muY,covXX,covXY,covYY]  //size 6*numGaus
               double* a_initialMeans = NULL, //[muX1,muY1,muX2,muY2,...] size 2*numGaus
               double* a_initialAmps = NULL,  // [amp1, amp2, ...] size numGaus
               double a_avgXSigma = -1,
               double a_avgYSigma = -1,
               int a_maxTrials = 500,
               double a_tolerance = 1.0E-7, // tolerance for change in log likelihood
               bool  a_forceZeroCovariance = false
               );
               //double* a_amplitudesToFill, // [amp1, amp2, ...] size numGaus
               //double* a_meansToFill, //[muX1,muY1,muX2,muY2,...] size 2*numGaus
               //double* a_covarMatrixToFill,//[covXX1,covXY1,covYY1,covXX2,covXY2,covYY2,...] size 3*numGaus

  double getChiSquared();
  bool   isTooManyGaus();  // not yet fully tested
  double getYield(int a_gausNumber); // stored final parameters are fractions of total, this converts to overall yield

  //This returns a TF1 of the gaussians (new external to class)
  TF1* getFunction();

  //this gets TH2D with result of fit to draw
  TH2D* getTGraphFor2D();

  TF1* getSingleGaus(int a_gausIndex);

  //This is for making the TF1 to plot for quality assurance
  double operator() (double *x, double *p);

protected:

  int m_dimensions;
  TH1D* m_dataHistogram;
  TH2D* m_dataHistogram2D;
  int m_numGaus;
  vector<double> m_parameters;



};




#endif
