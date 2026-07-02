#ifndef _SPLINEFITTER_
#define _SPLINEFITTER_


#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TVirtualFitter.h"
#include <TCanvas.h>
#include <TString.h>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"


#include "../headers/HistogramUtilities.h"

using namespace std;



// This class uses a series of second order polynomials to fit an arbitrary graph. for the use of N parabolas, there are 2N+2 free parameters including the x point locations where
// the continuity and first derivative constraints are applied. You can specify the number of points to fit and/or put in your own x positions. If you use findAdditionalSeedingPoints,
// a third derivative graph is used to find the number of points that are more than 1.5sigma away from 0 to find kinks in the second derivave. This places the points close to where they
// need to be to characterize the changing curvature of your funciton. You can get integrals and error if you want too.

// To see what this code can do, try:     testSplineFitter(NULL, 3);




class SplineFitter{

public:

  SplineFitter();

  void setGraphToFit(TGraphAsymmErrors* a_graph);
  //These two are automatically set by setGraphToFit, so these would be overrides after you set the graph
 // void setRegularzationParameter(double a_par){m_regularizationParForClosePoints = a_par;};
  void setMinimumPointSeparation(double a_par){m_minimumDistanceForRegularization = a_par;};
  void analyzeGraph(int a_numSplinePoints, TF1* &a_fitResult, TGraph* &a_splinePoints, double &a_chiSqr, double &a_numDegreeOfFreedom, bool a_fixXPoints = false, bool a_findAdditionalSeedingPoints = false, double* a_forcedXLocs = NULL, bool a_fixedEndDerivative = false, double a_fixedEndDerivativeValue = 0.0);
  //used to get errors from fits
  void getIntegralAndError(double a_lowBound, double a_highBound, double &a_integral, double &a_integralError, bool a_divideByWidth = false);
  TMatrixDSym* getCovMatrix(){return m_covarainceMatrixDSym;};
  //uses
  void analyzeGraph_WithBestNumPoints(int a_lowNumSplinePoints, int a_highNumSplinePoints, TF1* &a_fitResult, TGraph* &a_splinePoints, TGraph* &a_bestNumPointsGraph, double &a_chiSqr, double &a_numDegreeOfFreedom);

  // Testing Functions
  void testSplineFunction(); // generates points and makes sure the spline funciton has the right behavior
  void testSplineFitter(TGraphAsymmErrors* a_graph = NULL, int a_numPoints = 4, bool a_fixXPoints = false); // fits a sample graph to show that this class works properly. If you don't hand it a pointer, there is a default graph to fit in the class. 
  void testSplineFitterBestNumPoints(TGraphAsymmErrors* a_graph = NULL);

  //If you need to get an error on the integral, you need to use the class structure here
  void reload(TF1* a_previousFit, TMatrixDSym* a_covMatrix);
  TF1* getFitFunciton(){return m_fitResult;};

private:

  TGraphAsymmErrors* m_graphToFit;
  double m_paramSpace[40][4]; // x y y' y"
  int m_numParabolas;
  int m_numPoints;
  bool m_drawing;
  double m_regularizationParForClosePoints;
  double m_minimumDistanceForRegularization;
  vector<double> m_additionalSeedingPoints;
  TF1* m_fitResult; // storage of fit result from analyzeGraph
  vector<double> m_finalParameters; // storage of final parameters from analyzeGraph

  double m_currentParameters[82];
  double m_covarainceMatrix[1700]; // only populated after analyzeGraph
  TMatrixDSym* m_covarainceMatrixDSym;

  //evaluates a second order polynomial at a_x given certian parameters
  double parabEval(double a_x, double a_xLoc, double a_yLoc, double a_ySlope, double a_ySecDer);
  //double parabEvalSimple(double a_x, double &a_A, double &a_B, double &a_C);
  //Operator used for TF1 making
  double splineOperator(double* a_x, double* a_params);
  double splineOperatorWithParams(double* a_x, double* a_params);
  //Function to Minimize
  double splineChiSqr(const double* a_params);

  //Calculates continuity and first derivative constraints to fill unnecessary parameters. It updates m_currentParameters
  void calculateConstraintVariables(const double* a_params);
  void printCurrentParamSpace();

  //uses the third derivative to populate m_additionalSeedingPoints with canidates to seed x values
  void findSeedingPoints();

  //Make a graph that holds the points of continuity of constraints (points of spline stitching)
  TGraph* makeSplinePointGraph();

  bool isSorted(const double* a_params);
  void sortParamSpace();

  //Makes an example graph to fit in testing
  TGraphAsymmErrors* makeTestingGraph();
  //void propagateDerivativeToABC();

  ROOT::Math::Minimizer* m_minimizer;
  ROOT::Math::Functor* m_chiSqrdFunctor;


};

#endif