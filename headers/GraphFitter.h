#ifndef _GRAPH_FITTER_ 
#define _GRAPH_FITTER_

/*#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>


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
*/

#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TVirtualFitter.h"

#include "../headers/HistogramUtilities.h"
#include "../headers/PhysMath.h"

using namespace std;


//############################
//This class is for fitting a TGraphAsymmErrors with a chi squared that includes the integral 
//  of the function over the "bin". There is no simple ROOT way to do this at this time...
//##############################


class GraphFitter{
public:
  GraphFitter();

  void loadGraphAsym(TGraphAsymmErrors* a_graph);
  void loadGraph(TGraphErrors* a_graph); // converts it to an asym graph before fitting
  void loadTF1(TF1* a_funct, bool a_isSimpleFunctString = true); // if the function is connected to a class member function, it isn't simple
  void setParametersToBeFixed(vector<int> a_paramsToFix){m_paramsToFix = a_paramsToFix;};
  void setFitRange(double a_lowRange, double a_highRange);
  void minimize();
  int  getStatus(){return m_status;};
  TF1* getResultFunction(){return m_fitResult;};
  vector<double> getResultFunctionParameters();
  TMatrixDSym* getCovMat(){return m_covMat;};
  


  //for the minimizer only (must be public)
  double chiSqrFunction(const double* a_param);

protected:

  TGraphAsymmErrors* m_graphToFit;
  TF1*   m_fitFunction;
  TF1*   m_fitResult;
  int    m_numParameters;
  double m_lowFitRange;
  double m_highFitRange;
  int    m_firstPointInRange;
  int    m_lastPointInRange;
  int    m_status;
  double chiSqr;
  double chiSqrPerDOF;
  TMatrixDSym*   m_covMat;
  vector<int>    m_paramsToFix;
  vector<double> m_final_parameters;





};




#endif
