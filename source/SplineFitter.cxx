#ifndef _SPLINEFITTER_CXX_
#define _SPLINEFITTER_CXX_

#include "../headers/SplineFitter.h"

SplineFitter::SplineFitter(){
  m_numPoints = 0; 
  m_graphToFit = NULL;
  m_numParabolas = 0;
  m_minimizer = NULL;
  m_chiSqrdFunctor = NULL;
  m_drawing = false;
  m_fitResult = NULL;
  for(int iii = 0; iii < 20; iii++){
    m_paramSpace[iii][0] = 0.0;
    m_paramSpace[iii][1] = 0.0;
    m_paramSpace[iii][2] = 0.0;
    m_paramSpace[iii][3] = 0.0;
  }

  m_regularizationParForClosePoints = 1.0;
  m_minimumDistanceForRegularization = 1.0;

}


//============================================================================================================



void SplineFitter::setGraphToFit(TGraphAsymmErrors* a_graph){
  m_graphToFit = a_graph;
  m_minimumDistanceForRegularization = 1e99;

  int numPoints = m_graphToFit->GetN();
  for(int pointIndex = 1; pointIndex < numPoints; pointIndex++){
    for(int pointIndexTwo = pointIndex + 1; pointIndexTwo < numPoints; pointIndexTwo++){
      double dist = fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndexTwo]);
      if(dist < m_minimumDistanceForRegularization) m_minimumDistanceForRegularization = dist;
    }
  }

  m_minimumDistanceForRegularization = m_minimumDistanceForRegularization * 5.0;


}



//============================================================================================================

double SplineFitter::parabEval(double a_x, double a_xLoc, double a_yLoc, double a_ySlope, double a_ySecDer){
    return 0.5*a_ySecDer*a_x*a_x + (a_ySlope - a_ySecDer*a_xLoc)*a_x + a_yLoc + 0.5*a_ySecDer*a_xLoc*a_xLoc - a_ySlope*a_xLoc;
}

//============================================================================================================

/*
double SplineFitter::parabEvalSimple(double &a_x, double &a_A, double &a_B, double &a_C){
    return 0.5*a_A*a_x*a_x + a_B*a_x + a_C;
}

void SplineFitter::paramABCtoDerivative(double &a_xLoc, double &a_A, double &a_B, double &a_C, double &a_yLoc, double &a_ySlope, double &a_ySecDer){
    a_yLoc = 0.5*a_A*a_xLoc*a_xLoc + a_B*a_xLoc + a_C;
    a_ySlope = a_A*a_xLoc + a_B;
    a_ySecDer = a_A;
}

void SplineFitter::paramDerivativeToABC(double &a_xLoc, double &a_yLoc, double &a_ySlope, double &a_ySecDer, double &a_A, double &a_B, double &a_C){
    a_C = a_yLoc - a_ySlope*a_xLoc + 0.5*a_ySecDer*a_xLoc*a_xLoc; // y - y'x + 0.5y"x^2
    a_B = a_ySlope - a_ySecDer*a_xLoc; // y'-y"x
    a_A = a_ySecDer;
}
//============================================================================================================


void SplineFitter::propagateDerivativeToABC(){
    for(int iii = 0; iii < m_numParabolas; iii++){
        m_paramSpaceSimple[iii][0] = m_paramSpace[iii][0];
        m_paramSpaceSimple[iii][3] = m_paramSpace[iii][1] - m_paramSpace[iii][2]*m_paramSpace[iii][0] + 0.5*m_paramSpace[iii][3]*m_paramSpace[iii][0]*m_paramSpace[iii][0]; // y - y'x + 0.5y"x^2
        m_paramSpaceSimple[iii][2] = m_paramSpace[iii][2] - m_paramSpace[iii][3]*m_paramSpace[iii][0]; // y'-y"x
        m_paramSpaceSimple[iii][1] = m_paramSpace[iii][3];
    }
}
*/
////============================================================================================================


void SplineFitter::calculateConstraintVariables(const double* a_params){

  for(int pointIndex = 0; pointIndex < m_numPoints; pointIndex++){
    m_paramSpace[pointIndex][0] = a_params[2*pointIndex];
    m_paramSpace[pointIndex][1] = a_params[2*pointIndex + 1];
    m_currentParameters[2*pointIndex] = a_params[2*pointIndex];
    m_currentParameters[2*pointIndex+1] = a_params[2*pointIndex+1];
  }
  m_paramSpace[m_numPoints-1][2] = a_params[2*m_numPoints]; // y'_n
  m_currentParameters[2*m_numPoints] = a_params[2*m_numPoints];
  //m_currentParameters has been set, It doesn't have to be incorrect order, but the param space does

  if(!isSorted(a_params)){
    cout << " Note: SplineFitter - The points have become unsorted..."<< endl;
    sortParamSpace();
  }
  //calculate y"n for m_numPoints-1
  m_paramSpace[m_numPoints-1][3] = 2.0* m_paramSpace[m_numPoints-1][2]/(m_paramSpace[m_numPoints-1][0] - m_paramSpace[m_numPoints-2][0]) -2.0*(m_paramSpace[m_numPoints-1][1] - m_paramSpace[m_numPoints-2][1])/pow(m_paramSpace[m_numPoints-1][0] - m_paramSpace[m_numPoints-2][0],2.0);
  //calculate y'(n-1) then y"n for points 0 through m_numPoints-2
  for(int pointIndex = m_numPoints - 2; pointIndex > 0; pointIndex--){
    m_paramSpace[pointIndex][2] = m_paramSpace[pointIndex + 1][2] + m_paramSpace[pointIndex+1][3]* (m_paramSpace[pointIndex][0] - m_paramSpace[pointIndex + 1][0]); // y2' = y3' + y3"(x2-x3) 
    m_paramSpace[pointIndex][3] = 2.0* m_paramSpace[pointIndex][2]/(m_paramSpace[pointIndex][0] - m_paramSpace[pointIndex - 1][0]) 
                                   -2.0* (m_paramSpace[pointIndex][1] - m_paramSpace[pointIndex - 1][1])/pow(m_paramSpace[pointIndex][0] - m_paramSpace[pointIndex - 1][0],2.0);
                                  // y2" = 2*y2'/(x2 - x1) - 2*(y2-y1)/(x2 - x1)^2
  }




}



bool SplineFitter::isSorted(const double* a_params){
  for(int pointIndex = 0; pointIndex < m_numPoints - 1; pointIndex++){
    if(a_params[2*pointIndex] > a_params[2*(pointIndex+1)]) return false;
  }
  return true;
}

//============================================================================================================



double SplineFitter::splineChiSqr(const double* a_params){

  // x0 is fixed because y0 is enough to constrain y"0
  // parameter order        x0 y0 x1 y1 x2 y2 ... y'n


  bool newParameters = false;
  for(int iii = 0; iii < 2*m_numPoints + 1; iii++){
    if(m_currentParameters[iii] != a_params[iii]){
      newParameters = true;
      //cout << " new Parameters =  " << iii << endl;
      break;
    }
  }

  //change the current parameters if they have been updated
  if(newParameters){
    #ifdef _SPLINEFITTER_DEBUG_
      cout << " Parameters Given: ";
      for(int iii = 0; iii < 2*m_numPoints + 1; iii++) cout << "   " << a_params[iii];
      cout << endl;
    #endif 
    calculateConstraintVariables(a_params);

    #ifdef _SPLINEFITTER_DEBUG_
      printCurrentParamSpace();
    #endif

  }

  double chiSqr = 0.0;
  for(int graphPointIndex = 0; graphPointIndex < m_graphToFit->GetN(); graphPointIndex++){
    double X_ofPoint = m_graphToFit->GetX()[graphPointIndex];
    int pointIndex = -1;
    for(pointIndex = 1; pointIndex < m_numPoints; pointIndex++){
      if(X_ofPoint > m_paramSpace[pointIndex-1][0] && X_ofPoint < m_paramSpace[pointIndex][0]) break;
    }
    if(X_ofPoint < m_paramSpace[0][0]) pointIndex = 1;
    if(X_ofPoint > m_paramSpace[m_numPoints - 1][0]) pointIndex = m_numPoints -1;

    double Y_ofPoint = m_graphToFit->GetY()[graphPointIndex];
    double expected = parabEval(X_ofPoint,m_paramSpace[pointIndex][0],m_paramSpace[pointIndex][1],m_paramSpace[pointIndex][2],m_paramSpace[pointIndex][3]);
    double errorY = 0.0;
    if(expected > Y_ofPoint) errorY = m_graphToFit->GetErrorYhigh(graphPointIndex);
    else errorY = m_graphToFit->GetErrorYlow(graphPointIndex);
    //cout << "  yPoint=" << Y_ofPoint << " expected=" << expected << "   error=" << errorY << endl;
    if(errorY != 0.0) chiSqr += (Y_ofPoint - expected)*(Y_ofPoint - expected)/ (errorY*errorY);  
  }
  if(m_regularizationParForClosePoints > 0.0){
    for(int pointIndex = 1; pointIndex < m_numPoints; pointIndex++){
      for(int pointIndexTwo = pointIndex + 1; pointIndexTwo < m_numPoints; pointIndexTwo++){

        double dist = fabs(m_paramSpace[pointIndex][0] - m_paramSpace[pointIndexTwo][0]);

        /*
        int pointNearest = HistogramUtilities::getPointClosestToVal(m_graphToFit,m_paramSpace[pointIndex][0]);

        if(pointNearest == m_graphToFit->GetN()){
          double ratio = fabs(m_graphToFit->GetX()[pointNearest] - m_graphToFit->GetX()[pointNearest-1])/dist;
          chiSqr += m_regularizationParForClosePoints/ratio;
        }else{
          chiSqr += m_regularizationParForClosePoints*dist/fabs(m_graphToFit->GetX()[pointNearest + 1] - m_graphToFit->GetX()[pointNearest]);
        }
        */

        if(dist < m_minimumDistanceForRegularization){
          chiSqr += (m_minimumDistanceForRegularization/dist - 1.0)*m_regularizationParForClosePoints;
        }
      }
    }
  }
  #ifdef _SPLINEFITTER_DEBUG_
    //if(newParameters) 
    cout << "ChiSquared = " << chiSqr << endl;
  #endif

  return chiSqr;
}

//============================================================================================================

double SplineFitter::splineOperator(double* a_x, double* a_params){
  int pointIndex = -1;
  for(pointIndex = 1; pointIndex < m_numPoints; pointIndex++){
    if(a_x[0] > m_paramSpace[pointIndex-1][0] && a_x[0] < m_paramSpace[pointIndex][0]) break;
  }
  if(a_x[0] < m_paramSpace[0][0]) pointIndex = 1;
  if(a_x[0] > m_paramSpace[m_numPoints - 1][0]) pointIndex = m_numPoints -1;

  return parabEval(a_x[0],m_paramSpace[pointIndex ][0],m_paramSpace[pointIndex][1],m_paramSpace[pointIndex][2],m_paramSpace[pointIndex][3]);
}

//============================================================================================================

double SplineFitter::splineOperatorWithParams(double* a_x, double* a_params){

  // x0 is fixed because y0 is enough to constrain y"0
  // parameter order        x0 y0 x1 y1 x2 y2 ... y'n

  bool newParameters = false;
  for(int iii = 0; iii < 2*m_numPoints + 1; iii++){
    if(m_currentParameters[iii] != a_params[iii]){
      newParameters = true;
      break;
    }
  }

  //change the current parameters if they have been updated
  if(newParameters){
    calculateConstraintVariables(a_params);
  }

  int thePointIndex = -1;
  for(int pointIndex = 1; pointIndex < m_numPoints; pointIndex++){
    if(a_x[0] > m_paramSpace[pointIndex-1][0] && a_x[0] < m_paramSpace[pointIndex][0]){
      thePointIndex = pointIndex;
      break;
    }
  }
  if(a_x[0] < m_paramSpace[0][0]) thePointIndex = 1;
  if(a_x[0] > m_paramSpace[m_numPoints - 1][0]) thePointIndex = m_numPoints -1;

  return parabEval(a_x[0],m_paramSpace[thePointIndex][0],m_paramSpace[thePointIndex][1],m_paramSpace[thePointIndex][2],m_paramSpace[thePointIndex][3]);
}

//============================================================================================================

void SplineFitter::analyzeGraph(int a_numSplinePoints,TF1* &a_fitResult, TGraph* &a_splinePoints, double &a_chiSqr, double &a_numDegreeOfFreedom, bool a_fixXPoints, bool a_findAdditionalSeedingPoints, double* a_forcedXLocs, bool a_fixedEndDerivative, double a_fixedEndDerivativeValue){
  #ifdef _SPLINEFITTER_DEBUG_
    cout << "Starting SplineFitter::analyzeGraph(" << a_numSplinePoints << " , " << a_fitResult << " , " << a_splinePoints << " , " << a_chiSqr << " , "
         << a_numDegreeOfFreedom << " , " << a_fixXPoints << " , " << a_findAdditionalSeedingPoints << " , " << a_forcedXLocs << " , " << a_fixedEndDerivative << " , " << a_fixedEndDerivativeValue << " ) " << endl;
    cout << "  Graph ptr: " << m_graphToFit << endl;
    cout << "  Graph name: " << m_graphToFit->GetName() << endl;
    cout << "  Graph Points: " << m_graphToFit->GetN() << endl;
  #endif
  m_graphToFit->Sort();
  m_additionalSeedingPoints.clear();
  if(a_findAdditionalSeedingPoints){
    #ifdef _SPLINEFITTER_DEBUG_
      cout << "  Finding additional seeding points... "  << endl;
    #endif
    findSeedingPoints();
  }


  m_numPoints = a_numSplinePoints;
  m_numParabolas = m_numPoints - 1;

  double lowRange = m_graphToFit->GetX()[0];
  double highRange = m_graphToFit->GetX()[m_graphToFit->GetN()-1];
  double range = highRange - lowRange;
  double width = range / ((double) m_numParabolas);

  double lowYVal = HistogramUtilities::findSmallestYVal(m_graphToFit);
  double highYVal = HistogramUtilities::findLargestYVal(m_graphToFit);
  double rangeY = highYVal - lowYVal;


  #ifdef _SPLINEFITTER_DEBUG_
    cout << "  Number of Points: " << m_numPoints << endl;
    cout << "  Number of Parabolas: " << m_numParabolas << endl;
    cout << "  Low Range : " << lowRange << "   High Range : " << highRange <<  "  Range: " << range << "   Width: " << width << endl;
    cout << "  Low Value : " << lowYVal << "  High Value : " << highYVal <<  "   Range: " << rangeY << endl;
  #endif


  for(int pointIndex = 0; pointIndex < m_numPoints; pointIndex++){
    if(a_forcedXLocs == NULL) m_paramSpace[pointIndex][0] = lowRange + ((double) pointIndex)*width;
    else m_paramSpace[pointIndex][0] = a_forcedXLocs[pointIndex];
    m_paramSpace[pointIndex][1] = m_graphToFit->Eval(m_paramSpace[pointIndex][0]);
    m_currentParameters[2*pointIndex] =  m_paramSpace[pointIndex][0];
    m_currentParameters[2*pointIndex + 1] =  m_paramSpace[pointIndex][1];
  }

  if(a_findAdditionalSeedingPoints){
    int numAdditionalPoints = (int) m_additionalSeedingPoints.size();
    if(m_numPoints + numAdditionalPoints > 30){
      cout << "SplineFitter::analyzeGraph() - Too many points (" << m_numPoints + numAdditionalPoints << ") for minimization, restricting to 30" << endl;
      numAdditionalPoints = 30 - m_numPoints;
    }
    for(unsigned int iii = 0; iii < numAdditionalPoints; iii++){
      m_paramSpace[(m_numPoints+iii)][0] = m_additionalSeedingPoints[iii];
      m_paramSpace[(m_numPoints+iii)][1] = m_graphToFit->Eval(m_paramSpace[(m_numPoints+iii)][0]);
      m_currentParameters[2*(m_numPoints+iii)] =  m_paramSpace[(m_numPoints+iii)][0];
      m_currentParameters[2*(m_numPoints+iii) + 1] =  m_paramSpace[(m_numPoints+iii)][1];
    }
  
    m_numPoints = a_numSplinePoints + numAdditionalPoints;
    m_numParabolas = m_numPoints - 1;
  }

  //initialize slope at last point
  if(a_fixedEndDerivative){
    #ifdef _SPLINEFITTER_DEBUG_
      cout << "FIXING END DERIVATIVE" << endl;
    #endif
      m_paramSpace[m_numPoints - 1][2] = a_fixedEndDerivativeValue;
  }else{  
    #ifdef _SPLINEFITTER_DEBUG_
      cout << " NOT FIXING END DERIVATIVE" << endl;
    #endif
    m_paramSpace[m_numPoints - 1][2] = ( m_paramSpace[m_numPoints - 1][1] - m_graphToFit->Eval(m_paramSpace[m_numPoints - 1][0] - 0.1*width))/(0.1*width);    
  }
  m_currentParameters[2*m_numPoints] = m_paramSpace[m_numPoints - 1][2];

  calculateConstraintVariables(&m_currentParameters[0]); // sets m_currentParameters
  #ifdef _SPLINEFITTER_DEBUG_
    cout << "Initial Paramter Space:" << endl;
    printCurrentParamSpace();
  #endif

  // reset current parameters if sorted
  /*for(int pointIndex = 0; pointIndex < m_numPoints; pointIndex++){
    m_currentParameters[2*pointIndex] =  m_paramSpace[pointIndex][0];
    m_currentParameters[2*pointIndex + 1] =  m_paramSpace[pointIndex][1];
  }
  m_currentParameters[2*m_numPoints] = m_paramSpace[m_numPoints - 1][2];
  */

  #ifdef _SPLINEFITTER_DEBUG_
    cout << "  Setting Up Minimizer... " << endl;
  #endif
  
  //setup minimizer
  m_minimizer = ROOT::Math::Factory::CreateMinimizer();
  m_chiSqrdFunctor = new ROOT::Math::Functor(this,&SplineFitter::splineChiSqr,(2*m_numPoints) + 1);
  #ifdef _SPLINEFITTER_DEBUG_
    cout << "  Created Fnctor with NDim " << m_chiSqrdFunctor->NDim() << "  functor: " << &SplineFitter::splineChiSqr  << endl;
  #endif
  //fix the first x location
  m_minimizer->SetFunction(*m_chiSqrdFunctor);
  m_minimizer->SetVariable(0, Form("x_%02d",0), m_currentParameters[0], width);
  m_minimizer->FixVariable(0);
  m_minimizer->SetVariable(1, Form("y_%02d", 0), m_currentParameters[1], rangeY*0.1);
  m_minimizer->SetVariableLimits(1, lowYVal - rangeY*0.4, highYVal + rangeY*0.4); 

  //set parameters for the middle points
  for(int pointIndex = 1; pointIndex < m_numPoints - 1; pointIndex++){
    m_minimizer->SetVariable(2*pointIndex, Form("x_%02d",pointIndex), m_currentParameters[2*pointIndex] , width);
    if(a_fixXPoints) m_minimizer->FixVariable(2*pointIndex);
    else             m_minimizer->SetVariableLimits(2*pointIndex, m_currentParameters[0], m_currentParameters[2*(m_numPoints-1)]);
    m_minimizer->SetVariable(2*pointIndex+1, Form("y_%02d", pointIndex), m_currentParameters[2*pointIndex + 1], rangeY*0.1);
    m_minimizer->SetVariableLimits(2*pointIndex+1, lowYVal - rangeY*0.1, highYVal + rangeY*0.1);  
  }
  // fix the last x location
  m_minimizer->SetVariable(2*(m_numPoints-1), Form("x_%02d",m_numPoints - 1), m_currentParameters[2*(m_numPoints -1)]*1.000000001 , width);
  m_minimizer->FixVariable(2*(m_numPoints-1));
  m_minimizer->SetVariable(2*(m_numPoints-1)+1, Form("y_%02d", m_numPoints - 1), m_currentParameters[2*(m_numPoints-1)+1], rangeY*0.1);
  m_minimizer->SetVariableLimits(2*(m_numPoints-1)+1, lowYVal - rangeY*0.1, highYVal + rangeY*0.1);  

  //set the first derivative at the right side
  m_minimizer->SetVariable(2*m_numPoints, Form("dydx_%02d",m_numPoints - 1) ,m_currentParameters[2*m_numPoints], 100.0*rangeY/range);
  if(a_fixedEndDerivative){
    m_minimizer->FixVariable(2*m_numPoints);
  }else{  
    m_minimizer->SetVariableLimits(2*m_numPoints,m_currentParameters[2*m_numPoints] - 100.0*rangeY/range, m_currentParameters[2*m_numPoints] + 100.0*rangeY/range);
  }

  a_numDegreeOfFreedom = (double) m_graphToFit->GetN();
  a_numDegreeOfFreedom = a_numDegreeOfFreedom - ((double)m_minimizer->NFree());


  m_regularizationParForClosePoints = a_numDegreeOfFreedom;

  #ifdef _SPLINEFITTER_DEBUG_
    m_minimizer->SetPrintLevel(5);
  #else
    m_minimizer->SetPrintLevel(1);
  #endif
  m_minimizer->SetMaxIterations(3000);
  m_minimizer->SetMaxFunctionCalls(200000);
  m_minimizer->SetValidError(true);

  #ifdef _SPLINEFITTER_DEBUG_
    cout << "  Parameters Set. Ready to minimize... " << endl;
  #endif


  m_minimizer->Minimize();
  int status = m_minimizer->Status();
  #ifdef _SPLINEFITTER_DEBUG_
    cout << "  STATUS: " << status << endl;
  #endif
  const double* finalParameters = m_minimizer->X();
  m_finalParameters.clear();
  #ifdef _SPLINEFITTER_DEBUG_
    cout << " Final Parameters: ";
  #endif
    for(int iii = 0; iii < 2*m_numPoints + 1; iii++){
      m_finalParameters.push_back(finalParameters[iii]);
      #ifdef _SPLINEFITTER_DEBUG_
        cout << "    " << finalParameters[iii];
      #endif
    }
  #ifdef _SPLINEFITTER_DEBUG_
    cout << endl;
  #endif
    
  calculateConstraintVariables(&m_finalParameters[0]);

  #ifdef _SPLINEFITTER_DEBUG_
    printCurrentParamSpace();
  #endif

  a_chiSqr = splineChiSqr(finalParameters);


  //calculateConstraintVariables(&m_currentParameters[0]);
  //a_fitResult = new TF1("finalSplineFit",this,&SplineFitter::splineOperator,lowRange - 0.1*width,highRange+0.1*width,0,"SplineFitter","splineOperator");
  a_fitResult = new TF1("finalSplineFit",this,&SplineFitter::splineOperatorWithParams,lowRange - 0.1*width,highRange+0.1*width,(2*m_numPoints) + 1,"SplineFitter","splineOperatorWithParams");
  for(int iii = 0; iii < 2*m_numPoints + 1; iii++){
    a_fitResult->SetParameter(iii,m_finalParameters[iii]);
  }

  m_fitResult = a_fitResult;
  cout << "  FitResultPtr: " << a_fitResult << endl;
  a_splinePoints = makeSplinePointGraph();

  cout << " getting covariance matrix " << endl;
  m_minimizer->GetCovMatrix(&m_covarainceMatrix[0]);
  cout << "Making Covariance Matrix" << endl;
  m_covarainceMatrixDSym = new TMatrixDSym((2*m_numPoints) + 1, &m_covarainceMatrix[0]);
  cout << " Covariance Matrix Okay?" << endl;
  m_covarainceMatrixDSym->Print();

}



void SplineFitter::analyzeGraph_WithBestNumPoints(int a_lowNumSplinePoints, int a_highNumSplinePoints, TF1* &a_fitResult, TGraph* &a_splinePoints, TGraph* &a_bestNumPointsGraph, double &a_chiSqr, double &a_numDegreeOfFreedom){

  if(a_highNumSplinePoints - a_lowNumSplinePoints > 20){
    cout << "SplineFitter::analyzeGraph_WithBestNumPoints - You specified too wide a range of points... try something less than 20" << endl;
    return;
  }


  double trialNumPoints[20];
  double chiSqrPlusDof[20];

  int counter = 0;
  for(int points = a_lowNumSplinePoints; points <= a_highNumSplinePoints; points++){
    trialNumPoints[counter] = points;
    TF1* fitResult = NULL;
    TGraph* graphForPoints = NULL;
    double chiSquared = 0.0;
    double numDof = 0.0;
    analyzeGraph(points, fitResult, graphForPoints, chiSquared, numDof,false, false, NULL);
    chiSqrPlusDof[counter] = chiSquared + 2.0*numDof;
    counter++;

    TCanvas* canv = new TCanvas();
    m_graphToFit->Draw("AP");
    fitResult->Draw("same");
    graphForPoints->Draw("PSame");


  }

  a_bestNumPointsGraph = new TGraph(counter, trialNumPoints,chiSqrPlusDof);
  a_bestNumPointsGraph->SetName("BestNumPointsGraph");

  double smallestNumPoints = trialNumPoints[0];
  double smallestValue = chiSqrPlusDof[0];
  for(int iii = 1; iii < counter; iii++){
    if(chiSqrPlusDof[iii] < smallestValue){
      smallestValue = chiSqrPlusDof[iii];
      smallestNumPoints = trialNumPoints[0];
    }
  }

  analyzeGraph(((int) (smallestNumPoints + 0.2)), a_fitResult, a_splinePoints, a_chiSqr, a_numDegreeOfFreedom,false, false, NULL);
    
  return;

}






void SplineFitter::getIntegralAndError(double a_lowBound, double a_highBound, double &a_integral, double &a_integralError, bool a_divideByWidth){

  a_integral = m_fitResult->Integral(a_lowBound, a_highBound);
  a_integralError = m_fitResult->IntegralError(a_lowBound, a_highBound, &m_finalParameters[0], m_covarainceMatrix, 1.E-2);
  if(a_divideByWidth){
    a_integral = a_integral/(a_highBound - a_lowBound);
    a_integralError = a_integralError/(a_highBound - a_lowBound);
  }

  #ifdef _SPLINEFITTER_DEBUG_
    cout << " SplineFitter::getIntegralAndError(" <<  a_lowBound << ", " << a_highBound << ", " << a_integral << ", " << a_integralError << ", " << a_divideByWidth << endl;
  #endif


  return;
}


void SplineFitter::printCurrentParamSpace(){
  cout << "Parameter Space:" << endl;
  cout << "\tx\ty\ty\'\ty\"" << endl;
  for(int pointIndex = 0; pointIndex < m_numPoints; pointIndex++){

    cout << "\t" << m_paramSpace[pointIndex][0] 
         << "\t" << m_paramSpace[pointIndex][1] 
         << "\t" << m_paramSpace[pointIndex][2] 
         << "\t" << m_paramSpace[pointIndex][3]
         << endl;
  }
  

}



TGraph* SplineFitter::makeSplinePointGraph(){
    vector<double> x_vect;
    vector<double> y_vect;
    for(int iii = 0; iii < m_numPoints; iii++){
        x_vect.push_back(m_paramSpace[iii][0]);
        y_vect.push_back(m_paramSpace[iii][1]);
    }
    TGraph* graph = new TGraph(m_numPoints, &x_vect[0], &y_vect[0]);
    graph->SetName("graphSpline");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue+1);
    graph->SetLineColor(kBlue-1);
    graph->SetMarkerSize(0.5);
    return graph;
}





void SplineFitter::sortParamSpace(){

  bool sorted = false;
  while(!sorted){
    sorted = true;
    for(int pointIndex = 1; pointIndex < m_numPoints - 1; pointIndex++){
      if(m_paramSpace[pointIndex][0] > m_paramSpace[pointIndex+1][0]){
        sorted = false;
        double one = m_paramSpace[pointIndex+1][0];
        double two = m_paramSpace[pointIndex+1][1];
        //double three= m_paramSpace[pointIndex+1][2];  //if unsorted, it must be new variables anyway
        //double four = m_paramSpace[pointIndex+1][3];
         
        m_paramSpace[pointIndex+1][0] = m_paramSpace[pointIndex][0];
        m_paramSpace[pointIndex+1][1] = m_paramSpace[pointIndex][1];
        //m_paramSpace[pointIndex+1][2] = m_paramSpace[pointIndex][2];
        //m_paramSpace[pointIndex+1][3] = m_paramSpace[pointIndex][3];

        m_paramSpace[pointIndex][0] = one;
        m_paramSpace[pointIndex][1] = two;
        //m_paramSpace[pointIndex][2] = three;
        //m_paramSpace[pointIndex][3] = four;

      }
    } 
  }

  m_paramSpace[m_numPoints-1][2] = m_currentParameters[2*m_numPoints];


  //cout << "-------Sorted---------" << endl;
   
}




void SplineFitter::findSeedingPoints(){

  TGraphErrors* derivGraph =  HistogramUtilities::makeDerivativeGraph(m_graphToFit);
  TGraphErrors* secondDerivGraph = HistogramUtilities::makeDerivativeGraph(derivGraph);
  TGraphErrors* thirdDerivGraph = HistogramUtilities::makeDerivativeGraph(secondDerivGraph);

  derivGraph->SetTitle("1st Derivative");
  secondDerivGraph->SetTitle("2nd Derivative");
  thirdDerivGraph->SetTitle("3rd Derivative");

  TCanvas* canv = new TCanvas();
  canv->cd();
  gPad->Divide(1,4);
  canv->cd(1);
  m_graphToFit->Draw();
  canv->cd(2);
  derivGraph->Draw();
  canv->cd(3);
  secondDerivGraph->Draw();
  canv->cd(4);
  thirdDerivGraph->Draw();  

  m_additionalSeedingPoints.clear();
  for(int pointIndex = 0; pointIndex < thirdDerivGraph->GetN(); pointIndex++){

    if(fabs(thirdDerivGraph->GetY()[pointIndex]) > 1.5*thirdDerivGraph->GetEY()[pointIndex]){
      m_additionalSeedingPoints.push_back(thirdDerivGraph->GetX()[pointIndex]);
    }

  }
  return;



}









void SplineFitter::reload(TF1* a_previousFit, TMatrixDSym* a_covMatrix){
  if(!a_previousFit){
    cout << " SplineFitter::reload handed a null TF1 pointer..." << endl;
    return;
  }

  if(!a_covMatrix){
    cout << " SplineFitter::reload handed a null TMatrixDSym pointer..." << endl;
    return;
  }

  int numParameters = a_previousFit->GetNpar();
  m_numPoints=(numParameters-1)/2;
  m_numParabolas = m_numPoints - 1;
  double lowRange = 0;
  double highRange = 0;
  a_previousFit->GetRange(lowRange,highRange);

  m_finalParameters.clear();

  for(int iii = 0; iii < numParameters; iii++) m_finalParameters.push_back(0.0);

  #ifdef _SPLINEFITTER_DEBUG_
    cout << " SplineFitter Reload:  params: " << numParameters << " points: " << m_numPoints << "  parabolas: " 
         << m_numParabolas << " lowrange: " << lowRange << " highRange: " << highRange << endl;
    cout << " Fit Name: " << a_previousFit->GetName() << endl;
  #endif

  string previousFitName = a_previousFit->GetName();


  m_fitResult = new TF1(Form("%s_recreated",previousFitName.c_str()),this,&SplineFitter::splineOperatorWithParams,lowRange,highRange,numParameters,"SplineFitter","splineOperatorWithParams");
  //a_fitResult = new TF1("finalSplineFit",this,&SplineFitter::splineOperatorWithParams,lowRange - 0.1*width,highRange+0.1*width,(2*m_numPoints) + 1,"SplineFitter","splineOperatorWithParams");

  //cout << " m_fitResult : " << m_fitResult << endl;

  for(int paramIndex = 0; paramIndex < numParameters; paramIndex++){
    m_finalParameters[paramIndex] = a_previousFit->GetParameter(paramIndex);
    #ifdef _SPLINEFITTER_DEBUG_
      cout << " parameter " << paramIndex << " : " << m_finalParameters[paramIndex] << endl;
    #endif
    m_fitResult->SetParameter(paramIndex,m_finalParameters[paramIndex]);
  }
  
  //cout << "Now Constraining Variables" << endl;
  calculateConstraintVariables(&m_finalParameters[0]);
  //cout << " Done Constraining Variables" << endl;

  int matrixLength = numParameters*numParameters;
  double* matrixPtr = a_covMatrix->GetMatrixArray();
  for(int iii =0; iii < matrixLength; iii++){
    m_covarainceMatrix[iii] = matrixPtr[iii];
  }
  m_covarainceMatrixDSym = new TMatrixDSym((2*m_numPoints) + 1, &m_covarainceMatrix[0]);
  //everything should be loaded to try and calculate integral errors
  #ifdef _SPLINEFITTER_DEBUG_
    cout << " Spline Fit Funct Ptr: " << m_fitResult << endl;
    printCurrentParamSpace();
  #endif


  return;
}



















//###########################################
//                  TESTERS
//############################################



void SplineFitter::testSplineFunction(){

  m_numPoints = 4;
  m_numParabolas = 3;

  const double params[9] = {0.0,1.0, 0.2,2, 0.4,3, 0.6,1, -20.0};
  calculateConstraintVariables(params);
  TF1* fitFunction = new TF1("testSpline",this,&SplineFitter::splineOperator,-0.1,1.0,0,"SplineFitter","splineOperator");
  fitFunction->SetNpx(300);
  TGraph* graph = makeSplinePointGraph();
  TCanvas* canv = new TCanvas();
  fitFunction->Draw();
  graph->Draw("PSame");

}



void SplineFitter::testSplineFitter(TGraphAsymmErrors* a_graph, int a_numPoints, bool a_fixXPoints){
  
  TGraphAsymmErrors* graph = NULL;
  if(a_graph) graph = a_graph;
  else graph = makeTestingGraph();
  setGraphToFit(graph);


  TF1* fitFunction = NULL;
  TGraph* splinePoints = NULL;
  double chiSqr = 0.0;
  double ndof = 0.0;
  analyzeGraph(a_numPoints,fitFunction, splinePoints,chiSqr,ndof,a_fixXPoints,false,NULL);
  fitFunction->SetNpx(300);
  TCanvas* canv = new TCanvas();
  TH1I* frame = new TH1I("frame","TPC Efficiency Proton Spline Fit",1,0,5);
  frame->Draw();
  frame->GetYaxis()->SetRangeUser(-1,1);
  graph->Draw("PSame");
  fitFunction->Draw("same");
  splinePoints->Draw("PSame");

}

void SplineFitter::testSplineFitterBestNumPoints(TGraphAsymmErrors* a_graph){
  
  TGraphAsymmErrors* graph = NULL;
  if(a_graph) graph = a_graph;
  else graph = makeTestingGraph();
  setGraphToFit(graph);


  TF1* fitFunction = NULL;
  TGraph* splinePoints = NULL;
  TGraph* bestPointsGraph = NULL;
  double chiSqr = 0.0;
  double ndof = 0.0;
  analyzeGraph_WithBestNumPoints(3,7,fitFunction, splinePoints, bestPointsGraph, chiSqr, ndof);

  
  fitFunction->SetNpx(300);
  TCanvas* canv = new TCanvas();
  //TH1I* frame = new TH1I("frame","TPC Efficiency Proton Spline Fit",20,a_lowNumSplinePoints,a_highNumSplinePoints);
  //frame->Draw();
  bestPointsGraph->Draw("AP");


}

TGraphAsymmErrors* SplineFitter::makeTestingGraph(){

  vector<double> x(200,0.0);
  vector<double> y(200,0.0);
  vector<double> exl(200,0.0);
  vector<double> eyl(200,0.0);
  vector<double> exh(200,0.0);
  vector<double> eyh(200,0.0);

  x[0]=0.01; y[0]=0; exl[0]=0.01; exh[0]=0.01; eyl[0]=0; eyh[0]=0.000177003;
  x[1]=0.03; y[1]=0.34353; exl[1]=0.01; exh[1]=0.01; eyl[1]=0.00778816; eyh[1]=0.00784468;
  x[2]=0.05; y[2]=0.897678; exl[2]=0.01; exh[2]=0.01; eyl[2]=0.00640153; eyh[2]=0.00617335;
  x[3]=0.07; y[3]=0.927896; exl[3]=0.01; exh[3]=0.01; eyl[3]=0.00646188; eyh[3]=0.00612428;
  x[4]=0.09; y[4]=0.939057; exl[4]=0.01; exh[4]=0.01; eyl[4]=0.00449188; eyh[4]=0.0042945;
  x[5]=0.11; y[5]=0.935041; exl[5]=0.01; exh[5]=0.01; eyl[5]=0.00485981; eyh[5]=0.00464421;
  x[6]=0.13; y[6]=0.937867; exl[6]=0.01; exh[6]=0.01; eyl[6]=0.00548704; eyh[6]=0.00520101;
  x[7]=0.15; y[7]=0.93524; exl[7]=0.01; exh[7]=0.01; eyl[7]=0.00587879; eyh[7]=0.00556518;
  x[8]=0.17; y[8]=0.943654; exl[8]=0.01; exh[8]=0.01; eyl[8]=0.00567529; eyh[8]=0.00533809;
  x[9]=0.19; y[9]=0.951572; exl[9]=0.01; exh[9]=0.01; eyl[9]=0.00557999; eyh[9]=0.00520059;
  x[10]=0.21; y[10]=0.931719; exl[10]=0.01; exh[10]=0.01; eyl[10]=0.00646821; eyh[10]=0.00611043;
  x[11]=0.23; y[11]=0.928807; exl[11]=0.01; exh[11]=0.01; eyl[11]=0.00679615; eyh[11]=0.00641879;
  x[12]=0.25; y[12]=0.922918; exl[12]=0.01; exh[12]=0.01; eyl[12]=0.00719653; eyh[12]=0.00680801;
  x[13]=0.27; y[13]=0.939415; exl[13]=0.01; exh[13]=0.01; eyl[13]=0.00650613; eyh[13]=0.00609752;
  x[14]=0.29; y[14]=0.925626; exl[14]=0.01; exh[14]=0.01; eyl[14]=0.00733463; eyh[14]=0.00691625;
  x[15]=0.31; y[15]=0.921391; exl[15]=0.01; exh[15]=0.01; eyl[15]=0.00761692; eyh[15]=0.00719178;
  x[16]=0.33; y[16]=0.913858; exl[16]=0.01; exh[16]=0.01; eyl[16]=0.00789065; eyh[16]=0.00747691;
  x[17]=0.35; y[17]=0.929467; exl[17]=0.01; exh[17]=0.01; eyl[17]=0.0073985; eyh[17]=0.00694918;
  x[18]=0.37; y[18]=0.918548; exl[18]=0.01; exh[18]=0.01; eyl[18]=0.00799815; eyh[18]=0.00754765;
  x[19]=0.39; y[19]=0.908585; exl[19]=0.01; exh[19]=0.01; eyl[19]=0.00834676; eyh[19]=0.00791334;
  x[20]=0.41; y[20]=0.900621; exl[20]=0.01; exh[20]=0.01; eyl[20]=0.00915287; eyh[20]=0.00867858;
  x[21]=0.43; y[21]=0.904049; exl[21]=0.01; exh[21]=0.01; eyl[21]=0.00898003; eyh[21]=0.00850544;
  x[22]=0.45; y[22]=0.904844; exl[22]=0.01; exh[22]=0.01; eyl[22]=0.00886833; eyh[22]=0.00840104;
  x[23]=0.47; y[23]=0.910941; exl[23]=0.01; exh[23]=0.01; eyl[23]=0.00853255; eyh[23]=0.00806742;
  x[24]=0.49; y[24]=0.897179; exl[24]=0.01; exh[24]=0.01; eyl[24]=0.00940677; eyh[24]=0.00892461;
  x[25]=0.51; y[25]=0.903743; exl[25]=0.01; exh[25]=0.01; eyl[25]=0.00904966; eyh[25]=0.00856952;
  x[26]=0.53; y[26]=0.896768; exl[26]=0.01; exh[26]=0.01; eyl[26]=0.00935744; eyh[26]=0.00888225;
  x[27]=0.55; y[27]=0.893996; exl[27]=0.01; exh[27]=0.01; eyl[27]=0.00967885; eyh[27]=0.00918577;
  x[28]=0.57; y[28]=0.876242; exl[28]=0.01; exh[28]=0.01; eyl[28]=0.0101268; eyh[28]=0.00967341;
  x[29]=0.59; y[29]=0.895735; exl[29]=0.01; exh[29]=0.01; eyl[29]=0.00966277; eyh[29]=0.00916235;
  x[30]=0.61; y[30]=0.877306; exl[30]=0.01; exh[30]=0.01; eyl[30]=0.0101996; eyh[30]=0.00973536;
  x[31]=0.63; y[31]=0.896718; exl[31]=0.01; exh[31]=0.01; eyl[31]=0.00971424; eyh[31]=0.00920338;
  x[32]=0.65; y[32]=0.887774; exl[32]=0.01; exh[32]=0.01; eyl[32]=0.00977364; eyh[32]=0.00930167;
  x[33]=0.67; y[33]=0.88189; exl[33]=0.01; exh[33]=0.01; eyl[33]=0.0103786; eyh[33]=0.00987727;
  x[34]=0.69; y[34]=0.890578; exl[34]=0.01; exh[34]=0.01; eyl[34]=0.0102037; eyh[34]=0.00967582;
  x[35]=0.71; y[35]=0.87666; exl[35]=0.01; exh[35]=0.01; eyl[35]=0.0103694; eyh[35]=0.00989271;
  x[36]=0.73; y[36]=0.876238; exl[36]=0.01; exh[36]=0.01; eyl[36]=0.0106128; eyh[36]=0.010116;
  x[37]=0.75; y[37]=0.876248; exl[37]=0.01; exh[37]=0.01; eyl[37]=0.0106557; eyh[37]=0.0101549;
  x[38]=0.77; y[38]=0.898566; exl[38]=0.01; exh[38]=0.01; eyl[38]=0.00994005; eyh[38]=0.00939526;
  x[39]=0.79; y[39]=0.877734; exl[39]=0.01; exh[39]=0.01; eyl[39]=0.0105813; eyh[39]=0.0100805;
  x[40]=0.81; y[40]=0.869223; exl[40]=0.01; exh[40]=0.01; eyl[40]=0.0108164; eyh[40]=0.0103322;
  x[41]=0.83; y[41]=0.876877; exl[41]=0.01; exh[41]=0.01; eyl[41]=0.0106497; eyh[41]=0.0101465;
  x[42]=0.85; y[42]=0.872383; exl[42]=0.01; exh[42]=0.01; eyl[42]=0.0107852; eyh[42]=0.0102901;
  x[43]=0.87; y[43]=0.89172; exl[43]=0.01; exh[43]=0.01; eyl[43]=0.010405; eyh[43]=0.00985027;
  x[44]=0.89; y[44]=0.868177; exl[44]=0.01; exh[44]=0.01; eyl[44]=0.011111; eyh[44]=0.0106053;
  x[45]=0.91; y[45]=0.868621; exl[45]=0.01; exh[45]=0.01; eyl[45]=0.0113997; eyh[45]=0.010866;
  x[46]=0.93; y[46]=0.873563; exl[46]=0.01; exh[46]=0.01; eyl[46]=0.0110053; eyh[46]=0.0104847;
  x[47]=0.95; y[47]=0.866733; exl[47]=0.01; exh[47]=0.01; eyl[47]=0.0110049; eyh[47]=0.0105148;
  x[48]=0.97; y[48]=0.882664; exl[48]=0.01; exh[48]=0.01; eyl[48]=0.0107356; eyh[48]=0.0101961;
  x[49]=0.99; y[49]=0.865112; exl[49]=0.01; exh[49]=0.01; eyl[49]=0.0111274; eyh[49]=0.0106335;
  x[50]=1.01; y[50]=0.867444; exl[50]=0.01; exh[50]=0.01; eyl[50]=0.0113037; eyh[50]=0.0107841;
  x[51]=1.03; y[51]=0.856838; exl[51]=0.01; exh[51]=0.01; eyl[51]=0.0117028; eyh[51]=0.0111945;
  x[52]=1.05; y[52]=0.871001; exl[52]=0.01; exh[52]=0.01; eyl[52]=0.0110253; eyh[52]=0.0105147;
  x[53]=1.07; y[53]=0.877654; exl[53]=0.01; exh[53]=0.01; eyl[53]=0.0106768; eyh[53]=0.0101675;
  x[54]=1.09; y[54]=0.854839; exl[54]=0.01; exh[54]=0.01; eyl[54]=0.0118061; eyh[54]=0.0112974;
  x[55]=1.11; y[55]=0.869427; exl[55]=0.01; exh[55]=0.01; eyl[55]=0.011241; eyh[55]=0.010718;
  x[56]=1.13; y[56]=0.887605; exl[56]=0.01; exh[56]=0.01; eyl[56]=0.0105114; eyh[56]=0.0099683;
  x[57]=1.15; y[57]=0.854893; exl[57]=0.01; exh[57]=0.01; eyl[57]=0.0120791; eyh[57]=0.0115469;
  x[58]=1.17; y[58]=0.860414; exl[58]=0.01; exh[58]=0.01; eyl[58]=0.0117072; eyh[58]=0.0111831;
  x[59]=1.19; y[59]=0.856223; exl[59]=0.01; exh[59]=0.01; eyl[59]=0.0117484; eyh[59]=0.0112388;
  x[60]=1.21; y[60]=0.848574; exl[60]=0.01; exh[60]=0.01; eyl[60]=0.0114698; eyh[60]=0.0110128;
  x[61]=1.23; y[61]=0.866805; exl[61]=0.01; exh[61]=0.01; eyl[61]=0.0112168; eyh[61]=0.0107078;
  x[62]=1.25; y[62]=0.878988; exl[62]=0.01; exh[62]=0.01; eyl[62]=0.0110976; eyh[62]=0.0105416;
  x[63]=1.27; y[63]=0.871019; exl[63]=0.01; exh[63]=0.01; eyl[63]=0.0137677; eyh[63]=0.0129803;
  x[64]=1.29; y[64]=0.876471; exl[64]=0.01; exh[64]=0.01; eyl[64]=0.0185715; eyh[64]=0.0170975;
  x[65]=1.31; y[65]=0.859756; exl[65]=0.01; exh[65]=0.01; eyl[65]=0.0198846; eyh[65]=0.0184253;
  x[66]=1.33; y[66]=0.875; exl[66]=0.01; exh[66]=0.01; eyl[66]=0.0192548; eyh[66]=0.0176952;
  x[67]=1.35; y[67]=0.838906; exl[67]=0.01; exh[67]=0.01; eyl[67]=0.0209272; eyh[67]=0.0195573;
  x[68]=1.37; y[68]=0.880259; exl[68]=0.01; exh[68]=0.01; eyl[68]=0.0192768; eyh[68]=0.0176389;
  x[69]=1.39; y[69]=0.835664; exl[69]=0.01; exh[69]=0.01; eyl[69]=0.0226584; eyh[69]=0.0210985;
  x[70]=1.41; y[70]=0.869301; exl[70]=0.01; exh[70]=0.01; eyl[70]=0.0193159; eyh[70]=0.0178221;
  x[71]=1.43; y[71]=0.833333; exl[71]=0.01; exh[71]=0.01; eyl[71]=0.0219976; eyh[71]=0.0205495;
  x[72]=1.45; y[72]=0.80339; exl[72]=0.01; exh[72]=0.01; eyl[72]=0.0237797; eyh[72]=0.0224133;
  x[73]=1.47; y[73]=0.864615; exl[73]=0.01; exh[73]=0.01; eyl[73]=0.019708; eyh[73]=0.0182152;
  x[74]=1.49; y[74]=0.848993; exl[74]=0.01; exh[74]=0.01; eyl[74]=0.0214935; eyh[74]=0.0199363;
  x[75]=1.51; y[75]=0.818815; exl[75]=0.01; exh[75]=0.01; eyl[75]=0.0234336; eyh[75]=0.0219577;
  x[76]=1.53; y[76]=0.867508; exl[76]=0.01; exh[76]=0.01; eyl[76]=0.0197965; eyh[76]=0.0182539;
  x[77]=1.55; y[77]=0.869416; exl[77]=0.01; exh[77]=0.01; eyl[77]=0.0205777; eyh[77]=0.0188889;
  x[78]=1.57; y[78]=0.867069; exl[78]=0.01; exh[78]=0.01; eyl[78]=0.0193835; eyh[78]=0.0179077;
  x[79]=1.59; y[79]=0.869091; exl[79]=0.01; exh[79]=0.01; eyl[79]=0.0212112; eyh[79]=0.0194262;
  x[80]=1.61; y[80]=0.860068; exl[80]=0.01; exh[80]=0.01; eyl[80]=0.0210614; eyh[80]=0.019427;
  x[81]=1.63; y[81]=0.863222; exl[81]=0.01; exh[81]=0.01; eyl[81]=0.0196618; eyh[81]=0.0181928;
  x[82]=1.65; y[82]=0.854305; exl[82]=0.01; exh[82]=0.01; eyl[82]=0.0210576; eyh[82]=0.0194974;
  x[83]=1.67; y[83]=0.84918; exl[83]=0.01; exh[83]=0.01; eyl[83]=0.0212275; eyh[83]=0.0197051;
  x[84]=1.69; y[84]=0.83612; exl[84]=0.01; exh[84]=0.01; eyl[84]=0.022123; eyh[84]=0.0206287;
  x[85]=1.71; y[85]=0.893204; exl[85]=0.01; exh[85]=0.01; eyl[85]=0.0184128; eyh[85]=0.0167183;
  x[86]=1.73; y[86]=0.818471; exl[86]=0.01; exh[86]=0.01; eyl[86]=0.0223925; eyh[86]=0.0210443;
  x[87]=1.75; y[87]=0.869427; exl[87]=0.01; exh[87]=0.01; eyl[87]=0.0197811; eyh[87]=0.0182156;
  x[88]=1.77; y[88]=0.864769; exl[88]=0.01; exh[88]=0.01; eyl[88]=0.0212408; eyh[88]=0.0195145;
  x[89]=1.79; y[89]=0.882353; exl[89]=0.01; exh[89]=0.01; eyl[89]=0.0192394; eyh[89]=0.0175762;
  x[90]=1.81; y[90]=0.882562; exl[90]=0.01; exh[90]=0.01; eyl[90]=0.0200983; eyh[90]=0.0182866;
  x[91]=1.83; y[91]=0.850498; exl[91]=0.01; exh[91]=0.01; eyl[91]=0.0213017; eyh[91]=0.0197532;
  x[92]=1.85; y[92]=0.865385; exl[92]=0.01; exh[92]=0.01; eyl[92]=0.0200842; eyh[92]=0.0185261;
  x[93]=1.87; y[93]=0.899642; exl[93]=0.01; exh[93]=0.01; eyl[93]=0.0189408; eyh[93]=0.0170332;
  x[94]=1.89; y[94]=0.819079; exl[94]=0.01; exh[94]=0.01; eyl[94]=0.02274; eyh[94]=0.021345;
  x[95]=1.91; y[95]=0.853333; exl[95]=0.01; exh[95]=0.01; eyl[95]=0.0211836; eyh[95]=0.0196174;
  x[96]=1.93; y[96]=0.896226; exl[96]=0.01; exh[96]=0.01; eyl[96]=0.0179291; eyh[96]=0.0162696;
  x[97]=1.95; y[97]=0.852632; exl[97]=0.01; exh[97]=0.01; eyl[97]=0.0217921; eyh[97]=0.020147;
  x[98]=1.97; y[98]=0.848592; exl[98]=0.01; exh[98]=0.01; eyl[98]=0.0220562; eyh[98]=0.0204245;
  x[99]=1.99; y[99]=0.866242; exl[99]=0.01; exh[99]=0.01; eyl[99]=0.0199681; eyh[99]=0.0184163;
  x[100]=2.01; y[100]=0.887789; exl[100]=0.01; exh[100]=0.01; eyl[100]=0.0189762; eyh[100]=0.0172724;
  x[101]=2.03; y[101]=0.860465; exl[101]=0.01; exh[101]=0.01; eyl[101]=0.0207471; eyh[101]=0.0191543;
  x[102]=2.05; y[102]=0.843866; exl[102]=0.01; exh[102]=0.01; eyl[102]=0.0229458; eyh[102]=0.021247;
  x[103]=2.07; y[103]=0.863309; exl[103]=0.01; exh[103]=0.01; eyl[103]=0.0214481; eyh[103]=0.0197102;
  x[104]=2.09; y[104]=0.843972; exl[104]=0.01; exh[104]=0.01; eyl[104]=0.022388; eyh[104]=0.0207667;
  x[105]=2.11; y[105]=0.914286; exl[105]=0.01; exh[105]=0.01; eyl[105]=0.0177251; eyh[105]=0.0157529;
  x[106]=2.13; y[106]=0.814103; exl[106]=0.01; exh[106]=0.01; eyl[106]=0.0226574; eyh[106]=0.0213193;
  x[107]=2.15; y[107]=0.853659; exl[107]=0.01; exh[107]=0.01; eyl[107]=0.0216558; eyh[107]=0.0200174;
  x[108]=2.17; y[108]=0.872093; exl[108]=0.01; exh[108]=0.01; eyl[108]=0.0217294; eyh[108]=0.0198116;
  x[109]=2.19; y[109]=0.847751; exl[109]=0.01; exh[109]=0.01; eyl[109]=0.021904; eyh[109]=0.0203043;
  x[110]=2.21; y[110]=0.834507; exl[110]=0.01; exh[110]=0.01; eyl[110]=0.0227992; eyh[110]=0.0212339;
  x[111]=2.23; y[111]=0.822289; exl[111]=0.01; exh[111]=0.01; eyl[111]=0.0215955; eyh[111]=0.0203048;
  x[112]=2.25; y[112]=0.871245; exl[112]=0.01; exh[112]=0.01; eyl[112]=0.0229724; eyh[112]=0.0208547;
  x[113]=2.27; y[113]=0.849315; exl[113]=0.01; exh[113]=0.01; eyl[113]=0.0217028; eyh[113]=0.0201123;
  x[114]=2.29; y[114]=0.822742; exl[114]=0.01; exh[114]=0.01; eyl[114]=0.0227663; eyh[114]=0.0213317;
  x[115]=2.31; y[115]=0.845283; exl[115]=0.01; exh[115]=0.01; eyl[115]=0.0230454; eyh[115]=0.0213139;
  x[116]=2.33; y[116]=0.869565; exl[116]=0.01; exh[116]=0.01; eyl[116]=0.0202811; eyh[116]=0.0186368;
  x[117]=2.35; y[117]=0.831715; exl[117]=0.01; exh[117]=0.01; eyl[117]=0.0219655; eyh[117]=0.0205383;
  x[118]=2.37; y[118]=0.833333; exl[118]=0.01; exh[118]=0.01; eyl[118]=0.0231972; eyh[118]=0.0215924;
  x[119]=2.39; y[119]=0.885993; exl[119]=0.01; exh[119]=0.01; eyl[119]=0.0189672; eyh[119]=0.0172934;
  x[120]=2.41; y[120]=0.843648; exl[120]=0.01; exh[120]=0.01; eyl[120]=0.0214453; eyh[120]=0.0199569;
  x[121]=2.43; y[121]=0.836299; exl[121]=0.01; exh[121]=0.01; eyl[121]=0.0228325; eyh[121]=0.021242;
  x[122]=2.45; y[122]=0.836299; exl[122]=0.01; exh[122]=0.01; eyl[122]=0.0228325; eyh[122]=0.021242;
  x[123]=2.47; y[123]=0.83165; exl[123]=0.01; exh[123]=0.01; eyl[123]=0.0224206; eyh[123]=0.0209364;
  x[124]=2.49; y[124]=0.844444; exl[124]=0.01; exh[124]=0.01; eyl[124]=0.0211223; eyh[124]=0.0196682;
  x[125]=2.51; y[125]=0.822148; exl[125]=0.01; exh[125]=0.01; eyl[125]=0.022833; eyh[125]=0.0213964;
  x[126]=2.53; y[126]=0.861818; exl[126]=0.01; exh[126]=0.01; eyl[126]=0.0216593; eyh[126]=0.0199099;
  x[127]=2.55; y[127]=0.908772; exl[127]=0.01; exh[127]=0.01; eyl[127]=0.0180161; eyh[127]=0.0161049;
  x[128]=2.57; y[128]=0.839721; exl[128]=0.01; exh[128]=0.01; eyl[128]=0.0224094; eyh[128]=0.020836;
  x[129]=2.59; y[129]=0.876761; exl[129]=0.01; exh[129]=0.01; eyl[129]=0.0203721; eyh[129]=0.0186071;
  x[130]=2.61; y[130]=0.846154; exl[130]=0.01; exh[130]=0.01; eyl[130]=0.0216076; eyh[130]=0.0200684;
  x[131]=2.63; y[131]=0.87372; exl[131]=0.01; exh[131]=0.01; eyl[131]=0.0202374; eyh[131]=0.0185404;
  x[132]=2.65; y[132]=0.863971; exl[132]=0.01; exh[132]=0.01; eyl[132]=0.0216516; eyh[132]=0.0198723;
  x[133]=2.67; y[133]=0.856667; exl[133]=0.01; exh[133]=0.01; eyl[133]=0.0209984; eyh[133]=0.0194173;
  x[134]=2.69; y[134]=0.867596; exl[134]=0.01; exh[134]=0.01; eyl[134]=0.020838; eyh[134]=0.0191343;
  x[135]=2.71; y[135]=0.868056; exl[135]=0.01; exh[135]=0.01; eyl[135]=0.0207723; eyh[135]=0.0190724;
  x[136]=2.73; y[136]=0.884244; exl[136]=0.01; exh[136]=0.01; eyl[136]=0.0189548; eyh[136]=0.0173101;
  x[137]=2.75; y[137]=0.834437; exl[137]=0.01; exh[137]=0.01; eyl[137]=0.0220929; eyh[137]=0.0206207;
  x[138]=2.77; y[138]=0.847656; exl[138]=0.01; exh[138]=0.01; eyl[138]=0.0233253; eyh[138]=0.0215208;
  x[139]=2.79; y[139]=0.811258; exl[139]=0.01; exh[139]=0.01; eyl[139]=0.0231628; eyh[139]=0.0217932;
  x[140]=2.81; y[140]=0.846405; exl[140]=0.01; exh[140]=0.01; eyl[140]=0.0213382; eyh[140]=0.0198328;
  x[141]=2.83; y[141]=0.853821; exl[141]=0.01; exh[141]=0.01; eyl[141]=0.0211204; eyh[141]=0.0195572;
  x[142]=2.85; y[142]=0.833333; exl[142]=0.01; exh[142]=0.01; eyl[142]=0.0217791; eyh[142]=0.0203587;
  x[143]=2.87; y[143]=0.869258; exl[143]=0.01; exh[143]=0.01; eyl[143]=0.0208873; eyh[143]=0.0191518;
  x[144]=2.89; y[144]=0.852201; exl[144]=0.01; exh[144]=0.01; eyl[144]=0.0206161; eyh[144]=0.019143;
  x[145]=2.91; y[145]=0.833333; exl[145]=0.01; exh[145]=0.01; eyl[145]=0.0234613; eyh[145]=0.021821;
  x[146]=2.93; y[146]=0.839465; exl[146]=0.01; exh[146]=0.01; eyl[146]=0.0219545; eyh[146]=0.0204452;
  x[147]=2.95; y[147]=0.851613; exl[147]=0.01; exh[147]=0.01; eyl[147]=0.0209206; eyh[147]=0.0194122;
  x[148]=2.97; y[148]=0.855019; exl[148]=0.01; exh[148]=0.01; eyl[148]=0.0223143; eyh[148]=0.0205599;
  x[149]=2.99; y[149]=0.895082; exl[149]=0.01; exh[149]=0.01; eyl[149]=0.018406; eyh[149]=0.016681;
  x[150]=3.01; y[150]=0.859649; exl[150]=0.01; exh[150]=0.01; eyl[150]=0.0213899; eyh[150]=0.0197118;
  x[151]=3.03; y[151]=0.866667; exl[151]=0.01; exh[151]=0.01; eyl[151]=0.0222169; eyh[151]=0.0203053;
  x[152]=3.05; y[152]=0.868327; exl[152]=0.01; exh[152]=0.01; eyl[152]=0.0210223; eyh[152]=0.0192789;
  x[153]=3.07; y[153]=0.847973; exl[153]=0.01; exh[153]=0.01; eyl[153]=0.0216232; eyh[153]=0.0200602;
  x[154]=3.09; y[154]=0.808307; exl[154]=0.01; exh[154]=0.01; eyl[154]=0.0228668; eyh[154]=0.0215577;
  x[155]=3.11; y[155]=0.847059; exl[155]=0.01; exh[155]=0.01; eyl[155]=0.0234068; eyh[155]=0.0215984;
  x[156]=3.13; y[156]=0.807971; exl[156]=0.01; exh[156]=0.01; eyl[156]=0.0244037; eyh[156]=0.0229216;
  x[157]=3.15; y[157]=0.844595; exl[157]=0.01; exh[157]=0.01; eyl[157]=0.0218028; eyh[157]=0.020255;
  x[158]=3.17; y[158]=0.872807; exl[158]=0.01; exh[158]=0.01; eyl[158]=0.0231242; eyh[158]=0.0209509;
  x[159]=3.19; y[159]=0.861111; exl[159]=0.01; exh[159]=0.01; eyl[159]=0.0303967; eyh[159]=0.0270746;
  x[160]=3.21; y[160]=0.849398; exl[160]=0.01; exh[160]=0.01; eyl[160]=0.0290771; eyh[160]=0.0262873;
  x[161]=3.23; y[161]=0.834437; exl[161]=0.01; exh[161]=0.01; eyl[161]=0.0316058; eyh[161]=0.028674;
  x[162]=3.25; y[162]=0.87013; exl[162]=0.01; exh[162]=0.01; eyl[162]=0.0286193; eyh[162]=0.0254321;
  x[163]=3.27; y[163]=0.844444; exl[163]=0.01; exh[163]=0.01; eyl[163]=0.0327642; eyh[163]=0.0293887;
  x[164]=3.29; y[164]=0.85; exl[164]=0.01; exh[164]=0.01; eyl[164]=0.0295963; eyh[164]=0.0266976;
  x[165]=3.31; y[165]=0.838323; exl[165]=0.01; exh[165]=0.01; eyl[165]=0.0297433; eyh[165]=0.0270591;
  x[166]=3.33; y[166]=0.8375; exl[166]=0.01; exh[166]=0.01; eyl[166]=0.0304669; eyh[166]=0.0276731;
  x[167]=3.35; y[167]=0.865385; exl[167]=0.01; exh[167]=0.01; eyl[167]=0.028812; eyh[167]=0.0257066;
  x[168]=3.37; y[168]=0.840659; exl[168]=0.01; exh[168]=0.01; eyl[168]=0.0282965; eyh[168]=0.0258148;
  x[169]=3.39; y[169]=0.846667; exl[169]=0.01; exh[169]=0.01; eyl[169]=0.0308539; eyh[169]=0.0277933;
  x[170]=3.41; y[170]=0.858065; exl[170]=0.01; exh[170]=0.01; eyl[170]=0.0294845; eyh[170]=0.026423;
  x[171]=3.43; y[171]=0.82; exl[171]=0.01; exh[171]=0.01; eyl[171]=0.0326567; eyh[171]=0.0298343;
  x[172]=3.45; y[172]=0.879433; exl[172]=0.01; exh[172]=0.01; eyl[172]=0.0291483; eyh[172]=0.0255797;
  x[173]=3.47; y[173]=0.867089; exl[173]=0.01; exh[173]=0.01; eyl[173]=0.028484; eyh[173]=0.0254031;
  x[174]=3.49; y[174]=0.816327; exl[174]=0.01; exh[174]=0.01; eyl[174]=0.0332293; eyh[174]=0.0303832;
  x[175]=3.51; y[175]=0.846154; exl[175]=0.01; exh[175]=0.01; eyl[175]=0.030267; eyh[175]=0.0273277;
  x[176]=3.53; y[176]=0.835366; exl[176]=0.01; exh[176]=0.01; eyl[176]=0.03022; eyh[176]=0.0275112;
  x[177]=3.55; y[177]=0.888889; exl[177]=0.01; exh[177]=0.01; eyl[177]=0.0270586; eyh[177]=0.0236838;
  x[178]=3.57; y[178]=0.852761; exl[178]=0.01; exh[178]=0.01; eyl[178]=0.029112; eyh[178]=0.0262435;
  x[179]=3.59; y[179]=0.852564; exl[179]=0.01; exh[179]=0.01; eyl[179]=0.0298004; eyh[179]=0.0268058;
  x[180]=3.61; y[180]=0.895349; exl[180]=0.01; exh[180]=0.01; eyl[180]=0.0248482; eyh[180]=0.021793;
  x[181]=3.63; y[181]=0.854167; exl[181]=0.01; exh[181]=0.01; eyl[181]=0.0309473; eyh[181]=0.0276902;
  x[182]=3.65; y[182]=0.797468; exl[182]=0.01; exh[182]=0.01; eyl[182]=0.0330863; eyh[182]=0.0305957;
  x[183]=3.67; y[183]=0.851351; exl[183]=0.01; exh[183]=0.01; eyl[183]=0.0307218; eyh[183]=0.0275777;
  x[184]=3.69; y[184]=0.880795; exl[184]=0.01; exh[184]=0.01; eyl[184]=0.0279925; eyh[184]=0.0246465;
  x[185]=3.71; y[185]=0.869048; exl[185]=0.01; exh[185]=0.01; eyl[185]=0.0274291; eyh[185]=0.0245146;
  x[186]=3.73; y[186]=0.850932; exl[186]=0.01; exh[186]=0.01; eyl[186]=0.0294331; eyh[186]=0.0265445;
  x[187]=3.75; y[187]=0.827381; exl[187]=0.01; exh[187]=0.01; eyl[187]=0.0303515; eyh[187]=0.0277703;
  x[188]=3.77; y[188]=0.861635; exl[188]=0.01; exh[188]=0.01; eyl[188]=0.0288205; eyh[188]=0.0258052;
  x[189]=3.79; y[189]=0.893333; exl[189]=0.01; exh[189]=0.01; eyl[189]=0.0269158; eyh[189]=0.0234332;
  x[190]=3.81; y[190]=0.832215; exl[190]=0.01; exh[190]=0.01; eyl[190]=0.0319757; eyh[190]=0.0290249;
  x[191]=3.83; y[191]=0.86; exl[191]=0.01; exh[191]=0.01; eyl[191]=0.0298422; eyh[191]=0.0266619;
  x[192]=3.85; y[192]=0.857988; exl[192]=0.01; exh[192]=0.01; eyl[192]=0.0281888; eyh[192]=0.0253798;
  x[193]=3.87; y[193]=0.870748; exl[193]=0.01; exh[193]=0.01; eyl[193]=0.0292743; eyh[193]=0.0259307;
  x[194]=3.89; y[194]=0.859756; exl[194]=0.01; exh[194]=0.01; eyl[194]=0.0285009; eyh[194]=0.0255923;
  x[195]=3.91; y[195]=0.852071; exl[195]=0.01; exh[195]=0.01; eyl[195]=0.0286183; eyh[195]=0.0258564;
  x[196]=3.93; y[196]=0.844444; exl[196]=0.01; exh[196]=0.01; eyl[196]=0.0327642; eyh[196]=0.0293887;
  x[197]=3.95; y[197]=0.812903; exl[197]=0.01; exh[197]=0.01; eyl[197]=0.0325372; eyh[197]=0.0298662;
  x[198]=3.97; y[198]=0.831169; exl[198]=0.01; exh[198]=0.01; eyl[198]=0.0315023; eyh[198]=0.0286555;
  x[199]=3.99; y[199]=0.879747; exl[199]=0.01; exh[199]=0.01; eyl[199]=0.0274237; eyh[199]=0.0242341;
  
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(200,&x[0],&y[0],&exl[0],&exh[0],&eyl[0],&eyh[0]);
  graph->SetName("protonMidRapTPCEff");
  graph->SetMarkerStyle(21);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  
  return graph;

}






#endif




