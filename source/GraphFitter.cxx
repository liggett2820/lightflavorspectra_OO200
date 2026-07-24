#ifndef _GRAPH_FITTER_CXX_
#define _GRAPH_FITTER_CXX_

#include "../headers/GraphFitter.h"

GraphFitter::GraphFitter(){
  m_graphToFit = NULL;
  m_fitFunction = NULL;
  m_fitResult = NULL;
  m_numParameters = 0;
  m_status = 0;
  chiSqr = 0.0;
  chiSqrPerDOF = 0.0;
  m_covMat = NULL;
  

  m_lowFitRange = -999;
  m_highFitRange = 999;
  m_firstPointInRange = -1;
  m_lastPointInRange = -1;
}

void GraphFitter::loadGraphAsym(TGraphAsymmErrors* a_graph){
  if(!a_graph || a_graph->GetN() <= 0){
    cout << "ERROR: GraphFitter::loadGraph, you loaded a NULL or graph with zero points" << endl;
    return;
  }
  m_graphToFit = a_graph;
  m_firstPointInRange = 0;
  m_lastPointInRange  = m_graphToFit->GetN();
  m_lowFitRange       = m_graphToFit->GetX()[m_firstPointInRange] - m_graphToFit->GetEXlow()[m_firstPointInRange];
  m_highFitRange      = m_graphToFit->GetX()[m_lastPointInRange]  + m_graphToFit->GetEXhigh()[m_lastPointInRange];
};


void GraphFitter::loadGraph(TGraphErrors* a_graph){
	if(!a_graph || a_graph->GetN() <= 0){
    cout << "ERROR: GraphFitter::loadGraph, you loaded a NULL or graph with zero points" << endl;
    return;
  }
  m_graphToFit = HistogramUtilities::convertNoAsymToAsym(a_graph);
  m_firstPointInRange = 0;
  m_lastPointInRange  = m_graphToFit->GetN();
  m_lowFitRange       = m_graphToFit->GetX()[m_firstPointInRange] - m_graphToFit->GetEXlow()[m_firstPointInRange];
  m_highFitRange      = m_graphToFit->GetX()[m_lastPointInRange]  + m_graphToFit->GetEXhigh()[m_lastPointInRange];
}

void GraphFitter::loadTF1(TF1* a_funct, bool a_isSimpleFunctString){
  if(a_isSimpleFunctString){
	  m_fitFunction = (TF1*) a_funct->Clone();
	  m_fitFunction->SetName("GraphFitter_InternalCopyFunction");
  }else{
    m_fitFunction = a_funct;
  }
  m_numParameters = a_funct->GetNpar();

};


vector<double> GraphFitter::getResultFunctionParameters(){
  //vector<double> final_params;
  //for(int iii = 0; iii < m_numParameters; iii++){
  //  final_params.push_back(m_fitResult->GetParameter(iii));
  //}
  return m_final_parameters;
}


void GraphFitter::setFitRange(double a_lowRange, double a_highRange){
  if(!m_graphToFit){
    cout << "ERROR: You must loadGraph before calling GraphFitter::setFitRange" << endl;
    return;
  }

  m_firstPointInRange = HistogramUtilities::getPointClosestToXVal(m_graphToFit,a_lowRange);
  m_lastPointInRange  = HistogramUtilities::getPointClosestToXVal(m_graphToFit,a_highRange);
  m_lowFitRange       = m_graphToFit->GetX()[m_firstPointInRange] - m_graphToFit->GetEXlow()[m_firstPointInRange];
  m_highFitRange      = m_graphToFit->GetX()[m_lastPointInRange]  + m_graphToFit->GetEXhigh()[m_lastPointInRange];

  if(m_firstPointInRange > m_lastPointInRange){
    cout << "WARNING: GraphFitter::setFitRange Your graph is not sorted. I am sorting it for you." << endl;
    m_graphToFit->Sort();
    m_firstPointInRange = HistogramUtilities::getPointClosestToXVal(m_graphToFit,a_lowRange);
    m_lastPointInRange  = HistogramUtilities::getPointClosestToXVal(m_graphToFit,a_highRange);
    m_lowFitRange       = m_graphToFit->GetX()[m_firstPointInRange] - m_graphToFit->GetEXlow()[m_firstPointInRange];
    m_highFitRange      = m_graphToFit->GetX()[m_lastPointInRange] + m_graphToFit->GetEXhigh()[m_lastPointInRange];
  }

  if(m_firstPointInRange >= m_lastPointInRange){
    cout << "ERROR: GraphFitter::setFitRange, too few points to fit..." << endl;
    return;
  }
}



double GraphFitter::chiSqrFunction(const double* a_param){
  if(std::isnan(a_param[0])) return a_param[0];
  #ifdef _GRAPH_FITTER_DEBUG_
    cout << "GraphFitter::chiSqrFunction(";
  #endif
  for(int paramIndex = 0; paramIndex < m_numParameters; paramIndex++){
  	m_fitFunction->SetParameter(paramIndex, a_param[paramIndex]);
    #ifdef _GRAPH_FITTER_DEBUG_
      cout << a_param[paramIndex] << ", ";
    #endif
  }
  #ifdef _GRAPH_FITTER_DEBUG_
    cout << ")";
  #endif  
  double chiSqr = 0.0;
  for(int pointIndex = 0; pointIndex < m_graphToFit->GetN(); pointIndex++){
  	double data_Y       = m_graphToFit->GetY()[pointIndex];
  	double data_lowX    = m_graphToFit->GetX()[pointIndex] - m_graphToFit->GetEXlow()[pointIndex];
  	double data_highX   = m_graphToFit->GetX()[pointIndex] + m_graphToFit->GetEXhigh()[pointIndex];
  	double fit_integral = m_fitFunction->Integral(data_lowX,data_highX)/(data_highX - data_lowX);
  	double error        = m_graphToFit->GetEYhigh()[pointIndex];
  	if(fit_integral < data_Y) error = m_graphToFit->GetEYlow()[pointIndex];
    //if(error == 0.0) continue;
    double chiSqr_to_add = pow((fit_integral - data_Y)/error,2.0);
  	if(std::isfinite(chiSqr_to_add)) chiSqr += chiSqr_to_add;
    #ifdef _GRAPH_FITTER_DEBUG_
      cout << "    " << pointIndex << ": [(" << fit_integral << " - " << data_Y << ")/" << error << "]^2 -> x^2 = " << chiSqr << endl;
    #endif  
  }
  return chiSqr;
}


void GraphFitter::minimize(){
  m_final_parameters.clear();
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer();
  ROOT::Math::Functor* chiSqrdFunctor = new ROOT::Math::Functor(this,&GraphFitter::chiSqrFunction,m_numParameters);
  #ifdef _GRAPH_FITTER_DEBUG_
  cout << "created chi squared functor with NDim " << chiSqrdFunctor->NDim() << "  functor: " << &GraphFitter::chiSqrFunction  << endl;
  #endif
  minimizer->SetFunction(*chiSqrdFunctor);
  for(int paramIndex = 0; paramIndex < m_numParameters; paramIndex++){
    bool fixTheParameter = false;
    for(unsigned int iii = 0; iii < m_paramsToFix.size(); iii++){
      if(m_paramsToFix[iii] == paramIndex){
        fixTheParameter = true;
        break;
      }
    }
    if(fixTheParameter){
      #ifdef _GRAPH_FITTER_DEBUG_
        cout << "Fixing Parameter " << paramIndex << "  " <<  m_fitFunction->GetParName(paramIndex) << " to " << m_fitFunction->GetParameter(paramIndex) << endl;
      #endif
      minimizer->SetFixedVariable(paramIndex, m_fitFunction->GetParName(paramIndex), m_fitFunction->GetParameter(paramIndex));
    }else{
  	  double lowerLimit = 0.0;
  	  double upperLimit = 1.0;
  	  m_fitFunction->GetParLimits(paramIndex,lowerLimit,upperLimit);
  	  double delta = 0.1*(upperLimit - lowerLimit);
  	  minimizer->SetVariable(paramIndex, m_fitFunction->GetParName(paramIndex), m_fitFunction->GetParameter(paramIndex) , delta);
      minimizer->SetVariableLimits(paramIndex, lowerLimit, upperLimit);
      #ifdef _GRAPH_FITTER_DEBUG_
        cout << "  Minimizer Parameter " << paramIndex << " = " << m_fitFunction->GetParameter(paramIndex) << " with range [" << lowerLimit << "," << upperLimit << "]" << " and step size = " << delta << endl; 
      #endif
    }
  }

  //for(unsigned int iii = 0; iii < m_paramsToFix.size(); iii++){
  //	minimizer->FixVariable(m_paramsToFix[iii]);
  //}
  
  #ifdef _GRAPH_FITTER_DEBUG_
    minimizer->SetPrintLevel(5);
  #else
    minimizer->SetPrintLevel(0);
  #endif
  minimizer->SetMaxIterations(3000);
  minimizer->SetMaxFunctionCalls(100000);
  minimizer->SetValidError(true);
  minimizer->Minimize();
  m_status = minimizer->Status();
  if(m_status != 0){
    cout << "ERROR: GraphFitter::minimize() Graph Fit Unsuccessful with Status = " << m_status << endl;
    m_graphToFit->Print();
    m_fitFunction->Print();
    return;
  }

  minimizer->Hesse();
  
  
  const double* finalParams = minimizer->X();
  m_fitResult = (TF1*) m_fitFunction->Clone();
  for(int paramIndex = 0; paramIndex < m_numParameters; paramIndex++){
    m_final_parameters.push_back(finalParams[paramIndex]);
    m_fitResult->SetParameter(paramIndex, finalParams[paramIndex]);
  }

  double covMat[m_numParameters][m_numParameters];
  minimizer->GetCovMatrix(&covMat[0][0]);
  m_covMat = new TMatrixDSym(m_numParameters,&covMat[0][0]);

  #ifdef _GRAPH_FITTER_DEBUG_
    cout << " GraphFitter::minimize    numParameters: " << m_numParameters << "  size of matrix: " << m_covMat->GetNrows() << endl;
    cout << "   matrix array filled: " << endl;
    for(int col = 0; col < m_numParameters; col++){
      cout << "     ";
      for(int row = 0; row < m_numParameters; row++){
        cout << covMat[row][col] << " ";
      }
      cout << endl;
    }
    cout << "   resulting array filled: " << endl;
    for(int col = 0; col < m_numParameters; col++){
      cout << "     ";
      for(int row = 0; row < m_numParameters; row++){
        cout << covMat[row][col] << " ";
      }
      cout << endl;
    }
  #endif
}









#endif