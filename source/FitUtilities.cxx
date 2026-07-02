

#include "../headers/FitUtilities.h"


using namespace FitUtilities;

TF1* spectraFitToScale_FitUtil;

double scaleSpectraByInvFactor(double *x, double *par){
  double scale = 6.2831853072 * (x[0] + par[0]);

  return spectraFitToScale_FitUtil->Eval(x[0]) * scale;
}


double scalingOperator(double a_x, double a_mass, TF1* a_function){
  double scale = 6.2831853072 * (a_x + a_mass);
  return a_function->Eval(a_x) * scale;
}

void FitUtilities::bootstrapIntegral(double &a_integral, double &a_integralError,
              TF1* a_fitFunction, TFitResult* a_fitResult, double a_mass,
              double a_lowIntegralBound, double a_highIntegralBound){
  int numTrials = 10000;


  if(!a_fitFunction || !a_fitResult){
    cout << "NULL Pointer passed to bootstrapIntegral " << endl;
    cout << "Fit Function: " << a_fitFunction << "  Fit Result: " << a_fitResult << endl;
    a_integral = 0.0;
    a_integralError = 0.0;
    return;
  }

  double lowBound = 0;
  double highBound = 0;
  a_fitFunction->GetRange(lowBound, highBound);
  if(!(a_lowIntegralBound == -999 && a_highIntegralBound == -999)){
    lowBound = a_lowIntegralBound;
    highBound = a_highIntegralBound;
  }

  TF1* funct = (TF1*) a_fitFunction->Clone();
  funct->SetName("initalFunctCopy");
  funct->SetTitle("initalFunctCopy");


  //Decompose the Matrix
  TMatrixDSym covMatrix = a_fitResult->GetCovarianceMatrix();
  cout << "Covaraince Matrix: " << endl;
  int numParams = covMatrix.GetNcols();
  cout << "NPar: " << funct->GetNpar() << " matrix size: " << numParams << endl;
  covMatrix.Print();

  TDecompChol cholDecomp(covMatrix);
  TMatrixD upperTriangular = cholDecomp.GetU();
  cout << "Upper Triangular Matrix: " << endl;
  upperTriangular.Print();
  
  TRandom3 randGen;
  randGen.SetSeed(time(NULL));
  TVectorD baseParameters(numParams);
  for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++){
    baseParameters[paramIndex] = funct->GetParameter(paramIndex);
  }
  

  string functName = a_fitFunction->GetName();
  string functString =  functName + " * 6.2831853072 * (x + " + std::to_string(a_mass) + " )";
  cout << "Function String: " << functString << endl;
  TF1* mTScaledFunct = new TF1("mTScaledFunct",functString.c_str(),lowBound,highBound);

  double baseIntegral = mTScaledFunct->Integral(lowBound,highBound);

  double sumOfIntegrals = 0.0;
  double varianceSum = 0.0;
  vector<double> integralsArray;
  
  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    //setup bootstraping
    TVectorD gausDistArray(numParams);
    for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++){
      gausDistArray[paramIndex] = randGen.Gaus();
    }

    TVectorD parameters = upperTriangular*gausDistArray + baseParameters;

    for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++){
      mTScaledFunct->SetParameter(paramIndex, parameters[paramIndex]);
    }
    double integral = mTScaledFunct->Integral(lowBound,highBound);
    sumOfIntegrals += integral;
    integralsArray.push_back( integral );
    //varianceSum +=  (integral - baseIntegral) * (integral - baseIntegral);
  }

  double avgIntegral = sumOfIntegrals/ ((double) numTrials);
  

  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    varianceSum +=  (integralsArray[trialIndex] - avgIntegral) * (integralsArray[trialIndex] - avgIntegral);
  }
  double stdDev = sqrt( varianceSum / ((double) (numTrials - 1)) );

  cout << "Base Integral: " << baseIntegral << "    Avg Integral: " << avgIntegral  << "  StdDev: " << stdDev << endl;

  a_integral = baseIntegral;
  a_integralError = stdDev;
  return;
}


void FitUtilities::bootstrapIntegralMultiParameter(double &a_integral, double &a_avgIntegral, double &a_integralError,TF1* a_fitFunction, TFitResult* a_fitResult, double a_mass, int a_numTrials,double a_lowIntegralBound, double a_highIntegralBound){


  int numTrials = a_numTrials;


  if(!a_fitFunction){
    cout << "NULL Pointer passed to bootstrapIntegral " << endl;
    cout << "Fit Function: " << a_fitFunction << "  Fit Result: " << a_fitResult << endl;
    a_integral = 0.0;
    a_integralError = 0.0;
    return;
  }

  double lowBound = 0;
  double highBound = 0;
  a_fitFunction->GetRange(lowBound, highBound);
  if(!(a_lowIntegralBound == -999 && a_highIntegralBound == -999)){
    lowBound = a_lowIntegralBound;
    highBound = a_highIntegralBound;
  }

  TF1* funct = (TF1*) a_fitFunction->Clone();
  funct->SetName("initalFunctCopy");
  funct->SetTitle("initalFunctCopy");


  //Decompose the Matrix
  vector< int > liveParams;
  TMatrixDSym covMatrix = a_fitResult->GetCovarianceMatrix();
  TMatrixDSym redMatrix = PhysMath::reduceCovMatrix(covMatrix);

  cout << "Live Parameters: ";
  for(int iii = 0; iii < covMatrix.GetNrows(); iii++){
    if(covMatrix[iii][iii] != 0.0){
      liveParams.push_back(iii);
      cout << "  " << iii;
    }
  }
  cout << endl;

  cout << "Covaraince Matrix: " << endl;
  int numTotParams = covMatrix.GetNcols();
  int numLiveParams = redMatrix.GetNcols();
  cout << "NTotPar: " << funct->GetNpar() << " matrix size: " << numTotParams << endl;
  covMatrix.Print();

  cout << "Reduced Matrix: " << endl;
  redMatrix.Print();

  TDecompChol cholDecomp(redMatrix);
  cholDecomp.Decompose();
  TMatrixD lowerTriangular = cholDecomp.GetU(); // U is upper
  lowerTriangular.T();
  cout << "Lower Triangular Matrix: " << endl;
  lowerTriangular.Print();
  
  TRandom3 randGen;
  randGen.SetSeed(time(NULL));
  TVectorD baseParameters(numTotParams);
  for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++){
    baseParameters[paramIndex] = funct->GetParameter(paramIndex);
  }
    

  //string functName = a_fitFunction->GetName();
  //string functString =  functName + " * 6.2831853072 * (x + " + std::to_string(a_mass) + " )";
  //cout << "Function String: " << functString << endl;
  //TF1* mTScaledFunct = new TF1("mTScaledFunct",functString.c_str(),lowBound,highBound);
  spectraFitToScale_FitUtil = funct;
  //TF1* mTScaledFunct = new TF1(Form("%s_Scaled",spectraFitToScale_FitUtil->GetName()), scaleSpectraByInvFactor,0.0,10,1);
  //mTScaledFunct->SetParameter(0,a_mass);

  //double baseIntegral = mTScaledFunct->Integral(lowBound,highBound);
  double baseIntegral = PhysMath::integralWithMtScale(funct,a_mass, 0.01, 0.0, 30.0);



  double sumOfIntegrals = 0.0;
  double varianceSum = 0.0;
  vector<double> integralsArray;
  
  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    //setup bootstraping
    TVectorD gausDistArray(numTotParams);
    for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++){
      gausDistArray[paramIndex] = randGen.Gaus();
    }

    TVectorD parameters = baseParameters;

    for(unsigned int liveParamIndex = 0; liveParamIndex < numLiveParams; liveParamIndex++){
      int liveParam = liveParams[liveParamIndex];
      parameters[liveParam] = 0;
      for(unsigned int liveParamIndex2 = 0; liveParamIndex2 < numLiveParams; liveParamIndex2++){
        int liveParam2 = liveParams[liveParamIndex2];
        parameters[liveParam] += lowerTriangular[liveParamIndex][liveParamIndex2]*gausDistArray[liveParamIndex2];
      }
      parameters[liveParam] += baseParameters[liveParam];
    }

    /*cout << "  New Parameters: ";
    for(int iii = 0; iii < numTotParams; iii++){
      cout << " " << parameters[iii];
    } 
    cout << endl;
    */



    for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++){
      funct->SetParameter(paramIndex, parameters[paramIndex]);
    }

    double integral = PhysMath::integralWithMtScale(funct,a_mass, 0.001, 0.0, 35.0);;
    sumOfIntegrals += integral;
    integralsArray.push_back( integral );
    //varianceSum +=  (integral - baseIntegral) * (integral - baseIntegral);
  }

  double avgIntegral = sumOfIntegrals/ ((double) numTrials);
  

  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    varianceSum +=  (integralsArray[trialIndex] - avgIntegral) * (integralsArray[trialIndex] - avgIntegral);
  }
  double stdDev = sqrt( varianceSum / ((double) (numTrials - 1)) );

  cout << "Base Integral: " << baseIntegral << "    Avg Integral: " << avgIntegral  << "  StdDev: " << stdDev << endl;

  a_integral = baseIntegral;
  a_avgIntegral = avgIntegral;
  a_integralError = stdDev;
  return;
}






vector<double> FitUtilities::bootstrapIntegralSingleParameter(double &a_integral, double &a_avgIntegral, double &a_integralError,TF1* a_fitFunction, 
     TFitResult* a_fitResult, double a_mass, int a_parameterToVary,int a_numTrials, double a_paramError,double a_lowIntegralBound, double a_highIntegralBound, bool a_extrapolateToInf){
  int numTrials = a_numTrials;


  if(!a_fitFunction){
    cout << "NULL Pointer passed to bootstrapIntegral " << endl;
    cout << "Fit Function: " << a_fitFunction << "  Fit Result: " << a_fitResult << endl;
    a_integral = 0.0;
    a_integralError = 0.0;
    vector<double> emptyVect;
    return emptyVect;
  }

  double lowBound = 0;
  double highBound = 0;
  a_fitFunction->GetRange(lowBound, highBound);
  if(!(a_lowIntegralBound == -999 && a_highIntegralBound == -999)){
    lowBound = a_lowIntegralBound;
    highBound = a_highIntegralBound;
  }

  TString *fitFunctName = new TString(Form("%s",a_fitFunction->GetName()));
  cout << "fitFunctName: " << fitFunctName->Data() << endl;

  //Make a copy of the function
  TF1* funct;
  if( fitFunctName->Contains("bw",TString::kIgnoreCase) || fitFunctName->Contains("blast",TString::kIgnoreCase)){
    cout << "FitUtilities::bootstrapIntegralSingleParameter blastWave found" << endl;
    funct = new TF1("initalFunctCopy",PhysMath::blastWaveModelFit,0,10,6);
    funct->SetParameter(0,a_fitFunction->GetParameter(0));
    funct->SetParameter(1,a_fitFunction->GetParameter(1));
    funct->SetParameter(2,a_fitFunction->GetParameter(2));
    funct->FixParameter(3,a_fitFunction->GetParameter(3));
    funct->FixParameter(4,a_fitFunction->GetParameter(4));
    funct->FixParameter(5,a_fitFunction->GetParameter(5));
    cout << "set bw pars" << endl;
  }else{
    funct = (TF1*) a_fitFunction->Clone();
    funct->SetName("initalFunctCopy");
  }
  int numParams = funct->GetNpar();

  //determine the parameter error, given or grabbed from TFitResult
  double paramError = a_paramError;
  cout << "passed error: " << a_paramError << endl;
  if(a_fitResult){
    paramError = sqrt(a_fitResult->GetCovarianceMatrix()[a_parameterToVary][a_parameterToVary]);
  }

  cout << "paramError: " << paramError << endl;

  TRandom3 randGen;
  randGen.SetSeed(time(NULL));
  TVectorD baseParameters(numParams);
  for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++) baseParameters[paramIndex] = funct->GetParameter(paramIndex);
    //changed by Ben Aug 16 2021 to use passed bounds
//  a_integral = PhysMath::integralWithMtScale(funct,a_mass, 0.0005, 0.0, 5.0,true);
  //cout << "integtal1" << endl;
  a_integral = PhysMath::integralWithMtScale(funct,a_mass, 0.0005, lowBound, highBound,a_extrapolateToInf);
  //cout << "integral1 val: " << a_integral << endl;

  double sumOfIntegrals = 0.0;
  double sumOfParameter = 0.0;
  
  vector<double> integralsArray; // stores the integrals computed
  vector<double> parameterVaryArray; // stores the random parameter used for bootstrap
  
  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    //setup bootstraping
    TVectorD gausDistArray(numParams);
    TVectorD parameters(numParams);
    for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++){
      gausDistArray[paramIndex] = randGen.Gaus();
      parameters[paramIndex] = baseParameters[paramIndex];
    }
    
    //Randomly sample the parameter of interest
    parameters[a_parameterToVary] = randGen.Gaus(baseParameters[a_parameterToVary], paramError);
    sumOfParameter += parameters[a_parameterToVary];
    parameterVaryArray.push_back(parameters[a_parameterToVary]);
    for(int paramIndex = 0; paramIndex < funct->GetNpar(); paramIndex++) funct->SetParameter(paramIndex, parameters[paramIndex]);
    //changed by Ben Aug 16 2021 to use passed bounds
//    double integral = PhysMath::integralWithMtScale(funct,a_mass, 0.0005, 0.0, 5.0,true);
     // cout << "trial " << trialIndex << " integral: ";
    double integral = PhysMath::integralWithMtScale(funct,a_mass, 0.0005, lowBound, highBound,a_extrapolateToInf);
    //cout << integral << endl;
    sumOfIntegrals += integral;
    integralsArray.push_back( integral );
  }

  a_avgIntegral = sumOfIntegrals/ ((double) numTrials);
  double avgParam = sumOfParameter/ ((double) numTrials);

  double varianceParamSum = 0.0;
  double varianceSum = 0.0;
  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    varianceParamSum += (parameterVaryArray[trialIndex] - avgParam) * (parameterVaryArray[trialIndex] - avgParam);  
    varianceSum +=  (integralsArray[trialIndex] - a_avgIntegral) * (integralsArray[trialIndex] - a_avgIntegral);
  }
  a_integralError = sqrt( varianceSum / ((double) (numTrials - 1)) );
  double calcParamError = sqrt( varianceParamSum / ((double) (numTrials - 1)) );

  cout << "  Base Integral: " << a_integral << "  Avg Integral: " << a_avgIntegral  << "  StdDev: " << a_integralError << endl;
  cout << "  Base Parameter: " << baseParameters[a_parameterToVary] << " Avg Parameter: " << avgParam <<  "  Given Err: " << a_paramError << "  Calc Err: " << calcParamError << endl;
  return integralsArray;
}








void FitUtilities::testBootstrapper(){
  TF1* expFunct = new TF1("expFunct","[0]/([1]*6.2831853072*(x + 1.0))*exp(-x/[1])",0.0,50.0);
  expFunct->SetParNames("A","T");
  expFunct->SetNpx(300);
  expFunct->SetParameter(0,500.0);
  expFunct->SetParameter(1,5.0);


  TF1* gausFunct = new TF1("gausFunct","[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",-20.0,20.0);
  gausFunct->SetParNames("N","#mu","#sig");
  gausFunct->SetNpx(300);


  //exp test
  double integral = 0;
  double avgIntegral = 0;
  double integralError = 0;
  vector<double> integralList = bootstrapIntegralSingleParameter(integral,avgIntegral, integralError, expFunct, NULL, 1.0, 0, 20000, 5.0);

  double smallest = *min_element(integralList.begin(), integralList.end());
  double biggest = *max_element(integralList.begin(), integralList.end());
  TH1D* histo = new TH1D("paramHisto","Sampled Integrals; Counts;Integral Value",100,0.98*smallest, 1.02*biggest);
  for(int iii = 0; iii < (int) integralList.size(); iii++){
    histo->Fill(integralList[iii]);
  }

  TCanvas* canv = new TCanvas("c1");
  histo->Sumw2();
  histo->Draw("E");

  return;
}



