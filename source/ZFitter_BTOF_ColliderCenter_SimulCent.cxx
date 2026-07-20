// THIS IS INSERTED INTO ZFITTER.cxx // Just separate file becuase it's long
//============================================================================================================
//#########      BTOF FITTING   #####################
//============================================================================================================

#include "../headers/ZFitter.h"

// only one of these can be defined at a time
//#define _GAUSSIAN_BTOF_SIMUL_CENT_
//#define _SKEW_NORMAL_BTOF_SIMUL_CENT_
#define _STUDENT_T_BTOF_SIMUL_CENT_   // needs _BOOST_ defined in makefile_toggles.h boost library for integral mode
//#define _PEARSON_IV_BTOF_SIMUL_CENT_ // no integral support yet
//#define _DOUBLE_GAUS_BTOF_SIMUL_CENT_

//#define _CREATE_TEXT_FIT_FILE_

#define _ZFITTER_BTOF_COLLIDERCENTER_INTEGRAL_OPTION_

// extrapolate_BTOF_Fit_Parameters is copied verbatim from the original repo's
// source/ZFitter_BTOF_ColliderCenter.cxx -- a file this port otherwise excludes
// because its own fitBTOFColliderCenter driver is superseded by fitBTOF_SimulCent
// (defined below in this file). fitBTOF_SimulCent itself calls
// extrapolate_BTOF_Fit_Parameters directly (multiple live, uncommented call sites),
// so it's a genuine live dependency that had to be pulled out of the otherwise-
// excluded file rather than left behind. Verified via full-repo grep that this is the
// only active (non-commented) definition of this function anywhere in the original
// repo.
bool ZFitter::extrapolate_BTOF_Fit_Parameters(int a_particleIndexToExtrapolate, double a_lowMomFitRange, double a_highMomFitRange, int a_functionVersion, int a_centIndex, bool a_meanIsInvBeta, int a_specificParameterIndex, int a_convertInvBetaToM2){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to bool ZFitter::extrapolate_BTOF_Fit_Parameters(int " << a_particleIndexToExtrapolate << ", double " << a_lowMomFitRange << ", double " << a_highMomFitRange << ", int " << a_functionVersion << ", int " << a_centIndex << ", bool " << a_meanIsInvBeta << ", int " << a_specificParameterIndex << ")" << endl;
  #endif
  //####################################################################################################################################################################################################################################
  //###################    Extrapolate Gaus / Student-T / SkewNormal Parameters [lambda / sigma / m / nu]                                   ###################################################################################################################################
  //####################################################################################################################################################################################################################################
  bool   success  = true;
  double lowMtM0  = PhysMath::mTm0_from_y_p(m_currentRapidity,a_lowMomFitRange,m_partInfo->GetParticleMass(m_currentPartIndex));
  double highMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity,a_highMomFitRange,m_partInfo->GetParticleMass(m_currentPartIndex));
  if(lowMtM0 < 0.0){
    cout << " WARNING: the rapidity is too extreme for the extrapolation procedure you have given setting low mTm0 to 0" << endl;
    lowMtM0 = 0.0;
  }
  if(highMtM0 < 0.0){
    cout << " ERROR: the rapidity is too extreme for the extrapolation procedure, becuase high mTm0 is negative... (" << highMtM0 << " ) so continuing on without extrapolation" << endl;
    return false;
  }
  double midMtM0 = 0.5*(lowMtM0 + highMtM0);
  #ifdef _ZFITTER_DEBUG_
    cout << "Extrapolating " << m_partInfo->GetParticleName(a_particleIndexToExtrapolate) << " in " << m_partInfo->GetParticleName(m_currentPartIndex)  << "_space with mTm0 fit range = [ " << lowMtM0 << ",  " << highMtM0 << "]" <<  endl;
    cout << "      momentum fit range = [ " << a_lowMomFitRange << ", " << a_highMomFitRange << "]" << endl;
  #endif

  //make the histograms to extrapolate
  TH1D*          param1DHistos[5]   = {nullptr,nullptr,nullptr,nullptr,nullptr};
  TGraphErrors*  param1DGraphs[5]   = {nullptr,nullptr,nullptr,nullptr,nullptr};
  TF1*           extrapFunctions[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  double numberOfParams = 5;
  string distributionName = "Gaussian";
  string paramNames[5]    = {"N","mu","sigma","ignore","ignore"};
  string paramSymbols[5]  = {"N","#mu","#sigma","ignore","ignore"};
  if(a_functionVersion == 1) {paramNames[3] = "gamma"; paramSymbols[3] = "#gamma"; distributionName = "SkewNormal";};
  if(a_functionVersion == 2) {paramNames[3] = "nu";    paramSymbols[3] = "#nu";    distributionName = "StudentT";};
  if(a_functionVersion == 3) {paramNames[3] = "m";     paramSymbols[3] = "m";      paramNames[4] = "nu";        paramSymbols[4] = "#nu";        distributionName = "PearsonIV";};
  if(a_functionVersion == 4) {paramNames[3] = "r";     paramSymbols[3] = "r";      paramNames[4] = "sigma2";    paramSymbols[4] = "#sigma_{2}"; distributionName = "DoubleGauss";};
  string functionNames[5];
  int baseIndex = 1;
  if(a_specificParameterIndex > 0) baseIndex = a_specificParameterIndex;
  if(a_centIndex < 0){
    for(int jjj = 1; jjj < numberOfParams; jjj++){
      if(a_specificParameterIndex >= 0 && jjj != a_specificParameterIndex) continue;
      param1DHistos[jjj] = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][jjj][m_currentPlusMinusIndex]->ProjectionY(Form("param_proj_%s_%s_%s_Rap%02d",
        paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),m_currentRapBin),
        m_currentRapBin,m_currentRapBin);
      functionNames[jjj] = Form("bTOF_param_parab_%s_%s_%s_Rap%02d",paramNames[jjj].c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
        m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),m_currentRapBin);
      param1DHistos[jjj]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),paramSymbols[jjj].c_str()));
      HistogramUtilities::discardNanHistoData(param1DHistos[jjj]);
      param1DGraphs[jjj] = HistogramUtilities::makeGraph(param1DHistos[jjj],false,-0.1, 8.0, true);
    }
  }else{
    for(int jjj = 1; jjj < numberOfParams; jjj++){
      if(a_specificParameterIndex >= 0 && jjj != a_specificParameterIndex) continue;
      param1DHistos[jjj] = m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][a_centIndex][jjj][m_currentPlusMinusIndex]->ProjectionY(Form("param_proj_%s_%s_%s_Cent%02d_Rap%02d",
         paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),a_centIndex,m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
      functionNames[jjj] = Form("bTOF_param_parab_%s_%s_%s_Cent%02d_Rap%02d",paramNames[jjj].c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
        m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),a_centIndex,m_currentRapBin);
      param1DHistos[jjj]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),paramSymbols[jjj].c_str()));
      HistogramUtilities::discardNanHistoData(param1DHistos[jjj]);
      param1DGraphs[jjj] = HistogramUtilities::makeGraph(param1DHistos[jjj],false,-0.1, 8.0, true);
    }
  }

  if(a_specificParameterIndex >= 0 && HistogramUtilities::isEmpty(param1DHistos[a_specificParameterIndex])){
    cout << "ERROR: extrapolate_BTOF_Fit_Parameters. The histogram you are trying to fit is empty!" << endl;
    delete param1DHistos[a_specificParameterIndex];
    delete param1DGraphs[a_specificParameterIndex];
    delete extrapFunctions[a_specificParameterIndex];
    return false;
  }

  #ifdef _ZFITTER_DEBUG_
    if(a_specificParameterIndex >= 0) HistogramUtilities::printHistoInfo(param1DHistos[a_specificParameterIndex]);
  #endif



  if(a_meanIsInvBeta && a_convertInvBetaToM2){ // convert invBeta to m2 becasue you fit with a line
    if(param1DHistos[1]){
      for(int binX = 1; binX <= param1DHistos[1]->GetNbinsX(); binX++){
        double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,param1DHistos[1]->GetBinCenter(binX),m_partInfo->GetParticleMass(m_currentPartIndex));
        double m2       = PhysMath::massSquared_signed_fromInvBeta(momentum, param1DHistos[1]->GetBinContent(binX));
        double m2_Error = PhysMath::massSquaredWidth_from_invBetaWidth(momentum, param1DHistos[1]->GetBinContent(binX), param1DHistos[1]->GetBinError(binX));
        HistogramUtilities::setBinValues(param1DHistos[1],binX, m2, m2_Error);
      }
      string variableSymbol = paramSymbols[1] + "_{m^{2}}";
      param1DHistos[1]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),variableSymbol.c_str()));
    }
    for(int point = 0; param1DGraphs[1] && point < param1DGraphs[1]->GetN(); point++){
      double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,param1DGraphs[1]->GetX()[point],m_partInfo->GetParticleMass(m_currentPartIndex));
      double m2       = PhysMath::massSquared_signed_fromInvBeta(momentum, param1DGraphs[1]->GetY()[point]);
      double m2_Error = PhysMath::massSquaredWidth_from_invBetaWidth(momentum, param1DGraphs[1]->GetY()[point], param1DGraphs[1]->GetEY()[point]);
      param1DGraphs[1]->SetPoint(point,param1DGraphs[1]->GetX()[point],m2);
      param1DGraphs[1]->SetPointError(point,param1DGraphs[1]->GetEX()[point],m2_Error);
    }
  }

  //Make the graphs that ignore zero bins



  int lowBin  = param1DHistos[baseIndex]->FindBin(lowMtM0);
  int midBin  = param1DHistos[baseIndex]->FindBin(midMtM0);
  int highBin = param1DHistos[baseIndex]->FindBin(highMtM0);

  //  Find bins with actual fit content
  bool properLowBinFound = (param1DHistos[baseIndex]->GetBinError(lowBin) > 0.0);
  bool properHighBinFound = (param1DHistos[baseIndex]->GetBinError(highBin) > 0.0);
  if(!properLowBinFound){
    for(int binX = lowBin; binX < highBin; binX++){
      if(param1DHistos[baseIndex]->GetBinError(binX) > 0.0){
        properLowBinFound = true;
        lowBin = binX;
        break;
      }
    }
  }
  if(!properHighBinFound){
    properHighBinFound = false;
    for(int binX = highBin; binX > lowBin; binX--){
      if(param1DHistos[baseIndex]->GetBinError(binX) > 0.0){
        properHighBinFound = true;
        highBin = binX;
        break;
      }
    }
  }

  if(!properLowBinFound || !properHighBinFound){
    cout << " ERROR: The Extrapolation will not work. Could not find bins with actual fit data..." << endl;
    success = false;
  }
  if(midBin <= lowBin || midBin >=highBin){
    midBin = (lowBin+highBin)/2;
    cout << " WARNING: Mid-Bin Didn't work, so putting it to half the range" << endl;
  }
  if(highBin - lowBin < 2){
    cout << " ERROR: The extrapolation will fail because you only have two bins to extrapolate..." << endl;
    success = false;
  }
  if(!success){
    delete param1DHistos[a_specificParameterIndex];
    delete param1DGraphs[a_specificParameterIndex];
    delete extrapFunctions[a_specificParameterIndex];
    return false;
  }

  lowMtM0  = param1DHistos[baseIndex]->GetBinCenter(lowBin)  - 0.25*param1DHistos[baseIndex]->GetBinWidth(lowBin);
  highMtM0 = param1DHistos[baseIndex]->GetBinCenter(highBin) + 0.25*param1DHistos[baseIndex]->GetBinWidth(highBin);
  midMtM0  = 0.5*(lowMtM0 + highMtM0);
  midBin   = param1DHistos[baseIndex]->FindBin(midMtM0);
  midMtM0  = param1DHistos[baseIndex]->GetBinCenter(midBin);

  #ifdef _ZFITTER_DEBUG_
    cout << " Resized the fit range... lowMtM0: " << lowMtM0 << " midMtM0: " << midMtM0 << " highMtM0: " << highMtM0 << endl;
    cout <<  "                   LowBin: " << lowBin << " midBin: " << midBin << " highBin: " << highBin << endl;
    cout << "    Ptrs: " << param1DHistos[0] << "  " << param1DHistos[1] << "  " << param1DHistos[2] << "  " << param1DHistos[3] << "  " << param1DHistos[4] << endl;
  #endif

  bool fitSuccess = false;
  // ########    PARAMETER LOOP FOR MEAN WIDTH AND MAYBE NDF ###########
  for(int jjj = 1; jjj < numberOfParams; jjj++){
    if(!param1DHistos[jjj]) continue;
    double amplitudes[3];
    amplitudes[0] = param1DHistos[jjj]->GetBinContent(lowBin);
    amplitudes[1] = param1DHistos[jjj]->GetBinContent(midBin);
    amplitudes[2] = param1DHistos[jjj]->GetBinContent(highBin);
    if(amplitudes[0] == 0.0 && amplitudes[1] == 0.0 && amplitudes[2] == 0.0){
      #ifdef _ZFITTER_DEBUG_
        cout << " ERROR: all the amplitudes at the three points are zero... continuing" << endl;
      #endif
      continue;
    }


    if(a_meanIsInvBeta && !a_convertInvBetaToM2 && jjj == 1){
      extrapFunctions[jjj] = PhysMath::getInvBetaCalibrationFunction_byMtM0(functionNames[jjj].c_str(), m_partInfo->GetParticleMass(m_currentPartIndex),
                                m_partInfo->GetParticleMass(a_particleIndexToExtrapolate), m_currentRapidity, 0.0);
    }else{
      double slope_1  = (amplitudes[1] - amplitudes[0])/(midMtM0 - lowMtM0);
      double slope_2  = (amplitudes[2] - amplitudes[1])/(highMtM0 - midMtM0);
      double ySlope_0 = 0.5* (slope_1 + slope_2);
      double maxVal   = HistogramUtilities::getMaxInRange(param1DHistos[jjj], lowMtM0, highMtM0);
      double minVal   = HistogramUtilities::getMinInRange(param1DHistos[jjj], lowMtM0, highMtM0);
      double distanceScale = maxVal - minVal;
      double slopeScale = distanceScale / (highMtM0 - lowMtM0);

      #ifdef _ZFITTER_DEBUG_
        cout << " Initializing linear funciton with y0=" << amplitudes[0] << " y'= " << ySlope_0 << endl;
      #endif

      extrapFunctions[jjj] = PhysMath::linearFunct(functionNames[jjj].c_str(), amplitudes[0], ySlope_0,lowMtM0);
      if(a_functionVersion == 3 && jjj > 2){
        extrapFunctions[jjj]->SetParLimits(0,amplitudes[0] - 10.0*distanceScale,amplitudes[0] + 10.0*distanceScale);
        extrapFunctions[jjj]->FixParameter(1,0.0);
      }else{
        //extrapFunctions[jjj]->SetParLimits(0,amplitudes[0] - 10.0*distanceScale,amplitudes[0] + 10.0*distanceScale);
        extrapFunctions[jjj]->FixParameter(0,amplitudes[0]);
        extrapFunctions[jjj]->SetParLimits(1,-10.0*slopeScale,10.0*slopeScale);
      }
    }
    extrapFunctions[jjj]->SetRange(lowMtM0,highMtM0);


    TF1* initialFunct = (TF1*) extrapFunctions[jjj]->Clone();
    initialFunct->SetLineColor(kYellow);
    initialFunct->SetName("initialFunct");

    string fitOptions = "";
    #ifdef _ZFITTER_DEBUG_
      fitOptions += "";
    #else
      fitOptions += " Q";
    #endif
    int  status    = param1DGraphs[jjj]->Fit(functionNames[jjj].c_str(),fitOptions.c_str(),"",lowMtM0,highMtM0);
    TF1* fitResult = param1DGraphs[jjj]->GetFunction(functionNames[jjj].c_str());

    if(jjj == 3){
      if(!fitResult) fitResult = (TF1*) initialFunct->Clone();
      fitResult->SetRange(highMtM0,10.0);
      fitResult->SetParameter(0,amplitudes[2]);
      fitResult->SetParameter(1,0.0);
      fitResult->SetParameter(2,highMtM0);
    }

    #ifdef _ZFITTER_DEBUG_
      cout << "Fit Completed with status = " << status << " and fitResult: " << fitResult << endl;
    #endif

    //Push the results to the histograms
    //double mu = -999.999;
    if(fitResult){
      fitSuccess = true;
      for(int binX = highBin+1; binX <= param1DHistos[jjj]->GetNbinsX(); binX++){
        double center       = param1DHistos[jjj]->GetBinCenter(binX);
        double valueToPush  = fitResult->Eval(center);
        double errorToPush  = 1.0E-13;
        if(jjj != 1) errorToPush = sqrt( pow(fitResult->GetParError(1)*center,2) + pow(fitResult->GetParameter(1)*fitResult->GetParError(2),2) + pow(fitResult->GetParError(0),2));
        if(errorToPush <= 0.0) errorToPush = 1.0E-13;
        if(a_meanIsInvBeta && a_convertInvBetaToM2 && jjj == 1){ // convert m2 back to invBeta
          double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,param1DHistos[1]->GetBinCenter(binX),m_partInfo->GetParticleMass(m_currentPartIndex));
          errorToPush = PhysMath::invBetaWidth_from_MassSquaredWidth(valueToPush, errorToPush, momentum);
          valueToPush  = PhysMath::invBeta_from_MassSquaredSigned(valueToPush, momentum);
        }
        if(a_centIndex < 0){
          HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][jjj][m_currentPlusMinusIndex],m_currentRapBin, binX,valueToPush,errorToPush);
        }else{
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][a_centIndex][jjj][m_currentPlusMinusIndex],m_currentRapBin, binX,valueToPush,errorToPush);
        }
      }
    }else{
      cout << " ERROR: " << m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data() << " PARAMETER " << jjj << " EXTRAPOLATION FAILED!!!!!!!!! " << endl;
      #ifdef _ZFITTER_DEBUG_
        HistogramUtilities::printHistoInfo(param1DHistos[jjj]);
      #endif
    }


    TGraph* parabPoint = nullptr;
    if(jjj != 1){
      parabPoint = new TGraph();
      if(fitResult) parabPoint->SetPoint(0,fitResult->GetParameter(2),fitResult->GetParameter(0));
      parabPoint->SetMarkerStyle(20);
      parabPoint->SetMarkerColor(kRed);
    }


    // DRAW THE FITS
    TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
    fittingCanvas->SetWindowSize(1300,700);
    fittingCanvas->SetCanvasSize(1200,600);
    fittingCanvas->cd();
    gStyle->SetOptFit(0111);
    double lowDrawRange = lowMtM0 - 3.0*(highMtM0 - lowMtM0);
    if(lowDrawRange < 0.0) lowDrawRange = 0.0;
    double highDrawRange =  highMtM0 + 2.0*(highMtM0 - lowMtM0);
    //HistogramUtilities::setAxisRanges(param1DHistos[jjj], lowDrawRange,highDrawRange,false);
    TH1D* frame =  HistogramUtilities::generateGraphFrame(param1DGraphs[jjj], lowDrawRange, highDrawRange, false, false, param1DHistos[jjj]->GetXaxis()->GetTitle(), param1DHistos[jjj]->GetYaxis()->GetTitle());


    /*TLine* extraLine = NULL;
    if(mu != 999.999 && fitResult && mu < highMtM0 && mu > lowMtM0){
      extraLine = new TLine(mu,fitResult->Eval(mu), highMtM0, fitResult->Eval(mu));
      extraLine->SetLineColor(kRed);
    }*/
    //param1DHistos[jjj]->SetStats(true);
    //param1DHistos[jjj]->Draw("E");
    frame->Draw("E");
    #ifndef _BEFORE_ROOT_6p24_
    param1DGraphs[jjj]->SetStats(true);
    #endif
    param1DGraphs[jjj]->Draw("P");
    initialFunct->Draw("SAME");
    if(parabPoint) parabPoint->Draw("P");
    //if(extraLine) extraLine->Draw("same");
    //gPad->SetLogy();
    gPad->SetRightMargin(0.4);
    fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Rap_%02d_%s_Param%02d%s.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin,
      m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),jjj,(a_centIndex < 0 ? "" : Form("_Cent%02d",a_centIndex))));
    delete frame;

    TH1D* frame2 =  HistogramUtilities::generateGraphFrame(param1DGraphs[jjj], 0, highDrawRange, false, false, param1DHistos[jjj]->GetXaxis()->GetTitle(), param1DHistos[jjj]->GetYaxis()->GetTitle());
    frame2->Draw("E");
    #ifndef _BEFORE_ROOT_6p24_
    param1DGraphs[jjj]->SetStats(true);
    #endif
    param1DGraphs[jjj]->Draw("P");
    initialFunct->Draw("SAME");
    if(parabPoint) parabPoint->Draw("P");
    fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Rap_%02d_%s_Param%02d%s_FullRange.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin,
      m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),jjj,(a_centIndex < 0 ? "" : Form("_Cent%02d",a_centIndex))));

    //cout << " param1DHistos[jjj] " << param1DHistos[jjj] << endl;
    delete param1DHistos[jjj];
    delete param1DGraphs[jjj];
    //cout << " initialFunct " << initialFunct << endl;
    delete initialFunct;
    //cout << " extrapFunctions[jjj] " << extrapFunctions[jjj] << endl;
    delete extrapFunctions[jjj];
    //cout << " parabPoint " << parabPoint << endl;
    if(parabPoint) delete parabPoint;
    //cout << " fittingCanvas " << fittingCanvas << endl;
    delete frame2;
    delete fittingCanvas;
    //cout << " fitResult " << fitResult << endl; if(fitResult) delete fitResult;

  } // jjj parameter loop


  #ifdef _ZFITTER_DEBUG_
    cout << "End call to bool ZFitter::extrapolate_BTOF_Fit_Parameters(int " << a_particleIndexToExtrapolate << ", double " << a_lowMomFitRange << ", double " << a_highMomFitRange << ", int " << a_functionVersion << ", int " << a_centIndex << ", bool " << a_meanIsInvBeta << ", int " << a_specificParameterIndex << ")" << endl;
  #endif
  return fitSuccess;
}

double ZFitter::simultaneous_centrality_ZTOF_chisqr(const double* a_param){
  //mu_pi sigma_pi nu_pi p_pi    mu_K sigma_K nu_K p_K   mu_pro sigma_pro nu_pro p_pro  then For Each Cent: A_pi A_K A_pro ... A_pi A_K A_pro
  //   0    1        2     3       4     5      6    7     8       9         10    11
  double chisqr = 0.0;
  #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
    double alpha_pi  = a_param[2]/sqrt(1.0 - a_param[2]*a_param[2]);
    double alpha_K   = a_param[6]/sqrt(1.0 - a_param[6]*a_param[6]);
    double alpha_pro = a_param[10]/sqrt(1.0 - a_param[10]*a_param[10]);
  #endif
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    for(int binX = m_currentLowFitBin[centIndex]; binX <= m_currentHighFitBin[centIndex]; binX++){
      double content = m_currentHistosToFit_ByCent_InvBeta[centIndex]->GetBinContent(binX);
      #ifdef _ZFITTER_BTOF_COLLIDERCENTER_INTEGRAL_OPTION_
        double lowEdge  = m_currentHistosToFit_ByCent_InvBeta[centIndex]->GetBinLowEdge(binX);
        double highEdge = m_currentHistosToFit_ByCent_InvBeta[centIndex]->GetBinLowEdge(binX+1);
        if(content > 0.0){
          #ifdef _GAUSSIAN_BTOF_SIMUL_CENT_
            double prediction = PhysMath::gaus_finite_integral(a_param[12+(3*centIndex)],a_param[0],a_param[1], lowEdge, highEdge)
                              + PhysMath::gaus_finite_integral(a_param[13+(3*centIndex)],a_param[4],a_param[5], lowEdge, highEdge)
                              + PhysMath::gaus_finite_integral(a_param[14+(3*centIndex)],a_param[8],a_param[9], lowEdge, highEdge);
          #endif
          #ifdef _SKEW_NORMAL_BTOF_SIMUL_CENT_
            double prediction = PhysMath::skewNormal_finite_integral(a_param[12 + (3*centIndex)],a_param[0],a_param[1], alpha_pi, lowEdge,highEdge)
                              + PhysMath::skewNormal_finite_integral(a_param[13 + (3*centIndex)],a_param[4],a_param[5], alpha_K,  lowEdge,highEdge)
                              + PhysMath::skewNormal_finite_integral(a_param[14 + (3*centIndex)],a_param[8],a_param[9], alpha_K,  lowEdge,highEdge);
          #endif

          #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
            //double prediction = PhysMath::student_T_finite_integral(a_param[12 + (3*centIndex)],a_param[0],a_param[1],a_param[2] ,lowEdge,highEdge)
            //                  + PhysMath::student_T_finite_integral(a_param[13 + (3*centIndex)],a_param[4],a_param[5],a_param[6] ,lowEdge,highEdge)
            //                  + PhysMath::student_T_finite_integral(a_param[14 + (3*centIndex)],a_param[8],a_param[9],a_param[10] ,lowEdge,highEdge);
            double prediction = a_param[12 + (3*centIndex)]*PhysMath::student_T_BinContent(a_param[0],a_param[1],a_param[2] ,lowEdge,highEdge)
                              + a_param[13 + (3*centIndex)]*PhysMath::student_T_BinContent(a_param[4],a_param[5],a_param[6] ,lowEdge,highEdge)
                              + a_param[14 + (3*centIndex)]*PhysMath::student_T_BinContent(a_param[8],a_param[9],a_param[10] ,lowEdge,highEdge);
            //double prediction = a_param[12 + (3*centIndex)]*PhysMath::m2_student_T_BinContent(a_param[0],a_param[1],a_param[2] ,lowEdge,highEdge, m_currentMomOfBin)
            //                  + a_param[13 + (3*centIndex)]*PhysMath::m2_student_T_BinContent(a_param[4],a_param[5],a_param[6] ,lowEdge,highEdge, m_currentMomOfBin)
            //                  + a_param[14 + (3*centIndex)]*PhysMath::m2_student_T_BinContent(a_param[8],a_param[9],a_param[10] ,lowEdge,highEdge, m_currentMomOfBin);

          #endif

          #ifdef _PEARSON_IV_BTOF_SIMUL_CENT_
            double prediction =   a_param[12 + (3*centIndex)]*PhysMath::pearsonIV(a_param[0],a_param[1],a_param[2],a_param[3],center,true)
                                + a_param[13 + (3*centIndex)]*PhysMath::pearsonIV(a_param[4],a_param[5],a_param[6],a_param[7],center,true)
                                + a_param[14 + (3*centIndex)]*PhysMath::pearsonIV(a_param[8],a_param[9],a_param[10],a_param[11],center,true);
          #endif

          #ifdef _DOUBLE_GAUS_BTOF_SIMUL_CENT_
            double prediction = PhysMath::doubleGauss(a_param[12 + (3*centIndex)],a_param[0],a_param[1],a_param[2],a_param[3], center)
                              + PhysMath::doubleGauss(a_param[13 + (3*centIndex)],a_param[4],a_param[5],a_param[6],a_param[7], center)
                              + PhysMath::doubleGauss(a_param[14 + (3*centIndex)],a_param[8],a_param[9],a_param[10],a_param[11], center);
          #endif

          chisqr += pow((prediction - content)/m_currentHistosToFit_ByCent_InvBeta[centIndex]->GetBinError(binX),2);
        }
      #else
        double center  = m_currentHistosToFit_ByCent_InvBeta[centIndex]->GetBinCenter(binX);
        if(content > 0.0){
          #ifdef _GAUSSIAN_BTOF_SIMUL_CENT_
            double prediction = a_param[12 + (3*centIndex)]*PhysMath::gaus(a_param[0],a_param[1],center)
                              + a_param[13 + (3*centIndex)]*PhysMath::gaus(a_param[4],a_param[5],center)
                              + a_param[14 + (3*centIndex)]*PhysMath::gaus(a_param[8],a_param[9],center);
          #endif
          #ifdef _SKEW_NORMAL_BTOF_SIMUL_CENT_
            double prediction =   a_param[12 + (3*centIndex)]*PhysMath::skewNormal(a_param[0],a_param[1],a_param[2] ,center)
                                + a_param[13 + (3*centIndex)]*PhysMath::skewNormal(a_param[4],a_param[5],a_param[6] ,center)
                                + a_param[14 + (3*centIndex)]*PhysMath::skewNormal(a_param[8],a_param[9],a_param[10],center);
          #endif
          #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
            double prediction =   a_param[12 + (3*centIndex)]*PhysMath::student_T(a_param[0],a_param[1],a_param[2] ,center)
                                + a_param[13 + (3*centIndex)]*PhysMath::student_T(a_param[4],a_param[5],a_param[6] ,center)
                                + a_param[14 + (3*centIndex)]*PhysMath::student_T(a_param[8],a_param[9],a_param[10],center);
            //double prediction = a_param[12 + (3*centIndex)]*PhysMath::m2_student_T(a_param[0],a_param[1],a_param[2] ,center, m_currentMomOfBin)
            //                  + a_param[13 + (3*centIndex)]*PhysMath::m2_student_T(a_param[4],a_param[5],a_param[6] ,center, m_currentMomOfBin)
            //                  + a_param[14 + (3*centIndex)]*PhysMath::m2_student_T(a_param[8],a_param[9],a_param[10],center, m_currentMomOfBin);
          #endif

          #ifdef _PEARSON_IV_BTOF_SIMUL_CENT_
            double prediction =   a_param[12 + (3*centIndex)]*PhysMath::pearsonIV(a_param[0],a_param[1],a_param[2],a_param[3],center,true)
                                + a_param[13 + (3*centIndex)]*PhysMath::pearsonIV(a_param[4],a_param[5],a_param[6],a_param[7],center,true)
                                + a_param[14 + (3*centIndex)]*PhysMath::pearsonIV(a_param[8],a_param[9],a_param[10],a_param[11],center,true);
          #endif
          #ifdef _DOUBLE_GAUS_BTOF_SIMUL_CENT_
            double prediction = PhysMath::doubleGauss(a_param[12 + (3*centIndex)],a_param[0],a_param[1],a_param[2],a_param[3], center)
                              + PhysMath::doubleGauss(a_param[13 + (3*centIndex)],a_param[4],a_param[5],a_param[6],a_param[7], center)
                              + PhysMath::doubleGauss(a_param[14 + (3*centIndex)],a_param[8],a_param[9],a_param[10],a_param[11], center);
          #endif
          chisqr += pow((prediction - content)/m_currentHistosToFit_ByCent_InvBeta[centIndex]->GetBinError(binX),2);
        }
      #endif
    }
  }
  //#ifdef _ZFITTER_DEBUG_
  //  cout << "simultaneous_centrality_ZTOF_chisqr: " << chisqr << endl;
  //#endif
  return chisqr;
}





void ZFitter::fitBTOF_SimulCent(int a_partIndexSpace, int a_charge){
  #ifdef _ZFITTER_DEBUG_
    cout << " Call to ZFitter::fitBTOF_SimulCent(int" <<  a_partIndexSpace << ", int " << a_charge << " )" << endl;
  #endif

  TVirtualFitter::SetMaxIterations(10000000);
  bool usePreviousGoodParams_Step03 = true;
  double lowPearsonIV_m_bound = 0.001;


  double momToStartPreviousGoodParameters_ForKaon = 0.3;
  double momToStartPreviousGoodParameters_ForProton = 0.4;
  double momToIntroduceProton     = 1.1; // only for pi/K
  double pionLowMomParamExtrap    = 1.0; // fits from this to the next with lines
  double pionHighMomParamExtrap   = 1.65; // fits to this with lines
  double kaonLowMomParamExtrap    = 1.0; // fits from this to the next with lines
  double kaonHighMomParamExtrap   = 1.65; // fits to this with lines
  double kaonMinusLowMomParamExtrap    = 1.0; // fits from this to the next with lines
  double kaonMinusHighMomParamExtrap   = 1.65; // fits to this with lines
  double protonLowMomParamExtrap  = 1.5; // fits from this to the next with lines
  double protonHighMomParamExtrap = 2.4; // fits to this with lines
  double antiprotonLowMomParamExtrap  = 1.5; // fits from this to the next with lines
  double antiprotonHighMomParamExtrap = 2.4; // fits to this with lines
  double momToDoPionKaonMatrixAmpInit = 1.9;
  double momToFixKaonStudentTNu       = 1.7;
  //double rapidityToNotDoParameterExtrapolation = 0.75;
  double minStudentTNuVal = 1.001; // integral doesn't matter if not particle of interest
  double minStudentTNuVal_POI = 3.0; // integral diverges as nu -> 1.0
  bool   normalizePearsonIV = true; // add toggles to simul chi^2 if you change this
  bool   usePositiveChargeAsSeeds = true;
  bool   fixPOI_ToPositiveValues = (a_partIndexSpace == 2 && a_charge == -1);
  bool   doHesse = false;
  double invBetaResolution = 0.02;// 0.004 //0.00843;
  bool   drawInvBetaFits = false; // doesn't affect mass^2 drawing

  //### PION SPACE FIT RANGE INFO ###
  double upToMomToIgnoreKaons = 0.5;
  double highFitRangeWhenIngoringKaonsInPion = 0.12;
  double upToMomToIgnoreMuons = 0.4;
  double highMomForMuonIgnore = 0.85; // the smooth transition to the low cut
  double m2ForMuonCut = 0.9*pow(m_partInfo->GetParticleMass(8),2.0) + 0.1*pow(m_partInfo->GetParticleMass(0),2.0);
  double basicPionM2LowCut = -0.05;
  TF1* lowerPionFitRangeForPionLinFunct = PhysMath::linearFunct("lowPiForPiFunctRange", upToMomToIgnoreMuons, m2ForMuonCut, highMomForMuonIgnore, basicPionM2LowCut);

  //### KAON SPACE FIT RANGE INFO ###
  //double lowFitRangeForPionsInKaonSpace = pow(m_partInfo->GetParticleMass(0),2.0) + 0.25*(pow(m_partInfo->GetParticleMass(1),2.0) - pow(m_partInfo->GetParticleMass(0),2.0));
  double lowFitRangeForPionsInKaonSpace = pow(m_partInfo->GetParticleMass(0),2.0); // right pi tail is longer than left
  double reducedIsolatedKaon_lowFitRange  = 0.21;
  double reducedIsolatedKaon_highFitRange = 0.3;
  double highMomToIgnorePionsInKaonSpace  = 0.55;
  double lowerPionRangeMom = 0.85; // for kaon space, how much of left pion to include
  double upperPionRangeMom = 1.25;
  double highFitRangeForProtonsInKaonSpace = pow(m_partInfo->GetParticleMass(2),2.0);
  TF1* lowerPionFitRangeForKaonLinFunct = PhysMath::linearFunct("lowPiForKaonFunctRange", lowerPionRangeMom, lowFitRangeForPionsInKaonSpace, upperPionRangeMom, -0.1);

  double lowerProtonRangeMom = 2.0; // for kaon space, how much of left Proton to include
  double upperProtonRangeMom = 2.7;

  TF1* upperProtonFitRangeForKaonLinFunct = PhysMath::linearFunct("highProtonFunctRange", lowerProtonRangeMom, highFitRangeForProtonsInKaonSpace, upperProtonRangeMom, 1.1);


  //### PROTON SPACE FIT RANGE INFO ###
  double upToMomToIgnorePionKaonsForProtons = 0.9;
  double lowFitRangeForKaonsInProtonSpace   = 0.1;
  double lowFitRangeForPionsInProtonSpace   = -0.3;
  //use same low limits as Kaons


  //mass^2 width
  double mom_pion_array[5]       = {0.23,0.38,0.82,1.18,1.49};
  double m2width_pion_array[5]   = {0.0051,0.011,0.041,0.065,0.085};
  double mom_kaon_array[5]       = {0.42,0.67,0.97,1.42,2.0};
  double m2width_kaon_array[5]   = {0.0085,0.013,0.025,0.0425,0.0747};
  double mom_proton_array[8]     = {0.35,0.67,1.13,1.58,2.07,2.59,3.27,3.8};
  double m2width_proton_array[8] = {0.207,0.01833,0.027,0.041,0.061,0.0846,0.115,0.1404};

  TGraph* m2width_pion_graph   = new TGraph(5,mom_pion_array,  m2width_pion_array);
  TGraph* m2width_kaon_graph   = new TGraph(5,mom_kaon_array,  m2width_kaon_array);
  TGraph* m2width_proton_graph = new TGraph(8,mom_proton_array,m2width_proton_array);


  //particle charge rapidity param lowMom highMom        param = -999 means all of them   rapidity = -999 means all of them    param 0 is amplitude, 1 is mean
  // This first rules take priority over the last rules
  vector< vector<double> > extrapolationMatrix_pion = {

    // dont extrapolate amplitudes
    {0, 1, -999, 0, 999, 999},
    {1, 1, -999, 0, 999, 999},
    {2, 1, -999, 0, 999, 999},
    {0,-1, -999, 0, 999, 999},
    {1,-1, -999, 0, 999, 999},
    {2,-1, -999, 0, 999, 999},
    // dont extrapolate sigma
    {0, 1, -999, 2, 999, 999},
    {1, 1, -999, 2, 999, 999},
    {2, 1, -999, 2, 999, 999},
    {0,-1, -999, 2, 999, 999},
    {1,-1, -999, 2, 999, 999},
    {2,-1, -999, 2, 999, 999},
    // dont extrapolate the ignore parameter
    {0, 1, -999, 4, 999, 999},
    {1, 1, -999, 4, 999, 999},
    {2, 1, -999, 4, 999, 999},
    {0,-1, -999, 4, 999, 999},
    {1,-1, -999, 4, 999, 999},
    {2,-1, -999, 4, 999, 999},


    // put general rules here
    {0, 1, -999, -999, 1.3, 1.6},
    {1, 1, -999, -999, 1.3, 1.6},
    {2, 1, -999, -999, 1.7, 2.25},
    {0,-1, -999, -999, 1.3, 1.6},
    {1,-1, -999, -999, 1.3, 1.6},
    {2,-1, -999, -999, 1.7, 2.25}

  };


  vector< vector<double> > extrapolationMatrix_kaon = {

    // dont extrapolate amplitudes
    {0, 1, -999, 0, 999, 999},
    {1, 1, -999, 0, 999, 999},
    {2, 1, -999, 0, 999, 999},
    {0,-1, -999, 0, 999, 999},
    {1,-1, -999, 0, 999, 999},
    {2,-1, -999, 0, 999, 999},
    // dont extrapolate sigma
    {0, 1, -999, 2, 999, 999},
    {1, 1, -999, 2, 999, 999},
    {2, 1, -999, 2, 999, 999},
    {0,-1, -999, 2, 999, 999},
    {1,-1, -999, 2, 999, 999},
    {2,-1, -999, 2, 999, 999},
    // dont extrapolate the ignore parameter
    {0, 1, -999, 4, 999, 999},
    {1, 1, -999, 4, 999, 999},
    {2, 1, -999, 4, 999, 999},
    {0,-1, -999, 4, 999, 999},
    {1,-1, -999, 4, 999, 999},
    {2,-1, -999, 4, 999, 999},


    // put general rules here
    {0, 1, -999, -999, 1.3, 1.6},
    {1, 1, -999, -999, 1.3, 1.6},
    {2, 1, -999, -999, 1.7, 2.25},
    {0,-1, -999, -999, 1.3, 1.6},
    {1,-1, -999, -999, 1.3, 1.6},
    {2,-1, -999, -999, 1.7, 2.25}

  };



  vector< vector<double> > extrapolationMatrix_proton = {


    // -------------  specific rules  ---------------------
    // Pion Mean ------------------------------
    //{0, 1, 0.0, 1,  ,  },
    //{0, 1, 0.1, 1,  ,  },
    //{0, 1, 0.2, 1,  ,  },
    //{0, 1, 0.3, 1,  ,  },
    //{0, 1, 0.4, 1,  ,  },
    //{0, 1, 0.5, 1,  ,  },
    //{0, 1, 0.6, 1,  ,  },
    {0, 1, 0.7, 1, 1.5, 2.1},
    {0, 1, 0.8, 1, 1.72, 2.0},
    {0, 1, 0.9, 1, 1.85, 2.2},

    // Pion Sigma ------------------------------
    {0, 1, 0.0, 2, 1.8, 2.1},
    {0, 1, 0.1, 2, 1.8, 2.1},
    {0, 1, 0.2, 2, 1.8, 2.1},
    {0, 1, 0.3, 2, 1.8, 2.1},
    {0, 1, 0.4, 2, 1.8, 2.1},
    {0, 1, 0.5, 2, 1.8, 2.1},
    {0, 1, 0.6, 2, 1.8, 2.1},
    {0, 1, 0.7, 2, 1.8, 2.1},
    {0, 1, 0.8, 2, 2.0, 2.25},
    {0, 1, 0.9, 2, 1.9, 2.25},

    // Pion Nu ------------------------------
    {0, 1, 0.0, 3, 0, 1.6},
    {0, 1, 0.1, 3, 0, 1.6},
    {0, 1, 0.2, 3, 0, 1.6},
    {0, 1, 0.3, 3, 0, 1.6},
    {0, 1, 0.4, 3, 0, 1.6},
    {0, 1, 0.5, 3, 0, 1.6},
    {0, 1, 0.6, 3, 0, 1.6},
    {0, 1, 0.7, 3, 0, 1.5},
    {0, 1, 0.8, 3, 0, 1.5},
    {0, 1, 0.9, 3, 0, 1.9},

    // Kaon Mean ------------------------------
    //{1, 1, 0.0, 1,  ,  },
    //{1, 1, 0.1, 1,  ,  },
    //{1, 1, 0.2, 1,  ,  },
    //{1, 1, 0.3, 1,  ,  },
    //{1, 1, 0.4, 1,  ,  },
    //{1, 1, 0.5, 1,  ,  },
    //{1, 1, 0.6, 1,  ,  },
    {1, 1, 0.7, 1, 1.5, 2.0},
    {1, 1, 0.8, 1, 1.72, 2.0},
    {1, 1, 0.9, 1, 1.85, 2.0},

    // Kaon Sigma ------------------------------
    {1, 1, 0.0, 2, 1.7, 2.0},
    {1, 1, 0.1, 2, 1.7, 2.0},
    {1, 1, 0.2, 2, 1.7, 2.0},
    {1, 1, 0.3, 2, 1.7, 2.0},
    {1, 1, 0.4, 2, 1.7, 2.0},
    {1, 1, 0.5, 2, 1.7, 2.0},
    {1, 1, 0.6, 2, 1.7, 2.0},
    {1, 1, 0.7, 2, 1.72, 2.0},
    {1, 1, 0.8, 2, 1.72, 2.0},
    {1, 1, 0.9, 2, 1.85, 2.0},

    // Kaon Nu ------------------------------
    {1, 1, 0.0, 3, 0, 1.6},
    {1, 1, 0.1, 3, 0, 1.6},
    {1, 1, 0.2, 3, 0, 1.6},
    {1, 1, 0.3, 3, 0, 1.6},
    {1, 1, 0.4, 3, 0, 1.6},
    {1, 1, 0.5, 3, 0, 1.6},
    {1, 1, 0.6, 3, 0, 1.6},
    {1, 1, 0.7, 3, 0, 1.5},
    {1, 1, 0.8, 3, 0, 1.5},
    {1, 1, 0.9, 3, 0, 1.9},


    // Proton Mean ------------------------------
    //{2, 1, 0.0, 1,  ,  },
    //{2, 1, 0.1, 1,  ,  },
    //{2, 1, 0.2, 1,  ,  },
    //{2, 1, 0.3, 1,  ,  },
    //{2, 1, 0.4, 1,  ,  },
    //{2, 1, 0.5, 1,  ,  },
    //{2, 1, 0.6, 1,  ,  },
    //{2, 1, 0.7, 1, 1.5, 2.0},
    //{2, 1, 0.8, 1, 1.72, 2.0},
    //{2, 1, 0.9, 1, 1.85, 2.0},

    // Proton Sigma ------------------------------
    //{2, 1, 0.0, 2,  1.8, 2.0},
    //{2, 1, 0.1, 2,  1.8, 2.0},
    //{2, 1, 0.2, 2,  1.8, 2.0},
    //{2, 1, 0.3, 2,  1.8, 2.0},
    //{2, 1, 0.4, 2,  1.8, 2.0},
    //{2, 1, 0.5, 2,  1.8, 2.0},
    //{2, 1, 0.6, 2,  1.8, 2.0},
    //{2, 1, 0.7, 2, 1.72, 2.0},
    //{2, 1, 0.8, 2, 1.72, 2.0},
    //{2, 1, 0.9, 2, 1.85, 2.0},

    // Proton Nu ------------------------------
    //{2, 1, 0.0, 3, 0, 1.5},
    //{2, 1, 0.1, 3, 0, 1.5},
    //{2, 1, 0.2, 3, 0, 1.5},
    //{2, 1, 0.3, 3, 0, 1.5},
    //{2, 1, 0.4, 3, 0, 1.5},
    //{2, 1, 0.5, 3, 0, 1.5},
    //{2, 1, 0.6, 3, 0, 1.5},
    //{2, 1, 0.7, 3, 0, 1.5},
    //{2, 1, 0.8, 3, 0, 1.5},
    //{2, 1, 0.9, 3, 0, 1.9},


    //--------------------------   general rules  ----------------------------

    // dont extrapolate amplitudes
    {0, 1, -999, 0, 999, 999},
    {1, 1, -999, 0, 999, 999},
    {2, 1, -999, 0, 999, 999},
    {0,-1, -999, 0, 999, 999},
    {1,-1, -999, 0, 999, 999},
    {2,-1, -999, 0, 999, 999},
    // dont extrapolate sigma
    {0, 1, -999, 2, 999, 999},
    {1, 1, -999, 2, 999, 999},
    {2, 1, -999, 2, 999, 999},
    {0,-1, -999, 2, 999, 999},
    {1,-1, -999, 2, 999, 999},
    {2,-1, -999, 2, 999, 999},
    // dont extrapolate the ignore parameter
    {0, 1, -999, 4, 999, 999},
    {1, 1, -999, 4, 999, 999},
    {2, 1, -999, 4, 999, 999},
    {0,-1, -999, 4, 999, 999},
    {1,-1, -999, 4, 999, 999},
    {2,-1, -999, 4, 999, 999},



    {0, 1, -999, -999, 1.3, 1.6},
    {1, 1, -999, -999, 1.3, 1.6},
    {2, 1, -999, -999, 1.7, 2.25},
    {0,-1, -999, -999, 1.3, 1.6},
    {1,-1, -999, -999, 1.3, 1.6},
    {2,-1, -999, -999, 1.7, 2.25}

  };


  vector< vector<double> > extrapolationMatrix;
  if(m_currentPartIndex == 0) extrapolationMatrix = extrapolationMatrix_pion;
  if(m_currentPartIndex == 1) extrapolationMatrix = extrapolationMatrix_kaon;
  if(m_currentPartIndex == 2) extrapolationMatrix = extrapolationMatrix_proton;


  //particle rapidity param momentum      param = -999 means all of them   rapidity = -999 means all of them     param 0 is amplitude, 1 is mean
  vector< vector<double> > fixToPositiveValuesRuleMatrix = {
    // dont constrain amplitudes
    {0, -999, 0, 999},
    {1, -999, 0, 999},
    {2, -999, 0, 999},
    // dont constrain sigmas
    {0, -999, 2, 999},
    {1, -999, 2, 999},
    {2, -999, 2, 999},
    // dont constrain the ignore parameter
    {0, -999, 4, 999},
    {1, -999, 4, 999},
    {2, -999, 4, 999},

    // put general rules here
    {0, -999, -999, 999},
    {1, -999, -999, 1.0},
    {2, -999, -999, 1.0},

  };


  double rapidityToSwitchToM2Extrapolation[3][3] = {
    {0.85,0.95,2.05},
    {0.65,2.05,2.05},
    {0.55,0.75,2.05}
  };


  //##################################    DON"T EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING (ASK MATT FIRST PLEASE) OR MATT HAS GRADUATED #################################



  // 0 = gaussian    1 = skewNormal    2 = student-t  3 = pearsonIV 4 = doubleGaus
  #ifdef _GAUSSIAN_BTOF_SIMUL_CENT_
    int pionFunctionVersion   = 0;
    int kaonFunctionVersion   = 0;
    int protonFunctionVersion = 0;
  #endif
  #ifdef _SKEW_NORMAL_BTOF_SIMUL_CENT_
    int pionFunctionVersion   = 1;
    int kaonFunctionVersion   = 1;
    int protonFunctionVersion = 1;
  #endif
  #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
    int pionFunctionVersion   = 2;
    int kaonFunctionVersion   = 2;
    int protonFunctionVersion = 2;
  #endif
  #ifdef _PEARSON_IV_BTOF_SIMUL_CENT_
    int pionFunctionVersion   = 3;
    int kaonFunctionVersion   = 3;
    int protonFunctionVersion = 3;
  #endif
  #ifdef _DOUBLE_GAUS_BTOF_SIMUL_CENT_
    int pionFunctionVersion   = 4;
    int kaonFunctionVersion   = 4;
    int protonFunctionVersion = 4;
  #endif






  // MINIMIZER SETUP
  /*ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer();
  ROOT::Math::Functor* chiSqrdFunctor = new ROOT::Math::Functor(this,&ZFitter::simultaneous_centrality_ZTOF_chisqr,m_numCentralities*3 + 12);
  minimizer->SetFunction(*chiSqrdFunctor);
  #ifdef _ZFITTER_DEBUG_
    minimizer->SetPrintLevel(6);
  #else
    minimizer->SetPrintLevel(1);
  #endif
  minimizer->SetMaxIterations(3000);
  minimizer->SetMaxFunctionCalls(200000);
  minimizer->SetValidError(true);*/

  m_currentPartIndex  = a_partIndexSpace;
  m_currentPartCharge = a_charge;
  if(a_charge > 0) m_currentPlusMinusIndex = 0;
  else             m_currentPlusMinusIndex = 1;
  m_inAllCentMode = false;
  m_currentCentIndex = 0; // this is needed to find data for extrapolation procedure
  TH3D* histo3DPM = NULL;
  if(m_currentPartCharge == 1){ // Positive
    histo3DPM = m_rap_mTm0_invBetaBTOF_Plus[m_currentPartIndex][0];
  }else{ // Negative
    histo3DPM = m_rap_mTm0_invBetaBTOF_Minus[m_currentPartIndex][0];
  }

  #ifdef _ZFITTER_DEBUG_
    cout << " using ptr: " << histo3DPM << " for bining information " << endl;
  #endif

  int numHorBins  = histo3DPM->GetNbinsX();
  int numVertBins = histo3DPM->GetNbinsY();
  int rapZeroBin  = histo3DPM->GetXaxis()->FindBin(0.0); // this shouldn't be midrapdity
  if(rapZeroBin > numHorBins) rapZeroBin = numHorBins;
  if(rapZeroBin < 1) rapZeroBin = 1;
  bool rightBinExists = true;
  bool leftBinExists = true;

  if(m_cuts->getEventConfig() == "FixedTarget"){
    rapZeroBin = histo3DPM->GetNbinsX()/2;
  }


  #ifdef _ZFITTER_DEBUG_
    cout << "using y=0 bin: " << rapZeroBin << endl;
  #endif

  double fixed_pion_mass     = m_partInfo->GetParticleMass(0);
  double fixed_kaon_mass     = m_partInfo->GetParticleMass(1);
  double fixed_proton_mass   = m_partInfo->GetParticleMass(2);

  double fixed_pion_mass_sqr     =  pow(m_partInfo->GetParticleMass(0),2.0);
  double fixed_kaon_mass_sqr     =  pow(m_partInfo->GetParticleMass(1),2.0);
  double fixed_proton_mass_sqr   =  pow(m_partInfo->GetParticleMass(2),2.0);






  //###########################   RAPIDITY BIN LOOP   #######################################
  // Fixed 2026-07: was `rightBinExists && leftBinExists` (AND) -- same bug and fix as
  // ZFitter_ProtonTPC_SimulCent_SingleLoop.cxx (see that file's comment for the full
  // explanation). Note the `FixedTarget` special case just above (forcing rapZeroBin to
  // the middle bin) is a workaround for a DIFFERENT manifestation of this same
  // underlying edge-case problem in fixed-target mode -- it doesn't apply here (this
  // repo runs "ColliderCenter" mode) and doesn't fix the loop-logic bug itself, which is
  // what this OR change addresses generally.
  for(int HorBinCounter = 0; rightBinExists || leftBinExists; HorBinCounter++){
    if(HorBinCounter % 2 == 0){
      m_currentRapBin = rapZeroBin - (HorBinCounter/2);
      if(m_currentRapBin < 1){
        leftBinExists = false;
        continue;
      }
    }else{
      m_currentRapBin = rapZeroBin + ((HorBinCounter+1)/2);
      if(m_currentRapBin > numHorBins){
        rightBinExists = false;
        continue;
      }
    }
    if(m_rapBinsToTest.size() > 0){
      bool doThisRapBin = false;
      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(m_currentRapBin == m_rapBinsToTest[iii]){
            doThisRapBin = true;
        }
      }
      if(!doThisRapBin){
        continue;
      }
    }

    m_currentRapidity = histo3DPM->GetXaxis()->GetBinCenter(m_currentRapBin);

    #ifdef _ZFITTER_DEBUG_
      cout << " Working on RapBin:  " << m_currentRapBin << " y = " << m_currentRapidity << endl;
    #endif


    #ifdef _CREATE_TEXT_FIT_FILE_
      string outputTextFileName = Form("%s/%s_m2_RapBin_%02d.txt",
          m_textFileDir.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin);
      ofstream outputTextFile(outputTextFileName.c_str());
      if(outputTextFile.is_open()){
        #ifdef _ZFITTER_DEBUG_
          cout << outputTextFileName << " opened sucessfuly!" << endl;
        #endif
      }else{
        cout << "ERROR: " << outputTextFileName << " wasn't opened properly!!!!!!" << endl;
      }
      outputTextFile.close();
    #endif


    /*TH1D* rebinZBTOFHistoStorage[1000][16];
    for(int iii = 0; iii < 1000; iii++){
      for(int jjj = 0; jjj < 16; jjj++){
        rebinZBTOFHistoStorage[iii][jjj] = NULL;
      }
    }*/

    //##################   Construct Extrapolation Ranges for parameters from the broad matrix as well as momentum to fix to positive fit values
    bool   hasBeenExtrapolatedMatrix[3][2][5]; // part charge param
    double extrapolationMatrix_reduced[3][2][5][2]; // part charge param lowMom/highMom
    double fixToPositiveValuesRuleMatrix_reduced[3][5]; // part param

    for(unsigned int iii = 0; iii < 3; iii++){
      for(unsigned int jjj = 0; jjj < 2; jjj++){
        for(unsigned int kkk = 0; kkk < 5; kkk++){
          hasBeenExtrapolatedMatrix[iii][jjj][kkk] = false;
          extrapolationMatrix_reduced[iii][jjj][kkk][0] = 999;
          extrapolationMatrix_reduced[iii][jjj][kkk][1] = 999;
          fixToPositiveValuesRuleMatrix_reduced[iii][kkk] = 999;
        }
      }
    }

    //particle charge rapidity param lowMom highMom
    for(int extrapIndex = extrapolationMatrix.size() - 1; extrapIndex >= 0; extrapIndex--){
      if(extrapolationMatrix[extrapIndex][2] < -990 || fabs(fabs(m_currentRapidity) - extrapolationMatrix[extrapIndex][2]) < 0.01){
        int partIndex = ((int)(extrapolationMatrix[extrapIndex][0] + 0.01));
        int pmIndex = 0;
        if(extrapolationMatrix[extrapIndex][1] < 0.0) pmIndex = 1;
        int paramIndexRead = ((int)(extrapolationMatrix[extrapIndex][3] + 0.01));
        bool allParams = (paramIndexRead < -990);
        for(unsigned int paramIndex = 0; paramIndex < 5; paramIndex++){
          if(allParams || (paramIndexRead == paramIndex)){
            extrapolationMatrix_reduced[partIndex][pmIndex][paramIndex][0] = extrapolationMatrix[extrapIndex][4];
            extrapolationMatrix_reduced[partIndex][pmIndex][paramIndex][1] = extrapolationMatrix[extrapIndex][5];
          }
        }
      }
    }
    //particle rapidity param momentum
    for(int index = fixToPositiveValuesRuleMatrix.size() - 1; index >= 0; index--){
      if(fixToPositiveValuesRuleMatrix[index][1] < -990 || fabs(fabs(m_currentRapidity) - fixToPositiveValuesRuleMatrix[index][1]) < 0.01){
        int partIndex = ((int)(fixToPositiveValuesRuleMatrix[index][0] + 0.01));
        int paramIndexRead = ((int)(fixToPositiveValuesRuleMatrix[index][2] + 0.01));
        bool allParams = (paramIndexRead < -990);
        for(unsigned int paramIndex = 0; paramIndex < 5; paramIndex++){
          if(allParams || (paramIndexRead == paramIndex)){
            fixToPositiveValuesRuleMatrix_reduced[partIndex][paramIndex] = fixToPositiveValuesRuleMatrix[index][3];
          }
        }
      }
    }

    #ifdef _ZFITTER_DEBUG_
      cout << "BTOF Parameter Extrapolation Matrix: " << endl;
      for(unsigned int iii = 0; iii < 3; iii++){
        for(unsigned int jjj = 0; jjj < 2; jjj++){
          int charge = 1;
          if(jjj == 1) charge = -1;
          for(unsigned int kkk = 0; kkk < 5; kkk++){
            cout << "   " << m_partInfo->GetParticleName(iii,charge) << " parameter " << kkk << "  p_Tot = [" << extrapolationMatrix_reduced[iii][jjj][kkk][0] << ", " << extrapolationMatrix_reduced[iii][jjj][kkk][1] << "]";
            double mTm0_low  = PhysMath::mTm0_from_y_p(m_currentRapidity, extrapolationMatrix_reduced[iii][jjj][kkk][0], m_partInfo->GetParticleMass(m_currentPartIndex));
            double mTm0_high = PhysMath::mTm0_from_y_p(m_currentRapidity, extrapolationMatrix_reduced[iii][jjj][kkk][1], m_partInfo->GetParticleMass(m_currentPartIndex));
            cout << "  mT-m0 = [" << mTm0_low << ", " << mTm0_high << "]" << endl;
          }
        }
      }
      cout << "BTOF Positive Parameter Constraint Matrix: " << endl;
      for(unsigned int iii = 0; iii < 3; iii++){
        for(unsigned int kkk = 0; kkk < 5; kkk++){
          cout << "   " << m_partInfo->GetParticleName(iii) << " parameter " << kkk << "  p_Tot = " << fixToPositiveValuesRuleMatrix_reduced[iii][kkk];
          double mTm0  = PhysMath::mTm0_from_y_p(m_currentRapidity, fixToPositiveValuesRuleMatrix_reduced[iii][kkk], m_partInfo->GetParticleMass(m_currentPartIndex));
          cout << "   mT-m0 = " << mTm0 << endl;
        }
      }
    #endif

    const unsigned int prevParamSize = 15;
    double previousGoodParameters[prevParamSize] = {-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999};
    int firstMtM0BinFit       = -1;
    int lastMtM0BinFit        = -1;
    int firstProtonMtM0BinFit = -1;

    bool previousFitToggle_Pion   = false;
    bool previousFitToggle_Kaon   = false;
    bool previousFitToggle_Proton = false;

    double previousCentralityIntegratedYields_Pion = 0;
    double previousCentralityIntegratedYields_Kaon = 0;
    double previousCentralityIntegratedYields_Proton = 0;

    double centIntThreshold = 200.0;
    bool centralityIntegratedYieldsWereAboveThreshold_Pion = false;
    bool centralityIntegratedYieldsWereAboveThreshold_Kaon = false;
    bool centralityIntegratedYieldsWereAboveThreshold_Proton = false;

    //double fixedKaonStudentTNuValue = -999; //

    //###################################################    mT-m0 bin loop  #########################################################################

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentMtM0Bin  = VertBinIndex;
      m_currentMtM0     = histo3DPM->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      double eta        = PhysMath::eta_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      double theta      = 0.5*TMath::Pi() - PhysMath::theta_from_eta(eta);

      double fixed_pion_invBeta     =  PhysMath::invBeta_from_MassSquaredSigned(fixed_pion_mass_sqr,   m_currentMomOfBin);
      double fixed_kaon_invBeta     =  PhysMath::invBeta_from_MassSquaredSigned(fixed_kaon_mass_sqr,   m_currentMomOfBin);
      double fixed_proton_invBeta   =  PhysMath::invBeta_from_MassSquaredSigned(fixed_proton_mass_sqr, m_currentMomOfBin);


      #ifdef _ZFITTER_DEBUG_
        cout << "\n\n\n###########################################################################" << endl;
        cout << "Working on ZBTOF STEP 3: " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data()
             << " RapBin: " << m_currentRapBin << " MtM0Bin: " << m_currentMtM0Bin << endl;
      #endif

      //################ SETUP THE HISTOGRAM TO FIT ################
      // this histogram is just for setup, see _byCent[] below for the actual fit histograms
      m_currentHistoToFit = histo3DPM->ProjectionZ(
         Form("btof%s_cent%02d_%02d_%02d",m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),0,m_currentRapBin,m_currentMtM0Bin),
         m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
      m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(m_currentHistoToFit);
      if(m_currentIntegralofHistogram < 30){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping this bin because of low counts: " << m_currentIntegralofHistogram << " < 30"  << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }
      if(m_convertInvBetaToMassSquared){
        TH1D* invBetaToM2Histo = HistogramUtilities::convertZTOFToMassSqr(m_currentHistoToFit,m_currentMomOfBin,m_cuts->getMassBTOF(m_currentPartIndex),false,false);
        delete m_currentHistoToFit;
        m_currentHistoToFit = invBetaToM2Histo;
      }

      TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(m_currentHistoToFit,80,false,true);
      delete m_currentHistoToFit;
      m_currentHistoToFit = rebinHisto;
      HistogramUtilities::removeHighCountCorruptData(m_currentHistoToFit);
      m_currentHistoToFit->Scale(1.0,"width"); // make differential
      double lowRangeData  = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistoToFit));
      double highRangeData = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistoToFit)+1);

      //############## EXTRAPOLATION OF PARAMETERS  ##########################################
      if(m_currentPartCharge > 0 || (!usePositiveChargeAsSeeds && m_currentPartCharge < 0)){
        for(int partIndex = 0; partIndex < 3; partIndex++){
          int functVer = pionFunctionVersion;
          if(partIndex == 1) functVer = kaonFunctionVersion;
          if(partIndex == 2) functVer = protonFunctionVersion;
          for(int paramIndex = 0; paramIndex < 5; paramIndex++){
            if(!hasBeenExtrapolatedMatrix[partIndex][m_currentPlusMinusIndex][paramIndex] && m_currentMomOfBin > extrapolationMatrix_reduced[partIndex][m_currentPlusMinusIndex][paramIndex][1]){
              bool useM2ExtrapForMean = false;
              if(rapidityToSwitchToM2Extrapolation[m_currentPartIndex][partIndex] < fabs(m_currentRapidity)) useM2ExtrapForMean = true;
              hasBeenExtrapolatedMatrix[partIndex][m_currentPlusMinusIndex][paramIndex] = extrapolate_BTOF_Fit_Parameters(partIndex,
                    extrapolationMatrix_reduced[partIndex][m_currentPlusMinusIndex][paramIndex][0], extrapolationMatrix_reduced[partIndex][m_currentPlusMinusIndex][paramIndex][1],
                    functVer, -1, true, paramIndex,useM2ExtrapForMean); // -1 is all cent and true is invBeta->m2 mapping
            }
          }
        }
      }

      /*
      if(m_currentMomOfBin > pionHighMomParamExtrap && !pionHasBeenExtrapolated && fabs(m_currentRapidity) < rapidityToNotDoParameterExtrapolation){
      }
      if(m_currentMomOfBin > kaonHighMomParamExtrap && !kaonHasBeenExtrapolated && fabs(m_currentRapidity) < rapidityToNotDoParameterExtrapolation){
        kaonHasBeenExtrapolated = extrapolate_BTOF_Fit_Parameters(1, kaonLowMomParamExtrap, kaonHighMomParamExtrap, kaonFunctionVersion, -1, true); // -1 is all cent and true is invBeta->m2 mapping
      }
      if( !protonHasBeenExtrapolated && m_currentPartCharge == -1 && m_currentMomOfBin > antiprotonHighMomParamExtrap && fabs(m_currentRapidity) < rapidityToNotDoParameterExtrapolation){
        protonHasBeenExtrapolated = extrapolate_BTOF_Fit_Parameters(2, antiprotonLowMomParamExtrap, antiprotonHighMomParamExtrap, protonFunctionVersion, -1, true); // -1 is all cent and true is invBeta->m2 mapping
      }
      if( !protonHasBeenExtrapolated && m_currentPartCharge == 1 && m_currentMomOfBin > protonHighMomParamExtrap &&  fabs(m_currentRapidity) < rapidityToNotDoParameterExtrapolation){
        protonHasBeenExtrapolated = extrapolate_BTOF_Fit_Parameters(2, protonLowMomParamExtrap, protonHighMomParamExtrap, protonFunctionVersion, -1, true); // -1 is all cent and true is invBeta->m2 mapping
      }
      */

      //#################      FIND THE MEAN OF PI / K / E
      double pionPredict     = fixed_pion_invBeta;
      double kaonPredict     = fixed_kaon_invBeta;
      double protonPredict   = fixed_proton_invBeta;
      //double m2_width_pion   = PhysMath::m2_width_prediction(invBetaResolution,m_currentMomOfBin,fixed_pion_mass);
      //double m2_width_kaon   = PhysMath::m2_width_prediction(invBetaResolution,m_currentMomOfBin,fixed_kaon_mass);
      //double m2_width_proton = PhysMath::m2_width_prediction(invBetaResolution,m_currentMomOfBin,fixed_proton_mass);
      double m2_width_pion   = m2width_pion_graph->Eval(m_currentMomOfBin);
      double m2_width_kaon   = m2width_kaon_graph->Eval(m_currentMomOfBin);
      double m2_width_proton = m2width_proton_graph->Eval(m_currentMomOfBin);


      double lowRange  = -0.8;
      double highRange = 1.4;
      bool doPion   = true;
      bool doKaon   = true;
      bool doProton = true;
      TH1D* truncHisto    = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
      TH1D* orignalHisto  = m_currentHistoToFit;
      m_currentHistoToFit = truncHisto;
      delete orignalHisto;
      //rebinZBTOFHistoStorage[m_currentMtM0Bin] = m_currentHistoToFit;

      // From the Cent00 Histo, make the other Cent Histos with the same range
      for(int centIndex = 0; centIndex < 16; centIndex++){
        m_currentHistosToFit_ByCent[centIndex] = nullptr;
        m_currentHistosToFit_ByCent_InvBeta[centIndex] = nullptr;
      }
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        m_currentHistosToFit_ByCent[centIndex] = grabSpecificHistogram(m_currentPartIndex, m_currentPartCharge, centIndex, 1, m_currentRapBin, m_currentMtM0Bin, 60, false, -999, -999, true,false);
        if(!m_currentHistosToFit_ByCent[centIndex]){
          #ifdef _ZFITTER_DEBUG_
            cout << " ERROR: Could not grabSpecificHistogram... Cent " << centIndex << " probably doesn't exist    ptr: " << m_currentHistosToFit_ByCent[centIndex] << endl;
          #endif
          continue;
        }
        double counts = HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[centIndex], -1.0, true) * m_centEvents[m_currentPartIndex]->GetBinContent(centIndex+1);
        #ifdef _ZFITTER_DEBUG_
          cout << " Cent Index: " << centIndex << "  Counts: " << counts << "     N_Evt: " << m_centEvents[m_currentPartIndex]->GetBinContent(centIndex+1) << endl;
        #endif
        if(counts < 300){ // if there are too few counts for 80 bins, try again with 20-80 bins
          delete m_currentHistosToFit_ByCent[centIndex];
          int numNewBins = counts / 20.0;
          if(numNewBins < 20) numNewBins = 20;
          if(numNewBins > 60) numNewBins = 60;
          m_currentHistosToFit_ByCent[centIndex] = grabSpecificHistogram(m_currentPartIndex, m_currentPartCharge, centIndex, 1, m_currentRapBin, m_currentMtM0Bin, numNewBins, false, -999, -999, true,false);
        }

        if(m_convertInvBetaToMassSquared){
          TH1D* invBetaToM2Histo = HistogramUtilities::convertZTOFToMassSqr(m_currentHistosToFit_ByCent[centIndex],m_currentMomOfBin,m_cuts->getMassBTOF(m_currentPartIndex),true,true);
          delete m_currentHistosToFit_ByCent[centIndex];
          m_currentHistosToFit_ByCent[centIndex] = invBetaToM2Histo;

          m_currentHistosToFit_ByCent_InvBeta[centIndex] = HistogramUtilities::convertMassSqrToInvBeta(invBetaToM2Histo,m_currentMomOfBin,true,true);

        }

        //truncHisto    = HistogramUtilities::truncateHistogram(m_currentHistosToFit_ByCent[centIndex],lowRange,highRange);
        //orignalHisto  = m_currentHistosToFit_ByCent[centIndex];
        //m_currentHistosToFit_ByCent[centIndex] = truncHisto;
        //delete orignalHisto;
        //rebinZBTOFHistoStorage[centIndex][m_currentMtM0Bin] = m_currentHistosToFit_ByCent[centIndex];
        m_currentHistosToFit_ByCent[centIndex]->SetTitle(
          Form("%s Cent %d y=[%1.2f,%1.2f]  m_{T}-m_{0}=[%2.3f,%2.3f] p_{T}=[%2.3f,%2.3f] p=%2.3f #eta=%2.3f #theta=%2.3f N=%.2e; m^{2}; #frac{1}{N_{events}} #frac{dN}{dm^{2}}",
            m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data(),
            centIndex,
            histo3DPM->GetXaxis()->GetBinLowEdge(m_currentRapBin),
            histo3DPM->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
            histo3DPM->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
            histo3DPM->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
            PhysMath::pT_from_mTm0(histo3DPM->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),m_partInfo->GetParticleMass(m_currentPartIndex)),
            PhysMath::pT_from_mTm0(histo3DPM->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),m_partInfo->GetParticleMass(m_currentPartIndex)),
            m_currentMomOfBin, eta, theta, counts
            ));
      }




      #ifdef _ZFITTER_DEBUG_
        cout << " Step03 initial predictions | " << "pi: " << pionPredict << " k: " << kaonPredict << " p: " << protonPredict << endl;
        cout << "   truncated from " << lowRange << " to " << highRange << endl;
      #endif

      double pionAmpPredict     = 1.0;
      double pionSigmaInit      = 1.0;
      double pionSkewInit       = 1.0;
      double pionNuInit         = 0.0;

      double kaonAmpPredict     = 1.0;
      double kaonSigmaInit      = 1.0;
      double kaonSkewInit       = 1.0;  // 3rd param of whatever distribution
      double kaonNuInit         = 0.0;  // 4th param of whatever distribution

      double protonAmpPredict   = 1.0;
      double protonSigmaInit    = 1.0;
      double protonSkewInit     = 1.0;
      double protonNuInit       = 0.0;

      //    pion         koan         proton   Cent Amplitudes (pi/k/p)
      //  0 1 2 3       4 5 6 7     8 9 10 11    12  13  14 ...

      //######################    LOAD PREVIOUS GOOD PARAMETERS AS SEEDS   ###########################
      bool useThePrevPar = usePreviousGoodParams_Step03;
      if(m_currentPartIndex == 1 && m_currentMomOfBin < momToStartPreviousGoodParameters_ForKaon)   useThePrevPar = false;
      if(m_currentPartIndex == 2 && m_currentMomOfBin < momToStartPreviousGoodParameters_ForProton) useThePrevPar = false;
      if(previousGoodParameters[4*a_partIndexSpace + 1] <= 0.0) useThePrevPar = false;

      bool usedPionPrevParams = false;
      bool usedKaonPrevParams = false;
      bool usedProtonPrevParams = false;

      if(useThePrevPar){ // use the previous fit parameters as seeds

        if(std::isfinite(previousGoodParameters[0]) && std::isfinite(previousGoodParameters[1]) &&
          std::isfinite(previousGoodParameters[2]) && std::isfinite(previousGoodParameters[3])){
          //previousGoodParameters[0] < previousGoodParameters[4] &&
          //previousGoodParameters[0] < previousGoodParameters[8] &&
          //fabs(previousGoodParameters[0] - fixed_pion_invBeta) < 1.5*invBetaResolution){

          pionPredict    = PhysMath::invBeta_from_MassSquaredSigned(previousGoodParameters[0], m_currentMomOfBin);
          pionSigmaInit  = previousGoodParameters[1];
          pionSkewInit   = previousGoodParameters[2];
          pionNuInit     = previousGoodParameters[3];
          usedPionPrevParams = true;
        }
        if(std::isfinite(previousGoodParameters[4]) && std::isfinite(previousGoodParameters[5]) &&
           std::isfinite(previousGoodParameters[6]) && std::isfinite(previousGoodParameters[7])){
           // &&
           //previousGoodParameters[0] < previousGoodParameters[4] &&
           //previousGoodParameters[4] < previousGoodParameters[8] &&
           //fabs(previousGoodParameters[4] - fixed_kaon_invBeta) < 1.5*invBetaResolution){
          kaonPredict    = PhysMath::invBeta_from_MassSquaredSigned(previousGoodParameters[4], m_currentMomOfBin);
          kaonSigmaInit  = previousGoodParameters[5];
          kaonSkewInit   = previousGoodParameters[6];
          kaonNuInit     = previousGoodParameters[7];
          usedKaonPrevParams = true;
        }
        if(std::isfinite(previousGoodParameters[8]) && std::isfinite(previousGoodParameters[9]) &&
           std::isfinite(previousGoodParameters[10]) && std::isfinite(previousGoodParameters[11])){
           // &&
           //previousGoodParameters[0] < previousGoodParameters[8] &&
           //previousGoodParameters[4] < previousGoodParameters[8] &&
           //fabs(previousGoodParameters[8] - fixed_proton_invBeta) < 1.5*invBetaResolution){
          protonPredict    = PhysMath::invBeta_from_MassSquaredSigned(previousGoodParameters[8], m_currentMomOfBin);
          protonSigmaInit  = previousGoodParameters[9];
          protonSkewInit   = previousGoodParameters[10];
          protonNuInit     = previousGoodParameters[11];
          usedProtonPrevParams = true;
        }

        #ifdef _ZFITTER_DEBUG_
          cout << "  Tried Using previous good parameters..." << endl;
          cout << "      Pion:   " << (usedPionPrevParams ? " used  " : "NOT USED  ")   << previousGoodParameters[0]  << "\t" << previousGoodParameters[1]  << "\t" << previousGoodParameters[2]  << "\t" << previousGoodParameters[3] << "\t" << previousGoodParameters[4] << endl;
          cout << "      Kaon:   " << (usedKaonPrevParams ? " used  " : "NOT USED  ")   << previousGoodParameters[5]  << "\t" << previousGoodParameters[6]  << "\t" << previousGoodParameters[7]  << "\t" << previousGoodParameters[8] << "\t" << previousGoodParameters[9] << endl;
          cout << "      Proton: " << (usedProtonPrevParams ? " used  " : "NOT USED  ") << previousGoodParameters[10] << "\t" << previousGoodParameters[11] << "\t" << previousGoodParameters[12] << "\t" << previousGoodParameters[13] << "\t" << previousGoodParameters[14] << endl;
        #endif

      }else{
        #ifdef _ZFITTER_DEBUG_
          cout << "NOT Using Previous Good Parameters..." << endl;
        #endif
      }

      // Even if you have good parameters, you might be able to fit a particle you couldn't before, so try again
      if(!previousFitToggle_Pion || !usedPionPrevParams){
        int maxPionBin = HistogramUtilities::getMaxBinInRange(m_currentHistoToFit, 0.0, 0.1);
        double pion_m2 = m_currentHistoToFit->GetBinCenter(maxPionBin);
        pionPredict    = PhysMath::invBeta_from_MassSquaredSigned(pion_m2, m_currentMomOfBin);
        double pionData[3];
        //HistogramUtilities::getGausData(&pionData[0], m_currentHistoToFit, true,false, pion_m2 - 1.5*m2_width_pion, pion_m2 + 1.5*m2_width_pion);
        HistogramUtilities::getGausData(&pionData[0], m_currentHistoToFit, true,false, -0.2, 0.1);

        #ifdef _ZFITTER_DEBUG_
          cout << " Pion Gaus Data: " << pionData[0] << " " << pionData[1] << " " << pionData[2] << endl;
        #endif
        if(std::isfinite(pionData[1]) && std::isfinite(pionData[2])){
          pionPredict      = PhysMath::invBeta_from_MassSquaredSigned(pionData[1], m_currentMomOfBin);
          pionSigmaInit    = PhysMath::invBetaWidth_from_MassSquaredWidth(pionData[1], pionData[2], m_currentMomOfBin);
        }
        if(pionSigmaInit > 1.5*invBetaResolution || pionSigmaInit < 0.5*invBetaResolution){
          pionSigmaInit = invBetaResolution;
        }

        //if(std::isfinite(pionData[1]) && std::isfinite(pionData[2])){
        //  pionPredict        = pionData[1];
        //  pionSigmaInit      = pionData[2];
        //}
        //if(pionSigmaInit > 1.2*m2_width_pion || pionSigmaInit > 0.8*m2_width_pion){
        //  pionSigmaInit = m2_width_pion;
        //}
        if(pionFunctionVersion == 1)      pionSkewInit = 0.0;
        else if(pionFunctionVersion == 2) pionSkewInit = 5.0;
        else if(pionFunctionVersion == 3){
          pionSkewInit = 5.0;
          pionNuInit = 0.1;
        }else if(pionFunctionVersion == 4){ //doubleGaus
          pionSkewInit = 0.5; // frac
          pionNuInit = 1.01*pionSigmaInit;
        }
      }

      if(!previousFitToggle_Kaon || !usedKaonPrevParams){
        double kaonData[3];
        //HistogramUtilities::getGausData(&kaonData[0], m_currentHistoToFit, true,false, fixed_kaon_mass_sqr - 1.5*m2_width_kaon, fixed_kaon_mass_sqr + 1.5*m2_width_kaon);
        HistogramUtilities::getGausData(&kaonData[0], m_currentHistoToFit, true,false, 0.1, 0.5);

        #ifdef _ZFITTER_DEBUG_
          cout << " Kaon Gaus Data: " << kaonData[0] << " " << kaonData[1] << " " << kaonData[2] << endl;
        #endif
        if(std::isfinite(kaonData[1]) && std::isfinite(kaonData[2])){
          kaonPredict      = PhysMath::invBeta_from_MassSquaredSigned(kaonData[1], m_currentMomOfBin);
          kaonSigmaInit    = PhysMath::invBetaWidth_from_MassSquaredWidth(kaonData[1], kaonData[2], m_currentMomOfBin);
        }
        if(kaonSigmaInit > 1.5*invBetaResolution || kaonSigmaInit < 0.5*invBetaResolution){
          kaonSigmaInit = invBetaResolution;
        }

        //if(std::isfinite(kaonData[1]) && std::isfinite(kaonData[2])){
        //  kaonPredict   = kaonData[1];
        //  kaonSigmaInit = kaonData[2];
        //}
        //if(kaonSigmaInit > 1.2*m2_width_kaon || kaonSigmaInit > 0.8*m2_width_kaon){
        //  kaonSigmaInit = m2_width_kaon;
        //}
        if(kaonFunctionVersion == 1)      kaonSkewInit = 0.0;
        else if(kaonFunctionVersion == 2) kaonSkewInit = 5.0;
        else if(kaonFunctionVersion == 3){
          kaonSkewInit = 5.0;
          kaonNuInit = 0.1;
        }else if(kaonFunctionVersion == 4){ //doubleGaus
          kaonSkewInit = 0.5; // frac
          kaonNuInit = 1.01*pionSigmaInit;
        }

      }

      if(!previousFitToggle_Proton || !usedProtonPrevParams){
        double protonData[3];
        //HistogramUtilities::getGausData(&protonData[0], m_currentHistoToFit, true,false, fixed_proton_mass_sqr - 1.5*m2_width_proton, fixed_proton_mass_sqr + 1.5*m2_width_proton);
        HistogramUtilities::getGausData(&protonData[0], m_currentHistoToFit, true,false, 0.5, 1.5);
        //if(protonData[0] < 1.0) doProton = false;
        #ifdef _ZFITTER_DEBUG_
          cout << " Proton Gaus Data: " << protonData[0] << " " << protonData[1] << " " << protonData[2] << endl;
        #endif
        if(std::isfinite(protonData[1]) && std::isfinite(protonData[2])){
          protonPredict      = PhysMath::invBeta_from_MassSquaredSigned(protonData[1], m_currentMomOfBin);
          protonSigmaInit    = PhysMath::invBetaWidth_from_MassSquaredWidth(protonData[1], protonData[2], m_currentMomOfBin);
          cout << " Converted Gaus Data: " << protonData[0] << " " << protonPredict << " " << protonSigmaInit << endl;
        }
        if(protonSigmaInit > 1.5*invBetaResolution || protonSigmaInit < 0.5*invBetaResolution){
          protonSigmaInit = invBetaResolution;
        }
        //if(protonSigmaInit > 1.2*m2_width_proton || protonSigmaInit < 0.8*m2_width_proton){
        //  protonSigmaInit = m2_width_proton;
        //}

        if(protonFunctionVersion == 1)       protonSkewInit = 0.0;
        else if(protonFunctionVersion == 2)  protonSkewInit = 5.0;
        else if(protonFunctionVersion == 3){
          protonSkewInit = 5.0;
          protonNuInit   = 0.1;
        }else if(protonFunctionVersion == 4){ //doubleGaus
          protonSkewInit = 0.5; // frac
          protonNuInit   = 1.01*pionSigmaInit;
        }
      }

      // If extrapolated, use the extrapolation data (positive seeds will override)
      if(hasBeenExtrapolatedMatrix[0][m_currentPlusMinusIndex][1])  pionPredict     = m_Fit_Data_ZbTOF[m_currentPartIndex][0][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[0][m_currentPlusMinusIndex][2])  pionSigmaInit   = m_Fit_Data_ZbTOF[m_currentPartIndex][0][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[0][m_currentPlusMinusIndex][3])  pionSkewInit    = m_Fit_Data_ZbTOF[m_currentPartIndex][0][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[0][m_currentPlusMinusIndex][4])  pionNuInit      = m_Fit_Data_ZbTOF[m_currentPartIndex][0][4][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[1][m_currentPlusMinusIndex][1])  kaonPredict     = m_Fit_Data_ZbTOF[m_currentPartIndex][1][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[1][m_currentPlusMinusIndex][2])  kaonSigmaInit   = m_Fit_Data_ZbTOF[m_currentPartIndex][1][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[1][m_currentPlusMinusIndex][3])  kaonSkewInit    = m_Fit_Data_ZbTOF[m_currentPartIndex][1][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[1][m_currentPlusMinusIndex][4])  kaonNuInit      = m_Fit_Data_ZbTOF[m_currentPartIndex][1][4][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[2][m_currentPlusMinusIndex][1])  protonPredict   = m_Fit_Data_ZbTOF[m_currentPartIndex][2][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[2][m_currentPlusMinusIndex][2])  protonSigmaInit = m_Fit_Data_ZbTOF[m_currentPartIndex][2][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[2][m_currentPlusMinusIndex][3])  protonSkewInit  = m_Fit_Data_ZbTOF[m_currentPartIndex][2][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(hasBeenExtrapolatedMatrix[2][m_currentPlusMinusIndex][4])  protonNuInit    = m_Fit_Data_ZbTOF[m_currentPartIndex][2][4][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);

      #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
        if(pionSkewInit < 1.0001) pionSkewInit = 1.0001;
        if(kaonSkewInit < 1.0001) kaonSkewInit = 1.0001;
        if(protonSkewInit < 1.0001) protonSkewInit = 1.0001;
        if(m_currentPartIndex == 0 && pionSkewInit < 2.0)   pionSkewInit = 2.0;
        if(m_currentPartIndex == 1 && kaonSkewInit < 2.0)   kaonSkewInit = 2.0;
        if(m_currentPartIndex == 2 && protonSkewInit < 2.0) protonSkewInit = 2.0;
      #endif

      // Use the positive seeds if available
      if(usePositiveChargeAsSeeds && m_currentPartCharge < 0){
        #ifdef _ZFITTER_DEBUG_
          cout << "Getting Parameters from Positive Fits... " << endl;
        #endif

        double pionMeanTemp      = 0;
        double pionSigmaTemp     = 0;
        double pionP3Temp        = 0;
        double pionP4Temp        = 0;
        double kaonMeanTemp      = 0;
        double kaonSigmaTemp     = 0;
        double kaonP3Temp        = 0;
        double kaonP4Temp        = 0;
        double protonMeanTemp    = 0;
        double protonSigmaTemp   = 0;
        double protonP3Temp      = 0;
        double protonP4Temp      = 0;
        double pionMeanTemp_Err    = 0;
        double pionSigmaTemp_Err   = 0;
        double pionP3Temp_Err      = 0;
        double pionP4Temp_Err      = 0;
        double kaonMeanTemp_Err    = 0;
        double kaonSigmaTemp_Err   = 0;
        double kaonP3Temp_Err      = 0;
        double kaonP4Temp_Err      = 0;
        double protonMeanTemp_Err  = 0;
        double protonSigmaTemp_Err = 0;
        double protonP3Temp_Err    = 0;
        double protonP4Temp_Err    = 0;

        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][0][1][0],m_currentRapBin, m_currentMtM0, pionMeanTemp,   pionMeanTemp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][0][2][0],m_currentRapBin, m_currentMtM0, pionSigmaTemp,  pionSigmaTemp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][0][3][0],m_currentRapBin, m_currentMtM0, pionP3Temp,     pionP3Temp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][0][4][0],m_currentRapBin, m_currentMtM0, pionP4Temp,     pionP4Temp_Err);

        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][1][1][0],m_currentRapBin, m_currentMtM0, kaonMeanTemp,   kaonMeanTemp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][1][2][0],m_currentRapBin, m_currentMtM0, kaonSigmaTemp,  kaonSigmaTemp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][1][3][0],m_currentRapBin, m_currentMtM0, kaonP3Temp,     kaonP3Temp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][1][4][0],m_currentRapBin, m_currentMtM0, kaonP4Temp,     kaonP4Temp_Err);

        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][2][1][0],m_currentRapBin, m_currentMtM0, protonMeanTemp,  protonMeanTemp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][2][2][0],m_currentRapBin, m_currentMtM0, protonSigmaTemp, protonSigmaTemp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][2][3][0],m_currentRapBin, m_currentMtM0, protonP3Temp,    protonP3Temp_Err);
        HistogramUtilities::interpolate2D(m_Fit_Data_ZbTOF[m_currentPartIndex][2][4][0],m_currentRapBin, m_currentMtM0, protonP4Temp,    protonP4Temp_Err);



        bool pionWasFit = false;
        if(pionMeanTemp_Err > 0.0 || pionSigmaTemp_Err > 0.0) pionWasFit = true;
        bool kaonWasFit = false;
        if(kaonMeanTemp_Err > 0.0 || kaonSigmaTemp_Err > 0.0) kaonWasFit = true;
        bool protonWasFit = false;
        if(protonMeanTemp_Err > 0.0 || protonSigmaTemp_Err > 0.0) protonWasFit = true;

        #ifdef _ZFITTER_DEBUG_
          cout << " Before Positive Parameters Introduced:  " << endl;
          cout << "     Pion     mu:" << pionPredict   << "  sig: " << pionSigmaInit   << "  nu: " << pionSkewInit << "  p: " << pionNuInit  << endl;
          cout << "     Kaon     mu:" << kaonPredict << "  sig: " << kaonSigmaInit << "  nu: " << kaonSkewInit << "  p: " << kaonNuInit  << endl;
          cout << "     Proton   mu:" << protonPredict << "  sig: " << protonSigmaInit << "  nu: " << protonSkewInit << "  p: " << protonNuInit  << endl;
          cout << " After Positive Parameters Introduced:  " << endl;
        #endif

        if(pionSigmaTemp > 0.0 && pionWasFit){
          pionPredict   = pionMeanTemp;
          pionSigmaInit = pionSigmaTemp;
          pionSkewInit  = pionP3Temp;
          pionNuInit    = pionP4Temp;

          #ifdef _ZFITTER_DEBUG_
            cout << "     Pion     mu:" << pionPredict   << "  sig: " << pionSigmaInit   << "  nu: " << pionSkewInit << "  p: " << pionNuInit  << endl;
          #endif

        }
        if(kaonSigmaTemp > 0.0 && kaonWasFit){
          kaonPredict   = kaonMeanTemp;
          kaonSigmaInit = kaonSigmaTemp;
          kaonSkewInit  = kaonP3Temp;
          kaonNuInit    = kaonP4Temp;

          #ifdef _ZFITTER_DEBUG_
            cout << "     Kaon     mu:" << kaonPredict << "  sig: " << kaonSigmaInit << "  nu: " << kaonSkewInit << "  p: " << kaonNuInit  << endl;
          #endif

        }
        if(protonSigmaTemp > 0.0 && protonWasFit){
          protonPredict   = protonMeanTemp;
          protonSigmaInit = protonSigmaTemp;
          protonSkewInit  = protonP3Temp;
          protonNuInit    = protonP4Temp;

          #ifdef _ZFITTER_DEBUG_
            cout << "     Proton   mu:" << protonPredict << "  sig: " << protonSigmaInit << "  nu: " << protonSkewInit << "  p: " << protonNuInit  << endl;
          #endif

        }

      }



      if(!std::isfinite(pionSigmaInit)   || pionSigmaInit == 0.0 || HistogramUtilities::outsideHistogramRange(m_currentHistoToFit,fixed_pion_mass_sqr)){
        cout << "ERROR: Pion initialization failed... Not fitting..." << pionPredict << " " << pionSigmaInit << " " << pionSkewInit << " " << pionNuInit << endl;
        doPion = false;
      }
      if(!std::isfinite(kaonSigmaInit)   || kaonSigmaInit == 0.0 || HistogramUtilities::outsideHistogramRange(m_currentHistoToFit,fixed_kaon_mass_sqr)){
        cout << "ERROR: Kaon initialization failed... Not fitting..."  << kaonPredict << " " << kaonSigmaInit << " " << kaonSkewInit << " " << kaonNuInit << endl;
        doKaon = false;
      }
      if(!std::isfinite(protonSigmaInit) || protonSigmaInit == 0.0 || HistogramUtilities::outsideHistogramRange(m_currentHistoToFit,fixed_proton_mass_sqr)){
        cout << "ERROR: Proton initialization failed... Not fitting..." << protonPredict << " " << protonSigmaInit << " " << protonSkewInit << " " << protonNuInit << endl;
        doProton = false;
      }

      if(HistogramUtilities::outsideHistogramRange(m_currentHistoToFit,fixed_pion_mass_sqr - 0.8*m2_width_pion)){
        cout << "WARNING: Pion Not Within Histogram Range... Not fitting..." << pionPredict << " " << pionSigmaInit << " " << pionSkewInit << " " << pionNuInit << endl;
        doPion = false;
      }
      if(HistogramUtilities::outsideHistogramRange(m_currentHistoToFit,fixed_kaon_mass_sqr - 0.8*m2_width_kaon)){
        cout << "WARNING: Kaon Not Within Histogram Range... Not fitting..."  << kaonPredict << " " << kaonSigmaInit << " " << kaonSkewInit << " " << kaonNuInit << endl;
        doKaon = false;
      }
      if(HistogramUtilities::outsideHistogramRange(m_currentHistoToFit,fixed_proton_mass_sqr + 0.8*m2_width_proton)){
        cout << "WARNING: Proton Not Within Histogram Range... Not fitting..." << protonPredict << " " << protonSigmaInit << " " << protonSkewInit << " " << protonNuInit << endl;
        doProton = false;
      }

      //double maxValInHisto = HistogramUtilities::getMaxInRange(m_currentHistoToFit);


      //DETERMINE THE AMPLITUDES
      double amplitudes[16][3];
      TGraph* initialPoints[16][3];
      for(int iii = 0; iii < 16; iii++){
        for(int jjj = 0; jjj < 3; jjj++){
          initialPoints[iii][jjj] = nullptr;
        }
      }
      double defaultAmplitudes[16];
      double defaultMaxAmplitude[16];
      double wideRangeAmplitudes[16][3];
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        for(int jjj = 0; jjj < 3; jjj++){
          initialPoints[centIndex][jjj] = new TGraph();
          initialPoints[centIndex][jjj]->SetMarkerStyle(20);
          initialPoints[centIndex][jjj]->SetMarkerColor(m_partInfo->GetParticleColor(jjj));
        }
        TH1D* histo = m_currentHistosToFit_ByCent[centIndex];
        defaultMaxAmplitude[centIndex] = HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[centIndex],-1,true);
        defaultAmplitudes[centIndex] = 0.3*defaultMaxAmplitude[centIndex];
        wideRangeAmplitudes[centIndex][0] = HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[centIndex], -1, true, -0.5, 0.1);
        wideRangeAmplitudes[centIndex][1] = HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[centIndex], -1, true, 0.1, 0.5);
        wideRangeAmplitudes[centIndex][2] = HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[centIndex], -1, true, 0.5, 1.5);

        if(wideRangeAmplitudes[centIndex][0] <= 0.0) wideRangeAmplitudes[centIndex][0] = defaultAmplitudes[centIndex];
        if(wideRangeAmplitudes[centIndex][1] <= 0.0) wideRangeAmplitudes[centIndex][1] = defaultAmplitudes[centIndex];
        if(wideRangeAmplitudes[centIndex][2] <= 0.0) wideRangeAmplitudes[centIndex][2] = defaultAmplitudes[centIndex];

        double pi_funct_at_pi = 0.0;
        double pi_funct_at_K  = 0.0;
        double K_funct_at_K   = 0.0;
        double K_funct_at_pi  = 0.0;
        int    pionBin        = -1;
        int    kaonBin        = -1;
        int    protonBin      = -1;
        #ifdef _GAUSSIAN_BTOF_SIMUL_CENT_
          pionBin = histo->FindBin(pionPredict);
          kaonBin = histo->FindBin(fixed_kaon_mass_sqr);
          pi_funct_at_pi = PhysMath::gaus(pionPredict,pionSigmaInit,pionPredict);
          pi_funct_at_K  = PhysMath::gaus(pionPredict,pionSigmaInit,kaonPredict);
          pionAmpPredict =  histo->GetBinContent(pionBin)/pi_funct_at_pi;
          K_funct_at_K   = PhysMath::gaus(kaonPredict,kaonSigmaInit,kaonPredict);
          K_funct_at_pi  = PhysMath::gaus(kaonPredict,kaonSigmaInit,pionPredict);
          kaonAmpPredict =  histo->GetBinContent(kaonBin)/K_funct_at_K;
        #endif
        #ifdef _SKEW_NORMAL_BTOF_SIMUL_CENT_
          pionBin = histo->FindBin(fixed_pion_mass_sqr);
          kaonBin = histo->FindBin(fixed_kaon_mass_sqr);
          pi_funct_at_pi = PhysMath::skewNormal(pionPredict,pionSigmaInit,pionSkewInit,fixed_pion_mass_sqr);
          pi_funct_at_K  = PhysMath::skewNormal(pionPredict,pionSigmaInit,pionSkewInit,fixed_kaon_mass_sqr);
          pionAmpPredict =  histo->GetBinContent(pionBin)/pi_funct_at_pi;
          K_funct_at_K   = PhysMath::skewNormal(kaonPredict,kaonSigmaInit,kaonSkewInit,fixed_kaon_mass_sqr);
          K_funct_at_pi  = PhysMath::skewNormal(kaonPredict,kaonSigmaInit,kaonSkewInit,fixed_pion_mass_sqr);
          kaonAmpPredict =  histo->GetBinContent(kaonBin)/K_funct_at_K;
        #endif
        #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
          double pion_m2 = PhysMath::massSquared_signed_fromInvBeta(m_currentMomOfBin,pionPredict);
          double kaon_m2 = PhysMath::massSquared_signed_fromInvBeta(m_currentMomOfBin,kaonPredict);
          pionBin = histo->FindBin(pion_m2);
          kaonBin = histo->FindBin(kaon_m2);
          pi_funct_at_pi = PhysMath::m2_student_T(pionPredict,pionSigmaInit,pionSkewInit,pion_m2,m_currentMomOfBin);
          pi_funct_at_K  = PhysMath::m2_student_T(pionPredict,pionSigmaInit,pionSkewInit,kaon_m2,m_currentMomOfBin);
          pionAmpPredict = histo->GetBinContent(pionBin)/pi_funct_at_pi;
          K_funct_at_K   = PhysMath::m2_student_T(kaonPredict,kaonSigmaInit,kaonSkewInit,kaon_m2,m_currentMomOfBin);
          K_funct_at_pi  = PhysMath::m2_student_T(kaonPredict,kaonSigmaInit,kaonSkewInit,pion_m2,m_currentMomOfBin);
          kaonAmpPredict =  histo->GetBinContent(kaonBin)/K_funct_at_K;
        #endif
        #ifdef _PEARSON_IV_BTOF_SIMUL_CENT_
          pionBin = histo->FindBin(fixed_pion_mass_sqr);
          kaonBin = histo->FindBin(fixed_kaon_mass_sqr);
          pi_funct_at_pi = PhysMath::pearsonIV(pionPredict,pionSigmaInit,pionSkewInit,pionNuInit,fixed_pion_mass_sqr,normalizePearsonIV);
          pi_funct_at_K  = PhysMath::pearsonIV(pionPredict,pionSigmaInit,pionSkewInit,pionNuInit,fixed_kaon_mass_sqr,normalizePearsonIV);
          pionAmpPredict =  histo->GetBinContent(pionBin)/pi_funct_at_pi;
          K_funct_at_K   = PhysMath::pearsonIV(kaonPredict,kaonSigmaInit,kaonSkewInit,kaonNuInit,fixed_kaon_mass_sqr,normalizePearsonIV);
          K_funct_at_pi  = PhysMath::pearsonIV(kaonPredict,kaonSigmaInit,kaonSkewInit,kaonNuInit,fixed_pion_mass_sqr,normalizePearsonIV);
          kaonAmpPredict =  histo->GetBinContent(kaonBin)/K_funct_at_K;
        #endif
        #ifdef _DOUBLE_GAUS_BTOF_SIMUL_CENT_
          pionBin = histo->FindBin(fixed_pion_mass_sqr);
          kaonBin = histo->FindBin(fixed_kaon_mass_sqr);
          pi_funct_at_pi = PhysMath::doubleGauss(1.0,pionPredict,pionSigmaInit,pionSkewInit,pionNuInit,pionPredict);
          pi_funct_at_K  = PhysMath::doubleGauss(1.0,pionPredict,pionSigmaInit,pionSkewInit,pionNuInit,kaonPredict);
          pionAmpPredict =  histo->GetBinContent(pionBin)/pi_funct_at_pi;
          K_funct_at_K   = PhysMath::doubleGauss(1.0,kaonPredict,kaonSigmaInit,kaonSkewInit,kaonNuInit,pionPredict);
          K_funct_at_pi  = PhysMath::doubleGauss(1.0,kaonPredict,kaonSigmaInit,kaonSkewInit,kaonNuInit,kaonPredict);
          kaonAmpPredict =  histo->GetBinContent(kaonBin)/K_funct_at_K;
        #endif



        //   MATRIX AMPLITUDE METHOD FOR PI / K
        if(momToDoPionKaonMatrixAmpInit < m_currentMomOfBin && doPion && doKaon){
          double det = pi_funct_at_pi*K_funct_at_K - pi_funct_at_K*K_funct_at_pi;
          double pionVal = histo->GetBinContent(pionBin);
          double kaonVal = histo->GetBinContent(kaonBin);
          #ifdef _ZFITTER_DEBUG_
            cout << "Using matrix amplitude initialization for pion and kaons: " << endl;
            cout << "   Initial Amplitudes   Pi: " << pionAmpPredict << "  K: " << kaonAmpPredict << endl;
          #endif
          double newPionAmp = -1;
          double newKaonAmp = -1;
          if(det > 0){
            newPionAmp = (K_funct_at_K*pionVal - pi_funct_at_K*kaonVal)/det;
            newKaonAmp = (pi_funct_at_pi*kaonVal - K_funct_at_pi*pionVal)/det;
          }
          if(newPionAmp > 0.0 && newKaonAmp > 0.0){
            pionAmpPredict = newPionAmp;
            kaonAmpPredict = newKaonAmp;
            #ifdef _ZFITTER_DEBUG_
              cout << "   Matrix   Amplitudes   Pi: " << pionAmpPredict << "  K: " << kaonAmpPredict << "    (if these are the same as above, matrix failed)" <<  endl;
            #endif
          }
        }
        double proton_m2 = PhysMath::massSquared_signed_fromInvBeta(m_currentMomOfBin,protonPredict);
        protonBin = histo->FindBin(proton_m2);
        #ifdef _ZFITTER_DEBUG_
          cout << " Proton M2: " << proton_m2 << " from invBeta " << protonPredict << "  at bin " << protonBin << " with center " << histo->GetBinCenter(protonBin) << endl;
        #endif
        if(protonFunctionVersion == 0)       protonAmpPredict =  histo->GetBinContent(protonBin)/PhysMath::gaus(protonPredict,protonSigmaInit,protonPredict);
        else if(protonFunctionVersion == 1)  protonAmpPredict =  histo->GetBinContent(protonBin)/PhysMath::skewNormal(protonPredict,protonSigmaInit,protonSkewInit,protonPredict);
        else if(protonFunctionVersion == 2){
          double value = PhysMath::m2_student_T(protonPredict,protonSigmaInit,protonSkewInit,proton_m2,m_currentMomOfBin);
          protonAmpPredict =  histo->GetBinContent(protonBin)/value;
        }else if(protonFunctionVersion == 3)  protonAmpPredict =  histo->GetBinContent(protonBin)/PhysMath::pearsonIV(protonPredict,protonSigmaInit,protonSkewInit,protonNuInit,protonPredict,normalizePearsonIV);
        else if(protonFunctionVersion == 4)   protonAmpPredict =  histo->GetBinContent(protonBin)/PhysMath::doubleGauss(1.0,protonPredict,protonSigmaInit,protonSkewInit,protonNuInit,protonPredict);
        if(pionAmpPredict <= 0.0   || !std::isfinite(pionAmpPredict) ){
          pionAmpPredict = 0.0;
          //cout << "ERROR: Pion initialization failed... Not fitting..." << pionAmpPredict << " " << pionPredict << " " << pionSigmaInit << " " << pionSkewInit << " " << pionNuInit << endl;
          //doPion = false;
        }
        if(kaonAmpPredict <= 0.0   || !std::isfinite(kaonAmpPredict) ){
          kaonAmpPredict = 0.0;
          //cout << "ERROR: Kaon initialization failed... Not fitting..." << kaonAmpPredict << " " << kaonPredict << " " << kaonSigmaInit << " " << kaonSkewInit << " " << kaonNuInit << endl;
          //doKaon = false;
        }
        if(protonAmpPredict <= 0.0 || !std::isfinite(protonAmpPredict) ){
          protonAmpPredict = 0.0;
          //cout << "ERROR: Proton initialization failed... Not fitting..." << protonAmpPredict << " " << protonPredict << " " << protonSigmaInit << " " << protonSkewInit << " " << protonNuInit << endl;
          //doProton = false;
        }

        amplitudes[centIndex][0] = pionAmpPredict;
        amplitudes[centIndex][1] = kaonAmpPredict;
        amplitudes[centIndex][2] = protonAmpPredict;

        initialPoints[centIndex][0]->SetPoint(0,pion_m2,histo->GetBinContent(pionBin));
        initialPoints[centIndex][1]->SetPoint(0,kaon_m2,histo->GetBinContent(kaonBin));
        initialPoints[centIndex][2]->SetPoint(0,proton_m2,histo->GetBinContent(protonBin));



      }// end cent loop for amplitudes


      bool pionAmpsAllZero     = true;
      bool kaonAmpsAllZero     = true;
      bool protonAmpsAllZero   = true;
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        if(amplitudes[centIndex][0] > 0.0) pionAmpsAllZero     = false;
        if(amplitudes[centIndex][1] > 0.0) kaonAmpsAllZero     = false;
        if(amplitudes[centIndex][2] > 0.0) protonAmpsAllZero   = false;
      }

      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        if(!pionAmpsAllZero && amplitudes[centIndex][0] <= 0.0)   amplitudes[centIndex][0] = wideRangeAmplitudes[centIndex][0];
        if(!kaonAmpsAllZero && amplitudes[centIndex][1] <= 0.0)   amplitudes[centIndex][1] = wideRangeAmplitudes[centIndex][1];
        if(!protonAmpsAllZero && amplitudes[centIndex][2] <= 0.0) amplitudes[centIndex][2] = wideRangeAmplitudes[centIndex][2];
      }



      if(pionAmpsAllZero){
        #ifdef _ZFITTER_DEBUG_
          cout << "WARNGING: Pion Amplitudes all Zero!" << endl;
        #endif
        doPion     = false;

      }
      if(kaonAmpsAllZero){
        #ifdef _ZFITTER_DEBUG_
          cout << "WARNGING: Kaon Amplitudes all Zero!" << endl;
        #endif
        doKaon     = false;

      }
      if(protonAmpsAllZero){
        #ifdef _ZFITTER_DEBUG_
          cout << "WARNGING: Proton Amplitudes all Zero!" << endl;
        #endif
        doProton   = false;
      }




      if(m_currentPartIndex == 0 && m_currentMomOfBin < upToMomToIgnoreKaons) doKaon = false;
      if(m_currentPartIndex != 2 && m_currentMomOfBin < momToIntroduceProton) doProton = false;
      if(m_currentPartIndex == 1 && m_currentMomOfBin < highMomToIgnorePionsInKaonSpace){
        doPion   = false;
        doProton = false;
      }
      if(m_currentPartIndex == 2 && m_currentMomOfBin < upToMomToIgnorePionKaonsForProtons){
        doKaon = false;
        doPion = false;
      }

      if(!doPion && !doKaon && !doProton){
        cout << "ERROR: Nothing to fit here..." << endl;
        HistogramUtilities::printHistoInfo(m_currentHistoToFit);
      }


      //##############   CREATE THE FIT FUNCTION  ##################
      string paramNames[12] = {"#mu_{#pi}","#sigma_{#pi}","Ignore_{#pi}","Ignore_{#pi}","#mu_{K}","#sigma_{K}","Ignore_{K}","Ignore_{K}","#mu_{p}","#sigma_{p}","Ignore_{p}","Ignore_{p}"};
      if(pionFunctionVersion == 1){          paramNames[2]="#gamma_{#pi}";}
      if(pionFunctionVersion == 2){          paramNames[2]="#nu_{#pi}";}
      if(pionFunctionVersion == 3){          paramNames[2]="m_{#pi}"; paramNames[3]="#nu_{#pi}";}
      if(pionFunctionVersion == 4){          paramNames[2]="r_{#pi}"; paramNames[3]="#sigma_{2,#pi}";    }
      if(kaonFunctionVersion == 1){          paramNames[6]="#gamma_{K}";}
      if(kaonFunctionVersion == 2){          paramNames[6]="#nu_{K}";}
      if(kaonFunctionVersion == 3){          paramNames[6]="m_{K}"; paramNames[7]="#nu_{K}";}
      if(kaonFunctionVersion == 4){          paramNames[6]="r_{K}"; paramNames[7]="#sigma_{2,K}";        }
      if(protonFunctionVersion == 1){        paramNames[10]="#gamma_{p}";}
      if(protonFunctionVersion == 2){        paramNames[10]="#nu_{p}";}
      if(protonFunctionVersion == 3){        paramNames[10]="m_{p}"; paramNames[11]="#nu_{p}";}
      if(protonFunctionVersion == 4){        paramNames[10]="r_{p}"; paramNames[11]="#sigma_{2,p}";          }

      double parameters_simul_cent[70];
      double parameters_simul_cent_errors[70];
      parameters_simul_cent[0]  = pionPredict;
      parameters_simul_cent[1]  = pionSigmaInit;
      parameters_simul_cent[2]  = pionSkewInit;
      parameters_simul_cent[3]  = pionNuInit;
      parameters_simul_cent[4]  = kaonPredict;
      parameters_simul_cent[5]  = kaonSigmaInit;
      parameters_simul_cent[6]  = kaonSkewInit;
      parameters_simul_cent[7]  = kaonNuInit;
      parameters_simul_cent[8]  = protonPredict;
      parameters_simul_cent[9]  = protonSigmaInit;
      parameters_simul_cent[10] = protonSkewInit;
      parameters_simul_cent[11] = protonNuInit;

      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        parameters_simul_cent[12 + centIndex*3] = amplitudes[centIndex][0];
        parameters_simul_cent[13 + centIndex*3] = amplitudes[centIndex][1];
        parameters_simul_cent[14 + centIndex*3] = amplitudes[centIndex][2];
      }

      #ifdef _ZFITTER_DEBUG_

        cout << "doFitPion  : " << (doPion ? "TRUE" : "FALSE") << endl;
        cout << "doFitKaon  : " << (doKaon ? "TRUE" : "FALSE") << endl;
        cout << "doFitProton  : " << (doProton ? "TRUE" : "FALSE") << endl;

        cout << "###########  PARAMETERS FOR THE FIT ###########" << endl;
        cout << paramNames[0] << "  :  " << parameters_simul_cent[0] << endl;
        cout << paramNames[1] << "  :  " << parameters_simul_cent[1] << endl;
        cout << paramNames[2] << "  :  " << parameters_simul_cent[2] << endl;
        cout << paramNames[3] << "  :  " << parameters_simul_cent[3] << endl;
        cout << paramNames[4] << "  :  " << parameters_simul_cent[4] << endl;
        cout << paramNames[5] << "  :  " << parameters_simul_cent[5] << endl;
        cout << paramNames[6] << "  :  " << parameters_simul_cent[6] << endl;
        cout << paramNames[7] << "  :  " << parameters_simul_cent[7] << endl;
        cout << paramNames[8] << "  :  " << parameters_simul_cent[8] << endl;
        cout << paramNames[9] << "  :  " << parameters_simul_cent[9] << endl;
        cout << paramNames[10] << "  :  " << parameters_simul_cent[10] << endl;
        cout << paramNames[11] << "  :  " << parameters_simul_cent[11] << endl;
        for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
          cout << " A_Pi_" << centIndex << "  : " << parameters_simul_cent[12 + centIndex*3] << endl;
          cout << " A_K_" << centIndex << "   : " << parameters_simul_cent[13 + centIndex*3] << endl;
          cout << " A_P_" << centIndex << "   : " << parameters_simul_cent[14 + centIndex*3] << endl;
        }
        cout << "################################" << endl;

      #endif


      if(!doPion && !doKaon && !doProton){
        cout << "WARNING: Not Fitting Any Particle... Continuing to next mTm0 Bin..." << endl;
      }



      //Assign Variables And Ranges To Minimizer
      // MINIMIZER SETUP
      ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer();
      ROOT::Math::Functor* chiSqrdFunctor = new ROOT::Math::Functor(this,&ZFitter::simultaneous_centrality_ZTOF_chisqr,m_numCentralities*3 + 12);
      minimizer->SetFunction(*chiSqrdFunctor);
      #ifdef _ZFITTER_DEBUG_
        minimizer->SetPrintLevel(1);
      #else
        minimizer->SetPrintLevel(0);
      #endif
      minimizer->SetMaxIterations(1000);
      minimizer->SetMaxFunctionCalls(20000);
      if(doHesse) minimizer->SetValidError(true);

      //minimizer->Clear();
      for(int paramIndex = 0; paramIndex < 12; paramIndex++){
        minimizer->SetVariable(paramIndex, paramNames[paramIndex] , parameters_simul_cent[paramIndex], 0.01);
      }

      // funciton versions:  0 = gaussian    1 = skewNormal    2 = student-t  3 = pearsonIV
      #ifdef _GAUSSIAN_BTOF_SIMUL_CENT_
        //PION
        minimizer->SetVariable(0, paramNames[0] ,pionPredict, 0.001);
        minimizer->SetVariableLimits(0, pionPredict - 2.0*pionSigmaInit, pionPredict + 2.0*pionSigmaInit);
        minimizer->SetVariable(1, paramNames[1] ,pionSigmaInit, 0.01);
        minimizer->SetVariableLimits(1, 0.0, 50.0*pionSigmaInit);
        minimizer->SetFixedVariable(2, paramNames[2] ,0);
        minimizer->SetFixedVariable(3, paramNames[3] ,0);
        //KAON
        minimizer->SetVariable(4, paramNames[4] ,kaonPredict, 0.001);
        minimizer->SetVariableLimits(4, kaonPredict - 2.0*kaonSigmaInit, kaonPredict + 2.0*kaonSigmaInit);
        minimizer->SetVariable(5, paramNames[5] ,kaonSigmaInit, 0.01);
        minimizer->SetVariableLimits(5, 0.0, 50.0*kaonSigmaInit);
        minimizer->SetFixedVariable(6, paramNames[6] , 0);
        minimizer->SetFixedVariable(7, paramNames[7] , 0);
        //PROTON
        minimizer->SetVariable(8, paramNames[8] ,protonPredict, 0.001);
        minimizer->SetVariableLimits(8, protonPredict - 2.0*protonSigmaInit, protonPredict + 2.0*protonSigmaInit);
        minimizer->SetVariable(9, paramNames[9] ,protonSigmaInit, 0.01);
        minimizer->SetVariableLimits(9, 0.0, 50.0*protonSigmaInit);
        minimizer->SetFixedVariable(10, paramNames[10] , 0);
        minimizer->SetFixedVariable(11, paramNames[11] , 0);
      #endif

      #ifdef _SKEW_NORMAL_BTOF_SIMUL_CENT_
        //PION
        minimizer->SetVariable(0, paramNames[0] ,pionPredict, 0.001);
        minimizer->SetVariableLimits(0, pionPredict - 2.0*pionSigmaInit, pionPredict + 2.0*pionSigmaInit);
        minimizer->SetVariable(1, paramNames[1] ,pionSigmaInit, 0.01);
        minimizer->SetVariableLimits(1, 0.0, 50.0*pionSigmaInit);
        minimizer->SetVariable(2, paramNames[2] , pionSkewInit , 0.01);
        minimizer->SetVariableLimits(2, -50, 50);
        minimizer->SetFixedVariable(3, paramNames[3] , 0);
        //KAON
        minimizer->SetVariable(4, paramNames[4] ,kaonPredict, 0.001);
        minimizer->SetVariableLimits(4, kaonPredict - 2.0*kaonSigmaInit, kaonPredict + 2.0*kaonSigmaInit);
        minimizer->SetVariable(5, paramNames[5] ,kaonSigmaInit, 0.01);
        minimizer->SetVariableLimits(5, 0.0, 50.0*kaonSigmaInit);
        minimizer->SetVariable(6, paramNames[6] , kaonSkewInit , 0.01);
        minimizer->SetVariableLimits(6, -50, 50);
        minimizer->SetFixedVariable(7, paramNames[7] , 0);
        //PROTON
        minimizer->SetVariable(8, paramNames[8] ,protonPredict, 0.001);
        minimizer->SetVariableLimits(8, protonPredict - 2.0*protonSigmaInit, protonPredict + 2.0*protonSigmaInit);
        minimizer->SetVariable(9, paramNames[9] ,protonSigmaInit, 0.01);
        minimizer->SetVariableLimits(9, 0.0, 50.0*protonSigmaInit);
        minimizer->SetVariable(10, paramNames[10] , protonSkewInit , 0.01);
        minimizer->SetVariableLimits(10, -50, 50);
        minimizer->SetFixedVariable(11, paramNames[11] , 0);
      #endif

      #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
        //### PION MEAN
        minimizer->SetVariable(0, paramNames[0] ,pionPredict, 0.001);
        //minimizer->SetVariableLimits(0, pionPredict - 2.0*pionSigmaInit, pionPredict + 2.0*pionSigmaInit); //m2
        if(hasBeenExtrapolatedMatrix[0][m_currentPlusMinusIndex][1]) minimizer->SetFixedVariable(0,  paramNames[0], pionPredict);
        else                                                         minimizer->SetVariableLimits(0, pionPredict - 0.1, pionPredict + 0.1);
        //### PION SIGMA
        minimizer->SetVariable(1, paramNames[1] ,pionSigmaInit, 0.001);
        if(hasBeenExtrapolatedMatrix[0][m_currentPlusMinusIndex][2]) minimizer->SetFixedVariable(1, paramNames[1], pionSigmaInit);
        else{
          if(previousCentralityIntegratedYields_Pion > centIntThreshold){
            minimizer->SetVariableLimits(1, 0.5*pionSigmaInit, 2.0*pionSigmaInit);
          }else{
            minimizer->SetVariableLimits(1, 0.00001, 0.01);
          }
        }
        //### PION NU
        minimizer->SetVariable(2, paramNames[2] , pionSkewInit , 0.01);
        bool makePionGaussian = false;
        if(pionSkewInit > 200 && centralityIntegratedYieldsWereAboveThreshold_Pion){
          pionSkewInit = 250;
          makePionGaussian = true;
        }
        if(hasBeenExtrapolatedMatrix[0][m_currentPlusMinusIndex][3] || makePionGaussian) minimizer->SetFixedVariable(2, paramNames[2], pionSkewInit);
        else{
          if(m_currentPartIndex == 0)   minimizer->SetVariableLimits(2, minStudentTNuVal_POI, 250);
          else                          minimizer->SetVariableLimits(2, minStudentTNuVal, 250);
        }

        minimizer->SetFixedVariable(3, paramNames[3] , 0); // ignored parameter

        if(m_currentPartCharge < 0){
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[0][1]) minimizer->SetFixedVariable(0, paramNames[0] , pionPredict);
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[0][2]) minimizer->SetFixedVariable(1, paramNames[1] , pionSigmaInit);
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[0][3]) minimizer->SetFixedVariable(2, paramNames[2] , pionSkewInit);
          //if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[2][4]) minimizer->SetFixedVariable(11, paramNames[11] , protonNuInit);
        }



        //KAON
        minimizer->SetVariable(4, paramNames[4] ,kaonPredict, 0.001);
        if(hasBeenExtrapolatedMatrix[1][m_currentPlusMinusIndex][1]) minimizer->SetFixedVariable(4, paramNames[4] ,kaonPredict);
        else                                                         minimizer->SetVariableLimits(4, kaonPredict - 0.1, kaonPredict + 0.1);
        minimizer->SetVariable(5, paramNames[5] ,kaonSigmaInit, 0.01);
        if(hasBeenExtrapolatedMatrix[1][m_currentPlusMinusIndex][2]) minimizer->SetFixedVariable(5, paramNames[5] ,kaonSigmaInit);
        else{
          if(previousCentralityIntegratedYields_Kaon > centIntThreshold){
            minimizer->SetVariableLimits(5, 0.5*kaonSigmaInit, 2.0*kaonSigmaInit);
          }else{
            minimizer->SetVariableLimits(5, 0.00001, 0.01);
          }
        }
        minimizer->SetVariable(6, paramNames[6] , kaonSkewInit , 0.01);
        bool makeKaonGaussian = false;
        if(kaonSkewInit > 200 && centralityIntegratedYieldsWereAboveThreshold_Kaon){
          kaonSkewInit = 250;
          makeKaonGaussian = true;
        }
        if(hasBeenExtrapolatedMatrix[1][m_currentPlusMinusIndex][3] || makeKaonGaussian) minimizer->SetFixedVariable(6, paramNames[6] , kaonSkewInit);
        else{
          if(m_currentPartIndex == 1)   minimizer->SetVariableLimits(6, minStudentTNuVal_POI, 250);
          else                          minimizer->SetVariableLimits(6, minStudentTNuVal, 250);
        }
        minimizer->SetFixedVariable(7, paramNames[7] , 0);// ignored parameter

        /*if(m_currentMomOfBin > momToFixKaonStudentTNu){
          double storedKaonNu = m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][1][0][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin - 1);
          if(fixedKaonStudentTNuValue < 1.0 && storedKaonNu > 1.0){
            fixedKaonStudentTNuValue = storedKaonNu;
          }
          if(fixedKaonStudentTNuValue >= 1.0){
            minimizer->SetVariable(6, paramNames[6] , fixedKaonStudentTNuValue , 0.01);
            minimizer->FixVariable(6);
          }
        }*/

        if(m_currentPartCharge < 0){
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[1][1]) minimizer->SetFixedVariable(4, paramNames[4] , kaonPredict);
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[1][2]) minimizer->SetFixedVariable(5, paramNames[5] , kaonSigmaInit);
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[1][3]) minimizer->SetFixedVariable(6, paramNames[6] , kaonSkewInit);
        }


        //PROTON
        minimizer->SetVariable(8, paramNames[8] ,protonPredict, 0.001);
        if(hasBeenExtrapolatedMatrix[2][m_currentPlusMinusIndex][1]) minimizer->SetFixedVariable(8, paramNames[8] ,protonPredict);
        else                                                         minimizer->SetVariableLimits(8, protonPredict - 0.1, protonPredict + 0.1);
        minimizer->SetVariable(9, paramNames[9] ,protonSigmaInit, 0.001);
        if(hasBeenExtrapolatedMatrix[2][m_currentPlusMinusIndex][2]) minimizer->SetFixedVariable(9, paramNames[9] ,protonSigmaInit);
        else{
          if(previousCentralityIntegratedYields_Proton > centIntThreshold){
            minimizer->SetVariableLimits(9, 0.5*protonSigmaInit, 2.0*protonSigmaInit);
          }else{
            minimizer->SetVariableLimits(9, 0.00001, 0.01);
          }
        }
        minimizer->SetVariable(10, paramNames[10] , protonSkewInit, 0.01);
        bool makeProtonGaussian = false;
        if(protonSkewInit > 200 && centralityIntegratedYieldsWereAboveThreshold_Proton){
          protonSkewInit = 250;
          makeProtonGaussian = true;
        }
        if(hasBeenExtrapolatedMatrix[2][m_currentPlusMinusIndex][3] || makeProtonGaussian) minimizer->SetFixedVariable(10, paramNames[10] , protonSkewInit);
        else{
          if(m_currentPartIndex == 2)   minimizer->SetVariableLimits(10, minStudentTNuVal_POI, 250);
          else                          minimizer->SetVariableLimits(10, minStudentTNuVal, 250);
        }
        minimizer->SetFixedVariable(11, paramNames[11] , 0); // ignored parameter

        if(m_currentPartCharge < 0){
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[2][1]) minimizer->SetFixedVariable(8,  paramNames[8] , protonPredict);
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[2][2]) minimizer->SetFixedVariable(9,  paramNames[9] , protonSigmaInit);
          if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[2][3]) minimizer->SetFixedVariable(10, paramNames[10] , protonSkewInit);
          //if(m_currentMomOfBin > fixToPositiveValuesRuleMatrix_reduced[2][4]) minimizer->SetFixedVariable(11, paramNames[11] , protonNuInit);
        }


      #endif

      #ifdef _PEARSON_IV_BTOF_SIMUL_CENT_
        //PION
        minimizer->SetVariable(0, paramNames[0] ,pionPredict, 0.001);
        minimizer->SetVariableLimits(0, pionPredict - 0.1, pionPredict + 0.1);
        minimizer->SetVariable(1, paramNames[1] ,pionSigmaInit, 0.01);
        minimizer->SetVariableLimits(1, 0.0, 10.0*pionSigmaInit);
        minimizer->SetVariable(2, paramNames[2] , pionSkewInit , 0.01);
        minimizer->SetVariableLimits(2, 0.01, 350);
        minimizer->SetVariable(3, paramNames[3] , 0 , 0.01);
        minimizer->SetVariableLimits(3, -350, 350);
        //KAON
        minimizer->SetVariable(4, paramNames[4] ,kaonPredict, 0.001);
        minimizer->SetVariableLimits(4, kaonPredict - 0.1, kaonPredict + 0.1);
        minimizer->SetVariable(5, paramNames[5] ,kaonSigmaInit, 0.01);
        minimizer->SetVariableLimits(5, 0.0, 10.0*kaonSigmaInit);
        minimizer->SetVariable(6, paramNames[6] , kaonSkewInit , 0.01);
        minimizer->SetVariableLimits(6, 0.01, 350);
        minimizer->SetVariable(7, paramNames[7] , 0 , 0.01);
        minimizer->SetVariableLimits(7, -350, 350);
        //PROTON
        minimizer->SetVariable(8, paramNames[8] ,protonPredict, 0.001);
        minimizer->SetVariableLimits(8, protonPredict - 0.1, protonPredict + 0.1);
        minimizer->SetVariable(9, paramNames[9] ,protonSigmaInit, 0.01);
        minimizer->SetVariableLimits(9, 0.0, 10.0*protonSigmaInit);
        minimizer->SetVariable(10, paramNames[10] , protonSkewInit , 0.01);
        minimizer->SetVariableLimits(10, 0.01, 350);
        minimizer->SetVariable(11, paramNames[11] , 0 , 0.01);
        minimizer->SetVariableLimits(11, -350, 350);
      #endif

      #ifdef _DOUBLE_GAUS_BTOF_SIMUL_CENT_
        //PION
        minimizer->SetVariable(0, paramNames[0] ,pionPredict, 0.001);
        minimizer->SetVariableLimits(0, pionPredict - 0.1, pionPredict + 0.1);
        minimizer->SetVariable(1, paramNames[1] ,pionSigmaInit, 0.01);
        minimizer->SetVariableLimits(1, 0.0, 10.0*pionSigmaInit);
        if(m_currentPartIndex == 0){
          minimizer->SetVariable(2, paramNames[2] , pionSkewInit , 0.01);
          minimizer->SetVariableLimits(2, 0.0, 1.0);
        }else{
          minimizer->SetVariable(2, paramNames[2] , pionSkewInit , 0.01);
          minimizer->SetVariableLimits(2, -10.0, 10.0);
        }
        minimizer->SetVariable(3, paramNames[3] , pionNuInit, 0.01);
        minimizer->SetVariableLimits(3, 0.0, 10.0*pionSigmaInit);
        //KAON
        minimizer->SetVariable(4, paramNames[4] ,kaonPredict, 0.001);
        minimizer->SetVariableLimits(4, kaonPredict - 0.1, kaonPredict + 0.1);
        minimizer->SetVariable(5, paramNames[5] ,kaonSigmaInit, 0.01);
        minimizer->SetVariableLimits(5, 0.0, 10.0*kaonSigmaInit);
        if(m_currentPartIndex == 1){
          minimizer->SetVariable(6, paramNames[6] , pionSkewInit , 0.01);
          minimizer->SetVariableLimits(6, 0.0, 1.0);
        }else{
          minimizer->SetVariable(6, paramNames[6] , pionSkewInit , 0.01);
          minimizer->SetVariableLimits(6, -10.0, 10.0);
        }
        minimizer->SetVariable(7, paramNames[7] , kaonNuInit, 0.01);
        minimizer->SetVariableLimits(7, 0.0, 10.0*kaonSigmaInit);
        //PROTON
        minimizer->SetVariable(8, paramNames[8] ,protonPredict, 0.001);
        minimizer->SetVariableLimits(8, protonPredict - 0.1, protonPredict + 0.1);
        minimizer->SetVariable(9, paramNames[9] ,protonSigmaInit, 0.01);
        minimizer->SetVariableLimits(9, 0.0, 10.0*protonSigmaInit);
        if(m_currentPartIndex == 2){
          minimizer->SetVariable(10, paramNames[10] , pionSkewInit , 0.01);
          minimizer->SetVariableLimits(10, 0.0, 1.0);
        }else{
          minimizer->SetVariable(10, paramNames[10] , pionSkewInit , 0.01);
          minimizer->SetVariableLimits(10, -10.0, 10.0);
        }
        minimizer->SetVariable(11, paramNames[11] , protonNuInit, 0.01);
        minimizer->SetVariableLimits(11, 0.0, 10.0*protonSigmaInit);
      #endif

      //###################          SET AMPLITUDE PARAMETERS      ##############################################

      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->SetVariable(12 + centIndex*3, Form("N_{#pi}_%d",centIndex), parameters_simul_cent[12 + centIndex*3], 0.01*parameters_simul_cent[12 + centIndex*3]);
        if(parameters_simul_cent[12 + centIndex*3] == 0.0)  minimizer->FixVariable(12 + centIndex*3);
        else                                                minimizer->SetVariableLimits(12 + centIndex*3, 0.0, 1.25*defaultMaxAmplitude[centIndex]);
        minimizer->SetVariable(13 + centIndex*3, Form("N_{K}_%d",centIndex), parameters_simul_cent[13 + centIndex*3], 0.01*parameters_simul_cent[13 + centIndex*3]);
        if(parameters_simul_cent[13 + centIndex*3] == 0.0)  minimizer->FixVariable(13 + centIndex*3);
        else                                                minimizer->SetVariableLimits(13 + centIndex*3, 0.0, 1.25*defaultMaxAmplitude[centIndex]);
        minimizer->SetVariable(14 + centIndex*3, Form("N_{p}_%d",centIndex), parameters_simul_cent[14 + centIndex*3], 0.01*parameters_simul_cent[14 + centIndex*3]);
        if(parameters_simul_cent[14 + centIndex*3] == 0.0)  minimizer->FixVariable(14 + centIndex*3);
        else                                                minimizer->SetVariableLimits(14 + centIndex*3, 0.0, 1.25*defaultMaxAmplitude[centIndex]);
      }




      if(!doPion){
        minimizer->SetVariable(0, paramNames[0] , -3  , 0.001);
        minimizer->SetVariableLimits(0, -3.01, -2.99);
        minimizer->SetVariable(1, paramNames[1] , 1.0, 0.01);
        minimizer->SetVariableLimits(1, 0.99, 1.01);
        if(pionFunctionVersion == 4){
          minimizer->SetVariable(2, paramNames[2] , 0.0, 0.01);
          minimizer->SetVariable(3, paramNames[3] , 1.0 , 0.01);
        }else{
          minimizer->SetVariable(2, paramNames[2] , 3.0, 0.01);
          minimizer->SetVariableLimits(2, 2.9, 3.1);
          minimizer->SetVariable(3, paramNames[3] , 0.0 , 0.01);
          minimizer->SetVariableLimits(3, -0.1, 0.1);
        }
        minimizer->FixVariable(0);
        minimizer->FixVariable(1);
        minimizer->FixVariable(2);
        minimizer->FixVariable(3);
        for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
          minimizer->SetVariable(12 + centIndex*3, Form("N_{#pi}_%d",centIndex),0.0, 0.01);
          minimizer->SetVariableLimits(12 + centIndex*3, -0.01, 0.01);
          minimizer->FixVariable(12 + centIndex*3);
        }
        #ifdef _ZFITTER_DEBUG_
          cout << " Killing Pion" << endl;
        #endif
      }

      if(!doKaon){
        minimizer->SetVariable(4, paramNames[4] , -2  , 0.001);
        minimizer->SetVariableLimits(4, -2.01, -1.99);
        minimizer->SetVariable(5, paramNames[5] , 1.0, 0.01);
        minimizer->SetVariableLimits(5, 0.99, 1.01);
        if(kaonFunctionVersion == 4){
          minimizer->SetVariable(6, paramNames[6] , 0.0, 0.01);
          minimizer->SetVariable(7, paramNames[7] , 1.0 , 0.01);
        }else{
          minimizer->SetVariable(6, paramNames[6] , 3.0, 0.01);
          minimizer->SetVariableLimits(6, 2.99, 3.01);
          minimizer->SetVariable(7, paramNames[7] , 0  , 0.01);
          minimizer->SetVariableLimits(7, -0.01, 0.01);
        }
        minimizer->FixVariable(4);
        minimizer->FixVariable(5);
        minimizer->FixVariable(6);
        minimizer->FixVariable(7);
        for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
          minimizer->SetVariable(13 + centIndex*3, Form("N_{K}_%d",centIndex), 0.0, 0.01);
          minimizer->SetVariableLimits(13 + centIndex*3, -0.001, 0.001);
          minimizer->FixVariable(13 + centIndex*3);
        }
        #ifdef _ZFITTER_DEBUG_
          cout << " Killing Kaon" << endl;
        #endif
      }

      if(!doProton){
        minimizer->SetVariable(8, paramNames[8] , -3 , 0.001);
        minimizer->SetVariableLimits(8, -3.01, -2.99);
        minimizer->SetVariable(9, paramNames[9] , 1.0, 0.01);
        minimizer->SetVariableLimits(9, 0.99, 1.01);

        if(protonFunctionVersion == 4){
          minimizer->SetVariable(10, paramNames[10] , 0.0, 0.01);
          minimizer->SetVariable(11, paramNames[11] , 1.0 , 0.01);
        }else{
          minimizer->SetVariable(10, paramNames[10] , 3.0, 0.01);
          minimizer->SetVariableLimits(10, 2.99, 3.01);
          minimizer->SetVariable(11, paramNames[11] , 0  , 0.01);
          minimizer->SetVariableLimits(9, -0.01, 0.01);

        }
        minimizer->FixVariable(8);
        minimizer->FixVariable(9);
        minimizer->FixVariable(10);
        minimizer->FixVariable(11);
        for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
          minimizer->SetVariable(14 + centIndex*3, Form("N_{p}_%d",centIndex), 0.0, 0.01);
          minimizer->SetVariableLimits(14 + centIndex*3, -0.01, 0.01);
          minimizer->FixVariable(14 + centIndex*3);
        }
        #ifdef _ZFITTER_DEBUG_
          cout << " Killing Proton" << endl;
        #endif
      }

      if(!doPion && !doKaon && !doProton){
        cout << " WARNING: Nothing to fit... something wrong? : " << endl;
        HistogramUtilities::printHistoInfo(m_currentHistoToFit);
        continue;
      }



      //####### CREATE INITIAL SEED FUNCTION FOR DEBUGGING  ########
      TF1* initialFunct[16];
      TF1* initialFunct_invBeta[16];
      for(int centIndex = 0; centIndex < 16; centIndex++){
        initialFunct[centIndex] = NULL;
        initialFunct_invBeta[centIndex] = NULL;
      }
      if(m_drawInitialSeedsToFits){
        for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
          double pion_amp   = parameters_simul_cent[12 + centIndex*3];
          double kaon_amp   = parameters_simul_cent[13 + centIndex*3];
          double proton_amp = parameters_simul_cent[14 + centIndex*3];
          if(!doPion)   pion_amp     = 0.0;
          if(!doKaon)   kaon_amp     = 0.0;
          if(!doProton) proton_amp   = 0.0;

          #ifdef _GAUSSIAN_BTOF_SIMUL_CENT_
            initialFunct[centIndex] = PhysMath::getThreeGausFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_init",
              m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
              parameters_simul_cent,pion_amp,kaon_amp,proton_amp);
          #endif
          #ifdef _SKEW_NORMAL_BTOF_SIMUL_CENT_
            initialFunct[centIndex] = PhysMath::getThreeSkewFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_init",
              m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
              parameters_simul_cent,pion_amp,kaon_amp,proton_amp);
          #endif
          #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
            initialFunct[centIndex] = PhysMath::getThreeM2StudentTFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_init",
              m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
              parameters_simul_cent,pion_amp,kaon_amp,proton_amp,m_currentMomOfBin);
            initialFunct_invBeta[centIndex] = PhysMath::getThreeStudentTFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_init_invBeta",
              m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
              parameters_simul_cent,pion_amp,kaon_amp,proton_amp);
            double value = PhysMath::m2_student_T(protonPredict,protonSigmaInit,protonSkewInit,fixed_proton_mass_sqr,m_currentMomOfBin);
            double protonAmpPredict =  m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(fixed_proton_mass_sqr))/value;
            double functVal = initialFunct[centIndex]->Eval(fixed_proton_mass_sqr);
            cout << "Init Check: " << functVal << " =? " << protonAmpPredict*value << endl;
            initialFunct_invBeta[centIndex]->SetLineColor(kYellow);
            initialFunct_invBeta[centIndex]->SetRange(0.1,15);
          #endif
          #ifdef _PEARSON_IV_BTOF_SIMUL_CENT_
            initialFunct[centIndex] = PhysMath::getThreePearsonIVFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_init",
              m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
              parameters_simul_cent, pion_amp,kaon_amp,proton_amp);
          #endif
          #ifdef _DOUBLE_GAUS_BTOF_SIMUL_CENT_
            initialFunct[centIndex] = PhysMath::getThreeDoubleGaussFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_init",
              m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
              parameters_simul_cent, pion_amp,kaon_amp,proton_amp);
          #endif
          initialFunct[centIndex]->SetLineColor(kYellow);
          initialFunct[centIndex]->SetRange(-0.3,1.5);

        } //cent loop
      }


      //################   REDUCE HISTOGRAM RANGES BASED ON WHICH PARTICLES FIT ################
      //setup fit ranges
      m_currentLowFitRange = -0.25;
      m_currentHighFitRange = 1.5;
      for(int centIndex = 0;centIndex < m_numCentralities; centIndex++){
        m_currentLowFitBin[centIndex] = 1;
        m_currentHighFitBin[centIndex] = m_currentHistosToFit_ByCent[centIndex]->GetNbinsX();
      }

      //double lowFitRange = lowRange;
      //double highFitRange = highRange;
      //Fix the fit ranges for PION
      if(a_partIndexSpace == 0 && doPion){
        if(m_currentMomOfBin < upToMomToIgnoreMuons)      m_currentLowFitRange = m2ForMuonCut;
        else if(m_currentMomOfBin < highMomForMuonIgnore) m_currentLowFitRange = lowerPionFitRangeForPionLinFunct->Eval(m_currentMomOfBin);
        else                                              m_currentLowFitRange = basicPionM2LowCut;
        if(doProton){
          if(m_currentMomOfBin < lowerProtonRangeMom){
            m_currentHighFitRange = highFitRangeForProtonsInKaonSpace;
          }else if(m_currentMomOfBin < upperProtonRangeMom){
            m_currentHighFitRange = upperProtonFitRangeForKaonLinFunct->Eval(m_currentMomOfBin);
          }else{
            m_currentHighFitRange = 1.1;
          }
        }else{
          m_currentHighFitRange = 0.4;
        }
        if(!doKaon){
          m_currentHighFitRange = highFitRangeWhenIngoringKaonsInPion;
        }
      }

      //Fix the fit ranges for KAON
      if(a_partIndexSpace == 1 && doKaon){
        if(doPion){
          if(m_currentMomOfBin < lowerPionRangeMom){
            m_currentLowFitRange = lowFitRangeForPionsInKaonSpace;
          }else if(m_currentMomOfBin < upperPionRangeMom){
            m_currentLowFitRange = lowerPionFitRangeForKaonLinFunct->Eval(m_currentMomOfBin);
          }else{
            m_currentLowFitRange = -0.1;
          }
          if(m_currentMomOfBin > 1.1){
            m_currentLowFitRange = -0.55;
            if(doProton){
              m_currentHighFitRange = 1.5;
            }
          }
        }else{
          if(m_currentMomOfBin < highMomToIgnorePionsInKaonSpace){
            m_currentLowFitRange  = reducedIsolatedKaon_lowFitRange;
            m_currentHighFitRange = reducedIsolatedKaon_highFitRange;
          }else{
            m_currentLowFitRange = 0.15;
          }
        }
        if(doProton){
          if(m_currentMomOfBin < lowerProtonRangeMom){
            m_currentHighFitRange = highFitRangeForProtonsInKaonSpace;
          }else if(m_currentMomOfBin < upperProtonRangeMom){
            m_currentHighFitRange = upperProtonFitRangeForKaonLinFunct->Eval(m_currentMomOfBin);
          }else{
            m_currentHighFitRange = 1.1;
          }
        }else{
          m_currentHighFitRange = reducedIsolatedKaon_highFitRange;
        }
      }


      //Fix the fit ranges for proton
      if(a_partIndexSpace == 2 && doProton){
        m_currentLowFitRange = 0.5;
        if(doKaon)   m_currentLowFitRange = lowFitRangeForKaonsInProtonSpace;
        if(doPion)   m_currentLowFitRange = lowFitRangeForPionsInProtonSpace;
        m_currentHighFitRange = 1.5;
      }

      #ifdef _ZFITTER_DEBUG_
        cout << " LowFitRange: " << m_currentLowFitRange << "  highFitRange: " << m_currentHighFitRange << endl;
      #endif
      for(int centIndex = 0;centIndex < m_numCentralities; centIndex++){
        int binLow  =  m_currentHistosToFit_ByCent[centIndex]->FindBin(m_currentLowFitRange);
        if(binLow > 1) m_currentLowFitBin[centIndex] = binLow;
        int binHigh =  m_currentHistosToFit_ByCent[centIndex]->FindBin(m_currentHighFitRange);
        if(binHigh  <  m_currentHistosToFit_ByCent[centIndex]->GetNbinsX()) m_currentHighFitBin[centIndex] = binHigh;
        #ifdef _ZFITTER_DEBUG_
          cout << "    Cent: " << centIndex << "  bin " << m_currentLowFitBin[centIndex] << " to " << m_currentHighFitBin[centIndex] << endl;
        #endif
      }

      //// Modify the Histogram Ranges
      //for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      //  truncHisto    = HistogramUtilities::truncateHistogram(m_currentHistosToFit_ByCent[centIndex],lowFitRange,highFitRange);
      //  orignalHisto  = m_currentHistosToFit_ByCent[centIndex];
      //  m_currentHistosToFit_ByCent[centIndex] = truncHisto;
      //  delete orignalHisto;
      //  HistogramUtilities::removeHighCountCorruptData(m_currentHistosToFit_ByCent[centIndex]);
      //  //m_rebinZBTOFHistoStorage[centIndex][m_currentMtM0Bin] = m_currentHistosToFit_ByCent[centIndex];
      //}

      previousFitToggle_Pion   = doPion;
      previousFitToggle_Kaon   = doKaon;
      previousFitToggle_Proton = doProton;



      //Aggregate initial parameters and limits
      #ifdef _ZFITTER_DEBUG_
        cout << "Grabbing parameter limits" << endl;
      #endif
      double initialParams[55];
      double initialParamLimits[55][2];
      for(int iii = 0; iii < 9 + 3*m_numCentralities; iii++){
        ROOT::Fit::ParameterSettings* setting = new ROOT::Fit::ParameterSettings();
        minimizer->GetVariableSettings(iii, *setting);
        initialParamLimits[iii][0] = setting->LowerLimit();
        initialParamLimits[iii][1] = setting->UpperLimit();
        initialParams[iii]         = setting->Value();
        #ifdef _ZFITTER_DEBUG_
          cout << " Param " << iii << "  " << initialParams[iii] << "   [ " << initialParamLimits[iii][0] << ", " << initialParamLimits[iii][1] << " ] " << endl;
        #endif
        delete setting;
      }

      #ifdef _ZFITTER_DEBUG_
        cout << "Minimizing" << endl;
      #endif
      minimizer->Minimize();
      int status_Minimizer = minimizer->Status();
      const double* finalParams = minimizer->X();
      if(doHesse && (status_Minimizer == 0 || status_Minimizer == 4)) minimizer->Hesse();
      #ifdef _ZFITTER_DEBUG_
        cout << "simultaneous centrality minimizer status: " << status_Minimizer << endl;
      #endif

      // 2026-07: chi^2/ndf fit-quality diagnostic -- BTOF counterpart of the same fix
      // in ZFitter_PionKaonTPC_SimulCent_SingleLoop.cxx (see that file's Minimize()
      // block for the full explanation). Andrew asked for this to be saved per fit so
      // PresentZFitterSpectra.C's a_maxChiSqrNdf cut has real data to read -- the
      // m_ChiSqr_Cent_ZbTOF histograms were booked (ZFitter.cxx) but never filled in
      // this live joint/simultaneous-centrality fit path, confirmed via
      // macros/CheckChiSqrValues.C on a real ZFitter run (chi^2/ndf read back as
      // exactly 0 for every BTOF bin).
      //
      // simultaneous_centrality_ZTOF_chisqr (the minimizer's objective function,
      // defined earlier in this file) is the sum-of-squared-residuals chi^2 passed to
      // ROOT::Math::Minimizer via the Functor, so minimizer->MinValue() after
      // Minimize() is exactly the chi^2 at the best-fit point -- same reasoning as the
      // TPC fix, no separate recomputation needed. ndf = (data points with content>0
      // within [m_currentLowFitBin,m_currentHighFitBin] for
      // m_currentHistosToFit_ByCent_InvBeta[centIndex], summed over ALL linked
      // centralities -- mirroring exactly what the chisqr functor sums over) minus
      // (minimizer->NFree()).
      //
      // CAVEAT: same as the TPC fix -- centralities are fit JOINTLY here (pion+kaon+
      // proton templates x all centralities in one simultaneous minimization), so this
      // is ONE chi^2/ndf value per (rapidity bin, mT-m0 bin), not an independent
      // per-centrality figure. It's written into every centrality's m_ChiSqr_Cent_ZbTOF
      // slot below (in the centIndex loop further down) so PresentZFitterSpectra.C's
      // existing per-centrality-indexed read/cut works unchanged, but a bad fit
      // anywhere shows up identically in every centrality's reported value.
      double jointFitChiSqr = minimizer->MinValue();
      int jointFitNDataPoints = 0;
      for(int centIndexForNdf = 0; centIndexForNdf < m_numCentralities; centIndexForNdf++){
        for(int binX = m_currentLowFitBin[centIndexForNdf]; binX <= m_currentHighFitBin[centIndexForNdf]; binX++){
          if(m_currentHistosToFit_ByCent_InvBeta[centIndexForNdf]->GetBinContent(binX) > 0.0) jointFitNDataPoints++;
        }
      }
      int jointFitNDF = jointFitNDataPoints - (int) minimizer->NFree();
      // ndf<=0 (underconstrained fit) -> write -1 rather than divide by a non-positive
      // number; suppressBadFitQualityBins() in PresentZFitterSpectra.C already treats
      // chiSqrNdf<=0 as a failed fit.
      double jointFitChiSqrNDF = (jointFitNDF > 0) ? (jointFitChiSqr / jointFitNDF) : -1.0;
      #ifdef _ZFITTER_DEBUG_
        cout << "simultaneous centrality joint chi^2/ndf: " << jointFitChiSqr << " / " << jointFitNDF
             << " = " << jointFitChiSqrNDF << "  (nDataPoints=" << jointFitNDataPoints
             << ", nFree=" << minimizer->NFree() << ")" << endl;
      #endif

      bool okayFitParams = true;
      bool okayPionNu   = true;
      bool okayKaonNu   = true;
      bool okayProtonNu = true;
      bool peaksInOrder = true;
      for(int iii = 0; iii < 12; iii++){
        parameters_simul_cent[iii] = finalParams[iii];
        parameters_simul_cent_errors[iii] = sqrt(minimizer->CovMatrix(iii,iii));
        if(!std::isfinite(parameters_simul_cent[iii]))        okayFitParams = false;
        if(!std::isfinite(parameters_simul_cent_errors[iii])) okayFitParams = false;
        if(iii%4 == 3 && parameters_simul_cent[iii] > 200) okayFitParams = false;
      }

      if((doPion && doKaon && parameters_simul_cent[0] > parameters_simul_cent[4]) ||
         (doKaon && doProton && parameters_simul_cent[4] > parameters_simul_cent[8]) ||
         (doPion && doProton && parameters_simul_cent[0] > parameters_simul_cent[8])
        ){
        peaksInOrder = false;
        okayFitParams = false;
      }

      bool pionKaonSwap = false;
      if(doPion && doKaon && parameters_simul_cent[0] > parameters_simul_cent[4]){
        pionKaonSwap = true;
        parameters_simul_cent[0] = finalParams[4];
        parameters_simul_cent[1] = finalParams[5];
        parameters_simul_cent[2] = finalParams[6];
        parameters_simul_cent[3] = finalParams[7];
        parameters_simul_cent_errors[0] = sqrt(minimizer->CovMatrix(4,4));
        parameters_simul_cent_errors[1] = sqrt(minimizer->CovMatrix(5,5));
        parameters_simul_cent_errors[2] = sqrt(minimizer->CovMatrix(6,6));
        parameters_simul_cent_errors[3] = sqrt(minimizer->CovMatrix(7,7));
        parameters_simul_cent[4] = finalParams[0];
        parameters_simul_cent[5] = finalParams[1];
        parameters_simul_cent[6] = finalParams[2];
        parameters_simul_cent[7] = finalParams[3];
        parameters_simul_cent_errors[4] = sqrt(minimizer->CovMatrix(0,0));
        parameters_simul_cent_errors[5] = sqrt(minimizer->CovMatrix(1,1));
        parameters_simul_cent_errors[6] = sqrt(minimizer->CovMatrix(2,2));
        parameters_simul_cent_errors[7] = sqrt(minimizer->CovMatrix(3,3));
      }

      double yieldsCentralityIntegrated[3] = {0,0,0};
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        if(pionKaonSwap){
          parameters_simul_cent[12+3*centIndex]       = finalParams[13+3*centIndex];
          parameters_simul_cent_errors[12+3*centIndex] = sqrt(minimizer->CovMatrix(13+3*centIndex, 13+3*centIndex));
          parameters_simul_cent[13+3*centIndex]       = finalParams[12+3*centIndex];
          parameters_simul_cent_errors[13+3*centIndex] = sqrt(minimizer->CovMatrix(12+3*centIndex, 12+3*centIndex));
        }else{
          parameters_simul_cent[12+3*centIndex]       = finalParams[12+3*centIndex];
          parameters_simul_cent_errors[12+3*centIndex] = sqrt(minimizer->CovMatrix(12+3*centIndex, 12+3*centIndex));
          parameters_simul_cent[13+3*centIndex]       = finalParams[13+3*centIndex];
          parameters_simul_cent_errors[13+3*centIndex] = sqrt(minimizer->CovMatrix(13+3*centIndex, 13+3*centIndex));
        }

        parameters_simul_cent[14+3*centIndex]       = finalParams[14+3*centIndex];
        parameters_simul_cent_errors[14+3*centIndex] = sqrt(minimizer->CovMatrix(14+3*centIndex, 14+3*centIndex));
        if(!std::isfinite(parameters_simul_cent[12+3*centIndex]))        okayFitParams = false;
        if(!std::isfinite(parameters_simul_cent_errors[12+3*centIndex])) okayFitParams = false;
        if(!std::isfinite(parameters_simul_cent[13+3*centIndex]))        okayFitParams = false;
        if(!std::isfinite(parameters_simul_cent_errors[13+3*centIndex])) okayFitParams = false;
        if(!std::isfinite(parameters_simul_cent[14+3*centIndex]))        okayFitParams = false;
        if(!std::isfinite(parameters_simul_cent_errors[14+3*centIndex])) okayFitParams = false;

        previousCentralityIntegratedYields_Pion   += parameters_simul_cent[12+3*centIndex]*m_centEvents[m_currentPartIndex]->GetBinContent(1 + centIndex);
        previousCentralityIntegratedYields_Kaon   += parameters_simul_cent[13+3*centIndex]*m_centEvents[m_currentPartIndex]->GetBinContent(1 + centIndex);
        previousCentralityIntegratedYields_Proton += parameters_simul_cent[14+3*centIndex]*m_centEvents[m_currentPartIndex]->GetBinContent(1 + centIndex);
      }

      #ifdef _ZFITTER_DEBUG_
        cout << "Centrality Integrated Yields: " << endl;
        cout << "   PION:   " << previousCentralityIntegratedYields_Pion << endl;
        cout << "   KAON:   " << previousCentralityIntegratedYields_Kaon << endl;
        cout << "   PROTON: " << previousCentralityIntegratedYields_Proton << endl;
      #endif

      if(previousCentralityIntegratedYields_Pion   > centIntThreshold) centralityIntegratedYieldsWereAboveThreshold_Pion = true;
      if(previousCentralityIntegratedYields_Kaon   > centIntThreshold) centralityIntegratedYieldsWereAboveThreshold_Kaon = true;
      if(previousCentralityIntegratedYields_Proton > centIntThreshold) centralityIntegratedYieldsWereAboveThreshold_Proton = true;



      /*
      if(!okayFitParams){    // second round of fitting if first failed
        #ifdef _ZFITTER_DEBUG_
          cout << "WARNING: Initial Fit has nan/inf results or nu is too high!   Trying again! " << endl;
        #endif
        for(int iii = 0; iii < 12; iii++){
          if(iii%4 == 3 && parameters_simul_cent[iii] > 200) minimizer->SetFixedVariable(iii,paramNames[iii],200);
        }

        if(!peaksInOrder){ // force refit if they got swapped!
          minimizer->SetVariable(0,paramNames[0],fixed_pion_mass_sqr,0.01);
          minimizer->SetVariable(4,paramNames[4],fixed_kaon_mass_sqr,0.01);
          minimizer->SetVariable(8,paramNames[8],fixed_proton_mass_sqr,0.01);
        }


        minimizer->Minimize();
        status_Minimizer = minimizer->Status();
        const double* finalParamsRoundTwo = minimizer->X();
        //if(status_Minimizer == 0)
        minimizer->Hesse();
        cout << "simultaneous centrality minimizer status second round: " << status_Minimizer << endl;
        okayFitParams = true;
        for(int iii = 0; iii < 12; iii++){
          parameters_simul_cent[iii] = finalParamsRoundTwo[iii];
          parameters_simul_cent_errors[iii] = sqrt(minimizer->CovMatrix(iii,iii));
          if(!std::isfinite(parameters_simul_cent[iii]))        okayFitParams = false;
          if(!std::isfinite(parameters_simul_cent_errors[iii])) okayFitParams = false;
        }
        for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
          parameters_simul_cent[12+3*centIndex]       = finalParamsRoundTwo[12+3*centIndex];
          parameters_simul_cent_errors[12+3*centIndex] = sqrt(minimizer->CovMatrix(12+3*centIndex, 12+3*centIndex));
          parameters_simul_cent[13+3*centIndex]       = finalParamsRoundTwo[13+3*centIndex];
          parameters_simul_cent_errors[13+3*centIndex] = sqrt(minimizer->CovMatrix(13+3*centIndex, 13+3*centIndex));
          parameters_simul_cent[14+3*centIndex]       = finalParamsRoundTwo[14+3*centIndex];
          parameters_simul_cent_errors[14+3*centIndex] = sqrt(minimizer->CovMatrix(14+3*centIndex, 14+3*centIndex));
          if(!std::isfinite(parameters_simul_cent[12+3*centIndex]))        okayFitParams = false;
          if(!std::isfinite(parameters_simul_cent_errors[12+3*centIndex])) okayFitParams = false;
          if(!std::isfinite(parameters_simul_cent[13+3*centIndex]))        okayFitParams = false;
          if(!std::isfinite(parameters_simul_cent_errors[13+3*centIndex])) okayFitParams = false;
          if(!std::isfinite(parameters_simul_cent[14+3*centIndex]))        okayFitParams = false;
          if(!std::isfinite(parameters_simul_cent_errors[14+3*centIndex])) okayFitParams = false;
        }
      }*/






      if(firstMtM0BinFit == -1) firstMtM0BinFit = m_currentMtM0Bin;
      lastMtM0BinFit = m_currentMtM0Bin;

      #ifdef _ZFITTER_DEBUG_
        cout << "Loading the previous good parameters for next fit: " << endl;
      #endif
      for(int paramIndex = 0; paramIndex < prevParamSize; paramIndex++){
        if(paramIndex%4 == 0) previousGoodParameters[paramIndex] = PhysMath::massSquared_signed_fromInvBeta(m_currentMomOfBin, parameters_simul_cent[paramIndex]);
        else previousGoodParameters[paramIndex] = parameters_simul_cent[paramIndex];
        #ifdef _ZFITTER_DEBUG_
          if(paramIndex % 4 == 0) cout << endl;
          cout << "  " << previousGoodParameters[paramIndex];
        #endif
      }

      TF1* combinedFunct_CentStorage[16];
      for(int iii = 0; iii < 16; iii++) combinedFunct_CentStorage[iii] = NULL;

      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        double numberOfCentEvents = m_centEvents[m_currentPartIndex]->GetBinContent(centIndex+1);

        m_currentHistosToFit_ByCent[centIndex]->SetMarkerStyle(20);
        m_currentHistosToFit_ByCent[centIndex]->SetMarkerColor(kBlack);
        m_currentHistosToFit_ByCent[centIndex]->SetMarkerSize(0.5);
        m_currentHistosToFit_ByCent[centIndex]->SetLineWidth(2);
        double highestValue = HistogramUtilities::getMaxInRange(m_currentHistosToFit_ByCent[centIndex],m_currentLowFitRange,m_currentHighFitRange);
        double lowestValue  = HistogramUtilities::getMinInRange(m_currentHistosToFit_ByCent[centIndex],m_currentLowFitRange,m_currentHighFitRange);
        TF1* combinedFunct = NULL;
        TF1* pionFunct     = NULL;
        TF1* kaonFunct     = NULL;
        TF1* protonFunct   = NULL;
        TF1* combinedFunct_invBeta = NULL;
        TF1* pionFunct_invBeta     = NULL;
        TF1* kaonFunct_invBeta     = NULL;
        TF1* protonFunct_invBeta   = NULL;
        TLine* lowRangeLine  = new TLine(m_currentLowFitRange, lowestValue,m_currentLowFitRange, highestValue);
        TLine* highRangeLine = new TLine(m_currentHighFitRange,lowestValue,m_currentHighFitRange,highestValue);
        lowRangeLine->SetLineColor(kViolet);
        highRangeLine->SetLineColor(kViolet);

        #ifdef _GAUSSIAN_BTOF_SIMUL_CENT_
          combinedFunct = PhysMath::getThreeGausFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_combined",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent,
            parameters_simul_cent[12 + centIndex*3],parameters_simul_cent[13 + centIndex*3],parameters_simul_cent[14 + centIndex*3]);
          pionFunct  = PhysMath::getGaussianFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Pion",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[12 + centIndex*3], parameters_simul_cent[0], parameters_simul_cent[1]);
          kaonFunct  = PhysMath::getGaussianFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Kaon",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[13 + centIndex*3], parameters_simul_cent[4], parameters_simul_cent[5]);
          protonFunct = PhysMath::getGaussianFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Proton",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[14 + centIndex*3], parameters_simul_cent[8], parameters_simul_cent[9]);
        #endif
        #ifdef _SKEW_NORMAL_BTOF_SIMUL_CENT_
          combinedFunct = PhysMath::getThreeSkewFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_combined",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent,
            parameters_simul_cent[12 + centIndex*3],parameters_simul_cent[13 + centIndex*3],parameters_simul_cent[14 + centIndex*3]);
            pionFunct    = PhysMath::skewNormalFunct( parameters_simul_cent[12 + centIndex*3], parameters_simul_cent[0], parameters_simul_cent[1], parameters_simul_cent[2]);
            kaonFunct    = PhysMath::skewNormalFunct( parameters_simul_cent[13 + centIndex*3], parameters_simul_cent[4], parameters_simul_cent[5], parameters_simul_cent[6]);
            protonFunct  = PhysMath::skewNormalFunct( parameters_simul_cent[14 + centIndex*3], parameters_simul_cent[8], parameters_simul_cent[9], parameters_simul_cent[10]);
        #endif
        #ifdef _STUDENT_T_BTOF_SIMUL_CENT_
          combinedFunct = PhysMath::getThreeM2StudentTFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_combined",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent,
            parameters_simul_cent[12 + centIndex*3],parameters_simul_cent[13 + centIndex*3],parameters_simul_cent[14 + centIndex*3],m_currentMomOfBin);
          pionFunct  = PhysMath::getM2StudentTFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Pion",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[12 + centIndex*3], parameters_simul_cent[0], parameters_simul_cent[1], parameters_simul_cent[2],m_currentMomOfBin);
          kaonFunct  = PhysMath::getM2StudentTFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Kaon",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[13 + centIndex*3], parameters_simul_cent[4], parameters_simul_cent[5], parameters_simul_cent[6],m_currentMomOfBin);
          protonFunct = PhysMath::getM2StudentTFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Proton",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[14 + centIndex*3], parameters_simul_cent[8], parameters_simul_cent[9], parameters_simul_cent[10],m_currentMomOfBin);
          combinedFunct_invBeta = PhysMath::getThreeStudentTFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_combined_invBeta",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent,
            parameters_simul_cent[12 + centIndex*3],parameters_simul_cent[13 + centIndex*3],parameters_simul_cent[14 + centIndex*3]);
          pionFunct_invBeta  = PhysMath::getStudentTFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Pion_invBeta",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[12 + centIndex*3], parameters_simul_cent[0], parameters_simul_cent[1], parameters_simul_cent[2]);
          kaonFunct_invBeta  = PhysMath::getStudentTFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Kaon_invBeta",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[13 + centIndex*3], parameters_simul_cent[4], parameters_simul_cent[5], parameters_simul_cent[6]);
          protonFunct_invBeta = PhysMath::getStudentTFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Proton_invBeta",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[14 + centIndex*3], parameters_simul_cent[8], parameters_simul_cent[9], parameters_simul_cent[10]);
        #endif
        #ifdef _PEARSON_IV_BTOF_SIMUL_CENT_
          combinedFunct = PhysMath::getThreePearsonIVFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_combined",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent,
            parameters_simul_cent[12 + centIndex*3],parameters_simul_cent[13 + centIndex*3],parameters_simul_cent[14 + centIndex*3]);
            pionFunct    = PhysMath::pearsonIVFunct( parameters_simul_cent[12 + centIndex*3], parameters_simul_cent[0], parameters_simul_cent[1], parameters_simul_cent[2], parameters_simul_cent[3]);
            kaonFunct    = PhysMath::pearsonIVFunct( parameters_simul_cent[13 + centIndex*3], parameters_simul_cent[4], parameters_simul_cent[5], parameters_simul_cent[6], parameters_simul_cent[7]);
            protonFunct  = PhysMath::pearsonIVFunct( parameters_simul_cent[14 + centIndex*3], parameters_simul_cent[8], parameters_simul_cent[9], parameters_simul_cent[10], parameters_simul_cent[11]);
        #endif
        #ifdef _DOUBLE_GAUS_BTOF_SIMUL_CENT_
          combinedFunct = PhysMath::getThreeDoubleGaussFunct(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_combined",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent,
            parameters_simul_cent[12 + centIndex*3],parameters_simul_cent[13 + centIndex*3],parameters_simul_cent[14 + centIndex*3]);
          pionFunct  = PhysMath::getDoubleGaussianFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Pion",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[12 + centIndex*3], parameters_simul_cent[0], parameters_simul_cent[1], parameters_simul_cent[2], parameters_simul_cent[3]);
          kaonFunct  = PhysMath::getDoubleGaussianFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Kaon",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[13 + centIndex*3], parameters_simul_cent[4], parameters_simul_cent[5], parameters_simul_cent[6], parameters_simul_cent[7]);
          protonFunct = PhysMath::getDoubleGaussianFunction(Form("zbtof_step3_fit_%d_%d_%02d_%03d_%04d_Proton",
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
            parameters_simul_cent[14 + centIndex*3], parameters_simul_cent[8], parameters_simul_cent[9], parameters_simul_cent[10], parameters_simul_cent[11]);
        #endif

        combinedFunct_CentStorage[centIndex] = combinedFunct;
        if(pionFunct){
          pionFunct->SetLineColor(m_partInfo->GetParticleColor(0));
          pionFunct->SetRange(lowRange,highRange);
          pionFunct_invBeta->SetLineColor(m_partInfo->GetParticleColor(0));
          pionFunct_invBeta->SetRange(lowRange,highRange);
        }
        if(kaonFunct){
          kaonFunct->SetLineColor(m_partInfo->GetParticleColor(1));
          kaonFunct->SetRange(lowRange,highRange);
          kaonFunct_invBeta->SetLineColor(m_partInfo->GetParticleColor(1));
          kaonFunct_invBeta->SetRange(lowRange,highRange);
        }
        if(protonFunct){
          protonFunct->SetLineColor(m_partInfo->GetParticleColor(2));
          protonFunct->SetRange(lowRange,highRange);
          protonFunct_invBeta->SetLineColor(m_partInfo->GetParticleColor(2));
          protonFunct_invBeta->SetRange(lowRange,highRange);
        }
        double rangeExtension = 0.01;
        double lowHorDrawRange  = lowRange - rangeExtension;
        double highHorDrawRange =  highRange + rangeExtension;
        HistogramUtilities::setAxisRanges(m_currentHistosToFit_ByCent[centIndex], lowHorDrawRange, highHorDrawRange,true);
        m_currentHistosToFit_ByCent[centIndex]->SetStats(kFALSE);


        TCanvas* fittingCanvas = new TCanvas(Form("fittingCanvas_SimulBTOF_Cent%02d",centIndex),"Fitting Canvas");
        fittingCanvas->SetWindowSize(1300,700);
        fittingCanvas->SetCanvasSize(1200,600);
        fittingCanvas->Draw();
        TPad* topPad    = new TPad(Form("topPad_SimulBTOF_Cent%02d",centIndex),   "topPad_SimulBTOF",   0.0,0.3,1.0,1.0);
        TPad* bottomPad = new TPad(Form("bottomPad_SimulBTOF_Cent%02d",centIndex),"bottomPad_SimulBTOF",0.0,0.0,1.0,0.3);
        topPad->Draw();
        bottomPad->Draw();
        fittingCanvas->cd();
        topPad->cd();
        topPad->SetLogy(true);
        topPad->SetRightMargin(0.4);
        bottomPad->SetLogy(false);
        bottomPad->SetRightMargin(0.4);
        gStyle->SetOptFit(0111);

        m_currentHistosToFit_ByCent[centIndex]->Draw("E");
        lowRangeLine->Draw("same");
        highRangeLine->Draw("same");
        if(initialFunct[centIndex]) initialFunct[centIndex]->Draw("same");
        if(combinedFunct) combinedFunct->Draw("same");
        if(pionFunct)     pionFunct->Draw("same");
        if(kaonFunct)     kaonFunct->Draw("same");
        if(protonFunct)   protonFunct->Draw("same");

        if(initialPoints[centIndex][0]) initialPoints[centIndex][0]->Draw("P");
        if(initialPoints[centIndex][1]) initialPoints[centIndex][1]->Draw("P");
        if(initialPoints[centIndex][2]) initialPoints[centIndex][2]->Draw("P");


        TPaveText* fitParsTxt = new TPaveText(0.65,0.05,0.95,0.90,"NDC");
        if(parameters_simul_cent[12 + centIndex*3] != 0.0) fitParsTxt->AddText(Form("N_{#pi} \t %.2e #pm %.2e",parameters_simul_cent[12 + centIndex*3],parameters_simul_cent_errors[12 + centIndex*3]));
        if(parameters_simul_cent[13 + centIndex*3] != 0.0) fitParsTxt->AddText(Form("N_{K} \t %.2e #pm %.2e",parameters_simul_cent[13 + centIndex*3],parameters_simul_cent_errors[13 + centIndex*3]));
        if(parameters_simul_cent[14 + centIndex*3] != 0.0) fitParsTxt->AddText(Form("N_{p} \t %.2e #pm %.2e",parameters_simul_cent[14 + centIndex*3],parameters_simul_cent_errors[14 + centIndex*3]));
        for(int iii = 0; iii < 12; iii++){
          int index = iii/4;
          if( ((index == 0 && doPion) || (index == 1 && doKaon) || (index == 2 && doProton)) && paramNames[iii].find("Ignore") == std::string::npos){
            fitParsTxt->AddText(Form("%s \t %.2e #pm %.2e %s",paramNames[iii].c_str(),parameters_simul_cent[iii],parameters_simul_cent_errors[iii],(minimizer->IsFixedVariable(iii) ? "(fixed)" : Form(" #pm %.2e", parameters_simul_cent_errors[iii]))));
          }
        }
        fitParsTxt->SetFillStyle(0);
        fitParsTxt->SetBorderSize(5);
        fitParsTxt->Draw("same");

        bottomPad->cd();
        TH1D* pullHisto = HistogramUtilities::makeFitPullHistogram(m_currentHistosToFit_ByCent[centIndex], combinedFunct, true);
        pullHisto->SetMarkerStyle(20);
        pullHisto->SetMarkerSize(0.5);
        pullHisto->SetMarkerColor(kBlack);
        pullHisto->SetLineColor(kBlack);
        pullHisto->SetStats(false);

        TPaveText* fitTxt = new TPaveText(0.65,0.05,0.95,0.95,"NDC");
        if(initialParams[12 + centIndex*3] != 0.0) fitTxt->AddText(Form("N_{#pi} \t %.2e  [ %.2e, %.2e ]",initialParams[12 + centIndex*3],initialParamLimits[12 + centIndex*3][0],initialParamLimits[12 + centIndex*3][1]));
        if(initialParams[13 + centIndex*3] != 0.0) fitTxt->AddText(Form("N_{K} \t %.2e  [ %.2e, %.2e ]",  initialParams[13 + centIndex*3],initialParamLimits[13 + centIndex*3][0],initialParamLimits[13 + centIndex*3][1]));
        if(initialParams[14 + centIndex*3] != 0.0) fitTxt->AddText(Form("N_{p} \t %.2e  [ %.2e, %.2e ]",  initialParams[14 + centIndex*3],initialParamLimits[14 + centIndex*3][0],initialParamLimits[14 + centIndex*3][1]));
        for(int iii = 0; iii < 12; iii++){
          if(!minimizer->IsFixedVariable(iii) && paramNames[iii].find("Ignore") == std::string::npos){
            double percent = 100.0*(parameters_simul_cent[iii] - initialParamLimits[iii][0])/(initialParamLimits[iii][1] - initialParamLimits[iii][0]);
            bool atLimit = (percent > 99.9 ||  percent < 0.1);
            string limitStr = "";
            if(atLimit) limitStr = "AT LIMIT";
            fitTxt->AddText(Form("%s \t %f  [ %f, %f ]  %1.2f %% %s",paramNames[iii].c_str(),initialParams[iii],initialParamLimits[iii][0],initialParamLimits[iii][1], percent,limitStr.c_str() ));
          }
        }
        fitTxt->SetFillStyle(0);
        fitTxt->SetBorderSize(5);


        //TH1F* frame_pull = bottomPad->DrawFrame(m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin(lowHorDrawRange)),-15.0,
        //       m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( highHorDrawRange)),15.0);
        pullHisto->SetTitle(Form("; BTOF m^{2}; Standardized Pull"));
        pullHisto->Draw("E");
        pullHisto->GetXaxis()->SetRangeUser(lowHorDrawRange,highHorDrawRange);
        pullHisto->GetYaxis()->SetRangeUser(-15,15);
        TLine* oneLine = new TLine(lowHorDrawRange,0.0,highHorDrawRange,0.0);
        oneLine->SetLineStyle(7);
        oneLine->SetLineColor(kBlack);
        oneLine->Draw("same");
        fitTxt->Draw("same");
        gPad->Update();
        gSystem->ProcessEvents();
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
             m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex, m_currentRapBin, m_currentMtM0Bin));


        //gPad->SetLogy(false);
        //gPad->Update();
        if(m_saveNoLogImages){
          topPad->SetLogy(false);
          topPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
               m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex, m_currentRapBin, m_currentMtM0Bin));
        }









        if(drawInvBetaFits){
          double lowHorDrawRange_InvBeta  = PhysMath::invBeta_from_MassSquaredSigned(lowHorDrawRange,m_currentMomOfBin);
          double highHorDrawRange_InvBeta = PhysMath::invBeta_from_MassSquaredSigned(highHorDrawRange,m_currentMomOfBin);

          //HistogramUtilities::setAxisRanges(m_currentHistosToFit_ByCent_InvBeta[centIndex], lowHorDrawRange_InvBeta, highHorDrawRange_InvBeta,true);
          m_currentHistosToFit_ByCent_InvBeta[centIndex]->SetStats(kFALSE);
          TCanvas* fittingCanvas_InvBeta = new TCanvas(Form("fittingCanvas_SimulBTOF_Cent%02d_InvBeta",centIndex),"Fitting Canvas");
          fittingCanvas_InvBeta->SetWindowSize(1300,700);
          fittingCanvas_InvBeta->SetCanvasSize(1200,600);
          fittingCanvas_InvBeta->Draw();
          TPad* topPad_InvBeta    = new TPad(Form("topPad_SimulBTOF_Cent%02d_InvBeta",centIndex),   "topPad_SimulBTOF",   0.0,0.3,1.0,1.0);
          TPad* bottomPad_InvBeta = new TPad(Form("bottomPad_SimulBTOF_Cent%02d_InvBeta",centIndex),"bottomPad_SimulBTOF",0.0,0.0,1.0,0.3);
          topPad_InvBeta->Draw();
          bottomPad_InvBeta->Draw();
          fittingCanvas_InvBeta->cd();
          topPad_InvBeta->cd();
          topPad_InvBeta->SetLogy(true);
          topPad_InvBeta->SetRightMargin(0.4);
          bottomPad_InvBeta->SetLogy(false);
          bottomPad_InvBeta->SetRightMargin(0.4);
          gStyle->SetOptFit(0111);

          m_currentHistosToFit_ByCent_InvBeta[centIndex]->Draw("E");
          if(initialFunct_invBeta[centIndex]) initialFunct_invBeta[centIndex]->Draw("same");
          if(combinedFunct_invBeta) combinedFunct_invBeta->Draw("same");
          if(pionFunct_invBeta)     pionFunct_invBeta->Draw("same");
          if(kaonFunct_invBeta)     kaonFunct_invBeta->Draw("same");
          if(protonFunct_invBeta)   protonFunct_invBeta->Draw("same");
          fitParsTxt->Draw("same");

          bottomPad_InvBeta->cd();
          TH1D* pullHisto_InvBeta = HistogramUtilities::makeFitPullHistogram(m_currentHistosToFit_ByCent_InvBeta[centIndex], combinedFunct_invBeta, true);
          pullHisto_InvBeta->SetMarkerStyle(20);
          pullHisto_InvBeta->SetMarkerSize(0.5);
          pullHisto_InvBeta->SetMarkerColor(kBlack);
          pullHisto_InvBeta->SetLineColor(kBlack);
          pullHisto_InvBeta->SetStats(false);

          //TH1F* frame_pull = bottomPad->DrawFrame(m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin(lowHorDrawRange)),-15.0,
          //       m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( highHorDrawRange)),15.0);
          pullHisto_InvBeta->SetTitle(Form("; BTOF 1/#beta; Standardized Pull"));
          pullHisto_InvBeta->Draw("E");
          //pullHisto_InvBeta->GetXaxis()->SetRangeUser(lowHorDrawRange_InvBeta,highHorDrawRange_InvBeta);
          pullHisto_InvBeta->GetYaxis()->SetRangeUser(-15,15);
          TLine* oneLine_InvBeta = new TLine(lowHorDrawRange_InvBeta,0.0,highHorDrawRange_InvBeta,0.0);
          oneLine_InvBeta->SetLineStyle(7);
          oneLine_InvBeta->SetLineColor(kBlack);
          oneLine_InvBeta->Draw("same");
          fitTxt->Draw("same");
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas_InvBeta->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_InvBeta.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
               m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex, m_currentRapBin, m_currentMtM0Bin));
          if(m_saveNoLogImages){
            topPad_InvBeta->SetLogy(false);
            topPad_InvBeta->Update();
            gSystem->ProcessEvents();
            fittingCanvas_InvBeta->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_InvBeta.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
                 m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex, m_currentRapBin, m_currentMtM0Bin));
          }

          if(pionFunct_invBeta)    delete pionFunct_invBeta;
          if(kaonFunct_invBeta)    delete kaonFunct_invBeta;
          if(protonFunct_invBeta)  delete protonFunct_invBeta;
          if(pullHisto_InvBeta)    delete pullHisto_InvBeta;
          if(oneLine_InvBeta)      delete oneLine_InvBeta;
          if(topPad_InvBeta)       delete topPad_InvBeta;
          if(bottomPad_InvBeta)    delete bottomPad_InvBeta;
          if(fittingCanvas_InvBeta)   delete fittingCanvas_InvBeta;
        }

        //---------------------------------  push the fit parameters to the gobal storage histograms  ---------------------------------------
        if(centIndex == 0){
          if(doPion){
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][0][1][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[0],parameters_simul_cent_errors[0]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][0][2][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[1],parameters_simul_cent_errors[1]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][0][3][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[2],parameters_simul_cent_errors[2]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][0][4][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[3],parameters_simul_cent_errors[3]);
          }
          if(doKaon){
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][1][1][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[4],parameters_simul_cent_errors[4]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][1][2][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[5],parameters_simul_cent_errors[5]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][1][3][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[6],parameters_simul_cent_errors[6]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][1][4][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[7],parameters_simul_cent_errors[7]);
          }
          if(doProton){
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][2][1][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[8],parameters_simul_cent_errors[8]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][2][2][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[9],parameters_simul_cent_errors[9]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][2][3][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[10],parameters_simul_cent_errors[10]);
            HistogramUtilities::setBinValues(m_Fit_Data_ZbTOF[m_currentPartIndex][2][4][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[11],parameters_simul_cent_errors[11]);
          }
        }

        if(doPion){
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[12+3*centIndex] * numberOfCentEvents,parameters_simul_cent_errors[12+3*centIndex]*numberOfCentEvents);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][0][centIndex][1][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[0],parameters_simul_cent_errors[0]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][0][centIndex][2][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[1],parameters_simul_cent_errors[1]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][0][centIndex][3][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[2],parameters_simul_cent_errors[2]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][0][centIndex][4][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[3],parameters_simul_cent_errors[3]);
        }
        if(doKaon){
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[13+3*centIndex] * numberOfCentEvents,parameters_simul_cent_errors[13+3*centIndex]*numberOfCentEvents);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][1][centIndex][1][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[4],parameters_simul_cent_errors[4]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][1][centIndex][2][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[5],parameters_simul_cent_errors[5]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][1][centIndex][3][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[6],parameters_simul_cent_errors[6]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][1][centIndex][4][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[7],parameters_simul_cent_errors[7]);
        }
        if(doProton){
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][2][centIndex][0][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[14+3*centIndex] * numberOfCentEvents,parameters_simul_cent_errors[14+3*centIndex]*numberOfCentEvents);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][2][centIndex][1][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[8],parameters_simul_cent_errors[8]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][2][centIndex][2][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[9],parameters_simul_cent_errors[9]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][2][centIndex][3][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[10],parameters_simul_cent_errors[10]);
          HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][2][centIndex][4][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[11],parameters_simul_cent_errors[11]);
        }

        // 2026-07: shared joint-fit chi^2/ndf (computed once above, right after
        // Minimize()) written into every centrality's slot -- see the CAVEAT comment
        // there. This is what CheckChiSqrValues.C / PresentZFitterSpectra.C's
        // a_maxChiSqrNdf cut reads back.
        m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, jointFitChiSqrNDF);

        #ifdef _CREATE_TEXT_FIT_FILE_
          //if(!outputTextFile.is_open()) cout << "file is no longer open..." << endl;
          ofstream outputTextFile_local(outputTextFileName.c_str(),ios::app);
          outputTextFile_local << m_currentPartIndex << " " << m_currentPlusMinusIndex << " " << 1 << " " << centIndex << " " << m_currentRapBin << " " << m_currentMtM0Bin << " ";
          double spectraValue_forText = 0;
          double spectraValueErr_forText = 0;
        #endif
        if(m_currentPartIndex == 0){
          HistogramUtilities::setBinValues(m_Spectra_Cent_ZbTOF[m_currentPartIndex][centIndex][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[12+3*centIndex] * numberOfCentEvents,parameters_simul_cent_errors[12+3*centIndex]*numberOfCentEvents);
          #ifdef _CREATE_TEXT_FIT_FILE_
            spectraValue_forText = parameters_simul_cent[12+3*centIndex] * numberOfCentEvents;
            spectraValueErr_forText = parameters_simul_cent_errors[12+3*centIndex] * numberOfCentEvents;
          #endif
        }
        if(m_currentPartIndex == 1){
          HistogramUtilities::setBinValues(m_Spectra_Cent_ZbTOF[m_currentPartIndex][centIndex][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[13+3*centIndex] * numberOfCentEvents,parameters_simul_cent_errors[13+3*centIndex]*numberOfCentEvents);
          #ifdef _CREATE_TEXT_FIT_FILE_
            spectraValue_forText = parameters_simul_cent[13+3*centIndex] * numberOfCentEvents;
            spectraValueErr_forText = parameters_simul_cent_errors[13+3*centIndex] * numberOfCentEvents;
          #endif
        }
        if(m_currentPartIndex == 2){
          HistogramUtilities::setBinValues(m_Spectra_Cent_ZbTOF[m_currentPartIndex][centIndex][m_currentPlusMinusIndex],m_currentRapBin,m_currentMtM0Bin, parameters_simul_cent[14+3*centIndex] * numberOfCentEvents,parameters_simul_cent_errors[14+3*centIndex]*numberOfCentEvents);
          #ifdef _CREATE_TEXT_FIT_FILE_
            spectraValue_forText = parameters_simul_cent[14+3*centIndex] * numberOfCentEvents;
            spectraValueErr_forText = parameters_simul_cent_errors[14+3*centIndex] * numberOfCentEvents;
          #endif
        }

        #ifdef _CREATE_TEXT_FIT_FILE_
          if(!std::isfinite(spectraValue_forText))    spectraValue_forText = 0.0;
          if(!std::isfinite(spectraValueErr_forText)) spectraValueErr_forText = 999.999;
          outputTextFile_local << Form("%.10e %.10e",spectraValue_forText,spectraValueErr_forText) << endl;
          outputTextFile_local.close();
        #endif

        if(pionFunct)            delete pionFunct;
        if(kaonFunct)            delete kaonFunct;
        if(protonFunct)          delete protonFunct;
        if(lowRangeLine)         delete lowRangeLine;
        if(highRangeLine)        delete highRangeLine;
        if(pullHisto)            delete pullHisto;
        if(oneLine)              delete oneLine;
        if(fitParsTxt)           delete fitParsTxt;
        //delete frame_pull;
        if(topPad)                  delete topPad;
        if(bottomPad)               delete bottomPad;
        if(fittingCanvas)           delete fittingCanvas;
        if(initialFunct[centIndex]) delete initialFunct[centIndex];
        if(initialFunct_invBeta[centIndex]) delete initialFunct_invBeta[centIndex];
      } // end centrality loop


      if( !(status_Minimizer == 0 || status_Minimizer == 4) ){
        #ifdef _ZFITTER_DEBUG_
          cout << " Fit falied with satus: " << status_Minimizer  << endl;
        #endif
        for(int iii = 0; iii < prevParamSize; iii++) previousGoodParameters[iii] = -999;
        previousCentralityIntegratedYields_Pion = 0;
        previousCentralityIntegratedYields_Kaon = 0;
        previousCentralityIntegratedYields_Proton = 0;
        previousFitToggle_Pion   = false;
        previousFitToggle_Kaon   = false;
        previousFitToggle_Proton = false;

      }


      //DRAW IMAGE WITH ALL CENTRALITIES
      TCanvas* fittingCanvas = new TCanvas("fittingCanvas_SimulBTOF","Fitting Canvas");
      fittingCanvas->SetWindowSize(1300,700);
      fittingCanvas->SetCanvasSize(1200,600);
      fittingCanvas->Draw();
      fittingCanvas->cd();
      gPad->SetLogy();
      gPad->SetRightMargin(0.45);
      gStyle->SetOptFit(0111);



      //bottomPad->cd();
      //bottomPad->Clear();
      //topPad->cd();
      //topPad->SetLogy(true);
      m_currentHistosToFit_ByCent[0]->SetStats(kFALSE);
      m_currentHistosToFit_ByCent[0]->Draw("E");
      combinedFunct_CentStorage[0]->Draw("same");
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        #ifdef _ZFITTER_DEBUG_
          cout << " Ptr Check: " << m_currentHistosToFit_ByCent[centIndex] << "  " << combinedFunct_CentStorage[centIndex] << endl;
        #endif
        if(centIndex == 0){
          m_currentHistosToFit_ByCent[centIndex]->Draw("E");
        }else{
          m_currentHistosToFit_ByCent[centIndex]->Draw("E SAME");
        }
        combinedFunct_CentStorage[centIndex]->Draw("same");

      }

      //------------------   Make Text Boxes ---------------------------
      TPaveText* fitParsTxtSimul = new TPaveText(0.65,0.05,0.95,0.6,"NDC");
      for(int iii = 0; iii < 12; iii++){
        fitParsTxtSimul->AddText(Form("%s \t %.2e#pm%.2e %s",paramNames[iii].c_str(),parameters_simul_cent[iii],parameters_simul_cent_errors[iii],(minimizer->IsFixedVariable(iii) ? "(fixed)" : "")));
      }
      fitParsTxtSimul->SetFillStyle(0);
      fitParsTxtSimul->SetBorderSize(5);
      fitParsTxtSimul->Draw("same");

      TPaveText* fitTxtSimul = new TPaveText(0.65,0.60,0.95,0.97,"NDC");
      for(int iii = 0; iii < 12; iii++){
        fitTxtSimul->AddText(Form("%s \t %f  [ %f, %f ]",paramNames[iii].c_str(),initialParams[iii],initialParamLimits[iii][0],initialParamLimits[iii][1]));
      }
      fitTxtSimul->SetFillStyle(0);
      fitTxtSimul->SetBorderSize(5);

      fitTxtSimul->Draw("same");
      gPad->Update();
      gSystem->ProcessEvents();
      fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/SimulCent_RapBin_%02d_mTm0Bin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
           m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin));





      #ifdef _ZFITTER_DEBUG_
        cout << "minimizer: "     << minimizer << endl;
        cout << "chiSqrdFunctor: "<< chiSqrdFunctor << endl;
        cout << "fittingCanvas: " << fittingCanvas << endl;
      #endif
      delete fitTxtSimul;
      delete fitParsTxtSimul;
      delete fittingCanvas;

      //delete chiSqrdFunctor;
      delete minimizer;
      for(int iii = 0; iii < 16; iii++){
        #ifdef _ZFITTER_DEBUG_
          cout << " Cent " << iii << " " << combinedFunct_CentStorage[iii] << " " << m_currentHistosToFit_ByCent[iii] << endl;
        #endif
        if(combinedFunct_CentStorage[iii])   delete combinedFunct_CentStorage[iii];
        if(m_currentHistosToFit_ByCent[iii]) delete m_currentHistosToFit_ByCent[iii];
        if(m_currentHistosToFit_ByCent_InvBeta[iii]) delete m_currentHistosToFit_ByCent_InvBeta[iii];
      }

      for(int iii = 0; iii < 16; iii++){
        for(int jjj = 0; jjj < 3; jjj++){
          if(initialPoints[iii][jjj]) delete initialPoints[iii][jjj];
        }
      }


    }// End mTm0 Loop

    //outputTextFile.close();
    #ifdef _ZFITTER_DEBUG_
      cout << "Finished the mTm0 loop for BTOF Simultaneous Centrality Fitting" << endl;
    #endif

    /*for(int iii = 0; iii < 1000; iii++){
      for(int jjj = 0; jjj < 16; jjj++){
        if(rebinZBTOFHistoStorage[iii][jjj]) delete rebinZBTOFHistoStorage[iii][jjj];
      }
    }*/
  }// End Rapidity Loop

  //delete minimizer;   //commented to see if this is causing the wierd crashes
  //delete chiSqrdFunctor;
  //delete topPad;
  //delete bottomPad;
  //delete fittingCanvas;
  #ifdef _ZFITTER_DEBUG_
    cout << "Finished the Rapidity loop for BTOF Simultaneous Centrality Fitting" << endl;
    cout << "End of Function Call" << endl;
  #endif
  return;
}




//============================================================================================================##################################################################################################################
//============================================================================================================##################################################################################################################

/*
void ZFitter::extrapolate_BTOF_Fit_Parameters_PionKaon(int a_particleIndexToExtrapolate, double a_lowMomFitRange, double a_highMomFitRange, int a_functionVersion){
  //####################################################################################################################################################################################################################################
  //###################    Extrapolate Gaus / Student-T / SkewNormal Parameters [lambda / sigma / m / nu]                                   ###################################################################################################################################
  //####################################################################################################################################################################################################################################
  double lowMtM0  = PhysMath::mTm0_from_y_p(m_currentRapidity,a_lowMomFitRange,m_partInfo->GetParticleMass(m_currentPartIndex));
  double highMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity,a_highMomFitRange,m_partInfo->GetParticleMass(m_currentPartIndex));
  if(lowMtM0 < 0.0){
    cout << " WARNING: the rapidity is too extreme for the extrapolation procedure you have given setting low mTm0 to 0" << endl;
    lowMtM0 = 0.0;
  }
  if(highMtM0 < 0.0){
    cout << " WARNING: the rapidity is too extreme for the extrapolation procedure, becuase high mTm0 is negative... (" << highMtM0 << " ) so continuing on without extrapolation" << endl;
    return;
  }
  double midMtM0 = 0.5*(lowMtM0 + highMtM0);
  #ifdef _ZFITTER_DEBUG_
    cout << "  extrapolating " << m_partInfo->GetParticleName(a_particleIndexToExtrapolate) << " in " << m_partInfo->GetParticleName(m_currentPartIndex)  << "_space with mTm0 fit range from " << lowMtM0 << "  to  " << highMtM0 << endl;
    cout << "      momentum range was " << a_lowMomFitRange << "  to  " << a_highMomFitRange << endl;
  #endif

  //make the histograms to extrapolate
  TH1D* param1DHistos[5];
  TF1* extrapFunctions[5];
  double numberOfParams = 5;
  string distributionName = "Gaussian";
  string paramNames[5] = {"N","mu","sigma","ignore","ignore"};
  string paramSymbols[5] = {"N","#lambda","#sigma","ignore","ignore"};
  if(a_functionVersion == 1) {paramNames[3] = "gamma"; paramSymbols[3] = "#gamma"; distributionName = "SkewNormal";};
  if(a_functionVersion == 2) {paramNames[3] = "nu";    paramSymbols[3] = "#nu"; distributionName = "StudentT";};
  if(a_functionVersion == 3) {paramNames[3] = "m";    paramSymbols[3] = "m";paramNames[4] = "nu";    paramSymbols[4] = "#nu"; distributionName = "PearsonIV";};

  string functionNames[5];
  if(m_inAllCentMode){
    for(int jjj = 1; jjj < numberOfParams; jjj++){
      param1DHistos[jjj] = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][jjj][m_currentPlusMinusIndex]->ProjectionY(Form("param_proj_%s_%s_%s_Rap%02d",
         paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
      functionNames[jjj] = Form("bTOF_param_parab_%s_%s_%s_Rap%02d",paramNames[jjj].c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
        m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),m_currentRapBin);
      param1DHistos[jjj]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),paramSymbols[jjj].c_str()));
    }
  }else{
    for(int jjj = 1; jjj < numberOfParams; jjj++){
      param1DHistos[jjj] = m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][m_currentCentIndex][jjj][m_currentPlusMinusIndex]->ProjectionY(Form("param_proj_%s_%s_%s_Cent%02d_Rap%02d",
         paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
      functionNames[jjj] = Form("bTOF_param_parab_%s_%s_%s_Cent%02d_Rap%02d",paramNames[jjj].c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
        m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin);
      param1DHistos[jjj]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),paramSymbols[jjj].c_str()));
    }
  }

  int lowBin = param1DHistos[1]->FindBin(lowMtM0);
  int midBin = param1DHistos[1]->FindBin(midMtM0);
  int highBin =param1DHistos[1]->FindBin(highMtM0);

  bool properLowBinFound = true;
  if(param1DHistos[1]->GetBinError(lowBin) == 0.0){
    properLowBinFound = false;
    for(int binX = lowBin; binX < highBin; binX++){
      if(param1DHistos[1]->GetBinError(binX) > 0.0){
        properLowBinFound = true;
        lowBin = binX;
        break;
      }
    }
  }

  bool properHighBinFound = true;
  if(param1DHistos[1]->GetBinError(highBin) == 0.0){
    properHighBinFound = false;
    for(int binX = highBin; binX > lowBin; binX--){
      if(param1DHistos[1]->GetBinError(binX) > 0.0){
        properHighBinFound = true;
        highBin = binX;
        break;
      }
    }
  }

  if(!properLowBinFound || !properHighBinFound){
    cout << " ERROR: The Extrapolation will not work. Could not find bins with actual fit data..." << endl;
  }
  if(midBin <= lowBin || midBin >=highBin){
    midBin = (lowBin+highBin)/2;
    cout << " WARNING: Mid-Bin Didn't work, so putting it to half the range" << endl;
  }
  if(highBin - lowBin < 2){
    cout << " ERROR: The extrapolation will fail because you only have two bins to extrapolate..." << endl;
  }
  lowMtM0  =  param1DHistos[1]->GetBinCenter(lowBin)  - 0.25*param1DHistos[1]->GetBinWidth(lowBin);
  highMtM0 = param1DHistos[1]->GetBinCenter(highBin) + 0.25*param1DHistos[1]->GetBinWidth(highBin);
  midMtM0  = 0.5*(lowMtM0 + highMtM0);
  midBin   = param1DHistos[1]->FindBin(midMtM0);
  midMtM0  = param1DHistos[1]->GetBinCenter(midBin);

  #ifdef _ZFITTER_DEBUG_
    cout << " Resized the fit range... lowMtM0: " << lowMtM0 << " midMtM0: " << midMtM0 << " highMtM0: " << highMtM0 << endl;
    cout <<  "                   LowBin: " << lowBin << " midBin: " << midBin << " highBin: " << highBin << endl;
  #endif


  // ########    PARAMETER LOOP FOR MEAN WIDTH AND MAYBE NDF ###########
  for(int jjj = 1; jjj < numberOfParams; jjj++){
    double parabolaParams[3];
    double amplitudes[3];
    amplitudes[0] = param1DHistos[jjj]->GetBinContent(lowBin);
    amplitudes[1] = param1DHistos[jjj]->GetBinContent(midBin);
    amplitudes[2] = param1DHistos[jjj]->GetBinContent(highBin);
    if(amplitudes[0] == 0.0 && amplitudes[1] == 0.0 &&amplitudes[2] == 0.0){
      #ifdef _ZFITTER_DEBUG_
        cout << " all the amplitudes at the three points are zero... continuing" << endl;
      #endif
      continue;
    }

    double slope_1 = (amplitudes[1] - amplitudes[0])/(midMtM0 - lowMtM0);
    double slope_2 = (amplitudes[2] - amplitudes[1])/(highMtM0 - midMtM0);
    double ySlope_0 = 0.5* (slope_1 + slope_2);
    double maxVal = HistogramUtilities::getMaxInRange(param1DHistos[jjj], lowMtM0, highMtM0);
    double minVal = HistogramUtilities::getMinInRange(param1DHistos[jjj], lowMtM0, highMtM0);
    double distanceScale = maxVal - minVal;
    double slopeScale = distanceScale / (highMtM0 - lowMtM0);

    #ifdef _ZFITTER_DEBUG_
      cout << " Initializing linear funciton with y0=" << amplitudes[0] << " y'= " << ySlope_0 << endl;
    #endif

    extrapFunctions[jjj] = PhysMath::linearFunct(functionNames[jjj].c_str(), amplitudes[0], ySlope_0,lowMtM0);
    extrapFunctions[jjj]->SetRange(lowMtM0,highMtM0);

    if(a_functionVersion == 3 && jjj > 2){
      extrapFunctions[jjj]->SetParLimits(0,amplitudes[0] - 10.0*distanceScale,amplitudes[0] + 10.0*distanceScale);
      extrapFunctions[jjj]->FixParameter(1,0.0);
    }else{
      //extrapFunctions[jjj]->SetParLimits(0,amplitudes[0] - 10.0*distanceScale,amplitudes[0] + 10.0*distanceScale);
      extrapFunctions[jjj]->FixParameter(0,amplitudes[0]);
      extrapFunctions[jjj]->SetParLimits(1,-10.0*slopeScale,10.0*slopeScale);
    }


    TF1* initialFunct = (TF1*) extrapFunctions[jjj]->Clone();
    initialFunct->SetLineColor(kYellow);
    initialFunct->SetName("initialFunct");

    string fitOptions = "";
    #ifdef _ZFITTER_DEBUG_
      fitOptions += " V";
    #else
      fitOptions += " Q";
    #endif

    int status = param1DHistos[jjj]->Fit(functionNames[jjj].c_str(),fitOptions.c_str(),"",lowMtM0,highMtM0);
    TF1* fitResult = param1DHistos[jjj]->GetFunction(functionNames[jjj].c_str());

    //Push the results to the histograms
    //double mu = -999.999;
    if(fitResult){
      for(int binX = lowBin+1; binX <= param1DHistos[jjj]->GetNbinsX(); binX++){
        double center = param1DHistos[jjj]->GetBinCenter(binX);
        double valueToPush =  fitResult->Eval(center);
        if(m_inAllCentMode){
          m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][jjj][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, binX,valueToPush);
          m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][jjj][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin, binX, 1.0E-13);
        }else{
          m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][m_currentCentIndex][jjj][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, binX,valueToPush);
          m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][a_particleIndexToExtrapolate][m_currentCentIndex][jjj][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin, binX, 1.0E-13);
        }
      }
    }else{
      #ifdef _ZFITTER_DEBUG_
        cout << " ERROR: PION PARAMETER " << jjj << " EXTRAPOLATION FAILED!!!!!!!!! " << endl;
        HistogramUtilities::printHistoInfo(param1DHistos[jjj]);
      #endif
    }

    HistogramUtilities::discardNanHistoData(param1DHistos[jjj]);
    TGraph* parabPoint = new TGraph();
    if(fitResult) parabPoint->SetPoint(0,midMtM0,fitResult->GetParameter(0));
    parabPoint->SetMarkerStyle(20);
    parabPoint->SetMarkerColor(kRed);


    // DRAW THE FITS
    TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
    fittingCanvas->SetWindowSize(1300,700);
    fittingCanvas->SetCanvasSize(1200,600);
    fittingCanvas->cd();
    gStyle->SetOptFit(0111);
    HistogramUtilities::setAxisRanges(param1DHistos[jjj], lowMtM0,highMtM0,false);
    /*TLine* extraLine = NULL;
    if(mu != 999.999 && fitResult && mu < highMtM0 && mu > lowMtM0){
      extraLine = new TLine(mu,fitResult->Eval(mu), highMtM0, fitResult->Eval(mu));
      extraLine->SetLineColor(kRed);
    }* /
    param1DHistos[jjj]->SetStats(true);
    param1DHistos[jjj]->Draw("E");
    initialFunct->Draw("same");
    parabPoint->Draw("Psame");
    //if(extraLine) extraLine->Draw("same");
    //gPad->SetLogy();
    gPad->SetRightMargin(0.4);
    if(m_inAllCentMode){
      fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Param%02d_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),jjj,
      m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Param%02d_%s_in_%s_space_Cent%02d_Rap_%02d.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),jjj,
      m_partInfo->GetParticleName(a_particleIndexToExtrapolate,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin));
    }

    cout << " param1DHistos[jjj] " << param1DHistos[jjj] << endl;   delete param1DHistos[jjj];
    cout << " initialFunct " << initialFunct << endl;   delete initialFunct;
    cout << " extrapFunctions[jjj] " << extrapFunctions[jjj] << endl;   delete extrapFunctions[jjj];
    cout << " parabPoint " << parabPoint << endl;   delete parabPoint;
    cout << " fittingCanvas " << fittingCanvas << endl;   delete fittingCanvas;
    //cout << " fitResult " << fitResult << endl; if(fitResult) delete fitResult;

  } // jjj parameter loop
}

*/
