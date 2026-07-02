// THIS IS INSERTED INTO ZFITTER.cxx // Just separate file becuase it's long
//============================================================================================================
//#########      TPC FITTING FOR PION KAON ONLY   #####################
//============================================================================================================
#ifndef _ZFITTER_PIONKAONTPC_SIMULCENT_SINGLELOOP_
#define _ZFITTER_PIONKAONTPC_SIMULCENT_SINGLELOOP_

#include "../headers/ZFitter.h"

// The two minimizer chi-squared callbacks below (simultaneous_centrality_ZTPC_
// chisqr_E_K_Ratio_Linked / _E_Pi_Ratio_Linked) are copied verbatim from the
// original repo's source/ZFitter_PionKaonTPC_SimulCent.cxx -- a file this port
// otherwise excludes because its own fitTPCPionKaon_SimulCent/_Part2 driver methods
// are superseded by fitTPCPionKaon_SimulCent_RapidityLooper below. However, the
// RapidityLooper driver itself calls BOTH of these functors (selecting between them
// at runtime via the doElectronPionRatio toggle -- see the `if(doElectronPionRatio)`
// branch a bit further down in this file), so they are genuinely live dependencies,
// not dead code, and had to be pulled out of the otherwise-excluded file rather than
// left behind. Verified via full-repo grep that neither function is defined anywhere
// else that this port also excludes/keeps, so there is exactly one copy of each.
double ZFitter::simultaneous_centrality_ZTPC_chisqr_E_K_Ratio_Linked(const double* a_param){
  //    mu_pi sigma_pi gamma_pi    mu_K sigma_K    mu_e sigma_e then For Each Cent: A_pi A_k A_e ... A_pi A_k A_e
  double chisqr = 0.0;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    for(int binX = m_currentLowFitBin[centIndex]; binX <= m_currentHighFitBin[centIndex]; binX++){
      double center  = m_currentHistosToFit_ByCent[centIndex]->GetBinCenter(binX);
      double content = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(binX);
      if(content > 0.0){
        double prediction =   a_param[7 + (3*centIndex)]*PhysMath::skewNormal(a_param[0],a_param[1],a_param[2],center)
                            + a_param[8 + (3*centIndex)]*PhysMath::gaus(a_param[3],a_param[4],center)
                            + a_param[8 + (3*centIndex)]*a_param[9 + (3*centIndex)]*PhysMath::gaus(a_param[5],a_param[6],center);
        chisqr += pow((prediction - content)/m_currentHistosToFit_ByCent[centIndex]->GetBinError(binX),2);
      }
    }
  }
  //#ifdef _ZFITTER_DEBUG_
  //  cout << "simultaneous_centrality_ZTPC_chisqr: " << chisqr << endl;
  //#endif
  return chisqr;
}

double ZFitter::simultaneous_centrality_ZTPC_chisqr_E_Pi_Ratio_Linked(const double* a_param){
  //    mu_pi sigma_pi gamma_pi    mu_K sigma_K    mu_e sigma_e then For Each Cent: A_pi A_d A_e ... A_pi A_k A_e
  double chisqr = 0.0;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    for(int binX = m_currentLowFitBin[centIndex]; binX <= m_currentHighFitBin[centIndex]; binX++){
      double center  = m_currentHistosToFit_ByCent[centIndex]->GetBinCenter(binX);
      double content = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(binX);
      if(content > 0.0){
        double prediction =   a_param[7 + (3*centIndex)] * PhysMath::skewNormal(a_param[0],a_param[1],a_param[2],center)
                            + a_param[8 + (3*centIndex)] * PhysMath::gaus(a_param[3],a_param[4],center)
                            + a_param[7 + (3*centIndex)] * a_param[9 + (3*centIndex)] * PhysMath::gaus(a_param[5],a_param[6],center);
        chisqr += pow((prediction - content)/m_currentHistosToFit_ByCent[centIndex]->GetBinError(binX),2);
      }
    }
  }
  //#ifdef _ZFITTER_DEBUG_
  //  cout << "simultaneous_centrality_ZTPC_chisqr: " << chisqr << endl;
  //#endif
  return chisqr;
}

void ZFitter::fitTPCPionKaon_SimulCent_RapidityLooper(int a_partIndexSpace, int a_charge){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::fitTPCPionKaon_SimulCent_RapidityLooper(" << a_partIndexSpace << ", " << a_charge << ")" << endl;
  #endif
  m_currentCentIndex  = -1;
  m_inAllCentMode     = true;
  m_currentPartIndex  = a_partIndexSpace;
  m_currentPartCharge = a_charge;
  if(a_charge > 0) m_currentPlusMinusIndex = 0;
  else             m_currentPlusMinusIndex = 1;

  TH3D* histo3DPM = NULL;
  if(m_currentPartCharge == 1){ // Positive
    histo3DPM = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex];
  }else{ // Negative
    histo3DPM = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex];
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
  bool leftBinExists  = true;

  #ifdef _ZFITTER_DEBUG_
    cout << "using zero bin: " << rapZeroBin << endl;
  #endif

  //###########################   RAPIDITY BIN LOOP   #######################################
  for(int HorBinCounter = 0; rightBinExists && leftBinExists; HorBinCounter++){
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

    #ifdef _ZFITTER_DEBUG_
      cout << " Doing RapBin:  " << m_currentRapBin<< endl;
    #endif

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
    fitTPCPionKaon_SimulCent_ByRapidity();
  }
}





//============================================================================================================================================================================


bool ZFitter::calculateAverageElectronWidth(){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::calculateAverageElectronWidth()" << endl;
  #endif
  TH2D*  electronSTDDevHisto2D = m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex];
  double electronStdDevAverage = HistogramUtilities::meanOfZValsFromProj(electronSTDDevHisto2D, m_currentRapBin, 0, 1);
  #ifdef _ZFITTER_DEBUG_
    cout << "Electron StdDev average: " << electronStdDevAverage << " from mTm0 = [" << 0 << ", " << 1 << "]" << endl;
  #endif
  //for(int binY = 1; binY <= electronSTDDevHisto2D->GetNbinsY(); binY++){
  //  m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, binY,electronStdDevAverage);
  //  m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin, binY, 1.0E-6);
  //}// bin loop
  //return;

  TH1D* electronSTDDevHisto = electronSTDDevHisto2D->ProjectionY(Form("elec_param_proj_%s_%s_%s_Rap%02d",
     "StdDev", m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(3,m_currentPartCharge).Data(),m_currentRapBin),
     m_currentRapBin,m_currentRapBin);
  electronSTDDevHisto = (TH1D*) electronSTDDevHisto->Clone();
  string electronExtrapFunctName = Form("TPC_param_parab_%s_%s_%s_Rap%02d","StdDev",
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(3,m_currentPartCharge).Data(),m_currentRapBin);
  electronSTDDevHisto->SetTitle("Electron Z_{TPC} #sigma");
  electronSTDDevHisto->GetYaxis()->SetTitle("#sigma");
  electronSTDDevHisto->GetXaxis()->SetTitle("m_{T}-m_{0}");
  HistogramUtilities::discardNanHistoData(electronSTDDevHisto,true);

  #ifdef _ZFITTER_DEBUG_
    HistogramUtilities::printHistoInfo(electronSTDDevHisto);
  #endif

  //#################   Determine if there is a low pT region fit  ################
  bool electronLowPtWasFit = false;
  bool firstContentFound   = false;
  bool firstHoleFound      = false;

  double  mTm0OfFirstContent = -1;
  double  mTm0OfSecondContent = -1;

  double weights_lowPt  = 0.0;
  double average_lowPt  = 0.0;
  double weights_highPt = 0.0;
  double average_highPt = 0.0;

  int numLowPtPoints     = 0;
  int numHighPtPoints    = 0;
  for(int bin = 1; bin <= electronSTDDevHisto->GetNbinsX(); bin++){
    double content = electronSTDDevHisto->GetBinContent(bin);
    if(firstContentFound){
      if(firstHoleFound){
        if(content > 0 && std::isfinite(content)){
          electronLowPtWasFit = true;
          numHighPtPoints++;
          double weight  = 1.0/(electronSTDDevHisto->GetBinError(bin)*electronSTDDevHisto->GetBinError(bin));
          double content = electronSTDDevHisto->GetBinContent(bin);
          weights_highPt += weight;
          average_highPt += content*weight;
          if(mTm0OfSecondContent < 0.0) mTm0OfSecondContent = electronSTDDevHisto->GetBinCenter(bin);
        }
      }else{
        if(content <= 0 && std::isfinite(content)){
          firstHoleFound = true;
        }else{
          numLowPtPoints++;
          double weight  = 1.0/(electronSTDDevHisto->GetBinError(bin)*electronSTDDevHisto->GetBinError(bin));
          double content = electronSTDDevHisto->GetBinContent(bin);
          weights_lowPt += weight;
          average_lowPt += content*weight;
        }
      }
    }else{
      if(content > 0 && std::isfinite(content)){
        firstContentFound = true;
        mTm0OfFirstContent = electronSTDDevHisto->GetBinCenter(bin);
        numLowPtPoints++;
        double weight  = 1.0/(electronSTDDevHisto->GetBinError(bin)*electronSTDDevHisto->GetBinError(bin));
        double content = electronSTDDevHisto->GetBinContent(bin);
        weights_lowPt += weight;
        average_lowPt += content*weight;
      }
    }
  }

  if(!electronLowPtWasFit){
    numHighPtPoints = numLowPtPoints;
    numLowPtPoints = 0;
    weights_highPt = weights_lowPt;
    average_highPt = average_lowPt;
    weights_lowPt = 0.0;
    average_lowPt = 0.0;
    mTm0OfFirstContent = 0.0;
  }

  average_lowPt  = average_lowPt/weights_lowPt;
  average_highPt = average_highPt/weights_highPt;


  int lastFitSigmaBin = HistogramUtilities::lastBinWithContent(electronSTDDevHisto);
  double lastFitMtM0 = electronSTDDevHisto->GetBinLowEdge(lastFitSigmaBin);
  double lastFitSigmaWidth = 0.0;
  //######### modify the std-dev histogram by the spread of electrons in the bin ######
  for(int bin = 1; bin <= electronSTDDevHisto->GetNbinsX(); bin++){
    double content = electronSTDDevHisto->GetBinContent(bin);
    if(content == 0.0) continue;
    double lowMtM0  = electronSTDDevHisto->GetBinLowEdge(bin);
    double highMtM0 = electronSTDDevHisto->GetBinLowEdge(bin+1);
    double lowMomentum  = PhysMath::pT_from_mTm0(lowMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
    double highMomentum = PhysMath::pT_from_mTm0(lowMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
    double width = fabs(log(m_bichselCurves[3]->Eval(highMomentum)/m_bichselCurves[3]->Eval(lowMomentum)));
    if(lastFitSigmaBin == bin) lastFitSigmaWidth = width;
    electronSTDDevHisto->SetBinContent(bin, content - width);
  }




  #ifdef _ZFITTER_DEBUG_
    cout << "  Two Electron Regions? " << (electronLowPtWasFit ? "YES" : "NO") << endl;
    cout << "  Num Low Pt Points  : " << numLowPtPoints << endl;
    cout << "  Num High Pt Points : " << numHighPtPoints << endl;
  #endif

  if(numHighPtPoints < 3){
    cout << "WARNING: calculateAverageElectronWidth() can't determing widths due to insufficient fit bins:  " << numHighPtPoints << endl;
    return false;
  }

  //TF1* extrapFunct = PhysMath::getExponentialFunction(electronExtrapFunctName.c_str(), 0.04, -0.2, 0, 0.067);
  //TF1* extrapFunct = new TF1(electronExtrapFunctName.c_str(),"[0]*exp(-x/[1])/[1] + [2] + [3]*x",0,5.0);
  //extrapFunct->SetNpx(500);
  //extrapFunct->SetLineColor(kRed);
  //extrapFunct->SetParNames("N","T","y_0","a");
  //extrapFunct->SetRange(0,1);
  //extrapFunct->SetParameter(0,0.0);
  //extrapFunct->SetParameter(1,20.0);
  //extrapFunct->SetParameter(2,0.067);
  //extrapFunct->SetParameter(3,0.0);
  //extrapFunct->SetParLimits(0,0.0,1.0);
  //extrapFunct->SetParLimits(1,0.0001,200);
  //extrapFunct->SetParLimits(2,-0.1,0.1);
  //extrapFunct->SetParLimits(3,-0.1,0.5);

  //TF1* extrapFunct = new TF1(electronExtrapFunctName.c_str(),"(x < ([3]-[1])/([0]-[2]))*([0]*x + [1]) + (x >= ([3]-[1])/([0]-[2]))*([2]*x + [3])",0,5.0);
  //extrapFunct->SetNpx(500);
  //extrapFunct->SetLineColor(kRed);
  //extrapFunct->SetParNames("A","B","C","D");
  //extrapFunct->SetParameter(0,-0.01);
  //extrapFunct->SetParameter(1,0.067);
  //extrapFunct->SetParameter(2,0.0);
  //extrapFunct->SetParameter(3,0.067);

  double crossXInit = 0.5*(mTm0OfSecondContent + mTm0OfFirstContent);
  double ampInit    = (average_lowPt - average_highPt)/pow(mTm0OfFirstContent - crossXInit,2);
  TF1* extrapFunct = new TF1(electronExtrapFunctName.c_str(),"(x < [1])*[0]*pow(x - [1],2) + (x > [1])*[3]*(x - [1]) + [2]",0,5.0);
  extrapFunct->SetNpx(500);
  extrapFunct->SetLineColor(kRed);
  extrapFunct->SetParNames("A","x_0","y_0","m");
  extrapFunct->SetParameter(0,ampInit);
  extrapFunct->SetParLimits(0,-1.0,20.0);
  extrapFunct->SetParameter(1,crossXInit);
  extrapFunct->SetParLimits(1,0.0,0.5);
  extrapFunct->SetParameter(2,average_highPt);
  extrapFunct->SetParLimits(2,0.0,0.3);
  extrapFunct->SetParameter(3,0.001);
  extrapFunct->SetParLimits(3,-1.0,5.0);

  string initialFunctName = electronExtrapFunctName + "Initial";
  TF1* initialFunct = (TF1*) extrapFunct->Clone();
  initialFunct->SetLineColor(kYellow+1);
  initialFunct->SetName(initialFunctName.c_str());

  //if low pT isn't there, you should just use a line
  if(!electronLowPtWasFit){
    extrapFunct->FixParameter(0,0.0);
    extrapFunct->FixParameter(1,0.0);
  }

  if(numHighPtPoints <= 3){
    extrapFunct->FixParameter(3,0.0);
  }


  string fitOptions = "S";
  #ifdef _ZFITTER_DEBUG_
    fitOptions += "";
  #else
    fitOptions += "Q";
  #endif
  #ifdef _ZFITTER_CRASH_BYPASS_
    TF1* fitResult = NULL;
    bool doBypass = true;
    if(doBypass){
      fitResult = (TF1*) initialFunct->Clone();
    }else{
      TFitResultPtr resultPtr = electronSTDDevHisto->Fit(extrapFunct,fitOptions.c_str(),"",0,1);
      fitResult = electronSTDDevHisto->GetFunction(electronExtrapFunctName.c_str());
    }
  #else
    TFitResultPtr resultPtr = electronSTDDevHisto->Fit(extrapFunct,fitOptions.c_str(),"",0,1);
    //electronSTDDevHisto->Fit("pol2","","",0,1);
    TF1* fitResult = electronSTDDevHisto->GetFunction(electronExtrapFunctName.c_str());
  #endif

  if(fitResult){
    double lastFitSigma = fitResult->Eval(lastFitMtM0) + lastFitSigmaWidth;
    //######### modify the std-dev histogram by the spread of electrons in the bin ######
    for(int bin = 1; bin <= electronSTDDevHisto->GetNbinsX(); bin++){
      double content = electronSTDDevHisto->GetBinContent(bin);
      //if(content == 0.0) continue;
      double lowMtM0  = electronSTDDevHisto->GetBinLowEdge(bin);
      double highMtM0 = electronSTDDevHisto->GetBinLowEdge(bin+1);
      double lowMomentum  = PhysMath::pT_from_mTm0(lowMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      double highMomentum = PhysMath::pT_from_mTm0(lowMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      double width = fabs(log(m_bichselCurves[3]->Eval(highMomentum)/m_bichselCurves[3]->Eval(lowMomentum)));

      double mTm0 = electronSTDDevHisto2D->GetYaxis()->GetBinCenter(bin);
      double valueToPush = fitResult->Eval(mTm0) + width;
      if(mTm0 > lastFitMtM0 && valueToPush > lastFitSigma) valueToPush = lastFitSigma;
      if(valueToPush < 0.065) valueToPush = 0.065;
      if(valueToPush > 0.3) valueToPush = 0.3;
      m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, bin,valueToPush);
      m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin, bin, 1.0E-6);
    }
    /*
    for(int binY = 1; binY <= electronSTDDevHisto2D->GetNbinsY(); binY++){
      double mTm0 = electronSTDDevHisto2D->GetYaxis()->GetBinCenter(binY);
      m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, binY,fitResult->Eval(mTm0));
      m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin, binY, 1.0E-6);
    }// bin loop
    */
  }

  // DRAW THE FITS
  TCanvas* fittingCanvas = new TCanvas(Form("fittingCanvasStdElec%s_%02d", m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin),"Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  HistogramUtilities::setAxisRanges(electronSTDDevHisto, 0,0.5,true);
  electronSTDDevHisto->SetStats(true);
  electronSTDDevHisto->Draw("E");
  initialFunct->Draw("same");
  gPad->SetRightMargin(0.4);

  fittingCanvas->SaveAs(Form("%s/%s/%s/TPCParamExtrapolations/electronStdDev_in_%s_space_Rap_%02d.png", m_imagePreDir.c_str(),
      m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin));


  #ifndef _SAVE_POINTERS_
    delete extrapFunct;
    delete initialFunct;
    //if(fitResult)     delete fitResult;
    delete electronSTDDevHisto;
    delete fittingCanvas;
  #endif

  return true;
}


//============================================================================================================================================================================



bool ZFitter::calculateElectronRatios(){

  //############     FIT THE ELECTRON AMPLITUDE (or e/k ratio) WITH AN EXPONENTIAL  ###############
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::calculateElectronRatios()" << endl;
    cout << "Fitting e/pi and e/k ratio with an exponential for electron extrapolation...." << endl;
  #endif
  int lowBin  = HistogramUtilities::firstYBinWithContent(m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex], m_currentRapBin);
  int highBin = HistogramUtilities::lastYBinWithContent(m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex], m_currentRapBin);
  double lowMtM0  = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex]->GetYaxis()->GetBinLowEdge(lowBin);
  double highMtM0 = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex]->GetYaxis()->GetBinLowEdge(highBin);
  double midMtM0 = 0.5*(lowMtM0 + highMtM0);
  int midBin  = m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->GetYaxis()->FindBin(midMtM0);
  bool isEmpty = HistogramUtilities::isYSliceEmpty(m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex], m_currentRapBin);
  #ifdef _ZFITTER_DEBUG_
    cout << "  mTm0 Bin = [" << lowBin << ", " << midBin << ", " << highBin << "]" << endl;
    cout << "  mTm0 Val = [" << lowMtM0 << ", " << midMtM0 << ", " << highMtM0 << "]" << endl;
    cout << "  Is Histo Slice Empty?    " << (isEmpty ? "YES" : "NO") << endl;
  #endif

  if(highBin - lowBin < 2 || isEmpty){
    cout << "WARNING: calculateElectronRatios can't complete due to too few bins... continueing..." << endl;
    return false;
  }

  TF1*  elecPionRatioLowMtM0ExtrapLine = NULL;

  TH1D* electrondNdMtHisto_ByCent[16];
  TH1D* electronPionRatioHisto_ByCent[16];
  TH1D* electronKaonRatioHisto_ByCent[16];

  string electronExtrapFunctName;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    electrondNdMtHisto_ByCent[centIndex] = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][0][m_currentPlusMinusIndex]->ProjectionY(Form("param_proj_%s_%s_%s_Rap%02d_cent%02d",
       "dNdmT", m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(3,m_currentPartCharge).Data(),m_currentRapBin,centIndex),
       m_currentRapBin,m_currentRapBin);
    electronExtrapFunctName = Form("TPC_param_parab_%s_%s_%s_Rap%02d_Cent%02d","dNdmT",
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(3,m_currentPartCharge).Data(),m_currentRapBin,centIndex);
    electrondNdMtHisto_ByCent[centIndex]->GetYaxis()->SetTitle("Electron dN/dmT; m_{T}-m_{0};dN/dm_{T}");

    TH1D* structHisto = (TH1D*) electrondNdMtHisto_ByCent[centIndex]->Clone();
    HistogramUtilities::zeroBins(structHisto);
    electronPionRatioHisto_ByCent[centIndex] = (TH1D*) structHisto->Clone();
    electronPionRatioHisto_ByCent[centIndex]->GetYaxis()->SetTitle("N_{e}/N_{#pi}");
    electronKaonRatioHisto_ByCent[centIndex] = (TH1D*) structHisto->Clone();
    electronKaonRatioHisto_ByCent[centIndex]->GetYaxis()->SetTitle("N_{e}/N_{K}");
    #ifndef _SAVE_POINTERS_
      delete structHisto;
    #endif

    for(int mTm0Bin = 1; mTm0Bin <= electrondNdMtHisto_ByCent[centIndex]->GetNbinsX(); mTm0Bin++){
      double mTm0           = electrondNdMtHisto_ByCent[centIndex]->GetBinCenter(mTm0Bin);
      double electronAmp    = electrondNdMtHisto_ByCent[centIndex]->GetBinContent(mTm0Bin);
      double electronAmp_Err= electrondNdMtHisto_ByCent[centIndex]->GetBinError(mTm0Bin);
      double kaonAmp        = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,mTm0Bin);
      double kaonAmp_Err    = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,mTm0Bin);
      double pionAmp        = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,mTm0Bin);
      double pionAmp_Err    = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,mTm0Bin);

      double percentError_KR = sqrt(electronAmp_Err*electronAmp_Err/(electronAmp*electronAmp) + kaonAmp_Err*kaonAmp_Err/(kaonAmp*kaonAmp));
      double percentError_PR = sqrt(electronAmp_Err*electronAmp_Err/(electronAmp*electronAmp) + pionAmp_Err*pionAmp_Err/(pionAmp*pionAmp));
      double KR = electronAmp/kaonAmp;
      double PR = electronAmp/pionAmp;
      if(!std::isfinite(PR) || !std::isfinite(percentError_PR)){
        PR = 0.0;
        percentError_PR = 0.0;
      }
      if(!std::isfinite(KR) || !std::isfinite(percentError_KR)){
        KR = 0.0;
        percentError_KR = 0.0;
      }
      electronPionRatioHisto_ByCent[centIndex]->SetBinContent(mTm0Bin,PR);
      electronKaonRatioHisto_ByCent[centIndex]->SetBinContent(mTm0Bin,KR);
      electronPionRatioHisto_ByCent[centIndex]->SetBinError(mTm0Bin,PR*percentError_PR);
      electronKaonRatioHisto_ByCent[centIndex]->SetBinError(mTm0Bin,KR*percentError_KR);
      #ifdef _ZFITTER_DEBUG_
        if(std::isfinite(percentError_PR) && PR*KR != 0.0){
          cout << " For mTm0 bin : " << mTm0Bin << " mTm0: " << mTm0 << "   E: " << electronAmp << " Pi: " << pionAmp << " K: " << kaonAmp << " ratios: " << PR << " " << KR << endl;
        }
      #endif
    }
    electronPionRatioHisto_ByCent[centIndex]->SetName(Form("electronPionRatio_cent%02d_rap%03d",centIndex,m_currentRapBin));
    electronKaonRatioHisto_ByCent[centIndex]->SetName(Form("electronKaonRatio_cent%02d_rap%03d",centIndex,m_currentRapBin));
    HistogramUtilities::discardNanHistoData(electronPionRatioHisto_ByCent[centIndex],true);
    HistogramUtilities::discardNanHistoData(electronKaonRatioHisto_ByCent[centIndex],true);

  } // end pseudo-cent loop for simulataneous cent


  //############ Now Fit the e/K and e/pi ratio histograms ################
  TH1D* electronRatioHistos[3]    = {NULL,NULL,NULL};
  TH1D* electronRatioFitHistos[3] = {NULL,NULL,NULL};

  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    electronRatioHistos[0] = electrondNdMtHisto_ByCent[centIndex];
    electronRatioHistos[1] = electronPionRatioHisto_ByCent[centIndex];
    electronRatioHistos[2] = electronKaonRatioHisto_ByCent[centIndex];
    #ifndef _SAVE_POINTERS_
      if(m_electronRatioFitHistos_ByCent[centIndex][0]) delete m_electronRatioFitHistos_ByCent[centIndex][0];
      if(m_electronRatioFitHistos_ByCent[centIndex][1]) delete m_electronRatioFitHistos_ByCent[centIndex][1];
      if(m_electronRatioFitHistos_ByCent[centIndex][2]) delete m_electronRatioFitHistos_ByCent[centIndex][2];
    #endif
    m_electronRatioFitHistos_ByCent[centIndex][0] = NULL;
    m_electronRatioFitHistos_ByCent[centIndex][1] = NULL;
    m_electronRatioFitHistos_ByCent[centIndex][2] = NULL;
    electronRatioHistos[0]->Scale(1.0,"width"); // fit is dN/dmT
    #ifdef _ZFITTER_DEBUG_
      cout << "  extrapolating ELECTRON dN/dmT" << " with mTm0 fit range from " << lowMtM0 << "  to  " << highMtM0 << endl;
    #endif

    for(int electronRatioIndex = 1; electronRatioIndex < 3; electronRatioIndex++){
      m_currentHistoToFit = electronRatioHistos[electronRatioIndex];
      bool properLowBinFound = true;
      if(m_currentHistoToFit->GetBinError(lowBin) == 0.0){
        properLowBinFound = false;
        for(int binX = lowBin; binX < highBin; binX++){
          if(m_currentHistoToFit->GetBinError(binX) > 0.0){
            properLowBinFound = true;
            lowBin = binX;
            break;
          }
        }
      }

      bool properHighBinFound = true;
      if(m_currentHistoToFit->GetBinError(highBin) == 0.0){
        properHighBinFound = false;
        for(int binX = highBin; binX > lowBin; binX--){
          if(m_currentHistoToFit->GetBinError(binX) > 0.0){
            properHighBinFound = true;
            highBin = binX;
            break;
          }
        }
      }

      if(!properLowBinFound || !properHighBinFound){
        cout << " The Extrapolation will not work. Could not find bins with actual fit data..." << endl;
      }
      if(midBin <= lowBin || midBin >=highBin){
        midBin = (lowBin+highBin)/2;
      }
      if(highBin - lowBin < 2){
        cout << " The extrapolation will fail because you only have two bins to extrapolate..." << endl;
      }
      lowMtM0  = m_currentHistoToFit->GetBinCenter(lowBin)  - 0.25*m_currentHistoToFit->GetBinWidth(lowBin);
      highMtM0 = m_currentHistoToFit->GetBinCenter(highBin) + 0.25*m_currentHistoToFit->GetBinWidth(highBin);
      midMtM0  = 0.5*(lowMtM0 + highMtM0);
      midBin   = m_currentHistoToFit->FindBin(midMtM0);
      midMtM0  = m_currentHistoToFit->GetBinCenter(midBin);

      double amplitudes[3];
      amplitudes[0] = m_currentHistoToFit->GetBinContent(lowBin);
      amplitudes[1] = m_currentHistoToFit->GetBinContent(midBin);
      amplitudes[2] = m_currentHistoToFit->GetBinContent(highBin);
      double bSeed   = log(amplitudes[2]/amplitudes[0])/(highMtM0 - lowMtM0);
      double ampSeed = amplitudes[0]*exp(-bSeed*lowMtM0);

      #ifdef _ZFITTER_DEBUG_
        cout << " bSeed: " << bSeed << "    amp seed: " << ampSeed << endl;
      #endif

      TF1* ampExtractFunct = PhysMath::getExponentialFunction(electronExtrapFunctName.c_str(), ampSeed, bSeed, 0, 0);
      ampExtractFunct->SetRange(lowMtM0,highMtM0);
      ampExtractFunct->SetParLimits(0,0.0,10.0*ampSeed);
      if(bSeed < 0.0){
        ampExtractFunct->SetParLimits(1,15.0*bSeed,-10.0*bSeed);
      }else{
        ampExtractFunct->SetParLimits(1,-10.0*bSeed,15.0*bSeed);
      }
      ampExtractFunct->FixParameter(2,0.0);
      ampExtractFunct->SetParLimits(3,-5.0,5.0);

      double lowMtM0FitRange  = lowMtM0;
      double highMtM0FitRange = highMtM0;
      //If Kaon, only fit the last 6 bins for extrapolation
      if(m_currentPartIndex == 1 && highBin - lowBin > 6){
        lowMtM0FitRange = m_currentHistoToFit->GetBinLowEdge(highBin - 6);
      }



      string fitOptions = "I";
      #ifdef _ZFITTER_DEBUG_
        fitOptions += " V";
      #else
        fitOptions += " Q";
      #endif
      int  status    = m_currentHistoToFit->Fit(electronExtrapFunctName.c_str(),fitOptions.c_str(),"",lowMtM0FitRange,highMtM0FitRange);
      TF1* fitResult = m_currentHistoToFit->GetFunction(electronExtrapFunctName.c_str());

      //Push the results to the histograms
      if(fitResult){
        if(electronRatioIndex == 1){
          double y_val = fitResult->Eval(lowMtM0);
          double slope = fitResult->GetParameter(0)*fitResult->GetParameter(1)*exp(fitResult->GetParameter(1)*lowMtM0);
          elecPionRatioLowMtM0ExtrapLine = PhysMath::linearFunct("elecPionRatioLowMtM0ExtrapLine", fitResult->Eval(lowMtM0), slope, lowMtM0);
          #ifdef _ZFITTER_DEBUG_
            cout << " Electron line low mTm0 to " << y_val << " + " << slope << "*(mT-m0)" << endl;
          #endif
          elecPionRatioLowMtM0ExtrapLine->SetRange(0.0,lowMtM0);
        }
        electronRatioFitHistos[electronRatioIndex] = (TH1D*) m_currentHistoToFit->Clone();
        electronRatioFitHistos[electronRatioIndex]->SetName(Form("elecExtrap_cent%02d_rap%02d_%02d",centIndex,m_currentRapBin,electronRatioIndex));
        int startBinToModify = 1;
        int stopBinToModify = m_currentHistoToFit->GetNbinsX();
        if(m_currentPartIndex == 1){ //If a kaon, only push back from end of fits
          startBinToModify = highBin;
        }
        for(int binX = 1; binX <= stopBinToModify; binX++){
          double content      = m_currentHistoToFit->GetBinContent(binX);
          if(!(content <= 0.0 || binX >= startBinToModify)) continue;
          double center       = m_currentHistoToFit->GetBinCenter(binX);
          double lowBinEdge   = m_currentHistoToFit->GetBinLowEdge(binX);
          double highBinEdge  = m_currentHistoToFit->GetBinLowEdge(binX+1);
          double valueToPush  = fitResult->Integral(lowBinEdge,highBinEdge);
          double valueError   = fitResult->IntegralError(lowBinEdge,highBinEdge);  // cant be done with multithreading
          if(electronRatioIndex != 0){
            valueToPush = valueToPush / m_currentHistoToFit->GetBinWidth(binX);
            valueError  = valueError  / m_currentHistoToFit->GetBinWidth(binX);
          }
          if(valueError <= 0.0) valueError = 1E-13;
           electronRatioFitHistos[electronRatioIndex]->SetBinContent(binX, valueToPush);
           electronRatioFitHistos[electronRatioIndex]->SetBinError(binX, valueError);
          //if(electronRatioIndex == 2){ // only push kaon ratios
           /* if(m_inAllCentMode){
              cout << " Electron Bin " << binX << " Step 4 Reading " << m_Fit_Data_ZTPC[m_currentPartIndex][3][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, binX) << endl;
              cout << " Electron Bin " << binX << " Step 4 Writing " << valueToPush << endl;
              m_Fit_Data_ZTPC[m_currentPartIndex][3][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, binX, valueToPush);
              m_Fit_Data_ZTPC[m_currentPartIndex][3][0][m_currentPlusMinusIndex]->SetBinError(  m_currentRapBin, binX, valueError);
            }else{
              cout << " The previous electron amp val was " << m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][m_currentCentIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, binX) << endl;
              m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][m_currentCentIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, binX,valueToPush);
              m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][m_currentCentIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin, binX,   0.0);
            }*/
          //}
          //cout << "  Electron dN/dmT fit is pushing back " << valueToPush << " +/- " << valueError << "   with TotalHistoIntegral: " << m_rebinZTPCHistoIntegralStorage[binX] << "   prev val: " <<  endl;
        }//end bin loop
        m_electronRatioFitHistos_ByCent[centIndex][electronRatioIndex] = electronRatioFitHistos[electronRatioIndex];
      }else{
        #ifdef _ZFITTER_DEBUG_
          cout << " ERROR: EXTRAPOLATION FAILED!!!!!!!!! ZFitter::calculateElectronRatios() " << endl;
          HistogramUtilities::printHistoInfo(m_currentHistoToFit);
        #endif
        return false;
      }
      HistogramUtilities::discardNanHistoData(m_currentHistoToFit);


      // DRAW THE FITS
      TCanvas* fittingCanvas = new TCanvas("fittingCanvasOut","Fitting Canvas");
      fittingCanvas->SetWindowSize(1300,700);
      fittingCanvas->SetCanvasSize(1200,600);


      fittingCanvas->cd();
      gStyle->SetOptFit(0111);
      //electrondNdMtHisto->GetXaxis()->SetRangeUser(lowMtM0,highMtM0);
      HistogramUtilities::setAxisRanges(m_currentHistoToFit, lowMtM0,highMtM0,true);
      m_currentHistoToFit->SetStats(true);
      m_currentHistoToFit->Draw("E");
      //gPad->SetLogy();
      gPad->SetRightMargin(0.4);


      fittingCanvas->SaveAs(Form("%s/%s/%s/TPCParamExtrapolations/%s_in_%s_space_Cent%02d_Rap_%02d.png",m_imagePreDir.c_str(),
        m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), electronRatioIndex == 1 ? "electronPionRatio" : "electronKaonRatio",
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex,m_currentRapBin));

      #ifdef _ZFITTER_DEBUG_
        cout << " ampExtractFunct " << ampExtractFunct << endl;
        cout << " fittingCanvas " << fittingCanvas << endl;
        cout << " fitResult " << fitResult     << endl;
      #endif

      #ifndef _SAVE_POINTERS_
        delete ampExtractFunct;
        delete fitResult;
        delete fittingCanvas;
      #endif
    }// end ratio index loop
  }// end centrality loop

  return true;
}


//============================================================================================================================================================================


bool ZFitter::calculateElectronMeanExtrapolation(){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::calculateElectronMeanExtrapolation()" << endl;
  #endif




  TH1D* electronMeanHisto = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->ProjectionY(Form("param_proj_%s_%s_%s_Rap%02d",
    "mu", m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(3,m_currentPartCharge).Data(),m_currentRapBin),
    m_currentRapBin,m_currentRapBin);
  string functionName = Form("TPC_param_parab_%s_%s_%s_Rap%02d","mu",
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(3,m_currentPartCharge).Data(),m_currentRapBin);
  electronMeanHisto->GetYaxis()->SetTitle("Gaussian #mu - Prediction");

  //Readjust the mean by the electron bichsel prediction (fit difference between fit result and prediction)
  for(int mTm0Bin = 1; mTm0Bin <= electronMeanHisto->GetNbinsX(); mTm0Bin++){
    if(electronMeanHisto->GetBinError(mTm0Bin) > 0.0){
      double mTm0 = electronMeanHisto->GetBinCenter(mTm0Bin);
      double electronZTPC  = electronMeanHisto->GetBinContent(mTm0Bin);
      double momentumOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity,mTm0,m_partInfo->GetParticleMass(m_currentPartIndex));
      double prediction    = log(m_bichselCurves[3]->Eval(momentumOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(momentumOfBin));
      electronMeanHisto->SetBinContent(mTm0Bin,electronZTPC - prediction );
    }
  }
  int firstBin    = HistogramUtilities::firstBinWithContent(electronMeanHisto);
  int lastBin     = HistogramUtilities::lastBinWithContent(electronMeanHisto);
  int binDiff     = lastBin - firstBin + 1;
  double lowMtM0  = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->GetYaxis()->GetBinLowEdge(firstBin);
  double highMtM0 = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->GetYaxis()->GetBinLowEdge(lastBin);
  double midMtM0  = 0.5*(lowMtM0 + highMtM0);
  if(lowMtM0 <= 0.0){
    cout << " WARNING: the rapidity is too extreme for the electron extrapolation procedure you have given setting low mTm0 to 0" << endl;
    lowMtM0 = 0.0;
  }
  if(highMtM0 <= 0.0){
    cout << " ERROR: the rapidity is too extreme for the electron mean extrapolation procedure, becuase high mTm0 is negative... (" << highMtM0 << " ) so returning without extrapolation" << endl;
    return false;
  }

  //#################   Determine if there is a low pT region fit  ################
  bool electronLowPtWasFit = false;
  bool firstContentFound   = false;
  bool firstHoleFound      = false;

  double  mTm0OfFirstContent     = -1;
  int     binOfFirstContentBegin = -1;
  int     binOfFirstContentEnd   = -1;
  double  mTm0OfSecondContent    = -1;
  int     binOfSecondContentBegin = -1;
  int     binOfSecondContentEnd   = -1;


  double weights_lowPt  = 0.0;
  double average_lowPt  = 0.0;
  double weights_highPt = 0.0;
  double average_highPt = 0.0;

  int numLowPtPoints     = 0;
  int numHighPtPoints    = 0;
  for(int bin = 1; bin <= electronMeanHisto->GetNbinsX(); bin++){
    double content = electronMeanHisto->GetBinContent(bin);
    double error   = electronMeanHisto->GetBinError(bin);
    if(firstContentFound){
      if(firstHoleFound){
        if(error > 0.0 && std::isfinite(error)){
          electronLowPtWasFit = true;
          numHighPtPoints++;
          double weight  = 1.0/(electronMeanHisto->GetBinError(bin)*electronMeanHisto->GetBinError(bin));
          double content = electronMeanHisto->GetBinContent(bin);
          weights_highPt += weight;
          average_highPt += content*weight;
          if(mTm0OfSecondContent < 0.0){
            binOfSecondContentBegin = bin;
            mTm0OfSecondContent = electronMeanHisto->GetBinCenter(bin);
          }
          binOfSecondContentEnd = bin;
        }
      }else{
        if(error <= 0.0 && std::isfinite(error)){
          firstHoleFound = true;
          binOfFirstContentEnd = bin - 1;
        }else{
          numLowPtPoints++;
          double weight  = 1.0/(electronMeanHisto->GetBinError(bin)*electronMeanHisto->GetBinError(bin));
          double content = electronMeanHisto->GetBinContent(bin);
          weights_lowPt += weight;
          average_lowPt += content*weight;
        }
      }
    }else{
      if(error > 0.0 && std::isfinite(error)){
        firstContentFound      = true;
        binOfFirstContentBegin = bin;
        mTm0OfFirstContent = electronMeanHisto->GetBinCenter(bin);
        numLowPtPoints++;
        double weight  = 1.0/(electronMeanHisto->GetBinError(bin)*electronMeanHisto->GetBinError(bin));
        double content = electronMeanHisto->GetBinContent(bin);
        weights_lowPt += weight;
        average_lowPt += content*weight;
      }
    }
  }

  #ifdef _ZFITTER_DEBUG_
    cout << "electronLowPtWasFit : " << electronLowPtWasFit << endl;
    cout << "firstContentFound : " << firstContentFound << endl;
    cout << "firstHoleFound : " << firstHoleFound << endl;
    cout << "mTm0OfFirstContent : " << mTm0OfFirstContent << endl;
    cout << "binOfFirstContentBegin : " << binOfFirstContentBegin << endl;
    cout << "binOfFirstContentEnd : " << binOfFirstContentEnd << endl;
    cout << "mTm0OfSecondContent : " << mTm0OfSecondContent << endl;
    cout << "binOfSecondContentBegin : " << binOfSecondContentBegin << endl;
    cout << "binOfSecondContentEnd : " << binOfSecondContentEnd << endl;
    cout << "weights_lowPt : " << weights_lowPt << endl;
    cout << "average_lowPt : " << average_lowPt << endl;
    cout << "weights_highPt : " << weights_highPt << endl;
    cout << "average_highPt : " << average_highPt << endl;
    cout << "numLowPtPoints : " << numLowPtPoints << endl;
    cout << "numHighPtPoints : " << numHighPtPoints << endl;
  #endif


  if(!electronLowPtWasFit){
    numHighPtPoints = numLowPtPoints;
    numLowPtPoints  = 0;
    weights_highPt  = weights_lowPt;
    average_highPt  = average_lowPt;
    binOfSecondContentBegin = binOfFirstContentBegin;
    binOfSecondContentEnd   = binOfFirstContentEnd;
    binOfFirstContentBegin  = 0;
    binOfFirstContentEnd    = 0;
    weights_lowPt = 0.0;
    average_lowPt = 0.0;
    mTm0OfFirstContent = 0.0;
    #ifdef _ZFITTER_DEBUG_
      cout << "Swapping Electron Regions because Low Pt wasn't fit..." << endl;
    #endif
  }

  average_lowPt  = average_lowPt/weights_lowPt;
  average_highPt = average_highPt/weights_highPt;


  lowMtM0  = electronMeanHisto->GetBinCenter(binOfSecondContentBegin);
  highMtM0 = electronMeanHisto->GetBinCenter(binOfSecondContentEnd);
  midMtM0  = 0.5*(lowMtM0 + highMtM0);

  int lowBin  = binOfSecondContentBegin;
  int midBin  = electronMeanHisto->FindBin(midMtM0);
  int highBin = binOfSecondContentEnd;

  //if(binDiff > 5){
  //  lowMtM0  = electronMeanHisto->GetBinCenter(lastBin - 5);
  //  midMtM0  = electronMeanHisto->GetBinCenter(lastBin - 3);
  //  highMtM0 = electronMeanHisto->GetBinCenter(lastBin);
  //}

  //int lowBin  = electronMeanHisto->FindBin(lowMtM0);
  //int midBin  = electronMeanHisto->FindBin(midMtM0);
  //int highBin = electronMeanHisto->FindBin(highMtM0);




  #ifdef _ZFITTER_DEBUG_
    cout << "Low Bin: " << lowBin << "  MtM0: " << lowMtM0 << endl;
    cout << "Mid Bin: " << midBin << "  MtM0: " << midMtM0 << endl;
    cout << "High Bin: " << highBin << "  MtM0: " << highMtM0 << endl;
    cout << "Extrapolating electron" << " with mTm0 fit range from " << lowMtM0 << "  to  " << highMtM0 << endl;
  #endif

  /*bool properLowBinFound = true;
  if(electronMeanHisto[0]->GetBinError(lowBin) == 0.0){
    properLowBinFound = false;
    for(int binX = lowBin; binX < highBin; binX++){
      if(electronMeanHisto[0]->GetBinError(binX) > 0.0){
        properLowBinFound = true;
        lowBin = binX;
        break;
      }
    }
  }
  bool properHighBinFound = true;
  if(electronMeanHisto[0]->GetBinError(highBin) == 0.0){
    properHighBinFound = false;
    for(int binX = highBin; binX > lowBin; binX--){
      if(electronMeanHisto[0]->GetBinError(binX) > 0.0){
        properHighBinFound = true;
        highBin = binX;
        break;
      }
    }
  }

  if(!properLowBinFound || !properHighBinFound){
    cout << " ERROR: The electron Extrapolation will not work. Could not find bins with actual fit data..." << endl;
  }*/
  if(midBin <= lowBin || midBin >=highBin){
    midBin = (lowBin+highBin)/2;
    cout << " WARNING: Mid-Bin Didn't work, so putting it to half the range" << endl;
  }
  if(highBin - lowBin < 2){
    cout << " ERROR: The electron extrapolation will fail because you only have two bins to extrapolate..." << endl;
    return false;
  }
  lowMtM0  = electronMeanHisto->GetBinCenter(lowBin)  - 0.25*electronMeanHisto->GetBinWidth(lowBin);
  highMtM0 = electronMeanHisto->GetBinCenter(highBin) + 0.25*electronMeanHisto->GetBinWidth(highBin);
  midMtM0  = 0.5*(lowMtM0 + highMtM0);
  midBin   = electronMeanHisto->FindBin(midMtM0);
  midMtM0  = electronMeanHisto->GetBinCenter(midBin);


  #ifdef _ZFITTER_DEBUG_
    cout << " Resized the fit range... lowMtM0: " << lowMtM0 << " midMtM0: " << midMtM0 << " highMtM0: " << highMtM0 << endl;
    cout <<  "                   LowBin: " << lowBin << " midBin: " << midBin << " highBin: " << highBin << endl;
  #endif


  double amplitudes[3];
  amplitudes[0] = electronMeanHisto->GetBinContent(lowBin);
  amplitudes[1] = electronMeanHisto->GetBinContent(midBin);
  amplitudes[2] = electronMeanHisto->GetBinContent(highBin);
  if(amplitudes[0] == 0.0 && amplitudes[1] == 0.0 && amplitudes[2] == 0.0){
    cout << " ERROR: calculateElectronMeanExtrapolation all the amplitudes at the three points are zero... Returning..." << endl;
    return false;
  }

  double slope_1 = (amplitudes[1] - amplitudes[0])/(midMtM0 - lowMtM0);
  double slope_2 = (amplitudes[2] - amplitudes[1])/(highMtM0 - midMtM0);
  double ySlope_0 = 0.5 * (slope_1 + slope_2);
  double maxVal = HistogramUtilities::getMaxInRange(electronMeanHisto, lowMtM0, highMtM0);
  double minVal = HistogramUtilities::getMinInRange(electronMeanHisto, lowMtM0, highMtM0);
  double distanceScale = maxVal - minVal;
  double slopeScale = distanceScale / (highMtM0 - lowMtM0);

  #ifdef _ZFITTER_DEBUG_
    cout << " Initializing linear funciton with y0=" << amplitudes[1] << " y'= " << ySlope_0 << endl;
  #endif

  TF1* extrapFunction = PhysMath::linearFunct(functionName.c_str(), amplitudes[1], ySlope_0,midMtM0);
  extrapFunction->SetRange(lowMtM0,highMtM0);
  extrapFunction->SetParLimits(0,amplitudes[1] - 10.0*distanceScale,amplitudes[1] + 10.0*distanceScale);
  extrapFunction->SetParLimits(1,-10.0*slopeScale, 10.0*slopeScale);

  TF1* initialFunct = (TF1*) extrapFunction->Clone();
  initialFunct->SetLineColor(kYellow);
  initialFunct->SetName("initialFunct");

  int status     = electronMeanHisto->Fit(functionName.c_str(),"","",lowMtM0,highMtM0);
  TF1* fitResult = electronMeanHisto->GetFunction(functionName.c_str());

  //Push the results to the histograms
  if(fitResult){
    for(int binX = binOfFirstContentEnd + 1; binX <= electronMeanHisto->GetNbinsX(); binX++){
      double mTm0          = electronMeanHisto->GetBinCenter(binX);
      double momentumOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity,mTm0,m_partInfo->GetParticleMass(m_currentPartIndex));
      double prediction    = log(m_bichselCurves[3]->Eval(momentumOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(momentumOfBin));
      double valueToPush = fitResult->Eval(mTm0)  + prediction;
      m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin, binX,valueToPush);
      m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin, binX, 1.0E-13);
    }
  }else{
    #ifdef _ZFITTER_DEBUG_
      cout << "\n ERROR: ELECTRON MEAN PARAMETER EXTRAPOLATION FAILED!!!!!!!!! \n " << endl;
      HistogramUtilities::printHistoInfo(electronMeanHisto);
    #endif
    return false;
  }

  HistogramUtilities::discardNanHistoData(electronMeanHisto);
  TGraph* parabPoint = NULL;
  if(fitResult){
    parabPoint = new TGraph();
    parabPoint->SetPoint(0,midMtM0,fitResult->GetParameter(0));
    parabPoint->SetMarkerStyle(20);
    parabPoint->SetMarkerColor(kRed);
  }


  // DRAW THE FITS
  TCanvas* fittingCanvas = new TCanvas(Form("fittingCanvas%02d",1),"Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  HistogramUtilities::setAxisRanges(electronMeanHisto, 0,0.5,false);
  electronMeanHisto->SetStats(true);
  electronMeanHisto->Draw("E");

  initialFunct->Draw("same");
  if(parabPoint) parabPoint->Draw("Psame");
  gPad->SetRightMargin(0.4);
  fittingCanvas->SaveAs(Form("%s/%s/%s/TPCParamExtrapolations/Param%02d_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
          m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),1,
          m_partInfo->GetParticleName(3,m_currentPartCharge).Data(),
          m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));


  #ifdef _ZFITTER_DEBUG_
    cout << " electronMeanHisto[jjj] " << electronMeanHisto << endl;
    cout << " initialFunct " << initialFunct << endl;
    cout << " extrapFunctions[jjj] " << extrapFunction << endl;
    cout << " parabPoint " << parabPoint << endl;
  #endif

  #ifndef _SAVE_POINTERS_
    delete electronMeanHisto;
    if(initialFunct) delete initialFunct;
    delete extrapFunction;
    if(parabPoint) delete parabPoint;
    delete fittingCanvas;
  #endif

  return true;
}


//============================================================================================================================================================================


void ZFitter::calculateMeanShifts(double &a_delta_Pi,double &a_delta_K,double &a_delta_E){

  int partIndexIndexArray[3] = {0,1,3};

  for(int partIndexIndex = 0; partIndexIndex < 3; partIndexIndex++){
    int partIndex = partIndexIndexArray[partIndexIndex];
    TH2D* histo = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex];
    double weights = 0.0;
    double mean    = 0.0;

    for(int binY = 1; binY <= histo->GetNbinsY(); binY++){
      double mTm0      = histo->GetYaxis()->GetBinCenter(binY);
      double momentum  = PhysMath::pTot_from_y_mTm0(m_currentRapidity,mTm0,m_partInfo->GetParticleMass(m_currentPartIndex));
      double binCount  = histo->GetBinContent(m_currentRapBin,binY) - log(m_bichselCurves[partIndex]->Eval(momentum)/m_bichselCurves[m_currentPartIndex]->Eval(momentum));
      double binWeight = 1.0/(histo->GetBinError(m_currentRapBin,binY)*histo->GetBinError(m_currentRapBin,binY));
      if(!std::isfinite(binWeight) || !std::isfinite(binCount)) continue;
      mean    += binWeight*binCount;
      weights += binWeight;
    }
    mean = mean/weights;

    if(partIndex == 0) a_delta_Pi = mean;
    if(partIndex == 1) a_delta_K  = mean;
    if(partIndex == 3) a_delta_E  = mean;
  }

  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::calculateMeanShifts()" << endl;
    cout << "       Delta Mu Pi : " << a_delta_Pi << endl;
    cout << "       Delta Mu K  : " << a_delta_K << endl;
    cout << "       Delta Mu E  : " << a_delta_E << endl;
  #endif

  return;
}

//============================================================================================================================================================================


bool ZFitter::calculatePionGammaLowPtExtrapolation(){
  TH2D* histo = m_Fit_Data_ZTPC[m_currentPartIndex][0][3][m_currentPlusMinusIndex]; ;
  int firstBin = HistogramUtilities::firstYBinWithContent(histo, m_currentRapBin);
  //int lastBin  = HistogramUtilities::lastYBinWithContent(histo, m_currentRapBin);
  if(firstBin < 1 || firstBin >= histo->GetNbinsY()){
    cout << "ERROR: ZFitter::calculatePionGammaLowPtExtrapolation didn't have any filled bins to extrapolate..." << endl;
    return false;
  }
  double firstGamma = histo->GetBinContent(m_currentRapBin,firstBin);
  for(int binY = 1; binY < firstBin; binY++){
    histo->SetBinContent(m_currentRapBin,binY, firstGamma);
    histo->SetBinError(m_currentRapBin,binY, 1E-6);
  }
  #ifdef _ZFITTER_DEBUG_
    cout << "First good pion gamma found at bin " << firstBin << " with a value of " << firstGamma << "  " <<  endl;
    TH1D* proj = histo->ProjectionY("proj",m_currentRapBin,m_currentRapBin);
    HistogramUtilities::printHistoInfo(proj);
  #endif

  return true;
}



//============================================================================================================================================================================







void ZFitter::fitTPCPionKaon_SimulCent_ByRapidity(){
  #ifdef _ZFITTER_DEBUG_
    cout << " Call to ZFitter::fitTPCPionKaon_SimulCent_ByRapidity()" << endl;
  #endif

  TVirtualFitter::SetMaxIterations(1000000);

  bool   usePositiveChargeAsSeeds = true;
  bool   usePreviousGoodParams    = true;
  double kaonIsolatedStopMomentum    = 0.25;
  double kaonIsolatedLowMomForExtrap = 0.28;
  double sigmaGuess       = log(1.07); //0.067 ish for the ~7% resolution of the tpc
  double electronStartMom = 0.30; // for easy fiting
  double electronStopMom  = 0.40; // for easy fitting
  double electronStartMom_afterStdFix = 0.27; // for easy fiting
  double electronStopMom_afterStdFix  = 0.43; // for easy fitting
  double electronEarlyStopMom = 0.11; //electron region on left of pion peak

  #ifndef _ZFITTER_PHYSMATH_BYPASS_
  TF1* electronStopMomFunct = PhysMath::parabolaFunct("electronStartMomFunct",-0.6,0.45,0,0.40,0.6,0.45);
  #endif


  double electronStartFlatStdDevExtrap = 0.38; // fits from this to the next with a constant
  double electronStopFlatStdDevExtrap  = 0.40; // fits to this with a constant

  double electronStartDnDmtExtrap = 0.30; // fits from this to the next with an exponential  // also start of last fit constraint amplitude
  double electronStopDnDmtExtrap  = 0.4; // fits to this with an exponential
  double electronLowMom_FreeMeanAndWidthBelow = 0.11; // electron region on left of pion peak

  //double pionHighMomParamExtrap = 0.55; // fits to this with lines
  //double kaonHighMomParamExtrap = 0.55; // fits to this with lines

  double momentumToKillElectrons = 0.605;


  double pionSkewSeed = 0.75;

  double maximumMomentumToFit = 0.75;


  //double pionSeedParams[4] = {-0.1155,0.07511,24.85,-14.67};
  //double pionSeedParams[4] = {-0.1155,0.513,24.85,-14.67};
  double pionSeedParams[4] = {-0.1155,0.513,0.4613,-0.2854};




  //##################################    DON"T EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING (ASK MATT FIRST PLEASE) #################################

  TH3D* histo3DPM = NULL;
  if(m_currentPartCharge == 1){ // Positive
    histo3DPM = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex];
  }else{ // Negative
    histo3DPM = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex];
  }
  #ifndef _ZFITTER_PHYSMATH_BYPASS_
    electronStopMom = electronStopMomFunct->Eval(m_currentRapidity);
  #else
    electronStopMom = 0.4;
  #endif
  if(electronStopMom > 0.40) electronStopMom = 0.40;

  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int histogramIndex = 0; histogramIndex < 1000; histogramIndex++){
      #ifndef _SAVE_POINTERS_
        if(m_rebinZTPCHistoStorage_ByCent[centIndex][histogramIndex]) delete m_rebinZTPCHistoStorage_ByCent[centIndex][histogramIndex];
      #endif
      m_rebinZTPCHistoStorage_ByCent[centIndex][histogramIndex] = NULL;
    }
  }

  int fitting_round = 0;
  vector<string> roundTag;
  vector<double> mTm0IndexVector;
  vector<int>    mTm0BinsFitPerRound;
  vector<unsigned int>    indexOfRoundChange;

  //##############   STEP 0  #############
  double step0_lowlowMtM0  = PhysMath::mTm0_from_y_p(m_currentRapidity, electronEarlyStopMom, m_partInfo->GetParticleMass(m_currentPartIndex));
  double step0_lowMtM0     = PhysMath::mTm0_from_y_p(m_currentRapidity, electronStartMom, m_partInfo->GetParticleMass(m_currentPartIndex));
  double step0_highMtM0    = PhysMath::mTm0_from_y_p(m_currentRapidity, electronStopMom, m_partInfo->GetParticleMass(m_currentPartIndex));
  int    step0_lowlowBin   = histo3DPM->GetYaxis()->FindBin(step0_lowlowMtM0);
  int    step0_lowBin      = histo3DPM->GetYaxis()->FindBin(step0_lowMtM0);
  int    step0_highBin     = histo3DPM->GetYaxis()->FindBin(step0_highMtM0);
  if(step0_lowlowMtM0 > histo3DPM->GetYaxis()->GetBinCenter(step0_lowlowBin)) step0_lowlowBin++;
  if(step0_lowMtM0 > histo3DPM->GetYaxis()->GetBinCenter(step0_lowBin)) step0_lowBin++;
  if(step0_highMtM0 < histo3DPM->GetYaxis()->GetBinCenter(step0_highBin)) step0_highBin--;
  if(step0_lowBin <= 0 && step0_highBin <=0){
    step0_lowBin  = -1;
    step0_highBin = -2;
  }else if(step0_lowBin <= 0){
    step0_lowBin = 1;
  }
  for(int aaa = 1; m_currentPartIndex == 0 && aaa <= step0_lowlowBin; aaa++) mTm0IndexVector.push_back(aaa);
  for(int bbb = step0_lowBin; bbb <= step0_highBin; bbb++) mTm0IndexVector.push_back(bbb);
  roundTag.push_back("electronRegion");
  mTm0BinsFitPerRound.push_back(0);
  indexOfRoundChange.push_back(mTm0IndexVector.size());
  cout << "Step 0 " << roundTag[0] << " Bin Ranges [1," << step0_lowlowBin << "]" << "   [" << step0_lowBin << "," << step0_highBin << "]" << endl;


  //##############   STEP 1  #############
  double step1_lowMtM0  = PhysMath::mTm0_from_y_p(m_currentRapidity, electronStartMom_afterStdFix, m_partInfo->GetParticleMass(m_currentPartIndex));
  double step1_highMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity, electronStopMom_afterStdFix, m_partInfo->GetParticleMass(m_currentPartIndex));
  int    step1_lowBin   = histo3DPM->GetYaxis()->FindBin(step1_lowMtM0);
  int    step1_highBin  = histo3DPM->GetYaxis()->FindBin(step1_highMtM0);
  if(step1_lowMtM0 > histo3DPM->GetYaxis()->GetBinCenter(step1_lowBin)) step1_lowBin++;
  if(step0_highMtM0 > histo3DPM->GetYaxis()->GetBinCenter(step1_highBin)) step1_highBin--;
  if(step1_lowBin <= 0 && step1_highBin <=0){
    step1_lowBin  = -1;
    step1_highBin = -2;
  }else if(step1_lowBin <= 0){
    step1_lowBin = 1;
  }
  for(int ccc = 1; m_currentPartIndex == 0 && ccc <= step0_lowlowBin; ccc++) mTm0IndexVector.push_back(ccc);
  for(int ddd = step1_lowBin; ddd <= step1_highBin; ddd++) mTm0IndexVector.push_back(ddd);
  roundTag.push_back("electronRegionStdFix");
  mTm0BinsFitPerRound.push_back(0);
  indexOfRoundChange.push_back(mTm0IndexVector.size());
  cout << "Step 1 " << roundTag[1] << " Bin Ranges [1," << step0_lowlowBin << "]" << "   [" << step1_lowBin << "," << step1_highBin << "]" << endl;


  //##############   STEP 2  #############
  double step2_lowBin   = 1;
  //double step2_highMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity, pionHighMomParamExtrap,  m_partInfo->GetParticleMass(m_currentPartIndex));
  //if(pionHighMomParamExtrap <= kaonHighMomParamExtrap){
  //  step2_highMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity, kaonHighMomParamExtrap,  m_partInfo->GetParticleMass(m_currentPartIndex));
  //}
  //int    step2_highBin  = histo3DPM->GetYaxis()->FindBin(step2_highMtM0);
  int    step2_highBin  = histo3DPM->GetNbinsY();
  for(int fff = step2_lowBin; fff <= step2_highBin; fff++) mTm0IndexVector.push_back(fff);
  roundTag.push_back("electronRatioFixed");
  mTm0BinsFitPerRound.push_back(0);
  indexOfRoundChange.push_back(mTm0IndexVector.size());
  cout << "Step 2 " << roundTag[2] << " Bin Ranges [" << step2_lowBin << "," << step2_highBin << "]" << endl;


  #ifdef _ZFITTER_DEBUG_
    cout << "------------------   MtM0 Index Index Array ------------------" << endl;
    for(unsigned int ggg = 0; ggg < mTm0IndexVector.size(); ggg++) cout << mTm0IndexVector[ggg] << "  ";
    cout << endl << "-----------------------------------------------------------------" << endl;
  #endif


  double pionMeanShift   = 0;
  double kaonMeanShift   = 0;
  double protonMeanShift = 0;

  int numberOfWidthCalculations = 0;

  double previousGoodParameters[15] = {-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999};
  for(int paramIndex = 0; paramIndex < 11; paramIndex++) previousGoodParameters[paramIndex] = -999;

  bool electronWidthSuccess  = false;
  bool electronRatioSuccess  = false;
  bool electronMeanSuccess   = false;
  bool pionGammaSuccess      = false;


  for(unsigned int VertBinIndexIndex = 0; VertBinIndexIndex < mTm0IndexVector.size(); VertBinIndexIndex++){
    // What to do between trasitions of rounds of fitting
    if(VertBinIndexIndex == 0){
      cout << "########  STARTING ROUND 0: " << roundTag[fitting_round] << endl;
    }

    m_currentMtM0Bin  = mTm0IndexVector[VertBinIndexIndex];
    m_currentMtM0     = histo3DPM->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
    m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
    if(!std::isfinite(m_currentMomOfBin) || m_currentMomOfBin <= 0.0) continue;
    if(m_currentMomOfBin > maximumMomentumToFit) continue;

    //cout << "BIN CHECK: " << m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) << endl;

    if(VertBinIndexIndex > 0 && mTm0IndexVector[VertBinIndexIndex] <= mTm0IndexVector[VertBinIndexIndex - 1]){
      #ifdef _ZFITTER_DEBUG_
        cout << "Fitting Round Was " << fitting_round << " Is Now " << fitting_round + 1 << endl;
      #endif

      if(roundTag[fitting_round] == "electronRegion"){
        electronWidthSuccess = calculateAverageElectronWidth();
        if(!electronWidthSuccess) cout << "WARNING: Electron Width Calculation Failed..." << endl;
        else                     cout << "Electron Width Calculation Successful..." << endl;
        numberOfWidthCalculations++;
        calculateMeanShifts(pionMeanShift,kaonMeanShift,protonMeanShift);
      }
      if(roundTag[fitting_round] == "electronRegionStdFix"){
        if(!electronWidthSuccess && numberOfWidthCalculations < 2){
          electronWidthSuccess = calculateAverageElectronWidth();
          VertBinIndexIndex = indexOfRoundChange[0] - 1;
          fitting_round--;
          numberOfWidthCalculations++;
          continue; //reset to round
        }
        electronRatioSuccess = calculateElectronRatios();
        if(!electronRatioSuccess) cout << "WARNING: Electron Ratio Calculation Failed..." << endl;
        else                      cout << "Electron Ratio Calculation Successful..." << endl;
        calculateMeanShifts(pionMeanShift,kaonMeanShift,protonMeanShift);
        electronMeanSuccess = calculateElectronMeanExtrapolation();
        if(!electronMeanSuccess) cout << "WARNING: Electron Mean Calculation Failed..." << endl;
        else                     cout << "Electron Mean Calculation Successful..." << endl;

        pionGammaSuccess = calculatePionGammaLowPtExtrapolation();
        if(!electronMeanSuccess) cout << "WARNING: Pion Gamma Calculation Failed..." << endl;
        else                     cout << "Pion Gamma Calculation Successful..." << endl;


      }
      //clear the previous good paramters because the bins aren't consecutive at round change
      for(int paramIndex = 0; paramIndex < 11; paramIndex++) previousGoodParameters[paramIndex] = -999;

      fitting_round++;
      cout << "########  STARTING ROUND " << fitting_round << ": " << roundTag[fitting_round] << endl;

    }

    bool isConsecutiveBin = false;
    if(VertBinIndexIndex > 0 && mTm0IndexVector[VertBinIndexIndex] == mTm0IndexVector[VertBinIndexIndex - 1] + 1) isConsecutiveBin = true;




    cout << "\n\n\n###########################################################################" << endl;
    cout << "ZTPC Round " << fitting_round << "   " << roundTag[fitting_round] << " : " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data()
         << " RapBin: " << m_currentRapBin << " MtM0Bin: " << m_currentMtM0Bin << "  p= " << m_currentMomOfBin <<  endl;

    //################ GENERATE THE HISTOGRAMS TO FIT ################
    bool currentBinWasFitBefore = false;
    if(m_rebinZTPCHistoStorage_ByCent[0][m_currentMtM0Bin]) currentBinWasFitBefore = true;
      #ifdef _ZFITTER_DEBUG_
        cout << "  Fit Before? " << (currentBinWasFitBefore ? "yes" : "no") << "  with ptr: " << m_rebinZTPCHistoStorage_ByCent[0][m_currentMtM0Bin] <<  endl;
      #endif
    if(!currentBinWasFitBefore){
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        m_currentHistosToFit_ByCent[centIndex] = grabSpecificHistogram(m_currentPartIndex, m_currentPartCharge, centIndex, 0, m_currentRapBin, m_currentMtM0Bin, 80, false, -999, -999, true,false);
        m_rebinZTPCHistoStorage_ByCent[centIndex][m_currentMtM0Bin] = m_currentHistosToFit_ByCent[centIndex];
        m_currentHistosToFit_ByCent[centIndex]->SetTitle(Form("%s Cent %d y=[%1.2f,%1.2f] m_{T}-m_{0}=[%1.4f,%1.4f] p=%1.4f Simultaneous Cent Fit ; Z_{TPC %s}; #frac{1}{N_{events}} #frac{dN}{dZ_{TPC}}",
                    m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data(),
                    centIndex,
                    histo3DPM->GetXaxis()->GetBinLowEdge(m_currentRapBin),
                    histo3DPM->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
                    histo3DPM->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
                    histo3DPM->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
                    m_currentMomOfBin,
                    m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data()));
        HistogramUtilities::removeHighCountCorruptData(m_currentHistosToFit_ByCent[centIndex]);
        if(m_mTm0_ZTPC_by_Rapidity[m_currentPartIndex][centIndex][m_currentRapBin-1][m_currentPlusMinusIndex]){
          m_mTm0_ZTPC_by_Rapidity[m_currentPartIndex][centIndex][m_currentRapBin-1][m_currentPlusMinusIndex]->SetYProjection(m_currentMtM0Bin-1,m_currentHistosToFit_ByCent[centIndex],true);
        }else{
          cout << "WARNING: mTm0_ZTPC_by_Rapidity histogram not found for " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " Cent " << centIndex << " RapBin " << m_currentRapBin << " PlusMinus " << m_currentPlusMinusIndex << endl;
        }
      }
    }else{
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        m_currentHistosToFit_ByCent[centIndex] = m_rebinZTPCHistoStorage_ByCent[centIndex][m_currentMtM0Bin];
      }
    }
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      #ifdef _ZFITTER_DEBUG_
        cout << " Current Histo To Fit Cent " << centIndex << " Ptr: " << m_currentHistosToFit_ByCent[centIndex] << endl;
      #endif
      m_currentLowFitBin[centIndex]  = 1;
      m_currentHighFitBin[centIndex] = m_currentHistosToFit_ByCent[centIndex]->GetNbinsX();
    }


    double countsOfHisto = HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[0],-1,true) * m_centEvents[m_currentPartIndex]->GetBinContent(1);
    #ifdef _ZFITTER_DEBUG_
      cout <<  "mTm0BinsFitPerRound[0] = " << mTm0BinsFitPerRound[0] << endl;
      cout <<  "mTm0BinsFitPerRound[1] = " << mTm0BinsFitPerRound[1] << endl;
      cout <<  "mTm0BinsFitPerRound[2] = " << mTm0BinsFitPerRound[2] << endl;
      cout << "electronWidthSuccess: " << (electronWidthSuccess ? "TRUE" : "FALSE") << endl;
      cout << "electronRatioSuccess: " << (electronRatioSuccess ? "TRUE" : "FALSE") << endl;
      cout << "electronMeanSuccess: " << (electronMeanSuccess ? "TRUE" : "FALSE") << endl;
      cout << " Cent 0 has " << HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[0],-1,true)*m_centEvents[m_currentPartIndex]->GetBinContent(1) << " counts " << endl;
    #endif
    if(countsOfHisto < 600){
      cout << "   Too few counts in Cent 0 histogram (" << countsOfHisto << ") ... Continuing... " << endl;
      continue;
    }


    //################################    SETUP PION KAON ELECTRON PARAMETERS FOR FIT  ##################################

    double pionBichselPredict      =  log(m_bichselCurves[0]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + pionMeanShift;
    double kaonBichselPredict      =  log(m_bichselCurves[1]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + kaonMeanShift;
    double electronBichselPredict  =  log(m_bichselCurves[3]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + protonMeanShift;
    double pionMeanInit            =  pionBichselPredict;
    double kaonMeanInit            =  kaonBichselPredict;
    double electronMeanInit        =  electronBichselPredict;
    double pionSigmaInit           =  sigmaGuess;
    double pionSkewInit            =  pionSkewSeed;
    double kaonSigmaInit           =  sigmaGuess;
    double electronSigmaInit       =  sigmaGuess;

    //If the data has been fit before, use the previous parameters
    if(currentBinWasFitBefore){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Parameters from previous fits... " << endl;
      #endif
      pionMeanInit        = m_Fit_Data_ZTPC[m_currentPartIndex][0][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      pionSigmaInit       = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      pionSkewInit        = m_Fit_Data_ZTPC[m_currentPartIndex][0][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      kaonMeanInit        = m_Fit_Data_ZTPC[m_currentPartIndex][1][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      kaonSigmaInit       = m_Fit_Data_ZTPC[m_currentPartIndex][1][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      electronMeanInit    = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      electronSigmaInit   = m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
    } else if(usePositiveChargeAsSeeds && m_currentPartCharge < 0 ){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Parameters from Positive Fits... " << endl;
      #endif
      //double pionMeanTemp        = m_Fit_Data_ZTPC[m_currentPartIndex][0][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      //double pionSigmaTemp       = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double kaonMeanTemp          = m_Fit_Data_ZTPC[m_currentPartIndex][1][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double kaonSigmaTemp         = m_Fit_Data_ZTPC[m_currentPartIndex][1][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double electronMeanTemp      = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double electronSigmaTemp     = m_Fit_Data_ZTPC[m_currentPartIndex][3][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);

      //if(pionSigmaTemp > 0.0)   pionMeanInit = pionMeanTemp;
      //if(pionSigmaTemp > 0.0)   pionSigmaInit = pionSigmaTemp;
      if(kaonSigmaTemp > 0.0)     kaonMeanInit      = kaonMeanTemp;
      if(kaonSigmaTemp > 0.0)     kaonSigmaInit     = kaonSigmaTemp;
      if(electronSigmaTemp > 0.0) electronMeanInit  = electronMeanTemp;
      if(electronSigmaTemp > 0.0) electronSigmaInit = electronSigmaTemp;

      #ifdef _ZFITTER_DEBUG_
        cout << "       Read Kaon Params...       mu: " << kaonMeanInit << "  sig: " << kaonSigmaInit << endl;
        cout << "       Read Electron Params...   mu: " << electronMeanInit << "  sig: " << electronSigmaInit << endl;
      #endif
    }


    /*else if(roundTag[fitting_round] != "electronRegion" && roundTag[fitting_round] != "electronRegionStdFix"){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Parameters from gaussian calculation... " << endl;
      #endif
      double pionGausData[3];
      HistogramUtilities::getGausData(pionGausData, m_currentHistosToFit_ByCent[0],true,false, pionMeanInit - 2.0*sigmaGuess, pionMeanInit + 2.0*sigmaGuess);
      if(fabs(pionGausData[1] - pionMeanInit) < 1.2*sigmaGuess) pionMeanInit = pionGausData[1];
      double kaonGausData[3];
      HistogramUtilities::getGausData(kaonGausData, m_currentHistosToFit_ByCent[0],true,false, kaonMeanInit - 2.0*sigmaGuess, kaonMeanInit + 2.0*sigmaGuess);
      if(fabs(kaonGausData[1] - kaonMeanInit) < 1.2*sigmaGuess) kaonMeanInit = kaonGausData[1];
      double electronGausData[3];
      HistogramUtilities::getGausData(electronGausData, m_currentHistosToFit_ByCent[0],true,false, electronMeanInit - 2.0*sigmaGuess, electronMeanInit + 2.0*sigmaGuess);
      if(fabs(electronGausData[1] - electronMeanInit) < 1.2*sigmaGuess) electronMeanInit = electronGausData[1];
    }*/




    //######################    LOAD PREVIOUS GOOD PARAMETERS AS SEEDS   ###########################
    if(usePreviousGoodParams && !currentBinWasFitBefore && isConsecutiveBin){ // use the previous fit parameters as seeds
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Parameters from Previous Good Parameters... " << endl;
      #endif
      // Dont do in specific cent, the all cent would be better seeds
      if(previousGoodParameters[1] > 0 && std::isfinite(previousGoodParameters[0]) && std::isfinite(previousGoodParameters[1]) && std::isfinite(previousGoodParameters[2])){
        if(fabs(previousGoodParameters[0]) < 0.45) pionMeanInit        = previousGoodParameters[0] + log(m_bichselCurves[0]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        if(previousGoodParameters[1] < 0.2)        pionSigmaInit       = previousGoodParameters[1];
        if(fabs(previousGoodParameters[2]) > 0.95) pionSkewInit = 0.85;
        else                                       pionSkewInit = previousGoodParameters[2];
      }
      if(previousGoodParameters[4] > 0 && std::isfinite(previousGoodParameters[3]) && std::isfinite(previousGoodParameters[4])){
        if(fabs(previousGoodParameters[3]) < 0.45) kaonMeanInit        = previousGoodParameters[3] + log(m_bichselCurves[1]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        if(previousGoodParameters[4] < 0.2)        kaonSigmaInit       = previousGoodParameters[4];
      }
      if(previousGoodParameters[6] > 0 && std::isfinite(previousGoodParameters[5]) && std::isfinite(previousGoodParameters[6])){
        if(fabs(previousGoodParameters[5]) < 0.45) electronMeanInit    = previousGoodParameters[5] + log(m_bichselCurves[3]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        if(previousGoodParameters[6] < 0.2)        electronSigmaInit   = previousGoodParameters[6];
      }
    }
    if(roundTag[fitting_round] != "electronRegion" && roundTag[fitting_round] != "electronRegionStdFix" && electronMeanSuccess){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Electron Mean From Param Storage Histogram..." << endl;
      #endif
      electronMeanInit = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
    }
    if(roundTag[fitting_round] != "electronRegion" && electronWidthSuccess){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Electron StdDev From Param Storage Histogram..." << endl;
      #endif
      electronSigmaInit = m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
    }



    double lowRangeData  = m_currentHistosToFit_ByCent[0]->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistosToFit_ByCent[0]));
    double highRangeData = m_currentHistosToFit_ByCent[0]->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistosToFit_ByCent[0])+1);

    //#########################   DETERMINE WHICH PARTICLES AND PARAMETERS SHOULD BE FIT  ###################################
    bool doFitPion            = true; // these are defaults to be changed below under certain conditions
    bool doFitKaon            = true;
    bool doFitElectron        = true;
    bool fixElectronAmplitude = false;
    bool fixKaonAmplitude     = false;
    bool fixPionMean          = false;
    bool fixKaonMean          = false;
    bool fixElectronMean      = false;
    bool fixPionWidths        = false;
    bool fixKaonWidths        = false;
    bool fixElectronWidths    = false;
    bool fixPionGamma         = false;

    //################    DETERMINE FIT RANGES AND PARTICLE TOGGLES  ########################
    double lowRange  = pionMeanInit - 3.5*sigmaGuess;
    double highRange = kaonMeanInit + 3.5*sigmaGuess;
    if(pionBichselPredict > electronBichselPredict){
      lowRange   = lowRangeData;
      //lowRange  = electronBichselPredict - 3.5*sigmaGuess;
    }
    if(roundTag[fitting_round] == "electronRegion"){
      if(m_currentPartIndex == 0){
        if(m_currentMomOfBin > electronStartMom){
          lowRange      = pionBichselPredict - 0.4*sigmaGuess; // the pions will be fit in full in the next round
          highRange     = kaonBichselPredict + 6.0*sigmaGuess;
        }else{
          lowRange      = lowRangeData;
          highRange     = pionBichselPredict + 4.0*sigmaGuess;
        }
      }else{
        lowRange      = pionBichselPredict - 0.4*sigmaGuess; // the pions will be fit in full in the next round
        highRange     = kaonBichselPredict + 6.0*sigmaGuess;
      }
    }
    if(roundTag[fitting_round] == "electronRegionStdFix"){
      //lowRange      = kaonBichselPredict - 3.0*sigmaGuess;
      highRange     = kaonBichselPredict + 6.0*sigmaGuess;
    }

    if(roundTag[fitting_round] == "electronRatioFixed"){
      if(m_currentPartIndex == 1 && m_currentMomOfBin < kaonIsolatedStopMomentum){
        doFitElectron = false;
        doFitPion     = false;
        lowRange      = kaonBichselPredict - 3.0*sigmaGuess;
        highRange     = kaonBichselPredict + 4.5*sigmaGuess;
      }

      if(m_currentPartIndex == 0 && m_currentMomOfBin < kaonIsolatedStopMomentum){
        doFitKaon = false;
        //lowRange = kaonBichselPredict - 3.5*sigmaGuess;
        lowRange   = lowRangeData;
        if(pionBichselPredict > electronBichselPredict){
          //lowRange  = electronBichselPredict - 3.5*sigmaGuess;
          lowRange   = lowRangeData;
          highRange = pionBichselPredict     + 5.5*sigmaGuess;
        }else{
          lowRange  = pionBichselPredict     - 3.5*sigmaGuess;
          highRange = electronBichselPredict + 3.5*sigmaGuess;
        }
      }
    }

    if(m_currentMomOfBin > momentumToKillElectrons){
      doFitElectron = false;
      highRange = kaonMeanInit + 1.5*kaonSigmaInit;
    }

    if(lowRange < lowRangeData)   lowRange  = lowRangeData;
    if(highRange > highRangeData) highRange = highRangeData;
    m_currentLowFitRange  = lowRange;
    m_currentHighFitRange = highRange;

    if(m_currentHistosToFit_ByCent[0]->FindBin(m_currentHighFitRange) - m_currentHistosToFit_ByCent[0]->FindBin(m_currentLowFitRange) < 8 ){
      cout << "WARNING: Range is too small Z_TPC = [" << m_currentLowFitRange << " , " << m_currentHighFitRange << "]"
           << "  with bins = [" << m_currentHistosToFit_ByCent[0]->FindBin(m_currentLowFitRange) << " , " << m_currentHistosToFit_ByCent[0]->FindBin(m_currentHighFitRange)
           << "]   Continuing to next bin..." << endl;
      continue;
    }

    #ifdef _ZFITTER_DEBUG_
      cout << " LowFitRange: " << m_currentLowFitRange << "  highFitRange: " << m_currentHighFitRange << endl;
    #endif
    for(int centIndex = 0;centIndex < m_numCentralities; centIndex++){
      int binLow  =  m_currentHistosToFit_ByCent[centIndex]->FindBin(m_currentLowFitRange);
      if(binLow > 1 && binLow < m_currentHistosToFit_ByCent[centIndex]->GetNbinsX()) m_currentLowFitBin[centIndex] = binLow;
      else binLow = 1;
      int binHigh =  m_currentHistosToFit_ByCent[centIndex]->FindBin(m_currentHighFitRange);
      if(binHigh > 1  && binHigh  <  m_currentHistosToFit_ByCent[centIndex]->GetNbinsX()) m_currentHighFitBin[centIndex] = binHigh;
      else binHigh = m_currentHistosToFit_ByCent[centIndex]->GetNbinsX();
      #ifdef _ZFITTER_DEBUG_
        cout << "    Cent: " << centIndex << "  bin " << m_currentLowFitBin[centIndex] << " to " << m_currentHighFitBin[centIndex] << endl;
      #endif
    }


    if(m_currentPartIndex == 0 && m_currentMomOfBin < electronEarlyStopMom) doFitKaon = false;
    if(m_currentPartIndex == 0 && electronBichselPredict < pionBichselPredict) doFitKaon = false;

    if(roundTag[fitting_round] != "electronRegion" && roundTag[fitting_round] != "electronRegionStdFix"){

      if(electronLowMom_FreeMeanAndWidthBelow < m_currentMomOfBin){
        if(electronMeanSuccess)  fixElectronMean   = true;
        if(electronWidthSuccess) fixElectronWidths = true;
      }

      //Start Fixing the e/k or e/pi ratio after 75% of the mTm0 range
      if(m_currentMomOfBin - electronStartDnDmtExtrap > 0.75*(electronStopDnDmtExtrap - electronStartDnDmtExtrap)){
        if(electronRatioSuccess){
          if(m_currentPartIndex == 1){
            if(m_currentMomOfBin > electronStopMom_afterStdFix) fixElectronAmplitude = true;
          }else{
            fixElectronAmplitude = true;
          }
        }
        if(electronMeanSuccess) {
          electronMeanInit     = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          fixElectronMean      = true;
        }
        if(electronWidthSuccess){
          electronSigmaInit    = m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          fixElectronWidths    = true;
        }
      }

      if(m_currentPartIndex == 0 && m_currentMomOfBin < electronStartDnDmtExtrap){
        if(electronRatioSuccess) fixElectronAmplitude  = true;
      }

      if(m_currentPartIndex == 0 && pionGammaSuccess && m_currentMomOfBin < electronStartMom_afterStdFix){
        pionSkewInit = m_Fit_Data_ZTPC[m_currentPartIndex][0][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        fixPionGamma = true;
      }

    }

    if(roundTag[fitting_round] != "electronRegion" && electronWidthSuccess) fixElectronWidths = true;

    #ifdef _ZFITTER_DEBUG_
      cout << "####### This Fit will Have The Following Toggles: " << endl;
      cout << "doFitPion:            " << (doFitPion ? "TRUE" : "FALSE") << endl;
      cout << "doFitKaon:            " << (doFitKaon ? "TRUE" : "FALSE") << endl;
      cout << "doFitElectron:        " << (doFitElectron ? "TRUE" : "FALSE") << endl;
      cout << "fixKaonAmplitude:     " << (fixKaonAmplitude ? "TRUE" : "FALSE") << endl;
      cout << "fixElectronAmplitude: " << (fixElectronAmplitude ? "TRUE" : "FALSE") << endl;
      cout << "fixPionMean:          " << (fixPionMean ? "TRUE" : "FALSE") << endl;
      cout << "fixKaonMean:          " << (fixKaonMean ? "TRUE" : "FALSE") << endl;
      cout << "fixElectronMean:      " << (fixElectronMean ? "TRUE" : "FALSE") << endl;
      cout << "fixPionWidths:        " << (fixPionWidths ? "TRUE" : "FALSE") << endl;
      cout << "fixKaonWidths:        " << (fixKaonWidths ? "TRUE" : "FALSE") << endl;
      cout << "fixElectronWidths:    " << (fixElectronWidths ? "TRUE" : "FALSE") << endl;
      cout << "fixPionGamma:         " << (fixPionGamma ? "TRUE" : "FALSE") << endl;
    #endif



    //################  SETUP THE FIT FUNCTION ###################
    bool doElectronPionRatio = false;
    if(m_currentMomOfBin < electronEarlyStopMom)    doElectronPionRatio = true;
    //if(roundTag[fitting_round] == "electronRegion") doElectronPionRatio = true;
    if(m_currentPartIndex == 0) doElectronPionRatio = true;
    if(!doFitKaon) doElectronPionRatio = true;

    //    mu_pi sigma_pi gamma_pi  mu_k sigma_k  mu_e sigma_e  then For Each Cent: A_pi A_K e ... A_pi A_K e
    double parameters_simul_cent[75];
    for(int iii = 0; iii < 75; iii++) parameters_simul_cent[iii] = 0;
    parameters_simul_cent[0] = pionMeanInit;
    parameters_simul_cent[1] = pionSigmaInit;
    parameters_simul_cent[2] = pionSkewInit;
    parameters_simul_cent[3] = kaonMeanInit;
    parameters_simul_cent[4] = kaonSigmaInit;
    parameters_simul_cent[5] = electronMeanInit;
    parameters_simul_cent[6] = electronSigmaInit;
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      if(currentBinWasFitBefore){
        parameters_simul_cent[7 + centIndex*3]  = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        parameters_simul_cent[8 + centIndex*3]  = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        parameters_simul_cent[9 + centIndex*3]  = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{
        parameters_simul_cent[7 + centIndex*3] = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(pionMeanInit))    /PhysMath::skewNormal(pionMeanInit,pionSigmaInit,pionSkewInit,pionMeanInit);
        parameters_simul_cent[8 + centIndex*3] = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(kaonMeanInit))    /PhysMath::gaus(kaonMeanInit,kaonSigmaInit,kaonMeanInit);
        parameters_simul_cent[9 + centIndex*3] = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(electronMeanInit))/PhysMath::gaus(electronMeanInit,electronSigmaInit,electronMeanInit);
      }

      if(  roundTag[fitting_round] != "electronRegion" && roundTag[fitting_round] != "electronRegionStdFix" && electronRatioSuccess){
        if(doElectronPionRatio){
          parameters_simul_cent[9 + centIndex*3]  = m_electronRatioFitHistos_ByCent[centIndex][1]->GetBinContent(m_currentMtM0Bin);
          if(parameters_simul_cent[9 + centIndex*3] <= 0.0){
            parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[7 + centIndex*3];
          }
        }else{
          parameters_simul_cent[9 + centIndex*3]  = m_electronRatioFitHistos_ByCent[centIndex][2]->GetBinContent(m_currentMtM0Bin);
          if(parameters_simul_cent[9 + centIndex*3] <= 0.0){
            parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[8 + centIndex*3];
          }
        }
      }else{
        if(doElectronPionRatio){
          parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[7 + centIndex*3];
        }else{
          parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[8 + centIndex*3];
        }
      }
      if(!std::isfinite(parameters_simul_cent[9 + centIndex*3])){
        parameters_simul_cent[9 + centIndex*3] = 0.0;
      }
    }


    bool pionAmpsAllZero     = true;
    bool kaonAmpsAllZero     = true;
    bool electronAmpsAllZero = true;
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      if(parameters_simul_cent[7 + centIndex*3] > 0.0) pionAmpsAllZero     = false;
      if(parameters_simul_cent[8 + centIndex*3] > 0.0) kaonAmpsAllZero     = false;
      if(parameters_simul_cent[9 + centIndex*3] > 0.0) electronAmpsAllZero = false;
    }


    if(pionAmpsAllZero)     doFitPion     = false;
    if(kaonAmpsAllZero)     doFitKaon     = false;
    if(electronAmpsAllZero) doFitElectron = false;


    #ifdef _ZFITTER_DEBUG_
      cout << " Pi_Mean: " << parameters_simul_cent[0] << endl;
      cout << " Pi_Sig : " << parameters_simul_cent[1] << endl;
      cout << " Pi_Gam : " << parameters_simul_cent[2] << endl;
      cout << " K_Mean : " << parameters_simul_cent[3] << endl;
      cout << " K_Sig  : " << parameters_simul_cent[4] << endl;
      cout << " E_Mean : " << parameters_simul_cent[5] << endl;
      cout << " E_Sig  : " << parameters_simul_cent[6] << endl;
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        cout << " N_Pi_" << centIndex << " :" << parameters_simul_cent[7 + centIndex*3] << endl;
        cout << " N_K__" << centIndex << " :" << parameters_simul_cent[8 + centIndex*3] << endl;
        if(doElectronPionRatio){
          cout << " N_e/N_pi" << centIndex << " :" << parameters_simul_cent[9 + centIndex*3] << endl;
        }else{
          cout << " N_e/N_K" << centIndex << " :" << parameters_simul_cent[9 + centIndex*3] << endl;
        }
      }
    #endif


    bool initialParametersOkay = true;
    for(int iii = 0; iii < 7 + 3*m_numCentralities; iii++){
      if(!std::isfinite(parameters_simul_cent[iii])){
        initialParametersOkay = false;
        cout << "ERROR: Parameter " << iii << " = " << parameters_simul_cent[iii] << "     FAILED INITIALIZATION!!!!!" << endl;
        break;
      }
    }

    if(!initialParametersOkay){
      cout << "ERROR: Initial Parameters are not okay, probably infinite or NAN... Continuing... " << endl;
      continue;
    }



    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer();
    ROOT::Math::Functor* chiSqrdFunctor = NULL;
    if(doElectronPionRatio){
      chiSqrdFunctor = new ROOT::Math::Functor(this,&ZFitter::simultaneous_centrality_ZTPC_chisqr_E_Pi_Ratio_Linked,m_numCentralities*3 + 7);
    }else{
      chiSqrdFunctor = new ROOT::Math::Functor(this,&ZFitter::simultaneous_centrality_ZTPC_chisqr_E_K_Ratio_Linked,m_numCentralities*3 + 7);
    }
    minimizer->SetFunction(*chiSqrdFunctor);

    int printLevel = 0;
    #ifdef _ZFITTER_DEBUG_
      printLevel = 1;
    #endif
    minimizer->SetPrintLevel(printLevel);
    minimizer->SetMaxIterations(3000);
    minimizer->SetMaxFunctionCalls(200000);
    minimizer->SetValidError(true);

    //#################################   BASIC VARIABLE SETUP   ####################################
    string paramNames[55];
    double initialParams[55];
    double initialParamLimits[55][2];

    for(int iii = 0; iii < 55; iii++){
      initialParams[iii] = 0.0;
      initialParamLimits[iii][0] = 0.0;
      initialParamLimits[iii][1] = 0.0;
    }

    //PION PARAMETERS
    paramNames[0] = "#mu_{#pi}";
    initialParams[0] = pionMeanInit;
    minimizer->SetVariable(0, "#mu_{#pi}" ,pionMeanInit, 0.001);
    if(roundTag[fitting_round] == "electronRegion"){
      initialParamLimits[0][0] = pionMeanInit - 0.35;
      initialParamLimits[0][1] = pionMeanInit + 0.35;
      minimizer->SetVariableLimits(0, pionMeanInit - 0.35, pionMeanInit + 0.35);
    }else{
      initialParamLimits[0][0] = pionMeanInit - 0.15;
      initialParamLimits[0][1] = pionMeanInit + 0.15;
      minimizer->SetVariableLimits(0, pionMeanInit - 0.15, pionMeanInit + 0.15);
    }

    paramNames[1] = "#sigma_{#pi}";
    if(pionSigmaInit < 0.04) pionSigmaInit = 0.04;
    initialParams[1] = pionSigmaInit;
    minimizer->SetVariable(1, "#sigma_{#pi}" ,pionSigmaInit, 0.01);
    initialParamLimits[1][0] = 0.04;
    initialParamLimits[1][1] = 2.5*pionSigmaInit;
    minimizer->SetVariableLimits(1, 0.04, 2.5*pionSigmaInit);

    paramNames[2] = "#gamma_{#pi}";
    if(pionSkewInit > 0.99) pionSkewInit = 0.99;
    initialParams[2] = pionSigmaInit;
    minimizer->SetVariable(2, "#gamma_{#pi}" , pionSkewInit , 0.01);
    initialParamLimits[2][0] = 0.0;
    initialParamLimits[2][1] = 0.99;
    minimizer->SetVariableLimits(2, 0.0,0.99);


    //KAON PARAMETERS
    paramNames[3] = "#mu_{K}";
    initialParams[3] = kaonMeanInit;
    minimizer->SetVariable(3, "#mu_{K}" ,kaonMeanInit, 0.001);
    if(roundTag[fitting_round] == "electronRatioFixed" && m_currentPartIndex == 0 && m_currentMomOfBin > 0.3){
      initialParamLimits[3][0] = kaonMeanInit - 0.03;
      initialParamLimits[3][1] = kaonMeanInit + 0.03;
      minimizer->SetVariableLimits(3, kaonMeanInit - 0.03, kaonMeanInit + 0.03);
    }else{
      initialParamLimits[3][0] = kaonMeanInit - 0.15;
      initialParamLimits[3][1] = kaonMeanInit + 0.15;
      minimizer->SetVariableLimits(3, kaonMeanInit - 0.15, kaonMeanInit + 0.15);
    }
    //minimizer->SetVariableLimits(3, initialParamLimits[3][0], initialParamLimits[3][1]);

    paramNames[4] = "#sigma_{K}";
    if(kaonSigmaInit < 0.04) kaonSigmaInit = 0.04;
    if(kaonSigmaInit > 0.12) kaonSigmaInit = 0.12;
    initialParams[4] = kaonSigmaInit;
    if(m_currentPartIndex == 1 && m_currentMomOfBin < 0.3){
      initialParamLimits[4][0] = 0.04;
      initialParamLimits[4][1] = 0.12;
    }else{
      initialParamLimits[4][0] = 0.04;
      initialParamLimits[4][1] = 0.12;
    }
    minimizer->SetVariable(4, "#sigma_{K}" ,kaonSigmaInit, 0.01);
    minimizer->SetVariableLimits(4, initialParamLimits[4][0], initialParamLimits[4][1]);

    //ELECTRON PARAMETERS
    paramNames[5] = "#mu_{e}";
    initialParams[5] = electronMeanInit;
    minimizer->SetVariable(5, "#mu_{e}" ,electronMeanInit, 0.001);
    if(m_currentPartIndex == 0 && pionBichselPredict > electronBichselPredict){
      initialParamLimits[5][0] = electronMeanInit - 0.35;
      initialParamLimits[5][1] = electronMeanInit + 0.15;
      minimizer->SetVariableLimits(5, electronMeanInit - 0.35, electronMeanInit + 0.15);
    }else{
      initialParamLimits[5][0] = electronMeanInit - 0.15;
      initialParamLimits[5][1] = electronMeanInit + 0.15;
      minimizer->SetVariableLimits(5, electronMeanInit - 0.15, electronMeanInit + 0.15);
    }
    if(m_currentMomOfBin > 0.3){
      if(electronSigmaInit < 0.04) electronSigmaInit = 0.04;
      if(electronSigmaInit > 0.1) electronSigmaInit = 0.1;
      initialParamLimits[6][0] = 0.04;
      initialParamLimits[6][1] = 0.1;
    }else{
      if(electronSigmaInit < 0.04) electronSigmaInit = 0.04;
      if(electronSigmaInit > 0.2) electronSigmaInit = 0.2;
      initialParamLimits[6][0] = 0.04;
      initialParamLimits[6][1] = 0.2;

    }

    paramNames[6] = "#sigma_{e}";
    initialParams[6] = electronSigmaInit;
    minimizer->SetVariable(6, "#sigma_{e}" ,electronSigmaInit, 0.001);
    minimizer->SetVariableLimits(6, initialParamLimits[6][0], initialParamLimits[6][1]);







    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

      paramNames[7+centIndex*3] = Form("N_{#pi}_%d",centIndex);
      if(parameters_simul_cent[7 + centIndex*3] <= 0.0){
        initialParams[7+centIndex*3] = 0.0;
        initialParamLimits[7+centIndex*3][0] = 0.0;
        initialParamLimits[7+centIndex*3][1] = 1.0;
        minimizer->SetVariable(7 + centIndex*3, Form("N_{#pi}_%d",centIndex), 0.0, 0.0001);
        minimizer->SetVariableLimits(7 + centIndex*3, 0.0, 1.0);
      }else{
        initialParams[7+centIndex*3] = parameters_simul_cent[7 + centIndex*3];
        initialParamLimits[7+centIndex*3][0] = 0.0001*parameters_simul_cent[7 + centIndex*3];
        initialParamLimits[7+centIndex*3][1] = 100.0*parameters_simul_cent[7 + centIndex*3];
        minimizer->SetVariable(7 + centIndex*3, Form("N_{#pi}_%d",centIndex), parameters_simul_cent[7 + centIndex*3], 0.01*parameters_simul_cent[7 + centIndex*3]);
        minimizer->SetVariableLimits(7 + centIndex*3, 0.0001*parameters_simul_cent[7 + centIndex*3], 100.0*parameters_simul_cent[7 + centIndex*3]);
      }

      paramNames[8+centIndex*3] = Form("N_{K}_%d",centIndex);
      if(parameters_simul_cent[8 + centIndex*3] <= 0.0){
        initialParams[8+centIndex*3] = 0.0;
        initialParamLimits[8+centIndex*3][0] = 0.0;
        initialParamLimits[8+centIndex*3][1] = 1.0;
        minimizer->SetVariable(8 + centIndex*3, Form("N_{K}_%d",centIndex), 0.0, 0.0001);
        minimizer->SetVariableLimits(8 + centIndex*3, 0.0, 1.0);
      }else{
        initialParams[8+centIndex*3] = parameters_simul_cent[8+centIndex*3];
        initialParamLimits[8+centIndex*3][0] = 0.0001*parameters_simul_cent[8 + centIndex*3];
        initialParamLimits[8+centIndex*3][1] = 100.0*parameters_simul_cent[8 + centIndex*3];
        minimizer->SetVariable(8 + centIndex*3, Form("N_{K}_%d",centIndex), parameters_simul_cent[8 + centIndex*3], 0.01*parameters_simul_cent[8 + centIndex*3]);
        minimizer->SetVariableLimits(8 + centIndex*3, 0.0001*parameters_simul_cent[8 + centIndex*3], 100.0*parameters_simul_cent[8 + centIndex*3]);
      }
      //double electronPionRatioInit = m_electronRatioFitHistos_ByCent[centIndex][1]->GetBinContent(m_currentMtM0Bin);
      //double electronKaonRatioInit = m_electronRatioFitHistos_ByCent[centIndex][2]->GetBinContent(m_currentMtM0Bin);

      if(doElectronPionRatio){
        paramNames[9+centIndex*3] = Form("N_{e}/N_{#pi}_%d",centIndex);
        if(parameters_simul_cent[9 + centIndex*3] <= 0.0){
          initialParams[9+centIndex*3] = 0.0;
          initialParamLimits[9+centIndex*3][0] = 0.0;
          initialParamLimits[9+centIndex*3][1] = 15.0;
          minimizer->SetVariable(9 + centIndex*3, Form("N_{e}/N_{#pi}_%d",centIndex), 0.0, 0.0001);
          minimizer->SetVariableLimits(9 + centIndex*3,0.0, 15.0);
        }else{
          initialParams[9+centIndex*3] = parameters_simul_cent[8+centIndex*3];
          initialParamLimits[9+centIndex*3][0] = 0.0001*parameters_simul_cent[9 + centIndex*3];
          initialParamLimits[9+centIndex*3][1] = 100.0*parameters_simul_cent[9 + centIndex*3];
          minimizer->SetVariable(9 + centIndex*3, Form("N_{e}/N_{#pi}_%d",centIndex), parameters_simul_cent[9 + centIndex*3], 0.01*parameters_simul_cent[9 + centIndex*3]);
          minimizer->SetVariableLimits(9 + centIndex*3, 0.0001*parameters_simul_cent[9 + centIndex*3], 100.0*parameters_simul_cent[9 + centIndex*3]);
        }
      }else{
        paramNames[9+centIndex*3] = Form("N_{e}/N_{K}_%d",centIndex);
        if(parameters_simul_cent[9 + centIndex*3] <= 0.0){
          initialParams[9+centIndex*3] = 0.0;
          initialParamLimits[9+centIndex*3][0] = 0.0;
          initialParamLimits[9+centIndex*3][1] = 15.0;
          minimizer->SetVariable(9 + centIndex*3, Form("N_{e}/N_{K}_%d",centIndex), 0.0, 0.0001);
          minimizer->SetVariableLimits(9 + centIndex*3, 0.0, 15.0);
        }else{
          initialParams[9+centIndex*3] = parameters_simul_cent[8+centIndex*3];
          initialParamLimits[9+centIndex*3][0] = 0.0001*parameters_simul_cent[9 + centIndex*3];
          initialParamLimits[9+centIndex*3][1] = 100.0*parameters_simul_cent[9 + centIndex*3];
          minimizer->SetVariable(9 + centIndex*3, Form("N_{e}/N_{K}_%d",centIndex), parameters_simul_cent[9 + centIndex*3], 0.01*parameters_simul_cent[9 + centIndex*3]);
          minimizer->SetVariableLimits(9 + centIndex*3, 0.0001*parameters_simul_cent[9 + centIndex*3], 100.0*parameters_simul_cent[9 + centIndex*3]);
        }
      }
    }


    //-----------------------      FIXING PARAMETERS     ------------------------------
    //PION
    if(fixPionMean)       minimizer->FixVariable(0);
    if(fixPionWidths)     minimizer->FixVariable(1);
    if(fixPionGamma)      minimizer->FixVariable(2);
    if(!doFitPion){
      minimizer->SetFixedVariable(0, "#mu_{#pi}" ,0.0);
      minimizer->SetFixedVariable(1, "#sigma_{#pi}" ,1.0);
      minimizer->SetFixedVariable(2, "#gamma_{#pi}" , 0.1);
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->SetFixedVariable(7 + centIndex*3, Form("N_{#pi}_%d",centIndex), 0.0);
        if(doElectronPionRatio){
          minimizer->SetFixedVariable(9 + centIndex*3, Form("N_{e}/N_{#pi}_%d",centIndex), 0.0);
        }
      }
    }

    //KAON
    if(fixKaonAmplitude){
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->FixVariable(8 + centIndex*3);
      }
    }
    if(fixKaonMean)       minimizer->FixVariable(3);
    if(fixKaonWidths)     minimizer->FixVariable(4);
    if(!doFitKaon){
      minimizer->SetFixedVariable(3,"#mu_{K}",   0.0);
      minimizer->SetFixedVariable(4,"#sigma_{K}",1.0);
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->SetFixedVariable(8 + centIndex*3, Form("N_{K}_%d",centIndex), 0.0);
        if(!doElectronPionRatio){
          minimizer->SetFixedVariable(9 + centIndex*3, Form("N_{e}/N_{K}_%d",centIndex), 0.0);
        }
      }
    }

    //ELECTRON
    if(fixElectronAmplitude){   // ratio e/pi or e/K
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->FixVariable(9 + centIndex*3);
      }
    }
    if(fixElectronMean)   minimizer->FixVariable(5);
    if(fixElectronWidths) minimizer->FixVariable(6);
    if(!doFitElectron){
      minimizer->SetFixedVariable(5,"#mu_{e}",   0.0);
      minimizer->SetFixedVariable(6,"#sigma_{e}",1.0);
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->SetFixedVariable(9 + centIndex*3, Form("N_{e}_%d",centIndex), 0.0);
      }
    }


    //Aggregate initial parameters and limits
    #ifdef _ZFITTER_DEBUG_
      cout << "Grabbing parameter limits" << endl;
    #endif

    /*
    //string paramNames[55];

    //double initialParams[55];
    //double initialParamLimits[55][2];
    for(int iii = 0; iii < 7 + 3*m_numCentralities; iii++){
      ROOT::Fit::ParameterSettings* setting = new ROOT::Fit::ParameterSettings();
      minimizer->GetVariableSettings(iii, *setting);
      //initialParamLimits[iii][0] = setting->LowerLimit();
      //initialParamLimits[iii][1] = setting->UpperLimit();
      //initialParams[iii]         = setting->Value();
      paramNames[iii]            = setting->Name();

      if(setting->Value() != initialParams[iii]){
        cout << "WARNING: Parameter " << iii << " " << setting->Name() << " has a different value than expected" << endl;
        cout << "  Expected: " << initialParams[iii] << "  Found: " << setting->Value()
             << "  with limits [ " << setting->LowerLimit() << ", " << setting->UpperLimit() << " ]" << endl;
      }
      #ifdef _ZFITTER_DEBUG_
        cout << " Param " << iii << "  " << initialParams[iii] << "   [ " << initialParamLimits[iii][0] << ", " << initialParamLimits[iii][1] << " ] " << endl;
      #endif
      delete setting;
    }*/

    #ifdef _ZFITTER_DEBUG_
      cout << "PARAMETERS BEFORE PARAM FIXING --------------------" << endl;
      for(int iii = 0; iii < 7 + 3*m_numCentralities; iii++){
        cout << " Param " << iii << "  " << paramNames[iii] << "  " << initialParams[iii] << "   [ " << initialParamLimits[iii][0] << ", " << initialParamLimits[iii][1] << " ] " << endl;
      }
    #endif




    //###############################   DO THE MINIMIZATION   #######################################
    minimizer->Minimize();
    int status_Minimizer = minimizer->Status();
    const double* finalParams = minimizer->X();
    if(status_Minimizer == 0 || status_Minimizer == 4) minimizer->Hesse(); // HESSE CAN MAKE free() invalid pointer errors for failed fits BE CAREFUL!!!!!------------------------------------------------------------------------------------------------------------------------------
    m_FitStatus_ZTPC[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status_Minimizer);
    #ifdef _ZFITTER_DEBUG_
      cout << "simultaneous centrality minimizer status: " << status_Minimizer << endl;
    #endif



    double simulParams[55];
    double simulParamErrors[55];
    for(int iii = 0; iii < 7; iii++){
      simulParams[iii]      = finalParams[iii];
      simulParamErrors[iii] = sqrt(minimizer->CovMatrix(iii,iii));
    }
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      simulParams[7+3*centIndex]      = finalParams[7+3*centIndex];
      simulParamErrors[7+3*centIndex] = sqrt(minimizer->CovMatrix(7+3*centIndex, 7+3*centIndex));
      simulParams[8+3*centIndex]      = finalParams[8+3*centIndex];
      simulParamErrors[8+3*centIndex] = sqrt(minimizer->CovMatrix(8+3*centIndex, 8+3*centIndex));
      if(doElectronPionRatio){
        simulParams[9+3*centIndex]      = finalParams[9+3*centIndex] * finalParams[7+3*centIndex];
        simulParamErrors[9+3*centIndex] = simulParams[9+3*centIndex] * sqrt(minimizer->CovMatrix(7+3*centIndex, 7+3*centIndex)/(finalParams[7+3*centIndex] * finalParams[7+3*centIndex])
                                              + minimizer->CovMatrix(9+3*centIndex, 9+3*centIndex)/(finalParams[9+3*centIndex] * finalParams[9+3*centIndex]));
      }else{
        simulParams[9+3*centIndex]      = finalParams[9+3*centIndex] * finalParams[8+3*centIndex];
        simulParamErrors[9+3*centIndex] = simulParams[9+3*centIndex] * sqrt(minimizer->CovMatrix(8+3*centIndex, 8+3*centIndex)/(finalParams[8+3*centIndex] * finalParams[8+3*centIndex])
                                              + minimizer->CovMatrix(9+3*centIndex, 9+3*centIndex)/(finalParams[9+3*centIndex] * finalParams[9+3*centIndex]));
      }
    }

    if(status_Minimizer == 0 || status_Minimizer == 4){
      if(doFitPion){
        previousGoodParameters[0] = simulParams[0] - log(m_bichselCurves[0]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        previousGoodParameters[1] = simulParams[1];
        previousGoodParameters[2] = simulParams[2];
      }else{
        previousGoodParameters[0] = -999;
        previousGoodParameters[1] = -999;
        previousGoodParameters[2] = -999;
      }
      if(doFitKaon){
        previousGoodParameters[3] = simulParams[3] - log(m_bichselCurves[1]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        previousGoodParameters[4] = simulParams[4];
      }else{
        previousGoodParameters[3] = -999;
        previousGoodParameters[4] = -999;
      }
      if(doFitElectron){
        previousGoodParameters[5] = simulParams[5] - log(m_bichselCurves[3]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        previousGoodParameters[6] = simulParams[6];
      }else{
        previousGoodParameters[5] = -999;
        previousGoodParameters[6] = -999;
      }
    }else{
      for(int iii = 0; iii < 11; iii++) previousGoodParameters[iii] = -999;
    }

    TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
    fittingCanvas->SetWindowSize(1300,700);
    fittingCanvas->SetCanvasSize(1200,600);
    fittingCanvas->Draw();
    TPad* topPad = new TPad("topPad_PiKTPC","topPad_PiKTPC",0.0,0.3,1.0,1.0);
    TPad* bottomPad = new TPad("bottomPad_PiKTPC","bottomPad_PiKTPC",0.0,0.0,1.0,0.3);
    topPad->Draw();
    bottomPad->Draw();
    fittingCanvas->cd();
    topPad->cd();
    topPad->SetLogy(true);
    topPad->SetRightMargin(0.4);
    bottomPad->SetLogy(false);
    bottomPad->SetRightMargin(0.4);
    gStyle->SetOptFit(0111);

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      double numberOfCentEvents = m_centEvents[m_currentPartIndex]->GetBinContent(centIndex+1);
      m_currentHistosToFit_ByCent[centIndex]->SetMarkerStyle(20);
      m_currentHistosToFit_ByCent[centIndex]->SetMarkerColor(kBlack);
      m_currentHistosToFit_ByCent[centIndex]->SetMarkerSize(0.5);
      m_currentHistosToFit_ByCent[centIndex]->SetLineWidth(2);
      m_currentHistosToFit_ByCent[centIndex]->SetLineColor(kBlack);


      #ifndef _ZFITTER_PHYSMATH_BYPASS_
        TF1* combinedFunct = PhysMath::getSkewAndTwoGausFunct(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d",fitting_round,
          m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
          simulParams[7 + centIndex*3], simulParams[0],simulParams[1],simulParams[2],
          simulParams[8 + centIndex*3], simulParams[3], simulParams[4],
          simulParams[9 + centIndex*3], simulParams[5], simulParams[6]);
        TF1* pionFunct  = PhysMath::skewNormalFunct( simulParams[7 + centIndex*3], simulParams[0], simulParams[1], simulParams[2]);
        TF1* kaonFunct  = PhysMath::getGaussianFunction(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d_Kaon",fitting_round,
          m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin), simulParams[8 + centIndex*3], simulParams[3], simulParams[4]);
        TF1* electronFunct = PhysMath::getGaussianFunction(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d_Electron",fitting_round,
          m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin), simulParams[9 + centIndex*3], simulParams[5], simulParams[6]);
        pionFunct->SetLineColor(m_partInfo->GetParticleColor(0));
        pionFunct->SetRange(lowRange,highRange);
        electronFunct->SetLineColor(m_partInfo->GetParticleColor(3));
        electronFunct->SetRange(lowRange,highRange);
        kaonFunct->SetLineColor(m_partInfo->GetParticleColor(1));
        kaonFunct->SetRange(lowRange,highRange);
        TF1* pionFunct_init = NULL;
        TF1* koanFunct_init = NULL;
        TF1* electronFunct_init = NULL;
        if(m_drawInitialSeedsToFits){
          pionFunct_init  = PhysMath::skewNormalFunct( parameters_simul_cent[7 + 3*centIndex], parameters_simul_cent[0], parameters_simul_cent[1], parameters_simul_cent[2]);
          koanFunct_init  = PhysMath::getGaussianFunction(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d_Kaon_init",fitting_round,
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin), parameters_simul_cent[8 + 3*centIndex],parameters_simul_cent[3], parameters_simul_cent[4]);
          double electronN = parameters_simul_cent[9 + 3*centIndex] * parameters_simul_cent[8 + 3*centIndex];
          if(doElectronPionRatio) electronN = parameters_simul_cent[9 + 3*centIndex] * parameters_simul_cent[7 + 3*centIndex];
          electronFunct_init = PhysMath::getGaussianFunction(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d_Electron_init",fitting_round,
            m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin), electronN, parameters_simul_cent[5], parameters_simul_cent[6]);
          pionFunct_init->SetLineColor(kYellow+2);
          pionFunct_init->SetRange(lowRange,highRange);
          electronFunct_init->SetLineColor(kOrange);
          electronFunct_init->SetRange(lowRange,highRange);
          koanFunct_init->SetLineColor(kYellow+1);
          koanFunct_init->SetRange(lowRange,highRange);
        }
      #endif
      double highestValue  = HistogramUtilities::getMaxInRange(m_currentHistosToFit_ByCent[centIndex],m_currentLowFitRange,m_currentHighFitRange);
      double lowestValue   = HistogramUtilities::getMinInRange(m_currentHistosToFit_ByCent[centIndex],m_currentLowFitRange,m_currentHighFitRange);
      TLine* lowRangeLine  = new TLine(m_currentLowFitRange,lowestValue,m_currentLowFitRange,highestValue);
      TLine* highRangeLine = new TLine(m_currentHighFitRange,lowestValue,m_currentHighFitRange,highestValue);
      lowRangeLine->SetLineColor(kViolet);
      highRangeLine->SetLineColor(kViolet);



      topPad->cd();
      HistogramUtilities::setAxisRanges(m_currentHistosToFit_ByCent[centIndex], lowRange - 0.5*sigmaGuess, highRange + 0.5*sigmaGuess,true);
      topPad->SetLogy(true);
      m_currentHistosToFit_ByCent[centIndex]->SetStats(kFALSE);
      m_currentHistosToFit_ByCent[centIndex]->Draw("E");
      #ifndef _ZFITTER_PHYSMATH_BYPASS_
        if(pionFunct_init && doFitPion) pionFunct_init->Draw("same");
        if(koanFunct_init && doFitKaon) koanFunct_init->Draw("same");
        if(electronFunct_init && doFitElectron) electronFunct_init->Draw("same");
        combinedFunct->Draw("same");
        pionFunct->Draw("same");
        electronFunct->Draw("same");
        kaonFunct->Draw("same");
      #endif
      lowRangeLine->Draw("same");
      highRangeLine->Draw("same");




      TPaveText* fitParsTxt = new TPaveText(0.65,0.05,0.95,0.90,"NDC");

      fitParsTxt->AddText(Form("N_{#pi} \t %e#pm%e %s",simulParams[7 + centIndex*3],simulParamErrors[7 + centIndex*3],minimizer->IsFixedVariable(7 + 3*centIndex) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#mu_{#pi} \t %f#pm%f %s",simulParams[0],simulParamErrors[0],minimizer->IsFixedVariable(0) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#sigma_{#pi} \t %f#pm%f %s",simulParams[1],simulParamErrors[1],minimizer->IsFixedVariable(1) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#gamma_{#pi} \t %f#pm%f %s",simulParams[2],simulParamErrors[2],minimizer->IsFixedVariable(2) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("N_{K} \t %e#pm%e %s",simulParams[8 + centIndex*3], simulParamErrors[8 + centIndex*3],minimizer->IsFixedVariable(8 + 3*centIndex) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#mu_{K} \t %f#pm%f %s",simulParams[3],simulParamErrors[3],minimizer->IsFixedVariable(3) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#sigma_{K} \t %f#pm%f %s",simulParams[4],simulParamErrors[4],minimizer->IsFixedVariable(4) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("N_{e} \t %e#pm%e %s",simulParams[9 + centIndex*3], simulParamErrors[9 + centIndex*3],minimizer->IsFixedVariable(9 + centIndex*3) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#mu_{e} \t %f#pm%f %s",simulParams[5],simulParamErrors[5],minimizer->IsFixedVariable(5) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#sigma_{e} \t %f#pm%f %s",simulParams[6],simulParamErrors[6],minimizer->IsFixedVariable(6) ? "(fixed)" : "       "));
      fitParsTxt->SetFillStyle(0);
      fitParsTxt->SetBorderSize(5);
      fitParsTxt->Draw("same");





      bottomPad->cd();
      #ifndef _ZFITTER_PHYSMATH_BYPASS_
      TH1D* pullHisto = HistogramUtilities::makeFitPullHistogram(m_currentHistosToFit_ByCent[centIndex], combinedFunct, true);
      pullHisto->SetMarkerStyle(20);
      pullHisto->SetMarkerSize(0.5);
      pullHisto->SetMarkerColor(kBlack);
      pullHisto->SetLineColor(kBlack);
      #endif

      TPaveText* fitTxt = new TPaveText(0.65,0.05,0.95,0.95,"NDC");
      if(initialParams[7 + 3*centIndex] != 0.0) fitTxt->AddText(Form("N_{#pi} \t %.2e  [ %.2e, %.2e ]",initialParams[7 + 3*centIndex],initialParamLimits[7 + 3*centIndex][0],initialParamLimits[7 + 3*centIndex][1]));
      if(initialParams[8 + 3*centIndex] != 0.0) fitTxt->AddText(Form("N_{K} \t %.2e  [ %.2e, %.2e ]",  initialParams[8 + 3*centIndex],initialParamLimits[8 + 3*centIndex][0],initialParamLimits[8 + 3*centIndex][1]));
      if(initialParams[9 + 3*centIndex] != 0.0) fitTxt->AddText(Form("N_{e} \t %.2e  [ %.2e, %.2e ]",  initialParams[9 + 3*centIndex],initialParamLimits[9 + 3*centIndex][0],initialParamLimits[9 + 3*centIndex][1]));
      for(int iii = 0; iii < 7; iii++){
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



      TH1F* frame_pull = bottomPad->DrawFrame(m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( lowRange - 0.5*sigmaGuess)),-15.0,
             m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( highRange + 0.5*sigmaGuess)),15.0);
      frame_pull->SetTitle(Form("; Z_{TPC}; Standardized Pull"));
      #ifndef _ZFITTER_PHYSMATH_BYPASS_
      pullHisto->Draw("Same");
      #endif
      TLine* oneLine = new TLine(lowRange - 0.5*sigmaGuess,0.0,highRange + 0.5*sigmaGuess,0.0);
      oneLine->SetLineStyle(7);
      oneLine->SetLineColor(kBlack);
      oneLine->Draw("same");
      fitTxt->Draw("same");

      gPad->Update();
      gSystem->ProcessEvents();
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_STEP_%02d_mTm0Bin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
             m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex, m_currentRapBin, fitting_round, m_currentMtM0Bin));

      if(m_saveNoLogImages){
        topPad->SetLogy(false);
        topPad->Update();
        gSystem->ProcessEvents();
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_STEP_%02d_mTm0Bin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
               m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex,m_currentRapBin,fitting_round, m_currentMtM0Bin));
      }
      if(centIndex == 0){ // only write the all cent data once
        //PION
        m_Fit_Data_ZTPC[m_currentPartIndex][0][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[0]);
        m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[1]);
        m_Fit_Data_ZTPC[m_currentPartIndex][0][3][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[2]);
        m_Fit_Data_ZTPC[m_currentPartIndex][0][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[0]);
        m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[1]);
        m_Fit_Data_ZTPC[m_currentPartIndex][0][3][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[2]);
        //KAON
        m_Fit_Data_ZTPC[m_currentPartIndex][1][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[3]);
        m_Fit_Data_ZTPC[m_currentPartIndex][1][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[4]);
        m_Fit_Data_ZTPC[m_currentPartIndex][1][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[3]);
        m_Fit_Data_ZTPC[m_currentPartIndex][1][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[4]);
        //ELECTRON
        m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[5]);
        m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[6]);
        m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[5]);
        m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[6]);
      }
      //PION
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[7 + 3*centIndex]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[0]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[1]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][3][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[2]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,   simulParamErrors[7 + 3*centIndex]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[0]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[1]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][3][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[2]);
      //KAON
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[8 + 3*centIndex]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[3]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[4]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[8 + 3*centIndex]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[3]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[4]);
      //ELECTRON
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[9 + 3*centIndex]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[5]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[6]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[9 + 3*centIndex]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[5]);
      m_Fit_Data_Cent_ZTPC[m_currentPartIndex][3][centIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, simulParamErrors[6]);
      if(m_currentPartIndex == 0){
        m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[7 + 3*centIndex]*numberOfCentEvents);
        m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,   simulParamErrors[7 + 3*centIndex]*numberOfCentEvents);
      }
      if(m_currentPartIndex == 1){
        m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, simulParams[8 + 3*centIndex]*numberOfCentEvents);
        m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,   simulParamErrors[8 + 3*centIndex]*numberOfCentEvents);
      }
      #ifndef _SAVE_POINTERS_
        #ifndef _ZFITTER_PHYSMATH_BYPASS_
          delete combinedFunct;
          delete pionFunct;
          delete electronFunct;
          delete kaonFunct;
          if(pionFunct_init) delete pionFunct_init;
          if(koanFunct_init) delete koanFunct_init;
          if(electronFunct_init) delete electronFunct_init;
          delete pullHisto;
        #endif

        delete lowRangeLine;
        delete highRangeLine;
        delete oneLine;
        delete fitParsTxt;
        delete fitTxt;
      #endif
    } // end centrality loop
    //minimizer->Clear();
    #ifndef _SAVE_POINTERS_
      delete topPad;
      delete bottomPad;
      delete fittingCanvas;
      delete chiSqrdFunctor;
      delete minimizer;
    #endif



    mTm0BinsFitPerRound[fitting_round] = mTm0BinsFitPerRound[fitting_round] + 1;

  }// End mTm0 Loop
  for(int centIndex = 0; centIndex < 16; centIndex++){
    draw_ZVar_mTm0_data_with_fit_data(m_currentPartIndex,0,centIndex,m_currentRapBin-1,m_currentPartCharge);
  }


}



#endif
