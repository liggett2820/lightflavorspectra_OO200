// THIS IS INSERTED INTO ZFITTER.cxx // Just separate file becuase it's long
//============================================================================================================
//#########      TPC FITTING FOR PROTON ONLY   #####################
//============================================================================================================
#ifndef _ZFITTER_PROTONTPC_SIMULCENT_SINGLELOOP_
#define _ZFITTER_PROTONTPC_SIMULCENT_SINGLELOOP_

double ZFitter::simultaneous_centrality_ZTPC_skewNormalAndTwoGaus_Proton_chisqr(const double* a_param){
  //    mu_pi sigma_pi  mu_K sigma_K    mu_p sigma_p gamma_p then For Each Cent: A_pi A_K A_P ... A_pi A_K A_P 
  double chisqr = 0.0;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    for(int binX = m_currentLowFitBin[centIndex]; binX <= m_currentHighFitBin[centIndex]; binX++){
      double content = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(binX);
      if(content > 0.0){
          double center  = m_currentHistosToFit_ByCent[centIndex]->GetBinCenter(binX);
          double prediction =   a_param[7 + 3*centIndex]*PhysMath::gaus(a_param[0],a_param[1],center)
                              + a_param[8 + 3*centIndex]*PhysMath::gaus(a_param[2],a_param[3],center)
                              + a_param[9 + 3*centIndex]*PhysMath::skewNormal(a_param[4],a_param[5],a_param[6],center);
          chisqr += pow((prediction - content)/m_currentHistosToFit_ByCent[centIndex]->GetBinError(binX),2);
      }// content > 0
    }//bin loop
  }// cent loop
  return chisqr;
}

double ZFitter::simultaneous_centrality_ZTPC_ThreeSkewNormal_Proton_chisqr(const double* a_param){
  //    mu_pi sigma_pi gamma_pi   mu_K sigma_K gamma_pi    mu_p sigma_p gamma_p then For Each Cent: A_pi A_K A_P ... A_pi A_K A_P 
  double chisqr = 0.0;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    for(int binX = m_currentLowFitBin[centIndex]; binX <= m_currentHighFitBin[centIndex]; binX++){
      double content = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(binX);
      if(content > 0.0){
          //double center  = m_currentHistosToFit_ByCent[centIndex]->GetBinCenter(binX);
          //double prediction =   a_param[9  + 3*centIndex]*PhysMath::skewNormal(a_param[0],a_param[1],a_param[2],center)
          //                    + a_param[10 + 3*centIndex]*PhysMath::skewNormal(a_param[3],a_param[4],a_param[5],center)
          //                    + a_param[11 + 3*centIndex]*PhysMath::skewNormal(a_param[6],a_param[7],a_param[8],center);
          double low  =  m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(binX);
          double high =  m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(binX+1);
          double prediction =   a_param[9  + 3*centIndex]*PhysMath::skewNormalBinContent(a_param[0],a_param[1],a_param[2],low,high)
                              + a_param[10 + 3*centIndex]*PhysMath::skewNormalBinContent(a_param[3],a_param[4],a_param[5],low,high)
                              + a_param[11 + 3*centIndex]*PhysMath::skewNormalBinContent(a_param[6],a_param[7],a_param[8],low,high);
          chisqr += pow((prediction - content)/m_currentHistosToFit_ByCent[centIndex]->GetBinError(binX),2);
      }// content > 0
    }//bin loop
  }// cent loop
  return chisqr;
}




//============================================================================================================================================================================


void ZFitter::fitTPCProton_SimulCent_RapidityLooper(int a_charge){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::fitTPCProton_SimulCent_RapidityLooper(" << a_charge << ")" << endl;
  #endif
  m_currentCentIndex  = -1;
  m_inAllCentMode     = true;
  m_currentPartIndex  = 2;
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
  // Fixed 2026-07: was `rightBinExists && leftBinExists` (AND). rapZeroBin is
  // FindBin(0.0) on the ACTUAL loaded data histogram's axis (immutable, baked in by
  // etof's PicoBinner at yield-file creation time) -- for this O+O 200 GeV proton yield
  // file, y=0 lands in the LAST bin of a 31-bin axis (-3.03 to 0.07), i.e. right at the
  // edge. With AND, the very first step in the direction with no room (right, in this
  // case) sets that side's *Exists flag false, and the very next loop-condition check
  // kills the ENTIRE loop -- even though the other direction (left) still has ~14 more
  // populated, unexplored bins (confirmed via InspectRawYieldRapBins.C: 15 populated
  // bins total in the raw yield data, y=-1.38 to 0.02). OR lets the loop keep walking
  // whichever side still has room; the side that ran out just keeps re-setting its own
  // (already-false) flag every other iteration, a harmless no-op, until the other side
  // also exhausts itself and the loop ends naturally. Confirmed byte-identical to
  // etof's copy of this file before this fix (0 diff lines) -- this is a real,
  // pre-existing ZFitter algorithm bug, not something introduced by porting. It never
  // manifested in etof's own collider-mode analyses because their rapidity-zero bin
  // was never adjacent to an axis edge; etof's BTOF equivalent even has a special-cased
  // `if(getEventConfig()=="FixedTarget") rapZeroBin = GetNbinsX()/2;` workaround for a
  // different manifestation of this same edge case in fixed-target mode -- confirming
  // this exact problem was known before, just never fixed at the loop-logic level.
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
    fitTPCProton_SimulCent_ByRapidity();
  }
}


//============================================================================================================================================================================


void ZFitter::calculateProtonMeanShifts(double &a_delta_Pi,double &a_delta_K,double &a_delta_P){

  for(int partIndex = 0; partIndex < 3; partIndex++){
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
      if(std::isfinite(mean)){
      if(partIndex == 0) a_delta_Pi = mean;
      if(partIndex == 1) a_delta_K  = mean;
      if(partIndex == 2) a_delta_P  = mean;
    }
  }

  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::calculateProtonMeanShifts()" << endl;
    cout << "       Delta Mu Pi : " << a_delta_Pi << endl;
    cout << "       Delta Mu K  : " << a_delta_K << endl;
    cout << "       Delta Mu P  : " << a_delta_P << endl;
  #endif

  return;
}

//============================================================================================================================================================================


bool ZFitter::calculateProtonSkewExtrapolation(double &a_skew){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::calculateProtonSkewExtrapolation()" << endl;
  #endif

  TH2D* histo = m_Fit_Data_ZTPC[m_currentPartIndex][2][3][m_currentPlusMinusIndex];
  double weights = 0.0;
  double mean    = 0.0;

  int binsUsed = 0;
  for(int binY = histo->GetNbinsY(); binY > 0; binY--){
    if(binsUsed > 5) break;
    double mTm0      = histo->GetYaxis()->GetBinCenter(binY);
    double momentum  = PhysMath::pTot_from_y_mTm0(m_currentRapidity,mTm0,m_partInfo->GetParticleMass(m_currentPartIndex));
    double binCount  = histo->GetBinContent(m_currentRapBin,binY) - log(m_bichselCurves[2]->Eval(momentum)/m_bichselCurves[m_currentPartIndex]->Eval(momentum));
    double binWeight = 1.0/(histo->GetBinError(m_currentRapBin,binY)*histo->GetBinError(m_currentRapBin,binY));
    if(!std::isfinite(binWeight) || !std::isfinite(binCount)) continue;
    binsUsed++;
    mean    += binWeight*binCount;
    weights += binWeight;
  }
  a_skew = mean/weights;
  #ifdef _ZFITTER_DEBUG_
    cout << "       Skew = " << a_skew << " from mu: " << mean << " and weights: " << weights << endl;
  #endif

  bool goodEstimate = true;
  if(!std::isfinite(a_skew) || fabs(a_skew) > 0.99){
    a_skew = 0.0;
    goodEstimate = false;
  }

  return goodEstimate;
}

//============================================================================================================================================================================







void ZFitter::fitTPCProton_SimulCent_ByRapidity(){
  #ifdef _ZFITTER_DEBUG_
    cout << " Call to ZFitter::fitTPCProton_SimulCent_ByRapidity()" << endl;
    cout << "    m_currentPartIndex:  " << m_currentPartIndex << endl;
    cout << "    m_currentPartCharge: " << m_currentPartCharge << endl;
    cout << "    m_currentRapidity:   " << m_currentRapidity << endl;
    cout << "    Mass:                " << m_partInfo->GetParticleMass(m_currentPartIndex) << endl;
  #endif

  TVirtualFitter::SetMaxIterations(1000000);

  bool   usePositiveChargeAsSeeds    = true;
  bool   usePreviousGoodParams       = true; 
  bool   fixAntiProtonValuesToProton = true;
  bool   kaonBackgroundMode          = (m_currentPartCharge == -1); // moves the low fit range so kaons are fit as a background instead of a whole peak
  double protonIsolatedStopMomentum  = 0.25;
  double kaonIsolatedLowMomForExtrap = 0.28;
  double sigmaGuess       = log(1.07); //0.067 ish for the ~7% resolution of the tpc

  //double pionHighMomParamExtrap = 0.55; // fits to this with lines
  //double kaonHighMomParamExtrap = 0.55; // fits to this with lines

  double momentumToKillprotons = 0.605;
  double pionSkewSeed = 0.5;
  double maximumMomentumToFit = 2.0; //matrix is used below instead
  double startOfEasyRegion    = 0.15;
  double stopOfEasyRegion     = 0.4;
  double momIsolatedProtonStop = 0.6;
  double maxMomForGausDataCalc = 0.25; // only isolated step
  double momentumToKillPions   = 0.8;
  double momToFixProtonSkew    = 0.92;

  //double pionSeedParams[4] = {-0.1155,0.07511,24.85,-14.67};
  //double pionSeedParams[4] = {-0.1155,0.513,24.85,-14.67};
  double pionSeedParams[4] = {-0.1155,0.513,0.4613,-0.2854};

  TH3D* histo3DPM = NULL;
  if(m_currentPartCharge == 1){ // Positive
    histo3DPM = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex];
  }else{ // Negative
    histo3DPM = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex];
  }
  
  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int histogramIndex = 0; histogramIndex < 1000; histogramIndex++){
      #ifndef _SAVE_POINTERS_
        if(m_rebinZTPCHistoStorage_ByCent[centIndex][histogramIndex]) delete m_rebinZTPCHistoStorage_ByCent[centIndex][histogramIndex];
      #endif
      m_rebinZTPCHistoStorage_ByCent[centIndex][histogramIndex] = NULL;
    }
  }

  int                   fitting_round = 0;
  vector<string>        roundTag;
  vector<double>        mTm0IndexVector;
  vector<int>           mTm0BinsFitPerRound;
  vector<unsigned int>  indexOfRoundChange;

  //##############   STEP 0  #############
  double step0_lowMtM0     = PhysMath::mTm0_from_y_p(m_currentRapidity, startOfEasyRegion, m_partInfo->GetParticleMass(m_currentPartIndex));
  double step0_highMtM0    = PhysMath::mTm0_from_y_p(m_currentRapidity, stopOfEasyRegion, m_partInfo->GetParticleMass(m_currentPartIndex));
  int    step0_lowBin      = histo3DPM->GetYaxis()->FindBin(step0_lowMtM0);
  if(step0_lowBin <= 0) step0_lowBin = 1;
  int    step0_highBin     = histo3DPM->GetYaxis()->FindBin(step0_highMtM0);
  if(step0_highBin <= 0) step0_highBin = 1;
  roundTag.push_back("calibRegion");
  mTm0BinsFitPerRound.push_back(0);
  cout << "Step 0 " << roundTag[0] << " Bin Ranges [" << step0_lowBin << "," << step0_highBin << "]" << endl; 
  if(step0_lowBin == step0_highBin && step0_highBin == 1){
    fitting_round++;
  }else{
    for(int bbb = step0_lowBin; bbb <= step0_highBin; bbb++) mTm0IndexVector.push_back(bbb);
  }
  indexOfRoundChange.push_back(mTm0IndexVector.size());


  //##############   STEP 1  #############
  double step1_highMtM0    = PhysMath::mTm0_from_y_p(m_currentRapidity, momIsolatedProtonStop, m_partInfo->GetParticleMass(m_currentPartIndex));
  int    step1_lowBin      = 1;
  int    step1_highBin     = histo3DPM->GetYaxis()->FindBin(step0_highMtM0);
  if(step1_highBin <= 0) step1_highBin = 1;
  roundTag.push_back("isolatedRegion");
  mTm0BinsFitPerRound.push_back(0);
  cout << "Step 1 " << roundTag[1] << " Bin Ranges [" << step1_lowBin << "," << step1_highBin << "]" << endl; 
  if(step1_lowBin == step1_highBin && step1_highBin == 1){
    fitting_round++;
  }else{
    for(int bbb = step1_lowBin; bbb <= step1_highBin; bbb++) mTm0IndexVector.push_back(bbb);
  }
  indexOfRoundChange.push_back(mTm0IndexVector.size());


  //##############   STEP 2  #############
  double step2_highMtM0    = PhysMath::mTm0_from_y_p(m_currentRapidity, momToFixProtonSkew, m_partInfo->GetParticleMass(m_currentPartIndex));
  double step2_lowBin   = step1_highBin - 2; //= 1;
  if(step2_lowBin <= 0) step2_lowBin = 1;
  int    step2_highBin     = histo3DPM->GetYaxis()->FindBin(step2_highMtM0);
  if(step2_highBin <= 0) step2_highBin = 1;
  roundTag.push_back("mixingRegion");
  mTm0BinsFitPerRound.push_back(0);
  cout << "Step 2 " << roundTag[2] << " Bin Ranges [" << step2_lowBin << "," << step2_highBin << "]" << endl; 
  if(step2_lowBin == step2_highBin && step2_highBin == 1){
    fitting_round++;
  }else{
    for(int fff = step2_lowBin; fff <= step2_highBin; fff++) mTm0IndexVector.push_back(fff);
  }
  indexOfRoundChange.push_back(mTm0IndexVector.size());

  //##############   STEP 3  #############
  double step3_lowBin   = step2_highBin - 1;
  if(step3_lowBin < 1)  step3_lowBin = 1;
  int    step3_highBin  = histo3DPM->GetNbinsY();
  for(int fff = step3_lowBin ; fff <= step3_highBin; fff++) mTm0IndexVector.push_back(fff);
  roundTag.push_back("overlapRegion");
  mTm0BinsFitPerRound.push_back(0);
  indexOfRoundChange.push_back(mTm0IndexVector.size());
  cout << "Step 3 " << roundTag[3] << " Bin Ranges [" << step3_lowBin << "," << step3_highBin << "]" << endl; 


  #ifdef _ZFITTER_DEBUG_
    cout << "------------------   MtM0 Index Index Array ------------------" << endl;
    for(unsigned int ggg = 0; ggg < mTm0IndexVector.size(); ggg++) cout << mTm0IndexVector[ggg] << "  ";
    cout << endl << "-----------------------------------------------------------------" << endl;
  #endif


  double pionMeanShift   = 0;
  double kaonMeanShift   = 0;
  double protonMeanShift = 0;


  double previousGoodParameters[15] = {-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999};
  for(int paramIndex = 0; paramIndex < 11; paramIndex++) previousGoodParameters[paramIndex] = -999;

  bool protonSkewSuccess  = false;
  double fixedProtonSkew = 0.0;

  
  double theStopMomentum = 50.0;
  for(int iii = 0; iii < m_colliderStopTable.size(); iii++){
    if( m_colliderStopTable[iii][0] == 2 && m_colliderStopTable[iii][1] == m_currentPartCharge && fabs(fabs(m_currentRapidity) - m_colliderStopTable[iii][2]) < 0.01){
      theStopMomentum = m_colliderStopTable[iii][6];
      break;
    }
  }

  #ifdef _ZFITTER_DEBUG_
    cout << "  Table Momentum Stop: " << theStopMomentum << endl;
  #endif





  for(unsigned int VertBinIndexIndex = 0; VertBinIndexIndex < mTm0IndexVector.size(); VertBinIndexIndex++){
    // What to do between trasitions of rounds of fitting
    if(VertBinIndexIndex == 0){
      cout << "########  STARTING ROUND " << fitting_round << " : " << roundTag[fitting_round] << endl;
    }

    m_currentMtM0Bin  = mTm0IndexVector[VertBinIndexIndex];
    m_currentMtM0     = histo3DPM->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
    m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
    if(!std::isfinite(m_currentMomOfBin) || m_currentMomOfBin <= 0.0) continue;
    if(m_currentMomOfBin > maximumMomentumToFit) continue;
    if(m_currentMomOfBin > theStopMomentum) continue;


    if(VertBinIndexIndex > 0 && mTm0IndexVector[VertBinIndexIndex] <= mTm0IndexVector[VertBinIndexIndex - 1]){
      #ifdef _ZFITTER_DEBUG_
        cout << "Fitting Round Was " << fitting_round << " Is Now " << fitting_round + 1 << endl;
      #endif

      if(roundTag[fitting_round] == "calibRegion"){
        calculateProtonMeanShifts(pionMeanShift,kaonMeanShift,protonMeanShift);
      }

      if(roundTag[fitting_round] == "isolatedRegion"){
        calculateProtonMeanShifts(pionMeanShift,kaonMeanShift,protonMeanShift);
      }

      if(roundTag[fitting_round] == "mixingRegion"){
        protonSkewSuccess = calculateProtonSkewExtrapolation(fixedProtonSkew);
        if(!protonSkewSuccess)   cout << "WARNING: Proton Skew Calculation Failed...   Skew = " << fixedProtonSkew << endl;
        else                     cout << "Proton Skew Calculation Successful...   Skew = " << fixedProtonSkew << endl;
      }

      //clear the previous good paramters because the bins aren't consecutive at round change
      for(int paramIndex = 0; paramIndex < 11; paramIndex++) previousGoodParameters[paramIndex] = -999;

      fitting_round++;
      cout << "########  STARTING ROUND " << fitting_round << " : " << roundTag[fitting_round] << endl;

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
      cout << "mTm0BinsFitPerRound[0] = " << mTm0BinsFitPerRound[0] << endl;
      cout << "mTm0BinsFitPerRound[1] = " << mTm0BinsFitPerRound[1] << endl;
      cout << "mTm0BinsFitPerRound[2] = " << mTm0BinsFitPerRound[2] << endl;

      cout << "protonSkewSuccess: " << (protonSkewSuccess ? "TRUE" : "FALSE") << endl;
      //cout << "protonRatioSuccess: " << (protonRatioSuccess ? "TRUE" : "FALSE") << endl;
      //cout << "protonMeanSuccess: " << (protonMeanSuccess ? "TRUE" : "FALSE") << endl;
      cout << " Cent 0 has " << HistogramUtilities::sumOfContent(m_currentHistosToFit_ByCent[0],-1,true)*m_centEvents[m_currentPartIndex]->GetBinContent(1) << " counts " << endl;
    #endif  
    if(countsOfHisto < 600){
      cout << "   Too few counts in Cent 0 histogram (" << countsOfHisto << ") ... Continuing... " << endl;
      continue;   
    }


    //################################    SETUP PION KAON proton PARAMETERS FOR FIT  ##################################

    double pionBichselPredict      =  log(m_bichselCurves[0]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + pionMeanShift;
    double kaonBichselPredict      =  log(m_bichselCurves[1]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + kaonMeanShift;
    double protonBichselPredict    =  log(m_bichselCurves[2]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + protonMeanShift;
    double pionMeanInit      = pionBichselPredict;
    double kaonMeanInit      = kaonBichselPredict;
    double protonMeanInit    = protonBichselPredict;
    double pionSigmaInit     = sigmaGuess;
    double kaonSigmaInit     = sigmaGuess;
    double protonSigmaInit   = sigmaGuess;
    double pionSkewInit      = pionSkewSeed;
    double kaonSkewInit      = pionSkewSeed;
    double protonSkewInit    = pionSkewSeed;






    if((roundTag[fitting_round] == "calibRegion" || roundTag[fitting_round] == "isolatedRegion") && m_currentMomOfBin < maxMomForGausDataCalc){
      double protonGausData[3];
      HistogramUtilities::getGausData(protonGausData, m_currentHistosToFit_ByCent[0],true,true, protonMeanInit - 6.0*sigmaGuess, protonMeanInit + 6.0*sigmaGuess);
      #ifdef _ZFITTER_DEBUG_
        cout << "Attempting Gaussian Parameter Calculation....  N: " << protonGausData[0] << "    mu: " << protonGausData[1] << "     sig: " << protonGausData[2] << endl;
      #endif
      if(std::isfinite(protonGausData[1]) && std::isfinite(protonGausData[2]) && protonGausData[2] < 0.4){
        protonMeanInit  = protonGausData[1];
        protonSigmaInit = protonGausData[2];
      }
    }

    /*else if(roundTag[fitting_round] != "protonRegion" && roundTag[fitting_round] != "protonRegionStdFix"){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Parameters from gaussian calculation... " << endl;
      #endif
      double pionGausData[3];
      HistogramUtilities::getGausData(pionGausData, m_currentHistosToFit_ByCent[0],true,false, pionMeanInit - 2.0*sigmaGuess, pionMeanInit + 2.0*sigmaGuess);
      if(fabs(pionGausData[1] - pionMeanInit) < 1.2*sigmaGuess) pionMeanInit = pionGausData[1];
      double kaonGausData[3];
      HistogramUtilities::getGausData(kaonGausData, m_currentHistosToFit_ByCent[0],true,false, kaonMeanInit - 2.0*sigmaGuess, kaonMeanInit + 2.0*sigmaGuess);
      if(fabs(kaonGausData[1] - kaonMeanInit) < 1.2*sigmaGuess) kaonMeanInit = kaonGausData[1];
      double protonGausData[3];
      HistogramUtilities::getGausData(protonGausData, m_currentHistosToFit_ByCent[0],true,false, protonMeanInit - 2.0*sigmaGuess, protonMeanInit + 2.0*sigmaGuess);
      if(fabs(protonGausData[1] - protonMeanInit) < 1.2*sigmaGuess) protonMeanInit = protonGausData[1];
    }*/




    //######################    LOAD PREVIOUS GOOD PARAMETERS AS SEEDS   ###########################
    if(usePreviousGoodParams && !currentBinWasFitBefore && isConsecutiveBin){ // use the previous fit parameters as seeds
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Parameters from Previous Good Parameters... " << endl;
      #endif 
      // Dont do in specific cent, the all cent would be better seeds
      if(previousGoodParameters[1] > 0 && std::isfinite(previousGoodParameters[0]) && std::isfinite(previousGoodParameters[1])){
        if(fabs(previousGoodParameters[0]) < 0.45) pionMeanInit        = previousGoodParameters[0] + log(m_bichselCurves[0]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        if(previousGoodParameters[1] < 0.2)        pionSigmaInit       = previousGoodParameters[1];
        if(fabs(previousGoodParameters[2]) > 0.95) pionSkewInit = 0.85;
        else                                       pionSkewInit = previousGoodParameters[2];
      }
      if(previousGoodParameters[4] > 0 && std::isfinite(previousGoodParameters[3]) && std::isfinite(previousGoodParameters[4])){
        if(roundTag[fitting_round] == "overlapRegion" && m_currentMomOfBin > momentumToKillPions){
          kaonMeanInit        = previousGoodParameters[3] + log(m_bichselCurves[1]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          kaonSigmaInit       = previousGoodParameters[4];
          if(fabs(previousGoodParameters[5]) > 0.95) kaonSkewInit = 0.85;
          else                                       kaonSkewInit = previousGoodParameters[5];
        }else{
          if(fabs(previousGoodParameters[3]) < 0.45) kaonMeanInit        = previousGoodParameters[3] + log(m_bichselCurves[1]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          if(previousGoodParameters[4] < 0.2)        kaonSigmaInit       = previousGoodParameters[4];  
          if(fabs(previousGoodParameters[5]) > 0.95) kaonSkewInit = 0.85;
          else                                       kaonSkewInit = previousGoodParameters[5];         
        }
      }
      if(previousGoodParameters[7] > 0 && std::isfinite(previousGoodParameters[6]) && std::isfinite(previousGoodParameters[7] && std::isfinite(previousGoodParameters[8]))){
        if(fabs(previousGoodParameters[6]) < 0.45) protonMeanInit    = previousGoodParameters[6] + log(m_bichselCurves[2]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        if(previousGoodParameters[7] < 0.2)        protonSigmaInit   = previousGoodParameters[7];
        if(fabs(previousGoodParameters[8]) > 0.95) protonSkewInit    = 0.85;
        else                                       protonSkewInit    = previousGoodParameters[8];
      }
    }
    /*if(roundTag[fitting_round] != "protonRegion" && roundTag[fitting_round] != "protonRegionStdFix" && protonMeanSuccess){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting proton Mean From Param Storage Histogram..." << endl;
      #endif
      protonMeanInit = m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
    }
    if(roundTag[fitting_round] != "protonRegion" && protonWidthSuccess){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting proton StdDev From Param Storage Histogram..." << endl;
      #endif
      protonSigmaInit = m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
    }*/






    if(usePositiveChargeAsSeeds && m_currentPartCharge < 0){
      #ifdef _ZFITTER_DEBUG_
        cout << "Getting Parameters from Positive Fits... " << endl;
      #endif

      /*
      double pionMeanTemp        = m_Fit_Data_ZTPC[m_currentPartIndex][0][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double pionSigmaTemp       = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double kaonMeanTemp        = m_Fit_Data_ZTPC[m_currentPartIndex][1][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double kaonSigmaTemp       = m_Fit_Data_ZTPC[m_currentPartIndex][1][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double protonMeanTemp      = m_Fit_Data_ZTPC[m_currentPartIndex][2][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double protonSigmaTemp     = m_Fit_Data_ZTPC[m_currentPartIndex][2][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double protonSkewTemp      = m_Fit_Data_ZTPC[m_currentPartIndex][2][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);  
      double pionMeanTemp_Err        = m_Fit_Data_ZTPC[m_currentPartIndex][0][1][0]->GetBinError(m_currentRapBin,m_currentMtM0Bin);
      double pionSigmaTemp_Err       = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][0]->GetBinError(m_currentRapBin,m_currentMtM0Bin);
      double kaonMeanTemp_Err        = m_Fit_Data_ZTPC[m_currentPartIndex][1][1][0]->GetBinError(m_currentRapBin,m_currentMtM0Bin);
      double kaonSigmaTemp_Err       = m_Fit_Data_ZTPC[m_currentPartIndex][1][2][0]->GetBinError(m_currentRapBin,m_currentMtM0Bin);
      double protonMeanTemp_Err      = m_Fit_Data_ZTPC[m_currentPartIndex][2][1][0]->GetBinError(m_currentRapBin,m_currentMtM0Bin);
      double protonSigmaTemp_Err     = m_Fit_Data_ZTPC[m_currentPartIndex][2][2][0]->GetBinError(m_currentRapBin,m_currentMtM0Bin);
      double protonSkewTemp_Err      = m_Fit_Data_ZTPC[m_currentPartIndex][2][3][0]->GetBinError(m_currentRapBin,m_currentMtM0Bin);       
      */
      
      double pionMeanTemp      = 0;       
      double pionSigmaTemp     = 0;
      double pionSkewTemp      = 0;      
      double kaonMeanTemp      = 0;       
      double kaonSigmaTemp     = 0; 
      double kaonSkewTemp      = 0;     
      double protonMeanTemp    = 0;     
      double protonSigmaTemp   = 0;    
      double protonSkewTemp    = 0;     
      double pionMeanTemp_Err  = 0;   
      double pionSigmaTemp_Err = 0;  
      double pionSkewTemp_Err = 0;
      double kaonMeanTemp_Err  = 0;   
      double kaonSigmaTemp_Err = 0;  
      double kaonSkewTemp_Err = 0; 
      double protonMeanTemp_Err  = 0; 
      double protonSigmaTemp_Err = 0;
      double protonSkewTemp_Err  = 0; 

      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][0][1][0],m_currentRapBin, m_currentMtM0, pionMeanTemp,   pionMeanTemp_Err);
      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][0][2][0],m_currentRapBin, m_currentMtM0, pionSigmaTemp,  pionSigmaTemp_Err);
      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][0][3][0],m_currentRapBin, m_currentMtM0, protonSkewTemp, pionSkewTemp_Err);

      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][1][1][0],m_currentRapBin, m_currentMtM0, kaonMeanTemp,   kaonMeanTemp_Err);
      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][1][2][0],m_currentRapBin, m_currentMtM0, kaonSigmaTemp,  kaonSigmaTemp_Err);
      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][1][3][0],m_currentRapBin, m_currentMtM0, protonSkewTemp, kaonSkewTemp_Err);

      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][2][1][0],m_currentRapBin, m_currentMtM0, protonMeanTemp,  protonMeanTemp_Err);
      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][2][2][0],m_currentRapBin, m_currentMtM0, protonSigmaTemp, protonSigmaTemp_Err);
      HistogramUtilities::interpolate2D(m_Fit_Data_ZTPC[m_currentPartIndex][2][3][0],m_currentRapBin, m_currentMtM0, protonSkewTemp,  protonSkewTemp_Err);


      bool pionWasFit = false;
      if(pionMeanTemp_Err > 0.0 || pionSigmaTemp_Err > 0.0) pionWasFit = true;
      bool kaonWasFit = false;
      if(kaonMeanTemp_Err > 0.0 || kaonSigmaTemp_Err > 0.0) kaonWasFit = true;
      bool protonWasFit = false;
      if(protonMeanTemp_Err > 0.0 || protonSigmaTemp_Err > 0.0 || protonSkewTemp_Err > 0.0) protonWasFit = true;      

      if(pionSigmaTemp > 0.0 && pionWasFit){
        pionMeanInit  = pionMeanTemp;
        pionSigmaInit = pionSigmaTemp;
        pionSkewInit  = pionSkewInit;
      }
      if(kaonSigmaTemp > 0.0 && kaonWasFit){
        kaonMeanInit  = kaonMeanTemp;
        kaonSigmaInit = kaonSigmaTemp;
        kaonSkewInit  = kaonSkewInit;
      }
      if(protonSigmaTemp > 0.0 && protonWasFit){
        protonMeanInit  = protonMeanTemp;
        protonSigmaInit = protonSigmaTemp;
        protonSkewInit  = protonSkewTemp;
        #ifdef _ZFITTER_DEBUG_
          cout << "Using PROTON+ parameters:" << endl;
          cout << "   P_Mean : " << protonMeanInit << endl;
          cout << "   P_Sig  : " << protonSigmaInit << endl;
          cout << "   P_Gam : "  << protonSkewInit << endl;
        #endif
      }


    }

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
      kaonSkewInit        = m_Fit_Data_ZTPC[m_currentPartIndex][1][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      protonMeanInit      = m_Fit_Data_ZTPC[m_currentPartIndex][2][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      protonSigmaInit     = m_Fit_Data_ZTPC[m_currentPartIndex][2][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      protonSkewInit      = m_Fit_Data_ZTPC[m_currentPartIndex][2][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
    }





    double lowRange  = pionBichselPredict - 3.5*sigmaGuess;
    double highRange = protonBichselPredict + 3.5*sigmaGuess;
    double lowRangeData  = m_currentHistosToFit_ByCent[0]->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistosToFit_ByCent[0]));
    double highRangeData = m_currentHistosToFit_ByCent[0]->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistosToFit_ByCent[0])+1);

    //#########################   DETERMINE WHICH PARTICLES AND PARAMETERS SHOULD BE FIT  ###################################
    bool doFitPion     = true;
    bool doFitKaon     = true;
    bool doFitProton   = true;
    
    bool fixProtonAmplitude = false;
    bool fixKaonAmplitude   = false;

    bool fixPionMean     = false;
    bool fixKaonMean     = false;
    bool fixProtonMean   = false;

    bool fixPionWidths     = false;
    bool fixKaonWidths     = false;
    bool fixProtonWidths   = false;

    bool fixPionSkew       = false;
    bool fixKaonSkew       = false;
    bool fixProtonSkew     = false;

    if( usePositiveChargeAsSeeds && m_currentPartCharge < 0 && fixAntiProtonValuesToProton ){
      fixProtonMean   = true;
      fixProtonWidths = true;
      fixProtonSkew   = true;
    }


    if(roundTag[fitting_round] == "calibRegion"){
      lowRange      = protonMeanInit - 5.0*sigmaGuess;
      highRange     = protonMeanInit + 5.0*sigmaGuess;
      doFitPion     = false;
      doFitKaon     = false;
    }

    if(roundTag[fitting_round] == "isolatedRegion"){
      lowRange      = protonMeanInit - 5.0*sigmaGuess;
      highRange     = protonMeanInit + 5.0*sigmaGuess;
      doFitPion     = false;
      doFitKaon     = false;
    }

    if(roundTag[fitting_round] == "mixingRegion"){
      highRange     = protonMeanInit + 5.0*sigmaGuess;
      if(m_currentMomOfBin < momIsolatedProtonStop){
        lowRange      = protonMeanInit - 5.0*sigmaGuess;
        highRange     = protonMeanInit + 5.0*sigmaGuess;
        doFitPion     = false;
        doFitKaon     = false;
      }else{
        if(kaonBackgroundMode){ // Treat kaon as background and only fit right hand side of peak
          if(m_currentMomOfBin < 0.77){
            lowRange      = protonMeanInit - 3.0*sigmaGuess;
            highRange     = protonMeanInit + 5.5*sigmaGuess;
            doFitPion     = false;
            doFitKaon     = false;
          }else{
            if(m_currentMomOfBin < momentumToKillPions){
              lowRange      = protonMeanInit - 2.0*sigmaGuess;
              highRange     = protonMeanInit + 5.5*sigmaGuess;
              doFitPion     = false;
              doFitKaon     = false;
            }else{
              doFitPion     = false;
              lowRange      = protonMeanInit - 5.0*sigmaGuess;
             // double tempLowRange  = kaonBichselPredict + protonMeanShift + 1.5*sigmaGuess - sigmaGuess*(m_currentMomOfBin - 0.8)/0.15;
              double tempLowRange  = protonMeanInit - 5.0*sigmaGuess - sigmaGuess*(m_currentMomOfBin - 0.8)/0.15;
              //if(tempLowRange < pionBichselPredict - 0.5*sigmaGuess) tempLowRange = pionBichselPredict - 0.5*sigmaGuess;
              if(tempLowRange < lowRange) lowRange = tempLowRange;
            }
          }
        }else{
          lowRange      = pionBichselPredict -2.0*sigmaGuess;
          highRange     = protonMeanInit + 3.0*sigmaGuess;
        }
      }

      /*if(m_currentMomOfBin > momentumToKillPions){
        doFitPion     = false;
        lowRange      = protonMeanInit - 5.0*sigmaGuess;
        double tempLowRange  = kaonBichselPredict + protonMeanShift + 1.5*sigmaGuess;
        if(tempLowRange < lowRange) lowRange = tempLowRange;
      }else{ // region to fit pion and kaon
        double halfWay = 0.5*(pionBichselPredict - 3.5*sigmaGuess + protonMeanInit);
        int maxBin = HistogramUtilities::getMaxBinInRange(m_currentHistosToFit_ByCent[0], pionBichselPredict - 3.5*sigmaGuess, halfWay);
        lowRange = m_currentHistosToFit_ByCent[0]->GetBinCenter(maxBin) + sigmaGuess;
        #ifdef _ZFITTER_DEBUG_
          cout << "Re-adjusting limits in mixing region...   halfway=" << halfWay << "  maxBin=" << maxBin << "   lowRange = " << lowRange << endl;
        #endif
      }*/
          
      if(mTm0BinsFitPerRound[0] < 2 && mTm0BinsFitPerRound[1] < 2){
        protonSkewInit = 0.0;
        fixProtonSkew = true;
      }
    }

    if(roundTag[fitting_round] == "overlapRegion"){
      if(kaonBackgroundMode){ // Treat kaon as background and only fit right hand side of peak
        if(m_currentMomOfBin > momentumToKillPions){
          doFitPion     = false;
          //lowRange      = protonMeanInit - 5.0*sigmaGuess;
          //double tempLowRange_0  = kaonMeanInit + sigmaGuess;
          //double tempLowRange_1  = kaonBichselPredict + protonMeanShift + sigmaGuess;
          //if(tempLowRange_0 > lowRange) lowRange = tempLowRange_0;
          //if(tempLowRange_1 > lowRange) lowRange = tempLowRange_1;
          //if(lowRange > protonMeanInit - 3.0*sigmaGuess) lowRange = protonMeanInit - 3.0*sigmaGuess;
          //if(lowRange < protonMeanInit - 8.0*sigmaGuess) lowRange = protonMeanInit - 8.0*sigmaGuess;
          lowRange      = protonMeanInit - 5.0*sigmaGuess - sigmaGuess*(m_currentMomOfBin - 0.8)/0.15;
        }
      }else{
          lowRange      = pionBichselPredict - 2.0*sigmaGuess;
          highRange     = protonMeanInit + 3.0*sigmaGuess;
      }

      if(m_currentMomOfBin > momToFixProtonSkew && protonSkewSuccess){
        protonSkewInit = fixedProtonSkew;
        fixProtonSkew = true;
      }
    }



    //Set Ranges If Off Histogram Edges
    if(lowRange < lowRangeData)   lowRange  = lowRangeData;
    if(highRange > highRangeData) highRange = highRangeData;
    m_currentLowFitRange  = lowRange;
    m_currentHighFitRange = highRange;
    #ifdef _ZFITTER_DEBUG_
      cout << " LowFitRange: " << m_currentLowFitRange << "  highFitRange: " << m_currentHighFitRange << endl;
    #endif
    for(int centIndex = 0;centIndex < m_numCentralities; centIndex++){
      int binLow  =  m_currentHistosToFit_ByCent[centIndex]->FindBin(m_currentLowFitRange);
      if(binLow > 1 && binLow < m_currentHistosToFit_ByCent[centIndex]->GetNbinsX()) m_currentLowFitBin[centIndex] = binLow;
      else  m_currentLowFitBin[centIndex] = 1;
      int binHigh =  m_currentHistosToFit_ByCent[centIndex]->FindBin(m_currentHighFitRange);
      if(binHigh > 1  && binHigh  <  m_currentHistosToFit_ByCent[centIndex]->GetNbinsX()) m_currentHighFitBin[centIndex] = binHigh;
      else m_currentHighFitBin[centIndex] = m_currentHistosToFit_ByCent[centIndex]->GetNbinsX();
      #ifdef _ZFITTER_DEBUG_
        cout << "    Cent: " << centIndex << "  bin " << m_currentLowFitBin[centIndex] << " to " << m_currentHighFitBin[centIndex] << endl;
      #endif
    }


    if(m_currentHistosToFit_ByCent[0]->FindBin(m_currentLowFitRange) == m_currentHistosToFit_ByCent[0]->FindBin(m_currentHighFitRange)){
      cout << "WARNING: Range is too small Z_TPC = [" << m_currentLowFitRange << " , " << m_currentHighFitRange << "]" 
           << "  with bins = [" << m_currentHistosToFit_ByCent[0]->FindBin(m_currentLowFitRange) << " , " << m_currentHistosToFit_ByCent[0]->FindBin(m_currentHighFitRange)
           << "]   Continuing to next bin..." << endl;
      continue;
    }

    #ifdef _ZFITTER_DEBUG_
      cout << "####### This Fit will Have The Following Toggles: " << endl;
      cout << "doFitPion:            " << (doFitPion ? "TRUE" : "FALSE") << endl;
      cout << "doFitKaon:            " << (doFitKaon ? "TRUE" : "FALSE") << endl;
      cout << "doFitProton:          " << (doFitProton ? "TRUE" : "FALSE") << endl;
      cout << "fixKaonAmplitude:     " << (fixKaonAmplitude ? "TRUE" : "FALSE") << endl;
      cout << "fixProtonAmplitude:   " << (fixProtonAmplitude ? "TRUE" : "FALSE") << endl;
      cout << "fixPionMean:          " << (fixPionMean ? "TRUE" : "FALSE") << endl;
      cout << "fixKaonMean:          " << (fixKaonMean ? "TRUE" : "FALSE") << endl;
      cout << "fixProtonMean:        " << (fixProtonMean ? "TRUE" : "FALSE") << endl;
      cout << "fixPionWidths:        " << (fixPionWidths ? "TRUE" : "FALSE") << endl;
      cout << "fixKaonWidths:        " << (fixKaonWidths ? "TRUE" : "FALSE") << endl;
      cout << "fixProtonWidths:      " << (fixProtonWidths ? "TRUE" : "FALSE") << endl;
      cout << "fixPionSkew:          " << (fixPionSkew ? "TRUE" : "FALSE") << endl;
      cout << "fixKaonSkew:          " << (fixKaonSkew ? "TRUE" : "FALSE") << endl;
      cout << "fixProtonSkew:        " << (fixProtonSkew ? "TRUE" : "FALSE") << endl;
    #endif



    //################  SETUP THE FIT FUNCTION ###################
    bool doKaonPionRatio = false;
    if(roundTag[fitting_round] != "isolatedRegion") doKaonPionRatio = false; //change to true? --------------------------------------------------------------------------------------------------------------------------------------------------------
    if(!doFitKaon) doKaonPionRatio = false;

    //    mu_pi sigma_pi  mu_k sigma_k  mu_p sigma_p gamma_p then For Each Cent: A_pi A_K A_p ... A_pi A_K A_p 
    double parameters_simul_cent[75];
    for(int iii = 0; iii < 75; iii++) parameters_simul_cent[iii] = 0;
    parameters_simul_cent[0] = pionMeanInit;
    parameters_simul_cent[1] = pionSigmaInit;
    parameters_simul_cent[2] = pionSkewInit;
    parameters_simul_cent[3] = kaonMeanInit;
    parameters_simul_cent[4] = kaonSigmaInit;
    parameters_simul_cent[5] = kaonSkewInit;
    parameters_simul_cent[6] = protonMeanInit;
    parameters_simul_cent[7] = protonSigmaInit;
    parameters_simul_cent[8] = protonSkewInit;

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      if(currentBinWasFitBefore){
        parameters_simul_cent[9  + centIndex*3]  = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        parameters_simul_cent[10 + centIndex*3]  = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        parameters_simul_cent[11 + centIndex*3]  = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][2][centIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{
        parameters_simul_cent[9 + centIndex*3] = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(pionMeanInit)) / PhysMath::skewNormal(pionMeanInit,pionSigmaInit,pionSkewInit,pionMeanInit);
        if(roundTag[fitting_round] == "overlapRegion" && m_currentMomOfBin > momentumToKillPions){
          parameters_simul_cent[10 + centIndex*3] = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(m_currentLowFitRange)) / PhysMath::skewNormal(kaonMeanInit,kaonSigmaInit,kaonSkewInit,m_currentLowFitRange);          
        }else{
          parameters_simul_cent[10 + centIndex*3] = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(kaonMeanInit)) / PhysMath::skewNormal(kaonMeanInit,kaonSigmaInit,kaonSkewInit,kaonMeanInit);                    
        }
        parameters_simul_cent[11 + centIndex*3] = m_currentHistosToFit_ByCent[centIndex]->GetBinContent(m_currentHistosToFit_ByCent[centIndex]->FindBin(protonMeanInit)) / PhysMath::skewNormal(protonMeanInit,protonSigmaInit,protonSkewInit,protonMeanInit);         
      }
      /*
      if(  roundTag[fitting_round] != "protonRegion" && roundTag[fitting_round] != "protonRegionStdFix" && protonRatioSuccess){
        if(doprotonPionRatio){
          parameters_simul_cent[9 + centIndex*3]  = m_protonRatioFitHistos_ByCent[centIndex][1]->GetBinContent(m_currentMtM0Bin);
          if(parameters_simul_cent[9 + centIndex*3] <= 0.0){
            parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[7 + centIndex*3];
          }
        }else{
          parameters_simul_cent[9 + centIndex*3]  = m_protonRatioFitHistos_ByCent[centIndex][2]->GetBinContent(m_currentMtM0Bin);
          if(parameters_simul_cent[9 + centIndex*3] <= 0.0){
            parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[8 + centIndex*3];
          }
        }       
      }else{
        if(doKaonPionRatio){
          parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[7 + centIndex*3];
        }else{
          parameters_simul_cent[9 + centIndex*3]  = parameters_simul_cent[9 + centIndex*3]/parameters_simul_cent[8 + centIndex*3];
        } 
      }
      if(!std::isfinite(parameters_simul_cent[9 + centIndex*3])){
        parameters_simul_cent[9 + centIndex*3] = 0.0;
      }*/
    }  


    bool pionAmpsAllZero     = true;
    bool kaonAmpsAllZero     = true;
    bool protonAmpsAllZero   = true;
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      if(parameters_simul_cent[9  + centIndex*3] > 0.0) pionAmpsAllZero     = false;
      if(parameters_simul_cent[10 + centIndex*3] > 0.0) kaonAmpsAllZero     = false;
      if(parameters_simul_cent[11 + centIndex*3] > 0.0) protonAmpsAllZero   = false;
    }  

    if(pionAmpsAllZero)     doFitPion     = false;
    if(kaonAmpsAllZero)     doFitKaon     = false;
    if(protonAmpsAllZero)   doFitProton   = false;
    
    #ifdef _ZFITTER_DEBUG_
      cout << " Pi_Mean: " << parameters_simul_cent[0] << endl;
      cout << " Pi_Sig : " << parameters_simul_cent[1] << endl;
      cout << " Pi_Gam : " << parameters_simul_cent[2] << endl;
      cout << " K_Mean : " << parameters_simul_cent[3] << endl;
      cout << " K_Sig  : " << parameters_simul_cent[4] << endl;
      cout << " K_Gam  : " << parameters_simul_cent[5] << endl;
      cout << " P_Mean : " << parameters_simul_cent[6] << endl;
      cout << " P_Sig  : " << parameters_simul_cent[7] << endl;
      cout << " P_Gam  : " << parameters_simul_cent[8] << endl;
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        cout << " N_Pi_" << centIndex << " :" << parameters_simul_cent[9  + centIndex*3] << endl;
        cout << " N_K_" << centIndex << "  :" << parameters_simul_cent[10 + centIndex*3] << endl;
        cout << " N_P_" << centIndex << "  :" << parameters_simul_cent[11 + centIndex*3] << endl;
      }  
    #endif


    bool initialParametersOkay = true;
    for(int iii = 0; iii < 11 + 3*m_numCentralities; iii++){
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
    //if(m_currentMomOfBin < momToLinkPiKSigma){
      //chiSqrdFunctor = new ROOT::Math::Functor(this, &ZFitter::simultaneous_centrality_ZTPC_skewNormalAndTwoGaus_Proton_chisqr,(m_numCentralities*3) + 11);
      chiSqrdFunctor = new ROOT::Math::Functor(this, &ZFitter::simultaneous_centrality_ZTPC_ThreeSkewNormal_Proton_chisqr,(m_numCentralities*3) + 9);

    //}else{
    //  linkingPiKSigma = true;
    //  chiSqrdFunctor = new ROOT::Math::Functor(this, &ZFitter::simultaneous_centrality_ZTPC_skewNormalAndTwoGaus_sigmaLinked_chisqr,(m_numCentralities*3) + 7);
    //}
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
    minimizer->SetVariable(0, "#mu_{#pi}" ,pionMeanInit, 0.001);
    minimizer->SetVariableLimits(0, pionMeanInit - 0.2, pionMeanInit + 0.15);
    minimizer->SetVariable(1, "#sigma_{#pi}" ,pionSigmaInit, 0.001);
    if(m_currentPartCharge == 1){
      minimizer->SetVariableLimits(1, 0.04, 0.15);
    }else{
      minimizer->SetVariableLimits(1, 0.04, 0.25);
    }
    minimizer->SetVariable(2, "#gamma_{#pi}" ,pionSkewInit, 0.01);
    minimizer->SetVariableLimits(2, -0.01, 0.9);

    minimizer->SetVariable(3, "#mu_{K}" ,kaonMeanInit, 0.001);
    double midwayKaonPion = 0.1*protonBichselPredict + 0.9*kaonBichselPredict;
    minimizer->SetVariableLimits(3, kaonMeanInit - 0.15, midwayKaonPion);
    if(kaonBackgroundMode){
      if(roundTag[fitting_round] == "overlapRegion"){
        if(m_currentMomOfBin > momentumToKillPions){
          minimizer->SetVariableLimits(3, pionBichselPredict - sigmaGuess, midwayKaonPion);
        }else{
          minimizer->SetVariableLimits(3, kaonMeanInit - 0.15, midwayKaonPion);
        }
      }
    }else{
      minimizer->SetVariableLimits(3, kaonMeanInit - 0.2, midwayKaonPion);
    }

    minimizer->SetVariable(4, "#sigma_{K}" ,kaonSigmaInit, 0.01);
    minimizer->SetVariableLimits(4, 0.04, 0.09);
    if(kaonBackgroundMode && roundTag[fitting_round] == "overlapRegion" && m_currentMomOfBin > momentumToKillPions){
      minimizer->SetVariableLimits(4, 0.04, 0.15);
    }
    minimizer->SetVariable(5, "#gamma_{K}" ,kaonSkewInit, 0.01);
    minimizer->SetVariableLimits(5, -0.01, 0.85);

    minimizer->SetVariable(6, "#mu_{p}" ,protonMeanInit, 0.001);
    if(roundTag[fitting_round] == "calibRegion"){
      minimizer->SetVariableLimits(6, protonMeanInit - 0.75, protonMeanInit + 0.75);
    }else if(roundTag[fitting_round] == "isolatedRegion"){
      minimizer->SetVariableLimits(6, protonMeanInit - 0.35, protonMeanInit + 0.35);
    }else{
      minimizer->SetVariableLimits(6, protonMeanInit - 0.15, protonMeanInit + 0.15);
    }
    minimizer->SetVariable(7, "#sigma_{p}" ,protonSigmaInit, 0.01);
    minimizer->SetVariableLimits(7, 0.04, 0.13);
    minimizer->SetVariable(8, "#gamma_{p}" , protonSkewInit , 0.01);
    minimizer->SetVariableLimits(8, -0.8,0.9);

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      minimizer->SetVariable(9 + centIndex*3, Form("N_{#pi}_%d",centIndex), parameters_simul_cent[9 + centIndex*3], 0.01*parameters_simul_cent[9 + centIndex*3]);
      minimizer->SetVariableLimits(9 + centIndex*3, 0.0001*parameters_simul_cent[9 + centIndex*3], 100.0*parameters_simul_cent[9 + centIndex*3]);
      if(parameters_simul_cent[9 + centIndex*3] <= 0.0){
        minimizer->SetVariable(9 + centIndex*3, Form("N_{#pi}_%d",centIndex), 0.0, 0.0001);
        minimizer->SetVariableLimits(9 + centIndex*3, 0.0, 1.0);
      }
      minimizer->SetVariable(10 + centIndex*3, Form("N_{K}_%d",centIndex), parameters_simul_cent[10 + centIndex*3], 0.01*parameters_simul_cent[10 + centIndex*3]);
      minimizer->SetVariableLimits(10 + centIndex*3, 0.0001*parameters_simul_cent[10 + centIndex*3], 100.0*parameters_simul_cent[10 + centIndex*3]);
      if(parameters_simul_cent[10 + centIndex*3] <= 0.0){
        minimizer->SetVariable(10 + centIndex*3, Form("N_{K}_%d",centIndex), 0.0, 0.0001);
        minimizer->SetVariableLimits(10 + centIndex*3, 0.0, 1.0);        
      }
      //double protonPionRatioInit = m_protonRatioFitHistos_ByCent[centIndex][1]->GetBinContent(m_currentMtM0Bin);
      //double protonKaonRatioInit = m_protonRatioFitHistos_ByCent[centIndex][2]->GetBinContent(m_currentMtM0Bin);

      minimizer->SetVariable(11 + centIndex*3, Form("N_{p}_%d",centIndex), parameters_simul_cent[11 + centIndex*3], 0.01*parameters_simul_cent[11 + centIndex*3]);
      minimizer->SetVariableLimits(11 + centIndex*3, 0.0001*parameters_simul_cent[11 + centIndex*3], 100.0*parameters_simul_cent[11 + centIndex*3]);
      if(parameters_simul_cent[11 + centIndex*3] <= 0.0){
        minimizer->SetVariable(11 + centIndex*3, Form("N_{p}_%d",centIndex), 0.0, 0.0001);
        minimizer->SetVariableLimits(11 + centIndex*3, 0.0, 1.0);        
      }
      
    }  


    //-----------------------      FIXING PARAMETERS     ------------------------------
    //PION
    if(fixPionMean)       minimizer->FixVariable(0);
    if(fixPionWidths)     minimizer->FixVariable(1);
    if(fixPionSkew)       minimizer->FixVariable(2);
    if(!doFitPion){
      minimizer->SetFixedVariable(0, "#mu_{#pi}" ,0.0);
      minimizer->SetFixedVariable(1, "#sigma_{#pi}" ,1.0);
      minimizer->SetFixedVariable(2, "#gamma_{#pi}" ,0.0);
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->SetFixedVariable(9 + centIndex*3, Form("N_{#pi}_%d",centIndex), 0.0);
      }  
    }

    //KAON
    if(fixKaonAmplitude){
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->FixVariable(10 + centIndex*3);
      }
    }
    if(fixKaonMean)       minimizer->FixVariable(3);
    if(fixKaonWidths)     minimizer->FixVariable(4);
    if(fixKaonSkew)       minimizer->FixVariable(5);
    if(!doFitKaon){
      minimizer->SetFixedVariable(3, "#mu_{K}",    0.0);
      minimizer->SetFixedVariable(4, "#sigma_{K}", 1.0);
      minimizer->SetFixedVariable(5, "#gamma_{K}", 0.0);
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->SetFixedVariable(10 + centIndex*3, Form("N_{K}_%d",centIndex), 0.0);
      }  
    }

    //PROTON
    if(fixProtonAmplitude){   
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->FixVariable(11 + centIndex*3);
      }
    }      
    if(fixProtonMean)   minimizer->FixVariable(6);
    if(fixProtonWidths) minimizer->FixVariable(7);
    if(fixProtonSkew)   minimizer->FixVariable(8);   
    if(!doFitProton){
      minimizer->SetFixedVariable(6,"#mu_{p}",   0.0);
      minimizer->SetFixedVariable(7,"#sigma_{p}",1.0);
      minimizer->SetFixedVariable(8,"#gamma_{p}" , 0.1);
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        minimizer->SetFixedVariable(11 + centIndex*3, Form("N_{p}_%d",centIndex), 0.0); 
      }  
    }


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
    double simulParamLimits[55][2];
    for(int iii = 0; iii < 9 + 3*m_numCentralities; iii++){
      simulParams[iii]      = finalParams[iii];
      simulParamErrors[iii] = sqrt(minimizer->CovMatrix(iii,iii));
      //const ROOT::Math::IMultiGenFunction::VariableLimits& limits = minimizer->GetVariableLimits(iii);
      //simulParamLimits[iii][0] = limits.Min();
      //simulParamLimits[iii][1] = limits.Max();

      ROOT::Fit::ParameterSettings* setting = new ROOT::Fit::ParameterSettings();
      minimizer->GetVariableSettings(iii, *setting);
      simulParamLimits[iii][0] = setting->LowerLimit();
      simulParamLimits[iii][1] = setting->UpperLimit();

      //simulParamLimits[iii][0] = minimizer->GetVariableLowerLimit();
      //simulParamLimits[iii][1] = minimizer->GetVariableUpperLimit();

    }

    /*
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      simulParams[9+3*centIndex]      = finalParams[9+3*centIndex];
      simulParamErrors[9+3*centIndex] = sqrt(minimizer->CovMatrix(9+3*centIndex, 9+3*centIndex));
      simulParams[10+3*centIndex]      = finalParams[10+3*centIndex];
      simulParamErrors[10+3*centIndex] = sqrt(minimizer->CovMatrix(10+3*centIndex, 10+3*centIndex));
      simulParams[11+3*centIndex]      = finalParams[11+3*centIndex];
      simulParamErrors[11+3*centIndex] = sqrt(minimizer->CovMatrix(11+3*centIndex, 11+3*centIndex));
      / *if(doprotonPionRatio){
        simulParams[9+3*centIndex]      = finalParams[9+3*centIndex] * finalParams[7+3*centIndex];
        simulParamErrors[9+3*centIndex] = simulParams[9+3*centIndex] * sqrt(minimizer->CovMatrix(7+3*centIndex, 7+3*centIndex)/(finalParams[7+3*centIndex] * finalParams[7+3*centIndex])
                                              + minimizer->CovMatrix(9+3*centIndex, 9+3*centIndex)/(finalParams[9+3*centIndex] * finalParams[9+3*centIndex]));
      }else{
        simulParams[9+3*centIndex]      = finalParams[9+3*centIndex] * finalParams[8+3*centIndex];
        simulParamErrors[9+3*centIndex] = simulParams[9+3*centIndex] * sqrt(minimizer->CovMatrix(8+3*centIndex, 8+3*centIndex)/(finalParams[8+3*centIndex] * finalParams[8+3*centIndex])
                                              + minimizer->CovMatrix(9+3*centIndex, 9+3*centIndex)/(finalParams[9+3*centIndex] * finalParams[9+3*centIndex]));
      }* /
    }*/
    
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
        previousGoodParameters[5] = simulParams[5]; 
      }else{
        previousGoodParameters[3] = -999;
        previousGoodParameters[4] = -999;
        previousGoodParameters[5] = -999;
      }
      if(doFitProton){
        previousGoodParameters[6] = simulParams[6] - log(m_bichselCurves[2]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        previousGoodParameters[7] = simulParams[7];
        previousGoodParameters[8] = simulParams[8];    
      }else{
        previousGoodParameters[6] = -999;
        previousGoodParameters[7] = -999;
        previousGoodParameters[8] = -999;
      }       
    }else{
      for(int iii = 0; iii < 11; iii++) previousGoodParameters[iii] = -999;
    }

    TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
    fittingCanvas->SetWindowSize(1300,700);
    fittingCanvas->SetCanvasSize(1200,600);
    fittingCanvas->Draw();
    TPad* topPad = new TPad("topPad_P_TPC","topPad_P_TPC",0.0,0.3,1.0,1.0);
    TPad* bottomPad = new TPad("bottomPad_P_TPC","bottomPad_P_TPC",0.0,0.0,1.0,0.3);
    topPad->Draw();
    bottomPad->Draw();
    fittingCanvas->cd();
    topPad->cd();
    topPad->SetLogy(true);
    topPad->SetRightMargin(0.4);
    bottomPad->SetLogy(false);
    bottomPad->SetRightMargin(0.4);
    gStyle->SetOptFit(0111);

    TF1* combinedFunct_CentStorage[16];
    for(int iii = 0; iii < 16; iii++) combinedFunct_CentStorage[iii] = NULL;

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      double numberOfCentEvents = m_centEvents[m_currentPartIndex]->GetBinContent(centIndex+1);
      m_currentHistosToFit_ByCent[centIndex]->SetMarkerStyle(20);
      m_currentHistosToFit_ByCent[centIndex]->SetMarkerColor(kBlack);
      m_currentHistosToFit_ByCent[centIndex]->SetMarkerSize(0.5);
      m_currentHistosToFit_ByCent[centIndex]->SetLineWidth(2);
      m_currentHistosToFit_ByCent[centIndex]->SetLineColor(kBlack);


      #ifndef _ZFITTER_PHYSMATH_BYPASS_
        //TF1* combinedFunct = PhysMath::getSkewAndTwoGausFunct(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d",fitting_round,
        //  m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin),
        //  simulParams[9 + centIndex*3], simulParams[4],simulParams[5],simulParams[6], 
        //  simulParams[8 + centIndex*3], simulParams[2], simulParams[3], 
        //  simulParams[7 + centIndex*3], simulParams[0], simulParams[1]);
        //TF1* kaonFunct  = PhysMath::getGaussianFunction(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d_Kaon",fitting_round,
        //  m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin), simulParams[8 + centIndex*3], simulParams[2], simulParams[3]);
        //TF1* pionFunct = PhysMath::getGaussianFunction(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d_proton",fitting_round,
        //  m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin), simulParams[7 + centIndex*3], simulParams[0], simulParams[1]);
        TF1* combinedFunct = PhysMath::getThreeSkewFunct(Form("ztpc_step%d_fit_%d_%d_%02d_%03d,%04d",fitting_round,
                 m_currentPartIndex, m_currentPlusMinusIndex,centIndex,m_currentRapBin,m_currentMtM0Bin), 
                 simulParams[9  + centIndex*3], simulParams[0], simulParams[1], simulParams[2],
                 simulParams[10 + centIndex*3], simulParams[3], simulParams[4], simulParams[5],
                 simulParams[11 + centIndex*3], simulParams[6], simulParams[7], simulParams[8]);
        TF1* pionFunct    = PhysMath::skewNormalFunct( simulParams[9 + centIndex*3], simulParams[0], simulParams[1], simulParams[2]);
        TF1* kaonFunct    = PhysMath::skewNormalFunct( simulParams[10 + centIndex*3], simulParams[3], simulParams[4], simulParams[5]);
        TF1* protonFunct  = PhysMath::skewNormalFunct( simulParams[11 + centIndex*3], simulParams[6], simulParams[7], simulParams[8]);
        combinedFunct_CentStorage[centIndex] = combinedFunct;
        pionFunct->SetLineColor(m_partInfo->GetParticleColor(0));
        pionFunct->SetRange(lowRange,highRange);
        kaonFunct->SetLineColor(m_partInfo->GetParticleColor(1));
        kaonFunct->SetRange(lowRange,highRange);
        protonFunct->SetLineColor(m_partInfo->GetParticleColor(2));
        protonFunct->SetRange(lowRange,highRange);
        TF1* pionFunct_init = NULL;
        TF1* kaonFunct_init = NULL;
        TF1* protonFunct_init = NULL;
        if(m_drawInitialSeedsToFits){
          protonFunct_init    = PhysMath::skewNormalFunct( parameters_simul_cent[9 + centIndex*3], parameters_simul_cent[0], parameters_simul_cent[1], parameters_simul_cent[2]);
          kaonFunct_init      = PhysMath::skewNormalFunct( parameters_simul_cent[10 + centIndex*3], parameters_simul_cent[3], parameters_simul_cent[4], parameters_simul_cent[5]);
          pionFunct_init      = PhysMath::skewNormalFunct( parameters_simul_cent[11 + centIndex*3], parameters_simul_cent[6], parameters_simul_cent[7], parameters_simul_cent[8]);
          pionFunct_init->SetLineColor(kYellow+2);
          pionFunct_init->SetRange(lowRange,highRange);
          kaonFunct_init->SetLineColor(kYellow+1);
          kaonFunct_init->SetRange(lowRange,highRange);
          protonFunct_init->SetLineColor(kOrange);
          protonFunct_init->SetRange(lowRange,highRange);
          if(!doFitPion)   pionFunct_init->SetParameter(0,0.0);
          if(!doFitKaon)   kaonFunct_init->SetParameter(0,0.0);
          if(!doFitProton) protonFunct_init->SetParameter(0,0.0);

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
        if(pionFunct_init) pionFunct_init->Draw("same");
        if(kaonFunct_init) kaonFunct_init->Draw("same");
        if(protonFunct_init) protonFunct_init->Draw("same");
        combinedFunct->Draw("same");
        pionFunct->Draw("same");
        protonFunct->Draw("same");
        kaonFunct->Draw("same");
      #endif
      lowRangeLine->Draw("same");
      highRangeLine->Draw("same");




      TPaveText* fitParsTxt = new TPaveText(0.65,0.05,0.95,0.95,"NDC");
      fitParsTxt->AddText(Form("N_{#pi} \t %e#pm%e %s",simulParams[9 + centIndex*3],simulParamErrors[9 + centIndex*3],minimizer->IsFixedVariable(9 + 3*centIndex) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#mu_{#pi} \t %f#pm%f %s",simulParams[0],simulParamErrors[0],minimizer->IsFixedVariable(0) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#sigma_{#pi} \t %f#pm%f %s",simulParams[1],simulParamErrors[1],minimizer->IsFixedVariable(1) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#gamma_{#pi} \t %f#pm%f %s",simulParams[2],simulParamErrors[2],minimizer->IsFixedVariable(2) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("N_{K} \t %e#pm%e %s",simulParams[10 + centIndex*3], simulParamErrors[10 + centIndex*3],minimizer->IsFixedVariable(10 + 3*centIndex) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#mu_{K} \t %f#pm%f %s",simulParams[3],simulParamErrors[3],minimizer->IsFixedVariable(3) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#sigma_{K} \t %f#pm%f %s",simulParams[4],simulParamErrors[4],minimizer->IsFixedVariable(4) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#gamma_{K} \t %f#pm%f %s",simulParams[5],simulParamErrors[5],minimizer->IsFixedVariable(5) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("N_{p} \t %e#pm%e %s",simulParams[11 + centIndex*3], simulParamErrors[11 + centIndex*3],minimizer->IsFixedVariable(11 + centIndex*3) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#mu_{p} \t %f#pm%f %s",simulParams[6],simulParamErrors[6],minimizer->IsFixedVariable(6) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#sigma_{p} \t %f#pm%f %s",simulParams[7],simulParamErrors[7],minimizer->IsFixedVariable(7) ? "(fixed)" : "       "));
      fitParsTxt->AddText(Form("#gamma_{p} \t %f#pm%f %s",simulParams[8],simulParamErrors[8],minimizer->IsFixedVariable(8) ? "(fixed)" : "       "));
      fitParsTxt->SetFillStyle(0);
      fitParsTxt->SetBorderSize(5);
      fitParsTxt->Draw("same");





      bottomPad->cd();
      TH1D* pullHisto = nullptr;
      #ifndef _ZFITTER_PHYSMATH_BYPASS_
      pullHisto = HistogramUtilities::makeFitPullHistogram(m_currentHistosToFit_ByCent[centIndex], combinedFunct, true);
      pullHisto->SetMarkerStyle(20);
      pullHisto->SetMarkerSize(0.5);
      pullHisto->SetMarkerColor(kBlack);
      pullHisto->SetLineColor(kBlack);
      pullHisto->GetXaxis()->SetTitle(Form("Z_{TPC}"));
      pullHisto->GetYaxis()->SetTitle(Form("Standardized Pull"));
      pullHisto->GetXaxis()->SetRangeUser(m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( lowRange - 0.5*sigmaGuess)), m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( highRange + 0.5*sigmaGuess)+1));
      pullHisto->GetYaxis()->SetRangeUser(-15,15);
      #endif

      TPaveText* fitTxt = new TPaveText(0.65,0.03,0.95,0.97,"NDC");
      fitTxt->AddText(Form("N_{#pi} \t %e  [%e,%e]",     parameters_simul_cent[9 + centIndex*3],simulParamLimits[9 + centIndex*3][0],simulParamLimits[9 + centIndex*3][1]));
      fitTxt->AddText(Form("#mu_{#pi} \t %f  [%f,%f]",   parameters_simul_cent[0],  simulParamLimits[0][0], simulParamLimits[0][1]));
      fitTxt->AddText(Form("#sigma_{#pi} \t %f  [%f,%f]",parameters_simul_cent[1],  simulParamLimits[1][0], simulParamLimits[1][1]));
      fitTxt->AddText(Form("#gamma_{#pi} \t %f  [%f,%f]",parameters_simul_cent[2],  simulParamLimits[2][0], simulParamLimits[2][1]));
      fitTxt->AddText(Form("N_{K} \t %e  [%e,%e]",       parameters_simul_cent[10 + centIndex*3], simulParamLimits[10 + centIndex*3][0],simulParamLimits[10 + centIndex*3][1]));
      fitTxt->AddText(Form("#mu_{K} \t %f  [%f,%f]",     parameters_simul_cent[3],  simulParamLimits[3][0], simulParamLimits[3][1]));
      fitTxt->AddText(Form("#sigma_{K} \t %f  [%f,%f]",  parameters_simul_cent[4],  simulParamLimits[4][0], simulParamLimits[4][1]));
      fitTxt->AddText(Form("#gamma_{K} \t %f  [%f,%f]",  parameters_simul_cent[5],  simulParamLimits[5][0], simulParamLimits[5][1]));
      fitTxt->AddText(Form("N_{p} \t %e  [%e,%e]",       parameters_simul_cent[11 + centIndex*3], simulParamLimits[11 + centIndex*3][0],simulParamLimits[11 + centIndex*3][1]));
      fitTxt->AddText(Form("#mu_{p} \t %f  [%f,%f]",     parameters_simul_cent[6],  simulParamLimits[6][0], simulParamLimits[6][1]));
      fitTxt->AddText(Form("#sigma_{p} \t %f  [%f,%f]",  parameters_simul_cent[7],  simulParamLimits[7][0], simulParamLimits[7][1]));
      fitTxt->AddText(Form("#gamma_{p} \t %f  [%f,%f]",  parameters_simul_cent[8],  simulParamLimits[8][0], simulParamLimits[8][1]));
      fitTxt->SetFillStyle(0);
      fitTxt->SetBorderSize(5);
      

      //TH1F* frame_pull = bottomPad->DrawFrame(m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( lowRange - 0.5*sigmaGuess)),-15.0,
      //       m_currentHistosToFit_ByCent[centIndex]->GetBinLowEdge(m_currentHistosToFit_ByCent[centIndex]->FindBin( highRange + 0.5*sigmaGuess)),15.0);
      //frame_pull->SetTitle(Form("; Z_{TPC}; Standardized Pull"));
      if(pullHisto) pullHisto->Draw("E");
      TLine* oneLine = new TLine(lowRange - 0.5*sigmaGuess,0.0,highRange + 0.5*sigmaGuess,0.0);
      oneLine->SetLineStyle(7);
      oneLine->SetLineColor(kBlack);
      oneLine->Draw("same");
      fitTxt->Draw("same");
      gPad->Update();
      gSystem->ProcessEvents();
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_STEP_%02d_mTm0Bin_%02d_%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(), 
             m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex,m_currentRapBin,fitting_round, m_currentMtM0Bin,roundTag[fitting_round].c_str()));

      if(m_saveNoLogImages){
        topPad->SetLogy(false);
        topPad->Update();
        gSystem->ProcessEvents();
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_STEP_%02d_mTm0Bin_%02d_%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(), 
               m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex, m_currentRapBin,fitting_round, m_currentMtM0Bin,roundTag[fitting_round].c_str()));
      }
      if(centIndex == 0){ // only write the all cent data once
        //PION
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][0][1][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[0], simulParamErrors[0]);
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][0][2][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[1], simulParamErrors[1]);
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][0][3][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[2], simulParamErrors[2]);
        //KAON
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][1][1][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[3], simulParamErrors[3]);
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][1][2][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[4], simulParamErrors[4]);
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][1][3][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[5], simulParamErrors[5]);
        //PROTON
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][2][1][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[6], simulParamErrors[6]);
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][2][2][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[7], simulParamErrors[7]);
        HistogramUtilities::setBinValues(m_Fit_Data_ZTPC[m_currentPartIndex][2][3][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[8], simulParamErrors[8]);
      }
      //PION
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][0][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[9 + 3*centIndex], simulParamErrors[9 + 3*centIndex]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][1][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[0], simulParamErrors[0]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][2][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[1], simulParamErrors[1]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][0][centIndex][3][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[2], simulParamErrors[2]);     
      //KAON
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][0][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[10 + 3*centIndex], simulParamErrors[10 + 3*centIndex]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][1][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[3], simulParamErrors[3]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][2][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[4], simulParamErrors[4]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][1][centIndex][3][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[5], simulParamErrors[5]); 
      //PROTON
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][2][centIndex][0][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[11 + 3*centIndex], simulParamErrors[11 + 3*centIndex]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][2][centIndex][1][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[6], simulParamErrors[6]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][2][centIndex][2][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[7], simulParamErrors[7]);
      HistogramUtilities::setBinValues(m_Fit_Data_Cent_ZTPC[m_currentPartIndex][2][centIndex][3][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[8], simulParamErrors[8]); 

      HistogramUtilities::setBinValues(m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex], m_currentRapBin, m_currentMtM0Bin, simulParams[11 + 3*centIndex]*numberOfCentEvents, simulParamErrors[11 + 3*centIndex]*numberOfCentEvents);


      #ifndef _SAVE_POINTERS_
        #ifndef _ZFITTER_PHYSMATH_BYPASS_
          //delete combinedFunct;
          delete pionFunct;
          delete protonFunct;
          delete kaonFunct;
          if(pionFunct_init)   delete pionFunct_init;
          if(kaonFunct_init)   delete kaonFunct_init;
          if(protonFunct_init) delete protonFunct_init;
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



    //DRAW IMAGE WITH ALL CENTRALITIES
    TCanvas* fittingCanvasSimul = new TCanvas("fittingCanvasSimul_SimulTPC","Fitting Canvas");
    fittingCanvasSimul->SetWindowSize(1300,700);
    fittingCanvasSimul->SetCanvasSize(1200,600);
    fittingCanvasSimul->Draw();
    fittingCanvasSimul->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.3);
    gStyle->SetOptFit(0111);
    m_currentHistosToFit_ByCent[0]->SetStats(kFALSE);
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      #ifdef _ZFITTER_DEBUG_
        cout << " Ptr Check: " << m_currentHistosToFit_ByCent[centIndex] << endl;
      #endif
      if(centIndex%2==1){
        m_currentHistosToFit_ByCent[centIndex]->SetLineColor(kBlue);
        m_currentHistosToFit_ByCent[centIndex]->SetMarkerColor(kBlue);
        combinedFunct_CentStorage[centIndex]->SetLineColor(kViolet+1);
      }
      if(centIndex == 0){
        m_currentHistosToFit_ByCent[centIndex]->Draw("E");
      }else{
        m_currentHistosToFit_ByCent[centIndex]->Draw("E SAME");
      }
    }
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      #ifdef _ZFITTER_DEBUG_
        cout << " Ptr Check: " << combinedFunct_CentStorage[centIndex] << endl;
      #endif
      combinedFunct_CentStorage[centIndex]->Draw("same");
    }

    TPaveText* fitParsTxtSimul = new TPaveText(0.7,0.05,0.95,0.95,"NDC");
    fitParsTxtSimul->AddText(Form("#mu_{#pi} \t %f#pm%f %s",simulParams[0],simulParamErrors[0],minimizer->IsFixedVariable(0) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#sigma_{#pi} \t %f#pm%f %s",simulParams[1],simulParamErrors[1],minimizer->IsFixedVariable(1) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#gamma_{#pi} \t %f#pm%f %s",simulParams[2],simulParamErrors[2],minimizer->IsFixedVariable(2) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#mu_{K} \t %f#pm%f %s",simulParams[3],simulParamErrors[3],minimizer->IsFixedVariable(3) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#sigma_{K} \t %f#pm%f %s",simulParams[4],simulParamErrors[4],minimizer->IsFixedVariable(4) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#gamma_{K} \t %f#pm%f %s",simulParams[5],simulParamErrors[5],minimizer->IsFixedVariable(5) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#mu_{p} \t %f#pm%f %s",simulParams[6],simulParamErrors[6],minimizer->IsFixedVariable(6) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#sigma_{p} \t %f#pm%f %s",simulParams[7],simulParamErrors[7],minimizer->IsFixedVariable(7) ? "(fixed)" : "       "));
    fitParsTxtSimul->AddText(Form("#gamma_{p} \t %f#pm%f %s",simulParams[8],simulParamErrors[8],minimizer->IsFixedVariable(8) ? "(fixed)" : "       "));
    fitParsTxtSimul->SetFillStyle(0);
    fitParsTxtSimul->SetBorderSize(5);
    fitParsTxtSimul->Draw("same");

    gPad->Update();
    gSystem->ProcessEvents();
    fittingCanvasSimul->SaveAs(Form("%s/%s/%s/dEdxFits/SimulCent_RapBin_%02d_STEP_%02d_mTm0Bin_%02d_%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(), 
         m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, fitting_round, m_currentMtM0Bin,roundTag[fitting_round].c_str()));


    #ifndef _SAVE_POINTERS_
      for(int iii = 0; iii < 16; iii++){
        if(combinedFunct_CentStorage[iii]) delete combinedFunct_CentStorage[iii];
      }
      delete fitParsTxtSimul;
      delete fittingCanvasSimul;
      delete topPad;
      delete bottomPad;
      delete fittingCanvas;
      delete chiSqrdFunctor;
      delete minimizer;
    #endif


    mTm0BinsFitPerRound[fitting_round] = mTm0BinsFitPerRound[fitting_round] + 1; 

  }// End mTm0 Loop
}



#endif