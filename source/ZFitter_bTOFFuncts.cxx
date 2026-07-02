// THIS IS INSERTED INTO ZFITTER.cxx // Just separate file becuase it's long
//============================================================================================================
//#####      INVERSE BETA FITTING    #################
//============================================================================================================


//=======================================================================================================================
void ZFitter::initBTOFFunction(TF1* a_funct){
  int paramIndexCounter = 0;
  double integral = m_currentHistoToFit->Integral();
  for(unsigned int iii = 0; iii <  m_currentParameters.size(); iii++){




    if(m_paramToggles[iii]){
         //cout << " Funct Set Par " <<  paramIndexCounter << " to " <<  m_currentParameters[iii] << endl;
      a_funct->SetParameter(paramIndexCounter, m_currentParameters[iii]);
      a_funct->SetParName(paramIndexCounter, m_currentParamNames[paramIndexCounter].c_str());
      if(iii%3==0){
        if(m_currentParameters[iii] > 0.0 && m_currentParameters[iii] < 2.0){
          a_funct->SetParLimits(paramIndexCounter,0.0,1.1);
          //a_funct->SetParLimits(paramIndexCounter,0.0,1.3*integral);
        }else{
          cout << "ERROR:  bad amp init..." << endl;
          a_funct->SetParLimits(paramIndexCounter,0.0,1.1);
          a_funct->SetParameter(paramIndexCounter,0.5);
          //a_funct->SetParLimits(paramIndexCounter,0.0, 1.3*integral);
          //a_funct->SetParameter(paramIndexCounter, 0.5*integral);
        }
      }
      if(iii%3==1) a_funct->SetParLimits(paramIndexCounter,m_currentParameters[iii] - 0.01, m_currentParameters[iii] + 0.01);
      if(iii%3==2) a_funct->SetParLimits(paramIndexCounter,0.1*m_currentParameters[iii], m_currentParameters[iii]*10.0);

      if(m_massSquaredMode){
        if(iii%3==1) a_funct->SetParLimits(paramIndexCounter,m_currentParameters[iii] - 0.05, m_currentParameters[iii] + 0.05);
        if(iii%3==2) a_funct->SetParLimits(paramIndexCounter,0.1*m_currentParameters[iii], m_currentParameters[iii]*10.0);
      }
      paramIndexCounter++;
    }
  }
}

//=======================================================================================================================================================================
//=======================================================================================================================================================================
void ZFitter::initBTOFParams(){
  vector<double> params;
  vector<string> paramNames;
  vector<bool> paramToggles;
  m_currentFreeParameters = 0;
  //m_currentNumPartToFit = 0;
  double sigmaGuess = 0.007;
  //if(m_massSquaredMode) sigmaGuess = 0.05;

  //matrix inversion for amplitudes
  /*vector< double > matrixElements;
  vector< double > countsVector;
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_currentParticleTags[partIndex] == 0) continue;
    m_currentNumPartToFit++;
    double xPoint = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin);
    int predictBin = m_currentHistoToFit->FindBin(xPoint);
    double binContent = m_currentHistoToFit->GetBinContent(predictBin);
    countsVector.push_back(binContent);
    for(int partIndex2 = 0; partIndex2 < m_numParticles; partIndex2++){

      if(m_currentParticleTags[partIndex2] == 1 || m_currentParticleTags[partIndex2] == 3){
        double mu = m_InvBetaBTOFFitData[m_currentPartIndex][partIndex2[1]]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        double sig = m_InvBetaBTOFFitData[m_currentPartIndex][partIndex2[2]]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        double prob = gaus(mu,sig,xPoint);
        if(prob < 1E-5) matrixElements.push_back( 0.0 );
        else matrixElements.push_back( prob );
      }
      if(m_currentParticleTags[partIndex2] == 2){
        double mean = m_invBetaCurves[partIndex2]->Eval(m_currentMomOfBin);
        double stdDev = sigmaGuess;
        double prob = gaus(mean,stdDev,xPoint);
        if(prob < 1E-5) matrixElements.push_back( 0.0 );
        else matrixElements.push_back( prob );
      }
    }
  }


  cout << "Matrix STDVectorSize: " << matrixElements.size() << " NumPart: " << m_currentNumPartToFit << endl;
  TMatrixD mat(m_currentNumPartToFit, m_currentNumPartToFit, &matrixElements[0]);
  mat.Print();
  cout << "Det: " << mat.Determinant() << endl;

  TDecompLU lu(mat);
  TMatrixD matInv(m_currentNumPartToFit, m_currentNumPartToFit);
  if (!lu.Decompose()) {
      cout << "Decomposition failed, matrix singular ?" << endl;
      //cout << "condition number = " <<  mat.GetCondition() << endl;
   } else {
      lu.Invert(matInv);
   }

  matInv.Print();
  cout << "Det: " << matInv.Determinant() << endl;

  //TMatrixD unit = mat * matInv;
  //unit.Print();

  TVectorD countVect(m_currentNumPartToFit);
  for(int iii = 0; iii < m_currentNumPartToFit; iii++){
    countVect(iii) = countsVector[iii];
  }
  countVect.Print();

  TVectorD ampVect = mat * countVect;
  ampVect.Print();
  */



  //cout << "Initializing  Parameters" << endl;
  int partIndexCounter = 0;
  if(!m_separatePM_bToF_allCent){
    m_currentPlusMinusIndex = 0;
  }else{
    m_currentPlusMinusIndex = 0;
    if(m_currentPartCharge == -1) m_currentPlusMinusIndex = 1;
  }

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){


    if(m_currentParticleTags[partIndex] == 0) continue;

    if(m_massSquaredMode){
      //sigmaGuess = PhysMath::m2_width_prediction(35.0E-11,m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(partIndex),1.0);
      double gausDataArray[3];
      HistogramUtilities::getGausData(&gausDataArray[0], m_currentHistoToFit, true, false, m_cuts->getLowMassSquaredCut(partIndex), m_cuts->getLowMassSquaredCut(partIndex));
      sigmaGuess = gausDataArray[2];
      //sigParGuess = sqrt(sigmaGuess*sigmaGuess*(ndfParamGuess - 2.0)/ndfParamGuess);
    }

    //#####    TAG 1 #######    if fitting by centrality and already have widths and means
    if(m_currentParticleTags[partIndex] == 1){
      double mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      double invBetaPredictBin = m_currentHistoToFit->FindBin(mean);
      double sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      //double ampAtPredict = ampVect[partIndexCounter];
      double gausProb = gaus(mean, sig, mean);
      double ampAtPredict = 0.0;
      if(gausProb > 0.0){
        //ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / gausProb;
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
        //  previousAmpContributionFromOtherParts += m_currentIntegralofHistogram * params[paramIndex]*gaus(params[paramIndex+1],params[paramIndex+2],mean);
        //}
        //ampAtPredict = (m_currentHistoToFit->GetBinContent(invBetaPredictBin) - previousAmpContributionFromOtherParts)  / (gausProb * m_currentIntegralofHistogram);
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (gausProb*m_currentIntegralofHistogram);
        if(ampAtPredict <= 0.0) ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;

      }else{
        cout << "ERROR: PREDICTION OF AMPLITUDE FAILED!!!!   Amp = NAN " << endl;
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }

      params.push_back(ampAtPredict);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));
      params.push_back(mean);
      paramToggles.push_back(false);
      params.push_back(sig);
      paramToggles.push_back(false);

      partIndexCounter++;
      //m_currentNumPartToFit++;
      m_currentFreeParameters+=1;

      if(std::isnan(mean) || std::isnan(sig) || std::isnan(ampAtPredict) || sig==0){
        cout << "ERROR: Tag = 1 : " << m_partInfo->GetParticleName(partIndex) << "  N: " << ampAtPredict << "  Mean: " << mean << " Sig: " << sig << endl;
        unsigned int testIndex = 3*(partIndexCounter-1);// KILL the particle
        paramToggles[testIndex] = false; paramToggles[testIndex+1] = false; paramToggles[testIndex+2] = false;
        params[testIndex] = 0.0; params[testIndex+1] = 0.0; params[testIndex+2] = 1.0;
      }
      #ifdef _ZFITTER_DEBUG_
        cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 1 " << "  N:" << ampAtPredict << "  mu:" << mean << "(fixed)  sig:" << sig << "(fixed) " <<  endl;
      #endif
    }







     //#####    TAG 2 #######    Need to fit all parameters
    else if(m_currentParticleTags[partIndex] == 2){

      double mean = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
      if(m_massSquaredMode) mean = pow(m_partInfo->GetParticleMass(partIndex),2.0);
      int invBetaPredictBin = m_currentHistoToFit->FindBin(mean);
      double sig = sigmaGuess;
      if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
        mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
        sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      }else{
        bool gotOtherParticleSeed = false;
        if(m_useFitInfoForOppositeCharge && m_seedFromOppositeChargeDirection[partIndex] > 0 && m_currentPartCharge < 0.0 && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0){
            mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            gotOtherParticleSeed = true;
        }
        if(m_useFitInfoForOppositeCharge && m_seedFromOppositeChargeDirection[partIndex] < 0 && m_currentPartCharge > 0.0 && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0){
            mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            gotOtherParticleSeed = true;
        }

        if(!gotOtherParticleSeed && m_useMirroredRapidityAsSeed){
          double currentCenter = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
          double mirrorBin = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->FindBin(-currentCenter);
          if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->GetBinCenter(m_currentRapBin) > 0.0001
            && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin) > 0.0){
            mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin);
            sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin);
          }
        }
      }

      double gausProb = gaus(mean, sig, mean);
      double ampAtPredict = 0.0;
      if(gausProb > 0.0){
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
        //  previousAmpContributionFromOtherParts += params[paramIndex]*gaus(params[paramIndex+1],params[paramIndex+2],mean);
        //}
        //ampAtPredict = (m_currentHistoToFit->GetBinContent(invBetaPredictBin) - previousAmpContributionFromOtherParts)  / (gausProb * m_currentIntegralofHistogram);
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (gausProb*m_currentIntegralofHistogram);
        if(ampAtPredict <= 0.0) ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }else{
        cout << "ERROR: PREDICTION OF AMPLITUDE FAILED!!!!   Amp = NAN " << endl;
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }


      params.push_back(ampAtPredict);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(mean);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(sig);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));


      //m_currentNumPartToFit++;
      partIndexCounter++;
      m_currentFreeParameters+=3;

      if(std::isnan(mean) || std::isnan(sig) || std::isnan(ampAtPredict) || sig==0){
        cout << "ERROR: Tag = 2 : " << m_partInfo->GetParticleName(partIndex) << "  N: " << ampAtPredict << "  Mean: " << mean << " Sig: " << sig << endl;
        unsigned int testIndex = 3*(partIndexCounter-1); // KILL the particle
        paramToggles[testIndex] = false; paramToggles[testIndex+1] = false; paramToggles[testIndex+2] = false;
        params[testIndex] = 0.0; params[testIndex+1] = 0.0; params[testIndex+2] = 1.0;
      }
      #ifdef _ZFITTER_DEBUG_
        cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 2 " << "  N:" << ampAtPredict << "  mu:" << mean << "  sig:" << sig << " " <<  endl;
      #endif

    }
    else if(m_currentParticleTags[partIndex] == 3){

      if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) == 0.0){
        cout << "ERROR: particle tag is 3, but no fit data stored in AmpHisto!" << endl;
      }
      params.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) / m_currentIntegralofHistogram);
      paramToggles.push_back(false);
      params.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      paramToggles.push_back(false);
      params.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      paramToggles.push_back(false);
      partIndexCounter++;
      //m_currentNumPartToFit++;

      #ifdef _ZFITTER_DEBUG_
        cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 3 "
        << "  N:" << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)  / m_currentIntegralofHistogram
        << "(fixed)  mu:" << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)
        << "(fixed)  sig:" << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) << " (fixed)" <<  endl;
      #endif
    }
    else if(m_currentParticleTags[partIndex] == 9){

      double mean = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
      if(m_massSquaredMode) mean = pow(m_partInfo->GetParticleMass(partIndex),2.0);

      int invBetaPredictBin = m_currentHistoToFit->FindBin(mean);
      double sig = sigmaGuess;

      if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
        mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
        sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      }

      double gausProb = gaus(mean, sig, mean);
      double ampAtPredict = 0.0;
      if(gausProb > 0.0){
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
        //  previousAmpContributionFromOtherParts += params[paramIndex]*gaus(params[paramIndex+1],params[paramIndex+2],mean);
        //}
        //ampAtPredict = (m_currentHistoToFit->GetBinContent(invBetaPredictBin) - previousAmpContributionFromOtherParts)  / (gausProb * m_currentIntegralofHistogram);
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (gausProb*m_currentIntegralofHistogram);
        if(ampAtPredict <= 0.0) ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5/m_currentIntegralofHistogram;
      }else{
        cout << "ERROR: PREDICTION OF AMPLITUDE FAILED!!!!   Amp = NAN " << endl;
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5/m_currentIntegralofHistogram;
      }

      params.push_back(ampAtPredict);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(mean);
      paramToggles.push_back(false);
      //paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(sig);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));



      //m_currentNumPartToFit++;
      partIndexCounter++;
      m_currentFreeParameters+=2;

      if(std::isnan(mean) || std::isnan(sig) || std::isnan(ampAtPredict)  || sig==0 ){
        cout << "ERROR: Tag = 2 : " << m_partInfo->GetParticleName(partIndex) << "  N: " << ampAtPredict << "  Mean: " << mean << " Sig: " << sig  << endl;
        unsigned int testIndex = 4*(partIndexCounter-1);// KILL the particle
        paramToggles[testIndex] = false; paramToggles[testIndex+1] = false; paramToggles[testIndex+2] = false;
        params[testIndex] = 0.0; params[testIndex+1] = 0.0; params[testIndex+2] = 1.0;
      }
      #ifdef _ZFITTER_DEBUG_
        cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 2 " << "  N:" << ampAtPredict << "  mu:" << mean << "  sig:" << sig   <<  endl;
      #endif

    }

  }//end particle loop


  m_currentParameters = params;
  m_paramToggles = paramToggles;
  m_currentParamNames = paramNames;
  //printCurrentParams();


}

//============================================================================================================
//============================================================================================================

void ZFitter::initBTOFFunction_student(TF1* a_funct){
  int paramIndexCounter = 0;
  //double integral = m_currentHistoToFit->Integral();




  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
    int partIndex = m_currentLiveParticleIndicies[iii/4];
    if(m_paramToggles[iii]){
      //cout << " Funct Set Par " <<  paramIndexCounter << " to " <<  m_currentParameters[iii] << endl;
      a_funct->SetParameter(paramIndexCounter, m_currentParameters[iii]);
      a_funct->SetParName(paramIndexCounter, m_currentParamNames[paramIndexCounter].c_str());
      if(iii%4==0){
        if(m_currentParameters[iii] > 0.0 && m_currentParameters[iii] < 2.0){
          a_funct->SetParLimits(paramIndexCounter,0.0,1.3);
          //a_funct->SetParLimits(paramIndexCounter,0.0,1.3*integral);
        // cout << "  Limits on param " << paramIndexCounter << " are 0 to " << 3.0*m_currentParameters[iii] << endl;
        }else{
          a_funct->SetParameter(paramIndexCounter,0.5);
          a_funct->SetParLimits(paramIndexCounter,0.0,1.3);
          //a_funct->SetParameter(paramIndexCounter,0.5);
          //a_funct->SetParLimits(paramIndexCounter,0.0, 1.3*integral);
          //a_funct->SetParameter(paramIndexCounter, 0.5*integral);
          //cout << "  Limits on param " << paramIndexCounter << " are 0 to 1000" << endl;
        }
      }

      if(!m_inAllCentMode && m_allowSpecCentFitsToVary_BTOF){
        if(iii%4==1){
          if(!m_massSquaredMode){
            a_funct->SetParLimits(paramIndexCounter,m_currentParameters[iii] - 0.02, m_currentParameters[iii] + 0.02);
          }else{
            a_funct->SetParLimits(paramIndexCounter,0.5*(m_currentParameters[iii] + m_cuts->getLowMassSquaredCut(partIndex)),0.5*(m_currentParameters[iii] + m_cuts->getHighMassSquaredCut(partIndex)));
          }
        }
        else if(iii%4==2) {
          a_funct->SetParLimits(paramIndexCounter,0.01*m_currentParameters[iii], m_currentParameters[iii]*100.0);
        }
        else if(iii%4==3) {
          a_funct->SetParLimits(paramIndexCounter,m_NdofStudentTMinimum, 150.0); // if you set the lower limit less than 4ish, the tails contribute more to the integral. You'll get low DOF at high pT and the spectra will seem to have a higher temperature
          //YOU CAN"T CHANGE THE UPPER LIMIT TO OVER 350! You will get NAN as output in TMath::Student()!!!!!!!!!!!!!!!!!!!
        }
      }else{
        if(iii%4==1){
          if(!m_massSquaredMode){
            a_funct->SetParLimits(paramIndexCounter,m_currentParameters[iii] - 0.02, m_currentParameters[iii] + 0.02);
          }else{
            a_funct->SetParLimits(paramIndexCounter,0.5*(m_currentParameters[iii] + m_cuts->getLowMassSquaredCut(partIndex)),0.5*(m_currentParameters[iii] + m_cuts->getHighMassSquaredCut(partIndex)));
          }        }
        else if(iii%4==2) {
          a_funct->SetParLimits(paramIndexCounter,0.01*m_currentParameters[iii], m_currentParameters[iii]*100.0);
        }
        else if(iii%4==3) {
          a_funct->SetParLimits(paramIndexCounter,m_NdofStudentTMinimum, 150.0); // if you set the lower limit less than 4ish, the tails contribute more to the integral.
          //You'll get low DOF at high pT and the spectra will seem to have a higher temperature
          //YOU CAN"T CHANGE THE UPPER LIMIT TO OVER 350! You will get NAN as output in TMath::Student()!!!!!!!!!!!!!!!!!!!
        }
      }

      paramIndexCounter++;
    }
  }
}

//=======================================================================================================================

void ZFitter::initBTOFParams_student(){
  vector<double> params;
  vector<string> paramNames;
  vector<bool> paramToggles;
  m_currentLiveParticleIndicies.clear();
  m_currentFreeParameters = 0;
  //m_currentNumPartToFit = 0;
  double sigmaGuess = 0.011;
  double ndfParamGuess = 15.0;
  double sigParGuess = sqrt(sigmaGuess*sigmaGuess*(ndfParamGuess - 2.0)/ndfParamGuess); // convert gaussian std dev to sigma of student T
  double sigParGuessBackup = sigParGuess;

  //matrix inversion for amplitudes
  /*vector< double > matrixElements;
  vector< double > countsVector;
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_currentParticleTags[partIndex] == 0) continue;
    m_currentNumPartToFit++;
    double xPoint = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin);
    int predictBin = m_currentHistoToFit->FindBin(xPoint);
    double binContent = m_currentHistoToFit->GetBinContent(predictBin);
    countsVector.push_back(binContent);
    for(int partIndex2 = 0; partIndex2 < m_numParticles; partIndex2++){

      if(m_currentParticleTags[partIndex2] == 1 || m_currentParticleTags[partIndex2] == 3){
        double mu = m_InvBetaBTOFFitData[m_currentPartIndex][partIndex2[1]]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        double sig = m_InvBetaBTOFFitData[m_currentPartIndex][partIndex2[2]]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        double prob = gaus(mu,sig,xPoint);
        if(prob < 1E-5) matrixElements.push_back( 0.0 );
        else matrixElements.push_back( prob );
      }
      if(m_currentParticleTags[partIndex2] == 2){
        double mean = m_invBetaCurves[partIndex2]->Eval(m_currentMomOfBin);
        double stdDev = sigmaGuess;
        double prob = gaus(mean,stdDev,xPoint);
        if(prob < 1E-5) matrixElements.push_back( 0.0 );
        else matrixElements.push_back( prob );
      }
    }
  }


  cout << "Matrix STDVectorSize: " << matrixElements.size() << " NumPart: " << m_currentNumPartToFit << endl;
  TMatrixD mat(m_currentNumPartToFit, m_currentNumPartToFit, &matrixElements[0]);
  mat.Print();
  cout << "Det: " << mat.Determinant() << endl;

  TDecompLU lu(mat);
  TMatrixD matInv(m_currentNumPartToFit, m_currentNumPartToFit);
  if (!lu.Decompose()) {
      cout << "Decomposition failed, matrix singular ?" << endl;
      //cout << "condition number = " <<  mat.GetCondition() << endl;
   } else {
      lu.Invert(matInv);
   }

  matInv.Print();
  cout << "Det: " << matInv.Determinant() << endl;

  //TMatrixD unit = mat * matInv;
  //unit.Print();

  TVectorD countVect(m_currentNumPartToFit);
  for(int iii = 0; iii < m_currentNumPartToFit; iii++){
    countVect(iii) = countsVector[iii];
  }
  countVect.Print();

  TVectorD ampVect = mat * countVect;
  ampVect.Print();
  */
  //cout << "Initializing  Parameters" << endl;

  int partIndexCounter = 0;
  if(!m_separatePM_bToF_allCent){
    m_currentPlusMinusIndex = 0;
  }else{
    m_currentPlusMinusIndex = 0;
    if(m_currentPartCharge == -1) m_currentPlusMinusIndex = 1;
  }


  //HistogramUtilities::printHistoInfo(m_currentHistoToFit);


  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_currentParticleTags[partIndex] == 0) continue;

    m_currentLiveParticleIndicies.push_back(partIndex);

    double gausDataArray[3];
    if(m_massSquaredMode){
      //sigmaGuess = PhysMath::m2_width_prediction(35.0E-11,m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(partIndex),1.0);
      //cout << " Setup Tag for Part : " << partIndex << " with m^2 from " << m_cuts->getLowMassSquaredCut(partIndex) << " to " << m_cuts->getHighMassSquaredCut(partIndex) << endl;
      HistogramUtilities::getGausData(&gausDataArray[0], m_currentHistoToFit, true, false, m_cuts->getLowMassSquaredCut(partIndex), m_cuts->getHighMassSquaredCut(partIndex)); //bools are is differential and use errors
      sigmaGuess = 0.75*gausDataArray[2];
      //if(sigmaGuess > 0.1) sigmaGuess = 0.1;
      if(sigmaGuess <= 0.0){
        sigParGuess = sigParGuessBackup;
      }else{
        sigParGuess = sqrt(sigmaGuess*sigmaGuess*(ndfParamGuess - 2.0)/ndfParamGuess);
      }
      #ifdef _ZFITTER_DEBUG_
        cout << "    with initial guess N: " << gausDataArray[0] << "  mu: " << gausDataArray[1] << "   stdDev: " << gausDataArray[2] << "   sigPar: " << sigParGuess << endl;
      #endif
    }
    #ifdef _ZFITTER_DEBUG_
      cout << " Setup Tag for Part : " << partIndex << " with m^2 from " << m_cuts->getLowMassSquaredCut(partIndex) << " to " << m_cuts->getHighMassSquaredCut(partIndex) << endl;
    #endif

    //##### TAG 1 #####   if fitting buy centrality and already have widths and means
    if(m_currentParticleTags[partIndex] == 1){
      double mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      int invBetaPredictBin = m_currentHistoToFit->FindBin(mean);
      double sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      double ndf = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      if(ndf < 1.1) ndf = 1.1;

      double sTProb = studentT(mean, sig, ndf, mean);
      double ampAtPredict = 0.0;
      if(sTProb > 0.0){
        //ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (sTProb * m_currentIntegralofHistogram);
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 4){
        //  previousAmpContributionFromOtherParts += m_currentIntegralofHistogram * params[paramIndex]*studentT(params[paramIndex+1],params[paramIndex+2],params[paramIndex+3],mean);
        //}
        //ampAtPredict = (m_currentHistoToFit->GetBinContent(invBetaPredictBin) - previousAmpContributionFromOtherParts)  / (sTProb * m_currentIntegralofHistogram);
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (sTProb*m_currentIntegralofHistogram);
        if(ampAtPredict <= 0.0) ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }else{
        cout << "ERROR: PREDICTION OF AMPLITUDE FAILED!!!!   Amp = NAN " << endl;
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }

      params.push_back(ampAtPredict);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(mean);
      paramToggles.push_back(false);

      params.push_back(sig);
      paramToggles.push_back(false);

      params.push_back(ndf);
      paramToggles.push_back(false);

      partIndexCounter++;
      //m_currentNumPartToFit++;
      m_currentFreeParameters+=1;

      if(std::isnan(mean) || std::isnan(sig) || std::isnan(ampAtPredict) || std::isnan(ndf) || sig==0 || ndf <= 1.0){
        cout << "ERROR: Tag = 1 : " << m_partInfo->GetParticleName(partIndex) << "  N: " << ampAtPredict << "  Mean: " << mean << "(fixed) Sig: " << sig << "(fixed)  NDF: " << ndf  << " (fixed) " << endl;
        unsigned int testIndex = 4*(partIndexCounter-1);// KILL the particle
        paramToggles[testIndex] = false; paramToggles[testIndex+1] = false; paramToggles[testIndex+2] = false; paramToggles[testIndex+3] = false;
        params[testIndex] = 0.0; params[testIndex+1] = 0.0; params[testIndex+2] = 1.0; params[testIndex+3] = 2.0;
      }
      #ifdef _ZFITTER_DEBUG_
        cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 1 " << "  N:" << ampAtPredict << "  mu:" << mean << "(fixed)  sig:" << sig << "(fixed) ndof: " << ndf << " (fixed)"  <<  endl;
      #endif
    }

    else if(m_currentParticleTags[partIndex] == 2){

      double mean = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
      if(m_massSquaredMode){
        double simpleMass2Predict = pow(m_partInfo->GetParticleMass(partIndex),2.0);
        mean = gausDataArray[1];
        if(fabs(simpleMass2Predict - mean)/simpleMass2Predict > 0.1) mean = simpleMass2Predict;
        if(fabs(simpleMass2Predict - mean) > 0.1) mean = simpleMass2Predict;
      }




      int invBetaPredictBin = m_currentHistoToFit->FindBin(mean);
      double sig = sigParGuess;
      double ndf = ndfParamGuess;


      if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
        mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
        sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
        ndf = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);

        #ifdef _ZFITTER_DEBUG_
          cout << "  Fit Data Exists"  <<  endl;
        #endif

      }else{

        bool gotOtherParticleSeed = false;
        if(m_useFitInfoForOppositeCharge && m_seedFromOppositeChargeDirection[partIndex] > 0 && m_currentPartCharge < 0.0 && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0){
            mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            ndf = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][0]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
            gotOtherParticleSeed = true;
            #ifdef _ZFITTER_DEBUG_
              cout << "  Fit Data Exists for Positive Particles, using that as a seed..."  <<  endl;
            #endif
        }
        if(m_useFitInfoForOppositeCharge && m_seedFromOppositeChargeDirection[partIndex] < 0 && m_currentPartCharge > 0.0 && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0){
            mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            ndf = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][1]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
            gotOtherParticleSeed = true;
            #ifdef _ZFITTER_DEBUG_
              cout << "  Fit Data Exists for Negative Particles, using that as a seed..."  <<  endl;
            #endif
        }

        if(!gotOtherParticleSeed && m_useMirroredRapidityAsSeed){
          double currentCenter = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
          int mirrorBin = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->FindBin(-currentCenter);
          if(mirrorBin > 0 && mirrorBin <= m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetNbinsX()               // mirror bin is within range of histogram
            && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->GetBinCenter(m_currentRapBin) > 0.0001    // we are looking at a positive rapidity (nonZero)
            && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin) > 0.0){     // the mirrored bin has been fit
            if( (m_inAllCentMode && m_FitStatus_ZbTOF[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) == 0 ) // the other fit was good
              ||(!m_inAllCentMode && m_FitStatus_Cent_ZbTOF[m_currentPartIndex][m_currentCentIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) == 0  ) ){
              mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin);
              sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin);
              ndf = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->GetBinContent(mirrorBin, m_currentMtM0Bin);
              #ifdef _ZFITTER_DEBUG_
                cout << "  Fit Data Exists for the mirrored rapidity (" << -currentCenter << "), using that as a seed..."  <<  endl;
              #endif
            }

          }
        }

        if(m_massSquaredMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin - 1) > 0.0){
          mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin - 1);
          sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin - 1);
          ndf = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin - 1);
        }
      }






      if(ndf < 1.0000001) ndf = 1.0000001;

      double sTProb = studentT(mean, sig, ndf, mean);
      double ampAtPredict = 0.0;
      if(sTProb > 0.0){
        //ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / sTProb;
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 4){
        //  previousAmpContributionFromOtherParts += m_currentIntegralofHistogram * params[paramIndex]*studentT(params[paramIndex+1],params[paramIndex+2],params[paramIndex+3],mean);
        //}
        //ampAtPredict = (m_currentHistoToFit->GetBinContent(invBetaPredictBin) - previousAmpContributionFromOtherParts)  / (sTProb * m_currentIntegralofHistogram);
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (sTProb*m_currentIntegralofHistogram);
        if(ampAtPredict <= 0.0) ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }else{
        cout << "ERROR: PREDICTION OF AMPLITUDE FAILED!!!!   Amp = NAN " << endl;
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }

      params.push_back(ampAtPredict);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(mean);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(sig);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(ndf);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#nu_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));


      //m_currentNumPartToFit++;
      partIndexCounter++;
      m_currentFreeParameters+=4;

      if(std::isnan(mean) || std::isnan(sig) || std::isnan(ampAtPredict) || std::isnan(ndf) || sig==0 || ndf <= 1.0){
        cout << "ERROR: Tag = 2 : " << m_partInfo->GetParticleName(partIndex) << "  N: " << ampAtPredict << "  Mean: " << mean << " Sig: " << sig << "  NDF: " << ndf << endl;
        unsigned int testIndex = 4*(partIndexCounter-1);// KILL the particle
        paramToggles[testIndex] = false; paramToggles[testIndex+1] = false; paramToggles[testIndex+2] = false; paramToggles[testIndex+3] = false;
        params[testIndex] = 0.0; params[testIndex+1] = 0.0; params[testIndex+2] = 1.0; params[testIndex+3] = 2.0;
      }
      #ifdef _ZFITTER_DEBUG_
        cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 2 " << "  N:" << ampAtPredict << "  mu:" << mean << "  sig:" << sig << " ndof: " << ndf << " "  <<  endl;
      #endif

    }
    else if(m_currentParticleTags[partIndex] == 3){

      if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) == 0.0){
        cout << "ERROR: particle tag is 3, but no fit data stored in AmpHisto!" << endl;
      }
      params.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) / m_currentIntegralofHistogram);
      paramToggles.push_back(false);

      params.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      paramToggles.push_back(false);

      params.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      paramToggles.push_back(false);

      params.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      paramToggles.push_back(false);


      partIndexCounter++;
      //m_currentNumPartToFit++;
      #ifdef _ZFITTER_DEBUG_
      cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 3 "
      << "  N:" << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)
      << "  mu:" << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)
      << "  sig:" << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)
      << "  ndf:" << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)
      <<  endl;
      #endif
    }
    else if(m_currentParticleTags[partIndex] == 9){

      double mean = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
      if(m_massSquaredMode){
        double simpleMass2Predict = pow(m_partInfo->GetParticleMass(partIndex),2.0);
        mean = gausDataArray[1];
        if(fabs(simpleMass2Predict - mean)/simpleMass2Predict > 0.1) mean = simpleMass2Predict;
      }
      int invBetaPredictBin = m_currentHistoToFit->FindBin(mean);
      double sig = sigParGuess;
      double ndf = ndfParamGuess;

      if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
        mean = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
        sig = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
        ndf = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
      }

      if(ndf < 1.1) ndf = 1.1;
      double sTProb = studentT(mean, sig, ndf, mean);
      double ampAtPredict = 0.0;
      if(sTProb > 0.0){
        //ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (sTProb * m_currentIntegralofHistogram);
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 4){
        //  previousAmpContributionFromOtherParts += m_currentIntegralofHistogram * params[paramIndex]*studentT(params[paramIndex+1],params[paramIndex+2],params[paramIndex+3],mean);
        //}
        //ampAtPredict = (m_currentHistoToFit->GetBinContent(invBetaPredictBin) - previousAmpContributionFromOtherParts)  / (sTProb * m_currentIntegralofHistogram);
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) / (sTProb*m_currentIntegralofHistogram);
        if(ampAtPredict <= 0.0) ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }else{
        cout << "ERROR: PREDICTION OF AMPLITUDE FAILED!!!!   Amp = NAN " << endl;
        ampAtPredict = m_currentHistoToFit->GetBinContent(invBetaPredictBin) * 0.5 / m_currentIntegralofHistogram;
      }

      params.push_back(ampAtPredict);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(mean);
      paramToggles.push_back(false);
      //paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(sig);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(ndf);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#nu_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));


      //m_currentNumPartToFit++;
      partIndexCounter++;
      m_currentFreeParameters+=3;

      if(std::isnan(mean) || std::isnan(sig) || std::isnan(ampAtPredict) || std::isnan(ndf) || sig==0 || ndf <= 1.0){
        cout << "ERROR: Tag = 2 : " << m_partInfo->GetParticleName(partIndex) << "  N: " << ampAtPredict << "  Mean: " << mean << " Sig: " << sig << "  NDF: " << ndf << endl;
        unsigned int testIndex = 4*(partIndexCounter-1);// KILL the particle
        paramToggles[testIndex] = false; paramToggles[testIndex+1] = false; paramToggles[testIndex+2] = false; paramToggles[testIndex+3] = false;
        params[testIndex] = 0.0; params[testIndex+1] = 0.0; params[testIndex+2] = 1.0; params[testIndex+3] = 2.0;
      }
      #ifdef _ZFITTER_DEBUG_
        cout << "     " << m_partInfo->GetParticleName(partIndex,0) << " Tag: 2 " << "  N:" << ampAtPredict << "  mu:" << mean << "  sig:" << sig << " ndof: " << ndf << " "  <<  endl;
      #endif

    }
  }//end particle loop
  m_currentParameters = params;
  m_paramToggles = paramToggles;
  m_currentParamNames = paramNames;
  #ifdef _ZFITTER_DEBUG_
    cout << "Parameters after function initialization: " << endl;
    printCurrentParams();
    cout << "parameter names: " << endl;
    for(unsigned int iii = 0; iii < m_currentParamNames.size(); iii++) cout << "  " << m_currentParamNames[iii];
    cout << endl << "  parameter toggles : " << endl;
    for(unsigned int iii = 0; iii < m_paramToggles.size(); iii++) cout << "  " << m_paramToggles[iii];
  #endif
}



//==================================================================================================================================
//==================================================================================================================================

vector<double> ZFitter::fitIsolation_BTOF(int a_particleIndex){

  double ndfParamGuess = 8.0; // >2

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas_iso_tof","Fitting Canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);

  double lowRange = 0;
  double highRange = 0;

  if(a_particleIndex == 0){
    lowRange = 0.8;
    highRange = 0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(0), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(1), m_currentMomOfBin));
    if(!m_initialParticleToggles[1]){
      highRange = 0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(0), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(2), m_currentMomOfBin));
    }
  }else if(a_particleIndex == 2){
    lowRange = 0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(2), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(1), m_currentMomOfBin));
    highRange = 0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(2), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(4), m_currentMomOfBin));
    cout << "   proton range: " << lowRange << "  " << highRange << endl;
  }else if(a_particleIndex == 4){
    lowRange = 0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(4), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(2), m_currentMomOfBin));
    highRange =0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(4), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(5), m_currentMomOfBin));
  }else{
    lowRange = 0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(a_particleIndex), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(a_particleIndex-1), m_currentMomOfBin));
    highRange = 0.5*(PhysMath::invBeta(m_partInfo->GetParticleMass(a_particleIndex), m_currentMomOfBin) + PhysMath::invBeta(m_partInfo->GetParticleMass(a_particleIndex+1), m_currentMomOfBin));
  }


  lowRange = lowRange - PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), m_currentMomOfBin);
  highRange = highRange - PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), m_currentMomOfBin);

  TH1D* zoomHisto = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);

  double gausData[3];
  HistogramUtilities::getGausData(gausData, zoomHisto, true, true, lowRange, highRange);
  double counts = HistogramUtilities::sumOfContent(zoomHisto, -1, true);

  //double integral = HistogramUtilities::integral(zoomHisto);
  //double mean = zoomHisto->GetMean();
  //double stdDev = zoomHisto->GetStdDev();


  double integral = gausData[0];
  double mean = gausData[1];
  double stdDev = gausData[2];

  vector<double> returnParams;
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);

  if(integral == 0 || std::isnan(mean) || std::isnan(stdDev)) return returnParams;


  TString imageName = TString::Format("%s_rap_%02d_mTm0_%02d",
   m_partInfo->GetParticleName(a_particleIndex,0).Data(),m_currentRapBin,m_currentMtM0Bin);
  zoomHisto->SetName(imageName.Data());
  zoomHisto->SetTitle(Form("BTOF Particles for y=[%1.2f,%1.2f] m_{T}-m{0}=[%2.3f,%2.3f] (GeV/c) p=%2.3f : zoom on %s",
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
    m_currentMomOfBin,
    m_partInfo->GetParticleName(a_particleIndex,0).Data()
    ));
  zoomHisto->GetXaxis()->SetTitle("#beta^{-1}");
  zoomHisto->GetYaxis()->SetTitle("dN/d#beta^{-1}");
  if(m_massSquaredMode){
    m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
    m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
  }
  zoomHisto->SetStats(kTRUE);
  //cout << "mean: " << mean << "  stdDev: " << stdDev << "  mu: " << mu << "   sig: " << sig << endl;
  //fitFunction->SetParameter(0,(counts)*zoomHistoTemp->GetBinWidth(2));

  TF1* initialFunct = NULL;
  if(m_useStudentTDistributionsForTOF){

    double sigParGuess = sqrt(stdDev*stdDev*(ndfParamGuess - 2.0)/ndfParamGuess);

    m_student_T_Funct->SetRange(lowRange,highRange);
    m_student_T_Funct->SetParameter(0,integral);
    m_student_T_Funct->SetParLimits(0,0.1*integral, 10.0*integral);
    m_student_T_Funct->SetParameter(1,mean);
    m_student_T_Funct->SetParLimits(1, mean - 0.25, mean + 0.25);
    m_student_T_Funct->SetParameter(2,sigParGuess);
    m_student_T_Funct->SetParLimits(2,0.01*sigParGuess,100.0*sigParGuess);
    m_student_T_Funct->SetParameter(3,ndfParamGuess);
    m_student_T_Funct->SetParLimits(3,1.0,50.0);
    m_student_T_Funct->SetParNames("N","#mu","#sigma","#nu");
    initialFunct = (TF1*) m_student_T_Funct->Clone();
    initialFunct->SetLineColor(kViolet);
  }else{
    m_gausFunct->SetRange(lowRange,highRange);
    m_gausFunct->SetParameter(0,integral);
    m_gausFunct->SetParLimits(0,0.1*integral, 5.0*integral);
    m_gausFunct->SetParameter(1,mean);
    m_gausFunct->SetParLimits(1, mean - 0.25, mean + 0.25);
    m_gausFunct->SetParameter(2,stdDev);
    m_gausFunct->SetParLimits(2,0.001*stdDev,2.0*stdDev);
    m_gausFunct->SetParNames("N","#mu","#sigma");
    initialFunct = (TF1*) m_gausFunct->Clone();
    initialFunct->SetLineColor(kViolet);
  }

  zoomHisto->SetStats(kTRUE);
  gStyle->SetOptFit(0111);
  //m_logNormFunct->Draw("same");


  //FITTING HERE
  int status = -1;
  string fitOptions = "EQ";
  if(m_doIntegralOfBinForFits) fitOptions += "I";
  m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(zoomHisto,-1,true);
  if(m_useStudentTDistributionsForTOF)  status = zoomHisto->Fit("m_student_T_Funct",fitOptions.c_str());
  else status = zoomHisto->Fit("m_gausFunct",fitOptions.c_str());

  returnParams[8] = status;
  if(status == 0){
    TF1* fitResult = NULL;
    if(m_useStudentTDistributionsForTOF)  fitResult = zoomHisto->GetFunction("m_student_T_Funct");
    else fitResult = zoomHisto->GetFunction("m_gausFunct");
    if(fitResult->GetParameter(0) != 0.0){
      returnParams[0] = fitResult->GetParameter(0);
      returnParams[1] = fitResult->GetParameter(1);
      returnParams[2] = fitResult->GetParameter(2);
      returnParams[4] = fitResult->GetParError(0);
      returnParams[5] = fitResult->GetParError(1);
      returnParams[6] = fitResult->GetParError(2);
      if(m_useStudentTDistributionsForTOF){
        returnParams[3] = fitResult->GetParameter(3);
        returnParams[7] = fitResult->GetParError(3);
      }
      cout << "Zoom Params: " << returnParams[0] << " " <<returnParams[1] << " " << returnParams[2];
      if(m_useStudentTDistributionsForTOF) cout << " " << returnParams[3];
      cout << endl;
    }
    //draw result
    fittingCanvas->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits) initialFunct->Draw("same");
    fitResult->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();
    cout << "     ";
    fittingCanvas->SaveAs(Form("%s/%s/%s/BTOF_Zoomed/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
     m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),imageName.Data()));
  }else{
    cout << "Bad Fit for zoom " << m_partInfo->GetParticleName(a_particleIndex,0) << " Eta Bin " << m_currentRapBin << " pT Bin " << m_currentMtM0Bin << endl;
    fittingCanvas->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits) initialFunct->Draw("same");
    if(m_useStudentTDistributionsForTOF) m_student_T_Funct->Draw("same");
    else m_gausFunct->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();
    fittingCanvas->SaveAs(Form("%s/%s/%s/BTOF_Zoomed/BAD_%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),imageName.Data()));
  }


  delete zoomHisto;
  delete initialFunct;
  delete fittingCanvas;
  return returnParams;

}



//==================================================================================================================================

vector<double> ZFitter::fit_BTOF_withBackground(int a_particleIndex, vector<double> a_fixedParameters, vector<double> a_momentum_startUsingBkgd, bool a_useStudentT, bool a_useExpBackground, bool a_fitBackgroundOnRight){
  cout << "   in  ZFitter::fit_BTOF_withBackground(vector<double> a_fixedParameters)" << endl;
  double ndfParamGuess = 8.0; // >2 : input to student t, number of degrees of freedom
  double sigmaGuess = 0.01;
  double fracOfWindowForBackground = 0.2; //Using this fraction of histogram range from left to seed the background parameters

  //Setup Return Structure
  vector<double> returnParams;
  returnParams.push_back(-999); // parameters
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999); //parameter errors
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999); //status
  returnParams.push_back(-999); //chi^2


  TCanvas* fittingCanvas = new TCanvas("fittingCanvas_bkg_tof","Fitting Canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);


  TH1D* zoomHisto = m_currentHistoToFit;
  double lowRange = HistogramUtilities::getLowEdge(zoomHisto);
  double highRange = HistogramUtilities::getHighEdge(zoomHisto);
  /*
  if(m_massSquaredMode){
    lowRange = pow(m_partInfo->GetParticleMass(a_particleIndex)-0.1,2.0);
    highRange = pow(m_partInfo->GetParticleMass(a_particleIndex)+0.2,2.0);
  }*/
  //########   Create the Fit Function  #############
  TF1* fit_funct = NULL;
  if(a_useExpBackground){
    if(a_useStudentT){
      fit_funct = new TF1("StudentT_withPowerLawFunction",
      "[0]*TMath::Student((x - [1])/[2], [3])/[2]   + [6] +  [4]*exp(-[5]*(x - [7])) + [8]",lowRange,highRange);
      fit_funct->SetNpx(500);
      fit_funct->SetLineColor(kBlue);
      fit_funct->SetParNames("N","#mu","#sigma","#nu","A","b","ignore","x_{0}","y_{0}");
      fit_funct->FixParameter(6,0.0);

    }else{
      fit_funct = new TF1("Gaus_withPowerLawFunction",
      "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]   +[3] +[6] +  [4]*exp(-[5]*(x - [7])) + [8]",lowRange,highRange);
      fit_funct->SetNpx(500);
      fit_funct->SetLineColor(kBlue);
      fit_funct->SetParNames("N","#mu","#sigma","ignore","A","b","ignore2","x_{0}","y_{0}");
      fit_funct->FixParameter(3,0.0);
      fit_funct->FixParameter(6,0.0);
    }
  }else{
    if(a_useStudentT){
      fit_funct = new TF1("StudentT_withPowerLawFunction",
      "[0]*TMath::Student((x - [1])/[2], [3])/[2]   +  [4]*exp(-[5]*exp([6]*(x - [7]))) + [8]",lowRange,highRange);
      fit_funct->SetNpx(500);
      fit_funct->SetLineColor(kBlue);
      fit_funct->SetParNames("N","#mu","#sigma","#nu","A","b","b2","x_{0}","y_{0}");
    }else{
      fit_funct = new TF1("Gaus_withPowerLawFunction",
      "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]   +[3]+  [4]*exp(-[5]*exp([6]*(x - [7]))) + [8]",lowRange,highRange);
      fit_funct->SetNpx(500);
      fit_funct->SetLineColor(kBlue);
      fit_funct->SetParNames("N","#mu","#sigma","ignore","A","b","b2","x_{0}","y_{0}");
      fit_funct->FixParameter(3,0.0);
    }
  }






  //#############   Seed Background Parameters   #################
  bool ingoreBackground = false;

  double aAvg = 0.0;  // Average Amplitude Background Parameter
  double bAvg = 0.0;
  double bAvgHigh = 0.0;
  int bin1 = HistogramUtilities::firstBinWithContent(zoomHisto);
  double firstPoint = zoomHisto->GetBinLowEdge(bin1);
  double lastPoint = HistogramUtilities::getLowEdge(zoomHisto) + fracOfWindowForBackground * (HistogramUtilities::getHighEdge(zoomHisto) - HistogramUtilities::getLowEdge(zoomHisto));
  if(m_massSquaredMode){
    firstPoint = lowRange;
    //lastPoint = zoomHisto->FindBin(highRange);
    lastPoint = pow(m_partInfo->GetParticleMass(a_particleIndex)-0.01,2.0);

    if(a_fitBackgroundOnRight){
      cout << "background on right" << endl;
      lastPoint = zoomHisto->GetBinLowEdge(zoomHisto->FindBin(highRange)-1);
    }

  }
  int lastBinForBackground = zoomHisto->FindBin(lastPoint);
  double midpoint = (firstPoint + lastPoint)/2.0;
  if(a_fitBackgroundOnRight){
      cout << "midpoint on right" << endl;
      midpoint = pow(m_partInfo->GetParticleMass(a_particleIndex)+0.075,2.0);
  }
  int midBin = zoomHisto->FindBin(midpoint);



  cout << "Background points:" << endl;
  cout << "  lowPoint: " << firstPoint << " ( bin " << bin1 << ")" << endl;
  cout << "  midPoint: " << midpoint << " ( bin " << midBin << ")" << endl;
  cout << "  highPoint: " << lastPoint << " ( bin " << lastBinForBackground << ")" << endl;

  if(midBin <= bin1 || midBin >= lastBinForBackground ){
    cout << " ERROR: lastBinForBackground: " << lastBinForBackground << "    is too low. Point: " << lastPoint << "  with low limit " << HistogramUtilities::getLowEdge(zoomHisto) << endl;
  }



  if(a_useExpBackground){
    if(a_fixedParameters.size() > 5  && a_fixedParameters[5] > 0.0){
      bAvg = a_fixedParameters[5];

      double x1_avg = 0.5*(firstPoint + midpoint);
      double y1_avg = 0.0;
      double numBinsAdded = 0.0;
      for(int bin = bin1; bin < midBin; bin++){
        y1_avg += log(zoomHisto->GetBinContent(bin));
        numBinsAdded+= 1.0;
      }
      y1_avg = exp(y1_avg/numBinsAdded);

      aAvg = y1_avg*exp(bAvg*x1_avg);

    }else if(a_fitBackgroundOnRight){

      cout << "fitting on right" << endl;

      double x1_avg = 0.5*(firstPoint + midpoint);
      cout << "x1_avg: " << x1_avg << endl;
      double x2_avg = 0.5*(lastPoint + midpoint);
      cout << "x2_avg: " << x2_avg << endl;

      cout << "low content: " << zoomHisto->GetBinContent(bin1) << "   log: " << log(zoomHisto->GetBinContent(bin1)) << endl;
      cout << "mid content: " << zoomHisto->GetBinContent(midBin) << "   log: " << log(zoomHisto->GetBinContent(midBin)) << endl;
      cout << "high content: " << zoomHisto->GetBinContent(lastBinForBackground) << "   log: " << log(zoomHisto->GetBinContent(lastBinForBackground)) << endl;

      double y1_avg = exp(0.5*(log(zoomHisto->GetBinContent(bin1)) + log(zoomHisto->GetBinContent(midBin))));
      double y2_avg = exp(0.5*(log(zoomHisto->GetBinContent(midBin)) + log(zoomHisto->GetBinContent(lastBinForBackground))));

      cout << "y1_avg: " << y1_avg << endl;
      cout << "y2_avg: " << y2_avg << endl;

      bAvg = log(y1_avg/y2_avg)/(x2_avg - x1_avg);
      aAvg = 0.5*(y1_avg*exp(bAvg * x1_avg) + y2_avg*exp(bAvg * x2_avg));


      cout << "bAvg: " << bAvg << endl;
      cout << "aAvg: " << aAvg << endl;

    }else{


      // divide the lower region where the background is, into 2 parts to seed background parameters
      double x1_avg = 0.5*(firstPoint + midpoint);
      cout << "x1_avg: " << x1_avg << endl;
      double x2_avg = 0.5*(lastPoint + midpoint);
      cout << "x2_avg: " << x2_avg << endl;


      double y1_avg = 0.0;
      double numBinsAdded = 0.0;
      for(int bin = bin1; bin < midBin; bin++){
        y1_avg += log(zoomHisto->GetBinContent(bin));
        numBinsAdded+= 1.0;
      }
      y1_avg = exp(y1_avg/numBinsAdded);

      cout << "y1_avg: " << y1_avg << endl;

      double y2_avg = 0.0;
      numBinsAdded = 0.0;
      for(int bin = midBin; bin <= lastBinForBackground; bin++){
        if(zoomHisto->GetBinContent(bin) <= 0) continue;
        cout << "bin " << bin << " content: " << zoomHisto->GetBinContent(bin) << endl;
        cout << "log(" << zoomHisto->GetBinContent(bin) << "): " << log(zoomHisto->GetBinContent(bin)) << endl;
        y2_avg += log(zoomHisto->GetBinContent(bin));
        cout << "current y2_avg: " << y2_avg << endl;
        numBinsAdded+= 1.0;
      }
      y2_avg = exp(y2_avg/numBinsAdded);

      cout << "y2_avg: " << y2_avg << endl;

      bAvg = log(y1_avg/y2_avg)/(x2_avg - x1_avg);
      aAvg = y1_avg*exp(bAvg * x1_avg);


      cout << "bAvg: " << bAvg << endl;
      cout << "aAvg: " << aAvg << endl;


    /*double weights = 0.0;
    for(int bin = bin1; bin < lastBinForBackground; bin++){
      double x1 = zoomHisto->GetBinCenter(bin);
      double y1 = zoomHisto->GetBinContent(bin);
      double y1_err = zoomHisto->GetBinError(bin);
      double x2 = zoomHisto->GetBinCenter(bin+1);
      double y2 = zoomHisto->GetBinContent(bin+1);
      double y2_err = zoomHisto->GetBinError(bin+1);
      if(y1 <= 0.0 || y2 <= 0.0){
        continue;
      }

      double bCalc = log(y1/y2)/(x2-x1);
      double weight = (y1_err*y1_err/(y1*y1) + y2_err*y2_err/(y2*y2)) / ((x2-x1) * (x2-x1)); // sqr of error
      if(!std::isnan(weight) && !std::isnan(bCalc)){
        weights += weight;
        bAvg += bCalc*weight;
      }


      #ifdef _ZFITTER_DEBUG_
      cout << " (x1,y1) = (" << x1 << ", "  << y1  << ")  (x2,y2) = (" << x2 << ", "  << y2 << ")   bCalc = " << bCalc << "  weight= " << weight << endl;
      #endif
    }
    bAvg = bAvg/weights;


    weights = 0.0;
    for(int bin = bin1; bin < lastBinForBackground; bin++){
      double x1 = zoomHisto->GetBinCenter(bin);
      double y1 = zoomHisto->GetBinContent(bin);
      double y1_err = zoomHisto->GetBinError(bin);
      if(y1 <= 0.0){
        continue;
      }
      double weight = y1_err*y1_err*exp(2.0*bAvg*x1);
      aAvg += y1*exp(bAvg*x1)*weight;
      weights += weight;
    }
    aAvg = aAvg/weights;
    */


    }
  }else{

  }







  if(m_currentMomOfBin < a_momentum_startUsingBkgd.at(a_particleIndex)){
      fit_funct->FixParameter(4,0.0); // A
      fit_funct->FixParameter(5,0.0); // b
      fit_funct->FixParameter(6,0.0);
      fit_funct->FixParameter(7,0.0); // x_0
      fit_funct->FixParameter(8,0.0); // y_0
  }else if(bAvg < 0.0){ // || numBinsAdded < 2.0){ // retry background seed
    cout << " Error: b parameter is negative, should be positive, trying a different seed" << endl;

    double lowBinsCenter = zoomHisto->GetBinLowEdge(3);
    double avgLowBins = ( zoomHisto->GetBinContent(1) + zoomHisto->GetBinContent(2)
                        + zoomHisto->GetBinContent(3) + zoomHisto->GetBinContent(4) )/4.0;
    int lastContentBin = HistogramUtilities::lastBinWithContent(zoomHisto);
    double highBinsCenter = zoomHisto->GetBinLowEdge(lastContentBin-1);
    double avgHighBins = ( zoomHisto->GetBinContent(lastContentBin-3) + zoomHisto->GetBinContent(lastContentBin-2)
                         + zoomHisto->GetBinContent(lastContentBin-1) + zoomHisto->GetBinContent(lastContentBin) )/4.0;
    bAvg = log(avgLowBins/avgHighBins)/(highBinsCenter - lowBinsCenter);


    if(bAvg < 0.0 || aAvg < 0.0){
      cout << "  Error: Seed for background exponent failed : bAvg=" << bAvg << endl;
      cout << " Background will be silenced for this bin" << endl;
      fit_funct->FixParameter(4,0.0); // A
      fit_funct->FixParameter(5,1.0); // b
      if(!a_useExpBackground) fit_funct->FixParameter(6,1.0);
      fit_funct->FixParameter(7,0.0); // x_0
      fit_funct->FixParameter(8,0.0); // y_0
      ingoreBackground = true;
    }else{

      aAvg = avgLowBins*exp(bAvg*lowBinsCenter);
      cout << "  Seeding Second Try: A= " << aAvg << "  b:" << bAvg << endl;
      fit_funct->SetParameter(4,aAvg); // Amp
      fit_funct->SetParLimits(4,0.0,1.0E7 * aAvg);
      fit_funct->SetParameter(5,bAvg); // Exponential factor
//      fit_funct->SetParLimits(5,1E-4*bAvg,5*bAvg);
      fit_funct->SetParLimits(5,1E-4*bAvg,10);
      if(!a_useExpBackground) fit_funct->SetParameter(6,1.0);
      fit_funct->FixParameter(7,0.0); // x_0
      fit_funct->FixParameter(8,0.0); // y_0
    }
  }else{ // good background seed
    fit_funct->SetParameter(4,aAvg); // Amp
    fit_funct->SetParLimits(4,0.0,1.0E7 * aAvg);
    fit_funct->SetParameter(5,bAvg); // Exponential factor
//    fit_funct->SetParLimits(5,1E-4*bAvg,1E4*bAvg);
    fit_funct->SetParLimits(5,1E-4*bAvg,10);
    if(!a_useExpBackground) fit_funct->SetParameter(6,1.0);
    fit_funct->FixParameter(7,0.0); // x_0
    fit_funct->FixParameter(8,0.0); // y_0

    cout << "     Background Limits:  A = " << aAvg << " on [" << 0.0 << ", " << 1.0E7 * aAvg << "]" << endl;
    cout << "                         b = " << bAvg << " on [" << 1E-4*bAvg << ", " << 1E4*bAvg << "]" << endl;

  }




  //###########     Seed the Source Peak Fit Parameters  ###########
  double integral = HistogramUtilities::sumOfContent(zoomHisto, -1, true);
  double mean = 0.0; // center of peak / only for m_currentPartIndex
  if(m_massSquaredMode) mean = pow(m_partInfo->GetParticleMass(a_particleIndex),2.0 );
  double stdDev = zoomHisto->GetStdDev()*0.8;
  if(m_massSquaredMode) stdDev = 0.05;

  /*double gausData[8];
  double lowTightCut = HistogramUtilities::getLowEdge(zoomHisto) + 0.25*(HistogramUtilities::getHighEdge(zoomHisto) - HistogramUtilities::getLowEdge(zoomHisto));
  double highTightCut = HistogramUtilities::getHighEdge(zoomHisto) - 0.25*(HistogramUtilities::getHighEdge(zoomHisto) - HistogramUtilities::getLowEdge(zoomHisto));
  HistogramUtilities::getGausData(&gausData[0], zoomHisto, true, true, lowTightCut , highTightCut);

  cout << " gaus data : " << gausData[0] << " : " << gausData[1] << " : " << gausData[2] << endl;
  */

  double stdDevCalc = 0.0;
  double newAmp = 0.0;
  double sigmaParCalc = 0.0;
  if(a_fixedParameters.size() > 2 && a_fixedParameters[2] > 0.0){
    cout << " Using Fixed Parameters to Seed Amplitude" << endl;
    if(a_useStudentT){
      newAmp = zoomHisto->GetBinContent(zoomHisto->FindBin(a_fixedParameters[1]))   * a_fixedParameters[2] / TMath::Student(0.0, a_fixedParameters[3]);
    }else{
      newAmp = zoomHisto->GetBinContent(zoomHisto->FindBin(a_fixedParameters[1])) / gaus(a_fixedParameters[1], a_fixedParameters[2], a_fixedParameters[1]);
    }

    mean = a_fixedParameters[1];
    stdDevCalc = 1.0;
    sigmaParCalc = a_fixedParameters[2];
  }else{



    //vector< double > calcData =  HistogramUtilities::calcGausData(zoomHisto, 0.0, 12, 3);
    //newAmp = calcData[0];
    //mean = calcData[1];
    //sigmaParCalc = calcData[2];


    int peakBin = zoomHisto->FindBin(mean);
    cout << "peakBin: " << peakBin << " value: " << mean << endl;
    //double width = HistogramUtilities::getHighEdge(zoomHisto) - zoomHisto->GetBinLowEdge(1);
    double width = stdDev;
    double x3 = zoomHisto->GetBinCenter(peakBin);  // point at peak
    int testBin1 = zoomHisto->FindBin(x3 - 0.08*width);
    int testBin2 = zoomHisto->FindBin(x3 + 0.12*width);
    if(testBin1 >= peakBin) testBin1 = peakBin - 1;
    if(testBin2 <= peakBin) testBin2 = peakBin + 1;
    //int testBin1 = peakBin - 6;
    //int testBin2 = peakBin + 6;

    double x4 = zoomHisto->GetBinCenter(testBin1);
    double x5 = zoomHisto->GetBinCenter(testBin2);
    double y3 = HistogramUtilities::getBinContentAvg(zoomHisto, peakBin, 2);  // point at peak
    double y4 = HistogramUtilities::getBinContentAvg(zoomHisto, testBin1, 2);
    double y5 = HistogramUtilities::getBinContentAvg(zoomHisto, testBin2, 2);

    if(!ingoreBackground){ // subtract off the background
      double new_y3 = y3 - aAvg*exp(-bAvg*x3);
      double new_y4 = y4 - aAvg*exp(-bAvg*x4);
      double new_y5 = y5 - aAvg*exp(-bAvg*x5);

      if( new_y3 > 0.0 && new_y4 > 0.0 && new_y5 > 0.0){
        y3 = new_y3;
        y4 = new_y4;
        y5 = new_y5;
      }
    }

    cout << "Points for Gaussin seed: " <<  " p1 " << x3 << "  " << y3 << "    p2 " << x4 << "  " << y4 << "     p3 " << x5 << "  " << y5 << endl;

    PhysMath::get3PointGausParams(newAmp, mean, sigmaParCalc,x3,y3,x4,y4,x5,y5);
    if(a_useStudentT){ //convert stdDev to sigma of studentT
      sigmaParCalc = sqrt(sigmaParCalc*sigmaParCalc*(ndfParamGuess - 2.0)/ndfParamGuess);
    }

    if(m_massSquaredMode){
      mean = pow(m_partInfo->GetParticleMass(a_particleIndex),2.0 );
      if(fabs(sigmaParCalc) > 0.2) sigmaParCalc = stdDev;
      if(a_particleIndex == 1 && fabs(sigmaParCalc) > 0.1) sigmaParCalc = stdDev;
      if(newAmp < 0) newAmp = HistogramUtilities::getBinContentAvg(zoomHisto, peakBin, 2)/(sigmaParCalc*2.5066283);
   }


    /*
    int peakBin = zoomHisto->FindBin(mean);
    double x3 = zoomHisto->GetBinCenter(peakBin);  // point at peak
    double y3 = zoomHisto->GetBinContent(peakBin);
    int testBin = zoomHisto->FindBin(x3 + 1.5*gausData[2]);
    if(testBin > zoomHisto->GetNbinsX()){
      cout << "   Test bin is off the edge of the histogram for sigma seed, using last bin" << endl;
      testBin = zoomHisto->GetNbinsX();
    }
    double x4 = zoomHisto->GetBinCenter(testBin); //point 1 stdev right of peak
    double y4 = zoomHisto->GetBinContent(testBin);

    if(!ingoreBackground){ // subtract off the background
      y3 = y3 - aAvg*exp(-bAvg*x3);
      y4 = y4 - aAvg*exp(-bAvg*x4);
    }

    stdDevCalc = sqrt(0.5*(x4-x3)*(x4-x3)/log(y3/y4));
    cout << " first internal=" << 0.5*(x4-x3)*(x4-x3)/log(y3/y4) << endl;
    cout << "   stdDevCalc= " << stdDevCalc << endl;
    if(y3 == 0.0 || y4==0.0 || y3 == y4 || y3/y4 < 0.0){
      cout << "    stdDevCalc = " << stdDevCalc << "   from (x3,y3)=(" << x3 << ", " << y3 << ")   and (x4,y4)=(" << x4 << ", " << y4 << ")" << endl;
      cout << "      setting it to guess " << endl;
      stdDevCalc = sigmaGuess;
    }
    sigmaParCalc = sqrt(stdDevCalc*stdDevCalc*(ndfParamGuess - 2.0)/ndfParamGuess);
    cout << "      internal calc=" << stdDevCalc*stdDevCalc*(ndfParamGuess - 2.0)/ndfParamGuess << endl;
    cout << "   sigmaParCalc= " << sigmaParCalc << endl;

    newAmp = y3*sigmaParCalc/ TMath::Student(0.0, ndfParamGuess);
    */
  }



  cout << " Starting Parameters:   N:" << newAmp << "  mu:" << mean << "   sig:" << sigmaParCalc;
  if(a_useStudentT) cout << "   nu:" << ndfParamGuess << endl;
  cout << "                        A:" << aAvg << "  b:" << bAvg << "  x_0: 0.0   y_0:0.0" << endl;


  if(newAmp == 0 || std::isnan(mean) || std::isnan(sigmaParCalc)){
    cout << " Error in ZFitter::fitIsolation_BTOF_withBackground: integral = "<< newAmp << " mean = "<< mean << "  stdDev = "<< sigmaParCalc << endl;
    //delete fittingCanvas;
    return returnParams;
  }

  TString imageName = TString::Format("%s_withBackground_rap_%02d_mTm0_%02d",
      m_partInfo->GetParticleName(a_particleIndex,0).Data(),m_currentRapBin,m_currentMtM0Bin);
  zoomHisto->SetName(imageName.Data());
  zoomHisto->SetTitle(Form("BTOF Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] (GeV/c) p=%2.3f : zoom on %s",
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
    m_currentMomOfBin,
    m_partInfo->GetParticleName(m_currentPartIndex,0).Data()
    ));
  if(m_massSquaredMode){
    zoomHisto->GetXaxis()->SetTitle("m^{2}");
    zoomHisto->GetYaxis()->SetTitle("dN/d(m^{2})}");
  }else{
    zoomHisto->GetXaxis()->SetTitle("#beta^{-1}");
    zoomHisto->GetYaxis()->SetTitle("dN/d#beta^{-1}");
  }

  zoomHisto->SetStats(kTRUE);
  gStyle->SetOptFit(0111);




  //##########  Setup Initial Seed Function to Draw for Debugging  #############
  fit_funct->SetParameter(0,newAmp);
  fit_funct->SetParLimits(0,0.01*newAmp, 100.0*newAmp);
  fit_funct->SetParameter(1,mean);
  fit_funct->SetParLimits(1, mean - 0.15, mean + 0.15);
  if(m_massSquaredMode)fit_funct->SetParLimits(1, mean - 0.01, mean + 0.01);
  fit_funct->SetParameter(2,sigmaParCalc);
  fit_funct->SetParLimits(2,0.1*sigmaParCalc,5*sigmaParCalc);
  if(a_useStudentT){
    fit_funct->SetParameter(3,ndfParamGuess);
    fit_funct->SetParLimits(3,2.0,1000.0);
  }

 //Fix the parameters in argument
  if(a_fixedParameters.size() > 0){
    for(int jjj = 0; jjj < 9; jjj++){
      if(a_fixedParameters[jjj] != -999 && !std::isnan(a_fixedParameters[jjj]) && !std::isinf(a_fixedParameters[jjj])){
        fit_funct->FixParameter(jjj, a_fixedParameters[jjj]);
        cout << "    Fixing Parameter " << jjj << " : " << a_fixedParameters[jjj] << endl;
      }
    }
  }

  TF1* initialFunct = NULL;
  initialFunct = (TF1*) fit_funct->Clone();
  initialFunct->SetLineColor(kYellow);

  // Create
  TF1* source_funct_init = NULL;
  if(a_useStudentT){
    source_funct_init = new TF1("StudentT_Source_Init", "[0]*TMath::Student((x - [1])/[2], [3])/[2]",lowRange,highRange);
    source_funct_init->SetNpx(200);
    source_funct_init->SetLineColor(kYellow+1);
    source_funct_init->SetLineStyle(2);
    source_funct_init->SetParNames("N","#mu","#sigma","#nu");
    source_funct_init->SetParameter(0,fit_funct->GetParameter(0));
    source_funct_init->SetParameter(1,fit_funct->GetParameter(1));
    source_funct_init->SetParameter(2,fit_funct->GetParameter(2));
    source_funct_init->SetParameter(3,fit_funct->GetParameter(3));
  }else{
    source_funct_init = new TF1("gaus_Source_Init", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
    source_funct_init->SetNpx(200);
    source_funct_init->SetLineColor(kYellow+1);
    source_funct_init->SetLineStyle(2);
    source_funct_init->SetParNames("N","#mu","#sigma");
    source_funct_init->SetParameter(0,fit_funct->GetParameter(0));
    source_funct_init->SetParameter(1,fit_funct->GetParameter(1));
    source_funct_init->SetParameter(2,fit_funct->GetParameter(2));
  }

  // Create Red Background Function For Drawing
  TF1* background_funct_init = new TF1("Background_PowerLawFunction_Init", "[0]*exp(-[1]*(x - [3])) + [4] + [2]",lowRange,highRange);
  if(!a_useExpBackground) background_funct_init = new TF1("Background_PowerLawPowerLawFunction_Init", "[0]*exp(-[1]*exp([2]*(x - [3]))) + [4]",lowRange,highRange);
  background_funct_init->SetNpx(200);
  background_funct_init->SetLineColor(kYellow+2);
  background_funct_init->SetLineStyle(2);
  background_funct_init->SetParNames("A","b","ignore","x_{0}","y_{0}");
  if(!a_useExpBackground)background_funct_init->SetParNames("A","b","b2","x_{0}","y_{0}");
  background_funct_init->SetParameter(0,fit_funct->GetParameter(4));
  background_funct_init->SetParameter(1,fit_funct->GetParameter(5));
  background_funct_init->SetParameter(2,fit_funct->GetParameter(6));
  background_funct_init->SetParameter(3,fit_funct->GetParameter(7));
  background_funct_init->SetParameter(4,fit_funct->GetParameter(8));


  //###############   FIT   THE   CURVE   #####################
  cout << "      Starting the fit in funct" << endl;
  int status = -999;
  if(a_useStudentT){
  //  status = zoomHisto->Fit("StudentT_withPowerLawFunction","ES");
    status = zoomHisto->Fit(fit_funct,"ES");
  }else{
    //status = zoomHisto->Fit("Gaus_withPowerLawFunction","ES");
    status = zoomHisto->Fit(fit_funct,"ES");
  }


  if(status == 0){
    TF1* fitResult = NULL;
      if(a_useStudentT){
        fitResult = zoomHisto->GetFunction("StudentT_withPowerLawFunction");
      }else{
        fitResult = zoomHisto->GetFunction("Gaus_withPowerLawFunction");
      }
    if(fitResult->GetParameter(0) != 0.0){
      returnParams[0] = fitResult->GetParameter(0);
      returnParams[1] = fitResult->GetParameter(1);
      returnParams[2] = fitResult->GetParameter(2);
      returnParams[3] = fitResult->GetParameter(3);
      returnParams[4] = fitResult->GetParameter(4);
      returnParams[5] = fitResult->GetParameter(5);
      returnParams[6] = fitResult->GetParameter(6);
      returnParams[7] = fitResult->GetParameter(7);
      returnParams[8] = fitResult->GetParameter(8);

      returnParams[9] =  fitResult->GetParError(0);
      returnParams[10] =  fitResult->GetParError(1);
      returnParams[11] = fitResult->GetParError(2);
      returnParams[12] = fitResult->GetParError(3);
      returnParams[13] = fitResult->GetParError(4);
      returnParams[14] = fitResult->GetParError(5);
      returnParams[15] = fitResult->GetParError(6);
      returnParams[16] = fitResult->GetParError(7);
      returnParams[17] = fitResult->GetParError(9);
      returnParams[18] = status;
      returnParams[19] = fitResult->GetChisquare() / fitResult->GetNDF();

    }else{
      cout << " Error in ZFitter::fitIsolation_BTOF_withBackground: fit status = "<< status << endl;
    }


    // Create Violet Source Peak For Drawinga
    TF1* source_funct = NULL;

    if(a_useStudentT){
      source_funct = new TF1("StudentT_Source", "[0]*TMath::Student((x - [1])/[2], [3])/[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma","#nu");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
      source_funct->SetParameter(3,returnParams[3]);
    }else{
      source_funct = new TF1("gaus_Source", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
    }


    // Create Red Background Function For Drawing
    TF1* background_funct = new TF1("Background_PowerLawFunction", "[0]*exp(-[1]*(x - [3])) + [2]+[4]",lowRange,highRange);
    if(!a_useExpBackground) background_funct = new TF1("Background_PowerLawPowerLawFunction", "[0]*exp(-[1]*pow([2]*(x - [3]))) + [4]",lowRange,highRange);
    background_funct->SetNpx(200);
    background_funct->SetLineColor(kRed);
    background_funct->SetLineStyle(2);
    background_funct->SetParNames("A","b","ignore","x_{0}","y_{0}");
    if(!a_useExpBackground)background_funct->SetParNames("A","b","b2","x_{0}","y_{0}");
    background_funct->SetParameter(0,returnParams[4]);
    background_funct->SetParameter(1,returnParams[5]);
    background_funct->SetParameter(2,returnParams[6]);
    background_funct->SetParameter(3,returnParams[7]);
    background_funct->SetParameter(4,returnParams[8]);



    //draw result
    fittingCanvas->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits){
      initialFunct->Draw("same");
      source_funct_init->Draw("same");
      background_funct_init->Draw("same");
    }
    fitResult->Draw("same");
    background_funct->Draw("same");
    source_funct->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();

    cout << "     ";
    if(m_inAllCentMode){
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),imageName.Data()));
      }
    }else{
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(), m_currentCentIndex,imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(), m_currentCentIndex,imageName.Data()));
      }
    }

    delete source_funct;
    delete background_funct;
  }else{  //////   BAD FIT / NOT OPTIMAL



    TF1* fitResult = zoomHisto->GetFunction("StudentT_withPowerLawFunction");
    cout << "Error in ZFitter::fit_BTOF_withBackground: fit status: " << status << "   " << m_partInfo->GetParticleName(m_currentPartIndex,0)
         << " Rap Bin " << m_currentRapBin << " mTm0 Bin " << m_currentMtM0Bin << endl;

    cout << "  Pushing back parameters anyway..." << endl;
    if(fitResult){
      cout << " Pushing back parameters in Fit Result : " << fitResult << endl;
      returnParams[0] = fitResult->GetParameter(0);
      returnParams[1] = fitResult->GetParameter(1);
      returnParams[2] = fitResult->GetParameter(2);
      returnParams[3] = fitResult->GetParameter(3);
      returnParams[4] = fitResult->GetParameter(4);
      returnParams[5] = fitResult->GetParameter(5);
      if(std::isnan(fit_funct->GetParameter(4)) || std::isinf(fit_funct->GetParameter(4)) || std::isnan(fit_funct->GetParameter(5)) || std::isinf(fit_funct->GetParameter(5))){
        returnParams[4] = 0.0;
        returnParams[5] = 0.0;
      }
      returnParams[6] = fitResult->GetParameter(6);
      returnParams[7] = fitResult->GetParameter(7);
      returnParams[8] = fitResult->GetParameter(8);

      returnParams[9] = fitResult->GetParError(0);
      returnParams[10] = fitResult->GetParError(1);
      returnParams[11] = fitResult->GetParError(2);
      returnParams[12] = fitResult->GetParError(3);
      returnParams[13] = fitResult->GetParError(4);
      returnParams[14] = fitResult->GetParError(5);
      returnParams[15] = fitResult->GetParError(6);
      returnParams[16] = fitResult->GetParError(7);
      returnParams[17] = fitResult->GetParError(8);
      returnParams[18] = status;
    }else{
      cout << " Pushing back parameters in Fit Funciton : " << fit_funct << endl;
      returnParams[0] = fit_funct->GetParameter(0);
      returnParams[1] = fit_funct->GetParameter(1);
      returnParams[2] = fit_funct->GetParameter(2);
      returnParams[3] = fit_funct->GetParameter(3);
      returnParams[4] = fit_funct->GetParameter(4);
      returnParams[5] = fit_funct->GetParameter(5);
      if(std::isnan(fit_funct->GetParameter(4)) || std::isinf(fit_funct->GetParameter(4)) || std::isnan(fit_funct->GetParameter(5)) || std::isinf(fit_funct->GetParameter(5))){
        returnParams[4] = 0.0;
        returnParams[5] = 0.0;
      }
      returnParams[6] = fit_funct->GetParameter(6);
      returnParams[7] = fit_funct->GetParameter(7);
      returnParams[8] = fit_funct->GetParameter(8);

      returnParams[9] = fit_funct->GetParError(0);
      returnParams[10] = fit_funct->GetParError(1);
      returnParams[11] = fit_funct->GetParError(2);
      returnParams[12] = fit_funct->GetParError(3);
      returnParams[13] = fit_funct->GetParError(4);
      returnParams[14] = fit_funct->GetParError(5);
      returnParams[15] = fit_funct->GetParError(6);
      returnParams[16] = fit_funct->GetParError(7);
      returnParams[17] = fit_funct->GetParError(8);
      returnParams[18] = status;
    }

    // Create Violet Source Peak For Drawinga
    TF1* source_funct = NULL;

    if(a_useStudentT){
      source_funct = new TF1("StudentT_Source", "[0]*TMath::Student((x - [1])/[2], [3])/[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma","#nu");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
      source_funct->SetParameter(3,returnParams[3]);
    }else{
      source_funct = new TF1("gaus_Source", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
    }


    // Create Red Background Function For Drawing
    TF1* background_funct = new TF1("Background_PowerLawFunction", "[0]*exp(-[1]*(x - [3])) +[2]+ [4]",lowRange,highRange);
    if(!a_useExpBackground) background_funct = new TF1("Background_PowerLawPowerLawFunction", "[0]*exp(-[1]*([2]*(x - [3]))) + [4]",lowRange,highRange);
    background_funct->SetNpx(200);
    background_funct->SetLineColor(kRed);
    background_funct->SetLineStyle(2);
    background_funct->SetParNames("A","b","ignore","x_{0}","y_{0}");
    if(!a_useExpBackground) background_funct->SetParNames("A","b","b2","x_{0}","y_{0}");
    background_funct->SetParameter(0,returnParams[4]);
    background_funct->SetParameter(1,returnParams[5]);
    background_funct->SetParameter(2,returnParams[6]);
    background_funct->SetParameter(3,returnParams[7]);
    background_funct->SetParameter(4,returnParams[8]);






    cout << "Drawing Now" << endl;
    fittingCanvas->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    if(zoomHisto->GetMinimum() < 0){
      zoomHisto->SetMinimum(1e-3);
      if(zoomHisto->GetMaximum() < 0){
        zoomHisto->SetMaximum(3e6);
      }
    }
    zoomHisto->Draw("E");
    initialFunct->Draw("same");
    source_funct_init->Draw("same");
    background_funct_init->Draw("same");
    if(returnParams[0] < 0.0 || std::isnan(returnParams[1]) || returnParams[2] <= 0.0 || (a_useStudentT && returnParams[3] <= 1.0)
        || returnParams[4] < 0.0 || returnParams[5] < 0.0){
      cout << " Really bad fit, not drawing any funcitons" << endl;
      //return returnParams;
    }else{
      background_funct->Draw("same");
      source_funct->Draw("same");
    }

    //fitResult->SetLineColor(kOrange);
    //fitResult->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();

    cout << "     ";
    if(m_inAllCentMode){
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),imageName.Data()));
      }
    }else{
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(), m_currentCentIndex,imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(), m_currentCentIndex,imageName.Data()));
      }
    }

    delete source_funct;
    delete background_funct;

  }


  //delete zoomHisto;
  delete initialFunct;
  delete background_funct_init;
  delete source_funct_init;
  delete fittingCanvas;
  //delete zoomHisto;


  cout << "  Finished fit_BTOF_withBackground" << endl;


  return returnParams;

}
//==================================================================================================================================

vector<double> ZFitter::fitIsolation_BTOF_withBackground(vector<double> a_fixedParameters, bool a_useStudentT){
  cout << "   in  ZFitter::fitIsolation_BTOF_withBackground(vector<double> a_fixedParameters)" << endl;
  double ndfParamGuess = 8.0; // >2 : input to student t, number of degrees of freedom
  double sigmaGuess = 0.01;
  double fracOfWindowForBackground = 0.2; //Using this fraction of histogram range from left to seed the background parameters

  //Setup Return Structure
  vector<double> returnParams;
  returnParams.push_back(-999); // parameters
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999); //parameter errors
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999); //status
  returnParams.push_back(-999); //chi^2


  TCanvas* fittingCanvas = new TCanvas("fittingCanvas_iso_tof","Fitting Canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);


  TH1D* zoomHisto = m_currentHistoToFit;
  double lowRange = HistogramUtilities::getLowEdge(zoomHisto);
  double highRange = HistogramUtilities::getHighEdge(zoomHisto);

  //########   Create the Fit Function  #############
  TF1* fit_funct = NULL;
  if(a_useStudentT){
    fit_funct = new TF1("StudentT_withPowerLawFunction",
      "[0]*TMath::Student((x - [1])/[2], [3])/[2]   +  [4]*exp(-[5]*(x - [6])) + [7]",lowRange,highRange);
    fit_funct->SetNpx(500);
    fit_funct->SetLineColor(kBlue);
    fit_funct->SetParNames("N","#mu","#sigma","#nu","A","b","x_{0}","y_{0}");
  }else{
    fit_funct = new TF1("Gaus_withPowerLawFunction",
      "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]   +  [4]*exp(-[5]*(x - [6])) + [7]",lowRange,highRange);
    fit_funct->SetNpx(500);
    fit_funct->SetLineColor(kBlue);
    fit_funct->SetParNames("N","#mu","#sigma","ignore","A","b","x_{0}","y_{0}");
    fit_funct->FixParameter(3,0.0);
  }






  //#############   Seed Background Parameters   #################
  bool ingoreBackground = false;

  double aAvg = 0.0;  // Average Amplitude Background Parameter
  double bAvg = 0.0;
  int bin1 = HistogramUtilities::firstBinWithContent(zoomHisto);
  double firstPoint = zoomHisto->GetBinLowEdge(bin1);
  double lastPoint = HistogramUtilities::getLowEdge(zoomHisto) + fracOfWindowForBackground * (HistogramUtilities::getHighEdge(zoomHisto) - HistogramUtilities::getLowEdge(zoomHisto));
  int lastBinForBackground = zoomHisto->FindBin(lastPoint);
  double midpoint = (firstPoint + lastPoint)/2.0;
  int midBin = zoomHisto->FindBin(midpoint);

  if(midBin <= bin1 || midBin >= lastBinForBackground ){
    cout << " ERROR: lastBinForBackground: " << lastBinForBackground << "    is too low. Point: " << lastPoint << "  with low limit " << HistogramUtilities::getLowEdge(zoomHisto) << endl;
  }




  if(a_fixedParameters.size() > 5  && a_fixedParameters[5] > 0.0){
    bAvg = a_fixedParameters[5];

    double x1_avg = 0.5*(firstPoint + midpoint);
    double y1_avg = 0.0;
    double numBinsAdded = 0.0;
    for(int bin = bin1; bin < midBin; bin++){
      y1_avg += log(zoomHisto->GetBinContent(bin));
      numBinsAdded+= 1.0;
    }
    y1_avg = exp(y1_avg/numBinsAdded);

    aAvg = y1_avg*exp(bAvg*x1_avg);

  }else{


    // divide the lower region where the background is, into 2 parts to seed background parameters
    double x1_avg = 0.5*(firstPoint + midpoint);
    double x2_avg = 0.5*(lastPoint + midpoint);


    double y1_avg = 0.0;
    double numBinsAdded = 0.0;
    for(int bin = bin1; bin < midBin; bin++){
      y1_avg += log(zoomHisto->GetBinContent(bin));
      numBinsAdded+= 1.0;
    }
    y1_avg = exp(y1_avg/numBinsAdded);

    double y2_avg = 0.0;
    numBinsAdded = 0.0;
    for(int bin = midBin; bin <= lastBinForBackground; bin++){
      y2_avg += log(zoomHisto->GetBinContent(bin));
      numBinsAdded+= 1.0;
    }
    y2_avg = exp(y2_avg/numBinsAdded);

    bAvg = log(y1_avg/y2_avg)/(x2_avg - x1_avg);
    aAvg = y1_avg*exp(bAvg * x1_avg);



    /*double weights = 0.0;
    for(int bin = bin1; bin < lastBinForBackground; bin++){
      double x1 = zoomHisto->GetBinCenter(bin);
      double y1 = zoomHisto->GetBinContent(bin);
      double y1_err = zoomHisto->GetBinError(bin);
      double x2 = zoomHisto->GetBinCenter(bin+1);
      double y2 = zoomHisto->GetBinContent(bin+1);
      double y2_err = zoomHisto->GetBinError(bin+1);
      if(y1 <= 0.0 || y2 <= 0.0){
        continue;
      }

      double bCalc = log(y1/y2)/(x2-x1);
      double weight = (y1_err*y1_err/(y1*y1) + y2_err*y2_err/(y2*y2)) / ((x2-x1) * (x2-x1)); // sqr of error
      if(!std::isnan(weight) && !std::isnan(bCalc)){
        weights += weight;
        bAvg += bCalc*weight;
      }


      #ifdef _ZFITTER_DEBUG_
      cout << " (x1,y1) = (" << x1 << ", "  << y1  << ")  (x2,y2) = (" << x2 << ", "  << y2 << ")   bCalc = " << bCalc << "  weight= " << weight << endl;
      #endif
    }
    bAvg = bAvg/weights;


    weights = 0.0;
    for(int bin = bin1; bin < lastBinForBackground; bin++){
      double x1 = zoomHisto->GetBinCenter(bin);
      double y1 = zoomHisto->GetBinContent(bin);
      double y1_err = zoomHisto->GetBinError(bin);
      if(y1 <= 0.0){
        continue;
      }
      double weight = y1_err*y1_err*exp(2.0*bAvg*x1);
      aAvg += y1*exp(bAvg*x1)*weight;
      weights += weight;
    }
    aAvg = aAvg/weights;
    */

  }







  if(bAvg < 0.0){ // || numBinsAdded < 2.0){ // retry background seed
    cout << " Error: b parameter is negative, should be positive, trying a different seed" << endl;

    double lowBinsCenter = zoomHisto->GetBinLowEdge(3);
    double avgLowBins = ( zoomHisto->GetBinContent(1) + zoomHisto->GetBinContent(2)
                        + zoomHisto->GetBinContent(3) + zoomHisto->GetBinContent(4) )/4.0;
    int lastContentBin = HistogramUtilities::lastBinWithContent(zoomHisto);
    double highBinsCenter = zoomHisto->GetBinLowEdge(lastContentBin-1);
    double avgHighBins = ( zoomHisto->GetBinContent(lastContentBin-3) + zoomHisto->GetBinContent(lastContentBin-2)
                         + zoomHisto->GetBinContent(lastContentBin-1) + zoomHisto->GetBinContent(lastContentBin) )/4.0;
    bAvg = log(avgLowBins/avgHighBins)/(highBinsCenter - lowBinsCenter);


    if(bAvg < 0.0){
      cout << "  Error: Seed for background exponent failed : bAvg=" << bAvg << endl;
      cout << " Background will be silenced for this bin" << endl;
      fit_funct->FixParameter(4,0.0); // A
      fit_funct->FixParameter(5,1.0); // b
      fit_funct->FixParameter(6,0.0); // x_0
      fit_funct->FixParameter(7,0.0); // y_0
      ingoreBackground = true;
    }else{

      aAvg = avgLowBins*exp(bAvg*lowBinsCenter);
      cout << "  Seeding Second Try: A= " << aAvg << "  b:" << bAvg << endl;
      fit_funct->SetParameter(4,aAvg); // Amp
      fit_funct->SetParLimits(4,0.0,1.0E7 * aAvg);
      fit_funct->SetParameter(5,bAvg); // Exponential factor
      fit_funct->SetParLimits(5,1E-4*bAvg,1E4*bAvg);
      fit_funct->FixParameter(6,0.0); // x_0
      fit_funct->FixParameter(7,0.0); // y_0
    }
  }else{ // good background seed
    fit_funct->SetParameter(4,aAvg); // Amp
    fit_funct->SetParLimits(4,0.0,1.0E7 * aAvg);
    fit_funct->SetParameter(5,bAvg); // Exponential factor
    fit_funct->SetParLimits(5,1E-4*bAvg,1E4*bAvg);
    fit_funct->FixParameter(6,0.0); // x_0
    fit_funct->FixParameter(7,0.0); // y_0

    cout << "     Background Limits:  A = " << aAvg << " on [" << 0.0 << ", " << 1.0E7 * aAvg << "]" << endl;
    cout << "                         b = " << bAvg << " on [" << 1E-4*bAvg << ", " << 1E4*bAvg << "]" << endl;

  }




  //###########     Seed the Source Peak Fit Parameters  ###########
  double integral = HistogramUtilities::sumOfContent(zoomHisto, -1, true);
  double mean = 0.0; // center of peak / only for m_currentPartIndex
  if(m_massSquaredMode) mean = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0 );
  double stdDev = zoomHisto->GetStdDev()*0.8;

  /*double gausData[8];
  double lowTightCut = HistogramUtilities::getLowEdge(zoomHisto) + 0.25*(HistogramUtilities::getHighEdge(zoomHisto) - HistogramUtilities::getLowEdge(zoomHisto));
  double highTightCut = HistogramUtilities::getHighEdge(zoomHisto) - 0.25*(HistogramUtilities::getHighEdge(zoomHisto) - HistogramUtilities::getLowEdge(zoomHisto));
  HistogramUtilities::getGausData(&gausData[0], zoomHisto, true, true, lowTightCut , highTightCut);

  cout << " gaus data : " << gausData[0] << " : " << gausData[1] << " : " << gausData[2] << endl;
  */

  double stdDevCalc = 0.0;
  double newAmp = 0.0;
  double sigmaParCalc = 0.0;
  if(a_fixedParameters.size() > 4 && a_fixedParameters[2] > 0.0){
    cout << " Using Fixed Parameters to Seed Amplitude" << endl;
    if(a_useStudentT){
      newAmp = zoomHisto->GetBinContent(zoomHisto->FindBin(a_fixedParameters[1]))   * a_fixedParameters[2] / TMath::Student(0.0, a_fixedParameters[3]);
    }else{
      newAmp = zoomHisto->GetBinContent(zoomHisto->FindBin(a_fixedParameters[1])) / gaus(a_fixedParameters[1], a_fixedParameters[2], a_fixedParameters[1]);
    }

    mean = a_fixedParameters[1];
    stdDevCalc = 1.0;
    sigmaParCalc = a_fixedParameters[2];
  }else{



    //vector< double > calcData =  HistogramUtilities::calcGausData(zoomHisto, 0.0, 12, 3);
    //newAmp = calcData[0];
    //mean = calcData[1];
    //sigmaParCalc = calcData[2];


    int peakBin = zoomHisto->FindBin(mean);
    double width = HistogramUtilities::getHighEdge(zoomHisto) - zoomHisto->GetBinLowEdge(1);
    double x3 = zoomHisto->GetBinCenter(peakBin);  // point at peak
    int testBin1 = zoomHisto->FindBin(x3 - 0.08*width);
    int testBin2 = zoomHisto->FindBin(x3 + 0.12*width);
    //int testBin1 = peakBin - 6;
    //int testBin2 = peakBin + 6;

    double x4 = zoomHisto->GetBinCenter(testBin1);
    double x5 = zoomHisto->GetBinCenter(testBin2);
    double y3 = HistogramUtilities::getBinContentAvg(zoomHisto, peakBin, 2);  // point at peak
    double y4 = HistogramUtilities::getBinContentAvg(zoomHisto, testBin1, 2);
    double y5 = HistogramUtilities::getBinContentAvg(zoomHisto, testBin2, 2);

    if(!ingoreBackground){ // subtract off the background
      double new_y3 = y3 - aAvg*exp(-bAvg*x3);
      double new_y4 = y4 - aAvg*exp(-bAvg*x4);
      double new_y5 = y5 - aAvg*exp(-bAvg*x5);

      if( new_y3 > 0.0 && new_y4 > 0.0 && new_y5 > 0.0){
        y3 = new_y3;
        y4 = new_y4;
        y5 = new_y5;
      }
    }

    PhysMath::get3PointGausParams(newAmp, mean, sigmaParCalc,x3,y3,x4,y4,x5,y5);
    if(a_useStudentT){ //convert stdDev to sigma of studentT
      sigmaParCalc = sqrt(sigmaParCalc*sigmaParCalc*(ndfParamGuess - 2.0)/ndfParamGuess);
    }




    /*
    int peakBin = zoomHisto->FindBin(mean);
    double x3 = zoomHisto->GetBinCenter(peakBin);  // point at peak
    double y3 = zoomHisto->GetBinContent(peakBin);
    int testBin = zoomHisto->FindBin(x3 + 1.5*gausData[2]);
    if(testBin > zoomHisto->GetNbinsX()){
      cout << "   Test bin is off the edge of the histogram for sigma seed, using last bin" << endl;
      testBin = zoomHisto->GetNbinsX();
    }
    double x4 = zoomHisto->GetBinCenter(testBin); //point 1 stdev right of peak
    double y4 = zoomHisto->GetBinContent(testBin);

    if(!ingoreBackground){ // subtract off the background
      y3 = y3 - aAvg*exp(-bAvg*x3);
      y4 = y4 - aAvg*exp(-bAvg*x4);
    }

    stdDevCalc = sqrt(0.5*(x4-x3)*(x4-x3)/log(y3/y4));
    cout << " first internal=" << 0.5*(x4-x3)*(x4-x3)/log(y3/y4) << endl;
    cout << "   stdDevCalc= " << stdDevCalc << endl;
    if(y3 == 0.0 || y4==0.0 || y3 == y4 || y3/y4 < 0.0){
      cout << "    stdDevCalc = " << stdDevCalc << "   from (x3,y3)=(" << x3 << ", " << y3 << ")   and (x4,y4)=(" << x4 << ", " << y4 << ")" << endl;
      cout << "      setting it to guess " << endl;
      stdDevCalc = sigmaGuess;
    }
    sigmaParCalc = sqrt(stdDevCalc*stdDevCalc*(ndfParamGuess - 2.0)/ndfParamGuess);
    cout << "      internal calc=" << stdDevCalc*stdDevCalc*(ndfParamGuess - 2.0)/ndfParamGuess << endl;
    cout << "   sigmaParCalc= " << sigmaParCalc << endl;

    newAmp = y3*sigmaParCalc/ TMath::Student(0.0, ndfParamGuess);
    */
  }



  cout << " Starting Parameters:   N:" << newAmp << "  mu:" << mean << "   sig:" << sigmaParCalc;
  if(a_useStudentT) cout << "   nu:" << ndfParamGuess << endl;
  cout << "                        A:" << aAvg << "  b:" << bAvg << "  x_0: 0.0   y_0:0.0" << endl;


  if(newAmp == 0 || std::isnan(mean) || std::isnan(sigmaParCalc)){
    cout << " Error in ZFitter::fitIsolation_BTOF_withBackground: integral = "<< newAmp << " mean = "<< mean << "  stdDev = "<< sigmaParCalc << endl;
    //delete fittingCanvas;
    return returnParams;
  }

  TString imageName = TString::Format("%s_withBackground_rap_%02d_mTm0_%02d",
      m_partInfo->GetParticleName(m_currentPartIndex,0).Data(),m_currentRapBin,m_currentMtM0Bin);
  zoomHisto->SetName(imageName.Data());
  zoomHisto->SetTitle(Form("BTOF Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] (GeV/c) p=%2.3f : zoom on %s",
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
    m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
    m_currentMomOfBin,
    m_partInfo->GetParticleName(m_currentPartIndex,0).Data()
    ));
  if(m_massSquaredMode){
    zoomHisto->GetXaxis()->SetTitle("m^{2}");
    zoomHisto->GetYaxis()->SetTitle("dN/d(m^{2})}");
  }else{
    zoomHisto->GetXaxis()->SetTitle("#beta^{-1}");
    zoomHisto->GetYaxis()->SetTitle("dN/d#beta^{-1}");
  }

  zoomHisto->SetStats(kTRUE);
  gStyle->SetOptFit(0111);




  //##########  Setup Initial Seed Function to Draw for Debugging  #############
  fit_funct->SetParameter(0,newAmp);
  fit_funct->SetParLimits(0,0.01*newAmp, 10.0*newAmp);
  fit_funct->SetParameter(1,mean);
  fit_funct->SetParLimits(1, mean - 0.15, mean + 0.15);
  fit_funct->SetParameter(2,sigmaParCalc);
  fit_funct->SetParLimits(2,0.1*sigmaParCalc,10.0*sigmaParCalc);
  if(a_useStudentT){
    fit_funct->SetParameter(3,ndfParamGuess);
    fit_funct->SetParLimits(3,2.0,1000.0);
  }
  TF1* initialFunct = NULL;
  initialFunct = (TF1*) fit_funct->Clone();
  initialFunct->SetLineColor(kYellow);

  // Create
  TF1* source_funct_init = NULL;
  if(a_useStudentT){
    source_funct_init = new TF1("StudentT_Source_Init", "[0]*TMath::Student((x - [1])/[2], [3])/[2]",lowRange,highRange);
    source_funct_init->SetNpx(200);
    source_funct_init->SetLineColor(kYellow+1);
    source_funct_init->SetLineStyle(2);
    source_funct_init->SetParNames("N","#mu","#sigma","#nu");
    source_funct_init->SetParameter(0,fit_funct->GetParameter(0));
    source_funct_init->SetParameter(1,fit_funct->GetParameter(1));
    source_funct_init->SetParameter(2,fit_funct->GetParameter(2));
    source_funct_init->SetParameter(3,fit_funct->GetParameter(3));
  }else{
    source_funct_init = new TF1("gaus_Source_Init", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
    source_funct_init->SetNpx(200);
    source_funct_init->SetLineColor(kYellow+1);
    source_funct_init->SetLineStyle(2);
    source_funct_init->SetParNames("N","#mu","#sigma");
    source_funct_init->SetParameter(0,fit_funct->GetParameter(0));
    source_funct_init->SetParameter(1,fit_funct->GetParameter(1));
    source_funct_init->SetParameter(2,fit_funct->GetParameter(2));
  }

  // Create Red Background Function For Drawing
  TF1* background_funct_init = new TF1("Background_PowerLawFunction_Init", "[0]*exp(-[1]*(x - [2])) + [3]",lowRange,highRange);
  background_funct_init->SetNpx(200);
  background_funct_init->SetLineColor(kYellow+2);
  background_funct_init->SetLineStyle(2);
  background_funct_init->SetParNames("A","b","x_{0}","y_{0}");
  background_funct_init->SetParameter(0,fit_funct->GetParameter(4));
  background_funct_init->SetParameter(1,fit_funct->GetParameter(5));
  background_funct_init->SetParameter(2,fit_funct->GetParameter(6));
  background_funct_init->SetParameter(3,fit_funct->GetParameter(7));


  //Fix the parameters in argument
  if(a_fixedParameters.size() > 0){
    for(int jjj = 0; jjj < 8; jjj++){
      if(a_fixedParameters[jjj] != -999){
        fit_funct->FixParameter(jjj, a_fixedParameters[jjj]);
        cout << "    Fixing Parameter " << jjj << " : " << a_fixedParameters[jjj] << endl;
      }
    }
  }




  //###############   FIT   THE   CURVE   #####################
  cout << "      Starting the fit in funct" << endl;
  int status = -999;
  string fitOptions = "EQ";
  #ifdef _ZFITTER_DEBUG_
    fitOptions = "EV";
  #endif
  if(m_doIntegralOfBinForFits) fitOptions += "I";
  m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(zoomHisto,-1,true);
  if(a_useStudentT){
    status = zoomHisto->Fit("StudentT_withPowerLawFunction",fitOptions.c_str());
  }else{
    status = zoomHisto->Fit("Gaus_withPowerLawFunction",fitOptions.c_str());
  }


  if(status == 0){
    TF1* fitResult = NULL;
      if(a_useStudentT){
        fitResult = zoomHisto->GetFunction("StudentT_withPowerLawFunction");
      }else{
        fitResult = zoomHisto->GetFunction("Gaus_withPowerLawFunction");
      }
    if(fitResult->GetParameter(0) != 0.0){
      returnParams[0] = fitResult->GetParameter(0);
      returnParams[1] = fitResult->GetParameter(1);
      returnParams[2] = fitResult->GetParameter(2);
      returnParams[3] = fitResult->GetParameter(3);
      returnParams[4] = fitResult->GetParameter(4);
      returnParams[5] = fitResult->GetParameter(5);
      returnParams[6] = fitResult->GetParameter(6);
      returnParams[7] = fitResult->GetParameter(7);

      returnParams[8] =  fitResult->GetParError(0);
      returnParams[9] =  fitResult->GetParError(1);
      returnParams[10] = fitResult->GetParError(2);
      returnParams[11] = fitResult->GetParError(3);
      returnParams[12] = fitResult->GetParError(4);
      returnParams[13] = fitResult->GetParError(5);
      returnParams[14] = fitResult->GetParError(6);
      returnParams[15] = fitResult->GetParError(7);
      returnParams[16] = status;
      returnParams[17] = fitResult->GetChisquare() / fitResult->GetNDF();

    }else{
      cout << " Error in ZFitter::fitIsolation_BTOF_withBackground: fit status = "<< status << endl;
    }


    // Create Violet Source Peak For Drawinga
    TF1* source_funct = NULL;

    if(a_useStudentT){
      source_funct = new TF1("StudentT_Source", "[0]*TMath::Student((x - [1])/[2], [3])/[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma","#nu");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
      source_funct->SetParameter(3,returnParams[3]);
    }else{
      source_funct = new TF1("gaus_Source", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
    }


    // Create Red Background Function For Drawing
    TF1* background_funct = new TF1("Background_PowerLawFunction", "[0]*exp(-[1]*(x - [2])) + [3]",lowRange,highRange);
    background_funct->SetNpx(200);
    background_funct->SetLineColor(kRed);
    background_funct->SetLineStyle(2);
    background_funct->SetParNames("A","b","x_{0}","y_{0}");
    background_funct->SetParameter(0,returnParams[4]);
    background_funct->SetParameter(1,returnParams[5]);
    background_funct->SetParameter(2,returnParams[6]);
    background_funct->SetParameter(3,returnParams[7]);



    //draw result
    fittingCanvas->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits){
      initialFunct->Draw("same");
      source_funct_init->Draw("same");
      background_funct_init->Draw("same");
    }
    fitResult->Draw("same");
    background_funct->Draw("same");
    source_funct->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();

    cout << "     ";
    if(m_inAllCentMode){
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),imageName.Data()));
      }
    }else{
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(), m_currentCentIndex,imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(), m_currentCentIndex,imageName.Data()));
      }
    }

    delete source_funct;
    delete background_funct;
  }else{  //////   BAD FIT / NOT OPTIMAL



    TF1* fitResult = zoomHisto->GetFunction("StudentT_withPowerLawFunction");
    cout << "Error in ZFitter::fitIsolation_BTOF_withBackground: fit status: " << status << "   " << m_partInfo->GetParticleName(m_currentPartIndex,0)
         << " Rap Bin " << m_currentRapBin << " mTm0 Bin " << m_currentMtM0Bin << endl;

    cout << "  Pushing back parameters anyway..." << endl;
    if(fitResult){
      cout << " Pushing back parameters in Fit Result : " << fitResult << endl;
      returnParams[0] = fitResult->GetParameter(0);
      returnParams[1] = fitResult->GetParameter(1);
      returnParams[2] = fitResult->GetParameter(2);
      returnParams[3] = fitResult->GetParameter(3);
      returnParams[4] = fitResult->GetParameter(4);
      returnParams[5] = fitResult->GetParameter(5);
      returnParams[6] = fitResult->GetParameter(6);
      returnParams[7] = fitResult->GetParameter(7);

      returnParams[8] = fitResult->GetParError(0);
      returnParams[9] = fitResult->GetParError(1);
      returnParams[10] = fitResult->GetParError(2);
      returnParams[11] = fitResult->GetParError(3);
      returnParams[12] = fitResult->GetParError(4);
      returnParams[13] = fitResult->GetParError(5);
      returnParams[14] = fitResult->GetParError(6);
      returnParams[15] = fitResult->GetParError(7);
      returnParams[16] = status;
    }else{
      cout << " Pushing back parameters in Fit Funciton : " << fit_funct << endl;
      returnParams[0] = fit_funct->GetParameter(0);
      returnParams[1] = fit_funct->GetParameter(1);
      returnParams[2] = fit_funct->GetParameter(2);
      returnParams[3] = fit_funct->GetParameter(3);
      returnParams[4] = fit_funct->GetParameter(4);
      returnParams[5] = fit_funct->GetParameter(5);
      returnParams[6] = fit_funct->GetParameter(6);
      returnParams[7] = fit_funct->GetParameter(7);

      returnParams[8] = fit_funct->GetParError(0);
      returnParams[9] = fit_funct->GetParError(1);
      returnParams[10] = fit_funct->GetParError(2);
      returnParams[11] = fit_funct->GetParError(3);
      returnParams[12] = fit_funct->GetParError(4);
      returnParams[13] = fit_funct->GetParError(5);
      returnParams[14] = fit_funct->GetParError(6);
      returnParams[15] = fit_funct->GetParError(7);
      returnParams[16] = status;
    }

    // Create Violet Source Peak For Drawinga
    TF1* source_funct = NULL;

    if(a_useStudentT){
      source_funct = new TF1("StudentT_Source", "[0]*TMath::Student((x - [1])/[2], [3])/[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma","#nu");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
      source_funct->SetParameter(3,returnParams[3]);
    }else{
      source_funct = new TF1("gaus_Source", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
      source_funct->SetNpx(200);
      source_funct->SetLineColor(kViolet);
      source_funct->SetLineStyle(2);
      source_funct->SetParNames("N","#mu","#sigma");
      source_funct->SetParameter(0,returnParams[0]);
      source_funct->SetParameter(1,returnParams[1]);
      source_funct->SetParameter(2,returnParams[2]);
    }


    // Create Red Background Function For Drawing
    TF1* background_funct = new TF1("Background_PowerLawFunction", "[0]*exp(-[1]*(x - [2])) + [3]",lowRange,highRange);
    background_funct->SetNpx(200);
    background_funct->SetLineColor(kRed);
    background_funct->SetLineStyle(2);
    background_funct->SetParNames("A","b","x_{0}","y_{0}");
    background_funct->SetParameter(0,returnParams[4]);
    background_funct->SetParameter(1,returnParams[5]);
    background_funct->SetParameter(2,returnParams[6]);
    background_funct->SetParameter(3,returnParams[7]);






    cout << "Drawing Now" << endl;
    fittingCanvas->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    initialFunct->Draw("same");
    source_funct_init->Draw("same");
    background_funct_init->Draw("same");
    if(returnParams[0] < 0.0 || std::isnan(returnParams[1]) || returnParams[2] <= 0.0 || (a_useStudentT && returnParams[3] <= 1.0)
        || returnParams[4] <= 0.0 || returnParams[5] <= 0.0){
      cout << " Really bad fit, not drawing any funcitons" << endl;
    }else{
      background_funct->Draw("same");
      source_funct->Draw("same");
    }

    //fitResult->SetLineColor(kOrange);
    //fitResult->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();

    cout << "     ";
    if(m_inAllCentMode){
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),imageName.Data()));
      }
    }else{
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(), m_currentCentIndex,imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/%s__________BAD___FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(), m_currentCentIndex,imageName.Data()));
      }
    }

    delete source_funct;
    delete background_funct;

  }


  //delete zoomHisto;
  delete initialFunct;
  delete background_funct_init;
  delete source_funct_init;
  delete fittingCanvas;
  delete zoomHisto;


  cout << "  Finished fitIsolation_BTOF_withBackground" << endl;


  return returnParams;

}



//============================================================================================================
//============================================================================================================

void ZFitter::extrapolateMeanAndWidthParams_BTOF(int a_particleIndex){

  cout << " Starting extrapolateMeanAndWidthParams_BTOF for part " << a_particleIndex << endl;
  TCanvas* fittingCanvas = new TCanvas("FittingCanvasMeanAndWidth");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);
  //fittingCanvas->Divide(1,2);

  TF1* meanExtrap = new TF1("meanExtrap","[0]*x + [1]",0.005,5.0);
  TF1* stdDevExtrap = new TF1("stdDevExtrap","[0]*x + [1]",0.005,5.0);
  TF1* dofExtrap = new TF1("dofExtrap","[0]*x + [1]",0.005,5.0);

  if(m_extrap_fn_vers == 1){
    //    stdDevExtrap = new TF1("stdDevExtrap","(pow([0]*x+[1],2)+[2])*(x<[3])+(pow([0]*[3]+[1],2)+[2])*(x>=[3])",0.005,5.0);
    stdDevExtrap = new TF1("stdDevExtrap","([0]*pow(x-[2],2)+[1])*(x<[3])+([0]*pow([3]-[2],2)+[1])*(x>=[3])",0.005,5.0);
    //    meanExtrap = new TF1("meanExtrap","([0]*x+[1])*(x<[2])+([0]*[2]+[1])*(x>=[2])",0.005,5.0);
    //  meanExtrap = new TF1("meanExtrap","([0]*x+[1])*(x<[3])+([2]*x+[3]*([0]-[2])+[1])*(x>=[3]&&x<[4])+([2]*[4]+[3]*([0]-[2])+[1])*(x>=[4])",0.005,5.0);
  }

  if(m_extrap_fn_vers == 2){
    meanExtrap = new TF1("meanExtrap","[0]*exp(-[1]*x)",0.005,5.0);
    stdDevExtrap = new TF1("stdDevExtrap","[0]*exp([1]*x)+[2]",0.005,5.0);
  }


  //Get Mean Histogram
  TH1D* mean_histo = NULL;
  TH1D* stdDev_histo = NULL;
  TH1D* dof_histo = NULL;

  if(m_separatePM_bToF_allCent){
    if(m_currentPartCharge == -1){
      mean_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][1][1]->ProjectionY(Form("mean_proj_%s_%s_Rap%02d",
         m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
      stdDev_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][2][1]->ProjectionY(Form("stdDev_proj_%s_%s_Rap%02d",
         m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
      if(m_useStudentTDistributionsForTOF){
        dof_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][3][1]->ProjectionY(Form("dof_proj_%s_%s_Rap%02d",
           m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
           m_currentRapBin,m_currentRapBin);
      }
    }else{
      mean_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][1][0]->ProjectionY(Form("mean_proj_%s_%s_Rap%02d",
         m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
      stdDev_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][2][0]->ProjectionY(Form("stdDev_proj_%s_%s_Rap%02d",
         m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
      if(m_useStudentTDistributionsForTOF){
        dof_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][3][0]->ProjectionY(Form("dof_proj_%s_%s_Rap%02d",
          m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
          m_currentRapBin,m_currentRapBin);
      }
    }
  }else{
    mean_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->ProjectionY(Form("mean_proj_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,0).Data(),m_partInfo->GetParticleName(a_particleIndex,0).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);
    stdDev_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->ProjectionY(Form("stdDev_proj_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,0).Data(),m_partInfo->GetParticleName(a_particleIndex,0).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);
    if(m_useStudentTDistributionsForTOF){
      dof_histo = m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][3][m_currentPlusMinusIndex]->ProjectionY(Form("dof_proj_%s_%s_Rap%02d",
         m_partInfo->GetParticleName(m_currentPartIndex,0).Data(),m_partInfo->GetParticleName(a_particleIndex,0).Data(),m_currentRapBin),
         m_currentRapBin,m_currentRapBin);
    }
  }


  //############ Make Prediction/Fit Difference  ##################
  TH1D* mean_diff_histo = (TH1D*) mean_histo->Clone();
  mean_diff_histo->SetName("mean_diff_histo");
  mean_diff_histo->GetYaxis()->SetTitle("#beta^{-1} - #beta^{-1}_{predict}");
  mean_diff_histo->GetXaxis()->SetTitle(Form("m_{T} - m_{%s}",m_partInfo->GetParticleSymbol(m_currentPartIndex,0).Data()));
  for(int binX = 1; binX <= mean_histo->GetNbinsX(); binX++){
    double center = mean_histo->GetBinCenter(binX);
    double content = mean_histo->GetBinContent(binX);
    double error =  mean_histo->GetBinError(binX);
    double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,center,m_partInfo->GetParticleMass(m_currentPartIndex));
    // ln(fit/base)-ln(predict/base)
    double newContent = content - PhysMath::invBeta(m_partInfo->GetParticleMass(a_particleIndex), momentum) + PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), momentum);
    mean_diff_histo->SetBinContent(binX,newContent);
    mean_diff_histo->SetBinError(binX,error);
  }
  cout << " Mean Diff Histo Created " << endl;


  int minReqBins = 2;
  int numBinsWithContent = 0;
  int endBin = m_currentMtM0Bin;
  double endMtM0 = mean_diff_histo->GetBinCenter(endBin - 1) + 0.001;
  //endBin = m_currentMtM0Bin;
  //endMtM0 = mean_diff_histo->GetBinCenter(endBin - 1) + 0.001;

  int firstBinWithData = 999;
  for(int bin = endBin - 1; bin > 0 ; bin--){
    if(stdDev_histo->GetBinContent(bin) > 0.0){
      numBinsWithContent++;
      firstBinWithData = bin;
    }else{
      if(numBinsWithContent < minReqBins){
        numBinsWithContent = -999;
      }
    }
  }

  if(numBinsWithContent < minReqBins){
    cout << "End bin : " << endBin << " by not enough bins to constrain: " << numBinsWithContent <<  " < " << minReqBins << endl;
    cout << " You may have m_bTOFSeparatedStops high compared to the m_partBTOF_MomStartExtrapolation_Matrix_Plus" << endl;

    delete meanExtrap;
    delete stdDevExtrap;
    delete dofExtrap;
    delete mean_histo;
    delete stdDev_histo;
    delete dof_histo;
    delete mean_diff_histo;
    delete fittingCanvas;
    return;
    // Zero the rest of the fits?
  }




  //#########  FIT THE MEAN HISTOGRAM ###############
  double firstBin = 1;
  if(m_extrap_fn_vers == 2) firstBin = firstBinWithData;
  if(endBin - firstBinWithData > 4) firstBin = endBin - 4;
  //if(endBin - firstBinWithData > 5) firstBin = endBin - 5;
  //if(endBin - firstBinWithData > 6) firstBin = endBin - 6;
  //if(endBin - firstBinWithData > 7) firstBin = endBin - 7;
  //if(endBin - firstBinWithData > 8) firstBin = endBin - 8;
  //if(endBin - firstBinWithData > 9) firstBin = endBin - 9;
  //if(endBin - firstBinWithData > 10) firstBin = endBin - 10;
  double x1 = mean_diff_histo->GetBinCenter(firstBin);
  double x2 = mean_diff_histo->GetBinCenter(endBin - 1);
  double y1 = mean_diff_histo->GetBinContent(firstBin);
  double y2 = mean_diff_histo->GetBinContent(endBin - 1);
    if(m_extrap_fn_vers == 2){
    double yHold[20];
    double yMean = 0.0;
    double ySigma = 0.0;
    for(int iii=firstBin; iii<endBin-1; iii++){
      yHold[iii] = (mean_diff_histo->GetBinContent(iii));
      yMean += yHold[iii];
    }
    yMean = yMean/(endBin-1 - firstBin);
    for(int iii=firstBin; iii<endBin-1; iii++){
      ySigma += pow(yHold[iii] - yMean,2);
    }
    ySigma = sqrt(ySigma/(endBin-1 - firstBin));

    int badBin = -1;

    cout << "yMean: " << yMean << endl;
    cout << "ySigma: " << ySigma << endl;

    cout << "lowLimit: " << yMean-ySigma << "   highLimit: " << yMean+ySigma << endl;

    for(int iii=firstBin; iii<endBin-1; iii++){
      cout << "xValue: " << mean_diff_histo->GetBinCenter(iii) << "   yValue: " << yHold[iii] << endl;
      if( (yHold[iii] > (yMean+ySigma)) || (yHold[iii]<(yMean-ySigma)) ){
        cout << "badBin!" << endl;
        badBin = iii;
        break;
      }
    }
    if(badBin != -1){
      if(badBin == firstBin){
        x1 = mean_diff_histo->GetBinCenter(firstBin + 1);
        y1 = mean_diff_histo->GetBinContent(firstBin + 1);
      }else if(badBin == endBin - 1){
          x2 = mean_diff_histo->GetBinCenter(endBin - 2);
          y2 = mean_diff_histo->GetBinContent(endBin - 2);
      }
    }
  }
  double slope = (y2-y1)/(x2-x1);
  double intercept = y2 - x2*slope;
  meanExtrap->SetParameter(0,slope);
  meanExtrap->SetParameter(1,intercept);
  if(slope > 0.0){
    meanExtrap->SetParLimits(0,-10.0*slope,10.0*slope);
  }else{
    meanExtrap->SetParLimits(0,10.0*slope,-10.0*slope);
  }
  if(intercept > 0.0){
    meanExtrap->SetParLimits(1,-5.0*intercept,10.0*intercept);
  }else{
    meanExtrap->SetParLimits(1,10.0*intercept,-5.0*intercept);
  }

  /* if(m_extrap_fn_vers == 1){
    meanExtrap->SetParameter(2,0.1);
    meanExtrap->SetParameter(3,endMtM0-0.1);
    meanExtrap->FixParameter(4,endMtM0);
    //meanExtrap->FixParameter(3,100);
    //meanExtrap->SetParLimits(2,0,10);
    meanExtrap->SetParNames("slope","intercept","slope2","shift1","shift2");
  }*/

  if(m_extrap_fn_vers == 2){
    meanExtrap->SetParameter(0,-1);
    meanExtrap->SetParLimits(0,-10.0,10.0);
    meanExtrap->SetParameter(1,3);
    meanExtrap->SetParLimits(1,0.0,20.0);
    meanExtrap->SetParNames("scale","decayConstant");
  }

  int mean_status = mean_diff_histo->Fit("meanExtrap","","",x1-0.0001,x2+0.0001);
  TF1* fitResult = mean_diff_histo->GetFunction("meanExtrap");

  //Push the results to the histograms
  if(fitResult){
    for(int binX = endBin; binX <= mean_histo->GetNbinsX(); binX++){
      double center = mean_histo->GetBinCenter(binX);
      double error =  mean_histo->GetBinError(binX);
      double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,center,m_partInfo->GetParticleMass(m_currentPartIndex));
      double newMean = fitResult->Eval(center) + PhysMath::invBeta(m_partInfo->GetParticleMass(a_particleIndex), momentum) - PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), momentum);
      m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,binX, newMean);
    }
  }else{
    cout << "Error: Mean Extrapolation Failed with Status: " << mean_status << "  : Exiting!" << endl;
    delete meanExtrap;
    delete stdDevExtrap;
    delete dofExtrap;
    delete mean_histo;
    delete stdDev_histo;
    delete dof_histo;
    delete mean_diff_histo;
    delete fittingCanvas;
    return;
  }



  //##########  Make Std Dev Difference Histogram ###############

  stdDev_histo->SetTitle(Form(";m_{T} - m_{%s}; #sigma",m_partInfo->GetParticleSymbol(m_currentPartIndex,0).Data()));
  cout << " now fit stdDev " << endl;
  //FIT THE STDDEV DIFF HISTOGRAM
  x1 = stdDev_histo->GetBinCenter(firstBin);
  x2 = stdDev_histo->GetBinCenter(endBin - 1);
  y1 = stdDev_histo->GetBinContent(firstBin);
  y2 = stdDev_histo->GetBinContent(endBin - 1);
  if(m_extrap_fn_vers == 2){
    double yHold[20];
    double yMean = 0.0;
    double ySigma = 0.0;
    for(int iii=firstBin; iii<endBin-1; iii++){
      yHold[iii] = (stdDev_histo->GetBinContent(iii));
      yMean += yHold[iii];
    }
    yMean = yMean/(endBin-1 - firstBin);
    for(int iii=firstBin; iii<endBin-1; iii++){
      ySigma += pow(yHold[iii] - yMean,2);
    }
    ySigma = sqrt(ySigma/(endBin-1 - firstBin));

    int badBin = -1;

    cout << "yMean: " << yMean << endl;
    cout << "ySigma: " << ySigma << endl;

    cout << "lowLimit: " << yMean-ySigma << "   highLimit: " << yMean+ySigma << endl;

    for(int iii=firstBin; iii<endBin-1; iii++){
      cout << "xValue: " << stdDev_histo->GetBinCenter(iii) << "   yValue: " << yHold[iii] << endl;
      if( (yHold[iii] > (yMean+ySigma)) || (yHold[iii]<(yMean-ySigma)) ){
        cout << "badBin!" << endl;
        badBin = iii;
        break;
      }
    }
    if(badBin != -1){
      if(badBin == firstBin){
        x1 = stdDev_histo->GetBinCenter(firstBin + 1);
        y1 = stdDev_histo->GetBinContent(firstBin + 1);
      }else if(badBin == endBin - 1){
          x2 = stdDev_histo->GetBinCenter(endBin - 2);
          y2 = stdDev_histo->GetBinContent(endBin - 2);
      }
    }
  }
  slope = (y2-y1)/(x2-x1);
  intercept = y2 - x2*slope;
  stdDevExtrap->SetParameter(0,slope);
  stdDevExtrap->SetParameter(1,intercept);
  if(slope > 0.0){
    stdDevExtrap->SetParLimits(0,-5.0*slope,10.0*slope);
  }else{
    stdDevExtrap->SetParLimits(0,10.0*slope,-5.0*slope);
  }
  if(intercept > 0.0){
    stdDevExtrap->SetParLimits(1,-5.0*intercept,10.0*intercept);
  }else{
    stdDevExtrap->SetParLimits(1,10.0*intercept,-5.0*intercept);
  }

  if(m_extrap_fn_vers == 1){
    stdDev_histo->GetXaxis()->SetRangeUser(0.2,endMtM0);
    if(endMtM0-0.4 > 0.2) stdDev_histo->GetXaxis()->SetRangeUser(endMtM0-0.4,endMtM0);
    int stdDev_min_bin = stdDev_histo->GetMinimumBin();
    stdDev_histo->GetXaxis()->SetRangeUser(stdDev_histo->GetXaxis()->GetBinLowEdge(1),stdDev_histo->GetXaxis()->GetBinLowEdge(stdDev_histo->GetNbinsX()+1));
    double stdDev_min = stdDev_histo->GetBinCenter(stdDev_min_bin);
    double stdDev_min_val = stdDev_histo->GetBinContent(stdDev_min_bin);
    stdDevExtrap->SetParameter(0,1);
    stdDevExtrap->SetParLimits(0,0.0,100.0);
    stdDevExtrap->SetParameter(1,0.01);
    stdDevExtrap->SetParLimits(1,0.0,100.0);
    stdDevExtrap->FixParameter(2,stdDev_min);
    stdDevExtrap->FixParameter(3,stdDev_min);
   // stdDevExtrap->SetParLimits(3,0.0,10.0);
    stdDevExtrap->SetParNames("a","c","#mu","shift");
  }

  if(m_extrap_fn_vers == 2){
    stdDevExtrap->SetParameter(0,-1);
    stdDevExtrap->SetParLimits(0,-10.0,10.0);
    stdDevExtrap->SetParameter(1,-3);
    stdDevExtrap->SetParLimits(1,-10.0,10.);
    stdDevExtrap->SetParameter(2,0.0);
    stdDevExtrap->SetParLimits(2,-0.05,0.05);
    stdDevExtrap->SetParNames("scale","decayConstant","verticalShift");
  }

  int std_status = stdDev_histo->Fit("stdDevExtrap","","",x1-0.0001,x2+0.0001);
  TF1* fitResult_std = stdDev_histo->GetFunction("stdDevExtrap");

  if(m_extrap_fn_vers == 1){
    cout << "part binning: " << m_currentPartIndex << endl;
    cout << "part interest: " << a_particleIndex << endl;
    cout << "stdDev value at mTm0=" << endMtM0+0.5 << ": " << stdDevExtrap->Eval(endMtM0 +0.5) << endl;
  }

  //Push the results to the histograms
  if(fitResult_std){
    for(int binX = endBin; binX <= stdDev_histo->GetNbinsX(); binX++){
      double center = stdDev_histo->GetBinCenter(binX);
      double fitSTDDev = fitResult_std->Eval(center);
      if(fitSTDDev <= 0.0) fitSTDDev = 1e-12;
      m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,binX, fitSTDDev);
    }
  }else{
    cout << "Error: Std Dev Extrapolation Failed with Status: " << std_status << "  : Exiting!" << endl;
    delete meanExtrap;
    delete stdDevExtrap;
    delete dofExtrap;
    delete mean_histo;
    delete stdDev_histo;
    delete dof_histo;
    delete mean_diff_histo;
    delete fittingCanvas;
    return;
  }


  if(m_useStudentTDistributionsForTOF){
    //############ Make Prediction/Fit Difference  ##################
    dof_histo->GetYaxis()->SetTitle("#nu");
    dof_histo->GetXaxis()->SetTitle(Form("m_{T} - m_{%s}",m_partInfo->GetParticleSymbol(m_currentPartIndex,0).Data()));
    //FIT THE STDDEV DIFF HISTOGRAM
    x1 = dof_histo->GetBinCenter(firstBin);
    x2 = dof_histo->GetBinCenter(endBin - 1);
    y1 = dof_histo->GetBinContent(firstBin);
    y2 = dof_histo->GetBinContent(endBin - 1);
    slope = (y2-y1)/(x2-x1);
    intercept = y2 - x2*slope;
    dofExtrap->SetParameter(0,slope);
    dofExtrap->SetParameter(1,intercept);
    if(slope > 0.0){
      dofExtrap->SetParLimits(0,-5.0*slope,10.0*slope);
    }else{
      dofExtrap->SetParLimits(0,10.0*slope,-5.0*slope);
    }
    if(intercept > 0.0){
      dofExtrap->SetParLimits(1,-5.0*intercept,10.0*intercept);
    }else{
      dofExtrap->SetParLimits(1,10.0*intercept,-5.0*intercept);
    }
    int dof_status = dof_histo->Fit("dofExtrap","","",x1 - 0.0001,x2 + 0.0001);
    TF1* fitResult_dof = dof_histo->GetFunction("dofExtrap");

    //Push the results to the histograms
    if(fitResult_dof){
      for(int binX = endBin; binX <= dof_histo->GetNbinsX(); binX++){
        double center = dof_histo->GetBinCenter(binX);
        double fitNu = fitResult_dof->Eval(center);
        if(fitNu <= 1.0) fitNu = 1.001;
        m_Fit_Data_ZbTOF[m_currentPartIndex][a_particleIndex][3][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,binX, fitNu);
      }
    }else{
      cout << "Error: DOF Extrapolation Failed with Status: " << dof_status << "  : Exiting!" << endl;
      delete meanExtrap;
      delete stdDevExtrap;
      delete dofExtrap;
      delete mean_histo;
      delete stdDev_histo;
      delete dof_histo;
      delete mean_diff_histo;
      delete fittingCanvas;
      return;
    }

    // DRAW THE FITS
    fittingCanvas->cd();
    gStyle->SetOptFit(0111);
    dof_histo->GetXaxis()->SetRangeUser(x1-0.1,x2+0.1);
    dof_histo->SetStats(true);
    dof_histo->Draw("E");
    //gPad->SetLogy();
    gPad->SetRightMargin(0.4);
    fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/DOFExtrap_%s_in_%s_space_Rap_%02d_MtM0Bin_%02d.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin,m_currentMtM0Bin));

  }










  // DRAW THE FITS
  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  //proj->GetXaxis()->SetRangeUser(0.0,x2+0.5);
  mean_diff_histo->SetStats(true);
  mean_diff_histo->GetXaxis()->SetRangeUser(x1-0.01,x2+0.01);
  mean_diff_histo->Draw("E");
  //gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/MeanExtrap_%s_in_%s_space_Rap_%02d_MtM0Bin_%02d.png",m_imagePreDir.c_str(),
    m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin,m_currentMtM0Bin));

  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  //proj->GetXaxis()->SetRangeUser(0.0,x2+0.5);
  stdDev_histo->SetStats(true);
  stdDev_histo->GetXaxis()->SetRangeUser(x1-0.01,x2+0.01);
  stdDev_histo->Draw("E");
  //gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  gPad->Update();
  //gPad->ProccessEvents();
  //if(fitResult) fitResult->Draw("same");
  fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/StdDevExtrap_%s_in_%s_space_Rap_%02d_MtM0Bin_%02d.png",m_imagePreDir.c_str(),
    m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin,m_currentMtM0Bin));


  delete meanExtrap;
  delete stdDevExtrap;
  delete dofExtrap;
  delete mean_histo;
  delete stdDev_histo;
  delete dof_histo;
  delete mean_diff_histo;
  delete fittingCanvas;
  //delete stdDev_diff_histo;
  //m_partHasBeenExtrapolated[a_particleIndex] = true;

  m_partHasBeenExtrapolatedCharge[a_particleIndex][m_currentPlusMinusIndex] = true;

  cout << "Particle Successfully Extrapolated" << endl;



  return;
}


//============================================================================================================
//============================================================================================================

void ZFitter::extrapolateMeanAndWidthParams_BTOF_m2(){

  //This fits all necessary parameters with student ndf

  #ifdef _ZFITTER_DEBUG_
  cout << " Starting postFit_BTOF_parameterFits_m2 for part " << m_currentPartIndex << endl;
  #endif

  TCanvas* fittingCanvas = new TCanvas("FittingCanvasMeanAndWidth");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);
  //fittingCanvas->Divide(1,2);


  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
    if(m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) continue;

    double lowMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity,m_BTOF_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex],m_partInfo->GetParticleMass(m_currentPartIndex));
    double highMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity,m_BTOF_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex],m_partInfo->GetParticleMass(m_currentPartIndex));
    if(lowMtM0 < 0.0){
      cout << " the rapidity is too extreme for the extrapolation procedure you have given setting low mTm0 to 0" << endl;
      lowMtM0 = 0.0;
    }
    if(highMtM0 < 0.0){
      cout << " the rapidity is too extreme for the extrapolation procedure, becuase high mTm0 is negative... (" << highMtM0 << " ) so continuing on without extrapolation" << endl;
      continue;
    }
    if(m_BTOF_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex] > 900){
      cout << " no extrapolation range so continuing " << endl;
      continue;
    }
    double midMtM0 = 0.5*(lowMtM0 + highMtM0);
    #ifdef _ZFITTER_DEBUG_
      cout << "  extrapolating particle " << m_partInfo->GetParticleName(partIndex) << " with mTm0 fit range from " << lowMtM0 << "  to  " << highMtM0 << endl;
      cout << "      momentum range was " << m_BTOF_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex] << "  to  " << m_BTOF_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex] << endl;
    #endif


    //make the histograms to extrapolate
    TH1D* param1DHistos[4];
    TF1* extrapFunctions[4];
    double numberOfParams = 3;
    string distributionName = "Gaussian";
    if(m_useStudentTDistributionsForTOF){
      distributionName = "Student-T";
      numberOfParams = 4;
    }

    int pmIndexOfHistoToGrab = 0;
    if(m_separatePM_bToF_allCent && m_currentPartCharge < 0) pmIndexOfHistoToGrab = 1;

    string paramNames[4] = {"N","mean","width","ndf"};
    string paramSymbols[4] = {"N","#mu","#sigma","#nu"};
    string functionNames[4];
    if(m_inAllCentMode){
      for(int jjj = 0; jjj < numberOfParams; jjj++){
        param1DHistos[jjj] = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][jjj][pmIndexOfHistoToGrab]->ProjectionY(Form("param_proj_%s_%s_%s_Rap%02d",
           paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),m_currentRapBin),
           m_currentRapBin,m_currentRapBin);
        functionNames[jjj] = Form("param_parab_%s_%s_%s_Rap%02d",paramNames[jjj].c_str(),
          m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
          m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),m_currentRapBin);
        param1DHistos[jjj]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),paramSymbols[jjj].c_str()));
      }
    }else{
      for(int jjj = 0; jjj < numberOfParams; jjj++){
        param1DHistos[jjj] = m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][partIndex][jjj][m_currentCentIndex][pmIndexOfHistoToGrab]->ProjectionY(Form("param_proj_%s_%s_%s_Cent%02d_Rap%02d",
           paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin),
           m_currentRapBin,m_currentRapBin);

        functionNames[jjj] = Form("param_parab_%s_%s_%s_Cent%02d_Rap%02d",paramNames[jjj].c_str(),
          m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
          m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin);
        param1DHistos[jjj]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),paramSymbols[jjj].c_str()));
      }
    }

    int lowBin = param1DHistos[0]->FindBin(lowMtM0);
    int midBin = param1DHistos[0]->FindBin(midMtM0);
    int highBin =param1DHistos[0]->FindBin(highMtM0);

    bool properLowBinFound = true;
    if(param1DHistos[0]->GetBinError(lowBin) == 0.0){
      properLowBinFound = false;
      for(int binX = lowBin; binX < highBin; binX++){
        if(param1DHistos[0]->GetBinError(binX) > 0.0){
          properLowBinFound = true;
          lowBin = binX;
          break;
        }
      }
    }

    bool properHighBinFound = true;
    if(param1DHistos[0]->GetBinError(highBin) == 0.0){
      properHighBinFound = false;
      for(int binX = highBin; binX > lowBin; binX--){
        if(param1DHistos[0]->GetBinError(binX) > 0.0){
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
    lowMtM0 =  param1DHistos[0]->GetBinCenter(lowBin) - 0.25*param1DHistos[0]->GetBinWidth(lowBin);
    highMtM0 = param1DHistos[0]->GetBinCenter(highBin) + 0.25*param1DHistos[0]->GetBinWidth(highBin);
    midMtM0 = 0.5*(lowMtM0 + highMtM0);
    midBin = param1DHistos[0]->FindBin(midMtM0);
    midMtM0 = param1DHistos[0]->GetBinCenter(midBin);


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
      if(amplitudes[0] == 0.0 && amplitudes[1] == 0.0 &&amplitudes[2] == 0.0) continue;

      PhysMath::get3PointParabParams(midMtM0, parabolaParams[0], parabolaParams[1], parabolaParams[2],
               lowMtM0 ,amplitudes[0],midMtM0 ,amplitudes[1],highMtM0,amplitudes[2]);
      #ifdef _ZFITTER_DEBUG_
        cout << " Initializing parabola funciton with y0=" << parabolaParams[0] << " y'= " << parabolaParams[1] << " y'' = " << parabolaParams[2] << endl;
      #endif

      extrapFunctions[jjj] = PhysMath::parabolaFunct(functionNames[jjj].c_str(), parabolaParams[0], parabolaParams[1], parabolaParams[2],midMtM0);
      extrapFunctions[jjj]->SetRange(lowMtM0,highMtM0);
      double distanceScale = fabs(amplitudes[1] - amplitudes[0]);
      extrapFunctions[jjj]->SetParLimits(0,parabolaParams[0] - 10.0*distanceScale,parabolaParams[0] + 10.0*distanceScale);
      if(parabolaParams[1] < 0.0){
        extrapFunctions[jjj]->SetParLimits(1,150.0*parabolaParams[1],-100.0*parabolaParams[1]);
      }else{
        extrapFunctions[jjj]->SetParLimits(1,-100.0*parabolaParams[1],150.0*parabolaParams[1]);
      }
      if(parabolaParams[2] < 0.0){
        extrapFunctions[jjj]->SetParLimits(2,150.0*parabolaParams[2],-100.0*parabolaParams[2]);
      }else{
        extrapFunctions[jjj]->SetParLimits(2,-100.0*parabolaParams[2],150.0*parabolaParams[2]);
      }


      TF1* initialFunct = (TF1*) extrapFunctions[jjj]->Clone();
      initialFunct->SetLineColor(kYellow);
      initialFunct->SetName("initialFunct");

      int status = param1DHistos[jjj]->Fit(functionNames[jjj].c_str(),"","",lowMtM0,highMtM0);
      TF1* fitResult = param1DHistos[jjj]->GetFunction(functionNames[jjj].c_str());

      //Push the results to the histograms
      double mu = -999.999;
      if(fitResult){
        for(int binX = param1DHistos[jjj]->FindBin(lowMtM0); binX <= param1DHistos[jjj]->GetNbinsX(); binX++){
          double center = param1DHistos[jjj]->GetBinCenter(binX);
          double valueToPush =  fitResult->Eval(center);
          if(jjj == 3){
            mu = midMtM0 - (fitResult->GetParameter(1)/fitResult->GetParameter(2));
            if(center > mu && mu > lowMtM0) valueToPush = fitResult->Eval(mu);

          }
          if(m_inAllCentMode){
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][jjj][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX,valueToPush);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][jjj][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX, 0.0);
          }else{
            m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][partIndex][jjj][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX, valueToPush);
            m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][partIndex][jjj][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX, 0.0);
          }
        }
      }else{
        #ifdef _ZFITTER_DEBUG_
          cout << " ERROR: EXTRAPOLATION FAILED!!!!!!!!! " << endl;
        #endif
      }

      #ifdef _ZFITTER_DEBUG_
        HistogramUtilities::printHistoInfo(param1DHistos[jjj]);
      #endif

      HistogramUtilities::discardNanHistoData(param1DHistos[jjj]);
      TGraph* parabPoint = new TGraph();
      if(fitResult) parabPoint->SetPoint(0,midMtM0,fitResult->GetParameter(0));
      parabPoint->SetMarkerStyle(20);
      parabPoint->SetMarkerColor(kRed);


      // DRAW THE FITS
      fittingCanvas->cd();
      gStyle->SetOptFit(0111);
      HistogramUtilities::setAxisRanges(param1DHistos[jjj], lowMtM0,highMtM0);
      TLine* extraLine = NULL;
      if(mu != 999.999 && fitResult && mu < highMtM0 && mu > lowMtM0){
        extraLine = new TLine(mu,fitResult->Eval(mu), highMtM0, fitResult->Eval(mu));
        extraLine->SetLineColor(kRed);
      }
      param1DHistos[jjj]->SetStats(true);
      param1DHistos[jjj]->Draw("E");
      initialFunct->Draw("same");
      parabPoint->Draw("Psame");
      if(extraLine) extraLine->Draw("same");
      //gPad->SetLogy();
      gPad->SetRightMargin(0.4);
      if(m_inAllCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Param%02d_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
        m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),jjj,
        m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Param%02d_%s_in_%s_space_Cent%02d_Rap_%02d.png",m_imagePreDir.c_str(),
        m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),jjj,
        m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin));
      }


      delete param1DHistos[jjj];
      delete initialFunct;
      delete extrapFunctions[jjj];
      delete parabPoint;
      if(extraLine) delete extraLine;



    } // end parameter loop



    //############     FIT THE AMPLITUDE WITH AN EXPONENTIAL  ###############
    param1DHistos[0]->Scale(1.0,"width");


    double amplitudes[3];
    amplitudes[0] = param1DHistos[0]->GetBinContent(lowBin);
    amplitudes[1] = param1DHistos[0]->GetBinContent(midBin);
    amplitudes[2] = param1DHistos[0]->GetBinContent(highBin);


    double bSeed = log(amplitudes[2]/amplitudes[0])/(highMtM0 - lowMtM0);
    double ampSeed = amplitudes[0]*exp(-bSeed*lowMtM0);

    #ifdef _ZFITTER_DEBUG_
      cout << " bSeed: " << bSeed << "    amp seed: " << ampSeed << endl;
    #endif

    TF1* ampExtractFunct = PhysMath::getExponentialFunction(functionNames[0].c_str(), ampSeed, bSeed, 0, 0);
    ampExtractFunct->SetRange(lowMtM0,highMtM0);
    ampExtractFunct->SetParLimits(0,0.0,10.0*ampSeed);
    if(bSeed < 0.0){
      ampExtractFunct->SetParLimits(1,15.0*bSeed,-10.0*bSeed);
    }else{
      ampExtractFunct->SetParLimits(1,-10.0*bSeed,15.0*bSeed);
    }
    ampExtractFunct->FixParameter(2,0.0);
    ampExtractFunct->FixParameter(3,0.0);



    int status = param1DHistos[0]->Fit(functionNames[0].c_str(),"","",lowMtM0,highMtM0);
    TF1* fitResult = param1DHistos[0]->GetFunction(functionNames[0].c_str());

    //Push the results to the histograms
    if(fitResult){
      for(int binX = param1DHistos[0]->FindBin(lowMtM0); binX <= param1DHistos[0]->GetNbinsX(); binX++){
        double center = param1DHistos[0]->GetBinCenter(binX);
        if(m_inAllCentMode){
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX, fitResult->Eval(center) * param1DHistos[0]->GetBinWidth(binX));
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX, 0.0);

        }else{
          m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][partIndex][0][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX, fitResult->Eval(center) * param1DHistos[0]->GetBinWidth(binX));
          m_Fit_Data_Cent_ZbTOF[m_currentPartIndex][partIndex][0][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX,0.0);
        }
      }
    }else{
      #ifdef _ZFITTER_DEBUG_
        cout << " ERROR: EXTRAPOLATION FAILED!!!!!!!!! " << endl;
      #endif
    }

    #ifdef _ZFITTER_DEBUG_
      HistogramUtilities::printHistoInfo(param1DHistos[0]);
    #endif

    HistogramUtilities::discardNanHistoData(param1DHistos[0]);


    // DRAW THE FITS
    fittingCanvas->cd();
    gStyle->SetOptFit(0111);
    //param1DHistos[0]->GetXaxis()->SetRangeUser(lowMtM0,highMtM0);
    HistogramUtilities::setAxisRanges(param1DHistos[0], lowMtM0,highMtM0);
    param1DHistos[0]->SetStats(true);
    param1DHistos[0]->Draw("E");
    //gPad->SetLogy();
    gPad->SetRightMargin(0.4);
    if(m_inAllCentMode){
      fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Param%02d_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),0,
      m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/%s/BTOFParamExtrapolations/Param%02d_%s_in_%s_space_Cent%02d_Rap_%02d.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),0,
      m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin));
    }





     m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex] = true;
  } // end particle loop



  return;
}









//============================================================================================================
//============================================================================================================



void ZFitter::fitBTOF(bool a_saveImages){
  m_inAllCentMode = true;
  cout << "####################################################################################################" << endl;
  cout << "                        FITTING Z_BTOF ALL CENTRALITIES  " << endl;
  cout << "####################################################################################################" << endl;

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);
  fittingCanvas->cd();
  gPad->SetRightMargin(0.4);
  gPad->SetLogy();

  m_currentCentIndex = -1;
  int numHorBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsY();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    for(int iii  = 0; iii < 9; iii++){
      //m_partHasBeenExtrapolated[iii] = false;
      m_partHasBeenExtrapolatedCharge[iii][0] = false;
      m_partHasBeenExtrapolatedCharge[iii][1] = false;
    }

    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    // make testing do the same bin as dEdx so
    int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    //if(m_testing && (dEdxBin%(numDeDxBins/m_divTesting) != 0 || m_currentRapidity == m_cuts->getYCM())) continue;
    //if(m_testing && (dEdxBin%(numDeDxBins/m_divTesting) != 0 ? (float)m_currentRapidity != m_cuts->getYCM() : 0 )) continue;
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }
      if(!doThisRapBin) continue;
    }

    if(m_rapBinsToTest.size() > 0){
      doThisRapBin = false;
      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }
      if(!doThisRapBin) continue;
    }
    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      if(m_testing && VertBinIndex%m_divTestingMtM0!=0) continue;

      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      //cout << " :E" << HorBinIndex << "P" << VertBinIndex << ": " << endl;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }





      for (int pmIndex=0; pmIndex <= 1; pmIndex++){
        m_currentPlusMinusIndex = pmIndex;
        if( !m_separatePM_bToF_allCent ){
          if(m_currentPlusMinusIndex == 1) continue;
          m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->ProjectionZ(
              Form("BTOF_%02d_%02d_PlusMinusTogether",m_currentRapBin, m_currentMtM0Bin),
              m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
          m_currentPartCharge = 0;
        }else{
          if(m_seedFromOppositeChargeDirection[m_currentPartIndex] < 0){
            if(pmIndex == 0) m_currentPlusMinusIndex = 1;
            if(pmIndex == 1) m_currentPlusMinusIndex = 0;
          }
          if(m_currentPlusMinusIndex == 0){
            m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Plus_allCent[m_currentPartIndex]->ProjectionZ(
              Form("BTOF_%02d_%02d_Plus",m_currentRapBin, m_currentMtM0Bin),
              m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
            m_currentPartCharge = 1;
          }else{
            m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Minus_allCent[m_currentPartIndex]->ProjectionZ(
              Form("BTOF_%02d_%02d_Minus",m_currentRapBin, m_currentMtM0Bin),
              m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
            m_currentPartCharge = -1;
          }
        }

        bool isBeyondAllBTOFExtrapRanges = true;
        for(int partIndex = 0; partIndex < 9; partIndex++){
          if(m_BTOF_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex] < 900.0   // there was actually a setting put into the matrix for this particle
             &&m_currentMomOfBin < m_BTOF_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]){
            isBeyondAllBTOFExtrapRanges = false;
          }
        }
        if(isBeyondAllBTOFExtrapRanges){
          #ifdef _ZFITTER_DEBUG_
            cout << " the momentum is beyond all high ranges for the parameter extrapolation, so continuing" << endl;
          #endif
          continue;
        }

        cout << "###########################################################################" << endl;
        cout << "Working on All Cent ZBTOF " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
        m_currentNumPartToFit = 0;
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
        double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
        if((countsInHisto < 100 && !m_testing) || countsInHisto < 40){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
          #endif
          continue;
        }
        //HistogramUtilities::removeHighCountCorruptData(m_currentHistoToFit);
        m_currentIntegralofHistogram = countsInHisto;
        /*TH1D* tempHistoPtr = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,5.0);
        delete tempHistoPtr;
        */

        m_currentHistoToFit->Scale(1.0,"width");  //make differential

        double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
        double lowFunctRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highFunctRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);



        //set histogram titles
        m_currentHistoToFit->SetTitle(Form("BTOF Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
          m_currentMomOfBin
          ));

        if(m_massSquaredMode){
          m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
        }else{
          m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
        }


        //################# SET UP PARTICLE TAGS ####################
        bool atLeastOneToFit = false;
        //##########################################################
        //##########    FIXED TARGET PARTICLE TAGS   ###############
        //##########################################################
        if(m_cuts->getEventConfig() == "FixedTarget"){
        //cout << "setting particle tags" << endl;
        //old particle tags
        /*
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            m_currentParticleTags[partIndex] = 0;
            if(partIndex == 3 || partIndex == 8){
              m_currentParticleTags[partIndex] = 0;
              continue;
            }
            double bToFStop = m_cuts->getMomStops_bToF(partIndex,m_currentPartCharge);
            cout << "partIndex: " << partIndex << "        bToF Stop: " << bToFStop << "        m_currentMomOfBin: " << m_currentMomOfBin << endl;
            if(bToFStop < m_currentMomOfBin){
              m_currentParticleTags[partIndex] = 0;
              continue;
            }

            //Extrapolate Means and Widths if going under another peak
            if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
              if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);

              if(m_partHasBeenExtrapolated[partIndex]){
                m_currentParticleTags[partIndex] = 1;
                continue;
              }
            }

          //see if there is data at prediction
            if(!m_initialParticleToggles[partIndex]){
              m_currentParticleTags[partIndex] = 0;
            }else{
              double invBetaPredict = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin) - m_invBetaCurves[m_currentPartIndex]->Eval(m_currentMomOfBin);
              int invBetaBin = m_currentHistoToFit->FindBin(invBetaPredict);
              double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*
              m_currentHistoToFit->GetBinWidth(invBetaBin);
              if(countsInThatBin > 2 && invBetaPredict < highRange && invBetaPredict > lowRange){
                bool zoomedFitSuccess = false;

                */
              /*if(m_currentMomOfBin < m_bTOFSeparatedStops[partIndex]){
                vector<double> fitData = fitIsolation_BTOF(partIndex);
                if(fitData[0] > 0){
                  cout << "   Zoom Fit Success: N: " << fitData[0] << " mu:" << fitData[1] << " sig:" << fitData[2] << endl;
                  m_InvBetaBTOFFitData[m_currentPartIndex][partIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
                  m_InvBetaBTOFFitData[m_currentPartIndex][partIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
                  m_InvBetaBTOFFitData[m_currentPartIndex][partIndex][2]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
                  m_currentParticleTags[partIndex] = 3;
                  zoomedFitSuccess = true;
                  continue;
                }
              }*/
            /*
                if(!zoomedFitSuccess){
                m_currentParticleTags[partIndex] = 2; //fit amps, means and widths
                atLeastOneToFit = true;
              }

            }else{
              m_currentParticleTags[partIndex] = 0;
            }
            */
            /*if(partIndex == 3 &&
              m_electronParams[2] - m_invBetaCurves->Eval(m_currentMomOfBin) > m_currentHistoToFit->GetBinLowEdge(1)){
              m_soloInvBetaAmpHistos[partIndex]->SetBinContent(HorBinIndex,VertBinIndex, 1.0);
              m_soloInvBetaMeanHistos[partIndex]->SetBinContent(HorBinIndex,VertBinIndex,
                m_electronParams[2]- m_invBetaCurves->Eval(m_currentMomOfBin));
              m_soloInvBetaStdDevHistos[partIndex]->SetBinContent(HorBinIndex,VertBinIndex, m_electronParams[3]);
              m_currentParticleTags[partIndex] = 1;
            }*/
         // }//end else of init toggles
        //}//end particle loop

          //Matt's particle tags
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            m_currentParticleTags[partIndex] = 0;
            #ifdef _ZFITTER_DEBUG_
              cout << "partIndex " << partIndex << ": " << endl;
            #endif
            if(partIndex == 3 || partIndex == 8)  continue;
            if( (m_currentPartCharge == 1 && m_partBTOFMomStopMatrix_Plus[m_currentPartIndex][partIndex] < m_currentMomOfBin) || (m_currentPartCharge == -1 && m_partBTOFMomStopMatrix_Minus[m_currentPartIndex][partIndex] < m_currentMomOfBin) ){

              #ifdef _ZFITTER_DEBUG_
                cout << "momentum of bin past stop: " << m_currentMomOfBin << ">";
                if(m_currentPartCharge == 1) cout << m_partBTOFMomStopMatrix_Plus[m_currentPartIndex][partIndex] << endl;
                else if(m_currentPartCharge == -1 ) cout << m_partBTOFMomStopMatrix_Minus[m_currentPartIndex][partIndex] << endl;
              #endif

              continue;
            }
            if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex][m_currentPlusMinusIndex]
               && m_currentPartIndex != partIndex){

              #ifdef _ZFITTER_DEBUG_
                cout << "momentum of bin " << m_currentMomOfBin << " is lower than bToF separated stops " << m_bTOFSeparatedStops[m_currentPartIndex][m_currentPlusMinusIndex] << endl;
              #endif

              continue;
          }
            if( (m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
                 || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex])
                 || (m_currentPartCharge == 0 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) ){

            #ifdef _ZFITTER_DEBUG_
              cout << "particle " << partIndex << " is not set to be fit in particle " << m_currentPartIndex << " space" << endl;
            #endif

              continue;
            }

          //Extrapolate Means and Widths if going under another peak
          if( (m_currentPartCharge == 1 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
            || (m_currentPartCharge == 0 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
            || (m_currentPartCharge == -1 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex]) ){
            if(!m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]){
              extrapolateMeanAndWidthParams_BTOF(partIndex);
            }

            if(m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]){
              m_currentParticleTags[partIndex] = 1;
              continue;
            }
          }

          //see if there is data at prediction
            double invBetaPredict = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin)
                                    - PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), m_currentMomOfBin);
            if(m_massSquaredMode) invBetaPredict = pow(m_partInfo->GetParticleMass(partIndex),2.0);
            int invBetaBin = m_currentHistoToFit->FindBin(invBetaPredict);
            double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin);
            if(countsInThatBin > 4 && invBetaPredict < highRange && invBetaPredict > lowRange){
            m_currentParticleTags[partIndex] = 2; //fit amps, means and widths
            atLeastOneToFit = true;
          }

        }//end particle loop

        //#######   FROM SEPARATE BTOF PLUS MINUS FIT CODE - WHICH WAS DELETED ####################
        /*
        bool atLeastOneToFit = false;
        //cout << "setting particle tags" << endl;

        for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
          m_currentParticleTags[partIndex] = 0;
          if(partIndex == 3 || partIndex == 8)  continue;
          if(m_cuts->getMomStops_bToF(partIndex,m_currentPartCharge) < m_currentMomOfBin){
            continue;
          }
          if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex]
             && m_currentPartIndex != partIndex){
            continue;
          }
          if( (m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) ){
            continue;
          }

            //Extrapolate Means and Widths if going under another peak
          if( (m_currentPartCharge == 1 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]) || (m_currentPartCharge == -1 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex]) ) {
            if(!m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]){
              extrapolateMeanAndWidthParams_BTOF(partIndex);
            }

            if(m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]){
              m_currentParticleTags[partIndex] = 1;
              continue;
            }
          }

            //see if there is data at prediction
          double invBetaPredict = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin)
          - PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), m_currentMomOfBin);
          int invBetaBin = m_currentHistoToFit->FindBin(invBetaPredict);
          double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)
          *m_currentHistoToFit->GetBinWidth(invBetaBin);
          if(countsInThatBin > 2 && invBetaPredict < highRange && invBetaPredict > lowRange){
              m_currentParticleTags[partIndex] = 2; //fit amps, means and widths
              atLeastOneToFit = true;
            }

          }//end particle loop
          */



      }else if(m_cuts->getEventConfig() == "ColliderCenter"){

        //##########################################################
        //##########    COLLIDER CENTER PARTICLE TAGS   ############
        //##########################################################
        for(int partIndex = 0; partIndex < m_numParticles; partIndex++){

          m_currentParticleTags[partIndex] = 0;
          if(partIndex == 3 || partIndex == 8)  continue;
          //if(m_cuts->getMomStops_bToF(partIndex,m_currentPartCharge) < m_currentMomOfBin) continue;
          if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex][m_currentPlusMinusIndex]
               && m_currentPartIndex != partIndex) continue;
          if( (m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
                 || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex])
                 || (m_currentPartCharge == 0 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) ) continue;
          if(m_currentPartCharge == 1 && m_partBTOFMomStopMatrix_Plus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;
          if(m_currentPartCharge == -1 && m_partBTOFMomStopMatrix_Minus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;


          //Extrapolate Means and Widths if going under another peak
          if( (m_currentPartCharge == 1 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
            || (m_currentPartCharge == 0 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
            || (m_currentPartCharge == -1 && m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex]) ){
            if(!m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]){
              extrapolateMeanAndWidthParams_BTOF(partIndex);
            }

            if(m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]){
              m_currentParticleTags[partIndex] = 1;
              continue;
            }
          }

          //see if there is data at prediction
          double invBetaPredict = m_invBetaCurves[partIndex]->Eval(m_currentMomOfBin)
                                  - PhysMath::invBeta(m_partInfo->GetParticleMass(m_currentPartIndex), m_currentMomOfBin);
          if(m_massSquaredMode) invBetaPredict = pow(m_partInfo->GetParticleMass(partIndex),2.0);
          int invBetaBin = m_currentHistoToFit->FindBin(invBetaPredict);
          double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)
                                   *m_currentHistoToFit->GetBinWidth(invBetaBin);
          if(countsInThatBin > 2 && invBetaPredict < highRange && invBetaPredict > lowRange){
            m_currentParticleTags[partIndex] = 2; //fit amps, means and widths
            atLeastOneToFit = true;
          }

        }//end particle loop
      } // END COLLIDER CENTER TAGS
      //########## END PARTICLE TAG SETUP ##############











       m_currentNumPartToFit = 0;
       for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
         if(m_currentParticleTags[partIndex] != 0){
          m_currentNumPartToFit++;
         }
       }


      if(!atLeastOneToFit){
        cout << "NO PARTICLES TO FIT" << endl;
        continue;
      }

      if(m_currentParticleTags[m_currentPartIndex] == 0){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping this bin beacuse main particle tag = 0 (not fit)" << endl;
        #endif
        //Drawing Failure
        fittingCanvas->cd();
        m_currentHistoToFit->Draw("E");
        DrawPredictionLines_TOF(fittingCanvas);
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/Rap_%02d_mTm0_%02d______NOT_FIT.png",
        m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), HorBinIndex, VertBinIndex));
        delete m_currentHistoToFit;
        continue;
      }


      //############ Truncate the Histogram to only the particles bing fit #############
      getRangeZTOF(lowRange, highRange);
      if(m_currentNumPartToFit == 1){
        double x1 = - 0.011;
        double x2 = 0.0;
        double x3 = 0.011;
        if(m_massSquaredMode){
          x2 = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0);
          x1 = 0.5*(x2 + m_cuts->getLowMassSquaredCut(m_currentPartIndex));
          //if(x1<0.0) x1 = 0.0;
          x3 = 0.5*(x2+m_cuts->getHighMassSquaredCut(m_currentPartIndex));
        }
        /*   Doesn't work with variable bins, go old school...
        double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(x1), 2);
        double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(x2), 2);
        double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(x3), 2);
        */
        double y1  = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(x1));
        double y2  = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(x2));
        double y3  = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(x3) );
        double N_3p, mean_3p, sig_3p = 0.0;
        PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,x1,y1,x2,y2,x3,y3);
        #ifdef _ZFITTER_DEBUG_
          cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
          cout << " sample points (" << x1  << "," <<  y1 <<  ") ("  << x2  << ","  << y2  << ") ("  << x3  << ","  << y3 <<  ")" << endl;
        #endif
        if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
          lowRange  = mean_3p  - 3.0*sig_3p;
          highRange  = mean_3p + 3.0*sig_3p;
          #ifdef _ZFITTER_DEBUG_
            cout << "     After: " << lowRange << " to " << highRange << endl;
          #endif
          if(m_currentPartIndex == 0 && m_currentMomOfBin < 0.58){
            //if pion at low momentum, avoid the muon and electron peaks
            lowRange  = mean_3p  - sig_3p;
            highRange  = mean_3p + 3.0*sig_3p;
          }
        }
        #ifdef _ZFITTER_DEBUG_
          else{
            cout << " not changing fit range by 3pt gaus, initializtion was wrong somewhere" << endl;
          }
        #endif
      }


      m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
      m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
      m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
      m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);

      lowFunctRange = lowRange;
      highFunctRange = highRange;
      TH1D* ptrToDelete = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
      delete ptrToDelete;
      if(!m_currentHistoToFit){
        cout << "  ERROR: No Ptr after truncation of histogram" << endl;
        continue;
      }

      double truncatedTotalCounts = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
      double possibleNewBins = truncatedTotalCounts/6.0;
      if(possibleNewBins > 80.0) possibleNewBins = 80.0;
      TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(m_currentHistoToFit,(int) possibleNewBins,true,true);
      delete m_currentHistoToFit;
      m_currentHistoToFit = rebinHisto;
      HistogramUtilities::removeHighCountCorruptData(m_currentHistoToFit);


      m_currentHistoToFit->SetMarkerStyle(20);
      m_currentHistoToFit->SetMarkerColor(kBlack);






      //#########  INITIALIZE PARAMETERS #############
      if(m_useStudentTDistributionsForTOF)  initBTOFParams_student();
      else initBTOFParams();

      if(m_currentNumPartToFit == 0 || areBTOFParamsRidiculous()){
        #ifdef _ZFITTER_DEBUG_
          if(m_currentNumPartToFit == 0) cout << "Skipping this bin beacuse there are no particles to fit" << endl;
          else cout << "\n\n\n\n\n\n####!!!!!!!!!#######   BTOF parameters are ridiculous #####!!!!!!!!!!########### \n\n\n\n\n" << endl;
        #endif
        //Drawing Failure
        fittingCanvas->cd();
        m_currentHistoToFit->Draw("E");
        DrawPredictionLines_TOF(fittingCanvas);
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/Rap_%02d_mTm0_%02d______NOT_FIT.png",
        m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), HorBinIndex, VertBinIndex));
        delete m_currentHistoToFit;
        continue;
      }

      //cout << "Fitting " << m_currentNumPartToFit << " Particles" << "   With " << m_currentFreeParameters << " Parameters" << endl;
      #ifdef _ZFITTER_DEBUG_
        cout << "Fitting " << m_currentNumPartToFit << " Particles" << "   With " << m_currentFreeParameters << " Parameters" << endl;
        cout << "Initial Parameters: " << endl;
        printCurrentParams();
      #endif      //printCurrentParams();
      //do the fit


      //######   DEFINE THE FIT FUNCTION #########
      TF1* fitFunct = NULL;
      if(m_useStudentTDistributionsForTOF)  fitFunct = new TF1("fitFunctInvBeta",this,&ZFitter::ZTOFOperator_student,
        lowFunctRange,highFunctRange,m_currentFreeParameters,"ZFitter","ZTOFOperator_student");
      else fitFunct = new TF1("fitFunctInvBeta",this,&ZFitter::ZTOFOperator,lowFunctRange,highFunctRange,m_currentFreeParameters,"ZFitter","ZTOFOperator");
      fitFunct->SetNpx(500);

      //########### INITIALIZE FUNCTION WITH INITIAL PARAMETERS ##################
      if(m_useStudentTDistributionsForTOF)  initBTOFFunction_student(fitFunct);
      else initBTOFFunction(fitFunct);
      TF1* fitFunctInitial = (TF1*) fitFunct->Clone();
      fitFunctInitial->SetLineColor(kYellow);
      if(m_drawInitialSeedsToFits) fitFunctInitial->SetName("fitFunctInitial");




      //############  EXECUTE THE FIT ##################
      string printSetting = "EQ";
      #ifdef _ZFITTER_DEBUG_
        printSetting = "EV";
      #endif
      if(m_doIntegralOfBinForFits) printSetting += "I";
      m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
      int status = m_currentHistoToFit->Fit("fitFunctInvBeta",printSetting.c_str(),"",lowFunctRange,highFunctRange);//add Q to quiet
      TF1* fitResult = m_currentHistoToFit->GetFunction("fitFunctInvBeta");


      m_FitStatus_ZbTOF[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status+5);


      #ifdef  _ZFITTER_DEBUG_
        cout << "Final Status: " << status << endl;
      #endif

      if(status!=0){
        cout << "\n\n !!! BAD FIT !!! Status = " << status << " \n\n";
        if(status < 0) continue;
      }

      if(!fitResult || (status !=0 && !m_storeFitsNotStatusZero)){
        #ifdef  _ZFITTER_DEBUG_
          cout << "Fit failed, see parameters:" << endl;
          printCurrentParams();
        #endif

        //################   DRAW THE FIT IMAGE ###################
        if(a_saveImages){
          fittingCanvas->cd();
          m_currentHistoToFit->Draw("E");
          //gPad->SetLogx();
          gPad->SetLogy();
          //fitResult->SetLineColor(kRed);
          if(m_drawInitialSeedsToFits) fitFunctInitial->Draw("same");
          //fitResult->Draw("same");
          vector<TF1*> singleFuncts;
          if(m_useStudentTDistributionsForTOF) singleFuncts = drawIndividualFits_Student();
          else singleFuncts = drawIndividualFits_Gaus();

          DrawPredictionLines_TOF(fittingCanvas);

          gPad->Update();
          gPad->SetRightMargin(0.4);
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/RapBin_%02d_mTm0Bin_%02d_________BAD_FIT.png",
            m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), HorBinIndex, VertBinIndex));

          gPad->SetLogy(false);
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/NoLog_RapBin_%02d_mTm0Bin_%02d_________BAD_FIT.png",
            m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), HorBinIndex, VertBinIndex));

          for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
          for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
          m_TLinePtrs.clear();

        }
        delete m_currentHistoToFit;
        //delete fitResult;
        delete fitFunctInitial;
        delete fitFunct;
        m_previousGoodBTOFParams.clear();

        continue;

      }




      vector<double> ampErrors;
      vector<double> meanErrors;
      vector<double> stdDevErrors;
      vector<double> ndfErrors;

      int paramIndexCounter = 0;
      //cout << "Final Parameters::";
      if(fitResult){
        if(m_useStudentTDistributionsForTOF){
          for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
            if(iii%4==0) cout << endl;
            if(m_paramToggles[iii]){
              m_currentParameters[iii] = fitResult->GetParameter(paramIndexCounter);
              if(paramIndexCounter%4==0) ampErrors.push_back(fitResult->GetParError(paramIndexCounter) );
              if(paramIndexCounter%4==1) meanErrors.push_back(fitResult->GetParError(paramIndexCounter));
              if(paramIndexCounter%4==2) stdDevErrors.push_back(fitResult->GetParError(paramIndexCounter));
              if(paramIndexCounter%4==3) ndfErrors.push_back(fitResult->GetParError(paramIndexCounter));
              paramIndexCounter++;
            }
            //cout << m_currentParameters[iii] << "  ";
          }
        }else{
          for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
            if(iii%3==0) cout << endl;
            if(m_paramToggles[iii]){
              m_currentParameters[iii] = fitResult->GetParameter(paramIndexCounter);
              if(paramIndexCounter%3==0) ampErrors.push_back(fitResult->GetParError(paramIndexCounter));
              if(paramIndexCounter%3==1) meanErrors.push_back(fitResult->GetParError(paramIndexCounter));
              if(paramIndexCounter%3==2) stdDevErrors.push_back(fitResult->GetParError(paramIndexCounter));
              paramIndexCounter++;
            }
            //cout << m_currentParameters[iii] << "  ";
          }
        }
      }

      if(status == 0){
        m_previousGoodBTOFParams.clear();
        for(int paramIndex = 0; paramIndex < m_currentParameters.size(); paramIndex++){
          m_previousGoodBTOFParams.push_back(m_currentParameters[paramIndex]);
        }
      }

      //cout << endl;

      m_ChiSqr_ZbTOF[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitResult->GetChisquare() / fitResult->GetNDF());
      #ifdef  _ZFITTER_DEBUG_
        cout << "   NumBins: " << m_currentHistoToFit->GetNbinsX() << " Params: " << fitResult->GetNumberFreeParameters()
             << "  ndof: " <<  fitResult->GetNDF() << endl;
      #endif

      //push to particle fits
      int paramIndex = 0; //the particles not included in fit
      for(int partIndex = 0;status == 0 && partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] == 0) continue;
        if(m_currentParticleTags[partIndex] ==1){
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
          //m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ampErrors[paramIndex/3]);
          if(m_useStudentTDistributionsForTOF){
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ampErrors[paramIndex/4] * m_currentIntegralofHistogram);
            paramIndex+=4;
          }
          else{
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ampErrors[paramIndex/3] * m_currentIntegralofHistogram);
            paramIndex+=3;
          }
        }else if(m_currentParticleTags[partIndex] == 2){
          //if(m_currentParameters[paramIndex+2] <=0.0 ) cout << "Pushing Back Zero Sig" << endl;
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+1]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+2]);
          if(m_useStudentTDistributionsForTOF){
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+3]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ampErrors[paramIndex/4] * m_currentIntegralofHistogram);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, meanErrors[paramIndex/4]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, stdDevErrors[paramIndex/4]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ndfErrors[paramIndex/4]);
            paramIndex+=4;
          }else{
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ampErrors[paramIndex/3] * m_currentIntegralofHistogram);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, meanErrors[paramIndex/3]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, stdDevErrors[paramIndex/3]);
            paramIndex+=3;
          }
        }else{
          //if(m_currentParameters[paramIndex+2] <=0.0 ) cout << "Pushing Back Zero Sig" << endl;
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+1]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+2]);
          if(m_useStudentTDistributionsForTOF){
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+3]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ampErrors[paramIndex/4] * m_currentIntegralofHistogram);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, meanErrors[paramIndex/4]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, stdDevErrors[paramIndex/4]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ndfErrors[paramIndex/4]);
            paramIndex+=4;
          }else{
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, ampErrors[paramIndex/3] * m_currentIntegralofHistogram);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, meanErrors[paramIndex/3]);
            m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, stdDevErrors[paramIndex/3]);
            paramIndex+=3;
          }
        }
      }

      //################   DRAW THE FIT IMAGE ###################
      if(a_saveImages){
        fittingCanvas->cd();
        m_currentHistoToFit->Draw("E");
        //gPad->SetLogx();
        gPad->SetLogy();
        fitResult->SetLineColor(kRed);
        if(m_drawInitialSeedsToFits) fitFunctInitial->Draw("same");
        fitResult->Draw("same");
        vector<TF1*> singleFuncts;
        if(m_useStudentTDistributionsForTOF) singleFuncts = drawIndividualFits_Student();
        else singleFuncts = drawIndividualFits_Gaus();

        DrawPredictionLines_TOF(fittingCanvas);

        gPad->Update();
        gPad->SetRightMargin(0.4);
        gSystem->ProcessEvents();
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/RapBin_%02d_mTm0Bin_%02d.png",
          m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), HorBinIndex, VertBinIndex));

        gPad->SetLogy(false);
        gPad->Update();
        gSystem->ProcessEvents();
        fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits/NoLog_RapBin_%02d_mTm0Bin_%02d.png",
          m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), HorBinIndex, VertBinIndex));

        for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
        for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
        m_TLinePtrs.clear();
      }
      delete m_currentHistoToFit;
      //delete fitResult;
      delete fitFunctInitial;
      delete fitFunct;









    } // end plus minus loop

  }//end pt loop
  m_currentPlusMinusIndex = 0;
  m_currentPartCharge = 1;
  extrapolateMeanAndWidthParams_BTOF_m2();
  m_currentPlusMinusIndex = 1;
  m_currentPartCharge = -1;
  extrapolateMeanAndWidthParams_BTOF_m2();

  } //end rap loop



  ////draw the results
  //for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
  //  drawPartdEdxFits(partIndex);
  //}
  delete fittingCanvas;
}


void ZFitter::fitBTOF_withBackground(int a_particleIndex, vector<double> a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT , bool a_useExpBackground, bool a_fitBackgroundOnRight, double a_lowMassOverRide, double a_highMassOverRide){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING BTOF YIELD ALL CENTRALITIES - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;


  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = -1;
  m_inAllCentMode = true;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if(!m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]) cout << " Error in ZFitter::fitBTOF_withBackground, can't find m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  TH2D* base_histogram = dynamic_cast<TH2D*>(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->Clone());
  HistogramUtilities::zeroBins2D(base_histogram);
  base_histogram->SetName("BaseHistogram_ForReference");
  base_histogram->SetTitle("BaseHistogram_ForReference;y;m_{T}-m_{0}");

  cout << "  Starting to form Base Histograms" << endl;
  for(int partIndex=0; partIndex<9; partIndex++){

    if(m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]){
      //Background Amplitude
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetName(Form("%s Space %s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      //
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetName(Form("%s Space %s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b2",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetName(Form("%s Space %s BTOF Background b2",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetName(Form("%s Space %s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetName(Form("%s Space %s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
    }

    if(m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]){
      //Background Amplitude
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetName(Form("%s Space %s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      //
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetName(Form("%s Space %s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b2",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetName(Form("%s Space %s BTOF Background b2",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetName(Form("%s Space %s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetName(Form("%s Space %s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
    }
  }//end loop over partIndex

  delete base_histogram;

  m_currentPartCharge = 1;
  m_currentPlusMinusIndex = 0;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsY();
  cout << "numHorBins (rapidity): " << numHorBins << endl;
  cout << "numVertBins (mTm0): " << numVertBins << endl;
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;

      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      for(int partIndex=0; partIndex<9; partIndex++){
        if(partIndex == 3) continue;
        cout << "working in part " << m_currentPartIndex << " binning on partIndex " << partIndex << endl;
        if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;

        if(m_partBTOFMomStopMatrix_Plus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;
        m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Plus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


        double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
        cout << " low: " << lowRange << " high:" << highRange << endl;
        if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
        else{
          /*
          lowRange = pow(m_partInfo->GetParticleMass(partIndex)-0.1,2.0);
          if(lowRange < HistogramUtilities::getLowEdge(m_currentHistoToFit)) lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
          highRange = pow(m_partInfo->GetParticleMass(partIndex)+0.1,2.0);
          if(highRange > HistogramUtilities::getHighEdge(m_currentHistoToFit)) highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);*/

          getRangeZTOF_mass2(lowRange,highRange,partIndex);
        }
        cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
        if(highRange < lowRange){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }


        cout << "m_ZbTOF_LowRange["<< partIndex<<"]["<<m_currentPlusMinusIndex<<"]: " << m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex] << endl;


        /*
        m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
        m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
        m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        */
        if(m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]){
          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        }
        TH1D* origHistToFit = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
        delete origHistToFit;

        double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
        if(countsInHisto < 40){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }
        //TH1D* tempHistoPtr = m_currentHistoToFit;
        //m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
        //delete tempHistoPtr;

        cout << "m_currentHistoToFit: " << m_currentHistoToFit << endl;

        m_currentHistoToFit->Scale(1.0,"width");
        m_currentIntegralofHistogram = countsInHisto;
        m_currentHistoToFit->Scale(1.0/m_currentIntegralofHistogram);




        //set histogram titles
        m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
        partName.c_str(),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
        m_currentMomOfBin
        ));
        m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
        if(m_massSquaredMode){
          m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
        }


        //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
        //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;
        //Extrapolate Means and Widths if going under another peak
        //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
        //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
        //  if(m_partHasBeenExtrapolated[partIndex]){
        //    m_currentParticleTags[partIndex] = 1;
        //    continue;
        //  }
        //}


        double predictLocation = 0.0;
        cout << "m_currentHistoToFit: " << m_currentHistoToFit << endl;
        if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(partIndex),2.0 );
        int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
        double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin)*m_currentIntegralofHistogram;
        cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
        cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
        if(countsInThatBin < 1.0) continue;

        cout << " Starting the fit " << endl;
        vector<double> emptyVect;
        vector<double> fitData = fit_BTOF_withBackground(partIndex,emptyVect,a_momentum_startUsingBkgd,a_useStudentT,a_useExpBackground,a_fitBackgroundOnRight);
        cout << " out of fit" << endl;
        if(fitData[0] > 0){
          cout << "fit data size " << fitData.size() << endl;
          cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0] << endl;
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[12]);
          cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetTitle() << endl;
          cout << "here now doing backgournd data" << endl;

          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][0][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]*m_currentIntegralofHistogram);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]*m_currentIntegralofHistogram);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][1][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][2][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][3][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[16]);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][4][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[8]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[17]);
        }else{
          cout << "     Error:  fit_BTOF_withBackground : FAILED for Plus" << endl;
        }
      }//end partIndex loop
    }//end mtm0 loop
  } //end rap loop









  // #####################      FIT    MINUS       ######################


  if(!a_fitMinus) return;
  m_currentPartCharge = -1;
  m_currentPlusMinusIndex = 1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;

      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      for(int partIndex=0; partIndex<9; partIndex++){
        if(partIndex == 3) continue;
        if(!m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) continue;
        if(m_partBTOFMomStopMatrix_Minus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;
        m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Minus_allCent[m_currentPartIndex]->ProjectionZ(
          Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
          m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


        double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
        cout << " low: " << lowRange << " high:" << highRange << endl;
        if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
        else{
          /*
          lowRange = pow(m_partInfo->GetParticleMass(partIndex)-0.1,2.0);
          if(lowRange < HistogramUtilities::getLowEdge(m_currentHistoToFit)) lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
          highRange = pow(m_partInfo->GetParticleMass(partIndex)+0.1,2.0);
          if(highRange > HistogramUtilities::getHighEdge(m_currentHistoToFit)) highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);*/
          getRangeZTOF_mass2(lowRange,highRange,partIndex);

        }
        cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
        if(highRange < lowRange){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }

        /*
        m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
        m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
        m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        */

        if(m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]){
          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        }

        TH1D* origHistToFit = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
        delete origHistToFit;

        double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
        if(countsInHisto < 40){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }
        cout << " skipping because counts = " << countsInHisto << endl;
        //TH1D* tempHistoPtr = m_currentHistoToFit;
        //m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
        //delete tempHistoPtr;
        m_currentHistoToFit->Scale(1.0,"width");
        m_currentIntegralofHistogram = countsInHisto;
        m_currentHistoToFit->Scale(1.0/m_currentIntegralofHistogram);




        //set histogram titles
        m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
          partName.c_str(),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
          m_currentMomOfBin
          ));
        m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
        if(m_massSquaredMode){
          m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
        }

        //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
        //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;

        //Extrapolate Means and Widths if going under another peak
        //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
        //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
        //  if(m_partHasBeenExtrapolated[partIndex]){
        //    m_currentParticleTags[partIndex] = 1;
        //    continue;
        //  }
        //}


        double predictLocation = 0.0;
        if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(partIndex),2.0 );
        int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
        double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin)*m_currentIntegralofHistogram;
        cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
        cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
        if(countsInThatBin < 1.0) continue;

        cout << " Starting the fit " << endl;
        vector<double> constrainVect;
        if(a_constrainMinusToPlus){
          constrainVect.push_back(-999); // free the amplitude
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(-999); // free the scale of the background
          constrainVect.push_back(-999); // free the slope of the background
          //constrainVect.push_back(m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(-999);
          constrainVect.push_back(-999); // free the scale of the linear background
        }
        vector<double> fitData = fit_BTOF_withBackground(partIndex,constrainVect,a_momentum_startUsingBkgd,a_useStudentT,a_useExpBackground,a_fitBackgroundOnRight);
        cout << " out of fit" << endl;
        if(fitData[0] > 0){
          cout << "fit data size " << fitData.size() << endl;
          cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1] << endl;
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[12]);
          cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->GetTitle() << endl;
          cout << "here now doing backgournd data" << endl;
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]*m_currentIntegralofHistogram);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]*m_currentIntegralofHistogram);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[16]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[8]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[17]);
        }else{
          cout << "     Error:  fit_BTOF_withBackground : FAILED for Minus" << endl;
        }
        constrainVect.clear();
        #ifdef _CPP11_
        constrainVect.shrink_to_fit();
        #endif
      }//end partIndex loop
    }//end mtm0 loop
  } //end rap loop

}



void ZFitter::fitBTOF_NonMixing(int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT , double a_lowMassOverRide, double a_highMassOverRide){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING BTOF YIELD ALL CENTRALITIES - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;

  if(a_particleIndex < 0 || a_particleIndex == 3  || a_particleIndex > 7 ){
    cout << " Error in ZFitter::fitBTOF_NonMixing  with invalid argument, a_particleIndex" << endl;
    return;
  }
  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = -1;
  m_inAllCentMode = true;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if(!m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]) cout << " Error in ZFitter::fitBTOF_NonMixing, can't find m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  TH2D* base_histogram = dynamic_cast<TH2D*>(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->Clone());
  HistogramUtilities::zeroBins2D(base_histogram);
  base_histogram->SetName("BaseHistogram_ForReference");
  base_histogram->SetTitle("BaseHistogram_ForReference;y;m_{T}-m_{0}");

  cout << "  Starting to form Base Histograms" << endl;
  //Background Amplitude
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetName(Form("%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetName(Form("%s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetName(Form("%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetName(Form("%s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  //
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetName(Form("%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetName(Form("%s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetName(Form("%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetName(Form("%s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetName(Form("%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetName(Form("%s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetName(Form("%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetName(Form("%s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetName(Form("%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetName(Form("%s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetName(Form("%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetName(Form("%s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));


  delete base_histogram;

  m_currentPartCharge = 1;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsY();
  cout << "numHorBins (rapidity): " << numHorBins << endl;
  cout << "numVertBins (mTm0): " << numVertBins << endl;
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;

      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Plus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
      else{
        lowRange = pow(a_lowMassOverRide,2.0);
        highRange = pow(a_highMassOverRide,2.0);
      }
      cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }


      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);

      TH1D* origHistToFit = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
      delete origHistToFit;

      double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
      if(countsInHisto < 40){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");




      //set histogram titles
      m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
        partName.c_str(),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
        m_currentMomOfBin
        ));
      m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
      m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
      if(m_massSquaredMode){
        m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
      }


      //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
      //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;
      //Extrapolate Means and Widths if going under another peak
      //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
      //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
      //  if(m_partHasBeenExtrapolated[partIndex]){
      //    m_currentParticleTags[partIndex] = 1;
      //    continue;
      //  }
      //}

      double predictLocation = 0.0;
      if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0 );
      int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin);
      cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> emptyVect;
      vector<double> fitData = fitIsolation_BTOF_withBackground(emptyVect,a_useStudentT);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0] << endl;
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
      }else{
        cout << "     Error:  fitIsolation_BTOF_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop









  // #####################      FIT    MINUS       ######################


  if(!a_fitMinus) return;
  m_currentPartCharge = -1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;

      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Minus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
      else{
        lowRange = pow(a_lowMassOverRide,2.0);
        highRange = pow(a_highMassOverRide,2.0);
      }
      cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }

      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);


      TH1D* origHistToFit = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
      delete origHistToFit;

      double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
      if(countsInHisto < 40){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }
      cout << " skipping because counts = " << countsInHisto << endl;
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");




      //set histogram titles
      m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
        partName.c_str(),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
        m_currentMomOfBin
        ));
      m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
      m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
      if(m_massSquaredMode){
        m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
      }

      //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
      //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;

      //Extrapolate Means and Widths if going under another peak
      //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
      //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
      //  if(m_partHasBeenExtrapolated[partIndex]){
      //    m_currentParticleTags[partIndex] = 1;
      //    continue;
      //  }
      //}



      double predictLocation = 0.0;
      if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0 );
      int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin);
      cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;
      if(a_constrainMinusToPlus){
        constrainVect.push_back(-999); // free the amplitude
        constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999); // free the scale of the background
        constrainVect.push_back(-999); // free the slope of the background
        //constrainVect.push_back(m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999);
        constrainVect.push_back(-999); // free the scale of the linear background
      }
      vector<double> fitData = fitIsolation_BTOF_withBackground(constrainVect,a_useStudentT);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1] << endl;
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
      }else{
        cout << "     Error:  fitIsolation_BTOF_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop

}


void ZFitter::fitBTOF_NonMixing_Liggett(int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT , double a_lowMassOverRide, double a_highMassOverRide){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING BTOF YIELD ALL CENTRALITIES - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;

  if(a_particleIndex < 0 || a_particleIndex == 3  || a_particleIndex > 7 ){
    cout << " Error in ZFitter::fitBTOF_NonMixing_Liggett  with invalid argument, a_particleIndex" << endl;
    return;
  }
  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = -1;
  m_inAllCentMode = true;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if(!m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]) cout << " Error in ZFitter::fitBTOF_NonMixing_Liggett, can't find m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  TH2D* base_histogram = dynamic_cast<TH2D*>(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->Clone());
  HistogramUtilities::zeroBins2D(base_histogram);
  base_histogram->SetName("BaseHistogram_ForReference");
  base_histogram->SetTitle("BaseHistogram_ForReference;y;m_{T}-m_{0}");

  cout << "  Starting to form Base Histograms" << endl;
  //Background Amplitude
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetName(Form("%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetName(Form("%s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetName(Form("%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetName(Form("%s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  //
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetName(Form("%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetName(Form("%s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetName(Form("%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetName(Form("%s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetName(Form("%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetName(Form("%s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetName(Form("%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetName(Form("%s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetName(Form("%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetName(Form("%s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetName(Form("%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetName(Form("%s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));


  delete base_histogram;

  m_currentPartCharge = 1;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsY();
  cout << "numHorBins (rapidity): " << numHorBins << endl;
  cout << "numVertBins (mTm0): " << numVertBins << endl;
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;

      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Plus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
      else{
        highRange = pow(a_lowMassOverRide,2.0);
        lowRange = pow(a_highMassOverRide,2.0);
      }
      cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }

      if(!m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]){
          cout << "m_ZbTOF_LowRange not defined" << endl;
      }
      cout << "m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);" << endl;
      cout << m_currentPartIndex << " " << m_currentPlusMinusIndex << " " << m_currentRapBin << " "<< m_currentMtM0Bin << " " << lowRange << endl;

      /*double list[2] = {10,20}; // list is a double*
      double* ptr_to_first_element = &list[0];
      double second_element = ptr_to_first_element[1]; //dereferences double*

      TH2D*** doodad = m_ZbTOF_LowRange;
      TH2D* thingy = m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex];
      thingy->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);*/


      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
      cout << " m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);" << endl;
      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
      cout << "m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);" << endl;
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
      cout << "m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);" << endl;
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
      cout << "TH1D* origHistToFit = m_currentHistoToFit;" << endl;

      TH1D* origHistToFit = m_currentHistoToFit;
      cout << "m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);" << endl;
      m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
      cout << "delete origHistToFit;" << endl;
      delete origHistToFit;

      cout << "double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);" << endl;
      double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
      if(countsInHisto < 40){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }
      cout << "countsInHisto ok!" << endl;
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");




      //set histogram titles
      m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
        partName.c_str(),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
        m_currentMomOfBin
        ));
      m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
      m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
      if(m_massSquaredMode){
        m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
      }


      //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
      //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;
      //Extrapolate Means and Widths if going under another peak
      //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
      //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
      //  if(m_partHasBeenExtrapolated[partIndex]){
      //    m_currentParticleTags[partIndex] = 1;
      //    continue;
      //  }
      //}

      double predictLocation = 0.0;
      if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0 );
      int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin);
      cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> emptyVect;
      vector<double> fitData = fitIsolation_BTOF_withBackground(emptyVect,a_useStudentT);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0] << endl;
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
      }else{
        cout << "     Error:  fitIsolation_BTOF_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop









  // #####################      FIT    MINUS       ######################


  if(!a_fitMinus) return;
  m_currentPartCharge = -1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;

      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Minus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
      else{
        lowRange = pow(a_lowMassOverRide,2.0);
        highRange = pow(a_highMassOverRide,2.0);
      }
      cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }

      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
      m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
      m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);


      TH1D* origHistToFit = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
      delete origHistToFit;

      double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
      if(countsInHisto < 40){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
      }
      cout << " skipping because counts = " << countsInHisto << endl;
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");




      //set histogram titles
      m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
        partName.c_str(),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
        m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
        m_currentMomOfBin
        ));
      m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
      m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
      if(m_massSquaredMode){
        m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
      }

      //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
      //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;

      //Extrapolate Means and Widths if going under another peak
      //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
      //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
      //  if(m_partHasBeenExtrapolated[partIndex]){
      //    m_currentParticleTags[partIndex] = 1;
      //    continue;
      //  }
      //}



      double predictLocation = 0.0;
      if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0 );
      int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin);
      cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;
      if(a_constrainMinusToPlus){
        constrainVect.push_back(-999); // free the amplitude
        constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999); // free the scale of the background
        constrainVect.push_back(-999); // free the slope of the background
        //constrainVect.push_back(m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999);
        constrainVect.push_back(-999); // free the scale of the linear background
      }
      vector<double> fitData = fitIsolation_BTOF_withBackground(constrainVect,a_useStudentT);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1] << endl;
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][1]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
      }else{
        cout << "     Error:  fitIsolation_BTOF_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop

}




//========================================================================================================================================================================================================================
//========================================================================================================================================================================================================================







//fit as function of centrality using results from fitInvBeta
void ZFitter::fitBTOFCent(int a_centIndex, bool a_saveImages){

  cout << "Starting ZFitter::fitBTOFCent" << endl;

  m_inAllCentMode = false;
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);


  int numHorBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsY();
  //int numCentralities = m_numParticles; // cout << "Warning hardcoded numCentralities = m_numParticles" << endl;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    if(a_centIndex !=-1 && centIndex != a_centIndex) continue;
    m_currentCentIndex = centIndex;
    cout << "m_currentCentIndex: " <<  m_currentCentIndex << endl;
    for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

      m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
      cout << "m_currentRapidity: " << m_currentRapidity << endl;
      // make testing do the same bin as dEdx so
      int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
      int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
      //cout << "Rapidity: " << m_currentRapidity << endl;
      //if(m_testing && (dEdxBin%(numDeDxBins/m_divTesting) != 0 ? (float)m_currentRapidity != m_cuts->getYCM() : 0 )) continue;

      bool doThisRapBin = false;
      if(m_testing){
        if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
          doThisRapBin = false;
        }else{
          doThisRapBin = true;
        }
        if(!doThisRapBin) continue;
      }

      if(m_rapBinsToTest.size() > 0){
        doThisRapBin = false;
        for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
          if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
        }
        if(!doThisRapBin) continue;
      }


      //if(m_testing && !doThisRapBin) continue;
      for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
        if(m_testing && VertBinIndex%m_divTestingMtM0!=0) continue;

        m_currentRapBin = HorBinIndex;
        m_currentMtM0Bin = VertBinIndex;
        //cout << " :E" << HorBinIndex << "P" << VertBinIndex << ": " << endl;
        m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
        //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
        cout << "m_currentMtM0: " << m_currentMtM0 << endl;

        //#### PLUS MINUS HISTO LOOP #########
        for(int pmIndex= 0; pmIndex< 2; pmIndex++){

          TH3D* histo3DPM = NULL;
          string partChargeNames[2] = {"Positive","Negative"};
          string partChargeNamesSmall[2] = {"Plus","Minus"};
          if(pmIndex== 0){
            histo3DPM = m_rap_mTm0_invBetaBTOF_Plus[m_currentPartIndex][centIndex];
            m_currentPartCharge = 1;
          }else{
            histo3DPM = m_rap_mTm0_invBetaBTOF_Minus[m_currentPartIndex][centIndex];
            m_currentPartCharge = -1;
          }
          m_currentPlusMinusIndex = pmIndex;
          cout << "\n\n\n###########################################################################" << endl;
          cout << "Working on Specific Centrality ZBTOF " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " Cent: " << centIndex << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << endl;


          //################   SETUP HISTOGRAM TO FIT #####################
          m_currentHistoToFit = histo3DPM->ProjectionZ(
          Form("BTOF_%s_%02d_%02d",m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin,m_currentMtM0Bin), m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
          HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
          double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
          if((countsInHisto < 100 && !m_testing) || (countsInHisto < 20)) {
            cout << "not enough counts in hist" << endl;
            delete m_currentHistoToFit;
            #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts: " << countsInHisto << endl;
            #endif
            continue;
          }
          //HistogramUtilities::removeHighCountCorruptData(m_currentHistoToFit);
          m_currentIntegralofHistogram = countsInHisto;

          //TH1D* tempHistoPtr = m_currentHistoToFit;
          //m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,5.0);
          //delete tempHistoPtr;
          m_currentHistoToFit->Scale(1.0,"width");

          double lowRange = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistoToFit));
          double highRange = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistoToFit)+1 );
          double lowFunctRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
          double highFunctRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);




          m_currentHistoToFit->SetTitle(Form("BTOF %s Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] (GeV/c) p_{tot}=%2.3f",
            partChargeNames[pmIndex].c_str(),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
            m_currentMomOfBin
            ));
          m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
          if(m_massSquaredMode){
            m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
            m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
          }

          //###################### SETUP PARTICLE TAGS  ######################
          bool atLeastOne = false;
          if(m_cuts->getEventConfig() == "FixedTarget"){
              //##########################################################
              //##########    FIXED TARGET PARTICLE TAGS   ############
              //##########################################################

            //old particle tags
            /*
            for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
              double bToFStop = m_cuts->getMomStops_bToF(partIndex,m_currentPartCharge);
              if(bToFStop < m_currentMomOfBin){
                m_currentParticleTags[partIndex] = 0;
                continue;
              }
              //see if there is data at prediction
              if(!m_initialParticleToggles[partIndex]){
                m_currentParticleTags[partIndex] = 0;
              }else{
                if(!m_separatePM_bToF_allCent){
                  if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) != 0.0){
                    m_currentParticleTags[partIndex] = 1;
                    atLeastOne = true;
                  }else{
                    m_currentParticleTags[partIndex] = 0; // should have All Cent Fit already, just doing amplitudes in this function
                  }
                }else if(m_separatePM_bToF_allCent){
                  if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][pmIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) != 0.0){
                    m_currentParticleTags[partIndex] = 1;
                    atLeastOne = true;
                  }else{
                    m_currentParticleTags[partIndex] = 0; // should have All Cent Fit already, just doing amplitudes in this function
                  }
                }
              }//end else of init tags
            }//end particle loop
            */
            for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
              m_currentParticleTags[partIndex] = 0;
              if(partIndex == 3 || partIndex == 8) continue;
              if((m_currentPartCharge == 1 && m_partBTOFMomStopMatrix_Plus[m_currentPartIndex][partIndex] < m_currentMomOfBin) || (m_currentPartCharge == -1 && m_partBTOFMomStopMatrix_Minus[m_currentPartIndex][partIndex] < m_currentMomOfBin)) continue;
              if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex][m_currentPlusMinusIndex]
               && m_currentPartIndex != partIndex) continue;
                if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
                  || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) ) continue;

                  if(!m_separatePM_bToF_allCent){
                    if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
                      m_currentParticleTags[partIndex] = 1;
                      atLeastOne = true;
                    }
                  }else if(m_separatePM_bToF_allCent){
                    if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][pmIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
                      m_currentParticleTags[partIndex] = 1;
                      atLeastOne = true;
                    }
                  }
            }//end particle loop
          }else if(m_cuts->getEventConfig() == "ColliderCenter"){
              //##########################################################
              //##########    COLLIDER CENTER PARTICLE TAGS   ############
              //##########################################################
            cout << "num particles = " << m_numParticles << endl;
            for(int partIndex = 0; partIndex < m_numParticles; partIndex++){

              m_currentParticleTags[partIndex] = 0;
              if(partIndex == 3 || partIndex == 8) continue;
              //if(m_cuts->getMomStops_bToF(partIndex,m_currentPartCharge) < m_currentMomOfBin) continue;
              if(m_currentPartCharge == 1 && m_partBTOFMomStopMatrix_Plus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;
              if(m_currentPartCharge == -1 && m_partBTOFMomStopMatrix_Minus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;
              if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex][m_currentPlusMinusIndex]
                 && m_currentPartIndex != partIndex) continue;
              if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
                || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) ) continue;

              if(m_allowSpecCentFitsToVary_BTOF && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
                if( (m_currentPartCharge == 1 && m_currentMomOfBin < m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
                  || (m_currentPartCharge == 0 && m_currentMomOfBin < m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
                  || (m_currentPartCharge == -1 && m_currentMomOfBin < m_partBTOF_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex]) ){

                  m_currentParticleTags[partIndex] = 9;
                  atLeastOne = true;
                  continue;
                }
              }




              if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
                  m_currentParticleTags[partIndex] = 1;
                  atLeastOne = true;
              }

            }//end particle loop
          }
          // END PARTILE TAG SETUP



          m_currentNumPartToFit = 0;
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            if(m_currentParticleTags[partIndex] != 0) m_currentNumPartToFit = m_currentNumPartToFit + 1;
          }



          if(!atLeastOne){
            #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin because there are no particles to fit" << endl;
            #endif
            delete m_currentHistoToFit;
            continue;
          }

          if(m_currentParticleTags[m_currentPartIndex] == 0){
            #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin because the main particle will not be fit (tag = 0)" << endl;
            #endif
            delete m_currentHistoToFit;
            continue;
          }



          if(m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex] && m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin) == 0.0){
            getRangeZTOF(lowRange, highRange);

            if(m_currentNumPartToFit == 1){
              double N_3p, mean_3p, sig_3p = 0.0;
              if(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
                mean_3p = m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
                sig_3p = m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin);
                //cout << " window determined by old fit data " << "   sigBefore: " << sig_3p << endl;
                if(m_useStudentTDistributionsForTOF){
                  sig_3p = PhysMath::studentT_var_to_stdDev(
                    m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin),
                    m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin));
                }
                //cout << "      sigAfter: " << sig_3p << endl;
              }else{
                double x1 = - 0.011;
                double x2 = 0.0;
                double x3 = 0.011;
                if(m_massSquaredMode){
                  x2 = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0);
                  x1 = 0.75*x2;
                  if(x1<0.0) x1 = 0.0;
                  x3 = 1.2*x2;
                }
                double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(x1), 2);
                double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(x2), 2);
                double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(x3), 2);
                PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,x1,y1,x2,y2,x3,y3);
              }

              lowRange  = mean_3p  - 3.0*sig_3p;
              highRange  = mean_3p + 3.0*sig_3p;
              if(m_currentPartIndex == 0 && m_currentMomOfBin < 0.58){
                //if pion at low momentum, avoid the muon and electron peaks
                lowRange  = mean_3p  - sig_3p;
                highRange  = mean_3p + 3.0*sig_3p;
              }
              //cout << "   LowRange: " << lowRange << "   high: " << highRange << endl;
            }

            if(m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]){
              m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
              m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,1.0);
              m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
              m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,1.0);
            }
          }else{
            lowRange = m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            highRange = m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }

          lowFunctRange = lowRange;
          highFunctRange = highRange;
          TH1D* origHistToFit = m_currentHistoToFit;
          m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
          delete origHistToFit;
          if(!m_currentHistoToFit) continue;

          double truncatedTotalCounts = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
          double possibleNewBins = truncatedTotalCounts/6.0;
          if(possibleNewBins > 80.0) possibleNewBins = 80.0;
          TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(m_currentHistoToFit,(int) possibleNewBins,true,true);
          delete m_currentHistoToFit;
          m_currentHistoToFit = rebinHisto;
          HistogramUtilities::removeHighCountCorruptData(m_currentHistoToFit);



          m_currentHistoToFit->SetMarkerStyle(20);
          m_currentHistoToFit->SetMarkerColor(kBlack);


          //initialize parameters
          if(m_useStudentTDistributionsForTOF)  initBTOFParams_student();
          else initBTOFParams();



          if(m_currentNumPartToFit ==0){
            #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin because there are no particles to fit" << endl;
            #endif
            delete m_currentHistoToFit;
            continue;
          }





          //################   DEFINE FIT FUNCTION ################
          TF1* fitFunct = NULL;
          if(m_useStudentTDistributionsForTOF)  fitFunct = new TF1("fitFunctBTOF",this,&ZFitter::ZTOFOperator_student,
            lowFunctRange,highFunctRange,m_currentFreeParameters,"ZFitter","ZTOFOperator_student");
            else fitFunct = new TF1("fitFunctBTOF",this,&ZFitter::ZTOFOperator,lowFunctRange,highFunctRange,m_currentFreeParameters,"ZFitter","ZTOFOperator");
          fitFunct->SetNpx(500);



          if(m_useStudentTDistributionsForTOF)  initBTOFFunction_student(fitFunct);
          else initBTOFFunction(fitFunct);



          TF1* fitFunctInitial = (TF1*) fitFunct->Clone();
          fitFunctInitial->SetLineColor(kYellow);
          fitFunctInitial->SetName("fitFunctInitialPlus");


          //#################    FIT THE HISTOGRAM  ###################
          string fitOptionsZBTOFCent = "EQ";
          #ifdef _ZFITTER_DEBUG_
          fitOptionsZBTOFCent = "EV";
          #endif
          if(m_doIntegralOfBinForFits) fitOptionsZBTOFCent += "I";
          m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
          int status = m_currentHistoToFit->Fit("fitFunctBTOF",fitOptionsZBTOFCent.c_str(),"",lowFunctRange,highFunctRange);//add Q to quiet
          TF1* fitResult = m_currentHistoToFit->GetFunction("fitFunctBTOF");

          if(m_currentPartCharge == 1){
            m_FitStatus_Cent_ZbTOF[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status);
          }else{
            m_FitStatus_Cent_ZbTOF[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status);
          }

          #ifdef _ZFITTER_DEBUG_
          cout << "Final Status: " << status << endl;
          #endif

          if(status!=0){
            cout << "\n\n ABNORMAL FIT: status=" << status << "\n\n";
            if(status < 0) continue;
          }

          if(!fitResult || (status != 0 && !m_storeFitsNotStatusZero)){
            #ifdef _ZFITTER_DEBUG_
            cout << " NO FIT RESULT: See Parameters Below:" << endl;
            printCurrentParams();
            #endif

            //################   DRAW THE FIT IMAGE ###################
            if(a_saveImages){
              fittingCanvas->cd();
              m_currentHistoToFit->Draw("E");
              //gPad->SetLogx();
              gPad->SetLogy();
              //fitResult->SetLineColor(kRed);
              if(m_drawInitialSeedsToFits) fitFunctInitial->Draw("same");
              //fitResult->Draw("same");
              vector<TF1*> singleFuncts;
              if(m_useStudentTDistributionsForTOF) singleFuncts = drawIndividualFits_Student();
              else singleFuncts = drawIndividualFits_Gaus();

              DrawPredictionLines_TOF(fittingCanvas);

              gPad->Update();
              gPad->SetRightMargin(0.4);
              gSystem->ProcessEvents();
              fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_________BAD_FIT.png",
                m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), centIndex, HorBinIndex, VertBinIndex));

              gPad->SetLogy(false);
              gPad->Update();
              gSystem->ProcessEvents();
              fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_________BAD_FIT.png",
                m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), centIndex, HorBinIndex, VertBinIndex));

              for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
              for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
              m_TLinePtrs.clear();

            }
            delete m_currentHistoToFit;
            //delete fitResult;
            delete fitFunctInitial;
            delete fitFunct;

            continue;


          }



          vector< double > errorsAmp;
          int paramIndexCounter = 0;
          for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
            if(m_paramToggles[iii] && fitResult){
              m_currentParameters[iii] = fitResult->GetParameter(paramIndexCounter);
              errorsAmp.push_back(fitResult->GetParError(paramIndexCounter));
              paramIndexCounter++;
            }
          }



          if(fitResult){
            double chiSqrNDF = fitResult->GetChisquare() / fitResult->GetNDF();
            if(m_currentPartCharge > 0)  m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][centIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, chiSqrNDF);
            else                         m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][centIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, chiSqrNDF);
          }
          int divNum = -1;
          if(m_useStudentTDistributionsForTOF) divNum = 4;
          else divNum = 3;
          //push to particle fits
          int paramIndex = 0; //the particles not included in fit
          if(m_currentPartCharge == 1){
            for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
              if(m_currentParticleTags[partIndex] !=0 && errorsAmp.size() != 0){
                if(partIndex == m_currentPartIndex){
                  m_Spectra_Cent_ZbTOF[partIndex][centIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                    m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
                  m_Spectra_Cent_ZbTOF[partIndex][centIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                    errorsAmp[paramIndex/divNum] * m_currentIntegralofHistogram);
                }
                m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
                if(errorsAmp[paramIndex/divNum] > 0.0 ) m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                  errorsAmp[paramIndex/divNum] * m_currentIntegralofHistogram);

                  paramIndex+=divNum;
              }
            }
          }else{ //negative
            for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
              if(m_currentParticleTags[partIndex] !=0 && errorsAmp.size() != 0){
                if(partIndex == m_currentPartIndex){
                  m_Spectra_Cent_ZbTOF[partIndex][centIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                    m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
                  m_Spectra_Cent_ZbTOF[partIndex][centIndex][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                    errorsAmp[paramIndex/divNum] * m_currentIntegralofHistogram);
                }
                m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
                if(errorsAmp[paramIndex/divNum] > 0.0 ) m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                  errorsAmp[paramIndex/divNum] * m_currentIntegralofHistogram);

                  paramIndex+=divNum;
              }
            }
          }




          //##############     SAVE THE FIT IMAGE ############################
          if(a_saveImages){
            fittingCanvas->cd();
            m_currentHistoToFit->SetStats(true);
            m_currentHistoToFit->Draw("E");
            m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange, highRange);
            if(m_drawInitialSeedsToFits) fitFunctInitial->Draw("same");
            gPad->SetLogy();
            gPad->SetRightMargin(0.4);
            if(fitResult){
              fitResult->SetLineColor(kRed);
              fitResult->Draw("same");
            }
            vector<TF1*> singleFuncts;
            if(m_useStudentTDistributionsForTOF) singleFuncts = drawIndividualFits_Student();
            else singleFuncts = drawIndividualFits_Gaus();

            DrawPredictionLines_TOF(fittingCanvas);
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_%s.png",
              m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
              centIndex, m_currentRapBin, m_currentMtM0Bin, partChargeNames[pmIndex].c_str()));

            gPad->SetLogy(false);
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/%s/invBetaBTOFFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_%s.png",
              m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
              centIndex, m_currentRapBin, m_currentMtM0Bin, partChargeNames[pmIndex].c_str()));
            for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
            for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
            m_TLinePtrs.clear();
          }

        delete m_currentHistoToFit;
        delete fitFunctInitial;
          //delete fitResult;
        delete fitFunct;

        }// end plus minus histo loop
      }//end pt loop
    }
  }//end centrality loop

  ////draw the results
  //for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
  //  drawPartdEdxFits(partIndex);
  //}
  delete fittingCanvas;

}



void ZFitter::fitBTOFCent_withBackground(int a_centIndex, int a_particleIndex, vector<double> a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT , bool a_useExpBackground, bool a_fitBackgroundOnRight, double a_lowMassOverRide, double a_highMassOverRide){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING BTOF YIELD SPECIFIC CENTRALITIES - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;


  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = a_centIndex;
  m_inAllCentMode = false;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if(!m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]) cout << " Error in ZFitter::fitBTOFCent_withBackground, can't find m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]" << endl;
  /*
  TH2D* base_histogram = dynamic_cast<TH2D*>(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->Clone());
  HistogramUtilities::zeroBins2D(base_histogram);
  base_histogram->SetName("BaseHistogram_ForReference");
  base_histogram->SetTitle("BaseHistogram_ForReference;y;m_{T}-m_{0}");

  cout << "  Starting to form Base Histograms" << endl;
  for(int partIndex=0; partIndex<9; partIndex++){

    if(m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]){
      //Background Amplitude
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetName(Form("%s Space %s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      //
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetName(Form("%s Space %s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b2",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetName(Form("%s Space %s BTOF Background b2",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetName(Form("%s Space %s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetName(Form("%s Space %s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data()));
    }

    if(m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]){
      //Background Amplitude
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetName(Form("%s Space %s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      //
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetName(Form("%s Space %s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_b2",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetName(Form("%s Space %s BTOF Background b2",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetName(Form("%s Space %s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));

      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetName(Form("%sSpace_%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
      m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetName(Form("%s Space %s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data()));
    }
  }//end loop over partIndex

  delete base_histogram;

*/

  m_currentPartCharge = 1;
  m_currentPlusMinusIndex = 0;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsY();
  cout << "numHorBins (rapidity): " << numHorBins << endl;
  cout << "numVertBins (mTm0): " << numVertBins << endl;
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;
      if(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) <= 0.0){
          cout << "Skipping Bin - the centrality integrated fit was not done for this bin" << endl;
          continue;
        }

      cout << "###########################################################################" << endl;
      cout << "Working on Specific Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      for(int partIndex=0; partIndex<9; partIndex++){
        if(partIndex == 3) continue;
        cout << "working in part " << m_currentPartIndex << " binning on partIndex " << partIndex << endl;
        if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;

        if(m_partBTOFMomStopMatrix_Plus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;
        m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Plus[m_currentPartIndex][m_currentCentIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


        double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
        cout << " low: " << lowRange << " high:" << highRange << endl;
        if(m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex] && m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin) == 0.001 ){

          lowRange = m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          highRange = m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);

        }else if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
        else{
          /*
          lowRange = pow(m_partInfo->GetParticleMass(partIndex)-0.1,2.0);
          if(lowRange < HistogramUtilities::getLowEdge(m_currentHistoToFit)) lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
          highRange = pow(m_partInfo->GetParticleMass(partIndex)+0.1,2.0);
          if(highRange > HistogramUtilities::getHighEdge(m_currentHistoToFit)) highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
          */
          getRangeZTOF_mass2(lowRange,highRange,partIndex);
        }
        cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
        if(highRange < lowRange){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " <<  highRange << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }


        if(m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]){
          cout << "m_ZbTOF_LowRange["<< m_currentPartIndex<<"]["<<m_currentPlusMinusIndex<<"]: " << m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex] << endl;

          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        }

        TH1D* origHistToFit = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
        delete origHistToFit;

        double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
        if(countsInHisto < 20){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }

        //TH1D* tempHistoPtr = m_currentHistoToFit;
        //m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
        //delete tempHistoPtr;

        cout << "m_currentHistoToFit: " << m_currentHistoToFit << endl;

        m_currentHistoToFit->Scale(1.0,"width");
        m_currentIntegralofHistogram = countsInHisto;
        m_currentHistoToFit->Scale(1.0/m_currentIntegralofHistogram);




        //set histogram titles
        m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
          partName.c_str(),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
          m_currentMomOfBin
          ));
        m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
        if(m_massSquaredMode){
          m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
        }


        //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
        //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;
        //Extrapolate Means and Widths if going under another peak
        //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
        //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
        //  if(m_partHasBeenExtrapolated[partIndex]){
        //    m_currentParticleTags[partIndex] = 1;
        //    continue;
        //  }
        //}


        double predictLocation = 0.0;
        cout << "m_currentHistoToFit: " << m_currentHistoToFit << endl;
        if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(partIndex),2.0 );
        int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
        double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin)*m_currentIntegralofHistogram;
        cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
        cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
        if(countsInThatBin < 1.0) continue;

        cout << " Starting the fit " << endl;
        vector<double> constrainVect;
        if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) >= 0){
          constrainVect.push_back(-999); // free the amplitude
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][0][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0] << endl;
          if(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) <= 0.0 || m_currentMomOfBin < a_momentum_startUsingBkgd.at(partIndex)){
            //There was no background fit in centrality integrated fits
            constrainVect.push_back(0.0); // fix the scale of the background
            constrainVect.push_back(0.0); // fix the slope of the background
            constrainVect.push_back(0.0); // fix the slope2 of the background
            constrainVect.push_back(0.0); // fix the x_0 of the background
            constrainVect.push_back(0.0); // fix the y_0 of the background
          }else{
            constrainVect.push_back(-999); // free the scale of the background
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd slope
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd slope2
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd x_0
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd y_0
          }
        }else continue;
        if(constrainVect[2] <= 0.0){
          cout << " Error: the constrain vector has an unrealistic constraint for sigma = " << constrainVect[2] << endl;
        }
        vector<double> fitData = fit_BTOF_withBackground(partIndex,constrainVect,a_momentum_startUsingBkgd,a_useStudentT,a_useExpBackground,a_fitBackgroundOnRight);
        cout << " out of fit" << endl;
        if(fitData[0] > 0){
          cout << "fit data size " << fitData.size() << endl;
          cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0] << endl;
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[12]);
          cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0]->GetTitle() << endl;
          cout << "here now doing backgournd data" << endl;

          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][0][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]*m_currentIntegralofHistogram);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]*m_currentIntegralofHistogram);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][1][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][2][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][3][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[16]);
          cout << "m_Background_FitData_partBinning_ZbTOF["<<m_currentPartIndex<<"]["<<partIndex<<"][4][0]: " << m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0] << endl;

          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[8]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[17]);


          //Nov 23
          if(fitData[19] != -999){
            if(m_currentPartCharge > 0)  m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][a_centIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[19]);
          }
            if(partIndex == m_currentPartIndex){
              m_Spectra_Cent_ZbTOF[partIndex][a_centIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  fitData[0]*m_currentIntegralofHistogram);
              m_Spectra_Cent_ZbTOF[partIndex][a_centIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                  fitData[9]*m_currentIntegralofHistogram);
            }
            cout << "m_InvBetaBTOFYieldsCentPlus["<<m_currentPartIndex<<"]["<<partIndex<<"]["<<a_centIndex<<"]: " << m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][a_centIndex]<< endl;
            m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][a_centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]*m_currentIntegralofHistogram);
            if(fitData[9] > 0.0 ) m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][a_centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[9]*m_currentIntegralofHistogram);
        }else{
          cout << "     Error:  fit_BTOF_withBackground : FAILED for Plus" << endl;
        }
      }//end partIndex loop
    }//end mtm0 loop
  } //end rap loop









  // #####################      FIT    MINUS       ######################


  if(!a_fitMinus) return;
  m_currentPartCharge = -1;
  m_currentPlusMinusIndex = 1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
    bool doThisRapBin = false;
    if(m_testing){
      if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
        doThisRapBin = false;
      }else{
        doThisRapBin = true;
      }

      for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
        if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
      }

    }
    if(m_testing && !doThisRapBin) continue;

    for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
      m_currentRapBin = HorBinIndex;
      m_currentMtM0Bin = VertBinIndex;
      m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_invBetaBTOF_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;
      if(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) <= 0.0){
          cout << "Skipping Bin - the centrality integrated fit was not done for this bin" << endl;
          continue;
        }
      cout << "###########################################################################" << endl;
      cout << "Working on Specific Cent ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      for(int partIndex=0; partIndex<9; partIndex++){
        if(partIndex == 3) continue;
        if(!m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) continue;
        if(m_partBTOFMomStopMatrix_Minus[m_currentPartIndex][partIndex] < m_currentMomOfBin) continue;
        m_currentHistoToFit = m_rap_mTm0_invBetaBTOF_Minus[m_currentPartIndex][m_currentCentIndex]->ProjectionZ(
          Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
          m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


        double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
        cout << " low: " << lowRange << " high:" << highRange << endl;
        if(m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex] && m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin) == 0.001 ){

          lowRange = m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          highRange = m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);

        }else if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
        else{
          /*
          lowRange = pow(m_partInfo->GetParticleMass(partIndex)-0.1,2.0);
          if(lowRange < HistogramUtilities::getLowEdge(m_currentHistoToFit)) lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
          highRange = pow(m_partInfo->GetParticleMass(partIndex)+0.1,2.0);
          if(highRange > HistogramUtilities::getHighEdge(m_currentHistoToFit)) highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
          */
          getRangeZTOF_mass2(lowRange,highRange,partIndex);
        }
        cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
        if(highRange < lowRange){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }

        if(m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]){
          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
          m_ZbTOF_LowRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
          m_ZbTOF_HighRange[partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
        }


        TH1D* origHistToFit = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
        delete origHistToFit;

        double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
        if(countsInHisto < 40){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }
        cout << " skipping because counts = " << countsInHisto << endl;
        TH1D* tempHistoPtr = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
        delete tempHistoPtr;
        m_currentHistoToFit->Scale(1.0,"width");
        m_currentIntegralofHistogram = countsInHisto;
        m_currentHistoToFit->Scale(1.0/m_currentIntegralofHistogram);




        // set histogram titles
        m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
          partName.c_str(),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
          m_currentMomOfBin
          ));
        m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
        if(m_massSquaredMode){
          m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
        }

        //if(!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
        //if(m_currentMomOfBin < m_bTOFSeparatedStops[m_currentPartIndex] && m_currentPartIndex != partIndex) continue;

        //Extrapolate Means and Widths if going under another peak
        //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
        //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
        //  if(m_partHasBeenExtrapolated[partIndex]){
        //    m_currentParticleTags[partIndex] = 1;
        //    continue;
        //  }
        //}


        double predictLocation = 0.0;
        if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(partIndex),2.0 );
        int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
        double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin)*m_currentIntegralofHistogram;
        cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
        cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
        if(countsInThatBin < 1.0) continue;

        cout << " Starting the fit " << endl;
        vector<double> constrainVect;
        if(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) >= 0){
          constrainVect.push_back(-999); // free the amplitude
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          if(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) <= 0.0){
            //There was no background fit in centrality integrated fits
            constrainVect.push_back(0.0); // fix the scale of the background
            constrainVect.push_back(0.0); // fix the slope of the background
            constrainVect.push_back(0.0); // fix the slope2 of the background
            constrainVect.push_back(0.0); // fix the x_0 of the background
            constrainVect.push_back(0.0); // fix the y_0 of the background
          }else{
            constrainVect.push_back(-999); // free the scale of the background
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd slope
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd slope2
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd x_0
            constrainVect.push_back(m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd y_0
          }
        }
        if(constrainVect[2] <= 0.0){
          cout << " Error: the constrain vector has an unrealistic constraint for sigma = " << constrainVect[2] << endl;
        }
        vector<double> fitData = fit_BTOF_withBackground(partIndex,constrainVect,a_momentum_startUsingBkgd,a_useStudentT,a_useExpBackground,a_fitBackgroundOnRight);
        cout << " out of fit" << endl;
        if(fitData[0] > 0){
          cout << "fit data size " << fitData.size() << endl;
          cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1] << endl;
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]*m_currentIntegralofHistogram);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[12]);
          cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1]->GetTitle() << endl;
          cout << "here now doing backgournd data" << endl;
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]*m_currentIntegralofHistogram);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][0][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]*m_currentIntegralofHistogram);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][1][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][2][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][3][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[16]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[8]);
          m_Background_FitData_partBinning_ZbTOF[m_currentPartIndex][partIndex][4][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[17]);

          if(fitData[19] != -999){
            if(m_currentPartCharge < 0)  m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][a_centIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[19]);
          }
              if(partIndex == m_currentPartIndex){
                m_Spectra_Cent_ZbTOF[partIndex][a_centIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]*m_currentIntegralofHistogram);
                m_Spectra_Cent_ZbTOF[partIndex][a_centIndex][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[9]*m_currentIntegralofHistogram);
              }
              m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][a_centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]*m_currentIntegralofHistogram);
              if(fitData[9] > 0.0 ) m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][a_centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[9]*m_currentIntegralofHistogram);
        }else{
          cout << "     Error:  fitBTOFCent_withBackground : FAILED for Minus" << endl;
        }
        constrainVect.clear();
        #ifdef _CPP11_
          constrainVect.shrink_to_fit();
        #endif
      }//end partIndex loop
    }//end mtm0 loop
  } //end rap loop

}



void ZFitter::fitBTOFCent_NonMixing(int a_centIndex, int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus, bool a_useStudentT , double a_lowMassOverRide, double a_highMassOverRide){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING BTOF YIELD SPECIFIC CENTRALITIES - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;

  if(a_particleIndex < 0 || a_particleIndex == 3  || a_particleIndex > 7 ){
    cout << " Error in ZFitter::fitBTOFCent_NonMixing  with invalid argument, a_particleIndex" << endl;
    return;
  }
  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = a_centIndex;
  m_inAllCentMode = false;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if(!m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]) cout << " Error in ZFitter::fitBTOF_NonMixing, can't find m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  /*TH2D* base_histogram = dynamic_cast<TH2D*>(m_InvBetaBTOFFitData[m_currentPartIndex][m_currentPartIndex][0]->Clone());
  HistogramUtilities::zeroBins2D(base_histogram);
  base_histogram->SetName("BaseHistogram_ForReference");
  base_histogram->SetTitle("BaseHistogram_ForReference;y;m_{T}-m_{0}");

  cout << "  Starting to form Base Histograms" << endl;
  //Background Amplitude
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetName(Form("%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->SetName(Form("%s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetName(Form("%s_BTOF_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][0][1]->SetName(Form("%s BTOF Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  //
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetName(Form("%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][0]->SetName(Form("%s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetName(Form("%s_BTOF_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][1][1]->SetName(Form("%s BTOF Background b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetName(Form("%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][0]->SetName(Form("%s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetName(Form("%s_BTOF_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][2][1]->SetName(Form("%s BTOF Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetName(Form("%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][0]->SetName(Form("%s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetName(Form("%s_BTOF_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZbTOF[m_currentPartIndex][3][1]->SetName(Form("%s BTOF Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  */

  for(int pmIndex= 0; pmIndex< 2; pmIndex++){
    if(!a_fitMinus && pmIndex == 1) continue;
    TH3D* histo3DPM = NULL;
    string partName;
    string partChargeNames[2] = {"Positive","Negative"};
    string partChargeNamesSmall[2] = {"Plus","Minus"};
    if(pmIndex== 0){
      histo3DPM = m_rap_mTm0_invBetaBTOF_Plus[m_currentPartIndex][m_currentCentIndex];
      partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
      m_currentPartCharge = 1;
    }else{
      histo3DPM = m_rap_mTm0_invBetaBTOF_Minus[m_currentPartIndex][m_currentCentIndex];
      partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
      m_currentPartCharge = -1;
    }

    if(!histo3DPM){
      cout << "No m_rap_mTm0_invBetaBTOF_Minus[m_currentPartIndex][m_currentCentIndex] for particle " << m_currentPartIndex << "   ptr:" << histo3DPM<< endl;
      break;
    }

    int numHorBins = histo3DPM->GetNbinsX();
    int numVertBins = histo3DPM->GetNbinsY();
    for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
      cout << " Starting HorBin " << HorBinIndex << endl;
      m_currentRapidity = histo3DPM->GetXaxis()->GetBinCenter(HorBinIndex);
      //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
      int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
      bool doThisRapBin = false;
      if(m_testing){
        if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
          doThisRapBin = false;
        }else{
          doThisRapBin = true;
        }

        for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
          if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
        }

      }
      if(m_testing && !doThisRapBin) continue;

      for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
        m_currentRapBin = HorBinIndex;
        m_currentMtM0Bin = VertBinIndex;
        m_currentMtM0 = histo3DPM->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
        //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
        if(m_cuts->getMomStops_bToF(m_currentPartIndex,m_currentPartCharge) < m_currentMomOfBin) continue;
        if(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][pmIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) <= 0.0){
          cout << "Skipping Bin - the centrality integrated fit was not done for this bin" << endl;
          continue;
        }




        cout << "###########################################################################" << endl;
        cout << "Working on Cent " << m_currentCentIndex << " ZBTOF " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;



        //#############    SETTING UP HISTOGRAM TO FIT #################################
        m_currentHistoToFit = histo3DPM->ProjectionZ(
          Form("%s_BTOF_%02d_%02d_%02d",partName.c_str(),m_currentCentIndex, m_currentRapBin, m_currentMtM0Bin),
          m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);

        double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);

        cout << " low: " << lowRange << " high:" << highRange << endl;
        if(!m_massSquaredMode) getRangeZTOF_Simple(lowRange, highRange, a_lowMassOverRide, a_highMassOverRide);
        else{
          lowRange = pow(a_lowMassOverRide,2.0);
          highRange = pow(a_highMassOverRide,2.0);
        }
        cout << "   by ztofRange low: " << lowRange << "  high:" << highRange << endl;
        if(highRange < lowRange){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping because Ztof Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }
        TH1D* origHistToFit = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
        delete origHistToFit;

        double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
        if(countsInHisto < 20){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts : " << countsInHisto << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }
        cout << " skipping because counts = " << countsInHisto << endl;
        TH1D* tempHistoPtr = m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
        delete tempHistoPtr;
        m_currentHistoToFit->Scale(1.0,"width");

        //set histogram titles
        m_currentHistoToFit->SetTitle(Form("BTOF %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p_{tot}=%2.3f",
          partName.c_str(),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
          m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
          m_currentMomOfBin
          ));
        if(m_massSquaredMode){
          m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/dm^{2}");
        }else{
          m_currentHistoToFit->GetXaxis()->SetTitle("#beta^{-1}");
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/d#beta^{-1}");
        }




        //Extrapolate Means and Widths if going under another peak
        //if(m_currentMomOfBin > m_partBTOF_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
        //  if(!m_partHasBeenExtrapolated[partIndex]) extrapolateMeanAndWidthParams_BTOF(partIndex);
        //  if(m_partHasBeenExtrapolated[partIndex]){
        //    m_currentParticleTags[partIndex] = 1;
        //    continue;
        //  }
        //}




        double predictLocation = 0.0;
        if(m_massSquaredMode) predictLocation = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0 );
        int invBetaBin = m_currentHistoToFit->FindBin(predictLocation);
        double countsInThatBin = m_currentHistoToFit->GetBinContent(invBetaBin)*m_currentHistoToFit->GetBinWidth(invBetaBin);
        cout << " InvBeta Bin = " << invBetaBin << "   counts in that bin: " << countsInThatBin << endl;
        cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
        if(countsInThatBin < 1.0) continue;

        cout << " Starting the fit " << endl;
        vector<double> constrainVect;
        if(a_constrainMinusToPlus){
          constrainVect.push_back(-999); // free the amplitude
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][pmIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][pmIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          constrainVect.push_back(m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][pmIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
          if(m_Background_FitData_ZbTOF[m_currentPartIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) <= 0.0){
            //There was no background fit in centrality integrated fits
            constrainVect.push_back(0.0); // fix the scale of the background
            constrainVect.push_back(1.0); // fix the slope of the background
            constrainVect.push_back(0.0); // fix the x_0 of the background
            constrainVect.push_back(0.0); // fix the y_0 of the background
          }else{
            constrainVect.push_back(-999); // free the scale of the background
            constrainVect.push_back(m_Background_FitData_ZbTOF[m_currentPartIndex][1][pmIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd slope
            constrainVect.push_back(m_Background_FitData_ZbTOF[m_currentPartIndex][2][pmIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd x_0
            constrainVect.push_back(m_Background_FitData_ZbTOF[m_currentPartIndex][3][pmIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); //bkgd y_0
          }
        }
        if(constrainVect[2] <= 0.0){
          cout << " Error: the constrain vector has an unrealistic constraint for sigma = " << constrainVect[2] << endl;
        }
        vector<double> fitData = fitIsolation_BTOF_withBackground(constrainVect,a_useStudentT);
        cout << " out of fit" << endl;
        if(fitData[0] > 0){
          cout << "fit data size " << fitData.size() << endl;
          cout << "   ptr data: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0] << endl;
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][pmIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][pmIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][pmIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][1][pmIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][pmIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][2][pmIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][pmIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
          m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][3][pmIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);

          cout << " filling histo named: " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->GetName() << "       " << m_Fit_Data_ZbTOF[m_currentPartIndex][m_currentPartIndex][0][0]->GetTitle() << endl;
          cout << "here now doing backgournd data" << endl;
          m_Background_FitData_ZbTOF[m_currentPartIndex][0][pmIndex]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
          m_Background_FitData_ZbTOF[m_currentPartIndex][0][pmIndex]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
          m_Background_FitData_ZbTOF[m_currentPartIndex][1][pmIndex]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
          m_Background_FitData_ZbTOF[m_currentPartIndex][1][pmIndex]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
          m_Background_FitData_ZbTOF[m_currentPartIndex][2][pmIndex]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
          m_Background_FitData_ZbTOF[m_currentPartIndex][2][pmIndex]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
          m_Background_FitData_ZbTOF[m_currentPartIndex][3][pmIndex]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
          m_Background_FitData_ZbTOF[m_currentPartIndex][3][pmIndex]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);


          m_Spectra_Cent_ZbTOF[m_currentPartIndex][m_currentCentIndex][pmIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,  fitData[0]);
          m_Spectra_Cent_ZbTOF[m_currentPartIndex][m_currentCentIndex][pmIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,    fitData[8]);
          m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][m_currentPartIndex][m_currentCentIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,  fitData[0]);
          m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][m_currentPartIndex][m_currentCentIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,    fitData[8]);


        }else{
          cout << "     Error:  fitIsolation_BTOF_withBackground : FAILED for " << partName << endl;
        }

      }//end mtm0 loop
      for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
      m_TLinePtrs.clear();
    } //end rap loop
  }





  return;
}










//========================================================================================================================================================================================================================
//========================================================================================================================================================================================================================






/*
//fit as function of centrality using results from fitInvBeta
void ZFitter::fitMassSquaredCent(int a_centIndex, bool a_saveImages){

  cout << "Starting ZFitter::fitMassSquaredCent" << endl;

  m_inAllCentMode = false;
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);


  int numHorBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetNbinsY();
  //int numCentralities = m_numParticles; // cout << "Warning hardcoded numCentralities = m_numParticles" << endl;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    if(a_centIndex !=-1 && centIndex != a_centIndex) continue;
    m_currentCentIndex = centIndex;
    cout << "m_currentCentIndex: " <<  m_currentCentIndex << endl;
    for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
      m_currentRapidity = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
      cout << "m_currentRapidity: " << m_currentRapidity << endl;
      // make testing do the same bin as dEdx so
      int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
      int dEdxBin = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
      //cout << "Rapidity: " << m_currentRapidity << endl;
      //if(m_testing && (dEdxBin%(numDeDxBins/m_divTesting) != 0 ? (float)m_currentRapidity != m_cuts->getYCM() : 0 )) continue;

      bool doThisRapBin = false;
      if(m_testing){
        if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(m_currentRapidity - m_cuts->getYCM()) < 0.01 )  ){
          doThisRapBin = false;
        }else{
          doThisRapBin = true;
        }

        for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
          if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
        }

      }
      if(m_testing && !doThisRapBin) continue;

      for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
        m_currentRapBin = HorBinIndex;
        m_currentMtM0Bin = VertBinIndex;
        m_currentMtM0 = m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));

        //#### PLUS MINUS HISTO LOOP #########
        for(int pmIndex= 0; pmIndex< 2; pmIndex++){
          TH3D* histo3DPM = NULL;
          string partChargeNames[2] = {"Positive","Negative"};
          string partChargeNamesSmall[2] = {"Plus","Minus"};
          if(pmIndex== 0){
            histo3DPM = m_rap_mTm0_invBetaBTOF_Plus[m_currentPartIndex][centIndex];
            m_currentPartCharge = 1;
          }else{
            histo3DPM = m_rap_mTm0_invBetaBTOF_Minus[m_currentPartIndex][centIndex];
            m_currentPartCharge = -1;
          }
          m_currentPlusMinusIndex = pmIndex;
          cout << "\n\n\n###########################################################################" << endl;
          cout << "Working on Specific Centrality ZBTOF " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " Cent: " << centIndex << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << endl;


          //################   SETUP HISTOGRAM TO FIT #####################
          m_currentHistoToFit = histo3DPM->ProjectionZ(
          Form("BTOF_%s_%02d_%02d",m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin,m_currentMtM0Bin), m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
          HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
          double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
          if((countsInHisto < 100 && !m_testing) || (countsInHisto < 20)) {
            cout << "not enough counts in hist" << endl;
            delete m_currentHistoToFit;
            #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin due to low counts: " << countsInHisto << endl;
            #endif
            continue;
          }
          m_currentIntegralofHistogram = countsInHisto;

          //TH1D* tempHistoPtr = m_currentHistoToFit;
          //m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,5.0);
          //delete tempHistoPtr;
          HistogramUtilities::makeDifferential(m_currentHistoToFit, 1.0/m_currentIntegralofHistogram);
          //m_currentHistoToFit->Scale(1.0,"width");

          double lowRange = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistoToFit));
          double highRange = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistoToFit)+1 );
          double lowFunctRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
          double highFunctRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);


          m_currentHistoToFit->SetTitle(Form("BTOF %s Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] (GeV/c) p_{tot}=%2.3f",
            partChargeNames[pmIndex].c_str(),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
            m_rap_mTm0_invBetaBTOF_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
            m_currentMomOfBin
            ));
          m_currentHistoToFit->GetXaxis()->SetTitle("m^{2}");
          m_currentHistoToFit->GetYaxis()->SetTitle(Form("#frac{1}{%E} dN/d(m^{2})",m_currentIntegralofHistogram));


          double massSqrPrediction = pow(m_partInfo->GetParticleMass(m_currentPartIndex),2.0);
          if(m_currentPartIndex == 6 || m_currentPartIndex == 7) massSqrPrediction = massSqrPrediction/2.0;
          int lowBinReturn = -1;
          int highBinReturn = -1;
          HistogramUtilities::getPeakRange(m_currentHistoToFit, m_currentHistoToFit->FindBin(massSqrPrediction),
            lowBinReturn, highBinReturn,0.05);
          #ifdef _ZFITTER_DEBUG_
            cout << " HistogramUtilities::getPeakRange gives a lowbin " << lowBinReturn << "  and highBin " << highBinReturn << endl;
          #endif
          lowRange = m_currentHistoToFit->GetBinLowEdge(lowBinReturn);
          highRange = m_currentHistoToFit->GetBinLowEdge(highBinReturn+1);

          //find the number of particles within this range
          vector< int > liveParticleIndecies;
          vector< double > liveParticleMassSqr;
          for(int partIndex = 0; partIndex < 8; partIndex++){
            if(partIndex == 3) continue;
            double massSqr = pow(m_partInfo->GetParticleMass(partIndex),2.0);
            if(partIndex == 6 || partIndex == 7) massSqr = massSqr/2.0;
            if(massSqr > lowRange && massSqr < highRange){
              liveParticleIndecies.push_back(partIndex);
              liveParticleMassSqr.push_back(massSqr);
            }
          }

          m_currentNumPartToFit = liveParticleIndecies.size();

          if(m_currentNumPartToFit == 0){
            #ifdef _ZFITTER_DEBUG_
            cout << "Skipping this bin because there are no particles to fit" << endl;
            #endif
            delete m_currentHistoToFit;
            continue;
          }



          if(m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin) == 0.0){
            m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
            m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,1.0);
            m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
            m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,1.0);
          }else{
            lowRange = m_ZbTOF_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            highRange = m_ZbTOF_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }

          lowFunctRange = lowRange;
          highFunctRange = highRange;
          m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange,highRange);
          m_currentHistoToFit->SetMarkerStyle(20);
          m_currentHistoToFit->SetMarkerColor(kBlack);

          double x1_loc = HistogramUtilities::getLowEdge(m_currentHistoToFit) + 2;
          double x2_loc = HistogramUtilities::getLowEdge(m_currentHistoToFit) + 0.1*(HistogramUtilities::getHighEdge(m_currentHistoToFit) - HistogramUtilities::getLowEdge(m_currentHistoToFit)) ;
          double bin1 = m_currentHistoToFit->FindBin(x1_loc);
          double bin2 = m_currentHistoToFit->FindBin(x2_loc);
          if(bin1 == bin2) cout << " bad init because bin1 = bin2" << endl;
          double x1 = m_currentHistoToFit->GetBinCenter(bin1);
          double x2 = m_currentHistoToFit->GetBinCenter(bin2);
          double y1 = m_currentHistoToFit->GetBinContent(bin1);
          double y2 = m_currentHistoToFit->GetBinContent(bin2);
          double expConst = log(y1 - y2)/(x2-x1);
          double expScale = y1*exp(expConst*x1);


          TF1* fitFunction = new TF1("fitFunction_m2","[0]*TMath::Student((x - [1])/[2], [3])/[2] + [3]*TMath::Exp(-[4]*x)",0,20);
          fitFunction->SetParameter(0,1.0);
          fitFunction->SetParameter(1,pow(m_partInfo->GetParticleMass(liveParticleIndecies[0]),2.0));
          fitFunction->SetParameter(2,studentT_var_from_stdDev(0.5, 50)); // stdDev nu
          fitFunction->SetParameter(3,50);
          fitFunction->SetParameter(4,expScale);
          fitFunction->SetParameter(5,expConst);
          fitFunction->SetParameter(0,m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(liveParticleMassSqr[0]))/fitFunction->Eval(liveParticleMassSqr[0]));
          fitFunction->SetParLimits(0,0.1*fitFunction->GetParameter(0),10.0*fitFunction->GetParameter(0));
          fitFunction->SetParLimits(1,fitFunction->GetParameter(1) - 0.05,fitFunction->GetParameter(1) + 0.05);
          fitFunction->SetParLimits(2,0.00,15.0);
          fitFunction->SetParLimits(3,1.0,200);
          fitFunction->SetParLimits(4,0.2*expScale,5.0*expScale);
          fitFunction->SetParLimits(5,0.2*expConst,5.0*expConst);

          TF1* background_funct = new TF1("[0]*TMath::Exp(-[1]*x)",0,20);
          background_funct->SetLineColor(kRed);
          TF1* source_funct = new TF1("fitFunction_m2","[0]*TMath::Student((x - [1])/[2], [3])/[2]",0,20);
          source_funct->SetLineColor(kBlue);


          / *
          }else if(m_currentNumPartToFit == 2){
            fitFunction = new TF1("fitFunction_m2",
              "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2] + [3]*0.3989422804 * TMath::Exp(-0.5*( x - [4] )*(x - [4])/([5]*[5])) /[5]",0,20);
            fitFunction->SetParameter(0,1.0);
            fitFunction->SetParameter(1,pow(m_partInfo->GetParticleMass(liveParticleIndecies[0]),2.0));
            fitFunction->SetParameter(2,0.7);
            fitFunction->SetParameter(3,0.0);
            fitFunction->SetParameter(4,pow(m_partInfo->GetParticleMass(liveParticleIndecies[1]),2.0));
            fitFunction->SetParameter(5,0.7);
            double amp_1 = fitFunction->Eval(liveParticleMassSqr[0]);
            amp_1 = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(liveParticleMassSqr[0]))/amp_1;
            fitFunction->SetParameter(0,0.0);
            fitFunction->SetParameter(3,1.0);
            double amp_2 = fitFunction->Eval(liveParticleMassSqr[1]);
            fitFunction->SetParameter(3,m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(liveParticleMassSqr[1]))/amp_2);
            fitFunction->SetParameter(0,amp_1);

            fitFunction->SetParLimits(0,0.1*fitFunction->GetParameter(0),10.0*fitFunction->GetParameter(0));
            fitFunction->SetParLimits(1,fitFunction->GetParameter(1) - 0.05,fitFunction->GetParameter(1) + 0.05);
            fitFunction->SetParLimits(2,0.01,10.0);
            fitFunction->SetParLimits(3,0.1*fitFunction->GetParameter(3),10.0*fitFunction->GetParameter(3));
            fitFunction->SetParLimits(4,fitFunction->GetParameter(4) - 0.05,fitFunction->GetParameter(4) + 0.05);
            fitFunction->SetParLimits(5,0.01,10.0);

          }else if(m_currentNumPartToFit == 3){
            fitFunction = new TF1("fitFunction_m2",
              "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2] + [3]*0.3989422804 * TMath::Exp(-0.5*( x - [4] )*(x - [4])/([5]*[5])) /[5] + [6]*0.3989422804 * TMath::Exp(-0.5*( x - [7] )*(x - [7])/([8]*[8])) /[8]",0,20);
            fitFunction->SetParameter(0,1.0);
            fitFunction->SetParameter(1,pow(m_partInfo->GetParticleMass(liveParticleIndecies[0]),2.0));
            fitFunction->SetParameter(2,0.7);
            fitFunction->SetParameter(3,0.0);
            fitFunction->SetParameter(4,pow(m_partInfo->GetParticleMass(liveParticleIndecies[1]),2.0));
            fitFunction->SetParameter(5,0.7);
            fitFunction->SetParameter(6,0.0);
            fitFunction->SetParameter(7,pow(m_partInfo->GetParticleMass(liveParticleIndecies[2]),2.0));
            fitFunction->SetParameter(8,0.7);
            double amp_1 = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(liveParticleMassSqr[0]))/fitFunction->Eval(liveParticleMassSqr[0]);
            fitFunction->SetParameter(0,0.0);
            fitFunction->SetParameter(3,1.0);
            double amp_2 = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(liveParticleMassSqr[1]))/fitFunction->Eval(liveParticleMassSqr[1]);
            fitFunction->SetParameter(3,1.0);
            fitFunction->SetParameter(6,0.0);
            double amp_3 = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(liveParticleMassSqr[2]))/fitFunction->Eval(liveParticleMassSqr[2]);
            fitFunction->SetParameter(0,amp_1);
            fitFunction->SetParameter(3,amp_2);
            fitFunction->SetParameter(6,amp_3);

            fitFunction->SetParLimits(0,0.1*fitFunction->GetParameter(0),10.0*fitFunction->GetParameter(0));
            fitFunction->SetParLimits(1,fitFunction->GetParameter(1) - 0.05,fitFunction->GetParameter(1) + 0.05);
            fitFunction->SetParLimits(2,0.01,10.0);
            fitFunction->SetParLimits(3,0.1*fitFunction->GetParameter(3),10.0*fitFunction->GetParameter(3));
            fitFunction->SetParLimits(4,fitFunction->GetParameter(4) - 0.05,fitFunction->GetParameter(4) + 0.05);
            fitFunction->SetParLimits(5,0.01,10.0);
            fitFunction->SetParLimits(6,0.1*fitFunction->GetParameter(6),10.0*fitFunction->GetParameter(6));
            fitFunction->SetParLimits(7,fitFunction->GetParameter(7) - 0.05,fitFunction->GetParameter(7) + 0.05);
            fitFunction->SetParLimits(8,0.01,10.0);

          }* /

          fitFunction->SetNpx(300);




          TF1* fitFunctInitial = (TF1*) fitFunction->Clone();
          fitFunctInitial->SetLineColor(kYellow);
          fitFunctInitial->SetName("fitFunctInitialPlus");


          //#################    FIT THE HISTOGRAM  ###################
          string fitOptionsZBTOFCent = "EQ";
          #ifdef _ZFITTER_DEBUG_
          fitOptionsZBTOFCent = "EV";
          #endif
          int status = m_currentHistoToFit->Fit(fitFunction,fitOptionsZBTOFCent.c_str(),"",lowFunctRange,highFunctRange);//add Q to quiet
          TF1* fitResult = m_currentHistoToFit->GetFunction("fitFunction_m2");

          if(m_currentPartCharge > 0){
            m_FitStatus_Cent_ZbTOF[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status);
          }else{
            m_FitStatus_Cent_ZbTOF[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status);
          }

          #ifdef _ZFITTER_DEBUG_
            cout << "Final Status: " << status << endl;
          #endif

          if(status!=0){
            cout << "\n\n ABNORMAL FIT: status=" << status << "\n\n";
          }
          #ifdef _ZFITTER_DEBUG_
          if(!fitResult){
            cout << " NO FIT RESULT: See Parameters Below:" << endl;
            printCurrentParams();
          }
          #endif


          if(fitResult){
            double chiSqrNDF = fitResult->GetChisquare() / fitResult->GetNDF();
            if(m_currentPartCharge > 0)  m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][centIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, chiSqrNDF);
            else                         m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][centIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, chiSqrNDF);
          }

          //push to particle fits
          int paramIndex = 0; //the particles not included in fit
          if(m_currentPartCharge == 1){
            for(int iii = 0; iii < liveParticleIndecies.size(); iii++){
              int partIndex = liveParticleIndecies[iii];
              if(partIndex == m_currentPartIndex){
                m_Spectra_Cent_ZbTOF[partIndex][centIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  fitResult->GetParameter(paramIndex) * m_currentIntegralofHistogram );
                m_Spectra_Cent_ZbTOF[partIndex][centIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                  fitResult->GetParError(paramIndex) * m_currentIntegralofHistogram);
              }
              m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                fitResult->GetParameter(paramIndex) * m_currentIntegralofHistogram);
              if(fitResult->GetParError(paramIndex) > 0.0 ) m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                fitResult->GetParError(paramIndex) * m_currentIntegralofHistogram);
                paramIndex += 3;
            }

          }else{ //negative
            for(int iii = 0; iii < liveParticleIndecies.size(); iii++){
              int partIndex = liveParticleIndecies[iii];
              if(partIndex == m_currentPartIndex){
                m_Spectra_Cent_ZbTOF[partIndex][centIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  fitResult->GetParameter(paramIndex) * m_currentIntegralofHistogram );
                m_Spectra_Cent_ZbTOF[partIndex][centIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                  fitResult->GetParError(paramIndex) * m_currentIntegralofHistogram);
              }
              m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                fitResult->GetParameter(paramIndex) * m_currentIntegralofHistogram);
              if(fitResult->GetParError(paramIndex) > 0.0 ) m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                fitResult->GetParError(paramIndex) * m_currentIntegralofHistogram);
                paramIndex += 3;
            }
          }




          //##############     SAVE THE FIT IMAGE ############################
          if(a_saveImages){
            fittingCanvas->cd();
            m_currentHistoToFit->SetStats(true);
            m_currentHistoToFit->Draw("E");
            m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange, highRange);
            if(m_drawInitialSeedsToFits) fitFunctInitial->Draw("same");
            gPad->SetLogy();
            gPad->SetRightMargin(0.4);
            if(fitResult){
              fitResult->SetLineColor(kRed);
              fitResult->Draw("same");
            }
            vector<TF1*> singleFuncts;
            if(m_useStudentTDistributionsForTOF) singleFuncts = drawIndividualFits_Student();
            else singleFuncts = drawIndividualFits_Gaus();

            DrawPredictionLines_TOF(fittingCanvas);
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/%s/massSqrBTOFFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_%s.png",
              m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
              centIndex, m_currentRapBin, m_currentMtM0Bin, partChargeNames[pmIndex].c_str()));

            gPad->SetLogy(false);
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/%s/massSqrBTOFFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_%s.png",
              m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
              centIndex, m_currentRapBin, m_currentMtM0Bin, partChargeNames[pmIndex].c_str()));
            for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
            for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
            m_TLinePtrs.clear();
          }

        delete m_currentHistoToFit;
        delete fitFunctInitial;
          //delete fitResult;
        delete fitFunction;

        }// end plus minus histo loop
      }//end pt loop
    }
  }//end centrality loop

  ////draw the results
  //for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
  //  drawPartdEdxFits(partIndex);
  //}
  delete fittingCanvas;

}

*/
