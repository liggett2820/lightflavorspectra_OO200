// THIS IS INSERTED INTO ZFITTER.cxx // Just separate file becuase it's long
//============================================================================================================
//#########      TPC FITTING (dE/dx)   #####################
//============================================================================================================

void ZFitter::initZTPCParams(){
  #ifdef _ZFITTER_DEBUG_
    cout << "    #####    Initiallizing TPC Parameters:   ####" << endl;
  #endif
  vector<double> params;
  vector<string> paramNames;
  vector<bool> paramToggles;
  m_currentFreeParameters = 0;
  //m_currentNumPartToFit = 0;
  double sigFractionGuess = 0.05;
  double tpcResolution = 0.07; // 7 percent

  //double ndfParamGuess = 15.0;
  //double sigParGuess = sqrt(sigmaGuess*sigmaGuess*(ndfParamGuess - 2.0)/ndfParamGuess); // convert gaussian std dev to sigma of student T

  #ifdef _ZFITTER_DEBUG_
    cout << "  Number of particles to fit: " << m_currentNumPartToFit << endl;
  #endif

  bool matrixMethodFailed = false;
  TVectorD ampVect(m_currentNumPartToFit);


  if(m_useMatrixInitializationForTPC){
    //matrix inversion for amplitudes
    cout << "Matrix Method Starting: " << endl;
    vector< double > matrixElements;
    vector< double > countsVector;
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      //this particle index loop uses the mean of this partIndex and finds the probabiliies for all the other particle probabilities for that point
      if(m_currentParticleTags[partIndex] == 0) continue;
      double xPoint = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      cout << "     First Point: " << xPoint;
      if(std::isnan(xPoint)) matrixMethodFailed = true;

      int predictBin = m_currentHistoToFit->FindBin(xPoint);
      double binContent = m_currentHistoToFit->GetBinContent(predictBin);
      cout << " BinCounts: " << binContent << endl;
      countsVector.push_back(binContent);
      for(int partIndex2 = 0; partIndex2 < m_numParticles; partIndex2++){

        if(m_currentParticleTags[partIndex2] == 1 || m_currentParticleTags[partIndex2] == 3){
          double mu = -1;
          double sig = -1;
          if(m_currentPartCharge == 1){
            mu = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex2][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            sig = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex2][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }else{
            mu = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex2][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            sig = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex2][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }
          double prob = gaus(mu,sig,xPoint);
          cout << "        Mu:" << mu << " sig:" << sig << " prob:" << prob << endl;
          //if(prob < 1E-8) matrixElements.push_back( 0.0 );
          //else
          matrixElements.push_back( prob );
        }else if(m_currentParticleTags[partIndex2] == 2){
          //double mean = m_bichselCurves[partIndex2]->Eval(m_currentMomOfBin);
          //if(mean == 0.0) matrixMethodFailed = true;
          //double stdDev = 0.07*mean;
          //double mu,sig;
          //PhysMath::gausToLogNorm(mean,stdDev,mu,sig);

          //double mu = 0;
          //double sig = 0;
          //if(m_currentRapBin - 1 >= 0 && m_currentRapBin -1 < m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[1].size()){
          //mu = m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[1][m_currentRapBin-1]->Eval(PhysMath::betagamma(m_currentMomOfBin, m_partInfo->GetParticleMass(partIndex)));
          //sig =  m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[2][m_currentRapBin-1]->Eval(PhysMath::betagamma(m_currentMomOfBin, m_partInfo->GetParticleMass(partIndex)));
          //}
          //cout << "mu: " << mu << " sig: " << sig << endl;
          //if(sig == 0.0 || std::isnan(mu) || std::isnan(sig)){


          // #####  THIS SECTION MUST BE IDENTICAL TO THE ONE 40 LINES BELOW #####
          double mean = log(m_bichselCurves[partIndex2]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double stdDev = log(1.0 + tpcResolution);
          //double stdDev = sigFractionGuess;// change 40 lines below too
          if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
            mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }else{
            if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
              mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }

            if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
              mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }

          }

          //############################# END OF SECTION #######################


          double prob = gaus(mean,stdDev,xPoint);
          cout << "        Mu:" << mean << " sig:" << stdDev << " prob:" << prob << endl;
          //if(prob < 1E-8) matrixElements.push_back( 0.0 );
          //else
          matrixElements.push_back( prob );
        }else if(m_currentParticleTags[partIndex2] == 4){
          //initialization for matrix method is not going to work for this perfectly, but that's okay
          //these are not seeds to the fit
          double mean = log(m_bichselCurves[partIndex2]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double stdDev = log(1.0 + tpcResolution);
          //double stdDev = sigFractionGuess;
          double prob = gaus(mean,stdDev,xPoint);
          cout << "        Mu:" << mean << " sig:" << stdDev << " prob:" << prob << endl;
          //if(prob < 1E-8) matrixElements.push_back( 0.0 );
          //else
          matrixElements.push_back( prob );

        }else if(m_currentParticleTags[partIndex2] == 5){
          //initialization for matrix method is not going to work for this perfectly, but that's okay
          //these are not seeds to the fit
          double mean = log(m_bichselCurves[partIndex2]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double stdDev = log(1.0 + tpcResolution);
          //double stdDev = sigFractionGuess;
          double prob = gaus(mean,stdDev,xPoint);
          cout << "        Mu:" << mean << " sig:" << stdDev << " prob:" << prob << endl;
        }else if(m_currentParticleTags[partIndex2] == 6){
          //initialization for matrix method is not going to work for this perfectly, but that's okay
          //these are not seeds to the fit
          double mean = log(m_bichselCurves[partIndex2]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double stdDev = log(1.0 + tpcResolution);
          //double stdDev = sigFractionGuess;
          double prob = gaus(mean,stdDev,xPoint);
          cout << "        Mu:" << mean << " sig:" << stdDev << " prob:" << prob << endl;
        }else if(m_currentParticleTags[partIndex2] == 7){
          //initialization for matrix method is not going to work for this perfectly, but that's okay
          //these are not seeds to the fit
          double mean = log(m_bichselCurves[partIndex2]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double stdDev = log(1.0 + tpcResolution);
          //double stdDev = sigFractionGuess;// change 40 lines below too
          if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
            mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }else{
            if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
              mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }
            if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
              mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }

          }
          double prob = gaus(mean,stdDev,xPoint);
          cout << "        Mu:" << mean << " sig:" << stdDev << " prob:" << prob << endl;
          //if(prob < 1E-8) matrixElements.push_back( 0.0 );
          //else
          matrixElements.push_back( prob );

        }else if(m_currentParticleTags[partIndex2] == 8){
          //initialization for matrix method is not going to work for this perfectly, but that's okay
          //these are not seeds to the fit
          double mean = log(m_bichselCurves[partIndex2]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double stdDev = log(1.0 + tpcResolution);
          //double stdDev = sigFractionGuess;// change 40 lines below too
          if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
            mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }else{
            if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
              mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }
            if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
              mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }

          }
          double prob = gaus(mean,stdDev,xPoint);
          cout << "        Mu:" << mean << " sig:" << stdDev << " prob:" << prob << endl;

          matrixElements.push_back( prob );

        }
      }
    }



    cout<< "Calculating initial amplitudes..." << endl;
    cout << "Matrix STDVectorSize: " << matrixElements.size() << " NumPart: " << m_currentNumPartToFit << endl;
    TMatrixD mat(m_currentNumPartToFit, m_currentNumPartToFit, &matrixElements[0]);
    double determ = mat.Determinant();
    cout << "MATRIX: Rows = Same dE/dx point Columns = Same Particle with Det:" << determ << endl;
    mat.Print();

    if(determ > 0.0 ){
      TDecompLU lu(mat);
      TVectorD countVect(m_currentNumPartToFit);
      for(int iii = 0; iii < m_currentNumPartToFit; iii++){
        countVect(iii) = countsVector[iii];
      }
      cout << "Counts at specific mean locations: " << endl;
      countVect.Print();

      //TVectorD ampVect = matInv * countVect;
      bool goodSolve = false;
      ampVect=lu.Solve(countVect,goodSolve);
      cout << "Initial Amplitude Vector (size="<< ampVect.GetNrows() << ") " << endl;
      ampVect.Print();

      cout << "This should be zero:" << endl;
      TVectorD zeroVect = (mat*ampVect)-countVect;
      zeroVect.Print();
    }else{
      matrixMethodFailed = true;
    }
  }

  if(matrixMethodFailed){
   cout << "Matrix Method Failed!!!!" << endl;
 }



 double counts = m_currentHistoToFit->Integral("width");

 int partIndexCounter = 0;
 for(int partIndex = 0; partIndex < m_numParticles; partIndex++){

  if(m_currentParticleTags[partIndex] == 0) continue;

  if(m_currentParticleTags[partIndex] == 1){
    //AMPLITUDE FREE
    double mu = -1;
    double sig = -1;
    mu = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
    //Take weighted average by (N_p*mu_p + N_m*mu_m)(N_p + N_m)
    if(m_weightedAverage_mean[partIndex]) mu = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
    sig = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
    if(m_weightedAverage_stdDev[partIndex]) sig = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));


    #ifdef _ZFITTER_DEBUG_
      cout << " setup part " << partIndex << " for tag 1:   Mu: " << mu << "  Sig: " << sig << endl;
    #endif
    double amp = 0.0;
    if(m_useMatrixInitializationForTPC && !matrixMethodFailed){
      amp = ampVect(partIndexCounter)/m_currentIntegralofHistogram;
    }else{
      int predictBin = m_currentHistoToFit->FindBin(mu);
      //double previousAmpContributionFromOtherParts = 0.0;
      //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
      //  previousAmpContributionFromOtherParts += params[paramIndex] * gaus(params[paramIndex+1], params[paramIndex+2], mu);
      //}
      //amp = (m_currentHistoToFit->GetBinContent(predictBin) - previousAmpContributionFromOtherParts)  / gaus(mu,sig,mu);
      amp = m_currentHistoToFit->GetBinContent(predictBin)  / (gaus(mu,sig,mu) * m_currentIntegralofHistogram) ;
      if(std::isnan(amp)){
        cout << "SOMETHING WRONG:: amp Infinite : " << amp << " : for particle " << partIndex << " : Setting amp to ";
        amp = counts / ((double) m_currentNumPartToFit);
        cout << amp << endl;
        cout << "      Mu: " << mu << " predictBin: " << predictBin << " content: " << m_currentHistoToFit->GetBinContent(predictBin)  << " gaussProb: " << gaus(mu,sig,mu) << endl;
      }
    }

    if(std::isnan(amp)){
      cout << "SOMETHING WRONG:: amp Infinite" << endl;
    }
      //cout << "      Amp= " << amp << endl;
      //cout << "      Mean= " << mean << endl;
      //cout << "      Std= " << stdDev << endl;
    params.push_back(amp);
    paramToggles.push_back(true);
    paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

    params.push_back(mu);
    paramToggles.push_back(false);

    params.push_back(sig);
    paramToggles.push_back(false);
    if(sig==0){
      cout << "Error: sig is zero for solo dEdx" << endl;
    }
      //m_currentNumPartToFit++;
    m_currentFreeParameters++;
    partIndexCounter++;
  }else if(m_currentParticleTags[partIndex] == 2){
      //ALL PARAMETERS FREE
      //cout << "  Fitting Part: " << partIndex << endl;
      //cout << "    No Solo dE/dx Data" << partIndex << endl;

      //double mean = m_bichselCurves[partIndex]->Eval(m_currentMomOfBin);
      //double stdDev = 0.07*mean;// change 40 lines above too
      //double mu,sig;
      //PhysMath::gausToLogNorm(mean,stdDev,mu,sig);
      //
        //double mu = 0;
        //double sig = 0;

        //mu = m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[1][m_currentRapBin-1]->Eval(PhysMath::betagamma(m_currentMomOfBin, m_partInfo->GetParticleMass(partIndex)));
        //sig =  m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[2][m_currentRapBin-1]->Eval(PhysMath::betagamma(m_currentMomOfBin, m_partInfo->GetParticleMass(partIndex)));
      //cout << "mu: " << mu << " sig: " << sig << endl;
      //if(std::isnan(mu) || std::isnan(sig)){

      //####### SECTION TO BE IDENTICAL TO ABOVE ###########
      double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double stdDev = log(1.0 + tpcResolution);

      bool use3PtGausAmp = false;
      double threePtGausAmp = 0.0;

      if(m_bichselShiftsComplete){
        if(m_currentPartCharge == 1){
          if(m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin) > 0.0){
            mean += m_bichselShifts_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
            stdDev = m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
          }else if(m_tpcResolution_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin) > 0.0){
            mean += m_bichselShifts_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
            //stdDev = m_tpcResolution_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
          }
        }else{
          if(m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin) > 0.0){
            mean += m_bichselShifts_Minus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
            stdDev = m_tpcResolution_Minus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
          }else if(m_tpcResolution_Minus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin) > 0.0){
            mean += m_bichselShifts_Minus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
            //stdDev = m_tpcResolution_Minus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
          }
        }
      }

      if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
        mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{  //use bTOF isolated fits as a seed if they are there
        if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }
        if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }
      }
      bool do3PtGausCalc = true;
      bool gotOtherParticleSeed = false;
      if(m_useFitInfoForOppositeCharge && m_seedFromOppositeChargeDirection[partIndex] > 0 && m_currentPartCharge < 0.0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          gotOtherParticleSeed = true;
          do3PtGausCalc = false;
      }
      if(m_useFitInfoForOppositeCharge && m_seedFromOppositeChargeDirection[partIndex] < 0 && m_currentPartCharge > 0.0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          gotOtherParticleSeed = true;
          do3PtGausCalc = false;
      }

      if(!gotOtherParticleSeed && m_useMirroredRapidityAsSeed){
        double currentCenter = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
        double mirrorBin = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->FindBin(-currentCenter);
        if(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetXaxis()->GetBinCenter(m_currentRapBin) > 0.0001
          && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin) > 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(mirrorBin,m_currentMtM0Bin);
          do3PtGausCalc = false;
        }
      }


      //####  Try to calculate the initialization parameters by 3 points  ###
      if(do3PtGausCalc){
        double x1, x2, x3, y1 , y2 , y3 = 0.0;
        if(m_currentNumPartToFit == 1){
          x1 = mean - stdDev;
          x2 = mean;
          x3 = mean + stdDev;
        }else{
          if(m_particleLeftRightSeedDirection[partIndex] == -1){
            x1 = mean - stdDev;
            x2 = mean - 0.5*stdDev;
            x3 = mean;
          }else if(m_particleLeftRightSeedDirection[partIndex] == 1){
            x1 = mean + stdDev;
            x2 = mean + 0.5*stdDev;
            x3 = mean;
          }else{
            x1 = mean - stdDev;
            x2 = mean;
            x3 = mean + stdDev;
          }

        }
        int bin1 = m_currentHistoToFit->FindBin(x1);
        int bin2 = m_currentHistoToFit->FindBin(x2);
        int bin3 = m_currentHistoToFit->FindBin(x3);

        #ifdef _ZFITTER_DEBUG_
        cout << " Current Sample Pts:  " << x1 << "   " << x2 << "   " << x3 << endl;
        cout << " Current Sample Bins: " << bin1 << "  " << bin2 << "  " << bin3 << endl;
        #endif

        if(bin1 == bin2 || bin1 == bin3 || bin2 == bin3){
          bin1 = bin2 - 4;
          bin3 = bin2 + 4;
          x1 = m_currentHistoToFit->GetBinCenter(bin1);
          x3 = m_currentHistoToFit->GetBinCenter(bin3);
        }
        if(bin1 > m_currentHistoToFit->GetNbinsX()-2 || bin2 > m_currentHistoToFit->GetNbinsX()-2 || bin3 > m_currentHistoToFit->GetNbinsX()-2
          || bin1 < 1 || bin2 < 1 || bin3 < 1){
          cout << "   Can't initialize with 3 point gaus because bins are off the edge of the histogram range" << endl;
        }else{
          y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, bin1, 2);
          y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, bin2, 2);
          y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, bin3, 2);
          double N_init, mu_init, sig_init = 0.0;
          PhysMath::get3PointGausParams(N_init, mu_init, sig_init,x1,y1,x2,y2,x3,y3);
          if(sig_init <= 0.0 || std::isnan(N_init) || std::isnan(mu_init) || std::isnan(sig_init) || fabs(mu_init-mean) > 0.08 ){
            //Initialization didn't work / don't overwrite previous initialization
            cout << " 3 point Gaussian failed - no change in parameters" << endl;
          }else{
            //initialization worked
            use3PtGausAmp = true;
            threePtGausAmp = N_init;
            mean = mu_init;
            stdDev = sig_init;
          }
        }
      }// end 3 point gaus calculation






      //#################### END OF SECTION ################

      double amp = 0.0;
      if(m_useMatrixInitializationForTPC && !matrixMethodFailed){
        amp = ampVect(partIndexCounter)/m_currentIntegralofHistogram;
      }else{
        int predictBin = m_currentHistoToFit->FindBin(mean);

        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
        //  previousAmpContributionFromOtherParts += params[paramIndex]*gaus(params[paramIndex+1], params[paramIndex+2], mean);
        //}
        //amp = (m_currentHistoToFit->GetBinContent(predictBin) - previousAmpContributionFromOtherParts)  / gaus(mean,stdDev,mean);
        amp = m_currentHistoToFit->GetBinContent(predictBin)  / (gaus(mean,stdDev,mean) * m_currentIntegralofHistogram);
        if(std::isnan(amp)){
          cout << "SOMETHING WRONG:: amp Infinite : " << amp << " forcing it to ";
          amp = counts / ((double) m_currentNumPartToFit);
          cout << amp << endl;
        }
      }

      if(use3PtGausAmp) amp = threePtGausAmp/ m_currentIntegralofHistogram;

      #ifdef _ZFITTER_DEBUG_
        cout << " Setup Particle " << partIndex << " for Tag=2:   Mean: " << mean << "  StdDev: " << stdDev << " Amp: " << amp << endl;
      #endif

      params.push_back(amp);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      params.push_back(mean);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      params.push_back(stdDev);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      //m_currentNumPartToFit++;
      m_currentFreeParameters+=3;
      partIndexCounter++;
    }else if(m_currentParticleTags[partIndex] == 3){
      // ALL PARAMETERS FIXED

      double mu = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(m_weightedAverage_mean[partIndex]) mu = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      double sig = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      if(m_weightedAverage_stdDev[partIndex]) sig = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));

      params.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) / m_currentIntegralofHistogram);
      params.push_back(mu);
      params.push_back(sig);

      paramToggles.push_back(false);
      paramToggles.push_back(false);
      paramToggles.push_back(false);
      partIndexCounter++;
    }else if(m_currentParticleTags[partIndex] == 4){
      params.push_back(0.0);
      params.push_back(0.0);
      params.push_back(0.0);
      paramToggles.push_back(false);
      paramToggles.push_back(false);
      paramToggles.push_back(false);
      partIndexCounter++;
    }else if(m_currentParticleTags[partIndex] == 5){
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));
      params.push_back(0.0);
      params.push_back(0.0);
      params.push_back(0.0);
      paramToggles.push_back(true);
      paramToggles.push_back(false);
      paramToggles.push_back(false);
      m_currentFreeParameters++;
      partIndexCounter++;
    }else if(m_currentParticleTags[partIndex] == 6){
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));
      params.push_back(0.0);
      params.push_back(0.0);
      params.push_back(0.0);
      paramToggles.push_back(true);
      paramToggles.push_back(false);
      paramToggles.push_back(false);
      partIndexCounter++;
      m_currentFreeParameters+=1;
    }else if(m_currentParticleTags[partIndex] == 7){
      double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double stdDev = log(1.0 + tpcResolution);
      if(m_bichselShiftsComplete){
        if(m_currentPartCharge == 1){
          if(m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin) != 0.0){
            mean += m_bichselShifts_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
            stdDev = m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
          }else{
            mean += m_bichselShifts_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
            stdDev = m_tpcResolution_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
          }
        }else{
          if(m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin) != 0.0){
            mean += m_bichselShifts_Minus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
            stdDev = m_tpcResolution_Minus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
          }else{
            mean += m_bichselShifts_Minus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
            stdDev = m_tpcResolution_Minus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
          }
        }
      }else{
        if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }else{  //use bTOF isolated fits as a seed if they are there
          if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
            mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }
          if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
            mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }
        }
      }

      params.push_back(0.0);
      paramToggles.push_back(false);

      params.push_back(mean);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      params.push_back(stdDev);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      //m_currentNumPartToFit++;
      m_currentFreeParameters+=2;
      partIndexCounter++;
    }else if(m_currentParticleTags[partIndex] == 8){
     double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
     double stdDev = log(1.0 + tpcResolution);
     //double stdDev = sigFractionGuess;// change 40 lines above too

      if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
        mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{  //use bTOF isolated fits as a seed if they are there
        if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

        if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

      }

      if(m_weightedAverage_mean[partIndex]) mean = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      if(m_weightedAverage_stdDev[partIndex]) stdDev = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));

      //#################### END OF SECTION ################

      params.push_back(0.0);
      paramToggles.push_back(false);

      params.push_back(mean);
      paramToggles.push_back(false);

      params.push_back(stdDev);
      paramToggles.push_back(false);

      //m_currentNumPartToFit++;
      partIndexCounter++;
    }else if(m_currentParticleTags[partIndex] == 9){
      //AMPLITUDE FREE
      double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double stdDev = sigFractionGuess;// change 40 lines above too

      if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
        mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{  //use bTOF isolated fits as a seed if they are there
        if(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }
      }

      if(m_weightedAverage_mean[partIndex]) mean = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));

      double amp = 0.0;
      if(m_useMatrixInitializationForTPC && !matrixMethodFailed){
        amp = ampVect(partIndexCounter)/m_currentIntegralofHistogram;
      }else{
        int predictBin = m_currentHistoToFit->FindBin(mean);
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
        //  previousAmpContributionFromOtherParts += params[paramIndex]*gaus(params[paramIndex+1],params[paramIndex+2],mean);
        //}
        //amp = (m_currentHistoToFit->GetBinContent(predictBin) - previousAmpContributionFromOtherParts)  / gaus(mean,stdDev,mean);
        amp = m_currentHistoToFit->GetBinContent(predictBin)  / (gaus(mean,stdDev,mean) * m_currentIntegralofHistogram);
        if(std::isnan(amp)){
          cout << "SOMETHING WRONG:: amp Infinite" << endl;
          amp = counts / (((double) m_currentNumPartToFit) * m_currentIntegralofHistogram);
        }
      }

      if(std::isnan(amp)){
        cout << "SOMETHING WRONG:: amp Infinite" << endl;
      }
      //cout << "      Amp= " << amp << endl;
      //cout << "      Mean= " << mean << endl;
      //cout << "      Std= " << stdDev << endl;
      params.push_back(amp);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      params.push_back(mean);
      paramToggles.push_back(false);

      params.push_back(stdDev);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));
      if(stdDev==0){
        cout << "Error: sig is zero for solo dEdx" << endl;
      }
      //m_currentNumPartToFit++;
      m_currentFreeParameters+=2;
      partIndexCounter++;

    }else if(m_currentParticleTags[partIndex] == 10){
      //AMPLITUDE FREE
      double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double stdDev = sigFractionGuess;// change 40 lines above too

      if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
        mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{  //use bTOF isolated fits as a seed if they are there
        if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }
        if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

      if(m_weightedAverage_mean[partIndex]) mean = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));

      }
      double amp = 0.0;
      if(m_useMatrixInitializationForTPC && !matrixMethodFailed){
        amp = ampVect(partIndexCounter)/m_currentIntegralofHistogram;
      }else{
        int predictBin = m_currentHistoToFit->FindBin(mean);
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
        //  previousAmpContributionFromOtherParts += params[paramIndex]*gaus(params[paramIndex+1],params[paramIndex+2],mean);
        //}
        //amp = (m_currentHistoToFit->GetBinContent(predictBin) - previousAmpContributionFromOtherParts)  / gaus(mean,stdDev,mean);
        amp = m_currentHistoToFit->GetBinContent(predictBin)  / (gaus(mean,stdDev,mean) * m_currentIntegralofHistogram);
        if(std::isnan(amp)){
          cout << "SOMETHING WRONG:: amp Infinite" << endl;
          amp = counts / (((double) m_currentNumPartToFit) * m_currentIntegralofHistogram);
        }
      }

      if(std::isnan(amp)){
        cout << "SOMETHING WRONG:: amp Infinite" << endl;
      }
      //cout << "      Amp= " << amp << endl;
      //cout << "      Mean= " << mean << endl;
      //cout << "      Std= " << stdDev << endl;
      params.push_back(amp);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      params.push_back(mean);
      paramToggles.push_back(false);

      params.push_back(stdDev);
      paramToggles.push_back(false);
      if(stdDev==0){
        cout << "Error: sig is zero for solo dEdx" << endl;
      }
      //m_currentNumPartToFit++;
      m_currentFreeParameters++;
      partIndexCounter++;

    }else if(m_currentParticleTags[partIndex] == 11){
      //AMPLITUDE FREE
      double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double stdDev = sigFractionGuess;// change 40 lines above too

      if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
        mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{  //use bTOF isolated fits as a seed if they are there
        if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

        if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

      }
      double amp = 0.0;
      if(m_useMatrixInitializationForTPC && !matrixMethodFailed){
        amp = ampVect(partIndexCounter)/m_currentIntegralofHistogram;
      }else{
        int predictBin = m_currentHistoToFit->FindBin(mean);
        //double previousAmpContributionFromOtherParts = 0.0;
        //for(unsigned int paramIndex = 0; paramIndex < params.size(); paramIndex += 3){
        //  previousAmpContributionFromOtherParts += params[paramIndex]*gaus(params[paramIndex+1],params[paramIndex+2],mean);
        //}
        //amp = (m_currentHistoToFit->GetBinContent(predictBin) - previousAmpContributionFromOtherParts)  / gaus(mean,stdDev,mean);
        amp = m_currentHistoToFit->GetBinContent(predictBin)  / (gaus(mean,stdDev,mean) * m_currentIntegralofHistogram);
        if(std::isnan(amp)){
          cout << "SOMETHING WRONG:: amp Infinite" << endl;
          amp = counts / (((double) m_currentNumPartToFit) * m_currentIntegralofHistogram);
        }
      }

      if(std::isnan(amp)){
        cout << "SOMETHING WRONG:: amp Infinite" << endl;
      }
      //cout << "      Amp= " << amp << endl;
      //cout << "      Mean= " << mean << endl;
      //cout << "      Std= " << stdDev << endl;
      params.push_back(amp);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      params.push_back(mean);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data()));

      params.push_back(stdDev);
      paramToggles.push_back(false);
      if(stdDev==0){
        cout << "Error: sig is zero for solo dEdx" << endl;
      }
      //m_currentNumPartToFit++;
      m_currentFreeParameters+=2;
      partIndexCounter++;

    }else if(m_currentParticleTags[partIndex] == 12){
     double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double stdDev = sigFractionGuess;// change 40 lines above too

      if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
        mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{  //use bTOF isolated fits as a seed if they are there
        if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

        if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

      }

      if(m_weightedAverage_mean[partIndex]) mean = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));

      //#################### END OF SECTION ################

      params.push_back(0.0);
      paramToggles.push_back(false);

      params.push_back(mean);
      paramToggles.push_back(false);

      params.push_back(stdDev);
      paramToggles.push_back(false);

      //m_currentNumPartToFit++;
      partIndexCounter++;
    }else if(m_currentParticleTags[partIndex] == 13){
     double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double stdDev = sigFractionGuess;// change 40 lines above too

      if(m_usePionMinusFitAsSeedForPionPlus && m_currentPartCharge == 1 && partIndex == 0 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
        mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      }else{  //use bTOF isolated fits as a seed if they are there
        if(m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

        if(m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
          mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          stdDev = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }

      }

      if(m_weightedAverage_mean[partIndex]) mean = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      if(m_weightedAverage_stdDev[partIndex]) stdDev = (m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)*m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin))/(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) + m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      //#################### END OF SECTION ################

      params.push_back(0.0);
      paramToggles.push_back(false);

      params.push_back(mean);
      paramToggles.push_back(false);

      params.push_back(stdDev);
      paramToggles.push_back(false);

      //m_currentNumPartToFit++;
      partIndexCounter++;
    }






    //Edit electron aplitude
    if(partIndex == 3 && m_currentParticleTags[3] != 0 && m_currentParticleTags[0] != 0 && m_electronStudyComplete ){
        cout << " Edited Electron amplitude from : " << params[params.size()-3] << endl;
        if(params[params.size()-3] > 1.5*params[0]*m_basicElectronPionRatio) params[params.size()-3] = params[0]*m_basicElectronPionRatio;
        cout << " Edited Electron amplitude to  : " << params[params.size()-3] << endl;
    }

    //Edit kaon aplitude
    if(partIndex == 1 && m_currentParticleTags[1] != 0 && m_currentParticleTags[0] != 0 ){
        cout << " Edited Kaon amplitude from : " << params[params.size()-3] << endl;
        if(params[params.size()-3] > 1.5*params[0]*m_basicKaonPionRatio) params[params.size()-3] = params[0]*m_basicKaonPionRatio;
        cout << " Edited Kaon amplitude to : " << params[params.size()-3] << endl;
    }
  }//end particle loop

  m_currentParameters.clear();
  m_paramToggles.clear();
  m_currentParamNames.clear();

  m_currentParameters = params;
  m_paramToggles = paramToggles;
  m_currentParamNames = paramNames;

  partIndexCounter = 0;
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){

    if(m_currentParticleTags[partIndex] == 4){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double ampOfMimic = m_currentParameters[partMimicLoc];
        double meanOfMimic = m_currentParameters[partMimicLoc+1];
        double stdDevOfMimic = m_currentParameters[partMimicLoc+2];

        //constrain the amplitude
        TH2D* tofSpectraMimic = NULL;
        TH2D* tofSpectra = NULL;
        if(m_inAllCentMode && !m_separatePM_bToF_allCent){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0]->GetEntries() != 0 && m_currentPartCharge == 1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1]->GetEntries() != 0 && m_currentPartCharge == -1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1];
        }else{
          if(m_currentPartCharge == 1){
            tofSpectraMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][m_currentCentIndex];
          }else{
            tofSpectraMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][m_currentCentIndex];
          }
        }
        double tofAmpMimic = tofSpectraMimic->GetBinContent(tofSpectraMimic->GetXaxis()->FindBin(m_currentRapidity), tofSpectraMimic->GetYaxis()->FindBin(m_currentMtM0));
        double tofAmp =  tofSpectra->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
        //if(tofAmpMimic == 0 || tofAmp == 0){
        //  cout << "Could not get TOF Amplitudes to constrain in dEdx function. Rap: " << m_currentRapidity << " mTm0: " << m_currentMtM0 << endl;
        //  m_currentConstrainedAmpRaios[partIndex] = 0.0;
        //  m_currentParameters[partLoc] = 0;
        //  m_currentParameters[partLoc+1] = 0;
        //  m_currentParameters[partLoc+2] = 0.025;
        //  m_currentParticleTags[partIndex] = 3;
        //  continue;
        //}
        if(partIndex == 1 && m_doKaonDecay){ //if kaon is mimicing other particle, need to account for kaon decays
          double pctKaonLeft = PhysMath::calc_pct_kaon_left(m_currentRapidity,m_currentMtM0,m_partInfo->GetParticleMass(1));
          cout << "original kaon amp: " << tofAmp << endl;
          tofAmp = tofAmp/pctKaonLeft;
          cout << "decay corrected tof amp: " << tofAmp << endl;
        }
        m_currentConstrainedAmpRaios[partIndex] = tofAmp/tofAmpMimic;
        double amp = ampOfMimic*m_currentConstrainedAmpRaios[partIndex];
        cout << "TofAmpMimic: " << tofAmpMimic << "  tofAmp: " << tofAmp << "   ratio: " << m_currentConstrainedAmpRaios[partIndex] << " newAmp: " << amp << endl;

        //constrain the mean
        double ZTPCShiftOfMimic = meanOfMimic -  log(m_bichselCurves[partIndexToMimic]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + ZTPCShiftOfMimic;
        cout << "ZTPCShiftOfMimic " << ZTPCShiftOfMimic << "  new Mean " << mean << endl;

        //constrain the width
        double stdDev = stdDevOfMimic;
        //double fracStdDev = stdDevOfMimic / (meanOfMimic + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)));
        //double stdDev = fracStdDev * (mean + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)));
        //cout << "Frac Std DeV " << fracStdDev << " new stdDev " << stdDev << endl;

        m_currentParameters[partLoc] = amp;
        m_currentParameters[partLoc+1] = mean;
        m_currentParameters[partLoc+2] = stdDev;

      }
    }// end of tag 4
    else if(m_currentParticleTags[partIndex] == 5){

      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{

        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double meanOfMimic = m_currentParameters[partMimicLoc+1];

        double mu = -1;
        double sig = -1;
        if(m_currentPartCharge == 1){
          mu = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          sig = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }else{
          mu = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          sig = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
        }
        double amp = 0.0;
        if(m_useMatrixInitializationForTPC && !matrixMethodFailed){
          amp = ampVect(partIndexCounter);
        }else{
          int predictBin = m_currentHistoToFit->FindBin(mu);
          amp = m_currentHistoToFit->GetBinContent(predictBin) / gaus(mu,sig,mu);
          if(std::isnan(amp)){
            cout << "SOMETHING WRONG:: amp Infinite" << endl;
            amp = counts / ((double) m_currentNumPartToFit);
          }
        }
        if(std::isnan(amp)){
          cout << "SOMETHING VERY WRONG:: amp Infinite" << endl;
        }

        //constrain the mean
        double ZTPCShiftOfMimic = meanOfMimic -  log(m_bichselCurves[partIndexToMimic]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        double mean = mu + ZTPCShiftOfMimic;
        cout << "ZTPCShiftOfMimic " << ZTPCShiftOfMimic << "  new Mean " << mean << endl;


        m_currentParameters[partLoc] = amp;
        m_currentParameters[partLoc+1] = mean;
        m_currentParameters[partLoc+2] = sig;
        m_currentMeansForTagFive[partIndex] = mu;

      }
    } //end of tag 5
    else if(m_currentParticleTags[partIndex] == 6){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double ampOfMimic = m_currentParameters[partMimicLoc];
        double meanOfMimic = m_currentParameters[partMimicLoc+1];
        double stdDevOfMimic = m_currentParameters[partMimicLoc+2];



        //constrain the mean
        double ZTPCShiftOfMimic = meanOfMimic -  log(m_bichselCurves[partIndexToMimic]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
        double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + ZTPCShiftOfMimic;
        cout << "ZTPCShiftOfMimic " << ZTPCShiftOfMimic << "  new Mean " << mean << endl;

        //constrain the width
        double stdDev = stdDevOfMimic;
        //double fracStdDev = stdDevOfMimic / (meanOfMimic + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)));
        //double stdDev = fracStdDev * (mean + log(m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)));
        //cout << "Frac Std DeV " << fracStdDev << " new stdDev " << stdDev << endl;


        double amp = 0.0;
        if(m_useMatrixInitializationForTPC && !matrixMethodFailed){
          amp = ampVect(partIndexCounter);
        }else{
          int predictBin = m_currentHistoToFit->FindBin(mean);
          amp = m_currentHistoToFit->GetBinContent(predictBin) / gaus(mean,stdDev,mean);
          if(std::isnan(amp)){
            cout << "SOMETHING WRONG:: amp Infinite" << endl;
            amp = counts / ((double) m_currentNumPartToFit);
          }
        }
        if(std::isnan(amp)){
          cout << "SOMETHING VERY WRONG:: amp Infinite" << endl;
        }
        m_currentParameters[partLoc] = amp;
        m_currentParameters[partLoc+1] = mean;
        m_currentParameters[partLoc+2] = stdDev;



      }
    } //end of tag 6
    else if(m_currentParticleTags[partIndex] == 7){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double ampOfMimic = m_currentParameters[partMimicLoc];

        //constrain the amplitude
        TH2D* tofSpectraMimic = NULL;
        TH2D* tofSpectra = NULL;
        if(m_inAllCentMode && !m_separatePM_bToF_allCent){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0]->GetEntries() != 0 && m_currentPartCharge == 1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1]->GetEntries() != 0 && m_currentPartCharge == -1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1];

        }else{
          if(m_currentPartCharge == 1){
            tofSpectraMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][m_currentCentIndex];
          }else{
            tofSpectraMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][m_currentCentIndex];
          }
        }
        double tofAmpMimic = tofSpectraMimic->GetBinContent(tofSpectraMimic->GetXaxis()->FindBin(m_currentRapidity), tofSpectraMimic->GetYaxis()->FindBin(m_currentMtM0));
        double tofAmp =  tofSpectra->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
        //if(tofAmpMimic == 0 || tofAmp == 0){
        //  cout << "Could not get TOF Amplitudes to constrain in dEdx function. Rap: " << m_currentRapidity << " mTm0: " << m_currentMtM0 << endl;
        //  m_currentConstrainedAmpRaios[partIndex] = 0.0;
        //  m_currentParameters[partLoc] = 0;
        //  m_currentParameters[partLoc+1] = 0;
        //  m_currentParameters[partLoc+2] = 0.025;
        //  m_currentParticleTags[partIndex] = 3;
        //  continue;
        //}
        if(partIndex == 1 && m_doKaonDecay){ //if kaon is mimicing other particle, need to account for kaon decays
          double pctKaonLeft = PhysMath::calc_pct_kaon_left(m_currentRapidity,m_currentMtM0,m_partInfo->GetParticleMass(1));
          cout << "original kaon amp: " << tofAmp << endl;
          tofAmp = tofAmp/pctKaonLeft;
          cout << "decay corrected tof amp: " << tofAmp << endl;
        }
        m_currentConstrainedAmpRaios[partIndex] = tofAmp/tofAmpMimic;
        double amp = ampOfMimic*m_currentConstrainedAmpRaios[partIndex];
        cout << "TofAmpMimic: " << tofAmpMimic << "  tofAmp: " << tofAmp << "   ratio: " << m_currentConstrainedAmpRaios[partIndex] << " newAmp: " << amp << endl;

        m_currentParameters[partLoc] = amp;
      }
    }// end of tag 7
    else if(m_currentParticleTags[partIndex] == 8){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double ampOfMimic = m_currentParameters[partMimicLoc];

        //constrain the amplitude
        TH2D* tofSpectraMimic = NULL;
        TH2D* tofSpectra = NULL;
        if(m_inAllCentMode && !m_separatePM_bToF_allCent){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0]->GetEntries() != 0 && m_currentPartCharge == 1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1]->GetEntries() != 0 && m_currentPartCharge == -1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1];
        }else{
          if(m_currentPartCharge == 1){
            tofSpectraMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][m_currentCentIndex];
          }else{
            tofSpectraMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][m_currentCentIndex];
          }
        }
        double tofAmpMimic = tofSpectraMimic->GetBinContent(tofSpectraMimic->GetXaxis()->FindBin(m_currentRapidity), tofSpectraMimic->GetYaxis()->FindBin(m_currentMtM0));
        double tofAmp =  tofSpectra->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
        //if(tofAmpMimic == 0 || tofAmp == 0){
        //  cout << "Could not get TOF Amplitudes to constrain in dEdx function. Rap: " << m_currentRapidity << " mTm0: " << m_currentMtM0 << endl;
        //  m_currentConstrainedAmpRaios[partIndex] = 0.0;
        //  m_currentParameters[partLoc] = 0;
        //  m_currentParameters[partLoc+1] = 0;
        //  m_currentParameters[partLoc+2] = 0.025;
        //  m_currentParticleTags[partIndex] = 3;
        //  continue;
        //}
        if(partIndex == 1 && m_doKaonDecay){ //if kaon is mimicing other particle, need to account for kaon decays
          double pctKaonLeft = PhysMath::calc_pct_kaon_left(m_currentRapidity,m_currentMtM0,m_partInfo->GetParticleMass(1));
          cout << "original kaon amp: " << tofAmp << endl;
          tofAmp = tofAmp/pctKaonLeft;
          cout << "decay corrected tof amp: " << tofAmp << endl;
        }
        m_currentConstrainedAmpRaios[partIndex] = tofAmp/tofAmpMimic;
        double amp = ampOfMimic*m_currentConstrainedAmpRaios[partIndex];
        cout << "TofAmpMimic: " << tofAmpMimic << "  tofAmp: " << tofAmp << "   ratio: " << m_currentConstrainedAmpRaios[partIndex] << " newAmp: " << amp << endl;

        m_currentParameters[partLoc] = amp;
      }
    }// end of tag 8
    else if(m_currentParticleTags[partIndex] == 10){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double stdDevOfMimic = m_currentParameters[partMimicLoc+2];


        double stdDev = stdDevOfMimic;
        m_currentParameters[partLoc+2] = stdDev;
      }
    }
    else if(m_currentParticleTags[partIndex] == 11){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double stdDevOfMimic = m_currentParameters[partMimicLoc+2];


        double stdDev = stdDevOfMimic;
        m_currentParameters[partLoc+2] = stdDev;
      }
    }
    else if(m_currentParticleTags[partIndex] == 12){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      //int partIndexToMimicMimic = m_particleMimicIndicies[partIndexToMimic];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double ampOfMimic = m_currentParameters[partMimicLoc];
        //if(!m_inAllCentMode) ampOfMimic = m_currentParameters[partMimicLoc]*m_currentIntegralofHistogram;
        double stdDevOfMimic = m_currentParameters[partMimicLoc+2];

        //constrain the amplitude
        TH2D* tofSpectraMimic = NULL;
        //TH2D* tofSpectraMimicMimic = NULL;
        TH2D* tofSpectra = NULL;
        if(m_inAllCentMode && !m_separatePM_bToF_allCent){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
          //if(partIndexToMimicMimic != partIndexToMimic) tofSpectraMimicMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0]->GetEntries() != 0 && m_currentPartCharge == 1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
          //if(partIndexToMimicMimic != partIndexToMimic && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0]->GetEntries() != 0) tofSpectraMimicMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1]->GetEntries() != 0 && m_currentPartCharge == -1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1];
          //if(partIndexToMimicMimic != partIndexToMimic && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0]->GetEntries() != 0) tofSpectraMimicMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0];
        }else{
          if(m_currentPartCharge == 1){
            tofSpectraMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][m_currentCentIndex];
            //if(partIndexToMimicMimic != partIndexToMimic) tofSpectraMimicMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndexToMimicMimic][m_currentCentIndex];
          }else{
            tofSpectraMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][m_currentCentIndex];
            //if(partIndexToMimicMimic != partIndexToMimic) tofSpectraMimicMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndexToMimicMimic][m_currentCentIndex];
          }
        }
        double tofAmpMimic = tofSpectraMimic->GetBinContent(tofSpectraMimic->GetXaxis()->FindBin(m_currentRapidity), tofSpectraMimic->GetYaxis()->FindBin(m_currentMtM0));
        double tofAmp =  tofSpectra->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
        //double tofAmpMimicMimic;
        //if(tofSpectraMimicMimic != NULL) tofAmpMimicMimic = tofSpectraMimicMimic->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
        //if(tofAmpMimic == 0 || tofAmp == 0){
        //  cout << "Could not get TOF Amplitudes to constrain in dEdx function. Rap: " << m_currentRapidity << " mTm0: " << m_currentMtM0 << endl;
        //  m_currentConstrainedAmpRaios[partIndex] = 0.0;
        //  m_currentParameters[partLoc] = 0;
        //  m_currentParameters[partLoc+1] = 0;
        //  m_currentParameters[partLoc+2] = 0.025;
        //  m_currentParticleTags[partIndex] = 3;
        //  continue;
        //}
        if(partIndex == 1 && m_doKaonDecay){ //if kaon is mimicing other particle, need to account for kaon decays
          double pctKaonLeft = PhysMath::calc_pct_kaon_left(m_currentRapidity,m_currentMtM0,m_partInfo->GetParticleMass(1));
          cout << "original kaon amp: " << tofAmp << endl;
          tofAmp = tofAmp/pctKaonLeft;
          cout << "decay corrected tof amp: " << tofAmp << endl;
        }
        m_currentConstrainedAmpRaios[partIndex] = tofAmp/tofAmpMimic;
        cout << "m_currentIntegralofHistogram: " << m_currentIntegralofHistogram << endl;
        cout << "ampOfMimic: " << ampOfMimic << endl;
        double amp = ampOfMimic*m_currentConstrainedAmpRaios[partIndex];
        cout << "TofAmpMimic: " << tofAmpMimic << "  tofAmp: " << tofAmp << "   ratio: " << m_currentConstrainedAmpRaios[partIndex] << " newAmp: " << amp << endl;

        double stdDev = stdDevOfMimic;

        //if(m_particleMimicIndicies[partIndexToMimic] != partIndexToMimic) m_currentConstrainedAmpRaios[partIndex] = m_currentConstrainedAmpRaios[partIndex]*tofAmpMimic/tofAmpMimicMimic;


        m_currentParameters[partLoc] = amp;//m_currentIntegralofHistogram;
        //if(!m_inAllCentMode) m_currentParameters[partLoc] = amp/m_currentIntegralofHistogram;
        m_currentParameters[partLoc+2] = stdDev;
      }
    }// end of tag 12
    else if(m_currentParticleTags[partIndex] == 13){
      int partIndexToMimic = m_particleMimicIndicies[partIndex];
      //int partIndexToMimicMimic = m_particleMimicIndicies[partIndexToMimic];
      if(m_currentParticleTags[partIndexToMimic] == 0){
        cout << "Something Very Wrong, cant mimic " << m_particleMimicIndicies[partIndex] << " for " << partIndex << endl;
        return;
      }else{
        int partLoc = getCurrParamIndex(partIndex);
        int partMimicLoc = getCurrParamIndex(partIndexToMimic);
        double ampOfMimic = m_currentParameters[partMimicLoc];
        //if(!m_inAllCentMode) ampOfMimic = m_currentParameters[partMimicLoc]*m_currentIntegralofHistogram;
        double stdDevOfMimic = m_currentParameters[partMimicLoc+2];

        //constrain the amplitude
        TH2D* tofSpectraMimic = NULL;
        //TH2D* tofSpectraMimicMimic = NULL;
        TH2D* tofSpectra = NULL;
        if(m_inAllCentMode && !m_separatePM_bToF_allCent){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
          //if(partIndexToMimicMimic != partIndexToMimic) tofSpectraMimicMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0]->GetEntries() != 0 && m_currentPartCharge == 1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][0];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][0];
          //if(partIndexToMimicMimic != partIndexToMimic && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0]->GetEntries() != 0) tofSpectraMimicMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0];
        }else if(m_inAllCentMode && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1]->GetEntries() != 0 && m_currentPartCharge == -1){
          tofSpectraMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimic][0][1];
          tofSpectra = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndex][0][1];
          //if(partIndexToMimicMimic != partIndexToMimic && m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0]->GetEntries() != 0) tofSpectraMimicMimic = m_Fit_Data_ZbTOF[m_currentPartIndex][partIndexToMimicMimic][0][0];
        }else{
          if(m_currentPartCharge == 1){
            tofSpectraMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndex][m_currentCentIndex];
            //if(partIndexToMimicMimic != partIndexToMimic) tofSpectraMimicMimic = m_InvBetaBTOFYieldsCentPlus[m_currentPartIndex][partIndexToMimicMimic][m_currentCentIndex];
          }else{
            tofSpectraMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndexToMimic][m_currentCentIndex];
            tofSpectra = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndex][m_currentCentIndex];
            //if(partIndexToMimicMimic != partIndexToMimic) tofSpectraMimicMimic = m_InvBetaBTOFYieldsCentMinus[m_currentPartIndex][partIndexToMimicMimic][m_currentCentIndex];
          }
        }
        double tofAmpMimic = tofSpectraMimic->GetBinContent(tofSpectraMimic->GetXaxis()->FindBin(m_currentRapidity), tofSpectraMimic->GetYaxis()->FindBin(m_currentMtM0));
        double tofAmp =  tofSpectra->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
        //double tofAmpMimicMimic;
        //if(tofSpectraMimicMimic != NULL) tofAmpMimicMimic = tofSpectraMimicMimic->GetBinContent(tofSpectra->GetXaxis()->FindBin(m_currentRapidity), tofSpectra->GetYaxis()->FindBin(m_currentMtM0));
        //if(tofAmpMimic == 0 || tofAmp == 0){
        //  cout << "Could not get TOF Amplitudes to constrain in dEdx function. Rap: " << m_currentRapidity << " mTm0: " << m_currentMtM0 << endl;
        //  m_currentConstrainedAmpRaios[partIndex] = 0.0;
        //  m_currentParameters[partLoc] = 0;
        //  m_currentParameters[partLoc+1] = 0;
        //  m_currentParameters[partLoc+2] = 0.025;
        //  m_currentParticleTags[partIndex] = 3;
        //  continue;
        //}
        if(partIndex == 1 && m_doKaonDecay){ //if kaon is mimicing other particle, need to account for kaon decays
          double pctKaonLeft = PhysMath::calc_pct_kaon_left(m_currentRapidity,m_currentMtM0,m_partInfo->GetParticleMass(1));
          cout << "original kaon amp: " << tofAmp << endl;
          tofAmp = tofAmp/pctKaonLeft;
          cout << "decay corrected tof amp: " << tofAmp << endl;
        }
        m_currentConstrainedAmpRaios[partIndex] = tofAmp/tofAmpMimic;
        cout << "m_currentIntegralofHistogram: " << m_currentIntegralofHistogram << endl;
        cout << "ampOfMimic: " << ampOfMimic << endl;
        double amp = ampOfMimic*m_currentConstrainedAmpRaios[partIndex];
        cout << "TofAmpMimic: " << tofAmpMimic << "  tofAmp: " << tofAmp << "   ratio: " << m_currentConstrainedAmpRaios[partIndex] << " newAmp: " << amp << endl;

        //double stdDev = stdDevOfMimic;

        //if(m_particleMimicIndicies[partIndexToMimic] != partIndexToMimic) m_currentConstrainedAmpRaios[partIndex] = m_currentConstrainedAmpRaios[partIndex]*tofAmpMimic/tofAmpMimicMimic;


        m_currentParameters[partLoc] = amp;//m_currentIntegralofHistogram;
        //if(!m_inAllCentMode) m_currentParameters[partLoc] = amp/m_currentIntegralofHistogram;
        //m_currentParameters[partLoc+2] = stdDev;
      }
    }// end of tag 12


    if(m_currentParticleTags[partIndex] != 0 ) partIndexCounter++;
  }



}



//============================================================================================================
//============================================================================================================

void ZFitter::initZTPCFunction(TF1* a_funct){
  int paramIndexCounter = 0;
  double integral = m_currentHistoToFit->Integral("width");

  double lowSigmaLimit = 3.0*HistogramUtilities::getSmallestBinWidth(m_currentHistoToFit);
  double highSigmaLimit =0.6*(HistogramUtilities::getHighEdge(m_currentHistoToFit) - HistogramUtilities::getLowEdge(m_currentHistoToFit));


  if(m_allowSpecCentFitsToVary && !m_inAllCentMode){
    for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
      if(m_paramToggles[iii]){
        //cout << " Funct Set Par " <<  paramIndexCounter << " to " <<  m_currentParameters[iii] << endl;
        a_funct->SetParameter(paramIndexCounter, m_currentParameters[iii]);
        a_funct->SetParName(paramIndexCounter, m_currentParamNames[paramIndexCounter].c_str());
        if(iii%3==0){
          if(m_currentParameters[iii] > 0.0 && m_currentParameters[iii] < 2.0){
            //a_funct->SetParLimits(paramIndexCounter,1.0,integral);// 2.0*m_currentParameters[iii]);
            a_funct->SetParLimits(paramIndexCounter,0.0,1.1);
            //a_funct->SetParLimits(paramIndexCounter,0.75*m_currentParameters[iii],1.25*m_currentParameters[iii]);// 2.0*m_currentParameters[iii]);
            //a_funct->SetParLimits(paramIndexCounter,0.1*m_currentParameters[iii],10.0*m_currentParameters[iii]);// 2.0*m_currentParameters[iii]);
            //cout << "  Limits on param " << paramIndexCounter << " are 0 to " << 1000.0*m_currentParameters[iii] << endl;
          }else{
            cout << "ERROR: bad amp init - amplitude is " << m_currentParameters[iii] << endl;
            a_funct->SetParameter(paramIndexCounter, 0.5);
            //a_funct->SetParLimits(paramIndexCounter,1.0, 2.0);
            a_funct->SetParLimits(paramIndexCounter,0.0,1.1);
            //cout << "  Limits on param " << paramIndexCounter << " are 0 to 100" << endl;
          }
        }
        if(iii%3==1) a_funct->SetParLimits(paramIndexCounter,m_currentParameters[iii] - 0.03, m_currentParameters[iii] + 0.03);
        if(iii%3==2) a_funct->SetParLimits(paramIndexCounter,0.97*m_currentParameters[iii],1.03*m_currentParameters[iii]);
        paramIndexCounter++;
      }
    }
  }else{
    for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
      if(m_paramToggles[iii]){
        //cout << " Funct Set Par " <<  paramIndexCounter << " to " <<  m_currentParameters[iii] << endl;
        a_funct->SetParameter(paramIndexCounter, m_currentParameters[iii]);
        a_funct->SetParName(paramIndexCounter, m_currentParamNames[paramIndexCounter].c_str());
        if(iii%3==0){
          double amp = m_currentParameters[iii];
          if(m_currentParameters[iii] > 0.0 && m_currentParameters[iii] < 2.0){
            //a_funct->SetParLimits(paramIndexCounter,1.0,integral);
            a_funct->SetParLimits(paramIndexCounter,0.0,1.1);
            //a_funct->SetParLimits(paramIndexCounter,0.75*m_currentParameters[iii],1.25*m_currentParameters[iii]);
            //a_funct->SetParLimits(paramIndexCounter,0.1*m_currentParameters[iii],10.0*m_currentParameters[iii]);// 2.0*m_currentParameters[iii]);
            //cout << "  Limits on param " << paramIndexCounter << " are 0 to " << 1000.0*m_currentParameters[iii] << endl;
          }else{
            cout << "ERROR: bad amp init - amplitude is " << m_currentParameters[iii] << endl;
            a_funct->SetParameter(paramIndexCounter, 0.5);
            a_funct->SetParLimits(paramIndexCounter,0.0,1.1);
            //cout << "  Limits on param " << paramIndexCounter << " are 0 to 100" << endl;
          }
        }
        if(iii%3==1) a_funct->SetParLimits(paramIndexCounter,m_currentParameters[iii] - 0.1, m_currentParameters[iii] + 0.3);
        if(iii%3==2) a_funct->SetParLimits(paramIndexCounter,0.001,0.6);
        paramIndexCounter++;
      }
    }
  }
}




//============================================================================================================
//============================================================================================================



#ifdef _MAC_OSX_
void ZFitter::initMinimizer_TPC(){
  int paramIndexCounter = 0;
  double integral = m_currentHistoToFit->Integral("width");
  if(m_allowSpecCentFitsToVary && !m_inAllCentMode){
    for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
      if(m_paramToggles[iii]){
        m_minimizer->SetVariable(paramIndexCounter,m_currentParamNames[paramIndexCounter].c_str(), m_currentParameters[iii], 1E-5);
        if(iii%3==0){
          double amp = m_currentParameters[iii];
          m_minimizer->SetVariableLimits(paramIndexCounter,1,integral);
           //cout << " Funct Set Par " <<  paramIndexCounter << " to " <<  m_currentParameters[iii] << endl;
        }
        if(iii%3==1) m_minimizer->SetVariableLimits(paramIndexCounter,m_currentParameters[iii] - 0.02, m_currentParameters[iii] + 0.02);
        if(iii%3==2) m_minimizer->SetVariableLimits(paramIndexCounter,0,1); // was current predict to 0.3
        paramIndexCounter++;
      }
    }
  } else{
    for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
      if(m_paramToggles[iii]){
        m_minimizer->SetVariable(paramIndexCounter,m_currentParamNames[paramIndexCounter].c_str(), m_currentParameters[iii], 1E-5);
        if(iii%3==0){
          double amp = m_currentParameters[iii];
          m_minimizer->SetVariableLimits(paramIndexCounter,1,integral);
           //cout << " Funct Set Par " <<  paramIndexCounter << " to " <<  m_currentParameters[iii] << endl;
        }
        if(iii%3==1) m_minimizer->SetVariableLimits(paramIndexCounter,m_currentParameters[iii] - 0.08, m_currentParameters[iii] + 0.13);
        if(iii%3==2) m_minimizer->SetVariableLimits(paramIndexCounter,0,1);
        paramIndexCounter++;
      }
    }
  }

}
#endif




//==================================================================================================================================
//==================================================================================================================================

vector<double> ZFitter::fitIsolation_dEdx(int a_particleIndex){

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas_iso","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  double lowRange = 0;
  double highRange = 0;

  /*if(a_particleIndex == 0){
    lowRange = 1.5;
    highRange = sqrt(1.2*m_bichselCurves[0]->Eval(m_currentMomOfBin) * m_bichselCurves[1]->Eval(m_currentMomOfBin));
    if(!m_initialParticleToggles[1]){
      highRange = sqrt(1.2*m_bichselCurves[0]->Eval(m_currentMomOfBin) * m_bichselCurves[2]->Eval(m_currentMomOfBin));
    }
  }else if(a_particleIndex == 2){
    lowRange = sqrt(0.8*m_bichselCurves[2]->Eval(m_currentMomOfBin) * m_bichselCurves[1]->Eval(m_currentMomOfBin));
    highRange = sqrt(1.2*m_bichselCurves[2]->Eval(m_currentMomOfBin) * m_bichselCurves[4]->Eval(m_currentMomOfBin));
    if(m_bichselCurves[a_particleIndex-1]->Eval(m_currentMomOfBin) == 0.0) lowRange = 0.7*m_bichselCurves[2]->Eval(m_currentMomOfBin);
    if(m_bichselCurves[a_particleIndex+1]->Eval(m_currentMomOfBin) == 0.0) highRange = 1.3*m_bichselCurves[2]->Eval(m_currentMomOfBin);

  }else if(a_particleIndex == 4){
    lowRange = sqrt(0.8*m_bichselCurves[4]->Eval(m_currentMomOfBin) * m_bichselCurves[2]->Eval(m_currentMomOfBin));
    highRange = sqrt(1.2*m_bichselCurves[4]->Eval(m_currentMomOfBin) * m_bichselCurves[5]->Eval(m_currentMomOfBin));
    if(m_bichselCurves[a_particleIndex-1]->Eval(m_currentMomOfBin) == 0.0) lowRange = 0.7*m_bichselCurves[4]->Eval(m_currentMomOfBin);
    if(m_bichselCurves[a_particleIndex+1]->Eval(m_currentMomOfBin) == 0.0) highRange = 1.3*m_bichselCurves[4]->Eval(m_currentMomOfBin);
  }else{
    lowRange = sqrt(0.8*m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin) * m_bichselCurves[a_particleIndex-1]->Eval(m_currentMomOfBin));
    highRange = sqrt(1.2*m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin) * m_bichselCurves[a_particleIndex+1]->Eval(m_currentMomOfBin));
    if(m_bichselCurves[a_particleIndex-1]->Eval(m_currentMomOfBin) == 0.0) lowRange = 0.7*m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin);
    if(m_bichselCurves[a_particleIndex+1]->Eval(m_currentMomOfBin) == 0.0) highRange = 1.3*m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin);
  }

  if(highRange <= lowRange){
    lowRange = log(m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin));
    highRange = log(m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin));
    //cerr <<  "ZFitter::fitIsolation_dEdx(int a_particleIndex) low and high range are incorrect  lowBound=" << lowRange << "  highBound=" << highRange << " with mom=" << m_currentMomOfBin << endl;
  }

  lowRange = log(lowRange/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  highRange = log(highRange/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  */
  double range[2];
  getRangeDeDx(&range[0],a_particleIndex);
  lowRange = log(range[0] / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  highRange = log(range[1] / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
  if(lowRange > highRange){
    cout << " Error: lowRange(" << lowRange << ") > highRange(" << highRange << ")" << endl;
    lowRange = log(0.9*m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
    highRange = log(1.1*m_bichselCurves[a_particleIndex]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
    cout << "      Forced to " << lowRange << "  and " << highRange << endl;
  }
  TH1D* zoomHisto = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);

  //zoomHisto->GetXaxis()->SetMoreLogLabels();
  double counts = HistogramUtilities::sumOfContent(zoomHisto, -1, true);
  //double counts = zoomHisto->Integral();

  //double integral = HistogramUtilities::integral(zoomHisto,-1);
  //double mu = zoomHisto->GetMean();
  //double sig = zoomHisto->GetStdDev();

  double gausData[3];
  HistogramUtilities::getGausData(gausData, zoomHisto, true, false, lowRange, highRange);
  double integral = gausData[0];
  double mu = gausData[1];
  double sig = gausData[2];

  vector<double> returnParams;
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);

  if(integral == 0 || std::isnan(mu) || std::isnan(sig)) return returnParams;

  //double mu = 0.0;
  //double sig = 0.0;
  //PhysMath::gausToLogNorm(mean,stdDev,mu,sig);

    TString imageName = TString::Format("%s_Rap_%02d_mTm0_%02d_Zoomed",
        m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin,m_currentMtM0Bin);

    zoomHisto->SetName(imageName.Data());
    zoomHisto->SetTitle(Form("TPC Particles for y=[%1.2f,%1.2f] m_{T}-m{0}=[%2.3f,%2.3f] (GeV/c) p=%2.3f : zoom on %s",
      m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
      m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
      m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
      m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
      m_currentMomOfBin,
      m_partInfo->GetParticleName(a_particleIndex,0).Data()
      ));
    zoomHisto->GetXaxis()->SetTitle("Z_{TPC}");
    zoomHisto->GetYaxis()->SetTitle("dN/dZ_{TPC}");
    zoomHisto->SetStats(kTRUE);
  //cout << "mean: " << mean << "  stdDev: " << stdDev << "  mu: " << mu << "   sig: " << sig << endl;
  //fitFunction->SetParameter(0,(counts)*zoomHistoTemp->GetBinWidth(2));
    m_gausFunct->SetRange(lowRange,highRange);
    m_gausFunct->SetParameter(0,counts);
    m_gausFunct->SetParLimits(0,0.1*counts, 3.0*counts);
    m_gausFunct->SetParameter(1,mu);
    m_gausFunct->SetParLimits(1, mu - 0.3, mu + 0.4);
    m_gausFunct->SetParameter(2,sig);
    m_gausFunct->SetParLimits(2,0.001*sig,10.0*sig);
    m_gausFunct->SetParNames("N","#mu","#sigma");
    TF1* initialFunct = (TF1*) m_gausFunct->Clone();
    initialFunct->SetLineColor(kViolet);
    zoomHisto->SetStats(kTRUE);
    gStyle->SetOptFit(0111);
  //m_logNormFunct->Draw("same");


  //FITTING HERE

    cout << "   starting zoom fit" << endl;
  string fitOptions = "EQ";
  #ifdef _ZFITTER_DEBUG_
    fitOptions = "EV";
  #endif
  if(m_doIntegralOfBinForFits) fitOptions += "I";

  m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(zoomHisto,-1,true);
  int status = zoomHisto->Fit("m_gausFunct",fitOptions.c_str(),"");//,ZTPCPredict-0.5, ZTPCPredict+0.5);//,-5.0, 5.0);
  returnParams[6] = status;
  if(status == 0){
    TF1* fitResult = zoomHisto->GetFunction("m_gausFunct");
    if(fitResult->GetParameter(0) != 0.0){
      returnParams[0] = fitResult->GetParameter(0);
      returnParams[1] = fitResult->GetParameter(1);
      returnParams[2] = fitResult->GetParameter(2);
      returnParams[3] = fitResult->GetParError(0);
      returnParams[4] = fitResult->GetParError(1);
      returnParams[5] = fitResult->GetParError(2);
    }
    //draw result
    fittingCanvas->cd();
    //gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits) initialFunct->Draw("same");
    fitResult->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();
    cout << "     ";
    string pmNames[3] = {"Minus","","Plus"};
    if(m_inAllCentMode){
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/Rap_%02d_mTm0_%02d_%s_Zoomed.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str()));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/Rap_%02d_mTm0_%02d_%s_Zoomed.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentCentIndex, m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str() ));
    }

    //fittingCanvas->SaveAs(Form("%s/%s/%s/ZTPC_Zoomed/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
    //      m_partInfo->GetParticleName(m_currentPartIndex,0).Data(),imageName.Data()));
  }else{
    cout << "Bad Fit for zoom " << m_partInfo->GetParticleName(a_particleIndex,0) << " Rap Bin " << m_currentRapBin << " mTm0 Bin " << m_currentMtM0Bin << endl;
    fittingCanvas->cd();
    //gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits) initialFunct->Draw("same");
    m_gausFunct->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();
    //fittingCanvas->SaveAs(Form("%s/%s/%s/ZTPC_Zoomed/BAD_%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
    //       m_partInfo->GetParticleName(m_currentPartIndex,0).Data(),imageName.Data()));
    string pmNames[3] = {"Minus","","Plus"};
    if(m_inAllCentMode){
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/Rap_%02d_mTm0_%02d_%s_Zoomed________BAD.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str()));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/Rap_%02d_mTm0_%02d_%s_Zoomed________BAD.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentCentIndex, m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str() ));
    }
  }


  delete zoomHisto;
  delete initialFunct;
  delete fittingCanvas;

  return returnParams;

}


vector<double> ZFitter::fitIsolation_ZTPC_Pearson(){

  #ifdef _ZFITTER_DEBUG_
    cout << "Starting ZFitter::fitIsolation_ZTPC_Pearson() " << endl;
  #endif

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas_iso","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  double lowRange = m_currentHistoToFit->GetBinLowEdge(1);
  double highRange = m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX()+1);

  if(m_inAllCentMode){// || m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin) == 0.0){


    int maxBin = HistogramUtilities::getMaxBinInRange(m_currentHistoToFit, -0.068, 0.068);
    int lowBinReturn = -1;
    int highBinReturn = -1;
    HistogramUtilities::getPeakRange(m_currentHistoToFit, maxBin, lowBinReturn, highBinReturn,0.05);
    lowRange = m_currentHistoToFit->GetBinLowEdge(lowBinReturn);
    highRange = m_currentHistoToFit->GetBinLowEdge(highBinReturn+1);
    #ifdef _ZFITTER_DEBUG_
      cout << "    maxBin = " << maxBin  << "    lowBinReturn: " << lowBinReturn << "      highBinReturn: " << highBinReturn << endl;
      cout << "    Histogram Utilities predicts a range of [" << lowRange << ", " << highRange << "] from bin [" << lowBinReturn << ", " << highBinReturn << "]" << endl;
    #endif
    double mean_3p = m_currentHistoToFit->GetBinCenter(maxBin);
    /*double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(-0.068), 2);
    double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.0), 2);
    double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.068), 2);
    double N_3p, mean_3p, sig_3p = 0.0;
    PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,-0.068,y1,0.0,y2,0.068,y3);
    cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
    if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
      lowRange  = mean_3p  - 3.0*sig_3p;
      highRange  = mean_3p + 3.0*sig_3p;
      cout << "     After: " << lowRange << " to " << highRange << endl;
    }else{
      cout << "   3pt gaus failed, just not zooming in" << endl;
    }*/

    if(0.0 < lowRange || 0.0 > highRange){
      lowRange = -log(1.14);
      highRange = log(1.14); // 2 sigma guess
    }

    double countsInWindow = m_currentHistoToFit->Integral(m_currentHistoToFit->FindBin(lowRange), m_currentHistoToFit->FindBin(highRange));

    //if contaminated to the left, find inflection point
    if(m_particleLeftRightSeedDirection[m_currentPartIndex] > 0 && countsInWindow > 60){
      bool particleIsContaminatingFromLeft = false;
      double binCenterOfInflection = -1;
      int meanBin = m_currentHistoToFit->FindBin(mean_3p);
      for(int binX = m_currentHistoToFit->FindBin(mean_3p) - 5; binX > 0; binX--){
        if(m_currentHistoToFit->GetBinContent(binX) >  m_currentHistoToFit->GetBinContent(binX + 1)){
          particleIsContaminatingFromLeft = true;
          binCenterOfInflection = m_currentHistoToFit->GetBinCenter(binX + 1);
          break;
        }
      }
      if(particleIsContaminatingFromLeft){
        lowRange = mean_3p - (mean_3p - binCenterOfInflection)*0.8;
        #ifdef _ZFITTER_DEBUG_
          cout << "    Contaminant Inflection found at  " << mean_3p - (mean_3p - binCenterOfInflection)  << "    cutting off at low : " << lowRange <<  endl;
        #endif
      }
    }

    // right contaminant
    if(m_particleLeftRightSeedDirection[m_currentPartIndex] < 0 && countsInWindow > 60){
      bool particleIsContaminatingFromRight = false;
      double binCenterOfInflection = -1;
      for(int binX = m_currentHistoToFit->FindBin(mean_3p) + 5; binX <= m_currentHistoToFit->GetNbinsX() ; binX++){
        if(m_currentHistoToFit->GetBinContent(binX) >  m_currentHistoToFit->GetBinContent(binX - 1)){
          particleIsContaminatingFromRight = true;
          binCenterOfInflection = m_currentHistoToFit->GetBinCenter(binX - 1);
          break;
        }
      }
      if(particleIsContaminatingFromRight){
        highRange = mean_3p + (binCenterOfInflection - mean_3p)*0.8;
        #ifdef _ZFITTER_DEBUG_
          cout << "    Contaminant Inflection found at  " << mean_3p + (binCenterOfInflection - mean_3p)  << "    cutting off at high : " << highRange <<  endl;
        #endif
      }
    }


    m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
    m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
    m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
    m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);

    //countsInWindow = m_currentHistoToFit->Integral(m_currentHistoToFit->FindBin(lowRange), m_currentHistoToFit->FindBin(highRange));

    #ifdef _ZFITTER_DEBUG_
      cout << "There are " << countsInWindow << " counts in the range" << endl;
    #endif

  }else{
    lowRange = m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
    highRange = m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
    #ifdef _ZFITTER_DEBUG_
      cout << "Loading previous histogram range data..." << endl;
    #endif
  }

  #ifdef _ZFITTER_DEBUG_
    //cout << "There are " << countsInWindow << " counts in the range" << endl;
    cout << "Range set to [" << lowRange << ", " << highRange << "]" << endl;
    cout << "Number of Bins before truncation: " << m_currentHistoToFit->GetNbinsX() << endl;
    HistogramUtilities::printHistoInfo(m_currentHistoToFit);

  #endif



  TH1D* zoomHisto = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
  lowRange = HistogramUtilities::getLowEdge(zoomHisto);
  highRange = HistogramUtilities::getHighEdge(zoomHisto);

  #ifdef _ZFITTER_DEBUG_
    cout << "Number of Bins after truncation: " << zoomHisto->GetNbinsX() << endl;
    HistogramUtilities::printHistoInfo(zoomHisto);

  #endif

  zoomHisto->SetMarkerStyle(20);
  zoomHisto->SetMarkerColor(kBlack);
  //double counts = HistogramUtilities::sumOfContent(zoomHisto, -1, true);

  //double integral = HistogramUtilities::integral(zoomHisto,-1);
  //double mu = zoomHisto->GetMean();
  //double sig = zoomHisto->GetStdDev();

  double gausData[3];
  HistogramUtilities::getGausData(gausData, zoomHisto, true, false, lowRange, highRange);
  double integral = gausData[0];
  double mu = gausData[1];
  double sig = gausData[2];

    #ifdef _ZFITTER_DEBUG_
      cout << "  Gaus Data :  Integral: " << integral << "  mu: " << mu << "  sig: " << sig  << endl;
    #endif




  vector<double> returnParams;
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);

  returnParams.push_back(-999); //errors
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);
  returnParams.push_back(-999);

  returnParams.push_back(-999); //status

  if(integral == 0 || std::isnan(mu) || std::isnan(sig) || sig < 0.0){

  #ifdef _ZFITTER_DEBUG_
    cout << "ERROR in Pearson Fit with gaus data... Integral: " << integral << "  mu: " << mu << "  sig: " << sig << "    so returning... " << endl;
  #endif

    return returnParams;
  }
  //double mu = 0.0;
  //double sig = 0.0;
  //PhysMath::gausToLogNorm(mean,stdDev,mu,sig);

  TString imageName = TString::Format("%s_Rap_%02d_mTm0_%02d_Zoomed",
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin,m_currentMtM0Bin);

  zoomHisto->SetName(imageName.Data());
  zoomHisto->SetTitle(Form("TPC Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] (GeV/c) p=%2.3f : zoom on %s : PearsonIV Fit",
    m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
    m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
    m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
    m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
    m_currentMomOfBin,
    m_partInfo->GetParticleName(m_currentPartIndex,0).Data()
    ));
  zoomHisto->GetXaxis()->SetTitle("Z_{TPC}");
  zoomHisto->GetYaxis()->SetTitle("dN/dZ_{TPC}");
  zoomHisto->SetStats(kTRUE);
  //cout << "mean: " << mean << "  stdDev: " << stdDev << "  mu: " << mu << "   sig: " << sig << endl;
  //fitFunction->SetParameter(0,(counts)*zoomHistoTemp->GetBinWidth(2));
  TF1* pearsonFunct = PhysMath::pearsonIVFunct(1.0, mu, sig, 60.0, 0.0);
  string functionName = Form("pearson_%02d_%02d_%s", m_currentRapBin, m_currentMtM0Bin,m_partInfo->GetParticleName(m_currentPartIndex, m_currentPartCharge).Data());
  if(m_inAllCentMode) functionName += "_allCent";
  else functionName += Form("_Cent%02d",m_currentCentIndex);
  pearsonFunct->SetName(functionName.c_str());
  pearsonFunct->SetRange(lowRange,highRange);


  if(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0){
    //the function has been fit before

    if(m_allowSpecCentFitsToVary){
      pearsonFunct->SetParameter(1,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->SetParameter(2,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->SetParameter(3,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->SetParameter(4,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));

      pearsonFunct->SetParLimits(1, pearsonFunct->GetParameter(1) - 0.3, pearsonFunct->GetParameter(1) + 0.4);
      pearsonFunct->SetParLimits(2,0.0,100.0*pearsonFunct->GetParameter(2));
      pearsonFunct->SetParLimits(3,0.0,500);
      pearsonFunct->SetParLimits(4,-60.0,60.0);

      //50.0*tanh(a_nu);

    }else{
      pearsonFunct->SetParLimits(1,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin),
                                   m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->SetParLimits(2,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin),
                                   m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->SetParLimits(3,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin),
                                   m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->SetParLimits(4,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin),
                                   m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));



      pearsonFunct->FixParameter(1,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->FixParameter(2,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->FixParameter(3,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      pearsonFunct->FixParameter(4,m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));

      #ifdef _ZFITTER_DEBUG_
        cout << "    Fixing Parameters " << endl;
        for(int jjj = 1; jjj < 5; jjj++){
          cout << "        Param " << jjj << " = " << pearsonFunct->GetParameter(jjj)  << endl;
        }
      #endif
    }


  }else{
    //the function has NOT been fit before
    pearsonFunct->SetParLimits(1, mu - 0.3, mu + 0.4);
    pearsonFunct->SetParLimits(2,0.0,100.0*sig);
    pearsonFunct->SetParLimits(3,0.0,500);
    pearsonFunct->SetParLimits(4,-60.0,60.0);
  }
  pearsonFunct->SetParLimits(0,0.0, 3.0);

  pearsonFunct->SetParNames("N","#lambda (location)","#alpha (spread)","m (tail)","#nu (skew)");
  TF1* initialFunct = (TF1*) pearsonFunct->Clone();
  string initialFunctName = pearsonFunct->GetName();
  initialFunctName += "_initFunct";
  initialFunct->SetName(initialFunctName.c_str());
  initialFunct->SetLineColor(kYellow);
  zoomHisto->SetStats(kTRUE);
  gStyle->SetOptFit(0111);

  double integralOfHistogram = HistogramUtilities::sumOfContent(zoomHisto,-1,true);
  #ifdef _ZFITTER_DEBUG_
    cout << "  Integral to scale to make histo amp = 1 :  " << integralOfHistogram << endl;
  #endif
  // root doesn't like big parameters, so scale the histogram to make the amplitude reasonable
  for(int bin = 1; bin <= zoomHisto->GetNbinsX(); bin++){
    double newContent = zoomHisto->GetBinContent(bin)/integralOfHistogram;
    double newError = zoomHisto->GetBinError(bin)/integralOfHistogram;
    zoomHisto->SetBinContent(bin, newContent);
    zoomHisto->SetBinError(bin, newError);
  }


  //FITTING HERE
  #ifdef _ZFITTER_DEBUG_
    string fitOptionsPearsonTPC = "ESV";  // Put in B? Matt 2022 Jan 17
    cout << "  Integral of Fit region: " <<  HistogramUtilities::sumOfContent(zoomHisto,-1,true) << endl;
    cout << " Starting the PearsonIV fit..." << endl;
  #else
    string fitOptionsPearsonTPC = "ESQ";
  #endif

  if(m_doIntegralOfBinForFits) fitOptionsPearsonTPC += "I";

  TFitResultPtr fitResultPtr = zoomHisto->Fit(pearsonFunct->GetName(),fitOptionsPearsonTPC.c_str(),"",lowRange,highRange);
  returnParams[10] = fitResultPtr->Status();
  int status = returnParams[10];
  #ifdef _ZFITTER_DEBUG_
    cout << " Ended with status : " << status << endl;
  #endif
  if(status == 0){
    TF1* fitResult = zoomHisto->GetFunction(pearsonFunct->GetName());
    if(fitResult->GetParameter(0) > 0.0){
      returnParams[0] = fitResult->GetParameter(0)*integralOfHistogram;
      returnParams[1] = fitResult->GetParameter(1);
      returnParams[2] = fitResult->GetParameter(2);
      returnParams[3] = fitResult->GetParameter(3);
      returnParams[4] = fitResult->GetParameter(4);

      returnParams[5] = fitResult->GetParError(0)*integralOfHistogram;
      returnParams[6] = fitResult->GetParError(1);
      returnParams[7] = fitResult->GetParError(2);
      returnParams[8] = fitResult->GetParError(3);
      returnParams[9] = fitResult->GetParError(4);

      //correct the integral and error parameters
      Int_t numberOfSamplePoints = 1000;
      double* x_array= new double[numberOfSamplePoints];
      double* w_array= new double[numberOfSamplePoints];
      fitResult->CalcGaussLegendreSamplingPoints(numberOfSamplePoints,x_array,w_array,1e-15);
      returnParams[0] = fitResult->IntegralFast(numberOfSamplePoints,x_array,w_array,-0.5,0.5)*integralOfHistogram;
      returnParams[5] = fitResult->IntegralError(-0.5,0.5,fitResultPtr->GetParams(), fitResultPtr->GetCovarianceMatrix().GetMatrixArray() )*integralOfHistogram;

      #ifdef _ZFITTER_DEBUG_
        cout << "   N from fit =  " << fitResult->GetParameter(0)*integralOfHistogram << " +/- " << fitResult->GetParError(0)*integralOfHistogram << endl;
        cout << "   N from int =  " << returnParams[0] << " +/- " << returnParams[5] << endl;
      #endif

    }
    //draw result
    fittingCanvas->cd();
    //gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits) initialFunct->Draw("same");
    fitResult->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();
    cout << "     ";
    string pmNames[3] = {"Minus","","Plus"};
    if(m_inAllCentMode){
      gPad->SetLogy(true);
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/RapBin_%02d_mTm0Bin_%02d_%s_Zoomed.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str()));
      if(m_saveNoLogImages){
        gPad->SetLogy(false);
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/NoLog_RapBin_%02d_mTm0Bin_%02d_%s_Zoomed.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str()));
      }
    }else{
      gPad->SetLogy(true);
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_%s_Zoomed.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentCentIndex, m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str() ));
      if(m_saveNoLogImages){
        gPad->SetLogy(false);
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_%s_Zoomed.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentCentIndex, m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str() ));
      }
    }
  }else{   // status != 0
    cout << "Bad Fit for zoom " << m_partInfo->GetParticleName(m_currentPartIndex,0) << " Rap Bin " << m_currentRapBin << " mTm0 Bin " << m_currentMtM0Bin << endl;
    returnParams[0] = pearsonFunct->GetParameter(0)*integralOfHistogram;
    returnParams[1] = pearsonFunct->GetParameter(1);
    returnParams[2] = pearsonFunct->GetParameter(2);
    returnParams[3] = pearsonFunct->GetParameter(3);
    returnParams[4] = pearsonFunct->GetParameter(4);

    returnParams[5] = pearsonFunct->GetParError(0)*integralOfHistogram;
    returnParams[6] = pearsonFunct->GetParError(1);
    returnParams[7] = pearsonFunct->GetParError(2);
    returnParams[8] = pearsonFunct->GetParError(3);
    returnParams[9] = pearsonFunct->GetParError(4);

    fittingCanvas->cd();
    gPad->SetRightMargin(0.45);
    zoomHisto->Draw("E");
    if(m_drawInitialSeedsToFits) initialFunct->Draw("same");
    //m_gausFunct->Draw("same");
    gPad->Update();
    gSystem->ProcessEvents();
    string pmNames[3] = {"Minus","","Plus"};
    if(m_inAllCentMode){
      gPad->SetLogy(true);
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/RapBin_%02d_mTm0Bin_%02d_%s_Zoomed__________BAD.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str()));
      if(m_saveNoLogImages){
      gPad->SetLogy(false);
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/NoLog_RapBin_%02d_mTm0Bin_%02d_%s_Zoomed__________BAD.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str()));
      }
    }else{
      gPad->SetLogy(true);
      fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_%s_Zoomed__________BAD.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentCentIndex, m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str() ));
      if(m_saveNoLogImages){
        gPad->SetLogy(false);
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_%s_Zoomed__________BAD.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
          m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentCentIndex, m_currentRapBin, m_currentMtM0Bin, pmNames[1+m_currentPartCharge].c_str() ));
      }
    }
  }




  delete zoomHisto;
  delete initialFunct;
  delete pearsonFunct;

  delete fittingCanvas;

  return returnParams;
}





//==================================================================================================================================
//==================================================================================================================================

vector<double> ZFitter::fitIsolation_ZTPC_withBackground(vector<double> a_fixedParameters){
  cout << "   in  ZFitter::fitIsolation_ZTPC_withBackground(vector<double> a_fixedParameters)" << endl;
  if(!m_currentHistoToFit){
    cout << "There is nothing to fit. The m_currentHistoToFit ptr = " << m_currentHistoToFit << endl;
  }
  double sigmaGuess = log(1.07); //7 percent resolution
  double fracOfWindowForBackground = 0.2; //Using this fraction of histogram range from left to seed the background parameters

  //Setup Return Structure
  vector<double> returnParams;
  returnParams.push_back(-999); //0 parameters
  returnParams.push_back(-999); //1
  returnParams.push_back(-999); //2
  returnParams.push_back(-999); //3
  returnParams.push_back(-999); //4
  returnParams.push_back(-999); //5
  returnParams.push_back(-999); //6
  returnParams.push_back(-999); //7
  returnParams.push_back(-999); //8 parameter errors
  returnParams.push_back(-999); //9
  returnParams.push_back(-999); //10
  returnParams.push_back(-999); //11
  returnParams.push_back(-999); //12
  returnParams.push_back(-999); //13
  returnParams.push_back(-999); //14
  returnParams.push_back(-999); //15
  returnParams.push_back(-999); //16 status
  returnParams.push_back(-999); //17 chi^2


  TCanvas* fittingCanvas = new TCanvas("fittingCanvas_iso_tof","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);


  TH1D* zoomHisto = m_currentHistoToFit;
  double lowRange = HistogramUtilities::getLowEdge(zoomHisto);
  double highRange = HistogramUtilities::getHighEdge(zoomHisto);

  //########   Create the Fit Function  #############
  TF1* fit_funct = fit_funct = new TF1("Gaus_withPowerLawFunction",
    "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]   +  [4]*exp(-[5]*(x - [6])) + [7]",lowRange,highRange);
  fit_funct->SetNpx(500);
  fit_funct->SetLineColor(kBlue);
  fit_funct->SetParNames("N","#mu","#sigma","ignore","A","b","x_{0}","y_{0}");
  fit_funct->FixParameter(3,0.0);







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
  }







  if(bAvg < 0.0 || std::isnan(bAvg) || std::isnan(aAvg)){ // retry background seed
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
      fit_funct->SetParLimits(4,0.0,1.0E3 * aAvg);
      fit_funct->SetParameter(5,bAvg); // Exponential factor
      fit_funct->SetParLimits(5,1E-4*bAvg,1E4*bAvg);
      fit_funct->FixParameter(6,0.0); // x_0
      fit_funct->FixParameter(7,0.0); // y_0
    }
  }else{ // good background seed
    fit_funct->SetParameter(4,aAvg); // Amp
    fit_funct->SetParLimits(4,0.0,1.0E3 * aAvg);
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
  //double stdDev = zoomHisto->GetStdDev();
  double stdDev = sigmaGuess;



  double stdDevCalc = 0.0;
  double newAmp = 0.0;
  double sigmaParCalc = 0.0;
  if(a_fixedParameters.size() > 4 && a_fixedParameters[2] > 0.0){
    cout << " Using Fixed Parameters to Seed Amplitude" << endl;

    newAmp = zoomHisto->GetBinContent(zoomHisto->FindBin(a_fixedParameters[1])) / gaus(a_fixedParameters[1], a_fixedParameters[2], a_fixedParameters[1]);
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
    if(newAmp < 0.0){
      cout << "    Seeding encountered issues: Using a single guess for the peak parameters: N mu sig" << endl;
      newAmp = zoomHisto->GetBinContent(zoomHisto->FindBin(0.0)) / gaus(0.0, log(1.07), 0.0);
      stdDevCalc = 1.0;
      sigmaParCalc = log(1.07);
    }

  }



  cout << " Starting Parameters:   N:" << newAmp << "  mu:" << mean << "   sig:" << sigmaParCalc << endl;
  cout << "                        A:" << aAvg << "  b:" << bAvg << "  x_0: 0.0   y_0:0.0" << endl;


  if(newAmp == 0 || std::isnan(mean) || std::isnan(sigmaParCalc)){
    cout << " Error in ZFitter::fitIsolation_ZTPC_withBackground: integral = "<< newAmp << " mean = "<< mean << "  stdDev = "<< sigmaParCalc << endl;
    //delete fittingCanvas;
    return returnParams;
  }

  TString imageName = TString::Format("%s_withBackground_rap_%02d_mTm0_%02d",
      m_partInfo->GetParticleName(m_currentPartIndex,0).Data(),m_currentRapBin,m_currentMtM0Bin);
  zoomHisto->SetName(imageName.Data());
  zoomHisto->SetTitle(Form("ZTPC Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] (GeV/c) p=%2.3f : zoom on %s",
    m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
    m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
    m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
    m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
    m_currentMomOfBin,
    m_partInfo->GetParticleName(m_currentPartIndex,0).Data()
    ));
  zoomHisto->GetXaxis()->SetTitle("Z_{TPC}");
  zoomHisto->GetYaxis()->SetTitle("dN/dZ_{TPC}");
  zoomHisto->SetStats(kTRUE);
  gStyle->SetOptFit(0111);




  //##########  Setup Initial Seed Function to Draw for Debugging  #############
  fit_funct->SetParameter(0,newAmp);
  fit_funct->SetParLimits(0,0.01*newAmp, 10.0*newAmp);
  fit_funct->SetParameter(1,mean);
  fit_funct->SetParLimits(1, mean - 0.15, mean + 0.15);
  fit_funct->SetParameter(2,sigmaParCalc);
  fit_funct->SetParLimits(2,0.1*sigmaParCalc,10.0*sigmaParCalc);

  TF1* initialFunct = NULL;
  initialFunct = (TF1*) fit_funct->Clone();
  initialFunct->SetLineColor(kYellow);

  // Create
  TF1* source_funct_init = NULL;

  source_funct_init = new TF1("gaus_Source_Init", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
  source_funct_init->SetNpx(200);
  source_funct_init->SetLineColor(kYellow+1);
  source_funct_init->SetLineStyle(2);
  source_funct_init->SetParNames("N","#mu","#sigma");
  source_funct_init->SetParameter(0,fit_funct->GetParameter(0));
  source_funct_init->SetParameter(1,fit_funct->GetParameter(1));
  source_funct_init->SetParameter(2,fit_funct->GetParameter(2));


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



  cout << "Printout of Funciton: " << endl;
  fit_funct->Print();
  cout << endl;
  //###############   FIT   THE   CURVE   #####################
  cout << "      Starting the fit in funct" << endl;
  string fitOptions = "EQ";
  #ifdef _ZFITTER_DEBUG_
    fitOptions = "EV";
  #endif
  if(m_doIntegralOfBinForFits) fitOptions += "I";
  m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(zoomHisto,-1,true);
  int status = zoomHisto->Fit("Gaus_withPowerLawFunction",fitOptions.c_str());
  cout << "      Finished the Fit" << endl;
  if(status == 0){
    TF1* fitResult = zoomHisto->GetFunction("Gaus_withPowerLawFunction");

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
      cout << " Error in ZFitter::fitIsolation_ZTPC_withBackground: fit status = "<< status << endl;
    }


    // Create Violet Source Peak For Drawinga
    TF1* source_funct = new TF1("gaus_Source", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
    source_funct->SetNpx(200);
    source_funct->SetLineColor(kViolet);
    source_funct->SetLineStyle(2);
    source_funct->SetParNames("N","#mu","#sigma");
    source_funct->SetParameter(0,returnParams[0]);
    source_funct->SetParameter(1,returnParams[1]);
    source_funct->SetParameter(2,returnParams[2]);



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
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),imageName.Data()));
      }
    }else{
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(), m_currentCentIndex,imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(), m_currentCentIndex,imageName.Data()));
      }
    }

    delete source_funct;
    delete background_funct;
  }else{  //////   BAD FIT / NOT OPTIMAL



    TF1* fitResult = zoomHisto->GetFunction("Gaus_withPowerLawFunction");
    cout << "Error in ZFitter::fitIsolation_ZTPC_withBackground: fit status: " << status << "   " << m_partInfo->GetParticleName(m_currentPartIndex,0)
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
      returnParams[17] = fitResult->GetChisquare() / fitResult->GetNDF();
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
      returnParams[17] = fitResult->GetChisquare() / fitResult->GetNDF();
    }

    // Create Violet Source Peak For Drawinga
    TF1* source_funct = new TF1("gaus_Source", "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",lowRange,highRange);
    source_funct->SetNpx(200);
    source_funct->SetLineColor(kViolet);
    source_funct->SetLineStyle(2);
    source_funct->SetParNames("N","#mu","#sigma");
    source_funct->SetParameter(0,returnParams[0]);
    source_funct->SetParameter(1,returnParams[1]);
    source_funct->SetParameter(2,returnParams[2]);



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
    if(returnParams[0] < 0.0 || std::isnan(returnParams[1]) || returnParams[2] <= 0.0
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
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(),imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,1).Data(),imageName.Data()));
      }
    }else{
      if(m_currentPartCharge == -1){
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
         m_partInfo->GetParticleName(m_currentPartIndex,-1).Data(), m_currentCentIndex,imageName.Data()));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
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


  cout << "  Finished fitIsolation_ZTPC_withBackground" << endl;


  return returnParams;

}







//============================================================================================================
//============================================================================================================

void ZFitter::extrapolateMeanAndWidthParams_TPC(int a_particleIndex){

  cout << " Starting extrapolateMeanAndWidthParams_TPC for part " << a_particleIndex << endl;
  TCanvas* fittingCanvas = new TCanvas("FittingCanvasMeanAndWidth");
  fittingCanvas->SetWindowSize(1300,800);
  fittingCanvas->SetCanvasSize(1200,700);
  //fittingCanvas->Divide(1,2);

  TF1* meanExtrap = new TF1("meanExtrap","[0]*x + [1]",0.005,5.0);
  TF1* stdDevExtrap = new TF1("stdDevExtrap","[0]*x + [1]",0.005,5.0);


  //Get Mean Histogram
  TH1D* mean_histo = NULL;
  TH1D* stdDev_histo = NULL;
  mean_histo = m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->ProjectionY(Form("mean_proj_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);
  stdDev_histo = m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->ProjectionY(Form("stdDev_proj_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);

  //Make Prediction/Fit Difference
  TH1D* mean_diff_histo = (TH1D*) mean_histo->Clone();
  mean_diff_histo->SetName("mean_diff_histo");
  mean_diff_histo->GetYaxis()->SetTitle("#mu - #mu_{predict}");
  mean_diff_histo->GetXaxis()->SetTitle("m_{T} - m_{0}");
  for(int binX = 1; binX <= mean_histo->GetNbinsX(); binX++){
    double center = mean_histo->GetBinCenter(binX);
    double content = mean_histo->GetBinContent(binX);
    double error =  mean_histo->GetBinError(binX);
    double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,center,m_partInfo->GetParticleMass(m_currentPartIndex));
    // ln(fit/base)-ln(predict/base)
    double newContent = content - log(m_bichselCurves[a_particleIndex]->Eval(momentum) / m_bichselCurves[m_currentPartIndex]->Eval(momentum) );
    double newError = error;
    mean_diff_histo->SetBinContent(binX,newContent);
    mean_diff_histo->SetBinError(binX,newError);
  }


  int minReqBins = 3;
  int numBinsWithContent = 0;
  int endBin = -1;
  double endMtM0 = -1;

  endBin = m_currentMtM0Bin;
  endMtM0 = mean_diff_histo->GetBinCenter(endBin - 1) + 0.001;

  int firstBinWithData = 999;
  for(int bin = endBin - 1; bin > 0 ; bin--){
    if(stdDev_histo->GetBinContent(bin) > 0.0){
      numBinsWithContent++;
      firstBinWithData = bin;
    }else if(numBinsWithContent < minReqBins){
      numBinsWithContent = -999;
    }
  }

  if(stdDev_histo->GetBinContent(endBin - 1) <= 0.0){
    cout << "Cannot extrapolate because the bin before this doesn't have sigma data" << endl;
  }

  if(numBinsWithContent < minReqBins){
    cout << "End bin : " << endBin << " by not enough bins to constrain: " << numBinsWithContent <<  " < " << minReqBins << endl;
    cout << " You may have m_dEdxSeparatedStops high compared to the m_partDeDx_MomStartExtrapolation_Matrix_Plus" << endl;
    delete mean_histo;
    delete stdDev_histo;
    delete mean_diff_histo;
    return;
    // Zero the rest of the fits?
  }

  cout << "fitting mean histo" << endl;
  //FIT THE MEAN HISTOGRAM
  double firstBin = 1;
  if(endBin - firstBinWithData > 4) firstBin = endBin - 4;
  if(endBin - firstBinWithData > 5) firstBin = endBin - 5;
  //if(endBin - firstBinWithData > 6) firstBin = endBin - 6;
  //if(endBin - firstBinWithData > 7) firstBin = endBin - 7;
  //if(endBin - firstBinWithData > 8) firstBin = endBin - 8;
  //if(endBin - firstBinWithData > 9) firstBin = endBin - 9;
  //if(endBin - firstBinWithData > 10) firstBin = endBin - 10;
  double x1 = mean_diff_histo->GetBinCenter(firstBin);
  double x2 = mean_diff_histo->GetBinCenter(endBin - 1);
  double y1 = mean_diff_histo->GetBinContent(firstBin);
  double y2 = mean_diff_histo->GetBinContent(endBin - 1);
  double slope = (y2-y1)/(x2-x1);
  double intercept = y2 - x2*slope;
  meanExtrap->SetParameter(0,slope);
  meanExtrap->SetParameter(1,intercept);
  if(slope > 0.0){
    meanExtrap->SetParLimits(0,-5.0*slope,10.0*slope);
  }else{
    meanExtrap->SetParLimits(0,10.0*slope,-5.0*slope);
  }

  if(intercept > 0.0){
    meanExtrap->SetParLimits(1,-5.0*intercept,10.0*intercept);
  }else{
    meanExtrap->SetParLimits(1,10.0*intercept,-5.0*intercept);
  }
  int mean_status = mean_diff_histo->Fit("meanExtrap","","",x1-0.0001,x2+0.0001);
  TF1* fitResult = mean_diff_histo->GetFunction("meanExtrap");

  //Push the results to the histograms
  if(fitResult){
    for(int binX = endBin; binX <= mean_histo->GetNbinsX(); binX++){
      double center = mean_histo->GetBinCenter(binX);
      double error =  mean_histo->GetBinError(binX);
      double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,center,m_partInfo->GetParticleMass(m_currentPartIndex));
      double newMean = fitResult->Eval(center) + log(m_bichselCurves[a_particleIndex]->Eval(momentum) / m_bichselCurves[m_currentPartIndex]->Eval(momentum) );
      m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,binX, newMean);
    }
  }else{
    cout << "Error: Mean Extrapolation Failed with Status: " << mean_status << "  : Exiting!" << endl;
    delete mean_histo;
    delete stdDev_histo;
    return;
  }



  //Make Std Dev Difference Histogram
  TH1D* stdDev_diff_histo = (TH1D*) stdDev_histo->Clone();
  stdDev_diff_histo->SetName("stdDev_diff_histo");
  stdDev_diff_histo->SetTitle(";m_{T} - m_{0}; TPC Resolution");
  for(int binX = 1; binX <= stdDev_histo->GetNbinsX(); binX++){
    double newContent = exp(stdDev_histo->GetBinContent(binX)) - 1.0;
    double newError =  exp(stdDev_histo->GetBinContent(binX)) * stdDev_histo->GetBinError(binX);
    stdDev_diff_histo->SetBinContent(binX,newContent);
    stdDev_diff_histo->SetBinError(binX,newError);
  }

  cout << " fit stdDev " << endl;
  //FIT THE STDDEV DIFF HISTOGRAM
  x1 = stdDev_diff_histo->GetBinCenter(firstBin);
  x2 = stdDev_diff_histo->GetBinCenter(endBin - 1);
  y1 = stdDev_diff_histo->GetBinContent(firstBin);
  y2 = stdDev_diff_histo->GetBinContent(endBin - 1);
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
  int std_status = stdDev_diff_histo->Fit("stdDevExtrap","","",x1-0.0001,x2+0.0001);
  TF1* fitResult_std = stdDev_diff_histo->GetFunction("stdDevExtrap");

  //Push the results to the histograms
  if(fitResult_std){
    for(int binX = endBin; binX <= stdDev_diff_histo->GetNbinsX(); binX++){
      double center = stdDev_diff_histo->GetBinCenter(binX);
      double fitSTDDev = log(fitResult_std->Eval(center) + 1.0);
      m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,binX, fitSTDDev);
    }
  }else{
    cout << "Error: StdDev Extrapolation Failed with Status: " << std_status << "  : Exiting!" << endl;
    delete mean_histo;
    delete stdDev_histo;
    delete stdDev_diff_histo;
    return;
  }

  // DRAW THE FITS
  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  //proj->GetXaxis()->SetRangeUser(0.0,x2+0.5);
  mean_diff_histo->SetStats(true);
  mean_diff_histo->Draw("E");
  //gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxParamExtrapolations/MeanExtrap_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
    m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));

  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  //proj->GetXaxis()->SetRangeUser(0.0,x2+0.5);
  stdDev_diff_histo->SetStats(true);
  stdDev_diff_histo->Draw("E");
  //gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  gPad->Update();
  //gPad->ProccessEvents();
  //if(fitResult) fitResult->Draw("same");
  fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxParamExtrapolations/StdDevExtrap_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
    m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));


  delete fittingCanvas;
  delete mean_diff_histo;
  delete mean_histo;
  delete stdDev_histo;
  delete stdDev_diff_histo;

  m_partHasBeenExtrapolatedCharge[a_particleIndex][0] = true;

  return;
}


//============================================================================================================
//============================================================================================================

void ZFitter::extrapolateMeanAndWidthParams_TPC_SecondOrder(){

  //This fits all necessary parameters with student ndf

  #ifdef _ZFITTER_DEBUG_
    cout << " Starting postFit_ZTPC_SecondOrder "  << endl;
  #endif

  TCanvas* fittingCanvas = new TCanvas("FittingCanvasMeanAndWidth");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);
  //fittingCanvas->Divide(1,2);


  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]) continue;
    if(m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) continue;

    double lowMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity,m_TPC_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex],m_partInfo->GetParticleMass(m_currentPartIndex));
    double highMtM0 = PhysMath::mTm0_from_y_p(m_currentRapidity,m_TPC_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex],m_partInfo->GetParticleMass(m_currentPartIndex));
    if(lowMtM0 < 0.0){
      cout << " the rapidity is too extreme for the extrapolation procedure you have given setting low mTm0 to 0" << endl;
      lowMtM0 = 0.0;
    }
    if(highMtM0 < 0.0){
      cout << " the rapidity is too extreme for the extrapolation procedure, becuase high mTm0 is negative... (" << highMtM0 << " ) so continuing on without extrapolation" << endl;
      continue;
    }
    if(m_TPC_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex] > 900){
      cout << " no extrapolation range so continuing " << endl;
      continue;
    }
    double midMtM0 = 0.5*(lowMtM0 + highMtM0);
    #ifdef _ZFITTER_DEBUG_
      cout << "  extrapolating particle " << m_partInfo->GetParticleName(partIndex) << " with mTm0 fit range from " << lowMtM0 << "  to  " << highMtM0 << endl;
      cout << "      momentum range was " << m_TPC_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex] << "  to  " << m_TPC_Extrapolation_Range_Matrix[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex] << endl;
    #endif


    //make the histograms to extrapolate
    TH1D* param1DHistos[4];
    TF1* extrapFunctions[4];
    double numberOfParams = 3;
    string distributionName = "Gaussian";


    int pmIndexOfHistoToGrab = 0;
    if(m_currentPartCharge < 0) pmIndexOfHistoToGrab = 1;

    string paramNames[4] = {"N","mean","width","ndf"};
    string paramSymbols[4] = {"N","#mu","#sigma","#nu"};
    string functionNames[4];
    if(m_inAllCentMode){
      for(int jjj = 0; jjj < numberOfParams; jjj++){
        param1DHistos[jjj] = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][jjj][pmIndexOfHistoToGrab]->ProjectionY(Form("param_proj_%s_%s_%s_Rap%02d",
           paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),m_currentRapBin),
           m_currentRapBin,m_currentRapBin);
        functionNames[jjj] = Form("TPC_param_parab_%s_%s_%s_Rap%02d",paramNames[jjj].c_str(),
          m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
          m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),m_currentRapBin);
        param1DHistos[jjj]->GetYaxis()->SetTitle(Form("%s %s",distributionName.c_str(),paramSymbols[jjj].c_str()));
      }
    }else{
      for(int jjj = 0; jjj < numberOfParams; jjj++){
        param1DHistos[jjj] = m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][jjj][m_currentCentIndex][pmIndexOfHistoToGrab]->ProjectionY(Form("param_proj_%s_%s_%s_Cent%02d_Rap%02d",
           paramNames[jjj].c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin),
           m_currentRapBin,m_currentRapBin);

        functionNames[jjj] = Form("TPC_param_parab_%s_%s_%s_Cent%02d_Rap%02d",paramNames[jjj].c_str(),
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
            m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][jjj][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX,valueToPush);
            m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][jjj][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX, 0.0);
          }else{
            m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][jjj][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX, valueToPush);
            m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][jjj][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX, 0.0);
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
        fittingCanvas->SaveAs(Form("%s/%s/%s/TPCParamExtrapolations/Param%02d_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
        m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),jjj,
        m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
        m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/%s/TPCParamExtrapolations/Param%02d_%s_in_%s_space_Cent%02d_Rap_%02d.png",m_imagePreDir.c_str(),
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
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX, fitResult->Eval(center) * param1DHistos[0]->GetBinWidth(binX));
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX, 0.0);

        }else{
          m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][0][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinContent(m_currentRapBin, binX, fitResult->Eval(center) * param1DHistos[0]->GetBinWidth(binX));
          m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][0][m_currentCentIndex][pmIndexOfHistoToGrab]->SetBinError(m_currentRapBin, binX,0.0);
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
      fittingCanvas->SaveAs(Form("%s/%s/%s/TPCParamExtrapolations/Param%02d_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),0,
      m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/%s/TPCParamExtrapolations/Param%02d_%s_in_%s_space_Cent%02d_Rap_%02d.png",m_imagePreDir.c_str(),
      m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),0,
      m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
      m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentCentIndex,m_currentRapBin));
    }





     m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex] = true;
  } // end particle loop



  return;
}















void ZFitter::extrapolateMeanAndWidthParams_FromSoloTPC_TPC(int a_particleIndex){
  cout << " Starting extrapolateMeanAndWidthParams_FromSolo_TPC for part " << a_particleIndex << endl;
  TCanvas* fittingCanvas = new TCanvas("FittingCanvasMeanAndWidth");
  fittingCanvas->SetWindowSize(1300,800);
  fittingCanvas->SetCanvasSize(1200,700);
  //fittingCanvas->Divide(1,2);

  TF1* meanExtrap = new TF1("meanExtrap","[0]*x + [1]",0.005,5.0);
  TF1* stdDevExtrap = new TF1("stdDevExtrap","[0]*x + [1]",0.005,5.0);


  //Get Mean Histogram
  TH1D* mean_histo = m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->ProjectionY(Form("mean_proj_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);
  mean_histo->SetTitle(Form("Difference from BTOF Isolated dEdx Means %s Rapdity Bin %02d; m_{T}-m_{0}; #mu_{free} - #mu_{BTOF Isolated}",m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));
  TH1D* stdDev_histo = m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->ProjectionY(Form("stdDev_proj_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);
  stdDev_histo->SetTitle(Form("Difference from BTOF Isolated dEdx StdDev %s Rapdity Bin %02d; m_{T}-m_{0}; #sigma_{free} - #sigma_{BTOF Isolated}",m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));

  TH1D* mean_histo_solo = m_Fit_Data_Solo_ZTPC[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->ProjectionY(Form("mean_proj_solo_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);
  TH1D* stdDev_histo_solo = m_Fit_Data_Solo_ZTPC[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->ProjectionY(Form("stdDev_proj_solo_%s_%s_Rap%02d",
       m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),m_currentRapBin),
       m_currentRapBin,m_currentRapBin);


  mean_histo->Print(); cout << "\n\n\n";
  stdDev_histo->Print(); cout << "\n\n\n";
  mean_histo_solo->Print(); cout << "\n\n\n";
  stdDev_histo_solo->Print(); cout << "\n\n\n";

  mean_histo->Add(mean_histo_solo, -1.0);
  stdDev_histo->Add(stdDev_histo_solo, -1.0);

  delete mean_histo_solo;
  delete stdDev_histo_solo;




  int minReqBins = 3;
  int numBinsWithContent = 0;
  int endBin = m_currentMtM0Bin;
  double endMtM0 = mean_histo->GetBinCenter(endBin - 1) + 0.001;

  if(endBin < 5){
    cout << "cant extrapolate with this low of mTm0 bin = " << endBin << endl;
    delete mean_histo;
    delete stdDev_histo;
    return;
  }

  /*int firstBinWithData = 999;
  for(int bin = endBin - 1; bin > 0 ; bin--){
    if(stdDev_histo->GetBinContent(bin) > 0.0){
      numBinsWithContent++;
      firstBinWithData = bin;
    }else if(numBinsWithContent < minReqBins){
      numBinsWithContent = -999;
    }
  }

  if(stdDev_histo->GetBinContent(endBin - 1) <= 0.0){
    cout << "Cannot extrapolate because the bin before this doesn't have sigma data" << endl;
  }

  if(numBinsWithContent < minReqBins){
    cout << "End bin : " << endBin << " by not enough bins to constrain: " << numBinsWithContent <<  " < " << minReqBins << endl;
    cout << " You may have m_dEdxSeparatedStops high compared to the m_partDeDx_MomStartExtrapolation_Matrix_Plus" << endl;
    delete mean_histo;
    delete stdDev_histo;
    return;
    // Zero the rest of the fits?
  }
  */


  cout << "fitting mean histo" << endl;
  //FIT THE MEAN HISTOGRAM
  double firstBin = endBin - 4;
  //if(endBin - firstBinWithData > 5) firstBin = endBin - 5;
  //if(endBin - firstBinWithData > 6) firstBin = endBin - 6;
  //if(endBin - firstBinWithData > 7) firstBin = endBin - 7;
  //if(endBin - firstBinWithData > 8) firstBin = endBin - 8;
  //if(endBin - firstBinWithData > 9) firstBin = endBin - 9;
  //if(endBin - firstBinWithData > 10) firstBin = endBin - 10;
  double x1 = mean_histo->GetBinCenter(firstBin);
  double x2 = mean_histo->GetBinCenter(endBin - 1);
  double y1 = mean_histo->GetBinContent(firstBin);
  double y2 = mean_histo->GetBinContent(endBin - 1);
  double slope = (y2-y1)/(x2-x1);
  double intercept = y2 - x2*slope;
  meanExtrap->SetParameter(0,slope);
  meanExtrap->SetParameter(1,intercept);
  if(slope > 0.0){
    meanExtrap->SetParLimits(0,-5.0*slope,10.0*slope);
  }else{
    meanExtrap->SetParLimits(0,10.0*slope,-5.0*slope);
  }

  if(intercept > 0.0){
    meanExtrap->SetParLimits(1,-5.0*intercept,10.0*intercept);
  }else{
    meanExtrap->SetParLimits(1,10.0*intercept,-5.0*intercept);
  }
  int mean_status = mean_histo->Fit("meanExtrap","","",x1-0.0001,x2+0.0001);
  TF1* fitResult = mean_histo->GetFunction("meanExtrap");

  //Push the results to the histograms
  if(fitResult){
    for(int binX = endBin; binX <= mean_histo->GetNbinsX(); binX++){
      double center = mean_histo->GetBinCenter(binX);
      //double error =  mean_histo->GetBinError(binX);
      //double momentum = PhysMath::pTot_from_y_mTm0(m_currentRapidity,center,m_partInfo->GetParticleMass(m_currentPartIndex));
      double newMean = fitResult->Eval(center) + m_Fit_Data_Solo_ZTPC[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,binX);
      m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,binX, newMean);
    }
  }else{
    cout << "Error: Mean Extrapolation Failed with Status: " << mean_status << "  : Exiting!" << endl;
    delete mean_histo;
    delete stdDev_histo;
    return;
  }




  //cout << " fit stdDev " << endl;
  //FIT THE STDDEV DIFF HISTOGRAM
  x1 = stdDev_histo->GetBinCenter(firstBin);
  x2 = stdDev_histo->GetBinCenter(endBin - 1);
  y1 = stdDev_histo->GetBinContent(firstBin);
  y2 = stdDev_histo->GetBinContent(endBin - 1);
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
  int std_status = stdDev_histo->Fit("stdDevExtrap","","",x1-0.0001,x2+0.0001);
  TF1* fitResult_std = stdDev_histo->GetFunction("stdDevExtrap");

  //Push the results to the histograms
  if(fitResult_std){
    for(int binX = endBin; binX <= stdDev_histo->GetNbinsX(); binX++){
      double center = stdDev_histo->GetBinCenter(binX);
      double newSig = fitResult_std->Eval(center) + m_Fit_Data_Solo_ZTPC[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,binX);
      m_Fit_Data_ZTPC[m_currentPartIndex][a_particleIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,binX, newSig);
    }
  }else{
    cout << "Error: StdDev Extrapolation Failed with Status: " << std_status << "  : Exiting!" << endl;
    delete mean_histo;
    delete stdDev_histo;
    return;
  }

  // DRAW THE FITS
  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  //proj->GetXaxis()->SetRangeUser(0.0,x2+0.5);
  mean_histo->SetStats(true);
  mean_histo->Draw("E");
  //gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxParamExtrapolations/SoloFit_MeanExtrap_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
    m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));

  fittingCanvas->cd();
  gStyle->SetOptFit(0111);
  //proj->GetXaxis()->SetRangeUser(0.0,x2+0.5);
  stdDev_histo->SetStats(true);
  stdDev_histo->Draw("E");
  //gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  gPad->Update();
  //gPad->ProccessEvents();
  //if(fitResult) fitResult->Draw("same");
  fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxParamExtrapolations/SoloFit_StdDevExtrap_%s_in_%s_space_Rap_%02d.png",m_imagePreDir.c_str(),
    m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(a_particleIndex,m_currentPartCharge).Data(),
    m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_currentRapBin));


  delete fittingCanvas;
  delete mean_histo;
  delete stdDev_histo;

  m_partHasBeenExtrapolated_Solo_Charge[a_particleIndex][m_currentPlusMinusIndex] = true;

  return;
}



//============================================================================================================
//============================================================================================================






void ZFitter::fitSoloTPC(int a_ZTPC_Rebin, bool a_saveImages){
  cout << "####################################################################################################" << endl;
  cout << "                        FITTING BTOF ISOLATED ZTPC   " << endl;
  cout << "####################################################################################################" << endl;
  //min counts for a projection
  double minContentToConsider = 30;
  // make the image directory
  //string imageDir = m_imagePreDir + "/" + m_imgDirName + "/ZTPC_Solo";

  //CREATE THE GENERAL TF1 for Fitting (Log Normal)

  TF1* fitFunction = new TF1("fitFunction",
    "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*( x - [1] )/([2]*[2])) /[2]",
    -10,10);

  if(m_useStudentTDistributionsForTPC){
    delete fitFunction;
    fitFunction = (TF1*) m_student_T_Funct->Clone();
    fitFunction->SetName("fitFunction");
  }

  TF1* fitFunctionTwo = new TF1("fitFunctionTwo",
    "0.3989422804 * [0]*TMath::Exp(-0.5*( x - [1] )*( x - [1] )/([2]*[2])) /([2] ) + 0.3989422804 * [3]*TMath::Exp(-0.5*( x - [4] )*( x - [4] )/([5]*[5])) /([5])"
    ,-10,10);

  TCanvas* fitCanvas = new TCanvas("fittingCanvas2");
  fitCanvas->SetWindowSize(1300,700);
  fitCanvas->SetCanvasSize(1200,600);


  //cout << " Checking Electron Params : " << m_dEdxFitDataMinus[0][3][2]->GetBinContent(14,2) << endl;


  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    m_currentPlusMinusIndex = pmIndex;
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){

      if(partIndex == 3 || partIndex == 8) continue;

      TH3D* current3DHisto = NULL;
      if(pmIndex == 0){
        current3DHisto = m_rap_mTm0_dEdx_Plus_allCentIsolated[m_currentPartIndex][partIndex];
        m_currentPartCharge = 1;
      }else{
        current3DHisto = m_rap_mTm0_dEdx_Minus_allCentIsolated[m_currentPartIndex][partIndex];
        m_currentPartCharge = -1;
      }
      int numHorBins = current3DHisto->GetNbinsX();
      int numVertBins = current3DHisto->GetNbinsY();

      if(m_cuts->getEventConfig() == "ColliderCenter" || m_cuts->getEventConfig() == "FixedTarget"){
          if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][partIndex])
              || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]) ){
            #ifdef _ZFITTER_DEBUG_
             cout << "Continuing becasue of fit matrix toggle that is off for " << m_partInfo->GetParticleName(partIndex) << " in " << m_partInfo->GetParticleName(m_currentPartIndex) << " space " << endl;
            #endif
            continue;
        }
      }


      for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
        float currentRapVal = current3DHisto->GetXaxis()->GetBinCenter(HorBinIndex);
        //if(m_testing && ( HorBinIndex%(numHorBins/m_divTesting) != 0 ? currentRapVal != m_cuts->getYCM() : 0 )) continue;

        cout << "  attempting HorBinIndex: " << HorBinIndex <<endl;
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
          cout  << " rap vect has size " << endl;
          doThisRapBin = false;
          for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
            if(HorBinIndex == m_rapBinsToTest[iii]){
              doThisRapBin = true;
              cout << " rapidity found" << endl;
            }
          }
          if(!doThisRapBin){
            continue;
          }
        }

        for(int binVertIndex = 1; binVertIndex <= numVertBins; binVertIndex++){
          if(m_testing && binVertIndex%m_divTestingMtM0!=0) continue;
          cout << " attempting binVertIndex: " << binVertIndex << endl;
          m_currentRapBin = HorBinIndex;
          m_currentMtM0Bin = binVertIndex;
          m_currentRapidity = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
          m_currentMtM0 = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
          m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));


          if(m_cuts->getEventConfig() == "ColliderCenter" || m_cuts->getEventConfig() == "FixedTarget"){
              if((m_currentPartCharge == 1 &&  m_currentMomOfBin > m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex])
                  || (m_currentPartCharge == -1 &&  m_currentMomOfBin > m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex]) ){
                #ifdef _ZFITTER_DEBUG_
                  cout << "Continuing becasue of m_partDeDxMomStopMatrix " << m_partInfo->GetParticleName(partIndex) << " cuttoff < momentum "
                       << m_currentMomOfBin << " < " << m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex] << " or " << m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex] << endl;
                #endif
                continue;
            }
          }

          cout << "\n\n\n###########################################################################" << endl;
          cout << "Working on ZTPC (BTOF Isolated) " << m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data() << " in " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data()
                << "  space : "  << " RapBin: " << HorBinIndex << " MtM0Bin: " << binVertIndex << " Rap: " << m_currentRapidity << " MtM0: " << m_currentMtM0 << " p_tot: " << m_currentMomOfBin << endl;


          TH1D* proj = current3DHisto->ProjectionZ(Form("dEdxSoloProj_%02d_%01d_%01d_%02d_%02d",m_currentPartIndex, partIndex, pmIndex, HorBinIndex , binVertIndex),
                        HorBinIndex, HorBinIndex, binVertIndex, binVertIndex);
          proj->Rebin(a_ZTPC_Rebin);
          //TH1D* proj = HistogramUtilities::rebinMinCount(orig_proj, a_ZTPC_Rebin);
          HistogramUtilities::makeSqrtNErrors(proj);

          if(HistogramUtilities::sumOfContent(proj, -1) <= 1 ){  // && !m_testing)  || (counts < 20)){
            delete proj;
            cout << "Continuing: Histogram Doesn't have Counts (before truncation)" << endl;
            continue;
          }

          double lowRange = 0;
          double highRange = 0;
          lowRange = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) - 0.8;
          highRange = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin)) + 0.8;
          double lowRangeData = proj->GetBinLowEdge(HistogramUtilities::firstBinWithContent(proj));
          double highRangeData = proj->GetBinLowEdge(HistogramUtilities::lastBinWithContent(proj)+1);

          //double integral = proj->Integral();
          double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double stdDev = 0.06;

          double x_low = mean - stdDev;
          //double y_low = proj->GetBinContent(proj->FindBin(x_low));
          double x_mid = mean;
          //double y_mid = proj->GetBinContent(proj->FindBin(x_mid));
          double x_high = mean + stdDev;
          //double y_high = proj->GetBinContent(proj->FindBin(x_high));
          double y_low  = HistogramUtilities::getBinContentAvg(proj, proj->FindBin(x_low), 2);
          double y_mid  = HistogramUtilities::getBinContentAvg(proj, proj->FindBin(x_mid), 2);
          double y_high  = HistogramUtilities::getBinContentAvg(proj, proj->FindBin(x_high), 2);
          double N_init, mu_init, sig_init = 0;
          PhysMath::get3PointGausParams(N_init, mu_init, sig_init, x_low , y_low , x_mid , y_mid , x_high, y_high);

          if(sig_init > 0){
            lowRangeData  = mu_init  - 3.0*sig_init;
            highRangeData  = mu_init + 3.0*sig_init;
            mean = mu_init;
            stdDev = sig_init;
          }else{
            cout << "ERROR: 3 point gaussian failed in btof isolated dEdx fitting for zooming in and seeding fit params" << endl;
          }




          m_currentHistoToFit = HistogramUtilities::truncateHistogram(proj, lowRangeData, highRangeData);
          delete proj;
          double counts = HistogramUtilities::sumOfContent(m_currentHistoToFit, -1);

          if( minContentToConsider > counts ){  // && !m_testing)  || (counts < 20)){
            delete m_currentHistoToFit;
            cout << "Continuing: too few counts for solo fit partIndex:" << partIndex << " counts:" << counts << " < "  << minContentToConsider <<  endl;
            continue;
          }

          //double mean = m_currentHistoToFit->GetMean();
          //double stdDev = m_currentHistoToFit->GetStdDev();
          //double mean = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          //double stdDev = 0.06; //log(1.07) 7% width
          m_currentHistoToFit->Scale(1.0,"width");

          double truncatedTotalCounts = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
          double possibleNewBins = truncatedTotalCounts/6.0;
          if(possibleNewBins > 40.0) possibleNewBins = 40.0;
          TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(m_currentHistoToFit,(int) possibleNewBins,true,true);
          delete m_currentHistoToFit;
          m_currentHistoToFit = rebinHisto;



          //double integral = m_currentHistoToFit->Integral();

          //double ZTPCPredict = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));


          if(mean < m_currentHistoToFit->GetBinLowEdge(1)
           || mean > m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1)){
            cout << "Continuing: ZTPC Prediction is not within histogram range" << endl;
            delete m_currentHistoToFit;
            continue;
          }

          TString imageName = TString::Format("%s_RapBin_%02d_Mtm0Bin_%02d",
            m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),
            HorBinIndex, binVertIndex);
          m_currentHistoToFit->SetName(imageName.Data());
          m_currentHistoToFit->SetTitle(Form("bTOF Isolated %s for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f]; Z_{TPC} of %s ; dN/dZ_{TPC}",
            m_partInfo->GetParticleSymbol(partIndex,m_currentPartCharge).Data(),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
            m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data()
            ));
          m_currentHistoToFit->SetStats(kTRUE);
          m_currentHistoToFit->SetMarkerStyle(20);
          m_currentHistoToFit->SetMarkerColor(kBlack);

            //cout << "mean: " << mean << "  stdDev: " << stdDev << "  mu: " << mu << "   sig: " << sig << endl;
            //fitFunction->SetParameter(0,(counts)*m_currentHistoToFitTemp->GetBinWidth(2));
          //if(sig_init > 0){
          //  fitFunction->SetParameter(0,counts);
          //  fitFunction->SetParLimits(0,0.1*counts, 10.0*counts);
          //}else{
          if(m_useStudentTDistributionsForTPC){
            double ndof_seed = 50;
            double sig = PhysMath::studentT_var_from_stdDev(stdDev, ndof_seed);
            fitFunction->SetParameter(0,counts);
            fitFunction->SetParLimits(0,0.1*counts, 10.0*counts);
            fitFunction->SetParameter(1, mean);
            fitFunction->SetParLimits(1, mean - 0.15, mean + 0.35);
            fitFunction->SetParameter(2,sig);
            fitFunction->SetParLimits(2,0.05*sig,3.0*sig);
            fitFunction->SetParameter(3,ndof_seed);
            fitFunction->SetParLimits(3,8,250);
            fitFunction->SetParNames("N","#mu","#sigma","#nu");
          }else{
            fitFunction->SetParameter(0,counts);
            fitFunction->SetParLimits(0,0.1*counts, 10.0*counts);
            fitFunction->SetParameter(1, mean);
            fitFunction->SetParLimits(1, mean - 0.15, mean + 0.35);
            fitFunction->SetParameter(2,stdDev);
            fitFunction->SetParLimits(2,0.05*stdDev,3.0*stdDev);
            fitFunction->SetParNames("N","#mu","#sigma");
          }



           //draw the fits
           //TLine* predictLine = new TLine(ZTPCPredict,ZTPCPredict,0,100000);
           //predictLine->SetVertical();
           //predictLine->SetLineColor(kRed);

          if(a_saveImages){
            fitCanvas->cd();

            //gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetRightMargin(0.4);
            m_currentHistoToFit->SetStats(kTRUE);
            gStyle->SetOptFit(0111);
              //m_currentHistoToFit->Draw("E");
              //m_currentHistoToFit->SetAxisRange(-2.0, 2.0);
              //m_currentHistoToFit->GetXaxis()->SetRangeUser(-2.0, 2.0);
            m_currentHistoToFit->Draw("E");
            TPaveStats *statsBox = (TPaveStats*)fitFunction->FindObject("stats");
            if(statsBox){
              statsBox->SetX1NDC(.15);
              statsBox->SetY1NDC(.75);
              statsBox->SetX2NDC(.4);
              statsBox->SetY2NDC(.9);
            }
          }

          string fitVerbose = "EQ";
          #ifdef _ZFITTER_DEBUG_
            fitVerbose = "EV";
          #endif
          if(m_doIntegralOfBinForFits) fitVerbose += "I";
          m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
          //FITTING HERE
          int status = m_currentHistoToFit->Fit("fitFunction",fitVerbose.c_str(),"I",mean - 3.0*stdDev,mean + 3.0*stdDev);//,ZTPCPredict-0.5, ZTPCPredict+0.5);//,-5.0, 5.0);
          if(status == 0){
            TF1* fitResult = m_currentHistoToFit->GetFunction("fitFunction");
            if(fitResult->GetParameter(2) != 0.0){
                m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(0));
                m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(1));
                m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(2));
                m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(0));
                m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(1));
                m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(2));

                m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(0));
                m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(1));
                m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(2));
                m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(0));
                m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(1));
                m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(2));

                if(m_useStudentTDistributionsForTPC){
                  m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(3));
                  m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(3));
                  m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinContent(HorBinIndex,binVertIndex,fitResult->GetParameter(3));
                  m_Fit_Data_Solo_ZTPC[m_currentPartIndex][partIndex][3][m_currentPlusMinusIndex]->SetBinError(HorBinIndex,binVertIndex,fitResult->GetParError(3));
                }

                double differenceShift = fitResult->GetParameter(1) - log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
                m_Fit_Data_Solo_ZTPCShifts[m_currentPartIndex][partIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,differenceShift);
                m_Fit_Data_Solo_ZTPCShifts[m_currentPartIndex][partIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResult->GetParError(1));

            }
            //draw result
            if(a_saveImages){
              fitCanvas->cd();
              fitFunction->SetNpx(300);
              fitFunction->Draw("same");
              DrawPredictionLines_dEdx(fitCanvas);
              //gPad->SetLogx();
              gPad->SetLogy();
              gPad->Update();
              gSystem->ProcessEvents();
              fitCanvas->SaveAs(Form("%s/%s/%s/ZTPC_BTOFIsolated/%s/%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
                m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),imageName.Data()));
              gPad->SetLogy(false);
              gPad->Update();
              gSystem->ProcessEvents();
              fitCanvas->SaveAs(Form("%s/%s/%s/ZTPC_BTOFIsolated/%s/NoLog_%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),
                m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),m_partInfo->GetParticleName(partIndex,m_currentPartCharge).Data(),imageName.Data()));

            }
          }else{
            #ifdef _ZFITTER_DEBUG_
              cout << "BTOF Isolated ZTPC Fitting Failed with status: " << status << " for " << m_partInfo->GetParticleName(partIndex) << " in " << m_partInfo->GetParticleName(m_currentPartIndex) << endl;
            #endif
          }

          delete m_currentHistoToFit;
          //delete predictLine;
          //delete fitFunction;
        } // end vert bin
      }// end hor bin
    } //end part
  } // end plus minus loop




  //####################
  cout << "NOT SOLO FITTING ELECTRONS" << endl;
  cout << "OKAY IF YOU RAN ELECTRON STUDY" << endl;
  return;
  //####################




  cout << "Fitting Electrons" << endl;
  //Fit Electrons and Pions
  if(m_currentPartIndex < 2){ // either pion of kaon
    TH3D* current3DHisto = m_rap_mTm0_dEdx_allCent_PionElectron[m_currentPartIndex];
    //TH3D* current3DHisto = m_rap_mTm0_dEdx_allCentIsolated[3];
    int numHorBins = current3DHisto->GetNbinsX();
    int numVertBins = current3DHisto->GetNbinsY();
    for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
      float currentRapVal = current3DHisto->GetXaxis()->GetBinCenter(HorBinIndex);
        bool doThisRapBin = false;
        if(m_testing){
          if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(currentRapVal - m_cuts->getYCM()) < 0.01 )  ){
            doThisRapBin = false;
          }else{
            doThisRapBin = true;
          }

          for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
            if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
          }

        }
        if(m_testing && !doThisRapBin) continue;
      for(int binVertIndex = 1; binVertIndex <= numVertBins; binVertIndex++){
        m_currentRapBin = HorBinIndex;
        m_currentMtM0Bin = binVertIndex;
        m_currentRapidity = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
        m_currentMtM0 = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));


        // NEWER ISOLATED ELECTRON FITTING
        cout << "\n\nElectrons - HorBin: " << HorBinIndex << "  VertBin: " << binVertIndex << " ";
        TH1D* proj = current3DHisto->ProjectionZ(Form("dEdxElecPionProj_%02d_%02d", m_currentRapBin , m_currentMtM0Bin),
          m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
        double counts = HistogramUtilities::sumOfContent(proj, -1);
        HistogramUtilities::makeSqrtNErrors(proj);

      if( 200  > counts ){  // && !m_testing)  || (counts < 20)){
        delete proj;
        continue;
      }
      HistogramUtilities::makeSqrtNErrors(proj);
      //proj->Sumw2();
      proj->Scale(1.0,"width");
      m_currentHistoToFit = proj;
      double integral = m_currentHistoToFit->Integral();
      double mean = m_currentHistoToFit->GetMean();
      double stdDev = m_currentHistoToFit->GetStdDev();
      double mu_pi = log(m_bichselCurves[0]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double mu_e = log(m_bichselCurves[3]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
      double sig_pi = 0.07*mu_pi;
      double sig_e = 0.07*mu_e;

      cout << "  pi: " << mu_pi << "  "  << sig_pi << "  e: " << mu_e << "  " << sig_e << endl;
      if(m_currentMomOfBin < m_isolatedZTPC_starts[3]) continue;
      TString imageName = TString::Format("PionElectron_RapBin_%02d_mTm0Bin_%02d",m_currentRapBin, m_currentMtM0Bin);
      m_currentHistoToFit->SetName(imageName.Data());
      m_currentHistoToFit->SetTitle(imageName.Data());

      //if pion fit data exists, use that instead
      if( (m_currentPartCharge ==1 && m_Fit_Data_ZTPC[m_currentPartIndex][0][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0 )
        || (m_currentPartCharge ==-1 && m_Fit_Data_ZTPC[m_currentPartIndex][0][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0 )){
        double predictShift = -999;
      if(m_currentPartCharge == 1){
        predictShift = m_Fit_Data_ZTPC[m_currentPartIndex][0][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) - mu_pi;
        sig_pi = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);

      }else{
        predictShift = m_Fit_Data_ZTPC[m_currentPartIndex][0][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) - mu_pi;
        sig_pi = m_Fit_Data_ZTPC[m_currentPartIndex][0][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);

      }
      mu_pi = mu_pi + predictShift;
      mu_e = mu_e + predictShift;
      sig_e = sig_pi;
    }
    if(mu_e < m_currentHistoToFit->GetBinLowEdge(1)  || mu_e > m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX() + 1)){
      delete m_currentHistoToFit;
      continue;
    }

    double N_pion = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(mu_pi) ) /
    gaus(mu_pi, sig_pi, mu_pi);
    double N_electron = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(mu_e) ) /
    gaus(mu_e, sig_e, mu_e);
    cout << "  Npi: " << N_pion << " Ne: " << N_electron << endl;
      //pion
    fitFunctionTwo->SetParameter(0,N_pion);
    fitFunctionTwo->SetParLimits(0,0.0, 1.1*counts);
    fitFunctionTwo->SetParameter(1,mu_pi);
    fitFunctionTwo->SetParLimits(1,mu_pi - 0.1, mu_pi - 0.1);
    fitFunctionTwo->SetParameter(2,sig_pi);
    fitFunctionTwo->SetParLimits(2,0.9*sig_pi,1.1*sig_pi);
      //electron
    fitFunctionTwo->SetParameter(3,N_electron);
    fitFunctionTwo->SetParLimits(3,5.0, 1.1*counts);
    fitFunctionTwo->SetParameter(4,mu_e);
    fitFunctionTwo->SetParLimits(4,mu_e - 0.1, mu_e - 0.1);
    fitFunctionTwo->SetParameter(5,sig_e);
    fitFunctionTwo->SetParLimits(5,0.9*sig_e,1.1*sig_e);
    fitFunctionTwo->SetParNames("N_{#pi}","#mu_{#pi}","#sigma_{#pi}","N_{e}","#mu_{e}","#sigma_{e}");

    string fitOptions = "EQ";
    #ifdef _ZFITTER_DEBUG_
      fitOptions = "EV";
    #endif
    if(m_doIntegralOfBinForFits) fitOptions += "I";
    m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
    int status = m_currentHistoToFit->Fit("fitFunctionTwo",fitOptions.c_str());
    if(status != 0){
      cout << "ERROR: " << imageName.Data() << " did not fit Pion and Electron peak correctly" << endl;
      delete m_currentHistoToFit;
      continue;
    }

    TF1* fitResultTwo = m_currentHistoToFit->GetFunction("fitFunctionTwo");
    if(fitResultTwo->GetParameter(0) != 0.0){
      double orderCheck = (m_bichselCurves[0]->Eval(m_currentMomOfBin) - m_bichselCurves[3]->Eval(m_currentMomOfBin))*( fitFunctionTwo->GetParameter(1) - fitFunctionTwo->GetParameter(4));
      if(orderCheck > 0.0){
          /*m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(0));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(1));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(2));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(0));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(1));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(2));
          */
        if(m_currentPartCharge == 1){
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(3));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(4));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(5));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(3));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(4));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(5));
        }else{
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(3));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(4));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(5));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(3));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(4));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(5));
        }

      }else{
        if(m_currentPartCharge == 1){
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(0));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(1));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(2));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(0));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(1));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(2));
        }else{
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(0));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(1));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(2));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(0));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(1));
          m_Fit_Data_ZTPC[m_currentPartIndex][3][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(2));
        }
        fitFunctionTwo->SetParNames("N_{e}","#mu_{e}","#sigma_{e}","N_{#pi}","#mu_{#pi}","#sigma_{#pi}");
          /*m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(3));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(4));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParameter(5));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(3));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(4));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitResultTwo->GetParError(5));
          */
      }

    }

      //draw result
    if(a_saveImages){
      fitCanvas->cd();
      proj->SetStats(true);
      proj->Draw("E");
      fitResultTwo->SetNpx(500);
      fitResultTwo->Draw("same");
        //gPad->SetLogx();
      gPad->SetLogy();
      gPad->SetRightMargin(0.5);
      gPad->Update();

      DrawPredictionLines_dEdx(fitCanvas);
      fitCanvas->Update();
      gSystem->ProcessEvents();
        //fitCanvas->SaveAs(Form("%s/%s/%s/%s.png",imageDir.c_str(),m_partInfo->GetParticleName(m_currentPartIndex).Data(),"PionElectron",imageName.Data()));

      fitCanvas->SaveAs(Form("%s/%s/%s/ZTPC_BTOFIsolated/PionElectron/%s.png",m_imagePreDir.c_str(),m_imgDirName.c_str(),m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),imageName.Data()));

    }

    delete m_currentHistoToFit;
      //delete fitResultTwo;



      /*   PREVIOUS PION AND ELECTRON FITTING

      //cout << "\n\nHorBin: " << HorBinIndex << "  VertBin: " << binVertIndex << " ";
      TH1D* proj = current3DHisto->ProjectionZ(Form("dEdxElecPionProj_%02d_%02d", HorBinIndex , binVertIndex),
           HorBinIndex, HorBinIndex, binVertIndex, binVertIndex);
      double counts = HistogramUtilities::sumOfContent(proj, -1);

      if( minContentToConsider > counts ){  // && !m_testing)  || (counts < 20)){
        delete proj;
        continue;
      }
      proj->Sumw2();
      double mean = proj->GetMean();
      double stdDev = proj->GetStdDev();
      double var = stdDev*stdDev;
      double mu = 2.0*log(mean)-0.5*log(var+mean*mean);
      double sig = log(1.0 + var/(mean*mean));
      proj->Scale(1.0,"width");
      double integral = proj->Integral();
      double eta = current3DHisto->GetXaxis()->GetBinCenter(HorBinIndex);
      double pT = current3DHisto->GetYaxis()->GetBinCenter(binVertIndex);
      double m_currentMomOfBin = pT*TMath::CosH(eta);

      //FIT PIONS FIRST
      double ZTPCPredict = m_bichselCurves[0]->Eval(m_currentMomOfBin);
      if(m_currentMomOfBin < m_isolatedZTPC_starts[0]) continue;
      if(ZTPCPredict < proj->GetBinLowEdge(1)
         || ZTPCPredict > proj->GetBinLowEdge(proj->GetNbinsX() + 1)) continue;
      TString imageName = TString::Format("PionElectron_EtaBin_%02d_PtBin_%02d",HorBinIndex, binVertIndex);
      proj->SetName(imageName.Data());
      proj->SetTitle(imageName.Data());
   //          proj->SetStats(kTRUE);

      //cout << "mean: " << mean << "  stdDev: " << stdDev << "  mu: " << mu << "   sig: " << sig << endl;
      //fitFunction->SetParameter(0,(counts)*projTemp->GetBinWidth(2));
      fitFunction->SetParameter(0,counts);
      fitFunction->SetParLimits(0,0.1*counts, 2.0*counts);
      fitFunction->SetParameter(1,mu);
      fitFunction->SetParLimits(1,0.5*mu, 2.0*mu);
      fitFunction->SetParameter(2,sig);
      fitFunction->SetParLimits(2,0.5*sig,2.0*sig);
      fitFunction->SetParNames("amp","#mu","#sigma");
      //FIT Pion Peak Only
      int status = -1;
      status = proj->Fit("fitFunction","EQ","");
      if(status != 0){
        cout << "ERROR: " << imageName.Data() << " did not fit pion peak correctly" << endl;
        continue;
      }
      //cout << "Pion Peak Fit Successful!" << endl;
      TF1* fitResult = proj->GetFunction("fitFunction");
      double meanOfPion = exp(fitResult->GetParameter(1)
           +0.5*fitResult->GetParameter(2)*fitResult->GetParameter(2));
      double predictShift = ZTPCPredict - meanOfPion;
      double electron_ztpc = m_bichselCurves[3]->Eval(m_currentMomOfBin) + predictShift;
      double var_e = 0.07*0.07*electron_ztpc;
      double mu_e = 2.0*log(electron_ztpc)-0.5*log(var_e+electron_ztpc*electron_ztpc);
      double sig_e = log(1.0 + var_e/(electron_ztpc*electron_ztpc));
        //DRAWING
   //       TLine* predictLine = new TLine(electron_ztpc,electron_ztpc,0,10000);
   //          predictLine->SetVertical();
   //          predictLine->SetLineColor(kRed);
      fitCanvas->cd();
      gPad->SetLogy();
   //          proj->SetStats(kTRUE);
      gStyle->SetOptFit(0111);
      proj->Draw("E");
   //          predictLine->Draw("same");

      m_currentHistoToFit = proj;
      DrawPredictionLines(fitCanvas);
      TPaveStats *statsBox = (TPaveStats*)fitFunction->FindObject("stats");
      if(statsBox){
        statsBox->SetX1NDC(.15);
        statsBox->SetY1NDC(.75);
        statsBox->SetX2NDC(.4);
        statsBox->SetY2NDC(.9);
      }
      //now fit BOTH pions and electrons together
      //pion
      fitFunctionTwo->SetParameter(0,fitResult->GetParameter(0));
      fitFunctionTwo->SetParLimits(0,0.95 * fitResult->GetParameter(0), 1.05*fitResult->GetParameter(0));
      fitFunctionTwo->SetParameter(1,fitResult->GetParameter(1));
      fitFunctionTwo->SetParLimits(1,0.95*fitResult->GetParameter(1), 1.05*fitResult->GetParameter(1));
      fitFunctionTwo->SetParameter(2,fitResult->GetParameter(2));
      fitFunctionTwo->SetParLimits(2,0.95*fitResult->GetParameter(2), 1.05*fitResult->GetParameter(2));

      //electron
      fitFunctionTwo->SetParameter(3,fitResult->GetParameter(0)/1000.0);
      fitFunctionTwo->SetParLimits(3,0.05*fitFunctionTwo->GetParameter(3), 20.0*fitFunctionTwo->GetParameter(3));
      fitFunctionTwo->SetParameter(4,mu_e);
      fitFunctionTwo->SetParLimits(4,0.95*mu_e, 1.05*mu_e);
      fitFunctionTwo->SetParameter(5,sig_e);
      fitFunctionTwo->SetParLimits(5,0.95*sig_e, 1.05*sig_e);
      fitFunctionTwo->SetParNames("pion_amp","pion_#mu","pion_#sigma",
                       "electron_amp","electron_#mu","electron_#sigma");
      status = proj->Fit("fitFunctionTwo","EQ");
      if(status != 0){
        cout << "ERROR: " << imageName.Data() << " did not fit Pion and Electron peak correctly" << endl;
        continue;
      }
      TF1* fitResultTwo = proj->GetFunction("fitFunctionTwo");
      if(fitResultTwo->GetParameter(0) != 0.0){
        if(fitResultTwo->GetParameter(0) > fitResultTwo->GetParameter(3)){
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(0));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(1));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(2));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(0));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(1));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(2));

          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(3));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(4));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(5));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(3));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(4));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(5));
        }else{
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(0));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(1));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(2));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(0));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(1));
          m_dEdxFitData[m_currentPartIndex][partIndex][][3]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(2));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(3));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(4));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinContent(HorBinIndex,binVertIndex,fitResultTwo->GetParameter(5));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(3));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(4));
          m_dEdxFitData[m_currentPartIndex][partIndex][][0]->SetBinError(HorBinIndex,binVertIndex,fitResultTwo->GetParError(5));
        }

      }

      //draw result
      fitCanvas->cd();
      fitResultTwo->SetNpx(500);
      fitResultTwo->Draw("same");
      gPad->Update();
      m_currentHistoToFit = proj;
      DrawPredictionLines(fitCanvas);
      fitCanvas->Update();

      gSystem->ProcessEvents();
      fitCanvas->SaveAs(Form("%s/%s.png",imageDir.c_str(),imageName.Data()));
      delete proj;*/
    } // end vert bin
  }// end hor bin
  } //end if pion or kaon
  delete fitCanvas;


  //inFile->Close();

}







//============================================================================================================
//============================================================================================================

void ZFitter::fitTPC(bool a_saveImages, bool a_useProminenceAllCent){
  cout << "####################################################################################################" << endl;
  cout << "                        FITTING ZTPC INTEGRATED CENTRALITIES (ALL ADDED TOGETHER AS ONE FIT)  " << endl;
  cout << "####################################################################################################" << endl;

  m_inAllCentMode = true;
  m_currentCentIndex = -1;
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  gPad->SetLogy();
  gPad->SetRightMargin(0.45);


  #ifdef _MAC_OSX_
    if(m_useRootMinimizer){
      m_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
      m_minimizer->SetPrintLevel(4);
      m_minimizer->SetMaxIterations(3000);
      m_minimizer->SetMaxFunctionCalls(100000);
      m_minimizer->SetValidError(true);
    }
  #endif


  TH3D* current3DHisto = NULL;
  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    m_currentPlusMinusIndex = pmIndex;
    if(m_seedFromOppositeChargeDirection[m_currentPartIndex] < 0){
      if(pmIndex == 0) m_currentPlusMinusIndex = 1;
      if(pmIndex == 1) m_currentPlusMinusIndex = 0;
    }
    if(m_currentPlusMinusIndex == 0){
      if(a_useProminenceAllCent){
        current3DHisto = m_rap_mTm0_dEdx_AllCentRead_Plus[m_currentPartIndex];
      }else{
        current3DHisto = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex];
      }
      m_currentPartCharge = 1;
    }else{
      if(a_useProminenceAllCent){
        current3DHisto = m_rap_mTm0_dEdx_AllCentRead_Minus[m_currentPartIndex];
      }else{
        current3DHisto = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex];
      }
      m_currentPartCharge = -1;
    }

    #ifdef _ZFITTER_DEBUG_
      cout << " PM Index = " << m_currentPlusMinusIndex << " currentPartCharge= " << m_currentPartCharge << " ptr3dHisto = " << current3DHisto << endl;
    #endif
    if(!current3DHisto){
      cout << " ERROR: For some reason, the 3D Histogram for ZTPC All Cent is Null: " << current3DHisto << "   for part " << m_currentPartIndex << "   a_useProminenceAllCent: " << a_useProminenceAllCent << endl;

    }
    int numHorBins = current3DHisto->GetNbinsX();
    int numVertBins = current3DHisto->GetNbinsY();
    cout << m_partInfo->GetParticleName(m_currentPartIndex).Data() << " has " << numHorBins << " rapidity bins and " << numVertBins << " mT-m0 bins" << endl;
    for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
      for(int iii = 0; iii < 9; iii++){
        m_partHasBeenExtrapolated_Solo_Charge[iii][m_currentPlusMinusIndex] = false;
        m_partHasBeenExtrapolatedCharge[iii][m_currentPlusMinusIndex] = false;
      }
      m_currentRapidity = current3DHisto->GetXaxis()->GetBinCenter(HorBinIndex);
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
          if(HorBinIndex == m_rapBinsToTest[iii]){
            doThisRapBin = true;
          }
        }
        if(!doThisRapBin){
          continue;
        }
      }


      for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
        if(m_testing && VertBinIndex%m_divTestingMtM0!=0) continue;

        cout << "\n\n###########################################################################" << endl;
        cout << "Working on All Cent ZTPC " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << endl;

        m_currentRapBin = HorBinIndex;
        m_currentMtM0Bin = VertBinIndex;
        //m_currentRapidity = current3DHisto->GetXaxis()->GetBinCenter(m_currentRapBin);
        m_currentMtM0 = current3DHisto->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
        cout << "y = " << m_currentRapidity << " mTm0 = " << m_currentMtM0 << " mom =  " << m_currentMomOfBin << endl;

        if(     (m_currentPartCharge == 1  && m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin)
            ||  (m_currentPartCharge == -1 && m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin)){
          #ifdef _ZFITTER_DEBUG_
            cout << " Skipping bin because mom outside CutClass Mom cuts: " << m_currentMomOfBin << " < "
             << m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][m_currentPartIndex] << "  or  " << m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][m_currentPartIndex] << endl;
          #endif
          continue;
        }


        //###############   SETUP HISTOGRAM TO FIT ######################
        m_currentHistoToFit = current3DHisto->ProjectionZ(Form("dEdxProj_%02d_%02d_%02d",m_currentPartCharge, m_currentRapBin , m_currentMtM0Bin),
         m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);

        double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);

        if(countsInHisto < 30){
          #ifdef _ZFITTER_DEBUG_
            cout << " Skipping bin due to too few counts: " << countsInHisto << " < " << "30" << endl;
          #endif
          delete m_currentHistoToFit;
          continue;
        }
        m_currentIntegralofHistogram = countsInHisto;

        //TH1D* tempHistoPtr = m_currentHistoToFit;
        //m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,5.0);
        //delete tempHistoPtr;

        m_currentHistoToFit->Scale(1.0,"width");
        double lowRangeData = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistoToFit));
        double highRangeData = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistoToFit)+1);

        double lowRange = lowRangeData;
        double highRange = highRangeData;




        //If isolated with no mixing, just fit one gausian...
        // make sure the conditions for running fitIsolation_dEdx are the same here and down in the centrality specific fitting below
        /*
        if(m_currentMomOfBin < m_dEdxSeparatedStops[m_currentPartIndex] && m_currentPartIndex != 0){ //fit with one gaussian if not close to other particles / can't do pions because of electron contamination
          vector<double> fitData = fitIsolation_dEdx(m_currentPartIndex);
        if(m_currentPartCharge == 1){
          m_dEdxFitStatusHistoPlus[m_currentPartIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[6] + 5);
        }else{
          m_dEdxFitStatusHistoMinus[m_currentPartIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[6] + 5);
        }
        if(fitData[0] > 0){
          cout << "   Zoom Fit Success: N: " << fitData[0] << " mu:" << fitData[1] << " sig:" << fitData[2] << endl;
          if(pmIndex == 0){
            m_dEdxFitDataPlus[m_currentPartIndex][m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
            m_dEdxFitDataPlus[m_currentPartIndex][m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
            m_dEdxFitDataPlus[m_currentPartIndex][m_currentPartIndex][2]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
            m_dEdxFitDataPlus[m_currentPartIndex][m_currentPartIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
            m_dEdxFitDataPlus[m_currentPartIndex][m_currentPartIndex][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[4]);
            m_dEdxFitDataPlus[m_currentPartIndex][m_currentPartIndex][2]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[5]);
          }else{
            m_dEdxFitDataMinus[m_currentPartIndex][m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
            m_dEdxFitDataMinus[m_currentPartIndex][m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
            m_dEdxFitDataMinus[m_currentPartIndex][m_currentPartIndex][2]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
            m_dEdxFitDataMinus[m_currentPartIndex][m_currentPartIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
            m_dEdxFitDataMinus[m_currentPartIndex][m_currentPartIndex][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[4]);
            m_dEdxFitDataMinus[m_currentPartIndex][m_currentPartIndex][2]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[5]);
          }

        }*/




      m_currentNumPartToFit = 0;
      //################## SETUP PARTICLE TAGS ##################
      if(m_cuts->getEventConfig() == "FixedTarget"){
        //##########################################################
        //##########    FIXED TARGET PARTICLE TAGS   ###############
        //##########################################################
        bool atLeastOneNotKnown = false; //dont fit if soloFitterGot All Particles
        for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
          //see if there is data at prediction
          cout << "partIndex: " << partIndex << " particleToggle: " << m_initialParticleToggles[partIndex] << endl;
          /*
          if(!m_initialParticleToggles[partIndex]){
            m_currentParticleTags[partIndex] = 0;
            continue;
          }

          double TPCStop = m_cuts->getMomStops_TPC(partIndex,m_currentPartCharge);
          if(TPCStop < m_currentMomOfBin){
            m_currentParticleTags[partIndex] = 0;
            continue;
          }
          */

          // Dont fit partIndex if not in the fit matrix    or     mom is above limit
          if(m_currentPartCharge == 1  && (!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]
            || m_currentMomOfBin > m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex])){
            #ifdef _ZFITTER_DEBUG_
             cout << " Not fitting partIndex =  " << partIndex << " because of fit matrix or dEdxMomStop Matrix" << endl;
            #endif
            m_currentParticleTags[partIndex] = 0;
            continue;
          }
          if(m_currentPartCharge == -1 && (!m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]
            || m_currentMomOfBin > m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex])){
            #ifdef _ZFITTER_DEBUG_
             cout << " Not fitting partIndex =  " << partIndex << " because of fit matrix or dEdxMomStop Matrix" << endl;
            #endif
            m_currentParticleTags[partIndex] = 0;
            continue;
          }

          if(partIndex != 3){ // this particle is not an electron
            bool zoomedFitSuccess = false;
            /*if(m_currentMomOfBin < m_dEdxSeparatedStops[partIndex]){
              vector<double> fitData;

              cout << " Starting fitIsolation_dEdx with mom:" << m_currentMomOfBin << " with partIndex=" << partIndex << endl;
              cout << m_currentHistoToFit << endl;
              fitData = fitIsolation_dEdx(partIndex);
              if(fitData[0] > 0){
                m_dEdxFitData[m_currentPartIndex][partIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
                m_dEdxFitData[m_currentPartIndex][partIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
                m_dEdxFitData[m_currentPartIndex][partIndex][2]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
                m_currentParticleTags[partIndex] = 3;
                zoomedFitSuccess = true;
                continue;
              }
            }*/
            if(!zoomedFitSuccess){
              if( (m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0)
                   || (m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0) ){
                m_currentParticleTags[partIndex] = 9;
              }else{
                double dEdxPredict = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
                int dEdxPredictBin = m_currentHistoToFit->FindBin(dEdxPredict);
                double countsAtPredict = m_currentHistoToFit->GetBinContent(dEdxPredictBin)
                * m_currentHistoToFit->GetBinWidth(dEdxPredictBin);
                if(countsAtPredict > 5 && dEdxPredict > lowRange && dEdxPredict < highRange){

                      m_currentParticleTags[partIndex] = 2; //change to 2 when ready
                      atLeastOneNotKnown = true;
                    }else{
                      cout << "part" << partIndex << " getting 0" << endl;
                      m_currentParticleTags[partIndex] = 0;
                    }
                  }
                }

              }
          else{ //######## how to fit electrons #############
            m_currentParticleTags[3] = 1;  //the electron study fixes the means and widths of electrons use 5 to set means with mimic method
            double TPCStop = m_cuts->getMomStops_TPC(3,m_currentPartCharge);
            if(TPCStop < m_currentMomOfBin){
              m_currentParticleTags[3] = 0;
              continue;
            }
            //if electrons not fit, kill them for future fits in specific centrality
            if(m_currentParticleTags[3] == 0){
             if(m_currentPartCharge >  0){
               m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
               m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
               m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
             }else{
               m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
               m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
               m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
             }
           }
         }
       }

       #ifdef _ZFITTER_DEBUG_
        for(int iii = 0; iii < 8; iii++){
         if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][iii]) || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][iii])) continue;
         cout << "   " << m_partInfo->GetParticleName(iii) << "  has initial fit tag: " << m_currentParticleTags[iii] << " : " << m_particleTagDescriptions[m_currentParticleTags[iii]] << endl;
       }
       #endif


      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] == 0) continue;

        if( partIndex == 3 ){
          if( m_currentParticleTags[partIndex] == 9){
            m_currentParticleTags[partIndex] = 10;
            atLeastOneNotKnown = true;
          }else if( m_currentParticleTags[partIndex] == 2 ){
            m_currentParticleTags[partIndex] = 11;
            atLeastOneNotKnown = true;
          }
          continue;
        }

        if(!m_allowCentIntegratedFitsToVary){
          //original particle tags befor July 13, 2021
          if( m_currentParticleTags[partIndex] == 9 ){
            double partMean;
            double partMimicMean;
            double partElecMean = log(m_bichselCurves[3]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
            double partWidth;
            if(m_currentParticleTags[m_particleMimicIndicies[partIndex]] != 0){
            if(m_currentPartCharge == 1){
              partMean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partWidth = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partMimicMean = m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partWidth += m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }else{
              partMean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partWidth = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partMimicMean = m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partWidth += m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }
          }
            partWidth = 3.0*partWidth/2.0;
            int dEdxMeanBin = m_currentHistoToFit->FindBin(partMean);
            double countsAtMean = m_currentHistoToFit->GetBinContent(dEdxMeanBin)
                  * m_currentHistoToFit->GetBinWidth(dEdxMeanBin);
            if(m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 0 || m_particleMimicIndicies[partIndex] == partIndex ){
              m_currentParticleTags[partIndex] = 1;//was 9
              atLeastOneNotKnown = true;
            }else if( m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 2 ){
              m_currentParticleTags[partIndex] = 1; //was 10
              atLeastOneNotKnown = true;
            }else if( (m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 9 || m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 10 || m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 12) && hasTOFRatioToConstrainDeDx(partIndex) && (fabs(partMean - partMimicMean) < 0.25 || (partIndex == 1 &&fabs(partMean - partElecMean) < 0.25)) ){
              m_currentParticleTags[partIndex] = 13; //was 12
              atLeastOneNotKnown = true;
            }else{
              cout << "part " << partIndex <<" Mean: " << partMean << endl;
              cout << "part mimic " << m_particleMimicIndicies[partIndex] <<" Mean: " << partMimicMean << endl;
              cout << "fabs(partMean - partMimicMean): " << fabs(partMean - partMimicMean) << endl;
              cout << "fabs(partMean - partElecMean): " << fabs(partMean - partElecMean) << endl;
              cout << "hasTOFRatioToConstrainDeDx(" << partIndex << "): " << hasTOFRatioToConstrainDeDx(partIndex)   << endl;
              cout << "m_currentParticleTags[partIndex]: " << m_currentParticleTags[partIndex] << endl;
              cout << "m_currentParticleTags[m_particleMimicIndicies[partIndex]]: " << m_currentParticleTags[m_particleMimicIndicies[partIndex]] << endl;
              m_currentParticleTags[partIndex] = 1; //was 10
              atLeastOneNotKnown = true;
            }
          }else if( m_currentParticleTags[partIndex] == 2 ){
            if( m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 0 || m_particleMimicIndicies[partIndex] == partIndex ){
              m_currentParticleTags[partIndex] = 2;
              atLeastOneNotKnown = true;
            }else{
              m_currentParticleTags[partIndex] = 2; //was 11
              atLeastOneNotKnown = true;
            }

          }
        }//end !m_allowCentIntegratedFitsToVary
        else{
          if( m_currentParticleTags[partIndex] != 9 ) continue;
          double partMean;
          double partMimicMean;
          double partElecMean = log(m_bichselCurves[3]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
          double partWidth;
          if(m_currentParticleTags[m_particleMimicIndicies[partIndex]] != 0){
          if(m_currentPartCharge == 1){
            partMean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            partWidth = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            partMimicMean = m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            partWidth += m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }else{
            partMean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            partWidth = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            cout << "m_particleMimicIndicies[" << partIndex << "]: " << m_particleMimicIndicies[partIndex] << endl;
            cout << "m_Fit_Data_ZTPC[" << m_currentPartIndex << "][" << m_particleMimicIndicies[partIndex] << "][1][1]: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][1] << endl;
            cout << "value: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) << endl;
            partMimicMean = m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            partWidth += m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }
        }
          partWidth = 3.0*partWidth/2.0;
          int dEdxMeanBin = m_currentHistoToFit->FindBin(partMean);
          double countsAtMean = m_currentHistoToFit->GetBinContent(dEdxMeanBin)
                * m_currentHistoToFit->GetBinWidth(dEdxMeanBin);
          if(m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 0 || m_particleMimicIndicies[partIndex] == partIndex ){
            m_currentParticleTags[partIndex] = 2;
            atLeastOneNotKnown = true;
          }else if((fabs(partMean - partMimicMean) < 0.25 || (partIndex == 1 &&fabs(partMean - partElecMean) < 0.25))){
            if(hasTOFRatioToConstrainDeDx(partIndex)){
              m_currentParticleTags[partIndex] = 13;
              atLeastOneNotKnown = true;
            }else{
              m_currentParticleTags[partIndex] = 1;
              atLeastOneNotKnown = true;
            }
          }else{
            m_currentParticleTags[partIndex] = 2;
            atLeastOneNotKnown = true;
          }
        }

        /*
        if(partIndex == 2){
          cout << "proton particle tag: " << m_currentParticleTags[partIndex] << endl;
          atLeastOneNotKnown = true;
        }
        if( m_currentParticleTags[partIndex] == 9 && m_currentParticleTags[m_particleMimicIndicies[partIndex]] != 0 &&m_currentParticleTags[m_particleMimicIndicies[partIndex]] != 110){
          double partMean;
          double partMimicMean;
          if(m_currentPartCharge == 1){
            partMean = m_dEdxFitDataPlus[m_currentPartIndex][partIndex][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            partMimicMean = m_dEdxFitDataPlus[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }else if(m_currentPartCharge == -1){
            partMean = m_dEdxFitDataMinus[m_currentPartIndex][partIndex][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            partMimicMean = m_dEdxFitDataMinus[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }
          if( m_particleMimicIndicies[partIndex] == partIndex && m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 9 && fabs(partMean - partMimicMean)<0.1 ){
            m_currentParticleTags[partIndex] = 12;
            atLeastOneNotKnown = true;
          }else{
            m_currentParticleTags[partIndex] = 10;
            atLeastOneNotKnown = true;
          }
        }else if(m_currentParticleTags[partIndex] == 9){
          if(partIndex == 0) m_currentParticleTags[partIndex] = 2;
          if(partIndex == 1){
            if(hasTOFRatioToConstrainDeDx(partIndex)) m_currentParticleTags[partIndex] = 12;
            else m_currentParticleTags[partIndex] = 10;
          }
        }
        */
        /*
        if( m_particleMimicIndicies[partIndex] == partIndex ){
          if((m_currentPartCharge == 1 && m_dEdxFitDataPlus[m_currentPartIndex][partIndex][2]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0) || (m_currentPartCharge == -1 && m_dEdxFitDataMinus[m_currentPartIndex][partIndex][2]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0)){
            m_currentParticleTags[partIndex] = 9;
            atLeastOneNotKnown = true;
          }else{
            m_currentParticleTags[partIndex] = 2;
            atLeastOneNotKnown = true;
          }
        }
        if(m_currentParticleTags[partIndex] == 9){
          if(m_currentParticleTags[m_particleMimicIndicies[partIndex]] != 9){
            m_currentParticleTags[partIndex] = 10;
            atLeastOneNotKnown = true;
          }else{
            double partMean;
            double partMimicMean;
            if(m_currentPartCharge == 1){
              partMean = m_dEdxFitDataPlus[m_currentPartIndex][partIndex][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partMimicMean = m_dEdxFitDataPlus[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }else if(m_currentPartCharge == -1){
              partMean = m_dEdxFitDataMinus[m_currentPartIndex][partIndex][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              partMimicMean = m_dEdxFitDataMinus[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            }
            if(fabs(partMean - partMimicMean)<0.1){
              m_currentParticleTags[partIndex] = 12;
              atLeastOneNotKnown = true;
            }else{
              m_currentParticleTags[partIndex] = 10;
              atLeastOneNotKnown = true;
            }
          }
        } */
      }







          //this particle will not have any free parameters and will mimic the particle defined in m_particleMimicIndicies
          //cout << "Check: PartIndex " << partIndex << " with tag " <<  m_currentParticleTags[partIndex] << endl;
          //cout << "      mom: " << m_currentMomOfBin << " > " << m_particleMimicLowMom[partIndex] << endl;
          // cout << "     mom: " << m_currentMomOfBin << " < " << m_particleMimicHighMom[partIndex] << endl;
          // cout << "     mimc Particle: " << m_particleMimicIndicies[partIndex] << " with tag " << m_currentParticleTags[m_particleMimicIndicies[partIndex]] << endl;
        //cout << "Current Particle Tags: ";
        //for(int partIndex = 0; partIndex < m_numParticles; partIndex++) cout << m_currentParticleTags[partIndex] << "  ";
        //cout << endl;

        //initialize parameters
        if(!atLeastOneNotKnown) continue; // no use fitting if all widths are known













      }else if(m_cuts->getEventConfig() == "ColliderCenter"){
        //##########################################################
        //##########    COLLIDER CENTER PARTICLE TAGS   ############
        //##########################################################
        //cout << "Setting up Collider Center Particle Tags" << endl;
        bool atLeastOneNotKnown = false; //dont fit if soloFitterGot All Particles
        for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
          #ifdef _ZFITTER_DEBUG_
             cout << "Setting up Collider Center Particle Tag for " << partIndex << endl;
          #endif
          m_currentParticleTags[partIndex] = 0; //INITAL ASSUMPTION IS DON'T FIT

          // Dont fit partIndex if cutclass says not to
          //if(m_cuts->getMomStops_TPC(partIndex,m_currentPartCharge) < m_currentMomOfBin) continue;

          // Dont fit partIndex if not in the fit matrix    or     mom is above limit
          if(m_currentPartCharge == 1  && (!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]
            || m_currentMomOfBin > m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex])){
            #ifdef _ZFITTER_DEBUG_
             cout << " Not fitting partIndex =  " << partIndex << " because of fit matrix or dEdxMomStop Matrix" << endl;
            #endif
            continue;
          }
          if(m_currentPartCharge == -1 && (!m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]
            || m_currentMomOfBin > m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex])){
            #ifdef _ZFITTER_DEBUG_
             cout << " Not fitting partIndex =  " << partIndex << " because of fit matrix or dEdxMomStop Matrix" << endl;
            #endif
            continue;
          }


          //Isolate the particle of interest if it's all alone in dE/dx space with little or no overlap
          //if(m_currentMomOfBin < m_dEdxSeparatedStops[m_currentPartIndex] && m_currentPartIndex != 0 && m_currentPartIndex != partIndex) continue;
          if(m_currentMomOfBin < m_dEdxSeparatedStops[m_currentPartIndex][m_currentPlusMinusIndex] && m_currentPartIndex != partIndex){
            #ifdef _ZFITTER_DEBUG_
             cout << " Not fitting " << m_partInfo->GetParticleName(partIndex) << " because " << m_partInfo->GetParticleName(m_currentPartIndex) << " is separated enough to be fit on it's own" << endl;
            #endif
            continue;
          }

          if(m_usePearsonIVinSingle && m_currentMomOfBin < m_PearsonZTPCMaxMomStop[m_currentPartIndex][m_currentPlusMinusIndex] && m_currentPartIndex == partIndex){
            #ifdef _ZFITTER_DEBUG_
             cout << " Starting fitIsolation_dEdx with mom:" << m_currentMomOfBin << " with partIndex=" << partIndex << endl;
             cout << "Histogram Pointer: " << m_currentHistoToFit << endl;
            #endif

            bool zoomedFitSuccess = false;
            vector<double> fitData = fitIsolation_ZTPC_Pearson();
            if(fitData[0] > 0.0){
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[4]);

              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[5]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[6]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[7]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
              m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);

              zoomedFitSuccess = true;
              m_currentParticleTags[partIndex] = 0;
              continue;
            }
            if(!zoomedFitSuccess){
              m_currentParticleTags[partIndex] = 0;
              continue;
            }
          }


          //Extrapolate Means and Widths if going under another peak
          if(m_allowDeDxParameterExtrapolation && m_currentPartCharge == 1 && m_currentMomOfBin > m_partDeDx_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]){
            if(!m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]) extrapolateMeanAndWidthParams_TPC(partIndex);
            m_currentParticleTags[partIndex] = 1;
            #ifdef _ZFITTER_DEBUG_
             cout << " Fitting with tag 1 because extrapolated " << endl;
            #endif
            continue;
          }

          if(m_allowDeDxParameterExtrapolation && m_currentPartCharge == -1 && m_currentMomOfBin > m_partDeDx_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex]){
            if(!m_partHasBeenExtrapolatedCharge[partIndex][m_currentPlusMinusIndex]) extrapolateMeanAndWidthParams_TPC(partIndex);
            m_currentParticleTags[partIndex] = 1;
            #ifdef _ZFITTER_DEBUG_
              cout << " Fitting with tag 1 because extrapolated " << endl;
            #endif
            continue;
          }

          #ifdef _ZFITTER_DEBUG_
             cout << " Should I extrapolate the solo fit data?" << endl;
             cout << "      allow cent to vary: " << m_allowCentIntegratedFitsToVary << endl;
             cout << "      am_calculateBTOFIsolatedParameterDifferencesToTotalFitso : " << m_calculateBTOFIsolatedParameterDifferencesToTotalFits << endl;
             cout << "      m_partHasBeenExtrapolated_Solo_Charge[partIndex][m_currentPlusMinusIndex] : " << m_partHasBeenExtrapolated_Solo_Charge[partIndex][m_currentPlusMinusIndex]  << endl;
             cout << "       Current Mom " << m_currentMomOfBin << endl;
             cout << "      dEdxEstmatrix  = " << m_partDeDx_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]   <<   "   or     " << m_partDeDx_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex] << endl;
          #endif

          if(m_allowCentIntegratedFitsToVary && m_calculateBTOFIsolatedParameterDifferencesToTotalFits
            &&!m_partHasBeenExtrapolated_Solo_Charge[partIndex][m_currentPlusMinusIndex]
            &&  ((m_currentPartCharge == -1    &&       m_currentMomOfBin > m_partDeDx_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex])
                || (m_currentPartCharge == 1  &&       m_currentMomOfBin > m_partDeDx_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex]))
            ){
              extrapolateMeanAndWidthParams_FromSoloTPC_TPC(partIndex);
          }

          if(m_partHasBeenExtrapolated_Solo_Charge[partIndex][m_currentPlusMinusIndex]){
              m_currentParticleTags[partIndex] = 1;
              continue;
          }


          if(partIndex != 3){ // this particle is not an electron

            //if(m_currentMomOfBin < m_dEdxSeparatedStops[m_currentPartIndex] && m_currentPartIndex != 0) continue;
            double ztpc_predict = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
            if(m_bichselShiftsComplete){
              if(m_currentPartCharge == 1){
                if(m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin) != 0.0){
                  ztpc_predict += m_bichselShifts_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
                }else{
                  ztpc_predict += m_bichselShifts_Plus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
                }
              }else{
                if(m_tpcResolution_Plus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin) != 0.0){
                  ztpc_predict += m_bichselShifts_Minus[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin);
                }else{
                  ztpc_predict += m_bichselShifts_Minus[m_currentPartIndex][m_currentPartIndex]->GetBinContent(m_currentRapBin);
                }
              }
            }

            int predictBin = m_currentHistoToFit->FindBin(ztpc_predict);
            double countsAtPredict = m_currentHistoToFit->GetBinContent(predictBin) * m_currentHistoToFit->GetBinWidth(predictBin);
            #ifdef _ZFITTER_DEBUG_
             cout << " Z_TPC Predition: " << log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin)/m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin))
                  << " shifted to " << ztpc_predict << " which is bin: " << predictBin << " that has " << countsAtPredict << " counts" << endl;
            #endif


            if(m_allowCentIntegratedFitsToVary && m_currentMomOfBin < m_dEdxMaximumMomForTotalFit[partIndex][m_currentPlusMinusIndex]
               && ztpc_predict > lowRange && ztpc_predict < highRange && countsAtPredict > 5.0){
                 m_currentParticleTags[partIndex] = 2;
                 atLeastOneNotKnown = true;
            }else{
              if(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0 ){
                m_currentParticleTags[partIndex] = 1;
              }else if(countsAtPredict > 5.0 && ztpc_predict > lowRange && ztpc_predict < highRange){
                m_currentParticleTags[partIndex] = 2;
                atLeastOneNotKnown = true;
              }
            }

            /*if(m_currentMomOfBin > m_dEdxMaximumMomForTotalFit[partIndex][m_currentPlusMinusIndex] &&
                  m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0
                  ){
              if(!m_allowCentIntegratedFitsToVary ){
                m_currentParticleTags[partIndex] = 1; //use tof isolated fit
              }else{
                if(ztpc_predict > lowRange && ztpc_predict < highRange && countsAtPredict > 5.0){
                  m_currentParticleTags[partIndex] = 2;
                  atLeastOneNotKnown = true;
                }
              }
            }else{
              if(countsAtPredict > 5.0 && ztpc_predict > lowRange && ztpc_predict < highRange){
                m_currentParticleTags[partIndex] = 2;
                atLeastOneNotKnown = true;
              }
            }*/


          }else{ //######## how to fit electrons #############
            #ifdef _ZFITTER_DEBUG_
              cout << "  electron stdDev : " << m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) << endl;
              cout << "  electron mean : " << m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) << endl;
            #endif


            m_currentParticleTags[3] = 1;  //the electron study fixes the means and widths of electrons use 5 to set means with mimic method

            if( (m_currentPartCharge == 1  && (!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex]
                   || m_currentMomOfBin > m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex]))
              ||  (m_currentPartCharge == -1  && (!m_partToFitMatrix_Minus[m_currentPartIndex][partIndex]
                    || m_currentMomOfBin > m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex]))
              ){
              m_currentParticleTags[3] = 0;
              //if electrons not fit, kill them for future fits in specific centrality
              m_Fit_Data_ZTPC[m_currentPartIndex][3][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
              m_Fit_Data_ZTPC[m_currentPartIndex][3][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);
              m_Fit_Data_ZTPC[m_currentPartIndex][3][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, 0);

              continue; //particle loop
            }

          }
        }//end particle loop

        /*for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
          if(m_currentParticleTags[partIndex] == 0) continue;
          //this particle will not have any free parameters and will mimic the particle defined in m_particleMimicIndicies
          //cout << "Check: PartIndex " << partIndex << " with tag " <<  m_currentParticleTags[partIndex] << endl;
          //cout << "      mom: " << m_currentMomOfBin << " > " << m_particleMimicLowMom[partIndex] << endl;
          // cout << "     mom: " << m_currentMomOfBin << " < " << m_particleMimicHighMom[partIndex] << endl;
          // cout << "     mimc Particle: " << m_particleMimicIndicies[partIndex] << " with tag " << m_currentParticleTags[m_particleMimicIndicies[partIndex]] << endl;
          if( m_particleMimicLowMom[partIndex] < m_currentMomOfBin
              && m_particleMimicHighMom[partIndex] > m_currentMomOfBin
              && m_currentParticleTags[m_particleMimicIndicies[partIndex]] != 0){
            if(hasTOFRatioToConstrainDeDx(partIndex) && m_fullParticleMimic[partIndex]){
              m_currentParticleTags[partIndex] = 4;
              atLeastOneNotKnown = true;
            }else if(hasTOFRatioToConstrainDeDx(partIndex) && !m_fullParticleMimic[partIndex]){
              m_currentParticleTags[partIndex] = 7;
              atLeastOneNotKnown = true;
            }else{
              m_currentParticleTags[partIndex] = 6;
              atLeastOneNotKnown = true;
            }
          }
        }*/
        //cout << "Current Particle Tags: ";
        //for(int partIndex = 0; partIndex < m_numParticles; partIndex++) cout << m_currentParticleTags[partIndex] << "  ";
        //cout << endl;

        //if(m_useRootMinimizer){
        //  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        //    if(m_currentParticleTags[partIndex] == 0) continue;
        //    m_currentParticleTags[partIndex] = 2;
        //  }
        //}

        if(!atLeastOneNotKnown){
          #ifdef _ZFITTER_DEBUG_
            cout << "Skipping bin because there are no particles to fit" << endl;
          #endif
          continue; // no use fitting if all widths are known
        }

      }
      //########## END PARTICLE TAG SETUP ##############
      m_currentNumPartToFit = 0;
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] != 0) m_currentNumPartToFit = m_currentNumPartToFit + 1;
      }




      if(m_currentParticleTags[m_currentPartIndex] == 0){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping bin because main particle will not be fit (tag = 0)" << endl;
        #endif
        if(a_saveImages){
          fittingCanvas->cd();
          m_currentHistoToFit->Draw("E");
          DrawPredictionLines_dEdx(fittingCanvas);
          fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/Rap_%02d_mTm0_%02d___Main_Part_Not_Fit.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
               m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin));
        }
        continue;
      }

      #ifdef _ZFITTER_DEBUG_
      for(int iii = 0; iii < 8; iii++){
        if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][iii])
              || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][iii]) )continue;
        cout << "   " << m_partInfo->GetParticleName(iii) << "  has fit tag: " << m_currentParticleTags[iii] << " : " << m_particleTagDescriptions[m_currentParticleTags[iii]] << endl;
      }
      #endif







      //###########     SET THE HISTOGRAM RANGE TO FIT    ###############
      getRangeZTPC(lowRange, highRange);
      if(lowRange >= highRange){
        lowRange = lowRangeData;
        highRange = highRangeData;
      }

      if(m_currentNumPartToFit == 1){
        double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(-0.068), 2);
        double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.0), 2);
        double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.068), 2);
        double N_3p, mean_3p, sig_3p = 0.0;
        PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,-0.068,y1,0.0,y2,0.068,y3);
        cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
        if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
          lowRange  = mean_3p  - 3.0*sig_3p;
          highRange  = mean_3p + 3.0*sig_3p;
          cout << "     After: " << lowRange << " to " << highRange << endl;
        }else{
          cout << "   3pt gaus failed, just not zooming in" << endl;
        }
      }

      cout << " m_useRLBiasedSingleParticleZoom: " << m_useRLBiasedSingleParticleZoom << "   m_currentPartIndex: " << m_currentPartIndex << endl;
      if(m_useRLBiasedSingleParticleZoom && m_currentPartIndex == 0){
        double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(-0.068), 2);
        double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.0), 2);
        double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.068), 2);
        double N_3p, mean_3p, sig_3p = 0.0;
        PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,-0.068,y1,0.0,y2,0.068,y3);
        cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
        if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
          if(m_currentMomOfBin < 0.2){ // fit right side of peak
            lowRange  = mean_3p  - 0.5*sig_3p;
            highRange  = mean_3p + 3.0*sig_3p;
          }else{ // fit left side of peak
            lowRange  = mean_3p  - 3.0*sig_3p;
            highRange  = mean_3p + 0.5*sig_3p;
          }
          cout << "     After: " << lowRange << " to " << highRange << endl;
          for(unsigned int jjj = 0; jjj < 9; jjj++){ // kill other particles, just fit one peak
            if(m_currentPartIndex != jjj) m_currentParticleTags[jjj] = 0;
          }
          m_currentNumPartToFit = 1;
        }else{
          cout << "   3pt gaus failed, just not zooming in" << endl;
        }
      }else if(m_currentPartIndex != 0 && m_useRLBiasedSingleParticleZoom){
        double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(-0.068), 2);
        double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.0), 2);
        double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.068), 2);
        double N_3p, mean_3p, sig_3p = 0.0;
        PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,-0.068,y1,0.0,y2,0.068,y3);
        cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
        if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
          lowRange  = mean_3p  - 0.5*sig_3p;
          highRange  = mean_3p + 3.0*sig_3p;
          cout << "     After: " << lowRange << " to " << highRange << endl;
          for(unsigned int jjj = 0; jjj < 9; jjj++){ // kill other particles, just fit one peak
            if(m_currentPartIndex != jjj) m_currentParticleTags[jjj] = 0;
          }
          m_currentNumPartToFit = 1;
        }else{
          cout << "   3pt gaus failed, just not zooming in" << endl;
        }
      }







      m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
      m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);
      m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
      m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,0.001);


      TH1D* origHistToFit = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
      delete origHistToFit;
      m_currentHistoToFit->SetTitle(Form("Z_{TPC} for y=[%1.2f,%1.2f] m_{T}-m_{%s}=[%2.3f,%2.3f] p=%2.3f",
        m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
        m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
        m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data(),
        m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
        m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
        m_currentMomOfBin
        ));
      m_currentHistoToFit->GetXaxis()->SetTitle(Form("Z_{TPC} of %s", m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data()));
      m_currentHistoToFit->GetYaxis()->SetTitle("dN/dZ_{TPC}");

      double truncatedTotalCounts = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
      double possibleNewBins = truncatedTotalCounts/6.0;
      if(possibleNewBins > 80.0) possibleNewBins = 80.0;
      TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(m_currentHistoToFit,(int) possibleNewBins,true,true);
      delete m_currentHistoToFit;
      m_currentHistoToFit = rebinHisto;


      m_currentHistoToFit->SetMarkerStyle(20);
      m_currentHistoToFit->SetMarkerColor(kBlack);

      TF1* fitFunct = NULL;
      TF1* fitFunctInitial = NULL;
      TF1* fitResult = NULL;

      initZTPCParams();
      if(m_currentNumPartToFit == 0){
        if(a_saveImages){
          fittingCanvas->cd();
          m_currentHistoToFit->Draw("E");
          if(m_drawInitialSeedsToFits && fitFunctInitial) fitFunctInitial->Draw("same");
          DrawPredictionLines_dEdx(fittingCanvas);
          fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/Rap_%02d_mTm0_%02d___No_Particle_to_Fit.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
               m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin));
        }
        continue;
      }

      if(areDeDxParamsRidiculous()){
        cout << "###########       ERROR      dEdx parameters are ridiculous ################ " << endl;
        printCurrentParams();
        if(a_saveImages){
          fittingCanvas->cd();
          m_currentHistoToFit->Draw("E");
          if(m_drawInitialSeedsToFits && fitFunctInitial) fitFunctInitial->Draw("same");
          DrawPredictionLines_dEdx(fittingCanvas);
          fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/Rap_%02d_mTm0_%02d___BAD_Params.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
               m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin));
        }
        continue;
      }
      //cout << "Fitting " << m_currentNumPartToFit << " Particles" << "   With " << m_currentFreeParameters << " Parameters" << endl;
      //printCurrentParams();
      //do the fit
      fitFunct = new TF1("fitFunct",this,&ZFitter::ZTPCOperator,lowRange,highRange,m_currentFreeParameters,"ZFitter","ZTPCOperator");
      fitFunct->SetNpx(1200);


      vector< double > finalParamsRootMinimizer;
      vector< double > finalParamsGausMixModel(3*m_currentNumPartToFit, 0.0);

       //#############   MINUIT2 ROOT MINIMIZATION #################
      if(m_useRootMinimizer){
        #ifdef _ZFITTER_DEBUG_
          cout << "Using Root Minimizer ROOT6" << endl;
        #endif

        #ifdef _MAC_OSX_
          ROOT::Math::Functor chiSquaredFunctor(this,&ZFitter::chiSquared_TPC,m_currentFreeParameters);
          m_minimizer->SetFunction(chiSquaredFunctor);
          initMinimizer_TPC();

          m_minimizer->Minimize();
          int status = m_minimizer->Status();
          //if(status != 1) continue;
          const double* finalParams = m_minimizer->X();
          for(int iii = 0; iii < m_currentFreeParameters; iii++){
            finalParamsRootMinimizer.push_back(finalParams[iii]);
          }

          fitResult = (TF1*) fitFunct->Clone();
          fitResult->SetLineColor(kRed);
          for(int iii = 0; iii < m_currentFreeParameters; iii++){
            fitResult->SetParameter(iii,finalParams[iii]);
            fitResult->SetParError(iii,sqrt(m_minimizer->CovMatrix(iii,iii)));
        }
     	#endif

      //#############   GAUSSIAN MIXTURE MODEL #################
      }else if(m_useGausMixModel){
        #ifdef _ZFITTER_DEBUG_
          cout << "Using Gaussian Mixture Model" << endl;
        #endif
          vector<double> initialMeans(m_currentNumPartToFit,0.0);
          vector<double> initialAmps(m_currentNumPartToFit,0.0);

          int liveParticleCounter = 0;
          double ampTotalTemp = 0.0;
          for(int partIndex = 0; partIndex < 9; partIndex++){
            if(m_currentParticleTags[partIndex] != 0){
              initialMeans[liveParticleCounter] = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin) / m_bichselCurves[m_currentPartIndex]->Eval(m_currentMomOfBin));
              initialAmps[liveParticleCounter] = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(initialMeans[liveParticleCounter]));
              ampTotalTemp += initialAmps[liveParticleCounter];
              liveParticleCounter++;
            }
          }

          for(int iii = 0; iii < m_currentNumPartToFit; iii++){
            initialAmps[liveParticleCounter] = initialAmps[liveParticleCounter] / ampTotalTemp;
          }

          m_gausMixModel->doGMM(&finalParamsGausMixModel[0], &initialMeans[0],&initialAmps[0],300,1e-7, 0.2);
          m_currentParameters = finalParamsGausMixModel;

      }else{
        initZTPCFunction(fitFunct);
        fitFunctInitial = (TF1*) fitFunct->Clone();
        fitFunctInitial->SetLineColor(kYellow);
        fitFunctInitial->SetName("fitFunctInitial");

        string fitOptionsZTPCAllCent = "EQ";
        #ifdef _ZFITTER_DEBUG_
        fitOptionsZTPCAllCent = "EV";
        #endif
        if(m_doIntegralOfBinForFits) fitOptionsZTPCAllCent += "I";
        m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
        int status = m_currentHistoToFit->Fit("fitFunct",fitOptionsZTPCAllCent.c_str(),"",lowRange,highRange); //add Q to quiet
        fitResult = m_currentHistoToFit->GetFunction("fitFunct");

        if(m_currentPartCharge == 1){
          m_FitStatus_ZTPC[m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status + 5);
          m_ChiSqr_ZTPC[m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,((double) fitResult->GetChisquare())/((double)fitResult->GetNDF()) );
        }else{
          m_FitStatus_ZTPC[m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status + 5);
          m_ChiSqr_ZTPC[m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,((double) fitResult->GetChisquare())/((double)fitResult->GetNDF()) );
        }

        if(status != 0){
          cout << "Status: " << status << endl;
          cout << "\n\n\n\n\n\n BAD FIT \n\n\n\n\n\n";
          if(a_saveImages){
            fittingCanvas->cd();
            m_currentHistoToFit->Draw("E");
            if(m_drawInitialSeedsToFits && fitFunctInitial) fitFunctInitial->Draw("same");
            DrawPredictionLines_dEdx(fittingCanvas);
            fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/Rap_%02d_mTm0_%02d___BAD_FIT.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
                 m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin));
          }
          if(!m_storeFitsNotStatusZero) continue;
        }
      }


      //int paramIndexCounter = 0;
      ////cout << "Final Parameters::";
      //for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
      //  if(iii%3==0) cout << endl;
      //  if(m_paramToggles[iii]){
      //    m_currentParameters[iii] = fitResult->GetParameter(paramIndexCounter);
      //    paramIndexCounter++;
      //  }
      //  //cout << m_currentParameters[iii] << "  ";
      //}
      ////cout << endl;



      vector< double > allParamError; // even not fit
      int paramIndexCounter = 0;
      int liveParticleCounter = 0;
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] == 1){
          m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          paramIndexCounter++;
          liveParticleCounter++;
        }else if(m_currentParticleTags[partIndex] == 2){
          m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter);
          m_currentParameters[3*liveParticleCounter+1] = fitResult->GetParameter(paramIndexCounter+1);
          m_currentParameters[3*liveParticleCounter+2] = fitResult->GetParameter(paramIndexCounter+2);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
          allParamError.push_back(fitResult->GetParError(paramIndexCounter+2));
          paramIndexCounter+=3;
          liveParticleCounter++;
        }else if(m_currentParticleTags[partIndex] == 3){
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          liveParticleCounter++;
          //no free parameters
        }else if(m_currentParticleTags[partIndex] == 4){
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          liveParticleCounter++;
          // no free paramters
        }else if(m_currentParticleTags[partIndex] == 5){
          m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          liveParticleCounter++;
          paramIndexCounter++;
          // no free paramters
        }else if(m_currentParticleTags[partIndex] == 6){
          m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          liveParticleCounter++;
          paramIndexCounter++;
          // no free paramters
        }else if(m_currentParticleTags[partIndex] == 7){
          m_currentParameters[3*liveParticleCounter+1] = fitResult->GetParameter(paramIndexCounter);
          m_currentParameters[3*liveParticleCounter+2] = fitResult->GetParameter(paramIndexCounter+1);
          allParamError.push_back(0.0);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
          paramIndexCounter+=2;
          liveParticleCounter++;
        }else if(m_currentParticleTags[partIndex] == 9){
          m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter);
          m_currentParameters[3*liveParticleCounter+2] = fitResult->GetParameter(paramIndexCounter+1);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(0.0);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
          paramIndexCounter+=2;
          liveParticleCounter++;
        }else if(m_currentParticleTags[partIndex] == 10){
          m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          paramIndexCounter++;
          liveParticleCounter++;
        }else if(m_currentParticleTags[partIndex] == 11){
          m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter);
          m_currentParameters[3*liveParticleCounter+1] = fitResult->GetParameter(paramIndexCounter+1);
          allParamError.push_back(fitResult->GetParError(paramIndexCounter));
          allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
          allParamError.push_back(0.0);
          paramIndexCounter+=2;
          liveParticleCounter++;
        }else if(m_currentParticleTags[partIndex] == 12){
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          liveParticleCounter++;
        }else if(m_currentParticleTags[partIndex] == 13){
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          allParamError.push_back(0.0);
          liveParticleCounter++;
        }
      }
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] == 4 || m_currentParticleTags[partIndex] == 7 || m_currentParticleTags[partIndex] == 12 || m_currentParticleTags[partIndex] == 13){
          allParamError[getCurrParamIndex(partIndex)] = allParamError[getCurrParamIndex(m_particleMimicIndicies[partIndex])]
          /m_currentParameters[getCurrParamIndex(m_particleMimicIndicies[partIndex])];
          // no free paramters
        }
      }
      cout << "Particle Tags: ";
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++) cout << m_currentParticleTags[partIndex] << "  ";
        cout << endl;
      printCurrentParams();
      cout << "Parameter Errors: " << endl;
      for(int iii = 0; iii < (int) allParamError.size(); iii++) cout << allParamError[iii] << "  ";
        cout << endl;



      //######## Calculate Chi Squared ##############
      double cNDF = 0.0;

      if(m_useRootMinimizer){
	#ifdef _MAC_OSX_
        cNDF = chiSquared_TPC(&finalParamsRootMinimizer[0]) / ( (double) (HistogramUtilities::numBinsWithContentAboveThreshold(m_currentHistoToFit,1.0) + m_currentFreeParameters));
	#endif
      }else if(m_useGausMixModel){
        cNDF = m_gausMixModel->getChiSquared() / ( (double) (HistogramUtilities::numBinsWithContentAboveThreshold(m_currentHistoToFit,1.0) + 3*m_currentNumPartToFit));
      #ifdef _MAC_OSX_
      }else if(m_useRootMinimizer){
        cNDF = chiSquared_TPC(&finalParamsRootMinimizer[0]) / ( (double) (HistogramUtilities::numBinsWithContentAboveThreshold(m_currentHistoToFit,1.0) + m_currentFreeParameters));
      #endif
      }else{
        cNDF = fitResult->GetChisquare() / fitResult->GetNDF();
      }
      m_ChiSqr_ZTPC[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, cNDF);

      //push to particle fits

      double tempParamsInOrder[3*m_numParticles];
      int paramIndex = 0; //the particles not included in fit
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] !=0){

          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+1]);
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+2]);
          if(allParamError[paramIndex] > 0.0 ) m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, allParamError[paramIndex]* m_currentIntegralofHistogram);
          if(allParamError[paramIndex + 1] > 0.0 ) m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, allParamError[paramIndex+1]);
          if(allParamError[paramIndex + 2] > 0.0 ) m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, allParamError[paramIndex+2]);

          paramIndex+=3;
        }
      }

      // recorrect the particle fits if pions and electrons switched
      /*if(m_currentParticleTags[0] !=0 && m_currentParticleTags[3] != 0){
        int pionIndex = getCurrParamIndex(0);
        int electronIndex = getCurrParamIndex(3);

        double orderCheck = (m_bichselCurves[0]->Eval(m_currentMomOfBin) - m_bichselCurves[3]->Eval(m_currentMomOfBin))*(m_currentParameters[pionIndex+1] - m_currentParameters[electronIndex+1]);
        if(orderCheck < 0.0){
          if(m_currentPlusMinusIndex == 0){
            m_dEdxFitDataPlus[m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[electronIndex]);
            m_dEdxFitDataPlus[m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[electronIndex+1]);
            m_dEdxFitDataPlus[m_currentPartIndex][0][2]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[electronIndex+2]);
            m_dEdxFitDataPlus[m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, m_currentParameters[pionIndex]);
            m_dEdxFitDataPlus[m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, m_currentParameters[pionIndex+1]);
            m_dEdxFitDataPlus[m_currentPartIndex][3][2]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, m_currentParameters[pionIndex+2]);
          }else{
            m_dEdxFitDataMinus[m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[electronIndex]);
            m_dEdxFitDataMinus[m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[electronIndex+1]);
            m_dEdxFitDataMinus[m_currentPartIndex][0][2]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[electronIndex+2]);
            m_dEdxFitDataMinus[m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, m_currentParameters[pionIndex]);
            m_dEdxFitDataMinus[m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, m_currentParameters[pionIndex+1]);
            m_dEdxFitDataMinus[m_currentPartIndex][3][2]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, m_currentParameters[pionIndex+2]);
          }

          fitResult->SetParName(pionIndex,Form("N_{%s}",m_partInfo->GetParticleSymbol(3,0).Data()));
          fitResult->SetParName(pionIndex+1,Form("#mu_{%s}",m_partInfo->GetParticleSymbol(3,0).Data()));
          fitResult->SetParName(pionIndex+2,Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(3,0).Data()));

          fitResult->SetParName(electronIndex,Form("N_{%s}",m_partInfo->GetParticleSymbol(0,0).Data()));
          fitResult->SetParName(electronIndex+1,Form("#mu_{%s}",m_partInfo->GetParticleSymbol(0,0).Data()));
          fitResult->SetParName(electronIndex+2,Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(0,0).Data()));

        }//else{
        //  cout << "ERROR: Unable to recorrect the switching of electrons and pions in fits!!!!!!!!!!" << endl;
        //}
      }*/

      printCurrentParams();

      //save fit
      if(a_saveImages){
        fittingCanvas->cd();
        //gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetRightMargin(0.45);

        m_currentHistoToFit->Draw("E");

        if(m_drawInitialSeedsToFits && fitFunctInitial) fitFunctInitial->Draw("same");
        if(fitResult) fitResult->SetLineColor(kRed);
        if(fitResult) fitResult->Draw("same");
        vector<TF1*> singleFuncts = drawIndividualFits_Gaus();
        DrawPredictionLines_dEdx(fittingCanvas);
        gPad->Update();
        gSystem->ProcessEvents();
        fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/RapBin_%02d_mTm0Bin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
				     m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin));

        if(m_saveNoLogImages){
          gPad->SetLogy(false);
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits/NoLog_RapBin_%02d_mTm0Bin_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),
             m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), m_currentRapBin, m_currentMtM0Bin));
        }

        for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
        for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
        m_TLinePtrs.clear();

      }


      //save 2D Parameter Files

    delete m_currentHistoToFit;
      //delete fitResult;
    delete fitFunctInitial;
    delete fitFunct;

    }//end pt loop
  }// end eta loop
  } // end plus minus
  delete fittingCanvas;
  return;
} //end of fitZTPC



//============================================================================================================
//============================================================================================================



void ZFitter::fitTPC_NonMixing(int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING ZTPC YIELD ALL CENTRALITIES - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;

  if(a_particleIndex < 0 || a_particleIndex == 3  || a_particleIndex > 7 ){
    cout << " Error in ZFitter::fitTPC_NonMixing  with invalid argument, a_particleIndex" << endl;
    return;
  }
  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = -1;
  m_inAllCentMode = true;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if((m_currentPartCharge == 1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]) || (m_currentPartCharge == -1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]) ){
    cout << " Error in ZFitter::fitZTPC_NonMixing, can't find m_dEdxFitData*****[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  }
  TH2D* base_histogram = dynamic_cast<TH2D*>(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->Clone());
  HistogramUtilities::zeroBins2D(base_histogram);
  base_histogram->SetName("BaseHistogram_ForReference");
  base_histogram->SetTitle("BaseHistogram_ForReference;y;m_{T}-m_{0}");

  cout << "  Starting to form Base Histograms" << endl;
  //Background Amplitude
  m_Background_FitData_ZTPC[m_currentPartIndex][0][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetName(Form("%s_TPC_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetName(Form("%s TPC Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][0][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetName(Form("%s_TPC_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetName(Form("%s TPC Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZTPC[m_currentPartIndex][1][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetName(Form("%s_TPC_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetName(Form("%s TPC Background b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][1][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetName(Form("%s_TPC_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetName(Form("%s TPC Background b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZTPC[m_currentPartIndex][2][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetName(Form("%s_TPC_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetName(Form("%s TPC Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][2][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetName(Form("%s_TPC_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetName(Form("%s TPC Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZTPC[m_currentPartIndex][3][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetName(Form("%s_TPC_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetName(Form("%s TPC Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][3][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetName(Form("%s_TPC_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetName(Form("%s TPC Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  delete base_histogram;



  // #####################      FIT    Plus       ######################


  m_currentPartCharge = 1;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsY();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_dEdx_Plus_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;





      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
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


      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> emptyVect;
      vector<double> fitData = fitIsolation_ZTPC_withBackground(emptyVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;

        cout << "   ptr data: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0] << endl;
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->GetName()
               << "       " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
        m_FitStatus_ZTPC[m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_ZTPC[m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );

      }else{
        cout << "     Error:  fitIsolation_ZTPC_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop









  // #####################      FIT    MINUS       ######################


  if(!a_fitMinus) return;
  m_currentPartCharge = -1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_dEdx_Minus_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;




      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
      cout << "   by ztpcRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztpc Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
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



      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;
      if(a_constrainMinusToPlus){
        constrainVect.push_back(-999); // free the amplitude
        constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999); // free the scale of the background
        constrainVect.push_back(-999); // free the slope of the background
        //constrainVect.push_back(m_InvBetaBTOF_Background_FitData[m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999);
        constrainVect.push_back(-999); // free the scale of the linear background
      }
      vector<double> fitData = fitIsolation_ZTPC_withBackground(constrainVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        cout << "   ptr data: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1] << endl;
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->GetName() << "       "
            << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
        m_FitStatus_ZTPC[m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_ZTPC[m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );
      }else{
        cout << "     Error:  fitIsolation_BTOF_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop

}
//higher count cut
void ZFitter::fitTPC_NonMixing_Liggett(int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus, bool a_constrainMinusToPlus){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING ZTPC YIELD ALL CENTRALITIES - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;

  if(a_particleIndex < 0 || a_particleIndex == 3  || a_particleIndex > 7 ){
    cout << " Error in ZFitter::fitTPC_NonMixing_Liggett  with invalid argument, a_particleIndex" << endl;
    return;
  }
  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = -1;
  m_inAllCentMode = true;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if((m_currentPartCharge == 1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]) || (m_currentPartCharge == -1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]) ){
    cout << " Error in ZFitter::fitZTPC_NonMixing_Liggett, can't find m_dEdxFitData*****[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  }
  TH2D* base_histogram = dynamic_cast<TH2D*>(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->Clone());
  HistogramUtilities::zeroBins2D(base_histogram);
  base_histogram->SetName("BaseHistogram_ForReference");
  base_histogram->SetTitle("BaseHistogram_ForReference;y;m_{T}-m_{0}");

  cout << "  Starting to form Base Histograms" << endl;
  //Background Amplitude
  m_Background_FitData_ZTPC[m_currentPartIndex][0][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetName(Form("%s_TPC_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetName(Form("%s TPC Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][0][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetName(Form("%s_TPC_FitData_Bkgd_Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetName(Form("%s TPC Background Amp",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZTPC[m_currentPartIndex][1][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetName(Form("%s_TPC_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetName(Form("%s TPC Background b",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][1][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetName(Form("%s_TPC_FitData_Bkgd_b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetName(Form("%s TPC Background b",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZTPC[m_currentPartIndex][2][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetName(Form("%s_TPC_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetName(Form("%s TPC Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][2][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetName(Form("%s_TPC_FitData_Bkgd_x0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetName(Form("%s TPC Background x_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  m_Background_FitData_ZTPC[m_currentPartIndex][3][0] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetName(Form("%s_TPC_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetName(Form("%s TPC Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][3][1] = dynamic_cast<TH2D*>(base_histogram->Clone());
  m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetName(Form("%s_TPC_FitData_Bkgd_y0",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));
  m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetName(Form("%s TPC Background y_{0}",m_partInfo->GetParticleName(m_currentPartIndex,-1).Data()));

  delete base_histogram;



  // #####################      FIT    Plus       ######################


  m_currentPartCharge = 1;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsY();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_dEdx_Plus_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;

      if (m_currentMomOfBin<3.7) continue;





      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
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
      if(countsInHisto < 5000){
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


      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> emptyVect;
      vector<double> fitData = fitIsolation_ZTPC_withBackground(emptyVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;

        cout << "   ptr data: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0] << endl;
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->GetName()
               << "       " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZTPC[m_currentPartIndex][0][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
        m_FitStatus_ZTPC[m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_ZTPC[m_currentPartIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );

      }else{
        cout << "     Error:  fitIsolation_ZTPC_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop









  // #####################      FIT    MINUS       ######################


  if(!a_fitMinus) return;
  m_currentPartCharge = -1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->GetXaxis()->FindBin( m_currentRapidity );
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]){
        cout << "No m_rap_mTm0_dEdx_Minus_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;




      cout << "###########################################################################" << endl;
      cout << "Working on All Cent ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Minus_allCent[m_currentPartIndex]->ProjectionZ(
        Form("%s_BTOF_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
      cout << "   by ztpcRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztpc Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
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



      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;
      if(a_constrainMinusToPlus){
        constrainVect.push_back(-999); // free the amplitude
        constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999); // free the scale of the background
        constrainVect.push_back(-999); // free the slope of the background
        //constrainVect.push_back(m_InvBetaBTOF_Background_FitData[m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
        constrainVect.push_back(-999);
        constrainVect.push_back(-999); // free the scale of the linear background
      }
      vector<double> fitData = fitIsolation_ZTPC_withBackground(constrainVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        cout << "   ptr data: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1] << endl;
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
        m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[10]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
        //m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[11]);
        cout << " filling histo named: " << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->GetName() << "       "
            << m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]->GetTitle() << endl;
        cout << "here now doing backgournd data" << endl;
        m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[4]);
        m_Background_FitData_ZTPC[m_currentPartIndex][0][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[12]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[5]);
        m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[13]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[6]);
        m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[14]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetBinContent(m_currentRapBin, m_currentMtM0Bin, fitData[7]);
        m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->SetBinError(m_currentRapBin, m_currentMtM0Bin, fitData[15]);
        m_FitStatus_ZTPC[m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_ZTPC[m_currentPartIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );
      }else{
        cout << "     Error:  fitIsolation_BTOF_withBackground : FAILED for Plus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop

}


//============================================================================================================
//============================================================================================================



void ZFitter::setFXTTPCMeans(double *FXTMeanOverwriteMinMom, double *FXTMeanOverwriteMaxMom){
  //  if(m_currentPartIndex  > 1) return;
  cout << "starting ZFitter::setFXTTPCMeans" << endl;
  m_inAllCentMode = false;
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){

    int numHorBins = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetNbinsX();
    int numVertBins = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetNbinsY();
    for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
      float currentRapVal = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
      double MinMtm0 = PhysMath::mTm0_from_y_p(currentRapVal,FXTMeanOverwriteMinMom[partIndex],m_partInfo->GetParticleMass(m_currentPartIndex));
      double MaxMtm0 = PhysMath::mTm0_from_y_p(currentRapVal,FXTMeanOverwriteMaxMom[partIndex],m_partInfo->GetParticleMass(m_currentPartIndex));

      double previousBinMomentum = 0.0;
      /*
      for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
        //if(m_testing && VertBinIndex > 12) continue;
        m_currentEtaBin = HorBinIndex;
        m_currentPtBin = VertBinIndex;
        cout << "Working on " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << endl;

        m_currentRapidity = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(m_currentEtaBin);
        m_currentMtM0 = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentPtBin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
        //m_currentHistoToFit = m_rap_mTm0_dEdx_allCent[m_currentPartIndex][0]->ProjectionZ("projectionZ", HorBinIndex, HorBinIndex, binVertIndex, binVertIndex);;
        //if(m_ZTPCStops < momOfBin) continue;
        if(partIndex == 0 && m_currentMtM0 > 1.0){
          if(m_dEdxFitDataMinus[m_currentPartIndex][partIndex][0]->GetBinContent(HorBinIndex,VertBinIndex) > 0.06){
            m_dEdxFitDataMinus[m_currentPartIndex][partIndex][0]->SetBinContent(HorBinIndex,VertBinIndex,m_dEdxFitDataMinus[m_currentPartIndex][partIndex][0]->GetBinContent(HorBinIndex,VertBinIndex-1));
          }
        }
      }

      TF1 *pimMeanFit;
      if(m_currentPartIndex == 0){
        pimMeanFit = new TF1("pimMeanFit","[0]",MinMtm0,MaxMtm0);
        pimMeanFit->SetParameter(0,0.03);
      }else if(m_currentPartIndex == 1){
        pimMeanFit = new TF1("pimMeanFit","[0]*exp([1]*x)+[2]",MinMtm0,MaxMtm0);
        pimMeanFit->SetParameters(-0.4,-1,0.2);
      }
      m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->Fit(pimMeanFit,"R");
      */

      //if(m_testing && ( HorBinIndex%(numHorBins/m_divTesting) != 0 || currentRapVal == m_cuts->getYCM())) continue;
      bool doThisRapBin = false;
      if(m_testing){
        if(!(HorBinIndex%(numHorBins/m_divTesting) == 0 || fabs(currentRapVal - m_cuts->getYCM()) < 0.01 )  ){
          doThisRapBin = false;
        }else{
          doThisRapBin = true;
        }

        for(unsigned int iii = 0; iii < m_rapBinsToTest.size(); iii++){
          if(HorBinIndex == m_rapBinsToTest[iii]) doThisRapBin = true;
        }

      }
      if(m_testing && !doThisRapBin) continue;

      previousBinMomentum = 0.0;
      for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
        //if(m_testing && VertBinIndex > 12) continue;
        m_currentRapBin = HorBinIndex;
        m_currentMtM0Bin = VertBinIndex;
        cout << "Working on " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << endl;

        m_currentRapidity = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(m_currentRapBin);
        m_currentMtM0 = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
        //m_currentHistoToFit = m_rap_mTm0_dEdx_allCent[m_currentPartIndex][0]->ProjectionZ("projectionZ", HorBinIndex, HorBinIndex, binVertIndex, binVertIndex);;
        //if(m_ZTPCStops < momOfBin) continue;
        if(partIndex == 0 && m_currentMtM0 > 1.0){
          if(m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(HorBinIndex,VertBinIndex) > 0.06){
            m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->SetBinContent(HorBinIndex,VertBinIndex,m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][0][1]->GetBinContent(HorBinIndex,VertBinIndex-1));
          }
        }
        if(partIndex == 0 && m_currentMomOfBin > FXTMeanOverwriteMinMom[partIndex] && m_currentMomOfBin < FXTMeanOverwriteMaxMom[partIndex]){
          //m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->SetBinContent(HorBinIndex,VertBinIndex,pimMeanFit->Eval(m_currentMtM0));
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->SetBinContent(HorBinIndex,VertBinIndex,m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(HorBinIndex,VertBinIndex));
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->SetBinContent(HorBinIndex,VertBinIndex,m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(HorBinIndex,VertBinIndex));
        }else if(partIndex == 1 && m_currentMomOfBin > FXTMeanOverwriteMinMom[partIndex] && m_currentMomOfBin < FXTMeanOverwriteMaxMom[partIndex]){
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->SetBinContent(HorBinIndex,VertBinIndex,m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(HorBinIndex,VertBinIndex));
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->SetBinContent(HorBinIndex,VertBinIndex,m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(HorBinIndex,VertBinIndex));
        }else if(partIndex == 1 && m_currentMomOfBin < 0.5){
          m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->SetBinContent(HorBinIndex,VertBinIndex,m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(HorBinIndex,VertBinIndex));
        }
      }//end mTm0 loop
    }//end rapidity loop
   }//end particle loop
}//end FXT





//============================================================================================================
//============================================================================================================

void ZFitter::fitTPCCent(int a_centIndex, bool a_saveImages){
  cout << "  Starting to fit TPC by CENTRALITY for Cent Index: " << a_centIndex << endl;
  m_inAllCentMode = false;
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);


  #ifdef _ZFITTER_DEBUG_
    cout << " using ptr: " << m_rap_mTm0_dEdx_allCent[m_currentPartIndex] << " for bining information " << endl;
  #endif

  int numHorBins = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetNbinsY();
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    if(a_centIndex !=-1 && centIndex != a_centIndex) continue;
    m_currentCentIndex = centIndex;

    if(!m_rap_mTm0_dEdx_Plus[m_currentPartIndex][centIndex]){
      cout << "No Histo Plus for cent: " << centIndex << endl;
      cout << "Continuing to next Cent" << endl;
      continue;
    }
    if(!m_rap_mTm0_dEdx_Minus[m_currentPartIndex][centIndex]){
      cout << "No Histo Minus for cent : " << centIndex << endl;
      cout << "Continuing to next Cent" << endl;
      continue;
    }


    for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
      m_currentRapidity = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
      m_currentRapBin = HorBinIndex;

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


      //if(m_testing && ( HorBinIndex%(numHorBins/m_divTesting) != 0 ? currentRapVal != m_cuts->getYCM() : 0 )) continue;
      m_electronExtrapolationComplete = false;

      if(m_useElectronExtrapolation && m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]){
        //zero the amplitude bins because a new electron extrap will take it's place
        #ifdef _ZFITTER_DEBUG_
          cout << " using ptr: " << m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0] << " to zero electron information " << endl;
        #endif
        for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
          m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]->SetBinContent(m_currentRapBin,VertBinIndex, 0.0);
        }
      }

      double previousBinMomentum = 0.0;
      for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
        if(m_testing && VertBinIndex%m_divTestingMtM0!=0) continue;
        m_currentMtM0Bin = VertBinIndex;
        m_currentMtM0 = m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
        m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
        //m_currentHistoToFit = m_rap_mTm0_dEdx_allCent[m_currentPartIndex][0]->ProjectionZ("projectionZ", HorBinIndex, HorBinIndex, binVertIndex, binVertIndex);;
        //        if(m_ZTPCStops < momOfBin) continue;

        if(m_useElectronExtrapolation && m_Fit_Data_ZTPC[m_currentPartIndex][3][0][0]){
          if(m_currentMomOfBin > m_electronExtrapolation_HighMom && previousBinMomentum < m_electronExtrapolation_HighMom){
            makeElectronExtrapolationFunction();
          }
          previousBinMomentum = m_currentMomOfBin;
        }





        //#### PLUS MINUS HISTO LOOP #########
        for(int pmIndex= 1; pmIndex >= 0 ; pmIndex--){

          TH3D* histo3DPM = NULL;
          string partChargeNames[2] = {"Positive","Negative"};
          string partChargeNamesSmall[2] = {"Plus","Minus"};
          if(pmIndex == 0){
            histo3DPM = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][centIndex];
            if(!histo3DPM) continue;
            m_currentPartCharge = 1;
          }else{
            histo3DPM = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][centIndex];
            if(!histo3DPM) continue;
            m_currentPartCharge = -1;
          }
          m_currentPlusMinusIndex = pmIndex;
          cout << "\n\n\n###########################################################################" << endl;
          cout << "Working on Specific Centrality ZTPC " << m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data() << " Cent: " << centIndex << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << endl;
          //################ SETUP THE HISTOGRAM TO FIT ################
          m_currentHistoToFit = histo3DPM->ProjectionZ(
             Form("dEdx%s_cent%02d_%02d_%02d",m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(),centIndex,m_currentRapBin,m_currentMtM0Bin),
             m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
          HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
          double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
          if(countsInHisto < 30){
            #ifdef _ZFITTER_DEBUG_
              cout << "Skipping this bin because of low counts: " << countsInHisto << endl;
            #endif
            delete m_currentHistoToFit;
            continue;
          }
          m_currentIntegralofHistogram = countsInHisto;

          //TH1D* tempHistoPtr = m_currentHistoToFit;
          //m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,5.0);
          //delete tempHistoPtr;
          m_currentHistoToFit->Scale(1.0,"width");
          m_currentHistoToFit->SetTitle(Form("TPC %s Particles for y=[%1.2f,%1.2f] m_{T}-m_{0}=[%2.3f,%2.3f] p=%2.3f",
            partChargeNames[pmIndex].c_str(),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
            m_currentMomOfBin
            ));
          m_currentHistoToFit->GetXaxis()->SetTitle(Form("Z_{TPC} of %s",m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data()));
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/d(Z_{TPC})");


          double lowRangeData = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistoToFit));
          double highRangeData = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistoToFit)+1);
          double lowRange = 0;
          double highRange = 0;


          #ifdef _ZFITTER_DEBUG_
            cout << "Histogram Info Before Fit: " << endl;
            HistogramUtilities::printHistoInfo(m_currentHistoToFit);
          #endif

          //#################################
          //##### set up particle tags ######
          //#################################
          if(m_cuts->getEventConfig() == "FixedTarget"){
            //##########################################################
            //##########    FIXED TARGET PARTICLE TAGS   ############
            //##########################################################
            for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
              cout << "partIndex: " << partIndex << " particleToggle: " << m_initialParticleToggles[partIndex] << endl;
              if((m_currentPartCharge == 1 && (!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex] || m_currentMomOfBin > m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex]) )
              || (m_currentPartCharge == -1 && (!m_partToFitMatrix_Minus[m_currentPartIndex][partIndex] || m_currentMomOfBin > m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex]) )){
                m_currentParticleTags[partIndex] = 0;
                continue;
            }
           /*
            //if(!m_initialParticleToggles[partIndex]) m_currentParticleTags[partIndex] = 0;
            double TPCStop = 0.0;
            if(m_currentPartCharge == 1) TPCStop = m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex];
            if(m_currentPartCharge == -1) TPCStop = m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex];
            if(TPCStop < m_currentMomOfBin){
              m_currentParticleTags[partIndex] = 0;
              continue;
            }*/

            if(partIndex == 3 && m_electronExtrapolationComplete){
                m_currentParticleTags[partIndex] = 3; // the function will constrain the amplitude
              }else if( (m_currentPartCharge == 1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0 )
                     || (m_currentPartCharge == -1 && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) != 0.0 ) ){
                m_currentParticleTags[partIndex] = 1;
              }else{
                m_currentParticleTags[partIndex] = 0;
              }



        }//end particle loop


      #ifdef _ZFITTER_DEBUG_
          for(int iii = 0; iii < 8; iii++){
            if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][iii])
              || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][iii]) ) continue;
            cout << "   " << m_partInfo->GetParticleName(iii) << "  has initial fit tag: " << m_currentParticleTags[iii] << " : " << m_particleTagDescriptions[m_currentParticleTags[iii]] << endl;
          }
      #endif

          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            if(m_currentParticleTags[partIndex] == 0) continue;
            if(partIndex == 3) continue;
            //this particle will not have any free parameters and will mimic the particle defined in m_particleMimicIndicies
            //cout << "Check: PartIndex " << partIndex << " with tag " <<  m_currentParticleTags[partIndex] << endl;
            //cout << "      mom: " << m_currentMomOfBin << " > " << m_particleMimicLowMom[partIndex] << endl;
            // cout << "     mom: " << m_currentMomOfBin << " < " << m_particleMimicHighMom[partIndex] << endl;
            // cout << "     mimc Particle: " << m_particleMimicIndicies[partIndex] << " with tag " << m_currentParticleTags[m_particleMimicIndicies[partIndex]] << endl;
            if( m_currentParticleTags[m_particleMimicIndicies[partIndex]] == 0 || m_particleMimicIndicies[partIndex] == partIndex){
              m_currentParticleTags[partIndex] = 1;
            }else if( m_currentParticleTags[partIndex] == 1 ){
              double partMean;
              double partMimicMean;
              if(m_currentPartCharge == 1){
                partMean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
                partMimicMean = m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              }else{
                partMean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
                partMimicMean = m_Fit_Data_ZTPC[m_currentPartIndex][m_particleMimicIndicies[partIndex]][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
              }
              if( hasTOFRatioToConstrainDeDx(partIndex) && fabs(partMean - partMimicMean) < 0.25 ){
                m_currentParticleTags[partIndex] = 12;
                if(!m_allowSpecCentFitsToVary) m_currentParticleTags[partIndex] = 13;
              }else{
                m_currentParticleTags[partIndex] = 10;
                if(!m_allowSpecCentFitsToVary) m_currentParticleTags[partIndex] = 1;
              }
            }
            /*
              if(hasTOFRatioToConstrainDeDx(partIndex)){
                m_currentParticleTags[partIndex] = 12;
              }else{
                m_currentParticleTags[partIndex] = 10;
              }
            }
            */

          }//end particle loop

      #ifdef _ZFITTER_DEBUG_
          for(int iii = 0; iii < 8; iii++){
            if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][iii])
              || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][iii]) ) continue;
            cout << "   " << m_partInfo->GetParticleName(iii) << "  has fit tag: " << m_currentParticleTags[iii] << " : " << m_particleTagDescriptions[m_currentParticleTags[iii]] << endl;
          cout << "mean weighted average? " << m_weightedAverage_mean[iii] << endl;
          cout << "stdDev weighted average? " << m_weightedAverage_stdDev[iii] << endl;
          }
      #endif
          }else if(m_cuts->getEventConfig() == "ColliderCenter"){















          //##########################################################
          //##########    COLLIDER CENTER PARTICLE TAGS   ############
          //##########################################################
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            m_currentParticleTags[partIndex] = 0;
            if((m_currentPartCharge == 1 && (!m_partToFitMatrix_Plus[m_currentPartIndex][partIndex] || m_currentMomOfBin > m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][partIndex]) )
              || (m_currentPartCharge == -1 && (!m_partToFitMatrix_Minus[m_currentPartIndex][partIndex] || m_currentMomOfBin > m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][partIndex]) )) continue;

            //Isolate the particle of interest if it's all alone in dE/dx space with little or no overlap
            if(m_currentMomOfBin < m_dEdxSeparatedStops[m_currentPartIndex][m_currentPlusMinusIndex] && m_currentPartIndex != 0
              && m_currentPartIndex != partIndex) continue;

            if(m_usePearsonIVinSingle && m_currentMomOfBin < m_PearsonZTPCMaxMomStop[m_currentPartIndex][m_currentPlusMinusIndex] && m_currentPartIndex == partIndex){
              #ifdef _ZFITTER_DEBUG_
               cout << " Fitting with PearsonIV with mom:" << m_currentMomOfBin << " with partIndex=" << partIndex << endl;
               cout << "Histogram Pointer: " << m_currentHistoToFit << endl;
              #endif

              bool zoomedFitSuccess = false;
              vector<double> fitData = fitIsolation_ZTPC_Pearson();
                #ifdef _ZFITTER_DEBUG_
                  cout << "    Fit Parameters returned: " << endl;
                  for(int jjj = 0; jjj < 5; jjj++){
                    cout << "        Param " << jjj << " = " << fitData[jjj] << " +/- " << fitData[jjj + 5] << endl;
                  }
                  cout << "       with status " << fitData[10] << endl;
                #endif
              if(fitData[0] > 0){
                /*m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[0]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[1]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[2]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[3]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[4]);

                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[5]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[6]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[7]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][3][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[8]);
                m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][4][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, fitData[9]);
                */


                //push the results to the spectra
                m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]);
                m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[5]);
                if(m_currentPartCharge > 0){
                  m_dEdxYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]);
                  m_dEdxYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[5]);
                }else{
                  m_dEdxYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]);
                  m_dEdxYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[5]);
                }


                zoomedFitSuccess = true;
                m_currentParticleTags[partIndex] = 0;

                continue;

              }
              if(!zoomedFitSuccess){
                /*if( m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0 ){
                  double mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
                  if(mean > lowRangeData && mean < highRangeData){
                    m_currentParticleTags[partIndex] = 1; // the function will constrain mu and sig
                    continue;
                  }
                }*/
                m_currentParticleTags[partIndex] = 0;
                continue;
              }
            }


            if(partIndex == 3){ // Electron
              if(m_electronExtrapolationComplete){
                 m_currentParticleTags[partIndex] = 3; // the function will constrain the amplitude
              }else if( m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0 ){
                double mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
                if(mean > lowRangeData && mean < highRangeData){
                  m_currentParticleTags[partIndex] = 1; // the function will fix mu and sig
                }
              }

            }else{ // NOT Electron
              if( m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin) > 0.0 ){
                double mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);

                if(m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin) == 0.0){
                  if(mean > lowRangeData && mean < highRangeData){
                    m_currentParticleTags[partIndex] = 1; // the function will constrain mu and sig
                  }
                }else{
                  if(mean > m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)
                     && mean < m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)){
                      m_currentParticleTags[partIndex] = 1; // the function will constrain mu and sig
                  }
                }


              }else{
                m_currentParticleTags[partIndex] = 0; // should have ZTPC Fit already, just doing amplitudes in this function
              }
            }

            if(m_allowSpecCentFitsToVary && m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][2][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) > 0.0){
              if( (m_currentPartCharge == 1 && m_currentMomOfBin < m_partDeDx_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
                || (m_currentPartCharge == 0 && m_currentMomOfBin < m_partDeDx_MomStartExtrapolation_Matrix_Plus[m_currentPartIndex][partIndex])
                || (m_currentPartCharge == -1 && m_currentMomOfBin < m_partDeDx_MomStartExtrapolation_Matrix_Minus[m_currentPartIndex][partIndex]) ){
                double mean = m_Fit_Data_ZTPC[m_currentPartIndex][partIndex][1][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
                if(mean > lowRangeData && mean < highRangeData){
                  m_currentParticleTags[partIndex] = 9;
                }
                continue;
              }
            }


          }//end particle loop


            /*for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
              if(m_currentParticleTags[partIndex] == 0) continue;
              //this particle will not have any free parameters and will mimic the particle defined in m_particleMimicIndicies
              //cout << "Check: PartIndex " << partIndex << " with tag " <<  m_currentParticleTags[partIndex] << endl;
              //cout << "      mom: " << m_currentMomOfBin << " > " << m_particleMimicLowMom[partIndex] << endl;
              // cout << "     mom: " << m_currentMomOfBin << " < " << m_particleMimicHighMom[partIndex] << endl;
              // cout << "     mimc Particle: " << m_particleMimicIndicies[partIndex] << " with tag " << m_currentParticleTags[m_particleMimicIndicies[partIndex]] << endl;
              if( m_particleMimicLowMom[partIndex] < m_currentMomOfBin
                  && m_particleMimicHighMom[partIndex] > m_currentMomOfBin
                  && m_currentParticleTags[m_particleMimicIndicies[partIndex]] != 0){
                if(hasTOFRatioToConstrainDeDx(partIndex)){
                  m_currentParticleTags[partIndex] = 4;
                }else{
                  m_currentParticleTags[partIndex] = 1;
                }
              }
            }*/
          } // end Particle Tags Collider Center
          m_currentNumPartToFit = 0;
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            if(m_currentParticleTags[partIndex] != 0) m_currentNumPartToFit = m_currentNumPartToFit + 1;
          }



          #ifdef _ZFITTER_DEBUG_
          for(int iii = 0; iii < 8; iii++){
            if((m_currentPartCharge == 1 && !m_partToFitMatrix_Plus[m_currentPartIndex][iii])
              || (m_currentPartCharge == -1 && !m_partToFitMatrix_Minus[m_currentPartIndex][iii]) ) continue;
            cout << "   " << m_partInfo->GetParticleName(iii) << "  has fit tag: " << m_currentParticleTags[iii] << " : " << m_particleTagDescriptions[m_currentParticleTags[iii]] << endl;
          }
          #endif

          if(m_currentParticleTags[m_currentPartIndex] == 0){
            #ifdef _ZFITTER_DEBUG_
              cout << "Skipping this bin becasue main particle will not be fit (tag = 0)" << endl;
            #endif
            continue;
          }





         //###########     SET THE HISTOGRAM RANGE TO FIT    ###############
         if(m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinError(m_currentRapBin,m_currentMtM0Bin) == 0.0){
           // there wasn't a previous range yet...
           getRangeZTPC(lowRange, highRange);
           if(lowRange >= highRange){
             lowRange = lowRangeData;
             highRange = highRangeData;
           }

           if(m_currentNumPartToFit == 1){
             double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(-0.068), 2);
             double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.0), 2);
             double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.068), 2);
             double N_3p, mean_3p, sig_3p = 0.0;
             PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,-0.068,y1,0.0,y2,0.068,y3);
             //cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
             if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
               if(m_currentPartIndex == 0  && m_currentMomOfBin < 0.2){ // fit right side of peak
                 lowRange  = mean_3p  - 0.5*sig_3p;
                 highRange  = mean_3p + 3.0*sig_3p;
               }else{
                 lowRange  = mean_3p  - 3.0*sig_3p;
                 highRange  = mean_3p + 3.0*sig_3p;
               }
               //cout << "     After: " << lowRange << " to " << highRange << endl;
             }else{
               //cout << "   3pt gaus failed, just not zooming in" << endl;
             }
           }



           //cout << " m_useRLBiasedSingleParticleZoom: " << m_useRLBiasedSingleParticleZoom << "   m_currentPartIndex: " << m_currentPartIndex << endl;
           if(m_useRLBiasedSingleParticleZoom && m_currentPartIndex == 0){
             double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(-0.068), 2);
             double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.0), 2);
             double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.068), 2);
             double N_3p, mean_3p, sig_3p = 0.0;
             PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,-0.068,y1,0.0,y2,0.068,y3);
             //cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
             if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
               if(m_currentMomOfBin < 0.2){ // fit right side of peak
                 lowRange  = mean_3p  - 0.5*sig_3p;
                 highRange  = mean_3p + 3.0*sig_3p;
               }else{ // fit left side of peak
                 lowRange  = mean_3p  - 3.0*sig_3p;
                 highRange  = mean_3p + 0.5*sig_3p;
               }
               //cout << "     After: " << lowRange << " to " << highRange << endl;
               for(unsigned int jjj = 0; jjj < 9; jjj++){ // kill other particles, just fit one peak
                 if(m_currentPartIndex != jjj) m_currentParticleTags[jjj] = 0;
               }
               m_currentNumPartToFit = 1;
             }else{
               //cout << "   3pt gaus failed, just not zooming in" << endl;
             }
           }else if(m_currentPartIndex != 0 && m_useRLBiasedSingleParticleZoom){
             double y1  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(-0.068), 2);
             double y2  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.0), 2);
             double y3  = HistogramUtilities::getBinContentAvg(m_currentHistoToFit, m_currentHistoToFit->FindBin(0.068), 2);
             double N_3p, mean_3p, sig_3p = 0.0;
             PhysMath::get3PointGausParams(N_3p, mean_3p, sig_3p,-0.068,y1,0.0,y2,0.068,y3);
             //cout << "3pt Range Calc: Before: " << lowRange << " to " << highRange << endl;
             if(sig_3p > 0.0 && N_3p > 0.0 && !std::isnan(mean_3p) && !std::isnan(sig_3p)){
               lowRange  = mean_3p  - 0.5*sig_3p;
               highRange  = mean_3p + 3.0*sig_3p;
               //cout << "     After: " << lowRange << " to " << highRange << endl;
               for(unsigned int jjj = 0; jjj < 9; jjj++){ // kill other particles, just fit one peak
                 if(m_currentPartIndex != jjj) m_currentParticleTags[jjj] = 0;
               }
               m_currentNumPartToFit = 1;
             }else{
               //cout << "   3pt gaus failed, just not zooming in" << endl;
             }
           }

            m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,lowRange);
            m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,1.0);
            m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,highRange);
            m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,1.0);
          }else{
            lowRange = m_ZTPC_LowRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
            highRange = m_ZTPC_HighRange[m_currentPartIndex][m_currentPlusMinusIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
          }

          #ifdef _ZFITTER_DEBUG_
            cout << " Setting Range to [" << lowRange << " , " << highRange << "]" << endl;
          #endif





          TH1D* origHistToFit = m_currentHistoToFit;
          m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRange,highRange);
          delete origHistToFit;
          m_currentHistoToFit->SetTitle(Form("Z_{TPC}  for y=[%1.2f,%1.2f] m_{T}-m_{%s}=[%2.3f,%2.3f] p=%2.3f",
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetXaxis()->GetBinLowEdge(m_currentRapBin+1),
            m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data(),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin),
            m_rap_mTm0_dEdx_allCent[m_currentPartIndex]->GetYaxis()->GetBinLowEdge(m_currentMtM0Bin+1),
            m_currentMomOfBin
            ));
          m_currentHistoToFit->GetXaxis()->SetTitle(Form("Z_{TPC} of %s", m_partInfo->GetParticleSymbol(m_currentPartIndex,m_currentPartCharge).Data()));
          m_currentHistoToFit->GetYaxis()->SetTitle("dN/dZ_{TPC}");

          double truncatedTotalCounts = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
          double possibleNewBins = truncatedTotalCounts/6.0;
          if(possibleNewBins > 80.0) possibleNewBins = 80.0;
          TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(m_currentHistoToFit,(int) possibleNewBins,true,true);
          delete m_currentHistoToFit;
          m_currentHistoToFit = rebinHisto;


          m_currentHistoToFit->SetMarkerStyle(20);
          m_currentHistoToFit->SetMarkerColor(kBlack);



          //initialize parameters
          initZTPCParams();
          cout << "parameter seeds:" << endl;
          printCurrentParams();
          if(m_currentNumPartToFit == 0){
            #ifdef _ZFITTER_DEBUG_
              cout << "Skipping this bin because there are no particles to fit" << endl;
            #endif
            continue;
          }
          //cout << "    " << m_currentNumPartToFit << " Particles" << "   With " << m_currentFreeParameters << " Parameters" << endl;
          //do the fit
          TF1* fitFunct = new TF1("fitFunctCent",this,&ZFitter::ZTPCOperator,lowRange,highRange,m_currentFreeParameters,"ZFitter","ZTPCOperator");
          fitFunct->SetNpx(300);
          initZTPCFunction(fitFunct);
          TF1* fitFunctInitial = (TF1*) fitFunct->Clone();
          fitFunctInitial->SetLineColor(kYellow);
          fitFunctInitial->SetName("fitFunctInitialPlus");
          fitFunctInitial->Draw("same");

          string fitOptionsZTPCCent = "EQ";
          #ifdef _ZFITTER_DEBUG_
            fitOptionsZTPCCent = "EV";
          #endif
           if(m_doIntegralOfBinForFits) fitOptionsZTPCCent += "I";

          m_currentIntegralofHistogram = HistogramUtilities::sumOfContent(m_currentHistoToFit,-1,true);
          int status = m_currentHistoToFit->Fit("fitFunctCent",fitOptionsZTPCCent.c_str(),"",lowRange,highRange);//add Q to quiet
          TF1* fitResult = m_currentHistoToFit->GetFunction("fitFunctCent");

          if(m_currentPartCharge > 0){
            m_FitStatus_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status + 5);
          }else{
            m_FitStatus_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,status + 5);
          }
          #ifdef _ZFITTER_DEBUG_
            cout << " Final Status: " << status << "     fitResultPtr: " << fitResult << endl;
          #endif

          if(status!=0){
            cout << "\n\n ABNORMAL FIT: status=" << status << "\n\n";

            if(a_saveImages){
              fittingCanvas->cd();
              m_currentHistoToFit->SetStats(true);
              m_currentHistoToFit->Draw("E");
              m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange, highRange);
              gPad->SetLogy();
              gPad->SetRightMargin(0.4);
              if(m_drawInitialSeedsToFits) fitFunctInitial->Draw("same");
              if(fitResult){
                fitResult->SetLineColor(kRed);
                fitResult->Draw("same");
              }
              vector<TF1*> singleFuncts = drawIndividualFits_Gaus();
              DrawPredictionLines_dEdx(fittingCanvas);
              gPad->Update();
              gSystem->ProcessEvents();
              fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_%s____BAD_FIT.png",
                m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), centIndex,
                m_currentRapBin, m_currentMtM0Bin, partChargeNamesSmall[pmIndex].c_str()));

              if(m_saveNoLogImages){
              gPad->SetLogy(false);
              gPad->Update();
              gSystem->ProcessEvents();
              fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_%s____BAD_FIT.png",
                m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), centIndex,
                m_currentRapBin, m_currentMtM0Bin, partChargeNamesSmall[pmIndex].c_str()));
              }
              for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
              for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
              m_TLinePtrs.clear();
            }
            if(!m_storeFitsNotStatusZero) continue;
          }

          if(!fitResult){
           cout << " NO FIT RESULT: See Parameters Below" << endl;
           printCurrentParams();
           continue;
         }
          //vector< double > errorsAmp;

          //int paramIndexCounter = 0;
          //for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
          //  if(m_paramToggles[iii] && fitResult){
          //    m_currentParameters[iii] = fitResult->GetParameter(paramIndexCounter);
          //    //errorsAmp.push_back(fitResult->GetParError(paramIndexCounter));
          //    paramIndexCounter++;
          //  }
          //}
          //
          vector< double > allParamError; // even not fit
          int paramIndexCounter = 0;
          int liveParticleCounter = 0;
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){

            if(m_currentParticleTags[partIndex] == 1){
              m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              allParamError.push_back(fitResult->GetParError(paramIndexCounter)*m_currentIntegralofHistogram);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              paramIndexCounter++;
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 2){
              m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              m_currentParameters[3*liveParticleCounter+1] = fitResult->GetParameter(paramIndexCounter+1);
              m_currentParameters[3*liveParticleCounter+2] = fitResult->GetParameter(paramIndexCounter+2);
              allParamError.push_back(fitResult->GetParError(paramIndexCounter)*m_currentIntegralofHistogram);
              allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
              allParamError.push_back(fitResult->GetParError(paramIndexCounter+2));
              paramIndexCounter+=3;
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 3){
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              liveParticleCounter++;
              //no free parameters
            }else if(m_currentParticleTags[partIndex] == 4){
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              liveParticleCounter++;

              // no free paramters
            }else if(m_currentParticleTags[partIndex] == 5){
              m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              allParamError.push_back(fitResult->GetParError(paramIndexCounter)*m_currentIntegralofHistogram);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              liveParticleCounter++;
              paramIndexCounter++;
              // no free paramters
            }else if(m_currentParticleTags[partIndex] == 6){
              m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              allParamError.push_back(fitResult->GetParError(paramIndexCounter)*m_currentIntegralofHistogram);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              liveParticleCounter++;
              paramIndexCounter++;
              // no free paramters
            }else if(m_currentParticleTags[partIndex] == 7){
              m_currentParameters[3*liveParticleCounter+1] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              m_currentParameters[3*liveParticleCounter+2] = fitResult->GetParameter(paramIndexCounter+1);
              allParamError.push_back(0.0);
              allParamError.push_back(fitResult->GetParError(paramIndexCounter));
              allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
              paramIndexCounter+=2;
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 8){
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 9){
              m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              m_currentParameters[3*liveParticleCounter+2] = fitResult->GetParameter(paramIndexCounter+1);
              allParamError.push_back(fitResult->GetParError(paramIndexCounter)*m_currentIntegralofHistogram);
              allParamError.push_back(0.0);
              allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
              paramIndexCounter+=2;
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 10){
              m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              allParamError.push_back(fitResult->GetParError(paramIndexCounter)*m_currentIntegralofHistogram);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              paramIndexCounter++;
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 11){
              m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentIntegralofHistogram;
              m_currentParameters[3*liveParticleCounter+1] = fitResult->GetParameter(paramIndexCounter+1);
              allParamError.push_back(fitResult->GetParError(paramIndexCounter)*m_currentIntegralofHistogram);
              allParamError.push_back(fitResult->GetParError(paramIndexCounter+1));
              allParamError.push_back(0.0);
              paramIndexCounter+=2;
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 12){
              //ben here
              m_currentParameters[3*liveParticleCounter] = m_currentParameters[3*liveParticleCounter]*m_currentIntegralofHistogram;
              //m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentConstrainedAmpRaios[partIndex]*m_currentIntegralofHistogram;
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              liveParticleCounter++;
            }else if(m_currentParticleTags[partIndex] == 13){
              //ben here
              m_currentParameters[3*liveParticleCounter] = m_currentParameters[3*liveParticleCounter]*m_currentIntegralofHistogram;
              //m_currentParameters[3*liveParticleCounter] = fitResult->GetParameter(paramIndexCounter)*m_currentConstrainedAmpRaios[partIndex]*m_currentIntegralofHistogram;
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              allParamError.push_back(0.0);
              liveParticleCounter++;
            }

          }

          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            if(m_currentParticleTags[partIndex] == 4 || m_currentParticleTags[partIndex] == 7 || m_currentParticleTags[partIndex] == 8 || m_currentParticleTags[partIndex] == 12 || m_currentParticleTags[partIndex] == 13){
              allParamError[getCurrParamIndex(partIndex)] = allParamError[getCurrParamIndex(m_particleMimicIndicies[partIndex])]
              /m_currentParameters[getCurrParamIndex(m_particleMimicIndicies[partIndex])];
              // no free paramters
            }
          }

          cout << "Particle Tags: ";
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++) cout << m_currentParticleTags[partIndex] << "  ";
            cout << endl;

          printCurrentParams();
          cout << "Parameter Errors: " << endl;
          for(int iii = 0; iii < (int) allParamError.size(); iii++) cout << allParamError[iii] << "  ";
            cout << endl;

          //cout << endl;
          if(fitResult){
            double chiSqrNDF = fitResult->GetChisquare() / fitResult->GetNDF();
            m_ChiSqr_Cent_ZbTOF[m_currentPartIndex][centIndex][pmIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, chiSqrNDF);
          }
          //push to particle fits
          int paramIndex = 0; //the particles not included in fit
          for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
            if(m_currentParticleTags[partIndex] !=0){// && errorsAmp.size() != 0){
              if(partIndex == m_currentPartIndex){
                m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  m_currentParameters[getCurrParamIndex(m_currentPartIndex)]);
                m_Spectra_Cent_ZTPC[m_currentPartIndex][centIndex][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                 allParamError[getCurrParamIndex(m_currentPartIndex)]);
              }
              if(m_currentPartCharge > 0){
                m_dEdxYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  m_currentParameters[getCurrParamIndex(partIndex)]);
                m_dEdxYieldsCentPlus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                 allParamError[getCurrParamIndex(partIndex)]);
              }else{
                m_dEdxYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,
                  m_currentParameters[getCurrParamIndex(partIndex)]);
                m_dEdxYieldsCentMinus[m_currentPartIndex][partIndex][centIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin,
                 allParamError[getCurrParamIndex(partIndex)]);
              }

              if(m_currentParticleTags[partIndex] !=0){

                m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][centIndex][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex] * m_currentIntegralofHistogram);
                m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][centIndex][1][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+1]);
                m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][centIndex][2][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[paramIndex+2]);
                if(allParamError[paramIndex] > 0.0 ) m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][centIndex][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, allParamError[paramIndex]* m_currentIntegralofHistogram);
                if(allParamError[paramIndex + 1] > 0.0 ) m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][centIndex][1][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, allParamError[paramIndex+1]);
                if(allParamError[paramIndex + 2] > 0.0 ) m_Fit_Data_Cent_ZTPC[m_currentPartIndex][partIndex][centIndex][2][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, allParamError[paramIndex+2]);

                paramIndex+=3;
              }


              // for electron extrapolation fiiting
              if(partIndex == 3 && m_currentParticleTags[partIndex] != 0){
                m_Fit_Data_ZTPC[m_currentPartIndex][3][0][m_currentPlusMinusIndex]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, m_currentParameters[getCurrParamIndex(3)]);
                m_Fit_Data_ZTPC[m_currentPartIndex][3][0][m_currentPlusMinusIndex]->SetBinError(m_currentRapBin,m_currentMtM0Bin, allParamError[getCurrParamIndex(3)]);
              }
            }
          }

          //rescale the amplitude because the single funct drawing requires it
          for(unsigned int paramIndex = 0; paramIndex < m_currentParameters.size(); paramIndex+=3){
              m_currentParameters[paramIndex] = m_currentParameters[paramIndex]/m_currentIntegralofHistogram;
          }

          //save fit
          if(a_saveImages){
            fittingCanvas->cd();
            m_currentHistoToFit->SetStats(true);
            m_currentHistoToFit->Draw("E");
            m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange, highRange);
            gPad->SetLogy();
            gPad->SetRightMargin(0.4);
            if(m_drawInitialSeedsToFits) fitFunctInitial->Draw("same");
            if(fitResult){
              fitResult->SetLineColor(kRed);
              fitResult->Draw("same");

            }
            vector<TF1*> singleFuncts = drawIndividualFits_Gaus();
            DrawPredictionLines_dEdx(fittingCanvas);
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/RapBin_%02d_mTm0Bin_%02d_%s.png",
              m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), centIndex,
              m_currentRapBin, m_currentMtM0Bin, partChargeNamesSmall[pmIndex].c_str()));

            if(m_saveNoLogImages){
            gPad->SetLogy(false);
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/%s/dEdxFits_Cent%02d/NoLog_RapBin_%02d_mTm0Bin_%02d_%s.png",
              m_imagePreDir.c_str(), m_imgDirName.c_str(), m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data(), centIndex,
              m_currentRapBin, m_currentMtM0Bin, partChargeNamesSmall[pmIndex].c_str()));
            }
            for(unsigned int iii = 0; iii < singleFuncts.size(); iii++) delete singleFuncts[iii];
            for(unsigned int iii = 0; iii < m_TLinePtrs.size(); iii++) delete m_TLinePtrs[iii];
            m_TLinePtrs.clear();
          }


        delete m_currentHistoToFit;
          //delete fitResult;
        delete fitFunctInitial;
        delete fitFunct;
        }// end plus minus histo loop


      }//end pt loop
    }
  }// end centrality loop

  delete fittingCanvas;


  return;

};






















//============================================================================================================
//============================================================================================================



void ZFitter::fitTPCCent_NonMixing(int a_centIndex, int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING ZTPC YIELD CENTRALITY BIN " << a_centIndex << " - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;

  if(a_particleIndex < 0 || a_particleIndex == 3  || a_particleIndex > 7 ){
    cout << " Error in ZFitter::fitTPC_NonMixing  with invalid argument, a_particleIndex" << endl;
    return;
  }
  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = a_centIndex;
  m_inAllCentMode = false;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if((m_currentPartCharge == 1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0])
      || (m_currentPartCharge == -1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]) ){
    cout << " Error in ZFitter::fitZTPC_NonMixing, can't find m_dEdxFitData*****[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  }




  // #####################      FIT    Plus       ######################


  m_currentPartCharge = 1;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetNbinsY();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->FindBin( m_currentRapidity );
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]){
        cout << "No m_rap_mTm0_dEdx_Plus for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;

      if(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) <= 0.0) continue;



      cout << "###########################################################################" << endl;
      cout << "Working on Plus Cent" << m_currentCentIndex << " ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->ProjectionZ(
        Form("%s_ZTPC_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
      cout << "   by ztpcRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because ZTPC Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
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
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");



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


      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;
      constrainVect.push_back(-999); // free the amplitude
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(0.0); // the dof - just ignore this

      constrainVect.push_back(-999); // free the scale of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the slope of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the scale of the linear background
      cout << "  Made Constraint Vector " << endl;
      vector<double> fitData = fitIsolation_ZTPC_withBackground(constrainVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]);
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[8]);
        m_FitStatus_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );

      }
    }//end mtm0 loop
  } //end rap loop





  // #####################      FIT    MINUS       ######################
  if(!a_fitMinus) return;
  m_currentPartCharge = -1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->FindBin( m_currentRapidity );
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]){
        cout << "No m_rap_mTm0_dEdx_Minus_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;
      if(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) <= 0.0) continue;



      cout << "###########################################################################" << endl;
      cout << "Working on Minus Cent" << m_currentCentIndex << " ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->ProjectionZ(
        Form("%s_ZTPC_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
      cout << "   by ztpcRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztpc Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
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
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");

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



      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;

      constrainVect.push_back(-999); // free the amplitude
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(-999); // the dof - just ignore this

      constrainVect.push_back(-999); // free the scale of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the slope of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the scale of the linear background


      vector<double> fitData = fitIsolation_ZTPC_withBackground(constrainVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]);
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[8]);
        m_FitStatus_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );
      }else{
        cout << "     Error:  fitIsolation_ZTPC_withBackground : FAILED for Minus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop

}

void ZFitter::fitTPCCent_NonMixing_Liggett(int a_centIndex, int a_particleIndex, double a_momentum_startUsingBkgd, bool a_fitMinus){

  cout << "####################################################################################################" << endl;
  cout << "                        FITTING ZTPC YIELD CENTRALITY BIN " << a_centIndex << " - SIMPLE PARTICLES  " << endl;
  cout << "####################################################################################################" << endl;

  if(a_particleIndex < 0 || a_particleIndex == 3  || a_particleIndex > 7 ){
    cout << " Error in ZFitter::fitTPC_NonMixing_Liggett  with invalid argument, a_particleIndex" << endl;
    return;
  }
  m_currentPartIndex = a_particleIndex;
  m_currentCentIndex = a_centIndex;
  m_inAllCentMode = false;


  //======================     CREATE DATA STORAGE FOR BACKGROUND FUNCTIONS ==========================

  if((m_currentPartCharge == 1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][0])
      || (m_currentPartCharge == -1 && !m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][0][1]) ){
    cout << " Error in ZFitter::fitZTPC_NonMixing_Liggett, can't find m_dEdxFitData*****[m_currentPartIndex][m_currentPartIndex][0]" << endl;
  }




  // #####################      FIT    Plus       ######################


  m_currentPartCharge = 1;
  string partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  int numHorBins = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetNbinsX();
  int numVertBins = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetNbinsY();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->FindBin( m_currentRapidity );
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]){
        cout << "No m_rap_mTm0_dEdx_Plus for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Plus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;

      if(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) <= 0.0) continue;



      cout << "###########################################################################" << endl;
      cout << "Working on Plus Cent" << m_currentCentIndex << " ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Plus[m_currentPartIndex][m_currentCentIndex]->ProjectionZ(
        Form("%s_ZTPC_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
      cout << "   by ztpcRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because ZTPC Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
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
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");



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


      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;
      constrainVect.push_back(-999); // free the amplitude
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(0.0); // the dof - just ignore this

      constrainVect.push_back(-999); // free the scale of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][1][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the slope of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][2][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][3][0]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the scale of the linear background
      cout << "  Made Constraint Vector " << endl;
      vector<double> fitData = fitIsolation_ZTPC_withBackground(constrainVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]);
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[8]);
        m_FitStatus_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );

      }
    }//end mtm0 loop
  } //

  cout << "end rap loop" << endl;





  // #####################      FIT    MINUS       ######################
  if(a_fitMinus){
  m_currentPartCharge = -1;
  partName = m_partInfo->GetParticleName(m_currentPartIndex,m_currentPartCharge).Data();
  for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){

    cout << " Starting HorBin " << HorBinIndex << endl;
    m_currentRapidity = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->GetBinCenter(HorBinIndex);
    //int numDeDxBins = m_rap_mTm0_dEdx_Plus_allCent[m_currentPartIndex]->GetNbinsX();
    int dEdxBin = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->GetXaxis()->FindBin( m_currentRapidity );
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
      m_currentMtM0 = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
      m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(m_currentRapidity, m_currentMtM0, m_partInfo->GetParticleMass(m_currentPartIndex));
      //cout << "Momentum of Bin" << m_currentMomOfBin << endl;
      if(!m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]){
        cout << "No m_rap_mTm0_dEdx_Minus_allCent for particle " << m_currentPartIndex << endl;
        break;
      }

      if(m_partDeDxMomStopMatrix_Minus[m_currentPartIndex][m_currentPartIndex] < m_currentMomOfBin) continue;
      if(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->GetBinContent(m_currentRapBin, m_currentMtM0Bin) <= 0.0) continue;



      cout << "###########################################################################" << endl;
      cout << "Working on Minus Cent" << m_currentCentIndex << " ZTPC " << partName << " RapBin: " << HorBinIndex << " MtM0Bin: " << VertBinIndex << " Mom: " << m_currentMomOfBin << endl;
      //#############    SETTING UP HISTOGRAM TO FIT #################################
      m_currentHistoToFit = m_rap_mTm0_dEdx_Minus[m_currentPartIndex][m_currentCentIndex]->ProjectionZ(
        Form("%s_ZTPC_%02d_%02d",partName.c_str(),m_currentRapBin, m_currentMtM0Bin),
        m_currentRapBin, m_currentRapBin, m_currentMtM0Bin, m_currentMtM0Bin);
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);


      double lowRange = HistogramUtilities::getLowEdge(m_currentHistoToFit);
      double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
      cout << " low: " << lowRange << " high:" << highRange << endl;
      getRangeZTPC_Simple(lowRange, highRange);
      cout << "   by ztpcRange low: " << lowRange << "  high:" << highRange << endl;
      if(highRange < lowRange){
        #ifdef _ZFITTER_DEBUG_
          cout << "Skipping because Ztpc Range isn't in the correct order LowCut: " << lowRange << "   HighCut: " << highRange << endl;
        #endif
        delete m_currentHistoToFit;
        continue;
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
      TH1D* tempHistoPtr = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,4.0);
      delete tempHistoPtr;
      m_currentHistoToFit->Scale(1.0,"width");

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



      int meanBin = m_currentHistoToFit->FindBin(0.0);
      double countsInThatBin = m_currentHistoToFit->GetBinContent(meanBin)*m_currentHistoToFit->GetBinWidth(meanBin);
      cout << " Mean Bin = " << meanBin << "   counts in that bin: " << countsInThatBin << endl;
      cout << " Is Empty ? " << HistogramUtilities::isEmpty(m_currentHistoToFit) << endl;
      if(countsInThatBin < 1.0) continue;

      cout << " Starting the fit " << endl;
      vector<double> constrainVect;

      constrainVect.push_back(-999); // free the amplitude
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Fit_Data_ZTPC[m_currentPartIndex][m_currentPartIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(-999); // the dof - just ignore this

      constrainVect.push_back(-999); // free the scale of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][1][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the slope of the background
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][2][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin));
      constrainVect.push_back(m_Background_FitData_ZTPC[m_currentPartIndex][3][1]->GetBinContent(m_currentRapBin,m_currentMtM0Bin)); // free the scale of the linear background


      vector<double> fitData = fitIsolation_ZTPC_withBackground(constrainVect);
      cout << " out of fit" << endl;
      if(fitData[0] > 0){
        cout << "fit data size " << fitData.size() << endl;
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[0]);
        m_Spectra_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinError(m_currentRapBin,m_currentMtM0Bin,fitData[8]);
        m_FitStatus_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin, fitData[16] + 5);
        m_ChiSqr_Cent_ZTPC[m_currentPartIndex][m_currentCentIndex][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,fitData[17] );
      }else{
        cout << "     Error:  fitIsolation_ZTPC_withBackground : FAILED for Minus" << endl;
      }

    }//end mtm0 loop
  } //end rap loop
}
}
