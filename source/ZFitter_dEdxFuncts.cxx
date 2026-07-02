// THIS IS INSERTED INTO ZFITTER.cxx // Just separate file becuase it's long


//#########################################################################################################################
//                                        GENERAL dE/dx vs Mom FITTING
//#########################################################################################################################

double ZFitter::dEdxFunct(double a_dEdx){
  double val = 0;
  double params[3];
  for(int iii = 0; iii < m_currentNumPartToFit; iii++){
    params[0] = m_currentParameters[3*iii];
    params[1] = TMath::Log(   m_currentParameters[3*iii + 1]/
                TMath::Sqrt(1+m_currentParameters[3*iii + 2]*m_currentParameters[3*iii + 2]/(m_currentParameters[3*iii + 1]*m_currentParameters[3*iii + 1]))  );
    params[2] = TMath::Sqrt( TMath::Log(1+m_currentParameters[3*iii + 2]*m_currentParameters[3*iii + 2]/(m_currentParameters[3*iii + 1]*m_currentParameters[3*iii + 1]))  );
    val += logNormAmp(&params[0], a_dEdx);
  }
  if(std::isnan(val)){
   cout << "dEdx = " << a_dEdx << " Is Nan" << endl;
   printCurrentParams();
   cout << params[0] << " " << params[1] << " " << params[2] << endl;
   cout << " as Current Parameters " << endl;
  }
  return val;
}


//============================================================================================================


double ZFitter::dEdxOperator(double* a_x,double* a_params){
  int paramIndexCounter = 0;
  for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
    if(m_paramToggles[iii]){
      m_currentParameters[iii] = a_params[paramIndexCounter];
      paramIndexCounter++;
    }
  }
  return dEdxFunct(a_x[0]);
}

//============================================================================================================


void ZFitter::initDEDXFunction(TF1* a_funct){
    int paramIndexCounter = 0;
    double integral = m_currentHistoToFit->Integral();
    for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
      if(m_paramToggles[iii]){
         //cout << " Funct Set Par " <<  paramIndexCounter << " to " <<  m_currentParameters[iii] << endl;   
         a_funct->SetParameter(paramIndexCounter, m_currentParameters[iii]);
         a_funct->SetParName(paramIndexCounter, m_currentParamNames[paramIndexCounter].c_str());
         if(iii%3==0){
            if(m_currentParameters[iii] > 1){
              a_funct->SetParLimits(paramIndexCounter,0.0,1.3*integral);// 2.0*m_currentParameters[iii]);
              //cout << "  Limits on param " << paramIndexCounter << " are 0 to " << 1000.0*m_currentParameters[iii] << endl;
            }else{
              a_funct->SetParameter(paramIndexCounter, 0.5*integral);
              a_funct->SetParLimits(paramIndexCounter,0.0, 1.3*integral);
              //cout << "  Limits on param " << paramIndexCounter << " are 0 to 100" << endl;
            }

         }
         if(iii%3==1) a_funct->SetParLimits(paramIndexCounter,m_currentParameters[iii] - m_currentParameters[iii]/20.0, m_currentParameters[iii] + m_currentParameters[iii]/20.0 + 2.0);
         if(iii%3==2) a_funct->SetParLimits(paramIndexCounter, m_currentParameters[iii]/30.0, m_currentParameters[iii-1]/5.0+0.1);
         paramIndexCounter++;
      }
    }
}


//============================================================================================================



void ZFitter::initdEdxMomParams(){
  vector<double> params;
  vector<string> paramNames;
  vector<bool> paramToggles;
  m_currentFreeParameters = 0;
  m_currentNumPartToFit = 0;

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    /*if(m_currentParticleTags[partIndex] == 1){
      //cout << "    Has Solo dE/dx Fit Data" << endl;
      double mean = m_soloMeanHistos[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      double stdDev = m_soloStdDevHistos[m_currentPartIndex][partIndex]->GetBinContent(m_currentRapBin,m_currentMtM0Bin);
      int predictBin = m_currentHistoToFit->FindBin(mean);
      double amp = m_currentHistoToFit->GetBinContent(predictBin) / gaus(mean,stdDev,mean);
      if(std::isnan(amp)){
        cout << "SOMETHING WRONG:: amp Infinite" << endl;
      }
      //cout << "      Amp= " << amp << endl;
      //cout << "      Mean= " << mean << endl;
      //cout << "      Std= " << stdDev << endl;
      params.push_back(amp);
      paramToggles.push_back(true);
      params.push_back(mean);
      paramToggles.push_back(false);
      params.push_back(stdDev);
      paramToggles.push_back(false);
      if(stdDev==0){
        cout << "Error: sig is zero for solo dEdx" << endl;
      }
      m_currentNumPartToFit++;
      m_currentFreeParameters++;
    }*/
    if(m_currentParticleTags[partIndex] == 2){
      //cout << "  Fitting Part: " << partIndex << endl;
      //cout << "    No Solo dE/dx Data" << partIndex << endl;

      double mean = m_bichselCurves[partIndex]->Eval(m_currentMomOfBin) + m_ZTPCshift;
      int dEdxPredictBin = m_currentHistoToFit->FindBin(mean);
      double stdDev = mean*0.08;
      double var = stdDev*stdDev;
      double mu = TMath::Log(   mean/TMath::Sqrt(1+var/(mean*mean))  );
      double sig = TMath::Sqrt( TMath::Log(1+var/(mean*mean))  );
      double ampAtPredict = m_currentHistoToFit->GetBinContent(dEdxPredictBin) / logNorm(mu,sig,mean);
      /*
      cout << "      Amp= " << ampAtPredict << endl;
      cout << "      mean= " << mean << endl;
      cout << "      Mu= " << mu << endl;
      cout << "      stdDev= " << stdDev << endl;
      cout << "      sig= " << sig << endl;
      */
      params.push_back(ampAtPredict);
      paramToggles.push_back(true);
      paramNames.push_back(Form("N_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(mean);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#mu_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      params.push_back(stdDev);
      paramToggles.push_back(true);
      paramNames.push_back(Form("#sigma_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));

      m_currentNumPartToFit++;
      m_currentFreeParameters+=3;


   }
    // if tag = 0 then don't fit



  }//end particle loop
  m_currentParameters.clear();
  m_paramToggles.clear();
  m_currentParamNames.clear();

  m_currentParameters = params;
  m_paramToggles = paramToggles;
  m_currentParamNames = paramNames;
}



//============================================================================================================

void ZFitter::fitDeDx_Solo(int a_rebinMom){
  //make the data structures
  m_dEdxMomFittingCompleted = true;

  gStyle->SetOptFit(0111);

  double fullDeDxStops[9]; //tell when to not fit with full dEdx sepctrum
  fullDeDxStops[0]=0.26;
  fullDeDxStops[1]=0.36;
  fullDeDxStops[2]=0.64;
  fullDeDxStops[3]=0.0;
  fullDeDxStops[4]=0.8;
  fullDeDxStops[5]=0.5;
  fullDeDxStops[6]=0.5;
  fullDeDxStops[7]=0.5;
  fullDeDxStops[8]=0.1;

  TF1* logNormFunct = new TF1("logNormFunct",
      "[0]*0.3989422804 * TMath::Exp(-0.5*( TMath::Log(x) - [1] )*( TMath::Log(x) - [1] )/([2]*[2])) /([2] * x)",
      0,200);
  logNormFunct->SetNpx(500);
  logNormFunct->SetLineColor(kGreen);

  TH2I* dEdxPlus = HistogramUtilities::reBinVar2D(m_mom_dEdx_Plus,a_rebinMom,1);
  TH2I* dEdxMinus = HistogramUtilities::reBinVar2D(m_mom_dEdx_Minus,a_rebinMom,1);

  TCanvas* fittingCanvas = new TCanvas("dEdx_Solo_Fitting");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);
  TH1D* binBase = dEdxPlus->ProjectionX("momBins");
  HistogramUtilities::zeroBins(binBase);
  
  //Make the Fit Result Histograms for Parameter Storage
  string paramNames[3] = {"N","#mu","#sigma"};
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(!m_initialParticleToggles[partIndex]) continue;
    for(int paramIndex = 0; paramIndex < 3; paramIndex++){
      m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][paramIndex]= (TH1D*) binBase->Clone();
      m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][paramIndex]->SetName(Form("%s_%s",m_partInfo->GetParticleName(partIndex,1).Data(), paramNames[paramIndex].c_str()));
      m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][paramIndex]->SetTitle(Form("%s %s; p_{tot} (GeV/c); %s",m_partInfo->GetParticleName(partIndex,1).Data(),paramNames[paramIndex].c_str(), paramNames[paramIndex].c_str()));
      m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][paramIndex]= (TH1D*) binBase->Clone();
      m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][paramIndex]->SetName(Form("%s_%s",m_partInfo->GetParticleName(partIndex,-1).Data(), paramNames[paramIndex].c_str()));
      m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][paramIndex]->SetTitle(Form("%s %s; p_{tot} (GeV/c); %s",m_partInfo->GetParticleName(partIndex,-1).Data(),paramNames[paramIndex].c_str(), paramNames[paramIndex].c_str()));

    }
  }

  //## FIT THE bTOF isolated distributions ##
  for(int histoIndex = 0; histoIndex < 2; histoIndex++){
    string histoNameStr = "plus";
    if(histoIndex == 1){
      histoNameStr = "minus";
    }
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      if(!m_initialParticleToggles[partIndex] || !m_mom_dEdx_BTOFIsolated_Plus[partIndex]) continue;
      if(partIndex == 3 || partIndex == 8) continue;
      TH2I* histoPM = HistogramUtilities::reBinVar2D(m_mom_dEdx_BTOFIsolated_Plus[partIndex],a_rebinMom,1);
      if(histoIndex == 1){
        histoPM = HistogramUtilities::reBinVar2D(m_mom_dEdx_BTOFIsolated_Minus[partIndex],a_rebinMom,1);
      }
      for(int bin = 1;histoPM && bin <= histoPM->GetNbinsX(); bin++){
        TH1D* tempThing = histoPM->ProjectionY(Form("proj_%02d_%02d",histoIndex,bin),bin,bin);
        double momentum = histoPM->GetXaxis()->GetBinCenter(bin);
        double dEdxPredict = m_bichselCurves[partIndex]->Eval(momentum);
        m_currentHistoToFit = HistogramUtilities::truncateHistogram(tempThing,0.7*dEdxPredict,1.3*dEdxPredict);
        delete tempThing;
        //m_currentHistoToFit = tempThing;
        if(!m_currentHistoToFit) continue; // if the truncation doesn't work...
        double counts = HistogramUtilities::sumOfContent(m_currentHistoToFit, -1);
        if( 100 > counts ){  // && !m_testing)  || (counts < 20)){
          delete m_currentHistoToFit;
          cout << "Continuing: too few counts for solo fit partIndex:" << partIndex << " counts:" << counts << " at mom: " << momentum << endl; 
          continue;
        }
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
        //m_currentHistoToFit->Sumw2();
        m_currentHistoToFit->Scale(1.0,"width");     
        
        double integral = m_currentHistoToFit->Integral();
        double mean = m_currentHistoToFit->GetMean();
        double stdDev = m_currentHistoToFit->GetStdDev();
        double mu = 0.0;
        double sig = 0.0;      
        PhysMath::gausToLogNorm(mean,stdDev,mu,sig);
        

        TString imageName = TString::Format("%s_%s_momBin_%03d",
                                        m_partInfo->GetParticleName(partIndex,0).Data(),histoNameStr.c_str(),bin);
        m_currentHistoToFit->SetName(imageName.Data());
        m_currentHistoToFit->SetTitle(Form("TOF Isolated %s %s for p_{tot}=[%2.3f,%2.3f] (GeV/c)",
              m_partInfo->GetParticleName(partIndex,0).Data(),
              histoNameStr.c_str(),
              histoPM->GetXaxis()->GetBinLowEdge(bin),
              histoPM->GetXaxis()->GetBinLowEdge(bin+1)
              ));
        m_currentHistoToFit->GetXaxis()->SetTitle("dE/dx");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/ddEdx");
        m_currentHistoToFit->SetStats(kTRUE);
        //cout << "mean: " << mean << "  stdDev: " << stdDev << "  mu: " << mu << "   sig: " << sig << endl;
        //fitFunction->SetParameter(0,(counts)*m_currentHistoToFitTemp->GetBinWidth(2));
        logNormFunct->SetParameter(0,counts);
        logNormFunct->SetParLimits(0,0.6*counts, 1.5*counts);
        logNormFunct->SetParameter(1,mu);
        logNormFunct->SetParLimits(1, 0.5*mu, 2.0*mu);
        logNormFunct->SetParameter(2,sig);
        logNormFunct->SetParLimits(2,0.05*sig,2.0*sig);
        logNormFunct->SetParNames("N","#mu","#sigma");

        TF1* initialFunct = (TF1*) logNormFunct->Clone();
        initialFunct->SetLineColor(kYellow); 
        m_currentHistoToFit->SetStats(kTRUE);
        gStyle->SetOptFit(0111);
      //        predictLine->Draw("same");
        //logNormFunct->Draw("same");
        TPaveStats *statsBox = (TPaveStats*)logNormFunct->FindObject("stats");
        if(statsBox){
          statsBox->SetX1NDC(.15);
          statsBox->SetY1NDC(.75);
          statsBox->SetX2NDC(.4);
          statsBox->SetY2NDC(.9);
        }

        TH1D* histoClone = (TH1D*) m_currentHistoToFit->Clone();
         
 
        //FITTING HERE
        int status = m_currentHistoToFit->Fit("logNormFunct","EQ","");//,ZTPCPredict-0.5, ZTPCPredict+0.5);//,-5.0, 5.0);
        if(status == 0){
          TF1* fitResult = m_currentHistoToFit->GetFunction("logNormFunct"); 
          if(fitResult->GetParameter(0) != 0.0){
            if(histoIndex == 0){
              m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][0]->SetBinContent(bin,fitResult->GetParameter(0));
              m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->SetBinContent(bin,fitResult->GetParameter(1));
              m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->SetBinContent(bin,fitResult->GetParameter(2));
              m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][0]->SetBinError(bin,fitResult->GetParError(0));
              m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->SetBinError(bin,fitResult->GetParError(1));
              m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->SetBinError(bin,fitResult->GetParError(2));
            }else{
              m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][0]->SetBinContent(bin,fitResult->GetParameter(0));
              m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->SetBinContent(bin,fitResult->GetParameter(1));
              m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->SetBinContent(bin,fitResult->GetParameter(2));
              m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][0]->SetBinError(bin,fitResult->GetParError(0));
              m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->SetBinError(bin,fitResult->GetParError(1));
              m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->SetBinError(bin,fitResult->GetParError(2));
            }
          }
          //draw result
          fittingCanvas->cd();
          //gPad->SetLogx();
          //gPad->SetLogy();
          gPad->SetRightMargin(0.4);
          m_currentHistoToFit->GetXaxis()->SetMoreLogLabels();
          m_currentHistoToFit->SetLineColor(kRed);
          //histoClone->Draw();
          m_currentHistoToFit->Draw("E");
          initialFunct->Draw("same");
          //logNormFunct->SetNpx(500);
          //logNormFunct->Draw("same");

          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/dEdxMomFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),imageName.Data()));
        }

        delete m_currentHistoToFit;
        delete histoClone;
  
      }//mombin
      delete histoPM;
    }//particle
  }//plus minus


  //Now that we have isolated fits, the btof efficiency doesn't go as low in momentum as we like, but luckily the dEdx is separated enough to fit on its own

  TH2I* histoPlus = HistogramUtilities::reBinVar2D(m_mom_dEdx_Plus,a_rebinMom,1);
  TH2I* histoMinus = HistogramUtilities::reBinVar2D(m_mom_dEdx_Minus,a_rebinMom,1);

  //## FIT THE dEdx distributions at low momentum without BTOF##
  for(int histoIndex = 0; histoIndex < 2; histoIndex++){
    TH2I* histoPM = histoPlus;
    string histoNameStr = "plus";
    if(histoIndex == 1){
      histoPM = histoMinus;
      histoNameStr = "minus";
    }

    for(int bin = 1; bin <= histoPM->GetNbinsX(); bin++){
      TH1D* fullHisto = histoPM->ProjectionY(Form("proj_%02d_%02d",histoIndex,bin),bin,bin);
      double momentum = histoPM->GetXaxis()->GetBinCenter(bin);
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(!m_initialParticleToggles[partIndex]) continue;
        if(partIndex == 3 || partIndex == 8) continue;
        if(fullDeDxStops[partIndex] < momentum) continue;

        double dEdxPredict = m_bichselCurves[partIndex]->Eval(momentum);
        m_currentHistoToFit = HistogramUtilities::truncateHistogram(fullHisto,0.85*dEdxPredict,1.15*dEdxPredict);
        
        double counts = HistogramUtilities::sumOfContent(m_currentHistoToFit, -1);
        if( 200 > counts ){  // && !m_testing)  || (counts < 20)){
          delete m_currentHistoToFit;
          cout << "Continuing: too few counts for solo fit partIndex:" << partIndex << " counts:" << counts << " at mom: " << momentum << endl; 
          continue;
        }
        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
        //m_currentHistoToFit->Sumw2();
        m_currentHistoToFit->Scale(1.0,"width");  
        double integral = m_currentHistoToFit->Integral();
        double mean = m_currentHistoToFit->GetMean();
        double stdDev = m_currentHistoToFit->GetStdDev();
        double mu = 0.0;
        double sig = 0.0;      
        PhysMath::gausToLogNorm(mean,stdDev,mu,sig);       

        

        TString imageName = TString::Format("%s_%s_fullmomBin_%03d",
                  m_partInfo->GetParticleName(partIndex,0).Data(),histoNameStr.c_str(),bin);
        m_currentHistoToFit->SetName(imageName.Data());
        m_currentHistoToFit->SetTitle(Form(
                     "dE/dx Zoomed for %s %s p_{tot}=[%2.3f,%2.3f] (GeV/c)",
                     m_partInfo->GetParticleName(partIndex,0).Data(),
                     histoNameStr.c_str(),
                     histoPM->GetXaxis()->GetBinLowEdge(bin),
                     histoPM->GetXaxis()->GetBinLowEdge(bin+1)
                     ));
        m_currentHistoToFit->GetXaxis()->SetTitle("dE/dx");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/ddEdx");
        m_currentHistoToFit->SetStats(kTRUE);
        //cout << "mean: " << mean << "  stdDev: " << stdDev << "  mu: " << mu << "   sig: " << sig << endl;
        //fitFunction->SetParameter(0,(counts)*m_currentHistoToFitTemp->GetBinWidth(2));
        logNormFunct->SetParameter(0,counts);
        logNormFunct->SetParLimits(0,0.6*counts, 1.5*counts);
        logNormFunct->SetParameter(1,mu);
        logNormFunct->SetParLimits(1, 0.5*mu, 2.0*mu);
        logNormFunct->SetParameter(2,sig);
        logNormFunct->SetParLimits(2,0.05*sig,2.0*sig);
        logNormFunct->SetParNames("N","#mu","#sigma");

        TF1* initialFunct = (TF1*) logNormFunct->Clone();
        initialFunct->SetLineColor(kViolet); 
        m_currentHistoToFit->SetStats(kTRUE);
        
        //predictLine->Draw("same");
        //logNormFunct->Draw("same");
 
        //FITTING HERE
        int status = m_currentHistoToFit->Fit("logNormFunct","EQ","");//,ZTPCPredict-0.5, ZTPCPredict+0.5);//,-5.0, 5.0);
        if(status == 0){
          TF1* fitResult = m_currentHistoToFit->GetFunction("logNormFunct"); 
          if(fitResult->GetParameter(0) != 0.0){
            if(histoIndex == 0){
            
              double oldSigErr = m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->GetBinError(bin);
              if(oldSigErr > 0.0){
                double weight1 = pow(m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->GetBinError(bin),-2.0);
                double weight2 = pow(fitResult->GetParError(2),-2.0);
                double avgSig = (weight1*m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->GetBinContent(bin)
                               +weight2*fitResult->GetParameter(2))/(weight1+weight2);
                double sigsig = (1.0/weight1+1.0/weight2)/pow(weight1+weight2,2.0);

                double Mweight1 = pow(m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->GetBinError(bin),-2.0);
                double Mweight2 = pow(fitResult->GetParError(1),-2.0);
                double avgMu = (Mweight1*m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->GetBinContent(bin)
                               +Mweight2*fitResult->GetParameter(1))/(Mweight1+Mweight2);
                double musig = (1.0/Mweight1+1.0/Mweight2)/pow(Mweight1+Mweight2,2.0);

 
              //if(oldSigErr > 0.0 && fitResult->GetParError(2) < oldSigErr){
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][0]->SetBinContent(bin,fitResult->GetParameter(0));
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->SetBinContent(bin,avgMu);
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->SetBinContent(bin,avgSig);
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][0]->SetBinError(bin,fitResult->GetParError(0));
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->SetBinError(bin,musig);
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->SetBinError(bin,sigsig);
              //}
              }else{
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][0]->SetBinContent(bin,fitResult->GetParameter(0));
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->SetBinContent(bin,fitResult->GetParameter(1));
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->SetBinContent(bin,fitResult->GetParameter(2));
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][0]->SetBinError(bin,fitResult->GetParError(0));
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][1]->SetBinError(bin,fitResult->GetParError(1));
                m_mom_dEdxFitData_BTOFIsolated_Plus[partIndex][2]->SetBinError(bin,fitResult->GetParError(2));
              }
            }else{
              double oldSigErr = m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->GetBinError(bin);
              if(oldSigErr > 0.0){
                double weight1 = pow(m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->GetBinError(bin),-2.0);
                double weight2 = pow(fitResult->GetParError(2),-2.0);
                double avgSig = (weight1*m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->GetBinContent(bin)
                               +weight2*fitResult->GetParameter(2))/(weight1+weight2);
                double sigsig = (1.0/weight1+1.0/weight2)/pow(weight1+weight2,2.0);

                double Mweight1 = pow(m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->GetBinError(bin),-2.0);
                double Mweight2 = pow(fitResult->GetParError(1),-2.0);
                double avgMu = (Mweight1*m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->GetBinContent(bin)
                               +Mweight2*fitResult->GetParameter(1))/(Mweight1+Mweight2);
                double musig = (1.0/Mweight1+1.0/Mweight2)/pow(Mweight1+Mweight2,2.0);


              //if(oldSigErr > 0.0 && fitResult->GetParError(2) < oldSigErr){
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][0]->SetBinContent(bin,fitResult->GetParameter(0));
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->SetBinContent(bin,avgMu);
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->SetBinContent(bin,avgSig);
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][0]->SetBinError(bin,fitResult->GetParError(0));
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->SetBinError(bin,musig);
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->SetBinError(bin,sigsig);
              //}
              }else{
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][0]->SetBinContent(bin,fitResult->GetParameter(0));
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->SetBinContent(bin,fitResult->GetParameter(1));
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->SetBinContent(bin,fitResult->GetParameter(2));
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][0]->SetBinError(bin,fitResult->GetParError(0));
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][1]->SetBinError(bin,fitResult->GetParError(1));
                m_mom_dEdxFitData_BTOFIsolated_Minus[partIndex][2]->SetBinError(bin,fitResult->GetParError(2));
              }
            }
          }
          //draw result
          fittingCanvas->cd();
          //gPad->SetLogx();
          //gPad->SetLogy();
          gPad->SetRightMargin(0.4);
          m_currentHistoToFit->GetXaxis()->SetMoreLogLabels();
          m_currentHistoToFit->SetLineColor(kRed);
          m_currentHistoToFit->Draw("E");
          initialFunct->Draw("same"); 
          //logNormFunct->Draw("same");

          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/dEdxMomFits/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),imageName.Data()));
        }

        delete m_currentHistoToFit;
  
      }//particle
      delete fullHisto;
    }//bin
  }//plus minus

  delete histoPlus;
  delete histoMinus;

  cout << "FINISHED WITH ISOLATED dE/dx vs Mom FITS" << endl;
  return;
}

//============================================================================================================

void ZFitter::electronStudy(double a_lowMomRange, double a_highMomRange, double a_beginDeDxSoloPion, double a_endDeDxSoloPion, double a_dEdxStop){
  cout << "#########   ELECTRON STUDY ##########" << endl;
  
  double lowMomRange = a_lowMomRange;// = 0.29;
  double highMomRange = a_highMomRange;// = 0.33; //0.361
  double beginDeDxSoloPion = a_beginDeDxSoloPion;// = 2.8;
  double endDeDxSoloPion = a_endDeDxSoloPion;// = 3.6;
  //double endDeDxKaon = 6.8; // this should be right past the peak, but not too far because of kaon smearing
  double dEdx_dipLoc = a_dEdxStop;// = 4.75;

  //cout << "making image dir" << endl;
  //gSystem->mkdir(Form("%s/%s/electronStudy",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  cout << "made image dir" << endl;
  TCanvas* fittingCanvas = new TCanvas("fitting_canvas");
  cout << "made fitting canvas" << endl;
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);
  cout << "resized canvas" << endl;
  
  double midMom = 0.5*(lowMomRange+highMomRange);

  if(!m_mom_dEdx_Minus || !m_mom_dEdx_Plus){
    cout << "no m_mom_dEdx_..." << endl;
    return;
  }

  TF1* electronFunct = NULL;
  TF1* positronFunct = NULL;

  TH1D* histo_plus = m_mom_dEdx_Plus->ProjectionY("plus_proj",m_mom_dEdx_Plus->FindBin(lowMomRange),m_mom_dEdx_Plus->FindBin(highMomRange));

  //TH1I* histo_plus_int = HistogramUtilities::projectionY(m_mom_dEdx_Plus,lowMomRange,highMomRange);
  //TH1D* histo_plus  = HistogramUtilities::convertInt_to_DoubleHisto(histo_plus_int);
  //delete histo_plus_int;
  HistogramUtilities::makeSqrtNErrors(histo_plus);
  histo_plus->Scale(1.0,"width");
  TH1D* oldhisto_plus = histo_plus;
  TH1D* trucHisto = HistogramUtilities::truncateHistogram(histo_plus, beginDeDxSoloPion, dEdx_dipLoc + 0.4);
  delete oldhisto_plus;
  histo_plus = trucHisto;
  histo_plus->SetTitle(Form("Positron Study p_{tot} = [%1.2f,%1.2f]; dE/dx (keV/cm); dN/d(dE/dx)",lowMomRange,highMomRange));


  TH1D* histo_minus = m_mom_dEdx_Plus->ProjectionY("minus_proj",m_mom_dEdx_Plus->FindBin(lowMomRange),m_mom_dEdx_Plus->FindBin(highMomRange));
  //TH1I* histo_minus_int = HistogramUtilities::projectionY(m_mom_dEdx_Minus,lowMomRange,highMomRange);
  //TH1D* histo_minus  = HistogramUtilities::convertInt_to_DoubleHisto(histo_minus_int);
  //delete histo_minus_int;
  HistogramUtilities::makeSqrtNErrors(histo_plus);
  histo_minus->Scale(1.0,"width");
  TH1D* oldhisto_minus = histo_minus;
  trucHisto = HistogramUtilities::truncateHistogram(histo_minus, beginDeDxSoloPion, dEdx_dipLoc+0.4);
  delete oldhisto_minus;
  histo_minus = trucHisto;
  histo_minus->SetTitle(Form("Electron Study p_{tot} = [%1.2f,%1.2f]; dE/dx (keV/cm); dN/d(dE/dx)",lowMomRange,highMomRange));

  for(int pm_index = 0; pm_index < 2; pm_index++){
    TH1D* histo_dEdx = NULL;
    if(pm_index == 0){
      histo_dEdx = histo_plus;
    }else{
      histo_dEdx = histo_minus;
    }
    
    //###### fit the pion part of the distribution first  ######
    //cout << "Trying to get guas params" << endl;
    double gausParams[3]={0,0,0};
    HistogramUtilities::getGausData(&gausParams[0], histo_dEdx, true, false,beginDeDxSoloPion,endDeDxSoloPion);
    //cout << " Gaus Params: " << gausParams[0] << "  " << gausParams[1] << "  " << gausParams[2] << endl;
    double pionPredict = m_bichselCurves[0]->Eval(midMom);
    double piMean = pionPredict;
    double piStdDev = gausParams[2];
    double piInt = histo_dEdx->Integral("width"); //counts in range
    double piMu = 0.0;
    double piSig = 0.0;      
    PhysMath::gausToLogNorm(piMean,piStdDev,piMu,piSig); 
    m_logNormFunct->SetRange(beginDeDxSoloPion,endDeDxSoloPion);
    m_logNormFunct->SetParameter(0,piInt);
    m_logNormFunct->SetParLimits(0,0.6*piInt, 20.0*piInt);
    m_logNormFunct->SetParameter(1,piMu);
    m_logNormFunct->SetParLimits(1, 0.5*piMu, 2.0*piMu);
    m_logNormFunct->SetParameter(2,piSig);
    m_logNormFunct->SetParLimits(2,0.05*piSig,2.0*piSig);
    m_logNormFunct->SetParNames("N_{#pi}","#mu_{#pi}","#sigma_{#pi}");

    TF1* initPionFunct = (TF1*) m_logNormFunct->Clone();
    initPionFunct->SetLineColor(kViolet);
  
  
    int status = histo_dEdx->Fit("m_logNormFunct","E","",beginDeDxSoloPion,endDeDxSoloPion);//,ZTPCPredict-0.5, ZTPCPredict+0.5);//,-5.0, 5.0);
    TF1* fitResult = NULL;
    if(status == 0){
      fitResult = histo_dEdx->GetFunction("m_logNormFunct");
      piInt = fitResult->GetParameter(0);
      piMu = fitResult->GetParameter(1);
      piSig = fitResult->GetParameter(2);
    }else{
      cerr << "\n\n\n\n PION FITTING FOR THE ELECTRON STUDY FAILED IN SIMPLE FIT with Status: " << status << endl;
    }

    /*fittingCanvas->cd();
    gStyle->SetOptFit(0111);
    histo_dEdx->Draw("E");
    gPad->SetLogy();
    gPad->SetRightMargin(0.4);
    histo_dEdx->SetStats(true);
    histo_dEdx->GetXaxis()->SetMoreLogLabels();
    gSystem->ProcessEvents();
    fittingCanvas->SaveAs(Form("%s/%s/electronStudy/electronStudy_pionFit.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
    */

    TF1* pionFunct = (TF1*) m_logNormFunct->Clone();
    pionFunct->SetLineColor(m_partInfo->GetParticleColor(0));
    pionFunct->SetRange(lowMomRange,5.0);

    //get initial kaon parameters
    /*double gausParamsKaon[3]={0,0,0};
    HistogramUtilities::getGausData(&gausParamsKaon[0], histo_dEdx, true, false,dEdx_dipLoc,endDeDxKaon+0.3);
    cout << " Gaus Params: " << gausParamsKaon[0] << "  " << gausParamsKaon[1] << "  " << gausParamsKaon[2] << endl;
    double kaonPredict = m_bichselCurves[1]->Eval(midMom);
    double kaonMean = kaonPredict;
    double kaonStdDev = gausParamsKaon[2];
    double kaonInt = gausParamsKaon[0]; //counts in range
    double kaonMu = 0.0;
    double kaonSig = 0.0;      
    PhysMath::gausToLogNorm(kaonMean,kaonStdDev,kaonMu,kaonSig); 
    */

    //histo_dEdx->GetXaxis()->SetRangeUser(pionPredict,kaonPredict);
    string firstPart ="[0]*0.3989422804 * TMath::Exp(-0.5*( TMath::Log(x) - [1] )*( TMath::Log(x) - [1] )/([2]*[2])) /([2] * x)";
    string secondPart = "+[3]*0.3989422804 * TMath::Exp(-0.5*( TMath::Log(x) - [4] )*( TMath::Log(x) - [4] )/([5]*[5])) /([5] * x)";
    //string thirdPart = "+[6]*0.3989422804 * TMath::Exp(-0.5*( TMath::Log(x) - [7] )*( TMath::Log(x) - [7] )/([8]*[8])) /([8] * x)";
    string totalParts = firstPart+secondPart;//+thirdPart;
    TF1* twoLogNormFunct = new TF1("twoLogNormFunct",totalParts.c_str(),beginDeDxSoloPion,dEdx_dipLoc+0.5);
    twoLogNormFunct->FixParameter(0,piInt);
    twoLogNormFunct->FixParameter(1,piMu);
    twoLogNormFunct->FixParameter(2,piSig);
    //twoLogNormFunct->FixParameter(6,0);
    //twoLogNormFunct->FixParameter(7,kaonMu);
    //twoLogNormFunct->FixParameter(8,kaonSig);
  
    /*twoLogNormFunct->SetParameter(6,kaonInt);
    twoLogNormFunct->SetParameter(7,kaonMu);
    twoLogNormFunct->SetParameter(8,kaonSig);
    twoLogNormFunct->SetParLimits(6,0.5* kaonInt,2.0*kaonInt);
    twoLogNormFunct->SetParLimits(7,0.5* kaonMu,2.0*kaonMu);
    twoLogNormFunct->SetParLimits(8,0.5* kaonSig,2.0*kaonSig);
    */
    double elecMean = m_bichselCurves[3]->Eval(midMom);
    double elecStdDev = 0.07*elecMean;
    double elecMu = 0;
    double elecSig = 0;
    PhysMath::gausToLogNorm(elecMean,elecStdDev,elecMu,elecSig);
  
    double elecInt = histo_dEdx->GetBinContent(histo_dEdx->FindBin(elecMean))/logNorm(elecMu,elecSig,elecMean);
    twoLogNormFunct->SetParameter(3,elecInt);
    twoLogNormFunct->SetParLimits(3,0.5* elecInt,20.0*elecInt);
    twoLogNormFunct->SetParameter(4,elecMu);
    twoLogNormFunct->SetParLimits(4,0.5* elecMu,2.0*elecMu);
    twoLogNormFunct->SetParameter(5,elecSig);
    twoLogNormFunct->SetParLimits(5,0.5* elecSig,2.0*elecSig);
    twoLogNormFunct->SetParNames("N_{#pi}","#mu_{#pi}","#sigma_{#pi}",
            "N_{e}","#mu_{e}","#sigma_{e}");//,"N_{K}","#mu_{K}","#sigma_{K}");
    TF1* initTotFunct = (TF1*) twoLogNormFunct->Clone();
    initTotFunct->SetLineColor(kViolet);
  
    status = histo_dEdx->Fit("twoLogNormFunct","E","",endDeDxSoloPion,dEdx_dipLoc);  //FIT THE THING
    fitResult = NULL;
    if(status == 0){
      fitResult = histo_dEdx->GetFunction("twoLogNormFunct");
      elecInt = fitResult->GetParameter(3);
      elecMu = fitResult->GetParameter(4);
      elecSig = fitResult->GetParameter(5);
      
    }else{
      cerr << "\n\n\nELECTRON FITTING FOR THE ELECTRON STUDY FAILED IN SIMPLE FIT with Status: " << status << endl;
    }
  
    TF1* elecFunct = (TF1*) m_logNormFunct->Clone();
    elecFunct->SetParameter(0,elecInt);
    elecFunct->SetParameter(1,elecMu);
    elecFunct->SetParameter(2,elecSig);
    elecFunct->SetLineColor(m_partInfo->GetParticleColor(3));
    elecFunct->SetRange(lowMomRange,10);
  
  
    fittingCanvas->cd();
    gStyle->SetOptFit(0111);
    gPad->SetLogy();
    gPad->SetRightMargin(0.4);
    histo_dEdx->SetStats(true);
    histo_dEdx->Draw("E");
    histo_dEdx->GetXaxis()->SetMoreLogLabels();
    //initTotFunct->Draw("same");
    pionFunct->Draw("same");
    elecFunct->Draw("same");
    if(fitResult){
      fitResult->SetRange(beginDeDxSoloPion,dEdx_dipLoc);
      fitResult->Draw("same");
    }else{
      cout << "ERROR: ELECTRON STUDY FAILED!!!!!" << endl;
    }

    if(pm_index == 0){
      fittingCanvas->SaveAs(Form("%s/%s/electronStudy/electronStudy_totalFit.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
      positronFunct = elecFunct;
    }
    if(pm_index == 1){
      fittingCanvas->SaveAs(Form("%s/%s/electronStudy/positronStudy_totalFit.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
      electronFunct = elecFunct;
      if(fitResult) m_basicElectronPionRatio = fitResult->GetParameter(3) / fitResult->GetParameter(0);
    }
  
    //delete histo_dEdx;
    delete pionFunct;
    delete twoLogNormFunct;
    delete initTotFunct;
    delete initPionFunct;

  }//end pm loop


  fittingCanvas->cd();
  //gStyle->SetOptFit(0111);
  gPad->SetLogy();
  gPad->SetRightMargin(0.4);
  histo_plus->SetStats(false);
  histo_plus->Draw("E");
  histo_plus->GetXaxis()->SetMoreLogLabels();
  histo_minus->SetStats(false);
  histo_minus->Draw("E SAME");
  //initTotFunct->Draw("same");
  electronFunct->Draw("same");
  positronFunct->Draw("same");

  fittingCanvas->SaveAs(Form("%s/%s/electronStudy/electron_positron.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
    
  
  m_dEdxParams_Electrons[0][0] = positronFunct->GetParameter(1);
  m_dEdxParams_Electrons[0][1] = positronFunct->GetParameter(2);
  m_dEdxParams_Electrons[1][0] = electronFunct->GetParameter(1); // mu
  m_dEdxParams_Electrons[1][1] = electronFunct->GetParameter(2); // sig

  PhysMath::logNormToGaus( m_dEdxParams_Electrons[0][0],  m_dEdxParams_Electrons[0][1], 
                           m_dEdxParams_Electrons[0][2],  m_dEdxParams_Electrons[0][3]); // mu sig mean stdDev
  PhysMath::logNormToGaus( m_dEdxParams_Electrons[1][0],  m_dEdxParams_Electrons[1][1], 
                           m_dEdxParams_Electrons[1][2],  m_dEdxParams_Electrons[1][3]);


  


  //Pushing the fit data to the rapidity and mTm0 bins
  TH3D* current3DHisto = NULL;
  for(int partIndex = 0; partIndex < 9; partIndex++){
    //cout << "Trying to Push electron fit to particle " << partIndex << endl;
    if(!m_particleLoaded[partIndex]) continue;
    //cout << "Pushing electron fit to particle " << partIndex << endl;
    current3DHisto = m_rap_mTm0_dEdx_Plus_allCent[partIndex];
    if(!current3DHisto) cout << " there is no dEdx all cent plus for partIndex =" << partIndex << endl;
    for(int pmIndex = 0; pmIndex < 2; pmIndex++){
      //cout << "Getting Bins for pmIndex" << pmIndex << endl;

      if(pmIndex == 1){
        if(!m_partToFitMatrix_Minus[partIndex][3]) continue;
      }else{
        if(!m_partToFitMatrix_Plus[partIndex][3]) continue;
      }
        



      int numHorBins  = current3DHisto->GetNbinsX();
      int numVertBins = current3DHisto->GetNbinsY();
      for(int HorBinIndex = 1; HorBinIndex <= numHorBins; HorBinIndex++){
        for(int VertBinIndex = 1; VertBinIndex <= numVertBins; VertBinIndex++){
          m_currentRapBin = HorBinIndex;
          m_currentMtM0Bin = VertBinIndex;
          double rap = current3DHisto->GetXaxis()->GetBinCenter(m_currentRapBin);
          double mTm0 = current3DHisto->GetYaxis()->GetBinCenter(m_currentMtM0Bin);
          #ifdef _ZFITTER_DEBUG_
            cout << " now doing rap: " << rap << " mtm0: " << mTm0 << endl;
            cout << "  Necessary Pointers: " << m_Fit_Data_ZTPC[partIndex][3][1][0] << "  " << m_Fit_Data_ZTPC[partIndex][3][2][0] << "  " << m_Fit_Data_ZTPC[partIndex][3][1][1] << "  " << m_Fit_Data_ZTPC[partIndex][3][2][1] << endl;

          #endif
          m_currentMomOfBin = PhysMath::pTot_from_y_mTm0(rap, mTm0, m_partInfo->GetParticleMass(partIndex));
          if(pmIndex ==0){

            double newElectronMean = m_dEdxParams_Electrons[0][0] - log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin));
            double newSig = m_dEdxParams_Electrons[0][1];
            //cout << m_dEdxFitDataPlus[partIndex][3][1]->GetTitle();

            if(m_bichselShiftsComplete && m_bichselShifts_Plus[partIndex][0]){
              cout << "     electron mean before: " << newElectronMean << endl;
              newElectronMean += m_bichselShifts_Plus[partIndex][0]->GetBinContent(HorBinIndex); // correct electron by pion shift
              cout << "     electron mean after BichselShifts: " << newElectronMean << endl;
            }
            //if(HorBinIndex == 14 && VertBinIndex < 10) cout << "setting part index " << partIndex << " bin " << m_currentRapBin << " " << m_currentMtM0Bin << " mean: " << newElectronMean << " and stdDev: " << newSig << endl;
            m_Fit_Data_ZTPC[partIndex][3][1][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,newElectronMean);
            m_Fit_Data_ZTPC[partIndex][3][2][0]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,newSig);
          }else{
            double newElectronMean = m_dEdxParams_Electrons[1][0] - log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin));
            double newSig = m_dEdxParams_Electrons[1][1];
            if(m_bichselShiftsComplete  && m_bichselShifts_Minus[partIndex][0]){
              newElectronMean += m_bichselShifts_Minus[partIndex][0]->GetBinContent(HorBinIndex); // correct electron by pion shift
            }
            //if(HorBinIndex == 14 && VertBinIndex < 10) cout << "setting part index " << partIndex << " bin " << m_currentRapBin << " " << m_currentMtM0Bin << " mean: " << newElectronMean << " and stdDev: " << newSig << endl;

            m_Fit_Data_ZTPC[partIndex][3][1][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,newElectronMean);
            m_Fit_Data_ZTPC[partIndex][3][2][1]->SetBinContent(m_currentRapBin,m_currentMtM0Bin,newSig);
          }
        }//vert
      }//hor
    }//pm
  }



  cout << "Final Positron Parameters:  Mu: " << m_dEdxParams_Electrons[0][0]<< " Sig: " <<m_dEdxParams_Electrons[0][1] << " Mean: " << m_dEdxParams_Electrons[0][2]<< " StdDev: " << m_dEdxParams_Electrons[0][3]<< endl;
  cout << "Final Electron Parameters:  Mu: " << m_dEdxParams_Electrons[1][0]<< " Sig: " <<m_dEdxParams_Electrons[1][1] << " Mean: " << m_dEdxParams_Electrons[1][2]<< " StdDev: " << m_dEdxParams_Electrons[1][3]<< endl;
  cout << " Checking Electron Params : " << m_Fit_Data_ZTPC[0][3][2][1]->GetBinContent(14,2) << endl;

  
  m_electronStudyComplete = true;
  cout << "FINISHED ELECTRON STUDY" << endl;

}

//============================================================================================================


void ZFitter::fitElectronDeDxByRapidity(int a_partIndex){
  #ifdef _ZFITTER_DEBUG_
    cout << "Call to ZFitter::fitElectronDeDxByRapidity()" << endl;
  #endif

  //cout << "making image dir" << endl;
  //gSystem->mkdir(Form("%s/%s/electronStudy",m_imagePreDir.c_str(), m_imgDirName.c_str()));
  TCanvas* fittingCanvas = new TCanvas("fitting_canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);
  #ifdef _ZFITTER_DEBUG_
    cout << "Ptrs: " << m_electron_rap_dEdx[a_partIndex][0]  << "  " << m_electron_rap_dEdx[a_partIndex][1] << endl;
  #endif
  

  if(!m_electron_rap_dEdx[a_partIndex][0] || !m_electron_rap_dEdx[a_partIndex][1]){
    cout << "no m_electron_rap_dEdx for part " << a_partIndex << ", so cant do electron fitting..." << endl;
    return;
  }

  // make the rapidity histograms for the fit data
  string paramNames[3] = {"N","mu","sig"};
  string paramSymbols[3] = {"N","#mu","#sigma"};
  TH1D* rapBase = m_electron_rap_dEdx[a_partIndex][0]->ProjectionX(Form("rapBase_%d",a_partIndex),1,1);
  HistogramUtilities::zeroBins(rapBase);
  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    for(int iii = 0; iii < 3; iii++){
      m_electronDeDxParamsByRapidity[a_partIndex][iii][pmIndex] = (TH1D*) rapBase->Clone();
      m_electronDeDxParamsByRapidity[a_partIndex][iii][pmIndex]->SetName(
            Form("%s_byRap_%s",(pmIndex == 0 ? "positron" : "electron"),paramNames[iii].c_str()));
      m_electronDeDxParamsByRapidity[a_partIndex][iii][pmIndex]->SetTitle(
            Form("%s Fit %s; y_{%s}; %s",(pmIndex == 0 ? "Positron" : "Electron"),
              paramNames[iii].c_str(),m_partInfo->GetParticleSymbol(a_partIndex).Data(),paramSymbols[iii].c_str()));
    }
  }
  


  for(int rapBin = 1; rapBin <= rapBase->GetNbinsX(); rapBin++){
    for(int pmIndex = 0; pmIndex < 2; pmIndex++){
      TH2I* histo2D = m_electron_rap_dEdx[a_partIndex][pmIndex];
      TH1D* histo = histo2D->ProjectionY("proj",rapBin,rapBin);
      double content = HistogramUtilities::sumOfContent(histo);
      if(content < 30) continue;
      HistogramUtilities::makeSqrtNErrors(histo);
      TH1D* rebinHisto = HistogramUtilities::rebinEqualByCDF(histo,30,false,false);
      m_currentHistoToFit = HistogramUtilities::logifyHistoHorizontalAxis(rebinHisto);
      m_currentHistoToFit->Scale(1.0,"width");
      delete histo;
      delete rebinHisto;
      m_currentHistoToFit->SetTitle(Form("%s y = [%1.2f,%1.2f]",(pmIndex==0 ? "Positron" : "Electron"),
        histo2D->GetXaxis()->GetBinLowEdge(rapBin),histo2D->GetXaxis()->GetBinLowEdge(rapBin+1)));

      double gausParams[3]={0,0,0};
      HistogramUtilities::getGausData(&gausParams[0], m_currentHistoToFit, true, false);
      //cout << " Gaus Params: " << gausParams[0] << "  " << gausParams[1] << "  " << gausParams[2] << endl;
      double elecMean = gausParams[1];
      double elecStdDev = gausParams[2];
      double elecInt = m_currentHistoToFit->Integral("width"); //counts in range
    
      TF1* fitFunct = PhysMath::getGaussianFunction("electronFitFunct", elecInt, elecMean, elecStdDev);
      fitFunct->SetParLimits(0,0.6*elecInt, 20.0*elecInt);
      fitFunct->SetParLimits(1, 0.5*elecMean, 2.0*elecMean);
      fitFunct->SetParLimits(2,0.05*elecStdDev,2.0*elecStdDev);
      fitFunct->SetParNames("N_{e}","#mu_{e}","#sigma_{e}");
      fitFunct->SetLineColor(kViolet);
  
  
      int status = m_currentHistoToFit->Fit("electronFitFunct","EI");
      TF1* fitResult = NULL;
      if(status == 0){
        fitResult = m_currentHistoToFit->GetFunction("electronFitFunct");
        #ifdef _ZFITTER_DEBUG_
          elecInt  = fitResult->GetParameter(0);
          elecMean = fitResult->GetParameter(1);
          elecStdDev  = fitResult->GetParameter(2);
          cout << " RapBin: " << rapBin << " N: " << elecInt << " mu: " << elecMean << " sig: " << elecStdDev << endl;
        #endif
        for(int iii = 0; iii < 3; iii++){
          m_electronDeDxParamsByRapidity[a_partIndex][iii][pmIndex]->SetBinContent(rapBin,fitResult->GetParameter(iii));
          m_electronDeDxParamsByRapidity[a_partIndex][iii][pmIndex]->SetBinError(rapBin,fitResult->GetParError(iii));
        }
      }else{
        cerr << "\n\n\n\n PION FITTING FOR THE ELECTRON STUDY FAILED IN SIMPLE FIT with Status: " << status << endl;
      }
      fittingCanvas->cd();
      gStyle->SetOptFit(0111);
      gPad->SetRightMargin(0.4);
      m_currentHistoToFit->SetStats(true);
      m_currentHistoToFit->Draw("E");
      if(fitResult) fitResult->Draw("same");
      if(pmIndex == 0){
        fittingCanvas->SaveAs(Form("%s/%s/electronStudy/electron_rapBin%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),rapBin));
      }
      if(pmIndex == 1){
        fittingCanvas->SaveAs(Form("%s/%s/electronStudy/positron_rapBin%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),rapBin));
      }
      delete m_currentHistoToFit;
    }
  }

  // fit the parameter trends as a function of rapidity
  for(int paramIndex = 0; paramIndex <= 3; paramIndex++){
    for(int pmIndex = 0; pmIndex < 2; pmIndex++){
      m_currentHistoToFit = m_electronDeDxParamsByRapidity[a_partIndex][paramIndex][pmIndex];
      int middleRapBin = m_currentHistoToFit->GetNbinsX()/2;
      double middleRap = m_currentHistoToFit->GetBinCenter(middleRapBin);
      double middleRapValue = m_currentHistoToFit->GetBinContent(middleRapBin);
      string functName = Form("electronFitFunct_%d_%d",paramIndex,pmIndex);
      TF1* fitFunct =NULL;
      if(paramIndex == 0){
        fitFunct = PhysMath::getGaussianFunction(functName.c_str(),1.0,0.0,0.5);
        double amp = middleRapValue/fitFunct->Eval(middleRapValue);
        fitFunct->SetParameter(0,amp);
        fitFunct->SetParLimits(0,0.01*amp,100.0*amp);
        fitFunct->SetParLimits(1,-1.0,1.0);
        fitFunct->SetParLimits(2,0.001,2.0);
      }else if(paramIndex == 1){
        fitFunct = new TF1(functName.c_str(),
          "[3] - ([0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2]))/[2])",-5,5);
        fitFunct->SetParNames("c","N","#mu","#sigma");
        fitFunct->SetParameter(0,middleRapValue);
        fitFunct->SetParameter(1,0.0);
        fitFunct->FixParameter(2,0.0);
        fitFunct->SetParameter(3,0.5);
      }else if(paramIndex == 2){
        fitFunct = new TF1(functName.c_str(),"[0]");
        fitFunct->SetParameter(0,middleRapValue);
        fitFunct->SetParNames("p_0");
      }
      fitFunct->SetLineColor(kViolet);
  
      #ifdef _ZFITTER_DEBUG_
        HistogramUtilities::printHistoInfo(m_currentHistoToFit);
      #endif
      int status = m_currentHistoToFit->Fit(functName.c_str(),"E");
      if(status == 0){
        m_electronDeDxParamFitsByRapidity[a_partIndex][paramIndex][pmIndex] = m_currentHistoToFit->GetFunction(functName.c_str());
      }else{
        cerr << "\n\n\n\n ELECTRON FITTING FOR THE ELECTRON STUDY FAILED IN PARAMETER FIT with Status: " << status << endl;
      }


      fittingCanvas->cd();
      gStyle->SetOptFit(0111);
      gPad->SetRightMargin(0.4);
      m_currentHistoToFit->SetStats(true);
      m_currentHistoToFit->Draw("E");
      if(m_electronDeDxParamFitsByRapidity[a_partIndex][paramIndex][pmIndex] && status == 0) m_electronDeDxParamFitsByRapidity[a_partIndex][paramIndex][pmIndex]->Draw("same");
      if(pmIndex == 0){
        fittingCanvas->SaveAs(Form("%s/%s/electronStudy/electron_%s_rapFit.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),paramNames[paramIndex].c_str()));
      }
      if(pmIndex == 1){
        fittingCanvas->SaveAs(Form("%s/%s/electronStudy/positron_%s_rapFit.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),paramNames[paramIndex].c_str()));
      }
    }
  }
  #ifdef _ZFITTER_DEBUG_
    cout << "FINISHED ELECTRON RAPIDITY STUDY" << endl;
  #endif
}



//============================================================================================================


void ZFitter::fitDeDxVsMom(){
  //make the data structures
  m_ZTPCshift = 0;
  bool useGMM = true;
  int numParticleBackup = m_numParticles;
  m_numParticles = 4;

  TCanvas* fittingCanvas = new TCanvas("dEdx_Mom_Fitting");
  TH1D* binBase = m_mom_dEdx_Plus->ProjectionX("momBins");
  HistogramUtilities::zeroBins(binBase);

  
  
  
  //Make the Fit Result Histograms for Parameter Storage
  string paramNames[3] = {"amp","mu","sig"};
  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int paramIndex = 0; paramIndex < 3; paramIndex++){
      m_mom_dEdxFitData[partIndex][paramIndex]= (TH1D*) binBase->Clone();
      m_mom_dEdxFitData[partIndex][paramIndex]->SetName(Form("%s_%s",m_partInfo->GetParticleName(partIndex,0).Data(), paramNames[paramIndex].c_str()));
      m_mom_dEdxFitData[partIndex][paramIndex]->SetTitle(Form("%s %s; P (GeV/c), %s",m_partInfo->GetParticleName(partIndex,0).Data(),paramNames[paramIndex].c_str(), paramNames[paramIndex].c_str()));
    }
  }

  for(int momBin = 1; momBin <= m_mom_dEdx_Plus->GetNbinsX(); momBin++){
    m_currentMomOfBin = m_mom_dEdx_Plus->GetXaxis()->GetBinCenter(momBin);
    if( m_currentMomOfBin < 0.21 || 0.44 < m_currentMomOfBin) continue; // just fit the electrons
    m_currentHistoToFit = m_mom_dEdx_Plus->ProjectionY(Form("dEdxProj_%04d",momBin),momBin,momBin);
    double lowRange = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::firstBinWithContent(m_currentHistoToFit));
    double highRange = m_currentHistoToFit->GetBinLowEdge(HistogramUtilities::lastBinWithContent(m_currentHistoToFit)+1);
    double countsInHisto = HistogramUtilities::sumOfContent(m_currentHistoToFit);
    if((countsInHisto < 200 && !m_testing) || (countsInHisto < 20)) continue;
    cout << "Fitting Mom Bin: " << momBin << endl;
    TH1D* tempHistoPtr = m_currentHistoToFit;
    m_currentHistoToFit = HistogramUtilities::rebinMinCount(m_currentHistoToFit,8.0);
    delete tempHistoPtr;
    //make histo differential
    //m_currentHistoToFit->Sumw2();
    //Sumw2 created with rebinMinCount
    HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
    m_currentHistoToFit->Scale(1.0,"width");


    // ########### SETUP PARTICLE TAGS ##################
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      //see if there is data at prediction
      if(!m_initialParticleToggles[partIndex]){
        m_currentParticleTags[partIndex] = 0;
        continue;
      }
    

      if(partIndex != 3 && partIndex != 8 ){
        double dEdxPredict = m_bichselCurves[partIndex]->Eval(m_currentMomOfBin) +  m_ZTPCshift;
        int    dEdxPredictBin = m_currentHistoToFit->FindBin(dEdxPredict);
        double countsAtPredict = m_currentHistoToFit->GetBinContent(dEdxPredictBin)
                                 * m_currentHistoToFit->GetBinWidth(dEdxPredictBin);
        if(countsAtPredict > 15 && dEdxPredict > lowRange && dEdxPredict < highRange){
          m_currentParticleTags[partIndex] = 2;
        }else{
          m_currentParticleTags[partIndex] = 0;
        }

      }else if(partIndex == 3){ // if electron
        if(m_currentMomOfBin < 0.45){
          m_currentParticleTags[partIndex] = 2;
        }    
        //m_currentParticleTags[partIndex] = 0;
      }else if(partIndex == 8){ // if electron
        m_currentParticleTags[partIndex] = 0;
      }
      
    }//end particle loop

    #ifdef _ZFITTER_DEBUG_
    for(int iii = 0; iii < 8; iii++){
      if(!m_partToFitMatrix_Plus[m_currentPartIndex][iii]) continue;
      cout << "   " << m_partInfo->GetParticleName(iii) << "  has fit tag: " << m_currentParticleTags[iii] << " : " << m_particleTagDescriptions[m_currentParticleTags[iii]] << endl;
    }
    #endif

    //initialize parameters
    initdEdxMomParams();
    cout << "Num Particles To Fit" << endl;
    if(m_currentNumPartToFit == 0) continue;
    cout << "Doing the Fit" << endl;

    double lowRedRange = 0;
    double highRedRange = 0;
    getRangeDeDx(lowRedRange, highRedRange);
    TH1D* origHistToFit = m_currentHistoToFit;
    m_currentHistoToFit = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowRedRange,highRedRange);
    delete origHistToFit;

    


    if(useGMM){
      
      TH1D* tempHistoThing = m_currentHistoToFit;
      m_currentHistoToFit = HistogramUtilities::convertFromLogBins(m_currentHistoToFit);
      delete tempHistoThing;

      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);

      GausMixModel model;
      vector<double> initialMeans(m_currentNumPartToFit,0.0);
      vector<double> initialAmps(m_currentNumPartToFit,0.0);
      model.setNumberOfGaussians(m_currentNumPartToFit);
      model.importDataAsHistogram(m_currentHistoToFit);
      vector<double> finalParams(3*m_currentNumPartToFit);
      int liveParticleCounter = 0;
      double ampTotalTemp = 0.0;
      for(int partIndex = 0; partIndex < 9; partIndex++){
        if(m_currentParticleTags[partIndex] != 0){
          initialMeans[liveParticleCounter] = log(m_bichselCurves[partIndex]->Eval(m_currentMomOfBin));
          cout << " " << liveParticleCounter << "  mean " << initialMeans[liveParticleCounter] << endl;
          initialAmps[liveParticleCounter] = m_currentHistoToFit->GetBinContent(m_currentHistoToFit->FindBin(initialMeans[liveParticleCounter]));
          ampTotalTemp += initialAmps[liveParticleCounter];
          liveParticleCounter++;
        }
      }

      for(int iii = 0; iii < m_currentNumPartToFit; iii++){
        initialAmps[iii] = initialAmps[iii] / ampTotalTemp;
        cout << " " << iii << "  amp " << initialAmps[iii] << endl;
      }
      
      //model.doGMM(&finalParams[0], &initialMeans[0],&initialAmps[0],300,1e-7, 0.02);
      model.doGMM(&finalParams[0], &initialMeans[0],NULL,300,1e-7, 0.02);
      m_currentParameters = finalParams;
      TF1* fitResult = model.getFunction();


      vector< TF1* > indFits;
      for(int iii = 0; iii < m_currentNumPartToFit; iii++){
        indFits.push_back(model.getSingleGaus(iii));
      }


  
      int paramIndex = 0; //the particles not included in fit
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] !=0){
          m_mom_dEdxFitData[partIndex][0]->SetBinContent(momBin,m_currentParameters[paramIndex]);
          m_mom_dEdxFitData[partIndex][1]->SetBinContent(momBin,m_currentParameters[paramIndex+1]);
          m_mom_dEdxFitData[partIndex][2]->SetBinContent(momBin,m_currentParameters[paramIndex+2]);
          paramIndex+=3;
        }
      }
  
      
      //save fit
      fittingCanvas->cd();
      gPad->SetLogy();
      //gPad->SetLogx();
      m_currentHistoToFit->SetStats(kTRUE);
      gStyle->SetOptFit(0111);
      m_currentHistoToFit->Draw("E");
      //m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange,highRange);
      fitResult->SetLineColor(kRed);
      fitResult->Draw("same");
      for(int iii = 0; iii < m_currentNumPartToFit; iii++){
        indFits[iii]->Draw("same");
      }
      //drawIndividualFits();
      gPad->Update();
      gSystem->ProcessEvents();
      fittingCanvas->SaveAs(Form("%s/%s/dEdxMomFits/MomBin_%02d_dEdxFit.png",m_imagePreDir.c_str(), m_imgDirName.c_str(), momBin));

      for(int iii = 0; iii < m_currentNumPartToFit; iii++){
        delete indFits[iii];
      }
      
      



    }else{
      //do the fit
      TF1* fitFunct = new TF1("fitFunct",this,&ZFitter::dEdxOperator,lowRange,highRange,m_currentFreeParameters,"ZFitter","dEdxOperator");
      fitFunct->SetNpx(600);
      initDEDXFunction(fitFunct);
      TF1* fitFunctInitial = (TF1*) fitFunct->Clone();
      fitFunctInitial->SetLineColor(kYellow);
      fitFunctInitial->SetName("fitFunctInitial");
      //fitFunctInitial->Draw("same");
  
      int status = m_currentHistoToFit->Fit("fitFunct","E","",lowRange,highRange);//add Q to quiet
      TF1* fitResult = m_currentHistoToFit->GetFunction("fitFunct");
  
      if(status!=0){
        cout << "Status: " << status << endl;
        cout << "\n\n\n\n\n\n BAD FIT \n\n\n\n\n\n";
      }
  
  
      vector< double > errors;
      int paramIndexCounter = 0;
      for(int iii = 0; iii < (int) m_currentParameters.size(); iii++){
        if(m_paramToggles[iii]){
          m_currentParameters[iii] = fitResult->GetParameter(paramIndexCounter);
          errors.push_back(fitResult->GetParError(paramIndexCounter));
          paramIndexCounter++;
        }
      }
  
      double cNDF = fitResult->GetChisquare() / fitResult->GetNDF();
      cout << "Chi-Squared / ndf = " << cNDF << endl;  
  
      int paramIndex = 0; //the particles not included in fit
      for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
        if(m_currentParticleTags[partIndex] !=0){
          m_mom_dEdxFitData[partIndex][0]->SetBinContent(momBin,m_currentParameters[paramIndex]);
          m_mom_dEdxFitData[partIndex][1]->SetBinContent(momBin,m_currentParameters[paramIndex+1]);
          m_mom_dEdxFitData[partIndex][2]->SetBinContent(momBin,m_currentParameters[paramIndex+2]);
          
          m_mom_dEdxFitData[partIndex][0]->SetBinError(momBin,errors[paramIndex]);    
          m_mom_dEdxFitData[partIndex][1]->SetBinError(momBin,errors[paramIndex+1]);
          m_mom_dEdxFitData[partIndex][2]->SetBinError(momBin,errors[paramIndex+2]);
          paramIndex+=3;
        }
      }
  
  
      //save fit
      fittingCanvas->cd();
      gPad->SetLogy();
      gPad->SetLogx();
      m_currentHistoToFit->SetStats(kTRUE);
      gStyle->SetOptFit(0111);
      m_currentHistoToFit->Draw("E");
      m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange,highRange);
      fitResult->SetLineColor(kRed);
      fitResult->Draw("same");
      fitFunctInitial->Draw("same");
      //drawIndividualFits();
      gPad->Update();
      gSystem->ProcessEvents();
      fittingCanvas->SaveAs(Form("%s/%s/dEdxMomFits/MomBin_%02d_dEdxFit.png",m_imagePreDir.c_str(), m_imgDirName.c_str(), momBin));
      }




  }//end of momentum loop
  

  //set back to usual number of particles
  m_numParticles = numParticleBackup;
  delete fittingCanvas;
}










void ZFitter::fitBetaGammaDeDx(bool a_isBTOF, double a_betaGammaFitStart, bool a_saveImages){
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);

  // loop over plus and then minus particles
  for(int histoIndex = 0; histoIndex < 2; histoIndex++){
    TH2I* histoPM = NULL;
    string histoNameStr = "";
    if(a_isBTOF){
      if(histoIndex==0){
        histoPM = m_betaGammaBTOF_dEdx_Plus;
        histoNameStr = "btof_plus";
      }else{
        histoPM = m_betaGammaBTOF_dEdx_Minus;
        histoNameStr = "btof_minus";
      }
    }else{
      if(histoIndex==0){
        histoPM = m_betaGammaETOF_dEdx_Plus;
        histoNameStr = "etof_plus";
      }else{
        histoPM = m_betaGammaETOF_dEdx_Minus;
        histoNameStr = "etof_minus";
      }
    }
    
    //loop over betaGamma Bins
    for(int betaGammaBin = 1; betaGammaBin <= histoPM->GetNbinsY(); betaGammaBin++){
      cout << "Fitting dEdx for betaGammaBin: " << betaGammaBin << endl;
      m_currentHistoToFit = histoPM->ProjectionY(
           Form("betagamma_dEdx_%02d",betaGammaBin),betaGammaBin, betaGammaBin);
      
      double counts = HistogramUtilities::sumOfContent(m_currentHistoToFit, -1);
      if( 500 > counts ){  // && !m_testing)  || (counts < 20)){
        delete m_currentHistoToFit;
        //cout << "Continuing: too few counts for solo fit partIndex:" << partIndex << " counts:" << counts << " at mom: " << momentum << endl; 
        continue;
      }
      HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
      //m_currentHistoToFit->Sumw2();
      m_currentHistoToFit->Scale(1.0,"width");  
      double integral = m_currentHistoToFit->Integral();
      double mean = m_currentHistoToFit->GetMean();
      double stdDev = m_currentHistoToFit->GetStdDev();
      double mu = 0.0;
      double sig = 0.0;      
      PhysMath::gausToLogNorm(mean,stdDev,mu,sig);
      TString imageName = TString::Format("dEdx_%s_betagammma_%02d",histoNameStr.c_str(),betaGammaBin);
      m_currentHistoToFit->SetName(imageName.Data());
      m_currentHistoToFit->SetTitle(Form(
                   "dE/dx for %s #beta#gamma=[%2.3f,%2.3f]",
                   histoNameStr.c_str(),
                   histoPM->GetYaxis()->GetBinLowEdge(betaGammaBin),
                   histoPM->GetYaxis()->GetBinLowEdge(betaGammaBin+1)
                   ));
      m_currentHistoToFit->GetXaxis()->SetTitle("dE/dx");
      m_currentHistoToFit->GetYaxis()->SetTitle("dN/ddEdx");
      m_currentHistoToFit->SetStats(kTRUE);
      m_logNormFunct->SetParameter(0,counts);
      m_logNormFunct->SetParLimits(0,0.6*counts, 1.5*counts);
      m_logNormFunct->SetParameter(1,mu);
      m_logNormFunct->SetParLimits(1, 0.5*mu, 2.0*mu);
      m_logNormFunct->SetParameter(2,sig);
      m_logNormFunct->SetParLimits(2,0.05*sig,2.0*sig);
      m_logNormFunct->SetParNames("N","#mu","#sigma");
      TF1* initialFunct = (TF1*) m_logNormFunct->Clone();
      initialFunct->SetLineColor(kViolet); 
      m_currentHistoToFit->SetStats(kTRUE);
      gStyle->SetOptFit(0111);
      //predictLine->Draw("same");
      //logNormFunct->Draw("same");

      //FITTING HERE
      int status = m_currentHistoToFit->Fit("m_logNormFunct","EQ","");//,ZTPCPredict-0.5, ZTPCPredict+0.5);//,-5.0, 5.0);
      if(status == 0){
        TF1* fitResult = m_currentHistoToFit->GetFunction("m_logNormFunct"); 
        if(fitResult->GetParameter(0) != 0.0){
          if(a_isBTOF){
            if(histoIndex==0){
              m_betaGammaBTOF_DeDx_Data_Plus[0]->SetBinContent(betaGammaBin,fitResult->GetParameter(0));
              m_betaGammaBTOF_DeDx_Data_Plus[1]->SetBinContent(betaGammaBin,fitResult->GetParameter(1));
              m_betaGammaBTOF_DeDx_Data_Plus[2]->SetBinContent(betaGammaBin,fitResult->GetParameter(2));
              m_betaGammaBTOF_DeDx_Data_Plus[0]->SetBinError(betaGammaBin,fitResult->GetParError(0));
              m_betaGammaBTOF_DeDx_Data_Plus[1]->SetBinError(betaGammaBin,fitResult->GetParError(1));
              m_betaGammaBTOF_DeDx_Data_Plus[2]->SetBinError(betaGammaBin,fitResult->GetParError(2));
            }else{
              m_betaGammaBTOF_DeDx_Data_Minus[0]->SetBinContent(betaGammaBin,fitResult->GetParameter(0));
              m_betaGammaBTOF_DeDx_Data_Minus[1]->SetBinContent(betaGammaBin,fitResult->GetParameter(1));
              m_betaGammaBTOF_DeDx_Data_Minus[2]->SetBinContent(betaGammaBin,fitResult->GetParameter(2));
              m_betaGammaBTOF_DeDx_Data_Minus[0]->SetBinError(betaGammaBin,fitResult->GetParError(0));
              m_betaGammaBTOF_DeDx_Data_Minus[1]->SetBinError(betaGammaBin,fitResult->GetParError(1));
              m_betaGammaBTOF_DeDx_Data_Minus[2]->SetBinError(betaGammaBin,fitResult->GetParError(2));
            }
          }else{ // is ETOF
            if(histoIndex==0){
              m_betaGammaETOF_DeDx_Data_Plus[0]->SetBinContent(betaGammaBin,fitResult->GetParameter(0));
              m_betaGammaETOF_DeDx_Data_Plus[1]->SetBinContent(betaGammaBin,fitResult->GetParameter(1));
              m_betaGammaETOF_DeDx_Data_Plus[2]->SetBinContent(betaGammaBin,fitResult->GetParameter(2));
              m_betaGammaETOF_DeDx_Data_Plus[0]->SetBinError(betaGammaBin,fitResult->GetParError(0));
              m_betaGammaETOF_DeDx_Data_Plus[1]->SetBinError(betaGammaBin,fitResult->GetParError(1));
              m_betaGammaETOF_DeDx_Data_Plus[2]->SetBinError(betaGammaBin,fitResult->GetParError(2));
            }else{
              m_betaGammaETOF_DeDx_Data_Minus[0]->SetBinContent(betaGammaBin,fitResult->GetParameter(0));
              m_betaGammaETOF_DeDx_Data_Minus[1]->SetBinContent(betaGammaBin,fitResult->GetParameter(1));
              m_betaGammaETOF_DeDx_Data_Minus[2]->SetBinContent(betaGammaBin,fitResult->GetParameter(2));
              m_betaGammaETOF_DeDx_Data_Minus[0]->SetBinError(betaGammaBin,fitResult->GetParError(0));
              m_betaGammaETOF_DeDx_Data_Minus[1]->SetBinError(betaGammaBin,fitResult->GetParError(1));
              m_betaGammaETOF_DeDx_Data_Minus[2]->SetBinError(betaGammaBin,fitResult->GetParError(2));
            }
          }
        }
        //draw result
        if(a_saveImages){
          fittingCanvas->cd();
          m_currentHistoToFit->Draw("E");
          initialFunct->Draw("same"); 
          m_logNormFunct->Draw("same");
          gPad->SetLogx();
          gPad->SetLogy();
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),imageName.Data()));
        }
      } // end status 0
      delete m_currentHistoToFit;
      delete initialFunct;
      
    }// BetaGamma Bin Loop

  }//plus minus
  
  cout << "Done LogNormal Fitting BetaGamma dEdx" << endl;

  // now fit the means with a bichsel curve
  //TF1* betaGamma_dEdx_Fit = new TF1("betaGamma_dEdx_Fit",this,&ZFitter::dEdxBetaGammaOperator,0.001,100,7,"ZFitter","dEdxBetaGammaOperator");
  TF1* betaGamma_dEdx_Fit = new TF1("betaGamma_dEdx_Fit",
    "[0] * TMath::Power(1.0 + TMath::Power(1.0/x,2),[1]) * (TMath::Power(TMath::Abs(TMath::Log([2]*TMath::Power(x,2))),[3])+ [4] * TMath::Power(1.0+TMath::Power(1.0/x,2),[5])) - [6]",
    0.001,100);

  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    if(a_isBTOF){
      if(pmIndex == 0) m_currentHistoToFit = m_betaGammaBTOF_DeDx_Data_Plus[1];
      if(pmIndex == 1) m_currentHistoToFit = m_betaGammaBTOF_DeDx_Data_Minus[1]; 
    }else{
      if(pmIndex == 0) m_currentHistoToFit = m_betaGammaETOF_DeDx_Data_Plus[1];
      if(pmIndex == 1) m_currentHistoToFit = m_betaGammaETOF_DeDx_Data_Minus[1]; 
    }

      betaGamma_dEdx_Fit->SetParameter(0,1.56767);//
      betaGamma_dEdx_Fit->SetParameter(1,0.195761);
      betaGamma_dEdx_Fit->SetParameter(2,10.1956); //
      betaGamma_dEdx_Fit->SetParameter(3,0.181460);
      betaGamma_dEdx_Fit->SetParameter(4,2.25868);// 1.63
      betaGamma_dEdx_Fit->SetParameter(5,0.0103825);//
      betaGamma_dEdx_Fit->SetParameter(6,4.70897);
      for(int iii = 0; iii < 7; iii++){
        betaGamma_dEdx_Fit->SetParLimits(iii,0.8*betaGamma_dEdx_Fit->GetParameter(iii),1.2*betaGamma_dEdx_Fit->GetParameter(iii));
      }
      betaGamma_dEdx_Fit->SetNpx(1200);

      int status = m_currentHistoToFit->Fit("betaGamma_dEdx_Fit","E","",a_betaGammaFitStart,10);//add Q to quiet
      TF1* fitResult = m_currentHistoToFit->GetFunction("betaGamma_dEdx_Fit");
      if(fitResult){
        fitResult->SetRange(0.001,100);
      }
      
      //if(status!=0){
      //  cout << "Status: " << status << endl;
      //  cout << "\n\n\n\n\n\n BAD betaGamma_dEdx_Fit \n\n\n\n\n\n";
      //  continue;
      //}

      fittingCanvas->cd();
      //gPad->SetLogx();
      gPad->SetRightMargin(0.4);
      gPad->SetLogy();
      m_currentHistoToFit->Draw("E");
      if(fitResult) fitResult->Draw("same");
    
      if(a_isBTOF){
        if(pmIndex == 0){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitBTOFPlus.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
          if(fitResult){
            m_betaGammaBTOF_dEdx_FitFunction[0] = (TF1*) fitResult->Clone();
            m_betaGammaBTOF_dEdx_FitFunction[0]->SetName("bichselFitBTOFPlus");
          }
        }
        if(pmIndex == 1){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitBTOFMinus.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
          if(fitResult){
            m_betaGammaBTOF_dEdx_FitFunction[1] = (TF1*) fitResult->Clone();
            m_betaGammaBTOF_dEdx_FitFunction[1]->SetName("bichselFitBTOFMinus");
          }
        }
      }else{ //is ETOF
        if(pmIndex == 0){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitETOFPlus.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
          if(fitResult){
            m_betaGammaETOF_dEdx_FitFunction[0] = (TF1*) fitResult->Clone();
            m_betaGammaETOF_dEdx_FitFunction[0]->SetName("bichselFitETOFPlus");
          }
        }
        if(pmIndex == 1){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitETOFMinus.png",m_imagePreDir.c_str(), m_imgDirName.c_str()));
          if(fitResult){
            m_betaGammaETOF_dEdx_FitFunction[1] = (TF1*) fitResult->Clone();
            m_betaGammaETOF_dEdx_FitFunction[1]->SetName("bichselFitETOFMinus");
          }
        }
      }

  }




  delete fittingCanvas;
  return;

}










void ZFitter::fitEtaBetaGammaDeDx(string a_storageFileName,bool a_isBTOF, double a_betaGammaFitStart, bool a_saveImages){
  double lowBGCut = 0.5;
  //TGraphSmooth* smoother = new TGraphSmooth("normal");
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(m_canvasParameters[0],m_canvasParameters[1]);
  fittingCanvas->SetCanvasSize(m_canvasParameters[2],m_canvasParameters[3]);

  // loop over plus and then minus particles
  for(int histoIndex = 0; histoIndex < 2; histoIndex++){
    TH3I* histoPM = NULL;
    string histoNameStr = "";
    if(a_isBTOF){
      if(histoIndex==0){
        histoPM = m_eta_betaGammaBtof_dEdx_Plus;
        histoNameStr = "btof_plus";
      }else{
        histoPM = m_eta_betaGammaBtof_dEdx_Minus;
        histoNameStr = "btof_minus";
      }
    }else{
      if(histoIndex==0){
        histoPM = m_eta_betaGammaEtof_dEdx_Plus;
        histoNameStr = "etof_plus";
      }else{
        histoPM = m_eta_betaGammaEtof_dEdx_Minus;
        histoNameStr = "etof_minus";
      }
    }
    
    //loop over eta bins
    for(int etaBin = 1; etaBin <= histoPM->GetNbinsX(); etaBin++){
      for(int betaGammaBin = 1; betaGammaBin <= histoPM->GetNbinsY(); betaGammaBin++){
        cout << "Fitting dEdx for EtaBin: " << etaBin << " and betaGammaBin: " << betaGammaBin << endl;
        m_currentHistoToFit = histoPM->ProjectionZ(
             Form("betagamma_dEdx_%02d_%02d_%02d",etaBin,betaGammaBin,histoIndex), 
             etaBin, etaBin, betaGammaBin, betaGammaBin);
        double betaGamma = histoPM->GetYaxis()->GetBinCenter(betaGammaBin);   
        double counts = HistogramUtilities::sumOfContent(m_currentHistoToFit, -1);
        if( 250 > counts ){  // && !m_testing)  || (counts < 20)){
          delete m_currentHistoToFit;
          //cout << "Continuing: too few counts for solo fit partIndex:" << partIndex << " counts:" << counts << " at mom: " << momentum << endl; 
          continue;
        }

        int firstBinWithContent = HistogramUtilities::firstBinWithContent(m_currentHistoToFit);
        int lastBinWithContent = HistogramUtilities::lastBinWithContent(m_currentHistoToFit);
        double lowEdge = m_currentHistoToFit->GetBinLowEdge(firstBinWithContent);
        double highEdge = m_currentHistoToFit->GetBinLowEdge(lastBinWithContent+1);

        TH1D* truncHisto = HistogramUtilities::truncateHistogram(m_currentHistoToFit,lowEdge,highEdge);
        delete m_currentHistoToFit;
        m_currentHistoToFit = HistogramUtilities::rebinEqualByCDF(truncHisto,30,false,false);
        delete truncHisto;
      


        HistogramUtilities::makeSqrtNErrors(m_currentHistoToFit);
        //m_currentHistoToFit->Sumw2();
        m_currentHistoToFit->Scale(1.0,"width");  
        double integral = m_currentHistoToFit->Integral();
        double mean = m_currentHistoToFit->GetMean();
        double stdDev = m_currentHistoToFit->GetStdDev();
        double mu = 0.0;
        double sig = 0.0;      
        PhysMath::gausToLogNorm(mean,stdDev,mu,sig);

        TString imageName = TString::Format("dEdx_%s_eta_%02d_betagammma_%02d",histoNameStr.c_str(),etaBin,betaGammaBin);
        m_currentHistoToFit->SetName(imageName.Data());
        m_currentHistoToFit->SetTitle(Form(
                     "dE/dx for %s #eta=[%2.3f,%2.3f] #beta#gamma=[%2.3f,%2.3f]",
                     histoNameStr.c_str(),
                     histoPM->GetXaxis()->GetBinLowEdge(etaBin),
                     histoPM->GetXaxis()->GetBinLowEdge(etaBin+1),
                     histoPM->GetYaxis()->GetBinLowEdge(betaGammaBin),
                     histoPM->GetYaxis()->GetBinLowEdge(betaGammaBin+1)
                     ));
        m_currentHistoToFit->GetXaxis()->SetTitle("dE/dx");
        m_currentHistoToFit->GetYaxis()->SetTitle("dN/ddEdx");
        m_currentHistoToFit->SetStats(kTRUE);

        m_logNormFunct->SetParameter(0,counts);
        m_logNormFunct->SetParLimits(0,0.6*counts, 1.5*counts);
        m_logNormFunct->SetParameter(1,mu);
        m_logNormFunct->SetParLimits(1, 0.5*mu, 2.0*mu);
        m_logNormFunct->SetParameter(2,sig);
        m_logNormFunct->SetParLimits(2,0.05*sig,2.0*sig);
        m_logNormFunct->SetParNames("N","#mu","#sigma");

        TF1* initialFunct = (TF1*) m_logNormFunct->Clone();
        initialFunct->SetLineColor(kViolet); 
        m_currentHistoToFit->SetStats(kTRUE);
        gStyle->SetOptFit(0111);
        //predictLine->Draw("same");
        //logNormFunct->Draw("same");

        double lowBound = m_currentHistoToFit->GetBinLowEdge(1);
        double highBound = m_currentHistoToFit->GetBinLowEdge(m_currentHistoToFit->GetNbinsX()+1);

        if(betaGamma < 1.46){
          lowBound = pow(10.0, -1.418*log10(betaGamma)+0.487); 
        }
        if(betaGamma > 1.63) highBound = 3.32;

        //FITTING HERE
        int status = m_currentHistoToFit->Fit("m_logNormFunct","EIQ","",lowBound,highBound);//,ZTPCPredict-0.5, ZTPCPredict+0.5);//,-5.0, 5.0);
        if(status == 0){
          TF1* fitResult = m_currentHistoToFit->GetFunction("m_logNormFunct"); 
          if(fitResult->GetParameter(0) != 0.0){

            if(a_isBTOF){
              if(histoIndex==0){
                m_eta_betaGammaBTOF_DeDx_Data_Plus[0]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(0));
                m_eta_betaGammaBTOF_DeDx_Data_Plus[1]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(1));
                m_eta_betaGammaBTOF_DeDx_Data_Plus[2]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(2));
                m_eta_betaGammaBTOF_DeDx_Data_Plus[0]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(0));
                m_eta_betaGammaBTOF_DeDx_Data_Plus[1]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(1));
                m_eta_betaGammaBTOF_DeDx_Data_Plus[2]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(2));
              }else{
                m_eta_betaGammaBTOF_DeDx_Data_Minus[0]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(0));
                m_eta_betaGammaBTOF_DeDx_Data_Minus[1]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(1));
                m_eta_betaGammaBTOF_DeDx_Data_Minus[2]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(2));
                m_eta_betaGammaBTOF_DeDx_Data_Minus[0]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(0));
                m_eta_betaGammaBTOF_DeDx_Data_Minus[1]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(1));
                m_eta_betaGammaBTOF_DeDx_Data_Minus[2]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(2));
              }
            }else{ // is ETOF
              if(histoIndex==0){
                m_eta_betaGammaETOF_DeDx_Data_Plus[0]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(0));
                m_eta_betaGammaETOF_DeDx_Data_Plus[1]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(1));
                m_eta_betaGammaETOF_DeDx_Data_Plus[2]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(2));
                m_eta_betaGammaETOF_DeDx_Data_Plus[0]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(0));
                m_eta_betaGammaETOF_DeDx_Data_Plus[1]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(1));
                m_eta_betaGammaETOF_DeDx_Data_Plus[2]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(2));
              }else{
                m_eta_betaGammaETOF_DeDx_Data_Minus[0]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(0));
                m_eta_betaGammaETOF_DeDx_Data_Minus[1]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(1));
                m_eta_betaGammaETOF_DeDx_Data_Minus[2]->SetBinContent(etaBin,betaGammaBin,fitResult->GetParameter(2));
                m_eta_betaGammaETOF_DeDx_Data_Minus[0]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(0));
                m_eta_betaGammaETOF_DeDx_Data_Minus[1]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(1));
                m_eta_betaGammaETOF_DeDx_Data_Minus[2]->SetBinError(etaBin,betaGammaBin,fitResult->GetParError(2));
              }
            }
          }
          //draw result
          if(a_saveImages){
            fittingCanvas->cd();
            m_currentHistoToFit->Draw("E");
            initialFunct->Draw("same"); 
            m_logNormFunct->Draw("same");
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetRightMargin(0.4);
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/%s.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),imageName.Data()));
          }
        } // end status 0

        delete m_currentHistoToFit;
        delete initialFunct;
        

      }// BetaGamma Bin Loop
    }// Eta Bin Loop


      
   cout << "Making the Smoothed Parameter Plots for BetaGamma DeDx fitting" << endl;

    if(a_isBTOF){
      if(histoIndex==0){
        vector<TGraph* > NGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaBTOF_DeDx_Data_Plus[0], "eta", "#eta");
        vector<TGraph* > MuGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaBTOF_DeDx_Data_Plus[1], "eta", "#eta");
        vector<TGraph* > SigGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaBTOF_DeDx_Data_Plus[2], "eta", "#eta");
        for(int horIndex = 0; horIndex < (int) NGraph.size(); horIndex++){
          cout << "Smoothing " << NGraph[horIndex]->GetName() << endl;
          if(NGraph[horIndex]->GetN() == 0){
            cout << "    horIndex: " << horIndex << " has no points " << endl;
            continue;
          }
          cout << "  BEFORE:" << endl;
          cout << "    Ptr: " << NGraph[horIndex] << "  NumPoints: " << NGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << MuGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << SigGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;

          HistogramUtilities::truncateTGraph(NGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(MuGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(SigGraph[horIndex],lowBGCut, 1e20);

          cout << "  TRUNCATION:" << endl;
          cout << "    Ptr: " << NGraph[horIndex] << "  NumPoints: " << NGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << MuGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << SigGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;

          TGraphSmooth* N_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Mu_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Sig_smoother = new TGraphSmooth("normal");

          NGraph[horIndex] = N_smoother->SmoothSuper(NGraph[horIndex]);
          MuGraph[horIndex] = Mu_smoother->SmoothSuper(MuGraph[horIndex]);
          SigGraph[horIndex] = Sig_smoother->SmoothSuper(SigGraph[horIndex]);
          //delete N_smoother;
          //delete Mu_smoother;
          //delete Sig_smoother;
          cout << "  SMOOTH:" << endl;
          cout << "    Ptr: " << NGraph[horIndex] << "  NumPoints: " << NGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << MuGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << SigGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;          




        } 
        vector< TGraph* > tempGraph;
        m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData.push_back(tempGraph);
        m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData.push_back(tempGraph);
        m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData.push_back(tempGraph);
        for(int iii = 0; iii < (int) NGraph.size(); iii++){
          m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[0].push_back(NGraph[iii]);
          m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[1].push_back(MuGraph[iii]);
          m_eta_betaGammaBTOF_DeDx_Plus_SmoothFitData[2].push_back(SigGraph[iii]);
        }
                  
      }else{
        vector<TGraph* > NGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaBTOF_DeDx_Data_Minus[0], "eta", "#eta");
        vector<TGraph* > MuGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaBTOF_DeDx_Data_Minus[1], "eta", "#eta");
        vector<TGraph* > SigGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaBTOF_DeDx_Data_Minus[2], "eta", "#eta");
        for(int horIndex = 0; horIndex < (int) NGraph.size(); horIndex++){
          if(NGraph[horIndex]->GetN() == 0){
            cout << "horIndex: " << horIndex << " has no points " << endl;
            continue;
          }
           cout << "  BEFORE:" << endl;
          cout << "    Ptr: " << NGraph[horIndex] << "  NumPoints: " << NGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << MuGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << SigGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;

          HistogramUtilities::truncateTGraph(NGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(MuGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(SigGraph[horIndex],lowBGCut, 1e20);

          cout << "  TRUNCATION:" << endl;
          cout << "    Ptr: " << NGraph[horIndex] << "  NumPoints: " << NGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << MuGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << SigGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;

          TGraphSmooth* N_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Mu_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Sig_smoother = new TGraphSmooth("normal");

          NGraph[horIndex] = N_smoother->SmoothSuper(NGraph[horIndex]);
          MuGraph[horIndex] = Mu_smoother->SmoothSuper(MuGraph[horIndex]);
          SigGraph[horIndex] = Sig_smoother->SmoothSuper(SigGraph[horIndex]);
          //delete N_smoother;
          //delete Mu_smoother;
          //delete Sig_smoother;
          cout << "  SMOOTH:" << endl;
          cout << "    Ptr: " << NGraph[horIndex] << "  NumPoints: " << NGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << MuGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;
          cout << "    Ptr: " << SigGraph[horIndex] << "  NumPoints: " << MuGraph[horIndex]->GetN() << endl;

     
        }   
        
        cout << "Sizes: " << NGraph.size() << " " << MuGraph.size() << " " << SigGraph.size() << endl;
        vector< TGraph* > tempGraph;
        m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData.push_back(tempGraph);
        m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData.push_back(tempGraph);
        m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData.push_back(tempGraph);
        for(int iii = 0; iii < (int) NGraph.size(); iii++){
          m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[0].push_back(NGraph[iii]);
          m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[1].push_back(MuGraph[iii]);
          m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData[2].push_back(SigGraph[iii]);
        }


         // m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData.push_back(NGraph);
         // m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData.push_back(MuGraph);
         // m_eta_betaGammaBTOF_DeDx_Minus_SmoothFitData.push_back(SigGraph);  
      }
    }else{ // is ETOF
      if(histoIndex==0){
        vector<TGraph* > NGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaETOF_DeDx_Data_Plus[0], "eta", "#eta");
        vector<TGraph* > MuGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaETOF_DeDx_Data_Plus[1], "eta", "#eta");
        vector<TGraph* > SigGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaETOF_DeDx_Data_Plus[2], "eta", "#eta");
        for(int horIndex = 0; horIndex < (int) NGraph.size(); horIndex++){
          if(NGraph[horIndex]->GetN() == 0){
            cout << "horIndex: " << horIndex << " has no points " << endl;
            continue;
          }
          HistogramUtilities::truncateTGraph(NGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(MuGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(SigGraph[horIndex],lowBGCut, 1e20);
          TGraphSmooth* N_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Mu_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Sig_smoother = new TGraphSmooth("normal");

          NGraph[horIndex] = N_smoother->SmoothSuper(NGraph[horIndex]);
          MuGraph[horIndex] = Mu_smoother->SmoothSuper(MuGraph[horIndex]);
          SigGraph[horIndex] = Sig_smoother->SmoothSuper(SigGraph[horIndex]);
          delete N_smoother;
          delete Mu_smoother;
          delete Sig_smoother;       

        }   
          m_eta_betaGammaETOF_DeDx_Plus_SmoothFitData.push_back(NGraph);
          m_eta_betaGammaETOF_DeDx_Plus_SmoothFitData.push_back(MuGraph);
          m_eta_betaGammaETOF_DeDx_Plus_SmoothFitData.push_back(SigGraph);                
      }else{
        vector<TGraph* > NGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaETOF_DeDx_Data_Minus[0], "eta", "#eta");
        vector<TGraph* > MuGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaETOF_DeDx_Data_Minus[1], "eta", "#eta");
        vector<TGraph* > SigGraph =HistogramUtilities::VerticalProjectionGraphsNoErrors(
                               m_eta_betaGammaETOF_DeDx_Data_Minus[2], "eta", "#eta");
        for(int horIndex = 0; horIndex < (int) NGraph.size(); horIndex++){
          if(NGraph[horIndex]->GetN() == 0){
            cout << "horIndex: " << horIndex << " has no points " << endl;
            continue;
          }
          HistogramUtilities::truncateTGraph(NGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(MuGraph[horIndex],lowBGCut, 1e20);
          HistogramUtilities::truncateTGraph(SigGraph[horIndex],lowBGCut, 1e20);
          TGraphSmooth* N_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Mu_smoother = new TGraphSmooth("normal");
          TGraphSmooth* Sig_smoother = new TGraphSmooth("normal");

          NGraph[horIndex] = N_smoother->SmoothSuper(NGraph[horIndex]);
          MuGraph[horIndex] = Mu_smoother->SmoothSuper(MuGraph[horIndex]);
          SigGraph[horIndex] = Sig_smoother->SmoothSuper(SigGraph[horIndex]);
          delete N_smoother;
          delete Mu_smoother;
          delete Sig_smoother;        
        }     
          m_eta_betaGammaETOF_DeDx_Minus_SmoothFitData.push_back(NGraph);
          m_eta_betaGammaETOF_DeDx_Minus_SmoothFitData.push_back(MuGraph);
          m_eta_betaGammaETOF_DeDx_Minus_SmoothFitData.push_back(SigGraph);
      }
    }

  }//plus minus

  cout << "Done Fitting Eta BetaGamma dEdx" << endl;
  
  // now fit the means with a bichsel curve
  //TF1* betaGamma_dEdx_Fit = new TF1("betaGamma_dEdx_Fit",this,&ZFitter::dEdxBetaGammaOperator,0.001,100,7,"ZFitter","dEdxBetaGammaOperator");



  TH2D* current2DHisto = NULL;
  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    if(a_isBTOF){
      if(pmIndex == 0) current2DHisto = m_eta_betaGammaBTOF_DeDx_Data_Plus[1];
      if(pmIndex == 1) current2DHisto = m_eta_betaGammaBTOF_DeDx_Data_Minus[1]; 
    }else{
      if(pmIndex == 0) current2DHisto = m_eta_betaGammaETOF_DeDx_Data_Plus[1];
      if(pmIndex == 1) current2DHisto = m_eta_betaGammaETOF_DeDx_Data_Minus[1]; 
    }

    for(int etaBin = 1; etaBin <= current2DHisto->GetNbinsX(); etaBin++){
      if(a_isBTOF) m_eta_betaGammaBTOF_dEdx_FitFunction[pmIndex].push_back(NULL);
      if(!a_isBTOF) m_eta_betaGammaETOF_dEdx_FitFunction[pmIndex].push_back(NULL);

      m_currentHistoToFit = current2DHisto->ProjectionY(Form("dEdxProj_%02d", etaBin),
           etaBin, etaBin);
      
      int firstBinWithContent = HistogramUtilities::firstBinWithContent(m_currentHistoToFit);
      int lastBinWithContent = HistogramUtilities::lastBinWithContent(m_currentHistoToFit);
      double lowEdge = m_currentHistoToFit->GetBinLowEdge(firstBinWithContent);
      double highEdge = m_currentHistoToFit->GetBinLowEdge(lastBinWithContent+1);


      /*
      TF1* betaGamma_dEdx_Fit = new TF1("betaGamma_dEdx_Fit",
        "[0] * pow(1.0 + pow(1.0/x,2),[1]) * (pow(TMath::Abs(log([2]*pow(x,2))),[3])+ [4] * pow(1.0+pow(1.0/x,2),[5])) - [6]",
        0.001,100);
      betaGamma_dEdx_Fit->SetParameter(0,1.56767);//
      betaGamma_dEdx_Fit->SetParameter(1,0.195761);
      betaGamma_dEdx_Fit->SetParameter(2,10.1956); //
      betaGamma_dEdx_Fit->SetParameter(3,0.181460);
      betaGamma_dEdx_Fit->SetParameter(4,2.25868);// 1.63
      betaGamma_dEdx_Fit->SetParameter(5,0.0103825);//
      betaGamma_dEdx_Fit->SetParameter(6,4.70897);

      betaGamma_dEdx_Fit->SetParLimits(0,0.0,10);//
      betaGamma_dEdx_Fit->SetParLimits(1,0.0,10);
      betaGamma_dEdx_Fit->SetParLimits(2,0.0,50); //
      betaGamma_dEdx_Fit->SetParLimits(3,0.0,10);
      betaGamma_dEdx_Fit->SetParLimits(4,0.0,20);// 1.63
      betaGamma_dEdx_Fit->SetParLimits(5,0.0,5);//
      betaGamma_dEdx_Fit->SetParLimits(6,-20,50);
      //for(int iii = 0; iii < 7; iii++){
      //  betaGamma_dEdx_Fit->SetParLimits(iii,0.8*betaGamma_dEdx_Fit->GetParameter(iii),1.2*betaGamma_dEdx_Fit->GetParameter(iii));
      //}

      */

      TF1* betaGamma_dEdx_Fit = PhysMath::getBichselFunctionBetaGamma(Form("betaGamma_dEdx_Fit_%02d_%02d",pmIndex,etaBin));
      //betaGamma_dEdx_Fit->SetParameter(0,0.076);// A
      //betaGamma_dEdx_Fit->SetParameter(1,10.632);// B
      //betaGamma_dEdx_Fit->SetParameter(2,1.34E-5); //C
      //betaGamma_dEdx_Fit->SetParameter(3,1.863); //D
      //betaGamma_dEdx_Fit->SetParameter(4,1.948);// E

      betaGamma_dEdx_Fit->SetParameter(0,4.36870e+00);
      betaGamma_dEdx_Fit->SetParameter(1,1.27450e+00);
      betaGamma_dEdx_Fit->SetParameter(2,9.99997e-02);
      betaGamma_dEdx_Fit->SetParameter(3,1.74160e-01);
      betaGamma_dEdx_Fit->SetParameter(4,2.95910e-02);

      betaGamma_dEdx_Fit->SetParLimits(0,0,10000);// A
      betaGamma_dEdx_Fit->SetParLimits(1,0,100);// B
      betaGamma_dEdx_Fit->SetParLimits(2,0,100); //C
      betaGamma_dEdx_Fit->SetParLimits(3,0,100); //D
      betaGamma_dEdx_Fit->SetParLimits(4,0,100);// E

      betaGamma_dEdx_Fit->SetNpx(1200);

      TF1* initFunct = (TF1*) betaGamma_dEdx_Fit->Clone();
      initFunct->SetName("initFunct");
      initFunct->SetLineColor(kYellow);

      int status = m_currentHistoToFit->Fit(betaGamma_dEdx_Fit,"W","",lowEdge,highEdge);//add Q to quiet
      TF1* fitResult = m_currentHistoToFit->GetFunction(Form("betaGamma_dEdx_Fit_%02d_%02d",pmIndex,etaBin));
      if(fitResult){
        fitResult->SetRange(0.001,100);
      }

      //if(status!=0){
      //  cout << "Status: " << status << endl;
      //  cout << "\n\n\n\n\n\n BAD betaGamma_dEdx_Fit \n\n\n\n\n\n";
      //  continue;
      //}
      m_currentHistoToFit->SetMarkerStyle(20);
      m_currentHistoToFit->SetMarkerColor(kBlack);
      m_currentHistoToFit->SetMarkerSize(0.3);




      fittingCanvas->cd();
      //gPad->SetLogx();
      gPad->SetRightMargin(0.4);
      gPad->SetLogy();
      gPad->SetLogx();
      m_currentHistoToFit->Draw("PE");
      m_currentHistoToFit->GetXaxis()->SetRangeUser(lowEdge,highEdge);
      initFunct->Draw("same");
      if(fitResult) fitResult->Draw("same");

    
      if(a_isBTOF){
        if(pmIndex == 0){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitBTOFPlus_EtaIndex_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),etaBin-1));
          if(fitResult){
            m_eta_betaGammaBTOF_dEdx_FitFunction[0][etaBin-1] = (TF1*) fitResult->Clone();
            m_eta_betaGammaBTOF_dEdx_FitFunction[0][etaBin-1]->SetName(Form("bichselFitBTOFPlus_EtaIndex_%02d",etaBin-1));
            
            /*TF1* tempFit = new TF1("exp_betaGamma_dEdx_Fit",
              "exp([0] * TMath::Power(1.0 + TMath::Power(1.0/x,2),[1]) * (TMath::Power(TMath::Abs(TMath::Log([2]*TMath::Power(x,2))),[3])+ [4] * TMath::Power(1.0+TMath::Power(1.0/x,2),[5])) - [6])",
              0.001,100);
              tempFit->SetParameter(0,betaGamma_dEdx_Fit->GetParameter(0));
              tempFit->SetParameter(1,betaGamma_dEdx_Fit->GetParameter(1));
              tempFit->SetParameter(2,betaGamma_dEdx_Fit->GetParameter(2));
              tempFit->SetParameter(3,betaGamma_dEdx_Fit->GetParameter(3));
              tempFit->SetParameter(4,betaGamma_dEdx_Fit->GetParameter(4));
              tempFit->SetParameter(5,betaGamma_dEdx_Fit->GetParameter(5));
              tempFit->SetParameter(6,betaGamma_dEdx_Fit->GetParameter(6));*/
 
              


          }
        }
        if(pmIndex == 1){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitBTOFMinus_EtaIndex_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),etaBin-1));
          if(fitResult){
            m_eta_betaGammaBTOF_dEdx_FitFunction[1][etaBin-1] = (TF1*) fitResult->Clone();
            m_eta_betaGammaBTOF_dEdx_FitFunction[1][etaBin-1]->SetName(Form("bichselFitBTOFMinus_EtaIndex_%02d",etaBin-1));
          }
        }
      }else{ //is ETOF
        if(pmIndex == 0){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitETOFPlus_EtaIndex_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),etaBin-1));
          if(fitResult){
            m_eta_betaGammaETOF_dEdx_FitFunction[0][etaBin-1] = (TF1*) fitResult->Clone();
            m_eta_betaGammaETOF_dEdx_FitFunction[0][etaBin-1]->SetName(Form("bichselFitETOFPlus_EtaIndex_%02d",etaBin-1));
          }
        }
        if(pmIndex == 1){
          fittingCanvas->SaveAs(Form("%s/%s/eta_betaGamma_dEdx/BichselFitETOFMinus_EtaIndex_%02d.png",m_imagePreDir.c_str(), m_imgDirName.c_str(),etaBin-1));
          if(fitResult){
            m_eta_betaGammaETOF_dEdx_FitFunction[1][etaBin-1] = (TF1*) fitResult->Clone();
            m_eta_betaGammaETOF_dEdx_FitFunction[1][etaBin-1]->SetName(Form("bichselFitETOFMinus_EtaIndex_%02d",etaBin-1));
          }
        }
      }


    }//eta bin loop

  }


  TFile file(a_storageFileName.c_str(),"RECREATE");
  for(int etaIndex = 0; etaIndex < m_eta_betaGammaETOF_dEdx_FitFunction[0].size(); etaIndex++){
    cout << " attempting to write : " << m_eta_betaGammaBTOF_dEdx_FitFunction[0][etaIndex] << " eta " << etaIndex << " Plus" << endl;
    HistogramUtilities::ConditionalWrite((TF1*) m_eta_betaGammaBTOF_dEdx_FitFunction[0][etaIndex]->Clone());
  }
  for(int etaIndex = 0; etaIndex < m_eta_betaGammaETOF_dEdx_FitFunction[1].size(); etaIndex++){
    cout << " attempting to write : " << m_eta_betaGammaBTOF_dEdx_FitFunction[1][etaIndex] << " eta " << etaIndex << " Minus" << endl;
    HistogramUtilities::ConditionalWrite((TF1*) m_eta_betaGammaBTOF_dEdx_FitFunction[1][etaIndex]->Clone());
  }
  //file.Close();
  delete fittingCanvas;
  return;

}

void ZFitter::loadEtaBetaGammaDeDxFile(string a_storageFileName, int a_lowIndexToLoad, int a_highIndexToLoad){
  cout << " ZFitter::loadEtaBetaGammaDeDxFile(" << a_storageFileName  << ")" << endl;
  TFile file(a_storageFileName.c_str(),"READ");
  cout << " Getting indicies from " << m_eta_betaGammaBtof_dEdx_Plus << endl;
  for(int iii = 0; iii < a_lowIndexToLoad; iii++){
    m_eta_betaGammaBTOF_dEdx_FitFunction[0].push_back((TF1*) NULL);
    m_eta_betaGammaBTOF_dEdx_FitFunction[1].push_back((TF1*) NULL);    
  }
  for(int etaIndex = a_lowIndexToLoad; etaIndex <= a_highIndexToLoad; etaIndex++){
    cout << " attempting to load etaBetaGamma with etaIndex: " << etaIndex <<  " with name " << Form("betagamma_FitData/bichselFitBTOFPlus_EtaIndex_%02d",etaIndex) << endl;
    TF1* plusFit = (TF1*) file.Get(Form("betagamma_FitData/bichselFitBTOFPlus_EtaIndex_%02d",etaIndex));
    TF1* minusFit = (TF1*) file.Get(Form("betagamma_FitData/bichselFitBTOFMinus_EtaIndex_%02d",etaIndex));
    if(plusFit) m_eta_betaGammaBTOF_dEdx_FitFunction[0].push_back((TF1*) plusFit->Clone());
    else        m_eta_betaGammaBTOF_dEdx_FitFunction[0].push_back((TF1*) NULL);
    if(minusFit) m_eta_betaGammaBTOF_dEdx_FitFunction[1].push_back((TF1*) minusFit->Clone());
    else         m_eta_betaGammaBTOF_dEdx_FitFunction[1].push_back((TF1*) NULL);
    cout << " pushing back: " << plusFit << " " << minusFit << endl;
  }
  for(int etaIndex = a_highIndexToLoad + 1; etaIndex < m_eta_betaGammaBtof_dEdx_Plus->GetNbinsX(); etaIndex++){
    m_eta_betaGammaBTOF_dEdx_FitFunction[0].push_back((TF1*) NULL);
    m_eta_betaGammaBTOF_dEdx_FitFunction[1].push_back((TF1*) NULL);  
  }



  for(int etaIndex = 0; etaIndex < m_eta_betaGammaBTOF_dEdx_FitFunction[0].size(); etaIndex++){
    cout <<  "Loaded EtaIndex: " << etaIndex << "  ptr: " << m_eta_betaGammaBTOF_dEdx_FitFunction[0][etaIndex] << endl;
  }
  for(int etaIndex = 0; etaIndex < m_eta_betaGammaBTOF_dEdx_FitFunction[1].size(); etaIndex++){
    cout <<  "Loaded EtaIndex: " << etaIndex << "  ptr: " << m_eta_betaGammaBTOF_dEdx_FitFunction[1][etaIndex] << endl;
  }

  //file.Close();
  return;
}






