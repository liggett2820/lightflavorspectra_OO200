#ifndef _HYBRID_FEEDDOWN_MAKER_CXX_
#define _HYBRID_FEEDDOWN_MAKER_CXX_

#include "../headers/HybridFeeddownMaker.h"

#ifdef _MULTI_THREADING_
  std::mutex mtx;
#endif


HybridFeeddownMaker::HybridFeeddownMaker(){
  m_partInfo = new ParticleInfo();
  m_feeddown_particle_names[0] = "Lambda";          m_feeddown_particle_masses[0] = 1.115683;
  m_feeddown_particle_names[1] = "Antilambda";      m_feeddown_particle_masses[1] = 1.115683;
  m_feeddown_particle_names[2] = "SigmaPlus";       m_feeddown_particle_masses[2] = 1.18937;
  m_feeddown_particle_names[3] = "AntiSigmaMinus";  m_feeddown_particle_masses[3] = 1.18937;
  m_feeddown_particle_names[4] = "XiMinus";         m_feeddown_particle_masses[4] = 1.32171;
  m_feeddown_particle_names[5] = "AntiXiPlus";      m_feeddown_particle_masses[5] = 1.32171;
  m_feeddown_particle_names[6] = "XiZero";          m_feeddown_particle_masses[6] = 1.31486;
  m_feeddown_particle_names[7] = "OmegaMinus";      m_feeddown_particle_masses[7] = 1.67245;
  m_feeddown_particle_names[8] = "AntiOmegaPlus";   m_feeddown_particle_masses[8] = 1.67245;
  m_feeddown_particle_names[9] = "SigmaZero";       m_feeddown_particle_masses[9] = 1.192642;
  m_feeddown_particle_names[10] = "AntiSigmaZero";  m_feeddown_particle_masses[10] = 1.192642;
  m_feeddown_particle_names[11] = "K0_S";           m_feeddown_particle_masses[11] = 0.497614;
  m_feeddown_particle_names[12] = "Eta";            m_feeddown_particle_masses[12] = 0.547862;


  m_feeddown_particle_symbols[0] = "#Lambda";
  m_feeddown_particle_symbols[1] = "#bar{#Lambda}";
  m_feeddown_particle_symbols[2] = "#Sigma^{+}";
  m_feeddown_particle_symbols[3] = "#Sigma^{-}";
  m_feeddown_particle_symbols[4] = "#Xi^{-}";
  m_feeddown_particle_symbols[5] = "#Xi^{+}";
  m_feeddown_particle_symbols[6] = "#Xi^{0}";
  m_feeddown_particle_symbols[7] = "#Omega^{-}";
  m_feeddown_particle_symbols[8] = "#Omega^{+}";
  m_feeddown_particle_symbols[9] = "#Sigma^{0}";
  m_feeddown_particle_symbols[10] = "#bar{#Sigma}^{0}";
  m_feeddown_particle_symbols[11] = "K^{0}_{S}";
  m_feeddown_particle_symbols[12] = "#eta";

  m_response_nTuple = new TNtuple("feeddownDecayTuple","Parent and Daughters for Feeddown","ID_feed:y_feed:mtm0_feed:ID_d1:y_d1:mtm0_d1:ID_d2:y_d2:mTm0_d2:ID_d3:y_d3:mTm0_d3:ID_d4:y_d4:mTm0_d4");

  m_response_structure    = NULL;
  m_currentDataGraphToFit = NULL;
  m_currentSimGraphToFit  = NULL;
  m_currentDataHistoToFit = NULL;
  m_currentSimHistoToFit  = NULL;
  m_matchingParametersHisto = NULL;

  m_currentFeedDownIndex = -1;
  m_doNTuple = false;
  m_symmeterizeOutputs = true;
  m_imageDir = "";
  m_use_DCA_cut = false;
  m_DCA_Cut = 999;

  for(int feedDownIndex = 0; feedDownIndex < 20; feedDownIndex++){

    m_urqmd_dNdy[feedDownIndex] = NULL;
    m_urqmd_dNdy_scaled[feedDownIndex] = NULL;
    m_mappedUrQMD_dNdy[feedDownIndex]  = NULL;

    for(int centIndex = 0; centIndex < 16; centIndex++){
      m_urqmd_spectra[feedDownIndex][centIndex]             = NULL;
      m_urqmd_yields[feedDownIndex][centIndex]              = NULL;
      m_measuredSpectra[feedDownIndex][centIndex]           = NULL;
      m_mappedUrQMDSpectra[feedDownIndex][centIndex]        = NULL;
      m_unmappedUrQMDSpectra[feedDownIndex][centIndex]      = NULL;
      m_resampled_feeddown_yields[feedDownIndex][centIndex] = NULL;
      m_resampled_feeddown_spectra[feedDownIndex][centIndex] = NULL;
      m_measuredSpectra_Fits[feedDownIndex][centIndex]      = NULL;
      m_resampled_feeddown_dNdy[feedDownIndex][centIndex]   = NULL;
      m_mappedUrQMD_dNdyByCent[feedDownIndex][centIndex]    = NULL;
      m_urqmd_dNdy_scaled_byCent[feedDownIndex][centIndex]    = NULL;
      for(int partIndex = 0; partIndex < 9; partIndex++){
        m_produced_feeddown_dNdy[partIndex][feedDownIndex][centIndex][0] = NULL;
        m_produced_feeddown_dNdy[partIndex][feedDownIndex][centIndex][1] = NULL;
        m_produced_feeddown_DCA[partIndex][feedDownIndex][centIndex][0]  = NULL;
        m_produced_feeddown_DCA[partIndex][feedDownIndex][centIndex][1]  = NULL;
        m_produced_feeddown_dNdy_DCA[partIndex][feedDownIndex][centIndex][0] = NULL;
        m_produced_feeddown_dNdy_DCA[partIndex][feedDownIndex][centIndex][1] = NULL;
        m_produced_feeddown_yields[partIndex][feedDownIndex][centIndex][0]   = NULL;
        m_produced_feeddown_spectra[partIndex][feedDownIndex][centIndex][0]  = NULL;
        m_produced_feeddown_yields[partIndex][feedDownIndex][centIndex][1]   = NULL;
        m_produced_feeddown_spectra[partIndex][feedDownIndex][centIndex][1]  = NULL;
        m_produced_total_feeddown_yields[partIndex][centIndex][0]   = NULL;
        m_produced_total_feeddown_yields[partIndex][centIndex][1]   = NULL;
        m_produced_total_feeddown_spectra[partIndex][centIndex][0]  = NULL;
        m_produced_total_feeddown_spectra[partIndex][centIndex][1]  = NULL;
        m_produced_total_feeddown_dNdy[partIndex][centIndex][0]     = NULL;
        m_produced_total_feeddown_dNdy[partIndex][centIndex][1]     = NULL;
        m_produced_total_feeddown_dNdy_DCA[partIndex][centIndex][0] = NULL;
        m_produced_total_feeddown_dNdy_DCA[partIndex][centIndex][1] = NULL;
        if(partIndex < 4){
          if(partIndex < 2)  m_matchingParameters[feedDownIndex][centIndex][partIndex] = 1.0;
          if(partIndex == 3 || partIndex == 2 )  m_matchingParameters[feedDownIndex][centIndex][partIndex] = 0.0;
        }

        for(int iii = 0; iii < 50; iii++){
          m_produced_total_feeddown_spectraFits[partIndex][centIndex][0][iii] = NULL;
          m_produced_total_feeddown_spectraFits[partIndex][centIndex][1][iii] = NULL;
          m_produced_total_feeddown_spectraFits[partIndex][centIndex][0][iii] = NULL;
          m_produced_total_feeddown_spectraFitsCovMat[partIndex][centIndex][0][iii] = NULL;
          m_produced_total_feeddown_spectraFits[partIndex][centIndex][1][iii] = NULL;
          m_produced_total_feeddown_spectraFitsCovMat[partIndex][centIndex][1][iii] = NULL;
        }
      }
    }
  }

  m_matchingParametersHisto = NULL;


}

//====================================================================================================================================================================================================================================================================================================================

void HybridFeeddownMaker::loadUrQMDSpectraFile(string a_filename, bool a_oldNames){


  cout << "HybridFeeddownMaker::loadUrQMDSpectraFile   Loading " << a_filename << endl;
  TFile* inFile = new TFile(a_filename.c_str(), "READ");
  cout << "    File Ptr: " << inFile << endl;
  inFile->cd();

  string particleNames[13] = {"LAMBDA","ANTILAMBDA","SIGMA_Plus","ANTISIGMA_Minus","XI_Minus","ANTIXI_Plus","XI_0","OMEGA_Minus","ANTIOMEGA_Plus"," SIGMA_0","ANTISIGMA_0","KAON_0_S","ETA"};
  for(int feedDownIndex = 0; feedDownIndex < 13; feedDownIndex++){
    m_urqmd_dNdy[feedDownIndex] = (TH2D*) inFile->Get(Form("%s/dNdy_%s_convToDouble",particleNames[feedDownIndex].c_str(),particleNames[feedDownIndex].c_str()));
    cout << " Looading " << particleNames[feedDownIndex]  << "   Ptr: " << m_urqmd_dNdy[feedDownIndex] << endl;
    if(m_urqmd_dNdy[feedDownIndex]){
      m_mappedUrQMD_dNdy[feedDownIndex] = (TH2D*) m_urqmd_dNdy[feedDownIndex]->Clone();
      m_mappedUrQMD_dNdy[feedDownIndex]->GetYaxis()->SetTitle("Centrality Bin of Measurment");
      HistogramUtilities::zeroBins2D(m_mappedUrQMD_dNdy[feedDownIndex]);
      m_urqmd_dNdy_scaled[feedDownIndex] = (TH2D*) m_urqmd_dNdy[feedDownIndex]->Clone();
      m_urqmd_dNdy_scaled[feedDownIndex]->GetYaxis()->SetTitle("Centrality Bin of Simulation");
      string newName = m_urqmd_dNdy[feedDownIndex]->GetName();
      newName += "_urqmd";
      string scaledName = m_urqmd_dNdy[feedDownIndex]->GetName();
      scaledName += "_urqmd_scaled";
      m_urqmd_dNdy[feedDownIndex]->SetName(newName.c_str());
      m_urqmd_dNdy_scaled[feedDownIndex]->SetName(scaledName.c_str());
    }
    for(int centIndex = 0; centIndex < 16; centIndex++){
      if(!a_oldNames) m_urqmd_spectra[feedDownIndex][centIndex] = (TH2D*) inFile->Get(Form("%s/yield_%s_Cent%02d_convToDouble",particleNames[feedDownIndex].c_str(),particleNames[feedDownIndex].c_str(),centIndex));
      else            m_urqmd_spectra[feedDownIndex][centIndex] = (TH2D*) inFile->Get(Form("%s/Spectra for %s Cent %d",particleNames[feedDownIndex].c_str(),particleNames[feedDownIndex].c_str(),centIndex));

      cout << "    Cent: " << centIndex << "   Ptr: " << m_urqmd_spectra[feedDownIndex][centIndex] << endl;
    }
  }

}

//====================================================================================================================================================================================================================================================================================================================

void HybridFeeddownMaker::setResponseStructure(int a_numRapBins, double a_lowRapRange, double a_highRapRange, double a_numMtM0Bins, double a_lowMtM0Range, double a_highMtM0Range){
  m_response_structure = new TH2D( "responseStructure", "responseStructure", a_numRapBins, a_lowRapRange, a_highRapRange,  a_numMtM0Bins, a_lowMtM0Range, a_highMtM0Range);
  m_response_structure_1D = new TH1D("responseStructure_1D", "responseStructure1D",a_numRapBins, a_lowRapRange, a_highRapRange);
}

//====================================================================================================================================================================================================================================================================================================================

double HybridFeeddownMaker::calcScalingFactor(double a_par0, double a_par1, double a_mTm0, double a_mTm0_width, double a_mass){
  //return a_par0; // [0]
  //return a_par0/a_par1; // [0]*binWidth_old/binWidth_new
  //return a_par0*PhysMath::mT_from_mTm0_m(a_mTm0,a_mass)/(a_par1* (a_par1 * a_mTm0 + a_mass));  // [0] * binWidth_old * mT_old/(binWidth_new * mT_new)
  //return a_par0*(a_mTm0/a_par1 + a_mass)/(a_par1 * (a_mTm0 + a_mass));

  // SEARCH FOR  "Matching Parameter 0 Scaling Safety"  if you ever change this, that code must be edited for proper function too
  return   a_par0*(a_mTm0 + a_mass)*a_mTm0_width/(a_par1*a_mTm0 + a_mass); // still needs to be divided by bin width later
};


//====================================================================================================================================================================================================================================================================================================================

void HybridFeeddownMaker::loadDataSpectra(int a_feeddownParticleId, int a_centIndex, TGraphAsymmErrors* a_graph, bool a_mTm0Space){
  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout << "Call to HybridFeeddownMaker::loadDataSpectra(" << a_feeddownParticleId << "," << a_centIndex << "," << a_graph << "," << a_mTm0Space << ")" << endl;
  #endif
  if(!a_graph){
    cout << "ERROR: HybridFeeddownMaker::loadDataSpectra( " << a_feeddownParticleId << ", " <<  a_centIndex << ", " << a_graph << ", " << a_mTm0Space << ")   passed a null pointer!!!!!! " << endl;
    cout << "MEASURED SPECTRA NOT LOADED!!!!" << endl;
    return;
  }
  int foundFeedDownIndex = a_feeddownParticleId;
  /*for(int feedDownIndex = 0; feedDownIndex < 12; feedDownIndex++){
    if(strcasecmp(a_particlename.c_str(), m_feeddown_particle_names[feedDownIndex].c_str()) == 0){
      foundFeedDownIndex = feedDownIndex;
      break;
    }
  }*/
  if(foundFeedDownIndex < 0){
    cout << "ERROR: HybridFeeddownMaker::loadDataSpectra : Can't find particle index : " << a_feeddownParticleId << endl;
    cout << "     You can try: ";
    for(int feedDownIndex = 0; feedDownIndex < 13; feedDownIndex++) cout << "  " << m_feeddown_particle_names[feedDownIndex];
    cout << endl;
    return;
  }



  if(!a_mTm0Space){
    #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
      cout << "    Converting pT spectra to mT-m0 spectra with mass: " << m_feeddown_particle_masses[foundFeedDownIndex] << endl;
    #endif
    m_measuredSpectra[foundFeedDownIndex][a_centIndex] = HistogramUtilities::convertSpectra_pT_to_mTm0(a_graph, m_feeddown_particle_masses[foundFeedDownIndex]);
  }else{
    m_measuredSpectra[foundFeedDownIndex][a_centIndex] = a_graph;
  }

  m_measuredSpectra[foundFeedDownIndex][a_centIndex]->SetName(Form("measuredMtM0Spectra_%s_Cent_%02d",m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));

  return;

}

//====================================================================================================================================================================================================================================================================================================================


void HybridFeeddownMaker::loadScaleAndStretchParameters(string a_filename){
  cout << "HybridFeeddownMaker::loadScaleAndStretchParameters   Loading " << a_filename;
  if(!m_urqmd_dNdy_scaled[0]){
    cout << "\n\n\n ERROR: Call void HybridFeeddownMaker::loadUrQMDSpectraFile(string a_filename, bool a_oldNames) before HybridFeeddownMaker::loadScaleAndStretchParameters(string a_filename) \n\n\n" << endl;
  }
  TFile* inFile = new TFile(a_filename.c_str(), "READ");
  cout << "    Ptr: " << inFile << endl;
  inFile->cd();
  //Make the Scale Stretch data histo
  m_matchingParametersHisto = (TH3D*) inFile->Get("scaleAndStretchParameterHisto");
  if(!m_matchingParametersHisto){
    cout << "ERROR: Unable to load the scaleAndStretchParameterHisto   : " << m_matchingParametersHisto << endl;
    exit(EXIT_FAILURE);
  }
  for(int feedDownIndex = 0; feedDownIndex < 20; feedDownIndex++){
    for(int centIndex = 0; centIndex < 16; centIndex++){
      for(int paramIndex = 0; paramIndex < 4; paramIndex++){
        m_matchingParameters[feedDownIndex][centIndex][paramIndex] = m_matchingParametersHisto->GetBinContent(feedDownIndex+1,centIndex+1,paramIndex+1);
      }
      if(m_urqmd_dNdy_scaled[feedDownIndex]){
        for(int rapBin = 1; rapBin <= m_urqmd_dNdy_scaled[feedDownIndex]->GetNbinsX(); rapBin++){
          //  "Matching Parameter 0 Scaling Safety"  if you ever change calcScalingFactor(), this must be edited
          double newContent = m_matchingParameters[feedDownIndex][centIndex][0] * m_urqmd_dNdy_scaled[feedDownIndex]->GetBinContent(rapBin,centIndex + 1);
          double newError   = m_matchingParameters[feedDownIndex][centIndex][0] * m_urqmd_dNdy_scaled[feedDownIndex]->GetBinError(rapBin,centIndex + 1);
          m_urqmd_dNdy_scaled[feedDownIndex]->SetBinContent(rapBin,centIndex + 1, newContent);
          m_urqmd_dNdy_scaled[feedDownIndex]->SetBinError(rapBin,centIndex + 1, newError);
        }
      }
    }
    HistogramUtilities::discardNanHistoData(m_urqmd_dNdy_scaled[feedDownIndex]);
  }


}

//====================================================================================================================================================================================================================================================================================================================


void HybridFeeddownMaker::findSimToDataScaling(int a_feeddownParticleId, int a_centIndex , vector< int > a_centralityBinsToUse, double a_lowRapidityOfData, double a_highRapidityOfData, vector<int> a_extraFeedDownIdsToInclude, bool a_isMostPeripheralCentralityMeasured){
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  fittingCanvas->cd();
  gPad->SetRightMargin(0.4);
  gPad->SetLogy();

  cout << "\n\n#### Starting findSimToDataScaling for Feeddown Id: " << a_feeddownParticleId << "  with mass = " << m_feeddown_particle_masses[a_feeddownParticleId] << endl;
  if(a_feeddownParticleId < 0 || a_feeddownParticleId > 12){
    cout << "ERROR: Bad input ID HybridFeeddownMaker::findSimToDataScaling(int a_feeddownParticleId)" << endl;
    return;
  }
  m_currentDataGraphToFit = m_measuredSpectra[a_feeddownParticleId][a_centIndex];
  m_currentDataHistoToFit = HistogramUtilities::makeHistogram_Asymm(m_currentDataGraphToFit);
  m_currentDataHistoToFit->SetLineColor(kYellow+2);
  m_currentSimHistoToFit_inDataBinning = (TH1D*) m_currentDataHistoToFit->Clone();
  m_currentSimHistoToFit_inDataBinning->SetName("m_currentSimHistoToFit_inDataBinning");
  m_currentSimHistoToFit_inDataBinning->SetLineColor(kGreen);
  m_currentDataHistoToFit_BinEdges = HistogramUtilities::getBinEdges(m_currentDataHistoToFit);

  cout << "Fitting Data Graph: " << m_measuredSpectra[a_feeddownParticleId][a_centIndex];
  cout << "  " << m_measuredSpectra[a_feeddownParticleId][a_centIndex]->GetName() << " with N: " << m_measuredSpectra[a_feeddownParticleId][a_centIndex]->GetN() << endl;

  m_currentFeedDownIndex = a_feeddownParticleId;







  //#####  Because these are spectra and we need to add up a few of them... things get tricky with error propagation
  m_currentSimHistoToFit = NULL;
  double dNdyNotMapped = 0.0;

  //Setup empty spectra histogram for urqmd
  TH2D* spectraHisto = m_urqmd_spectra[a_feeddownParticleId][0];
  m_currentSimHistoToFit = spectraHisto->ProjectionY(Form("simulationProjection_%s",m_feeddown_particle_names[a_feeddownParticleId].c_str()),1,1);
  HistogramUtilities::zeroBins(m_currentSimHistoToFit);
  double lowRapBin  = spectraHisto->GetXaxis()->FindBin(a_lowRapidityOfData);
  double highRapBin = spectraHisto->GetXaxis()->FindBin(a_highRapidityOfData);
  double lowRap  = spectraHisto->GetXaxis()->GetBinLowEdge(lowRapBin);
  double highRap = spectraHisto->GetXaxis()->GetBinLowEdge(highRapBin + 1);
  cout << "Setting the ranges: " << "  " << lowRapBin<< "  " << highRapBin<< "  " << lowRap<< "  " << highRap << endl;

  //Determine the number of events in each centrality bin
  double N_events = 0.0;
  vector<double> centEvents;
  for(unsigned int centVectIndex = 0; centVectIndex < a_centralityBinsToUse.size(); centVectIndex++){
    int centIndexSimulation = a_centralityBinsToUse[centVectIndex];
    cout << " Adding the Centrality: " << centIndexSimulation << " to the comparison" << endl;
    int midRapBin = m_urqmd_dNdy[a_feeddownParticleId]->GetXaxis()->FindBin(0.0);
    int binToCalcNEvt = midRapBin;
    if(m_urqmd_dNdy[a_feeddownParticleId]->GetBinContent(midRapBin,centIndexSimulation+1) <= 0.0){
      for(int binX = 1; binX <= m_urqmd_dNdy[a_feeddownParticleId]->GetNbinsX(); binX++){
        if(m_urqmd_dNdy[a_feeddownParticleId]->GetBinContent(binX,centIndexSimulation+1) > 0.0){
          binToCalcNEvt = binX;
          break;
        }
      }
    }
    if(m_urqmd_dNdy[a_feeddownParticleId]->GetBinContent(binToCalcNEvt,centIndexSimulation+1) <= 0.0){
      centEvents.push_back(0.0);
      continue; //there is no dNdy for any rapidity for this bin
    }
    double events_this_centrality = pow(m_urqmd_dNdy[a_feeddownParticleId]->GetBinContent(binToCalcNEvt,centIndexSimulation+1)/m_urqmd_dNdy[a_feeddownParticleId]->GetBinError(binToCalcNEvt,centIndexSimulation+1),2.0)
              *(1.0 + (1.0/(m_urqmd_dNdy[a_feeddownParticleId]->GetBinContent(binToCalcNEvt,centIndexSimulation+1)*m_urqmd_dNdy[a_feeddownParticleId]->GetXaxis()->GetBinWidth(binToCalcNEvt))) );
    centEvents.push_back(events_this_centrality);
    N_events += events_this_centrality;
    cout << "---- N_events: " << N_events << "   cent_events: " << events_this_centrality << endl;
    cout << " Mid Rap Bin: " << midRapBin << "  dN/dy: " << m_urqmd_dNdy[a_feeddownParticleId]->GetBinContent(binToCalcNEvt,centIndexSimulation+1) << " +/- "
         << m_urqmd_dNdy[a_feeddownParticleId]->GetBinError(binToCalcNEvt,centIndexSimulation+1)
         << "   rapidity bin width: " << m_urqmd_dNdy[a_feeddownParticleId]->GetXaxis()->GetBinWidth(binToCalcNEvt) << endl;

    for(int rapBin = 1; rapBin <= m_mappedUrQMD_dNdy[a_feeddownParticleId]->GetNbinsX(); rapBin++){
      // be careful with a_centIndex (measurment) and centIndexSimulation
      double newContent = m_mappedUrQMD_dNdy[a_feeddownParticleId]->GetBinContent(rapBin, a_centIndex + 1) + m_urqmd_dNdy[a_feeddownParticleId]->GetBinContent(rapBin,centIndexSimulation + 1);
      double newError   = PhysMath::addInQuadriture(m_mappedUrQMD_dNdy[a_feeddownParticleId]->GetBinContent(rapBin, a_centIndex + 1),m_urqmd_dNdy[a_feeddownParticleId]->GetBinError(rapBin,centIndexSimulation + 1));
      m_mappedUrQMD_dNdy[a_feeddownParticleId]->SetBinContent(rapBin,a_centIndex + 1, newContent);
      m_mappedUrQMD_dNdy[a_feeddownParticleId]->SetBinError(rapBin,a_centIndex + 1, newError);
    }

  } // end cent Index


  cout << " Beginning urqmd spectra graph ptr: " << m_urqmd_spectra[a_feeddownParticleId][0] << endl;
  for(int mTm0Bin = 1; mTm0Bin <= m_urqmd_spectra[a_feeddownParticleId][0]->GetNbinsY(); mTm0Bin++){
    double mTm0Width = m_urqmd_spectra[a_feeddownParticleId][0]->GetYaxis()->GetBinWidth(mTm0Bin);
    double mTm0      = m_urqmd_spectra[a_feeddownParticleId][0]->GetYaxis()->GetBinCenter(mTm0Bin);
    double mT        = mTm0 + m_feeddown_particle_masses[a_feeddownParticleId];

    double N_counts = 0.0;

    for(unsigned int centVectIndex = 0; centVectIndex < a_centralityBinsToUse.size(); centVectIndex++){
      int centIndexSimulation = a_centralityBinsToUse[centVectIndex];
      cout << " Adding the Centrality: " << centIndexSimulation << " to the comparison" << endl;
      spectraHisto = m_urqmd_spectra[a_feeddownParticleId][centIndexSimulation];

      for(int rapBin = lowRapBin; rapBin <= highRapBin; rapBin++){
        double content = spectraHisto->GetBinContent(rapBin,mTm0Bin);
        if(content == 0.0) continue;
        double error = spectraHisto->GetBinError(rapBin,mTm0Bin);
        double rapWidth = spectraHisto->GetXaxis()->GetBinWidth(rapBin);
        double factor = rapWidth * mTm0Width * 2.0 * 3.14159265359 * mT;
        N_counts += content* factor * centEvents[centVectIndex];

        //#### ADD THE EXTRA SPECTRA (ex. lambda + sigma0)
        for(int extraIndex = 0; extraIndex < a_extraFeedDownIdsToInclude.size(); extraIndex++){
          N_counts += m_urqmd_spectra[a_extraFeedDownIdsToInclude[extraIndex]][centIndexSimulation]->GetBinContent(rapBin,mTm0Bin) * factor * centEvents[centVectIndex];
        }

      }
    } // end cent Index

    dNdyNotMapped += N_counts;

    double newContent = N_counts/(N_events * 2.0 * 3.14159265359 * mT * (highRap - lowRap) * mTm0Width);
    double percentError = TMath::Sqrt(1.0/N_events + 1.0/N_counts);
    double newContentError = percentError * newContent;
    if(N_counts == 0.0 || !std::isfinite(newContent))      newContent = 0.0;
    if(N_events == 0.0 || !std::isfinite(newContentError)) newContentError = 0.0;
    cout << " Setting Spectra Projection  mTm0 = " << mTm0 << " val: " << newContent << " +/- " << newContentError << " N_counts: " << N_counts << "  N_events: " << N_events << "  mTm0Width: " << mTm0Width << " mT: " << mT << endl;
    m_currentSimHistoToFit->SetBinContent(mTm0Bin,newContent);
    m_currentSimHistoToFit->SetBinError(mTm0Bin,newContentError);

  }//end mTm0 Loop


  dNdyNotMapped = dNdyNotMapped / (N_events * (highRap - lowRap));








  //### Now we have made a simulation spectra with the same rapidity range
  m_currentSimGraphToFit = HistogramUtilities::convertNoAsymToAsym(HistogramUtilities::makeGraph(m_currentSimHistoToFit));
  m_currentSimGraphToFit->SetName(Form("spectra_urqmd_%s_Cent%02d_notMapped",m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
  m_currentSimGraphToFit->SetTitle(Form("Spectra from UrQMD for %s y=[%1.2f,%1.2f]",m_feeddown_particle_names[a_feeddownParticleId].c_str(),lowRap,highRap));
  m_unmappedUrQMDSpectra[a_feeddownParticleId][a_centIndex] = m_currentSimGraphToFit;
  m_currentSimHistoToFit_BinEdges = HistogramUtilities::getBinEdges(m_currentSimHistoToFit);

  double mTm0OfFirstData = m_currentDataGraphToFit->GetX()[0];
  double scale_seed      = m_currentDataGraphToFit->GetY()[0]/m_currentSimGraphToFit->Eval(mTm0OfFirstData);


  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout << "Minimizer Setup... Ptrs: " << m_currentSimGraphToFit << " " << m_currentDataGraphToFit << " " << m_currentSimHistoToFit << " " << m_currentDataHistoToFit << endl;
    cout << "  m_currentSimGraphToFit  N: " << m_currentSimGraphToFit->GetN() << endl;
    cout << "  m_currentDataGraphToFit N: " << m_currentDataGraphToFit->GetN() << endl;
    cout << "  m_currentSimHistoToFit  N: " << m_currentSimHistoToFit->GetNbinsX() << endl;
    cout << "  m_currentDataHistoToFit N: " << m_currentDataHistoToFit->GetNbinsX() << endl;
  #endif


  //### Create the minimizer
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer();
  ROOT::Math::Functor* m_chiSqrdFunctor = new ROOT::Math::Functor(this,&HybridFeeddownMaker::simToDataMatchingChiSqr_BinSplitting,2);
  minimizer->SetFunction(*m_chiSqrdFunctor);
  minimizer->SetVariable(0, "Scale", scale_seed,0.5);
  minimizer->SetVariableLimits(0,0.0, 200.0);
  minimizer->SetVariable(1, "Stretch", 1.0, 0.5);
  minimizer->SetVariableLimits(1, 0.001, 20.0);
  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    minimizer->SetPrintLevel(1);
  #else
    minimizer->SetPrintLevel(0);
  #endif
  minimizer->SetMaxIterations(300);
  minimizer->SetMaxFunctionCalls(10000);
  minimizer->SetValidError(true);
  minimizer->Minimize();
  int status = minimizer->Status();
  const double* finalParams = minimizer->X();
  int ndf = m_currentDataGraphToFit->GetN() - 2;
  double chiSqr = simToDataMatchingChiSqr_BinSplitting(finalParams);
  cout << "Simulation to Data Fit Status: " << status << "     X^2/ndf= " << chiSqr/ndf << endl;
  if(!(status == 0 || status == 4)){
    cout << "################################  ERROR:  STATUS OF FIT = " << status << "   IS BAD!!!!!  ################################" << endl;;
    return;
  }
  if(status == 4){
    cout << "WARNING:  Minimizer Fit Status = " << status << "  might be bad..." << endl;;
  }

  string combineNameList = " and";

  //need to push back scaling for the other centralities
  int highestCentIndex = -1;
  for(int centIndexIndex = 0; centIndexIndex < a_centralityBinsToUse.size(); centIndexIndex++){
    int centIndexSimulation = a_centralityBinsToUse[centIndexIndex];
    if(centIndexSimulation > highestCentIndex) highestCentIndex = centIndexSimulation;
    m_matchingParameters[a_feeddownParticleId][centIndexSimulation][0] = finalParams[0];
    m_matchingParameters[a_feeddownParticleId][centIndexSimulation][1] = finalParams[1];
    m_matchingParameters[a_feeddownParticleId][centIndexSimulation][2] = sqrt( minimizer->CovMatrix(0,0));
    m_matchingParameters[a_feeddownParticleId][centIndexSimulation][3] = sqrt( minimizer->CovMatrix(1,1));
    for(int extraIndexIndex = 0; extraIndexIndex < a_extraFeedDownIdsToInclude.size(); extraIndexIndex++){
      m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][0] = finalParams[0];
      m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][1] = finalParams[1];
      m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][2] = sqrt( minimizer->CovMatrix(0,0));
      m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][3] = sqrt( minimizer->CovMatrix(1,1));
      combineNameList += " ";
      combineNameList += m_feeddown_particle_names[a_extraFeedDownIdsToInclude[extraIndexIndex]].c_str();
    }
  }
  if(a_isMostPeripheralCentralityMeasured){  // if this is the highest centrality measured, push back the same scaling for future centralities
    for(int centIndexSimulation = highestCentIndex + 1; centIndexSimulation < 16; centIndexSimulation++){
      m_matchingParameters[a_feeddownParticleId][centIndexSimulation][0] = finalParams[0];
      m_matchingParameters[a_feeddownParticleId][centIndexSimulation][1] = finalParams[1];
      m_matchingParameters[a_feeddownParticleId][centIndexSimulation][2] = sqrt( minimizer->CovMatrix(0,0));
      m_matchingParameters[a_feeddownParticleId][centIndexSimulation][3] = sqrt( minimizer->CovMatrix(1,1));
      for(int extraIndexIndex = 0; extraIndexIndex < a_extraFeedDownIdsToInclude.size(); extraIndexIndex++){
        m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][0] = finalParams[0];
        m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][1] = finalParams[1];
        m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][2] = sqrt( minimizer->CovMatrix(0,0));
        m_matchingParameters[a_extraFeedDownIdsToInclude[extraIndexIndex]][centIndexSimulation][3] = sqrt( minimizer->CovMatrix(1,1));
      }
    }
  }


  //Zero the Bin Mapped Simulation
  HistogramUtilities::zeroBins(m_currentSimHistoToFit_inDataBinning);




  double dNdyMapped = 0.0;
  m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex] = (TGraphAsymmErrors*) m_currentSimGraphToFit->Clone(); // this already has the centralities and particles combined
  m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetName(Form("spectra_urqmd_%s_Cent%02d_mapped",m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
  m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetTitle(Form("Spectra from UrQMD for %s y=[%1.2f,%1.2f] (Scaled and Stretched to Data) %s",m_feeddown_particle_names[a_feeddownParticleId].c_str(),lowRap,highRap,combineNameList.c_str()));

  for(int pointIndex = 0; pointIndex < m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetN(); pointIndex++){
    double mTm0       = m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetX()[pointIndex];
    double mTm0_width = m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetEXlow()[pointIndex]
                        + m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetEXhigh()[pointIndex];
    double newMtM0        = finalParams[1] * m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetX()[pointIndex];
    double new_mTm0_width = finalParams[1] * mTm0_width;
    double scalingFactor  = calcScalingFactor(finalParams[0], finalParams[1], mTm0, mTm0_width, m_feeddown_particle_masses[m_currentFeedDownIndex]);
    double newValue       = scalingFactor * m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetY()[pointIndex];
    double newValueErr    = scalingFactor * m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetEYlow()[pointIndex];

    m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetPoint(pointIndex,newMtM0, newValue/new_mTm0_width);
    m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetPointError(pointIndex,m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetEXlow()[pointIndex] * finalParams[1],m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetEXhigh()[pointIndex] * finalParams[1],
                                                                           newValueErr/new_mTm0_width,newValueErr/new_mTm0_width);

    dNdyMapped += newValue * 2.0 * 3.14159265359 * (newMtM0 + m_feeddown_particle_masses[a_feeddownParticleId]);

    HistogramUtilities::addToHistogram(m_currentSimHistoToFit_inDataBinning,newValue,newValueErr,
                          (m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetX()[pointIndex]  - m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetEXlow()[pointIndex]),
                          (m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetX()[pointIndex]  + m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->GetEXhigh()[pointIndex]));
  }
  m_currentSimHistoToFit_inDataBinning->Scale(1.0,"width");


  //make a scaled version of the urqmd dN/dy
  double dNdyScale = dNdyMapped/dNdyNotMapped;
  cout << "UrQMD dN/dy: " << dNdyNotMapped << "  has been scaled to dN/dy: " << dNdyMapped << "   Ratio: " << dNdyScale << "     with params: " <<  finalParams[0] << " and " << finalParams[1] << endl;
  // Scale the dN/dy in the centIndex of measurment
  for(int rapBin = 1; rapBin <= m_mappedUrQMD_dNdy[a_feeddownParticleId]->GetNbinsX(); rapBin++){
    double newContent = dNdyScale * m_mappedUrQMD_dNdy[a_feeddownParticleId]->GetBinContent(rapBin,a_centIndex + 1);
    double newError   = dNdyScale * m_mappedUrQMD_dNdy[a_feeddownParticleId]->GetBinError(rapBin,a_centIndex + 1);
    m_mappedUrQMD_dNdy[a_feeddownParticleId]->SetBinContent(rapBin,a_centIndex + 1, newContent);
    m_mappedUrQMD_dNdy[a_feeddownParticleId]->SetBinError(rapBin,a_centIndex + 1, newError);
  }
  HistogramUtilities::discardNanHistoData(m_mappedUrQMD_dNdy[a_feeddownParticleId]);

  // Scale the dN/dy in the cent Index simulation
  for(int centIndexIndex = 0; centIndexIndex < a_centralityBinsToUse.size(); centIndexIndex++){
    int centIndexSimulation = a_centralityBinsToUse[centIndexIndex];
    for(int rapBin = 1; rapBin <= m_urqmd_dNdy_scaled[a_feeddownParticleId]->GetNbinsX(); rapBin++){
      double newContent = dNdyScale * m_urqmd_dNdy_scaled[a_feeddownParticleId]->GetBinContent(rapBin,centIndexSimulation + 1);
      double newError   = dNdyScale * m_urqmd_dNdy_scaled[a_feeddownParticleId]->GetBinError(rapBin,centIndexSimulation + 1);
      m_urqmd_dNdy_scaled[a_feeddownParticleId]->SetBinContent(rapBin,centIndexSimulation + 1, newContent);
      m_urqmd_dNdy_scaled[a_feeddownParticleId]->SetBinError(rapBin,centIndexSimulation + 1, newError);
    }
    //#### Do the same to the EXTRA SPECTRA (ex. lambda + sigma0)
    for(int extraIndex = 0; extraIndex < a_extraFeedDownIdsToInclude.size(); extraIndex++){
      for(int rapBin = 1; rapBin <= m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->GetNbinsX(); rapBin++){
        double newContent = dNdyScale * m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->GetBinContent(rapBin,centIndexSimulation + 1);
        double newError   = dNdyScale * m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->GetBinError(rapBin,centIndexSimulation + 1);
        m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->SetBinContent(rapBin,centIndexSimulation + 1, newContent);
        m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->SetBinError(rapBin,centIndexSimulation + 1, newError);
      }
    }
  }
  // copy the scaling parameters to the centralities not measured
  if(a_isMostPeripheralCentralityMeasured){
    for(int centIndexSimulation = highestCentIndex + 1; centIndexSimulation < 16; centIndexSimulation++){
      for(int rapBin = 1; rapBin <= m_urqmd_dNdy_scaled[a_feeddownParticleId]->GetNbinsX(); rapBin++){
        double newContent = dNdyScale * m_urqmd_dNdy_scaled[a_feeddownParticleId]->GetBinContent(rapBin,centIndexSimulation + 1);
        double newError   = dNdyScale * m_urqmd_dNdy_scaled[a_feeddownParticleId]->GetBinError(rapBin,centIndexSimulation + 1);
        m_urqmd_dNdy_scaled[a_feeddownParticleId]->SetBinContent(rapBin,centIndexSimulation + 1, newContent);
        m_urqmd_dNdy_scaled[a_feeddownParticleId]->SetBinError(rapBin,centIndexSimulation + 1, newError);
      }
      for(int extraIndex = 0; extraIndex < a_extraFeedDownIdsToInclude.size(); extraIndex++){
        for(int rapBin = 1; rapBin <= m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->GetNbinsX(); rapBin++){
          double newContent = dNdyScale * m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->GetBinContent(rapBin,centIndexSimulation + 1);
          double newError   = dNdyScale * m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->GetBinError(rapBin,centIndexSimulation + 1);
          m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->SetBinContent(rapBin,centIndexSimulation + 1, newContent);
          m_urqmd_dNdy_scaled[a_extraFeedDownIdsToInclude[extraIndex]]->SetBinError(rapBin,centIndexSimulation + 1, newError);
        }
      }
    }
  }
  HistogramUtilities::discardNanHistoData(m_urqmd_dNdy_scaled[a_feeddownParticleId]);


  m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetLineColor(kRed);
  m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerStyle(20);
  m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerSize(0.5);
  m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerColor(kRed);

  m_unmappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetLineColor(kViolet);
  m_unmappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerStyle(20);
  m_unmappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerSize(0.5);
  m_unmappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerColor(kViolet);

  m_measuredSpectra[a_feeddownParticleId][a_centIndex]->SetLineColor(kBlue);
  m_measuredSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerStyle(20);
  m_measuredSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerSize(0.5);
  m_measuredSpectra[a_feeddownParticleId][a_centIndex]->SetMarkerColor(kBlue);

  if(m_imageDir != ""){
    fittingCanvas->cd();
    double lowVertRange, highVertRange, lowHorRange, highHorRange;
    HistogramUtilities::verticalDrawRange(m_measuredSpectra[a_feeddownParticleId][a_centIndex], lowVertRange,highVertRange,-999,-999,true);
    HistogramUtilities::horizontalDrawRange(m_measuredSpectra[a_feeddownParticleId][a_centIndex], lowHorRange, highHorRange);
    lowHorRange = 0;
    highHorRange = 1.3*highHorRange;
    lowVertRange = 0.2*lowVertRange;
    highVertRange = 5.0*highVertRange;



    TH1I* baseHisto = new TH1I("base",Form("%s Cent %02d y_{data}=[%2.3f,%2.3f]; m_{T}-m_{0} [GeV/c^2]; Invariant Spectra",m_feeddown_particle_names[a_feeddownParticleId].c_str(), a_centIndex , a_lowRapidityOfData ,a_highRapidityOfData),10,lowHorRange,highHorRange);
    baseHisto->GetYaxis()->SetRangeUser(lowVertRange,highVertRange);
    baseHisto->Draw();

    TPaveText* fitParsTxt = new TPaveText(0.65,0.15,0.95,0.85,"NDC");
    fitParsTxt->AddText(Form("#Chi^{2}/ndf \t %e/%e = %e",chiSqr,ndf,chiSqr/ndf));
    fitParsTxt->AddText(Form("Y Scale \t %e#pm%e",finalParams[0],sqrt( minimizer->CovMatrix(0,0))));
    fitParsTxt->AddText(Form("X Stretch \t %e#pm%e",finalParams[1],sqrt( minimizer->CovMatrix(1,1))));
    fitParsTxt->SetFillStyle(0);
    fitParsTxt->SetBorderSize(5);
    fitParsTxt->Draw("same");

    m_measuredSpectra[a_feeddownParticleId][a_centIndex]->Draw("P");
    m_unmappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->Draw("P");
    m_mappedUrQMDSpectra[a_feeddownParticleId][a_centIndex]->Draw("P");
    m_currentSimHistoToFit_inDataBinning->SetLineWidth(2);
    m_currentSimHistoToFit_inDataBinning->Draw("SAME");
    m_currentDataHistoToFit->SetLineWidth(2);
    m_currentDataHistoToFit->Draw("SAME");
    fittingCanvas->SaveAs(Form("%s/%s_ScaleAndStretch_Cent%02d.png",
          m_imageDir.c_str(),m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
    delete baseHisto;
  }


  delete fittingCanvas;

}


//====================================================================================================================================================================================================================================================================================================================

#ifndef _BEFORE_ROOT_6p24_
void HybridFeeddownMaker::fitMeasurement(int a_feeddownParticleId,  int a_centIndex){


  if(!m_measuredSpectra[a_feeddownParticleId][a_centIndex]){
    cout << "ERROR: in fitMeasurement. Cant fit because feedDownIndex " << a_feeddownParticleId << " cent " << a_centIndex << " doesn't exist" << endl;
  }

  string bwName = Form("blastWave_%02d_%02d",a_feeddownParticleId,a_centIndex);

  TF1* blastFunct = PhysMath::blastWaveFunct(bwName.c_str(), 1.0, 0.5, 0.1, 0, 0, m_feeddown_particle_masses[a_feeddownParticleId]);
  blastFunct->SetParLimits(1,0.0,1.0);
  blastFunct->SetParLimits(2,0,1.0);
  blastFunct->FixParameter(3,0.5);
  blastFunct->FixParameter(4,1.0);
  blastFunct->FixParameter(5,m_partInfo->GetParticleMass(2));


  double X_ofPoint = m_measuredSpectra[a_feeddownParticleId][a_centIndex]->GetX()[0];
  double Y_ofPoint = m_measuredSpectra[a_feeddownParticleId][a_centIndex]->GetY()[0];
  double functEvalAtXPoint = blastFunct->Eval(X_ofPoint);
  double ampPredict = Y_ofPoint/functEvalAtXPoint;

  blastFunct->SetParameter(0,ampPredict);
  blastFunct->SetParLimits(0,0.01*ampPredict,10.0*ampPredict);

  int status = m_measuredSpectra[a_feeddownParticleId][a_centIndex]->Fit(bwName.c_str(),"EX0");
  m_measuredSpectra_Fits[a_feeddownParticleId][a_centIndex] = m_measuredSpectra[a_feeddownParticleId][a_centIndex]->GetFunction(bwName.c_str());
  if(status != 0){
    cout << "  There was a bad fit... figure out what went wrong..." << endl;
    m_measuredSpectra_Fits[30][30]->GetName();
  }

  return;


}
#endif


//====================================================================================================================================================================================================================================================================================================================

double HybridFeeddownMaker::simToDataMatchingChiSqr(const double* a_param){
  //no XErrorBar calc
  //cout << " Data Ptr: " << m_currentDataGraphToFit << " Sim Ptr: " << m_currentSimGraphToFit << endl;
  //cout << " Scale: " << a_param[0] << " Stretch: " << a_param[1] << endl;
  //m_currentDataGraphToFit->Print();
  //m_currentSimGraphToFit->Print();

  //#### Now, we need 1/N_evnt*dN/dmTdy to properly scale the UrQMD Simulation
  // as you increase the stretch, it needs to evaluate a new dmT and a new mT to recalculate yield

  double mass = m_feeddown_particle_masses[m_currentFeedDownIndex];

  double firstMtM0 = m_currentSimGraphToFit->GetX()[0];
  double lastMtM0 = m_currentSimGraphToFit->GetX()[m_currentSimGraphToFit->GetN()-1];

  double chiSqr = 0.0;
  for(int pointIndex = 0; pointIndex < m_currentDataGraphToFit->GetN(); pointIndex++){
    double mTm0 = m_currentDataGraphToFit->GetX()[pointIndex];
    double mT = mTm0 + mass;
    //double scalingFactor = a_param[0]*mT/(a_param[1] * (a_param[1] * mTm0 + mass)); // I know right
    double scalingFactor = a_param[1]*mT/(a_param[0] * (mTm0/a_param[1] + mass)); // This is an old version that doesn't use calcScalingFactor
    //double scalingFactor = calcScalingFactor(a_param[0], a_param[1], mTm0, mass);

    double specPoint = m_currentDataGraphToFit->GetY()[pointIndex] * scalingFactor;
    double specPointErr = m_currentDataGraphToFit->GetEYlow()[pointIndex] * scalingFactor;
    double yOther = 0.0;
    double yOtherErr = 0.0;
    double newDataMtM0 = mTm0/a_param[1];
    HistogramUtilities::evalWithError(m_currentSimGraphToFit, newDataMtM0,yOther,yOtherErr);
    if( newDataMtM0 < firstMtM0 ){
      yOther = m_currentSimGraphToFit->GetY()[0];
      yOtherErr = m_currentSimGraphToFit->GetEYlow()[0];
      cout << " Pushing before first point" << endl;
    }
    if( newDataMtM0 > lastMtM0 ){
      yOther = m_currentSimGraphToFit->GetY()[m_currentSimGraphToFit->GetN()-1];
      yOtherErr = m_currentSimGraphToFit->GetEYlow()[m_currentSimGraphToFit->GetN()-1];
      cout << " Pushing beyond last point" << endl;
    }
    if(yOther > specPoint) specPointErr = m_currentDataGraphToFit->GetEYhigh()[pointIndex]*scalingFactor;
    cout << "   yData = " << specPoint << " yDataErr = " << specPointErr << "  ySim: " << yOther << " ySimErr: " << yOtherErr << endl;
    chiSqr += pow(specPoint - yOther,2.0) / (specPointErr*specPointErr + yOtherErr*yOtherErr);
    //cout << " ChiSqrd: " << chiSqr << "  x: " << mTm0 << " y: " << specPoint << " +/- " << specPointErr << " y_urqmd: " << yOther << " +/- " << yOtherErr << endl;

  }
  cout << "## P1: " << a_param[0] << " P2: " << a_param[1] << " ChiSqr: " << chiSqr << endl;
  return chiSqr;
}

//====================================================================================================================================================================================================================================================================================================================

double HybridFeeddownMaker::simToDataMatchingChiSqr_BinSplitting(const double* a_param){
  //#### Now, we need 1/N_evnt*dN/dmTdy to properly scale the UrQMD Simulation
  // as you increase the stretch, it needs to evaluate a new dmT and a new mT to recalculate yield
  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout << "## P1: " << a_param[0] << " P2: " << a_param[1] << " Ptr: " << m_currentDataHistoToFit << endl;
  #endif

  int numDataPoints = m_currentDataGraphToFit->GetN();
  HistogramUtilities::zeroBins(m_currentSimHistoToFit_inDataBinning);
  double sim_content_in_data_binning[100];
  double sim_content_in_data_binning_sqrd_err[100];
  for(int iii = 0; iii < 100; iii++){
    sim_content_in_data_binning[iii] = 0.0;
    sim_content_in_data_binning_sqrd_err[iii] = 0.0;
  }

  double mass      = m_feeddown_particle_masses[m_currentFeedDownIndex];
  int numSimPoints = m_currentSimGraphToFit->GetN();
  //#############   Scale&Stretch The Simulation And Bin It into Data's Structure for Chi^2
  for(int sim_index = 0; sim_index < numSimPoints; sim_index++ ){
    double mTm0          = m_currentSimGraphToFit->GetX()[sim_index];
    double mTm0_width    = m_currentSimGraphToFit->GetEXlow()[sim_index] + m_currentSimGraphToFit->GetEXhigh()[sim_index];
    double scalingFactor = calcScalingFactor(a_param[0], a_param[1], mTm0, mTm0_width, mass);
    double specPoint         = scalingFactor * m_currentSimGraphToFit->GetY()[sim_index];
    double specPointLowErr   = scalingFactor * m_currentSimGraphToFit->GetEYlow()[sim_index];
    //double specPointHighErr  = scalingFactor * m_currentSimGraphToFit->GetEYhigh()[sim_index];
    double lowMtM0Edge   = a_param[1]*(mTm0 - m_currentSimGraphToFit->GetEXlow()[sim_index]);
    double highMtM0Edge  = a_param[1]*(mTm0 + m_currentSimGraphToFit->GetEXhigh()[sim_index]);

    HistogramUtilities::addToHistogram(m_currentSimHistoToFit_inDataBinning, specPoint,specPointLowErr,lowMtM0Edge,highMtM0Edge);



    /*int    firstDataIndex  = m_currentDataHistoToFit->FindBin(lowMtM0Edge) - 1; // for this method to work, the data can't have holes (zero content bins in the middle)
    int    lastDataIndex   = m_currentDataHistoToFit->FindBin(highMtM0Edge) - 1;
    if(firstDataIndex == lastDataIndex){
      sim_content_in_data_binning[firstDataIndex]          += specPoint;
      sim_content_in_data_binning_sqrd_err[firstDataIndex] += specPointLowErr*specPointLowErr;
    }else{
      double width       = highMtM0Edge - lowMtM0Edge;
      double partialFrac = (m_currentDataHistoToFit_BinEdges[firstDataIndex + 1] - lowMtM0Edge)/width;
      sim_content_in_data_binning[firstDataIndex]          +=  partialFrac*specPoint;
      sim_content_in_data_binning_sqrd_err[firstDataIndex] +=  partialFrac*specPointLowErr*specPointLowErr;
      for(int index = firstDataIndex+1; index < lastDataIndex; index++){
        partialFrac = (m_currentDataHistoToFit_BinEdges[index + 1] - m_currentDataHistoToFit_BinEdges[index])/width;
        sim_content_in_data_binning[index] += partialFrac*specPoint;
        sim_content_in_data_binning_sqrd_err[index] +=  partialFrac*specPointLowErr*specPointLowErr;
      }
      partialFrac = (highMtM0Edge - m_currentDataHistoToFit_BinEdges[lastDataIndex])/width;
      sim_content_in_data_binning[lastDataIndex] +=  partialFrac*specPoint;
      sim_content_in_data_binning_sqrd_err[lastDataIndex] =  partialFrac*specPointLowErr*specPointLowErr;
    }*/
  }
  m_currentSimHistoToFit_inDataBinning->Scale(1.0,"width");

  //###################   Calculate Chi Squared  ##########################
  double chiSqr = 0.0;
  for(int data_bin = 1; data_bin <= numDataPoints; data_bin++ ){
  //for(int data_index = 0; data_index < numDataPoints; data_index++ ){
    //double data          = m_currentDataGraphToFit->GetY()[data_index];
    //double sim           = sim_content_in_data_binning[data_index];
    //double sim_sqrd_err  = sim_content_in_data_binning_sqrd_err[data_index];
    //double data_err      = m_currentDataGraphToFit->GetEYlow()[data_index];
    //if(sim > data) data_err = m_currentDataGraphToFit->GetEYhigh()[data_index];

    double data          = m_currentDataHistoToFit->GetBinContent(data_bin);
    double sim           = m_currentSimHistoToFit_inDataBinning->GetBinContent(data_bin);
    double sim_err       = m_currentSimHistoToFit_inDataBinning->GetBinError(data_bin);
    double data_err      = m_currentDataHistoToFit->GetBinError(data_bin);
    if(std::isfinite(sim)){
      chiSqr += pow(data - sim,2.0) / (data_err*data_err + sim_err*sim_err);
    }
    #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
      cout << "  " << data_bin << " : " <<  data << " " << sim << "  " << (data_err*data_err + sim_err*sim_err) << endl;
    #endif
  }
  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout <<  "   ChiSqr: " << chiSqr << endl;
  #endif
  return chiSqr;

}



//====================================================================================================================================================================================================================================================================================================================

void HybridFeeddownMaker::generateResponseFunctions(int a_feeddownParticleId, int a_centIndex, int a_numTrials, double a_feeddownSpectraIndexToUse, double a_extraFactor, bool a_useFitFunct){
  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout << "HybridFeeddownMaker::generateResponseFunctions(int " <<  a_feeddownParticleId  << ", int " << a_centIndex << ", int " << a_numTrials << ", double " << a_feeddownSpectraIndexToUse << ", double " << a_extraFactor << ", bool " << a_useFitFunct << ")" << endl;
  #endif
  int feedDownSpectraIndexToUse = a_feeddownParticleId;
  if(a_feeddownSpectraIndexToUse >= 0) feedDownSpectraIndexToUse = a_feeddownSpectraIndexToUse;
  cout << " Using spectra index : " << feedDownSpectraIndexToUse << " for index :" << a_feeddownParticleId << " with an extra scale of " << a_extraFactor << endl;
  int partIndicies[4];
  int partCharges[4];
  partIndicies[0] = -1;
  partCharges[0]  =  0;
  partIndicies[1] = -1;
  partCharges[1]  =  0;
  partIndicies[2] = -1;
  partCharges[2]  =  0;
  partIndicies[3] = -1;
  partCharges[3]  =  0;

  double branchingRatio = 1.0;

  switch(a_feeddownParticleId){
    case 0 : //Lambda
      branchingRatio = 0.639;
      partIndicies[0] = 2; partCharges[0] = 1;
      partIndicies[1] = 0; partCharges[1] = -1;
      break;
    case 1 : //Antilambda
      branchingRatio = 0.639;
      partIndicies[0] = 2; partCharges[0] = -1;
      partIndicies[1] = 0; partCharges[1] = 1;
      break;
    case 2 : //Simga Plus
      branchingRatio = 0.5157;
      partIndicies[0] = 2; partCharges[0] = 1;
      break;
    case 3 : //AntiSigma Minus
      branchingRatio = 0.5157;
      partIndicies[0] = 2; partCharges[0] = -1;
      break;
    case 4 : //Xi Minus
      branchingRatio = 0.63827793;
      partIndicies[0] = 2; partCharges[0] = 1;
      partIndicies[1] = 0; partCharges[1] = -1;
      partIndicies[2] = 0; partCharges[2] = -1;
      break;
    case 5 : //AntiXi Plus
      branchingRatio = 0.63827793;
      partIndicies[0] = 2; partCharges[0] = -1;
      partIndicies[1] = 0; partCharges[1] = 1;
      partIndicies[2] = 0; partCharges[2] = 1;
      break;
    case 6 : //Xi Zero
      branchingRatio = 0.63596475;
      partIndicies[0] = 2; partCharges[0] = 1;
      partIndicies[1] = 0; partCharges[1] = -1;
      break;
    case 7 : //Omega Minus
      branchingRatio = 1.0; // branching is built-in
      partIndicies[0] = 2; partCharges[0] = 1;
      partIndicies[1] = 1; partCharges[1] = -1;
      partIndicies[2] = 0; partCharges[2] = -1;
      partIndicies[3] = 0; partCharges[3] = -1;
      break;
    case 8 : //AntiOmega Plus
      branchingRatio = 1.0; // branching is built-in
      partIndicies[0] = 2; partCharges[0] = -1;
      partIndicies[1] = 1; partCharges[1] = 1;
      partIndicies[2] = 0; partCharges[2] = 1;
      partIndicies[3] = 0; partCharges[3] = 1;
      break;
    case 9 : //Sigma Zero
      branchingRatio = 0.639;
      partIndicies[0] = 2; partCharges[0] = 1;
      partIndicies[1] = 0; partCharges[1] = -1;
      break;
    case 10 : //ANTISIGMA_0
      branchingRatio = 0.639; // for lambda proton decay
      partIndicies[0] = 2; partCharges[0] = -1;
      partIndicies[1] = 0; partCharges[1] = 1;
      break;
    case 11 : // 11: KAON_0_S
      branchingRatio = 0.6920;
      partIndicies[0] = 0; partCharges[0] = 1;
      partIndicies[1] = 0; partCharges[1] = -1;
      break;
    case 12 : //12: ETA
      branchingRatio = 1.0; // branching is built-in
      partIndicies[0] = 0; partCharges[0] = 1;
      partIndicies[1] = 0; partCharges[1] = -1;
      break;
  }


  TH2D* responseHistograms[4];
  TH1D* responseDnDy[4];
  TH1D* responseDCA[4];
  TH2D* responseDnDyDCA[4];
  for(int iii = 0; iii < 3; iii++){
    if(partIndicies[iii] < 0) continue;
    vector<double> rapBinEdges = HistogramUtilities::getBinEdges(m_response_structure_1D);
    vector<double> dcaBinEdges = HistogramUtilities::makeEvenBinEdges(100,0,15);
    #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
      cout << "Bin Edges for response structure: " << endl;
      for(unsigned int mmm = 0; mmm < rapBinEdges.size(); mmm++){
        cout << "  " << rapBinEdges[mmm];
      }
      cout << endl;
    #endif
    responseHistograms[iii] = (TH2D*) m_response_structure->Clone();
    responseHistograms[iii]->SetName(Form("daughterYield_%s_from_parent_%s_decay_Cent%02d", m_partInfo->GetParticleName(partIndicies[iii],partCharges[iii]).Data(), m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
    responseDnDy[iii] = (TH1D*) m_response_structure_1D->Clone();
    responseDnDy[iii]->SetName(Form("daughterDnDy_%s_from_parent_%s_decay_Cent%02d", m_partInfo->GetParticleName(partIndicies[iii],partCharges[iii]).Data(), m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
    responseDCA[iii] = new TH1D("",Form("DCA of %s from %s",m_partInfo->GetParticleSymbol(partIndicies[iii],partCharges[iii]).Data(),m_feeddown_particle_names[a_feeddownParticleId].c_str()),100,0.0,15.0);
    responseDCA[iii]->SetName(Form("daughterDCA_%s_from_parent_%s_decay_Cent%02d", m_partInfo->GetParticleName(partIndicies[iii],partCharges[iii]).Data(), m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
    responseDnDyDCA[iii] = new TH2D("","",((int) rapBinEdges.size())-1,&rapBinEdges[0],((int)dcaBinEdges.size())-1,&dcaBinEdges[0]);
    responseDnDyDCA[iii]->SetTitle(Form("dN/dy DCA of %s from %s;y;DCA [cm]",m_partInfo->GetParticleSymbol(partIndicies[iii],partCharges[iii]).Data(),m_feeddown_particle_names[a_feeddownParticleId].c_str()));
    responseDnDyDCA[iii]->SetName(Form("daughterDnDyDCA_%s_from_parent_%s_decay_Cent%02d", m_partInfo->GetParticleName(partIndicies[iii],partCharges[iii]).Data(), m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));

  }
  if(partIndicies[3] >= 0){ // for the two omega pion minus decays
    responseHistograms[3] = responseHistograms[2];
    responseDnDy[3]       = responseDnDy[2];
    responseDCA[3]        = responseDCA[2];
    responseDnDyDCA[3]    = responseDnDyDCA[2];
  }

  m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex] = (TH2D*) m_response_structure->Clone();
  m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex]->SetName(Form("resampledYield_from_parent_%s_decay_Cent%02d",  m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
  m_resampled_feeddown_spectra[a_feeddownParticleId][a_centIndex] = (TH2D*) m_response_structure->Clone();
  m_resampled_feeddown_spectra[a_feeddownParticleId][a_centIndex]->SetName(Form("resampledSpectra_from_parent_%s_decay_Cent%02d",  m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));


  TH2D* feeddownYieldHisto = (TH2D*) m_urqmd_spectra[feedDownSpectraIndexToUse][a_centIndex]->Clone();
  double totalProductionPerEvent = 0.0;
  //Convert the spectra to a simple yield
  for(int mTm0Bin = 1; mTm0Bin <= feeddownYieldHisto->GetNbinsY(); mTm0Bin++){
    double mTm0Width = feeddownYieldHisto->GetYaxis()->GetBinWidth(mTm0Bin);
    double mTm0      = feeddownYieldHisto->GetYaxis()->GetBinCenter(mTm0Bin);
    double mT        = mTm0 + m_feeddown_particle_masses[feedDownSpectraIndexToUse];
    for(int rapBin = 1; rapBin <= feeddownYieldHisto->GetNbinsX(); rapBin++){
      double content  = feeddownYieldHisto->GetBinContent(rapBin,mTm0Bin);
      if(content == 0.0) continue;
      double error    = feeddownYieldHisto->GetBinError(rapBin,mTm0Bin);
      double rapWidth = feeddownYieldHisto->GetXaxis()->GetBinWidth(rapBin);
      double factor   = rapWidth * mTm0Width * 2.0 * 3.14159265359 * mT;
      feeddownYieldHisto->SetBinContent(rapBin,mTm0Bin, content * factor);
      totalProductionPerEvent += content * factor;
      feeddownYieldHisto->SetBinError(rapBin,mTm0Bin,0);
    }
  }//end mTm0 Loop


  // make a scaled and stretched version of feeddownYield Histo
  string name  = Form("%s_%s_%d",feeddownYieldHisto->GetName(),"scaledAndStretched",a_feeddownSpectraIndexToUse);
  int numBinsX = feeddownYieldHisto->GetNbinsX();
  int numBinsY = feeddownYieldHisto->GetNbinsY();
  //double* edgesX = new double[numBinsX+1];
  //double* edgesY = new double[numBinsY+1];
  vector<double> edgesX(numBinsX+1,0);
  vector<double> edgesY(numBinsY+1,0);
  for(int binX = 1;binX <= numBinsX + 1; binX++){
    edgesX[binX-1] = feeddownYieldHisto->GetXaxis()->GetBinLowEdge(binX);
  }
  //#### DO THE STRETCHING
  for(int binY = 1;binY <= numBinsY + 1; binY++){
    edgesY[binY-1] = feeddownYieldHisto->GetYaxis()->GetBinLowEdge(binY) * m_matchingParameters[feedDownSpectraIndexToUse][a_centIndex][1];
  }

  //###### DO THE SCALING
  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout << "Using matching parameters : " << m_matchingParameters[feedDownSpectraIndexToUse][a_centIndex][0] << "   " << m_matchingParameters[feedDownSpectraIndexToUse][a_centIndex][1] << endl;
  #endif
  TH2D* feeddownYieldHisto_Scaled = new TH2D(name.c_str(),Form("%s_%s",feeddownYieldHisto->GetName()," Scaled and Stretched"),numBinsX,&edgesX[0],numBinsY,&edgesY[0]);
  for(int binX = 1; binX <= numBinsX ; binX++ ){
    for(int binY = 1; binY <= numBinsY ; binY++ ){
      double mTm0_original  = feeddownYieldHisto->GetYaxis()->GetBinCenter(binY);
      double mTm0_width     = feeddownYieldHisto->GetYaxis()->GetBinWidth(binY);
      double new_mTm0_width = m_matchingParameters[feedDownSpectraIndexToUse][a_centIndex][1] * mTm0_width;
      double scalingFactor  = calcScalingFactor(m_matchingParameters[feedDownSpectraIndexToUse][a_centIndex][0], m_matchingParameters[feedDownSpectraIndexToUse][a_centIndex][1], mTm0_original, mTm0_width, m_feeddown_particle_masses[a_feeddownParticleId]);
      feeddownYieldHisto_Scaled->SetBinContent(binX,binY,feeddownYieldHisto->GetBinContent(binX,binY) * scalingFactor / new_mTm0_width);
      feeddownYieldHisto_Scaled->SetBinError(binX,binY,feeddownYieldHisto->GetBinError(binX,binY) * scalingFactor / new_mTm0_width);
    }
  }
  if(a_feeddownParticleId == a_feeddownSpectraIndexToUse) HistogramUtilities::symmeterizeHistoX(feeddownYieldHisto_Scaled);


  //Each sample just gives a probability density, but doesn't change the overall scale. You need to scale each decay with a weight
  //This is different than
  //double decayWeight = HistogramUtilities::sumOfContent(feeddownYieldHisto_Scaled)/HistogramUtilities::sumOfContent(feeddownYieldHisto);
  //cout << " Decay Weight : " << decayWeight << " and scale parameter: " <<  m_matchingParameters[feedDownSpectraIndexToUse][a_centIndex][0] << endl;
  double decayWeight = 1.0;  // this decay weight is now just an overall scaling at the end

  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout << "Some Pointers: " << feeddownYieldHisto_Scaled << " " << m_feeddown_particle_masses[feedDownSpectraIndexToUse] << " " << m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex] << endl;
  #endif



  #ifdef _MULTI_THREADING_
    thread_local TRandom3* randObj = new TRandom3();
    randObj->SetSeed(time(NULL) + std::hash<std::thread::id>{}(std::this_thread::get_id()));
  #else
    TRandom3* randObj = new TRandom3();
    randObj->SetSeed(time(NULL));
  #endif


  ROOT::Math::Boost boostObj;
  TGenPhaseSpace generator;

  double parent_mass = m_feeddown_particle_masses[a_feeddownParticleId];
  double daughter_rapidity[4];
  double daughter_mTm0[4];
  double daughter_dca[4];
  int daughter_vertex_index[4];
  ROOT::Math::XYZVector* primaryVertex    = new ROOT::Math::XYZVector(0,0,0);
  ROOT::Math::XYZVector* secondaryVertex  = NULL;
  ROOT::Math::XYZVector* thirdVertex      = NULL;
  ROOT::Math::XYZVector* fourthVertex     = NULL;

  ROOT::Math::PxPyPzEVector* parent_vect = new ROOT::Math::PxPyPzEVector(0,0,0,0);
  ROOT::Math::PxPyPzEVector* daughter_vect[4];
  ROOT::Math::PxPyPzEVector* discard_vect; // for uncharged things like pi_0


  double parent_rapidity = 0;
  double parent_mTm0 = 0;
  double parent_pTot = 0;
  double parent_pX   = 0;
  double parent_pY   = 0;
  double parent_pT   = 0;
  double parent_pZ   = 0;
  double parent_energy = 0;

  //########################   TRIAL DECAY LOOP   ################################
  for(int trialIndex = 0; trialIndex < a_numTrials; trialIndex++){
    //if(randObj->Uniform() > branchingRatio) continue;
    secondaryVertex = NULL;
    thirdVertex     = NULL;
    fourthVertex    = NULL;
    discard_vect    = NULL;

    for(int iii = 0; iii < 4; iii++){
      daughter_rapidity[iii] = -999;
      daughter_mTm0[iii]     = -999;
      daughter_vect[iii]     = NULL;
    }

    if(a_useFitFunct){
      //sample the thing. Assume flat dNdy and go from y=[-1,1]
      #ifdef _BEFORE_ROOT_6p24_
        parent_mTm0 = m_measuredSpectra_Fits[feedDownSpectraIndexToUse][a_centIndex]->GetRandom();
      #else
        parent_mTm0 = m_measuredSpectra_Fits[feedDownSpectraIndexToUse][a_centIndex]->GetRandom(randObj);
      #endif
      parent_rapidity = randObj->Uniform(-1.0,1.0);
    }else{
      // randomly sample the scaled histogram to get parent kinematic info
      #ifdef _BEFORE_ROOT_6p24_
        feeddownYieldHisto_Scaled->GetRandom2(parent_rapidity,parent_mTm0);
      #else
        feeddownYieldHisto_Scaled->GetRandom2(parent_rapidity,parent_mTm0,randObj);
      #endif
    }

    parent_pT  = PhysMath::pT_from_mTm0(parent_mTm0,parent_mass);


    //parent_mTm0 = parent_mTm0 * m_matchingParameters[a_feeddownParticleId][a_centIndex][1]; //stretch the parent mtm0
    // IF YOU ARE ASSUMING THE SPECTRA ARE THE SAME, IT IS THE SAME IN pT SPACE, so map to eta pt and convert back to other particle
    if(a_feeddownParticleId >= 0 && feedDownSpectraIndexToUse != a_feeddownParticleId){
      parent_rapidity = PhysMath::rapFromEtaPt(m_feeddown_particle_masses[a_feeddownParticleId], PhysMath::eta_from_y_mTm0(parent_rapidity, parent_mTm0, m_feeddown_particle_masses[feedDownSpectraIndexToUse]), parent_pT);
      parent_mTm0     = PhysMath::mTm0(m_feeddown_particle_masses[a_feeddownParticleId], parent_pT);
    }
    m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex]->Fill(parent_rapidity,parent_mTm0);
    m_resampled_feeddown_spectra[a_feeddownParticleId][a_centIndex]->Fill(parent_rapidity,parent_mTm0);

    parent_pZ     = PhysMath::pZ_from_y_mTm0(parent_rapidity,parent_mTm0,parent_mass);
    parent_pTot   = sqrt(parent_pT*parent_pT + parent_pZ*parent_pZ);
    parent_energy = sqrt(parent_pTot*parent_pTot + parent_mass*parent_mass);

    randObj->Circle(parent_pX,parent_pY,parent_pT);
    parent_vect->SetPxPyPzE(parent_pX,parent_pY,parent_pZ,parent_energy);

    //double pionZero_rapidity = 0;
    //double pionZero_mTm0 = 0;
    daughter_vertex_index[0] = 2;
    daughter_vertex_index[1] = 2;
    daughter_vertex_index[2] = 2;
    daughter_vertex_index[3] = 2;

    switch(a_feeddownParticleId){
      case 0 : //Lambda
        PhysMath::lambdaDecay( &boostObj, randObj, primaryVertex, parent_vect,    daughter_vect[0], daughter_vect[1],secondaryVertex); // Lambda -> p + piMinus
        break;
      case 1 : //Antilambda
        PhysMath::lambdaDecay( &boostObj, randObj, primaryVertex, parent_vect,    daughter_vect[0], daughter_vect[1],secondaryVertex); // Lambda -> pbar + piPlus
        break;
      case 2 : //Simga Plus
        PhysMath::sigmaPlusDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1],secondaryVertex);
        break;
      case 3 : //AntiSigma Minus
        PhysMath::sigmaPlusDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1],secondaryVertex); // SigmaPlus -> p + pi0
        break;
      case 4 : //Xi Minus
        daughter_vertex_index[0] = 3;
        daughter_vertex_index[1] = 2;
        daughter_vertex_index[2] = 3;
        PhysMath::xiMinusDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1], daughter_vect[2],secondaryVertex,thirdVertex); // XiMinus -> p + piMinus + piMinus
        break;
      case 5 : //AntiXi Plus
        daughter_vertex_index[0] = 3;
        daughter_vertex_index[1] = 2;
        daughter_vertex_index[2] = 3;
        PhysMath::xiMinusDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1], daughter_vect[2],secondaryVertex,thirdVertex); // XiMinus -> p + piMinus + piMinus
        break;
      case 6 : //Xi Zero
        daughter_vertex_index[0] = 3;
        daughter_vertex_index[1] = 2;
        daughter_vertex_index[2] = 3;
        PhysMath::xiZeroDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], discard_vect, daughter_vect[2],secondaryVertex,thirdVertex); // XiZero -> p + pi0 + piMinus
        break;
      case 7 : //Omega Minus
        daughter_vertex_index[0] = 4;
        daughter_vertex_index[1] = 1;
        daughter_vertex_index[2] = 3;
        daughter_vertex_index[3] = 3;
        PhysMath::omegaMinusDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0],discard_vect, daughter_vect[1], daughter_vect[2], daughter_vect[3],secondaryVertex,thirdVertex,fourthVertex, daughter_vertex_index[0],daughter_vertex_index[1],daughter_vertex_index[2],daughter_vertex_index[3]); // XiZero -> p + pi0 + piMinus
        break;
      case 8 : //AntiOmega Plus
        daughter_vertex_index[0] = 4;
        daughter_vertex_index[1] = 1;
        daughter_vertex_index[2] = 3;
        daughter_vertex_index[3] = 3;
        PhysMath::omegaMinusDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0],discard_vect, daughter_vect[1], daughter_vect[2], daughter_vect[3],secondaryVertex,thirdVertex,fourthVertex, daughter_vertex_index[0],daughter_vertex_index[1],daughter_vertex_index[2],daughter_vertex_index[3]); // XiZero -> p + pi0 + piMinus
        break;
      case 9 : // Sigma Zero
        daughter_vertex_index[0] = 3;
        daughter_vertex_index[1] = 3;
        PhysMath::sigmaZeroDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1],secondaryVertex,thirdVertex); // Sigma Zero -> p + pi- + gamma
        break;
      case 10 : //ANTISIGMA_0
        daughter_vertex_index[0] = 3;
        daughter_vertex_index[1] = 3;
        PhysMath::sigmaZeroDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1],secondaryVertex,thirdVertex); // Anti-Sigma Zero -> pbar + pi+ + gamma
        break;
      case 11 : // 11: KAON_0_S
        PhysMath::kaonZeroShortDecay( &boostObj, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1],secondaryVertex);  // K0S -> pi+ + pi-
        break;
      case 12 : //12: ETA
        PhysMath::etaDecay( &generator, randObj, primaryVertex, parent_vect, daughter_vect[0], daughter_vect[1],secondaryVertex); // eta -> pi+ & pi-
        break;
    }

    /*
    switch(a_feeddownParticleId){
      case 0 : //Lambda
        PhysMath::lambdaDecay(     randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1]); // Lambda -> p + piMinus
        break;
      case 1 : //Antilambda
        PhysMath::lambdaDecay(     randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1]); // Lambda -> p + piMinus
        break;
      case 2 : //Simga Plus
        PhysMath::sigmaPlusDecay(  randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1]); // SigmaPlus -> p + pi0
        break;
      case 3 : //AntiSigma Minus
        PhysMath::sigmaPlusDecay(  randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1]); // SigmaPlus -> p + pi0
        break;
      case 4 : //Xi Minus
        PhysMath::xiMinusDecay(    randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1],daughter_rapidity[2],daughter_mTm0[2]); // XiMinus -> p + piMinus + piMinus
        break;
      case 5 : //AntiXi Plus
        PhysMath::xiMinusDecay(    randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1],daughter_rapidity[2],daughter_mTm0[2]); // XiMinus -> p + piMinus + piMinus
        break;
      case 6 : //Xi Zero
        PhysMath::xiZeroDecay(     randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], pionZero_rapidity   ,pionZero_mTm0   ,daughter_rapidity[2],daughter_mTm0[2]); // XiZero -> p + pi0 + piMinus
        break;
      case 7 : //Omega Minus
        PhysMath::omegaMinusDecay( randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1],daughter_rapidity[2],daughter_mTm0[2],daughter_rapidity[3],daughter_mTm0[3]);
        break;
      case 8 : //AntiOmega Plus
        PhysMath::omegaMinusDecay( randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1],daughter_rapidity[2],daughter_mTm0[2],daughter_rapidity[3],daughter_mTm0[3]);
        break;
      case 9 : // Sigma Zero
        PhysMath::sigmaZeroDecay(randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1]); // Sigma Zero -> p + pi- + gamma
        break;
      case 10 : //ANTISIGMA_0
        PhysMath::sigmaZeroDecay(randObj,  parent_rapidity, parent_mTm0, daughter_rapidity[0], daughter_mTm0[0], daughter_rapidity[1],daughter_mTm0[1]); // Anti-Sigma Zero -> pbar + pi+ + gamma
        break;
      case 11 : // 11: KAON_0_S
        PhysMath::kaonZeroShortDecay(randObj,  parent_rapidity, parent_mTm0,  daughter_rapidity[0],  daughter_mTm0[0],  daughter_rapidity[1],  daughter_mTm0[1]); // K0S -> pi+ + pi-
        break;
      case 12 : //12: ETA
        PhysMath::etaDecay(randObj,  parent_rapidity, parent_mTm0,  daughter_rapidity[0],  daughter_mTm0[0],  daughter_rapidity[1],  daughter_mTm0[1]); // eta -> pi+ & pi-
        break;

    }*/
    // omegas have other decay channels, so it may not give a result every time
    if(!secondaryVertex){
      if(thirdVertex)      delete thirdVertex;
      if(fourthVertex)     delete fourthVertex;
      if(daughter_vect[0]) delete daughter_vect[0];
      if(daughter_vect[1]) delete daughter_vect[1];
      if(daughter_vect[2]) delete daughter_vect[2];
      if(daughter_vect[3]) delete daughter_vect[3];
      if(discard_vect)     delete discard_vect;
      continue;
    }

    for(int daughter_index = 0; daughter_index < 4; daughter_index++){
      if(!daughter_vect[daughter_index]) continue;
      PhysMath::rapidity_mTm0(daughter_vect[daughter_index], daughter_rapidity[daughter_index], daughter_mTm0[daughter_index]);
      if(daughter_vertex_index[daughter_index] == 2)      daughter_dca[daughter_index] = dca_from_decay_helix(primaryVertex, secondaryVertex, daughter_vect[daughter_index], partCharges[daughter_index]);
      else if(daughter_vertex_index[daughter_index] == 3) daughter_dca[daughter_index] = dca_from_decay_helix(primaryVertex, thirdVertex,     daughter_vect[daughter_index], partCharges[daughter_index]);
      else if(daughter_vertex_index[daughter_index] == 4) daughter_dca[daughter_index] = dca_from_decay_helix(primaryVertex, fourthVertex,    daughter_vect[daughter_index], partCharges[daughter_index]);
    }

    if(m_doNTuple){
      #ifdef _MULTI_THREADING_
        {
          std::lock_guard<std::mutex> lock(mtx);
      #endif
        m_response_nTuple->Fill(a_feeddownParticleId,parent_rapidity,parent_mTm0,
                            partIndicies[0],daughter_rapidity[0],daughter_mTm0[0],
                            partIndicies[1],daughter_rapidity[1],daughter_mTm0[1],
                            partIndicies[2],daughter_rapidity[2],daughter_mTm0[2],
                            partIndicies[3],daughter_rapidity[3],daughter_mTm0[3]);
      #ifdef _MULTI_THREADING_
        }
      #endif
    }

    //#ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    //  cout << "Parent y: " << parent_rapidity << "  mTm0: " << parent_mTm0 << "\n"
    //  << "     D1 y: " << daughter_rapidity[0] << " mTm0: " << daughter_mTm0[0] << "   dca: " << daughter_dca[0] << "\n"
    //  << "     D2 y: " << daughter_rapidity[1] << " mTm0: " << daughter_mTm0[1] << "   dca: " << daughter_dca[1] << "\n"
    //  << "     D3 y: " << daughter_rapidity[2] << " mTm0: " << daughter_mTm0[2] << "   dca: " << daughter_dca[2] << "\n"
    //  << "     D4 y: " << daughter_rapidity[3] << " mTm0: " << daughter_mTm0[3] << "   dca: " << daughter_dca[3] << endl;
    //#endif

    for(int iii = 0; iii < 4; iii++){
      if(partIndicies[iii] < 0) continue;
      if(!m_use_DCA_cut || (m_use_DCA_cut && daughter_dca[iii] < m_DCA_Cut)){
        responseHistograms[iii]->Fill(daughter_rapidity[iii],daughter_mTm0[iii],decayWeight);
        responseDnDy[iii]->Fill(daughter_rapidity[iii],decayWeight);
      }
      responseDCA[iii]->Fill(daughter_dca[iii]);
      responseDnDyDCA[iii]->Fill(daughter_rapidity[iii],daughter_dca[iii],decayWeight);
    }

    if(secondaryVertex)  delete secondaryVertex;
    if(thirdVertex)      delete thirdVertex;
    if(fourthVertex)     delete fourthVertex;
    if(daughter_vect[0]) delete daughter_vect[0];
    if(daughter_vect[1]) delete daughter_vect[1];
    if(daughter_vect[2]) delete daughter_vect[2];
    if(daughter_vect[3]) delete daughter_vect[3];
    if(discard_vect)     delete discard_vect;
  }
  #ifdef _HYBRID_FEEDDOWN_MAKER_DEBUG_
    cout << " Finished Trials" << endl;
  #endif

  //Convert the resampled parent into an actual spectra
  TH2D* origHisto = m_resampled_feeddown_spectra[a_feeddownParticleId][a_centIndex];
  m_resampled_feeddown_spectra[a_feeddownParticleId][a_centIndex] = HistogramUtilities::convertFromYieldToSpectra2D(m_resampled_feeddown_spectra[a_feeddownParticleId][a_centIndex], a_numTrials, m_feeddown_particle_masses[a_feeddownParticleId]);
  delete origHisto;

  //Project the Yield to get a dN/dy
  m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex] = m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex]->ProjectionX(Form("resampled_feeddown_dNdy_%s_Cent%02d",
                        m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex),1,m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex]->GetNbinsY());
  m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]->Scale(1.0/((double) a_numTrials),"width");

  //Convert 2D Yield to differential
  HistogramUtilities::makeDifferential(m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex],1.0/((double) a_numTrials));

  // m_urqmd_dNdy_scaled is cent by dN/dy (is already differential (not in cent))
  m_urqmd_dNdy_scaled_byCent[a_feeddownParticleId][a_centIndex] = m_urqmd_dNdy_scaled[a_feeddownParticleId]->ProjectionX(Form("scaled_feeddown_dNdy_%s_Cent%02d",m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex),a_centIndex+1,a_centIndex+1);
  HistogramUtilities::discardNanHistoData(m_urqmd_dNdy_scaled_byCent[a_feeddownParticleId][a_centIndex]);

  // m_mappedUrQMD_dNdy is cent by dN/dy (is already differential (not in cent))
  // This is someone unnecessary, because this is the centrality index space of the measurment, not the simulation and current work
  m_mappedUrQMD_dNdyByCent[a_feeddownParticleId][a_centIndex] = m_mappedUrQMD_dNdy[a_feeddownParticleId]->ProjectionX(Form("scaled_feeddown_dNdy_%s_CentData%02d",m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex),a_centIndex+1,a_centIndex+1);
  HistogramUtilities::discardNanHistoData(m_mappedUrQMD_dNdyByCent[a_feeddownParticleId][a_centIndex]);
  HistogramUtilities::discardNanHistoData(m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]);

  //double dNdyFromTrials  = HistogramUtilities::sumOfContent(m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex],-1,true);
  //double dNdyFromScaling = HistogramUtilities::sumOfContent(m_urqmd_dNdy_scaled_byCent[a_feeddownParticleId][a_centIndex],-1,true);

  double dNdyFromTrials  = HistogramUtilities::integral(m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]);
  double dNdyFromScaling = HistogramUtilities::integral(m_urqmd_dNdy_scaled_byCent[a_feeddownParticleId][a_centIndex],
           HistogramUtilities::getLowEdge(m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]),
           HistogramUtilities::getHighEdge(m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]));
  double dNdyScaleFactor     = dNdyFromScaling/dNdyFromTrials;
  double daughterScaleFactor = branchingRatio*dNdyFromScaling/dNdyFromTrials;
  m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]->Scale(dNdyScaleFactor);
  m_resampled_feeddown_yields[a_feeddownParticleId][a_centIndex]->Scale(dNdyScaleFactor);
  m_resampled_feeddown_spectra[a_feeddownParticleId][a_centIndex]->Scale(dNdyScaleFactor);

  cout << "  dN/dy readjustment   Trials: " << dNdyFromTrials << "   UrQMDScaled: " << dNdyFromScaling << "    ratio: " << dNdyScaleFactor;
  double dNdyFromTrialsAfter  = HistogramUtilities::integral(m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]);
  cout << "  Trials After: " << dNdyFromTrialsAfter << endl;


  for(int iii = 0; iii < 4; iii++){
    if(partIndicies[iii] < 0) continue;
    int pmIndex = 1;
    if(partCharges[iii] > 0) pmIndex = 0;
    HistogramUtilities::makeSqrtNErrors(responseDnDy[iii]);
    //responseDnDy[iii]->Scale(((double) totalProductionPerEvent)/((double) a_numTrials)*a_extraFactor); // gets back to estimate of production per event
    m_produced_feeddown_DCA[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]  = responseDCA[iii];
    m_produced_feeddown_dNdy[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex] = responseDnDy[iii];
    m_produced_feeddown_dNdy[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetName(Form("daughterDnDy_%s_from_parent_%s_decay_Cent%02d", m_partInfo->GetParticleName(partIndicies[iii],partCharges[iii]).Data(), m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
    m_produced_feeddown_dNdy_DCA[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex] = responseDnDyDCA[iii];
    for(int rapBin = 1; rapBin <= responseDnDy[iii]->GetNbinsX(); rapBin++){
      double content = responseDnDy[iii]->GetBinContent(rapBin);
      if(content == 0.0) continue;
      double error = sqrt(content);
      double rapWidth = responseDnDy[iii]->GetBinWidth(rapBin);
      //double factor = ((double) totalProductionPerEvent)*a_extraFactor / ((double) a_numTrials);
      double factor = a_extraFactor / ((double) a_numTrials);
      m_produced_feeddown_dNdy[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinContent(rapBin, factor * content / rapWidth);
      m_produced_feeddown_dNdy[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinError(rapBin, factor * error / rapWidth);
    }
    m_produced_feeddown_dNdy[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->Scale(daughterScaleFactor);

  }

  for(int iii = 0; iii < 4; iii++){
    if(partIndicies[iii] < 0) continue;
    int pmIndex = 1;
    if(partCharges[iii] > 0) pmIndex = 0;
    HistogramUtilities::makeSqrtNErrors2D(responseHistograms[iii]);
    //responseHistograms[iii]->Scale(((double) totalProductionPerEvent) / ((double) a_numTrials) * a_extraFactor); // gets back to estimate of production per event
    m_produced_feeddown_yields[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]  = responseHistograms[iii];
    m_produced_feeddown_spectra[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex] = (TH2D*) responseHistograms[iii]->Clone();
    m_produced_feeddown_spectra[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetName(Form("daughterSpectra_%s_from_parent_%s_decay_Cent%02d", m_partInfo->GetParticleName(partIndicies[iii],partCharges[iii]).Data(), m_feeddown_particle_names[a_feeddownParticleId].c_str(),a_centIndex));
    for(int mTm0Bin = 1; mTm0Bin <= responseHistograms[iii]->GetNbinsY(); mTm0Bin++){
      double mTm0Width = responseHistograms[iii]->GetYaxis()->GetBinWidth(mTm0Bin);
      double mTm0 = responseHistograms[iii]->GetYaxis()->GetBinCenter(mTm0Bin);
      double mT = mTm0 + m_partInfo->GetParticleMass(partIndicies[iii]);
      for(int rapBin = 1; rapBin <= responseHistograms[iii]->GetNbinsX(); rapBin++){
        double content = responseHistograms[iii]->GetBinContent(rapBin,mTm0Bin);
        if(content == 0.0) continue;
        double error = responseHistograms[iii]->GetBinError(rapBin,mTm0Bin);
        double trials = (double) a_numTrials;
        double percentError = sqrt( (error*error) / (content*content) + 1.0/trials );
        double rapWidth = responseHistograms[iii]->GetXaxis()->GetBinWidth(rapBin);
        //double factor = ((double) totalProductionPerEvent)*a_extraFactor / (trials * rapWidth * mTm0Width * 2.0 * 3.14159265359 * mT);
        double factor = a_extraFactor / (trials * rapWidth * mTm0Width * 2.0 * 3.14159265359 * mT);
        //m_produced_feeddown_yields[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinContent(rapBin, mTm0Bin, ((double) totalProductionPerEvent) * a_extraFactor * content / trials);
        //m_produced_feeddown_yields[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinError(rapBin, mTm0Bin, percentError * ((double) totalProductionPerEvent) * a_extraFactor * content / trials);
        m_produced_feeddown_yields[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinContent(rapBin, mTm0Bin, a_extraFactor * content / (trials * rapWidth * mTm0Width));
        m_produced_feeddown_yields[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinError(rapBin, mTm0Bin, percentError * a_extraFactor * content / (trials * rapWidth * mTm0Width));
        m_produced_feeddown_spectra[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinContent(rapBin, mTm0Bin, content * factor);
        m_produced_feeddown_spectra[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->SetBinError(rapBin, mTm0Bin, percentError * content * factor);
      }
    }//end mTm0 Loop
    m_produced_feeddown_yields[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->Scale(daughterScaleFactor);
    m_produced_feeddown_spectra[partIndicies[iii]][a_feeddownParticleId][a_centIndex][pmIndex]->Scale(daughterScaleFactor);
  }

  HistogramUtilities::discardNanHistoData(m_mappedUrQMD_dNdyByCent[a_feeddownParticleId][a_centIndex]);
  HistogramUtilities::discardNanHistoData(m_urqmd_dNdy_scaled_byCent[a_feeddownParticleId][a_centIndex]);
  HistogramUtilities::discardNanHistoData(m_resampled_feeddown_dNdy[a_feeddownParticleId][a_centIndex]);




}


void HybridFeeddownMaker::combineFeedDownSpectra(){

  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    int charge = 1;
    if(pmIndex == 1) charge = -1;
    for(int partIndex = 0; partIndex < 9; partIndex++){
      cout << " combining spectra for pmIndex: " << pmIndex << " partIndex: " << partIndex << endl;
      for(int centIndex = 0; centIndex < 16; centIndex++){
        //find the first yield that exists for each particle
        bool feeddownSpectraFound = false;
        TH2D* ptrToFirstYield   = NULL;
        TH1D* ptrToFirstDnDy    = NULL;
        TH2D* ptrToFirstDnDyDCA = NULL;
        int feedDownIndex_ofFirst = -1;
        for(int feedDownIndex = 0; feedDownIndex < 13; feedDownIndex++){
          if(m_produced_feeddown_yields[partIndex][feedDownIndex][centIndex][pmIndex]){
            ptrToFirstYield   = m_produced_feeddown_yields[partIndex][feedDownIndex][centIndex][pmIndex];
            ptrToFirstDnDy    = m_produced_feeddown_dNdy[partIndex][feedDownIndex][centIndex][pmIndex];
            ptrToFirstDnDyDCA = m_produced_feeddown_dNdy_DCA[partIndex][feedDownIndex][centIndex][pmIndex];
            feeddownSpectraFound = true;
            feedDownIndex_ofFirst = feedDownIndex;
            break;
          }
        }
        if(!feeddownSpectraFound){
          cout << "  WARNING: No Spectra Found For Particle " << partIndex << " in centrality " << centIndex << endl;
          continue;
        }

        TH2D* completeYields = (TH2D*) ptrToFirstYield->Clone();
        completeYields->SetName(Form("totalFeedDownYields_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));
        completeYields->SetTitle(Form("Combined Feed-Down Yields (dN/dydmT) for %s Centrality %02d; y ; m_{T} - m_{0} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));

        TH1D* completeDnDy = (TH1D*) ptrToFirstDnDy->Clone();
        completeDnDy->SetName(Form("totalFeedDown_dNdy_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));
        completeDnDy->SetTitle(Form("Combined Feed-Down dN/dy for %s Centrality %02d; y; dN/dy",m_partInfo->GetParticleSymbol(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));

        TH2D* completeDnDyDCA = (TH2D*) ptrToFirstDnDyDCA->Clone();
        completeDnDyDCA->SetName(Form("totalFeedDown_dNdy_DCA_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));
        completeDnDyDCA->SetTitle(Form("Combined Feed-Down dN/dy vs DCA for %s Centrality %02d; y; DCA (cm); dN/dy",m_partInfo->GetParticleSymbol(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));



        for(int feedDownIndex = feedDownIndex_ofFirst+1; feedDownIndex < 13; feedDownIndex++){
          if(m_produced_feeddown_dNdy[partIndex][feedDownIndex][centIndex][pmIndex]){
            completeDnDy->Add(m_produced_feeddown_dNdy[partIndex][feedDownIndex][centIndex][pmIndex]);
          }
          if(m_produced_feeddown_dNdy_DCA[partIndex][feedDownIndex][centIndex][pmIndex]){
            completeDnDyDCA->Add(m_produced_feeddown_dNdy_DCA[partIndex][feedDownIndex][centIndex][pmIndex]);
          }
          if(m_produced_feeddown_yields[partIndex][feedDownIndex][centIndex][pmIndex]){
            completeYields->Add(m_produced_feeddown_yields[partIndex][feedDownIndex][centIndex][pmIndex]);
          }
        }

        TH2D* completeSpectra = (TH2D*) completeYields->Clone();
        completeSpectra->SetName(Form("totalFeedDownSpectra_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));
        completeSpectra->SetTitle(Form("Combined Feed-Down Spectra for %s Centrality %02d; y ; m_{T} - m_{0} [GeV/c^{2}]",m_partInfo->GetParticleSymbol(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex));


        //Convert From Yield to Spectra
        for(int mTm0Bin = 1; mTm0Bin <= completeYields->GetNbinsY(); mTm0Bin++){
          double mTm0Width = completeYields->GetYaxis()->GetBinWidth(mTm0Bin);
          double mT = completeYields->GetYaxis()->GetBinCenter(mTm0Bin) + m_partInfo->GetParticleMass(partIndex);
          for(int rapBin = 1; rapBin <= completeYields->GetNbinsX(); rapBin++){
            double content = completeYields->GetBinContent(rapBin,mTm0Bin);
            if(content == 0.0) continue;
            double error = completeYields->GetBinError(rapBin,mTm0Bin);
            double rapWidth = completeYields->GetXaxis()->GetBinWidth(rapBin);
            //double factor = rapWidth * mTm0Width * 2.0 * 3.14159265359 * mT;
            double factor =  2.0 * 3.14159265359 * mT;

            completeSpectra->SetBinContent(rapBin, mTm0Bin, content / factor);
            completeSpectra->SetBinError(rapBin, mTm0Bin, error / factor);

            //completeYields->SetBinContent(rapBin, mTm0Bin, content / (rapWidth * mTm0Width));
            //completeYields->SetBinError(rapBin, mTm0Bin, error / (rapWidth * mTm0Width));
          }
        }//end mTm0 Loop
        m_produced_total_feeddown_yields[partIndex][centIndex][pmIndex]   = completeYields;
        m_produced_total_feeddown_spectra[partIndex][centIndex][pmIndex]  = completeSpectra;
        m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]     = completeDnDy;
        m_produced_total_feeddown_dNdy_DCA[partIndex][centIndex][pmIndex] = completeDnDyDCA;
      }
    }
  }




}





void HybridFeeddownMaker::fitTotalFeedDownSpectra(int a_partIndex, int a_charge, int a_centIndex,string a_imageDir){

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  fittingCanvas->cd();
  gPad->SetRightMargin(0.4);
  gPad->SetLogy();

  int pmIndex = 0;
  if(a_charge == -1) pmIndex = 1;
  TH2D* spectra2D = m_produced_total_feeddown_spectra[a_partIndex][a_centIndex][pmIndex];

  if(!spectra2D){
    cout << "  ERROR: Trying to fit the Total Feeddown Spectra, but it doesn't exist for part : " << a_partIndex << " charge: " << a_charge << "  cent: " << a_centIndex << endl;
    return;
  }


  for(int rapBin = 1; rapBin <= spectra2D->GetNbinsX(); rapBin++){

    TH1D* projHisto = spectra2D->ProjectionY(Form("fdProj_%02d_%02d_%02d_%02d",a_partIndex,pmIndex,a_centIndex,rapBin),rapBin,rapBin);
    projHisto->SetTitle(Form("%s Feed Down Spectra Cent %02d yIndex %02d;m_{T}-m_{0}; Feed Down Spectra",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),a_centIndex,rapBin-1));
    string bw_name = Form("blast_%02d_%02d_%02d_%02d",a_partIndex,pmIndex,a_centIndex,rapBin);
    TF1* blastWave =  PhysMath::blastWaveFunct(bw_name.c_str(), 1, 0.3, 0.3, 0.5, 1.0, 0.93827);
    blastWave->FixParameter(4,1.0);
    blastWave->SetParLimits(1,0,1.0);
    blastWave->SetParLimits(2,0,1.0);
    blastWave->SetParLimits(3,0,1.0);
    blastWave->SetParLimits(5,0,10.0);

    double initAmp = projHisto->GetBinContent(3)/blastWave->Eval(projHisto->GetBinCenter(3));

    blastWave->SetParameter(0,initAmp);
    blastWave->SetParLimits(0,0.001*initAmp,1000.0*initAmp);


    TFitResultPtr resultPtr = projHisto->Fit(bw_name.c_str(),"S");
    TF1* fitResult = (TF1*) projHisto->GetFunction(bw_name.c_str())->Clone();
    fitResult->SetName(Form("totalFeedDownSpectraFit_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,rapBin-1));

    m_produced_total_feeddown_spectraFits[a_partIndex][a_centIndex][pmIndex][rapBin - 1] = fitResult; // sum of all the spectra from all parents  // part cent +/- rapIndex
    m_produced_total_feeddown_spectraFitsCovMat[a_partIndex][a_centIndex][pmIndex][rapBin - 1] = (TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();


    if(a_imageDir != ""){
      fittingCanvas->cd();
      projHisto->Draw("E");
      fitResult->Draw("Same");
      fittingCanvas->SaveAs(Form("%s/feeddownSpectraFit_%s_Cent%02d_yIndex%02d.png",
            a_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,rapBin-1));
    }


    delete projHisto;
  }//end rap loop








}




  // FEED DOWN ID NAMES: 0:Lambda 1:Antilambda 2:SigmaPlus 3:AntiSigmaMinus 4:XiMinus 5:AntiXiPlus 6:XiZero 7:OmegaMinus 8:AntiOmegaPlus  9: SigmaZero 10:ANTISIGMA_0  11: KAON_0_S 12: ETA



double HybridFeeddownMaker::getSigmaPlusOverLambdaRatio(vector<int> a_centIndicies){
  double sumOfSigmaPlus = 0.0;
  double sumOfLambda = 0.0;
  int midRapBin = m_urqmd_dNdy[2]->GetXaxis()->FindBin(0.0);
  for(int centIndexIndex = 0; centIndexIndex < a_centIndicies.size();centIndexIndex++){
    //for(int rapBin = 1; rapBin <= m_urqmd_dNdy[2]->GetNbinsX(); rapBin++){
    //}
    sumOfSigmaPlus += m_urqmd_dNdy[2]->GetBinContent(midRapBin,a_centIndicies[centIndexIndex]+1);
    sumOfLambda    += m_urqmd_dNdy[0]->GetBinContent(midRapBin,a_centIndicies[centIndexIndex]+1);

  }
  return sumOfSigmaPlus/sumOfLambda;
}
double HybridFeeddownMaker::getAntiSigmaMinusOverAntiLambdaRatio(vector<int> a_centIndicies){
  double sumOfAntiSigmaMinus = 0.0;
  double sumOfAntiLambda = 0.0;
  int midRapBin = m_urqmd_dNdy[3]->GetXaxis()->FindBin(0.0);
  for(int centIndexIndex = 0; centIndexIndex < a_centIndicies.size();centIndexIndex++){
    //for(int rapBin = 1; rapBin <= m_urqmd_dNdy[2]->GetNbinsX(); rapBin++){
    //}
    sumOfAntiSigmaMinus += m_urqmd_dNdy[3]->GetBinContent(midRapBin,a_centIndicies[centIndexIndex]+1);
    sumOfAntiLambda    += m_urqmd_dNdy[1]->GetBinContent(midRapBin,a_centIndicies[centIndexIndex]+1);

  }
  return sumOfAntiSigmaMinus/sumOfAntiLambda;
}

double HybridFeeddownMaker::getSigmaPlusOverLambdaRatio(int a_centIndex){
  int midRapBin = m_urqmd_dNdy[2]->GetXaxis()->FindBin(0.0);
  double ratio = m_urqmd_dNdy[2]->GetBinContent(midRapBin,a_centIndex+1)/m_urqmd_dNdy[0]->GetBinContent(midRapBin,a_centIndex+1);
  if(ratio <= 0.0){
    cout << "ERROR: HybridFeeddownMaker::getSigmaPlusOverLambdaRatio(" << a_centIndex << ")" << " found a ratio of " << ratio << endl;
    cout << "Necessary Pointers: " << m_urqmd_dNdy[2] << "   and    " << m_urqmd_dNdy[0] << endl;
    HistogramUtilities::printHistoInfo(m_urqmd_dNdy[0]);
    HistogramUtilities::printHistoInfo(m_urqmd_dNdy[2]);
    exit(EXIT_FAILURE);
  }
  return ratio;
}
double HybridFeeddownMaker::getAntiSigmaMinusOverAntiLambdaRatio(int a_centIndex){
  int midRapBin = m_urqmd_dNdy[3]->GetXaxis()->FindBin(0.0);
  double ratio = m_urqmd_dNdy[3]->GetBinContent(midRapBin,a_centIndex+1)/m_urqmd_dNdy[1]->GetBinContent(midRapBin,a_centIndex+1);
  if(ratio <= 0.0){
    cout << "ERROR: HybridFeeddownMaker::getAntiSigmaMinusOverAntiLambdaRatio(" << a_centIndex << ")" << " found a ratio of " << ratio << endl;
    exit(EXIT_FAILURE);
  }
  return ratio;
}

double HybridFeeddownMaker::dca_from_decay_helix(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, int a_daughter_charge){
  if(a_daughter_charge == 0) return PhysMath::dca_from_decay(a_primaryVertex,a_secondaryVertex,a_daughter_vect);
  TVector3 momentum_inGeV(a_daughter_vect->Px(),a_daughter_vect->Py(),a_daughter_vect->Pz());
  TVector3 decay_vertex(0.01*a_primaryVertex->x(),0.01*a_primaryVertex->y(),0.01*a_primaryVertex->z());
  TVector3 event_vertex(0.01*a_secondaryVertex->x(),0.01*a_secondaryVertex->y(),0.01*a_secondaryVertex->z());
  Helix helix(momentum_inGeV,decay_vertex,-0.5*1E-11,a_daughter_charge);
  return 100.0*fabs(helix.geometricSignedDistance(event_vertex));
}


void HybridFeeddownMaker::write(string a_outFileName){

  if(m_symmeterizeOutputs) symmeterizeTotalHistograms();

  //Make the Scale Stretch data histo
  m_matchingParametersHisto = new TH3D("scaleAndStretchParameterHisto","Scale and Stretch Parameters; FeedDownIndex; CentIndex; ParamIndex",20,-0.5,19.5,16,-0.5,15.5,4,-0.5,3.5);
  for(int feedDownIndex = 0; feedDownIndex < 20; feedDownIndex++){
    for(int centIndex = 0; centIndex < 16; centIndex++){
      for(int paramIndex = 0; paramIndex < 4; paramIndex++){
        m_matchingParametersHisto->SetBinContent(feedDownIndex+1,centIndex+1,paramIndex+1,m_matchingParameters[feedDownIndex][centIndex][paramIndex]);
      }
    }
  }

  TFile* outFile = new TFile(a_outFileName.c_str(),"RECREATE");
  outFile->cd();
  if(m_doNTuple) m_response_nTuple->Write();
  HistogramUtilities::ConditionalWrite(m_matchingParametersHisto);
  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    for(int partIndex = 0; partIndex < 9; partIndex++){
      for(int centIndex = 0; centIndex < 16; centIndex++){
        HistogramUtilities::ConditionalWrite(m_produced_total_feeddown_yields[partIndex][centIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_produced_total_feeddown_spectra[partIndex][centIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_produced_total_feeddown_dNdy_DCA[partIndex][centIndex][pmIndex]);

        if(m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]){
          double centralDnDy = 0.0;
          double square_error = 0.0;
          double lowBin    = m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->FindBin(-0.5);
          double lowRange  = m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->GetBinLowEdge(lowBin);
          double highBin   = m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->FindBin(0.5);
          double highRange = m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->GetBinLowEdge(highBin+1);
          for(int binX = lowBin; binX <= highBin; binX++){
            centralDnDy += m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->GetBinContent(binX) * m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->GetBinWidth(binX);
            square_error += pow(m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->GetBinError(binX) * m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]->GetBinWidth(binX),2.0);
          }
          double error = sqrt(square_error)/ (highRange - lowRange);
          centralDnDy = centralDnDy / (highRange - lowRange);
          cout << "  dN/dy y=[-0.5,0.5] for " << m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1) << " Cent : " << centIndex << "   dNdy: " << centralDnDy << "  +/- " << error << endl;
        }

      }
    }
  }


  outFile->mkdir("TotalFeedDownSpectraFits");
  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    for(int partIndex = 0; partIndex < 9; partIndex++){
      outFile->mkdir(Form("TotalFeedDownSpectraFits/%s",m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1).Data()));
      outFile->cd(Form("TotalFeedDownSpectraFits/%s",m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1).Data()));
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int yIndex = 0; yIndex < 50; yIndex++){
          HistogramUtilities::ConditionalWrite(m_produced_total_feeddown_spectraFits[partIndex][centIndex][pmIndex][yIndex]);
          if(m_produced_total_feeddown_spectraFitsCovMat[partIndex][centIndex][pmIndex][yIndex]){
            m_produced_total_feeddown_spectraFitsCovMat[partIndex][centIndex][pmIndex][yIndex]->Write(
                  Form("totalFeedDownSpectraFitCovMat_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex,yIndex));
          }
        }//rap loop
      }//cent loop
    }// particle loop
  }// plus minus loop



  outFile->cd();
  for(int feedDownIndex = 0; feedDownIndex < 13; feedDownIndex++){
    TGraphErrors* scaleGraph = new TGraphErrors();
    scaleGraph->SetName(Form("scaleParam_byCent_%s",m_feeddown_particle_names[feedDownIndex].c_str()));
    TGraphErrors* stretchGraph = new TGraphErrors();
    stretchGraph->SetName(Form("stretchParam_byCent_%s",m_feeddown_particle_names[feedDownIndex].c_str()));
    HistogramUtilities::ConditionalWrite(m_urqmd_dNdy[feedDownIndex]);
    HistogramUtilities::ConditionalWrite(m_urqmd_dNdy_scaled[feedDownIndex]);
    HistogramUtilities::ConditionalWrite(m_mappedUrQMD_dNdy[feedDownIndex]);

    for(int centIndex = 0; centIndex < 16; centIndex++){
      HistogramUtilities::ConditionalWrite(m_urqmd_spectra[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_urqmd_yields[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_measuredSpectra[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_mappedUrQMDSpectra[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_unmappedUrQMDSpectra[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_resampled_feeddown_yields[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_resampled_feeddown_dNdy[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_resampled_feeddown_spectra[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_measuredSpectra_Fits[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_urqmd_dNdy_scaled_byCent[feedDownIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_mappedUrQMD_dNdyByCent[feedDownIndex][centIndex]);
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        for(int partIndex = 0; partIndex < 9; partIndex++){
          HistogramUtilities::ConditionalWrite(m_produced_feeddown_dNdy[partIndex][feedDownIndex][centIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_produced_feeddown_dNdy_DCA[partIndex][feedDownIndex][centIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_produced_feeddown_DCA[partIndex][feedDownIndex][centIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_produced_feeddown_yields[partIndex][feedDownIndex][centIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_produced_feeddown_spectra[partIndex][feedDownIndex][centIndex][pmIndex]);
          //if(feedDownIndex == 0) HistogramUtilities::ConditionalWrite(m_produced_total_feeddown_spectra[partIndex][centIndex][pmIndex]);
        }
      }
      int scalePointIndex = scaleGraph->GetN();
      if(m_matchingParameters[feedDownIndex][centIndex][0] < -900) continue;
      scaleGraph->SetPoint(scalePointIndex, centIndex, m_matchingParameters[feedDownIndex][centIndex][0]);
      scaleGraph->SetPointError(scalePointIndex, 0.0, m_matchingParameters[feedDownIndex][centIndex][2]);
      stretchGraph->SetPoint(scalePointIndex,centIndex,m_matchingParameters[feedDownIndex][centIndex][1]);
      stretchGraph->SetPointError(scalePointIndex,0.0,m_matchingParameters[feedDownIndex][centIndex][3]);
    }
    scaleGraph->Write();
    stretchGraph->Write();
  }
}




void HybridFeeddownMaker::symmeterizeTotalHistograms(){

  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    for(int partIndex = 0; partIndex < 9; partIndex++){
      for(int centIndex = 0; centIndex < 16; centIndex++){
        HistogramUtilities::symmeterizeHistoX(m_produced_total_feeddown_yields[partIndex][centIndex][pmIndex]);
        HistogramUtilities::symmeterizeHistoX(m_produced_total_feeddown_spectra[partIndex][centIndex][pmIndex]);
        HistogramUtilities::symmeterizeHistoX(m_produced_total_feeddown_dNdy[partIndex][centIndex][pmIndex]);
        HistogramUtilities::symmeterizeHistoX(m_produced_total_feeddown_dNdy_DCA[partIndex][centIndex][pmIndex]);
      }
    }
  }
}





#endif
