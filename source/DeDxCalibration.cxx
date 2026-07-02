#ifndef _DEDX_CALIBRATION_CXX_
#define _DEDX_CALIBRATION_CXX_

#include "../headers/DeDxCalibration.h"


DeDxCalibration::DeDxCalibration(){
  m_initialCalibCurve      = NULL;
  m_initialCalibCurve_betaGamma = NULL;
  m_vZ_calib_struct        = NULL;
  m_eta_calib_struct       = NULL;
  m_betaGamma_calib_struct = NULL;
  m_betaGamma_calib_struct_fromTOF = NULL;
  m_currentHistoToFit      = NULL;
  m_imageDir = "";
  m_partInfo = new ParticleInfo();
  m_numBins_cdf_rebin = 15;
  m_testing = false;

  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int vZ_index = 0; vZ_index < g_DeDxCalib_MaxNumVzBins; vZ_index++){
      for(int eta_index = 0; eta_index < g_DeDxCalib_MaxNumEtaBins; eta_index++){
        m_mean_calibrationGraphs[partIndex][vZ_index][eta_index]     = NULL;
        m_meanPlusSigma_calibrationGraphs[partIndex][vZ_index][eta_index]     = NULL;
        m_sigma_calibrationGraphs[partIndex][vZ_index][eta_index]    = NULL;
        m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]  = NULL;
        m_meanPlusSigma_calibrationFunctions[partIndex][vZ_index][eta_index] = NULL;
        m_sigma_calibrationFunctions[partIndex][vZ_index][eta_index] = NULL;
        for(int bg_index = 0; bg_index < g_DeDxCalib_MaxNumBetaGammaBins; bg_index++){
          m_lnDeDx_calibrationHistograms[partIndex][vZ_index][eta_index][bg_index] = NULL;
          m_lnDeDx_calibrationHistograms_fromTOF[partIndex][vZ_index][eta_index][bg_index] = NULL;
        }
      }
    }
  }

  m_maxBetaGammaForTPCFitting[0] = 4.3;
  m_maxBetaGammaForTPCFitting[1] = 4.0;
  m_maxBetaGammaForTPCFitting[2] = 4.1;
  m_maxBetaGammaForTPCFitting[3] = 50.0;
  m_maxBetaGammaForTPCFitting[4] = 50.0;
  m_maxBetaGammaForTPCFitting[5] = 50.0;
  m_maxBetaGammaForTPCFitting[6] = 50.0;
  m_maxBetaGammaForTPCFitting[7] = 50.0;
  m_maxBetaGammaForTPCFitting[8] = 50.0;


}

//##############################################################################################################################################


void DeDxCalibration::loadPicoBinnerOutput(string a_filename, int a_partIndex, int a_betaGamma_rebin){
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  //std::unique_ptr<TFile> inFile( TFile::Open(a_filename.c_str(),"READ") );
  if(!inFile || inFile->IsZombie()){
    cout << "ERROR:  DeDxCalibration::loadPicoBinnerOutput.. couldn't find file " << a_filename << endl;
    return;
  }

  //TF1* readCalibFunct = (TF1*) inFile->Get("calibration_InvBetaGammaSqr_dEdx");
  TH1D* readParams = (TH1D*) inFile->Get("calibration_InvBetaGammaSqr_dEdx_params");
  if(readParams){ //if every particle was done at once in a single root file

    double scale = readParams->GetBinContent(6);
    if(!m_initialCalibCurve){
      m_initialCalibCurve = PhysMath::getBichselFunctionInvBetaGammaSqr("initalCalib_invBetaGammaSqr");
      m_initialCalibCurve->SetParameter(0,readParams->GetBinContent(1)/scale);// A = P1*P2
      m_initialCalibCurve->SetParameter(1,readParams->GetBinContent(2)/scale);// B = P1
      m_initialCalibCurve->SetParameter(2,readParams->GetBinContent(3)/scale);// C = P3
      m_initialCalibCurve->SetParameter(3,readParams->GetBinContent(4)/scale);// D = P4/2
      m_initialCalibCurve->SetParameter(4,readParams->GetBinContent(5)/scale);// E = P5/2
      m_initialCalibCurve->SetLineColor(kBlue);
      //m_initialCalibCurve->SetLineStyle(10);
      m_initialCalibCurve->SetRange(0.0001,30);
      //m_initialCalibCurve = readCalibFunct;
      m_initialCalibCurve_betaGamma = PhysMath::getBichselFunctionBetaGamma("initalCalib_betaGamma");
      m_initialCalibCurve_betaGamma->SetParameter(0,m_initialCalibCurve->GetParameter(0));// A = P1*P2
      m_initialCalibCurve_betaGamma->SetParameter(1,m_initialCalibCurve->GetParameter(1));// B = P1
      m_initialCalibCurve_betaGamma->SetParameter(2,m_initialCalibCurve->GetParameter(2));// C = P3
      m_initialCalibCurve_betaGamma->SetParameter(3,m_initialCalibCurve->GetParameter(3));// D = P4/2
      m_initialCalibCurve_betaGamma->SetParameter(4,m_initialCalibCurve->GetParameter(4));// E = P5/2
      m_initialCalibCurve_betaGamma->SetLineColor(kBlue);
      //m_initialCalibCurve_betaGamma->SetLineStyle(10);
      m_initialCalibCurve_betaGamma->SetRange(0.0001,30);
    }
  }else{
    cout << "Error: DeDxCalibration::loadPicoBinnerOutput  Can't Load calibration_InvBetaGammaSqr_dEdx"  << endl;
    return;
  }



  if(!m_vZ_calib_struct){
  	TH1I* vz_struct = (TH1I*) inFile->Get("Vz_Calibration/vZ_calib_struct");
    if(!vz_struct){
      cout << "ERROR: No Vz Calibration Structure can be read... you must have an old file or didn't define the VZ Calib in makefile_toggles.h " << endl;
      return;
    }
    //m_vZ_calib_struct = (TH1I*) vz_struct->Clone();
    m_vZ_calib_struct = vz_struct;
  }

  for(int vZ_index = 0; vZ_index < m_vZ_calib_struct->GetNbinsX(); vZ_index++){
    TH3I* histo3D_Int = (TH3I*) inFile->Get(Form("Vz_Calibration/calib_%s_VzIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index));
    if(!histo3D_Int){
      cout << "ERROR: Can't load the Vz 3D histogram for the calibration..." << endl;
      continue;
    }
    TH3I* histo3D_fromTOF_Int = (TH3I*) inFile->Get(Form("Vz_Calibration/calib_%s_VzIndex%02d_fromTOF",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index));


    TH3D* histo3D     = HistogramUtilities::convertInt_to_DoubleHisto3D(histo3D_Int);
    TH3D* histo3D_Rebin = histo3D;
    if(a_betaGamma_rebin > 1){
      histo3D_Rebin = HistogramUtilities::reBinVar3D(histo3D,1,a_betaGamma_rebin,1);
      delete histo3D;
    }
    TH3D* histo3D_fromTOF = HistogramUtilities::convertInt_to_DoubleHisto3D(histo3D_fromTOF_Int);


    // The read histogram has axis: eta betagamma ln(dedx)

    //Construct Axis Structures So eta_index and betagamma_index make sense
    if(!m_eta_calib_struct){
      vector< double > eta_edges = HistogramUtilities::getBinEdges_X(histo3D_Rebin); 
      m_eta_calib_struct = new TH1I("eta_calib_struct","Bin Structure for #eta; #eta; N/A",((int) eta_edges.size())-1,&eta_edges[0]);

      vector< double > bg_edges = HistogramUtilities::getBinEdges_Y(histo3D_Rebin); 
      m_betaGamma_calib_struct = new TH1I("betaGamma_calib_struct","Bin Structure for #beta#gamma; #beta#gamma; N/A",((int) bg_edges.size())-1,&bg_edges[0]);
      
      vector< double > bg_edges_fromTOF = HistogramUtilities::getBinEdges_Y(histo3D_fromTOF); 
      m_betaGamma_calib_struct_fromTOF = new TH1I("betaGamma_calib_struct_fromTOF","Bin Structure for #beta#gamma; #beta#gamma; N/A",((int) bg_edges_fromTOF.size())-1,&bg_edges_fromTOF[0]);
      cout << "setting up from TOF with " << m_betaGamma_calib_struct_fromTOF->GetNbinsX() << "bins" << endl;
    }

    //############   Make the 1D Projections for Future Fitting  ###################
    for(int eta_index = 0; eta_index < histo3D_Rebin->GetNbinsX(); eta_index++){
      for(int bg_index = 0; bg_index < histo3D_Rebin->GetNbinsY(); bg_index++){
        TH1D* histo = histo3D_Rebin->ProjectionZ(Form("proj_%s_vZ%02d_eta%02d_bg%03d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index,bg_index),eta_index+1,eta_index+1,bg_index+1,bg_index+1);
        double counts = HistogramUtilities::sumOfContent(histo);
        if(counts < 30){
          delete histo;
          continue;
        }
        HistogramUtilities::makeSqrtNErrors(histo);
        m_lnDeDx_calibrationHistograms[a_partIndex][vZ_index][eta_index][bg_index] = HistogramUtilities::rebinEqualByCDF(histo, m_numBins_cdf_rebin, false, true); // isDifferrential doLogCDF
        m_lnDeDx_calibrationHistograms[a_partIndex][vZ_index][eta_index][bg_index]->Scale(0.01,"width");
        m_lnDeDx_calibrationHistograms[a_partIndex][vZ_index][eta_index][bg_index]->SetTitle(
          Form("Calibration %s V_{Z}=[%4.2f,%4.2f] #eta=[%1.3f,%1.3f] #beta#gamma=[%1.3f,%1.3f]; ln(dE/dx) - c_{0}; dN/d(ln(dE/dx))",
          	  m_partInfo->GetParticleName(a_partIndex).Data(),
              m_vZ_calib_struct->GetBinLowEdge(vZ_index+1),m_vZ_calib_struct->GetBinLowEdge(vZ_index+2),
              m_eta_calib_struct->GetBinLowEdge(eta_index+1),m_eta_calib_struct->GetBinLowEdge(eta_index+2),
              m_betaGamma_calib_struct->GetBinLowEdge(bg_index+1),m_betaGamma_calib_struct->GetBinLowEdge(bg_index+2)));
        delete histo; // the rebin is saved
      }
      for(int bg_index = 0;bg_index < histo3D_fromTOF->GetNbinsY(); bg_index++){
        TH1D* histo = histo3D_fromTOF->ProjectionZ(Form("proj_%s_vZ%02d_eta%02d_bg%03d_fromTOF",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index,bg_index),eta_index+1,eta_index+1,bg_index+1,bg_index+1);
        double counts = HistogramUtilities::sumOfContent(histo);
        if(counts < 30){
          delete histo;
          continue;
        }
        HistogramUtilities::makeSqrtNErrors(histo);
        m_lnDeDx_calibrationHistograms_fromTOF[a_partIndex][vZ_index][eta_index][bg_index] = HistogramUtilities::rebinEqualByCDF(histo, m_numBins_cdf_rebin, false, true); // isDifferrential doLogCDF
        m_lnDeDx_calibrationHistograms_fromTOF[a_partIndex][vZ_index][eta_index][bg_index]->Scale(0.01,"width"); // so N doesn't hit that stupid 1E7 limit in root
        m_lnDeDx_calibrationHistograms_fromTOF[a_partIndex][vZ_index][eta_index][bg_index]->SetTitle(
          Form("Calibration %s V_{Z}=[%4.2f,%4.2f] #eta=[%1.3f,%1.3f] #beta#gamma=[%1.3f,%1.3f] from TOF; ln(dE/dx) - c_{0}; dN/d(ln(dE/dx))",
              m_partInfo->GetParticleName(a_partIndex).Data(),
              m_vZ_calib_struct->GetBinLowEdge(vZ_index+1),m_vZ_calib_struct->GetBinLowEdge(vZ_index+2),
              m_eta_calib_struct->GetBinLowEdge(eta_index+1),m_eta_calib_struct->GetBinLowEdge(eta_index+2),
              m_betaGamma_calib_struct_fromTOF->GetBinLowEdge(bg_index+1),m_betaGamma_calib_struct_fromTOF->GetBinLowEdge(bg_index+2)));
        delete histo; // the rebin is saved
      }      
    }// end eta loop

    delete histo3D_Rebin;
    delete histo3D_fromTOF;
  }// end Vz Index Loop

  //inFile->Close();

  
}

//##############################################################################################################################################

void DeDxCalibration::loadSaveFile(string a_filename){
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  m_initialCalibCurve           = (TF1*)  inFile->Get("calibration_InvBetaGammaSqr_dEdx");
  m_initialCalibCurve_betaGamma = (TF1*)  inFile->Get("initalCalib_betaGamma");
  m_vZ_calib_struct   = (TH1I*) inFile->Get("vZ_calib_struct");
  m_eta_calib_struct  = (TH1I*) inFile->Get("eta_calib_struct");
  m_betaGamma_calib_struct         = (TH1I*) inFile->Get("betaGamma_calib_struct");
  m_betaGamma_calib_struct_fromTOF = (TH1I*) inFile->Get("betaGamma_calib_struct_fromTOF");

  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int vZ_index = 0; vZ_index < g_DeDxCalib_MaxNumVzBins; vZ_index++){
      for(int eta_index = 0; eta_index < g_DeDxCalib_MaxNumEtaBins; eta_index++){
        m_mean_calibrationGraphs[partIndex][vZ_index][eta_index] = (TGraphAsymmErrors*) inFile->Get(
          Form("calibration_graphs/mean_calibration_graph_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(partIndex).Data(),vZ_index,eta_index));
        m_sigma_calibrationGraphs[partIndex][vZ_index][eta_index] = (TGraphAsymmErrors*) inFile->Get(
          Form("calibration_graphs/sigma_calibration_graph_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(partIndex).Data(),vZ_index,eta_index));
      }
    }
  }



}



void DeDxCalibration::setInitialCalibrationCurve(TF1* a_funct){
  if(!a_funct){
    cout << "ERROR: DeDxCalibration::setInitialCalibrationCurve, Bad function pointer supplied:" << a_funct << endl;
    return;
  }
  m_initialCalibCurve = (TF1*) a_funct->Clone();
  m_initialCalibCurve->SetName("initialCalibrationCurve");
}

//##############################################################################################################################################

void DeDxCalibration::setVzStructureHistogram(TH1I* a_histo){
  if(!a_histo){
    cout << "ERROR: DeDxCalibration::setInitialCalibrationCurve, Bad histogram pointer supplied:" << a_histo << endl;
    return;
  }
  m_vZ_calib_struct = a_histo;
}

//##############################################################################################################################################

void DeDxCalibration::setCdfRebinNumBins(int a_numBins){
  m_numBins_cdf_rebin = a_numBins;
}


//##############################################################################################################################################


void DeDxCalibration::setImageDir(string a_directory){
  m_imageDir = a_directory;
  gSystem->mkdir(m_imageDir.c_str());
}

//##############################################################################################################################################

void DeDxCalibration::setTestingToggle(bool a_toggle){
  m_testing = a_toggle;
}



void DeDxCalibration::extractMeansAndWidths(int a_partIndex, bool a_saveImages){
  TVirtualFitter::SetMaxIterations( 1000000 );
  if(m_imageDir.length() > 1){
    gSystem->mkdir(Form("%s/%s", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));
    gSystem->mkdir(Form("%s/%s/01_GaussianFits", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));
  }

  #ifdef _DEDX_CALIBRATION_DEBUG_
    cout << " eta structure check : " << m_eta_calib_struct << endl;
    cout << "       name: " << m_eta_calib_struct->GetName() << endl;
  #endif

  for(int vZ_index = 0; vZ_index < m_vZ_calib_struct->GetNbinsX(); vZ_index++){
    if(m_testing && vZ_index != 3) continue;
    for(int eta_index = 0; eta_index < m_eta_calib_struct->GetNbinsX(); eta_index++){
      //Make the TGraphs for the mean and sigma for this Vz_Eta combination
      if(m_testing && eta_index != 4) continue;
      if(!m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]){
        m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]  = new TGraphAsymmErrors();
        m_meanPlusSigma_calibrationGraphs[a_partIndex][vZ_index][eta_index] = new TGraphAsymmErrors();
        m_sigma_calibrationGraphs[a_partIndex][vZ_index][eta_index] = new TGraphAsymmErrors();
        m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetName(
          Form("mean_calibration_graph_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index));
        m_meanPlusSigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetName(
          Form("meanPlusSigma_calibration_graph_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index));
        m_sigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetName(
          Form("sigma_calibration_graph_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index));
        m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetTitle(
          Form("%s ln(dE/dx) #mu Calibration V_{Z}=[%4.2f,%4.2f] #eta=[%1.3f,%1.3f]",m_partInfo->GetParticleSymbol(a_partIndex).Data(),
            m_vZ_calib_struct->GetBinLowEdge(vZ_index+1),m_vZ_calib_struct->GetBinLowEdge(vZ_index+2),
            m_eta_calib_struct->GetBinLowEdge(eta_index+1),m_eta_calib_struct->GetBinLowEdge(eta_index+2)));
        m_meanPlusSigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetTitle(
          Form("%s ln(dE/dx) #mu + #sigma Calibration V_{Z}=[%4.2f,%4.2f] #eta=[%1.3f,%1.3f]",m_partInfo->GetParticleSymbol(a_partIndex).Data(),
            m_vZ_calib_struct->GetBinLowEdge(vZ_index+1),m_vZ_calib_struct->GetBinLowEdge(vZ_index+2),
            m_eta_calib_struct->GetBinLowEdge(eta_index+1),m_eta_calib_struct->GetBinLowEdge(eta_index+2)));        
        m_sigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetTitle(
          Form("%s ln(dE/dx) #sigma Calibration V_{Z}=[%4.2f,%4.2f] #eta=[%1.3f,%1.3f]",m_partInfo->GetParticleSymbol(a_partIndex).Data(),
            m_vZ_calib_struct->GetBinLowEdge(vZ_index+1),m_vZ_calib_struct->GetBinLowEdge(vZ_index+2),
            m_eta_calib_struct->GetBinLowEdge(eta_index+1),m_eta_calib_struct->GetBinLowEdge(eta_index+2)));
      }


      int mode = -1; // 0 is TPC and 1 is TOF
      vector<double> bg_index_index_vect;
      for(int bg_index = 0; bg_index < m_betaGamma_calib_struct->GetNbinsX(); bg_index++){
        bg_index_index_vect.push_back(bg_index);
      }
      for(int bg_index = 0; bg_index < m_betaGamma_calib_struct_fromTOF->GetNbinsX(); bg_index++){
        bg_index_index_vect.push_back(bg_index);
      }

      #ifdef _DEDX_CALIBRATION_DEBUG_
        cout << "bg_index_index_vect size = " << m_betaGamma_calib_struct->GetNbinsX() << " + " << m_betaGamma_calib_struct_fromTOF->GetNbinsX() << endl;
      #endif

      double previousGoodParameters[6] = {-999,-999,-999,-999,-999,-999};
      for(unsigned int bg_index_index = 0; bg_index_index < bg_index_index_vect.size(); bg_index_index++){
        int bg_index = bg_index_index_vect[bg_index_index];
        if(bg_index == 0){
          mode++;
        }
        if(mode == 0){
          m_currentHistoToFit = m_lnDeDx_calibrationHistograms[a_partIndex][vZ_index][eta_index][bg_index];
        }else if(mode == 1){
          m_currentHistoToFit = m_lnDeDx_calibrationHistograms_fromTOF[a_partIndex][vZ_index][eta_index][bg_index];
        }
        if(!m_currentHistoToFit) continue;
        double betaGamma = 0;
        if(mode == 0) betaGamma =  m_betaGamma_calib_struct->GetBinCenter(bg_index+1);
        if(mode == 1) betaGamma =  m_betaGamma_calib_struct_fromTOF->GetBinCenter(bg_index+1);
        if(mode == 0 && betaGamma > m_maxBetaGammaForTPCFitting[a_partIndex]) continue;

        #ifdef _DEDX_CALIBRATION_DEBUG_
          cout << "###########################################################################" << endl;
          cout << "    Fitting ln(dE/dx) for " << m_partInfo->GetParticleName(a_partIndex) 
               << "  vZ_index = " << vZ_index << " eta_index = " << eta_index << " bg_index = " << bg_index << endl;
          HistogramUtilities::printHistoInfo(m_currentHistoToFit);
        #endif



        bool   useTwoGaussians = false;
        double lowFitRange  = HistogramUtilities::getLowEdge(m_currentHistoToFit);
        double highFitRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);

        if(mode == 0){ // only contamination in pure dE/dx fiting
          if(a_partIndex == 0){
            if(betaGamma < 0.9)  lowFitRange   = -0.07;
            if(betaGamma > 1.1)  highFitRange  = 0.18;
            if(betaGamma > 1.5)  highFitRange  = 0.12;
            if(betaGamma > 4.5)  highFitRange = 0.07;
          }else if(a_partIndex == 1){
            if(betaGamma > 1.8)  lowFitRange = -0.2;
            if(betaGamma > 2.1)  lowFitRange = -0.12;
            if(betaGamma > 2.4)  lowFitRange = -0.08;
            //if(betaGamma > 3.0)  lowFitRange = -0.15;
            if(betaGamma > 3.4){
              lowFitRange = -1.0;
              useTwoGaussians = true;
            }
          }else if(a_partIndex == 2){
            if(betaGamma > 0.75)  lowFitRange = -0.18;
            if(betaGamma > 0.94)  lowFitRange = -0.12;
            if(betaGamma > 1.0)   lowFitRange = -0.075;
            //if(betaGamma > 0.41) highFitRange = 0.1;
            //if(betaGamma > 0.45) highFitRange = 0.07;
          }
        }

        double gausData[3];
        HistogramUtilities::getGausData(gausData, m_currentHistoToFit, true,true,-0.15,0.15);
        //cout << " Gaus Data: " << gausData[0] << " " << gausData[1] << " " << gausData[2] << endl;
        if(!std::isfinite(gausData[0])) gausData[0] = 1.0;
        if(!std::isfinite(gausData[1])) gausData[1] = 0.0;
        if(!std::isfinite(gausData[2])) gausData[2] = 0.07;
        string fitFunctName =  Form("gaus_%d_%02d_%02d_%03d",a_partIndex,vZ_index,eta_index,bg_index);
        TF1* fitFunct = NULL;
        if(useTwoGaussians){
          fitFunct =  PhysMath::getTwoGaussianFunction(fitFunctName.c_str(),gausData[0],gausData[1],gausData[2],0,0.0,0.068);
          if(previousGoodParameters[0] > 0){
            for(int iii = 0; iii < 6; iii++) fitFunct->SetParameter(iii,previousGoodParameters[iii]);
          }else{
            if(a_partIndex == 1){
              fitFunct->SetParameter(3,5.0*gausData[0]);
              fitFunct->SetParameter(4,-0.4);
              fitFunct->SetParameter(5,1.5*gausData[2]);
            }else if(a_partIndex == 2){
              fitFunct->SetParameter(3,gausData[0]);
              fitFunct->SetParameter(4,-0.4);
              fitFunct->SetParameter(5,1.5*gausData[2]);
            }
          }
        }else{
          fitFunct = PhysMath::getGaussianFunction(fitFunctName.c_str(),gausData[0],gausData[1],gausData[2]);
        }

        if(useTwoGaussians){
          fitFunct->SetParLimits(0, 0.01*gausData[0],30.0*gausData[0]);
          fitFunct->SetParLimits(1, gausData[1]-10.0*gausData[2],gausData[1]+10.0*gausData[2]);
          fitFunct->SetParLimits(2, 0.0, 10.0*gausData[2]);
          fitFunct->SetParLimits(3, 0.01*fitFunct->GetParameter(3),250.0*fitFunct->GetParameter(3));
          fitFunct->SetParLimits(4, fitFunct->GetParameter(4)-10.0*gausData[2],fitFunct->GetParameter(4)+10.0*gausData[2]);
          fitFunct->SetParLimits(5, 0.0, 10.0*gausData[2]);          
        }else{
          fitFunct->SetParLimits(0, 0.1*gausData[0],10.0*gausData[0]);
          fitFunct->SetParLimits(1, gausData[1]-3.0*gausData[2],gausData[1]+3.0*gausData[2]);
          fitFunct->SetParLimits(2, 0.0, 5.0*gausData[2]);
        }

        fitFunct->SetRange(-1.0,1.0);

        TF1* initialFunct = (TF1*) fitFunct->Clone();
        initialFunct->SetName("initialFunct");
        initialFunct->SetLineColor(kYellow+1);

        string fitOptions = "E I";
        #ifdef _DEDX_CALIBRATION_DEBUG_
           //fitOptions += "V";
           HistogramUtilities::printTF1Info(fitFunct);
           HistogramUtilities::printTF1Info(initialFunct);
           cout << "    ln(dE/dx) = [" << lowFitRange << "," << highFitRange << "]" << endl;
        #else
           fitOptions += " Q";
        #endif
        
        int status = m_currentHistoToFit->Fit(fitFunct,fitOptions.c_str(),"",lowFitRange,highFitRange);
        TF1*  fitResult = m_currentHistoToFit->GetFunction(fitFunctName.c_str());
        TH1D* pullHisto = NULL;
        if(fitResult){
          fitResult->SetRange(lowFitRange,highFitRange);
          double previousPrediction = log(m_initialCalibCurve->Eval(1.0/(betaGamma*betaGamma)));
          fitResult = m_currentHistoToFit->GetFunction(fitFunctName.c_str());

          double mean     = fitResult->GetParameter(1);
          double mean_err = fitResult->GetParError(1);
          double sigma    = fitResult->GetParameter(2);
          double sigma_err= fitResult->GetParError(2);

          if(useTwoGaussians){
            if(fabs(fitResult->GetParameter(4)) < fabs(fitResult->GetParameter(1))){
              mean     = fitResult->GetParameter(4);
              mean_err = fitResult->GetParError(4);
              sigma    = fitResult->GetParameter(5);
              sigma_err= fitResult->GetParError(5);              
            }
          }


          int pointIndex = m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]->GetN();
          m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPoint(pointIndex,betaGamma,mean + previousPrediction);
          m_meanPlusSigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPoint(pointIndex,betaGamma,mean + previousPrediction + sigma);
          m_sigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPoint(pointIndex,betaGamma,sigma);
          if(mode == 0){
            m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPointError(pointIndex,
              0.5*m_betaGamma_calib_struct->GetBinWidth(bg_index+1),0.5*m_betaGamma_calib_struct->GetBinWidth(bg_index+1),
              mean_err,mean_err);
            m_meanPlusSigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPointError(pointIndex,
              0.5*m_betaGamma_calib_struct->GetBinWidth(bg_index+1),0.5*m_betaGamma_calib_struct->GetBinWidth(bg_index+1),
              sqrt(mean_err*mean_err + sigma_err*sigma_err),sqrt(mean_err*mean_err + sigma_err*sigma_err));
            m_sigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPointError(pointIndex,
              0.5*m_betaGamma_calib_struct->GetBinWidth(bg_index+1),0.5*m_betaGamma_calib_struct->GetBinWidth(bg_index+1),
              sigma_err,sigma_err);
          }else if(mode == 1){
            m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPointError(pointIndex,
              0.5*m_betaGamma_calib_struct_fromTOF->GetBinWidth(bg_index+1),0.5*m_betaGamma_calib_struct_fromTOF->GetBinWidth(bg_index+1),
              mean_err,mean_err);
            m_meanPlusSigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPointError(pointIndex,
              0.5*m_betaGamma_calib_struct_fromTOF->GetBinWidth(bg_index+1),0.5*m_betaGamma_calib_struct_fromTOF->GetBinWidth(bg_index+1),
              sqrt(mean_err*mean_err + sigma_err*sigma_err),sqrt(mean_err*mean_err + sigma_err*sigma_err));
            m_sigma_calibrationGraphs[a_partIndex][vZ_index][eta_index]->SetPointError(pointIndex,
              0.5*m_betaGamma_calib_struct_fromTOF->GetBinWidth(bg_index+1),0.5*m_betaGamma_calib_struct_fromTOF->GetBinWidth(bg_index+1),
              sigma_err,sigma_err);            
          }

          pullHisto = HistogramUtilities::makeFitPullHistogram(m_currentHistoToFit, fitResult, false);
          pullHisto->SetMarkerStyle(20);
          pullHisto->SetMarkerSize(0.5);
          pullHisto->SetMarkerColor(kBlack);
          pullHisto->SetLineColor(kBlack);

          if(useTwoGaussians){
            for(int iii = 0; iii < 6; iii++) previousGoodParameters[iii] = fitResult->GetParameter(iii);
          }else{
            for(int iii = 0; iii < 6; iii++) previousGoodParameters[iii] = -999;
          }

          // status == 0
        }else{
            for(int iii = 0; iii < 6; iii++) previousGoodParameters[iii] = -999;
        }

        if(m_imageDir.length() > 1 && a_saveImages){
          TCanvas* fittingCanvas = new TCanvas(Form("fittingCanvas_%d_%02d_%02d_%03d",a_partIndex,vZ_index,eta_index,bg_index),"Fitting Canvas");
          fittingCanvas->SetWindowSize(1400,1250);
          fittingCanvas->SetCanvasSize(1350,1150);
          fittingCanvas->Draw();
          TPad* topPad_base    = new TPad(Form("top_%d_%02d_%02d_%03d",a_partIndex,vZ_index,eta_index,bg_index),"topPad",0.0,0.3,1.0,1.0);
          topPad_base->Draw();
          TPad* bottomPad_base = new TPad(Form("bottom_%d_%02d_%02d_%03d",a_partIndex,vZ_index,eta_index,bg_index),"bottomPad",0.0,0.0,1.0,0.3);
          bottomPad_base->Draw();


          double lowRange  = HistogramUtilities::getLowEdge(m_currentHistoToFit);
          double highRange = HistogramUtilities::getHighEdge(m_currentHistoToFit);
          topPad_base->cd();
          gPad->SetRightMargin(0.38);
          gPad->SetLogy(true);
          m_currentHistoToFit->Draw("E");
          m_currentHistoToFit->GetXaxis()->SetRangeUser(lowRange,highRange);
          initialFunct->Draw("SAME");
          if(fitResult) fitResult->Draw("SAME");
          gPad->Update();

          bottomPad_base->cd();
          gPad->SetRightMargin(0.38);
          TH1F *frame_pull_sys = bottomPad_base->DrawFrame(lowRange,-10.0,highRange,10.0);
          frame_pull_sys->SetTitle(Form(";ln(dE/dx) - c_{0}; Standardized Pulls"));
          if(pullHisto) pullHisto->Draw("PSame");
          TLine* oneLine = new TLine(lowRange,0.0,highRange,0.0);
          oneLine->SetLineStyle(7);
          oneLine->SetLineColor(kBlack);
          oneLine->Draw("same");
          gPad->Update();
          gSystem->ProcessEvents();
        
          fittingCanvas->SaveAs(Form("%s/%s/01_GaussianFits/dEdx_Fit_VzIndex_%02d_etaIndex_%02d%s_bgIndex_%03d%s.png",
              m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data(), vZ_index,eta_index,(mode==0 ? "" : "_fromTOF"),
              bg_index,(status==0 ? "" : "___FAILURE")));
          delete oneLine;
        }

        delete initialFunct;

        //if(fitFunct) delete fitFunct;
        //delete m_currentHistoToFit;
      
      } // BetaGamma Index Loop
    } // Eta Index Loop
  } // Vz Index Loop
}





//##############################################################################################################################################

void DeDxCalibration::calibrateMeans(int a_partIndex, bool a_doOneSigma){
  #ifdef _DEDX_CALIBRATION_DEBUG_
    cout << "Call to  DeDxCalibration::calibrateMeans(" << a_partIndex << ")" << endl;
  #endif
  bool useChrisFloresFunction = false;
  TVirtualFitter::SetMaxIterations( 1000000 );
  if(m_imageDir.length() > 1){
    if(a_doOneSigma){
      gSystem->mkdir(Form("%s/%s", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));
      gSystem->mkdir(Form("%s/%s/03_WidthCalibration", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));
    }else{
      gSystem->mkdir(Form("%s/%s", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));
      gSystem->mkdir(Form("%s/%s/02_MeanCalibration", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));      
    }

  }

  for(int vZ_index = 0; vZ_index < m_vZ_calib_struct->GetNbinsX(); vZ_index++){
    if(m_testing && vZ_index != 3) continue;
    for(int eta_index = 0; eta_index < m_eta_calib_struct->GetNbinsX(); eta_index++){
      if(m_testing && eta_index != 4) continue;
      m_currentGraphToFit = m_mean_calibrationGraphs[a_partIndex][vZ_index][eta_index];
      if(a_doOneSigma)  m_currentGraphToFit = m_meanPlusSigma_calibrationGraphs[a_partIndex][vZ_index][eta_index];
      if(!m_currentGraphToFit) continue;
      if(m_currentGraphToFit->GetN() < 5) continue;


      TGraphAsymmErrors* expGraph = new TGraphAsymmErrors();
      for(int pointIndex = 0; pointIndex < m_currentGraphToFit->GetN(); pointIndex++){
        double midpoint  = exp(m_currentGraphToFit->GetY()[pointIndex]);
        double lowpoint  = exp(m_currentGraphToFit->GetY()[pointIndex] - m_currentGraphToFit->GetEYlow()[pointIndex]);
        double highpoint = exp(m_currentGraphToFit->GetY()[pointIndex] + m_currentGraphToFit->GetEYhigh()[pointIndex]);
        expGraph->SetPoint(pointIndex,m_currentGraphToFit->GetX()[pointIndex],midpoint);
        expGraph->SetPointError(pointIndex,m_currentGraphToFit->GetEXlow()[pointIndex],m_currentGraphToFit->GetEXhigh()[pointIndex],
                                midpoint - lowpoint, highpoint - midpoint);
      }
      expGraph->SetTitle(m_currentGraphToFit->GetTitle());
      m_currentGraphToFit = expGraph;

      #ifdef _DEDX_CALIBRATION_DEBUG_
        cout << "###########################################################################" << endl;
        cout << "    Fitting Mean of ln(dE/dx) for " << m_partInfo->GetParticleName(a_partIndex) 
             << "  vZ_index = " << vZ_index << " eta_index = " << eta_index << endl;
        //HistogramUtilities::printGraphInfo(m_currentGraphToFit);
      #endif

      string fitFunctName =  Form("mean_%d_%02d_%02d",a_partIndex,vZ_index,eta_index);
      TF1* fitFunct = NULL;
      if(useChrisFloresFunction){
        fitFunct = PhysMath::getBichselFunctionBetaGamma_ChrisFlores(fitFunctName);
        /*fitFunct->SetParameter(0,);// A
        fitFunct->SetParameter(1,);// B
        fitFunct->SetParameter(2,);// C
        fitFunct->SetParameter(3,);// D
        fitFunct->SetParameter(4,);// E
        fitFunct->SetParameter(5,);// F
        fitFunct->SetParameter(6,);// G
        */
      }else{
        /*fitFunct = PhysMath::getBichselFunctionBetaGamma(fitFunctName);
        fitFunct->SetParameter(0,1.86476e+00);// A = P1*P2
        fitFunct->SetParameter(1,3.29351e-01);// B = P1
        fitFunct->SetParameter(2,1.88172e-02);// C = P3
        fitFunct->SetParameter(3,1.36804e+00);// D = P4/2
        fitFunct->SetParameter(4,1.09079e+00);// E = P5/2
        fitFunct->SetRange(0.0001,30.0);*/

        fitFunct = PhysMath::getBichselFunctionBetaGamma_Harasty(fitFunctName); 
      }

      TF1* initialFunct = (TF1*) fitFunct->Clone();
      string initialFunctName = Form("mean_%d_%02d_%02d_init",a_partIndex,vZ_index,eta_index);
      initialFunct->SetName(initialFunctName.c_str());
      initialFunct->SetLineColor(kYellow+1);

      double lowFitRange  = 0.0001;
      double highFitRange = 30.0;

      string fitOptions = "EX0 ";
      #ifdef _DEDX_CALIBRATION_DEBUG_
         fitOptions += "V";
         HistogramUtilities::printTF1Info(fitFunct);
      #else
         fitOptions += "Q";
      #endif
      
      //###########   FIRST ROUND FITTING   ###############
      int status_first = m_currentGraphToFit->Fit(fitFunct,fitOptions.c_str(),"",lowFitRange,highFitRange);
      TF1* fitFunct_Second = NULL;
      string fitFunct_SecondName = fitFunctName + "_Second";
      if(useChrisFloresFunction){
        fitFunct_Second = PhysMath::getBichselFunctionBetaGamma_ChrisFlores(fitFunct_SecondName);
        fitFunct_Second->SetParameter(0,fitFunct->GetParameter(0));
        fitFunct_Second->SetParameter(1,fitFunct->GetParameter(1));
        fitFunct_Second->SetParameter(2,fitFunct->GetParameter(2));
        fitFunct_Second->SetParameter(3,fitFunct->GetParameter(3));
        fitFunct_Second->SetParameter(4,fitFunct->GetParameter(4));
        fitFunct_Second->SetParameter(5,fitFunct->GetParameter(5));
        fitFunct_Second->SetParameter(6,fitFunct->GetParameter(6));
        fitFunct_Second->SetParLimits(0,0.01*fitFunct->GetParameter(0),100.0*fitFunct->GetParameter(0)); // A = P1*P2
        fitFunct_Second->SetParLimits(1,0.01*fitFunct->GetParameter(1),100.0*fitFunct->GetParameter(1)); // B = P1
        fitFunct_Second->SetParLimits(2,0.001*fitFunct->GetParameter(2),100.0*fitFunct->GetParameter(2));// C = P3
        fitFunct_Second->SetParLimits(3,0.01*fitFunct->GetParameter(3),100.0*fitFunct->GetParameter(3)); // D = P4/2
        fitFunct_Second->SetParLimits(4,0.01*fitFunct->GetParameter(4),100.0*fitFunct->GetParameter(4)); // E = P5/2
        fitFunct_Second->SetParLimits(5,0.01*fitFunct->GetParameter(5),100.0*fitFunct->GetParameter(5)); // D = P4/2
        fitFunct_Second->SetParLimits(6,0.01*fitFunct->GetParameter(6),100.0*fitFunct->GetParameter(6)); // E = P5/2
        fitFunct_Second->SetRange(0.0001,30.0);
        
      }else{
        /*fitFunct_Second = PhysMath::getBichselFunctionBetaGamma(fitFunct_SecondName);
        fitFunct_Second->SetParameter(0,fitFunct->GetParameter(0));// A = P1*P2
        fitFunct_Second->SetParameter(1,fitFunct->GetParameter(1));// B = P1
        fitFunct_Second->SetParameter(2,fitFunct->GetParameter(2));// C = P3
        fitFunct_Second->SetParameter(3,fitFunct->GetParameter(3));// D = P4/2
        fitFunct_Second->SetParameter(4,fitFunct->GetParameter(4));// E = P5/2
        fitFunct_Second->SetParLimits(0,0.01*fitFunct->GetParameter(0), 100.0*fitFunct->GetParameter(0)); // A = P1*P2
        fitFunct_Second->SetParLimits(1,0.01*fitFunct->GetParameter(1), 100.0*fitFunct->GetParameter(1)); // B = P1
        fitFunct_Second->SetParLimits(2,0.001*fitFunct->GetParameter(2),100.0*fitFunct->GetParameter(2)); // C = P3
        fitFunct_Second->SetParLimits(3,0.01*fitFunct->GetParameter(3), 100.0*fitFunct->GetParameter(3)); // D = P4/2
        fitFunct_Second->SetParLimits(4,0.01*fitFunct->GetParameter(4), 100.0*fitFunct->GetParameter(4)); // E = P5/2
        fitFunct_Second->SetRange(0.0001,30.0);*/

        fitFunct_Second = PhysMath::getBichselFunctionBetaGamma_Harasty(fitFunct_SecondName);
        for(int paramIndex = 0; paramIndex < 13; paramIndex++){
          double par = fitFunct->GetParameter(paramIndex);
          double par_low  = 0;
          double par_high = 0;
          fitFunct->GetParLimits(paramIndex, par_low,par_high);
          double par_width = par_high - par_low;
          fitFunct_Second->SetParameter(paramIndex,par);
          double new_par_low  = par_low;
          double new_par_high = par_high;
          if(par - par_low < 0.05*par_width){
            if(new_par_low > 0) new_par_low = 0.01*new_par_low;
            else        new_par_low = 100.0*new_par_low;
          }
          if(par_high - par < 0.05*par_width){
            if(new_par_high > 0) new_par_high = 100.0*new_par_high;
            else        new_par_high = 0.01*new_par_high;
          }          
          fitFunct_Second->SetParLimits(paramIndex,new_par_low,new_par_high);
        }

        fitFunct_Second->SetRange(0.0001,30.0);

      }      


      //#################   ROUND 2 FITTING ##################
      fitOptions += " S ";


      TFitResultPtr fitResultPtrPtr = m_currentGraphToFit->Fit(fitFunct_Second,fitOptions.c_str(),"",lowFitRange,highFitRange);
      TFitResult*   fitResultPtr = fitResultPtrPtr.Get();
      TF1*  fitResult              = NULL;
      TGraphAsymmErrors* pullGraph = NULL;
      TPaveText* fitParsTxt        = NULL;
      int status = -999;
      fitResult = m_currentGraphToFit->GetFunction(fitFunct_SecondName.c_str());
      if(fitResultPtr){
        status = fitResultPtr->Status();
        //m_mean_calibrationFunctions[a_partIndex][vZ_index][eta_index] = (TF1*) fitResult->Clone();

        TF1* funct_store = NULL;
        if(a_doOneSigma){
          m_meanPlusSigma_calibrationFunctions[a_partIndex][vZ_index][eta_index] = PhysMath::getLnBichselFunctionBetaGamma_Harasty(
            Form("mean_calibration_fit_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index));
          funct_store = m_meanPlusSigma_calibrationFunctions[a_partIndex][vZ_index][eta_index];
        }else{
          m_mean_calibrationFunctions[a_partIndex][vZ_index][eta_index] = PhysMath::getLnBichselFunctionBetaGamma_Harasty(
            Form("mean_calibration_fit_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index));
          funct_store = m_mean_calibrationFunctions[a_partIndex][vZ_index][eta_index];
        }
        //m_mean_calibrationFunctions[a_partIndex][vZ_index][eta_index] = PhysMath::getLnBichselFunctionBetaGamma(
        //    Form("mean_calibration_fit_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index));


        for(int paramIndex = 0; paramIndex < 13; paramIndex++){
          double par = fitFunct->GetParameter(paramIndex);
          fitFunct_Second->SetParameter(paramIndex,par);
        }
        funct_store->SetLineColor(kRed);
        funct_store->SetLineWidth(2);
        funct_store->SetRange(0.0001,30);

        pullGraph = HistogramUtilities::makeFitPullGraph(m_currentGraphToFit, fitResult, false);
        pullGraph->SetMarkerStyle(20);
        pullGraph->SetMarkerSize(0.5);
        pullGraph->SetMarkerColor(kBlack);
        pullGraph->SetLineColor(kBlack);

        fitParsTxt = new TPaveText(0.65,0.15,0.95,0.85,"NDC");
        fitParsTxt->AddText(Form("A  %.2e #pm %.2e",fitResult->GetParameter(0),fitResult->GetParError(0)));
        fitParsTxt->AddText(Form("B  %.2e #pm %.2e",fitResult->GetParameter(1),fitResult->GetParError(1)));
        fitParsTxt->AddText(Form("C  %.2e #pm %.2e",fitResult->GetParameter(2),fitResult->GetParError(2)));
        fitParsTxt->AddText(Form("D  %.2e #pm %.2e",fitResult->GetParameter(3),fitResult->GetParError(3)));
        fitParsTxt->AddText(Form("E  %.2e #pm %.2e",fitResult->GetParameter(4),fitResult->GetParError(4)));

        fitParsTxt->AddText(Form("F  %.2e #pm %.2e",fitResult->GetParameter(5),fitResult->GetParError(5)));
        fitParsTxt->AddText(Form("G  %.2e #pm %.2e",fitResult->GetParameter(6),fitResult->GetParError(6)));
        fitParsTxt->AddText(Form("H  %.2e #pm %.2e",fitResult->GetParameter(7),fitResult->GetParError(7)));
        fitParsTxt->AddText(Form("I  %.2e #pm %.2e",fitResult->GetParameter(8),fitResult->GetParError(8)));
        fitParsTxt->AddText(Form("J  %.2e #pm %.2e",fitResult->GetParameter(9),fitResult->GetParError(9)));
        fitParsTxt->AddText(Form("K  %.2e #pm %.2e",fitResult->GetParameter(10),fitResult->GetParError(10)));
        fitParsTxt->AddText(Form("L  %.2e #pm %.2e",fitResult->GetParameter(11),fitResult->GetParError(11)));

 

        fitParsTxt->SetFillStyle(0);
        fitParsTxt->SetBorderSize(5);
      }


      
      if(m_imageDir.length() > 1){
        TCanvas* fittingCanvas = new TCanvas(Form("fittingCanvas_%d_%02d_%02d_mean",a_partIndex,vZ_index,eta_index),"Fitting Canvas");
        fittingCanvas->SetWindowSize(1300,1250);
        fittingCanvas->SetCanvasSize(1200,1150);
        //fittingCanvas->Divide(1,2);
        fittingCanvas->Draw();
        TPad* topPad    = new TPad(Form("top_%d_%02d_%02d_mean",a_partIndex,vZ_index,eta_index),"topPad",0.0,0.3,1.0,1.0);
        topPad->Draw();
        TPad* bottomPad = new TPad(Form("bottom_%d_%02d_%02d_mean",a_partIndex,vZ_index,eta_index),"bottomPad",0.0,0.0,1.0,0.3);
        bottomPad->Draw();

        double lowRange  = 0.0;
        double highRange = 0.0;
        HistogramUtilities::horizontalDrawRange(m_currentGraphToFit, lowRange,highRange);

        double lowYRange  = 0;
        double highYRange = 0;
        HistogramUtilities::verticalDrawRange(m_currentGraphToFit, lowYRange,highYRange, -999, -999, true);

        topPad->cd();
        gPad->SetRightMargin(0.38);
        gPad->SetLogy(true);
        TH1D* topFrame = new TH1D("topFrame",Form("%s;#beta#gamma; #mu of ln(dE/dx)",m_currentGraphToFit->GetTitle()),1,lowRange,highRange);
        topFrame->SetBinContent(1,-1);
        topFrame->GetYaxis()->SetRangeUser(lowYRange,highYRange);
        topFrame->Draw("E");
        m_initialCalibCurve_betaGamma->Draw("SAME");
        initialFunct->Draw("SAME");
        m_currentGraphToFit->Draw("P");
        if(fitParsTxt) fitParsTxt->Draw("SAME");
        gPad->Update();
        //topPad->Draw();


        bottomPad->cd();
        //fittingCanvas->cd(2);
        gPad->SetRightMargin(0.38);
        gPad->SetLogy(false);
        TH1D* bottomFrame = new TH1D("bottomFrame",";#beta#gamma; Standardized Pulls",1,lowRange,highRange); 
        bottomFrame->SetBinContent(1,-20);
        bottomFrame->Draw("E");
        bottomFrame->GetYaxis()->SetRangeUser(-10,10);
        TLine* oneLine = new TLine(lowRange,0.0,highRange,0.0);
        oneLine->SetLineStyle(7);
        oneLine->SetLineColor(kBlack);
        oneLine->Draw("SAME");
        if(pullGraph) pullGraph->Draw("P");
        gPad->Update();
        gSystem->ProcessEvents();
        
        if(a_doOneSigma){
          fittingCanvas->SaveAs(Form("%s/%s/03_WidthCalibration/dEdx_MeanPlusSigma_Fit_VzIndex_%02d_etaIndex_%02d%s.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data(), vZ_index,eta_index, (status==0 ? "" : "___FAILURE")));
        }else{
          fittingCanvas->SaveAs(Form("%s/%s/02_MeanCalibration/dEdx_Fit_VzIndex_%02d_etaIndex_%02d%s.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data(), vZ_index,eta_index, (status==0 ? "" : "___FAILURE")));          
        }
        delete oneLine;
        //delete topFrame;
        //delete bottomFrame;
      }
      delete initialFunct;
      //if(fitFunct) delete fitFunct;
    } // Eta Index Loop
  } // Vz Index Loop




}


//##############################################################################################################################################

void DeDxCalibration::calibrateWidths(int a_partIndex){

  #ifdef _DEDX_CALIBRATION_DEBUG_
    cout << "Call to  DeDxCalibration::calibrateWidths(" << a_partIndex << ")" << endl;
  #endif
  if(m_imageDir.length() > 1){
    gSystem->mkdir(Form("%s/%s", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));
    gSystem->mkdir(Form("%s/%s/03_WidthCalibration", m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data()));
  }

  for(int vZ_index = 0; vZ_index < m_vZ_calib_struct->GetNbinsX(); vZ_index++){
    if(m_testing && vZ_index != 3) continue;
    for(int eta_index = 0; eta_index < m_eta_calib_struct->GetNbinsX(); eta_index++){
      if(m_testing && eta_index != 4) continue;
      m_currentGraphToFit = m_sigma_calibrationGraphs[a_partIndex][vZ_index][eta_index];
      if(!m_currentGraphToFit) continue;
      if(m_currentGraphToFit->GetN() < 5) continue;
      #ifdef _DEDX_CALIBRATION_DEBUG_
        cout << "###########################################################################" << endl;
        cout << "    Fitting Width of ln(dE/dx) for " << m_partInfo->GetParticleName(a_partIndex) 
             << "  vZ_index = " << vZ_index << " eta_index = " << eta_index << endl;
        HistogramUtilities::printGraphInfo(m_currentGraphToFit);
      #endif

      string fitFunctName =  Form("sigma_%d_%02d_%02d",a_partIndex,vZ_index,eta_index);
      //TF1* fitFunct = PhysMath::parabolaFunct(fitFunctName,0.067, 0.0, 0.0,1.0);
      //fitFunct->SetParLimits(0,0.0,0.2);
      //fitFunct->SetParLimits(1,-20,20);    
      //fitFunct->SetParLimits(2,-40,40); 
      //fitFunct->FixParameter(3,1.0);   

      /*
      TF1* fitFunct = new TF1(fitFunctName.c_str(), "[0]*pow([5]*x*x +[3]*x+[4],-[1]) + [2]", 0,15);
      fitFunct->SetParNames("Scale","Exponent","y_0","a","b");
      fitFunct->SetNpx(250);
      fitFunct->SetParameter(1,60);
      fitFunct->SetParameter(2,0.07);
      fitFunct->SetParameter(3,3.9E-2);
      fitFunct->SetParameter(4,0.9);
      fitFunct->SetParameter(5,0.0);
      //fitFunct->SetParLimits(0,0.0,1.0);
      fitFunct->SetParLimits(1,0,1000.0);
      fitFunct->SetParLimits(2,0.0,0.1);
      //fitFunct->SetParLimits(3,0.0,100.0);
      //fitFunct->SetParLimits(4,0.0,); // let it be free
      //fitFunct->SetParLimits(5,0.0,0.3);
      
      //do the seeding
      fitFunct->SetParameter(0,1.0);
      double xValToUse = HistogramUtilities::findSmallestXVal(m_currentGraphToFit);
      double ampToUse = m_currentGraphToFit->Eval(xValToUse) - fitFunct->GetParameter(2);
      if(ampToUse > 0.0){
        fitFunct->SetParameter(0,ampToUse/(fitFunct->Eval(xValToUse) - fitFunct->GetParameter(2)));
      }*/
      
      TF1* fitFunct = new TF1(fitFunctName.c_str(), "[0]", 0,15);
      fitFunct->SetParNames("y_0");
      fitFunct->SetNpx(250);
      fitFunct->SetParameter(0,0.07);


      TF1* initialFunct = (TF1*) fitFunct->Clone();
      initialFunct->SetName("initialFunct");
      initialFunct->SetLineColor(kYellow+1);

      double lowFitRange  = 0.0;
      double highFitRange = 30.0;

      string fitOptions = "EX0 ";
      #ifdef _DEDX_CALIBRATION_DEBUG_
         //fitOptions += "V";
         HistogramUtilities::printTF1Info(fitFunct);
         cout << "    ln(dE/dx) = [" << lowFitRange << "," << highFitRange << "]" << endl;
      #else
         fitOptions += "Q";
      #endif
        
      int status = m_currentGraphToFit->Fit(fitFunct,fitOptions.c_str(),"",lowFitRange,highFitRange);
      TF1*  fitResult = NULL;
      TGraphAsymmErrors* pullGraph = NULL;
      if(status == 0){
        fitResult = m_currentGraphToFit->GetFunction(fitFunctName.c_str());
        m_sigma_calibrationFunctions[a_partIndex][vZ_index][eta_index] = (TF1*) fitResult->Clone();
        m_sigma_calibrationFunctions[a_partIndex][vZ_index][eta_index]->SetName(
          Form("sigma_calibration_fit_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex).Data(),vZ_index,eta_index));

        pullGraph = HistogramUtilities::makeFitPullGraph(m_currentGraphToFit, fitResult, false);
        pullGraph->SetMarkerStyle(20);
        pullGraph->SetMarkerSize(0.5);
        pullGraph->SetMarkerColor(kBlack);
        pullGraph->SetLineColor(kBlack);
      } // status == 0

      if(m_imageDir.length() > 1){
        TCanvas* fittingCanvas = new TCanvas(Form("fittingCanvas_%d_%02d_%02d",a_partIndex,vZ_index,eta_index),"Fitting Canvas");
        fittingCanvas->SetWindowSize(1300,1250);
        fittingCanvas->SetCanvasSize(1200,1150);
        fittingCanvas->Draw();
        TPad* topPad_base    = new TPad(Form("top_%d_%02d_%02d",a_partIndex,vZ_index,eta_index),"topPad",0.0,0.3,1.0,1.0);
        TPad* bottomPad_base = new TPad(Form("bottom_%d_%02d_%02d",a_partIndex,vZ_index,eta_index),"bottomPad",0.0,0.0,1.0,0.3);
        gPad->SetRightMargin(0.38);
        gPad->SetRightMargin(0.38);
        topPad_base->Draw();
        bottomPad_base->Draw();


        double lowRange  = 0.0;
        double highRange = 0.0;
        HistogramUtilities::horizontalDrawRange(m_currentGraphToFit, lowRange,highRange);

        double lowYRange  = 0;
        double highYRange = 0;
        HistogramUtilities::verticalDrawRange(m_currentGraphToFit, lowYRange,highYRange);

        topPad_base->cd();
        TH1F *frame_pull = bottomPad_base->DrawFrame(lowRange,lowYRange,highRange,highYRange);
        frame_pull->SetTitle(Form(";#beta#gamma; #sigma of ln(dE/dx)"));         
        initialFunct->Draw("same"); 
        m_currentGraphToFit->Draw("P Same");
        if(fitResult) fitResult->Draw("same");

        bottomPad_base->cd();
        TH1F *frame_pull_bottom = bottomPad_base->DrawFrame(lowRange,-10.0,highRange,10.0);
        frame_pull_bottom->SetTitle(Form(";#beta#gamma; Standardized Pulls"));
        if(pullGraph) pullGraph->Draw("PSame");
        TLine* oneLine = new TLine(lowRange,0.0,highRange,0.0);
        oneLine->SetLineStyle(7);
        oneLine->SetLineColor(kBlack);
        oneLine->Draw("same");
        gSystem->ProcessEvents();
        
        fittingCanvas->SaveAs(Form("%s/%s/03_WidthCalibration/dEdx_Fit_VzIndex_%02d_etaIndex_%02d%s.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex).Data(), vZ_index,eta_index, (status==0 ? "" : "___FAILURE")));
        gSystem->ProcessEvents();
        delete oneLine;
      }
      delete initialFunct;
      //if(fitFunct) delete fitFunct;
    } // Eta Index Loop
  } // Vz Index Loop







}








//##############################################################################################################################################

void DeDxCalibration::makeCalibrationFile(string a_outputFileName){
  TFile* outFile = new TFile(a_outputFileName.c_str(),"RECREATE");
  HistogramUtilities::ConditionalWrite(m_initialCalibCurve);
  HistogramUtilities::ConditionalWrite(m_initialCalibCurve_betaGamma);
  HistogramUtilities::ConditionalWrite(m_vZ_calib_struct);
  HistogramUtilities::ConditionalWrite(m_eta_calib_struct);
  HistogramUtilities::ConditionalWrite(m_betaGamma_calib_struct);
  HistogramUtilities::ConditionalWrite(m_betaGamma_calib_struct_fromTOF);
  outFile->mkdir("calibration_functions");
  outFile->cd("calibration_functions");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int vZ_index = 0; vZ_index < g_DeDxCalib_MaxNumVzBins; vZ_index++){
      for(int eta_index = 0; eta_index < g_DeDxCalib_MaxNumEtaBins; eta_index++){
        HistogramUtilities::ConditionalWrite(m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]);
        if(m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]){
          string newName = m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]->GetName();
          newName += "_params";
          TH1D* param_histo = new TH1D(newName.c_str(),"Parameters",5,0,5);
          for(int iii = 0; iii < 5; iii++){
            param_histo->SetBinContent(iii+1,m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]->GetParameter(iii));
          }
          HistogramUtilities::ConditionalWrite(param_histo);
        }
        HistogramUtilities::ConditionalWrite(m_meanPlusSigma_calibrationFunctions[partIndex][vZ_index][eta_index]);
        HistogramUtilities::ConditionalWrite(m_sigma_calibrationFunctions[partIndex][vZ_index][eta_index]);
      }
    }
  }

  outFile->mkdir("calibration_graphs");
  outFile->cd("calibration_graphs");
  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int vZ_index = 0; vZ_index < g_DeDxCalib_MaxNumVzBins; vZ_index++){
      for(int eta_index = 0; eta_index < g_DeDxCalib_MaxNumEtaBins; eta_index++){
        HistogramUtilities::ConditionalWrite(m_mean_calibrationGraphs[partIndex][vZ_index][eta_index]);
        HistogramUtilities::ConditionalWrite(m_meanPlusSigma_calibrationGraphs[partIndex][vZ_index][eta_index]);
        HistogramUtilities::ConditionalWrite(m_sigma_calibrationGraphs[partIndex][vZ_index][eta_index]);
      }
    }
  }
}










#endif