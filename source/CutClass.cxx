#ifndef __CUT_CLASS_CXX__
#define __CUT_CLASS_CXX__

#include "../headers/CutClass.h"

CutClass::CutClass(){
  m_numCentralities = 0;
  m_useLogMTBinning      = false;
  m_variableMtM0BinEdges = false;
  m_cutOnTriggers        = false;
  m_useStRefMultCorr     = false;
  m_beamXCenter = 0;
  m_beamYCenter = 0;
  m_ZLowCut     = -400;
  m_ZHighCut    = 400;
  m_radialCut   = 100;
  m_DCA         = 100;
  m_pileUpCut   = 10000;
  m_pileUpLineCuts[0] = 1.0;
  m_pileUpLineCuts[1] = 100.0;
  m_pileUpLineCuts[2] = 0.0;
  m_pileUpLineCuts[3] = -10;
  m_pileUpLineCuts[4] = 0.0;
  m_pileUpLineCuts[5] = -10;

  m_tofMatchCut = -1; //Greater or Equal To

  m_useHui_dt_27GeVPileup = false;

  m_useBTOFBetaGammaCuts = false;
  m_usePhiCut    = false;
  m_lowPhiCut    = -30;
  m_highPhiCut   = -30;
  m_useT0        = false;
  m_nTofT0Cut    = -1;
  m_NHitsFit     = 0;
  m_NHitsDeDx    = 0;
  m_fitMaxRatio  = 0.5001;

  m_btofLocY = 5;
  m_btofLocZ = 5;

  m_etofDetlaXY_Center_X   = 0.0;
  m_etofDeltaXY_Center_Y   = 0.0;
  m_etofDeltaX_Ellipse_Cut = 20.0;
  m_etofDeltaY_Ellipse_Cut = 20.0;
  m_etof_max_cluster_size  = 99;
  m_etof_match_flag_cut    = -1;

  for(int iii = 0; iii < 9; iii++){
    m_invBetaCut_bTOF.push_back(0.03);
    m_invBetaCut_eTOF.push_back(0.03);
    m_nSigmaCuts.push_back(2.0);
  }
  m_eventConfig  = "ColliderCenter";
  m_vertexConfig = "Center";
  m_yCM        = 0.0;
  m_sqrts_NN   = 27;
  m_beamEnergy = 13.5;
  m_starver    = "SL19b";

  #ifdef _STREFMULTCORR_
    m_refMultCorr = CentralityMaker::instance()->getRefMultCorr();
    m_refMultCorr9Bins = true;
  #endif

   m_particleMasses[0] = 0.13957;    //Pion
   m_particleMasses[1] = 0.49367;    //Kaon
   m_particleMasses[2] = 0.93827;    //Proton
   m_particleMasses[3] = 0.00051099; //Electron
   m_particleMasses[4] = 1.87561;    //Deuteron
   m_particleMasses[5] = 2.80925;    //Triton
   m_particleMasses[6] = 2.8094;     //Helion
   m_particleMasses[7] = 3.72742;    //Alpha
   m_particleMasses[8] = 0.105658;    //Muon


  m_particleMasses_bTOF.push_back(0.13957);
  m_particleMasses_bTOF.push_back(0.49367);
  m_particleMasses_bTOF.push_back(0.93827);
  m_particleMasses_bTOF.push_back(0.00051099);
  m_particleMasses_bTOF.push_back(1.87561);
  m_particleMasses_bTOF.push_back(2.80925);
  m_particleMasses_bTOF.push_back(2.8094);
  m_particleMasses_bTOF.push_back(3.72742);
  m_particleMasses_bTOF.push_back(0.105658);

  m_particleMasses_eTOF.push_back(0.13957);
  m_particleMasses_eTOF.push_back(0.49367);
  m_particleMasses_eTOF.push_back(0.93827);
  m_particleMasses_eTOF.push_back(0.00051099);
  m_particleMasses_eTOF.push_back(1.87561);
  m_particleMasses_eTOF.push_back(2.80925);
  m_particleMasses_eTOF.push_back(2.8094);
  m_particleMasses_eTOF.push_back(3.72742);
  m_particleMasses_eTOF.push_back(0.105658);


  m_rejectedEvents = new TH1D("rejectedEventCuts","Cuts that rejected Events",9,0,9);
  m_rejectedEvents->GetXaxis()->SetBinLabel(1,"Z Vertex");
  m_rejectedEvents->GetXaxis()->SetBinLabel(2,"Radial Vertex");
  m_rejectedEvents->GetXaxis()->SetBinLabel(3,"bTOF Matched Tracks");
  m_rejectedEvents->GetXaxis()->SetBinLabel(4,"Num T0 Tracks");
  m_rejectedEvents->GetXaxis()->SetBinLabel(5,"Triggers");
  m_rejectedEvents->GetXaxis()->SetBinLabel(6,"Run Number");
  m_rejectedEvents->GetXaxis()->SetBinLabel(7,"Centrality Bin");
  m_rejectedEvents->GetXaxis()->SetBinLabel(8,"Pile-Up");
  m_rejectedEvents->GetXaxis()->SetBinLabel(9,"Passed");



  m_rejectedTracks = new TH1D("rejectedTrackCuts","Cuts that rejected Tracks",19,0,19);
  m_rejectedTracks->GetXaxis()->SetBinLabel(1,"Primary Track");
  m_rejectedTracks->GetXaxis()->SetBinLabel(2,"nHitsFit/nHitsMax");
  m_rejectedTracks->GetXaxis()->SetBinLabel(3,"nHitsFit");
  m_rejectedTracks->GetXaxis()->SetBinLabel(4,"nHitsDeDx");
  m_rejectedTracks->GetXaxis()->SetBinLabel(5,"DCA Cut");
  m_rejectedTracks->GetXaxis()->SetBinLabel(6,"Phi Cut");
  m_rejectedTracks->GetXaxis()->SetBinLabel(7,"Eta Cut");
  m_rejectedTracks->GetXaxis()->SetBinLabel(8,"Pt Cut");
  m_rejectedTracks->GetXaxis()->SetBinLabel(9,"BTOF Traits Exists");
  m_rejectedTracks->GetXaxis()->SetBinLabel(10,"BTOF bTOFMatchFlag==0");
  m_rejectedTracks->GetXaxis()->SetBinLabel(11,"BTOF LocalY");
  m_rejectedTracks->GetXaxis()->SetBinLabel(12,"BTOF LocalZ");
  m_rejectedTracks->GetXaxis()->SetBinLabel(13,"BTOF bTofBeta < 0");
  m_rejectedTracks->GetXaxis()->SetBinLabel(14,"BTOF BetaGamma_dEdx");
  m_rejectedTracks->GetXaxis()->SetBinLabel(15,"ETOF DeltaXY");
  m_rejectedTracks->GetXaxis()->SetBinLabel(16,"ETOF Match Flag");
  m_rejectedTracks->GetXaxis()->SetBinLabel(17,"ETOF Cluster Size");
  m_rejectedTracks->GetXaxis()->SetBinLabel(18,"ETOF BetaGamma_dEdx");
  m_rejectedTracks->GetXaxis()->SetBinLabel(19,"ETOF Beta");


  m_cutInformation = new TH1D("cutInformation","Cut Class Information",50,0,50);
  //m_cutInformation->GetXaxis()->SetBinLabel(0,"useT0");
  m_cutInformation->GetXaxis()->SetBinLabel(1,"Vx Center");
  m_cutInformation->GetXaxis()->SetBinLabel(2,"Vy Center");
  m_cutInformation->GetXaxis()->SetBinLabel(3,"Vz Lower");
  m_cutInformation->GetXaxis()->SetBinLabel(4,"Vz Higher");
  m_cutInformation->GetXaxis()->SetBinLabel(5,"Vr Cut");
  m_cutInformation->GetXaxis()->SetBinLabel(6,"bTOF Match Cut");
  m_cutInformation->GetXaxis()->SetBinLabel(7,"nTOF T0 Cut");
  m_cutInformation->GetXaxis()->SetBinLabel(8,"pile-up Cut");
  m_cutInformation->GetXaxis()->SetBinLabel(9,"Track DCA");
  m_cutInformation->GetXaxis()->SetBinLabel(10,"nHitsFit");
  m_cutInformation->GetXaxis()->SetBinLabel(11,"nHitsDeDx");
  m_cutInformation->GetXaxis()->SetBinLabel(12,"nHitsFit/nHitsMax");
  m_cutInformation->GetXaxis()->SetBinLabel(13,"bTOF Local Y");
  m_cutInformation->GetXaxis()->SetBinLabel(14,"bTOF Local Z");
  m_cutInformation->GetXaxis()->SetBinLabel(15,"eTOF Delta Cut");
  m_cutInformation->GetXaxis()->SetBinLabel(16,"eTOF DeltaXY X");
  m_cutInformation->GetXaxis()->SetBinLabel(17,"eTOF DeltaXY Y");
  m_cutInformation->GetXaxis()->SetBinLabel(18,"eTOF DeltaX Ellipse");
  m_cutInformation->GetXaxis()->SetBinLabel(19,"eTOF DeltaY Ellipse");
  m_cutInformation->GetXaxis()->SetBinLabel(20,"useT0Calc");
  m_cutInformation->GetXaxis()->SetBinLabel(21,"usePhiCut");
  m_cutInformation->GetXaxis()->SetBinLabel(22,"useEtaCut");
  m_cutInformation->GetXaxis()->SetBinLabel(23,"usePtCut");
  m_cutInformation->GetXaxis()->SetBinLabel(24,"useBTOFBetaGammaCuts");
  m_cutInformation->GetXaxis()->SetBinLabel(25,"useBTOFInvBetaGammaSqrCuts");
  m_cutInformation->GetXaxis()->SetBinLabel(26,"useETOFInvBetaGammaSqrCuts");
  m_cutInformation->GetXaxis()->SetBinLabel(27,"useModNSigma_InvBetaGammaSqrBTOFCuts");
  m_cutInformation->GetXaxis()->SetBinLabel(28,"useBTOFElectronCut");
  m_cutInformation->GetXaxis()->SetBinLabel(29,"useModNSigma_InvBetaGammaSqrETOFCuts");
  m_cutInformation->GetXaxis()->SetBinLabel(30,"etof_max_cluster_size");
  m_cutInformation->GetXaxis()->SetBinLabel(31,"useETOFElectronCut");


  m_cutInformation->GetXaxis()->SetBinLabel(32,"Trigger 1");
  m_cutInformation->GetXaxis()->SetBinLabel(33,"Trigger 2");
  m_cutInformation->GetXaxis()->SetBinLabel(34,"Trigger 3");
  m_cutInformation->GetXaxis()->SetBinLabel(35,"Trigger 4");
  m_cutInformation->GetXaxis()->SetBinLabel(36,"Trigger 5");
  m_cutInformation->GetXaxis()->SetBinLabel(37,"Trigger 6");
  m_cutInformation->GetXaxis()->SetBinLabel(38,"Trigger 7");
  m_cutInformation->GetXaxis()->SetBinLabel(39,"Trigger 8");
  m_cutInformation->GetXaxis()->SetBinLabel(40,"Trigger 9");
  m_cutInformation->GetXaxis()->SetBinLabel(41,"Trigger 10");









  for(int partIndex = 0; partIndex<9; partIndex++){
    m_mTm0_structs.push_back(NULL);
    m_tpc_rap_structs.push_back(NULL);
    m_btof_rap_structs.push_back(NULL);
    m_etof_rap_structs.push_back(NULL);
    m_btof_mTm0_structs.push_back(NULL);
    m_etof_mTm0_structs.push_back(NULL);
    for(int p2 = 0; p2 < 9; p2++){
      if(partIndex == p2){
        m_particleFitMatrix_Plus[partIndex][p2] = true;
        m_particleFitMatrix_Minus[partIndex][p2] = true;
      }else{
        m_particleFitMatrix_Plus[partIndex][p2] = false;
        m_particleFitMatrix_Minus[partIndex][p2] = false;
      }
    }
  }

  m_particleFitMatrix_Plus[0][0] = true;
  m_particleFitMatrix_Plus[0][1] = true;
  m_particleFitMatrix_Plus[0][2] = true;

  m_particleFitMatrix_Plus[1][0] = true;
  m_particleFitMatrix_Plus[1][1] = true;
  m_particleFitMatrix_Plus[1][2] = true;

  m_particleFitMatrix_Plus[2][0] = true;
  m_particleFitMatrix_Plus[2][1] = true;
  m_particleFitMatrix_Plus[2][2] = true;

  m_particleFitMatrix_Minus[0][0] = true;
  m_particleFitMatrix_Minus[0][1] = true;
  m_particleFitMatrix_Minus[0][2] = true;

  m_particleFitMatrix_Minus[1][0] = true;
  m_particleFitMatrix_Minus[1][1] = true;
  m_particleFitMatrix_Minus[1][2] = true;

  m_particleFitMatrix_Minus[2][0] = true;
  m_particleFitMatrix_Minus[2][1] = true;
  m_particleFitMatrix_Minus[2][2] = true;

  m_useModNSigma_InvBetaGammaSqrBTOFCuts = false;
  m_lowInvBetaGammaSqr_ToUseModNSigmaBTOFCut = 500.0;
  m_lowModNSigma_BTOF_Cut = -999;
  m_highModNSigma_BTOF_Cut = 999;

  m_useModNSigma_InvBetaGammaSqrETOFCuts = false;
  m_lowInvBetaGammaSqr_ToUseModNSigmaETOFCut = 500.0;
  m_lowModNSigma_ETOF_Cut = -999;
  m_highModNSigma_ETOF_Cut = 999;

  m_lowInvBetaGammaSqrBTOFCut = -999;
  m_lowInvBetaGammaSqrETOFCut = -999;

  m_useBTOFElectronCut = false;
  m_electron_lowInvBetaGammaSqrBTOFCut = -0.08;
  m_electron_highInvBetaGammaSqrBTOFCut = 0.08;
  m_electron_lowDeDxBTOFCut = 3.55;
  m_electron_highDeDxBTOFCut = 5.0;

  m_useETOFElectronCut = false;
  m_electron_lowInvBetaGammaSqrETOFCut = -0.08;
  m_electron_highInvBetaGammaSqrETOFCut = 0.08;
  m_electron_lowDeDxETOFCut = 3.55;
  m_electron_highDeDxETOFCut = 5.0;

  // create the calibration function
  m_invBetaGammaSqr_dEdx_Calib_funct = new TF1("cutClassCalibFunct","pow(1.0 + x,[3])*([0] - [1]*log([2] + pow(x,[4]))) - [1]",0.0,50);
  m_invBetaGammaSqr_dEdx_Calib_funct->SetParNames("A","B","C","D","E","F");
  m_invBetaGammaSqr_dEdx_Calib_funct->SetParameter(0,1.97740);// A = P1*P2
  m_invBetaGammaSqr_dEdx_Calib_funct->SetParameter(1,2.76679e-01);// B = P1
  m_invBetaGammaSqr_dEdx_Calib_funct->SetParameter(2,1.72493e-02); // C = P3
  m_invBetaGammaSqr_dEdx_Calib_funct->SetParameter(3,1.23514); // D = P4/2
  m_invBetaGammaSqr_dEdx_Calib_funct->SetParameter(4,1.12612); // E = P5/2
  m_invBetaGammaSqr_dEdx_Calib_funct->SetNpx(400);

}


void CutClass::makeCutInformationHistogram(){
  m_cutInformation->SetBinContent(1,m_beamXCenter);
  m_cutInformation->SetBinContent(2,m_beamYCenter);
  m_cutInformation->SetBinContent(3,m_ZLowCut);
  m_cutInformation->SetBinContent(4,m_ZHighCut);
  m_cutInformation->SetBinContent(5,m_radialCut);
  m_cutInformation->SetBinContent(6,m_tofMatchCut);
  m_cutInformation->SetBinContent(7,m_nTofT0Cut);
  m_cutInformation->SetBinContent(8,m_pileUpCut);
  m_cutInformation->SetBinContent(9,m_DCA);
  m_cutInformation->SetBinContent(10,m_NHitsFit);
  m_cutInformation->SetBinContent(11,m_NHitsDeDx);
  m_cutInformation->SetBinContent(12,m_fitMaxRatio);
  m_cutInformation->SetBinContent(13,m_btofLocY);
  m_cutInformation->SetBinContent(14,m_btofLocZ);
  m_cutInformation->SetBinContent(15,m_etofDeltaCut);
  m_cutInformation->SetBinContent(16,m_etofDetlaXY_Center_X);
  m_cutInformation->SetBinContent(17,m_etofDeltaXY_Center_Y);
  m_cutInformation->SetBinContent(18,m_etofDeltaX_Ellipse_Cut);
  m_cutInformation->SetBinContent(19,m_etofDeltaY_Ellipse_Cut);
  m_cutInformation->SetBinContent(20,m_useT0);
  m_cutInformation->SetBinContent(21,m_usePhiCut);
  m_cutInformation->SetBinContent(22,m_useEtaCut);
  m_cutInformation->SetBinContent(23,m_usePtCut);
  m_cutInformation->SetBinContent(24,m_useBTOFBetaGammaCuts);
  m_cutInformation->SetBinContent(25,m_useBTOFInvBetaGammaSqrCuts);
  m_cutInformation->SetBinContent(26,m_useETOFInvBetaGammaSqrCuts);
  m_cutInformation->SetBinContent(27,m_useModNSigma_InvBetaGammaSqrBTOFCuts);
  m_cutInformation->SetBinContent(28,m_useBTOFElectronCut);
  m_cutInformation->SetBinContent(29,m_useModNSigma_InvBetaGammaSqrETOFCuts);
  m_cutInformation->SetBinContent(30,m_etof_max_cluster_size);
  m_cutInformation->SetBinContent(31,m_useETOFElectronCut);

  for(unsigned int iii = 0; iii < m_triggers.size(); iii++){
    m_cutInformation->SetBinContent(32 + iii, m_triggers[iii]);
  }


}




void CutClass::setZRange(double a_zLow, double a_zHigh){
  m_ZLowCut = a_zLow;
  m_ZHighCut = a_zHigh;
  //cout << m_ZLowCut << " = " << a_zLow << endl << m_ZHighCut << " = " << a_zHigh << endl;
}

std::pair<double,double> CutClass::getZVertexCuts(){
  std::pair<double,double> zVtx;
  zVtx.first = m_ZLowCut;
  zVtx.second = m_ZHighCut;
  return zVtx;
}

std::pair<double,double> CutClass::getBeamSpotLocation(){
  std::pair<double,double> beamXY;
  beamXY.first = m_beamXCenter;
  beamXY.second = m_beamYCenter;
  return beamXY;
}


double CutClass::getLowMassSquaredCut(int a_partIndex){
  if(m_lowMassCut.at(a_partIndex) > 0.0) return pow(m_lowMassCut.at(a_partIndex),2.0);
  return -pow(m_lowMassCut.at(a_partIndex),2.0);
}

double CutClass::getHighMassSquaredCut(int a_partIndex){
  if(m_highMassCut.at(a_partIndex) > 0.0) return pow(m_highMassCut.at(a_partIndex),2.0);
  return -pow(m_highMassCut.at(a_partIndex),2.0);
}

/*
void CutClass::setBeamLocation(double a_x, double a_y){
  m_beamXCenter = a_x;
  m_beamYCenter = a_y;
  //cout << m_ZLowCut << " = " << a_zLow << endl << m_ZHighCut << " = " << a_zHigh << endl;
  }*/


void CutClass::setTriggers(int a_numTriggers, unsigned int* a_triggers){
  for(int iii = 0; iii < a_numTriggers; ++iii){
    m_triggers.push_back(a_triggers[iii]);
  }
  return;
}


void CutClass::setHLTTriggers(int a_numHLTTriggers, unsigned int* a_HLTTriggers){
  for(int iii = 0; iii < a_numHLTTriggers; ++iii){
    m_HLTTriggers.push_back(a_HLTTriggers[iii]);
  }
  return;
}
#ifdef _HAS_ETOF_
void CutClass::setEToFTriggers(int a_numEToFTriggers, unsigned int* a_EToFTriggers){
  for(int iii = 0; iii < a_numEToFTriggers; ++iii){
    m_EToFTriggers.push_back(a_EToFTriggers[iii]);
  }
  return;
}
#endif




#ifdef _STREFMULTCORR_

void CutClass::setStRefMultCorrEvent(StPicoEvent* a_event){
#ifdef _CUTCLASS_DEBUG_
  cout << "initializing StRefMultCorr with refMult:" << a_event->refMult()
  << "  Vz:" <<a_event->primaryVertex().z() << " ZDCx:" << a_event->ZDCx() << endl;;
#endif

  m_refMultCorr->initEvent(a_event->refMult(),
    a_event->primaryVertex().z(),
    a_event->ZDCx());
}


void CutClass::setStRefMultCorrEvent(int a_refMult, double a_Vz, double a_ZDCx){
#ifdef _CUTCLASS_DEBUG_
  cout << "initializing StRefMultCorr with refMult:" << a_refMult
  << "  Vz:" <<a_Vz << " ZDCx:" << a_ZDCx << endl;;
#endif

  m_refMultCorr->initEvent(a_refMult,a_Vz,a_ZDCx);
}

void CutClass::setStRefMultCorrCentralities(){
//  Centrality binning:
//  //     Bin       Centrality (16)   Centrality (9)
//  //     0            75-80%            70-80%
//  //     1            70-75%            60-70%
//  //     2            65-70%            50-60%
//  //     3            60-65%            40-50%
//  //     4            55-60%            30-40%
//  //     5            50-55%            20-30%
//  //     6            45-50%            10-20%
//  //     7            40-45%             5-10%
//  //     8            35-40%             0- 5%
//  //     9            30-35%
//  //    10            25-30%
//  //    11            20-25%
//  //    12            15-20%
//  //    13            10-15%
//  //    14             5-10%
//  //    15             0- 5%

  m_centPercents.clear();
  m_centEdges.clear();

  if(m_refMultCorr9Bins){
    m_numCentralities = 9;
    m_centPercents.push_back(0);
    m_centPercents.push_back(5);
    m_centPercents.push_back(10);
    m_centPercents.push_back(20);
    m_centPercents.push_back(30);
    m_centPercents.push_back(40);
    m_centPercents.push_back(50);
    m_centPercents.push_back(60);
    m_centPercents.push_back(70);
    m_centPercents.push_back(80);
  }else{
    m_numCentralities = 16;
    m_centPercents.push_back(0);
    m_centPercents.push_back(5);
    m_centPercents.push_back(10);
    m_centPercents.push_back(15);
    m_centPercents.push_back(20);
    m_centPercents.push_back(25);
    m_centPercents.push_back(30);
    m_centPercents.push_back(35);
    m_centPercents.push_back(40);
    m_centPercents.push_back(45);
    m_centPercents.push_back(50);
    m_centPercents.push_back(55);
    m_centPercents.push_back(60);
    m_centPercents.push_back(65);
    m_centPercents.push_back(70);
    m_centPercents.push_back(75);
    m_centPercents.push_back(80);
  }

#ifdef _CUTCLASS_DEBUG_
  cout << "Centralities: " << m_numCentralities << endl;
  for(int iii = 0; iii < m_centPercents.size(); iii++){
   cout << m_centPercents.at(iii) << "  ";
 }
 cout << endl;
#endif

}
#endif















void CutClass::setCentralities(int a_numCent, double* a_cent_edges, int* a_cent_percents){
  m_numCentralities = a_numCent;
  m_centPercents.clear();
  m_centEdges.clear();
  m_centPercents.push_back(0.0);
  if(a_cent_edges[0] != 0) m_centEdges.push_back(1000.0*a_cent_edges[0]);
  else m_centEdges.push_back(100000.0);
  for (int iii = 0; iii < a_numCent; iii++){
    m_centPercents.push_back(a_cent_percents[iii]);
    m_centEdges.push_back(a_cent_edges[iii]);
    //    cout << m_centPercents[iii] << " " << m_centEdges[iii] << endl;
  }
}


void CutClass::setYCM(Float_t a_yCM){
  if (m_eventConfig == "ColliderCenter"){
    m_yCM = 0.0;
    return;
  }
  else if (m_eventConfig == "FixedTarget"){
    if (a_yCM > 0) m_yCM = -a_yCM;
    else m_yCM = a_yCM;
  }
}


void CutClass::setEtaPtBinningInfo(  int a_numPtBins, double a_pTMin, double a_pTMax,
  int a_numEtaBins_TPC, double a_EtaMin_TPC, double a_EtaMax_TPC,
  int a_numEtaBins_bTOF, double a_EtaMin_bTOF, double a_EtaMax_bTOF,
  int a_numEtaBIns_eTOF, double a_EtaMin_eTOF, double a_EtaMax_eTOF){

  m_pT_struct = HistogramUtilities::make1DHistLogBins("pt_bin_struct","pt_bin_struct",a_numPtBins,a_pTMin,a_pTMax);
  m_eta_tpc_struct = new TH1D("tpc_bin_struct","tpc_bin_struct",a_numEtaBins_TPC,a_EtaMin_TPC,a_EtaMax_TPC);
  m_eta_btof_struct = new TH1D("btof_bin_struct","btof_bin_struct",a_numEtaBins_bTOF,a_EtaMin_bTOF,a_EtaMax_bTOF);
  m_eta_etof_struct = new TH1D("etof_bin_struct","etof_bin_struct",a_numEtaBIns_eTOF,a_EtaMin_eTOF,a_EtaMax_eTOF);
}


TH1D* CutClass::getEtaPtBinStructure(int a_tag){
  switch(a_tag){
    case 0: return m_pT_struct;
    case 1: return m_eta_tpc_struct;
    case 2: return m_eta_btof_struct;
    case 3: return m_eta_etof_struct;
    default : cout << "you did not select a proper tag in CutClass::getBinStructure... 0=pt, eta-> 1=tpc, 2=btof, 3=etof " << endl; return NULL;
  }
  return NULL;
}


void CutClass::setRapMtM0BinningInfo(
 vector<int> a_numMtM0Bins, vector<double> a_mTm0Mins, vector<double> a_mTm0Maxs,
 vector<int> a_numRapBins_tpc, vector<double> a_RapMins_tpc, vector<double> a_RapMaxs_tpc,
 vector<int> a_numRapBins_btof, vector<double> a_RapMins_btof, vector<double> a_RapMaxs_btof,
 vector<int> a_numRapBins_etof, vector<double> a_RapMins_etof, vector<double> a_RapMaxs_etof){

  for(int partIndex = 0; partIndex<9; partIndex++){
    m_mTm0_structs.push_back(NULL);
    m_tpc_rap_structs.push_back(NULL);
    m_btof_rap_structs.push_back(NULL);
    m_etof_rap_structs.push_back(NULL);

    if(m_useLogMTBinning){
      m_mTm0_structs[partIndex] = HistogramUtilities::make1DHistLogBins(Form("mTm0_bin_struct_part%02d",partIndex),"",a_numMtM0Bins[partIndex],a_mTm0Mins[partIndex],a_mTm0Maxs[partIndex]);
    }else{
      m_mTm0_structs[partIndex] = new TH1D(Form("mTm0_bin_struct_part%02d",partIndex),"",a_numMtM0Bins[partIndex],a_mTm0Mins[partIndex],a_mTm0Maxs[partIndex]);
    }
    m_tpc_rap_structs[partIndex]= new TH1D(Form("tpc_rap_bin_struct_part%02d",partIndex),"",a_numRapBins_tpc[partIndex],a_RapMins_tpc[partIndex],a_RapMaxs_tpc[partIndex]);
    m_btof_rap_structs[partIndex] = new TH1D(Form("btof_rap_bin_struct_part%02d",partIndex),"",a_numRapBins_btof[partIndex],a_RapMins_btof[partIndex],a_RapMaxs_btof[partIndex]);
    m_etof_rap_structs[partIndex] = new TH1D(Form("etof_rap_bin_struct_part%02d",partIndex),"",a_numRapBins_etof[partIndex],a_RapMins_etof[partIndex],a_RapMaxs_etof[partIndex]);
  }
  return;
}

void CutClass::setVariableRapMtM0BinningInfo(int a_partIndex,
                                   int a_numRapBins_tpc, double a_lowRap_tpc, double a_highRap_tpc,
                                   int a_numMtM0Bins_tpc, vector<double> a_mTm0BinEdges_tpc,
                                   int a_numRapBins_btof, double a_lowRap_btof, double a_highRap_btof,
                                   int a_numMtM0Bins_btof, vector<double> a_mTm0BinEdges_btof,
                                   int a_numRapBins_etof, double a_lowRap_etof, double a_highRap_etof,
                                   int a_numMtM0Bins_etof, vector<double> a_mTm0BinEdges_etof){

  m_variableMtM0BinEdges = true;
  m_tpc_rap_structs[a_partIndex] = new TH1D(Form("tpc_rap_bin_struct_part%02d",a_partIndex),"",a_numRapBins_tpc,a_lowRap_tpc,a_highRap_tpc);
  m_btof_rap_structs[a_partIndex]  = new TH1D(Form("btof_rap_bin_struct_part%02d",a_partIndex),"",a_numRapBins_btof,a_lowRap_btof,a_highRap_btof);
  m_etof_rap_structs[a_partIndex]  = new TH1D(Form("etof_rap_bin_struct_part%02d",a_partIndex),"",a_numRapBins_etof,a_lowRap_etof,a_highRap_etof);

  m_mTm0_structs[a_partIndex] = new TH1D(Form("mTm0_tpc_bin_struct_part%02d",a_partIndex),"",a_numMtM0Bins_tpc,&a_mTm0BinEdges_tpc[0]);
  m_btof_mTm0_structs[a_partIndex] = new TH1D(Form("mTm0_bTof_bin_struct_part%02d",a_partIndex),"",a_numMtM0Bins_btof,&a_mTm0BinEdges_btof[0]);
  if(a_mTm0BinEdges_etof.size() > 0) m_etof_mTm0_structs[a_partIndex] = new TH1D(Form("mTm0_etof_bin_struct_part%02d",a_partIndex),"",a_numMtM0Bins_etof,&a_mTm0BinEdges_etof[0]);

}

void CutClass::setVariableRapMtM0BinningInfo(int a_partIndex,
                                     int a_numRapBins_tpc, double a_lowRap_tpc, double a_highRap_tpc,
                                     int a_numMtM0Bins_tpc, double* a_mTm0BinEdges_tpc,
                                     int a_numRapBins_btof, double a_lowRap_btof, double a_highRap_btof,
                                     int a_numMtM0Bins_btof, double* a_mTm0BinEdges_btof,
                                     int a_numRapBins_etof, double a_lowRap_etof, double a_highRap_etof,
                                     int a_numMtM0Bins_etof, double* a_mTm0BinEdges_etof){
  m_variableMtM0BinEdges = true;
  m_tpc_rap_structs[a_partIndex] = new TH1D(Form("tpc_rap_bin_struct_part%02d",a_partIndex),"",a_numRapBins_tpc,a_lowRap_tpc,a_highRap_tpc);
  m_btof_rap_structs[a_partIndex]  = new TH1D(Form("btof_rap_bin_struct_part%02d",a_partIndex),"",a_numRapBins_btof,a_lowRap_btof,a_highRap_btof);
  m_etof_rap_structs[a_partIndex]  = new TH1D(Form("etof_rap_bin_struct_part%02d",a_partIndex),"",a_numRapBins_etof,a_lowRap_etof,a_highRap_etof);

  m_mTm0_structs[a_partIndex] = new TH1D(Form("mTm0_tpc_bin_struct_part%02d",a_partIndex),"",a_numMtM0Bins_tpc,a_mTm0BinEdges_tpc);
  m_btof_mTm0_structs[a_partIndex] = new TH1D(Form("mTm0_bTof_bin_struct_part%02d",a_partIndex),"",a_numMtM0Bins_btof,a_mTm0BinEdges_btof);
  if(a_numMtM0Bins_etof > 0) m_etof_mTm0_structs[a_partIndex] = new TH1D(Form("mTm0_etof_bin_struct_part%02d",a_partIndex),"",a_numMtM0Bins_etof,a_mTm0BinEdges_etof);

}

void CutClass::setParticleFitMatrix(bool a_particleFitMatrix[9][9], int a_charge){

  if(a_charge>0){
    for(int partBaseIndex=0; partBaseIndex<9; partBaseIndex++){
      for(int partIndex=0; partIndex<9; partIndex++){
        m_particleFitMatrix_Plus[partBaseIndex][partIndex] = a_particleFitMatrix[partBaseIndex][partIndex];
      }
    }
  }
  else if(a_charge<0){
    for(int partBaseIndex=0; partBaseIndex<9; partBaseIndex++){
      for(int partIndex=0; partIndex<9; partIndex++){
        m_particleFitMatrix_Minus[partBaseIndex][partIndex] = a_particleFitMatrix[partBaseIndex][partIndex];
      }
    }
  }

  return;
}

void CutClass::setBTOFBetaGammaCuts(bool a_useCuts, vector<double> a_lowLineXPts, vector<double> a_lowLineYPts,
				    vector<double> a_highLineXPts, vector<double> a_highLineYPts){
  for(int iii = 0; iii < a_lowLineXPts.size(); iii+=2){
    double slope = (log(a_lowLineYPts[iii]) - log(a_lowLineYPts[iii+1]))/(log(a_lowLineXPts[iii]) - log(a_lowLineXPts[iii+1]));
    double intercept = log(a_lowLineYPts[iii]) - slope * log(a_lowLineXPts[iii]);
    m_BTOFBetaGammaLowLineSlopes.push_back(slope);
    m_BTOFBetaGammaLowLineIntercepts.push_back(intercept);
  }

  for(int iii = 0; iii < a_highLineXPts.size(); iii++){
    m_BTOFBetaGammaHighLineXPoints.push_back(log(a_highLineXPts[iii]));
    m_BTOFBetaGammaHighLineYPoints.push_back(log(a_highLineYPts[iii]));
  }

  m_useBTOFBetaGammaCuts = a_useCuts;

}

void CutClass::setBTOFInvBetaGammaSqrCuts(bool a_useCuts, vector<double> a_lowLineXPts, vector<double> a_lowLineYPts,
                                                          vector<double> a_highLineXPts, vector<double> a_highLineYPts){
  for(int iii = 0; iii < a_lowLineXPts.size(); iii+=2){
    double slope = (a_lowLineYPts[iii] - a_lowLineYPts[iii+1])/(a_lowLineXPts[iii] - a_lowLineXPts[iii+1]);
    double intercept = a_lowLineYPts[iii] - slope * a_lowLineXPts[iii];
    m_BTOFInvBetaGammaSqrLowLineSlopes.push_back(slope);
    m_BTOFInvBetaGammaSqrLowLineIntercepts.push_back(intercept);
  }

  for(int iii = 0; iii < a_highLineXPts.size(); iii++){
    m_BTOFInvBetaGammaSqrHighLineXPoints.push_back(a_highLineXPts[iii]);
    m_BTOFInvBetaGammaSqrHighLineYPoints.push_back(a_highLineYPts[iii]);
  }
  m_useBTOFInvBetaGammaSqrCuts = a_useCuts;

}

void CutClass::setETOFInvBetaGammaSqrCuts(bool a_useCuts, vector<double> a_lowLineXPts, vector<double> a_lowLineYPts,
                                                          vector<double> a_highLineXPts, vector<double> a_highLineYPts){
  for(int iii = 0; iii < a_lowLineXPts.size(); iii+=2){
    double slope = (a_lowLineYPts[iii] - a_lowLineYPts[iii+1])/(a_lowLineXPts[iii] - a_lowLineXPts[iii+1]);
    double intercept = a_lowLineYPts[iii] - slope * a_lowLineXPts[iii];
    m_ETOFInvBetaGammaSqrLowLineSlopes.push_back(slope);
    m_ETOFInvBetaGammaSqrLowLineIntercepts.push_back(intercept);
  }

  for(int iii = 0; iii < a_highLineXPts.size(); iii++){
    m_ETOFInvBetaGammaSqrHighLineXPoints.push_back(a_highLineXPts[iii]);
    m_ETOFInvBetaGammaSqrHighLineYPoints.push_back(a_highLineYPts[iii]);
  }

  m_useETOFInvBetaGammaSqrCuts = a_useCuts;

}


#ifndef _MAC_OSX_
bool CutClass::isGoodBTof(StPicoBTofPidTraits* a_traits){
  if(!a_traits){
    m_rejectedTracks->Fill(8.5);
    return false;
  }
  if(a_traits->btofMatchFlag() == 0){
    m_rejectedTracks->Fill(9.5);
    return false;
  }
  if(fabs(a_traits->btofYLocal()) > m_btofLocY){
    m_rejectedTracks->Fill(10.5);
    return false;
  }
  if(fabs(a_traits->btofZLocal()) > m_btofLocZ){
    m_rejectedTracks->Fill(11.5);
    return false;
  }
  if(a_traits->btofBeta() <= 0){
    m_rejectedTracks->Fill(12.5);
    return false;
  }

  //cout << " In isGoodBTof" << endl;
  //this cuts out double pion and double kaon peaks
  /*if(a_dEdx != -999){
    if(m_useBTOFBetaGammaCuts && !isGoodTrack_GivenBTOFBetaGamma(a_dEdx,PhysMath::beta_to_betagamma(a_traits->btofBeta())) ){
      return false;
    }
    double invBGSqr = PhysMath::invBetaGammaSqr_signed_fromInvBeta(1.0/a_traits->btofBeta());
    if(invBGSqr < m_lowInvBetaGammaSqrBTOFCut){
      return false;
    }

    if(m_useBTOFInvBetaGammaSqrCuts && !isGoodTrack_GivenBTOFInvBetaGammaSqr(a_dEdx,invBGSqr) ){
      return false;
    }

  }*/
  return true;
}

#endif


#ifdef _STREFMULTCORR_
bool CutClass::isBadStRefMultCorrRunId(int a_runID){
  bool badRunBool = m_refMultCorr->isBadRun(a_runID);
  if(badRunBool){
    m_rejectedEvents->Fill(5.5);
  }
  return badRunBool;
};

#endif

bool CutClass::isGoodBTof(double a_localY, double a_localZ){//, double a_btofBeta, double a_dEdx){
  #ifdef _CUTCLASS_DEBUG_

  cout << "  isGoodBTOF?" << endl;
  cout << "       localY = " << a_localY << endl;
  cout << "       localZ = " << a_localZ << endl;
  //cout << "       beta   = " << a_btofBeta << endl;
  //cout << "       dEdx   = " << a_dEdx << endl;


  #endif

  /*if(a_btofMatchFlag != 0){

    #ifdef _CUTCLASS_DEBUG_
      cout << "  rejected btof track because of btofMatchFlag" << endl;
    #endif
    return false;
  }*/

  ///if(a_btofMatchFlag == 0) return false;
  if(fabs(a_localY) > m_btofLocY){
    m_rejectedTracks->Fill(10.5);
    #ifdef _CUTCLASS_DEBUG_
      cout << "  rejected btof track because of local y : " << fabs(a_localY) << " > " << m_btofLocY << endl;
    #endif
    return false;
  }
  if(fabs(a_localZ) > m_btofLocZ){
    m_rejectedTracks->Fill(11.5);
    #ifdef _CUTCLASS_DEBUG_
      cout << "  rejected btof track because of local z : " << fabs(a_localZ) << " > " << m_btofLocZ << endl;
    #endif
    return false;

  }

  /*if(a_beta <= 0){
    #ifdef _CUTCLASS_DEBUG_
      cout << "  rejected btof track because of bTOF Beta : " << a_beta << " <=0" << endl;
    #endif
    return false;

  }*/

  /*  //MuDst embedding files do not have dE/dx simulator
  //cout << " In isGoodBTof" << endl;
  //this cuts out double pion and double kaon peaks
  if(m_useBTOFBetaGammaCuts){
    if(!isGoodTrack_GivenBTOFBetaGamma(a_dEdx,PhysMath::beta_to_betagamma(a_btofBeta)) ){
      #ifdef _CUTCLASS_DEBUG_
        cout << "  rejected btof track because of betaGamma DeDx : " << PhysMath::beta_to_betagamma(a_btofBeta)  << "  and  " << a_dEdx << endl;
      #endif
      return false;

    }

  }*/
  return true;
}


bool CutClass::isGoodTrack_GivenBTOFBetaGamma(double a_dEdx, double a_betaGamma){

  for(int iii = 0; iii < (int) m_BTOFBetaGammaLowLineIntercepts.size(); iii++){
    if(log(a_dEdx) < m_BTOFBetaGammaLowLineSlopes[iii]*log(a_betaGamma) + m_BTOFBetaGammaLowLineIntercepts[iii]){
     m_rejectedTracks->Fill(13.5);
     return false;
    }
  }
  for(int iii = 0; iii < (int) m_BTOFBetaGammaHighLineYPoints.size(); iii+=2){
    if(log(a_betaGamma) > m_BTOFBetaGammaHighLineXPoints[iii]
        && log(a_betaGamma) < m_BTOFBetaGammaHighLineXPoints[iii + 1]){
      double slope = (m_BTOFBetaGammaHighLineYPoints[iii] - m_BTOFBetaGammaHighLineYPoints[iii+1])
	       /(m_BTOFBetaGammaHighLineXPoints[iii] - m_BTOFBetaGammaHighLineXPoints[iii+1]);
      double intercept = m_BTOFBetaGammaHighLineYPoints[iii] - slope * m_BTOFBetaGammaHighLineXPoints[iii];
      if(log(a_dEdx) > slope*log(a_betaGamma)+intercept){
        m_rejectedTracks->Fill(13.5);
        return false;
      }
    }
  }

  return true;
}

bool CutClass::isGoodTrack_GivenBTOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr){
  if(!m_useBTOFInvBetaGammaSqrCuts) return true;
  if(a_invBetaGammaSqr < m_lowInvBetaGammaSqrBTOFCut){
    m_rejectedTracks->Fill(13.5);
    return false;
  }
  for(int iii = 0; iii < (int) m_BTOFInvBetaGammaSqrLowLineIntercepts.size(); iii++){
    if(a_dEdx < m_BTOFInvBetaGammaSqrLowLineSlopes[iii]*a_invBetaGammaSqr + m_BTOFInvBetaGammaSqrLowLineIntercepts[iii]){
     m_rejectedTracks->Fill(13.5);
     return false;
    }
  }
  for(int iii = 0; iii < (int) m_BTOFInvBetaGammaSqrHighLineYPoints.size(); iii+=2){
    if(a_invBetaGammaSqr > m_BTOFInvBetaGammaSqrHighLineXPoints[iii]
        && a_invBetaGammaSqr < m_BTOFInvBetaGammaSqrHighLineXPoints[iii + 1]){
      double slope = (m_BTOFInvBetaGammaSqrHighLineYPoints[iii] - m_BTOFInvBetaGammaSqrHighLineYPoints[iii+1])
               /(m_BTOFInvBetaGammaSqrHighLineXPoints[iii] - m_BTOFInvBetaGammaSqrHighLineXPoints[iii+1]);
      double intercept = m_BTOFInvBetaGammaSqrHighLineYPoints[iii] - slope * m_BTOFInvBetaGammaSqrHighLineXPoints[iii];
      if(a_dEdx > slope*a_invBetaGammaSqr+intercept){
        m_rejectedTracks->Fill(13.5);
        return false;
      }
    }
  }
  if(m_useModNSigma_InvBetaGammaSqrBTOFCuts && a_invBetaGammaSqr > m_lowInvBetaGammaSqr_ToUseModNSigmaBTOFCut){
    double modNSigma = getModNSigma_bTOF(a_dEdx,a_invBetaGammaSqr);
    if(modNSigma < m_lowModNSigma_BTOF_Cut || modNSigma > m_highModNSigma_BTOF_Cut){
      m_rejectedTracks->Fill(13.5);
      return false;
    }
  }
  return true;
}

bool CutClass::isElectron_GivenBTOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr){
  if(m_useBTOFElectronCut){
    if(a_invBetaGammaSqr < m_electron_lowInvBetaGammaSqrBTOFCut || a_invBetaGammaSqr > m_electron_highInvBetaGammaSqrBTOFCut
       || a_dEdx < m_electron_lowDeDxBTOFCut || a_dEdx > m_electron_highDeDxBTOFCut){
      return false;
    }
    return true;
  }
  return false;
}

/*
bool CutClass::isGoodTrack_GivenETOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr){
  if(a_invBetaGammaSqr < m_lowInvBetaGammaSqrETOFCut){
    m_rejectedTracks->Fill(17.5);
    return false;
  }
  for(int iii = 0; iii < (int) m_ETOFInvBetaGammaSqrLowLineIntercepts.size(); iii++){
    if(a_dEdx < m_ETOFInvBetaGammaSqrLowLineSlopes[iii]*a_invBetaGammaSqr + m_ETOFInvBetaGammaSqrLowLineIntercepts[iii]){
     m_rejectedTracks->Fill(17.5);
     return false;
    }
  }
  for(int iii = 0; iii < (int) m_ETOFInvBetaGammaSqrHighLineYPoints.size(); iii+=2){
    if(a_invBetaGammaSqr > m_ETOFInvBetaGammaSqrHighLineXPoints[iii]
        && a_invBetaGammaSqr < m_ETOFInvBetaGammaSqrHighLineXPoints[iii + 1]){
      double slope = (m_ETOFInvBetaGammaSqrHighLineYPoints[iii] - m_ETOFInvBetaGammaSqrHighLineYPoints[iii+1])
               /(m_ETOFInvBetaGammaSqrHighLineXPoints[iii] - m_ETOFInvBetaGammaSqrHighLineXPoints[iii+1]);
      double intercept = m_ETOFInvBetaGammaSqrHighLineYPoints[iii] - slope * m_ETOFInvBetaGammaSqrHighLineXPoints[iii];
      if(a_dEdx > slope*a_invBetaGammaSqr+intercept){
        m_rejectedTracks->Fill(17.5);
        return false;
      }
    }
  }
  if(m_useModNSigma_InvBetaGammaSqrETOFCuts && a_invBetaGammaSqr > m_lowInvBetaGammaSqr_ToUseModNSigmaETOFCut){
    double modNSigma = getModNSigma_eTOF(a_dEdx,a_invBetaGammaSqr);
    if(modNSigma < m_lowModNSigma_ETOF_Cut || modNSigma > m_highModNSigma_ETOF_Cut){
      m_rejectedTracks->Fill(17.5);
      return false;
    }
  }
  return true;
}*/


bool CutClass::isGoodTrack_GivenETOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr){
  if(!m_useETOFInvBetaGammaSqrCuts) return true;
  if(a_invBetaGammaSqr < m_lowInvBetaGammaSqrETOFCut){
    m_rejectedTracks->Fill(17.5);
    return false;
  }
  for(int iii = 0; iii < (int) m_ETOFInvBetaGammaSqrLowLineIntercepts.size(); iii++){
    if(a_dEdx < m_ETOFInvBetaGammaSqrLowLineSlopes[iii]*a_invBetaGammaSqr + m_ETOFInvBetaGammaSqrLowLineIntercepts[iii]){
     m_rejectedTracks->Fill(17.5);
     return false;
    }
  }
  for(int iii = 0; iii < (int) m_ETOFInvBetaGammaSqrHighLineYPoints.size(); iii+=2){
    if(a_invBetaGammaSqr > m_ETOFInvBetaGammaSqrHighLineXPoints[iii]
        && a_invBetaGammaSqr < m_ETOFInvBetaGammaSqrHighLineXPoints[iii + 1]){
      double slope = (m_ETOFInvBetaGammaSqrHighLineYPoints[iii] - m_ETOFInvBetaGammaSqrHighLineYPoints[iii+1])
               /(m_ETOFInvBetaGammaSqrHighLineXPoints[iii] - m_ETOFInvBetaGammaSqrHighLineXPoints[iii+1]);
      double intercept = m_ETOFInvBetaGammaSqrHighLineYPoints[iii] - slope * m_ETOFInvBetaGammaSqrHighLineXPoints[iii];
      if(a_dEdx > slope*a_invBetaGammaSqr+intercept){
        m_rejectedTracks->Fill(17.5);
        return false;
      }
    }
  }
  if(m_useModNSigma_InvBetaGammaSqrETOFCuts && a_invBetaGammaSqr > m_lowInvBetaGammaSqr_ToUseModNSigmaETOFCut){
    double modNSigma = getModNSigma_bTOF(a_dEdx,a_invBetaGammaSqr);
    if(modNSigma < m_lowModNSigma_ETOF_Cut || modNSigma > m_highModNSigma_ETOF_Cut){
      m_rejectedTracks->Fill(17.5);
      return false;
    }
  }
  return true;
}


bool CutClass::isElectron_GivenETOFInvBetaGammaSqr(double a_dEdx, double a_invBetaGammaSqr){
  if(m_useETOFElectronCut){
    if(a_invBetaGammaSqr < m_electron_lowInvBetaGammaSqrETOFCut || a_invBetaGammaSqr > m_electron_highInvBetaGammaSqrETOFCut
       || a_dEdx < m_electron_lowDeDxETOFCut || a_dEdx > m_electron_highDeDxETOFCut){
      return false;
    }
    return true;
  }
  return false;
}




double CutClass::getModNSigma_bTOF(double a_dEdx, double a_invBetaGammaSqr){
  return log(a_dEdx/m_invBetaGammaSqr_dEdx_Calib_funct->Eval(a_invBetaGammaSqr))/log(1.07);
}
double CutClass::getModNSigma_eTOF(double a_dEdx, double a_invBetaGammaSqr){
  return log(a_dEdx/m_invBetaGammaSqr_dEdx_Calib_funct->Eval(a_invBetaGammaSqr))/log(1.07);
}

double CutClass::getDeDxPrediction(double a_p, double a_m){
  return m_invBetaGammaSqr_dEdx_Calib_funct->Eval(pow(a_m/a_p,2));
}
double CutClass::getDeDxPrediciton_byPart(double a_p, int a_partIndex){
  if(a_partIndex != 3) return m_invBetaGammaSqr_dEdx_Calib_funct->Eval(pow(m_particleMasses[a_partIndex]/a_p,2));
  return 3.98508;
}



TH1D* CutClass::getRapMtM0BinStructure(int a_particleIndex, int a_tag){
  switch(a_tag){
    case 0: return m_mTm0_structs[a_particleIndex];
    case 1: return m_tpc_rap_structs[a_particleIndex];
    case 2: return m_btof_rap_structs[a_particleIndex];
    case 3: return m_etof_rap_structs[a_particleIndex];
    case 4: return m_btof_mTm0_structs[a_particleIndex];
    case 5: return m_etof_mTm0_structs[a_particleIndex];
    default : cout << "you did not select a proper tag in CutClass::getBinStructure... 0=mtm0, rap-> 1=tpc, 2=btof, 3=etof 4=mTm0BTOF 5=mTm0ETOF" << endl; return NULL;
  }

  return NULL;
}



int CutClass::getNRapidityBins(int a_particleIndex, int a_tag){
  TH1D* rapStructHist = getRapMtM0BinStructure(a_particleIndex,a_tag);
  if(rapStructHist){
    return rapStructHist->GetNbinsX();
  }else{
    return 0;
  }
  //return getRapMtM0BinStructure(a_particleIndex,a_tag)->GetNbinsX();
}

double CutClass::getLowRapidity(int a_particleIndex, int a_tag){
  //TH1D *rapStructHist = m_tpc_rap_structs[a_particleIndex];
  //return rapStructHist->GetBinLowEdge(1);
  return getRapMtM0BinStructure(a_particleIndex,a_tag)->GetBinLowEdge(1);
}

double CutClass::getHighRapidity(int a_particleIndex, int a_tag){
  //TH1D *rapStructHist = m_tpc_rap_structs[a_particleIndex];
  //return rapStructHist->GetBinLowEdge(rapStructHist->GetNbinsX()+1);
  TH1D *rapStructHist = getRapMtM0BinStructure(a_particleIndex,a_tag);
  return rapStructHist->GetBinLowEdge(rapStructHist->GetNbinsX()+1);
}

double CutClass::getRapidityRangeHigh(int a_yIndex, int a_partIndex, int a_tag){
  return getRapMtM0BinStructure(a_partIndex,a_tag)->GetBinLowEdge(a_yIndex+2);
}


double CutClass::getRapidityRangeLow(int a_yIndex, int a_partIndex, int a_tag){
  return getRapMtM0BinStructure(a_partIndex,a_tag)->GetBinLowEdge(a_yIndex+1);
}

double CutClass::getRapidityRangeCenter(int a_yIndex, int a_partIndex, int a_tag){
  return getRapMtM0BinStructure(a_partIndex,a_tag)->GetBinCenter(a_yIndex+1);
}

int CutClass::getN_mTm0Bins(int a_particleIndex){
  TH1D *mTm0StructHist = m_tpc_rap_structs[a_particleIndex];

  return mTm0StructHist->GetNbinsY();
}

double CutClass::getLow_mTm0(int a_particleIndex){
  TH1D *mTm0StructHist = m_tpc_rap_structs[a_particleIndex];

  return mTm0StructHist->GetYaxis()->GetBinLowEdge(1);
}

double CutClass::getHigh_mTm0(int a_particleIndex){
  TH1D *mTm0StructHist = m_tpc_rap_structs[a_particleIndex];

  return mTm0StructHist->GetYaxis()->GetBinLowEdge(mTm0StructHist->GetNbinsY()+1);
}






#ifndef _MAC_OSX_  // no StPicoDst stuff --------------------------------------------------------------------------------------------------------------------------------------------------------------


bool CutClass::isHLTEvent(StPicoEvent* a_event){
  bool HLTfoundInTriggers = false;
  if(m_cutOnTriggers){
    vector <unsigned int> trigIDs = a_event->triggerIds();
    for(int reqTrigIndex = 0; reqTrigIndex < (int)m_HLTTriggers.size(); reqTrigIndex++){
      for(int eventTrigIndex = 0; eventTrigIndex < (int)trigIDs.size(); eventTrigIndex++){
	//    cout << trigIDs[eventTrigIndex] << " to " << m_triggers[reqTrigIndex] << endl;
        if(trigIDs[eventTrigIndex] == m_HLTTriggers[reqTrigIndex]) HLTfoundInTriggers = true;
        if(HLTfoundInTriggers) break;
      }
      if(HLTfoundInTriggers) break;
    }
  }
  return HLTfoundInTriggers;
}

#ifdef _HAS_ETOF_
bool CutClass::isEToFEvent(StPicoEvent* a_event){
  bool EToFfoundInTriggers = false;
  if(m_cutOnTriggers){
    vector <unsigned int> trigIDs = a_event->triggerIds();
    for(int reqTrigIndex = 0; reqTrigIndex < (int)m_EToFTriggers.size(); reqTrigIndex++){
      for(int eventTrigIndex = 0; eventTrigIndex < (int)trigIDs.size(); eventTrigIndex++){
	//    cout << trigIDs[eventTrigIndex] << " to " << m_triggers[reqTrigIndex] << endl;
        if(trigIDs[eventTrigIndex] == m_EToFTriggers[reqTrigIndex]) EToFfoundInTriggers = true;
        if(EToFfoundInTriggers) break;
      }
      if(EToFfoundInTriggers) break;
    }
  }
  return EToFfoundInTriggers;
}
#endif

bool CutClass::isGoodEvent(StPicoEvent* a_event){
  //Vertex
  double z =(double) a_event->primaryVertex().Z();
  double r = TMath::Sqrt( pow(a_event->primaryVertex().X()-m_beamXCenter,2) + pow(a_event->primaryVertex().Y()-m_beamYCenter,2) );

  #ifdef _CUTCLASS_DEBUG_
    cout << "isGoodEvent(): zVertex " << z << "=?["<< m_ZLowCut <<","<< m_ZHighCut <<"]" << endl;
    cout << "               rVertex " << r << " > " << m_radialCut << endl;
    cout << "               tofMatch " << a_event->nBTOFMatch() << " < " << m_tofMatchCut << endl;
    if(m_useT0){
      cout << "               nTofT0 " <<  a_event->nTofT0() << " < " <<  m_nTofT0Cut << endl;
    }
    cout << "               triggers:" << endl;
    for(int eventTrigIndex = 0; eventTrigIndex < (int)a_event->triggerIds().size(); eventTrigIndex++){
     cout << "     " << a_event->triggerIds()[eventTrigIndex];
   }
   cout << "\n               in Req Trigs: " << endl;
   for(int reqTrigIndex = 0; reqTrigIndex < (int)m_triggers.size(); reqTrigIndex++){
    cout << "    " << m_triggers[reqTrigIndex];
  }
  cout << endl;
  #endif

  if (z < m_ZLowCut){
     m_rejectedEvents->Fill(0.5);
     #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of z Vertex = " << z << " < " << m_ZLowCut << endl;
     #endif

     return false;
  }
  if (z > m_ZHighCut){
     m_rejectedEvents->Fill(0.5);
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of z Vertex = " << z << " > " << m_ZHighCut << endl;
    #endif
     return false;
  }
  if (r > m_radialCut){
     m_rejectedEvents->Fill(1.5);
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of radius  = " << r << " > " << m_radialCut << endl;
    #endif
     return false;
  }
  if (m_tofMatchCut > a_event->nBTOFMatch()){
     m_rejectedEvents->Fill(2.5);
     #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of tofMatch = " << a_event->nBTOFMatch() << " < " << m_tofMatchCut << endl;
     #endif

     return false;
  }
  if(m_useT0){
    if (a_event->nTofT0() < m_nTofT0Cut){
       m_rejectedEvents->Fill(3.5);
      #ifdef _CUTCLASS_DEBUG_
       cout << " Rejected because of nTOFT0= " << a_event->nTofT0() << " < " << m_nTofT0Cut << endl;
      #endif
       return false;
    }
  }
    // Trigger
  bool goodEvent = true;
  if(m_cutOnTriggers){
    bool foundInTriggers = false;
    vector <unsigned int> trigIDs = a_event->triggerIds();
    for(int reqTrigIndex = 0; reqTrigIndex < (int)m_triggers.size(); reqTrigIndex++){
        //foundInTriggers = false;
      for(int eventTrigIndex = 0; eventTrigIndex < (int)trigIDs.size(); eventTrigIndex++){
          //    cout << trigIDs[eventTrigIndex] << " to " << m_triggers[reqTrigIndex] << endl;
        if(trigIDs[eventTrigIndex] == m_triggers[reqTrigIndex]) foundInTriggers = true;
        if(foundInTriggers) break;
      }
      if(foundInTriggers) break;
    }

    if(!foundInTriggers) goodEvent = false;
  }
  #ifdef _CUTCLASS_DEBUG_
  if(!goodEvent) {
    #ifdef _CUTCLASS_DEBUG_
     cout << " Rejected because of Triggers " << endl;
    #endif
    m_rejectedEvents->Fill(4.5);

  }else{
    #ifdef _CUTCLASS_DEBUG_
      cout << "Passed Triggers" << endl;
    #endif
  }
  #endif
  m_rejectedEvents->Fill(8.5);
  return goodEvent;
}









bool CutClass::isGoodTrack(StPicoTrack* a_track, TVector3* a_EventVertex, bool a_ignoreDCA){
  #ifdef _CUTCLASS_DEBUG_
    cout << "    isGoodTrack(): Primary:  " << a_track->isPrimary()  << endl;
    cout << "               nHitsFit: " << a_track->nHitsFit() << "   nHitsMax: " << a_track->nHitsMax() << endl;
    cout << "               nHitsFit/nHitsMax : " << ((double) (a_track->nHitsFit())) / ((double)(a_track->nHitsMax())) << " < " << m_fitMaxRatio << endl;
    cout << "               nHitsDedx : " << a_track->nHitsDedx() << " < " << m_NHitsDeDx << endl;
    cout << "               nHitsFit : " << a_track->nHitsFit() << " < " << m_NHitsFit << endl;
    cout << "               gDCA : " << a_track->gDCA(*a_EventVertex).Mag() << " >? " << m_DCA << endl;
    if(m_usePhiCut) cout << "               phi : " << TMath::ATan2(a_track->pMom().Y(), a_track->pMom().X()) << " =! [" << m_lowPhiCut << ", " << m_highPhiCut << "]" << endl;
    if(m_useEtaCut) cout << "               eta : " << a_track->pMom().PseudoRapidity() << " != [" << m_lowEtaCut <<", " << m_highEtaCut << "]" << endl;
    if(m_usePtCut) cout << "               pT : " << a_track->pPt() << " < " << m_lowPtCut << endl;
  #endif

  if(!a_track->isPrimary()) {
      m_rejectedTracks->Fill(0.5);
      return false;
  }
  // cout << "Track is primary" << endl;
  if(  ((double) (a_track->nHitsFit())) / ((double)(a_track->nHitsMax())) < m_fitMaxRatio ) {
      m_rejectedTracks->Fill(1.5);
      return false;
  }
  //  cout << "Track has good enough hit fraction: " << ((double) (a_track->nHitsFit())) / ((double)(a_track->nHitsMax())) << " > " << m_fitMaxRatio << endl;
  if(a_track->nHitsFit() < m_NHitsFit) {
      m_rejectedTracks->Fill(2.5);
      return false;
  }
  //  cout << "Track has enough hits: " << a_track->nHitsFit() << " > " << m_NHitsFit << endl;
  if(a_track->nHitsDedx() < m_NHitsDeDx) {
      m_rejectedTracks->Fill(3.5);
      return false;
  }
  //  cout << "Track has enough dEdx hits: " << a_track->nHitsDedx() << " > " << m_NHitsDeDx << endl;
  if(!a_ignoreDCA && a_track->gDCA(*a_EventVertex).Mag() > m_DCA) {
      m_rejectedTracks->Fill(4.5);
      return false;
  }
  //  cout << "Track has good enough DCA: " << a_track->gDCA(*a_EventVertex).Mag() << " < " << m_DCA << endl;
  //  cout << "Track is good" << endl;
  if(m_usePhiCut){
    double phi = TMath::ATan2(a_track->pMom().Y(), a_track->pMom().X());
    if( m_lowPhiCut < phi && phi < m_highPhiCut ) {
        m_rejectedTracks->Fill(5.5);
        return false;
    }
  }

  if(m_useEtaCut){
    double eta = a_track->pMom().PseudoRapidity();
    if( m_lowEtaCut > eta || m_highEtaCut < eta){
      m_rejectedTracks->Fill(6.5);
      return false;
    }
  }

  if(m_usePtCut){
    if( a_track->pPt() < m_lowPtCut){
      m_rejectedTracks->Fill(7.5);
      return false;
    }
  }

  return true;

}






/*
bool CutClass::isGoodBTof(StPicoBTofPidTraits* a_traits, StPicoTrack* a_track){
  if(!a_traits) return false;
  if(a_traits->btofMatchFlag() == 0) return false;
  if(fabs(a_traits->btofYLocal()) > m_btofLocY) return false;
  if(fabs(a_traits->btofZLocal()) > m_btofLocZ) return false;
  if(a_traits->btofBeta() < 0) return false;

  //cout << " In isGoodBTof" << endl;
  //this cuts out double pion and double kaon peaks
  if(a_track && m_useBTOFBetaGammaCuts){
    double betaGamma = PhysMath::beta_to_betagamma(a_traits->btofBeta());
    double dEdx = a_track->dEdx();
    if(!isGoodTrack_GivenBTOFBetaGamma(a_track->dEdx(),PhysMath::beta_to_betagamma(a_traits->btofBeta())) ) return false;
  }

  return true;
}

bool CutClass::isGoodTrack_GivenBTOFBetaGamma(double a_dEdx, double a_betaGamma){

  for(int iii = 0; iii < (int) m_BTOFBetaGammaLowLineIntercepts.size(); iii++){
    if(log(a_dEdx) < m_BTOFBetaGammaLowLineSlopes[iii]*log(a_betaGamma) + m_BTOFBetaGammaLowLineIntercepts[iii]) return false;
  }
  for(int iii = 0; iii < (int) m_BTOFBetaGammaHighLineYPoints.size(); iii+=2){
    if(log(a_betaGamma) > m_BTOFBetaGammaHighLineXPoints[iii] && log(a_betaGamma) < m_BTOFBetaGammaHighLineXPoints[iii + 1]){
      double slope = (m_BTOFBetaGammaHighLineYPoints[iii] - m_BTOFBetaGammaHighLineYPoints[iii+1])
                        /(m_BTOFBetaGammaHighLineXPoints[iii] - m_BTOFBetaGammaHighLineXPoints[iii+1]);
      double intercept = m_BTOFBetaGammaHighLineYPoints[iii] - slope * m_BTOFBetaGammaHighLineXPoints[iii];
      if(log(a_dEdx) > slope*log(a_betaGamma)+intercept) return false;
    }
  }

  return true;
}
*/

#endif // mac osx not defined --------------------------------------------------------------------------------------------------------------------------------------------------------------

// Pulled out from the _MAC_OSX_ gate above -- see the matching comment in
// headers/CutClass.h. This function's body only touches plain member variables
// (m_etof_match_flag_cut, m_rejectedTracks, m_etof_max_cluster_size,
// m_etofDetlaXY_Center_X/Y, m_etofDeltaX_Ellipse_Cut/Y) plus the
// StPicoETofPidTraits/StPicoETofHit types it's declared with -- nothing here
// actually depends on _MAC_OSX_ being undefined. Logic is otherwise byte-for-byte
// identical to the original.
#ifdef _HAS_ETOF_
bool CutClass::isGoodETof(StPicoETofPidTraits* a_traits, StPicoETofHit* a_hit){
  if(m_etof_match_flag_cut < 0){
    if(a_traits->matchFlag() == 0){
      m_rejectedTracks->Fill(15.5);
      return false;
    }
  }else{
    if(a_traits->matchFlag() != m_etof_match_flag_cut){
      m_rejectedTracks->Fill(15.5);
      return false;
    }
  }
  if( a_hit->clusterSize() > m_etof_max_cluster_size){
    m_rejectedTracks->Fill(16.5);
    return false;
  }
  if(pow((a_traits->deltaX()-m_etofDetlaXY_Center_X)/m_etofDeltaX_Ellipse_Cut,2) + pow((a_traits->deltaY() - m_etofDeltaXY_Center_Y)/m_etofDeltaY_Ellipse_Cut,2) > 1.0){
    m_rejectedTracks->Fill(14.5);
    return false;
  }
  //if(a_traits->beta() > 1.5){
  //  m_rejectedTracks->Fill(18.5);
  //  return false;
  //}

  //double invBGSqr = PhysMath::invBetaGammaSqr_signed_fromInvBeta(1.0/a_traits->beta());
  //if(invBGSqr < m_lowInvBetaGammaSqrETOFCut) return false;
  //if(m_useETOFInvBetaGammaSqrCuts && !isGoodTrack_GivenETOFInvBetaGammaSqr(a_track->dEdx(),invBGSqr)){
  //    return false;
  //}

  return true;

}
#endif






bool CutClass::isGoodEvent(double a_x, double a_y, double a_z, int a_nBTOFMatch, int a_nTofT0){
  //Vertex
  double r = TMath::Sqrt( pow(a_x-m_beamXCenter,2) + pow(a_y-m_beamYCenter,2) );

  #ifdef _CUTCLASS_DEBUG_
    cout << "isGoodEvent(): zVertex " << a_z << "=?["<< m_ZLowCut <<","<< m_ZHighCut <<"]" << endl;
    cout << "               rVertex " << r << " > " << m_radialCut << endl;
    cout << "               tofMatch " << a_nBTOFMatch << " < " << m_tofMatchCut << endl;
    if(m_useT0){
      cout << "               nTofT0 " <<  a_nTofT0 << " < " <<  m_nTofT0Cut << endl;
    }
  #endif


  if (a_z < m_ZLowCut){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of z Vertex = " << a_z << " < " << m_ZLowCut << endl;
    #endif
    m_rejectedEvents->Fill(0.5);
    return false;
  }
  if (a_z > m_ZHighCut){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of z Vertex = " << a_z << " > " << m_ZHighCut << endl;
    #endif
    m_rejectedEvents->Fill(0.5);
    return false;
  }
  if (r > m_radialCut){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of radius  = " << r << " > " << m_radialCut << endl;
    #endif
    m_rejectedEvents->Fill(1.5);
    return false;
  }
  if (m_tofMatchCut > a_nBTOFMatch){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of tofMatch = " << a_nBTOFMatch << " < " << m_tofMatchCut << endl;
    #endif
    m_rejectedEvents->Fill(2.5);
    return false;
  }
  if(m_useT0){
    if(a_nTofT0 < m_nTofT0Cut){
      #ifdef _CUTCLASS_DEBUG_
        cout << " Rejected because of nTOFT0= " << a_nTofT0 << " < " << m_nTofT0Cut << endl;
      #endif
      m_rejectedEvents->Fill(3.5);
      return false;
    }
  }
  return true;

}



bool CutClass::isGoodTrack(bool a_isPrimary,int a_nHitsFit, int a_nHitsMax, int a_nHitsDedx, double a_gDCA, double a_phi){

  #ifdef _CUTCLASS_DEBUG_
    cout << "    isGoodTrack(): Primary:  " << a_isPrimary  << endl;
    cout << "               nHitsFit: " << a_nHitsFit << "   nHitsMax: " << a_nHitsMax << endl;
    cout << "               nHitsFit/nHitsMax : " << ((double) (a_nHitsFit)) / ((double)(a_nHitsMax)) << " < " << m_fitMaxRatio << endl;
    cout << "               nHitsDedx : " << a_nHitsDedx << " < " <<  m_NHitsDeDx << endl;
    cout << "               nHitsFit : " << a_nHitsFit << " < " << m_NHitsFit << endl;
    cout << "               gDCA : " << a_gDCA << " >? " << m_DCA << endl;
    if(m_usePhiCut) cout << "               phi : " << a_phi << " =! [" << m_lowPhiCut << ", " << m_highPhiCut << "]" << endl;
  #endif


  if(!a_isPrimary){
     #ifdef _CUTCLASS_DEBUG_
       cout << " Rejected because not primary " << endl;
     #endif
     m_rejectedTracks->Fill(0.5);
     return false;
  }
  // cout << "Track is primary" << endl;
  if(  ((double) (a_nHitsFit)) / ((double)(a_nHitsMax)) < m_fitMaxRatio ){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because fit/max ratio hits " << endl;
    #endif
   m_rejectedTracks->Fill(1.5);
   return false;
  }
  //  cout << "Track has good enough hit fraction: " << ((double) (a_track->nHitsFit())) / ((double)(a_track->nHitsMax())) << " > " << m_fitMaxRatio << endl;
  if(a_nHitsFit < m_NHitsFit){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because hits fit " << a_nHitsFit << " < " << m_NHitsFit << endl;
    #endif
    m_rejectedTracks->Fill(2.5);
    return false;
  }
  //  cout << "Track has enough hits: " << a_track->nHitsFit() << " > " << m_NHitsFit << endl;
  if(a_nHitsDedx < m_NHitsDeDx){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because hits fit dEdx" << a_nHitsDedx << " < " << m_NHitsDeDx << endl;
    #endif
    m_rejectedTracks->Fill(3.5);
    return false;
  }
  //  cout << "Track has enough dEdx hits: " << a_track->nHitsDedx() << " > " << m_NHitsDeDx << endl;
  if(a_gDCA > m_DCA){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of DCA " << a_gDCA << " < " << m_DCA << endl;
    #endif
    m_rejectedTracks->Fill(4.5);
    return false;
  }
  //  cout << "Track has good enough DCA: " << a_track->gDCA(*a_EventVertex).Mag() << " < " << m_DCA << endl;
  //  cout << "Track is good" << endl;
  if(m_usePhiCut){
    if( m_lowPhiCut < a_phi && a_phi < m_highPhiCut ){
     #ifdef _CUTCLASS_DEBUG_
       cout << " Rejected because phi cut  " << a_phi << " =? [" << m_lowPhiCut << ", " << m_highPhiCut << "]" << endl;
     #endif
     m_rejectedTracks->Fill(5.5);
     return false;
    }
  }

  #ifdef _CUTCLASS_DEBUG_
    cout << " GOOD TRACK from CutClass_Debug " << endl;
  #endif

  return true;
}




bool CutClass::isGoodTrack(int a_nHitsFit, double a_nHitsOverNPoss, int a_nHitsDedx, double a_gDCA, double a_phi){

  #ifdef _CUTCLASS_DEBUG_
    cout << "    isGoodTrack():  "  << endl;
    cout << "               nHitsFit: " << a_nHitsFit << "   nHitsRatio: " << a_nHitsOverNPoss << endl;
    cout << "               nHitsFit/nHitsMax : " << ((double) (a_nHitsFit)) / ((double)(a_nHitsMax)) << " < " << m_fitMaxRatio << endl;
    cout << "               nHitsDedx : " << a_nHitsDedx << " < " <<  m_NHitsDeDx << endl;
    cout << "               nHitsFit : " << a_nHitsFit << " < " << m_NHitsFit << endl;
    cout << "               gDCA : " << a_gDCA << " >? " << m_DCA << endl;
    if(m_usePhiCut) cout << "               phi : " << a_phi << " =! [" << m_lowPhiCut << ", " << m_highPhiCut << "]" << endl;
  #endif

  if(  a_nHitsOverNPoss < m_fitMaxRatio ){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because fit/max ratio hits " << endl;
    #endif
   m_rejectedTracks->Fill(1.5);
   return false;
  }
  //  cout << "Track has good enough hit fraction: " << ((double) (a_track->nHitsFit())) / ((double)(a_track->nHitsMax())) << " > " << m_fitMaxRatio << endl;
  if(a_nHitsFit < m_NHitsFit){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because hits fit " << a_nHitsFit << " < " << m_NHitsFit << endl;
    #endif
    m_rejectedTracks->Fill(2.5);
    return false;
  }
  //  cout << "Track has enough hits: " << a_track->nHitsFit() << " > " << m_NHitsFit << endl;
  if(a_nHitsDedx < m_NHitsDeDx){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because hits fit dEdx" << a_nHitsDedx << " < " << m_NHitsDeDx << endl;
    #endif
    m_rejectedTracks->Fill(3.5);
    return false;
  }
  //  cout << "Track has enough dEdx hits: " << a_track->nHitsDedx() << " > " << m_NHitsDeDx << endl;
  if(a_gDCA > m_DCA){
    #ifdef _CUTCLASS_DEBUG_
      cout << " Rejected because of DCA " << a_gDCA << " < " << m_DCA << endl;
    #endif
    m_rejectedTracks->Fill(4.5);
    return false;
  }
  //  cout << "Track has good enough DCA: " << a_track->gDCA(*a_EventVertex).Mag() << " < " << m_DCA << endl;
  //  cout << "Track is good" << endl;
  if(m_usePhiCut){
    if( m_lowPhiCut < a_phi && a_phi < m_highPhiCut ){
     #ifdef _CUTCLASS_DEBUG_
       cout << " Rejected because phi cut  " << a_phi << " =? [" << m_lowPhiCut << ", " << m_highPhiCut << "]" << endl;
     #endif
     m_rejectedTracks->Fill(5.5);
     return false;
    }
  }

  #ifdef _CUTCLASS_DEBUG_
    cout << " GOOD TRACK from CutClass_Debug " << endl;
  #endif

  return true;
  }




bool CutClass::isSpecificParticleByTOF(double a_tofMass, int a_partIndex){
  if(a_partIndex >= 0 && a_partIndex < (int) m_highMassCut.size()){
    if( a_tofMass > m_lowMassCut[a_partIndex] && a_tofMass < m_highMassCut[a_partIndex]){
      return true;
    }
  }
  return false;
}

bool CutClass::isSpecificParticleByBTOF(double a_invBeta,double a_p, int a_partIndex){
  #ifdef _CutCLASS_DEBUG_
    cout << "CutClass::isSpecificParticleByBTOF(" << a_invBeta << ", " << a_p << ", " << a_partIndex << ")" << endl;
    cout << "    with internal mass: " << m_particleMasses_bTOF[a_partIndex] << "  and cut: " << m_invBetaCut_bTOF[a_partIndex] << endl;
  #endif
  return PhysMath::tofSelectionCut(a_p,m_particleMasses_bTOF[a_partIndex],a_invBeta,m_invBetaCut_bTOF[a_partIndex]);
}

bool CutClass::isSpecificParticleByETOF(double a_invBeta,double a_p, int a_partIndex){
  return PhysMath::tofSelectionCut(a_p,m_particleMasses_eTOF[a_partIndex],a_invBeta,m_invBetaCut_eTOF[a_partIndex]);
}




bool CutClass::isPileUp(int a_centVariable, int a_tofmatch){

  if(m_useHui_dt_27GeVPileup){
    return (a_tofmatch < 2.7*a_centVariable - 162.0 || a_tofmatch > 3.69*a_centVariable + 300);
  }

  if(m_useStRefMultCorr){
   #ifdef _STREFMULTCORR_
     #ifdef _CUTCLASS_DEBUG_
     cout << "Using StRefMultCorr as pileup, isPileup: " << !m_refMultCorr->passnTofMatchRefmultCut(a_centVariable, a_tofmatch) << endl;
     #endif
     bool isPileUpBool = !m_refMultCorr->passnTofMatchRefmultCut(a_centVariable, a_tofmatch);
     if(isPileUpBool){
       m_rejectedEvents->Fill(7.5);
     }
     return isPileUpBool;
   #else
     cout << "ERROR: St refmultcorr defined but not in  CutClass::isPileUp" << endl;
   #endif
    return true; //shouldn't get here!
  }else{
    double xVar = (double) a_centVariable;
    double yVar = (double) a_tofmatch;
    double top     =  (m_pileUpLineCuts[0]*xVar) + m_pileUpLineCuts[1];
    double bottom  =  (m_pileUpLineCuts[2]*xVar) + m_pileUpLineCuts[3];
    double bottom2 =  (m_pileUpLineCuts[4]*xVar) + m_pileUpLineCuts[5];
    if( top < yVar || bottom > yVar || bottom2 > yVar || (double) a_centVariable > m_pileUpCut){
     m_rejectedEvents->Fill(7.5);
     return true;
    }
    return false;
  }
}


int CutClass::centralityIndex(int a_mult){
  //if(isPileUp(a_mult)){ return -1;
  //}else{
  if(!m_useStRefMultCorr){
    for(int iii = 0; iii < m_numCentralities+1; iii++){
      if(a_mult >= m_centEdges[iii]) return iii-1;
    }
    m_rejectedEvents->Fill(6.5);
    return -1;
  }else{
    #ifdef _STREFMULTCORR_
      int corrIndex = 0;
      int centIndex = 0;
      if(m_refMultCorr9Bins){
        centIndex = 8 - m_refMultCorr->getCentralityBin9();
        if(centIndex < 0 || centIndex > 8){
          centIndex = -1;
          m_rejectedEvents->Fill(6.5);
        }
      }else{
        centIndex = 15 - m_refMultCorr->getCentralityBin16();
        if(centIndex < 0 || centIndex > 15){
          centIndex = -1;
          m_rejectedEvents->Fill(6.5);
        }
      }
      #ifdef _CUTCLASS_DEBUG_
        cout <<  "StRefMultCorr Found CentIndex: " << centIndex << endl;
      #endif
      return centIndex;
    #else
      return -1;
    #endif
  }
}


void CutClass::getParticleFitMatrix(int a_charge, bool a_particleFitMatrix[9][9]){
  if(a_charge>0){
    for(int partIndex=0; partIndex<9; partIndex++){
      for(int partBaseIndex=0; partBaseIndex<9; partBaseIndex++){
        a_particleFitMatrix[partBaseIndex][partIndex] = m_particleFitMatrix_Plus[partBaseIndex][partIndex];
      }
    }
  }
  else if(a_charge<0){
    for(int partBaseIndex=0; partBaseIndex<9; partBaseIndex++){
      for(int partIndex=0; partIndex<9; partIndex++){
        a_particleFitMatrix[partBaseIndex][partIndex] = m_particleFitMatrix_Minus[partBaseIndex][partIndex];
      }
    }
  }

  return;
}

double CutClass::getMomStops_TPC(int a_partIndex, int a_charge){
  //cout << "plus size: " << m_momStops_TPC_Plus.size() << "  minus size: " << m_momStops_TPC_Minus.size() << " trying part: " << a_partIndex << endl;
  if(a_charge > 0) return m_momStops_TPC_Plus.at(a_partIndex);
  else if(a_charge < 0) return m_momStops_TPC_Minus.at(a_partIndex);
  return -999;
}

double CutClass::getMomStops_bToF(int a_partIndex, int a_charge){
  //cout << "plus size: " << m_momStops_TPC_Plus.size() << "  minus size: " << m_momStops_TPC_Minus.size() << " trying part: " << a_partIndex << endl;
  if(a_charge > 0){
    if(a_partIndex >= m_momStops_bToF_Plus.size()) cout << " No momentum stopn plus in CutClass for partIndex: " << a_partIndex << endl;
    return m_momStops_bToF_Plus.at(a_partIndex);

  }
  else if(a_charge < 0){
    if(a_partIndex >= m_momStops_bToF_Minus.size()) cout << " No momentum stop minus in CutClass for partIndex: " << a_partIndex << endl;
    return m_momStops_bToF_Minus.at(a_partIndex);

  }
  return -999;
}


double CutClass::getMomStart_eToF(int a_partIndex, int a_charge){
  //cout << "plus size: " << m_momStops_TPC_Plus.size() << "  minus size: " << m_momStops_TPC_Minus.size() << " trying part: " << a_partIndex << endl;
  if(a_charge > 0) return m_momStarts_eToF_Plus.at(a_partIndex);
  else if(a_charge < 0) return m_momStarts_eToF_Minus.at(a_partIndex);
  return -999;
}

double CutClass::getMomStop_eToF(int a_partIndex, int a_charge){
  //cout << "plus size: " << m_momStops_TPC_Plus.size() << "  minus size: " << m_momStops_TPC_Minus.size() << " trying part: " << a_partIndex << endl;
  if(a_charge > 0) return m_momStops_eToF_Plus.at(a_partIndex);
  else if(a_charge < 0) return m_momStops_eToF_Minus.at(a_partIndex);
  return -999;
}

#endif
