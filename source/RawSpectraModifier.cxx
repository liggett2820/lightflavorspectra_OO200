#include "../headers/RawSpectraModifier.h"

RawSpectraModifier::RawSpectraModifier(){
  m_BTOF_Eff_Calc = false;
  m_ETOF_Eff_Calc = false;
  m_numParticles  = 0;

  m_chrisCorrectionFileMode = false;
  m_pTLossMode    = false;
  m_combined16To9 = false;
  m_equalCountsRebinForTPC = -1;

  for(int iii = 0; iii < 9; iii++){
    m_particleLoaded[iii]  = false;
    m_numEtaBins_TPC[iii]  = -1;
    m_numPtBins_TPC[iii]   = -1;
    m_numEtaBins_BTOF[iii] = -1;
    m_numPtBins_BTOF[iii]  = -1;
    m_numEtaBins_ETOF[iii] = -1;
    m_numPtBins_ETOF[iii]  = -1;
    m_tofIndexOffset[iii]  = -1;
    m_etofIndexOffset[iii] = -1; 
    m_centEvents[iii]      = NULL;
    m_rap_mTm0_dEdx_allCent[iii]        = NULL;
    m_rap_mTm0_invBetaBTOF_allCent[iii] = NULL;
    m_rap_mTm0_invBetaETOF_allCent[iii] = NULL;
    m_picoBinner_rapidity_structure[iii][0]  = NULL;
    m_picoBinner_rapidity_structure[iii][1]  = NULL;
    m_picoBinner_rapidity_structure[iii][2]  = NULL;
    m_embedding_rapidity_structure[iii]      = NULL;        
    m_hybridFeedDown_rapidity_structure[iii] = NULL;   


  }



  for(int iii = 0; iii < 9; iii++){
    for(int jjj = 0; jjj < 16; jjj++){

      m_BTOFEffPlusRatio[iii][jjj] = NULL;
      m_BTOFEffMinusRatio[iii][jjj] = NULL;
      m_BTOFEffPlus[iii][jjj] = NULL;
      m_BTOFEffMinus[iii][jjj] = NULL;
      m_ETOFEffPlus[iii][jjj] = NULL;
      m_ETOFEffMinus[iii][jjj] = NULL;


      m_dEdxSpectraPlus[iii][jjj] = NULL;
      m_InvBetaBTOFSpectraPlus[iii][jjj] = NULL;
      m_InvBetaETOFSpectraPlus[iii][jjj] = NULL;

      m_dEdxSpectraMinus[iii][jjj] = NULL;
      m_InvBetaBTOFSpectraMinus[iii][jjj] = NULL;
      m_InvBetaETOFSpectraMinus[iii][jjj] = NULL;

      m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[iii][jjj] = NULL;
      m_rap_mTm0_NumTPCTracks_ForETOF_Plus[iii][jjj] = NULL;
      m_rap_mTm0_NumBTOFTracks_Plus[iii][jjj] = NULL;
      m_rap_mTm0_NumETOFTracks_Plus[iii][jjj] = NULL;
      m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[iii][jjj] = NULL;
      m_rap_mTm0_NumTPCTracks_ForETOF_Minus[iii][jjj] = NULL;
      m_rap_mTm0_NumBTOFTracks_Minus[iii][jjj] = NULL;
      m_rap_mTm0_NumETOFTracks_Minus[iii][jjj] = NULL;

      m_btofScaleHistoPlus[iii][jjj] = NULL;
      m_btofScaleHistoMinus[iii][jjj] = NULL;

      for(int kkk = 0; kkk < g_RawSpecMod_MaxNumRapBins; kkk++){
        m_feeddownSpectraDeDx[iii][jjj][kkk][0] = NULL;
        m_feeddownSpectraInvBeta[iii][jjj][kkk][0] = NULL;
        m_feeddownSpectraInvBetaETOF[iii][jjj][kkk][0] = NULL;
        m_feeddownSpectraDeDx[iii][jjj][kkk][1] = NULL;
        m_feeddownSpectraInvBeta[iii][jjj][kkk][1] = NULL;
        m_feeddownSpectraInvBetaETOF[iii][jjj][kkk][1] = NULL;


        m_rawSpectraPlusDeDx[iii][jjj][kkk] = NULL;
        m_rawSpectraMinusDeDx[iii][jjj][kkk] = NULL;
        m_rawSpectraPlusInvBeta[iii][jjj][kkk] = NULL;
        m_rawSpectraMinusInvBeta[iii][jjj][kkk] = NULL;
        m_rawSpectraPlusInvBetaETOF[iii][jjj][kkk] = NULL;
        m_rawSpectraMinusInvBetaETOF[iii][jjj][kkk] = NULL;

        m_spectraPlusDeDx[iii][jjj][kkk] = NULL;
        m_spectraMinusDeDx[iii][jjj][kkk] = NULL;
        m_spectraPlusInvBeta[iii][jjj][kkk] = NULL;
        m_spectraMinusInvBeta[iii][jjj][kkk] = NULL;
        m_spectraPlusInvBetaETOF[iii][jjj][kkk] = NULL;
        m_spectraMinusInvBetaETOF[iii][jjj][kkk] = NULL;

        m_inclusiveSpectraPlusDeDx[iii][jjj][kkk] = NULL;
        m_inclusiveSpectraMinusDeDx[iii][jjj][kkk] = NULL;
        m_inclusiveSpectraPlusInvBeta[iii][jjj][kkk] = NULL;
        m_inclusiveSpectraMinusInvBeta[iii][jjj][kkk] = NULL;
        m_inclusiveSpectraPlusInvBetaETOF[iii][jjj][kkk] = NULL;
        m_inclusiveSpectraMinusInvBetaETOF[iii][jjj][kkk] = NULL;

      }
      for(int kkk = 0; kkk < g_RawSpecMod_MaxNumRapBins; kkk++){
        m_spectraStruct[iii][jjj][kkk][0] = NULL;
        m_spectraStruct[iii][jjj][kkk][1] = NULL;
      }

    }
  }




  // Updated 2026-07-15 (was "SL19b"): macros/RunRawSpectraModifier.C comments out
  // mod->setStarVer(...) as unused, so this constructor default is what's actually in
  // effect at runtime for this stage -- it wasn't just an inert default. Changed to
  // match the SL24y starver now used everywhere else in this analysis (see the matching
  // change/comment in macros/SetCutClass.C's setStarver(...) call), so stage 3
  // (RawSpectraModifier) reads dE/dx off the same Bichsel calibration as stage 2
  // (ZFitter) instead of silently mixing SL19b into a pipeline that's otherwise SL24y
  // throughout. Requires submodule/ParticleInfo/PidFunctions_Data/SL24y_BichselCurves.root
  // to exist (generated on SDCC via CollectBichselFunctions.csh SL24y, see SetCutClass.C).
  m_partInfo = new ParticleInfo("SL24y",true);

  m_loadedETOF = false;
}

void RawSpectraModifier::setStarVer(string a_starVer){
  m_starVer = a_starVer;
  m_partInfo = new ParticleInfo(m_starVer,true);
  return;
}
void RawSpectraModifier::setDynamicRebinForJankSpectra(vector< vector<double> > a_dynamicRebin_Plus,  vector< vector<double> > a_dynamicRebin_Minus){
  m_dynamicRebin_Plus = a_dynamicRebin_Plus;
  m_dynamicRebin_Minus = a_dynamicRebin_Minus;
}




void RawSpectraModifier::loadDataFile(string a_filename, int a_partIndex, int a_numCent,  bool a_converted16to9, bool a_combine5and10CentBins, bool a_combine40and50_60and70_centBins, bool a_useETOF, bool a_populateJankSpectra,double a_tpcEffJank){
  cout << "Starting RawSpectraModifier::loadDataFile(" << a_filename << "," << a_partIndex << "," << a_numCent << ")" << endl;
  
  m_combined16To9 = a_converted16to9;
  if(a_partIndex + 1 > m_numParticles) m_numParticles = a_partIndex + 1;
  cout << "num particles = " << m_numParticles << endl;
  m_numCentralities = a_numCent;


  bool firstOneLoaded = true;
  for(int iii = 0; iii < 9; iii++){
    if(m_particleLoaded[iii]) firstOneLoaded = false;
  }
  m_particleLoaded[a_partIndex] = true;






  cout << "Now Loading TFile: " << a_filename << endl; 
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "Warning in RawSpectraModifier::loadDataHistograms.. couldn't find file: " << a_filename << endl;
    return;
  }

  if(!inFile->IsOpen()){
    cout << " ERROR: TFILE NOT OPEN" << endl;
  }

  TH1D* centEvents = (TH1D*) inFile->Get("centEvents");
  if(centEvents){ //if every particle was done at once in a single root file
    m_centEvents[0] = centEvents;
    m_centEvents[1] = centEvents;
    m_centEvents[2] = centEvents;
    m_centEvents[3] = centEvents;
    m_centEvents[4] = centEvents;
    m_centEvents[5] = centEvents;
    m_centEvents[6] = centEvents;
    m_centEvents[7] = centEvents;
    m_centEvents[8] = centEvents;
  }else{
    string centEvents_name = Form("centEvents%s", m_partInfo->GetParticleName(a_partIndex,0).Data());
    cout << "Part Name" << m_partInfo->GetParticleName(a_partIndex,0).Data() << endl;
    cout << " Trying " << centEvents_name << endl;
    centEvents = (TH1D*) inFile->Get(centEvents_name.c_str());
    if(centEvents){ //if every particle was done at once in a single root file
      m_centEvents[a_partIndex] = centEvents;
    }else{
      cout << "Error: Can't Load " << centEvents_name << "   ptr: " << centEvents << endl;
    }

    TH1D* dca = (TH1D*) inFile->Get("dca");
    cout << " could we grab dca? " << dca << endl;
    cout << dca->GetName() << endl;



  }



  if(a_converted16to9){
    if(m_centEvents[a_partIndex]->GetNbinsX() != 16){
      cout << "you cant use a_convert16to9 if you don't have 16 centrality bins... centbins = " << m_centEvents[a_partIndex]->GetNbinsX() << endl;
    }
    TH1D* newCentEvents = new TH1D(Form("centEvents%s_16to9", m_partInfo->GetParticleName(a_partIndex,0).Data()),
     Form("Centrality Events for %s",m_partInfo->GetParticleName(a_partIndex,0).Data()), 9, -0.5, 8.5);

    // 0-5% and 5-10% don't change 
    double newBinContent = m_centEvents[a_partIndex]->GetBinContent(1);
    double newBinError = sqrt(newBinContent);
    newCentEvents->SetBinContent(1,newBinContent);
    newCentEvents->SetBinError(1,newBinError);
    newBinContent = m_centEvents[a_partIndex]->GetBinContent(2);
    newBinError = sqrt(newBinContent);
    newCentEvents->SetBinContent(2,newBinContent);
    newCentEvents->SetBinError(2,newBinError);

    //make 10-20% = 10-15% + 15-20%
    for(int centIndex = 2; centIndex < 16; centIndex+=2){
      int newCentIndex = 1+centIndex/2;
      cout << "cent: " <<  newCentIndex << " = " << centIndex << " + " << centIndex+1 << endl;
      double newBinContent = m_centEvents[a_partIndex]->GetBinContent(centIndex+1) + m_centEvents[a_partIndex]->GetBinContent(centIndex+2);
      cout << "NewBinContent: " << newBinContent  << endl;
      double newBinError = sqrt(newBinContent);
      newCentEvents->SetBinContent(newCentIndex + 1,newBinContent);
      newCentEvents->SetBinError(newCentIndex + 1,newBinError);
    }
    
    //delete m_centEvents[a_partIndex];
    cout << "CentPtrBefor: " << newCentEvents;
    m_centEvents[a_partIndex] = newCentEvents;
    cout << "  after: " << m_centEvents[a_partIndex] << endl;

  }

  if(a_combine5and10CentBins){
    TH1D* newCentEvents = new TH1D(Form("centEvents%s_5and10", m_partInfo->GetParticleName(a_partIndex,0).Data()),
           Form("Centrality Events for %s",m_partInfo->GetParticleName(a_partIndex,0).Data()), 8, -0.5, 8.5);
    // 0-5% and 5-10% don't change 
    double newBinContent = m_centEvents[a_partIndex]->GetBinContent(1) + m_centEvents[a_partIndex]->GetBinContent(2);
    double newBinError = sqrt(newBinContent);
    newCentEvents->SetBinContent(1,newBinContent);
    newCentEvents->SetBinError(1,newBinError);
    
    //make 10-20% = 10-15% + 15-20%
    for(int centIndex = 2; centIndex < 16; centIndex++){
      int newCentIndex = centIndex - 1;
      double newBinContent = m_centEvents[a_partIndex]->GetBinContent(centIndex+1);
      double newBinError = sqrt(newBinContent);
      newCentEvents->SetBinContent(newCentIndex + 1,newBinContent);
      newCentEvents->SetBinError(newCentIndex + 1,newBinError);
    }
    
    //delete m_centEvents[a_partIndex];
    cout << "CentPtrBefor: " << newCentEvents;
    m_centEvents[a_partIndex] = newCentEvents;
    cout << "  after: " << m_centEvents[a_partIndex] << endl;

  }

  if(a_combine40and50_60and70_centBins){
    if(m_centEvents[a_partIndex]->GetNbinsX() != 9 && (m_centEvents[a_partIndex]->GetNbinsX() != 8 && !a_combine5and10CentBins)){
      cout << "you cant use a_combine40and50_60and70_centBins if you don't have 8 or 9 centrality bins... centbins = " << m_centEvents[a_partIndex]->GetNbinsX() << endl;
    }
    TH1D* newCentEvents = new TH1D(Form("centEvents%s_40to80", m_partInfo->GetParticleName(a_partIndex,0).Data()),
     Form("Centrality Events for %s",m_partInfo->GetParticleName(a_partIndex,0).Data()), m_centEvents[a_partIndex]->GetNbinsX()-2, -0.5, 8.5);

    for(int centIndex=0; centIndex<m_centEvents[a_partIndex]->GetNbinsX()-4; centIndex++){
      double newBinContent = m_centEvents[a_partIndex]->GetBinContent(centIndex+1);
      double newBinError = sqrt(newBinContent);
      newCentEvents->SetBinContent(centIndex + 1, newBinContent);
      newCentEvents->SetBinError(centIndex + 1, newBinError);
    }

    for(int centIndex=m_centEvents[a_partIndex]->GetNbinsX()-4; centIndex<m_centEvents[a_partIndex]->GetNbinsX(); centIndex+=2){
      int newCentIndex = m_centEvents[a_partIndex]->GetNbinsX()-4 +(centIndex-(m_centEvents[a_partIndex]->GetNbinsX()-4))/2;
      cout << "cent: " <<  newCentIndex << " = " << centIndex << " + " << centIndex+1 << endl;
      double newBinContent = m_centEvents[a_partIndex]->GetBinContent(centIndex+1) + m_centEvents[a_partIndex]->GetBinContent(centIndex+2);
      cout << "NewBinContent: " << newBinContent  << endl;
      double newBinError = sqrt(newBinContent);
      newCentEvents->SetBinContent(newCentIndex + 1,newBinContent);
      newCentEvents->SetBinError(newCentIndex + 1,newBinError);
    }
        
    //delete m_centEvents[a_partIndex];
    cout << "CentPtrBefor: " << newCentEvents;
    m_centEvents[a_partIndex] = newCentEvents;
    cout << "  after: " << m_centEvents[a_partIndex] << endl;

  }
  #ifdef _RAWSPECMOD_DEBUG_
    cout << "done loading centEvents" << endl;
  #endif



  string firstDir = "yields";
  inFile->cd();
  const char* partName      = m_partInfo->GetParticleName(a_partIndex).Data();
  const char* partNamePlus  = m_partInfo->GetParticleName(a_partIndex,1).Data();
  const char* partNameMinus = m_partInfo->GetParticleName(a_partIndex,-1).Data();
  m_rap_mTm0_dEdx_allCent[a_partIndex] = (TH3D*) inFile->Get( Form("/%s/%s_space/zTPCPlus_%s_Cent0", firstDir.c_str(),partName,partName) );
  m_rap_mTm0_invBetaBTOF_allCent[a_partIndex] = (TH3D*) inFile->Get( Form("/%s/%s_space/zBTOFPlus_%s_Cent0", firstDir.c_str(),partName,partName) ); 
  m_rap_mTm0_invBetaETOF_allCent[a_partIndex] = (TH3D*) inFile->Get( Form("/%s/%s_space/zETOFPlus_%s_Cent0", firstDir.c_str(),partName,partName) ); 
  
  if(m_rap_mTm0_invBetaETOF_allCent[a_partIndex] && a_useETOF) m_loadedETOF = true;
  cout << "ETOF Loaded: " << m_loadedETOF << endl;
  cout <<"got the all cent histos  " << m_rap_mTm0_dEdx_allCent[a_partIndex] << " " <<  m_rap_mTm0_invBetaBTOF_allCent[a_partIndex] 
  << " " << m_rap_mTm0_invBetaETOF_allCent[a_partIndex] << endl;
  m_numEtaBins_TPC[a_partIndex]= (int) m_rap_mTm0_dEdx_allCent[a_partIndex]->GetNbinsX();
  m_numEtaBins_BTOF[a_partIndex]= (int) m_rap_mTm0_invBetaBTOF_allCent[a_partIndex]->GetNbinsX();
  if(m_rap_mTm0_invBetaETOF_allCent[a_partIndex]) m_numEtaBins_ETOF[a_partIndex]= (int) m_rap_mTm0_invBetaETOF_allCent[a_partIndex]->GetNbinsX();

  m_numPtBins_BTOF[a_partIndex]= (int) m_rap_mTm0_invBetaBTOF_allCent[a_partIndex]->GetNbinsY();
  m_numPtBins_TPC[a_partIndex]= (int) m_rap_mTm0_dEdx_allCent[a_partIndex]->GetNbinsY(); 
  if(m_rap_mTm0_invBetaETOF_allCent[a_partIndex]) m_numPtBins_ETOF[a_partIndex]= (int) m_rap_mTm0_invBetaETOF_allCent[a_partIndex]->GetNbinsY();

  int midBinTOF = m_rap_mTm0_invBetaBTOF_allCent[a_partIndex]->GetNbinsX()/2; //to accommidate FXT
  int checkBinTPC_B = m_rap_mTm0_dEdx_allCent[a_partIndex]->GetXaxis()->FindBin(m_rap_mTm0_invBetaBTOF_allCent[a_partIndex]->GetXaxis()->GetBinCenter(midBinTOF));
  m_tofIndexOffset[a_partIndex] = checkBinTPC_B-midBinTOF;
  cout << "BTOF OFFSET: " << m_tofIndexOffset[a_partIndex] << endl;
  
  if(m_rap_mTm0_invBetaETOF_allCent[a_partIndex]){
    int checkBinTPC_E = m_rap_mTm0_dEdx_allCent[a_partIndex]->GetXaxis()->FindBin( m_rap_mTm0_invBetaETOF_allCent[a_partIndex]->GetXaxis()->GetBinCenter(1) );
    m_etofIndexOffset[a_partIndex] = checkBinTPC_E - 1; // because of first bin
    cout << "ETOF OFFSET: " << m_etofIndexOffset[a_partIndex] << endl;
  }

  cout << "Loading for PartIndex: " << a_partIndex << endl;
  inFile->cd();

  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    TH2I* DTBP_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nTPCTracksPlus_ForBTOF_%s_Cent%d",
      firstDir.c_str(), partName,partName,centIndex));
    TH2I* DTBM_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nTPCTracksMinus_ForBTOF_%s_Cent%d",
      firstDir.c_str(), partName,partName,centIndex));
    TH2I* DTEP_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nTPCTracksPlus_ForETOF_%s_Cent%d",
      firstDir.c_str(), partName,partName,centIndex));
    TH2I* DTEM_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nTPCTracksMinus_ForETOF_%s_Cent%d",
      firstDir.c_str(), partName,partName,centIndex));
    TH2I* TTP_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nBTOFTracksPlus_%s_Cent%d",
      firstDir.c_str(), partName, partName,centIndex));
    TH2I* TTM_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nBTOFTracksMinus_%s_Cent%d",
      firstDir.c_str(), partName, partName,centIndex));
    TH2I* ETP_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nETOFTracksPlus_%s_Cent%d",
      firstDir.c_str(), partName, partName,centIndex));
    TH2I* ETM_int = (TH2I*) inFile->Get(Form("/%s/%s_space/nETOFTracksMinus_%s_Cent%d",
      firstDir.c_str(), partName, partName,centIndex));

    TH2I* SPEC_TP_int = (TH2I*) inFile->Get(Form("/%s/%s_space/SpectraTPC_%s_Cent%d",
      firstDir.c_str(), partName, partNamePlus,centIndex));
    TH2I* SPEC_TM_int = (TH2I*) inFile->Get(Form("/%s/%s_space/SpectraTPC_%s_Cent%d",
      firstDir.c_str(), partName, partNameMinus,centIndex));
    TH2I* SPEC_BP_int = (TH2I*) inFile->Get(Form("/%s/%s_space/SpectraBTOF_%s_Cent%d",
      firstDir.c_str(), partName, partNamePlus,centIndex));
    TH2I* SPEC_BM_int = (TH2I*) inFile->Get(Form("/%s/%s_space/SpectraBTOF_%s_Cent%d",
      firstDir.c_str(), partName, partNameMinus,centIndex));
    TH2I* SPEC_EP_int = (TH2I*) inFile->Get(Form("/%s/%s_space/SpectraETOF_%s_Cent%d",
      firstDir.c_str(), partName, partNamePlus,centIndex));
    TH2I* SPEC_EM_int = (TH2I*) inFile->Get(Form("/%s/%s_space/SpectraETOF_%s_Cent%d",
      firstDir.c_str(), partName, partNameMinus,centIndex));

    TH2D* DTBP = HistogramUtilities::convertInt_to_DoubleHisto2D(DTBP_int);
    TH2D* DTBM = HistogramUtilities::convertInt_to_DoubleHisto2D(DTBM_int);
    TH2D* TTP = HistogramUtilities::convertInt_to_DoubleHisto2D(TTP_int);
    TH2D* TTM = HistogramUtilities::convertInt_to_DoubleHisto2D(TTM_int);
    delete DTBP_int;
    delete DTBM_int;
    delete TTP_int;
    delete TTM_int;

    TH2D* DTEP = NULL;  
    TH2D* DTEM = NULL;
    TH2D* ETP  = NULL;
    TH2D* ETM  = NULL;
    if(DTEP_int){
      DTEP = HistogramUtilities::convertInt_to_DoubleHisto2D(DTEP_int);
      DTEM = HistogramUtilities::convertInt_to_DoubleHisto2D(DTEM_int);
      ETP = HistogramUtilities::convertInt_to_DoubleHisto2D(ETP_int);
      ETM = HistogramUtilities::convertInt_to_DoubleHisto2D(ETM_int);
      delete DTEP_int;
      delete DTEM_int;
      delete ETP_int;
      delete ETM_int;
    }

    if(a_populateJankSpectra){

      vector< vector< double > > rebinInfo_Plus  = m_dynamicRebin_Plus;
      vector< vector< double > > rebinInfo_Minus = m_dynamicRebin_Minus;

      //Dynamic Rebin Information
      for(int rebinIndex = rebinInfo_Plus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int attemptPartIndex = (int) rebinInfo_Plus[rebinIndex][0];
        if(attemptPartIndex != a_partIndex){
          rebinInfo_Plus.erase(rebinInfo_Plus.begin() + rebinIndex);
        }
      }
      for(int rebinIndex = rebinInfo_Minus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int attemptPartIndex = (int) rebinInfo_Minus[rebinIndex][0];
        if(attemptPartIndex != a_partIndex){
          rebinInfo_Minus.erase(rebinInfo_Minus.begin() + rebinIndex);
        }
      }
      vector< vector< double > > rebinInfo_TPC_Plus  = m_dynamicRebin_Plus;
      vector< vector< double > > rebinInfo_bTOF_Plus = m_dynamicRebin_Plus;
      vector< vector< double > > rebinInfo_eTOF_Plus = m_dynamicRebin_Plus;

      for(int rebinIndex = rebinInfo_TPC_Plus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int detIndex = (int) rebinInfo_TPC_Plus[rebinIndex][1];
        if(detIndex != 0){
          rebinInfo_TPC_Plus.erase(rebinInfo_TPC_Plus.begin() + rebinIndex);
        }else{
          rebinInfo_TPC_Plus[rebinIndex].erase(rebinInfo_TPC_Plus[rebinIndex].begin()); //partIndex
          rebinInfo_TPC_Plus[rebinIndex].erase(rebinInfo_TPC_Plus[rebinIndex].begin()); //detIndex
        }
      }
 
      for(int rebinIndex = rebinInfo_bTOF_Plus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int detIndex = (int) rebinInfo_bTOF_Plus[rebinIndex][1];
        if(detIndex != 1){
          rebinInfo_bTOF_Plus.erase(rebinInfo_bTOF_Plus.begin() + rebinIndex);
        }else{
          rebinInfo_bTOF_Plus[rebinIndex].erase(rebinInfo_bTOF_Plus[rebinIndex].begin()); //partIndex
          rebinInfo_bTOF_Plus[rebinIndex].erase(rebinInfo_bTOF_Plus[rebinIndex].begin()); //detIndex
        }
      }
 
      for(int rebinIndex = rebinInfo_eTOF_Plus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int detIndex = (int) rebinInfo_eTOF_Plus[rebinIndex][1];
        if(detIndex != 2){
          rebinInfo_eTOF_Plus.erase(rebinInfo_eTOF_Plus.begin() + rebinIndex);
        }else{
          rebinInfo_eTOF_Plus[rebinIndex].erase(rebinInfo_eTOF_Plus[rebinIndex].begin()); //partIndex
          rebinInfo_eTOF_Plus[rebinIndex].erase(rebinInfo_eTOF_Plus[rebinIndex].begin()); //detIndex
        }
      }


      vector< vector< double > > rebinInfo_TPC_Minus  = m_dynamicRebin_Minus;
      vector< vector< double > > rebinInfo_bTOF_Minus = m_dynamicRebin_Minus;
      vector< vector< double > > rebinInfo_eTOF_Minus = m_dynamicRebin_Minus;

      for(int rebinIndex = rebinInfo_TPC_Minus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int detIndex = (int) rebinInfo_TPC_Minus[rebinIndex][1];
        if(detIndex != 0){
          rebinInfo_TPC_Minus.erase(rebinInfo_TPC_Minus.begin() + rebinIndex);
        }else{
          rebinInfo_TPC_Minus[rebinIndex].erase(rebinInfo_TPC_Minus[rebinIndex].begin()); //partIndex
          rebinInfo_TPC_Minus[rebinIndex].erase(rebinInfo_TPC_Minus[rebinIndex].begin()); //detIndex
        }
      }
 
      for(int rebinIndex = rebinInfo_bTOF_Minus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int detIndex = (int) rebinInfo_bTOF_Minus[rebinIndex][1];
        if(detIndex != 1){
          rebinInfo_bTOF_Minus.erase(rebinInfo_bTOF_Minus.begin() + rebinIndex);
        }else{
          rebinInfo_bTOF_Minus[rebinIndex].erase(rebinInfo_bTOF_Minus[rebinIndex].begin()); //partIndex
          rebinInfo_bTOF_Minus[rebinIndex].erase(rebinInfo_bTOF_Minus[rebinIndex].begin()); //detIndex
        }
      }
 
      for(int rebinIndex = rebinInfo_eTOF_Minus.size() - 1; rebinIndex >= 0; rebinIndex--){
        int detIndex = (int) rebinInfo_eTOF_Minus[rebinIndex][1];
        if(detIndex != 2){
          rebinInfo_eTOF_Minus.erase(rebinInfo_eTOF_Minus.begin() + rebinIndex);
        }else{
          rebinInfo_eTOF_Minus[rebinIndex].erase(rebinInfo_eTOF_Minus[rebinIndex].begin()); //partIndex
          rebinInfo_eTOF_Minus[rebinIndex].erase(rebinInfo_eTOF_Minus[rebinIndex].begin()); //detIndex
        }
      }


      TH2D* SPEC_TP = NULL; 
      TH2D* SPEC_TM = NULL; 
      TH2D* SPEC_BP = NULL; 
      TH2D* SPEC_BM = NULL; 
      TH2D* SPEC_EP = NULL; 
      TH2D* SPEC_EM = NULL; 
      if(SPEC_TP_int) SPEC_TP = HistogramUtilities::convertInt_to_DoubleHisto2D(SPEC_TP_int); 
      if(SPEC_TM_int) SPEC_TM = HistogramUtilities::convertInt_to_DoubleHisto2D(SPEC_TM_int); 
      if(SPEC_BP_int) SPEC_BP = HistogramUtilities::convertInt_to_DoubleHisto2D(SPEC_BP_int); 
      if(SPEC_BM_int) SPEC_BM = HistogramUtilities::convertInt_to_DoubleHisto2D(SPEC_BM_int); 
      if(SPEC_EP_int) SPEC_EP = HistogramUtilities::convertInt_to_DoubleHisto2D(SPEC_EP_int); 
      if(SPEC_EM_int) SPEC_EM = HistogramUtilities::convertInt_to_DoubleHisto2D(SPEC_EM_int); 
      if(SPEC_TP_int) delete SPEC_TP_int;
      if(SPEC_TM_int) delete SPEC_TM_int;
      if(SPEC_BP_int) delete SPEC_BP_int; 
      if(SPEC_BM_int) delete SPEC_BM_int; 
      if(SPEC_EP_int) delete SPEC_EP_int; 
      if(SPEC_EM_int) delete SPEC_EM_int; 

      vector< TH2D* > histosToDelete_double;
      if(rebinInfo_TPC_Plus.size() > 0 && SPEC_TP){
        histosToDelete_double.push_back(SPEC_TP);
        SPEC_TP = HistogramUtilities::reBinBar2D_specificY(SPEC_TP, rebinInfo_TPC_Plus);
      }
      if(rebinInfo_TPC_Minus.size() > 0 && SPEC_TM){
        histosToDelete_double.push_back(SPEC_TM);
        SPEC_TM = HistogramUtilities::reBinBar2D_specificY(SPEC_TM, rebinInfo_TPC_Minus);
      }
      if(rebinInfo_bTOF_Plus.size() > 0 && SPEC_BP){
        histosToDelete_double.push_back(SPEC_BP);
        SPEC_BP = HistogramUtilities::reBinBar2D_specificY(SPEC_BP, rebinInfo_bTOF_Plus);
      }
      if(rebinInfo_bTOF_Minus.size() > 0 && SPEC_BM){
        histosToDelete_double.push_back(SPEC_BM);
        SPEC_BM = HistogramUtilities::reBinBar2D_specificY(SPEC_BM, rebinInfo_bTOF_Minus);
      }
      if(rebinInfo_eTOF_Plus.size() > 0){
        histosToDelete_double.push_back(SPEC_EP);
        SPEC_EP = HistogramUtilities::reBinBar2D_specificY(SPEC_EP, rebinInfo_eTOF_Plus);
      }
      if(rebinInfo_eTOF_Minus.size() > 0){
        histosToDelete_double.push_back(SPEC_EM);
        SPEC_EM = HistogramUtilities::reBinBar2D_specificY(SPEC_EM, rebinInfo_eTOF_Minus);
      }
      for(int iii = 0; iii < histosToDelete_double.size(); iii++){
        delete histosToDelete_double[iii];
      }

      m_dEdxSpectraPlus[a_partIndex][centIndex]         = SPEC_TP;     
      m_dEdxSpectraMinus[a_partIndex][centIndex]        = SPEC_TM;  

      if(a_tpcEffJank != 1.0){
        for(int binX = 1; binX <= SPEC_TP->GetNbinsX(); binX++){
          for(int binY = 1; binY <= SPEC_TP->GetNbinsY(); binY++){
            SPEC_TP->SetBinContent(binX,binY,SPEC_TP->GetBinContent(binX,binY)/a_tpcEffJank);
            SPEC_TP->SetBinError(binX,binY,SPEC_TP->GetBinError(binX,binY)/a_tpcEffJank);
            SPEC_TM->SetBinContent(binX,binY,SPEC_TM->GetBinContent(binX,binY)/a_tpcEffJank);
            SPEC_TM->SetBinError(binX,binY,SPEC_TM->GetBinError(binX,binY)/a_tpcEffJank);
          }         
        }

      }


      m_InvBetaBTOFSpectraPlus[a_partIndex][centIndex]  = SPEC_BP; 
      m_InvBetaBTOFSpectraMinus[a_partIndex][centIndex] = SPEC_BM;
      m_InvBetaETOFSpectraPlus[a_partIndex][centIndex]  = SPEC_EP; 
      m_InvBetaETOFSpectraMinus[a_partIndex][centIndex] = SPEC_EM;

      if(centIndex == 0){

        if(m_dEdxSpectraPlus[a_partIndex][centIndex] && !m_picoBinner_rapidity_structure[a_partIndex][0]){
          vector< double > edges = HistogramUtilities::getBinEdges_X(m_dEdxSpectraPlus[a_partIndex][centIndex]);
          m_picoBinner_rapidity_structure[a_partIndex][0] = new TH1D(Form("picoBinner_%s_tpc_rap_struct",
               m_partInfo->GetParticleName(a_partIndex,0).Data()),"rapStruct;y",edges.size()-1,&edges[0]);
        }

        if(m_InvBetaBTOFSpectraPlus[a_partIndex][centIndex] && !m_picoBinner_rapidity_structure[a_partIndex][1]){
          vector< double > edges = HistogramUtilities::getBinEdges_X(m_InvBetaBTOFSpectraPlus[a_partIndex][centIndex]);
          m_picoBinner_rapidity_structure[a_partIndex][1] = new TH1D(Form("picoBinner_%s_btof_rap_struct",
               m_partInfo->GetParticleName(a_partIndex,0).Data()),"rapStruct;y",edges.size()-1,&edges[0]);
        }
        if(m_InvBetaETOFSpectraPlus[a_partIndex][centIndex] && !m_picoBinner_rapidity_structure[a_partIndex][2]){
          vector< double > edges = HistogramUtilities::getBinEdges_X(m_InvBetaETOFSpectraPlus[a_partIndex][centIndex]);
          m_picoBinner_rapidity_structure[a_partIndex][2] = new TH1D(Form("picoBinner_%s_etof_rap_struct",
               m_partInfo->GetParticleName(a_partIndex,0).Data()),"rapStruct;y",edges.size()-1,&edges[0]);
        }
      }




    } // end of populate jank


    m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[a_partIndex][centIndex] = DTBP;
    m_rap_mTm0_NumTPCTracks_ForETOF_Plus[a_partIndex][centIndex] = DTEP;
    m_rap_mTm0_NumBTOFTracks_Plus[a_partIndex][centIndex] = TTP;
    m_rap_mTm0_NumETOFTracks_Plus[a_partIndex][centIndex] = ETP;
    m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[a_partIndex][centIndex] = DTBM;
    m_rap_mTm0_NumTPCTracks_ForETOF_Minus[a_partIndex][centIndex] = DTEM;
    m_rap_mTm0_NumBTOFTracks_Minus[a_partIndex][centIndex] = TTM;
    m_rap_mTm0_NumETOFTracks_Minus[a_partIndex][centIndex] = ETM;
  }

  //##########     COMBINE 16 to 9 Centrality Bins  ################
  if(a_converted16to9){
    cout << "  Converting 16 to 9 Cent Bins..." << endl;
    if(!m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[a_partIndex][15]){
      cout << "  ERROR: You cant use the a_convert16to9CentBins unless you actually have 16 bins in the BTOF Efficiency file..." << endl;
      return;
    }

    for(int centIndex = 2; centIndex < 16; centIndex+=2){
      int newCentIndex = 1+centIndex/2;
      m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[a_partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[a_partIndex][centIndex+1]);
      m_rap_mTm0_NumBTOFTracks_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumBTOFTracks_Plus[a_partIndex][centIndex],m_rap_mTm0_NumBTOFTracks_Plus[a_partIndex][centIndex+1]);
      m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[a_partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[a_partIndex][centIndex+1]);
      m_rap_mTm0_NumBTOFTracks_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumBTOFTracks_Minus[a_partIndex][centIndex],m_rap_mTm0_NumBTOFTracks_Minus[a_partIndex][centIndex+1]);
      if(m_rap_mTm0_NumTPCTracks_ForETOF_Plus[a_partIndex][newCentIndex]){
        m_rap_mTm0_NumTPCTracks_ForETOF_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumTPCTracks_ForETOF_Plus[a_partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForETOF_Plus[a_partIndex][centIndex+1]);
        m_rap_mTm0_NumETOFTracks_Plus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumETOFTracks_Plus[a_partIndex][centIndex],m_rap_mTm0_NumETOFTracks_Plus[a_partIndex][centIndex+1]);
        m_rap_mTm0_NumTPCTracks_ForETOF_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumTPCTracks_ForETOF_Minus[a_partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForETOF_Minus[a_partIndex][centIndex+1]);
        m_rap_mTm0_NumETOFTracks_Minus[a_partIndex][newCentIndex]->Add(m_rap_mTm0_NumETOFTracks_Minus[a_partIndex][centIndex],m_rap_mTm0_NumETOFTracks_Minus[a_partIndex][centIndex+1]);
      }
    }

    for(int centIndex = 9; centIndex < 16; centIndex++){
      delete m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[a_partIndex][centIndex];
      delete m_rap_mTm0_NumBTOFTracks_Plus[a_partIndex][centIndex];
      delete m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[a_partIndex][centIndex];
      delete m_rap_mTm0_NumBTOFTracks_Minus[a_partIndex][centIndex];
      m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[a_partIndex][centIndex]   = NULL;
      m_rap_mTm0_NumBTOFTracks_Plus[a_partIndex][centIndex]          = NULL;
      m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[a_partIndex][centIndex]  = NULL;
      m_rap_mTm0_NumBTOFTracks_Minus[a_partIndex][centIndex]         = NULL;
    }

  }


  if(a_combine5and10CentBins || a_combine40and50_60and70_centBins){
    cout << "\n\n\n\n WARNING !!!!!!!!!!      BTOF EFFICIENCY GRAPHS HAVE NOT BEEN COMBINED with YOUR a_combine5and10CentBins a_combine40and50_60and70_centBins options \n\n\n" << endl;
  }




  cout << "Finished loading cent data" << endl;

  for(int partIndex = 0; partIndex < 9; partIndex++ ){
    for(int centIndex = 0; centIndex < 16; centIndex++ ){
      if(m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex]){
        cout << m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex]->GetName() << "       " << m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex]->GetTitle() << endl;
      }
    }
  }
}


void RawSpectraModifier::loadSpectraFile(string a_filename){
  cout << "Starting RawSpectraModifier::loadSpectraFile(" << a_filename <<")"<<endl;
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "ERROR: in RawSpectraModifier::loadSpectraFile.. couldn't find file: " << a_filename << endl;
    return;
  }

  if(!inFile->IsOpen()){
    cout << "ERROR: in RawSpectraModifier::loadSpectraFile.. file is not open: " << a_filename << endl;
    return;
  }
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){


    if(!m_particleLoaded[partIndex]){
      cout << "Particle " << partIndex << "'s data file from PicoBinner has not been loaded, so skipping" << endl;
      continue;
    }
    cout << " part " << partIndex << " picoBinner file has been loaded, so continuing to get spectra" << endl;
    const char* partName = m_partInfo->GetParticleName(partIndex).Data();

    /*m_picoBinner_rapidity_structure[partIndex][0] = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_tpc_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
    m_picoBinner_rapidity_structure[partIndex][1] = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_btof_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
    m_picoBinner_rapidity_structure[partIndex][2] = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_etof_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
    */

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

      // Fixed 2026-07-15: these Form() strings were missing the underscore before
      // Plus/Minus that ZFitter.cxx's makeSpectra() actually writes (e.g.
      // "dEdxSpectra_Pion_Plus_Cent00", not "dEdxSpectra_PionPlus_Cent00" -- see
      // ZFitter.cxx's SetName(Form("dEdxSpectra_%s_Plus_Cent%02d",...)) and the
      // InvBetaBTOFSpectra_/InvBetaETOFSpectra_ equivalents). TFile::Get() silently
      // returns NULL on a name mismatch rather than erroring, so every one of these
      // six lookups was quietly returning NULL and would have surfaced later as a
      // null-pointer use rather than a clear "not found" message. Added the missing
      // underscores so these actually match what gets written.
      m_dEdxSpectraPlus[partIndex][centIndex]         = (TH2D*) inFile->Get(Form("DeDx_FitData/%s/dEdxSpectra_%s_Plus_Cent%02d",partName,partName,centIndex));
      m_dEdxSpectraMinus[partIndex][centIndex]        = (TH2D*) inFile->Get(Form("DeDx_FitData/%s/dEdxSpectra_%s_Minus_Cent%02d",partName,partName,centIndex));
      m_InvBetaBTOFSpectraPlus[partIndex][centIndex]  = (TH2D*) inFile->Get(Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%s_Plus_Cent%02d",partName,partName,centIndex));
      m_InvBetaBTOFSpectraMinus[partIndex][centIndex] = (TH2D*) inFile->Get(Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%s_Minus_Cent%02d",partName,partName,centIndex));
      m_InvBetaETOFSpectraPlus[partIndex][centIndex]  = (TH2D*) inFile->Get(Form("ETOF_FitData/%s/InvBetaETOFSpectra_%s_Plus_Cent%02d",partName,partName,centIndex));
      m_InvBetaETOFSpectraMinus[partIndex][centIndex] = (TH2D*) inFile->Get(Form("ETOF_FitData/%s/InvBetaETOFSpectra_%s_Minus_Cent%02d",partName,partName,centIndex));

      if(centIndex == 0){
        if(m_dEdxSpectraPlus[partIndex][centIndex]){
          vector< double > edges = HistogramUtilities::getBinEdges_X(m_dEdxSpectraPlus[partIndex][centIndex]);
          m_picoBinner_rapidity_structure[partIndex][0] = new TH1D(Form("picoBinner_%s_tpc_rap_struct",
               m_partInfo->GetParticleName(partIndex,0).Data()),"rapStruct;y",edges.size()-1,&edges[0]);
        }
        if(m_InvBetaBTOFSpectraPlus[partIndex][centIndex]){
          vector< double > edges = HistogramUtilities::getBinEdges_X(m_InvBetaBTOFSpectraPlus[partIndex][centIndex]);
          m_picoBinner_rapidity_structure[partIndex][1] = new TH1D(Form("picoBinner_%s_btof_rap_struct",
               m_partInfo->GetParticleName(partIndex,0).Data()),"rapStruct;y",edges.size()-1,&edges[0]);
        }
        if(m_InvBetaETOFSpectraPlus[partIndex][centIndex]){
          vector< double > edges = HistogramUtilities::getBinEdges_X(m_InvBetaETOFSpectraPlus[partIndex][centIndex]);
          m_picoBinner_rapidity_structure[partIndex][2] = new TH1D(Form("picoBinner_%s_etof_rap_struct",
               m_partInfo->GetParticleName(partIndex,0).Data()),"rapStruct;y",edges.size()-1,&edges[0]);
        }
      }

      #ifdef _RAWSPECMOD_DEBUG_
        cout << " Spectra TH2D Ptrs: " << m_dEdxSpectraPlus[partIndex][centIndex]  << "  " << m_dEdxSpectraMinus[partIndex][centIndex] << "  " << m_InvBetaBTOFSpectraPlus[partIndex][centIndex]  
             << "  " << m_InvBetaBTOFSpectraMinus[partIndex][centIndex] << "  " << m_InvBetaETOFSpectraPlus[partIndex][centIndex]  << "  " << m_InvBetaETOFSpectraMinus[partIndex][centIndex] << endl;
      #endif

    } // end cent
  } // end part


  /*for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    vector< vector< TGraphErrors* > > tempCentRap;
    m_spectraPlusDeDx.push_back(tempCentRap);
    m_spectraMinusDeDx.push_back(tempCentRap);
    m_spectraPlusInvBeta.push_back(tempCentRap);
    m_spectraMinusInvBeta.push_back(tempCentRap);
    m_spectraPlusInvBetaETOF.push_back(tempCentRap);
    m_spectraMinusInvBetaETOF.push_back(tempCentRap);        

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      vector< TGraphErrors* > tempRap;
      m_spectraPlusDeDx[partIndex].push_back(tempRap);
      m_spectraMinusDeDx[partIndex].push_back(tempRap);
      m_spectraPlusInvBeta[partIndex].push_back(tempRap);
      m_spectraMinusInvBeta[partIndex].push_back(tempRap);
      m_spectraPlusInvBetaETOF[partIndex].push_back(tempRap);
      m_spectraMinusInvBetaETOF[partIndex].push_back(tempRap);

      //TPC
      for(int HorBinIndex = 1; HorBinIndex <= m_numEtaBins_TPC; HorBinIndex++){

        TGraphErrors* DP = (TGraphErrors*) inFile->Get( Form("RawSpectra_%s/rawSpectra_%s_Cent%02d_yIndex%02d",
              m_partInfo->GetParticleName(partIndex,1).Data(),
              m_partInfo->GetParticleName(partIndex,1).Data(),
              centIndex,HorBinIndex - 1)); 
          
        TGraphErrors* DM = (TGraphErrors*) inFile->Get( Form("RawSpectra_%s/rawSpectra_%s_Cent%02d_yIndex%02d",
              m_partInfo->GetParticleName(partIndex,-1).Data(),
              m_partInfo->GetParticleName(partIndex,-1).Data(),
              centIndex,HorBinIndex - 1));

        m_spectraPlusDeDx[partIndex][centIndex].push_back(DP);
        m_spectraMinusDeDx[partIndex][centIndex].push_back(DM);
        
      }//end rapbin loop
      
      //BTOF
      for(int HorBinIndex = 1; HorBinIndex <= m_numEtaBins_BTOF; HorBinIndex++){
        TGraphErrors* IP = (TGraphErrors*) inFile->Get( Form("RawSpectra_%s/rawSpectraTOF_%s_Cent%02d_yIndex%02d",
              m_partInfo->GetParticleName(partIndex,1).Data(),
              m_partInfo->GetParticleName(partIndex,1).Data(),
              centIndex,HorBinIndex - 1));

        TGraphErrors* IM = (TGraphErrors*) inFile->Get( Form("RawSpectra_%s/rawSpectraTOF_%s_Cent%02d_yIndex%02d",
              m_partInfo->GetParticleName(partIndex,-1).Data(),
              m_partInfo->GetParticleName(partIndex,-1).Data(),
              centIndex,HorBinIndex - 1));
        m_spectraPlusInvBeta[partIndex][centIndex].push_back(IP);
        m_spectraMinusInvBeta[partIndex][centIndex].push_back(IM);

      }//end rapbin loop

      //ETOF
      for(int HorBinIndex = 1; HorBinIndex <= m_numEtaBins_ETOF; HorBinIndex++){
        TGraphErrors* IP = (TGraphErrors*) inFile->Get( Form("RawSpectra_%s/rawSpectraETOF_%s_Cent%02d_yIndex%02d",
              m_partInfo->GetParticleName(partIndex,1).Data(),
              m_partInfo->GetParticleName(partIndex,1).Data(),
              centIndex,HorBinIndex - 1));

        TGraphErrors* IM = (TGraphErrors*) inFile->Get( Form("RawSpectra_%s/rawSpectraETOF_%s_Cent%02d_yIndex%02d",
              m_partInfo->GetParticleName(partIndex,-1).Data(),
              m_partInfo->GetParticleName(partIndex,-1).Data(),
              centIndex,HorBinIndex - 1));

        if(!IP){
          cout << "Warning: Plus ETOF Eta Bin " << HorBinIndex << " was read as NULL" << endl;
        }
        if(!IM){
          cout << "Warning: Minus ETOF Eta Bin " << HorBinIndex << " was read as NULL" << endl;
        }


        m_spectraPlusInvBetaETOF[partIndex][centIndex].push_back(IP);
        m_spectraMinusInvBetaETOF[partIndex][centIndex].push_back(IM);

      }//end rapbin loop

    } // end cent
  } // end part
  */

}

void RawSpectraModifier::loadAndApplyTPCEffAndEnergyLossAndBTOFEffFile(string a_filename, bool a_useTPCEffSystematic, bool a_useBTOFEffSystematic, bool a_useEnergyLossSystematic, int* m_offsets, bool  a_useCent16ForCent9,
  bool a_useTPCEffDataPoints ,bool a_useBTOFEffDataPoints , bool a_doEnergyLoss, bool a_doTPCEff, bool a_doBTOFEff, bool a_useEnergyLossErrorAsXErrorBars){
  bool doSplineFitting = false;
  int offsetsHARDCODE[9] = {0,0,0,0,0,0,0,0,0};
  if(m_offsets){
    for(int iii = 0; iii < 9; iii++){
      offsetsHARDCODE[iii] = m_offsets[iii];
    }
    cout << "WARNING: Offsets probably shouldn't be applied! This is done automatically now" << endl;
  }
  //else{
  //  cout << "HARDCODED offsets, CODE MUST BE CHANGED ERROR ERROR, or you have matching indecies in your correction file and spectra..." << endl;
  //}
  

  #ifdef _RAWSPECMOD_DEBUG_
    cout << "Starting RawSpectraModifier::loadTPCEffAndEnergyLossFile(" << a_filename << "," << a_useTPCEffSystematic<< "," << a_useBTOFEffSystematic<< "," << a_useEnergyLossSystematic << "," << m_offsets << ")" <<endl;
  #endif
  
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "ERROR: Warning in RawSpectraModifier::loadTPCEffAndEnergyLossFile couldn't find file: " << a_filename << endl;
    return;
  }

  bool enbed_needed = (a_doEnergyLoss || a_doTPCEff);

  //###################    make the TH2D of the efficiencies  #################################
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(!m_particleLoaded[partIndex]) continue;
    #ifdef _RAWSPECMOD_DEBUG_
      cout << "Histogram Structure = m_dEdxSpectraPlus[" << partIndex << "][0]: " << m_dEdxSpectraPlus[partIndex][0] << endl;
    #endif
    TH2D* histo_struct = (TH2D*) m_dEdxSpectraPlus[partIndex][0]->Clone();
    HistogramUtilities::zeroBins2D(histo_struct);
    histo_struct->SetName("struct");
    histo_struct->GetXaxis()->SetTitle("y");
    histo_struct->GetYaxis()->SetTitle(Form("m_{T}-m_{%s}",m_partInfo->GetParticleSymbol(partIndex,0).Data()));  

    // --------  READ THE RAPIDITY STRUCTURES  ------------------------
    TH1D* rap_tpc = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_tpc_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
    TH1D* rap_btof = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_btof_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
    TH1D* rap_etof = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_etof_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
    TH1D* rap_embed = (TH1D*) inFile->Get(Form("%s/embedding_%s_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
  
    if(rap_tpc && !m_picoBinner_rapidity_structure[partIndex][0])   m_picoBinner_rapidity_structure[partIndex][0] = (TH1D*) rap_tpc->Clone();
    if(rap_btof && !m_picoBinner_rapidity_structure[partIndex][1])  m_picoBinner_rapidity_structure[partIndex][1] = (TH1D*) rap_btof->Clone();
    if(rap_etof && !m_picoBinner_rapidity_structure[partIndex][2])  m_picoBinner_rapidity_structure[partIndex][2] = (TH1D*) rap_etof->Clone();
    if(rap_embed) m_embedding_rapidity_structure[partIndex]     = (TH1D*) rap_embed->Clone();

    #ifdef _RAWSPECMOD_DEBUG_
      cout << "  Rapidity Structure Pointers...  " << endl;
      cout << "    TPC: " << m_picoBinner_rapidity_structure[partIndex][0]  << endl;
      cout << "    BTOF: " << m_picoBinner_rapidity_structure[partIndex][1]  << endl;
      cout << "    ETOF: " << m_picoBinner_rapidity_structure[partIndex][2] << endl;
      cout << "    EMBED: " << m_embedding_rapidity_structure[partIndex] << endl;
      HistogramUtilities::printHistoInfo(m_picoBinner_rapidity_structure[partIndex][0]);
      HistogramUtilities::printHistoInfo(m_picoBinner_rapidity_structure[partIndex][1]);
      HistogramUtilities::printHistoInfo(m_picoBinner_rapidity_structure[partIndex][2]);
      HistogramUtilities::printHistoInfo(m_embedding_rapidity_structure[partIndex]);
    #endif


    if(!a_useCent16ForCent9){
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        m_TPCEfficiency_Nominal_Plus[partIndex][centIndex]     = (TH2D*) histo_struct->Clone();
        m_TPCEfficiency_Systematic_Plus[partIndex][centIndex]  = (TH2D*) histo_struct->Clone();
        m_TPCEfficiency_Nominal_Minus[partIndex][centIndex]    = (TH2D*) histo_struct->Clone();
        m_TPCEfficiency_Systematic_Minus[partIndex][centIndex] = (TH2D*) histo_struct->Clone();
        m_TPCEfficiency_Nominal_Plus[partIndex][centIndex]->SetName(Form("TPCEff_Nom_Plus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndex));
        m_TPCEfficiency_Nominal_Plus[partIndex][centIndex]->SetTitle(Form("Nominal TPC Efficiency for Positive %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
        m_TPCEfficiency_Systematic_Plus[partIndex][centIndex]->SetName(Form("TPCEff_Sys_Plus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndex));
        m_TPCEfficiency_Systematic_Plus[partIndex][centIndex]->SetTitle(Form("Systematic TPC Efficiency for Positive %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
        m_TPCEfficiency_Nominal_Minus[partIndex][centIndex]->SetName(Form("TPCEff_Nom_Minus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndex));
        m_TPCEfficiency_Nominal_Minus[partIndex][centIndex]->SetTitle(Form("Nominal TPC Efficiency for Negative %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
        m_TPCEfficiency_Systematic_Minus[partIndex][centIndex]->SetName(Form("TPCEff_Sys_Minus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndex));
        m_TPCEfficiency_Systematic_Minus[partIndex][centIndex]->SetTitle(Form("Systematic TPC Efficiency for Negative %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
      }//cent loop end
    }else{
      cout << "WARNING: You are using a 16 cent to 9 cent approximation for the efficiency" << endl;
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        int centIndexToGrab = centIndex;
        if(centIndex == 0) centIndexToGrab = 0;
        if(centIndex == 1) centIndexToGrab = 1;
        if(centIndex == 2) centIndexToGrab = 2;
        if(centIndex == 3) centIndexToGrab = 4;
        if(centIndex == 4) centIndexToGrab = 6;
        if(centIndex == 5) centIndexToGrab = 8;
        if(centIndex == 6) centIndexToGrab = 10;
        if(centIndex == 7) centIndexToGrab = 12;
        if(centIndex == 8) centIndexToGrab = 14;

        m_TPCEfficiency_Nominal_Plus[partIndex][centIndex] = (TH2D*) histo_struct->Clone();
        m_TPCEfficiency_Systematic_Plus[partIndex][centIndex] = (TH2D*) histo_struct->Clone();
        m_TPCEfficiency_Nominal_Minus[partIndex][centIndex] = (TH2D*) histo_struct->Clone();
        m_TPCEfficiency_Systematic_Minus[partIndex][centIndex] = (TH2D*) histo_struct->Clone();
  
        m_TPCEfficiency_Nominal_Plus[partIndex][centIndex]->SetName(Form("TPCEff_Nom_Plus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndexToGrab));
        m_TPCEfficiency_Nominal_Plus[partIndex][centIndex]->SetTitle(Form("Nominal TPC Efficiency for Positive %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
        m_TPCEfficiency_Systematic_Plus[partIndex][centIndex]->SetName(Form("TPCEff_Sys_Plus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndexToGrab));
        m_TPCEfficiency_Systematic_Plus[partIndex][centIndex]->SetTitle(Form("Systematic TPC Efficiency for Positive %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
        m_TPCEfficiency_Nominal_Minus[partIndex][centIndex]->SetName(Form("TPCEff_Nom_Minus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndexToGrab));
        m_TPCEfficiency_Nominal_Minus[partIndex][centIndex]->SetTitle(Form("Nominal TPC Efficiency for Negative %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
        m_TPCEfficiency_Systematic_Minus[partIndex][centIndex]->SetName(Form("TPCEff_Sys_Minus_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,0).Data(), centIndexToGrab));
        m_TPCEfficiency_Systematic_Minus[partIndex][centIndex]->SetTitle(Form("Systematic TPC Efficiency for Negative %ss; m_{T}^{emb}-m_{0}; TPC Efficiency",m_partInfo->GetParticleName(partIndex,0).Data()));
      }//cent loop end
    }

    m_EnergyLoss_Nominal_Plus[partIndex]     = (TH2D*) histo_struct->Clone();
    m_EnergyLoss_Systematic_Plus[partIndex]  = (TH2D*) histo_struct->Clone();
    m_EnergyLoss_Nominal_Minus[partIndex]    = (TH2D*) histo_struct->Clone();
    m_EnergyLoss_Systematic_Minus[partIndex] = (TH2D*) histo_struct->Clone();

    m_EnergyLoss_Nominal_Plus[partIndex]->SetName(Form("EnergyLoss_Nom_Plus_%s",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_EnergyLoss_Nominal_Plus[partIndex]->SetTitle(Form("Nominal Energy Loss Efficiency for Positive %ss; m_{T}-m_{0}; m_{T}^{emb}-m_{T}^{reco} ",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_EnergyLoss_Systematic_Plus[partIndex]->SetName(Form("EnergyLoss_Sys_Plus_%s",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_EnergyLoss_Systematic_Plus[partIndex]->SetTitle(Form("Systematic Energy Loss Efficiency for Positive %ss; m_{T}-m_{0}; m_{T}^{emb}-m_{T}^{reco} ",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_EnergyLoss_Nominal_Minus[partIndex]->SetName(Form("EnergyLoss_Nom_Minus_%s",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_EnergyLoss_Nominal_Minus[partIndex]->SetTitle(Form("Nominal Energy Loss Efficiency for Negative %ss; m_{T}-m_{0}; m_{T}^{emb}-m_{T}^{reco} ",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_EnergyLoss_Systematic_Minus[partIndex]->SetName(Form("EnergyLoss_Sys_Minus_%s",m_partInfo->GetParticleName(partIndex,0).Data()));
    m_EnergyLoss_Systematic_Minus[partIndex]->SetTitle(Form("Systematic Energy Loss Efficiency for Negative %ss; m_{T}-m_{0}; m_{T}^{emb}-m_{T}^{reco} ",m_partInfo->GetParticleName(partIndex,0).Data()));


  }//part

  #ifdef _RAWSPECMOD_DEBUG_
    cout << " Finished Making the TPC and Energy Loss Correction TH2Ds" << endl;
  #endif

  //###############     PARICLE LOOP        ####################
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){   
    if(!m_particleLoaded[partIndex] || !m_picoBinner_rapidity_structure[partIndex][0] || (enbed_needed &&!m_embedding_rapidity_structure[partIndex])){
      #ifdef _RAWSPECMOD_DEBUG_
        cout << "Trying Part Index " << partIndex << " but continuing because don't have all stuff" << endl;
      #endif
      continue;
    }
    double mass = m_partInfo->GetParticleMass(partIndex);

    //#### EnergyLoss and TPC Efficiency
    //for(int yIndex = 0; yIndex < m_dEdxSpectraPlus[partIndex][0]->GetNbinsX(); yIndex++){
    for(int y_Counter = 0; y_Counter < 101; y_Counter++){ // goes from -5 to 5
      double rapidity = -5.0 + 0.1*((double) y_Counter);
      int yBin_tpc  = m_picoBinner_rapidity_structure[partIndex][0]->FindBin(rapidity);
      int yBin_btof = m_picoBinner_rapidity_structure[partIndex][1]->FindBin(rapidity);
      int yBin_etof = 0;
      if(m_loadedETOF) yBin_etof = m_picoBinner_rapidity_structure[partIndex][2]->FindBin(rapidity);
      int yBin_embedding = 0;
      if(enbed_needed) yBin_embedding = m_embedding_rapidity_structure[partIndex]->FindBin(rapidity);

      int yIndex_tpc = yBin_tpc - 1;
      int yIndex_btof = yBin_btof - 1;
      int yIndex_etof = yBin_etof - 1;
      int yIndex_embedding = yBin_embedding - 1;
  
      if(   HistogramUtilities::outsideHistogramRange(m_picoBinner_rapidity_structure[partIndex][0],rapidity) 
         && HistogramUtilities::outsideHistogramRange(m_picoBinner_rapidity_structure[partIndex][1],rapidity)
         && (!enbed_needed || (enbed_needed && HistogramUtilities::outsideHistogramRange(m_embedding_rapidity_structure[partIndex],rapidity)))){
        #ifdef _RAWSPECMOD_DEBUG_
          cout << " Not Computing rapidity = " << rapidity << " becasue it is not in bounds of any histogram..." << endl;
          continue;
        #endif
      }


      //  ----------------  READ ENERGY LOSS FITS AND GRAPHS  ------------------------
      #ifdef _RAWSPECMOD_DEBUG_
        cout << "Trying to read: " << Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),yIndex_embedding) << endl;
      #endif

      TF1* energy_loss_nom_plus_fit = (TF1*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d", 
        m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),yIndex_embedding)   );
      TF1* energy_loss_sys_plus_fit = (TF1*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d_Sys", 
        m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),yIndex_embedding)  );

      TF1* energy_loss_nom_minus_fit = (TF1*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d", 
        m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),yIndex_embedding) );
      TF1* energy_loss_sys_minus_fit = (TF1*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d_Sys", 
        m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),yIndex_embedding) );

      TGraphErrors* energy_loss_nom_plus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d_Conf", 
        m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),yIndex_embedding)   );
      TGraphErrors* energy_loss_sys_plus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d_Sys_Conf", 
        m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),yIndex_embedding)  );
      TGraphErrors* energy_loss_nom_minus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d_Conf", 
        m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),yIndex_embedding) );
      TGraphErrors* energy_loss_sys_minus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EnergyLossFits/energyLossFit_%s_yIndex%d_Sys_Conf", 
        m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),yIndex_embedding) );

      TMatrixDSym* energy_loss_nom_plus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EnergyLossFits/energyLossCovMat_%s_AllCent_yIndex%d",
        m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),yIndex_embedding));
      TMatrixDSym* energy_loss_sys_plus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EnergyLossFits/energyLossCovMat_%s_AllCent_yIndex%d_Sys",
        m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),yIndex_embedding));
      TMatrixDSym* energy_loss_nom_minus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EnergyLossFits/energyLossCovMat_%s_AllCent_yIndex%d",
        m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),yIndex_embedding));
      TMatrixDSym* energy_loss_sys_minus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EnergyLossFits/energyLossCovMat_%s_AllCent_yIndex%d_Sys",
        m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),yIndex_embedding));

      #ifdef _RAWSPECMOD_DEBUG_
        cout << "Pointers For Energy Loss TF1:  " << energy_loss_nom_plus_fit << "  " << energy_loss_sys_plus_fit << "  " << energy_loss_nom_minus_fit << "  " << energy_loss_sys_minus_fit << endl;
        cout << "         For Conf TGraphs: " <<energy_loss_nom_plus_fit_conf << "  " << energy_loss_sys_plus_fit_conf << "  " << energy_loss_nom_minus_fit_conf << "  " << energy_loss_sys_minus_fit_conf << endl;
        cout << "         For Cov Matrix: " << energy_loss_nom_plus_fit_covMat << "  " << energy_loss_sys_plus_fit_covMat << "  " << energy_loss_nom_minus_fit_covMat << "  " << energy_loss_sys_minus_fit_covMat << endl;
      #endif


      //  ------- Modifiy the spectra in each centrality  ------------------------
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        int centIndexToGrab = centIndex;
        if(a_useCent16ForCent9){
          if(centIndex > 8) continue;
          if(centIndex == 0) centIndexToGrab = 0;
          if(centIndex == 1) centIndexToGrab = 1;
          if(centIndex == 2) centIndexToGrab = 2;
          if(centIndex == 3) centIndexToGrab = 4;
          if(centIndex == 4) centIndexToGrab = 6;
          if(centIndex == 5) centIndexToGrab = 8;
          if(centIndex == 6) centIndexToGrab = 10;
          if(centIndex == 7) centIndexToGrab = 12;
          if(centIndex == 8) centIndexToGrab = 14;
        }
        #ifdef _RAWSPECMOD_DEBUG_
          cout << "####  Working on TPC Eff / Energy Loss " << m_partInfo->GetParticleName(partIndex) << " Cent " << centIndex << "     grabbing Cent: " << centIndexToGrab << endl;
        #endif

        TGraphErrors* spectraGraphs[6] = {NULL,NULL,NULL,NULL,NULL,NULL};
        if(yIndex_tpc >= 0){
          spectraGraphs[0] = m_spectraPlusDeDx[partIndex][centIndex][yIndex_tpc];
          spectraGraphs[1] = m_spectraMinusDeDx[partIndex][centIndex][yIndex_tpc];          
        }
        if(yIndex_btof >= 0){
          spectraGraphs[2] = m_spectraPlusInvBeta[partIndex][centIndex][yIndex_btof];
          spectraGraphs[3] = m_spectraMinusInvBeta[partIndex][centIndex][yIndex_btof];
        }
        if(m_loadedETOF && yIndex_etof >= 0){
          spectraGraphs[4] = m_spectraPlusInvBetaETOF[partIndex][centIndex][yIndex_etof];
          spectraGraphs[5] = m_spectraMinusInvBetaETOF[partIndex][centIndex][yIndex_etof];
        }


        #ifdef _RAWSPECMOD_DEBUG_
          cout << "   yIndex  for     tpc: " << yIndex_tpc << "    bToF: " << yIndex_btof << "    eToF: " << yIndex_etof << "     emb: " << yIndex_embedding << endl;
          cout << "   Ptrs for tpc  : " << spectraGraphs[0] << " " << spectraGraphs[1] << endl;
          cout << "   Ptrs for bToF : " << spectraGraphs[2] << " " << spectraGraphs[3] << endl;
          cout << "   Ptrs for eToF : " << spectraGraphs[4] << " " << spectraGraphs[5] << endl;

          if(spectraGraphs[0]) cout << "     dEdx graph plus: "  << spectraGraphs[0]->GetName() << endl;
          if(spectraGraphs[1]) cout << "     dEdx graph minus: " << spectraGraphs[1]->GetName() << endl;
          if(spectraGraphs[2]) cout << "     bToF graph plus: "  << spectraGraphs[2]->GetName() << endl;
          if(spectraGraphs[3]) cout << "     bToF graph minus: " << spectraGraphs[3]->GetName() << endl;
          if(m_loadedETOF){
            if(spectraGraphs[4]) cout << "     eToF graph plus: "  << spectraGraphs[4]->GetName() << endl;
            if(spectraGraphs[5]) cout << "     bToF graph minus: " << spectraGraphs[5]->GetName() << endl;
          }
        #endif

        TF1*          shiftSpecFunct[2]       = {NULL,NULL};
        TGraphErrors* shiftConfSpecFunct[2]   = {NULL,NULL};
        TMatrixDSym*  shiftCovMatSpecFunct[2] = {NULL,NULL};
        if(!a_useEnergyLossSystematic){
          shiftSpecFunct[0]       = energy_loss_nom_plus_fit;
          shiftSpecFunct[1]       = energy_loss_nom_minus_fit;
          shiftConfSpecFunct[0]   = energy_loss_nom_plus_fit_conf;
          shiftConfSpecFunct[1]   = energy_loss_nom_minus_fit_conf;
          shiftCovMatSpecFunct[0] = energy_loss_nom_plus_fit_covMat;
          shiftCovMatSpecFunct[1] = energy_loss_nom_minus_fit_covMat;
        }else{
          shiftSpecFunct[0]       = energy_loss_sys_plus_fit;
          shiftSpecFunct[1]       = energy_loss_sys_minus_fit;
          shiftConfSpecFunct[0]   = energy_loss_sys_plus_fit_conf;
          shiftConfSpecFunct[1]   = energy_loss_sys_minus_fit_conf;
          shiftCovMatSpecFunct[0] = energy_loss_sys_plus_fit_covMat;
          shiftCovMatSpecFunct[1] = energy_loss_sys_minus_fit_covMat;
        }




        //####  READ IN TPC Efficiency ######
        TGraphAsymmErrors* tpc_eff_plus = (TGraphAsymmErrors*) inFile->Get( Form("%s/EfficiencyGraphs/tpcEfficiencyGraph_%s_Cent%02d_yIndex%02d", 
            m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_embedding)   );
        TGraphAsymmErrors* tpc_eff_minus = (TGraphAsymmErrors*) inFile->Get(  Form("%s/EfficiencyGraphs/tpcEfficiencyGraph_%s_Cent%02d_yIndex%02d", 
            m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_embedding)   );

        TF1* tpc_eff_nom_plus_fit = (TF1*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_embedding)   );
        TF1* tpc_eff_sys_plus_fit = (TF1*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_embedding)  );
        TF1* tpc_eff_nom_minus_fit = (TF1*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_embedding) );
        TF1* tpc_eff_sys_minus_fit = (TF1*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_embedding) );
        TF1* tpc_eff_nom_plus_fit_highPt = (TF1*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_HighPt_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_embedding)   );
        TF1* tpc_eff_nom_minus_fit_highPt = (TF1*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_HighPt_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_embedding) );

        TGraphErrors* tpc_eff_nom_plus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Conf", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_embedding)   );
        TGraphErrors* tpc_eff_sys_plus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys_Conf", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_embedding)  );
        TGraphErrors* tpc_eff_nom_minus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Conf", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_embedding) );
        TGraphErrors* tpc_eff_sys_minus_fit_conf = (TGraphErrors*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys_Conf", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_embedding) );

        TMatrixDSym* tpc_eff_nom_plus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyCovMat_%s_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_embedding));
        TMatrixDSym* tpc_eff_sys_plus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyCovMat_%s_Cent%d_yIndex%d_Sys",
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_embedding));
        TMatrixDSym* tpc_eff_nom_minus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyCovMat_%s_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_embedding));
        TMatrixDSym* tpc_eff_sys_minus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyCovMat_%s_Cent%d_yIndex%d_Sys",
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_embedding));
        TMatrixDSym* tpc_eff_nom_plus_fit_covMat_highPt = (TMatrixDSym*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyCovMat_%s_HighPt_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_embedding));
        TMatrixDSym* tpc_eff_nom_minus_fit_covMat_highPt = (TMatrixDSym*) inFile->Get( Form("%s/EfficiencyFits/tpcEfficiencyCovMat_%s_HighPt_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_embedding));

        //####  BTOF EFFICIENCY FIT DATA ######

        TGraphAsymmErrors* btof_eff_plus = (TGraphAsymmErrors*) inFile->Get(Form("%s/BTOFEfficiencyGraphs/btofFromDataEfficiencyGraph_%s_Cent%02d_yIndex%02d", 
            m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_btof)   );
        TGraphAsymmErrors* btof_eff_minus = (TGraphAsymmErrors*) inFile->Get( Form("%s/BTOFEfficiencyGraphs/btofFromDataEfficiencyGraph_%s_Cent%02d_yIndex%02d", 
            m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_btof)   );

        TF1* btof_eff_nom_plus_fit = (TF1*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_btof)   );
        TF1* btof_eff_sys_plus_fit = (TF1*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d_Sys", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_btof)  );
        TF1* btof_eff_nom_minus_fit = (TF1*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_btof) );
        TF1* btof_eff_sys_minus_fit = (TF1*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d_Sys", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_btof) );
        TF1* btof_eff_nom_plus_fit_highPt = (TF1*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyFit_DataDriven_%s_HighPt_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexToGrab,yIndex_btof)   );
        TF1* btof_eff_nom_minus_fit_highPt = (TF1*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyFit_DataDriven_%s_HighPt_Cent%d_yIndex%d", 
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexToGrab,yIndex_btof) );

        TMatrixDSym* btof_eff_nom_plus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyCovMat_DataDriven_%s_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_btof));
        TMatrixDSym* btof_eff_sys_plus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyCovMat_DataDriven_%s_Cent%d_yIndex%d_Sys",
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_btof));
        TMatrixDSym* btof_eff_nom_minus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyCovMat_DataDriven_%s_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_btof));
        TMatrixDSym* btof_eff_sys_minus_fit_covMat = (TMatrixDSym*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyCovMat_DataDriven_%s_Cent%d_yIndex%d_Sys",
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_btof));
        TMatrixDSym* btof_eff_nom_plus_fit_covMat_highPt = (TMatrixDSym*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyCovMat_DataDriven_HighPt_%s_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_btof));
        TMatrixDSym* btof_eff_nom_minus_fit_covMat_highPt = (TMatrixDSym*) inFile->Get( Form("%s/BTOFEfficiencyFits/btofEfficiencyCovMat_DataDriven_HighPt_%s_Cent%d_yIndex%d",
          m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_btof));



        #ifdef _RAWSPECMOD_DEBUG_
          cout << " TPC Eff TF1 Ptrs: " << tpc_eff_nom_plus_fit << " " << tpc_eff_sys_plus_fit << " " << tpc_eff_nom_minus_fit << " " << tpc_eff_sys_minus_fit << endl;
          if(m_chrisCorrectionFileMode){
            cout << "      Ptrs Conf: " << tpc_eff_nom_plus_fit_conf << " " << tpc_eff_sys_plus_fit_conf << " " << tpc_eff_nom_minus_fit_conf << " " << tpc_eff_sys_minus_fit_conf << endl;        
          }else{
            cout << "      Ptrs CovMat: " << tpc_eff_nom_plus_fit_covMat << " " << tpc_eff_sys_plus_fit_covMat << " " << tpc_eff_nom_minus_fit_covMat << " " << tpc_eff_sys_minus_fit_covMat << endl;        
          }

          cout << "################  EFFICIENCY PTR CHECK:" << endl;
          cout << "\ttpc_eff_plus: " << tpc_eff_plus << endl;
          cout << "\ttpc_eff_minus: " << tpc_eff_minus << endl;
          cout << "\ttpc_eff_nom_plus_fit: " << tpc_eff_nom_plus_fit << endl;
          cout << "\ttpc_eff_sys_plus_fit: " << tpc_eff_sys_plus_fit << endl;
          cout << "\ttpc_eff_nom_minus_fit: " << tpc_eff_nom_minus_fit << endl;
          cout << "\ttpc_eff_sys_minus_fit: " << tpc_eff_sys_minus_fit << endl;
          cout << "\ttpc_eff_nom_plus_fit_highPt: " << tpc_eff_nom_plus_fit_highPt << endl;
          cout << "\ttpc_eff_nom_minus_fit_highPt: " << tpc_eff_nom_minus_fit_highPt << endl;

          cout << "\ttpc_eff_nom_plus_fit_conf: " << tpc_eff_nom_plus_fit_conf << endl;
          cout << "\ttpc_eff_sys_plus_fit_conf: " << tpc_eff_sys_plus_fit_conf << endl;
          cout << "\ttpc_eff_nom_minus_fit_conf: " << tpc_eff_nom_minus_fit_conf << endl;
          cout << "\ttpc_eff_sys_minus_fit_conf: " << tpc_eff_sys_minus_fit_conf << endl;

          cout << "\ttpc_eff_nom_plus_fit_covMat: " << tpc_eff_nom_plus_fit_covMat << endl;
          cout << "\ttpc_eff_sys_plus_fit_covMat: " << tpc_eff_sys_plus_fit_covMat << endl;
          cout << "\ttpc_eff_nom_minus_fit_covMat: " << tpc_eff_nom_minus_fit_covMat << endl;
          cout << "\ttpc_eff_sys_minus_fit_covMat: " << tpc_eff_sys_minus_fit_covMat << endl;
          cout << "\ttpc_eff_nom_plus_fit_covMat_highPt: " << tpc_eff_nom_plus_fit_covMat_highPt << endl;
          cout << "\ttpc_eff_nom_minus_fit_covMat_highPt: " << tpc_eff_nom_minus_fit_covMat_highPt <<endl << endl;

          cout << "\tbtof_eff_plus: "  << btof_eff_plus << endl;
          cout << "\tbtof_eff_minus: " << btof_eff_minus << endl;

          cout << "\tbtof_eff_nom_plus_fit: " << btof_eff_nom_plus_fit << endl;
          cout << "\tbtof_eff_sys_plus_fit: " << btof_eff_sys_plus_fit << endl;
          cout << "\tbtof_eff_nom_minus_fit: " << btof_eff_nom_minus_fit << endl;
          cout << "\tbtof_eff_sys_minus_fit: " << btof_eff_sys_minus_fit << endl;
          cout << "\tbtof_eff_nom_plus_fit_highPt: " << btof_eff_nom_plus_fit_highPt << endl;
          cout << "\tbtof_eff_nom_minus_fit_highPt: " << btof_eff_nom_minus_fit_highPt << endl;

          cout << "\tbtof_eff_nom_plus_fit_covMat: " << btof_eff_nom_plus_fit_covMat << endl;
          cout << "\tbtof_eff_sys_plus_fit_covMat: " << btof_eff_sys_plus_fit_covMat << endl;
          cout << "\tbtof_eff_nom_minus_fit_covMat: " << btof_eff_nom_minus_fit_covMat << endl;
          cout << "\tbtof_eff_sys_minus_fit_covMat: " << btof_eff_sys_minus_fit_covMat << endl;
          cout << "\tbtof_eff_nom_plus_fit_covMat_highPt: " << btof_eff_nom_plus_fit_covMat_highPt << endl;
          cout << "\tbtof_eff_nom_minus_fit_covMat_highPt: " << btof_eff_nom_minus_fit_covMat_highPt << endl;

        #endif





        TF1* tpcEffFunct[2] = {NULL, NULL};;
        TGraphAsymmErrors* tpcEffGr[2] = {NULL, NULL};
        tpcEffGr[0] = tpc_eff_plus;
        tpcEffGr[1] = tpc_eff_minus;
        TGraphErrors* tpcEffConf[2] = {NULL, NULL};
        TMatrixDSym* tpcEffCovMat[2] = {NULL, NULL};
        TF1* tpcEffFunct_HighPt[2] = {NULL, NULL};
        TMatrixDSym* tpcEffCovMat_HighPt[2] = {NULL, NULL};
        TGraphAsymmErrors* btofEffGr[2] = {NULL, NULL};
        btofEffGr[0] = btof_eff_plus;
        btofEffGr[1] = btof_eff_minus;
        TF1* btofEffFunct[2] = {NULL, NULL};;
        TMatrixDSym* btofEffCovMat[2] = {NULL, NULL};
        TF1* btofEffFunct_HighPt[2] = {NULL, NULL};
        TMatrixDSym* btofEffCovMat_HighPt[2] = {NULL, NULL};
        if(!a_useTPCEffSystematic){
          tpcEffFunct[0] = tpc_eff_nom_plus_fit;
          tpcEffFunct[1] = tpc_eff_nom_minus_fit;
          tpcEffConf[0] = tpc_eff_nom_plus_fit_conf;
          tpcEffConf[1] = tpc_eff_nom_minus_fit_conf;
          tpcEffCovMat[0] = tpc_eff_nom_plus_fit_covMat;
          tpcEffCovMat[1] = tpc_eff_nom_minus_fit_covMat;
          tpcEffFunct_HighPt[0]  = tpc_eff_nom_plus_fit_highPt;
          tpcEffFunct_HighPt[1]  = tpc_eff_nom_minus_fit_highPt;
          tpcEffCovMat_HighPt[0] = tpc_eff_nom_plus_fit_covMat_highPt;
          tpcEffCovMat_HighPt[1] = tpc_eff_nom_minus_fit_covMat_highPt;
        }else{
          tpcEffFunct[0]  = tpc_eff_sys_plus_fit;
          tpcEffFunct[1]  = tpc_eff_sys_minus_fit;
          tpcEffConf[0]   = tpc_eff_sys_plus_fit_conf;
          tpcEffConf[1]   = tpc_eff_sys_minus_fit_conf;
          tpcEffCovMat[0] = tpc_eff_sys_plus_fit_covMat;
          tpcEffCovMat[1] = tpc_eff_sys_minus_fit_covMat;
        }
        if(!a_useBTOFEffSystematic){
          btofEffFunct[0]  = btof_eff_nom_plus_fit;
          btofEffFunct[1]  = btof_eff_nom_minus_fit;
          btofEffCovMat[0] = btof_eff_nom_plus_fit_covMat;
          btofEffCovMat[1] = btof_eff_nom_minus_fit_covMat;
          btofEffFunct_HighPt[0]  = btof_eff_nom_plus_fit_highPt;
          btofEffFunct_HighPt[1]  = btof_eff_nom_minus_fit_highPt;
          btofEffCovMat_HighPt[0] = btof_eff_nom_plus_fit_covMat_highPt;
          btofEffCovMat_HighPt[1] = btof_eff_nom_minus_fit_covMat_highPt;
        }else{
          btofEffFunct[0] = btof_eff_sys_plus_fit;
          btofEffFunct[1] = btof_eff_sys_minus_fit;
          btofEffCovMat[0] = btof_eff_sys_plus_fit_covMat;
          btofEffCovMat[1] = btof_eff_sys_minus_fit_covMat;
        }
        for(int detIndex = 0; detIndex < 3; detIndex++){
          string detNames[3] = {"TPC","BTOF","ETOF"};
          //if(detIndex == 1 && (yIndexModified_btof < 0 || yIndexModified_btof >= m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ){
          //  cout << "  WARNING: the yIndexModified_btof is out of the range : " << yIndexModified_btof << endl;
          //  continue;
          //}
          //if(detIndex == 2 && m_loadedETOF && (yIndexModified_etof < 0 || yIndexModified_etof >= m_InvBetaETOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ){
          //  cout << "   WARNING: the yIndexModified_etof is out of the range : " << yIndexModified_etof << endl;
          //  continue;
          //}

          for(int pmIndex = 0; pmIndex<2; pmIndex++){
            int charge = 1;
            if(pmIndex == 1) charge = -1;
            cout << "\n\n###################################################################################################################" << endl;
            cout << "  Now Modifying " << m_partInfo->GetParticleName(partIndex,charge) << "  Cent " << centIndex << " Detector: " << detNames[detIndex] << " y = " << rapidity << endl;
            cout << "###################################################################################################################" << endl;


            TGraphErrors* spectraGraph   = spectraGraphs[2*detIndex+pmIndex];
            TF1* shiftFunct              = shiftSpecFunct[pmIndex];
            TGraphErrors* shiftFunctConf = shiftConfSpecFunct[pmIndex];
            TMatrixDSym* shiftCovMat     = shiftCovMatSpecFunct[pmIndex];
            TF1* effFunct                = tpcEffFunct[pmIndex];
            TGraphAsymmErrors* effGr     = tpcEffGr[pmIndex];
            TGraphErrors* effGraphForcedErr = NULL;
            if(effGr) effGraphForcedErr = HistogramUtilities::convertAsymToNoAsym(effGr);
            TGraphErrors* effConf        = tpcEffConf[pmIndex];
            TMatrixDSym* effCovMat       = tpcEffCovMat[pmIndex];
            TF1* effFunct_HighPt                = tpcEffFunct_HighPt[pmIndex];
            TMatrixDSym* effCovMat_HighPt       = tpcEffCovMat_HighPt[pmIndex];
            SplineFitter* splineFit_energyLoss = NULL;
            SplineFitter* splineFit_tpcEff     = NULL;

            TGraphAsymmErrors* btofEffGr_spec     = btofEffGr[pmIndex];
            TGraphErrors* btofEffGraphForcedErr = NULL;
            if(btofEffGr_spec) btofEffGraphForcedErr = HistogramUtilities::convertAsymToNoAsym(btofEffGr_spec);
            TF1* btofEffFunct_spec                = btofEffFunct[pmIndex];
            TMatrixDSym* btofEffCovMat_spec       = btofEffCovMat[pmIndex];
            TF1* btofEffFunct_HighPt_spec                = btofEffFunct_HighPt[pmIndex];
            TMatrixDSym* btofEffCovMat_HighPt_spec       = btofEffCovMat_HighPt[pmIndex];

            TH2D* energyLossHistoToFill = NULL;
            if(a_useEnergyLossSystematic && pmIndex == 0)  energyLossHistoToFill = m_EnergyLoss_Systematic_Plus[partIndex];
            if(a_useEnergyLossSystematic && pmIndex == 1)  energyLossHistoToFill = m_EnergyLoss_Systematic_Minus[partIndex];
            if(!a_useEnergyLossSystematic && pmIndex == 0) energyLossHistoToFill = m_EnergyLoss_Nominal_Plus[partIndex];
            if(!a_useEnergyLossSystematic && pmIndex == 1) energyLossHistoToFill = m_EnergyLoss_Nominal_Minus[partIndex];
            if(!energyLossHistoToFill) cout << " ERROR: No energy loss histo when doing tpc efficiency integral" << endl;
            TH2D* tpcEffHistogramToFill = NULL;
            if(a_useTPCEffSystematic && pmIndex == 0)  tpcEffHistogramToFill = m_TPCEfficiency_Systematic_Plus[partIndex][centIndex];
            if(a_useTPCEffSystematic && pmIndex == 1)  tpcEffHistogramToFill = m_TPCEfficiency_Systematic_Minus[partIndex][centIndex];
            if(!a_useTPCEffSystematic && pmIndex == 0) tpcEffHistogramToFill = m_TPCEfficiency_Nominal_Plus[partIndex][centIndex];
            if(!a_useTPCEffSystematic && pmIndex == 1) tpcEffHistogramToFill = m_TPCEfficiency_Nominal_Minus[partIndex][centIndex];


            if(spectraGraph){
              cout << "\n\n######  Spectra to Modify: " << spectraGraph->GetName() << endl;
            }else{
              cout << "  WARNING: Spectra does not exist...  Part: " << partIndex << " pmIndex: " << pmIndex <<  " Cent: " << centIndex << " detIndex " << detIndex;
              if(detIndex == 0)      cout << "  yIndex: " << yIndex_tpc << endl;
              else if(detIndex == 1) cout << "  yIndex: " << yIndex_btof << endl;
              else if(detIndex == 2) cout << "  yIndex: " << yIndex_etof << endl;
              continue;
            }

            #ifdef _RAWSPECMOD_DEBUG_
              cout << "Low pT Covariance Matrix: " << endl;
              if(effCovMat) effCovMat->Print();
              cout << "High pT Covariance Matrix: " << endl;
              if(effCovMat_HighPt) effCovMat_HighPt->Print();
            #endif


            TH2D* referenceHistogram = NULL;
            if(detIndex == 0){
              if(pmIndex == 0) referenceHistogram = m_dEdxSpectraPlus[partIndex][centIndex];
              else             referenceHistogram = m_dEdxSpectraMinus[partIndex][centIndex];
            }else if(detIndex == 1){
              if(pmIndex == 0) referenceHistogram = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
              else             referenceHistogram = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];              
            }else{
              if(pmIndex == 0) referenceHistogram = m_InvBetaETOFSpectraPlus[partIndex][centIndex];
              else             referenceHistogram = m_InvBetaETOFSpectraMinus[partIndex][centIndex];                  
            }

            bool doEnergyLoss = true;
            if(!spectraGraph || (m_chrisCorrectionFileMode && !shiftFunctConf) || (!m_chrisCorrectionFileMode && !shiftCovMat)  || !shiftFunct){
            //if(partIndex == 0 || !spectraGraph || (m_chrisCorrectionFileMode && !shiftFunctConf) || !shiftFunct){
              cout << "  WARNING: Not applying pT Loss for Part: " << partIndex << " pmIndex: " << pmIndex <<  " Cent: " << centIndex << "  y: " << rapidity <<  " detIndex " << detIndex << endl;
              cout << "       Energy Loss ptrs: " << spectraGraph << "  "  << shiftFunct << " " << shiftCovMat << "  " << shiftFunctConf << "  need one of last two " << endl;
              doEnergyLoss = false;
            }

            bool doTPCEfficiency = true;
            if(!spectraGraph || (!a_useTPCEffDataPoints && !effFunct) || (a_useTPCEffDataPoints && !effGr) || (m_chrisCorrectionFileMode && !effConf) || (!m_chrisCorrectionFileMode && !a_useTPCEffDataPoints && !effCovMat)){
              cout << "  WARNING: Not applying TPC Eff for Part: " << partIndex << " pmIndex: " << pmIndex <<  "  Cent: " << centIndex << "  y: " << rapidity   << " detIndex " << detIndex  << " tot " << 2*detIndex+pmIndex << endl;
              cout << "        TPC Eff Ptrs: " << spectraGraph << "  " << effFunct << "  " << effGr <<  "  " << effCovMat << "  " << effConf << endl;// << " " << shiftFunctConf << endl;
              doTPCEfficiency = false;
            }

            bool doBTOFEfficiency = true;
            if(detIndex != 1 || !spectraGraph || (a_useBTOFEffDataPoints && (!btofEffGr_spec || btofEffGr_spec->GetN() < 2 )) || (!a_useBTOFEffDataPoints &&(!btofEffFunct_spec ||  !btofEffCovMat_spec))){
              if(detIndex == 1){
                cout << "  WARNING: Not applying BTOF Eff for Part: " << partIndex << " pmIndex: " << pmIndex <<  "  Cent: " << centIndex  << " y: " << rapidity << " detIndex " << detIndex  << " tot " << 2*detIndex+pmIndex << endl;
                cout << "        BTOF Eff Ptrs: " << spectraGraph << "  " << btofEffGr_spec << "   " << btofEffFunct_spec <<  "  " << btofEffCovMat_spec << "  "  << endl;// << " " << shiftFunctConf << endl;
              }
              doBTOFEfficiency = false;
            }

            // USER OVERRIDES
            if(!a_doEnergyLoss)  doEnergyLoss     = false;
            if(!a_doTPCEff)      doTPCEfficiency  = false;
            if(!a_doBTOFEff)     doBTOFEfficiency = false;

            #ifdef _RAWSPECMOD_DEBUG_
              if(doEnergyLoss){
                cout << "Applying Energy Loss for Part: " << partIndex << " pmIndex: " << pmIndex <<  " Cent: " << centIndex << " y: " << rapidity  << " detIndex " << detIndex << endl;
              }
              if(doTPCEfficiency){
                cout << "Applying TPC Eff for Part: " << partIndex << " pmIndex: " << pmIndex <<  " Cent: " << centIndex << " y: " << rapidity << " detIndex " << detIndex << endl;
              }
              if(doBTOFEfficiency){
                cout << "Applying BTOF Eff for Part: " << partIndex << " pmIndex: " << pmIndex <<  " Cent: " << centIndex << " y: " << rapidity << " detIndex " << detIndex << endl;
              }
            #endif




            if(!doEnergyLoss && !doTPCEfficiency && !doBTOFEfficiency){
              #ifdef _RAWSPECMOD_DEBUG_
                cout << "  WARNING: Not applying energy loss or tpc efficiency or btof efficiency to this particle" <<  endl;
              #endif
              continue;
            }
            //spline fit for tpc Eff Nominal
            if(doSplineFitting && !m_chrisCorrectionFileMode){
              if(!a_useEnergyLossSystematic){
                splineFit_energyLoss = new SplineFitter();
                #ifdef _RAWSPECMOD_DEBUG_
                  cout << " Using the spline fitter energy loss with ptr: " << splineFit_energyLoss << "  TF1: " << shiftFunct << "   covMat: " << shiftCovMat <<  endl;
                #endif
                splineFit_energyLoss->reload(shiftFunct, shiftCovMat);
                shiftFunct = splineFit_energyLoss->getFitFunciton();
              }
              if(!a_useTPCEffSystematic){
                splineFit_tpcEff = new SplineFitter();
                #ifdef _RAWSPECMOD_DEBUG_
                  cout << " Using the spline fitter tpc efficiency with ptr: " << splineFit_tpcEff << "  TF1: " << effFunct << "   covMat: " << effCovMat <<  endl;
                #endif
                splineFit_tpcEff->reload(effFunct, effCovMat);
                effFunct = splineFit_tpcEff->getFitFunciton();
              }
              //else{
                //don't do anything, the functions don't need to be reloaded
              //}
            }


            #ifdef _RAWSPECMOD_DEBUG_
              if(doEnergyLoss){
                cout << "     Working on detector index " << detIndex << "  pmIndex " << pmIndex << " tot " << 2*detIndex + pmIndex << endl;
                cout << "     ######### ENERGY LOSS TITLE CHECK ############" << endl;
                cout << "       ---------- " << spectraGraph->GetName() << endl;
                if(shiftFunct) cout << "       ---------- " << shiftFunct->GetName() << endl;
                cout << "       ---------- " << spectraGraph->GetTitle() << endl;
                if(shiftFunct) cout << "       ---------- " << shiftFunct->GetTitle() << endl;
                if(shiftFunctConf) cout << "       ---------- " << shiftFunctConf->GetName() << endl;
                if(m_chrisCorrectionFileMode) cout << "   CHRIS pT Loss Correction Mode is On!!!!" << endl;
              }
              if(doTPCEfficiency){
                cout << "     ######### TPC EFF TITLE CHECK ############" << endl;
                cout << "       ---------- " << spectraGraph->GetName() << endl;
                if(effFunct) cout << "       ---------- " << effFunct->GetName() << endl;
                if(effGr)   cout << "       ---------- " << effGr->GetName() << endl;
                if(effConf) cout << "       ---------- " << effConf->GetName() << endl;
                cout << "       ---------- " << spectraGraph->GetTitle() << endl;
                if(effFunct)cout << "       ---------- " << effFunct->GetTitle() << endl;
                if(effGr)   cout << "       ---------- " << effGr->GetTitle() << endl;
                if(effConf) cout << "       ---------- " << effConf->GetTitle() << endl;              
              }

              if(doBTOFEfficiency){
                cout << "     #########  BTOF EFF TITLE CHECK ############" << endl;
                cout << "       ---------- " << spectraGraph->GetName() << endl;
                if(btofEffFunct_spec) cout << "       ---------- " << btofEffFunct_spec->GetName() << endl;
                if(btofEffGr_spec)    cout << "       ---------- " << btofEffGr_spec->GetName() << endl;
                cout << "       ---------- " << spectraGraph->GetTitle() << endl;
                if(btofEffFunct_spec) cout << "       ---------- " << btofEffFunct_spec->GetTitle() << endl;
                if(btofEffGr_spec)    cout << "       ---------- " << btofEffGr_spec->GetTitle() << endl;
              }
            #endif


            bool didEnergyLoss   = false;
            bool didTPCEffiency  = false;
            bool didBTOFEffiency = false;

            SpectraStructure* spectraStruct = findAssociatedSpectraStructure(partIndex,centIndex, charge, rapidity);

            int numPointsActuallyApplied_EnergyLoss = 0;
            int numPointsActuallyApplied_TPCEff     = 0;
            int numPointsActuallyApplied_BTOFEff    = 0;

            //#############       LOOP OVER THE SPECTRA POINTS  ###############################
            for(int pointIndex = 0; pointIndex < spectraGraph->GetN(); pointIndex++){
              //Get the point info
              double mTm0  = spectraGraph->GetX()[pointIndex];
    
              //new embedding routine uses mT-m0 instead of pT
              double binOfGraphPoint = referenceHistogram->GetYaxis()->FindBin(mTm0);
              double mTm0Low  = referenceHistogram->GetYaxis()->GetBinLowEdge(binOfGraphPoint);
              double mTm0High = referenceHistogram->GetYaxis()->GetBinLowEdge(binOfGraphPoint+1);
              double corrected_mTm0     = mTm0;
              double corrected_mTm0Low  = mTm0Low;
              double corrected_mTm0High = mTm0High;
  

              //####################   ENERGY LOSS CORRECTION #############################
              if(doEnergyLoss){
                if(m_chrisCorrectionFileMode){
                  double pT_ofPoint = PhysMath::pT_from_mTm0(mTm0,mass);
                  double lowPt  = PhysMath::pT_from_mTm0(mTm0 - spectraGraph->GetEX()[pointIndex],mass);
                  double highPt = PhysMath::pT_from_mTm0(mTm0 + spectraGraph->GetEX()[pointIndex],mass);
                  double pTLoss = 0.0; //reco - emb
                  double pTLossError = 0.0; 
                  HistogramUtilities::evalWithError(shiftFunctConf, pT_ofPoint,pTLoss, pTLossError);
                  pTLoss = shiftFunct->Integral(lowPt,highPt)/(highPt - lowPt);
                  corrected_mTm0     = PhysMath::mTm0(mass,pT_ofPoint - pTLoss); // becasue reco - emb in fitting
                  double corrected_mTm0ErrLow  = PhysMath::mTm0(mass,pT_ofPoint - pTLossError - shiftFunct->Eval(pT_ofPoint - pTLossError));
                  double corrected_mTm0ErrHigh = PhysMath::mTm0(mass,pT_ofPoint + pTLossError - shiftFunct->Eval(pT_ofPoint + pTLossError));
                  corrected_mTm0Low  = PhysMath::mTm0(mass,lowPt - shiftFunct->Eval(lowPt));
                  corrected_mTm0High = PhysMath::mTm0(mass,highPt - shiftFunct->Eval(highPt));
                  double newXErr = 0.0;
                  if(corrected_mTm0ErrHigh - corrected_mTm0 > corrected_mTm0 - corrected_mTm0ErrLow){
                    newXErr = corrected_mTm0ErrHigh - corrected_mTm0;
                  }else{
                    newXErr = corrected_mTm0 - corrected_mTm0ErrLow;
                  }
                  double scale = (mTm0 + mass)/(corrected_mTm0 + mass); // becasue of 1/mT in the invariant yield
                  spectraGraph->SetPoint(pointIndex, corrected_mTm0, scale*spectraGraph->GetY()[pointIndex]);
                  spectraGraph->SetPointError(pointIndex,fabs(newXErr),scale*spectraGraph->GetEY()[pointIndex]);
                  energyLossHistoToFill->SetBinContent(energyLossHistoToFill->GetXaxis()->FindBin(rapidity),binOfGraphPoint,corrected_mTm0 - mTm0);
                  numPointsActuallyApplied_EnergyLoss++;
                  didEnergyLoss = true;
                }else{
                  double shift = 0; 
                  double shiftErr = 0;
                  if(m_pTLossMode){
                    double pTshift = 0; 
                    double pTshiftErr = 0;
                    double pT_ofPoint = PhysMath::pT_from_mTm0(mTm0,mass);
                    double pTLow  = PhysMath::pT_from_mTm0(mTm0Low,mass);
                    double pTHigh = PhysMath::pT_from_mTm0(mTm0High,mass);
                    if(a_useEnergyLossSystematic){
                      pTshift    = shiftFunct->Integral(pTLow,pTHigh)/(pTHigh - pTLow);
                      pTshiftErr = shiftFunct->IntegralError(pTLow,pTHigh, shiftFunct->GetParameters(), shiftCovMat->GetMatrixArray())/(pTHigh - pTLow);
                    }else{
                      //Use the spline fitter
                      if(doSplineFitting){
                        splineFit_energyLoss->getIntegralAndError(pTLow, pTHigh, pTshift, pTshiftErr, true);// bool = divide integral by bin width
                      }else{
                        pTshift    = shiftFunct->Integral(pTLow,pTHigh)/(pTHigh - pTLow);
                        pTshiftErr = shiftFunct->IntegralError(pTLow,pTHigh, shiftFunct->GetParameters(), shiftCovMat->GetMatrixArray())/(pTHigh - pTLow);
                      }
                    }
                    shift = PhysMath::mTm0(mass,pT_ofPoint + pTshift) - PhysMath::mTm0(mass,pT_ofPoint) ;
                    shiftErr = shift*pTshiftErr/pTshift;
                  
                    //spectraStruct->addEnergyLossCorrectionPoint(mTm0,shift, shiftErr, true); //bool is integrated value
                    //spectraStruct->addEnergyLossCorrectionPoint(mTm0,shiftFunct->Eval(mTm0), 0.0, true); //bool is integrated value

                    corrected_mTm0 = mTm0 + shift; // because we switched to embed - reco, the shift is added (Points move to higher mTm0)
                    corrected_mTm0Low  = PhysMath::mTm0(mass,pTLow + shiftFunct->Eval(pTLow));
                    corrected_mTm0High = PhysMath::mTm0(mass,pTHigh + shiftFunct->Eval(pTHigh));
                    bool doBinRecentering = true;
                    double scale = 1.0;
                    if(doBinRecentering){
                      if(a_useEnergyLossErrorAsXErrorBars){
                        scale = (mTm0 + mass)/(corrected_mTm0 + mass); // becasue of 1/mT in the invariant yield
                        spectraGraph->SetPoint(pointIndex, corrected_mTm0, scale*spectraGraph->GetY()[pointIndex]);
                        spectraGraph->SetPointError(pointIndex, shiftErr, scale*spectraGraph->GetEY()[pointIndex]);
                      }else{
                        corrected_mTm0 = 0.5*(corrected_mTm0High + corrected_mTm0Low);
                        scale = (mTm0 + mass)/(corrected_mTm0 + mass); // becasue of 1/mT in the invariant yield
                        spectraGraph->SetPoint(pointIndex, corrected_mTm0, scale*spectraGraph->GetY()[pointIndex]);
                        spectraGraph->SetPointError(pointIndex, 0.5*(corrected_mTm0High - corrected_mTm0Low), scale*spectraGraph->GetEY()[pointIndex]);
                      }
                      numPointsActuallyApplied_EnergyLoss++;
                      didEnergyLoss = true;
                      #ifdef _RAWSPECMOD_DEBUG_
                        cout << "   mTm0: " << mTm0 << " -> " << corrected_mTm0 << " scale: " << scale << " was [" << mTm0Low << "," << mTm0High << "]   now [" << corrected_mTm0Low << "," << corrected_mTm0High << "]" << endl;
                      #endif
                    }else{
                      scale = (mTm0 + mass)/(corrected_mTm0 + mass); // becasue of 1/mT in the invariant yield
                      spectraGraph->SetPoint(pointIndex, corrected_mTm0, scale*spectraGraph->GetY()[pointIndex]);
                      if(a_useEnergyLossErrorAsXErrorBars) spectraGraph->SetPointError(pointIndex, shiftErr, scale*spectraGraph->GetEY()[pointIndex]);
                      else spectraGraph->SetPointError(pointIndex, corrected_mTm0 - corrected_mTm0Low, scale*spectraGraph->GetEY()[pointIndex]);
                      energyLossHistoToFill->SetBinContent(energyLossHistoToFill->GetXaxis()->FindBin(rapidity),binOfGraphPoint,shift);
                      numPointsActuallyApplied_EnergyLoss++;
                      didEnergyLoss = true;
                      #ifdef _RAWSPECMOD_DEBUG_
                        cout << " EnergyLoss mTm0: " << mTm0 << " corrected to " << corrected_mTm0 << " by integral from " << mTm0Low << " to " << mTm0High << " with an extra scale " << scale << endl;
                        cout << "                     now bin from " << corrected_mTm0Low << " to " << corrected_mTm0High << endl;
                      #endif
                    }

                  }else{ // mTm0 mode
                    if(a_useEnergyLossSystematic){
                      shift    = shiftFunct->Integral(mTm0Low,mTm0High)/(mTm0High - mTm0Low);
                      shiftErr = shiftFunct->IntegralError(mTm0Low,mTm0High, shiftFunct->GetParameters(), shiftCovMat->GetMatrixArray())/(mTm0High - mTm0Low);
                    }else{
                      //Use the spline fitter
                      if(doSplineFitting){
                        splineFit_energyLoss->getIntegralAndError(mTm0Low, mTm0High, shift, shiftErr, true);// bool = divide integral by bin width
                      }else{
                        shift    = shiftFunct->Integral(mTm0Low,mTm0High)/(mTm0High - mTm0Low);
                        shiftErr = shiftFunct->IntegralError(mTm0Low,mTm0High, shiftFunct->GetParameters(), shiftCovMat->GetMatrixArray())/(mTm0High - mTm0Low);
                      }
                    }
                  
                    //spectraStruct->addEnergyLossCorrectionPoint(mTm0,shift, shiftErr, true); //bool is integrated value
                    //spectraStruct->addEnergyLossCorrectionPoint(mTm0,shiftFunct->Eval(mTm0), 0.0, true); //bool is integrated value

                    corrected_mTm0 = mTm0 + shift; // because we switched to embed - reco, the shift is added (Points move to higher mTm0)
                    corrected_mTm0Low = mTm0Low + shiftFunct->Eval(mTm0Low);
                    corrected_mTm0High = mTm0High + shiftFunct->Eval(mTm0High);
                    double scale = (mTm0 + mass)/(corrected_mTm0 + mass); // becasue of 1/mT in the invariant yield
                    spectraGraph->SetPoint(pointIndex, corrected_mTm0, scale*spectraGraph->GetY()[pointIndex]);
                    spectraGraph->SetPointError(pointIndex, shiftErr, scale*spectraGraph->GetEY()[pointIndex]);
                    energyLossHistoToFill->SetBinContent(energyLossHistoToFill->GetXaxis()->FindBin(rapidity),binOfGraphPoint,shift);
                    numPointsActuallyApplied_EnergyLoss++;
                    didEnergyLoss = true;
                    #ifdef _RAWSPECMOD_DEBUG_
                      cout << " EnergyLoss mTm0: " << mTm0 << " corrected to " << corrected_mTm0 << " by integral from " << mTm0Low << " to " << mTm0High << " with an extra scale " << scale << endl;
                      cout << "                     now bin from " << corrected_mTm0Low << " to " << corrected_mTm0High << endl;
                    #endif
                  }


                }
              }// if doEnergyLoss



              //####################   TPC CORRECTION #############################
              if(doTPCEfficiency){
                if(m_chrisCorrectionFileMode){// #######################################    LEGACY METHOD #################################################
                  double eff = 0.0;
                  double effErr = 0.0;
                  if(a_useTPCEffDataPoints){
                    HistogramUtilities::evalWithError(effGr,corrected_mTm0,eff,effErr);
                  }else{
                    //to set the effErr
                    HistogramUtilities::evalWithError(effConf,corrected_mTm0,eff,effErr); // okay aproximation to integral error (difficult to do without covariance matrix)
                    //to set the tpcEfficiency
                    eff = effFunct->Integral(corrected_mTm0Low,corrected_mTm0High)/(corrected_mTm0High-corrected_mTm0Low);
                  }
                  if( eff > 1.0 ){
                    eff = 1.0;
                    cout << "  WARNING: TPC Eff > 1 (" << eff << " > 1)   so setting to 1.0 " << endl;
                  }
                  else if( eff <= 0.0 ){ //removing point
                    spectraGraph->RemovePoint(pointIndex);
                    pointIndex--;
                    cout << "  WARNING: TPC Eff < 0 (" << eff << " < 0)   so removing point " << endl;
                    continue;
                  }

                  //access the tgraph points because yield has now changed
                  double newContent = spectraGraph->GetY()[pointIndex]/eff;
                  double newContentError = newContent*
                      TMath::Sqrt((spectraGraph->GetEY()[pointIndex]*spectraGraph->GetEY()[pointIndex])/(spectraGraph->GetY()[pointIndex]*spectraGraph->GetY()[pointIndex])
                                + (effErr*effErr)/(eff*eff));
                  spectraGraph->SetPoint(pointIndex, corrected_mTm0, newContent);
                  spectraGraph->SetPointError(pointIndex,spectraGraph->GetEX()[pointIndex],newContentError);
                  tpcEffHistogramToFill->SetBinContent(tpcEffHistogramToFill->GetXaxis()->FindBin(rapidity),binOfGraphPoint,eff);
                  tpcEffHistogramToFill->SetBinError(tpcEffHistogramToFill->GetXaxis()->FindBin(rapidity),binOfGraphPoint,effErr);
                  didTPCEffiency = true;
                  numPointsActuallyApplied_TPCEff++;

                }else{ // Not old chris data  // #########  TPC EFF  BETTER METHOD ##############
                  double eff = 0.0;
                  double effErr = 0.0;
                  
                  if(a_useTPCEffDataPoints){
                    HistogramUtilities::integrateGraph(effGraphForcedErr, corrected_mTm0Low, corrected_mTm0High, eff, effErr, true);
                    double pointEval = 0.0;
                    double pointEvalErr = 0.0;
                    HistogramUtilities::evalWithError(effGr,corrected_mTm0,pointEval,pointEvalErr);
                    #ifdef _RAWSPECMOD_DEBUG_
                      cout << "  Compare TPC Eff by Integral: " << eff << " +/- " << effErr << "     and by Point Eval: " << pointEval  << " +/- " << pointEvalErr << endl;
                    #endif
                  }else{
                    //DETERMINE IF THE POINT USES THE LOW PT OR HIGH PT FIT
                    TF1*         functToUseEff  = NULL;
                    TMatrixDSym* matrixToUseEff = NULL;
                    double lowPtFit_LowRange = 0.0;
                    double lowPtFit_highRange = 0.0;
                    effFunct->GetRange(lowPtFit_LowRange,lowPtFit_highRange);
                    if(a_useTPCEffSystematic){  // Systematic function includes the full range of the fit, so no need to do separate things for different slices
                      functToUseEff  = effFunct;
                      matrixToUseEff = effCovMat;
                      eff    = functToUseEff->Integral(corrected_mTm0Low,corrected_mTm0High)/(corrected_mTm0High - corrected_mTm0Low);
                      effErr = functToUseEff->IntegralError(corrected_mTm0Low,corrected_mTm0High, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/(corrected_mTm0High - corrected_mTm0Low);
                      #ifdef _RAWSPECMOD_DEBUG_
                        double pointEval = 0.0;
                        double pointEvalErr = 0.0;
                        HistogramUtilities::evalWithError(effGr,corrected_mTm0,pointEval,pointEvalErr);
                        cout << "      Eff Integral Method: " << eff << " +/- " << effErr << "   Point Method: " << pointEval << " +/- " << pointEvalErr << endl;
                        //effErr = pointEvalErr;
                      #endif
                      
                    }else{
                      if(corrected_mTm0 < lowPtFit_LowRange){
                        //Just Calculate Eff By The Data Points because the fit doesn't extend to low pT
                        cout << "graph ptr: " << effGraphForcedErr << endl;
                        HistogramUtilities::integrateGraph(effGraphForcedErr, corrected_mTm0Low, corrected_mTm0High, eff, effErr, true);
                        double pointEval = 0.0;
                        double pointEvalErr = 0.0;
                        HistogramUtilities::evalWithError(effGr,corrected_mTm0,pointEval,pointEvalErr);
                        cout << "  Compare Integral: " << eff << " +/- " << effErr << "       Point Eval: " << pointEval  << " +/- " << pointEvalErr << endl;

                      }else if(corrected_mTm0 < lowPtFit_highRange){
                        functToUseEff = effFunct;
                        matrixToUseEff = effCovMat;
                        eff    = functToUseEff->Integral(corrected_mTm0Low,corrected_mTm0High)/(corrected_mTm0High - corrected_mTm0Low);
                        effErr = functToUseEff->IntegralError(corrected_mTm0Low,corrected_mTm0High, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/(corrected_mTm0High - corrected_mTm0Low);
                      }else if(effFunct_HighPt){
                        functToUseEff = effFunct_HighPt;
                        matrixToUseEff = effCovMat_HighPt;
                        eff    = functToUseEff->Integral(corrected_mTm0Low,corrected_mTm0High)/(corrected_mTm0High - corrected_mTm0Low);
                        effErr = functToUseEff->IntegralError(corrected_mTm0Low,corrected_mTm0High, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/(corrected_mTm0High - corrected_mTm0Low);
                      }else{
                        //just do a eff of the last little bit of the fit (when the line doesn't exist)
                        functToUseEff = effFunct;
                        matrixToUseEff = effCovMat;
                        eff    = functToUseEff->Integral(lowPtFit_highRange-0.1,lowPtFit_highRange)/0.1;
                        effErr = functToUseEff->IntegralError(lowPtFit_highRange-0.1,lowPtFit_highRange, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/0.1;
                      }                    
                    }
                  }

                  //if(doSplineFitting){
                  //  splineFit_tpcEff->getIntegralAndError(corrected_mTm0Low, corrected_mTm0High, eff, effErr, true);// bool = divide integral by bin width
                  //}

                  if( eff > 1.0 ){
                    eff = 1.0;
                    cout << "  WARNING: TPC Eff > 1 (" << eff << " > 1)   so setting to 1.0 " << endl;
                  }
                  else if( eff <= 0.0 ){ //removing point
                    spectraGraph->RemovePoint(pointIndex);
                    pointIndex--;
                    cout << "  WARNING: TPC Eff < 0 (" << eff << " < 0)   so removing point " << endl;
                    continue;
                  }

                  //access the tgraph points because yield has now changed
                  double newContent = spectraGraph->GetY()[pointIndex]/eff;
                  double newContentError = newContent*
                      TMath::Sqrt((spectraGraph->GetEY()[pointIndex]*spectraGraph->GetEY()[pointIndex])/(spectraGraph->GetY()[pointIndex]*spectraGraph->GetY()[pointIndex])
                                + (effErr*effErr)/(eff*eff));
                  #ifdef _RAWSPECMOD_DEBUG_
                    cout << "       mT-m0: " << corrected_mTm0 << " N: " << spectraGraph->GetY()[pointIndex] << " +/- " << spectraGraph->GetEY()[pointIndex] << " -> " <<  newContent << " +/- " << newContentError << " by eff: " << eff << " +/- " << effErr <<  endl;
                  #endif
                  spectraGraph->SetPoint(pointIndex, spectraGraph->GetX()[pointIndex], newContent);
                  spectraGraph->SetPointError(pointIndex,spectraGraph->GetEX()[pointIndex],newContentError);
                  tpcEffHistogramToFill->SetBinContent(tpcEffHistogramToFill->GetXaxis()->FindBin(rapidity),binOfGraphPoint,eff);
                  tpcEffHistogramToFill->SetBinError(tpcEffHistogramToFill->GetXaxis()->FindBin(rapidity),binOfGraphPoint,effErr);
                  didTPCEffiency = true;
                  numPointsActuallyApplied_TPCEff++;
                }// Chris Correciton Mode
              }// if doTPCEfficiency


            //#######################     BTOF EFFICIENCY ########################################
            if(doBTOFEfficiency){ // THIS IS DATA DRIVEN, DON"T USE CORRECTED mT-m0!
              double eff = 0.0;
              double effErr = 0.0;

              if(a_useBTOFEffDataPoints){
                HistogramUtilities::integrateGraph(btofEffGraphForcedErr, mTm0Low, mTm0High, eff, effErr, true);
                double pointEval = 0.0;
                double pointEvalErr = 0.0;
                HistogramUtilities::evalWithError(btofEffGr_spec,mTm0,pointEval,pointEvalErr);
                #ifdef _RAWSPECMOD_DEBUG_
                  cout << "  Compare BTOF TGraph Integral : " << eff << " +/- " << effErr << "  vs Point Eval:  " << pointEval  << " +/- " << pointEvalErr << endl;
                #endif
              }else{
                //DETERMINE IF THE POINT USES THE LOW PT OR HIGH PT FIT
                //This is done with data, so reco not emb pT. This uses the original bounds of the mTm0 bin.
                TF1* functToUseEff = NULL;
                TMatrixDSym* matrixToUseEff = NULL;
                double lowPtFit_LowRange  = 0.0;
                double lowPtFit_highRange = 0.0;
                btofEffFunct_spec->GetRange(lowPtFit_LowRange,lowPtFit_highRange);
                string methodString = "";
                if(a_useBTOFEffSystematic){
                  if(mTm0 < lowPtFit_LowRange){
                    //HistogramUtilities::integrateGraph(btofEffGraphForcedErr, mTm0Low, mTm0High, eff, effErr, true);
                    //double pointEval = 0.0;
                    double pointEvalErr = 0.0;
                    HistogramUtilities::evalWithError(btofEffGr_spec,mTm0,eff,effErr);
                  }else{
                    functToUseEff  = btofEffFunct_spec;
                    matrixToUseEff = btofEffCovMat_spec;
                    eff    = functToUseEff->Integral(mTm0Low,mTm0High)/(mTm0High - mTm0Low);
                    effErr = functToUseEff->IntegralError(mTm0Low,mTm0High, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/(mTm0High - mTm0Low);
                    methodString = "sysetmatic function integral";
                  }
                }else{
                  if(mTm0 < lowPtFit_LowRange){
                    //Just Calculate Eff By The Data Points because the fit doesn't extend to low pT
                    HistogramUtilities::integrateGraph(btofEffGraphForcedErr, mTm0Low, mTm0High, eff, effErr, true);
                    double pointEval = 0.0;
                    double pointEvalErr = 0.0;
                    HistogramUtilities::evalWithError(btofEffGr_spec,mTm0,pointEval,pointEvalErr);
                    cout << "  Compare : " << eff << " +/- " << effErr << "   " << pointEval  << " +/- " << pointEvalErr << endl;
                    methodString = "nominal point eval";

                  }else if(mTm0 < lowPtFit_highRange){
                    functToUseEff = btofEffFunct_spec;
                    matrixToUseEff = btofEffCovMat_spec;
                    eff    = functToUseEff->Integral(mTm0Low,mTm0High)/(mTm0High - mTm0Low);
                    effErr = functToUseEff->IntegralError(mTm0Low,mTm0High, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/(mTm0High - mTm0Low);
                    methodString = "nominal legendre function integral";
                    if(!std::isfinite(effErr) || !std::isfinite(eff)){
                      cout << " ERROR: legendre funciton produce a horrible nan inf...    Eff: " << eff << " +/- " << effErr << endl;
                      cout << "     mTm0 = [" << mTm0Low << "," << mTm0High << "]     Point Funct Eval: " << functToUseEff->Eval(mTm0) << endl;
                      cout << " Function: " << endl;
                      functToUseEff->Print();
                      cout << " Covariance Matrix: " << endl;
                      matrixToUseEff->Print();

                    }
                  }else if(btofEffFunct_HighPt_spec){
                    functToUseEff = btofEffFunct_HighPt_spec;
                    matrixToUseEff = btofEffCovMat_HighPt_spec;
                    eff    = functToUseEff->Integral(mTm0Low,mTm0High)/(mTm0High - mTm0Low);
                    effErr = functToUseEff->IntegralError(mTm0Low,mTm0High, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/(mTm0High - mTm0Low);
                    methodString = "nominal linear function integral";
                  }else{
                    //just do a eff of the last little bit of the fit (when the line doesn't exist)
                    functToUseEff = btofEffFunct_spec;
                    matrixToUseEff = btofEffCovMat_spec;
                    eff    = functToUseEff->Integral(lowPtFit_highRange-0.1,lowPtFit_highRange)/0.1;
                    effErr = functToUseEff->IntegralError(lowPtFit_highRange-0.1,lowPtFit_highRange, functToUseEff->GetParameters(), matrixToUseEff->GetMatrixArray())/0.1;
                    methodString = "nominal legendre function integral even though high pT";
                  }
                }
                #ifdef _RAWSPECMOD_DEBUG_
                  cout << "  Given the ranges 0 -> " << lowPtFit_LowRange << " -> " << lowPtFit_highRange << "  with mTm0=" << mTm0 << "  and ptrs: " << btofEffFunct_spec << "   " << btofEffFunct_HighPt_spec << " chose " << functToUseEff << " by " << methodString << endl;
                #endif                
              }



              if(std::isnan(eff) || std::isnan(effErr)){
                cout << "ERROR: The BTOF EFF or its error is Nan : " << eff << "+/-" << effErr << endl;
                continue;
              }
              if( eff > 1.0 ){
                eff = 1.0;
                cout << "  WARNING: BTOF Eff > 1 (" << eff << " > 1)   so setting to 1.0 " << endl;
              }
              else if( eff <= 0.0 ){ //removing point
                spectraGraph->RemovePoint(pointIndex);
                pointIndex--;
                cout << "  WARNING: BTOF Eff < 0 (" << eff << " < 0)   so removing point mTm0:" << corrected_mTm0 << endl;
                continue;
              }

              //access the tgraph points because yield has now changed
              double newContent = spectraGraph->GetY()[pointIndex]/eff;
              double newContentError = newContent*
                  TMath::Sqrt((spectraGraph->GetEY()[pointIndex]*spectraGraph->GetEY()[pointIndex])/(spectraGraph->GetY()[pointIndex]*spectraGraph->GetY()[pointIndex])
                            + (effErr*effErr)/(eff*eff));
              #ifdef _RAWSPECMOD_DEBUG_
                cout << "  BTOF EFF Applied at mT-m0: " << corrected_mTm0 << " Yield: " << spectraGraph->GetY()[pointIndex] << " +/- " << spectraGraph->GetEY()[pointIndex] << " --> " << newContent << " +/- " << newContentError << "  by BTOF EFF: " << eff <<  " +/- " << effErr << endl;
              #endif
              spectraGraph->SetPoint(pointIndex, spectraGraph->GetX()[pointIndex], newContent);
              spectraGraph->SetPointError(pointIndex,spectraGraph->GetEX()[pointIndex],newContentError);
              didBTOFEffiency = true;
              numPointsActuallyApplied_BTOFEff++;
            }// if doBTOFEfficiency




            }// mTm0 point loop

            double fracApplied_EnergyLoss = ((double) numPointsActuallyApplied_EnergyLoss)/((double)spectraGraph->GetN());
            double fracApplied_TPCEff = ((double) numPointsActuallyApplied_TPCEff)/((double)spectraGraph->GetN());
            double fracApplied_BTOFEff = ((double) numPointsActuallyApplied_BTOFEff)/((double)spectraGraph->GetN());

            string correctionCode = " Corr:";
            if(didEnergyLoss && fracApplied_EnergyLoss > 0.5)   correctionCode = correctionCode + "E";
            if(didTPCEffiency && fracApplied_TPCEff > 0.5)      correctionCode = correctionCode + "T";
            if(didBTOFEffiency && fracApplied_BTOFEff > 0.5)    correctionCode = correctionCode + "B";
            string newTitle = spectraGraph->GetTitle();
            newTitle = newTitle + correctionCode;
            spectraGraph->SetTitle(newTitle.c_str());
            cout << " Changing Title from " << spectraGraph->GetTitle() << " to " << newTitle << endl;


          }//plusminus loop
        }//detector Index loop

        if(tpc_eff_nom_plus_fit)       delete tpc_eff_nom_plus_fit;
        if(tpc_eff_sys_plus_fit)       delete tpc_eff_sys_plus_fit;
        if(tpc_eff_nom_minus_fit)      delete tpc_eff_nom_minus_fit;
        if(tpc_eff_nom_plus_fit_conf)  delete tpc_eff_nom_plus_fit_conf;
        if(tpc_eff_sys_plus_fit_conf)  delete tpc_eff_sys_plus_fit_conf;
        if(tpc_eff_nom_minus_fit_conf) delete tpc_eff_nom_minus_fit_conf;
        if(tpc_eff_sys_minus_fit_conf) delete tpc_eff_sys_minus_fit_conf;

        if(btof_eff_nom_plus_fit)  delete btof_eff_nom_plus_fit;
        if(btof_eff_sys_plus_fit)  delete btof_eff_sys_plus_fit;
        if(btof_eff_nom_minus_fit) delete btof_eff_nom_minus_fit;

      }//cent loop
      if(energy_loss_nom_plus_fit)       delete energy_loss_nom_plus_fit;
      if(energy_loss_sys_plus_fit)       delete energy_loss_sys_plus_fit;
      if(energy_loss_nom_minus_fit)      delete energy_loss_nom_minus_fit;
      if(energy_loss_sys_minus_fit)      delete energy_loss_sys_minus_fit;
      if(energy_loss_nom_plus_fit_conf)  delete energy_loss_nom_plus_fit_conf;
      if(energy_loss_sys_plus_fit_conf)  delete energy_loss_sys_plus_fit_conf;
      if(energy_loss_nom_minus_fit_conf) delete energy_loss_nom_minus_fit_conf;
      if(energy_loss_sys_minus_fit_conf) delete energy_loss_sys_minus_fit_conf;
    }//yIndex loop
  }//part

  cout << "End Read Efficiency for TPC and Energy Loss" << endl;
  //inFile->Close();
}



//##################################################################################################################################################################################





void RawSpectraModifier::applyTPCEffConf(TGraphErrors* spectraGraph, TGraphAsymmErrors* effGr, TF1* effFunct){

  TFitResultPtr fitResultPtr = effGr->Fit(effFunct,"S","",0.0,2.0);
  //  double points[(int) spectraGraph->GetN()] = spectraGraph->GetX();
  double *points = spectraGraph->GetX();
  double effErrors[(int) spectraGraph->GetN()];
  fitResultPtr->GetConfidenceIntervals((int) spectraGraph->GetN(), 1, 1, points, effErrors, 0.683, false);


  cout << " #########  TITLE CHECK ############" << endl;
  cout << "         ---------- " << spectraGraph->GetName() << endl;
  cout << "         ---------- " << effFunct->GetName() << endl;

  for (int pointIndex=0; pointIndex<spectraGraph->GetN(); pointIndex++){
    //Get the point info
    double mTm0  = spectraGraph->GetX()[pointIndex];
    double yield = spectraGraph->GetY()[pointIndex];
    double xErr  = spectraGraph->GetEX()[pointIndex];
    double yErr  = spectraGraph->GetEY()[pointIndex];
    double eff = effFunct->Eval(mTm0);

    if( eff > 1.0 ) eff = 1.0;
    if( eff <= 0.0 ){ //removing point
      spectraGraph->RemovePoint(pointIndex);
      pointIndex--;
      continue;
    }
    /*if( eff < 0.1 ){ 
      //removing point
      spectraGraph->RemovePoint(pointIndex);
      pointIndex--;
      continue;
    }*/
    //eff = 1e-50;
    double effErr = effErrors[pointIndex];
    double newContent = yield/eff;
    double newContentError = newContent*TMath::Sqrt(
     (yErr*yErr)/(yield*yield) + (effErr*effErr)/(eff*eff));

    cout << "spectraGraph mT-m0: " << mTm0 << " old content: " << yield << " eff: " << eff <<  " new content: " << newContent << endl;

    spectraGraph->SetPoint(pointIndex, mTm0, newContent);
    spectraGraph->SetPointError(pointIndex,xErr,newContentError);

    //    if(detIndex == 0 && pmIndex == 0) cout << "m_spectraPlusDeDx mT-m0: " << m_spectraPlusDeDx[partIndex][centIndex][yIndexModified]->GetX()[pointIndex] << " old content: " << m_spectraPlusDeDx[partIndex][centIndex][yIndexModified]->GetY()[pointIndex] << endl;
    }// pt bin
  }



//##################################################################################################################################################################################



void RawSpectraModifier::loadAndApplyFeedDownCorrections(string a_filename, bool a_useSystematic, int a_backwardsCentNum, int* m_offsets){

    int offsetsHARDCODE[9] = {0,0,0,0,0,0,0,0,0};
    if(m_offsets){
      for(int iii = 0; iii < 9; iii++){
        offsetsHARDCODE[iii] = m_offsets[iii];
      }
    }


    cout << "HARDCODED offsets, CODE MUST BE CHANGED ERROR ERROR" << endl;


    cout << "Starting RawSpectraModifier::loadAndApplyFeedDownCorrections(" << a_filename <<")"<<endl;
    TFile* inFile = new TFile(a_filename.c_str(),"READ");
    if(!inFile){
      cout << "Warning in RawSpectraModifier::loadAndApplyFeedDownCorrections couldn't find file: " << a_filename << endl;
      return;
    }

    //make the TH2D of the efficiencies
    //#####  PARTICLE INDEX LOOP ##############
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){   
      if(!m_particleLoaded[partIndex]) continue;
      cout << "Getting Eff for Part " << partIndex << endl;

      m_hybridFeedDown_rapidity_structure[partIndex] = (TH1D*) inFile->Get(Form("%s/hybridFeedDown_%s_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));

      //###### CENTRALITY INDEX LOOP ####################
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        cout << "  workin on cent " << centIndex << endl;
        #ifdef _RAWSPECMOD_DEBUG_
          cout << " working with pointer : " << m_dEdxSpectraPlus[partIndex][centIndex] << endl;
        #endif
        //###### RAPICITY INDEX LOOP ####################
        for(int yIndex = 0; yIndex < m_dEdxSpectraPlus[partIndex][centIndex]->GetNbinsX(); yIndex ++){
          int rapBin = yIndex + 1;
          int yIndexModified = yIndex - offsetsHARDCODE[partIndex];
          if( yIndexModified < 0 || yIndexModified >= m_dEdxSpectraPlus[partIndex][0]->GetNbinsX()){
            cout << "  WARNING: the yIndexModified = " << yIndexModified << " is not contained in the spectra TH2D " << endl;
            continue;
          }
          double rapidityCenter_dEdx = m_dEdxSpectraPlus[partIndex][0]->GetXaxis()->GetBinCenter(yIndexModified+1);


          int centIndexModified = centIndex;
          if(a_backwardsCentNum != 0) centIndexModified = a_backwardsCentNum - centIndex;

          cout << "Trying " << Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitExp", 
            m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexModified,yIndex) << endl;
          TF1* feed_down_nom_plus_graph = (TF1*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitExp", 
            m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexModified,yIndex)   );
          TF1* feed_down_sys_plus_graph = (TF1*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitPower_Sys", 
            m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexModified,yIndex)  );
          TF1* feed_down_nom_minus_graph = (TF1*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitExp", 
            m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexModified,yIndex) );
          TF1* feed_down_sys_minus_graph = (TF1*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitPower_Sys", 
            m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexModified,yIndex) );

          TGraphErrors* feed_down_nom_plus_graph_conf = (TGraphErrors*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitExp_Conf", 
            m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexModified,yIndex)   );
          TGraphErrors* feed_down_sys_plus_graph_conf = (TGraphErrors*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitPower_Sys_Conf", 
            m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndexModified,yIndex)  );
          TGraphErrors* feed_down_nom_minus_graph_conf = (TGraphErrors*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitExp_Conf", 
            m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexModified,yIndex) );
          TGraphErrors* feed_down_sys_minus_graph_conf = (TGraphErrors*) inFile->Get( Form("%s/FeedDownBackgroundFits/feeddownBackground_%s_cent%02d_yIndex%02d_FitPower_Sys_Conf", 
            m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndexModified,yIndex) );



          TMatrixDSym* feed_down_nom_plus_covMat  = (TMatrixDSym*) inFile->Get(Form("%s/FeedDownBackgroundFits/FeedDownCovMat_%s_Cent%02d_yIndex%d",m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
          TMatrixDSym* feed_down_sys_plus_covMat  = (TMatrixDSym*) inFile->Get(Form("%s/FeedDownBackgroundFits/FeedDownCovMat_%s_Cent%02d_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
          TMatrixDSym* feed_down_nom_minus_covMat = (TMatrixDSym*) inFile->Get(Form("%s/FeedDownBackgroundFits/FeedDownCovMat_%s_Cent%02d_yIndex%d",m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex));
          TMatrixDSym* feed_down_sys_minus_covMat = (TMatrixDSym*) inFile->Get(Form("%s/FeedDownBackgroundFits/FeedDownCovMat_%s_Cent%02d_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex));


          cout << "Ptrs: " << feed_down_nom_plus_graph << " " << feed_down_sys_plus_graph << " " << feed_down_nom_minus_graph << " " << feed_down_sys_minus_graph << endl;
          cout << "Ptrs Conf: " << feed_down_nom_plus_graph_conf << " " << feed_down_sys_plus_graph_conf << " " << feed_down_nom_minus_graph_conf << " " << feed_down_sys_minus_graph_conf << endl;        
          cout << "Ptrs CovMat: " << feed_down_nom_plus_covMat << "  " << feed_down_sys_plus_covMat << "  " << feed_down_nom_minus_covMat << "  " << feed_down_sys_minus_covMat << endl;


          cout << "Applying Feed Down EFF For " << m_partInfo->GetParticleName(partIndex,0) << " Cent: " << centIndex << endl;

          TGraphErrors* spectraGraphs[6] = {NULL,NULL,NULL,NULL,NULL,NULL};
          spectraGraphs[0] = m_spectraPlusDeDx[partIndex][centIndex][yIndexModified];
          spectraGraphs[1] = m_spectraMinusDeDx[partIndex][centIndex][yIndexModified];
          cout << "yIndexModified = " << yIndexModified << endl;
          int yIndexModified_btof = m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetXaxis()->FindBin(rapidityCenter_dEdx) - 1;

          cout << "yIndexModified_btof = " << yIndexModified_btof << endl;

          spectraGraphs[2] = m_spectraPlusInvBeta[partIndex][centIndex][yIndexModified_btof];
          spectraGraphs[3] = m_spectraMinusInvBeta[partIndex][centIndex][yIndexModified_btof];
          int yIndexModified_etof = 0;
          if(m_loadedETOF){
            yIndexModified_etof = m_InvBetaETOFSpectraPlus[partIndex][centIndex]->GetXaxis()->FindBin(rapidityCenter_dEdx) - 1;
            cout << "yIndexModified_etof = " << yIndexModified_etof << endl;
            spectraGraphs[4] = m_spectraPlusInvBetaETOF[partIndex][centIndex][yIndexModified_etof];
            spectraGraphs[5] = m_spectraMinusInvBetaETOF[partIndex][centIndex][yIndexModified_etof];
          }

          // make a copies to make feeddown spectra (spectra * feeddownFraction)
          if(spectraGraphs[0]) m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][0] = (TGraphErrors*) spectraGraphs[0]->Clone();
          if(spectraGraphs[2]) m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][0] = (TGraphErrors*) spectraGraphs[2]->Clone();
          if(spectraGraphs[4]) m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][0] = (TGraphErrors*) spectraGraphs[4]->Clone();
          if(spectraGraphs[1]) m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][1] = (TGraphErrors*) spectraGraphs[1]->Clone();
          if(spectraGraphs[2]) m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][1] = (TGraphErrors*) spectraGraphs[3]->Clone();
          if(spectraGraphs[5]) m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][1] = (TGraphErrors*) spectraGraphs[4]->Clone();

          if(m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][0]){
            m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][0]->SetName(Form("feeddownSpectra_TPC_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
            m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][0]->SetTitle(Form("Feed-Down Spectra TPC"));
            HistogramUtilities::zeroTGraph(m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][0]);
          }
          if(m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][0]){
            m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][0]->SetName(Form("feeddownSpectra_bTOF_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
            m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][0]->SetTitle(Form("Feed-Down Spectra BTOF"));
            HistogramUtilities::zeroTGraph(m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][0]);
          }
          if(m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][0]){
            m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][0]->SetName(Form("feeddownSpectra_eTOF_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
            m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][0]->SetTitle(Form("Feed-Down Spectra TPC"));
            HistogramUtilities::zeroTGraph(m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][0]);
          }
          if(m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][1]){
            m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][1]->SetName(Form("feeddownSpectra_TPC_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
            m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][1]->SetTitle(Form("Feed-Down Spectra TPC"));
            HistogramUtilities::zeroTGraph(m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][1]);
          }
          if(m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][1]){
            m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][1]->SetName(Form("feeddownSpectra_bTOF_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
            m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][1]->SetTitle(Form("Feed-Down Spectra BTOF"));
            HistogramUtilities::zeroTGraph(m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][1]);
          }
          if(m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][1]){
            m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][1]->SetName(Form("feeddownSpectra_eTOF_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex));
            m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][1]->SetTitle(Form("Feed-Down Spectra ETOF"));
            HistogramUtilities::zeroTGraph(m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][1]);
          }

          TGraphErrors* feeddownSpectraGraphs[6];
          feeddownSpectraGraphs[0] = m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][0];
          feeddownSpectraGraphs[1] = m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][0];
          feeddownSpectraGraphs[2] = m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][0];
          feeddownSpectraGraphs[3] = m_feeddownSpectraDeDx[partIndex][centIndex][yIndex][1];
          feeddownSpectraGraphs[4] = m_feeddownSpectraInvBeta[partIndex][centIndex][yIndex][1];
          feeddownSpectraGraphs[5] = m_feeddownSpectraInvBetaETOF[partIndex][centIndex][yIndex][1];

          TF1* feedDownFunct[2] = {NULL,NULL};
          TGraphErrors* feedDownConf[2] = {NULL,NULL};
          TMatrixDSym* feedDownCovMat[2] = {NULL,NULL};
          if(!a_useSystematic){
            feedDownFunct[0] = feed_down_nom_plus_graph;
            feedDownFunct[1] = feed_down_nom_minus_graph;
            feedDownConf[0] = feed_down_nom_plus_graph_conf;
            feedDownConf[1] = feed_down_nom_minus_graph_conf;
            feedDownCovMat[0] = feed_down_nom_plus_covMat;
            feedDownCovMat[1] = feed_down_nom_minus_covMat;
          }else{
            feedDownFunct[0] = feed_down_sys_plus_graph;
            feedDownFunct[1] = feed_down_sys_minus_graph;
            feedDownConf[0] = feed_down_sys_plus_graph_conf;
            feedDownConf[1] = feed_down_sys_minus_graph_conf;
            feedDownCovMat[0] = feed_down_sys_plus_covMat;
            feedDownCovMat[1] = feed_down_sys_minus_covMat;
          }


          for(int detIndex = 0; detIndex < 3; detIndex++){
            if(detIndex == 1 && (yIndexModified_btof < 0 || yIndexModified_btof >= m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ) continue;
            if(detIndex == 2 && m_loadedETOF && (yIndexModified_etof < 0 || yIndexModified_etof >= m_InvBetaETOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ) continue;
            for(int pmIndex = 0; pmIndex<2; pmIndex++){
              TGraphErrors* spectraGraph = spectraGraphs[2*detIndex+pmIndex];
              TGraphErrors* feeddownSpectraGraph = feeddownSpectraGraphs[2*detIndex+pmIndex];
              TF1* funct = feedDownFunct[pmIndex];
              TGraphErrors* confidence = feedDownConf[pmIndex];
              TMatrixDSym* covMat = feedDownCovMat[pmIndex];

              if(!spectraGraph) continue;
              cout << "     Working on det index " << detIndex << "  pmIndex " << pmIndex << " tot " << 2*detIndex+pmIndex << endl;
              cout << "        Ptrs: " << spectraGraph << " " << funct << endl;// << " " << shiftFunctConf << endl;

            if(!funct){
              cout << "No feed down Efficiency Function" << endl;
              continue;
            }

            cout << " #########  TITLE CHECK ############" << endl;
            cout << "         ---------- " << spectraGraph->GetTitle() << endl;
            cout << "         ---------- " << funct->GetName() << endl;

            for (int pointIndex=0; pointIndex<spectraGraph->GetN(); pointIndex++){

              //Get the point info
              double mTm0  = spectraGraph->GetX()[pointIndex];
              double yield = spectraGraph->GetY()[pointIndex];
              double xErr  = spectraGraph->GetEX()[pointIndex];
              double yErr  = spectraGraph->GetEY()[pointIndex];
              

              if(m_chrisCorrectionFileMode){ // for reverse compatability

                double  correctionFactor = 1.0 - funct->Eval(mTm0);
                if(correctionFactor > 1.0) correctionFactor = 1.0;
                if(correctionFactor <= 0.0) correctionFactor = 1E-50;
                int pointNum = HistogramUtilities::getPointClosestToVal(confidence,mTm0);
                double corrError = confidence->GetErrorY(pointNum);
                double newContent = yield*correctionFactor;
                double newContentError = newContent*TMath::Sqrt(
                   (yErr*yErr)/(yield*yield) + (corrError*corrError)/(correctionFactor*correctionFactor));
  
                spectraGraph->SetPoint(pointIndex, mTm0, newContent);
                spectraGraph->SetPointError(pointIndex,xErr,newContentError);

              }else{ // Better new method

                cout << " \n\nERROR: The integral ranges for feed down are not completetly accurate! Note to update code structure\n\n" << endl;
                double lowMtM0 = mTm0 - xErr;
                double highMtM0 = mTm0 + xErr;
                double feeddownFraction    = funct->Integral(lowMtM0,highMtM0)/(highMtM0 - lowMtM0);
                double correctionFactorErr = funct->IntegralError(lowMtM0,highMtM0, funct->GetParameters(), covMat->GetMatrixArray())/(highMtM0 - lowMtM0);
                double correctionFactor = 1.0 - feeddownFraction;
                
                if(correctionFactor > 1.0) correctionFactor = 1.0;
                if(correctionFactor <= 0.0) correctionFactor = 1E-50;
                double newContent = yield*correctionFactor;
                double newContentError = newContent*TMath::Sqrt( (yErr*yErr)/(yield*yield) + (correctionFactorErr*correctionFactorErr)/(correctionFactor*correctionFactor));
                double feeddownContent = yield*feeddownFraction;
                double feeddownContentErr = feeddownContent*TMath::Sqrt( (yErr*yErr)/(yield*yield) + (correctionFactorErr*correctionFactorErr)/(feeddownFraction*feeddownFraction));
                spectraGraph->SetPoint(pointIndex, mTm0, newContent);
                spectraGraph->SetPointError(pointIndex,xErr,newContentError);
                
                feeddownSpectraGraph->SetPoint(pointIndex, mTm0, feeddownContent);
                feeddownSpectraGraph->SetPointError(pointIndex,xErr,feeddownContentErr);
              }// nominal


            }// pt bin 
          }//plusminus
        }//detector Index

        delete feed_down_nom_plus_graph;
        delete feed_down_sys_plus_graph;
        delete feed_down_nom_minus_graph;
        delete feed_down_nom_plus_graph_conf;
        delete feed_down_sys_plus_graph_conf;
        delete feed_down_nom_minus_graph_conf;
        delete feed_down_sys_minus_graph_conf;
      }//yIndex
    }//cent
  }//part

  cout << "End Apply FeedDownCorrections" << endl;
  //inFile->Close();

}



void RawSpectraModifier::loadAndApplyHybridFeedDownCorrections(string a_filename, bool a_useSinglePointInsteadOfIntegral){

  cout << "RawSpectraModifier::loadAndApplyHybridFeedDownCorrections(" << a_filename <<")"<<endl;
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "ERROR: in RawSpectraModifier::loadAndApplyHybridFeedDownCorrections couldn't find file: " << a_filename << endl;
    return;
  }

  string detNames[3] = {"TPC","BTOF","ETOF"};

  //#####  PARTICLE INDEX LOOP ##############
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){   
    if(!m_particleLoaded[partIndex]) continue;


    //if(!m_picoBinner_rapidity_structure[partIndex][0]){
      TH1D* rap_tpc = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_tpc_rap_struct",
               m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
      TH1D* rap_btof = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_btof_rap_struct",
               m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
      TH1D* rap_etof = (TH1D*) inFile->Get(Form("%s/picoBinner_%s_etof_rap_struct",
               m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
      TH1D* rap_embed = (TH1D*) inFile->Get(Form("%s/embedding_%s_rap_struct",
               m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
  
      if(rap_tpc)   m_picoBinner_rapidity_structure[partIndex][0] = (TH1D*) rap_tpc->Clone();
      if(rap_btof)  m_picoBinner_rapidity_structure[partIndex][1] = (TH1D*) rap_btof->Clone();
      if(rap_etof)  m_picoBinner_rapidity_structure[partIndex][2] = (TH1D*) rap_etof->Clone();
    //}

    #ifdef _RAWSPECMOD_DEBUG_
      cout << "  Rapidity Structure Pointers...  TPC: " << m_picoBinner_rapidity_structure[partIndex][0] 
           << "   BTOF: " << m_picoBinner_rapidity_structure[partIndex][1] 
           << "   ETOF: " << m_picoBinner_rapidity_structure[partIndex][2]
           << "    EMBED: " << m_embedding_rapidity_structure[partIndex] << endl;
    #endif

    TH1D* rap_fd = (TH1D*) inFile->Get(Form("%s/hybridFeedDown_%s_rap_struct",
             m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,0).Data()));
    if(rap_fd) m_hybridFeedDown_rapidity_structure[partIndex] = (TH1D*) rap_fd->Clone();

    if(!m_picoBinner_rapidity_structure[partIndex][0] || !m_hybridFeedDown_rapidity_structure[partIndex]){
      cout << " ERROR: No rapidity structures for partIndex: " << partIndex << "  continuing..." << endl;
      continue;
    }

    #ifdef _RAWSPECMOD_DEBUG_
      cout << "  Getting FeedDown Spectra for " << m_partInfo->GetParticleName(partIndex) << endl;
    #endif
    //###### CENTRALITY INDEX LOOP ####################
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      #ifdef _RAWSPECMOD_DEBUG_
        cout << "  workin' on cent " << centIndex << endl;
        cout << "  working with TPC Spectra pointer : " << m_dEdxSpectraPlus[partIndex][centIndex] << endl;
      #endif

      //TH2D* feedDownSpectra[2];
      //feedDownSpectra[0] = (TH2D*) inFile->Get(Form("totalFeedDownSpectra_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex));
      //feedDownSpectra[1] = (TH2D*) inFile->Get(Form("totalFeedDownSpectra_%s_Cent%02d",m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex));


      //#ifdef _RAWSPECMOD_DEBUG_
      //  cout << "FeedDown Spectra Ptrs: " << feedDownSpectra[0] << " " << feedDownSpectra[1] <<  endl;
      //#endif

      //###### RAPICITY INDEX LOOP ####################
      //for(int yIndex = 0; yIndex < m_dEdxSpectraPlus[partIndex][centIndex]->GetNbinsX(); yIndex++){
      for(int y_Counter = 0; y_Counter < 101; y_Counter++){ // goes from -5 to 5
        double rapidity = -5.0 + 0.1*((double) y_Counter);
        cout <<" ptr: " << m_picoBinner_rapidity_structure[partIndex][0] << endl;
        int yBin_tpc  = m_picoBinner_rapidity_structure[partIndex][0]->FindBin(rapidity);
        int yBin_btof = m_picoBinner_rapidity_structure[partIndex][1]->FindBin(rapidity);
        int yBin_etof = 0;
        if(m_loadedETOF) yBin_etof = m_picoBinner_rapidity_structure[partIndex][2]->FindBin(rapidity);
        int yBin_feeddown = m_hybridFeedDown_rapidity_structure[partIndex]->FindBin(rapidity);

        int yIndex_tpc = yBin_tpc - 1;
        int yIndex_btof = yBin_btof - 1;
        int yIndex_etof = yBin_etof - 1;
        int yIndex_feeddown = yBin_feeddown - 1;



        //int rapBin = yIndex + 1;
        //int yIndexModified = yIndex;
        //if( yIndexModified < 0 || yIndexModified >= m_dEdxSpectraPlus[partIndex][0]->GetNbinsX()){
        //  cout << "  WARNING: the yIndexModified = " << yIndexModified << " is not contained in the spectra TH2D " << endl;
        //  continue;
        //}
        //double rapidityCenter_dEdx = m_dEdxSpectraPlus[partIndex][0]->GetXaxis()->GetBinCenter(yIndexModified+1);

        #ifdef _RAWSPECMOD_DEBUG_
          cout << "Applying Feed Down Correction For " << m_partInfo->GetParticleName(partIndex,0) << " Cent: " << centIndex << " rapidity: " << rapidity << endl;
        #endif

        TGraphErrors* feedDownSpectra[2];
        string positiveFDName = Form("%s/FeedDownBackgroundGraphs/hybridFeedDown_%s_Cent%02d_yIndex%02d_graph",m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_feeddown);
        string negativeFDName = Form("%s/FeedDownBackgroundGraphs/hybridFeedDown_%s_Cent%02d_yIndex%02d_graph",m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_feeddown);
        feedDownSpectra[0] = (TGraphErrors*) inFile->Get(positiveFDName.c_str());
        feedDownSpectra[1] = (TGraphErrors*) inFile->Get(negativeFDName.c_str());

        TF1* feedDownSpectraFit[2];
        string positiveFDFitName = Form("%s/FeedDownBackgroundFits/FeedDownFit_Hybrid_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_feeddown);
        string negativeFDFitName = Form("%s/FeedDownBackgroundFits/FeedDownFit_Hybrid_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_feeddown);
        feedDownSpectraFit[0] = (TF1*) inFile->Get(positiveFDFitName.c_str());
        feedDownSpectraFit[1] = (TF1*) inFile->Get(negativeFDFitName.c_str());

        string positiveFDCovMatName = Form("%s/FeedDownBackgroundFits/FeedDownCovMat_Hybrid_%s_Cent%02d_yIndex%d",m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_feeddown);
        string negativeFDCovMatName = Form("%s/FeedDownBackgroundFits/FeedDownCovMat_Hybrid_%s_Cent%02d_yIndex%d",m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_feeddown);
        TMatrixDSym* feedDownSpectraFitCovMat[2];
        feedDownSpectraFitCovMat[0] = (TMatrixDSym*) inFile->Get(positiveFDCovMatName.c_str());
        feedDownSpectraFitCovMat[1] = (TMatrixDSym*) inFile->Get(negativeFDCovMatName.c_str());

        TF1* feedDownSpectraFit_lowPt[2];
        string positiveFDFitName_lowPt = Form("%s/FeedDownBackgroundFits/FeedDownFit_Hybrid_%s_Cent%d_yIndex%d_lowPt",m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_feeddown);
        string negativeFDFitName_lowPt = Form("%s/FeedDownBackgroundFits/FeedDownFit_Hybrid_%s_Cent%d_yIndex%d_lowPt",m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_feeddown);
        feedDownSpectraFit_lowPt[0] = (TF1*) inFile->Get(positiveFDFitName_lowPt.c_str());
        feedDownSpectraFit_lowPt[1] = (TF1*) inFile->Get(negativeFDFitName_lowPt.c_str());

        string positiveFDCovMatName_lowPt = Form("%s/FeedDownBackgroundFits/FeedDownCovMat_Hybrid_%s_Cent%02d_yIndex%d_lowPt",m_partInfo->GetParticleName(partIndex,1).Data(),m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,yIndex_feeddown);
        string negativeFDCovMatName_lowPt = Form("%s/FeedDownBackgroundFits/FeedDownCovMat_Hybrid_%s_Cent%02d_yIndex%d_lowPt",m_partInfo->GetParticleName(partIndex,-1).Data(),m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex,yIndex_feeddown);
        TMatrixDSym* feedDownSpectraFitCovMat_lowPt[2];
        feedDownSpectraFitCovMat_lowPt[0] = (TMatrixDSym*) inFile->Get(positiveFDCovMatName_lowPt.c_str());
        feedDownSpectraFitCovMat_lowPt[1] = (TMatrixDSym*) inFile->Get(negativeFDCovMatName_lowPt.c_str());


        TGraphErrors* spectraGraphs[6] = {NULL,NULL,NULL,NULL,NULL,NULL};
        spectraGraphs[0] = m_spectraPlusDeDx[partIndex][centIndex][yIndex_tpc];
        spectraGraphs[1] = m_spectraMinusDeDx[partIndex][centIndex][yIndex_tpc];
        spectraGraphs[2] = m_spectraPlusInvBeta[partIndex][centIndex][yIndex_btof];
        spectraGraphs[3] = m_spectraMinusInvBeta[partIndex][centIndex][yIndex_btof];
        if(m_loadedETOF){
          spectraGraphs[4] = m_spectraPlusInvBetaETOF[partIndex][centIndex][yIndex_etof];
          spectraGraphs[5] = m_spectraMinusInvBetaETOF[partIndex][centIndex][yIndex_etof];
        }

        if(   HistogramUtilities::outsideHistogramRange(m_picoBinner_rapidity_structure[partIndex][0],rapidity) 
           && HistogramUtilities::outsideHistogramRange(m_picoBinner_rapidity_structure[partIndex][1],rapidity)
           && HistogramUtilities::outsideHistogramRange(m_hybridFeedDown_rapidity_structure[partIndex],rapidity)){
          #ifdef _RAWSPECMOD_DEBUG_
            cout << " Not Computing rapidity = " << rapidity << " becasue it is not in bounds of any histogram..." << endl;
            continue;
          #endif
        }
        cout << " Rapidity Index TPC: " << yIndex_tpc << "   BTOF: " << yIndex_btof << "   ETOF: " << yIndex_etof << "   Feed-Down: " << yIndex_feeddown << endl;

        #ifdef _RAWSPECMOD_DEBUG_
          cout << "  Pointer Check: " << endl;
          cout << "    F-D FitName: " << positiveFDFitName << endl;
          cout << "    F-D FitName: " << negativeFDFitName << endl;
          cout << "    F-D Fit: " << feedDownSpectraFit[0] << " " << feedDownSpectraFit[1] << endl;
          cout << "    F-D CovMatName: " << positiveFDCovMatName << endl;
          cout << "    F-D CovMatName: " << negativeFDCovMatName << endl;
          cout << "    F-D CovMat: " << feedDownSpectraFitCovMat[0] << " " << feedDownSpectraFitCovMat[1] << endl;          
        #endif

        /*int rapIndex_feeddown = -1;
        if(feedDownSpectra[0]){
          rapIndex_feeddown = feedDownSpectra[0]->GetXaxis()->FindBin(rapidityCenter_dEdx) - 1;
        }else if(feedDownSpectra[1]){
          rapIndex_feeddown = feedDownSpectra[1]->GetXaxis()->FindBin(rapidityCenter_dEdx) - 1;
        }else{
          cout << "\n\n WARNING: There are no feed down spectra to subtract! Continuing to next rapidity! \n\n" << endl;
          continue;
        }

        TGraphErrors* feedDownGraph[2] = {NULL,NULL};
        if(feedDownSpectra[0]) feedDownGraph[0] = HistogramUtilities::makeGraph(feedDownSpectra[0],rapIndex_feeddown+1);
        if(feedDownSpectra[1]) feedDownGraph[1] = HistogramUtilities::makeGraph(feedDownSpectra[1],rapIndex_feeddown+1);
        */

        for(int detIndex = 0; detIndex < 3; detIndex++){
          if(!m_loadedETOF && detIndex == 2) continue;
          /*if(detIndex == 1 && (yIndexModified_btof < 0 || yIndexModified_btof >= m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ){
            cout << "  WARNING: yIndexModified_btof " << yIndexModified_btof << " is beyond the scope of the spectra histogram" << endl;
            continue;
          }
          if(detIndex == 2 && m_loadedETOF && (yIndexModified_etof < 0 || yIndexModified_etof >= m_InvBetaETOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ){
            cout << "  WARNING: yIndexModified_etof " << yIndexModified_etof << " is beyond the scope of the spectra histogram" << endl;
            continue;
          } */

          for(int pmIndex = 0; pmIndex<2; pmIndex++){
            cout << "\n\n###################################################################################################################" << endl;
            cout << "  Now Feed-Down Modifying " << m_partInfo->GetParticleName(partIndex,pmIndex == 0?1:-1) << "  Cent " << centIndex << " Detector: " << detNames[detIndex] << " y = " << rapidity << endl;
            cout << "###################################################################################################################" << endl;

            TGraphErrors* spectraGraph     = spectraGraphs[2*detIndex+pmIndex];
            TGraphErrors* feeddownSpectra  = feedDownSpectra[pmIndex];
            TF1* feedDownFit               = feedDownSpectraFit[pmIndex];
            TMatrixDSym* feedDownFitCovMat = feedDownSpectraFitCovMat[pmIndex];
            TF1* feedDownFit_lowPt               = feedDownSpectraFit_lowPt[pmIndex];
            TMatrixDSym* feedDownFitCovMat_lowPt = feedDownSpectraFitCovMat_lowPt[pmIndex];
            if(!spectraGraph || !feedDownFit || !feedDownFitCovMat){
              cout << " WARNING: No feeddown spectra graph for this rapidity (" << rapidity << ") - spectraGraph:" << spectraGraph << "feedDownFit:" << feedDownFit << "feedDownFitCovMat :" << feedDownFitCovMat << endl;
              continue;
            }
            #ifdef _RAWSPECMOD_DEBUG_
              cout << "     Working on det index " << detIndex << "  pmIndex " << pmIndex << " tot " << 2*detIndex+pmIndex;
              cout << "        Ptrs: " << spectraGraph  << endl;// << " " << shiftFunctConf << endl;
              cout << " #########  TITLE CHECK ############" << endl;
              cout << "         ---------- " << spectraGraph->GetTitle() << endl;
              cout << "         ---------- " << feedDownFit->GetTitle() << endl;
              cout << "         ---------- " << feeddownSpectra->GetTitle() << endl;
            #endif

            double lowFitRange = 0.0;
            double highFitRange = 20.0;
            feedDownFit->GetRange(lowFitRange,highFitRange);


            int numPointsActuallyApplied_FeedDown = 0;
            for(int pointIndex=0; pointIndex < spectraGraph->GetN(); pointIndex++){
              //Get the point info
              double mTm0  = spectraGraph->GetX()[pointIndex];
              double yield = spectraGraph->GetY()[pointIndex];
              double xErr  = spectraGraph->GetEX()[pointIndex];
              double yErr  = spectraGraph->GetEY()[pointIndex];
              
              double lowMtM0 = mTm0 - xErr;
              double highMtM0 = mTm0 + xErr;
              if(xErr == 0.0){
                cout << "ERROR: The xErr of the point of the spectra is zero. Feed Down Subtraction wont work..." << endl;
                if(pointIndex < spectraGraph->GetN()-1) xErr = 0.5*(spectraGraph->GetX()[pointIndex+1] - spectraGraph->GetX()[pointIndex]);
                else xErr = 0.5*(spectraGraph->GetX()[pointIndex] - spectraGraph->GetX()[pointIndex-1]);
              }


              double feedDown_SpectraVal    = 0.0;
              double feedDown_SpectraValErr = 0.0;
              #ifdef _RAWSPECMOD_DEBUG_
                string integralMode = "High pT Integral";
              #endif
              if(mTm0 > lowFitRange){ // lower than fit start (high statistics region)
                feedDown_SpectraVal    = feedDownFit->Integral(lowMtM0,highMtM0)/(highMtM0 - lowMtM0);
                feedDown_SpectraValErr = feedDownFit->IntegralError(lowMtM0,highMtM0, feedDownFit->GetParameters(), feedDownFitCovMat->GetMatrixArray())/(highMtM0 - lowMtM0);
              }else{  // in the fit region
                if(feedDownFit_lowPt){
                  feedDown_SpectraVal    = feedDownFit->Integral(lowMtM0,highMtM0)/(highMtM0 - lowMtM0);
                  feedDown_SpectraValErr = feedDownFit->IntegralError(lowMtM0,highMtM0, feedDownFit->GetParameters(), feedDownFitCovMat->GetMatrixArray())/(highMtM0 - lowMtM0);
                  #ifdef _RAWSPECMOD_DEBUG_
                    integralMode = "Low pT Integral";
                  #endif
                }else{
                  HistogramUtilities::evalWithError(feeddownSpectra,mTm0,feedDown_SpectraVal,feedDown_SpectraValErr);
                  //HistogramUtilities::integrateGraph(feeddownSpectraGraph, lowMtM0,highMtM0, feedDown_SpectraVal, feedDown_SpectraValErr, true);
                  #ifdef _RAWSPECMOD_DEBUG_
                    integralMode = "Low pT Point Eval";
                  #endif
                }

              }
              double newContent = yield - feedDown_SpectraVal;
              double newContentError = sqrt( (yErr*yErr) + (feedDown_SpectraValErr*feedDown_SpectraValErr));
              //if(a_useSinglePointInsteadOfIntegral){
              //  newContent = yield - eval;
              //  newContentError = sqrt( (yErr*yErr) + (evalErr*evalErr));
              //}

              #ifdef _RAWSPECMOD_DEBUG_
                cout << "   mTm0: " << mTm0 << " Frac FD: " << feedDown_SpectraVal/yield << " N_inc = " << yield << " +- " << yErr 
                           << "   N_fd = " <<  feedDown_SpectraVal << " +- " << feedDown_SpectraValErr 
                           << " ->  N_pri = " << newContent << " +/- " << newContentError << "  by " << integralMode << endl;
              #endif

              if(newContent < 0){
                //newContent = 0;
                //newContentError = 0; keep error for fitting
                cout << "ERROR: Hybrid FeedDown Subtraction yielded a negative spectra point for mTm0=" << mTm0 << " with N = " << newContent << "  so removing the point" << endl;
                spectraGraph->RemovePoint(pointIndex);
                pointIndex--;
                continue;
              }

              spectraGraph->SetPoint(pointIndex, mTm0, newContent);
              spectraGraph->SetPointError(pointIndex,xErr,newContentError);
              numPointsActuallyApplied_FeedDown++;
            }// mTm0 point loop 


            double fracApplied_FeedDown = ((double) numPointsActuallyApplied_FeedDown)/((double)spectraGraph->GetN());
            string correctionCode = " Corr:";
            if(fracApplied_FeedDown > 0.5)   correctionCode = correctionCode + "F";
            else{
              cout << "ERROR: Feed-Down correction was not applied to more than half of the points..." << endl;
            }
            string newTitle = spectraGraph->GetTitle();
            newTitle = newTitle + correctionCode;
            spectraGraph->SetTitle(newTitle.c_str());

          }//plusminus loop
        }//detector Index loop
      }//yIndex
    }//cent
  }//part

  cout << "End Apply FeedDownCorrections" << endl;
  //inFile->Close();


}








void RawSpectraModifier::loadAndApplyKnockoutCorrections(string a_filename, int a_backwardsCentNum, int m_offset){
  int partIndex = 2;
  //cout << "HARDCODED offsets, CODE MUST BE CHANGED ERROR ERROR" << endl;
  #ifdef _RAWSPECMOD_DEBUG_
    cout << "Starting RawSpectraModifier::loadAndApplyKnockoutCorrections(" << a_filename <<")"<<endl;
  #endif
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "Warning in RawSpectraModifier::loadAndApplyKnockoutCorrections couldn't find file: " << a_filename << endl;
    return;
  }

  //m_knockout_mTm0_structure = (TH1D*) inFile->Get(Form("%s/knockout_mTm0_struct",m_partInfo->GetParticleName(partIndex,1).Data()));

  //make the TH2D of the efficiencies
  
  if(!m_particleLoaded[partIndex]) return;
  cout << "Getting Eff for Part " << partIndex << endl;
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
    cout << "  workin on cent " << centIndex << endl;
    for(int yIndex = 0; yIndex < m_dEdxSpectraPlus[partIndex][centIndex]->GetNbinsX(); yIndex ++){
      int rapBin = yIndex + 1;
      int yIndexModified = yIndex - m_offset;
      if( yIndexModified < 0 || yIndexModified >= m_dEdxSpectraPlus[partIndex][0]->GetNbinsX()) continue;
      double rapidityCenter_dEdx = m_dEdxSpectraPlus[partIndex][0]->GetXaxis()->GetBinCenter(yIndexModified+1);
      
      int centIndexModified = centIndex;
      if(a_backwardsCentNum != 0) centIndexModified = a_backwardsCentNum - centIndex;
      cout << "Trying " << Form("%s/KnockoutBackgroundFits/knockoutProton_%d_%d_Fit", 
        m_partInfo->GetParticleName(partIndex,1).Data(),centIndexModified,yIndex) << endl;
      TF1* knockout_nom_plus_graph = (TF1*) inFile->Get( Form("%s/KnockoutBackgroundFits/knockoutProton_%d_%d_Fit", 
        m_partInfo->GetParticleName(partIndex,1).Data(),centIndexModified,yIndex));
      cout << "Ptrs: " << knockout_nom_plus_graph  << endl;
      cout << "Applying Feed Down EFF For " << m_partInfo->GetParticleName(partIndex,0) << " Cent: " << centIndex << endl;

      TGraphErrors* spectraGraphs[3] = {NULL,NULL,NULL};
      spectraGraphs[0] = m_spectraPlusDeDx[partIndex][centIndex][yIndexModified];
      cout << "yIndexModified = " << yIndexModified << endl;
      int yIndexModified_btof = m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetXaxis()->FindBin(rapidityCenter_dEdx) - 1;
      cout << "yIndexModified_btof = " << yIndexModified_btof << endl;
      spectraGraphs[1] = m_spectraPlusInvBeta[partIndex][centIndex][yIndexModified_btof];
      int yIndexModified_etof = 0;
      if(m_loadedETOF){
        yIndexModified_etof = m_InvBetaETOFSpectraPlus[partIndex][centIndex]->GetXaxis()->FindBin(rapidityCenter_dEdx) - 1;
        cout << "yIndexModified_etof = " << yIndexModified_etof << endl;
        spectraGraphs[2] = m_spectraPlusInvBetaETOF[partIndex][centIndex][yIndexModified_etof];
      }


      for(int detIndex = 0; detIndex < 3; detIndex++){
        if(detIndex == 1 && (yIndexModified_btof < 0 || yIndexModified_btof >= m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ) continue;
        if(detIndex == 2 && m_loadedETOF && (yIndexModified_etof < 0 || yIndexModified_etof >= m_InvBetaETOFSpectraPlus[partIndex][centIndex]->GetNbinsX() ) ) continue;
        
        TGraphErrors* spectraGraph = spectraGraphs[detIndex];
        TF1* funct = knockout_nom_plus_graph;
        if(!spectraGraph) continue;
        cout << "     Working on det index " << detIndex   << endl;
          cout << "        Ptrs: " << spectraGraph << " " << funct << endl;// << " " << shiftFunctConf << endl;
          if(!funct){
            cout << "No knockout Efficiency Function" << endl;
            continue;
          }
          cout << " #########  TITLE CHECK ############" << endl;
          cout << "         ---------- " << spectraGraph->GetTitle() << endl;
          cout << "         ---------- " << funct->GetName() << endl;

          for (int pointIndex=0; pointIndex<spectraGraph->GetN(); pointIndex++){

            //Get the point info
            double mTm0  = spectraGraph->GetX()[pointIndex];
            double yield = spectraGraph->GetY()[pointIndex];
            double xErr  = spectraGraph->GetEX()[pointIndex];
            double yErr  = spectraGraph->GetEY()[pointIndex];
            double  correctionFactor = 1.0 - funct->Eval(mTm0);

            if(correctionFactor > 1.0) correctionFactor = 1.0;
            if(correctionFactor <= 0.0) correctionFactor = 1E-50;

            double newContent = yield*correctionFactor;
            double newContentError = yErr*correctionFactor;
            spectraGraph->SetPoint(pointIndex, mTm0, newContent);
            spectraGraph->SetPointError(pointIndex,xErr,newContentError);
            
          }// pt bin 

      }//detector Index
      delete knockout_nom_plus_graph;

    }//yIndex
  }//cent


  cout << "End Apply FeedDownCorrections" << endl;
  //inFile->Close();

}




//Different detectors have different rapidity binning
//The bins are centered at the same values, but some have more or fewer bins, so the indicies are not the same
//Convert rapidity index from input detector to output detector
int RawSpectraModifier::convertRapidityIndex(int a_yIndexSource, int a_sourceDetectorIndex, int a_targetDetectorIndex, int a_partIndex){
  //Detector indices:
  //TPC - 0
  //TOF - 1
  //eToF - 2

  if( !m_particleLoaded[a_partIndex] ) return -999;

  if(a_sourceDetectorIndex == a_targetDetectorIndex) return a_yIndexSource;
  TH2D *sourceDetector = NULL;
  TH2D *targetDetector = NULL;

  switch (a_sourceDetectorIndex){
    case 0:
    sourceDetector = m_dEdxSpectraPlus[a_partIndex][0];
    break;
    case 1:
    sourceDetector = m_InvBetaBTOFSpectraPlus[a_partIndex][0];
    break;
    case 2:
    if(!m_InvBetaETOFSpectraPlus[a_partIndex][0]) return -999;
    sourceDetector = m_InvBetaETOFSpectraPlus[a_partIndex][0];
    break;
    default:
    break;
  }

  switch (a_targetDetectorIndex){
    case 0:
    targetDetector = m_dEdxSpectraPlus[a_partIndex][0];
    break;
    case 1:
    targetDetector = m_InvBetaBTOFSpectraPlus[a_partIndex][0];
    break;
    case 2:
    if(!m_InvBetaETOFSpectraPlus[a_partIndex][0]) return -999;
    targetDetector = m_InvBetaETOFSpectraPlus[a_partIndex][0];
    break;
    default:
    break;
  }

  int rapBin = a_yIndexSource +1;
  double rapidityCenter_source = sourceDetector->GetXaxis()->GetBinCenter(rapBin);
  int rapBin_target = targetDetector->GetXaxis()->FindBin(rapidityCenter_source);
  if(rapBin_target < 1 || rapBin_target > targetDetector->GetNbinsX()) return -999;
  int yIndex_target = rapBin_target - 1;

  return yIndex_target;
}



void RawSpectraModifier::convertSpectraToInvariant(  double a_efficiency, double a_effErr){
  #ifdef _RAWSPECMOD_DEBUG_
    cout << "Call to RawSpectraModifier::convertSpectraToInvariant()" << endl;
  #endif
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(!m_particleLoaded[partIndex]) continue;
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      TH2D* specPM[6] = {NULL,NULL,NULL,NULL,NULL,NULL};
      specPM[0] = m_dEdxSpectraPlus[partIndex][centIndex];
      specPM[1] = m_dEdxSpectraMinus[partIndex][centIndex];
      specPM[2] = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
      specPM[3] = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];
      if(m_loadedETOF){
        specPM[4] = m_InvBetaETOFSpectraPlus[partIndex][centIndex];
        specPM[5] = m_InvBetaETOFSpectraMinus[partIndex][centIndex];
      }
      double numberOfCentEvents = m_centEvents[partIndex]->GetBinContent(centIndex+1);
      #ifdef _RAWSPECMOD_DEBUG_
        cout << "#######################################################################################################" << endl;
        cout << "Converting to Invariant: " << m_partInfo->GetParticleName(partIndex,0) << "   Centrality Index: " << centIndex << endl;
      #endif
      for(int histoIndex = 0; histoIndex < 6; histoIndex++){
        if(!specPM[histoIndex]) continue;
        //############# CONVERT FROM NUMBER OF COUNTS IN EACH BIN TO INVARIANT YIELDS ###############
        for(int rapBin = 1; rapBin <= specPM[histoIndex]->GetNbinsX(); rapBin++){
          double rapidity       = specPM[histoIndex]->GetXaxis()->GetBinCenter(rapBin);
          double rapidity_width = specPM[histoIndex]->GetXaxis()->GetBinWidth(rapBin);
          for(int mTm0Bin = 1; mTm0Bin <= specPM[histoIndex]->GetNbinsY(); mTm0Bin++){
            double content      = specPM[histoIndex]->GetBinContent(rapBin, mTm0Bin);
            if(content == 0.0) continue;
            double mTm0      = specPM[histoIndex]->GetYaxis()->GetBinCenter(mTm0Bin);
            double mTm0Width = specPM[histoIndex]->GetYaxis()->GetBinWidth(mTm0Bin);
            double mT        = mTm0 + m_partInfo->GetParticleMass(partIndex);
            double scale = 1.0 /(rapidity_width * mTm0Width * 2.0 * TMath::Pi() *
                                 mT * numberOfCentEvents * a_efficiency); 
            double contentError = specPM[histoIndex]->GetBinError(rapBin, mTm0Bin);           
            double percentError = TMath::Sqrt(
                1.0/numberOfCentEvents
              +(contentError*contentError/(content*content))
              +(a_effErr*a_effErr/(a_efficiency*a_efficiency))
             );
            double newContent = content * scale;
            double newContentError = percentError * newContent;
            specPM[histoIndex]->SetBinContent(rapBin,mTm0Bin, newContent);
            specPM[histoIndex]->SetBinError(rapBin,mTm0Bin, newContentError);
            #ifdef _RAWSPECMOD_DEBUG_
              cout << "   y: " << rapidity << " mT-m0: " << mTm0 << " mT: " << mT << " dy: " << rapidity_width << " dmT: " << mTm0Width << " N_evt: " << numberOfCentEvents 
                   << "    " << content << " -> " << newContent << endl;
            #endif
          } //pt Bin
        }// eta Bin
      }//histo loop
    }//cent
  }//part
  #ifdef _RAWSPECMOD_DEBUG_
    cout << " Finished Converting Spectra to Invariant" << endl;
  #endif
}





void RawSpectraModifier::applyBTOFEff(bool a_useFit){
  cout << "\n\n\n\n\n\n\n\n############  APPLYING BTOF EFFICIENCY ###########\n\n\n\n\n" << endl;
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(!m_particleLoaded[partIndex]) continue;
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      cout << "##################################" << endl;
      cout << "Applying BTOF EFF For " << m_partInfo->GetParticleName(partIndex,0) << " Cent: " << centIndex << endl;
      TH2D* bTofSpecPM[2] = {NULL,NULL};
      bTofSpecPM[0] = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
      bTofSpecPM[1] = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];
      TH2D* effHistoPM[2] = {NULL,NULL};
      //if statement here
      if(a_useFit){
        effHistoPM[0] = m_BTOFEffPlus[partIndex][centIndex];
        effHistoPM[1] = m_BTOFEffMinus[partIndex][centIndex];
      }
      else{
        effHistoPM[0] = m_BTOFEffPlusRatio[partIndex][centIndex];
        effHistoPM[1] = m_BTOFEffMinusRatio[partIndex][centIndex];
      }
      for(int pmIndex = 0; pmIndex<2; pmIndex++){
        if(!bTofSpecPM[pmIndex]) continue;
        if(!effHistoPM[pmIndex]){
          cout << "Error: The BTOF Efficiency doesn't exist for this particle. Continuing..." << endl;
          continue;
        }
        //############# BTOF Eff Applicaton ###############
        for(int etaBin = 1; etaBin <= bTofSpecPM[pmIndex]->GetNbinsX(); etaBin++){
          double eta = bTofSpecPM[pmIndex]->GetXaxis()->GetBinCenter(etaBin);
          for(int pTBin = 1; pTBin <= bTofSpecPM[pmIndex]->GetNbinsY(); pTBin++){
            double amp_ofPoint = bTofSpecPM[pmIndex]->GetBinContent(etaBin,pTBin);
            double ampErr_ofPoint = bTofSpecPM[pmIndex]->GetBinError(etaBin,pTBin);
            double pT_ofPoint = bTofSpecPM[pmIndex]->GetYaxis()->GetBinCenter(pTBin);
            double tofEff = effHistoPM[pmIndex]->GetBinContent(effHistoPM[pmIndex]->GetXaxis()->FindBin(eta),effHistoPM[pmIndex]->GetYaxis()->FindBin(pT_ofPoint));
            double tofEffErr = effHistoPM[pmIndex]->GetBinError(effHistoPM[pmIndex]->GetXaxis()->FindBin(eta),effHistoPM[pmIndex]->GetYaxis()->FindBin(pT_ofPoint));
            double newContent = amp_ofPoint/tofEff;
            double newContentError = newContent*TMath::Sqrt(
             (ampErr_ofPoint*ampErr_ofPoint)/(amp_ofPoint*amp_ofPoint)
             +(tofEffErr*tofEffErr)/(tofEff*tofEff));
            if(tofEff > 1.0 || tofEff <= 0.0 || std::isnan(newContentError) || std::isnan(newContent)){
              if(!(amp_ofPoint == 0.0)){
                #ifdef _RAWSPECMOD_DEBUG_
                  cout << "Nan Content or tofEff>1" << endl;
                  cout << "Part: " << partIndex << " cent: " << centIndex << "  RapBin: " << etaBin << "mTm0Bin: " << pTBin << " tofEff: " << tofEff
                   << "     content: " << amp_ofPoint << "  err: " << ampErr_ofPoint
                   <<  "    newCont: " << newContent << "  newErr: " << newContentError;
                  cout << "    BTOF NOT APPLIED!!" << endl;
                #endif
              }
              bTofSpecPM[pmIndex]->SetBinContent(etaBin,pTBin, 0.0);
              bTofSpecPM[pmIndex]->SetBinError(etaBin,pTBin, 0);
            }else{
              #ifdef _RAWSPECMOD_DEBUG_
              cout << "TOF APP: eta: " << eta << " pT: " << pT_ofPoint <<  " tofEff: " << tofEff << " err: " << tofEffErr
               << "\tcontent:\t" << amp_ofPoint << "  +/- " << ampErr_ofPoint
               << "\tnewContent:\t" << newContent << " +/- " << newContentError;
              cout << "    BTOF applied Okay" << endl;
              #endif
              bTofSpecPM[pmIndex]->SetBinContent(etaBin,pTBin, newContent);
              bTofSpecPM[pmIndex]->SetBinError(etaBin,pTBin, newContentError);
            }            

          } 
        }
        
      }//plusminus
    }//cent
  }//part
  cout << "done applying BTOF Efficiency correction" << endl;
};




void RawSpectraModifier::matchBTOFtoTPC(double* a_momMatch, bool a_matchBinByBin, double* a_lowRapToFit, double* a_highRapToFit){


  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){   
    

    if(!m_particleLoaded[partIndex]) continue;

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      cout << "Matching BTOF For " << m_partInfo->GetParticleName(partIndex,0) << " Cent: " << centIndex << endl;

      TH2D* specPM[4] = {NULL,NULL,NULL,NULL};
      specPM[0] = m_dEdxSpectraPlus[partIndex][centIndex];
      specPM[1] = m_dEdxSpectraMinus[partIndex][centIndex];
      specPM[2] = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
      specPM[3] = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];

      bool killPlus = false;
      bool killMinus = false;


      TH1D* matchHistoPlus = specPM[2]->ProjectionX(Form("btofMatchHisto_%s_cent%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex));
      matchHistoPlus->SetTitle(Form("BTOF Scaling Factor to dE/dx for %s with Centrality Index %02d",m_partInfo->GetParticleSymbol(partIndex,1).Data(),centIndex));
      HistogramUtilities::zeroBins(matchHistoPlus);

      TH1D* matchHistoMinus = specPM[2]->ProjectionX(Form("btofMatchHisto_%s_cent%02d",m_partInfo->GetParticleName(partIndex,-1).Data(),centIndex));
      matchHistoMinus->SetTitle(Form("BTOF Scaling Factor to dE/dx for %s with Centrality Index %02d",m_partInfo->GetParticleSymbol(partIndex,-1).Data(),centIndex));
      HistogramUtilities::zeroBins(matchHistoMinus);      

      for(int yIndex = 0; yIndex < m_dEdxSpectraPlus[partIndex][centIndex]->GetNbinsX(); yIndex ++){
        int rapBin = yIndex + 1;
        double rapidityCenter_dEdx = m_dEdxSpectraPlus[partIndex][0]->GetXaxis()->GetBinCenter(yIndex+1);
        int rapBin_bTOF = m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetXaxis()->FindBin(rapidityCenter_dEdx);
        if(rapBin_bTOF < 1 || rapBin_bTOF > m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetNbinsX()) continue;
        int yIndex_btof = rapBin_bTOF - 1;

        double mTm0_fromMom = PhysMath::mTm0_from_y_p(rapidityCenter_dEdx, a_momMatch[partIndex], m_partInfo->GetParticleMass(partIndex));

        //if(rapidityCenter_dEdx < a_lowRapToFit[partIndex] || rapidityCenter_dEdx > a_highRapToFit[partIndex]) continue;

        TGraphErrors* spectraGraphs[4] = {NULL,NULL,NULL,NULL};
        spectraGraphs[0] = m_spectraPlusDeDx[partIndex][centIndex][yIndex];
        spectraGraphs[1] = m_spectraMinusDeDx[partIndex][centIndex][yIndex];
        spectraGraphs[2] = m_spectraPlusInvBeta[partIndex][centIndex][yIndex_btof];
        spectraGraphs[3] = m_spectraMinusInvBeta[partIndex][centIndex][yIndex_btof];

        
        for(int pmIndex = 0; pmIndex<2; pmIndex++){
          TGraphErrors* dEdxGraph = spectraGraphs[pmIndex];
          TGraphErrors* bTOFGraph = spectraGraphs[pmIndex+2];
          if(!dEdxGraph || !bTOFGraph) continue;
          if(bTOFGraph->GetN() <= 0 || dEdxGraph->GetN() <= 0) continue;

          int pointDeDx =  HistogramUtilities::getPointClosestToVal(dEdxGraph,mTm0_fromMom);
          int pointBTOF =  HistogramUtilities::getPointClosestToVal(bTOFGraph,mTm0_fromMom); 
          
          double mTm0_dEdx  = dEdxGraph->GetX()[pointDeDx];
          double yield_dEdx = dEdxGraph->GetY()[pointDeDx];
          double xErr_dEdx  = dEdxGraph->GetEX()[pointDeDx];
          double yErr_dEdx  = dEdxGraph->GetEY()[pointDeDx];

          double mTm0_bTOF  = bTOFGraph->GetX()[pointBTOF];
          double yield_bTOF = bTOFGraph->GetY()[pointBTOF];
          double xErr_bTOF  = bTOFGraph->GetEX()[pointBTOF];
          double yErr_bTOF  = bTOFGraph->GetEY()[pointBTOF];


          double matchingVal = yield_dEdx / yield_bTOF;
          double matchingErr = matchingVal * sqrt(yErr_bTOF*yErr_bTOF/(yield_bTOF*yield_bTOF) 
            + yErr_dEdx*yErr_dEdx/(yield_dEdx*yield_dEdx));
          if(pmIndex == 0){
            matchHistoPlus->SetBinContent(rapBin_bTOF,matchingVal);
            matchHistoPlus->SetBinError(rapBin_bTOF,matchingErr);
          }else{
            matchHistoMinus->SetBinContent(rapBin_bTOF,matchingVal);
            matchHistoMinus->SetBinError(rapBin_bTOF,matchingErr);
          }

          if(a_matchBinByBin){
            for (int pointIndex=0; pointIndex<bTOFGraph->GetN(); pointIndex++){

              //Get the point info
              double mTm0  = bTOFGraph->GetX()[pointIndex];
              double yield = bTOFGraph->GetY()[pointIndex];
              double xErr  = bTOFGraph->GetEX()[pointIndex];
              double yErr  = bTOFGraph->GetEY()[pointIndex];

              double newContent = yield*matchingVal;
              double newContentError = newContent*sqrt(
               (yErr*yErr)/(yield*yield) + (matchingErr*matchingErr)/(matchingVal*matchingVal));
              bTOFGraph->SetPoint(pointIndex, mTm0, newContent);
              bTOFGraph->SetPointError(pointIndex,xErr,newContentError);

            }// pt bin 
          }



        }//pmindex
      }//yIndex




      if(a_matchBinByBin){
        m_btofScaleHistoPlus[partIndex][centIndex] = matchHistoPlus;
        m_btofScaleHistoMinus[partIndex][centIndex] = matchHistoMinus;
        continue;
      }




      if(matchHistoPlus->GetMaximum() == 0.0) killPlus = true;
      if(matchHistoMinus->GetMaximum() == 0.0) killMinus = true;

      double fitScalePlus = -999;
      double fitScalePlusErr = -999;

      if(!killPlus){
        //Now we have the match histos, fit them with a constant
        TF1* constFunctPlus = new TF1("constFunctPlus", "[0]",-2.5,2.5);
        int statusPlus = matchHistoPlus->Fit("constFunctPlus","EQ","",a_lowRapToFit[partIndex], a_highRapToFit[partIndex]);
        TF1* fitResultPlus = matchHistoPlus->GetFunction("constFunctPlus");
        fitScalePlus = fitResultPlus->GetParameter(0);
        fitScalePlusErr = fitResultPlus->GetParError(0);
      }

      double fitScaleMinus = -999;
      double fitScaleMinusErr = -999;

      if(!killMinus){
        //Same for Minus
        TF1* constFunctMinus = new TF1("constFunctMinus", "[0]",-2.5,2.5);
        int statusMinus = matchHistoMinus->Fit("constFunctMinus","EQ","",a_lowRapToFit[partIndex], a_highRapToFit[partIndex]);
        TF1* fitResultMinus = matchHistoMinus->GetFunction("constFunctMinus");
        fitScaleMinus = fitResultMinus->GetParameter(0);
        fitScaleMinusErr = fitResultMinus->GetParError(0);
      }

      for(int yIndex = 0; yIndex < m_dEdxSpectraPlus[partIndex][centIndex]->GetNbinsX(); yIndex ++){
        int rapBin = yIndex + 1;
        double rapidityCenter_dEdx = m_dEdxSpectraPlus[partIndex][0]->GetXaxis()->GetBinCenter(yIndex+1);
        int rapBin_bTOF = m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetXaxis()->FindBin(rapidityCenter_dEdx);
        if(rapBin_bTOF < 1 || rapBin_bTOF > m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetNbinsX()) continue;
        int yIndex_btof = rapBin_bTOF - 1;

        TGraphErrors* spectraGraphs[2] = {NULL,NULL};
        spectraGraphs[0] = m_spectraPlusInvBeta[partIndex][centIndex][yIndex_btof];
        spectraGraphs[1] = m_spectraMinusInvBeta[partIndex][centIndex][yIndex_btof];
        
        for(int pmIndex = 0; pmIndex<2; pmIndex++){
          if(pmIndex == 0 && killPlus) continue;
          else if(pmIndex == 1 && killMinus) continue;
          TGraphErrors* bTOFGraph = spectraGraphs[pmIndex];
          if(!bTOFGraph) continue;
          if(bTOFGraph->GetN() <= 0) continue;

          double matchVal = fitScalePlus;
          double matchErr = fitScalePlusErr;
          if(pmIndex == 1){
            matchVal = fitScaleMinus;
            matchErr = fitScaleMinusErr;
          }

          for (int pointIndex=0; pointIndex<bTOFGraph->GetN(); pointIndex++){

            //Get the point info
            double mTm0  = bTOFGraph->GetX()[pointIndex];
            double yield = bTOFGraph->GetY()[pointIndex];
            double xErr  = bTOFGraph->GetEX()[pointIndex];
            double yErr  = bTOFGraph->GetEY()[pointIndex];

            double newContent = yield*matchVal;
            double newContentError = newContent*sqrt(
             (yErr*yErr)/(yield*yield) + (matchErr*matchErr)/(matchVal*matchVal));
            bTOFGraph->SetPoint(pointIndex, mTm0, newContent);
            bTOFGraph->SetPointError(pointIndex,xErr,newContentError);
            
          }// pt bin 

        }//pmindex
      }//yIndex

      m_btofScaleHistoPlus[partIndex][centIndex] = matchHistoPlus;
      m_btofScaleHistoMinus[partIndex][centIndex] = matchHistoMinus;


    }//cent
  }//part

  cout << "Finished TOF MATCHING Function" << endl;
}


void RawSpectraModifier::cleanSpectraFXT(){
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

      TH2D* specPM[2] = {NULL,NULL};
      specPM[0] = m_dEdxSpectraPlus[partIndex][centIndex];
      specPM[1] = m_dEdxSpectraMinus[partIndex][centIndex];

      for(int pmIndex = 0; pmIndex< 2; pmIndex++){
        //#### CLEAN TPC SPECTRA #######
        for(int rapBin = 1; rapBin <= specPM[pmIndex]->GetNbinsX(); rapBin++){
          TGraphErrors* spectraGraphs[2] = {m_spectraPlusDeDx[partIndex][centIndex][rapBin -1],m_spectraMinusDeDx[partIndex][centIndex][rapBin -1]};
          double rapCenter = specPM[pmIndex]->GetXaxis()->GetBinCenter(rapBin);
          bool point0 = true;
          for (int pointIndex=0; pointIndex<spectraGraphs[pmIndex]->GetN(); pointIndex++){
            double mTm0  = spectraGraphs[pmIndex]->GetX()[pointIndex];
            if(partIndex == 0 && pmIndex == 0){
              if(rapCenter >= -0.35 && mTm0 > 1.5){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }
              if(rapCenter == -1.65 && (mTm0 > 0.35 || mTm0 < 0.2)){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }
              if(rapCenter == -1.55 && mTm0 > 0.4){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }
              if(rapCenter == -1.45 && mTm0 > 0.95){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }              
              if(rapCenter == -0.55 && mTm0 > 1.6){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }   
              if(rapCenter == -0.15 && mTm0 > 1.2){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }   
              if(rapCenter == -0.05 && mTm0 > 1.0){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }   
            }//if pip
            if(partIndex == 0 && pmIndex == 1){
              if(rapCenter <= -1.75 && mTm0 < 0.4){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }
              if(rapCenter <= -1.55 && rapCenter >= -1.65 && point0){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
                point0 = false;
              }
              if(rapCenter >= -0.45 && mTm0 > 1.15){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }
              if(rapCenter >= -0.15 && mTm0 > 1.1){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }
            }

            point0 = false;
          }//end point loop
        }//End rap     
      }//End PM
    }//End Cent
  }//end partIndex
}

void RawSpectraModifier::cleanSpectra(double* a_momStart_tpc,
  double* a_momStart_btof, 
  double* a_momStart_etof,
  double* a_momStop_tpc,
  double* a_momStop_btof, 
  double* a_momStop_etof,
  double* a_lowEtaCutTPC,
  double* a_highEtaCutTPC,
  double* a_lowEtaCutBTOF,
  double* a_highEtaCutBTOF,
  double* a_minMtM0,
  double* a_maxMtM0,
  double a_maxErrorFraction,
  double* a_minRapidity,
  double* a_maxRapidity,
  bool a_removeByPreviousRatio,
  vector<vector<double>> a_btof_mtm0_start_matrix){

  #ifdef _RAWSPECMOD_DEBUG_
  cout << "Starting RawSpectraModifier::cleanSpectra(" << a_momStart_tpc << ", " << a_momStart_btof << ", " << 
    a_momStart_etof << ", " << a_momStop_tpc << ", " << a_momStop_btof << ", " <<  a_momStop_etof << ", " << a_lowEtaCutTPC << ", " 
    << a_highEtaCutTPC << ", " << a_lowEtaCutBTOF << ", " << a_highEtaCutBTOF << ", " << a_minMtM0 << ", " <<
    a_maxMtM0 << ", " << a_maxErrorFraction << ", " << a_minRapidity << ", " << a_maxRapidity << ", " << a_removeByPreviousRatio << " )" << endl;
  #endif

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    m_momStart_tpc[partIndex]  = a_momStart_tpc[partIndex];
    m_momStart_btof[partIndex] = a_momStart_btof[partIndex];
    m_momStart_etof[partIndex] = a_momStart_etof[partIndex];
    m_momStop_tpc[partIndex]   = a_momStop_tpc[partIndex];
    m_momStop_btof[partIndex]  = a_momStop_btof[partIndex];
    m_momStop_etof[partIndex]  = a_momStop_etof[partIndex];

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      TH2D* specPM[6] = {NULL,NULL,NULL,NULL,NULL,NULL};
      specPM[0] = m_dEdxSpectraPlus[partIndex][centIndex];
      specPM[1] = m_dEdxSpectraMinus[partIndex][centIndex];
      specPM[2] = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
      specPM[3] = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];
      if(m_loadedETOF){
        specPM[4] = m_InvBetaETOFSpectraPlus[partIndex][centIndex];
        specPM[5] = m_InvBetaETOFSpectraMinus[partIndex][centIndex];
      }
      for(int pmIndex = 0; pmIndex< 2; pmIndex++){
        //#### CLEAN TPC SPECTRA #######
        cout << "Cleaning TPC Spectra" << endl;
        for(int rapBin = 1; specPM[pmIndex] && rapBin <= specPM[pmIndex]->GetNbinsX(); rapBin++){
          TGraphErrors* spectraGraphs[2] = {m_spectraPlusDeDx[partIndex][centIndex][rapBin -1],m_spectraMinusDeDx[partIndex][centIndex][rapBin -1]};
          double rapCenter = specPM[pmIndex]->GetXaxis()->GetBinCenter(rapBin);

          cout << "\n\n###################################################################################################################" << endl;
          cout << "  Now Cleaning " << m_partInfo->GetParticleName(partIndex,pmIndex == 0?1:-1) << "  Cent " << centIndex << " Detector: " << "TPC" << " y = " << rapCenter << endl;
          cout << "###################################################################################################################" << endl;

          //If the energy loss correction pushed the mTm0 points high enough that they are out of order, delete everything before the inversion
          if(spectraGraphs[pmIndex]->GetN() > 1){
            int point_of_start_inversion = -1;
            for(int pointIndex = spectraGraphs[pmIndex]->GetN() - 2; pointIndex >= 0; pointIndex--){
              if(spectraGraphs[pmIndex]->GetX()[pointIndex] > spectraGraphs[pmIndex]->GetX()[pointIndex+1]){
                point_of_start_inversion = pointIndex;
                break;
              }
            }        
            if(point_of_start_inversion >= 0){
              for(int pointIndex = point_of_start_inversion; pointIndex >= 0; pointIndex--){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                #ifdef _RAWSPECMOD_DEBUG_ 
                 cout << "point " << pointIndex << " removed since energy loss correction over corrected this lowPt point (mTm0 not in ascending order)" << endl;
                #endif
              }                
            } 
          }

          double stop_tpc_mtm0 = PhysMath::mTm0(m_partInfo->GetParticleMass(partIndex),a_momStop_tpc[partIndex]);
          if(a_btof_mtm0_start_matrix.size() > 0){
            stop_tpc_mtm0 = 5.0;
            for(unsigned int iii = 0; iii < a_btof_mtm0_start_matrix.size(); iii++){
              if(a_btof_mtm0_start_matrix[iii][0] == partIndex && fabs(a_btof_mtm0_start_matrix[iii][1] - fabs(rapCenter)) < 0.005){
                stop_tpc_mtm0 = a_btof_mtm0_start_matrix[iii][2];
                break;
              }
            }
          }

          for (int pointIndex=0; pointIndex < spectraGraphs[pmIndex]->GetN(); pointIndex++){
            //Get the point info
            double mTm0  = spectraGraphs[pmIndex]->GetX()[pointIndex];
            double yield = spectraGraphs[pmIndex]->GetY()[pointIndex];
            double yieldPrev = 0.1;
            if(pointIndex != 0) yieldPrev = spectraGraphs[pmIndex]->GetY()[pointIndex-1];
            double xErr  = spectraGraphs[pmIndex]->GetEX()[pointIndex];
            double yErr  = spectraGraphs[pmIndex]->GetEY()[pointIndex];
            double momentum = PhysMath::pTot_from_y_mTm0(rapCenter, mTm0, m_partInfo->GetParticleMass(partIndex));   
            double eta      = PhysMath::eta_from_y_mTm0(rapCenter, mTm0, m_partInfo->GetParticleMass(partIndex));
            #ifdef _RAWSPECMOD_DEBUG_
              cout << "############################" << endl;
              cout << "spectra name: " << spectraGraphs[pmIndex]->GetName() << " mT-m0: " << mTm0 << " yield: " << yield << endl;
              cout << "rapidity: " << rapCenter  << " mTm0: " << mTm0 <<  "eta: " << eta << " momentum: " << momentum   << endl;
              cout << "yield: " << yield  << "    yieldPrev: " << yieldPrev  << "    yield/yieldPrev: " << yield/yieldPrev << endl;
            #endif
           /* if(momentum > a_momStop_tpc[partIndex] || momentum < a_momStart_tpc[partIndex]
             || std::isnan(mTm0) || std::isnan(yield) || std::isnan(xErr) || std::isnan(yErr)
             || eta < a_lowEtaCutTPC[partIndex] || eta > a_highEtaCutTPC[partIndex]
             || a_minMtM0[partIndex] > mTm0 || a_maxErrorFraction < yErr/yield || (pointIndex != 0 && (1.0*yield/yieldPrev) > 2)){*/

            //if(momentum > a_momStop_tpc[partIndex]){
            if(mTm0 > stop_tpc_mtm0){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since momentum was above stop" << endl;
             #endif
           }else if(momentum < a_momStart_tpc[partIndex]){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since mometum was below start" << endl;
             #endif
           }else if(std::isnan(mTm0)){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed mtm0 NAN" << endl;
             #endif
           }else if(std::isnan(yield)){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since yield NAN" << endl;
             #endif
           }else if(std::isnan(xErr)){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since mtm0 err NAN" << endl;
             #endif
           }else if(std::isnan(yErr)){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since yeild err NAN" << endl;
             #endif
           }else if(eta < a_lowEtaCutTPC[partIndex]){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since eta below start" << endl;
             #endif
           }else if(eta > a_highEtaCutTPC[partIndex]){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since eta above cut" << endl;
             #endif
           }else if(a_minMtM0[partIndex] > mTm0){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since mtm0 below min" << endl;
             #endif
           }else if(a_maxMtM0[partIndex] < mTm0){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since mtm0 above max" << endl;
             #endif
           }else if(a_minRapidity[partIndex] > rapCenter){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since rap below min" << endl;
             #endif
           }else if(a_maxRapidity[partIndex] < rapCenter){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since rap above max" << endl;
             #endif
           }else if(a_maxErrorFraction < yErr/yield){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed since error/yield too high" << endl;
             #endif
           }else if(a_removeByPreviousRatio && (pointIndex != 0 && (1.0*yield/yieldPrev) > 2)){
             spectraGraphs[pmIndex]->RemovePoint(pointIndex);
             pointIndex--;
             #ifdef _RAWSPECMOD_DEBUG_ 
               cout << "point removed yield/yieldPrev > 2" << endl;
             #endif
           }
          }// pt bin 
        }// eta Bin
        cout << "Cleaning TPC Spectra Successful" << endl;

        //#### CLEAN BTOF SPECTRA #######
        cout << "Cleaning BTOF Spectra" << endl;
        for(int rapBin = 1; specPM[pmIndex+2] && rapBin <= specPM[pmIndex+2]->GetNbinsX(); rapBin++){
          TGraphErrors* spectraGraphs[2] = {m_spectraPlusInvBeta[partIndex][centIndex][rapBin-1],
          m_spectraMinusInvBeta[partIndex][centIndex][rapBin-1]};
          double rapCenter = specPM[pmIndex+2]->GetXaxis()->GetBinCenter(rapBin);

          cout << "\n\n###################################################################################################################" << endl;
          cout << "  Now Cleaning " << m_partInfo->GetParticleName(partIndex,pmIndex == 0?1:-1) << "  Cent " << centIndex << " Detector: " << "BTOF" << " y = " << rapCenter << endl;
          cout << "###################################################################################################################" << endl;


          //If the energy loss correction pushed the mTm0 points high enough that they are out of order, delete everything before the inversion
          if(spectraGraphs[pmIndex]->GetN() > 1){
            int point_of_start_inversion = -1;
            for(int pointIndex = spectraGraphs[pmIndex]->GetN() - 2; pointIndex >= 0; pointIndex--){
              if(spectraGraphs[pmIndex]->GetX()[pointIndex] > spectraGraphs[pmIndex]->GetX()[pointIndex+1]){
                point_of_start_inversion = pointIndex;
                break;
              }
            }        
            if(point_of_start_inversion >= 0){
              for(int pointIndex = point_of_start_inversion; pointIndex >= 0; pointIndex--){
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                #ifdef _RAWSPECMOD_DEBUG_ 
                 cout << "point " << pointIndex << " removed since energy loss correction over corrected this lowPt point (mTm0 not in ascending order)" << endl;
                #endif
              }                
            } 
          }


          double start_btof_mtm0 = PhysMath::mTm0(m_partInfo->GetParticleMass(partIndex),a_momStart_btof[partIndex]);
          if(a_btof_mtm0_start_matrix.size() > 0){
            start_btof_mtm0 = 5.0;
            for(unsigned int iii = 0; iii < a_btof_mtm0_start_matrix.size(); iii++){
              if(a_btof_mtm0_start_matrix[iii][0] == partIndex && fabs(a_btof_mtm0_start_matrix[iii][1] - fabs(rapCenter)) < 0.005){
                start_btof_mtm0 = a_btof_mtm0_start_matrix[iii][2];
                break;
              }
            }
          }

          
          for (int pointIndex=0; pointIndex<spectraGraphs[pmIndex]->GetN(); pointIndex++){
          //Get the point info
            double mTm0  = spectraGraphs[pmIndex]->GetX()[pointIndex];
            double yield = spectraGraphs[pmIndex]->GetY()[pointIndex];
            double yieldPrev = 0.1;
            if(pointIndex != 0) yieldPrev = spectraGraphs[pmIndex]->GetY()[pointIndex-1];
            double xErr  = spectraGraphs[pmIndex]->GetEX()[pointIndex];
            double yErr  = spectraGraphs[pmIndex]->GetEY()[pointIndex];
            double momentum = PhysMath::pTot_from_y_mTm0(rapCenter, mTm0, m_partInfo->GetParticleMass(partIndex));   
            double eta      = PhysMath::eta_from_y_mTm0(rapCenter, mTm0, m_partInfo->GetParticleMass(partIndex));

            if(momentum > a_momStop_btof[partIndex] || mTm0 < start_btof_mtm0
             || std::isnan(mTm0) || std::isnan(yield) || std::isnan(xErr) || std::isnan(yErr)
             || eta < a_lowEtaCutBTOF[partIndex] || eta > a_highEtaCutBTOF[partIndex]
             || a_minMtM0[partIndex] > mTm0 || a_maxErrorFraction < yErr/yield || (a_removeByPreviousRatio && (pointIndex != 0 && (1.0*yield/yieldPrev) > 2) )
             || a_maxMtM0[partIndex] < mTm0 
             || a_minRapidity[partIndex] > rapCenter ||a_maxRapidity[partIndex] < rapCenter) {
              #ifdef _RAWSPECMOD_DEBUG_
                cout << "yield: " << yield <<  "    yieldPrev: " << yieldPrev <<  "    yield/yieldPrev: " << yield/yieldPrev << "\t has been Removed" << endl;
              #endif
                spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
            }
          }// pt bin 
        }// eta Bin
        cout << "Cleaning BTOF Spectra Successful" << endl;


        //#### CLEAN ETOF SPECTRA #######
        if(m_loadedETOF){
          for(int rapBin = 1; specPM[pmIndex+4] && rapBin <= specPM[pmIndex+4]->GetNbinsX(); rapBin++){
            TGraphErrors* spectraGraphs[2] = {m_spectraPlusInvBetaETOF[partIndex][centIndex][rapBin-1],
            m_spectraMinusInvBetaETOF[partIndex][centIndex][rapBin-1]};
            double rapCenter = specPM[pmIndex+4]->GetXaxis()->GetBinCenter(rapBin);

            cout << "\n\n###################################################################################################################" << endl;
            cout << "  Now Cleaning " << m_partInfo->GetParticleName(partIndex,pmIndex == 0?1:-1) << "  Cent " << centIndex << " Detector: " << "ETOF" << " y = " << rapCenter << endl;
            cout << "###################################################################################################################" << endl;
            //If the energy loss correction pushed the mTm0 points high enough that they are out of order, delete everything before the inversion
            if(spectraGraphs[pmIndex]->GetN() > 1){
              int point_of_start_inversion = -1;
              for(int pointIndex = spectraGraphs[pmIndex]->GetN() - 2; pointIndex >= 0; pointIndex--){
                if(spectraGraphs[pmIndex]->GetX()[pointIndex] > spectraGraphs[pmIndex]->GetX()[pointIndex+1]){
                  point_of_start_inversion = pointIndex;
                  break;
                }
              }        
              if(point_of_start_inversion >= 0){
                for(int pointIndex = point_of_start_inversion; pointIndex >= 0; pointIndex--){
                  spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                  #ifdef _RAWSPECMOD_DEBUG_ 
                   cout << "point " << pointIndex << " removed since energy loss correction over corrected this lowPt point (mTm0 not in ascending order)" << endl;
                  #endif
                }                
              } 
            }


            
            for (int pointIndex=0; pointIndex<spectraGraphs[pmIndex]->GetN(); pointIndex++){
            //Get the point info
              double mTm0  = spectraGraphs[pmIndex]->GetX()[pointIndex];
              double yield = spectraGraphs[pmIndex]->GetY()[pointIndex];
              double xErr  = spectraGraphs[pmIndex]->GetEX()[pointIndex];
              double yErr  = spectraGraphs[pmIndex]->GetEY()[pointIndex];
              double momentum = PhysMath::pTot_from_y_mTm0(rapCenter, mTm0, m_partInfo->GetParticleMass(partIndex));   
              if(momentum > a_momStop_etof[partIndex] || momentum < a_momStart_etof[partIndex]
               || std::isnan(mTm0) || std::isnan(yield) || std::isnan(xErr) || std::isnan(yErr)
               || a_minMtM0[partIndex] > mTm0 || a_maxErrorFraction < yErr/yield || a_maxMtM0[partIndex] < mTm0 
               || a_minRapidity[partIndex] > rapCenter ||a_maxRapidity[partIndex] < rapCenter){
                 spectraGraphs[pmIndex]->RemovePoint(pointIndex);
                 pointIndex--;
              }
            }// pt bin 
          }// eta Bin
        }// if etof


      }//plus minus loop

    }//cent loop
  }//part loop

  cout << " Finished Cleaning Spectra" << endl;
}// end of clean spectra


void RawSpectraModifier::trimSpectra(vector< vector <double> > a_trimData, vector< vector<double> > a_pointByPointData, double a_mTm0Tolerence){
  //rap spec to modify loop
  #ifdef _RAWSPECMOD_DEBUG_
          cout << "Trim Data Size: " << a_trimData.size() << endl;
  #endif
  for(int iii = 0; iii < (int) a_trimData.size(); iii++){
    if(a_trimData[iii].size() != 6){
      cout << "ERROR: Pt by Pt at index " << iii << " doesn't have 5 parameters..." << endl;
      continue;
    }
    //### interpret trim data ###
    int partIndex = (int) a_trimData[iii][0];
    int centIndex = (int) a_trimData[iii][1];
    double charge = a_trimData[iii][2];
    double rapVal = a_trimData[iii][3];
    double minMtM0 =   a_trimData[iii][4];
    double maxMtM0 =   a_trimData[iii][5];
    cout << "here" << endl;
    cout << "ptr: " << m_picoBinner_rapidity_structure[partIndex][0] << endl;
    if(!m_picoBinner_rapidity_structure[partIndex][0]){
      cout << "ERROR: You have specified a particle that doen't have a rapidity structure.... Part: " << partIndex << endl;
      continue;
    }
    int rapIndexTPC  = m_picoBinner_rapidity_structure[partIndex][0]->FindBin(rapVal) - 1;
    int rapIndexBTOF = m_picoBinner_rapidity_structure[partIndex][1]->FindBin(rapVal) - 1;
    int rapIndexETOF = 0;
    if(m_loadedETOF) rapIndexETOF = m_picoBinner_rapidity_structure[partIndex][2]->FindBin(rapVal) - 1;
    cout << "here again" << endl;
    #ifdef _RAWSPECMOD_DEBUG_
      cout << "\n\n###################################################################################################################" << endl;
      cout << "  Now Trimming " << m_partInfo->GetParticleName(partIndex,charge) << "  Cent " << centIndex  << " y = " << rapVal << " mTm0 = [" << minMtM0 << "," << maxMtM0 << "]" << endl;
      cout << "###################################################################################################################" << endl;
      cout << "            tpcIndex: " << rapIndexTPC << "  btofIndex: " << rapIndexBTOF << "  etofIndex: " << rapIndexETOF << endl;
    #endif
    //### Get Spectra In Question ###
    TGraphErrors* tpcSpectrum = NULL;
    TGraphErrors* btofSpectrum = NULL;
    TGraphErrors* etofSpectrum = NULL;
    if(charge < 0.0){
      if(rapIndexTPC >= 0) tpcSpectrum  = m_spectraMinusDeDx[partIndex][centIndex][rapIndexTPC];
      if(rapIndexBTOF >= 0) btofSpectrum = m_spectraMinusInvBeta[partIndex][centIndex][rapIndexBTOF];
      if(m_loadedETOF && rapIndexETOF) etofSpectrum = m_spectraMinusInvBetaETOF[partIndex][centIndex][rapIndexETOF];
    }else{
      if(rapIndexTPC >= 0) tpcSpectrum  = m_spectraPlusDeDx[partIndex][centIndex][rapIndexTPC];
      if(rapIndexBTOF >= 0) btofSpectrum = m_spectraPlusInvBeta[partIndex][centIndex][rapIndexBTOF];
      if(m_loadedETOF && rapIndexETOF) etofSpectrum = m_spectraPlusInvBetaETOF[partIndex][centIndex][rapIndexETOF];
    }

    //### Trim TPC Points ###
    for (int pointIndex = 0; tpcSpectrum && pointIndex < tpcSpectrum->GetN(); pointIndex++){
      double mTm0  = tpcSpectrum->GetX()[pointIndex];
      if(mTm0 < minMtM0 || mTm0 > maxMtM0){
        tpcSpectrum->RemovePoint(pointIndex);
        #ifdef _RAWSPECMOD_DEBUG_
          cout << "  Removed TPC: " << m_partInfo->GetParticleName(partIndex,charge).Data() << " Cent: " << centIndex << " Rap: " << rapVal << "  " << mTm0 << "!=[" << minMtM0 << "," << maxMtM0 << "]" << endl;
        #endif
        pointIndex--;
      }
    }

    //### Trim BTOF Points ###
    for (int pointIndex =0; btofSpectrum && pointIndex < btofSpectrum->GetN(); pointIndex++){
      double mTm0  = btofSpectrum->GetX()[pointIndex];
      if(mTm0 < minMtM0 || mTm0 > maxMtM0){
        btofSpectrum->RemovePoint(pointIndex);
        #ifdef _RAWSPECMOD_DEBUG_
          cout << "  Removed BTOF: " << m_partInfo->GetParticleName(partIndex,charge).Data() << " Cent: " << centIndex << " Rap: " << rapVal << "  " << mTm0 << "!=[" << minMtM0 << "," << maxMtM0 << "]" << endl;
        #endif
        pointIndex--;
      }
    }

    //### Trim ETOF Points ###
    for (int pointIndex =0; m_loadedETOF && etofSpectrum && pointIndex < etofSpectrum->GetN(); pointIndex++){
      double mTm0  = etofSpectrum->GetX()[pointIndex];
      if(mTm0 < minMtM0 || mTm0 > maxMtM0){
        etofSpectrum->RemovePoint(pointIndex);
        #ifdef _RAWSPECMOD_DEBUG_
          cout << "  Removed ETOF: " << m_partInfo->GetParticleName(partIndex,charge).Data() << " Cent: " << centIndex << " Rap: " << rapVal << "  " << mTm0 << "!=[" << minMtM0 << "," << maxMtM0 << "]" << endl;
        #endif
        pointIndex--;
      }
    }
  }// end trim loop


  //##############   specific point deletion loop ##############################
  #ifdef _RAWSPECMOD_DEBUG_
    cout << "Point By Point Data Size: " << a_pointByPointData.size() << endl;
  #endif
  for(int iii = 0; iii < (int) a_pointByPointData.size(); iii++){
    if(a_pointByPointData[iii].size() != 5){
      cout << "ERROR: Pt by Pt at index " << iii << " doesn't have 5 parameters..." << endl;
      continue;
    }
    //### interpret trim data ###
    int partIndex = (int) a_pointByPointData[iii][0];
    int centIndex = (int) a_pointByPointData[iii][1];
    double charge = a_pointByPointData[iii][2];
    double rapVal = a_pointByPointData[iii][3];
    double mTm0_to_delete =   a_pointByPointData[iii][4];

    #ifdef _RAWSPECMOD_DEBUG_
      cout << "\n\n###################################################################################################################" << endl;
      cout << "  Now Spot Cleaning " << m_partInfo->GetParticleName(partIndex,charge) << "  Cent " << centIndex  << " y = " << rapVal << " mTm0 to delete = " << mTm0_to_delete << endl;
      cout << "###################################################################################################################" << endl;
    #endif

    cout << "ptr: " << m_picoBinner_rapidity_structure[partIndex][0] << endl;
    if(!m_picoBinner_rapidity_structure[partIndex][0]){
      cout << "ERROR: You have specified a particle that doen't have a rapidity structure.... Part: " << partIndex << endl;
      continue;
    }

    int rapIndexTPC  = m_picoBinner_rapidity_structure[partIndex][0]->FindBin(rapVal) - 1;
    int rapIndexBTOF = m_picoBinner_rapidity_structure[partIndex][1]->FindBin(rapVal) - 1;
    int rapIndexETOF = 0;
    if(m_loadedETOF) rapIndexETOF= m_picoBinner_rapidity_structure[partIndex][2]->FindBin(rapVal) - 1;
    
    #ifdef _RAWSPECMOD_DEBUG_
      cout << "            tpcIndex: " << rapIndexTPC << "  btofIndex: " << rapIndexBTOF << "  etofIndex: " << rapIndexETOF << endl;
    #endif
    //### Get Spectra In Question ###
    TGraphErrors* tpcSpectrum = NULL;
    TGraphErrors* btofSpectrum = NULL;
    TGraphErrors* etofSpectrum = NULL;
    if(charge < 0.0){
      if(rapIndexTPC >= 0) tpcSpectrum  = m_spectraMinusDeDx[partIndex][centIndex][rapIndexTPC];
      if(rapIndexBTOF >= 0) btofSpectrum = m_spectraMinusInvBeta[partIndex][centIndex][rapIndexBTOF];
      if(m_loadedETOF && rapIndexETOF) etofSpectrum = m_spectraMinusInvBetaETOF[partIndex][centIndex][rapIndexETOF];
    }else{
      if(rapIndexTPC >= 0) tpcSpectrum  = m_spectraPlusDeDx[partIndex][centIndex][rapIndexTPC];
      if(rapIndexBTOF >= 0) btofSpectrum = m_spectraPlusInvBeta[partIndex][centIndex][rapIndexBTOF];
      if(m_loadedETOF && rapIndexETOF) etofSpectrum = m_spectraPlusInvBetaETOF[partIndex][centIndex][rapIndexETOF];
    }

    //### Trim TPC Points ###
    for (int pointIndex =0;tpcSpectrum && pointIndex < tpcSpectrum->GetN(); pointIndex++){
      double mTm0  = tpcSpectrum->GetX()[pointIndex];
      if( fabs(mTm0 - mTm0_to_delete) < a_mTm0Tolerence ){
        tpcSpectrum->RemovePoint(pointIndex);
        #ifdef _RAWSPECMOD_DEBUG_
          cout << "  Removed TPC: " << m_partInfo->GetParticleName(partIndex,charge).Data() << " Cent: " << centIndex << " Rap: " << rapVal << "  mTm0:" << mTm0 << "  in tolerence (" << a_mTm0Tolerence << ") of " << mTm0_to_delete << endl;
        #endif
        pointIndex--;
      }
    }

    //### Trim BTOF Points ###
    for (int pointIndex =0; btofSpectrum && pointIndex < btofSpectrum->GetN(); pointIndex++){
      double mTm0  = btofSpectrum->GetX()[pointIndex];
      if( fabs(mTm0 - mTm0_to_delete) < a_mTm0Tolerence ){
        btofSpectrum->RemovePoint(pointIndex);
        #ifdef _RAWSPECMOD_DEBUG_
          cout << "  Removed BTOF: " << m_partInfo->GetParticleName(partIndex,charge).Data() << " Cent: " << centIndex << " Rap: " << rapVal << "  mTm0:" << mTm0 << "  in tolerence (" << a_mTm0Tolerence << ") of " << mTm0_to_delete << endl;
        #endif
        pointIndex--;
      }
    }

    //### Trim ETOF Points ###
    for (int pointIndex =0; m_loadedETOF && etofSpectrum && pointIndex < etofSpectrum->GetN(); pointIndex++){
      double mTm0  = etofSpectrum->GetX()[pointIndex];
      if( fabs(mTm0 - mTm0_to_delete) < a_mTm0Tolerence ){
        etofSpectrum->RemovePoint(pointIndex);
        #ifdef _RAWSPECMOD_DEBUG_
          cout << "  Removed ETOF: " << m_partInfo->GetParticleName(partIndex,charge).Data() << " Cent: " << centIndex << " Rap: " << rapVal << "  mTm0:" << mTm0 << "  in tolerence (" << a_mTm0Tolerence << ") of " << mTm0_to_delete << endl;
        #endif
        pointIndex--;
      }
    }
  }// end point by point deletion loop
}





void RawSpectraModifier::writeMemberData(string a_filename, string a_option){
  #ifdef _RAWSPECMOD_DEBUG_
    cout << "Call to RawSpectraModifier::writeMemberData(string " << a_filename << ", string " << a_option << ")" << endl;
  #endif
  TFile* outFile = new TFile(a_filename.c_str(),a_option.c_str());
  if(!outFile->IsOpen()){
    cout << "ERROR: COULD NOT OPEN OUTPUT FILE! " << a_filename << endl;
  }
  #ifdef _RAWSPECMOD_DEBUG_
    else{
      cout << a_filename.c_str() << " successfully created and opened" << endl;
    }
  #endif
  outFile->cd();

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    if(m_centEvents[partIndex]) m_centEvents[partIndex]->Write();
    HistogramUtilities::ConditionalWrite(m_picoBinner_rapidity_structure[partIndex][0]);
    HistogramUtilities::ConditionalWrite(m_picoBinner_rapidity_structure[partIndex][1]);
    HistogramUtilities::ConditionalWrite(m_picoBinner_rapidity_structure[partIndex][2]);
    HistogramUtilities::ConditionalWrite(m_embedding_rapidity_structure[partIndex]);
    HistogramUtilities::ConditionalWrite(m_hybridFeedDown_rapidity_structure[partIndex]);

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      HistogramUtilities::ConditionalWrite(m_dEdxSpectraPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_dEdxSpectraMinus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_InvBetaBTOFSpectraPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_InvBetaBTOFSpectraMinus[partIndex][centIndex]);
      if(m_loadedETOF) HistogramUtilities::ConditionalWrite(m_InvBetaETOFSpectraPlus[partIndex][centIndex]);
      if(m_loadedETOF) HistogramUtilities::ConditionalWrite(m_InvBetaETOFSpectraMinus[partIndex][centIndex]);
    }
  }

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
     HistogramUtilities::ConditionalWrite(m_BTOFEffPlus[partIndex][centIndex]);
     HistogramUtilities::ConditionalWrite(m_ETOFEffPlus[partIndex][centIndex]);
     HistogramUtilities::ConditionalWrite(m_BTOFEffMinus[partIndex][centIndex]);   
     HistogramUtilities::ConditionalWrite(m_ETOFEffMinus[partIndex][centIndex]);
     HistogramUtilities::ConditionalWrite(m_BTOFEffPlusRatio[partIndex][centIndex]);
     HistogramUtilities::ConditionalWrite(m_BTOFEffMinusRatio[partIndex][centIndex]);
    } // end cent
  } // end part

  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int centIndex = 0; centIndex < 16; centIndex++){
      HistogramUtilities::ConditionalWrite(m_btofScaleHistoPlus[partIndex][centIndex]);
      HistogramUtilities::ConditionalWrite(m_btofScaleHistoMinus[partIndex][centIndex]);
    }
  }

  outFile->mkdir("Spectra");
  

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    outFile->cd("Spectra");
    outFile->mkdir(m_partInfo->GetParticleName(partIndex,1).Data());
    outFile->cd(Form("Spectra/%s",m_partInfo->GetParticleName(partIndex,1).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      for(int rapIndex = 0; rapIndex < 100; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_spectraStruct[partIndex][centIndex][rapIndex][0],Form("SpectraStructure_%s_Cent%02d_RapIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,rapIndex));
      }// end rap
    } // end cent
  } // end part

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    outFile->cd("Spectra");
    outFile->mkdir(m_partInfo->GetParticleName(partIndex,-1).Data());
    outFile->cd(Form("Spectra/%s",m_partInfo->GetParticleName(partIndex,-1).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      for(int rapIndex = 0; rapIndex < 100; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_spectraStruct[partIndex][centIndex][rapIndex][1],Form("SpectraStructure_%s_Cent%02d_RapIndex%02d",m_partInfo->GetParticleName(partIndex,1).Data(),centIndex,rapIndex));
      }// end rap
    } // end cent
  } // end part


  //for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
  //  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
  //    m_TPCEfficiency_Nominal_Plus[partIndex][centIndex]->Write();
  //    m_TPCEfficiency_Systematic_Plus[partIndex][centIndex]->Write();
  //    m_TPCEfficiency_Nominal_Minus[partIndex][centIndex]->Write();
  //    m_TPCEfficiency_Systematic_Minus[partIndex][centIndex]->Write();
  //  }
  //
  //  m_EnergyLoss_Nominal_Plus[partIndex]->Write();
  //  m_EnergyLoss_Systematic_Plus[partIndex]->Write();
  //  m_EnergyLoss_Nominal_Minus[partIndex]->Write();
  //  m_EnergyLoss_Systematic_Minus[partIndex]->Write();
  //}



  outFile->Close();
  cout << "Wrote Member Data Successfully" << endl;
}


void RawSpectraModifier::makeBTofEff(string a_imageDir,double* a_momMatch,double* a_momStop, bool a_combinePlusMinus, bool a_usePlusEffForMinus){
  m_BTOF_Eff_Calc = true;
  gStyle->SetNumberContours(99);
  gSystem->mkdir(a_imageDir.c_str());
  cout << "Just made " << a_imageDir.c_str() << endl;

  //TCanvas* imgCanvas = new TCanvas();
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  
  TCanvas* effFitCanvas = new TCanvas();
  effFitCanvas->Divide(2,1);

  int chargeArr[2] = {-1,1};

  TH2D* tpcHisto = NULL;
  TH2D* tofHisto = NULL;
  TVirtualFitter::SetMaxIterations(10000000);
  int chargeIndex = 0;
  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    int chargeVal = chargeArr[pmIndex];
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      if(!m_particleLoaded[partIndex]) continue;
      cout << "############### " << "Part " << partIndex << "   Charge: " << pmIndex << "  #############" << endl;
      if(partIndex == 3) continue;
      //vector<TH2D* > tempCent; // temp vector to be added to m_TOFEffPlus...
      //vector<TH2D* > tempCentRatio; // temp vector to be added to m_TOFEffPlus...
  
      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        //if(pmIndex == 0){
        if(a_combinePlusMinus){
          tpcHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[partIndex][centIndex]);
          tofHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumBTOFTracks_Plus[partIndex][centIndex],m_rap_mTm0_NumBTOFTracks_Minus[partIndex][centIndex]);
        }else{
          if( chargeVal == -1 && !a_usePlusEffForMinus){
            tpcHisto = m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[partIndex][centIndex];
            tofHisto = m_rap_mTm0_NumBTOFTracks_Minus[partIndex][centIndex];
          }else{
            tpcHisto = m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex];
            tofHisto = m_rap_mTm0_NumBTOFTracks_Plus[partIndex][centIndex];
          }
        }
  
        //reduce the TPC Tracks to the BTOF Histo's Space (They should already have the same binning, but just in case)
        TH2D* tpcReduced = (TH2D*) tofHisto->Clone();
        HistogramUtilities::zeroBins2D(tpcReduced);
        for(int binx = 1; binx <= tpcReduced->GetNbinsX(); binx++ ){
          int origBin = tpcHisto->GetXaxis()->FindBin(tpcReduced->GetXaxis()->GetBinCenter(binx));
          for(int biny = 1; biny <= tpcReduced->GetNbinsY(); biny++ ){
            tpcReduced->SetBinContent(binx, biny, tpcHisto->GetBinContent(origBin,biny));
          }
        }
        HistogramUtilities::makeSqrtNErrors2D(tpcReduced);
        HistogramUtilities::makeSqrtNErrors2D(tofHisto);
  
        TH2D* ratioHisto = (TH2D*) tpcReduced->Clone();
        ratioHisto->SetTitle(Form("BTOF Efficiency %s", m_partInfo->GetParticleName(partIndex,chargeVal).Data()));
        ratioHisto->SetName(Form("BTofEff_%s_Cent%02d", m_partInfo->GetParticleName(partIndex,chargeVal).Data(),centIndex));
        ratioHisto->Divide(tofHisto,tpcReduced,1.0,1.0,"B");//bayesian

        //clean up nonsence valued bins
        for(int horBin = 1; horBin <= ratioHisto->GetNbinsX() ; horBin++){
          for(int vertBin = 1; vertBin <= ratioHisto->GetNbinsY(); vertBin++){
            double content = ratioHisto->GetBinContent(horBin,vertBin);
            double error = ratioHisto->GetBinError(horBin,vertBin);
            if(std::isnan(content) || content < 0.0) ratioHisto->SetBinContent(horBin,vertBin,0);
            if(std::isnan(error)) ratioHisto->SetBinError(horBin,vertBin,1e50);
          }
        }
        if(chargeVal == -1){
          m_BTOFEffMinusRatio[partIndex][centIndex] = ratioHisto;
        }else{
          m_BTOFEffPlusRatio[partIndex][centIndex]  = ratioHisto;
        }
      }// end cent loop
  
    }//end part
  }

}


double RawSpectraModifier::btofEffOperator(double* a_x, const double* a_params){
  if(a_x[0] > m_btofLowMtM0Ignore && a_x[0] < m_btofHighMtM0Ignore){
    return m_currentBTOF_Eff_ToFit->GetBinContent(m_currentBTOF_Eff_ToFit->FindBin(a_x[0]));
  }
  return a_params[0]*exp(a_params[1]*pow(a_x[0],-a_params[3])) + a_params[2] + a_params[4]*a_x[0];
}




void RawSpectraModifier::fitBTofEff(string a_imageDir, double* a_momMatch, double* a_momStop, bool a_reduceFitRange, vector< int > a_fitTags ){
  m_BTOF_Eff_Calc = true;
  gStyle->SetNumberContours(99);
  gSystem->mkdir(a_imageDir.c_str());
  cout << "Just made " << a_imageDir.c_str() << endl;

  //m_BTOFEffPlus.clear();
  //m_BTOFEffMinus.clear();

  //TCanvas* imgCanvas = new TCanvas();
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  TCanvas* effFitCanvas = new TCanvas();
  effFitCanvas->Divide(2,1);

  int charge[2] = {1,-1};

  TH2D* tpcHisto = NULL;
  TH2D* tofHisto = NULL;
  TVirtualFitter::SetMaxIterations(10000000);
  for(int chargeIndex = 0; chargeIndex < 2; chargeIndex++){
    double chargeVal = charge[chargeIndex];
    for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
      if(!m_particleLoaded[partIndex]) continue;
      cout << "############### " << m_partInfo->GetParticleName(partIndex,chargeVal).Data() << "  #############" << endl;
      if(partIndex == 3) continue;

      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

        if(chargeVal == 1){
          if(a_fitTags[partIndex] == 1){
            tpcHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[partIndex][centIndex]);
            tofHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumBTOFTracks_Plus[partIndex][centIndex],m_rap_mTm0_NumBTOFTracks_Minus[partIndex][centIndex]);
          }else{
            tpcHisto = m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex];
            tofHisto = m_rap_mTm0_NumBTOFTracks_Plus[partIndex][centIndex];
          }
        }else{ // negative particles
          if(a_fitTags[partIndex] == 1){
            tpcHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumTPCTracks_ForBTOF_Plus[partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[partIndex][centIndex]);
            tofHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumBTOFTracks_Plus[partIndex][centIndex],m_rap_mTm0_NumBTOFTracks_Minus[partIndex][centIndex]);
          }else if(a_fitTags[partIndex] == 0){
            tpcHisto = m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[partIndex][centIndex];
            tofHisto = m_rap_mTm0_NumBTOFTracks_Minus[partIndex][centIndex];
          }else if(a_fitTags[partIndex] == 2){
            tpcHisto = m_rap_mTm0_NumTPCTracks_ForBTOF_Minus[partIndex][centIndex];
            tofHisto = m_rap_mTm0_NumBTOFTracks_Minus[partIndex][centIndex];

            m_BTOFEffPlus[partIndex][centIndex] = dynamic_cast<TH2D*>(m_BTOFEffPlus[partIndex][centIndex]->Clone() );
          }
        }

        HistogramUtilities::makeSqrtNErrors2D(tpcHisto);
        HistogramUtilities::makeSqrtNErrors2D(tofHisto);
        TH2D* ratioHisto = (TH2D*) tpcHisto->Clone();
        ratioHisto->SetTitle(Form("BTOF Efficiency %s", m_partInfo->GetParticleName(partIndex,chargeVal).Data()));
        ratioHisto->SetName(Form("BTofEff_%s", m_partInfo->GetParticleName(partIndex,chargeVal).Data()));
        ratioHisto->Divide(tofHisto,tpcHisto,1.0,1.0,"B");//bayesian

        //smooth histogram (filled with fit data)
        TH2D* effHistoFit = (TH2D*) ratioHisto->Clone();
        HistogramUtilities::zeroBins2D(effHistoFit);
        effHistoFit->SetName(Form("BTofEff_%s_Fit", m_partInfo->GetParticleName(partIndex,chargeVal).Data()));
        effHistoFit->SetTitle(Form("BTOF Efficiency %s Fit", m_partInfo->GetParticleName(partIndex,chargeVal).Data()));
        cout <<  effHistoFit->GetTitle() << endl;

        for(int horBin = 1; horBin <= ratioHisto->GetNbinsX() ; horBin++){
          for(int vertBin = 1; vertBin <= ratioHisto->GetNbinsY(); vertBin++){
            double content = ratioHisto->GetBinContent(horBin,vertBin);
            double error = ratioHisto->GetBinError(horBin,vertBin);
            if(std::isnan(content)) ratioHisto->SetBinContent(horBin,vertBin,0);
            if(std::isnan(error)) ratioHisto->SetBinError(horBin,vertBin,1e50);
          }
        } 

        for(int horBin = 1; horBin <= effHistoFit->GetNbinsX() ; horBin++){
          double rapidity = effHistoFit->GetXaxis()->GetBinCenter(horBin);
          double tofCounts = 0;
          for(int vertBin = 1; vertBin <= tofHisto->GetNbinsY(); vertBin++){
            tofCounts += tofHisto->GetBinContent( horBin, vertBin);
          }
          if(tofCounts < 100) continue;
          cout << "   WORKING horBin " << horBin << endl;

          //double mTm0Low = a_momMatch[partIndex] / TMath::CosH(tofHisto->GetXaxis()->GetBinCenter(horBin));
          //double mTm0High = a_momStop[partIndex] / TMath::CosH(tofHisto->GetXaxis()->GetBinCenter(horBin));

          double mTm0Low = PhysMath::mTm0_from_y_p(tofHisto->GetXaxis()->GetBinCenter(horBin), 
            a_momMatch[partIndex], m_partInfo->GetParticleMass(partIndex));
          double mTm0High = PhysMath::mTm0_from_y_p(tofHisto->GetXaxis()->GetBinCenter(horBin), 
            a_momStop[partIndex], m_partInfo->GetParticleMass(partIndex));

          mTm0Low = mTm0Low - 0.1;
          mTm0High = mTm0High + 0.1;

          m_currentBTOF_Eff_ToFit = ratioHisto->ProjectionY(Form("eff_%s_%02d_%02d",
            m_partInfo->GetParticleName(partIndex,chargeVal).Data(),centIndex,horBin),horBin,horBin,"e");

          double x1 = m_currentBTOF_Eff_ToFit->GetBinCenter( m_currentBTOF_Eff_ToFit->FindBin(mTm0Low));
          double y1 = m_currentBTOF_Eff_ToFit->GetBinContent( m_currentBTOF_Eff_ToFit->FindBin(mTm0Low));
          double x2 = m_currentBTOF_Eff_ToFit->GetBinCenter( m_currentBTOF_Eff_ToFit->FindBin(mTm0High));
          double y2 = m_currentBTOF_Eff_ToFit->GetBinContent( m_currentBTOF_Eff_ToFit->FindBin(mTm0High));

          double slope = (y2 - y1)/(x2 - x1);
          double intercept =y1 -slope*x1;

          double chindf[3] = {5000,5000,5000};
          TF1* fitList[3];
          TFitResultPtr fitResultPtr[3];
          gStyle->SetOptStat(0);
          gStyle->SetOptFit(1111);
          gStyle->SetOptStat(0);

          if(!a_reduceFitRange){
            mTm0Low = 0.0;
            mTm0High = 5.0;
          }


          //TF1* chrisFunct = new TF1("chrisFunct", "[0]*exp([1]*pow(x,-[3]))+[2] + [4]*x",mTm0Low,mTm0High);
          TF1* chrisFunct = new TF1("chrisFunct",this,&RawSpectraModifier::btofEffOperator,0.0,3.5,5,"RawSpectraModifier","ZTPCOperator");
          chrisFunct->SetNpx(1200);

          //chrisFunct->SetParameter(0,0.6);
          //chrisFunct->SetParameter(1,-4.0e-5);
          //chrisFunct->SetParameter(2,0.0);
          //chrisFunct->SetParameter(3,-0.02);
          //chrisFunct->SetParameter(4,0.0);
          chrisFunct->SetParameter(0,0.001);
          chrisFunct->SetParLimits(0,0.0,20.0);
          chrisFunct->SetParameter(1,-0.112);
          chrisFunct->SetParLimits(1,-5.0,0.0);
          chrisFunct->SetParameter(2,0.0);
          chrisFunct->SetParLimits(2,-5.0, 5.0);
          chrisFunct->SetParameter(3,2.5);
          chrisFunct->SetParLimits(3,0.0,20.0);
          chrisFunct->SetParameter(4,0.0);
          chrisFunct->SetParLimits(4,-0.5,0.5);
          fittingCanvas->cd();
          m_currentBTOF_Eff_ToFit->Draw("E");
          gPad->SetRightMargin(0.4);
          //fitResultPtr[1] = m_currentBTOF_Eff_ToFit->Fit("chrisFunct","SQ","",mTm0Low,mTm0High);
          fitResultPtr[1] = m_currentBTOF_Eff_ToFit->Fit("chrisFunct","SQ","",0.0,3.5);
          cout << "IsEmpty? " << fitResultPtr[1]->IsEmpty() << endl;
          if(fitResultPtr[1]->IsEmpty()) continue;          

          cout << "-- FINAL STATUS -- " << fitResultPtr[1]->Status() << endl;
          fitList[1] = m_currentBTOF_Eff_ToFit->GetFunction("chrisFunct");
          TF1* totalFunct = new TF1("totalFunct", "[0]*exp([1]*pow(x,-[3]))+[2] + [4]*x",0.0,3.5);
          totalFunct->SetParameter(0,fitList[1]->GetParameter(0));
          totalFunct->SetParameter(1,fitList[1]->GetParameter(1));
          totalFunct->SetParameter(2,fitList[1]->GetParameter(2));
          totalFunct->SetParameter(3,fitList[1]->GetParameter(3));
          totalFunct->SetParameter(4,fitList[1]->GetParameter(4));
          totalFunct->SetLineColor(kRed);

          if(fitList[1]){
            chindf[1] = fitList[1]->GetChisquare()/((double) fitList[1]->GetNDF());
            cout << "Eval: " << fitList[1]->Eval(0.9) << endl;
          }else{
            continue;
          }
          fitList[1]->SetLineColor(kGreen);
          fitList[1]->Draw("same");
          totalFunct->Draw("same");



          // Determine a range of particle contamination 
          int index_low = partIndex - 1;
          int index_high = partIndex + 1;
          if(partIndex == 2) index_high = 4;
          if(partIndex == 4) index_low = 2;
          double lowMomCont = 50.0;
          if(partIndex > 0) lowMomCont = m_partInfo->momentumOfContamination(partIndex, index_low, 2.0);
          double highMomCont = m_partInfo->momentumOfContamination(partIndex, index_high, 2.0);

          double lowIndexMtM0Cont = PhysMath::mTm0_from_y_p(rapidity, lowMomCont, m_partInfo->GetParticleMass(partIndex));
          double highIndexMtM0Cont = PhysMath::mTm0_from_y_p(rapidity, highMomCont, m_partInfo->GetParticleMass(partIndex));
          m_btofLowMtM0Ignore = 50.0;
          if(lowIndexMtM0Cont > highIndexMtM0Cont && highIndexMtM0Cont > 0.0 ){
            m_btofLowMtM0Ignore = highIndexMtM0Cont;
          }else{
            if(lowIndexMtM0Cont > 0.0) m_btofLowMtM0Ignore = lowIndexMtM0Cont;
          }

          m_btofHighMtM0Ignore = 1.2;

          cout << "Current Ignore region, mT-m0 != [" << m_btofLowMtM0Ignore << "," << m_btofHighMtM0Ignore << "]"<< endl;

          TLine* lowIndexLine = new TLine(lowIndexMtM0Cont,0.0,lowIndexMtM0Cont,1.0);
          TLine* highIndexLine = new TLine(highIndexMtM0Cont,0.0,highIndexMtM0Cont,1.0);
          lowIndexLine->SetLineColor(kViolet);
          highIndexLine->SetLineColor(kViolet+1);
          //lowIndexLine->Draw("same");
          //highIndexLine->Draw("same");

          TLine* lowRangeLine = new TLine(m_btofLowMtM0Ignore,0.0,m_btofLowMtM0Ignore,1.0);
          TLine* highRangeLine = new TLine(m_btofHighMtM0Ignore,0.0,m_btofHighMtM0Ignore,1.0);
          lowRangeLine->SetLineColor(kRed);
          highRangeLine->SetLineColor(kRed+1);
          lowRangeLine->Draw("same");
          highRangeLine->Draw("same");

          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/BTOF_%s_Rap_%02d_Cent_%02d_ChrisFunct.png",
            a_imageDir.c_str(),m_partInfo->GetParticleName(
              partIndex,chargeVal).Data(), horBin,centIndex));

          if(!(fitResultPtr[1].Get()))continue;
          if(fitResultPtr[1].Get()->IsEmpty()) continue;
          cout << "IsEmpty? " << fitResultPtr[1]->IsEmpty() << endl;
          if(fitResultPtr[1]->IsEmpty()) continue;
          int fitOfSmall = 1;
          double binCenters[(int) effHistoFit->GetNbinsY()];
          for(int vertBin = 1; vertBin <= effHistoFit->GetNbinsY(); vertBin++){
            binCenters[vertBin-1] = effHistoFit->GetYaxis()->GetBinCenter(vertBin);
          }
          double effErrors[(int) effHistoFit->GetNbinsY()];
          fitResultPtr[fitOfSmall]->GetConfidenceIntervals((int) effHistoFit->GetNbinsY(), 1, 1, binCenters, effErrors, 0.683, false);
          cout << "Size " << effHistoFit->GetNbinsY() << endl;
          for(int vertBin = 1; vertBin <= effHistoFit->GetNbinsY(); vertBin++){
            double mtcent = effHistoFit->GetYaxis()->GetBinCenter(vertBin);
              //int vertIndex = vertBin - 1;
            if(mtcent < mTm0Low || mTm0High < mtcent) continue;
              //cout << "Eval at " << effHistoFit->GetYaxis()->GetBinCenter(vertBin) << " between " << mTm0Low << " and " << mTm0High << endl;
              //cout << "Funct Ptr: " << fitList[fitOfSmall] << endl;
              //cout << "   par 0 " << fitList[fitOfSmall]->GetParameter(0) << endl;
            double effVal = fitList[fitOfSmall]->Eval(effHistoFit->GetYaxis()->GetBinCenter(vertBin));
            effHistoFit->SetBinContent(horBin,vertBin, effVal);
            if(effErrors[vertBin-1] > ratioHisto->GetBinError(horBin,vertBin)) effErrors[vertBin-1] = ratioHisto->GetBinError(horBin,vertBin);
            if(std::isnan(effErrors[vertBin-1])) effErrors[vertBin-1] = ratioHisto->GetBinError(horBin,vertBin);
            if(effErrors[vertBin-1] <= 0.0) {effErrors[vertBin-1] = ratioHisto->GetBinError(horBin,vertBin);}
            if(effErrors[vertBin-1] > 100.0) effErrors[vertBin-1] = 5.0;
 
            effHistoFit->SetBinError(horBin,vertBin,effErrors[vertBin-1]);
             //cout << vertBin << " " << effErrors[vertBin-1] << endl;
          }
        }
        effFitCanvas->cd(1);
        gPad->SetLogz();
        ratioHisto->GetZaxis()->SetRangeUser(0.001,1.0);
        ratioHisto->Draw("colz");
   
        effFitCanvas->cd(2);
        gPad->SetLogz();
        effHistoFit->GetZaxis()->SetRangeUser(0.001,1.0);
        effHistoFit->Draw("colz");
        effFitCanvas->SaveAs(Form("%s/BTOF_%s_Cent_%02d.png",
           a_imageDir.c_str(), m_partInfo->GetParticleName(partIndex,chargeVal).Data(),centIndex));
        if(chargeVal == -1){
          m_BTOFEffMinus[partIndex][centIndex] = effHistoFit;
        }else{
          m_BTOFEffPlus[partIndex][centIndex] = effHistoFit;
        }
      }// end cent loop
    }//end particle loop
  } // end charge loop

}





void RawSpectraModifier::fitETofEff(string a_imageDir,double* a_momMatch,double* a_momStop){
  m_ETOF_Eff_Calc = true;
  gStyle->SetNumberContours(99);
  gSystem->mkdir(a_imageDir.c_str());

  //TCanvas* imgCanvas = new TCanvas();
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  TCanvas* effFitCanvas = new TCanvas();
  effFitCanvas->Divide(2,1);

  int charge[2] = {-1,1};

  TH2D* tpcHisto = NULL;
  TH2D* tofHisto = NULL;
  TVirtualFitter::SetMaxIterations(10000000);
  int chargeIndex = 0;
  //for(int chargeIndex = 0; chargeIndex < 2; chargeIndex++){
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    cout << "############### " << "Part " << partIndex << "   Charge: " << chargeIndex << "  #############" << endl;
      //if(partIndex !=2) continue;
    if(partIndex == 3) continue;
      //vector<TH2D* > tempCent; // temp vector to be added to m_TOFEffPlus...

      for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
        //if(chargeIndex == 0){
        //  tpcHisto = m_rap_mTm0_NumTPCTracks_Minus[partIndex][centIndex];
        //  tofHisto = m_rap_mTm0_NumETOFTracks_Minus[partIndex][centIndex];
        //}else{
        //  tpcHisto = m_rap_mTm0_NumTPCTracks_Plus[partIndex][centIndex];
        //  tofHisto = m_rap_mTm0_NumETOFTracks_Plus[partIndex][centIndex];
        //}
        tpcHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumTPCTracks_ForETOF_Minus[partIndex][centIndex],m_rap_mTm0_NumTPCTracks_ForETOF_Plus[partIndex][centIndex]);
        tofHisto = HistogramUtilities::addTH2D(m_rap_mTm0_NumETOFTracks_Minus[partIndex][centIndex],m_rap_mTm0_NumETOFTracks_Plus[partIndex][centIndex]);
        //reduce the TPC Tracks to the ETOF Eta Space
        TH2D* tpcReduced = (TH2D*) tofHisto->Clone();
        HistogramUtilities::zeroBins2D(tpcReduced);
        for(int binx = 1; binx <= tpcReduced->GetNbinsX(); binx++ ){
          int origBin = tpcHisto->FindBin(tpcReduced->GetXaxis()->GetBinCenter(binx));
          for(int biny = 1; biny <= tpcReduced->GetNbinsY(); biny++ ){
            tpcReduced->SetBinContent(binx, biny, tpcHisto->GetBinContent(origBin,biny));
          }
        }
        HistogramUtilities::makeSqrtNErrors2D(tpcReduced);
        HistogramUtilities::makeSqrtNErrors2D(tofHisto);
        TH2D* ratioHisto = (TH2D*) tpcReduced->Clone();
        ratioHisto->SetTitle(Form("ETOF Efficiency %s", m_partInfo->GetParticleName(partIndex,charge[chargeIndex]).Data()));
        ratioHisto->SetName(Form("ETofEff_%s", m_partInfo->GetParticleName(partIndex,charge[chargeIndex]).Data()));
        ratioHisto->Divide(tofHisto,tpcReduced,1.0,1.0,"B");//bayesian


        TH2D* effHistoFit = (TH2D*) ratioHisto->Clone();
        HistogramUtilities::zeroBins2D(effHistoFit);
        effHistoFit->SetName(Form("ETofEff_%s_Fit", m_partInfo->GetParticleName(partIndex,charge[chargeIndex]).Data()));
        effHistoFit->SetTitle(Form("ETOF Efficiency %s Fit", m_partInfo->GetParticleName(partIndex,charge[chargeIndex]).Data()));

        for(int horIndex = 1; horIndex <= effHistoFit->GetNbinsX() ; horIndex++){
          double tofCounts = 0;
          for(int vertBin = 1; vertBin <= tofHisto->GetNbinsY(); vertBin++){
            tofCounts += tofHisto->GetBinContent( horIndex, vertBin);
          }
          if(tofCounts < 100) continue;
          cout << "   WORKING HorIndex " << horIndex << endl;

         //double mTm0Low = a_momMatch[partIndex] / TMath::CosH(tofHisto->GetXaxis()->GetBinCenter(horIndex));
         //double mTm0High = a_momStop[partIndex] / TMath::CosH(tofHisto->GetXaxis()->GetBinCenter(horIndex));

          double mTm0Low = PhysMath::mTm0_from_y_p(tofHisto->GetXaxis()->GetBinCenter(horIndex), 
            a_momMatch[partIndex], m_partInfo->GetParticleMass(partIndex));
          double mTm0High = PhysMath::mTm0_from_y_p(tofHisto->GetXaxis()->GetBinCenter(horIndex), 
            a_momStop[partIndex], m_partInfo->GetParticleMass(partIndex));

          mTm0Low = mTm0Low - 0.1;
          mTm0High = mTm0High + 0.1;

          TH1D* effCurve = ratioHisto->ProjectionY(Form("eff_%s_%02d_%02d",
            m_partInfo->GetParticleName(partIndex,charge[chargeIndex]).Data(),centIndex,horIndex),horIndex,horIndex,"e");

          double x1 = effCurve->GetBinCenter( effCurve->FindBin(mTm0Low));
          double y1 = effCurve->GetBinContent( effCurve->FindBin(mTm0Low));
          double x2 = effCurve->GetBinCenter( effCurve->FindBin(mTm0High));
          double y2 = effCurve->GetBinContent( effCurve->FindBin(mTm0High));

          double slope = (y2 - y1)/(x2 - x1);
          double intercept =y1 -slope*x1;

          double chindf[3] = {5000,5000,5000};
          TF1* fitList[3];
          TFitResultPtr fitResultPtr[3];
          gStyle->SetOptStat(0);
          gStyle->SetOptFit(1111);
          gStyle->SetOptStat(0);
           /*
          //TF1* funct = new TF1("errFit", "[0]*TMath::Erf((x-[1])/[2])",0,3.0);
          TF1* lineFunct = new TF1("lineFunct", "[0]*x+[1]",mTm0Low,mTm0High);
          lineFunct->SetParameter(0,slope);
          lineFunct->SetParLimits(0,slope-0.2,slope+0.2);
          lineFunct->SetParameter(1,intercept);
          lineFunct->SetParLimits(1,intercept-0.5,intercept+0.5);
          fittingCanvas->cd();
          effCurve->Draw("E");
          gPad->SetRightMargin(0.4);
          fitResultPtr[0] = effCurve->Fit("lineFunct","SE","",mTm0Low,mTm0High);
          cout << "-- FINAL STATUS -- " << fitResultPtr[0]->Status() << endl;
          fitList[0] = effCurve->GetFunction("lineFunct");          
          if(fitList[0]){
            chindf[0] = fitList[0]->GetChisquare()/((double) fitList[0]->GetNDF());
          }
          fitList[0]->SetLineColor(kGreen);
          fitList[0]->Draw("same");
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s_Rap_%02d_Cent_%02d_Linear.png",
              a_imageDir.c_str(),m_partInfo->GetParticleName(
              partIndex,charge[chargeIndex]).Data(), horIndex,centIndex)); 
          */



          TF1* chrisFunct = new TF1("chrisFunct", "[0]*exp([1]*pow(x,-[3]))+[2] + [4]*x",mTm0Low,mTm0High);
          chrisFunct->SetParameter(0,0.001);
          chrisFunct->SetParameter(1,-0.112);
          chrisFunct->SetParameter(2,intercept);
          chrisFunct->SetParameter(3,2.5);
          chrisFunct->SetParameter(4,slope);
          fittingCanvas->cd();
          effCurve->Draw("E");
          gPad->SetRightMargin(0.4);
          fitResultPtr[1] = effCurve->Fit("chrisFunct","SE","",mTm0Low,mTm0High);
          cout << "IsEmpty? " << fitResultPtr[1]->IsEmpty() << endl;
          if(fitResultPtr[1]->IsEmpty()) continue;          

          cout << "-- FINAL STATUS -- " << fitResultPtr[1]->Status() << endl;
          fitList[1] = effCurve->GetFunction("chrisFunct");
          if(fitList[1]){
            chindf[1] = fitList[1]->GetChisquare()/((double) fitList[1]->GetNDF());
            cout << "Eval: " << fitList[1]->Eval(0.9) << endl;
          }
          fitList[1]->SetLineColor(kGreen);
          fitList[1]->Draw("same");
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/ETOF_%s_Eta_%02d_Cent_%02d_Chris.png",
            a_imageDir.c_str(),m_partInfo->GetParticleName(
              partIndex,charge[chargeIndex]).Data(), horIndex,centIndex));
          
          /*
          TF1* polyFunct = new TF1("polyFunct", "[0]*log([1]*x+[2])+[3]",mTm0Low,mTm0High);
          polyFunct->SetParameter(0,0.5);
          polyFunct->SetParameter(1,1.0);
          polyFunct->SetParameter(2,0.5);
          polyFunct->SetParameter(3,intercept);
          fittingCanvas->cd();
          effCurve->Draw("E");
          gPad->SetRightMargin(0.4);
          fitResultPtr[2] = effCurve->Fit("polyFunct","SE","",mTm0Low,mTm0High);
          cout << "-- FINAL STATUS -- " << fitResultPtr[2]->Status() << endl;
          fitList[2] = effCurve->GetFunction("polyFunct");
          if(fitList[2]){
            chindf[2] = fitList[2]->GetChisquare()/((double) fitList[2]->GetNDF());
          }
          fitList[2]->SetLineColor(kGreen);
          fitList[2]->Draw("same");
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s_Rap_%02d_Cent_%02d_PolyLog.png",
              a_imageDir.c_str(),m_partInfo->GetParticleName(
              partIndex,charge[chargeIndex]).Data(), horIndex,centIndex));


          //status = fitResultPtr->Status();
          //if(!(status == 4000 || status == 0)){
          //  cout << "ERROR in Efficiency fits" << endl;
          //  cout << "Status code: " << status << endl;
                  //continue;
          //}
          //

          double smallestChi = 1E7;
          int fitOfSmall = -1;
          for(int fitIndex = 0; fitIndex < 3; fitIndex++){
            //chindf[fitIndex] = fitList[fitIndex]->GetChisquare()/((double) fitList[fitIndex]->GetNDF());
            cout << "Chi^2 for " << fitIndex << " = " << chindf[fitIndex] << endl;
            if(chindf[fitIndex] < smallestChi && fitList[fitIndex] && fitResultPtr[fitIndex]->Status() == 0){
              smallestChi = chindf[fitIndex];
              fitOfSmall = fitIndex;
            }
          }
          cout << "Smallest: " << fitOfSmall << " with " << smallestChi << endl;
          fitOfSmall = 2;
          cout << "Ptrs: " << fitResultPtr[0].Get() << " " << fitResultPtr[1].Get() << " " << fitResultPtr[2].Get() << endl;
          cout << "Functs: " << fitList[0] << " " << fitList[1] << " " << fitList[2] << endl;
          */

          if(!(fitResultPtr[1].Get()))continue;
          if(fitResultPtr[1].Get()->IsEmpty()) continue;
          cout << "IsEmpty? " << fitResultPtr[1]->IsEmpty() << endl;
          if(fitResultPtr[1]->IsEmpty()) continue;
          int fitOfSmall = 1;
          double binCenters[(int) effHistoFit->GetNbinsY()];
          for(int vertBin = 1; vertBin <= effHistoFit->GetNbinsY(); vertBin++){
            binCenters[vertBin-1] = effHistoFit->GetYaxis()->GetBinCenter(vertBin);
          }
          double effErrors[(int) effHistoFit->GetNbinsY()];
          fitResultPtr[fitOfSmall]->GetConfidenceIntervals((int) effHistoFit->GetNbinsY(), 1, 1, binCenters, effErrors, 0.683, false);
          cout << "Size " << effHistoFit->GetNbinsY() << endl;
          for(int vertBin = 1; vertBin <= effHistoFit->GetNbinsY(); vertBin++){
           double mtcent = effHistoFit->GetYaxis()->GetBinCenter(vertBin);
             //int vertIndex = vertBin - 1;
           if(mtcent < mTm0Low || mTm0High < mtcent) continue;
             //cout << "Eval at " << effHistoFit->GetYaxis()->GetBinCenter(vertBin) << " between " << mTm0Low << " and " << mTm0High << endl;
             //cout << "Funct Ptr: " << fitList[fitOfSmall] << endl;
             //cout << "   par 0 " << fitList[fitOfSmall]->GetParameter(0) << endl;
           double effVal = fitList[fitOfSmall]->Eval(effHistoFit->GetYaxis()->GetBinCenter(vertBin));
           effHistoFit->SetBinContent(horIndex,vertBin, effVal);
           if(effErrors[vertBin-1] > ratioHisto->GetBinError(horIndex,vertBin)) effErrors[vertBin-1] = ratioHisto->GetBinError(horIndex,vertBin);
           if(std::isnan(effErrors[vertBin-1])) effErrors[vertBin-1] = ratioHisto->GetBinError(horIndex,vertBin);
           if(effErrors[vertBin-1] <= 0.0) {effErrors[vertBin-1] = ratioHisto->GetBinError(horIndex,vertBin);}
           if(effErrors[vertBin-1] > 100.0) effErrors[vertBin-1] = 5.0;

           effHistoFit->SetBinError(horIndex,vertBin,effErrors[vertBin-1]);
             //cout << vertBin << " " << effErrors[vertBin-1] << endl;
         }


           /*for(int vertBin = 0; vertBin <= effHistoFit->GetNbinsY(); vertBin++){
             double lowEdge = effHistoFit->GetYaxis()->GetBinLowEdge(vertBin);
             double highEdge = effHistoFit->GetYaxis()->GetBinLowEdge(vertBin+1);
             double lowVal = fitResult->Eval(lowEdge);
             double highVal = fitResult->Eval(highEdge);
             effHistoFit->SetBinContent(horIndex,vertBin, (lowVal+highVal)*0.5);
             effHistoFit->SetBinError(horIndex,vertBin, ratioHisto->GetBinError(horIndex,vertBin));
           }*/


       }
       effFitCanvas->cd(1);
       gPad->SetLogz();
       ratioHisto->GetZaxis()->SetRangeUser(0.001,1.0);
       ratioHisto->Draw("colz");

       effFitCanvas->cd(2);
       gPad->SetLogz();
       effHistoFit->GetZaxis()->SetRangeUser(0.001,1.0);
       effHistoFit->Draw("colz");
       effFitCanvas->SaveAs(Form("%s/ETOF_%s_Cent_%02d.png",
        a_imageDir.c_str(), m_partInfo->GetParticleName(
          partIndex,charge[chargeIndex]).Data(),centIndex));

       //if(chargeVal == -1){
          m_ETOFEffMinus[partIndex][centIndex] = effHistoFit;
        //}else{
          m_ETOFEffPlus[partIndex][centIndex] = effHistoFit;
        //}

      }// end cent loop

    }//end particle loop
  //}//end charge loop


  }



//########################################################################################################################
void RawSpectraModifier::makeRapiditySpectraTGraphs(){
  #ifdef _RAWSPECMOD_DEBUG_
    cout << "Call to RawSpectraModifier::makeRapiditySpectraTGraphs..." << endl;
  #endif

  double simple16CentEdges[17] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80};
  double simple9CentEdges[10] = {0,5,10,20,30,40,50,60,70,80};

  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int centIndex = 0; centIndex < 16; centIndex++){
      for(int rapIndex = 0; rapIndex < g_RawSpecMod_MaxNumRapBins; rapIndex++){
        m_spectraPlusDeDx[partIndex][centIndex][rapIndex] = NULL;
        m_spectraMinusDeDx[partIndex][centIndex][rapIndex] = NULL;
        m_spectraPlusInvBeta[partIndex][centIndex][rapIndex] = NULL;
        m_spectraMinusInvBeta[partIndex][centIndex][rapIndex] = NULL;
        m_spectraPlusInvBetaETOF[partIndex][centIndex][rapIndex] = NULL;
        m_spectraMinusInvBetaETOF[partIndex][centIndex][rapIndex] = NULL;
        m_rawSpectraPlusDeDx[partIndex][centIndex][rapIndex] = NULL;
        m_rawSpectraMinusDeDx[partIndex][centIndex][rapIndex] = NULL;
        m_rawSpectraPlusInvBeta[partIndex][centIndex][rapIndex] = NULL;
        m_rawSpectraMinusInvBeta[partIndex][centIndex][rapIndex] = NULL;
        m_rawSpectraPlusInvBetaETOF[partIndex][centIndex][rapIndex] = NULL;
        m_rawSpectraMinusInvBetaETOF[partIndex][centIndex][rapIndex] = NULL;
        m_inclusiveSpectraPlusDeDx[partIndex][centIndex][rapIndex] = NULL;
        m_inclusiveSpectraMinusDeDx[partIndex][centIndex][rapIndex] = NULL;
        m_inclusiveSpectraPlusInvBeta[partIndex][centIndex][rapIndex] = NULL;
        m_inclusiveSpectraMinusInvBeta[partIndex][centIndex][rapIndex] = NULL;
        m_inclusiveSpectraPlusInvBetaETOF[partIndex][centIndex][rapIndex] = NULL;
        m_inclusiveSpectraMinusInvBetaETOF[partIndex][centIndex][rapIndex] = NULL;
      }
    }
  }

  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      int lowCentEdge = 0;
      int highCentEdge = 80;
      if(m_combined16To9 || m_numCentralities==9){
        lowCentEdge  = simple9CentEdges[centIndex];
        highCentEdge = simple9CentEdges[centIndex+1];
      }else{
        lowCentEdge = simple16CentEdges[centIndex];
        highCentEdge = simple16CentEdges[centIndex+1];
      }

      // TPC Spectra
      vector<TGraphErrors*> projections_dEdx_Plus = HistogramUtilities::VerticalProjectionGraphs(
        m_dEdxSpectraPlus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_dEdx_Minus = HistogramUtilities::VerticalProjectionGraphs(
        m_dEdxSpectraMinus[partIndex][centIndex], "y", "y", true);
      // bTOF Spectra
      vector<TGraphErrors*> projections_BTOF_Plus = HistogramUtilities::VerticalProjectionGraphs(
        m_InvBetaBTOFSpectraPlus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_BTOF_Minus = HistogramUtilities::VerticalProjectionGraphs(
        m_InvBetaBTOFSpectraMinus[partIndex][centIndex], "y", "y", true);
      // eTOF Spectra
      vector<TGraphErrors*> projections_ETOF_Plus;
      vector<TGraphErrors*> projections_ETOF_Minus;
      if(m_loadedETOF){
        projections_ETOF_Plus = HistogramUtilities::VerticalProjectionGraphs(
          m_InvBetaETOFSpectraPlus[partIndex][centIndex], "y", "y", true);
        projections_ETOF_Minus= HistogramUtilities::VerticalProjectionGraphs(
          m_InvBetaETOFSpectraMinus[partIndex][centIndex], "y", "y", true);
      }

      for(int rapIndex = 0; rapIndex < projections_dEdx_Plus.size(); rapIndex++){
        m_rawSpectraPlusDeDx[partIndex][centIndex][rapIndex]  = projections_dEdx_Plus[rapIndex];
        m_rawSpectraMinusDeDx[partIndex][centIndex][rapIndex] = projections_dEdx_Minus[rapIndex];
        m_rawSpectraPlusDeDx[partIndex][centIndex][rapIndex]->SetName( Form("raw_dEdxSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
        m_rawSpectraMinusDeDx[partIndex][centIndex][rapIndex]->SetName(Form("raw_dEdxSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));

        m_spectraPlusDeDx[partIndex][centIndex][rapIndex]     = (TGraphErrors*) projections_dEdx_Plus[rapIndex]->Clone();
        m_spectraMinusDeDx[partIndex][centIndex][rapIndex]    = (TGraphErrors*) projections_dEdx_Minus[rapIndex]->Clone();
        m_spectraPlusDeDx[partIndex][centIndex][rapIndex]->SetName(    Form("corrected_dEdxSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
        m_spectraMinusDeDx[partIndex][centIndex][rapIndex]->SetName(   Form("corrected_dEdxSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));

        #ifdef _RAWSPECMOD_DEBUG_
          if(fabs(m_dEdxSpectraPlus[partIndex][centIndex]->GetXaxis()->GetBinCenter(rapIndex+1)) < 0.001){
            cout << " Example Rapidity: 0.0     Num point raw: " << m_rawSpectraPlusDeDx[partIndex][centIndex][rapIndex]->GetN() << "     corrected: " << m_spectraPlusDeDx[partIndex][centIndex][rapIndex]->GetN() << endl;
          }
        #endif
      }
      for(int rapIndex = 0; rapIndex < projections_BTOF_Plus.size(); rapIndex++){
        m_rawSpectraPlusInvBeta[partIndex][centIndex][rapIndex]  = projections_BTOF_Plus[rapIndex];
        m_rawSpectraMinusInvBeta[partIndex][centIndex][rapIndex] = projections_BTOF_Minus[rapIndex];
        m_rawSpectraPlusInvBeta[partIndex][centIndex][rapIndex]->SetName(Form("raw_InvBetaBTOFSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
        m_rawSpectraMinusInvBeta[partIndex][centIndex][rapIndex]->SetName(Form("raw_InvBetaBTOFSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));

        m_spectraPlusInvBeta[partIndex][centIndex][rapIndex]     = (TGraphErrors*) projections_BTOF_Plus[rapIndex]->Clone();
        m_spectraMinusInvBeta[partIndex][centIndex][rapIndex]    = (TGraphErrors*) projections_BTOF_Minus[rapIndex]->Clone();
        m_spectraPlusInvBeta[partIndex][centIndex][rapIndex]->SetName(Form("corrected_InvBetaBTOFSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
        m_spectraMinusInvBeta[partIndex][centIndex][rapIndex]->SetName(Form("corrected_InvBetaBTOFSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
      }
      for(int rapIndex = 0;m_loadedETOF && rapIndex < projections_ETOF_Plus.size(); rapIndex++){
        m_rawSpectraPlusInvBetaETOF[partIndex][centIndex][rapIndex]  = projections_ETOF_Plus[rapIndex];
        m_rawSpectraMinusInvBetaETOF[partIndex][centIndex][rapIndex] = projections_ETOF_Minus[rapIndex];
        m_rawSpectraPlusInvBetaETOF[partIndex][centIndex][rapIndex]->SetName(Form("raw_InvBetaETOFSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
        m_rawSpectraMinusInvBetaETOF[partIndex][centIndex][rapIndex]->SetName(Form("raw_InvBetaETOFSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));

        m_spectraPlusInvBetaETOF[partIndex][centIndex][rapIndex]  = (TGraphErrors*) projections_ETOF_Plus[rapIndex]->Clone();
        m_spectraMinusInvBetaETOF[partIndex][centIndex][rapIndex] = (TGraphErrors*) projections_ETOF_Minus[rapIndex]->Clone();
        m_spectraPlusInvBetaETOF[partIndex][centIndex][rapIndex]->SetName(Form("corrected_InvBetaETOFSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
        m_spectraMinusInvBetaETOF[partIndex][centIndex][rapIndex]->SetName(Form("corrected_InvBetaETOFSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,rapIndex));
      }



      //############################   MAKE SPECTRA STRUCTURES #####################################
      vector<double> rapidityList_positiveParticles;
      vector<double> rapIndexList_positiveParticles;
      vector<double> rapidityList_negativeParticles;
      vector<double> rapIndexList_negativeParticles;

      //// //############################  TPC SPECTRA 
      for(int rapIndex = 0; rapIndex < projections_dEdx_Plus.size(); rapIndex++){
        double rapVal = m_dEdxSpectraPlus[partIndex][centIndex]->GetXaxis()->GetBinCenter(rapIndex+1);
        m_spectraStruct[partIndex][centIndex][rapIndex][0] = new SpectraStructure();
        m_spectraStruct[partIndex][centIndex][rapIndex][0]->Construct(partIndex, 1, rapIndex,
           rapVal,centIndex,lowCentEdge, highCentEdge);
        m_spectraStruct[partIndex][centIndex][rapIndex][0]->setRawTPCSpectra(projections_dEdx_Plus[rapIndex]);
        rapidityList_positiveParticles.push_back(rapVal);
        rapIndexList_positiveParticles.push_back(rapIndex);
      }
      for(int rapIndex = 0; rapIndex < projections_dEdx_Minus.size(); rapIndex++){
        double rapVal = m_dEdxSpectraMinus[partIndex][centIndex]->GetXaxis()->GetBinCenter(rapIndex+1);
        m_spectraStruct[partIndex][centIndex][rapIndex][1] = new SpectraStructure();
        m_spectraStruct[partIndex][centIndex][rapIndex][1]->Construct(partIndex, -1, rapIndex,
           rapVal,centIndex,lowCentEdge, highCentEdge);
        m_spectraStruct[partIndex][centIndex][rapIndex][1]->setRawTPCSpectra(projections_dEdx_Minus[rapIndex]);
        rapidityList_negativeParticles.push_back(rapVal);
        rapIndexList_negativeParticles.push_back(rapIndex);

      }
      //// //############################  BTOF SPECTRA PLUS

      for(int rapIndex = 0; rapIndex < projections_BTOF_Plus.size(); rapIndex++){
        double rapVal = m_InvBetaBTOFSpectraPlus[partIndex][centIndex]->GetXaxis()->GetBinCenter(rapIndex+1);
        int rapIndexOfAlreadyMade = -1;
        bool spectraStructAlreadyMade = false;
        for(int iii = 0; iii < rapidityList_positiveParticles.size(); iii++){
          if(fabs(rapVal - rapidityList_positiveParticles[iii]) < 0.005){
            spectraStructAlreadyMade = true;
            rapIndexOfAlreadyMade = rapIndexList_positiveParticles[iii];
          }
        }

        if(spectraStructAlreadyMade){
          m_spectraStruct[partIndex][centIndex][rapIndexOfAlreadyMade][0]->setRawBTOFSpectra(projections_BTOF_Plus[rapIndex]);
        }else{
          int nextNullRapIndex = -1;
          for(int iii = 0; iii < g_RawSpecMod_MaxNumRapBins; iii++){
            if(!m_spectraStruct[partIndex][centIndex][rapIndex][0]){
              nextNullRapIndex = iii; break;
            }
          }
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][0] = new SpectraStructure();
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][0]->Construct(partIndex, 1, nextNullRapIndex,
             rapVal,centIndex,lowCentEdge, highCentEdge);
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][0]->setRawBTOFSpectra(projections_BTOF_Plus[rapIndex]);
          rapidityList_positiveParticles.push_back(rapVal);
          rapIndexList_positiveParticles.push_back(nextNullRapIndex);

        }
      }

      //// //############################  BTOF SPECTRA MINUS


      for(int rapIndex = 0; rapIndex < projections_BTOF_Minus.size(); rapIndex++){
        double rapVal = m_InvBetaBTOFSpectraMinus[partIndex][centIndex]->GetXaxis()->GetBinCenter(rapIndex+1);
        int rapIndexOfAlreadyMade = -1;
        bool spectraStructAlreadyMade = false;
        for(int iii = 0; iii < rapidityList_negativeParticles.size(); iii++){
          if(fabs(rapVal - rapidityList_negativeParticles[iii]) < 0.005){
            spectraStructAlreadyMade = true;
            rapIndexOfAlreadyMade = rapIndexList_negativeParticles[iii];
          }
        }


        #ifdef _RAWSPECMOD_DEBUG_
          cout << " BTOF Spectra has corresponding TPC Spectra? : " << spectraStructAlreadyMade << "   rap: " << rapVal << endl;
          cout << "    already made rapIndex: " << rapIndexOfAlreadyMade << "    " << m_spectraStruct[partIndex][centIndex][rapIndexOfAlreadyMade][1] << endl;
        #endif

        if(spectraStructAlreadyMade){
          m_spectraStruct[partIndex][centIndex][rapIndexOfAlreadyMade][1]->setRawBTOFSpectra(projections_BTOF_Minus[rapIndex]);
        }else{
          int nextNullRapIndex = -1;
          for(int iii = 0; iii < g_RawSpecMod_MaxNumRapBins; iii++){
            if(!m_spectraStruct[partIndex][centIndex][rapIndex][1]){
              nextNullRapIndex = iii; break;
            }
          }
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1] = new SpectraStructure();
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1]->Construct(partIndex, -1, nextNullRapIndex,
             rapVal,centIndex,lowCentEdge, highCentEdge);
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1]->setRawBTOFSpectra(projections_BTOF_Minus[rapIndex]);
          rapidityList_negativeParticles.push_back(rapVal);
          rapIndexList_negativeParticles.push_back(nextNullRapIndex);

        }
      }



      for(int rapIndex = 0; rapIndex < projections_ETOF_Plus.size(); rapIndex++){
        double rapVal = m_InvBetaETOFSpectraPlus[partIndex][centIndex]->GetXaxis()->GetBinCenter(rapIndex+1);
        int rapIndexOfAlreadyMade = -1;
        bool spectraStructAlreadyMade = false;
        for(int iii = 0; iii < rapidityList_positiveParticles.size(); iii++){
          if(fabs(rapVal - rapidityList_positiveParticles[iii]) < 0.005){
            spectraStructAlreadyMade = true;
            rapIndexOfAlreadyMade = rapIndexList_positiveParticles[iii];
          }
        }

        if(spectraStructAlreadyMade){
          m_spectraStruct[partIndex][centIndex][rapIndexOfAlreadyMade][1]->setRawETOFSpectra(projections_ETOF_Plus[rapIndex]);
        }else{
          int nextNullRapIndex = -1;
          for(int iii = 0; iii < g_RawSpecMod_MaxNumRapBins; iii++){
            if(!m_spectraStruct[partIndex][centIndex][rapIndex][1]){
              nextNullRapIndex = iii; break;
            }
          }
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1] = new SpectraStructure();
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1]->Construct(partIndex, -1, nextNullRapIndex,
             rapVal,centIndex,lowCentEdge, highCentEdge);
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1]->setRawETOFSpectra(projections_ETOF_Plus[rapIndex]);
          rapidityList_positiveParticles.push_back(rapVal);
          rapIndexList_positiveParticles.push_back(nextNullRapIndex);

        }
      }


      for(int rapIndex = 0; rapIndex < projections_ETOF_Minus.size(); rapIndex++){
        double rapVal = m_InvBetaETOFSpectraMinus[partIndex][centIndex]->GetXaxis()->GetBinCenter(rapIndex+1);
        int rapIndexOfAlreadyMade = -1;
        bool spectraStructAlreadyMade = false;
        for(int iii = 0; iii < rapidityList_negativeParticles.size(); iii++){
          if(fabs(rapVal - rapidityList_negativeParticles[iii]) < 0.005){
            spectraStructAlreadyMade = true;
            rapIndexOfAlreadyMade = rapIndexList_negativeParticles[iii];
          }
        }

        if(spectraStructAlreadyMade){
          m_spectraStruct[partIndex][centIndex][rapIndexOfAlreadyMade][1]->setRawETOFSpectra(projections_ETOF_Minus[rapIndex]);
        }else{
          int nextNullRapIndex = -1;
          for(int iii = 0; iii < g_RawSpecMod_MaxNumRapBins; iii++){
            if(!m_spectraStruct[partIndex][centIndex][iii][1]){
              nextNullRapIndex = iii; break;
            }
          }
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1] = new SpectraStructure();
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1]->Construct(partIndex, -1, nextNullRapIndex,
             rapVal,centIndex,lowCentEdge, highCentEdge);
          m_spectraStruct[partIndex][centIndex][nextNullRapIndex][1]->setRawETOFSpectra(projections_ETOF_Minus[rapIndex]);
          rapidityList_negativeParticles.push_back(rapVal);
          rapIndexList_negativeParticles.push_back(nextNullRapIndex);
        }
      }
      cout <<" finished making the projections" << endl;
    } //centrality
  }//particle
  return;
}





//########################################################################################################################

void RawSpectraModifier::writeRapiditySpectra(string a_filename, string a_option){
  #ifdef _RAWSPECMOD_DEBUG_
    cout << "RawSpectraModifier::writeRapiditySpectra(" << a_filename << ", "<< a_option << ")" << endl;
  #endif
  TGraph2DErrors* graph2D_RapMtM0_Plus[9][16];
  TGraph2DErrors* graph2D_RapMtM0_Minus[9][16];
  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int centIndex = 0; centIndex < 16; centIndex++){
      graph2D_RapMtM0_Plus[partIndex][centIndex] = NULL;
      graph2D_RapMtM0_Minus[partIndex][centIndex] = NULL;
    }
  }

  //########### MAKE THE GRAPHS  ######################
  for(int partIndex = 0; partIndex < 9; partIndex++){
    if(partIndex == 3) continue;
    double mass = m_partInfo->GetParticleMass(partIndex);
    for(int centIndex = 0; centIndex < 16; centIndex++){
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        int charge = 0;
        if(pmIndex == 0) charge = 1;
        else             charge = -1;
        #ifdef _RAWSPECMOD_DEBUG_
          cout << " Part" << partIndex << "Cent" << centIndex << endl;
          cout << "    making TGraph2DErrors" << endl;
        #endif
        TGraph2DErrors* rawSpectrum2DRapMtM0 = new TGraph2DErrors();
        rawSpectrum2DRapMtM0->SetName(Form("rawSpectra2D_RapMtM0_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpectrum2DRapMtM0->SetTitle(Form("Yield for Centrality Index: %02d  ; y ; m_{T} - m_{0} ; Yield", centIndex));

        if(pmIndex == 0){
          graph2D_RapMtM0_Plus[partIndex][centIndex] = rawSpectrum2DRapMtM0;
        }else{
          graph2D_RapMtM0_Minus[partIndex][centIndex] = rawSpectrum2DRapMtM0;
        }

        //FILL THE Rapidity TGRAPH 2D
        TH2D* specHistos[3] = {NULL,NULL,NULL};
        if(pmIndex == 0){
          specHistos[0] = m_dEdxSpectraPlus[partIndex][centIndex];
          specHistos[1] = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
          if(m_loadedETOF) specHistos[2] = m_InvBetaETOFSpectraPlus[partIndex][centIndex];
        }else{
          specHistos[0] = m_dEdxSpectraMinus[partIndex][centIndex];
          specHistos[1] = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];
          if(m_loadedETOF) specHistos[2] = m_InvBetaETOFSpectraMinus[partIndex][centIndex];
        }

        for(int histoIndex = 0; histoIndex < 3; histoIndex++){
          if(!specHistos[histoIndex]) continue;
          for(int rapBin = 1; rapBin <= specHistos[histoIndex]->GetNbinsX(); rapBin++){
            double rapCenter = specHistos[histoIndex]->GetXaxis()->GetBinCenter(rapBin);
            TGraphErrors* specGraph[6] = {NULL,NULL,NULL, NULL,NULL,NULL};
            specGraph[0] = m_spectraPlusDeDx[partIndex][centIndex][rapBin - 1];
            specGraph[1] = m_spectraMinusDeDx[partIndex][centIndex][rapBin - 1];
            specGraph[2] = m_spectraPlusInvBeta[partIndex][centIndex][rapBin - 1];
            specGraph[3] = m_spectraMinusInvBeta[partIndex][centIndex][rapBin - 1];
            if(m_loadedETOF){
              specGraph[4] = m_spectraPlusInvBetaETOF[partIndex][centIndex][rapBin - 1];
              specGraph[5] = m_spectraMinusInvBetaETOF[partIndex][centIndex][rapBin - 1];
            }

            TGraphErrors* graphToAdd = specGraph[2*histoIndex + pmIndex];

            for (int pointIndex=0; pointIndex < graphToAdd->GetN(); pointIndex++){
              //Get the point info
              double mTm0  = graphToAdd->GetX()[pointIndex];
              double yield = graphToAdd->GetY()[pointIndex];
              double xErr  = graphToAdd->GetEX()[pointIndex];
              double yErr  = graphToAdd->GetEY()[pointIndex];
              int pointNum2D = rawSpectrum2DRapMtM0->GetN();
              rawSpectrum2DRapMtM0->SetPoint(pointNum2D, rapCenter, mTm0, yield);
              rawSpectrum2DRapMtM0->SetPointError(pointNum2D, 0.5*specHistos[histoIndex]->GetXaxis()->GetBinWidth(rapBin) ,xErr,yErr);
            }// point loop
          }//eta loop
        }//histoIndex

        cout << "Done Filling the full and 2d graphs" << endl;

      }// end plus minus loop
    } // end cent
  } // end part

  //write everything out
  TFile* outFile = new TFile(a_filename.c_str(),a_option.c_str());

  if(!outFile->IsOpen()){
    cout << "Error: " << a_filename.c_str() << " could not be opened correctly" << endl;
  }

  for(int partIndex = 0; partIndex < 9; partIndex++){
    //if(partIndex == 3) continue;
    string spectra2DdirPlus = Form("Spectra_%s_2D",m_partInfo->GetParticleName(partIndex,1).Data());
    string spectra2DdirMinus = Form("Spectra_%s_2D",m_partInfo->GetParticleName(partIndex,-1).Data());
    string spectraGraphDirPlus = Form("Spectra_%s",m_partInfo->GetParticleName(partIndex,1).Data());
    string spectraGraphDirMinus = Form("Spectra_%s",m_partInfo->GetParticleName(partIndex,-1).Data());
    #ifdef _RAWSPECMOD_DEBUG_
      cout << "Writing part Index: " << partIndex << endl;
    #endif
    outFile->cd();
    outFile->mkdir(spectra2DdirPlus.c_str());
    outFile->cd(spectra2DdirPlus.c_str());
    for(int centIndex = 0; centIndex < 16; centIndex++){
      HistogramUtilities::ConditionalWrite(graph2D_RapMtM0_Plus[partIndex][centIndex]);
    } // end cent

    outFile->cd();
    outFile->mkdir(spectra2DdirMinus.c_str());
    outFile->cd(spectra2DdirMinus.c_str());
    for(int centIndex = 0; centIndex < 16; centIndex++){
      HistogramUtilities::ConditionalWrite(graph2D_RapMtM0_Minus[partIndex][centIndex]);
    } // end cent
    #ifdef _RAWSPECMOD_DEBUG_
      cout << "Wrote the full histos and 2D Graphs" << endl;
    #endif
    outFile->cd();
    outFile->mkdir(spectraGraphDirPlus.c_str());
    outFile->cd(spectraGraphDirPlus.c_str());
    for(int centIndex = 0; centIndex < 16; centIndex++){
      for(int projIndex = 0; projIndex < g_RawSpecMod_MaxNumRapBins; projIndex++){
        HistogramUtilities::ConditionalWrite(m_rawSpectraPlusDeDx[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_inclusiveSpectraPlusDeDx[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_spectraPlusDeDx[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_rawSpectraPlusInvBeta[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_inclusiveSpectraPlusInvBeta[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_spectraPlusInvBeta[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_rawSpectraPlusInvBetaETOF[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_inclusiveSpectraPlusInvBetaETOF[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_spectraPlusInvBetaETOF[partIndex][centIndex][projIndex]);
      }
    } // end cent
    
    outFile->mkdir(spectraGraphDirMinus.c_str());
    outFile->cd(spectraGraphDirMinus.c_str());
    for(int centIndex = 0; centIndex < 16; centIndex++){
      for(int projIndex = 0; projIndex < g_RawSpecMod_MaxNumRapBins; projIndex++){
        HistogramUtilities::ConditionalWrite(m_rawSpectraMinusDeDx[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_inclusiveSpectraMinusDeDx[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_spectraMinusDeDx[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_rawSpectraMinusInvBeta[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_inclusiveSpectraMinusInvBeta[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_spectraMinusInvBeta[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_rawSpectraMinusInvBetaETOF[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_inclusiveSpectraMinusInvBetaETOF[partIndex][centIndex][projIndex]);
        HistogramUtilities::ConditionalWrite(m_spectraMinusInvBetaETOF[partIndex][centIndex][projIndex]);
      }
    } // end cent
  } // end part

  outFile->Close();
  cout << "Wrote Spectra Successfully" << endl;
  return;

}



void RawSpectraModifier::populateInclusiveSpectra(){

  #ifdef _RAWSPECMOD_DEBUG_
    cout << " Call to RawSpectraModifier::populateInclusiveSpectra()... " << endl;
  #endif

  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int centIndex = 0; centIndex < 16; centIndex++){
      for(int yIndex = 0; yIndex < g_RawSpecMod_MaxNumRapBins; yIndex++){
        if(m_spectraPlusDeDx[partIndex][centIndex][yIndex]){
          m_inclusiveSpectraPlusDeDx[partIndex][centIndex][yIndex] = (TGraphErrors*) m_spectraPlusDeDx[partIndex][centIndex][yIndex]->Clone();
          string inclusiveName = Form("inclusive_dEdxSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,yIndex);
          m_inclusiveSpectraPlusDeDx[partIndex][centIndex][yIndex]->SetName(inclusiveName.c_str());
        }
        if(m_spectraMinusDeDx[partIndex][centIndex][yIndex]){
          m_inclusiveSpectraMinusDeDx[partIndex][centIndex][yIndex] = (TGraphErrors*) m_spectraMinusDeDx[partIndex][centIndex][yIndex]->Clone();
          string inclusiveName = Form("inclusive_dEdxSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,yIndex);
          m_inclusiveSpectraMinusDeDx[partIndex][centIndex][yIndex]->SetName(inclusiveName.c_str());
        }
        if(m_spectraPlusInvBeta[partIndex][centIndex][yIndex]){
          m_inclusiveSpectraPlusInvBeta[partIndex][centIndex][yIndex] = (TGraphErrors*) m_spectraPlusInvBeta[partIndex][centIndex][yIndex]->Clone();
          string inclusiveName = Form("inclusive_InvBetaBTOFSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,yIndex);
          m_inclusiveSpectraPlusInvBeta[partIndex][centIndex][yIndex]->SetName(inclusiveName.c_str());
        }
        if(m_spectraMinusInvBeta[partIndex][centIndex][yIndex]){
          m_inclusiveSpectraMinusInvBeta[partIndex][centIndex][yIndex] = (TGraphErrors*) m_spectraMinusInvBeta[partIndex][centIndex][yIndex]->Clone();
          string inclusiveName = Form("inclusive_InvBetaBTOFSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,yIndex);
          m_inclusiveSpectraMinusInvBeta[partIndex][centIndex][yIndex]->SetName(inclusiveName.c_str());
        }
        if(m_spectraPlusInvBetaETOF[partIndex][centIndex][yIndex]){
          m_inclusiveSpectraPlusInvBetaETOF[partIndex][centIndex][yIndex] = (TGraphErrors*) m_spectraPlusInvBetaETOF[partIndex][centIndex][yIndex]->Clone();
          string inclusiveName = Form("inclusive_InvBetaETOFSpectra_%s_Plus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,yIndex);
          m_inclusiveSpectraPlusInvBetaETOF[partIndex][centIndex][yIndex]->SetName(inclusiveName.c_str());
        }
        if(m_spectraMinusInvBetaETOF[partIndex][centIndex][yIndex]){
          m_inclusiveSpectraMinusInvBetaETOF[partIndex][centIndex][yIndex] = (TGraphErrors*) m_spectraMinusInvBetaETOF[partIndex][centIndex][yIndex]->Clone();
          string inclusiveName = Form("inclusive_InvBetaETOFSpectra_%s_Minus_Cent%02d_yIndex_%02d",m_partInfo->GetParticleName(partIndex,0).Data(),centIndex,yIndex);
          m_inclusiveSpectraMinusInvBetaETOF[partIndex][centIndex][yIndex]->SetName(inclusiveName.c_str());
        }
      }//rapidity
    }//cent
  }//part
}


//########################################################################################################################



/*void RawSpectraModifier::makeRapiditySpectra(string a_filename, string a_option){

  TGraph2DErrors* graph2D_RapMtM0_Plus[m_numParticles][m_numCentralities];
  TGraph2DErrors* graph2D_RapMtM0_Minus[m_numParticles][m_numCentralities];

  TH2D* rawSpec_RapMtM0_FULL_Plus[m_numParticles][m_numCentralities];
  TH2D* rawSpec_RapMtM0_FULL_Minus[m_numParticles][m_numCentralities];  

  //########### MAKE THE GRAPHS  ######################


  //Now Make Spectra from Fits
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    //if(partIndex ==3) continue;
    double mass = m_partInfo->GetParticleMass(partIndex);

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        cout << "Setting up TGraph2Ds part=" << partIndex << "   Cent Index " <<  centIndex << endl;
        int charge = 0;
        if(pmIndex == 0){ 
          charge = 1;
        }else{
          charge = -1;
        }

        cout << " Part" << partIndex << "Cent" << centIndex << endl;
        cout << "    making TGraph2DErrors" << endl;

        TGraph2DErrors* rawSpectrum2DRapMtM0 = new TGraph2DErrors();
        rawSpectrum2DRapMtM0->SetName(Form("rawSpectra2D_RapMtM0_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpectrum2DRapMtM0->SetTitle(Form("Yield for Centrality Index: %02d  ; y ; m_{T} - m_{0} ; Yield", centIndex));

        TH2D* rawSpecRapMtM0FULL = (TH2D*) m_dEdxSpectraPlus[partIndex][centIndex]->Clone();
        rawSpecRapMtM0FULL->SetName(Form("rawSpectra_RapMtM0_%s_Tot_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpecRapMtM0FULL->SetTitle(
                   Form("%s Yield for Centrality Index: %02d  ; y ; m_{T} - m_{0} ; Yield", 
                          m_partInfo->GetParticleName(partIndex,charge).Data(), centIndex));

	HistogramUtilities::zeroBins2D(rawSpecRapMtM0FULL);

        if(pmIndex == 0){
            graph2D_RapMtM0_Plus[partIndex][centIndex] = rawSpectrum2DRapMtM0;
            rawSpec_RapMtM0_FULL_Plus[partIndex][centIndex] = rawSpecRapMtM0FULL;
        }else{
            graph2D_RapMtM0_Minus[partIndex][centIndex] = rawSpectrum2DRapMtM0;
            rawSpec_RapMtM0_FULL_Minus[partIndex][centIndex] = rawSpecRapMtM0FULL;
        }

        //FILL THE Rapidity TGRAPH 2D

        TH2D* specHistos[3] = {NULL,NULL,NULL};
        if(pmIndex == 0){
          specHistos[0] = m_dEdxSpectraPlus[partIndex][centIndex];
          specHistos[1] = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
          if(m_loadedETOF) specHistos[2] = m_InvBetaETOFSpectraPlus[partIndex][centIndex];
        }else{
          specHistos[0] = m_dEdxSpectraMinus[partIndex][centIndex];
          specHistos[1] = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];
          if(m_loadedETOF) specHistos[2] = m_InvBetaETOFSpectraMinus[partIndex][centIndex];
        }

        for(int histoIndex = 0; histoIndex < 3; histoIndex++){
          if(!specHistos[histoIndex]) continue;
          for(int rapBin = 1; rapBin <= specHistos[histoIndex]->GetNbinsX(); rapBin++){
            double rapCenter = specHistos[histoIndex]->GetXaxis()->GetBinCenter(rapBin);
            for(int mTm0Bin = 1; mTm0Bin <= specHistos[histoIndex]->GetNbinsY(); mTm0Bin++){
              double mTm0Center = specHistos[histoIndex]->GetYaxis()->GetBinCenter(mTm0Bin);
              double content = specHistos[histoIndex]->GetBinContent(rapBin,mTm0Bin);
              double contentError = specHistos[histoIndex]->GetBinError(rapBin,mTm0Bin);
              if(content != 0.0 && !std::isnan(content) && !std::isnan(contentError) && content > 0.0){
                int pointNum2D = rawSpectrum2DRapMtM0->GetN();
                rawSpectrum2DRapMtM0->SetPoint(pointNum2D, rapCenter, mTm0Center, content);
                rawSpectrum2DRapMtM0->SetPointError(pointNum2D, 0.0,0.0,contentError);

                int horBin = rapBin;
                int vertBin = mTm0Bin;
                if(histoIndex==1) horBin+=m_tofIndexOffset[partIndex];
                if(histoIndex==2) horBin+=m_etofIndexOffset[partIndex];

                rawSpecRapMtM0FULL->SetBinContent(horBin,vertBin,content);
                rawSpecRapMtM0FULL->SetBinError(horBin,vertBin,contentError);
              }
            }//ptLoop
          }//eta loop
        }//histoIndex

        cout << "Done Filling the full and 2d graphs" << endl;

      }// end plus minus loop

    } // end cent

  } // end part

  //write everything out

  TFile* outFile = new TFile(a_filename.c_str(),a_option.c_str());

  if(!outFile->IsOpen()){
    cout << "Error: " << a_filename.c_str() << " could not be opened correctly" << endl;
  }
  //vector< string > initDir;


  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    //if(partIndex == 3) continue;
    cout << "Writing part Index: " << partIndex << endl;
    outFile->cd();
    outFile->mkdir(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,1).Data()));
    outFile->cd(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,1).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      graph2D_RapMtM0_Plus[partIndex][centIndex]->Write();
      rawSpec_RapMtM0_FULL_Plus[partIndex][centIndex]->Write();
    } // end cent

    outFile->cd();
    outFile->mkdir(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,-1).Data()));
    outFile->cd(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,-1).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      graph2D_RapMtM0_Minus[partIndex][centIndex]->Write();
      rawSpec_RapMtM0_FULL_Minus[partIndex][centIndex]->Write();      
    } // end cent

    cout << "Wrote the full histos and 2D Graphs" << endl;

    outFile->cd();
    outFile->mkdir(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,1).Data()));
    outFile->cd(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,1).Data()));

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

      vector<TGraphErrors*> projections_dEdx = HistogramUtilities::VerticalProjectionGraphs(
                m_dEdxSpectraPlus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_BTOF = HistogramUtilities::VerticalProjectionGraphs(
                m_InvBetaBTOFSpectraPlus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_ETOF;
      if(m_loadedETOF){
        projections_ETOF = HistogramUtilities::VerticalProjectionGraphs(
                m_InvBetaETOFSpectraPlus[partIndex][centIndex], "y", "y", true);
      }
      vector<TGraphErrors*> projections_total = HistogramUtilities::VerticalProjectionGraphs(
                rawSpec_RapMtM0_FULL_Plus[partIndex][centIndex], "y", "y", true);

     cout <<" finished making the projections" << endl;

      for(int projIndex = 0; projIndex < (int) projections_dEdx.size(); projIndex++)  projections_dEdx[projIndex]->Write();
      for(int projIndex = 0; projIndex < (int) projections_BTOF.size(); projIndex++)  projections_BTOF[projIndex]->Write();
      if(m_loadedETOF){
        for(int projIndex = 0; projIndex < (int) projections_ETOF.size(); projIndex++)  projections_ETOF[projIndex]->Write();
      }
      for(int projIndex = 0; projIndex < (int) projections_total.size(); projIndex++)  projections_total[projIndex]->Write();

    } // end cent

    outFile->mkdir(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,-1).Data()));
    outFile->cd(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,-1).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

      vector<TGraphErrors*> projections_dEdx = HistogramUtilities::VerticalProjectionGraphs(
                m_dEdxSpectraMinus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_BTOF = HistogramUtilities::VerticalProjectionGraphs(
                m_InvBetaBTOFSpectraMinus[partIndex][centIndex], "y", "y", true);

      vector<TGraphErrors*> projections_ETOF;
      if(m_loadedETOF){
         projections_ETOF= HistogramUtilities::VerticalProjectionGraphs(
                m_InvBetaETOFSpectraMinus[partIndex][centIndex], "y", "y", true);
      }
      vector<TGraphErrors*> projections_total = HistogramUtilities::VerticalProjectionGraphs(
                rawSpec_RapMtM0_FULL_Plus[partIndex][centIndex], "y", "y", true);
      for(int projIndex = 0; projIndex < (int) projections_dEdx.size(); projIndex++)  projections_dEdx[projIndex]->Write();
      for(int projIndex = 0; projIndex < (int) projections_BTOF.size(); projIndex++)  projections_BTOF[projIndex]->Write();
      if(m_loadedETOF){
        for(int projIndex = 0; projIndex < (int) projections_ETOF.size(); projIndex++)  projections_ETOF[projIndex]->Write();
      }
      for(int projIndex = 0; projIndex < (int) projections_total.size(); projIndex++)  projections_total[projIndex]->Write();
      
    } // end cent

  } // end part

  outFile->Close();
  return;

}


*/

//########################################################################################################################
//########################################################################################################################










void RawSpectraModifier::makeRapiditySpectraFromEtaSpectra(string a_filename, string a_option){


  TGraph2DErrors* graph2D_EtaPt_Plus[m_numParticles][m_numCentralities];
  TGraph2DErrors* graph2D_EtaPt_Minus[m_numParticles][m_numCentralities];

  TGraph2DErrors* graph2D_RapMtM0_Plus[m_numParticles][m_numCentralities];
  TGraph2DErrors* graph2D_RapMtM0_Minus[m_numParticles][m_numCentralities];

  //these store the errors of the above rapidity mtm0 spectra so I can inerpolate the errors
  TGraph2D* graph2D_Error_RapMtM0_Plus[m_numParticles][m_numCentralities];
  TGraph2D* graph2D_Error_RapMtM0_Minus[m_numParticles][m_numCentralities];

  TH2D* rawSpec_RapMtM0_Interpolated_Plus[m_numParticles][m_numCentralities];
  TH2D* rawSpec_RapMtM0_Interpolated_Minus[m_numParticles][m_numCentralities];

  TGraph2DErrors* graph2D_RapMtM0_Interpolated_Plus[m_numParticles][m_numCentralities];
  TGraph2DErrors* graph2D_RapMtM0_Interpolated_Minus[m_numParticles][m_numCentralities];


  TH2D* rapHisto_tpc_plus[m_numParticles][m_numCentralities];
  TH2D* rapHisto_tpc_minus[m_numParticles][m_numCentralities];
  TH2D* rapHisto_btof_plus[m_numParticles][m_numCentralities];
  TH2D* rapHisto_btof_minus[m_numParticles][m_numCentralities];
  TH2D* rapHisto_etof_plus[m_numParticles][m_numCentralities];
  TH2D* rapHisto_etof_minus[m_numParticles][m_numCentralities];  




  //########### MAKE THE GRAPHS  ######################


  //Now Make Spectra from Fits
  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    //if(partIndex ==3) continue;
    double mass = m_partInfo->GetParticleMass(partIndex);

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      for(int pmIndex = 0; pmIndex < 2; pmIndex++){
        cout << "Setting up TGraph2Ds part=" << partIndex << "   Cent Index " <<  centIndex << endl;
        int charge = 0;
        if(pmIndex == 0){ 
          charge = 1;
        }else{
          charge = -1;
        }

        TH2D* rawSpecInterp = HistogramUtilities::make2DHistLogBins("binStructure", "binStructure; y ; m_{T} - m_{0}; Yield",
          false, 31, -2.05, 2.05,
          true, m_rap_mTm0_dEdx_allCent[partIndex]->GetNbinsY(), 
          PhysMath::mTm0(mass,m_rap_mTm0_dEdx_allCent[partIndex]->GetYaxis()->GetBinLowEdge(1)),
          PhysMath::mTm0(mass,m_rap_mTm0_dEdx_allCent[partIndex]->GetYaxis()->GetBinLowEdge(m_rap_mTm0_dEdx_allCent[partIndex]->GetNbinsY()+1))
          );

        TGraph2DErrors* rawSpectrum2DEtaPt = new TGraph2DErrors();
        rawSpectrum2DEtaPt->SetName("graph2D1");
        TGraph2DErrors* rawSpectrum2DRapMtM0 = new TGraph2DErrors();
        rawSpectrum2DRapMtM0->SetName("graph2D2");
        TGraph2D* rawSpectrum2DRapMtM0_Error = new TGraph2D();
        rawSpectrum2DRapMtM0_Error->SetName("rapErrors");
        TGraph2DErrors* rawSpectrum2DRapMtM0Interp = new TGraph2DErrors();
        rawSpectrum2DRapMtM0Interp->SetName("graph2D4");


        rawSpectrum2DEtaPt->SetName(Form("rawSpectra2D_EtaPt_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpectrum2DEtaPt->SetTitle(Form("Yield for Centrality Index: %02d ; #eta ; p_{T} ; Yield", centIndex));
        
        rawSpectrum2DRapMtM0->SetName(Form("rawSpectra2D_RapMtM0_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpectrum2DRapMtM0->SetTitle(Form("Yield for Centrality Index: %02d  ; y ; m_{T} - m_{0} ; Yield", centIndex));

        rawSpectrum2DRapMtM0_Error->SetName(Form("rawSpectra2D_RapMtM0_Error_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpectrum2DRapMtM0_Error->SetTitle(Form("Yield Error for Centrality Index: %02d  ; y ; m_{T} - m_{0} ; Yield", centIndex));
        
        rawSpectrum2DRapMtM0Interp->SetName(Form("rawSpectra2D_RapMtM0Interp_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpectrum2DRapMtM0Interp->SetTitle(Form("Yield for Centrality Index: %02d ; y ; m_{T} - m_{0} ; Yield", centIndex));
        
        rawSpecInterp->SetName(Form("rawSpectraHisto_RapMtM0Interp_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        rawSpecInterp->SetTitle(Form("Yield for Centrality Index: %02d ; y ; m_{T} - m_{0} ; Yield", centIndex));

        if(pmIndex == 0){
          graph2D_EtaPt_Plus[partIndex][centIndex] = rawSpectrum2DEtaPt;
          graph2D_RapMtM0_Plus[partIndex][centIndex] = rawSpectrum2DRapMtM0;
          graph2D_Error_RapMtM0_Plus[partIndex][centIndex] = rawSpectrum2DRapMtM0_Error;          
          rawSpec_RapMtM0_Interpolated_Plus[partIndex][centIndex] = rawSpecInterp;          
          graph2D_RapMtM0_Interpolated_Plus[partIndex][centIndex] = rawSpectrum2DRapMtM0Interp;      
        }else{
          graph2D_EtaPt_Minus[partIndex][centIndex] = rawSpectrum2DEtaPt;
          graph2D_RapMtM0_Minus[partIndex][centIndex] = rawSpectrum2DRapMtM0;
          graph2D_Error_RapMtM0_Minus[partIndex][centIndex] = rawSpectrum2DRapMtM0_Error;
          rawSpec_RapMtM0_Interpolated_Minus[partIndex][centIndex] = rawSpecInterp;
          graph2D_RapMtM0_Interpolated_Minus[partIndex][centIndex] = rawSpectrum2DRapMtM0Interp;
        }



        //FILL THE ETA PT TGRAPH 2D

        TH2D* specHistos[3] = {NULL,NULL,NULL};
        if(pmIndex == 0){
          specHistos[0] = m_dEdxSpectraPlus[partIndex][centIndex];
          specHistos[1] = m_InvBetaBTOFSpectraPlus[partIndex][centIndex];
          specHistos[2] = m_InvBetaETOFSpectraPlus[partIndex][centIndex];
        }else{
          specHistos[0] = m_dEdxSpectraMinus[partIndex][centIndex];
          specHistos[1] = m_InvBetaBTOFSpectraMinus[partIndex][centIndex];
          specHistos[2] = m_InvBetaETOFSpectraMinus[partIndex][centIndex];
        }

        for(int histoIndex = 0; histoIndex < 3; histoIndex++){
          for(int etaBin = 1; etaBin <= specHistos[histoIndex]->GetNbinsX(); etaBin++){
            double etaCenter = specHistos[histoIndex]->GetXaxis()->GetBinCenter(etaBin);
            for(int pTBin = 1; pTBin <= specHistos[histoIndex]->GetNbinsY(); pTBin++){
              double pTCenter = specHistos[histoIndex]->GetYaxis()->GetBinCenter(pTBin);
              double rapidity = PhysMath::rapFromEtaPt(mass, etaCenter, pTCenter);
              double mTm0 = PhysMath::mTm0(mass, pTCenter);
              double content = specHistos[histoIndex]->GetBinContent(etaBin,pTBin);
              double contentError = specHistos[histoIndex]->GetBinError(etaBin,pTBin);
              if(content != 0.0 && !std::isnan(content) && !std::isnan(contentError)){
                int pointNum2DEtaPt = rawSpectrum2DEtaPt->GetN();
                rawSpectrum2DEtaPt->SetPoint(pointNum2DEtaPt, etaCenter, pTCenter, content);
                rawSpectrum2DEtaPt->SetPointError(pointNum2DEtaPt, 0.0,0.0,contentError);

                int pointNum2D = rawSpectrum2DRapMtM0->GetN();                
                rawSpectrum2DRapMtM0->SetPoint(pointNum2D, rapidity, mTm0, PhysMath::eta_to_rap_factor(mass, etaCenter,pTCenter)*content);
                rawSpectrum2DRapMtM0->SetPointError(pointNum2D,0.0,0.0,PhysMath::eta_to_rap_factor(mass, etaCenter,pTCenter)*contentError);
                rawSpectrum2DRapMtM0_Error->SetPoint(pointNum2D, rapidity, mTm0, PhysMath::eta_to_rap_factor(mass, etaCenter,pTCenter)*contentError);
              }
            }//ptLoop
          }//eta loop
        }//histoIndex



        // now, use delauney interplation to get spectra in terms of just rapidity
        for(int binX = 1; binX <= rawSpecInterp->GetNbinsX(); binX++){
          double x_center = rawSpecInterp->GetXaxis()->GetBinCenter(binX);
          for(int binY = 1; binY <= rawSpecInterp->GetNbinsY(); binY++){
            double y_center = rawSpecInterp->GetYaxis()->GetBinCenter(binY);
            double content = rawSpectrum2DRapMtM0->Interpolate(x_center,y_center);
            double error = rawSpectrum2DRapMtM0_Error->Interpolate(x_center,y_center);
            if(!std::isnan(content) && !std::isnan(error) && fabs(error/content) < 1.0){  
              rawSpecInterp->SetBinContent(binX,binY,content);
              rawSpecInterp->SetBinError(binX,binY,error);
            }
          }
        }

        // make a histogram from the 
        HistogramUtilities::makeGraph2D(rawSpecInterp,rawSpectrum2DRapMtM0Interp);



        //now extimate which points came from which detectors, due to interpolation
        

        TH2D* tpcHisto = (TH2D*) rawSpecInterp->Clone();
        tpcHisto->SetName(Form("rawSpectra_RapMtM0_TPC_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        TH2D* btofHisto = (TH2D*) rawSpecInterp->Clone();
        btofHisto->SetName(Form("rawSpectra_RapMtM0_BTOF_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));
        TH2D* etofHisto = (TH2D*) rawSpecInterp->Clone();
        etofHisto->SetName(Form("rawSpectraHisto_RapMtM0_ETOF_%s_Cent%02d",
          m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex));

        //#### TPC SPECTRA #######  zero the yield if its not in the momentum range previously specified
        for(int rapBin = 1; rapBin <= tpcHisto->GetNbinsX(); rapBin++){
          double rapCenter = tpcHisto->GetXaxis()->GetBinCenter(rapBin);
          for(int mTm0Bin = 1; mTm0Bin <= tpcHisto->GetNbinsY(); mTm0Bin++){
            double mTm0Center = tpcHisto->GetYaxis()->GetBinCenter(mTm0Bin);
            double momentum = PhysMath::pTot_from_y_mTm0(rapCenter, mTm0Center, mass);
            if(momentum > m_momStop_tpc[partIndex] || momentum < m_momStart_tpc[partIndex]){
             tpcHisto->SetBinContent(rapBin,mTm0Bin,0.0);
             tpcHisto->SetBinError(rapBin,mTm0Bin,0.0);
           }
          } //pt Bin
        }// eta Bin

        //#### BTOF SPECTRA #######
        for(int rapBin = 1; rapBin <= btofHisto->GetNbinsX(); rapBin++){
          double rapCenter = btofHisto->GetXaxis()->GetBinCenter(rapBin);
          for(int mTm0Bin = 1; mTm0Bin <= btofHisto->GetNbinsY(); mTm0Bin++){
            double mTm0Center = btofHisto->GetYaxis()->GetBinCenter(mTm0Bin);
            double momentum = PhysMath::pTot_from_y_mTm0(rapCenter, mTm0Center, mass);
            if(momentum > m_momStop_btof[partIndex] || momentum < m_momStart_btof[partIndex]){
             btofHisto->SetBinContent(rapBin,mTm0Bin,0.0);
             btofHisto->SetBinError(rapBin,mTm0Bin,0.0);
           }
          } //pt Bin
        }// eta Bin

        //#### ETOF SPECTRA #######
        for(int rapBin = 1; rapBin <= etofHisto->GetNbinsX(); rapBin++){
          double rapCenter = etofHisto->GetXaxis()->GetBinCenter(rapBin);
          for(int mTm0Bin = 1; mTm0Bin <= etofHisto->GetNbinsY(); mTm0Bin++){
            double mTm0Center = etofHisto->GetYaxis()->GetBinCenter(mTm0Bin);
            double momentum = PhysMath::pTot_from_y_mTm0(rapCenter, mTm0Center, mass);
            if(momentum > m_momStop_etof[partIndex] || momentum < m_momStart_etof[partIndex]){
             etofHisto->SetBinContent(rapBin,mTm0Bin,0.0);
             etofHisto->SetBinError(rapBin,mTm0Bin,0.0);
           }

          } //pt Bin
        }// eta Bin


        if(pmIndex == 0){
          rapHisto_tpc_plus[partIndex][centIndex] = tpcHisto;
          rapHisto_btof_plus[partIndex][centIndex] = btofHisto;
          rapHisto_etof_plus[partIndex][centIndex] = etofHisto;

        }else{
          rapHisto_tpc_minus[partIndex][centIndex] = tpcHisto;
          rapHisto_btof_minus[partIndex][centIndex] = btofHisto;
          rapHisto_etof_minus[partIndex][centIndex] = etofHisto;
        }

      }// end plus minus loop

    } // end cent

  } // end part

  //write everything out

  TFile* outFile = new TFile(a_filename.c_str(),a_option.c_str());

  //vector< string > initDir;


  for(int partIndex = 0; partIndex < m_numParticles; partIndex++){
    //if(partIndex == 3) continue;
    cout << "Writing part Index: " << partIndex << endl;
    outFile->mkdir(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,1).Data()));
    outFile->cd(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,1).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      graph2D_EtaPt_Plus[partIndex][centIndex]->Write();
      graph2D_RapMtM0_Plus[partIndex][centIndex]->Write();
      graph2D_Error_RapMtM0_Plus[partIndex][centIndex]->Write();
      rawSpec_RapMtM0_Interpolated_Plus[partIndex][centIndex]->Write();
      graph2D_RapMtM0_Interpolated_Plus[partIndex][centIndex]->Write();
    } // end cent

    outFile->mkdir(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,-1).Data()));
    outFile->cd(Form("RawSpectra_%s_2D",m_partInfo->GetParticleName(partIndex,-1).Data()));
    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){
      graph2D_EtaPt_Minus[partIndex][centIndex]->Write();
      graph2D_RapMtM0_Minus[partIndex][centIndex]->Write();
      graph2D_Error_RapMtM0_Minus[partIndex][centIndex]->Write();
      rawSpec_RapMtM0_Interpolated_Minus[partIndex][centIndex]->Write();
      graph2D_RapMtM0_Interpolated_Minus[partIndex][centIndex]->Write();
      
    } // end cent



    outFile->mkdir(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,1).Data()));
    outFile->cd(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,1).Data()));

    for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

      vector<TGraphErrors*> projections_dEdx = HistogramUtilities::VerticalProjectionGraphs(
        rapHisto_tpc_plus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_BTOF = HistogramUtilities::VerticalProjectionGraphs(
        rapHisto_btof_plus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_ETOF = HistogramUtilities::VerticalProjectionGraphs(
        rapHisto_etof_plus[partIndex][centIndex], "y", "y", true);
      vector<TGraphErrors*> projections_total = HistogramUtilities::VerticalProjectionGraphs(
        rawSpec_RapMtM0_Interpolated_Plus[partIndex][centIndex], "y", "y", true);
      for(int projIndex = 0; projIndex < (int) projections_dEdx.size(); projIndex++)  projections_dEdx[projIndex]->Write();
        for(int projIndex = 0; projIndex < (int) projections_BTOF.size(); projIndex++)  projections_BTOF[projIndex]->Write();
          for(int projIndex = 0; projIndex < (int) projections_ETOF.size(); projIndex++)  projections_ETOF[projIndex]->Write();
            for(int projIndex = 0; projIndex < (int) projections_total.size(); projIndex++)  projections_total[projIndex]->Write();

    } // end cent

  outFile->mkdir(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,-1).Data()));
  outFile->cd(Form("RawSpectra_%s",m_partInfo->GetParticleName(partIndex,-1).Data()));
  for(int centIndex = 0; centIndex < m_numCentralities; centIndex++){

    vector<TGraphErrors*> projections_dEdx = HistogramUtilities::VerticalProjectionGraphs(
      rapHisto_tpc_minus[partIndex][centIndex], "y", "y", true);
    vector<TGraphErrors*> projections_BTOF = HistogramUtilities::VerticalProjectionGraphs(
      rapHisto_btof_minus[partIndex][centIndex], "y", "y", true);
    vector<TGraphErrors*> projections_ETOF = HistogramUtilities::VerticalProjectionGraphs(
      rapHisto_etof_minus[partIndex][centIndex], "y", "y", true);
    vector<TGraphErrors*> projections_total = HistogramUtilities::VerticalProjectionGraphs(
      rawSpec_RapMtM0_Interpolated_Minus[partIndex][centIndex], "y", "y", true);
    for(int projIndex = 0; projIndex < (int) projections_dEdx.size(); projIndex++)  projections_dEdx[projIndex]->Write();
      for(int projIndex = 0; projIndex < (int) projections_BTOF.size(); projIndex++)  projections_BTOF[projIndex]->Write();
        for(int projIndex = 0; projIndex < (int) projections_ETOF.size(); projIndex++)  projections_ETOF[projIndex]->Write();
          for(int projIndex = 0; projIndex < (int) projections_total.size(); projIndex++)  projections_total[projIndex]->Write();

    } // end cent

  } // end part

  outFile->Close();
  return;

}














SpectraStructure* RawSpectraModifier::findAssociatedSpectraStructure(int a_partIndex, int a_centIndex, int a_charge, double a_rapidity, double a_rapidityTolerence){
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;
  SpectraStructure* returnPointer = NULL;
  for(int index = 0; index < 100; index++){
   if(m_spectraStruct[a_partIndex][a_centIndex][index][pmIndex]
       && fabs(m_spectraStruct[a_partIndex][a_centIndex][index][pmIndex]->getRapidityValue() - a_rapidity) < a_rapidityTolerence){
     return m_spectraStruct[a_partIndex][a_centIndex][index][pmIndex];
   }
  }
  return returnPointer;
}












