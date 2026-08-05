#ifndef _EMBEDDING_BINNER_CXX_
#define _EMBEDDING_BINNER_CXX_

#include "../headers/EmbeddingBinner.h"

EmbeddingBinner::EmbeddingBinner(string a_particleName, CutClass* a_cuts){

  m_partInfo = new ParticleInfo();

  m_inFile = nullptr;
  m_tree   = nullptr;
  m_event  = nullptr;
  m_track  = nullptr;

  TString partStr = Form("%s",a_particleName.c_str());
  if( partStr.EqualTo("pip",TString::kIgnoreCase) ){
    m_partIndex = 0;
    m_particleCharge = 1;
    m_pidSymbol = "#pi^{+}";
    m_mass = 0.13957;
  }else if( partStr.EqualTo("pim",TString::kIgnoreCase) ){
    m_partIndex = 0;
    m_particleCharge = -1;
    m_pidSymbol = "#pi^{-}";
    m_mass = 0.13957;
  }else if( partStr.EqualTo("Kp",TString::kIgnoreCase) ){
    m_partIndex = 1;
    m_particleCharge = 1;
    m_pidSymbol = "K^{+}";
    m_mass = 0.49367;
  }else if( partStr.EqualTo("Km",TString::kIgnoreCase) ){
    m_partIndex = 1;
    m_particleCharge = -1;
    m_pidSymbol = "K^{-}";
    m_mass = 0.49367;
  }else if( partStr.EqualTo("pro",TString::kIgnoreCase) ){
    m_partIndex = 2;
    m_particleCharge = 1;
    m_pidSymbol = "p";
    m_mass = 0.93827;
  }else if( partStr.EqualTo("pbar",TString::kIgnoreCase) ){
    m_partIndex = 2;
    m_particleCharge = -1;
    m_pidSymbol = "#bar{p}";
    m_mass = 0.93827;
  }else if( partStr.EqualTo("deu",TString::kIgnoreCase) ){
    m_partIndex = 4;
    m_particleCharge = 2;
    m_pidSymbol = "d";
    m_mass = 1.87561;
  }else if( partStr.EqualTo("tri",TString::kIgnoreCase) ){
    m_partIndex = 5;
    m_particleCharge = 3;
    m_pidSymbol = "t";
    m_mass = 2.80925;
  }else if( partStr.EqualTo("hel",TString::kIgnoreCase) ){
    m_partIndex = 6;
    m_particleCharge = 2;
    m_pidSymbol = "{}^{3}He";
    m_mass = 2.8094;
  }else if( partStr.EqualTo("alp",TString::kIgnoreCase) ){
    m_partIndex = 7;
    m_particleCharge = 2;
    m_pidSymbol = "#alpha";
    m_mass = 3.72742;
  }

  m_cuts = a_cuts;
  m_nCentBins = m_cuts->getNCentBins();

  if(m_cuts->getEventConfig() == "ColliderCenter"){
    m_nEtaBins = 40;
    m_etaMin = -2.0;
    m_etaMax = 2.0;
  }else if(m_cuts->getEventConfig() == "FixedTarget"){
    m_nEtaBins = 26;
    m_etaMin = -2.5;
    m_etaMax = 0.1;
  }

  m_mTm0DivNum = 1;
  m_useCutClassMtM0Binning = false;


  for(int centIndex = 0; centIndex < 16; centIndex++){
    m_unweighted_eff_y_mTm0[centIndex] = nullptr;
    m_weighted_eff_y_mTm0[centIndex] = nullptr;
    m_eff_y_mTm0_noDCA[centIndex] = nullptr;
    m_dca_eff_y_mTm0[centIndex] = nullptr;
    m_eff_y_mTm0_DCA[centIndex] = nullptr;

  }
  m_emb_weights_y_mTm0_DCA = nullptr;
  m_data_weights_y_mTm0_DCA = nullptr; // this is read-in

}

EmbeddingBinner::~EmbeddingBinner(){
  HistogramUtilities::ConditionalDelete(m_centEvents); 
  HistogramUtilities::ConditionalDelete(m_xyHisto);
  HistogramUtilities::ConditionalDelete(m_zHisto);
  HistogramUtilities::ConditionalDelete(m_tofMatchVsRefMult); 
  HistogramUtilities::ConditionalDelete(m_tofTrayMultVsRefMult); 
  HistogramUtilities::ConditionalDelete(m_numPrimaryDataTracks);
  HistogramUtilities::ConditionalDelete(m_RefMult);
  HistogramUtilities::ConditionalDelete(m_RefMultTrue);
  HistogramUtilities::ConditionalDelete(m_embTrackHisto_byRefMult); 
  HistogramUtilities::ConditionalDelete(m_matchTrackHisto_byRefMult); 

  //CENTRALITY DEPENDANT HISTOGRAMS
  for(int centIndex = 0; centIndex < 16; centIndex++){
    HistogramUtilities::ConditionalDelete(m_embTrackHisto[centIndex]); 
    HistogramUtilities::ConditionalDelete(m_matchTrackHisto[centIndex]);
    HistogramUtilities::ConditionalDelete(m_embTrackHisto_pTEta[centIndex]);
    HistogramUtilities::ConditionalDelete(m_matchTrackHisto_pTEta[centIndex]); 
    HistogramUtilities::ConditionalDelete(m_pTLossHisto[centIndex]); 
    HistogramUtilities::ConditionalDelete(m_mTm0LossHisto[centIndex]);
    HistogramUtilities::ConditionalDelete(m_rapidityLossHisto[centIndex]); 
    HistogramUtilities::ConditionalDelete(m_btof_tracks_reco[centIndex]); 
    HistogramUtilities::ConditionalDelete(m_tpc_tracks_forBTOF_reco[centIndex]); 
    HistogramUtilities::ConditionalDelete(m_btof_tracks_emb[centIndex]); 
    HistogramUtilities::ConditionalDelete(m_tpc_tracks_forBTOF_emb[centIndex]); 

  }

  HistogramUtilities::ConditionalDelete(m_mTm0ReponseMatrixHisto);

  // the reconstructed tracks that were embedded
  HistogramUtilities::ConditionalDelete(m_track_NHits_reco);      
  HistogramUtilities::ConditionalDelete(m_track_NHitsFit_reco);   
  HistogramUtilities::ConditionalDelete(m_track_NHitsPoss_reco);  
  HistogramUtilities::ConditionalDelete(m_track_NHitsDedx_reco);  
  HistogramUtilities::ConditionalDelete(m_track_NHitsFitOverPoss_reco);   
  HistogramUtilities::ConditionalDelete(m_track_pT_reco);     
  HistogramUtilities::ConditionalDelete(m_track_mTm0_reco);   
  HistogramUtilities::ConditionalDelete(m_track_eta_reco);    
  HistogramUtilities::ConditionalDelete(m_track_rapidity_reco);   
  HistogramUtilities::ConditionalDelete(m_track_phi_reco);    
  HistogramUtilities::ConditionalDelete(m_track_pTot_reco);   
  HistogramUtilities::ConditionalDelete(m_track_pZ_reco);     
  HistogramUtilities::ConditionalDelete(m_track_gDCA_reco);   
  
  HistogramUtilities::ConditionalDelete(m_track_pT_emb); 
  HistogramUtilities::ConditionalDelete(m_track_mTm0_emb);   
  HistogramUtilities::ConditionalDelete(m_track_eta_emb); 
  HistogramUtilities::ConditionalDelete(m_track_rapidity_emb);   
  HistogramUtilities::ConditionalDelete(m_track_phi_emb); 
  HistogramUtilities::ConditionalDelete(m_track_pTot_emb);   
  HistogramUtilities::ConditionalDelete(m_track_pZ_emb); 

}



//_________________________________________________________________________________________________________________________________________________________________________________________________________________

void EmbeddingBinner::loadEmbeddingOutput(string a_filename){
  m_inFile = new TFile(a_filename.c_str(),"READ");
  if(!m_inFile || !m_inFile->IsOpen()){
    cout << "ERROR: EmbeddingBinner::loadEmbeddingOutput can't load " << a_filename << endl;
    return;
  }
  m_tree = (TTree*) m_inFile->Get("eventTree");
  if(!m_tree){
    cout << "ERROR: Couldn't load the tree! " << m_tree << endl;
  }
  m_event = new MattMcEvent();
  m_tree->SetBranchAddress("m_event",&m_event);
  m_track = NULL;

}

//_________________________________________________________________________________________________________________________________________________________________________________________________________________


void EmbeddingBinner::loadDataWeightHistograms(string a_filename, int a_mTm0Rebin){
  if(m_partIndex >= 3){
    cout << "ERROR: Particle Index > 2 not supported yet in EmbeddingBinner::loadDataWeightHistogrmas" << endl;
  }
  m_inFile = new TFile(a_filename.c_str(),"READ");
  if(!m_inFile || !m_inFile->IsOpen()){
    cout << "ERROR: EmbeddingBinner::loadDataWeightHistogrmas can't load " << a_filename << endl;
    return;
  }
  m_data_weights_y_mTm0_DCA = (TH3I*) m_inFile->Get(Form("y_mTm0_DCA_%s",m_partInfo->GetParticleName(m_partIndex,m_particleCharge).Data()));
  if(!m_data_weights_y_mTm0_DCA){
    cout << "ERROR: Couldn't load the DCA histogram! " << m_data_weights_y_mTm0_DCA << endl;
    return;
  }

  if(a_mTm0Rebin != 1){
    TH3I* storage_histo = m_data_weights_y_mTm0_DCA;
    m_data_weights_y_mTm0_DCA = HistogramUtilities::reBinVar3D(m_data_weights_y_mTm0_DCA,1,a_mTm0Rebin,1);
    delete storage_histo;
  }

  vector<double> rapidity_edges = HistogramUtilities::getBinEdges_X(m_data_weights_y_mTm0_DCA);
  vector<double> mTm0_edges     = HistogramUtilities::getBinEdges_Y(m_data_weights_y_mTm0_DCA);
  m_dca_eff_y_mTm0_data = new TH2D("dca_eff_y_mTm0_data","DCA Cut Efficiency from Data; y; m_{T}-m_{0} [GeV/c^{2}]",((int)rapidity_edges.size())-1,&rapidity_edges[0],((int)mTm0_edges.size())-1,&mTm0_edges[0]);

  double low_fit_range = 1.5;
  double high_fit_range = 2.9;

  // loop throught the dca histograms and fit them with exponentials
  for(int binX = 1; binX <= m_data_weights_y_mTm0_DCA->GetNbinsX(); binX++){
    for(int binY = 1; binY <= m_data_weights_y_mTm0_DCA->GetNbinsY(); binY++){
      TH1D* proj_histo = m_data_weights_y_mTm0_DCA->ProjectionZ(Form("proj_z%d%d",binX,binY), binX, binX,binY,binY);
      double sum_of_content_in_range = HistogramUtilities::sumOfContent(proj_histo,-1,false,low_fit_range,high_fit_range);
      if(sum_of_content_in_range < 4000){
        delete proj_histo;
        continue;
      }
      HistogramUtilities::makeDifferential(proj_histo);
      TF1* function = PhysMath::getFitSafeExponentialFunction(Form("funct_%d%d",binX,binY), 1.0, 1.0);
      double inv_slope_seed = (low_fit_range - high_fit_range)/log(proj_histo->GetBinContent(proj_histo->FindBin(high_fit_range))/proj_histo->GetBinContent(proj_histo->FindBin(low_fit_range)));
      if(!std::isfinite(inv_slope_seed)){
        inv_slope_seed = 1.0;
      }
      function->SetParameter(1,inv_slope_seed);
      double scale_seed = proj_histo->GetBinContent(proj_histo->FindBin(low_fit_range)) / function->Eval(low_fit_range);
      if(!std::isfinite(scale_seed) || scale_seed <= 0.0) scale_seed = 1.0;
      function->SetParameter(0,scale_seed);
      function->SetParLimits(0,0.01*scale_seed,100.0*scale_seed);
      function->SetParLimits(1,0.001*inv_slope_seed,100.0*inv_slope_seed);
      TFitResultPtr result_ptr = proj_histo->Fit(function, "S","",low_fit_range,high_fit_range);
      TFitResult* result = result_ptr.Get();
      if(!result){
        delete proj_histo;
        continue;
      }

      TF1* funct_result = proj_histo->GetFunction(function->GetName());
      TMatrixDSym* cov_matrix = (TMatrixDSym*) result->GetCovarianceMatrix().Clone();
      //double highest_edge = proj_histo->GetBinLowEdge(HistogramUtilities::lastBinWithContent(proj_histo)+1);
      double highest_edge = proj_histo->GetBinLowEdge(proj_histo->FindBin(high_fit_range) + 1);
      double extra_integral = function->Integral(highest_edge, 50.0);///(50.0 - highest_edge);
      double extra_integral_error = function->IntegralError(highest_edge, 50.0,nullptr,cov_matrix->GetMatrixArray());
      double extra_ammount = extra_integral/(50.0 - highest_edge);
      double extra_ammount_error = extra_integral_error/(50.0 - highest_edge);

      double integral_below_cut_error = 0.0;
      double integral_above_cut_error = 0.0;      
      double integral_below_cut = HistogramUtilities::integral(proj_histo,0,m_cuts->getDCACut(),integral_below_cut_error);
      double integral_above_cut = HistogramUtilities::integral(proj_histo,m_cuts->getDCACut(),highest_edge,integral_above_cut_error);

      double efficiency = integral_below_cut/(integral_below_cut + integral_above_cut + extra_integral);

      double denom_error = sqrt(pow(integral_below_cut_error,2) + pow(integral_above_cut_error,2) + pow(extra_integral_error,2));
      double efficiency_error = efficiency * sqrt(pow(integral_below_cut_error/integral_below_cut,2) + pow(denom_error/(integral_below_cut + integral_above_cut + extra_integral),2));
      if(funct_result->GetParameter(1) > 0.0 && std::isfinite(efficiency) && std::isfinite(efficiency_error)){
        HistogramUtilities::setBinValues(m_dca_eff_y_mTm0_data,binX,binY,efficiency,efficiency_error);
      }

      cout << binX << "  " << binY << endl;
      cout << "   integral_below_cut " << integral_below_cut << " +/- " << integral_below_cut_error << endl;
      cout << "   integral_above_cut " << integral_above_cut << " +/- " << integral_above_cut_error << endl;
      cout << "   extra_integral     " << extra_integral << " +/- " << extra_integral_error << endl;
      cout << "   efficiency:        " << efficiency <<  "+/-" << efficiency_error << "       simple: " << integral_below_cut/(integral_below_cut + integral_above_cut) << endl;


      // images for testing
      
      TCanvas* canv = new TCanvas();
      canv->SetWindowSize(1000,700);
      canv->SetCanvasSize(1000,700);
      gPad->SetRightMargin(0.4);
      gPad->SetLogy();
      proj_histo->GetXaxis()->SetTitle("m_{T}-m_{0} [GeV/c^{2}]");
      proj_histo->GetYaxis()->SetTitle("Number or Tracks");
      proj_histo->Draw("E");
      canv->SaveAs(Form("/home/matthewharasty/ResearchData/dca_expo_fits/dca_%d_%d_%03d_%03d.png",m_partIndex,m_particleCharge,binX,binY));
      delete proj_histo;
      

    }
  }


}



//_________________________________________________________________________________________________________________________________________________________________________________________________________________


void EmbeddingBinner::closeInputFile(){
  if(m_inFile){
    if(m_inFile->IsOpen()) m_inFile->Close();
    delete m_inFile;
    m_inFile = nullptr;
  }
}

//_________________________________________________________________________________________________________________________________________________________________________________________________________________


void EmbeddingBinner::makeHistograms(){

  #ifdef _EMBEDDING_BINNER_DEBUG_
    cout << "Call to EmbeddingBinner::makeHistograms()" << endl;
  #endif

  m_centEvents = new TH1I("centEvents","Number of Events for Each Centrality Bin; Centrality Bin; Number of Events",
                        16, -0.5, ((double) 16) - 0.5);
 
  //## EVENT LEVEL Histograms###
  m_xyHisto = new TH2I("xyHisto","X vs Y Vertex Location",1000,-5,5,1000,-5,5);
  m_zHisto = NULL;
  if(m_cuts->getEventConfig() == "ColliderCenter"){
    m_zHisto = new TH1I("zVertex","Z Vertex Location",1000,-250,250);
  }
  if(m_cuts->getEventConfig() == "FixedTarget" ){
    m_zHisto  = new TH1I("zVertex","Z Vertex Location",1000,160,240);
  }
  m_tofMatchVsRefMult    = new TH2I("tofMatchVsRefMult","TOF Matches Against RefMult True",1000,0,1000,1000,0,1000);
  m_tofTrayMultVsRefMult = new TH2I("tofTrayMultVsRefMult","TOF Tray Mult Against RefMult True",1000,0,1000,500,0,5000);
  m_numPrimaryDataTracks = new TH1I("numPriamryDataTracks","Number of Real Primary Data Tracks; Number of real primary tracks; Number of events",250,0,2500);
  m_RefMult = new TH1I("refMult_withEmb", "refMult", 1000, 0., 1000.);
  m_RefMultTrue = new TH1I("refMult_withoutEmb", "refMultTrue", 1000, 0., 1000.);

  vector<double> rapidity_binEdges = HistogramUtilities::getBinEdges(m_cuts->getRapMtM0BinStructure(m_partIndex,1));
  vector<double> mTm0_binEdges;
  if(m_useCutClassMtM0Binning){
    mTm0_binEdges = HistogramUtilities::getSmallerBinEdges(m_cuts->getRapMtM0BinStructure(m_partIndex,0), m_mTm0DivNum);
  }else{
    mTm0_binEdges = HistogramUtilities::makeEvenBinEdges(200,0.0,4.0); 
  }

  vector<double> mTm0_binEdges_CutClass = HistogramUtilities::getBinEdges(m_cuts->getRapMtM0BinStructure(m_partIndex,0));
  vector<double> mTm0_binEdges_CutClass_BTOF = HistogramUtilities::getBinEdges(m_cuts->getRapMtM0BinStructure(m_partIndex,4));


  m_vZ_structure = new TH1I("vZ_structure", "Event V_{Z} Structure; V_{Z} cm; Events", g_EmbBin_MaxVzBins, -70, 70); // must change max Vz Bins if you are extending this...


  vector<double>  lossEdges    = HistogramUtilities::makeEvenBinEdges(300,-0.1,0.1);
  vector<double>  pTEdges      = HistogramUtilities::makeEvenBinEdges(150,0.0,3.0);
  vector<double>  refMultEdges = HistogramUtilities::makeEvenBinEdges(511,-0.5,510.5); // must be integer, max of 510 is from 19.6 so you can get efficiencies for lower collider energies
  vector<double>  binEdges_dca = HistogramUtilities::makeEvenBinEdges(32,0,3.657);
  m_nRapidityBins = rapidity_binEdges.size() - 1;
  m_nMtM0Bins     = mTm0_binEdges.size() - 1 ;

  m_embTrackHisto_byRefMult = new TH3I("embTrackHisto_byRefMult",
               Form("EmbTracks Rapidity Vs. m_{T}-m_{%s} %s;y_{%s}^{Emb}; m_{T}^{Emb}-m_{%s} (GeV/c^{2}); refMult",
                 m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                 m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0],511,&refMultEdges[0]);
  m_matchTrackHisto_byRefMult  = new TH3I("matchTrackHisto_byRefMult",
               Form("EmbTracks Rapidity Vs. m_{T}-m_{%s} %s;y_{%s}^{Emb}; m_{T}^{Emb}-m_{%s} (GeV/c^{2}); refMult",
                 m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                 m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0],511,&refMultEdges[0]);

  //CENTRALITY DEPENDANT HISTOGRAMS
  for(int centIndex = 0; centIndex < 16; centIndex++){
     m_embTrackHisto[centIndex] = new TH2I(Form("embTrackHisto_Cent%d",centIndex),
                     Form("EmbTracks Rapidity Vs. m_{T}-m_{%s} %s;y_{%s}^{Emb}; m_{T}^{Emb}-m_{%s} (GeV/c^{2})",
                          m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                     m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);

     m_matchTrackHisto[centIndex] = new TH2I(Form("matchTrackHisto_Cent%d",centIndex),
                       Form("MatchedTracks Rapidity Vs. m_{T}-m_{%s} %s;y_{%s}^{Emb};(m_{T}^{Emb}-m_{%s} (GeV/c^{2})",
                        m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                       m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);

     m_embTrackHisto_pTEta[centIndex] = new TH2I(Form("embTrackHisto_pTEta_Cent%d",centIndex),
                     Form("EmbTracks #eta Vs. p_{T} %s;#eta^{Emb}; p_{T}^{Emb} (GeV/c)",
                          m_pidSymbol.c_str()),
                     m_nEtaBins,m_etaMin,m_etaMax,250,0,2.5);

     m_matchTrackHisto_pTEta[centIndex] = new TH2I(Form("matchTrackHisto_pTEta_Cent%d",centIndex),
                       Form("MatchedTracks #eta Vs. p_{T} %s;#eta^{Emb};p_{T}^{Emb} (GeV/c)",
                        m_pidSymbol.c_str()),
                       m_nEtaBins,m_etaMin,m_etaMax,250,0,2.5);

     m_pTLossHisto[centIndex] = new TH3I(Form("pTLossHisto_Cent%d",centIndex),
                      Form("p_{T} Loss, %s;y_{%s}^{Reco};p_{T}^{Reco} (GeV/c); p_{T}^{Emb}-p_{T}^{Reco} (GeV/c)",m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                      m_nRapidityBins,&rapidity_binEdges[0],pTEdges.size()-1,&pTEdges[0],lossEdges.size()-1,&lossEdges[0]);

     m_mTm0LossHisto[centIndex] = new TH3I(Form("mTm0LossHisto_Cent%d",centIndex),
                       Form("m_{T} Loss, %s;y_{%s}^{Reco};m_{T}^{Reco} - m_{0} (GeV/c^{2});m_{T}^{Emb}-m_{T}^{Reco} (GeV/c^{2})",m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                       m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0],lossEdges.size()-1,&lossEdges[0]);

     m_rapidityLossHisto[centIndex] = new TH3I(Form("rapidityLossHisto_Cent%d",centIndex),
                         Form("Rapidity Loss %s;y_{%s}^{Reco};m_{T}^{Reco} - m_{0} (GeV/c^{2}); y_{%s}^{Emb}-y_{%s}^{Reco} ",
                         m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                         m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0],lossEdges.size()-1,&lossEdges[0]);


     m_btof_tracks_reco[centIndex] = new TH2I(Form("btof_tracks_reco_Cent%d",centIndex),
                                         Form("BTOF Matched Tracks %s;y_{%s}^{Reco};m_{T}^{Reco}-m_{%s} (GeV/c^{2})",
                                              m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                                         m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);

     m_tpc_tracks_forBTOF_reco[centIndex] = new TH2I(Form("tpc_tracks_reco_Cent%d",centIndex),
                                           Form("TPC Matched Tracks %s;y_{%s}^{Reco};m_{T}^{Reco}-m_{%s} (GeV/c^{2})",
                                              m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                                           m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);

     m_btof_tracks_emb[centIndex] = new TH2I(Form("btof_tracks_emb_Cent%d",centIndex),
                                         Form("BTOF Matched Tracks %s;y_{%s}^{Emb};m_{T}^{Emb}-m_{%s} (GeV/c^{2})",
                                              m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                                         m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);

     m_tpc_tracks_forBTOF_emb[centIndex] = new TH2I(Form("tpc_tracks_emb_Cent%d",centIndex),
                                           Form("TPC Matched Tracks %s;y_{%s}^{Emb};m_{T}^{Emb}-m_{%s} (GeV/c^{2})",
                                              m_pidSymbol.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                                           m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);


    m_unweighted_eff_y_mTm0[centIndex] = new TEfficiency(Form("unweighted_eff_y_mTm0_Cent%d",centIndex),
                    Form("Unweighted Centrality %d;y;m_{T}-m_{0} [GeV/c^{2}]",centIndex),m_nRapidityBins,&rapidity_binEdges[0], m_nMtM0Bins,&mTm0_binEdges[0]);
    m_weighted_eff_y_mTm0[centIndex] = new TEfficiency(Form("weighted_eff_y_mTm0_Cent%d",centIndex),
                    Form("Weighted Centrality %d;y;m_{T}-m_{0} [GeV/c^{2}]",centIndex),m_nRapidityBins,&rapidity_binEdges[0], m_nMtM0Bins,&mTm0_binEdges[0]);
    m_weighted_eff_y_mTm0[centIndex]->SetUseWeightedEvents();
    m_eff_y_mTm0_noDCA[centIndex] = new TEfficiency(Form("eff_y_mTm0_noDCA_Cent%d",centIndex),
                    Form("Efficiency without DCA Centrality %d;y;m_{T}-m_{0} [GeV/c^{2}]",centIndex),m_nRapidityBins,&rapidity_binEdges[0], m_nMtM0Bins,&mTm0_binEdges[0]);
    m_weighted_eff_y_mTm0_noDCA[centIndex] = new TEfficiency(Form("weighted_eff_y_mTm0_noDCA_Cent%d",centIndex),
                    Form("Weighted Efficiency without DCA Centrality %d;y;m_{T}-m_{0} [GeV/c^{2}]",centIndex),m_nRapidityBins,&rapidity_binEdges[0], m_nMtM0Bins,&mTm0_binEdges[0]);
    m_dca_eff_y_mTm0[centIndex] = new TEfficiency(Form("dca_eff_y_mTm0_Cent%d",centIndex),
                    Form("DCA Cut Efficiency Centrality %d;y;m_{T}-m_{0} [GeV/c^{2}]",centIndex),m_nRapidityBins,&rapidity_binEdges[0], m_nMtM0Bins,&mTm0_binEdges[0]);
    m_eff_y_mTm0_DCA[centIndex] = new TEfficiency(Form("eff_y_mTm0_dca_Cent%d",centIndex),
                    Form("Efficiency Centrality %d;y;m_{T}-m_{0} [GeV/c^{2}];gDCA [cm]",centIndex),m_nRapidityBins,&rapidity_binEdges[0], m_nMtM0Bins,&mTm0_binEdges[0],((int)binEdges_dca.size())-1,&binEdges_dca[0]);
  }
  
  m_emb_weights_y_mTm0_DCA = new TH3I("emb_weights_y_mTm0_DCA","Embedding;y;m_{T}-m_{0} [GeV/c^{2}]; DCA [cm]",
                    m_nRapidityBins,&rapidity_binEdges[0], m_nMtM0Bins,&mTm0_binEdges[0], ((int)binEdges_dca.size())-1,&binEdges_dca[0]);

  for(int vZ_index = 0; vZ_index < g_EmbBin_MaxVzBins; vZ_index++){
    string vZInfo = Form("V_{Z}=[%f,%f]",m_vZ_structure->GetBinLowEdge(vZ_index+1),m_vZ_structure->GetBinLowEdge(vZ_index+2));
    m_embTrackHisto_byRefMult_byVz[vZ_index] = new TH3I(Form("embTrackHisto_byRefMult_Vz%03d",vZ_index),
                 Form("EmbTracks Rapidity Vs. m_{T}-m_{%s} %s %s;y_{%s}^{Emb}; m_{T}^{Emb}-m_{%s} (GeV/c^{2}); refMult",
                   m_pidSymbol.c_str(),m_pidSymbol.c_str(),vZInfo.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                   m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0],511,&refMultEdges[0]);
    m_matchTrackHisto_byRefMult_byVz[vZ_index]  = new TH3I(Form("matchTrackHisto_byRefMult_Vz%03d",vZ_index),
                 Form("EmbTracks Rapidity Vs. m_{T}-m_{%s} %s %s;y_{%s}^{Emb}; m_{T}^{Emb}-m_{%s} (GeV/c^{2}); refMult",
                   m_pidSymbol.c_str(),m_pidSymbol.c_str(),vZInfo.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                   m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0],511,&refMultEdges[0]);

    //CENTRALITY DEPENDANT HISTOGRAMS
    for(int centIndex = 0; centIndex < 16; centIndex++){
       m_embTrackHisto_byVz[centIndex][vZ_index] = new TH2I(Form("embTrackHisto_Cent%d_Vz%03d",centIndex,vZ_index),
                       Form("EmbTracks Rapidity Vs. m_{T}-m_{%s} %s %s;y_{%s}^{Emb}; m_{T}^{Emb}-m_{%s} (GeV/c^{2})",
                            m_pidSymbol.c_str(),m_pidSymbol.c_str(),vZInfo.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                       m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);

       m_matchTrackHisto_byVz[centIndex][vZ_index] = new TH2I(Form("matchTrackHisto_Cent%d_Vz%03d",centIndex,vZ_index),
                         Form("MatchedTracks Rapidity Vs. m_{T}-m_{%s} %s %s;y_{%s}^{Emb};(m_{T}^{Emb}-m_{%s} (GeV/c^{2})",
                          m_pidSymbol.c_str(),m_pidSymbol.c_str(),vZInfo.c_str(),m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                         m_nRapidityBins,&rapidity_binEdges[0],m_nMtM0Bins,&mTm0_binEdges[0]);

       m_embTrackHisto_pTEta_byVz[centIndex][vZ_index] = new TH2I(Form("embTrackHisto_pTEta_Cent%d_Vz%03d",centIndex,vZ_index),
                       Form("EmbTracks #eta Vs. p_{T} %s %s;#eta^{Emb}; p_{T}^{Emb} (GeV/c)",
                            m_pidSymbol.c_str(),vZInfo.c_str()),
                       m_nEtaBins,m_etaMin,m_etaMax,250,0,2.5);

       m_matchTrackHisto_pTEta_byVz[centIndex][vZ_index] = new TH2I(Form("matchTrackHisto_pTEta_Cent%d_Vz%03d",centIndex,vZ_index),
                         Form("MatchedTracks #eta Vs. p_{T} %s %s;#eta^{Emb};p_{T}^{Emb} (GeV/c)",
                          m_pidSymbol.c_str(),vZInfo.c_str()),
                         m_nEtaBins,m_etaMin,m_etaMax,250,0,2.5);
    }
  }













  //double mTm0BinWidth_Response = 0.025;
  //double mTm0Max_Response = 
  vector <double> binEdges_mTm0_repsonse = HistogramUtilities::makeEvenBinEdges(250,0,4.5);
  vector <double> binEdges_mTm0_repsonse_emb = HistogramUtilities::makeEvenBinEdges(550,0,4.5);
  vector <double> binEdges_pT_repsonse = HistogramUtilities::makeEvenBinEdges(50,0,2.0);
  vector <double> binEdges_pT_loss_repsonse = HistogramUtilities::makeEvenBinEdges(400,-0.1,0.1);


  m_mTm0ReponseMatrixHisto  = new TH3I("mTm0ResponseMatrix",
                       Form("Response Matrix %s; y_{%s}^{Reco}; p_{T}^{Reco} (GeV/c); p_{T}^{Emb} (GeV/c)",m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                       m_nRapidityBins,&rapidity_binEdges[0], ((int)binEdges_mTm0_repsonse.size()) - 1, &binEdges_mTm0_repsonse[0], ((int)binEdges_mTm0_repsonse.size())-1,&binEdges_mTm0_repsonse[0]);

  /*m_mTm0ReponseMatrixHistoBTOF = new TH3I("mTm0ResponseMatrixBTOF",
                       Form("Response Matrix %s; y_{%s}^{Reco}; m_{T}^{Reco} - m_{0} (GeV/c^{2}); m_{T}^{Emb} - m_{0} (GeV/c^{2})",m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                       m_nRapidityBins,&rapidity_binEdges[0],((int)mTm0_binEdges_CutClass_BTOF.size()) - 1,&mTm0_binEdges_CutClass_BTOF[0],((int)binEdges_mTm0_repsonse.size())-1,&binEdges_mTm0_repsonse[0]);
  */

  m_pTLossReponseMatrixHisto = new TH3I("pTLossResponseMatrix",
                       Form("Response Matrix %s; y_{%s}^{Reco}; p_{T}^{Reco} [GeV/c]; p_{T}^{Emb} - p_{T}^{Reco} (GeV/c)",m_pidSymbol.c_str(),m_pidSymbol.c_str()),
                       m_nRapidityBins,&rapidity_binEdges[0], ((int)binEdges_pT_repsonse.size()) - 1, &binEdges_pT_repsonse[0], ((int)binEdges_pT_loss_repsonse.size())-1,&binEdges_pT_loss_repsonse[0]);



  // the reconstructed tracks that were embedded
  m_track_NHits_reco      = new TH1I("track_NHits_reco","track_NHits_reco", 100, 0, 100);
  m_track_NHitsFit_reco   = new TH1I("track_NHitsFit_reco","track_NHitsFit_reco",100, 0, 100);
  m_track_NHitsPoss_reco  = new TH1I("track_NHitsPoss_reco","track_NHitsPoss_reco",100, 0, 100);
  m_track_NHitsDedx_reco  = new TH1I("track_NHitsDedx_reco","track_NHitsDedx_reco",100, 0, 100);
  m_track_NHitsFitOverPoss_reco   = new TH1I("track_NHitsFitOverPoss_reco","track_NHitsFitOverPoss_reco",100, 0, 1.0);
  m_track_pT_reco     = new TH1I("track_pT_reco","track_pT_reco",300,0,6.0);
  m_track_mTm0_reco   = new TH1I("track_mTm0_reco","track_mTm0_reco", 300 , 0 , 6.0);
  m_track_eta_reco    = new TH1I("track_eta_reco","track_eta_reco", 200, -4.0,4.0);
  m_track_rapidity_reco   = new TH1I("track_rapidity_reco","track_rapidity_reco",200, -4.0,4.0);
  m_track_phi_reco    = new TH1I("track_phi_reco","track_phi_reco", 100, -3.2,3.2);
  m_track_pTot_reco   = new TH1I("track_pTot_reco","track_pTot_reco", 325, 0, 6.5);
  m_track_pZ_reco     = new TH1I("track_pZ_reco","track_pZ_reco", 200, 0,6.0);
  m_track_gDCA_reco   = new TH1I("track_gDCA_reco","track_gDCA_reco",200 , 0,3.5);
  
  //embeded recoa (the simulated variables, not reconstructed variables)
  m_track_pT_emb = new TH1I("track_pT_emb","track_pT_emb",300,0,6.0);
  m_track_mTm0_emb   = new TH1I("track_mTm0_emb","track_mTm0_emb", 300 , 0 , 6.0);
  m_track_eta_emb = new TH1I("track_eta_emb","track_eta_emb", 200, -4.0,4.0);
  m_track_rapidity_emb   = new TH1I("track_rapidity_emb","track_rapidity_emb",200, -4.0,4.0);
  m_track_phi_emb = new TH1I("track_phi_emb","track_phi_emb", 100, -3.2,3.2);
  m_track_pTot_emb   = new TH1I("track_pTot_emb","track_pTot_emb", 325, 0, 6.5);
  m_track_pZ_emb = new TH1I("track_pZ_emb ","track_pZ_emb ", 200, 0,6.0);

}




//_________________________________________________________________________________________________________________________________________________________________________________________________________________



void EmbeddingBinner::fillHistograms(){
  #ifdef _EMBEDDING_BINNER_DEBUG_
    cout << "Call to EmbeddingBinner::fillHistograms()" << endl;
  #endif

  int numEvents = m_tree->GetEntries();
  for(int eventIndex = 0; eventIndex < numEvents; eventIndex++){
    m_tree->GetEntry(eventIndex);
    if(!m_event){
      cerr << "ERROR: No MattEvent!    Ptr: " << m_event << endl;
      continue;
    }



    //#############################################    EVENT LEVEL CUTS   ##########################################################

    //NOTE: StRefMultCorr cuts out bad runs in the making of the TTree

    int centralityIndex = m_event->centIndex();
    if(centralityIndex < 0 || centralityIndex > 15){
      #ifdef _EMBEDDING_BINNER_DEBUG_
        cout << "Bad Centrality - CentBin: " << centralityIndex << "    refMult_dataOnly = " << m_event->refMult_data_only() << endl;
      #endif
      continue;
    }
    if(m_cuts->isPileUp(m_event->refMult_data_only(),m_event->numTofMatch())){
      #ifdef _EMBEDDING_BINNER_DEBUG_
        cout << "  Event rejected becuase of pileup: refMult_dataOnly = " << m_event->refMult_data_only() << "   and bTof Match: " << m_event->numTofMatch() << endl; 
      #endif
      continue;
    }

    //##################   GOOD EVENT CUTS ######################
    if(!m_cuts->isGoodEvent(m_event->vertex_X(), m_event->vertex_Y(), m_event->vertex_Z(), m_event->numTofMatch(), m_event->numTofT0())){
      #ifdef _EMBEDDING_BINNER_DEBUG_
        cout << " Rejecting Event by isGoodEvent  " << endl;
      #endif
      continue;
    }

     
    //##################   TRIGGER CUTS ######################
    // N/A  (they only embed into mb)

    int vZ_index = m_vZ_structure->FindBin(m_event->vertex_Z()) - 1;
    bool valid_vZ_index = (vZ_index >= 0 && vZ_index < g_EmbBin_MaxVzBins );


    m_centEvents->Fill(centralityIndex);
    m_xyHisto->Fill(m_event->vertex_X(), m_event->vertex_Y());
    m_zHisto->Fill(m_event->vertex_Z());
    m_tofMatchVsRefMult->Fill(m_event->refMult_data_only(),m_event->numTofMatch());
    m_numPrimaryDataTracks->Fill(m_event->numPrimaryTracks());
    m_RefMultTrue->Fill(m_event->refMult_data_only());
    m_RefMult->Fill(m_event->refMult());


    
    int numTracks = m_event->numMcTracks();
    for(int trackIndex = 0; trackIndex < numTracks; trackIndex++){
      m_track = m_event->getTrack(trackIndex);

      if(!m_track){
        cout << "ERROR: Track not found! Index: " << trackIndex << endl;
      }

      double mTm0_emb = m_track->mTm0_emb(m_mass);
      double rapidity_emb = m_track->rapidity_emb(m_mass);

      m_track_pT_emb->Fill(      m_track->pT_emb());
      m_track_mTm0_emb->Fill(    mTm0_emb);
      m_track_eta_emb->Fill(     m_track->eta_emb());
      m_track_rapidity_emb->Fill(rapidity_emb);
      m_track_phi_emb->Fill(     m_track->phi_emb());
      m_track_pTot_emb->Fill(    m_track->pTot());
      m_track_pZ_emb->Fill(      m_track->pZ());

      m_embTrackHisto[centralityIndex]->Fill(       rapidity_emb,mTm0_emb);
      m_embTrackHisto_pTEta[centralityIndex]->Fill( m_track->eta_emb(),     m_track->pT_emb());
      m_embTrackHisto_byRefMult->Fill(              rapidity_emb,mTm0_emb,(m_cuts->getEventConfig() == "ColliderCenter" ? m_event->refMult_data_only() : m_event->numPrimaryTracks()));



      if(valid_vZ_index){
        m_embTrackHisto_byVz[centralityIndex][vZ_index]->Fill(       rapidity_emb,mTm0_emb);
        m_embTrackHisto_pTEta_byVz[centralityIndex][vZ_index]->Fill( m_track->eta_emb(),     m_track->pT_emb());
        m_embTrackHisto_byRefMult_byVz[vZ_index]->Fill(              rapidity_emb,mTm0_emb,(m_cuts->getEventConfig() == "ColliderCenter" ? m_event->refMult_data_only() : m_event->numPrimaryTracks()));
      }

      bool passedTrackCuts = false;
      //bool passedTrackCuts_withoutDCA = false;
      bool goodDCA = (m_track->gDCA() < m_cuts->getDCACut());
      double mTm0 = m_track->mTm0(m_mass);
      double rapidity = m_track->rapidity(m_mass);
      if(m_track->wasReconstructed()){
        passedTrackCuts = m_cuts->isGoodTrack(true,m_track->nHitsFit(), m_track->nHitsFit()/m_track->nHitsRatio(), m_track->nHitsDeDx(),m_track->gDCA(), m_track->phi());
        //passedTrackCuts_withoutDCA = m_cuts->isGoodTrack(true,m_track->nHitsFit(), m_track->nHitsFit()/m_track->nHitsRatio(), m_track->nHitsDeDx(),0.0, m_track->phi());
        if(passedTrackCuts){

          m_track_NHitsFit_reco->Fill(m_track->nHitsFit());
          m_track_NHitsPoss_reco->Fill(m_track->nHitsFit()/m_track->nHitsRatio());
          m_track_NHitsDedx_reco->Fill(m_track->nHitsDeDx());
          m_track_NHitsFitOverPoss_reco->Fill(m_track->nHitsRatio());
          m_track_pT_reco->Fill(m_track->pT());
          m_track_mTm0_reco->Fill(mTm0);
          m_track_eta_reco->Fill(m_track->eta());
          m_track_rapidity_reco->Fill(rapidity);
          m_track_phi_reco->Fill(m_track->phi());
          m_track_pTot_reco->Fill(m_track->pTot());
          m_track_pZ_reco->Fill(m_track->pZ());
          m_track_gDCA_reco->Fill(m_track->gDCA());

          //The filled tracks are still in embedded pt and eta (not using muTrack on purpose)
          m_matchTrackHisto[centralityIndex]->Fill(m_track->rapidity_emb(m_mass), m_track->mTm0_emb(m_mass)); //this should be emb (MC)
          m_matchTrackHisto_pTEta[centralityIndex]->Fill(m_track->eta_emb(),m_track->pT_emb());
          m_matchTrackHisto_byRefMult->Fill(m_track->rapidity_emb(m_mass), m_track->mTm0_emb(m_mass), (m_cuts->getEventConfig() == "ColliderCenter" ? m_event->refMult_data_only() : m_event->numPrimaryTracks()));
          if(valid_vZ_index){
            m_matchTrackHisto_byVz[centralityIndex][vZ_index]->Fill(m_track->rapidity_emb(m_mass), m_track->mTm0_emb(m_mass)); //this should be emb (MC)
            m_matchTrackHisto_pTEta_byVz[centralityIndex][vZ_index]->Fill(m_track->eta_emb(),m_track->pT_emb());
            m_matchTrackHisto_byRefMult_byVz[vZ_index]->Fill(m_track->rapidity_emb(m_mass), m_track->mTm0_emb(m_mass), (m_cuts->getEventConfig() == "ColliderCenter" ? m_event->refMult_data_only() : m_event->numPrimaryTracks()));
          }


          m_pTLossHisto[centralityIndex]->Fill(      rapidity, m_track->pT(), m_track->pT_emb() - m_track->pT()); // embed (true) - reco (measured)
          m_mTm0LossHisto[centralityIndex]->Fill(    rapidity, mTm0,      m_track->mTm0_emb(m_mass) - mTm0);
          m_mTm0ReponseMatrixHisto->Fill(rapidity, m_track->pT(), m_track->pT_emb());
          m_pTLossReponseMatrixHisto->Fill(rapidity, m_track->pT(), m_track->pT_emb() - m_track->pT());
          //m_mTm0ReponseMatrixHistoBTOF->Fill(rapidity, m_track->pT(m_mass), m_track->pT_emb(m_mass));
          m_rapidityLossHisto[centralityIndex]->Fill(rapidity, mTm0,     m_track->rapidity_emb(m_mass) - rapidity);
    
          m_tpc_tracks_forBTOF_reco[centralityIndex]->Fill(rapidity,mTm0);
          m_tpc_tracks_forBTOF_emb[centralityIndex] ->Fill(m_track->rapidity_emb(m_mass),m_track->mTm0_emb(m_mass));
          if(m_cuts->isGoodBTof(m_track->btofLocalY(), m_track->btofLocalZ())){
            m_btof_tracks_reco[centralityIndex]->Fill(rapidity,mTm0);
            m_btof_tracks_emb[centralityIndex]->Fill(m_track->rapidity_emb(m_mass),m_track->mTm0_emb(m_mass));
          }
        } // passed track cuts
      } // End Was Reconstructed

      m_unweighted_eff_y_mTm0[centralityIndex]->Fill(passedTrackCuts,rapidity,mTm0);
      double track_weight = 1.0;
      if(m_emb_weights_y_mTm0_DCA && m_data_weights_y_mTm0_DCA){
        track_weight = 1.0;
      }
      m_weighted_eff_y_mTm0[centralityIndex]->FillWeighted(passedTrackCuts,track_weight,rapidity,mTm0);

      if(goodDCA){
        m_eff_y_mTm0_noDCA[centralityIndex]->Fill(passedTrackCuts,rapidity,mTm0);
        m_weighted_eff_y_mTm0_noDCA[centralityIndex]->FillWeighted(passedTrackCuts,track_weight,rapidity,mTm0);
      }
      m_dca_eff_y_mTm0[centralityIndex]->Fill(goodDCA,rapidity,mTm0);
      m_eff_y_mTm0_DCA[centralityIndex]->Fill(passedTrackCuts,rapidity,mTm0,m_track->gDCA());



    }// End Track Loop
  }// end of event loop
}


//_________________________________________________________________________________________________________________________________________________________________________________________________________________


void EmbeddingBinner::fillWeightedHistograms(){
  #ifdef _EMBEDDING_BINNER_DEBUG_
    cout << "Call to EmbeddingBinner::fillWeightedHistograms()" << endl;
  #endif

  int numEvents = m_tree->GetEntries();
  for(int eventIndex = 0; eventIndex < numEvents; eventIndex++){
    m_tree->GetEntry(eventIndex);
    if(!m_event){
      cerr << "ERROR: No MattEvent!    Ptr: " << m_event << endl;
      continue;
    }

    //#############################################    EVENT LEVEL CUTS   ##########################################################

    //NOTE: StRefMultCorr cuts out bad runs in the making of the TTree

    int centralityIndex = m_event->centIndex();
    if(centralityIndex < 0 || centralityIndex > 15){
      #ifdef _EMBEDDING_BINNER_DEBUG_
        cout << "Bad Centrality - CentBin: " << centralityIndex << "    refMult_dataOnly = " << m_event->refMult_data_only() << endl;
      #endif
      continue;
    }
    if(m_cuts->isPileUp(m_event->refMult_data_only(),m_event->numTofMatch())){
      #ifdef _EMBEDDING_BINNER_DEBUG_
        cout << "  Event rejected becuase of pileup: refMult_dataOnly = " << m_event->refMult_data_only() << "   and bTof Match: " << m_event->numTofMatch() << endl; 
      #endif
      continue;
    }

    //##################   GOOD EVENT CUTS ######################
    if(!m_cuts->isGoodEvent(m_event->vertex_X(), m_event->vertex_Y(), m_event->vertex_Z(), m_event->numTofMatch(), m_event->numTofT0())){
      #ifdef _EMBEDDING_BINNER_DEBUG_
        cout << " Rejecting Event by isGoodEvent  " << endl;
      #endif
      continue;
    }

     
    //##################   TRIGGER CUTS ######################
    // N/A  (they only embed into mb)

    int numTracks = m_event->numMcTracks();
    for(int trackIndex = 0; trackIndex < numTracks; trackIndex++){
      m_track = m_event->getTrack(trackIndex);

      if(!m_track){
        cout << "ERROR: Track not found! Index: " << trackIndex << endl;
      }

      //bool passedTrackCuts = m_track->wasReconstructed() && m_cuts->isGoodTrack(true,m_track->nHitsFit(), m_track->nHitsFit()/m_track->nHitsRatio(), m_track->nHitsDeDx(),m_track->gDCA(), m_track->phi());
      bool passedTrackCuts_withoutDCA = m_track->wasReconstructed() && m_cuts->isGoodTrack(true,m_track->nHitsFit(), m_track->nHitsFit()/m_track->nHitsRatio(), m_track->nHitsDeDx(),0.0, m_track->phi());
      if(passedTrackCuts_withoutDCA) m_emb_weights_y_mTm0_DCA->Fill(m_track->rapidity(m_mass), m_track->mTm0(m_mass),m_track->gDCA());
    }// End Track Loop
  }// end of event loop
}

//_________________________________________________________________________________________________________________________________________________________________________________________________________________


void EmbeddingBinner::write(string a_outFileName){
  #ifdef _EMBEDDING_BINNER_DEBUG_
    cout << "Call to EmbeddingBinner::write(" << a_outFileName << ")" << endl;;
  #endif
  TFile* outFile = new TFile(a_outFileName.c_str(), "RECREATE");
  if(!outFile->IsOpen()){
    cout << "ERROR: Output File is broken... Can't make or can't read TFile " << endl;
  }

  HistogramUtilities::ConditionalWrite(m_cuts->getEventRejectionHistogram());
  HistogramUtilities::ConditionalWrite(m_cuts->getTrackRejectionHistogram());
  HistogramUtilities::ConditionalWrite(m_centEvents);
  HistogramUtilities::ConditionalWrite(m_xyHisto);
  HistogramUtilities::ConditionalWrite(m_zHisto);
  HistogramUtilities::ConditionalWrite(m_tofMatchVsRefMult);
  HistogramUtilities::ConditionalWrite(m_tofTrayMultVsRefMult);
  HistogramUtilities::ConditionalWrite(m_numPrimaryDataTracks);
  HistogramUtilities::ConditionalWrite(m_RefMult);
  HistogramUtilities::ConditionalWrite(m_RefMultTrue);
  
  // the reconstructed tracks 
  HistogramUtilities::ConditionalWrite(m_track_NHitsFit_reco);
  HistogramUtilities::ConditionalWrite(m_track_NHitsDedx_reco);
  HistogramUtilities::ConditionalWrite(m_track_NHitsFitOverPoss_reco);
  HistogramUtilities::ConditionalWrite(m_track_pT_reco);
  HistogramUtilities::ConditionalWrite(m_track_mTm0_reco);
  HistogramUtilities::ConditionalWrite(m_track_eta_reco);
  HistogramUtilities::ConditionalWrite(m_track_rapidity_reco);
  HistogramUtilities::ConditionalWrite(m_track_phi_reco);
  HistogramUtilities::ConditionalWrite(m_track_pTot_reco);
  HistogramUtilities::ConditionalWrite(m_track_pZ_reco);
  HistogramUtilities::ConditionalWrite(m_track_gDCA_reco);
  
  //embeded track histograms (the simulated variables, not reconstructed variables)
  HistogramUtilities::ConditionalWrite(m_track_pT_emb);
  HistogramUtilities::ConditionalWrite(m_track_mTm0_emb);
  HistogramUtilities::ConditionalWrite(m_track_eta_emb);
  HistogramUtilities::ConditionalWrite(m_track_rapidity_emb);
  HistogramUtilities::ConditionalWrite(m_track_phi_emb);
  HistogramUtilities::ConditionalWrite(m_track_pTot_emb);
  HistogramUtilities::ConditionalWrite(m_track_pZ_emb);

  outFile->cd();
  outFile->mkdir("EfficiencyHistograms");
  outFile->cd("EfficiencyHistograms");
  for(int centIndex = 0; centIndex < m_nCentBins; centIndex++){
    HistogramUtilities::ConditionalWrite(m_embTrackHisto[centIndex]);
    HistogramUtilities::ConditionalWrite(m_matchTrackHisto[centIndex]);
    HistogramUtilities::ConditionalWrite(m_embTrackHisto_pTEta[centIndex]);
    HistogramUtilities::ConditionalWrite(m_matchTrackHisto_pTEta[centIndex]);
    HistogramUtilities::ConditionalWrite(m_tpc_tracks_forBTOF_reco[centIndex]);
    HistogramUtilities::ConditionalWrite(m_btof_tracks_reco[centIndex]);
    HistogramUtilities::ConditionalWrite(m_btof_tracks_emb[centIndex]);
    HistogramUtilities::ConditionalWrite(m_tpc_tracks_forBTOF_emb[centIndex]);
  }  
  m_embTrackHisto_byRefMult->Write();
  m_matchTrackHisto_byRefMult->Write();
  m_dca_eff_y_mTm0_data->Write();
  for(int centIndex = 0; centIndex < m_nCentBins; centIndex++){
    HistogramUtilities::ConditionalWrite(m_unweighted_eff_y_mTm0[centIndex]);
    HistogramUtilities::ConditionalWrite(m_weighted_eff_y_mTm0[centIndex]);
    HistogramUtilities::ConditionalWrite(m_eff_y_mTm0_noDCA[centIndex]);
    HistogramUtilities::ConditionalWrite(m_weighted_eff_y_mTm0_noDCA[centIndex]);
    HistogramUtilities::ConditionalWrite(m_dca_eff_y_mTm0[centIndex]);
    HistogramUtilities::ConditionalWrite(m_eff_y_mTm0_DCA[centIndex]);
  }

  HistogramUtilities::ConditionalWrite(m_vZ_structure);
  for(int vZ_index = 0; vZ_index < g_EmbBin_MaxVzBins; vZ_index++){
    HistogramUtilities::ConditionalWrite(m_embTrackHisto_byRefMult_byVz[vZ_index]);
    HistogramUtilities::ConditionalWrite(m_matchTrackHisto_byRefMult_byVz[vZ_index]);
    for(int centIndex = 0; centIndex < 16; centIndex++){
       HistogramUtilities::ConditionalWrite(m_embTrackHisto_byVz[centIndex][vZ_index]);
       HistogramUtilities::ConditionalWrite(m_matchTrackHisto_byVz[centIndex][vZ_index]);
       HistogramUtilities::ConditionalWrite(m_embTrackHisto_pTEta_byVz[centIndex][vZ_index]); 
       HistogramUtilities::ConditionalWrite(m_matchTrackHisto_pTEta_byVz[centIndex][vZ_index]);
    }
  }

  outFile->cd();
  outFile->mkdir("EnergyLossHistograms");
  outFile->cd("EnergyLossHistograms");
  for(int centIndex=0; centIndex<m_nCentBins; centIndex++){
    HistogramUtilities::ConditionalWrite(m_pTLossHisto[centIndex]);
    HistogramUtilities::ConditionalWrite(m_mTm0LossHisto[centIndex]);
    HistogramUtilities::ConditionalWrite(m_rapidityLossHisto[centIndex]);
  }
  outFile->cd();
  HistogramUtilities::ConditionalWrite(m_mTm0ReponseMatrixHisto);
  HistogramUtilities::ConditionalWrite(m_pTLossReponseMatrixHisto);
  //HistogramUtilities::ConditionalWrite(m_mTm0ReponseMatrixHistoBTOF);

  if(outFile){
    outFile->Close();
  }
}

#endif