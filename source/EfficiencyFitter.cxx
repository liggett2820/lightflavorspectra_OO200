#include "../headers/EfficiencyFitter.h"
#ifndef _EFF_FITTER_CXX_
#define _EFF_FITTER_CXX_


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


EfficiencyFitter::EfficiencyFitter(){
  gStyle->SetOptFit(1);
  TVirtualFitter::SetMaxIterations(10000);

  m_partInfo = new ParticleInfo();
  //m_cuts = a_cuts;
  m_uniqueIdentifier = 0;
  m_ycm = 0.0;
  m_testingOneCent = false;
  m_testingMidRapidity = false;
  m_makeChrisConfidenceGraphs = false;
  m_pTLossMode = false;
  m_useIntegralChiSqr = true;
  m_useChebyshevSpacing = false;
  m_useEnergyLossTProfiles = true; //CHANGED THIS ONE 
  m_fixed_target_mode = false;

  m_tpc_eff_no_fit_num_bins = 80;
  m_btof_eff_no_fit_num_bins = 80;

  m_currentFitFunction = NULL;
  m_currentGraphToFit  = NULL;
  m_currentProtonDCA   = NULL;
  m_currentAntiProtonDCA = NULL;
  m_knockout_mTm0_structure = NULL;
  //m_embedding_vZ_structure = NULL;


  //   ##############  NULL ALL THE THINGS  ##############

  HistogramUtilities::nullptr_array_1D(m_centEvents,9);
  HistogramUtilities::nullptr_array_1D(m_embedding_rapidity_structure,9);
  HistogramUtilities::nullptr_array_1D(m_embedding_pseudorapidity_structure,9);
  HistogramUtilities::nullptr_array_1D(m_hybridFeedDown_rapidity_structure,9);
  HistogramUtilities::nullptr_array_2D(m_picoBinner_rapidity_structure,9,3);
  HistogramUtilities::nullptr_array_2D(m_mTm0ReponseMatrixHisto,9,2);
  HistogramUtilities::nullptr_array_2D(m_mTm0ReponseMatrixHisto_PDF,9,2);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Data_AllCent_Plus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Data_AllCent_Minus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Data_AllCent_Combined,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Emb_AllCent_Plus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Emb_AllCent_Minus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Emb_AllCent_Combined,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Data_ChiSqrd_AllCent_Plus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Data_ChiSqrd_AllCent_Minus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Data_ChiSqrd_AllCent_Combined,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Emb_ChiSqrd_AllCent_Plus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Emb_ChiSqrd_AllCent_Minus,9);
  HistogramUtilities::nullptr_array_1D(m_BTOFEff_Emb_ChiSqrd_AllCent_Combined,9);

  HistogramUtilities::nullptr_array_2D(m_BTOF_fromData_Efficiency_ChiSqrd,9,2);
  HistogramUtilities::nullptr_array_2D(m_BTOF_Efficiency_ChiSqrd,9,2);
  HistogramUtilities::nullptr_array_2D(m_TPC_Efficiency_ChiSqrd,9,2);
  HistogramUtilities::nullptr_array_2D(m_TPC_Efficiency_ChiSqrd_Sys,9,2);
  HistogramUtilities::nullptr_array_2D(m_EnergyLoss_ChiSqrd,9,2);
  HistogramUtilities::nullptr_array_2D(m_EnergyLoss_ChiSqrd_Sys,9,2);
  HistogramUtilities::nullptr_array_2D(m_TPC_Efficiency_ChiSqrd_Compare,9,2);

  for(unsigned int partIndex = 0; partIndex < 9; partIndex++){

    /*m_centEvents[partIndex] = NULL;
    m_picoBinner_rapidity_structure[partIndex][0]   = NULL;
    m_picoBinner_rapidity_structure[partIndex][1]   = NULL;
    m_picoBinner_rapidity_structure[partIndex][2]   = NULL;
    m_embedding_rapidity_structure[partIndex]       = NULL;
    m_embedding_pseudorapidity_structure[partIndex] = NULL;
    m_hybridFeedDown_rapidity_structure[partIndex]  = NULL;

    m_mTm0ReponseMatrixHisto[partIndex][0]     = NULL;
    m_mTm0ReponseMatrixHisto[partIndex][1]     = NULL;
    m_mTm0ReponseMatrixHisto_PDF[partIndex][0] = NULL;
    m_mTm0ReponseMatrixHisto_PDF[partIndex][1] = NULL;
    */

    //m_embTrackHisto2D_AllCent_Plus[partIndex] = NULL;
    //m_matchTrackHisto2D_AllCent_Plus[partIndex] = NULL;
    //m_effHisto2DCent_Plus[partIndex] = NULL;
    //m_pTLossHisto3D_AllCent_Plus[partIndex] = NULL;
    //m_embTrackHisto2D_AllCent_Minus[partIndex] = NULL;
    //m_matchTrackHisto2D_AllCent_Minus[partIndex] = NULL;
    //m_effHisto2DCent_Minus[partIndex] = NULL;
    //m_pTLossHisto3D_AllCent_Minus[partIndex] = NULL;

    //m_TPC_Eff_ChiSqrd_AllCent_Plus[partIndex] = NULL; // part
    //m_TPC_Eff_ChiSqrd_AllCent_Minus[partIndex] = NULL; // part
    /*m_BTOFEff_Data_AllCent_Plus[partIndex] = NULL;
    m_BTOFEff_Data_AllCent_Minus[partIndex] = NULL;
    m_BTOFEff_Data_AllCent_Combined[partIndex] = NULL;
    m_BTOFEff_Emb_AllCent_Plus[partIndex] = NULL;
    m_BTOFEff_Emb_AllCent_Minus[partIndex] = NULL;
    m_BTOFEff_Emb_AllCent_Combined[partIndex] = NULL;

    m_BTOFEff_Data_ChiSqrd_AllCent_Plus[partIndex] = NULL; // part rap
    m_BTOFEff_Data_ChiSqrd_AllCent_Minus[partIndex] = NULL; // part rap
    m_BTOFEff_Data_ChiSqrd_AllCent_Combined[partIndex] = NULL; // part rap
    m_BTOFEff_Emb_ChiSqrd_AllCent_Plus[partIndex] = NULL; // part rap
    m_BTOFEff_Emb_ChiSqrd_AllCent_Minus[partIndex] = NULL; // part rap
    m_BTOFEff_Emb_ChiSqrd_AllCent_Combined[partIndex] = NULL; // part rap


    m_BTOF_fromData_Efficiency_ChiSqrd[partIndex][0] = NULL;
    m_BTOF_fromData_Efficiency_ChiSqrd[partIndex][1] = NULL;
    m_BTOF_Efficiency_ChiSqrd[partIndex][0] = NULL;
    m_BTOF_Efficiency_ChiSqrd[partIndex][1] = NULL;
    m_TPC_Efficiency_ChiSqrd[partIndex][0] = NULL;
    m_TPC_Efficiency_ChiSqrd[partIndex][1] = NULL;
    m_TPC_Efficiency_ChiSqrd_Sys[partIndex][0] = NULL;
    m_TPC_Efficiency_ChiSqrd_Sys[partIndex][1] = NULL;
    m_EnergyLoss_ChiSqrd[partIndex][0] = NULL;
    m_EnergyLoss_ChiSqrd[partIndex][1] = NULL;
    m_EnergyLoss_ChiSqrd_Sys[partIndex][0] = NULL;
    m_EnergyLoss_ChiSqrd_Sys[partIndex][1] = NULL;
    m_TPC_Efficiency_ChiSqrd_Compare[partIndex][0] = NULL;
    m_TPC_Efficiency_ChiSqrd_Compare[partIndex][1] = NULL;
    */

    m_particleLoaded_feeddown[partIndex][0] = false;
    m_particleLoaded_feeddown[partIndex][1] = false;
    m_polynomialOrderForTPCEff[partIndex][0] = 12;
    m_polynomialOrderForBTOFEff[partIndex][0] = 13;
    m_polynomialOrderForTPCEff[partIndex][1] = 12;
    m_polynomialOrderForBTOFEff[partIndex][1] = 13;

    for(unsigned int centIndex = 0; centIndex < 16; centIndex++){

      //m_embTrackHisto2D_Plus[partIndex][centIndex] = NULL; // particle, centrality
      //m_matchTrackHisto2D_Plus[partIndex][centIndex] = NULL;
      //m_effHisto2D_Plus[partIndex][centIndex] = NULL;
      //m_pTLossHisto3D_Plus[partIndex][centIndex] = NULL;
      //m_embTrackHisto2D_Minus[partIndex][centIndex] = NULL; // particle, centrality
      //m_matchTrackHisto2D_Minus[partIndex][centIndex] = NULL;
      //m_effHisto2D_Minus[partIndex][centIndex] = NULL;
      //m_pTLossHisto3D_Minus[partIndex][centIndex] = NULL;

      m_BTOFEff_Data_Plus[partIndex][centIndex] = NULL;
      m_BTOFEff_Data_Minus[partIndex][centIndex] = NULL;
      m_BTOFEff_Data_Combined[partIndex][centIndex] = NULL;

      m_BTOFEff_Emb_Plus[partIndex][centIndex] = NULL;
      m_BTOFEff_Emb_Minus[partIndex][centIndex] = NULL;
      m_BTOFEff_Emb_Combined[partIndex][centIndex] = NULL;

      m_tpcTrack_For_BTOFEff_Emb_Plus[partIndex][centIndex] = NULL;
      m_bTOFTrack_For_BTOFEff_Emb_Plus[partIndex][centIndex] = NULL;
      m_tpcTrack_For_BTOFEff_Emb_Minus[partIndex][centIndex] = NULL;
      m_bTOFTrack_For_BTOFEff_Emb_Minus[partIndex][centIndex] = NULL;

      m_tpcTrack_For_BTOFEff_Data_Plus[partIndex][centIndex] = NULL;
      m_bTOFTrack_For_BTOFEff_Data_Plus[partIndex][centIndex] = NULL;
      m_tpcTrack_For_BTOFEff_Data_Minus[partIndex][centIndex] = NULL;
      m_bTOFTrack_For_BTOFEff_Data_Minus[partIndex][centIndex] = NULL;


      for(unsigned int fitIndex = 0; fitIndex < 8; fitIndex++){
        for(unsigned int paramIndex = 0; paramIndex < 6; paramIndex++){
          m_BTOFEff_Data_FitParams_Plus[partIndex][centIndex][fitIndex][paramIndex] = NULL; //part cent fit param
          m_BTOFEff_Data_FitParams_Minus[partIndex][centIndex][fitIndex][paramIndex] = NULL;
          m_BTOFEff_Data_FitParams_Combined[partIndex][centIndex][fitIndex][paramIndex] = NULL;
          m_BTOFEff_Emb_FitParams_Plus[partIndex][centIndex][fitIndex][paramIndex] = NULL; //part cent fit param
          m_BTOFEff_Emb_FitParams_Minus[partIndex][centIndex][fitIndex][paramIndex] = NULL;
          m_BTOFEff_Emb_FitParams_Combined[partIndex][centIndex][fitIndex][paramIndex] = NULL;
        }
      }


    }

    for(unsigned int fitIndex = 0; fitIndex < 8; fitIndex++){
      for(unsigned int paramIndex = 0; paramIndex < 6; paramIndex++){
        m_BTOFEff_Data_FitParams_AllCent_Plus[partIndex][fitIndex][paramIndex] = NULL;
        m_BTOFEff_Data_FitParams_AllCent_Minus[partIndex][fitIndex][paramIndex] = NULL;
        m_BTOFEff_Data_FitParams_AllCent_Combined[partIndex][fitIndex][paramIndex] = NULL;
        m_BTOFEff_Emb_FitParams_AllCent_Plus[partIndex][fitIndex][paramIndex] = NULL;
        m_BTOFEff_Emb_FitParams_AllCent_Minus[partIndex][fitIndex][paramIndex] = NULL;
        m_BTOFEff_Emb_FitParams_AllCent_Combined[partIndex][fitIndex][paramIndex] = NULL;
      }
    }

  }



  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    for(int partIndex = 0; partIndex < 9; partIndex++){
      m_particleLoaded_feeddown[partIndex][pmIndex] = false;
      m_particleLoaded_embedding[partIndex][pmIndex] = false;
      m_particleLoaded_btofEff[partIndex][pmIndex] = false;

      m_TPC_Efficiency_ChiSqrd[partIndex][pmIndex] = NULL;
      m_TPC_Efficiency_ChiSqrd_Sys[partIndex][pmIndex] = NULL;
      m_EnergyLoss_ChiSqrd[partIndex][pmIndex] = NULL;
      m_EnergyLoss_ChiSqrd_Sys[partIndex][pmIndex] = NULL;
      m_FeedDown_ChiSqrd[partIndex][pmIndex] = NULL;
      m_FeedDown_ChiSqrd_Sys[partIndex][pmIndex] = NULL;

      //m_FeedDown_Hybrid_ChiSqrd[partIndex][pmIndex] = NULL; // part charge {+1,-1} TH2D is x=rapidity y= centrality 17th bin is AllCent
      //m_FeedDown_Hybrid_ChiSqrd_Sys[partIndex][pmIndex] = NULL; // part charge {+1,-1}

      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        m_FeedDown_primaryOverTotal_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_FeedDown_feeddownOverTotal_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_FeedDown_Fits_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_FeedDown_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;
        m_FeedDown_CovMat_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_FeedDown_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;

        //m_FeedDown_Hybrid_AllCent[partIndex][rapIndex][pmIndex] = NULL; //part rapidity charge {+1,-1}
        //m_FeedDown_Hybrid_Fits_AllCent[partIndex][rapIndex][pmIndex] = NULL; //part rapidity charge {+1,-1}
        //m_FeedDown_Hybrid_CovMat_AllCent[partIndex][rapIndex][pmIndex] = NULL; // part cent rap


        m_TPC_Efficiency_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Efficiency_AllCent_EtaPt[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_Fits_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_Fits_HighPt_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_FitErrGraphs_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_FitErrGraphs_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_CovMat_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_CovMat_HighPt_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_TPC_Eff_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;

        m_EnergyLoss_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_EnergyLossProfile_AllCent[partIndex][rapIndex][pmIndex] = NULL;    //part rapidity charge {+1,-1}
        m_EnergyLoss_Fits_AllCent[partIndex][rapIndex][pmIndex] = NULL; // particle rapdiity
        m_EnergyLoss_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL; // particle rapdiity
        m_EnergyLoss_FitErrGraphs_AllCent[partIndex][rapIndex][pmIndex] = NULL; // part cent rap
        m_EnergyLoss_FitErrGraphs_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL; // part cent rap
        m_EnergyLoss_CovMat_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_EnergyLoss_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;

        m_BTOF_Efficiency_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_Eff_Fits_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_Eff_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_Eff_CovMat_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_Eff_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;

        m_BTOF_fromData_Efficiency_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_fromData_Eff_Fits_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_fromData_Eff_Fits_HighPt_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_fromData_Eff_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_fromData_Eff_CovMat_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_fromData_Eff_CovMat_HighPt_AllCent[partIndex][rapIndex][pmIndex] = NULL;
        m_BTOF_fromData_Eff_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex] = NULL;

        for(int centIndex = 0; centIndex < 16; centIndex++){
          m_FeedDown_primaryOverTotal_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_FeedDown_feeddownOverTotal_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_FeedDown_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_FeedDown_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_FeedDown_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_FeedDown_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;

          m_FeedDown_Hybrid_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL; //part cent rapidity charge {+1,-1}
          m_FeedDown_Hybrid_ByCent_graphs[partIndex][centIndex][rapIndex][pmIndex] = NULL; //part cent rapidity charge {+1,-1}
          m_FeedDown_Hybrid_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL; //part cent rapidity charge {+1,-1}
          m_FeedDown_Hybrid_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL; // part cent rap
          m_FeedDown_Hybrid_Fits_ByCent_lowPt[partIndex][centIndex][rapIndex][pmIndex] = NULL; //part cent rapidity charge {+1,-1}
          m_FeedDown_Hybrid_CovMat_ByCent_lowPt[partIndex][centIndex][rapIndex][pmIndex] = NULL; // part cent rap

          m_TPC_Efficiency_ByCent[partIndex][centIndex][rapIndex][pmIndex]       = NULL;
          m_TPC_Efficiency_ByCent_EtaPt[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_Fits_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_FitErrGraphs_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_FitErrGraphs_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_CovMat_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_TPC_Eff_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;

          m_EnergyLoss_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_EnergyLossProfile_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL; //part cent rapidity charge {+1,-1}

          m_EnergyLoss_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL; // particle rapdiity
          m_EnergyLoss_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL; // particle rapdiity
          m_EnergyLoss_FitErrGraphs_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL; // part cent rap
          m_EnergyLoss_FitErrGraphs_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL; // part cent rap
          m_EnergyLoss_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_EnergyLoss_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;

          m_BTOF_Efficiency_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_Eff_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_Eff_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_Eff_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_Eff_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;

          m_BTOF_fromData_Efficiency_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_fromData_Eff_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_fromData_Eff_Fits_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_fromData_Eff_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_fromData_Eff_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_fromData_Eff_CovMat_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex] = NULL;
          m_BTOF_fromData_Eff_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex] = NULL;
        }
      }
    }
  }



  //for(int partIndex = 0; partIndex <9; partIndex++) m_BTOF_highMtM0ToLine[partIndex] = NULL;
  //m_BTOF_highMtM0ToLine[0] = PhysMath::linearFunct("pionBTOFLine",0,1.3,1.0,1.7);
  //m_BTOF_highMtM0ToLine[1] = PhysMath::linearFunct("kaonBTOFLine",0,1.3,1.0,1.7);
  //m_BTOF_highMtM0ToLine[2] = PhysMath::linearFunct("protonBTOFLine",0,1.3,1.0,1.7);
  //m_BTOF_highMtM0ToLine[4] = PhysMath::linearFunct("deuteronBTOFLine",0,1.3,1.0,1.7);
  //m_BTOF_highMtM0ToLine[5] = PhysMath::linearFunct("tritonBTOFLine",0,1.3,1.0,1.7);
  //
  //for(int partIndex = 0; partIndex <9; partIndex++) m_TPC_Eff_highMtM0ToLine[partIndex] = NULL;
  //m_TPC_Eff_highMtM0ToLine[0] = PhysMath::linearFunct("pionTPC_EffLine",0,1.3,1.0,1.7);
  //m_TPC_Eff_highMtM0ToLine[1] = PhysMath::linearFunct("kaonTPC_EffLine",0,1.3,1.0,1.7);
  //m_TPC_Eff_highMtM0ToLine[2] = PhysMath::linearFunct("protonTPC_EffLine",0,1.3,1.0,1.7);
  //m_TPC_Eff_highMtM0ToLine[4] = PhysMath::linearFunct("deuteronTPC_EffLine",0,1.3,1.0,1.7);
  //m_TPC_Eff_highMtM0ToLine[5] = PhysMath::linearFunct("tritonTPC_EffLine",0,1.3,1.0,1.7);

  for(int partIndex = 0; partIndex <9; partIndex++) m_BTOF_highMtM0ToLine[partIndex] = NULL;
  m_BTOF_highMtM0ToLine[0] = PhysMath::parabolaFunct("pionBTOFParab",-1,1.8,0,1.8,1,1.8);
  m_BTOF_highMtM0ToLine[1] = PhysMath::parabolaFunct("kaonBTOFParab",-1,2.8,0,1.5,1,2.8);
  m_BTOF_highMtM0ToLine[2] = PhysMath::parabolaFunct("protonBTOFParab",-1,2.8,0,1.5,1,2.8);
  m_BTOF_highMtM0ToLine[4] = PhysMath::parabolaFunct("deuteronBTOFParab",-1,2.8,0,1.5,1,2.8);
  m_BTOF_highMtM0ToLine[5] = PhysMath::parabolaFunct("tritonBTOFParab",-1,2.8,0,1.5,1,2.8);

  for(int partIndex = 0; partIndex <9; partIndex++) m_TPC_Eff_highMtM0ToLine[partIndex] = NULL;
  m_TPC_Eff_highMtM0ToLine[0] = PhysMath::parabolaFunct("pionTPC_EffParab",-1,1,0,0.6,1,1);
  m_TPC_Eff_highMtM0ToLine[1] = PhysMath::parabolaFunct("kaonTPC_EffParab",-1,1.6,0,2.3,1,1.6);
  m_TPC_Eff_highMtM0ToLine[2] = PhysMath::parabolaFunct("protonTPC_EffParab",-1,1.7,0,1,1,1.7);
  m_TPC_Eff_highMtM0ToLine[4] = PhysMath::parabolaFunct("deuteronTPC_EffParab",-1,1.7,0,1,1,1.7);
  m_TPC_Eff_highMtM0ToLine[5] = PhysMath::parabolaFunct("tritonTPC_EffParab",-1,1.7,0,1,1,1.7);

  //PROTON DCA INFORMATION
  m_proton_knockout_ChiSqrd = NULL;
  for(int centIndex = 0; centIndex < 16; centIndex++){
    m_proton_dca_ChiSqrd[centIndex] = NULL;
    m_proton_knockout_histogram_ByCent[centIndex] = NULL;

  }
  for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
    m_proton_knockout_AllCent[rapIndex] = NULL;
    m_proton_knockout_Fits_AllCent[rapIndex] = NULL;
    for(int mTm0Index = 0; mTm0Index < g_EffFit_MaxNumMtM0Bins; mTm0Index++){
      m_proton_dca_AllCent[rapIndex][mTm0Index][0] = NULL;
      m_proton_dca_AllCent[rapIndex][mTm0Index][1] = NULL;
      m_proton_dca_Fits_AllCent[rapIndex][mTm0Index] = NULL;
    }
    for(int centIndex = 0; centIndex < 16; centIndex++){
      m_proton_knockout_ByCent[centIndex][rapIndex] = NULL;
      m_proton_knockout_Fits_ByCent[centIndex][rapIndex] = NULL;
      m_proton_knockout_CovMat_ByCent[centIndex][rapIndex] = NULL;
      for(int mTm0Index = 0; mTm0Index < g_EffFit_MaxNumMtM0Bins; mTm0Index++){
        m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][0] = NULL;
        m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][1] = NULL;
        m_proton_dca_Fits_ByCent[centIndex][rapIndex][mTm0Index] = NULL;
      }
    }
  }

  m_spline  = NULL;
  m_spline5 = NULL;

}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________



void EfficiencyFitter::setLegendreOrderForTPCEff(int a_partIndex, int a_charge, int a_order){
  if(a_charge < 0) m_polynomialOrderForTPCEff[a_partIndex][1] = a_order;
  else m_polynomialOrderForTPCEff[a_partIndex][0] = a_order;
  if(a_order > 29) cout << " ERROR: EfficiencyFitter::setLegendreOrderForTPCEff, you can't select that high of an order for the polynomial (>29)" << endl;
}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void EfficiencyFitter::setLegendreOrderForBTOFEff(int a_partIndex, int a_charge, int a_order){
  if(a_charge < 0) m_polynomialOrderForBTOFEff[a_partIndex][1] = a_order;
  else m_polynomialOrderForBTOFEff[a_partIndex][0] = a_order;
  if(a_order > 29) cout << " ERROR: EfficiencyFitter::setLegendreOrderForBTOFEff, you can't select that high of an order for the polynomial (>29)" << endl;

}




//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________



void EfficiencyFitter::loadEmbeddingFiles(string a_inFileName, string a_inFileNameTwo, int a_partIndex, int a_charge, int a_mtm0Rebin, bool a_convert16to9CentBins, bool a_savePtLossGausImages, bool a_fitEnergyLossByCent, int a_energyLossMtM0Rebin, bool a_doEnergyLossFits, vector<double> a_refMultCuts){
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << " call to EfficiencyFitter::loadEmbeddingFiles(string " << a_inFileName << ",string" << a_inFileNameTwo << ",int " << a_partIndex << ", int "
         << a_charge << ", int " << a_mtm0Rebin << ", bool " << a_convert16to9CentBins << ", bool "
         << a_savePtLossGausImages << " , bool " << a_fitEnergyLossByCent << " , int " << a_energyLossMtM0Rebin
         << ", bool " << a_doEnergyLossFits << " , vector<double> " << " cuts... " << ")" <<endl;
  #endif
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;

  // Fitting Canvas for 1D Gaussian Fits of pT loss
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  gSystem->mkdir(Form("%s/%s",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  gSystem->mkdir(Form("%s/%s/EnergyLoss_AllCent",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  gSystem->mkdir(Form("%s/%s/EnergyLoss_ByCent",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  gSystem->mkdir(Form("%s/%s/EnergyLoss_AllCent/mTLossFits",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  cout << " Loading " << m_partInfo->GetParticleName(a_partIndex,a_charge).Data() << endl;

  //Open the Embedding File
  //string fullInFileName = a_inFileName_FirstHalf + m_partInfo->GetParticleName(a_partIndex,a_charge).Data() + a_inFileName_SecondHalf;
  TFile* inFile = new TFile(a_inFileName.c_str(),"READ");
  if(!inFile->IsOpen()){
  	cout << "ERROR: Embedding File is not open!!  " << a_inFileName << endl;
  }
  TFile* inFileTwo = NULL;
  bool secondFileLoaded = false;
  if(a_inFileNameTwo != ""){
    cout << "Loading negative charge file to add to plus: " << a_inFileNameTwo << endl;
    inFileTwo = new TFile(a_inFileNameTwo.c_str(),"READ");
    if(inFileTwo->IsOpen()){
      secondFileLoaded = true;
    }else{
      cout << "ERROR: second file failed to open... " << a_inFileNameTwo << endl;
      inFileTwo = NULL;
    }
  }
  inFile->cd();

  TH3I* embTrackHisto_byRefMult   = (TH3I*) inFile->Get("EfficiencyHistograms/embTrackHisto_byRefMult")->Clone();
  TH3I* matchTrackHisto_byRefMult = (TH3I*) inFile->Get("EfficiencyHistograms/matchTrackHisto_byRefMult")->Clone();
  bool loaded3DRefMult = (embTrackHisto_byRefMult && matchTrackHisto_byRefMult);
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << " 3D PTR CHECK: " << embTrackHisto_byRefMult << " " << matchTrackHisto_byRefMult << "   loaded bool: " << loaded3DRefMult << endl;
  #endif

  if(secondFileLoaded){
    inFileTwo->cd();
    TH3I* embTrackHisto_byRefMult_Two   = (TH3I*) inFileTwo->Get("EfficiencyHistograms/embTrackHisto_byRefMult")->Clone();
    TH3I* matchTrackHisto_byRefMult_Two = (TH3I*) inFileTwo->Get("EfficiencyHistograms/matchTrackHisto_byRefMult")->Clone();
    embTrackHisto_byRefMult->Add(embTrackHisto_byRefMult_Two);
    matchTrackHisto_byRefMult->Add(matchTrackHisto_byRefMult_Two);
    delete embTrackHisto_byRefMult_Two;
    delete matchTrackHisto_byRefMult_Two;
  }



  TH2D* embTracksAllCent         = NULL;
  TH2D* matchTracksAllCent       = NULL;
  TH2D* embTracksAllCent_EtaPt   = NULL;
  TH2D* matchTracksAllCent_EtaPt = NULL;
  TH3D* energyLossAllCent        = NULL;
  TH2D* tpc_forBTOFAllCent       = NULL;
  TH2D* btof_forBTOFAllCent      = NULL;

  TH2D* embeddedTrackHistogramRebin[16]; // cent
  TH2D* matchTrackHistogramRebin[16]; // cent
  TH2D* embeddedTrackHistogramRebin_EtaPt[16]; // cent
  TH2D* matchTrackHistogramRebin_EtaPt[16]; // cent
  TH3D* energyLossHistogramRebin[16]; // cent
  TH2D* tpc_forBTOFHistogramRebin[16];
  TH2D* btof_forBTOFHistogramRebin[16];



  vector< TH2I* > embTrackHisto_byRefMult_Vector;
  vector< TH2I* > matchTrackHisto_byRefMult_Vector;
  if(loaded3DRefMult){
    embTrackHisto_byRefMult_Vector   = HistogramUtilities::projeciton2DSlices_Z(embTrackHisto_byRefMult,a_refMultCuts);
    matchTrackHisto_byRefMult_Vector = HistogramUtilities::projeciton2DSlices_Z(matchTrackHisto_byRefMult,a_refMultCuts);
  }

  m_mTm0ReponseMatrixHisto[a_partIndex][pmIndex] = (TH3I*) inFile->Get("mTm0ResponseMatrix");
  if(m_mTm0ReponseMatrixHisto[a_partIndex][pmIndex]){
    m_mTm0ReponseMatrixHisto[a_partIndex][pmIndex]->SetName(Form("mTm0ResponseMatrix_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex] = HistogramUtilities::convertInt_to_DoubleHisto3D(m_mTm0ReponseMatrixHisto[a_partIndex][pmIndex]);
    m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->SetName(Form("mTm0ResponseMatrix_%s_PDF",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  }

  //m_mTm0ReponseMatrixHisto_BTOF[a_partIndex][pmIndex] = (TH3I*) inFile->Get("mTm0ResponseMatrixBTOF");
  //m_mTm0ReponseMatrixHisto_BTOF[a_partIndex][pmIndex]->SetName(Form("mTm0ResponseMatrixBTOF_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  //m_mTm0ReponseMatrixHisto_BTOF_PDF[a_partIndex][pmIndex] = HistogramUtilities::convertInt_to_DoubleHisto3D(m_mTm0ReponseMatrixHisto_BTOF[a_partIndex][pmIndex]);
  //m_mTm0ReponseMatrixHisto_BTOF_PDF[a_partIndex][pmIndex]->SetName(Form("mTm0ResponseMatrixBTOF_%s_PDF",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));

  //--- Centrality Index Loop  ------
  int numLoadedCentBins = 0;
  for(int centIndex = 0; centIndex < 16; centIndex++){
    inFile->cd();
    TH2I* test = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex));
    if(!test){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Stopping Loading at Cent " << centIndex << endl;
      #endif
      break;
    }

    numLoadedCentBins++;

    TH2I* embeddedTrackHistogram_Int     = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex))->Clone();
    TH2I* matchTrackHistogram_Int        = (TH2I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_Cent%d",centIndex))->Clone();
    TH2I* embeddedTrackHistogram_EtaPt_Int     = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_pTEta_Cent%d",centIndex))->Clone();
    TH2I* matchTrackHistogram_EtaPt_Int        = (TH2I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_pTEta_Cent%d",centIndex))->Clone();
    TH3I* mTLossHistogram_Int            = NULL;
    if(m_pTLossMode) mTLossHistogram_Int = (TH3I*) inFile->Get(Form("EnergyLossHistograms/pTLossHisto_Cent%d",centIndex))->Clone();
    else             mTLossHistogram_Int = (TH3I*) inFile->Get(Form("EnergyLossHistograms/mTm0LossHisto_Cent%d",centIndex))->Clone();
    TH2I* tpc_forBTOF_Histogram_Int      = (TH2I*) inFile->Get(Form("EfficiencyHistograms/tpc_tracks_reco_Cent%d",centIndex))->Clone();
    TH2I* btof_forBTOF_Histogram_Int     = (TH2I*) inFile->Get(Form("EfficiencyHistograms/btof_tracks_reco_Cent%d",centIndex))->Clone();

    if(secondFileLoaded){
      inFileTwo->cd();
      TH2I* embeddedTrackHistogram_Int_Two     = (TH2I*) inFileTwo->Get(Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex))->Clone();
      TH2I* matchTrackHistogram_Int_Two        = (TH2I*) inFileTwo->Get(Form("EfficiencyHistograms/matchTrackHisto_Cent%d",centIndex))->Clone();
      TH2I* embeddedTrackHistogram_EtaPt_Int_Two     = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_pTEta_Cent%d",centIndex))->Clone();
      TH2I* matchTrackHistogram_EtaPt_Int_Two        = (TH2I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_pTEta_Cent%d",centIndex))->Clone();
      TH3I* mTLossHistogram_Int_Two            = NULL;
      if(m_pTLossMode) mTLossHistogram_Int_Two = (TH3I*) inFileTwo->Get(Form("EnergyLossHistograms/pTLossHisto_Cent%d",centIndex))->Clone();
      else             mTLossHistogram_Int_Two = (TH3I*) inFileTwo->Get(Form("EnergyLossHistograms/mTm0LossHisto_Cent%d",centIndex))->Clone();
      TH2I* tpc_forBTOF_Histogram_Int_Two      = (TH2I*) inFileTwo->Get(Form("EfficiencyHistograms/tpc_tracks_reco_Cent%d",centIndex))->Clone();
      TH2I* btof_forBTOF_Histogram_Int_Two     = (TH2I*) inFileTwo->Get(Form("EfficiencyHistograms/btof_tracks_reco_Cent%d",centIndex))->Clone();

      embeddedTrackHistogram_Int->Add(embeddedTrackHistogram_Int_Two);
      matchTrackHistogram_Int->Add(matchTrackHistogram_Int_Two);
      embeddedTrackHistogram_EtaPt_Int->Add(embeddedTrackHistogram_EtaPt_Int_Two);
      matchTrackHistogram_EtaPt_Int->Add(matchTrackHistogram_EtaPt_Int_Two);
      mTLossHistogram_Int->Add(mTLossHistogram_Int_Two);
      tpc_forBTOF_Histogram_Int->Add(tpc_forBTOF_Histogram_Int_Two);
      btof_forBTOF_Histogram_Int->Add(btof_forBTOF_Histogram_Int_Two);

      delete embeddedTrackHistogram_Int_Two;
      delete matchTrackHistogram_Int_Two;
      delete embeddedTrackHistogram_EtaPt_Int_Two;
      delete matchTrackHistogram_EtaPt_Int_Two;
      delete mTLossHistogram_Int_Two;
      delete tpc_forBTOF_Histogram_Int_Two;
      delete btof_forBTOF_Histogram_Int_Two;
    }

    //##########    STORE THE RAPIDITY STRUCTURE #####################
    if(!m_embedding_rapidity_structure[a_partIndex]){
      vector<double> rap_bin_edges = HistogramUtilities::getBinEdges_X(embeddedTrackHistogram_Int);
      m_embedding_rapidity_structure[a_partIndex] = new TH1D(Form("embedding_%s_rap_struct",m_partInfo->GetParticleName(a_partIndex,0).Data()),
        Form("Efficiency Rapidity Structure from Embedding for %s; y_{%s}",
              m_partInfo->GetParticleSymbol(a_partIndex,0).Data(),m_partInfo->GetParticleSymbol(a_partIndex,0).Data()),
              (int)(rap_bin_edges.size() - 1),&rap_bin_edges[0]);
    }

    if(! m_embedding_pseudorapidity_structure[a_partIndex] ){
      vector<double> eta_bin_edges = HistogramUtilities::getBinEdges_X(embeddedTrackHistogram_EtaPt_Int);
      m_embedding_pseudorapidity_structure[a_partIndex] = new TH1D(Form("embedding_%s_eta_struct",m_partInfo->GetParticleName(a_partIndex,0).Data()),
        Form("Efficiency Pseudorapidity Structure from Embedding for %s; #eta_{%s}",
              m_partInfo->GetParticleSymbol(a_partIndex,0).Data(),m_partInfo->GetParticleSymbol(a_partIndex,0).Data()),
              (int)(eta_bin_edges.size() - 1),&eta_bin_edges[0]);
    }


    //This overwrites the efficiency histograms with refMult integrated histograms when embedding for a particular energy is not available
    if(loaded3DRefMult && a_refMultCuts.size() > 0 && centIndex  < embTrackHisto_byRefMult_Vector.size() && centIndex < matchTrackHisto_byRefMult_Vector.size()){
      embeddedTrackHistogram_Int = embTrackHisto_byRefMult_Vector.at(centIndex);    // modified from centIndex + 1 to just centIndex
      matchTrackHistogram_Int    = matchTrackHisto_byRefMult_Vector.at(centIndex);  //
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  Using refMult cuts to make emb/match histograms, instead of the main _CentXX histograms... (not an error, just make sure this is what you want)" << endl;
        cout << "     For CENT " << centIndex << " : " << endl;
        cout << "            " << embeddedTrackHistogram_Int->GetTitle() << endl;
        cout << "            " << matchTrackHistogram_Int->GetTitle() << endl;
      #endif
    }

    if(!embeddedTrackHistogram_Int){
        if(centIndex == 0){
          cout << " ERROR: Something is not loading right... embeddedTrackHistogram is NULL for Cent00" << endl;
        }
        continue;
    }

    TH2D* embeddedTrackHistogram       = HistogramUtilities::convertInt_to_DoubleHisto2D(embeddedTrackHistogram_Int);
    TH2D* matchTrackHistogram          = HistogramUtilities::convertInt_to_DoubleHisto2D(matchTrackHistogram_Int);
    TH2D* embeddedTrackHistogram_EtaPt = HistogramUtilities::convertInt_to_DoubleHisto2D(embeddedTrackHistogram_EtaPt_Int);
    TH2D* matchTrackHistogram_EtaPt    = HistogramUtilities::convertInt_to_DoubleHisto2D(matchTrackHistogram_EtaPt_Int);
    TH3D* mTLossHistogram              = HistogramUtilities::convertInt_to_DoubleHisto3D(mTLossHistogram_Int);
    TH2D* tpc_forBTOF_Histogram        = HistogramUtilities::convertInt_to_DoubleHisto2D(tpc_forBTOF_Histogram_Int);
    TH2D* btof_forBTOF_Histogram       = HistogramUtilities::convertInt_to_DoubleHisto2D(btof_forBTOF_Histogram_Int);

    gSystem->mkdir(Form("%s/%s/EnergyLoss_ByCent/Cent%02d",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex));
    gSystem->mkdir(Form("%s/%s/EnergyLoss_ByCent/Cent%02d/mTLossFits",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex));

    // REBIN THE HISTOGRAMS IN mT-m0
    embeddedTrackHistogramRebin[centIndex]       = HistogramUtilities::reBinVar2D(embeddedTrackHistogram,1,a_mtm0Rebin);
    matchTrackHistogramRebin[centIndex]          = HistogramUtilities::reBinVar2D(matchTrackHistogram,1,a_mtm0Rebin);
    embeddedTrackHistogramRebin_EtaPt[centIndex] = HistogramUtilities::reBinVar2D(embeddedTrackHistogram_EtaPt,1,a_mtm0Rebin);
    matchTrackHistogramRebin_EtaPt[centIndex]    = HistogramUtilities::reBinVar2D(matchTrackHistogram_EtaPt,1,a_mtm0Rebin);
    energyLossHistogramRebin[centIndex]          = HistogramUtilities::reBinVar3D(mTLossHistogram,1,a_energyLossMtM0Rebin,1);
    tpc_forBTOFHistogramRebin[centIndex]         = HistogramUtilities::reBinVar2D(tpc_forBTOF_Histogram,1,a_mtm0Rebin);
    btof_forBTOFHistogramRebin[centIndex]        = HistogramUtilities::reBinVar2D(btof_forBTOF_Histogram,1,a_mtm0Rebin);
    HistogramUtilities::makeSqrtNErrors2D(embeddedTrackHistogramRebin[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(matchTrackHistogramRebin[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(embeddedTrackHistogramRebin_EtaPt[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(matchTrackHistogramRebin_EtaPt[centIndex]);
    HistogramUtilities::makeSqrtNErrors3D(energyLossHistogramRebin[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(tpc_forBTOFHistogramRebin[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(btof_forBTOFHistogramRebin[centIndex]);

    if(a_mtm0Rebin != 1){
      delete embeddedTrackHistogram;
      delete matchTrackHistogram;
      delete embeddedTrackHistogram_EtaPt;
      delete matchTrackHistogram_EtaPt;
      delete mTLossHistogram;
      delete tpc_forBTOF_Histogram;
      delete btof_forBTOF_Histogram;
    }
    delete embeddedTrackHistogram_Int;
    delete matchTrackHistogram_Int;
    delete embeddedTrackHistogram_EtaPt_Int;
    delete matchTrackHistogram_EtaPt_Int;
    delete mTLossHistogram_Int;
    delete tpc_forBTOF_Histogram_Int;
    delete btof_forBTOF_Histogram_Int;


    //############      CREATE CENTRALITY INTEGRATED  HISTOS ##################
    if(centIndex == 0){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  Making centrality integrated histograms..." << endl;
      #endif
      embTracksAllCent = (TH2D*) embeddedTrackHistogramRebin[0]->Clone();
      embTracksAllCent->SetName(Form("EmbbedTrack_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      matchTracksAllCent = (TH2D*) matchTrackHistogramRebin[0]->Clone();
      matchTracksAllCent->SetName(Form("MatchTrack_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      embTracksAllCent_EtaPt = (TH2D*) embeddedTrackHistogramRebin_EtaPt[0]->Clone();
      embTracksAllCent_EtaPt->SetName(Form("EmbbedTrack_AllCent_%s_EtaPt",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      matchTracksAllCent_EtaPt = (TH2D*) matchTrackHistogramRebin_EtaPt[0]->Clone();
      matchTracksAllCent_EtaPt->SetName(Form("MatchTrack_AllCent_%s_EtaPt",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      energyLossAllCent = (TH3D*) energyLossHistogramRebin[0]->Clone();
      energyLossAllCent->SetName(Form("mTLoss_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      tpc_forBTOFAllCent = (TH2D*) tpc_forBTOFHistogramRebin[0]->Clone();
      tpc_forBTOFAllCent->SetName(Form("tpc_forBTOF_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      btof_forBTOFAllCent = (TH2D*) btof_forBTOFHistogramRebin[0]->Clone();
      btof_forBTOFAllCent->SetName(Form("btof_forBTOF_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Name: " << embTracksAllCent->GetName()    << " Ptr: "  << embTracksAllCent << endl;
        cout << " Name: " << matchTracksAllCent->GetName()  << " Ptr: "  << matchTracksAllCent << endl;
        cout << " Name: " << embTracksAllCent_EtaPt->GetName()    << " Ptr: "  << embTracksAllCent_EtaPt << endl;
        cout << " Name: " << matchTracksAllCent_EtaPt->GetName()  << " Ptr: "  << matchTracksAllCent_EtaPt << endl;
        cout << " Name: " << energyLossAllCent->GetName()   << " Ptr: "  << energyLossAllCent << endl;
      #endif

    }else{
      embTracksAllCent->Add(embeddedTrackHistogramRebin[centIndex]);
      matchTracksAllCent->Add(matchTrackHistogramRebin[centIndex]);
      embTracksAllCent_EtaPt->Add(embeddedTrackHistogramRebin_EtaPt[centIndex]);
      matchTracksAllCent_EtaPt->Add(matchTrackHistogramRebin_EtaPt[centIndex]);
      energyLossAllCent->Add(energyLossHistogramRebin[centIndex]);
      tpc_forBTOFAllCent->Add(tpc_forBTOFHistogramRebin[centIndex]);
      btof_forBTOFAllCent->Add(btof_forBTOFHistogramRebin[centIndex]);
    }


    //#### MAKE CHI^2 HISTOGRAMS ####
    if(!m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]){
      m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex] = new TH2D(
        Form("energyLoss_%s_ChiSqrd",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()),
        Form("#Chi^{2}/dof for Nominal Energy Loss Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()),
        embTracksAllCent->GetNbinsX(),embTracksAllCent->GetXaxis()->GetBinLowEdge(1),embTracksAllCent->GetXaxis()->GetBinLowEdge(embTracksAllCent->GetNbinsX()+1),
        17,-0.5,16.5);
      for(int iii = 0; iii < 16; iii++){
        m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(iii+1, Form("Cent %02d",iii));
      }
      m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(17, "All Cent");
      m_EnergyLoss_ChiSqrd_Sys[a_partIndex][pmIndex] = (TH2D*) m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->Clone();
      m_EnergyLoss_ChiSqrd_Sys[a_partIndex][pmIndex]->SetName(Form("energyLoss_%s_ChiSqrd_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      m_EnergyLoss_ChiSqrd_Sys[a_partIndex][pmIndex]->SetTitle(Form("#Chi^{2}/dof for Systematic Energy Loss Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()));
    }

    //#### MAKE CHI^2 HISTOGRAMS ####
    if(!m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]){
      m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex] = new TH2D(
        Form("tpcEff_%s_ChiSqrd",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()),
        Form("#Chi^{2}/dof for Nominal TPC Efficiency Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()),
        embTracksAllCent->GetNbinsX(),embTracksAllCent->GetXaxis()->GetBinLowEdge(1),embTracksAllCent->GetXaxis()->GetBinLowEdge(embTracksAllCent->GetNbinsX()+1),
        17,-0.5,16.5);
      for(int iii = 0; iii < numLoadedCentBins; iii++){
        m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(iii+1, Form("Cent %02d",iii));
      }
      m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(17, "All Cent");
      m_TPC_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex] = (TH2D*) m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]->Clone();
      m_TPC_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetName(Form("tpcEff_%s_ChiSqrd_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      m_TPC_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetTitle(Form("#Chi^{2}/dof for Systematic TPC Efficiency Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()));
    }

    //#### MAKE CHI^2 HISTOGRAMS ####
    if(!m_BTOF_Efficiency_ChiSqrd[a_partIndex][pmIndex]){
      m_BTOF_Efficiency_ChiSqrd[a_partIndex][pmIndex] = new TH2D(
        Form("BTOFEff_%s_ChiSqrd",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()),
        Form("#Chi^{2}/dof for Nominal BTOF Efficiency Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()),
        btof_forBTOFAllCent->GetNbinsX(),btof_forBTOFAllCent->GetXaxis()->GetBinLowEdge(1),btof_forBTOFAllCent->GetXaxis()->GetBinLowEdge(btof_forBTOFAllCent->GetNbinsX()+1),
        17,-0.5,16.5);
      for(int iii = 0; iii < numLoadedCentBins; iii++){
        m_BTOF_Efficiency_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(iii+1, Form("Cent %02d",iii));
      }
      m_BTOF_Efficiency_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(17, "All Cent");
      m_BTOF_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex] = (TH2D*) m_BTOF_Efficiency_ChiSqrd[a_partIndex][pmIndex]->Clone();
      m_BTOF_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetName(Form("BTOFEff_%s_ChiSqrd_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      m_BTOF_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetTitle(Form("#Chi^{2}/dof for Systematic BTOF Efficiency Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()));
    }
  }   //cent Index




  //##########     COMBINE 16 to 9 Centrality Bins  ################
  if(a_convert16to9CentBins){
    cout << "Converting 16 to 9 Cent Bins..." << endl;
    if(!embeddedTrackHistogramRebin[15]){
      cout << "  ERROR: You cant use the a_convert16to9CentBins unless you actually have 16 bins in the feed down file..." << endl;
      return;
    }

    for(int centIndex = 2; centIndex < 16; centIndex+=2){
      int newCentIndex = 1+centIndex/2;
      embeddedTrackHistogramRebin[newCentIndex]->Add(embeddedTrackHistogramRebin[centIndex],embeddedTrackHistogramRebin[centIndex+1]);
      matchTrackHistogramRebin[newCentIndex]->Add(matchTrackHistogramRebin[centIndex],matchTrackHistogramRebin[centIndex+1]);
      embeddedTrackHistogramRebin_EtaPt[newCentIndex]->Add(embeddedTrackHistogramRebin_EtaPt[centIndex],embeddedTrackHistogramRebin_EtaPt[centIndex+1]);
      matchTrackHistogramRebin_EtaPt[newCentIndex]->Add(matchTrackHistogramRebin_EtaPt[centIndex],matchTrackHistogramRebin_EtaPt[centIndex+1]);
      energyLossHistogramRebin[newCentIndex]->Add(energyLossHistogramRebin[centIndex],energyLossHistogramRebin[centIndex+1]);
      tpc_forBTOFHistogramRebin[newCentIndex]->Add(tpc_forBTOFHistogramRebin[centIndex],tpc_forBTOFHistogramRebin[centIndex+1]);
      btof_forBTOFHistogramRebin[newCentIndex]->Add(btof_forBTOFHistogramRebin[centIndex],btof_forBTOFHistogramRebin[centIndex+1]);
    }

    for(int centIndex = 9; centIndex < 16; centIndex++){
      delete embeddedTrackHistogramRebin[centIndex];
      delete matchTrackHistogramRebin[centIndex];
      delete embeddedTrackHistogramRebin_EtaPt[centIndex];
      delete matchTrackHistogramRebin_EtaPt[centIndex];
      delete energyLossHistogramRebin[centIndex];
      delete tpc_forBTOFHistogramRebin[centIndex];
      delete btof_forBTOFHistogramRebin[centIndex];
    }
  }


  //######################      Make The TPC Efficiency Projections       ########################
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "  Making Projections for TPC Efficiency..." << endl;
  #endif
  for(int rapBin = 1; rapBin <= embTracksAllCent->GetNbinsX(); rapBin++){
    int rapIndex = rapBin - 1;
    if(m_testingMidRapidity && rapBin != embTracksAllCent->GetXaxis()->FindBin(m_ycm)) continue;
    double rapidity = embTracksAllCent->GetXaxis()->GetBinCenter(rapBin);
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - rapidity) < 0.005){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating rapidity = " << rapidity << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }


    TH1D* embTrackHisto   = embTracksAllCent->ProjectionY(Form("proj_emb_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);
    TH1D* matchTrackHisto = matchTracksAllCent->ProjectionY(Form("proj_mat_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);
    TH1D* diff_histo = (TH1D*) embTrackHisto->Clone();
    TH1D* rebinDiffTrack = nullptr;
    for(int bin = 1; bin <= diff_histo->GetNbinsX(); bin++){
      double content = embTrackHisto->GetBinContent(bin) - matchTrackHisto->GetBinContent(bin);
      double error = sqrt(content);
      if(!std::isfinite(error) || error < 0.0) error = 0.0;
      diff_histo->SetBinContent(bin,content);
      diff_histo->SetBinError(bin,error);
    }
    if(m_useChebyshevSpacing){
      int firstBin = HistogramUtilities::firstBinWithContent(matchTrackHisto);
      int lastBin  = HistogramUtilities::lastBinWithContent(matchTrackHisto);
      double lowEdgeFirstBin = matchTrackHisto->GetBinLowEdge(firstBin);
      double highEdgeLastBin = matchTrackHisto->GetBinLowEdge(lastBin+1);
      TH1D* rebinMatchTrack =  HistogramUtilities::rebinHistoToChebyshevSpacing(matchTrackHisto, lowEdgeFirstBin, highEdgeLastBin,3,20);
      TH1D* rebinEmbTrack = (TH1D*) rebinMatchTrack->Clone();
      HistogramUtilities::zeroBins(rebinEmbTrack);
      HistogramUtilities::rebinHistoToAlternateBinStructure(embTrackHisto,rebinEmbTrack);
      delete embTrackHisto;
      delete matchTrackHisto;
      embTrackHisto = rebinEmbTrack;
      matchTrackHisto = rebinMatchTrack;
    }else{
      int firstBin = HistogramUtilities::firstBinWithContent(matchTrackHisto);
      int lastBin  = HistogramUtilities::lastBinWithContent(matchTrackHisto);
      double lowEdgeFirstBin = matchTrackHisto->GetBinLowEdge(firstBin);
      double highEdgeLastBin = matchTrackHisto->GetBinLowEdge(lastBin+1);

      TH1D* truncDiffHisto = HistogramUtilities::truncateHistogram(diff_histo,lowEdgeFirstBin,highEdgeLastBin);
      delete diff_histo;
      diff_histo = truncDiffHisto;
      rebinDiffTrack = HistogramUtilities::rebinEqualByCDF(diff_histo, m_tpc_eff_no_fit_num_bins,false,false); // rebins so more bins where content is changing, less bins where there is less change
      HistogramUtilities::zeroBins(rebinDiffTrack);
      TH1D* rebinMatchTrack = (TH1D*) rebinDiffTrack->Clone();
      TH1D* rebinEmbTrack   = (TH1D*) rebinDiffTrack->Clone();
      HistogramUtilities::rebinHistoToAlternateBinStructure(matchTrackHisto,rebinMatchTrack);
      HistogramUtilities::rebinHistoToAlternateBinStructure(embTrackHisto,rebinEmbTrack);
      delete embTrackHisto;
      delete matchTrackHisto;
      embTrackHisto   = rebinEmbTrack;
      matchTrackHisto = rebinMatchTrack;
      HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
      HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "  ## Rapidity Bin = " << rapBin << endl;
      cout << "  TH1D*'s : " << embTrackHisto  << " " << matchTrackHisto << endl;
      cout << "  NumBins : " << embTrackHisto->GetNbinsX()  << " " << matchTrackHisto->GetNbinsX() << endl;
    #endif

    // create the efficiency graph from the two rebined projection TH1Ds
    m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
    m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
    m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality Integrated y=[%1.2f,%1.2f]; m_{T}-m_{0}; TPC Efficiency",
               m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),embTracksAllCent->GetXaxis()->GetBinLowEdge(rapBin),embTracksAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
    m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->Divide(matchTrackHisto,embTrackHisto,"cl=0.68268949 b(1,1) mode");

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "    TPC Eff All Cent Ptr: " << m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] << "    " <<  m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->GetName() << endl;
    #endif

    delete embTrackHisto;
    delete matchTrackHisto;
    delete diff_histo;

    for(int centIndex = 0; centIndex < numLoadedCentBins; centIndex++){
      if(!embeddedTrackHistogramRebin[centIndex]){
        m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = NULL;
        continue;
      }
      if(m_testingOneCent && centIndex != 0) continue;

      TH1D* emb_histo_cent   = embeddedTrackHistogramRebin[centIndex]->ProjectionY(Form("proj_emb_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);
      TH1D* match_histo_cent = matchTrackHistogramRebin[centIndex]->ProjectionY(Form("proj_match_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);
      TH1D* diff_histo_cent = (TH1D*) emb_histo_cent->Clone();
      for(int bin = 1; bin <= diff_histo_cent->GetNbinsX(); bin++){
        double content = emb_histo_cent->GetBinContent(bin) - match_histo_cent->GetBinContent(bin);
        double error = sqrt(content);
        if(!std::isfinite(error) || error < 0.0) error = 0.0;
        diff_histo_cent->SetBinContent(bin,content);
        diff_histo_cent->SetBinError(bin,error);
      }

      if(m_useChebyshevSpacing){
        int firstBin = HistogramUtilities::firstBinWithContent(match_histo_cent);
        int lastBin  = HistogramUtilities::lastBinWithContent(match_histo_cent);
        double lowEdgeFirstBin = match_histo_cent->GetBinLowEdge(firstBin);
        double highEdgeLastBin = match_histo_cent->GetBinLowEdge(lastBin+1);
        TH1D* rebinMatchTrack =  HistogramUtilities::rebinHistoToChebyshevSpacing(match_histo_cent, lowEdgeFirstBin, highEdgeLastBin,3,20);
        TH1D* rebinEmbTrack = (TH1D*) rebinMatchTrack->Clone();
        HistogramUtilities::zeroBins(rebinEmbTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(emb_histo_cent,rebinEmbTrack);
        delete emb_histo_cent;
        delete match_histo_cent;
        emb_histo_cent   = rebinEmbTrack;
        match_histo_cent = rebinMatchTrack;
        HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
        HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);
      }else{
        //int firstBin = HistogramUtilities::firstBinWithContent(match_histo_cent);
        //int lastBin  = HistogramUtilities::lastBinWithContent(match_histo_cent);
        //double lowEdgeFirstBin = match_histo_cent->GetBinLowEdge(firstBin);
        //double highEdgeLastBin = match_histo_cent->GetBinLowEdge(lastBin+1);

        /*
        TH1D* truncMatchHisto = HistogramUtilities::truncateHistogram(match_histo_cent,lowEdgeFirstBin,highEdgeLastBin);
        delete match_histo_cent;
        match_histo_cent = truncMatchHisto;
        TH1D* rebinMatchTrack = HistogramUtilities::rebinEqualByCDF(match_histo_cent, m_tpc_eff_no_fit_num_bins,false,false); // rebins so more bins where content is changing, less bins where there is less change
        TH1D* rebinEmbTrack = (TH1D*) rebinMatchTrack->Clone();
        HistogramUtilities::zeroBins(rebinEmbTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(emb_histo_cent,rebinEmbTrack);
        delete emb_histo_cent;
        delete match_histo_cent;
        emb_histo_cent   = rebinEmbTrack;
        match_histo_cent = rebinMatchTrack;
        HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
        HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);
        */

        //TH1D* truncDiffHisto = HistogramUtilities::truncateHistogram(diff_histo_cent,lowEdgeFirstBin,highEdgeLastBin);
        //delete diff_histo_cent;
        //diff_histo_cent = truncDiffHisto;
        //TH1D* rebinDiffTrack = HistogramUtilities::rebinEqualByCDF(diff_histo_cent, m_tpc_eff_no_fit_num_bins,false,false); // rebins so more bins where content is changing, less bins where there is less change
        TH1D* rebinMatchTrack = (TH1D*) rebinDiffTrack->Clone();
        TH1D* rebinEmbTrack   = (TH1D*) rebinDiffTrack->Clone();
        //HistogramUtilities::zeroBins(rebinMatchTrack);
        //HistogramUtilities::zeroBins(rebinEmbTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(match_histo_cent,rebinMatchTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(emb_histo_cent,rebinEmbTrack);
        delete emb_histo_cent;
        delete match_histo_cent;
        emb_histo_cent   = rebinEmbTrack;
        match_histo_cent = rebinMatchTrack;
        HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
        HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);


      }
      delete diff_histo_cent;

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  TH1D*'s Cent : " << emb_histo_cent  << " " << match_histo_cent <<  endl;
        cout << "     NumBins :   " << emb_histo_cent->GetNbinsX()  << " " << match_histo_cent->GetNbinsX() << endl;
      #endif

      m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
      m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
      m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality %02d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
         embTracksAllCent->GetXaxis()->GetBinLowEdge(rapBin),embTracksAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
      m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->Divide(match_histo_cent,emb_histo_cent,"cl=0.68268949 b(1,1) mode");

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "    TPC Eff  Cent " << centIndex << " Ptr: " << m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] << "  " << m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->GetName() << endl;
      #endif
      delete emb_histo_cent;
      delete match_histo_cent;
      //delete total_histo_cent;
    } // end cent index loop
  }// end rap loop



  cout << "  Making Pseudorapidity Projections for TPC Efficiency..." << endl;
  for(int etaBin = 1; etaBin <= embTracksAllCent_EtaPt->GetNbinsX(); etaBin++){
    int etaIndex = etaBin - 1;
    if(m_testingMidRapidity && etaBin != embTracksAllCent_EtaPt->GetXaxis()->FindBin(m_ycm)) continue;
    double eta = embTracksAllCent_EtaPt->GetXaxis()->GetBinCenter(etaBin);
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - eta) < 0.05){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating eta = " << eta << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }


    TH1D* embTrackHisto   = embTracksAllCent_EtaPt->ProjectionY(Form("proj_emb_eta%02d_%02d_%d",etaIndex,a_partIndex,pmIndex),etaBin,etaBin);
    TH1D* matchTrackHisto = matchTracksAllCent_EtaPt->ProjectionY(Form("proj_mat_eta%02d_%02d_%d",etaIndex,a_partIndex,pmIndex),etaBin,etaBin);
    m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][etaIndex][pmIndex] = new TGraphAsymmErrors();
    m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][etaIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_AllCent_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),etaIndex));
    m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][etaIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality Integrated #eta=[%1.2f,%1.2f]; m_{T}-m_{0}; TPC Efficiency",
               m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),embTracksAllCent_EtaPt->GetXaxis()->GetBinLowEdge(etaBin),embTracksAllCent_EtaPt->GetXaxis()->GetBinLowEdge(etaBin+1)));
    m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][etaIndex][pmIndex]->Divide(matchTrackHisto,embTrackHisto,"cl=0.68268949 b(1,1) mode");

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "    TPC Eff All Cent Ptr: " << m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][etaIndex][pmIndex] << "    " <<  m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][etaIndex][pmIndex]->GetName() << endl;
    #endif

    delete embTrackHisto;
    delete matchTrackHisto;

    for(int centIndex = 0; centIndex < numLoadedCentBins; centIndex++){
      if(!embeddedTrackHistogramRebin_EtaPt[centIndex]){
        m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][etaIndex][pmIndex] = NULL;
        continue;
      }
      if(m_testingOneCent && centIndex != 0) continue;

      TH1D* emb_histo_cent   = embeddedTrackHistogramRebin_EtaPt[centIndex]->ProjectionY(Form("proj_emb_eta%02d_cent%02d_%02d_%d",etaIndex,centIndex,a_partIndex,pmIndex),etaBin,etaBin);
      TH1D* match_histo_cent = matchTrackHistogramRebin_EtaPt[centIndex]->ProjectionY(Form("proj_match_eta%02d_cent%02d_%02d_%d",etaIndex,centIndex,a_partIndex,pmIndex),etaBin,etaBin);

      int firstBin = HistogramUtilities::firstBinWithContent(match_histo_cent);
      int lastBin  = HistogramUtilities::lastBinWithContent(match_histo_cent);
      double lowEdgeFirstBin = match_histo_cent->GetBinLowEdge(firstBin);
      double highEdgeLastBin = match_histo_cent->GetBinLowEdge(lastBin+1);
      TH1D* truncMatchHisto = HistogramUtilities::truncateHistogram(match_histo_cent,lowEdgeFirstBin,highEdgeLastBin);
      delete match_histo_cent;
      match_histo_cent = truncMatchHisto;
      TH1D* rebinMatchTrack = HistogramUtilities::rebinEqualByCDF(match_histo_cent, 80,false,false); // rebins so more bins where content is changing, less bins where there is less change
      TH1D* rebinEmbTrack = (TH1D*) rebinMatchTrack->Clone();
      HistogramUtilities::zeroBins(rebinEmbTrack);
      HistogramUtilities::rebinHistoToAlternateBinStructure(emb_histo_cent,rebinEmbTrack);
      delete emb_histo_cent;
      delete match_histo_cent;
      emb_histo_cent   = rebinEmbTrack;
      match_histo_cent = rebinMatchTrack;
      HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
      HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);


      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  TH1D*'s Cent : " << emb_histo_cent  << " " << match_histo_cent <<  endl;
        cout << "     NumBins :   " << emb_histo_cent->GetNbinsX()  << " " << match_histo_cent->GetNbinsX() << endl;
      #endif

      m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][etaIndex][pmIndex] = new TGraphAsymmErrors();
      m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][etaIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_Cent%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,etaIndex));
      m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][etaIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality %02d #eta=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
         embTracksAllCent_EtaPt->GetXaxis()->GetBinLowEdge(etaBin),embTracksAllCent_EtaPt->GetXaxis()->GetBinLowEdge(etaBin+1)));
      m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][etaIndex][pmIndex]->Divide(match_histo_cent,emb_histo_cent,"cl=0.68268949 b(1,1) mode");

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "    TPC Eff  Cent " << centIndex << " Ptr: " << m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][etaIndex][pmIndex] << "  " << m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][etaIndex][pmIndex]->GetName() << endl;
      #endif
      delete emb_histo_cent;
      delete match_histo_cent;
      //delete total_histo_cent;
    } // end cent index loop
  }// end rap loop











  //######################      Make The BTOF Projections       ########################
  cout << "  Making Projections for BTOF Efficiency..." << endl;
  for(int rapBin = 1; rapBin <= tpc_forBTOFAllCent->GetNbinsX(); rapBin++){
    int rapIndex = rapBin - 1;
    if(m_testingMidRapidity && rapBin != tpc_forBTOFAllCent->GetXaxis()->FindBin(m_ycm)) continue;
    double rapidity = tpc_forBTOFAllCent->GetXaxis()->GetBinCenter(rapBin);
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - rapidity) < 0.005){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating rapidity = " << rapidity << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }


    TH1D* tpcTrackHisto = tpc_forBTOFAllCent->ProjectionY(Form("proj_tpc_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);
    TH1D* btofTrackHsito = btof_forBTOFAllCent->ProjectionY(Form("proj_btof_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "  ## Rapidity Bin = " << rapBin << endl;
      cout << "  TH1D*'s : " << tpcTrackHisto  << " " << btofTrackHsito << endl;
      cout << "  NumBins : " << tpcTrackHisto->GetNbinsX()  << " " << btofTrackHsito->GetNbinsX() << endl;

    #endif

    m_BTOF_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
    m_BTOF_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("btofEfficiencyGraph_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
    m_BTOF_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("BTOF Effiency %s Centrality Integrated y=[%1.2f,%1.2f]; m_{T}-m_{0}; BTOF Efficiency",
               m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin),tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
    m_BTOF_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->Divide(btofTrackHsito,tpcTrackHisto,"cl=0.68268949 b(1,1) mode");

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "    BTOF Eff All Cent Ptr: " << m_BTOF_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] << "    " <<  m_BTOF_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->GetName() << endl;
    #endif

    delete tpcTrackHisto;
    delete btofTrackHsito;

    for(int centIndex = 0; centIndex < numLoadedCentBins; centIndex++){
      if(!btof_forBTOFHistogramRebin[centIndex]){
        m_BTOF_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = NULL;
        continue;
      }
      if(m_testingOneCent && centIndex != 0) continue;

      TH1D* tpcTrackHistoCent = tpc_forBTOFHistogramRebin[centIndex]->ProjectionY(Form("proj_tpc_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);
      TH1D* btofTrackHistoCent = btof_forBTOFHistogramRebin[centIndex]->ProjectionY(Form("proj_btof_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  TH1D*'s Cent : " << tpcTrackHistoCent  << " " << btofTrackHistoCent <<  endl;
        cout << "     NumBins :   " << tpcTrackHistoCent->GetNbinsX()  << " " << btofTrackHistoCent->GetNbinsX() << endl;
      #endif

      m_BTOF_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
      m_BTOF_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("btofEfficiencyGraph_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
      m_BTOF_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("BTOF Effiency %s Centrality %02d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
         tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin),tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
      m_BTOF_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->Divide(btofTrackHistoCent,tpcTrackHistoCent,"cl=0.68268949 b(1,1) mode");

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "    TPC Eff  Cent " << centIndex << " Ptr: " << m_BTOF_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] << "  " << m_BTOF_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->GetName() << endl;
      #endif

      delete tpcTrackHistoCent;
      delete btofTrackHistoCent;
      //delete total_histo_cent;

    } // end cent index loop
  }// end rap loop



  //####  NOTE: If your statistics changes significantly from bin to bin (due to multiple embedding samples from different pT Eta ranges)
  //            Since y and pT reco are shifts of the embedding sample, for a bin on edge of statistics change, you wont have contributions
  //            from higher pT.

  TH2D* energyLossAllCent_Counts = HistogramUtilities::projection2D(energyLossAllCent, 1, energyLossAllCent->GetNbinsZ());







  /*
  //######################         Make the smooth reponse matrix pdf            ########################
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "  Making Projections for pT Response Matrix Efficiency..." << endl;
  #endif
  for(int rapBin = 1; rapBin <= energyLossAllCent->GetNbinsX(); rapBin++){
    int rapIndex = rapBin - 1;
    if(m_testingMidRapidity && rapBin != energyLossAllCent->GetXaxis()->FindBin(m_ycm)) continue;
    double rapidity = energyLossAllCent->GetXaxis()->GetBinCenter(rapBin);
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - rapidity) < 0.005){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating rapidity = " << rapidity << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }

    double previousGoodParameters[4] = {-999,-999,-999,-999};
    for(int mTm0Bin = 1; a_doEnergyLossFits && mTm0Bin <= energyLossAllCent->GetNbinsY(); mTm0Bin++){
      TH1D* mTLossHistogram = HistogramUtilities::projectionZ(m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex],rapBin,mTm0Bin,Form("_%d_%d_%lu",a_partIndex,pmIndex,getUniqueID()));
      double totalCounts = HistogramUtilities::sumOfContent(mTLossHistogram,-1,false);
      if(totalCounts < 10.0){
        delete mTLossHistogram;
        continue;
      }
      TH1D* histoToDelete = mTLossHistogram;
      mTLossHistogram =  HistogramUtilities::rebinEqualByCDF(mTLossHistogram, 20,false, false);
      delete histoToDelete;
      double mTm0_or_pT = m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetYaxis()->GetBinCenter(mTm0Bin);
      totalCounts = HistogramUtilities::sumOfContent(mTLossHistogram,-1,false);
      int firstBinWithContent = HistogramUtilities::firstBinWithContent(mTLossHistogram);
      int lastBinWithContent  = HistogramUtilities::lastBinWithContent(mTLossHistogram);
      double lowRangeForData  = mTLossHistogram->GetBinLowEdge(firstBinWithContent);
      double highRangeForData = mTLossHistogram->GetBinLowEdge(lastBinWithContent+1);

      if(m_pTLossMode){
        mTLossHistogram->SetTitle(Form("p_{T} response for %s Centrality Integrated y=[%1.2f,%1.2f] p_{T}^{reco}=[%1.2f,%1.2f]",
          m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),
          m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetXaxis()->GetBinLowEdge(rapBin),m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetXaxis()->GetBinLowEdge(rapBin+1),
          m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetYaxis()->GetBinLowEdge(mTm0Bin),m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetYaxis()->GetBinLowEdge(mTm0Bin+1)));
        mTLossHistogram->GetXaxis()->SetTitle("p_{T}^{emb} (GeV/c)");
        mTLossHistogram->GetYaxis()->SetTitle("dN_{Tracks}/dp_{T}^{emb} (c/GeV)");
      }else{
        mTLossHistogram->SetTitle(Form("m_{T} response for %s Centrality Integrated y=[%1.2f,%1.2f] m_{T}-m_{0}^{reco}=[%1.2f,%1.2f]",
          m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),
          m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetXaxis()->GetBinLowEdge(rapBin),m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetXaxis()->GetBinLowEdge(rapBin+1),
          m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetYaxis()->GetBinLowEdge(mTm0Bin),m_mTm0ReponseMatrixHisto_PDF[a_partIndex][pmIndex]->GetYaxis()->GetBinLowEdge(mTm0Bin+1)));
        mTLossHistogram->GetXaxis()->SetTitle("m_{T}^{emb} (GeV/c^{2})");
        mTLossHistogram->GetYaxis()->SetTitle("dN_{Tracks}/dm_{T}^{emb} (c^{2}/GeV)");
      }
      HistogramUtilities::makeSqrtNErrors(mTLossHistogram);
      HistogramUtilities::makeDifferential(mTLossHistogram);


      / * double gausData[3];
      if(previousGoodParameters[0] < 0.0){
        HistogramUtilities::getGausData(gausData, mTLossHistogram, true, false);
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "    Gaus Data: " << gausData[0] << "  " << gausData[1] << "  " << gausData[2] << endl;
        #endif
      }* /

      double amp_init  = 1.0;
      double mean_init = mTLossHistogram->GetMean();//gausData[1];
      double sig_init  = mTLossHistogram->GetStdDev();// gausData[2];
      double skew_init = 0.0;

      /*
      if(previousGoodParameters[0] > 0.0){
        amp_init  = 1.0;
        mean_init = previousGoodParameters[1];
        sig_init  = previousGoodParameters[2];
        skew_init = previousGoodParameters[3];
      }
      * /

      string gausFunctName = Form("response_matrix_funct%03d_%03d_%03d_%03d",rapIndex,a_partIndex,pmIndex,mTm0Bin);
      TF1* gausFunct = PhysMath::skewNormalFunct(1.0, mean_init, sig_init, skew_init);
      gausFunct->SetName(gausFunctName.c_str());
      gausFunct->SetRange(HistogramUtilities::getLowEdge(mTLossHistogram),HistogramUtilities::getHighEdge(mTLossHistogram));
      int highestBin = HistogramUtilities::getMaxBinInRange(mTLossHistogram);
      gausFunct->SetParameter(0,mTLossHistogram->GetBinContent(highestBin)/gausFunct->Eval(mTLossHistogram->GetBinCenter(highestBin)));
      gausFunct->SetParLimits(0,0.01*gausFunct->GetParameter(0), 10.0*gausFunct->GetParameter(0));
      gausFunct->SetParLimits(1,HistogramUtilities::getLowEdge(mTLossHistogram), HistogramUtilities::getHighEdge(mTLossHistogram));
      //gausFunct->SetParLimits(2,0.0, 150.0*gausFunct->GetParameter(2));
      gausFunct->SetParLimits(2,0.0, 3.0);
      gausFunct->SetParLimits(3,-1.0,1.0);
      if(m_pTLossMode){
        //gausFunct->SetParameter(0,totalCounts);
        //gausFunct->SetParLimits(0,0.0,100.0*totalCounts);
        gausFunct->SetParLimits(1,-0.5, 0.5);
        gausFunct->SetParLimits(2,0.0, 1.0);
      }

      TF1* initialFunct = (TF1*) gausFunct->Clone();
      initialFunct->SetLineColor(kYellow);
      initialFunct->SetName("initialFunct");


      string fitOptions = "I";
      #ifndef _EFFICIENCY_FITTER_DEBUG_
        fitOptions += " Q";
      #endif
      int status = mTLossHistogram->Fit(gausFunctName.c_str(),fitOptions.c_str());
      TF1* fitResult = mTLossHistogram->GetFunction(gausFunctName.c_str());
      if(!fitResult || status != 0){
        cout << "  ERROR: Energy Loss Calculation Failed with Fit Status: " << status << " and ptr: " << fitResult << endl;
        delete mTLossHistogram;
        delete gausFunct;
        delete initialFunct;
        previousGoodParameters[0] = -999;
        previousGoodParameters[1] = -999;
        previousGoodParameters[2] = -999;
        previousGoodParameters[3] = -999;
        continue;
      }

      previousGoodParameters[0] = fitResult->GetParameter(0);
      previousGoodParameters[1] = fitResult->GetParameter(1);
      previousGoodParameters[2] = fitResult->GetParameter(2);
      previousGoodParameters[3] = fitResult->GetParameter(3);

      if(a_savePtLossGausImages){
        fittingCanvas->cd();
        if(m_pTLossMode) mTLossHistogram->GetXaxis()->SetRangeUser(-0.035,0.035);
        gPad->SetRightMargin(0.45);
        mTLossHistogram->Draw("E");
        initialFunct->Draw("SAME");
        gPad->Update();
        gSystem->ProcessEvents();
        fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_AllCent/mTLossFits/ResponseFit_RapIndex%03d_mTm0Index%03d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex,mTm0Bin - 1));
      }
      delete gausFunct;
      delete mTLossHistogram;
      delete initialFunct;
    }
  }
  */














  //######################      Make The pT-Loss graphs       ########################
  cout << "  Making Projections for pT-Loss Efficiency..." << endl;
  vector< TH2D* > pT2DProjections;
  if(m_useEnergyLossTProfiles) pT2DProjections = HistogramUtilities::projeciton2DSlices_X(energyLossAllCent);
  for(int rapBin = 1; rapBin <= energyLossAllCent->GetNbinsX(); rapBin++){
    int rapIndex = rapBin - 1;
    if(m_testingMidRapidity && rapBin != energyLossAllCent->GetXaxis()->FindBin(m_ycm)) continue;
    double rapidity = energyLossAllCent->GetXaxis()->GetBinCenter(rapBin);
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - rapidity) < 0.005){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating rapidity = " << rapidity << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }


    if(m_useEnergyLossTProfiles){
      energyLossAllCent->GetXaxis()->SetRange(rapBin,rapBin);
      //TH2D* proj2DTemp = (TH2D*) energyLossAllCent->Project3D("zy");
      TH2D* proj2DTemp = pT2DProjections[rapBin - 1];
      for(int mTm0Bin = 1; mTm0Bin <= energyLossAllCent->GetNbinsY(); mTm0Bin++){
        m_EnergyLossProfile_AllCent[a_partIndex][rapIndex][pmIndex] = proj2DTemp->ProfileX(Form("energyLossProfile_%s_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex),mTm0Bin,mTm0Bin);
        if(m_pTLossMode){
          m_EnergyLossProfile_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("Energy Loss %s Centrality Integrated y=[%1.2f,%1.2f]; p_{T} (GeV/c); p_{T}^{Emb}-p_{T}^{Reco} (GeV/c)",
                 m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
        }else{
          m_EnergyLossProfile_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("Energy Loss %s Centrality Integrated y=[%1.2f,%1.2f]; m_{T} (GeV/c^{2}); m_{T}^{Emb}-m_{T}^{Reco} (GeV/c^{2})",
                 m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
        }
      }

    }else{

      m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
      m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("energyLossGraph_%s_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
      if(m_pTLossMode){
        m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("Energy Loss %s Centrality Integrated y=[%1.2f,%1.2f]; p_{T} (GeV/c; p_{T}^{Emb}-p_{T}^{Reco} (GeV/c)",
                 m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
      }else{
        m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("Energy Loss %s Centrality Integrated y=[%1.2f,%1.2f]; m_{T} (GeV/c^{2}); m_{T}^{Emb}-m_{T}^{Reco} (GeV/c^{2})",
                 m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
      }

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Creating energy loss graph: " << m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->GetName() << endl;
      #endif

      double previousStdDev = -1.0;
      for(int mTm0Bin = 1; a_doEnergyLossFits && mTm0Bin <= energyLossAllCent->GetNbinsY(); mTm0Bin++){
        TH1D* mTLossHistogram = HistogramUtilities::projectionZ(energyLossAllCent,rapBin,mTm0Bin,Form("_%d_%d_%lu",a_partIndex,pmIndex,getUniqueID()));
        double totalCounts = HistogramUtilities::sumOfContent(mTLossHistogram,-1,false);
        if(totalCounts < 10.0){
          delete mTLossHistogram;
          continue;
        }
        TH1D* histoToDelete = mTLossHistogram;
        mTLossHistogram =  HistogramUtilities::rebinEqualByCDF(mTLossHistogram, 15,false, false);
        delete histoToDelete;
        double mTm0_or_pT = energyLossAllCent->GetYaxis()->GetBinCenter(mTm0Bin);
        totalCounts = HistogramUtilities::sumOfContent(mTLossHistogram,-1,false);
        int firstBinWithContent = HistogramUtilities::firstBinWithContent(mTLossHistogram);
        int lastBinWithContent  = HistogramUtilities::lastBinWithContent(mTLossHistogram);
        double lowRangeForData  = mTLossHistogram->GetBinLowEdge(firstBinWithContent);
        double highRangeForData = mTLossHistogram->GetBinLowEdge(lastBinWithContent+1);

        if(m_pTLossMode){
          mTLossHistogram->SetTitle(Form("p_{T} loss for %s Centrality Integrated y=[%1.2f,%1.2f] p_{T}=[%1.2f,%1.2f]",
            m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),
            energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1),
            energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin),energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin+1)));
          mTLossHistogram->GetXaxis()->SetTitle("p_{T}^{emb} - p_{T}^{reco} (GeV/c)");
          mTLossHistogram->GetYaxis()->SetTitle("dN_{Tracks}/d(p_{T}^{emb} - p_{T}^{reco}) (c/GeV)");
        }else{
          mTLossHistogram->SetTitle(Form("m_{T} loss for %s Centrality Integrated y=[%1.2f,%1.2f] m_{T}-m_{0}=[%1.2f,%1.2f]",
            m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),
            energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1),
            energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin),energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin+1)));
          mTLossHistogram->GetXaxis()->SetTitle("m_{T}^{emb} - m_{T}^{reco} (GeV)");
          mTLossHistogram->GetYaxis()->SetTitle("dN_{Tracks}/d(m_{T}^{emb} - m_{T}^{reco}) (c^{2}/GeV)");
        }

        HistogramUtilities::makeSqrtNErrors(mTLossHistogram);
        HistogramUtilities::makeDifferential(mTLossHistogram);
        //mTLossHistogram->Print();

        double gausData[3];
        HistogramUtilities::getGausData(gausData, mTLossHistogram, true, false);
        /*if(totalCounts < 30.0){
          int numPoints = m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->GetN();
          m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->SetPoint(numPoints,energyLossAllCent->GetYaxis()->GetBinCenter(mTm0Bin),gausData[1]);
          double standardError = gausData[2]/sqrt(gausData[0]);

          m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->SetPointError(numPoints,
                 0.5*energyLossAllCent->GetYaxis()->GetBinWidth(mTm0Bin),0.5*energyLossAllCent->GetYaxis()->GetBinWidth(mTm0Bin),
                 standardError,standardError);
          delete mTLossHistogram;
          continue;
        }*/


        //cout << "    Gaus Data: " << gausData[0] << "  " << gausData[1] << "  " << gausData[2] << endl;
        string gausFunctName = Form("gausFunct%03d_%03d_%03d_%03d",rapIndex,a_partIndex,pmIndex,mTm0Bin);
        TF1* gausFunct = PhysMath::getGaussianFunction(gausFunctName,1.0, gausData[1], gausData[2]);
        gausFunct->SetRange(HistogramUtilities::getLowEdge(mTLossHistogram),HistogramUtilities::getHighEdge(mTLossHistogram));
        //TF1* gausFunct = new TF1(gausFunctName.c_str(),"[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",HistogramUtilities::getLowEdge(mTLossHistogram),HistogramUtilities::getHighEdge(mTLossHistogram));
        //gausFunct->SetParameter(0,1.0);
        //gausFunct->SetParameter(1,gausData[1]);
        //gausFunct->SetParameter(2,gausData[2]);
        gausFunct->SetParameter(0,mTLossHistogram->GetBinContent(mTLossHistogram->FindBin(gausData[1]))/gausFunct->Eval(gausData[1]));

        gausFunct->SetParLimits(0,0.01*gausFunct->GetParameter(0), 10.0*gausFunct->GetParameter(0));
        gausFunct->SetParLimits(1,HistogramUtilities::getLowEdge(mTLossHistogram), HistogramUtilities::getHighEdge(mTLossHistogram));
        gausFunct->SetParLimits(2,0.0, 150.0*gausFunct->GetParameter(2));
        if(m_pTLossMode){
          gausFunct->SetParameter(0,totalCounts);
          gausFunct->SetParLimits(0,0.0,100.0*totalCounts);
          gausFunct->SetParameter(0,0.0);
          gausFunct->SetParLimits(1,-0.5, 0.5);
          gausFunct->SetParLimits(2,0.0, 1.0);

        }
          //cout << "  Check: " << mTLossHistogram->GetName()   << "     " << mTLossHistogram << endl;
          //cout << "         " << gausFunct->GetName()   << "     " << gausFunct << endl;

        string fitOptions = "I E";
        //if(totalCounts < 500) fitOptions = fitOptions + "L";
        int status = mTLossHistogram->Fit(gausFunctName.c_str(),fitOptions.c_str());
        //cout << "Fit Status = " << status << endl;
        TF1* fitResult = mTLossHistogram->GetFunction(gausFunctName.c_str());
        //cout << "  result fit ptr: " << fitResult << endl;
        if(!fitResult || status != 0){
          cout << "  ERROR: Energy Loss Calculation Failed with Fit Status: " << status << " and ptr: " << fitResult << endl;
          delete mTLossHistogram;
          delete gausFunct;
          continue;
        }

        //cout << "   Fit::::   N: " << fitResult->GetParameter(0)*sqrt(2.0*TMath::Pi())*fitResult->GetParameter(2) << "  mu: " << fitResult->GetParameter(1) << "  sig: " << fitResult->GetParameter(2) << endl;
        //cout << "   Raw::::   N: " << gausData[0]*mTLossHistogram->GetBinWidth(1) << "  mu: " << gausData[1] << "  sig: " << gausData[2] << endl;


        //figure out if this bin should be used
        cout << "WARNING: Removed code that removes energy loss points on line 1073ish ... maybe put it back in at a later date" << endl;
        bool useThisPoint = true;
        /*if(fitResult->GetParameter(1) + 0.5*fitResult->GetParameter(2) < 0.0){
          useThisPoint = false;
          if(previousStdDev > 0.0){
            for(int pointIndex = 0;pointIndex < m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->GetN(); pointIndex++){
              double X_ofPoint = m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->GetX()[pointIndex];
              if(fabs(mTm0_or_pT - X_ofPoint) < previousStdDev){
                m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->RemovePoint(pointIndex);
                pointIndex--;
              }
            }
          }
        }else{
          previousStdDev = fitResult->GetParameter(2);
        }*/
        /*if(mTm0Bin > 10){
          double beforeAvg = energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin - 2) + energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin - 1);
          double afterAvg  = energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin + 1) + energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin + 2);
          double ratioAvg = afterAvg/beforeAvg;
          if(ratioAvg < 0.9 || ratioAvg > 1.1) useThisPoint = false;
        }*/

        if(useThisPoint){
          int numPoints = m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->GetN();
          m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->SetPoint(numPoints,embTracksAllCent->GetYaxis()->GetBinCenter(mTm0Bin),fitResult->GetParameter(1));
          m_EnergyLoss_AllCent[a_partIndex][rapIndex][pmIndex]->SetPointError(numPoints,
                 0.5*embTracksAllCent->GetYaxis()->GetBinWidth(mTm0Bin),0.5*embTracksAllCent->GetYaxis()->GetBinWidth(mTm0Bin),
                 fitResult->GetParError(1), fitResult->GetParError(1));
        }

        if(a_savePtLossGausImages){
          fittingCanvas->cd();
          if(m_pTLossMode) mTLossHistogram->GetXaxis()->SetRangeUser(-0.035,0.035);
          gPad->SetRightMargin(0.45);
          mTLossHistogram->Draw("E");
          gPad->Update();
          gSystem->ProcessEvents();
          fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_AllCent/mTLossFits/mTLossFit_RapIndex%03d_mTm0Index%03d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex,mTm0Bin - 1));
        }
        delete gausFunct;
        delete mTLossHistogram;
      }
    }

    for(int centIndex = 0; a_doEnergyLossFits && a_fitEnergyLossByCent && centIndex < numLoadedCentBins; centIndex++){
      if(!energyLossHistogramRebin[centIndex]){
        m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = NULL;
        continue;
      }
      if(m_testingOneCent && centIndex != 0) continue;


      if(m_useEnergyLossTProfiles){
        energyLossAllCent->GetXaxis()->SetRange(rapBin,rapBin);
        TH2D* proj2DTemp = (TH2D*) energyLossAllCent->Project3D("zy");
        for(int mTm0Bin = 1; mTm0Bin <= energyLossAllCent->GetNbinsY(); mTm0Bin++){
          m_EnergyLossProfile_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = proj2DTemp->ProfileY(Form("energyLossProfile_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex),mTm0Bin,mTm0Bin);
          if(m_pTLossMode){
            m_EnergyLossProfile_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("Energy Loss %s Centrality %02d y=[%1.2f,%1.2f]; p_{T} (GeV/c); p_{T}^{Emb}-p_{T}^{Reco} (GeV/c)",
                   m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
          }else{
            m_EnergyLossProfile_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("Energy Loss %s Centrality %02d y=[%1.2f,%1.2f]; m_{T} (GeV/c^{2}); m_{T}^{Emb}-m_{T}^{Reco} (GeV/c^{2})",
                   m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
          }
        }

      }else{

        m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
        m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("energyLossGraph_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
        m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("Energy Loss %s Centrality Index %02d y=[%1.2f,%1.2f]",
                 m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));

        for(int mTm0Bin = 1; mTm0Bin <= energyLossHistogramRebin[centIndex]->GetNbinsZ(); mTm0Bin++){
          //TH1D* mTLossHistogram = energyLossHistogramRebin[centIndex]->ProjectionZ(Form("proj_loss_rap%02d_part%02d_pm%d_mTm0%02d_Cent%02d",rapIndex,a_partIndex,pmIndex,mTm0Bin,centIndex),rapBin,rapBin,mTm0Bin,mTm0Bin);
          TH1D* mTLossHistogram = HistogramUtilities::projectionZ(energyLossHistogramRebin[centIndex],rapBin,mTm0Bin,Form("_%d_%d_%lu",a_partIndex,pmIndex,getUniqueID()));
          double totalCounts = HistogramUtilities::sumOfContent(mTLossHistogram,-1,false);
          if(totalCounts < 10.0){
            delete mTLossHistogram;
            continue;
          }
          TH1D* histoToDelete = mTLossHistogram;
          mTLossHistogram =  HistogramUtilities::rebinEqualByCDF(mTLossHistogram, 25,false, false);
          delete histoToDelete;
          totalCounts = HistogramUtilities::sumOfContent(mTLossHistogram,-1,false);

          if(m_pTLossMode){
            mTLossHistogram->SetTitle(Form("p_{T} loss for %s Centrality %02d y=[%1.2f,%1.2f] p_{T}=[%1.2f,%1.2f]",
              m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),centIndex,
              energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1),
              energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin),energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin+1)));
            mTLossHistogram->GetXaxis()->SetTitle("p_{T}^{emb} - p_{T}^{reco} (GeV/c)");
            mTLossHistogram->GetYaxis()->SetTitle("dN_{Tracks}/d(p_{T}^{emb} - p_{T}^{reco}) (GeV^{-1})");
            HistogramUtilities::makeSqrtNErrors(mTLossHistogram);
            HistogramUtilities::makeDifferential(mTLossHistogram);
          }else{
            mTLossHistogram->SetTitle(Form("p_{T} loss for %s Centrality %02d y=[%1.2f,%1.2f] m_{T}-m_{0}=[%1.2f,%1.2f]",
              m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),centIndex,
              energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin),energyLossAllCent->GetXaxis()->GetBinLowEdge(rapBin+1),
              energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin),energyLossAllCent->GetYaxis()->GetBinLowEdge(mTm0Bin+1)));
            mTLossHistogram->GetXaxis()->SetTitle("m_{T}^{emb} - m_{T}^{reco} (GeV/c^{2})");
            mTLossHistogram->GetYaxis()->SetTitle("dN_{Tracks}/d(m_{T}^{emb} - m_{T}^{reco}) (GeV^{-1})");
            HistogramUtilities::makeSqrtNErrors(mTLossHistogram);
            HistogramUtilities::makeDifferential(mTLossHistogram);
          }

          //mTLossHistogram->Print();


          double gausData[3];
          HistogramUtilities::getGausData(gausData, mTLossHistogram, true, false);

          if(totalCounts < 30.0){

            double standardError = gausData[2]/sqrt(gausData[0]);
            //cout << "    Gaus Data: " << gausData[0] << "  " << gausData[1] << "  " << gausData[2] << "  " << standardError << endl;
            int numPoints = m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->GetN();
            m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetPoint(numPoints,energyLossAllCent->GetYaxis()->GetBinCenter(mTm0Bin),gausData[1]);
            m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetPointError(numPoints,
                   0.5*energyLossAllCent->GetYaxis()->GetBinWidth(mTm0Bin),0.5*energyLossAllCent->GetYaxis()->GetBinWidth(mTm0Bin),
                   standardError,standardError);
            delete mTLossHistogram;
            continue;
          }

          //cout << "    Gaus Data: " << gausData[0] << "  " << gausData[1] << "  " << gausData[2] << endl;
          TF1* gausFunct = new TF1(Form("gausFunct%03d_%03d_%03d_%03d_%03d",rapIndex,a_partIndex,pmIndex,mTm0Bin,centIndex),"[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",HistogramUtilities::getLowEdge(mTLossHistogram),HistogramUtilities::getHighEdge(mTLossHistogram));
          gausFunct->SetParameter(0,1.0);
          gausFunct->SetParameter(1,gausData[1]);
          gausFunct->SetParameter(2,gausData[2]);
          gausFunct->SetParameter(0,mTLossHistogram->GetBinContent(mTLossHistogram->FindBin(gausData[1]))/gausFunct->Eval(gausData[1]));

          gausFunct->SetParLimits(0,0.5*gausFunct->GetParameter(0), 2.0*gausFunct->GetParameter(0));
          gausFunct->SetParLimits(1,HistogramUtilities::getLowEdge(mTLossHistogram), HistogramUtilities::getHighEdge(mTLossHistogram));
          gausFunct->SetParLimits(2,0.0, 1.5*gausFunct->GetParameter(2));


          //cout << "  Check: " << mTLossHistogram->GetName()   << "     " << mTLossHistogram << endl;
          //cout << "         " << gausFunct->GetName()   << "     " << gausFunct << endl;

          string fitOptions = "QI";
          if(totalCounts < 500) fitOptions = fitOptions + "L";

          int status = mTLossHistogram->Fit(Form("gausFunct%03d_%03d_%03d_%03d_%03d",rapIndex,a_partIndex,pmIndex,mTm0Bin,centIndex),fitOptions.c_str());
          //cout << "Fit Status = " << status << endl;
          TF1* fitResult = mTLossHistogram->GetFunction(Form("gausFunct%03d_%03d_%03d_%03d_%03d",rapIndex,a_partIndex,pmIndex,mTm0Bin,centIndex));
          //cout << "  result fit ptr: " << fitResult << endl;
          if(!fitResult){
            delete mTLossHistogram;
            delete gausFunct;
            cout << " ERROR: energy loss fit for spceific mT-m0 bin (" << mTm0Bin << ") Failed!" << endl;
            continue;
          }
          //cout << "   Fit::::   N: " << fitResult->GetParameter(0)*sqrt(2.0*TMath::Pi())*fitResult->GetParameter(2) << "  mu: " << fitResult->GetParameter(1) << "  sig: " << fitResult->GetParameter(2) << endl;
          //cout << "   Raw::::   N: " << gausData[0]*mTLossHistogram->GetBinWidth(1) << "  mu: " << gausData[1] << "  sig: " << gausData[2] << endl;


          //figure out if this bin should be used
          bool useThisPoint = true;
          /*if(mTm0Bin > 5){
            double beforeAvg = energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin - 2) + energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin - 1);
            double afterAvg  = energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin + 1) + energyLossAllCent_Counts->GetBinContent(rapBin,mTm0Bin + 2);
            double ratioAvg = afterAvg/beforeAvg;
            if(ratioAvg < 0.9 || ratioAvg > 1.1) useThisPoint = false;
          }*/

          //if(mTm0Bin > energyLossAllCent_Counts->GetNbinsY() - 5 )useThisPoint = false;

          if(useThisPoint){
            int numPoints = m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->GetN();
            m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetPoint(numPoints,energyLossAllCent->GetYaxis()->GetBinCenter(mTm0Bin),fitResult->GetParameter(1));
            m_EnergyLoss_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetPointError(numPoints,
                   0.5*energyLossAllCent->GetYaxis()->GetBinWidth(mTm0Bin),0.5*energyLossAllCent->GetYaxis()->GetBinWidth(mTm0Bin),
                   fitResult->GetParError(1),fitResult->GetParError(1));
          }

          if(a_savePtLossGausImages){
            fittingCanvas->cd();
            gPad->SetRightMargin(0.45);
            mTLossHistogram->Draw("E");
            gPad->Update();
            gSystem->ProcessEvents();
            fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_ByCent/Cent%02d/mTLossFits/mTLossFit_RapIndex%03d_mTm0Index%03d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex,mTm0Bin - 1));
          }


          delete gausFunct;
          delete mTLossHistogram;

        }
      }

    } // end cent index loop

  }// end rap loop

  delete fittingCanvas;

}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination(int a_partIndex, int a_charge, vector<string> a_inFileNames,vector<vector<double>> a_VzRanges, string a_picobinner_file_for_Vz_data_histo, vector<double> a_refMultCuts){
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Call to EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination" << endl;
  #endif


  bool doRefMultCuts = (a_refMultCuts.size() > 0);
  int pmIndex = ((a_charge > 0) ? 0 : 1);

  // delete everything to make sure this function works
  for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
    HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]);
    HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex]);
    m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] = nullptr;
    m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex] = nullptr;
    for(int centIndex = 0; centIndex < 16; centIndex++){
      HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]);
      HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]);
      m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = nullptr;
      m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex] = nullptr;
    }
  }



  TFile* inFile_PicoBinner = new TFile(a_picobinner_file_for_Vz_data_histo.c_str(),"READ");
  if(!inFile_PicoBinner || !inFile_PicoBinner->IsOpen()){
    cout << "ERROR:  EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination.. couldn't find file: " << a_picobinner_file_for_Vz_data_histo << endl;
    return;
  }
  TH1I* histo_Vz_data_Int = (TH1I*) inFile_PicoBinner->Get("zVertex")->Clone();
  if(!histo_Vz_data_Int){
    cout << "Couldn't read Z Vertex Histo! from " << a_picobinner_file_for_Vz_data_histo << endl;
  }

  //inFile_PicoBinner->Close();
  TH1D* histo_Vz_data = HistogramUtilities::convertInt_to_DoubleHisto(histo_Vz_data_Int);

  const unsigned int numFiles = a_inFileNames.size();
  const unsigned int numVzBins = 28; // from the other
  TH1I* histo_Vz_structure = nullptr;

  TH2D* embTracksAllCentAllFiles[numVzBins];
  TH2D* matchTracksAllCentAllFiles[numVzBins];
  TH2D* embTracksAllCent_EtaPtAllFiles[numVzBins];
  TH2D* matchTracksAllCent_EtaPtAllFiles[numVzBins];

  TH2D* embTracksAllCent[numVzBins][numFiles];
  TH2D* matchTracksAllCent[numVzBins][numFiles];
  TH2D* embTracksAllCent_EtaPt[numVzBins][numFiles];
  TH2D* matchTracksAllCent_EtaPt[numVzBins][numFiles];

  TH2D* embeddedTrackHistogram[16][numVzBins][numFiles];
  TH2D* matchTrackHistogram[16][numVzBins][numFiles];
  TH2D* embeddedTrackHistogram_EtaPt[16][numVzBins][numFiles];
  TH2D* matchTrackHistogram_EtaPt[16][numVzBins][numFiles];

  TH3D* embTrackHisto_byRefMult[numVzBins][numFiles];
  TH3D* matchTrackHisto_byRefMult[numVzBins][numFiles];

  //TGraphAsymmErrors* graph_eff_AllCent[g_EffFit_MaxNumRapBins][numVzBins][numFiles];
  //TGraphAsymmErrors* graph_eff_EtaPt_AllCent[g_EffFit_MaxNumRapBins][numVzBins][numFiles];
  //TGraphAsymmErrors* graph_eff[16][g_EffFit_MaxNumRapBins][numVzBins][numFiles];
  //TGraphAsymmErrors* graph_eff_EtaPt[16][g_EffFit_MaxNumRapBins][numVzBins][numFiles];

  HistogramUtilities::nullptr_array_1D(embTracksAllCentAllFiles, numVzBins);
  HistogramUtilities::nullptr_array_1D(matchTracksAllCentAllFiles, numVzBins);
  HistogramUtilities::nullptr_array_1D(embTracksAllCent_EtaPtAllFiles, numVzBins);
  HistogramUtilities::nullptr_array_1D(matchTracksAllCent_EtaPtAllFiles, numVzBins);


  for(unsigned int vZ_index = 0; vZ_index < numVzBins; vZ_index++){
    for(unsigned int file_index = 0; file_index < numFiles; file_index++){
      embTracksAllCent[vZ_index][file_index] = nullptr;
      matchTracksAllCent[vZ_index][file_index] = nullptr;
      embTracksAllCent_EtaPt[vZ_index][file_index] = nullptr;
      matchTracksAllCent_EtaPt[vZ_index][file_index] = nullptr;
      embTrackHisto_byRefMult[vZ_index][file_index] = nullptr;
      matchTrackHisto_byRefMult[vZ_index][file_index] = nullptr;
    }
  }

  for(unsigned int centIndex = 0; centIndex < 16; centIndex++){
    for(unsigned int vZ_index = 0; vZ_index < numVzBins; vZ_index++){
      for(unsigned int file_index = 0; file_index < numFiles; file_index++){
        embeddedTrackHistogram[centIndex][vZ_index][file_index] = nullptr;
        matchTrackHistogram[centIndex][vZ_index][file_index] = nullptr;
        embeddedTrackHistogram_EtaPt[centIndex][vZ_index][file_index] = nullptr;
        matchTrackHistogram_EtaPt[centIndex][vZ_index][file_index] = nullptr;
      }
    }
  }

  /*for(unsigned int rap_index = 0; rap_index < g_EffFit_MaxNumRapBins; rap_index++){
    for(unsigned int file_index = 0; file_index < numFiles; file_index++){
      for(unsigned int vZ_index = 0; vZ_index < numVzBins; vZ_index++){
        for(unsigned int centIndex = 0; centIndex < 16; centIndex++){
          graph_eff[centIndex][rap_index][vZ_index][file_index] = nullptr;
          graph_eff_EtaPt[centIndex][rap_index][vZ_index][file_index] = nullptr;
        }
        graph_eff_AllCent[rap_index][vZ_index][file_index] = nullptr;
        graph_eff_EtaPt_AllCent[rap_index][vZ_index][file_index] = nullptr;
      }
    }
  }*/



  //################   READ IN THE DATA  ##################################
  for(unsigned int file_index = 0; file_index < a_inFileNames.size(); file_index++){
    TFile* inFile = new TFile(a_inFileNames[file_index].c_str(),"READ");
    if(!inFile || !inFile->IsOpen()){
      cout << "ERROR:  EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination.. couldn't find file: " << a_inFileNames[file_index] << endl;
      return;
    }

    if(file_index == 0){
      histo_Vz_structure = (TH1I*) inFile->Get("EfficiencyHistograms/vZ_structure")->Clone();
      //m_embedding_vZ_structure = histo_Vz_structure;
    }

    for(unsigned int vZ_index = 0; vZ_index < numVzBins; vZ_index++){
      if(doRefMultCuts){
        TH3I* embTrackHisto_byRefMult_Int          = (TH3I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_byRefMult_Vz%03d",vZ_index))->Clone();
        TH3I* matchTrackHisto_byRefMult_Int        = (TH3I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_byRefMult_Vz%03d",vZ_index))->Clone();
        if(!embTrackHisto_byRefMult_Int || !matchTrackHisto_byRefMult_Int){
          cout << "ERROR: Couln't read embTrackHisto_byRefMult in " << a_inFileNames[file_index] << "..." << endl;
          inFile->Close();
          continue;
        }
        embTrackHisto_byRefMult[vZ_index][file_index]        = HistogramUtilities::convertInt_to_DoubleHisto3D(embTrackHisto_byRefMult_Int);
        matchTrackHisto_byRefMult[vZ_index][file_index]      = HistogramUtilities::convertInt_to_DoubleHisto3D(matchTrackHisto_byRefMult_Int);
        delete embTrackHisto_byRefMult_Int;
        delete matchTrackHisto_byRefMult_Int;
      }

      for(int centIndex = 0; centIndex < 16; centIndex++){
        TH2I* test           = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex));
        if(!test){
          cout << "WARNING: couldn't read " << Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex) << " in " << a_inFileNames[file_index] << "..." << endl;
          continue;
        }
        TH2I* embeddedTrackHistogram_Int           = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_Cent%d_Vz%03d",centIndex,vZ_index))->Clone();
        TH2I* matchTrackHistogram_Int              = (TH2I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_Cent%d_Vz%03d",centIndex,vZ_index))->Clone();
        TH2I* embeddedTrackHistogram_EtaPt_Int     = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_pTEta_Cent%d_Vz%03d",centIndex,vZ_index))->Clone();
        TH2I* matchTrackHistogram_EtaPt_Int        = (TH2I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_pTEta_Cent%d_Vz%03d",centIndex,vZ_index))->Clone();



        embeddedTrackHistogram[centIndex][vZ_index][file_index]           = HistogramUtilities::convertInt_to_DoubleHisto2D(embeddedTrackHistogram_Int);
        matchTrackHistogram[centIndex][vZ_index][file_index]              = HistogramUtilities::convertInt_to_DoubleHisto2D(matchTrackHistogram_Int);
        embeddedTrackHistogram_EtaPt[centIndex][vZ_index][file_index]     = HistogramUtilities::convertInt_to_DoubleHisto2D(embeddedTrackHistogram_EtaPt_Int);
        matchTrackHistogram_EtaPt[centIndex][vZ_index][file_index]        = HistogramUtilities::convertInt_to_DoubleHisto2D(matchTrackHistogram_EtaPt_Int);

        delete embeddedTrackHistogram_Int;
        delete matchTrackHistogram_Int;
        delete embeddedTrackHistogram_EtaPt_Int;
        delete matchTrackHistogram_EtaPt_Int;

        if(centIndex == 0){
          embTracksAllCent[vZ_index][file_index]         = (TH2D*) embeddedTrackHistogram[centIndex][vZ_index][file_index]->Clone();
          matchTracksAllCent[vZ_index][file_index]       = (TH2D*) matchTrackHistogram[centIndex][vZ_index][file_index]->Clone();
          embTracksAllCent_EtaPt[vZ_index][file_index]   = (TH2D*) embeddedTrackHistogram_EtaPt[centIndex][vZ_index][file_index]->Clone();
          matchTracksAllCent_EtaPt[vZ_index][file_index] = (TH2D*) matchTrackHistogram_EtaPt[centIndex][vZ_index][file_index]->Clone();
        }else{
          embTracksAllCent[vZ_index][file_index]->Add(embeddedTrackHistogram[centIndex][vZ_index][file_index]);
          matchTracksAllCent[vZ_index][file_index]->Add(matchTrackHistogram[centIndex][vZ_index][file_index]);
          embTracksAllCent_EtaPt[vZ_index][file_index]->Add(embeddedTrackHistogram_EtaPt[centIndex][vZ_index][file_index]);
          matchTracksAllCent_EtaPt[vZ_index][file_index]->Add(matchTrackHistogram_EtaPt[centIndex][vZ_index][file_index]);
        }
      }// cent index

      if(file_index == 0){
        embTracksAllCentAllFiles[vZ_index] = (TH2D*) embTracksAllCent[vZ_index][file_index]->Clone();
        matchTracksAllCentAllFiles[vZ_index] = (TH2D*) matchTracksAllCent[vZ_index][file_index]->Clone();
        embTracksAllCent_EtaPtAllFiles[vZ_index] = (TH2D*) embTracksAllCent_EtaPt[vZ_index][file_index]->Clone();
        matchTracksAllCent_EtaPtAllFiles[vZ_index] = (TH2D*) matchTracksAllCent_EtaPt[vZ_index][file_index]->Clone();
      }else{
        embTracksAllCentAllFiles[vZ_index]->Add(embTracksAllCent[vZ_index][file_index]);
        matchTracksAllCentAllFiles[vZ_index]->Add(matchTracksAllCent[vZ_index][file_index]);
        embTracksAllCent_EtaPtAllFiles[vZ_index]->Add(embTracksAllCent_EtaPt[vZ_index][file_index]);
        matchTracksAllCent_EtaPtAllFiles[vZ_index]->Add(matchTracksAllCent_EtaPt[vZ_index][file_index]);
      }
    }// vZ_index loop
  }// file_index loop



  vector< double > vZ_edges;
  for(unsigned int file_index = 0; file_index < a_VzRanges.size(); file_index++){
    for(unsigned int vZ_index = 0; vZ_index < a_VzRanges[file_index].size(); vZ_index++){
      vZ_edges.push_back(a_VzRanges[file_index][vZ_index]);
    }
  }
  std::sort(vZ_edges.begin(), vZ_edges.end());

  double smallest_Vz = vZ_edges[0];
  double largest_Vz = vZ_edges[vZ_edges.size()-1];

  for(double location = -210.0; location < 215.0; location += 10.0){
    if(smallest_Vz < location && location < largest_Vz){
      bool found = false;
      for(unsigned int index = 0; index < vZ_edges.size(); index++){
        if(fabs(vZ_edges[index] - location) < 5.0 ){
          found = true;
        }
      }
      vZ_edges.push_back(location);
    }
  }

  std::sort(vZ_edges.begin(), vZ_edges.end());
  vZ_edges.erase( unique( vZ_edges.begin(), vZ_edges.end() ), vZ_edges.end() );
  const unsigned int numVzRanges = vZ_edges.size() - 1;
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Vz Edges: ";
    for(unsigned int index = 0; index < vZ_edges.size(); index++){
      cout << vZ_edges[index] << "   ";
    }
    cout << endl;
    cout << "    " << numVzRanges << "  Vz Ranges " << endl;
  #endif


  TH1D* histo_Vz_rebin_structure = new TH1D("histo_Vz_rebin_structure","",((int) numVzRanges),&vZ_edges[0]);


  // Combine Histograms by Vz Range and Generate Efficiencies for each range

  TH2D* embTracksAllCent_AllVzRanges         = nullptr;
  TH2D* matchTracksAllCent_AllVzRanges       = nullptr;
  TH2D* embTracksAllCent_EtaPt_AllVzRanges   = nullptr;
  TH2D* matchTracksAllCent_EtaPt_AllVzRanges = nullptr;

  TH2D* embTracksAllCent_VzRanges[numVzRanges];
  TH2D* matchTracksAllCent_VzRanges[numVzRanges];
  TH2D* embTracksAllCent_EtaPt_VzRanges[numVzRanges];
  TH2D* matchTracksAllCent_EtaPt_VzRanges[numVzRanges];

  TH2D* embeddedTrackHistogram_VzRanges[16][numVzRanges];
  TH2D* matchTrackHistogram_VzRanges[16][numVzRanges];
  TH2D* embeddedTrackHistogram_EtaPt_VzRanges[16][numVzRanges];
  TH2D* matchTrackHistogram_EtaPt_VzRanges[16][numVzRanges];

  TH3D* embTrackHisto_byRefMult_VzRanges[numVzRanges];
  TH3D* matchTrackHisto_byRefMult_VzRanges[numVzRanges];

  TGraphAsymmErrors* graph_eff_AllCent_VzRanges[g_EffFit_MaxNumRapBins][numVzRanges];
  TGraphAsymmErrors* graph_eff_EtaPt_AllCent_VzRanges[g_EffFit_MaxNumRapBins][numVzRanges];
  TGraphAsymmErrors* graph_eff_VzRanges[16][g_EffFit_MaxNumRapBins][numVzRanges];
  TGraphAsymmErrors* graph_eff_EtaPt_VzRanges[16][g_EffFit_MaxNumRapBins][numVzRanges];

  HistogramUtilities::nullptr_array_1D(embTracksAllCent_VzRanges, numVzRanges);
  HistogramUtilities::nullptr_array_1D(matchTracksAllCent_VzRanges, numVzRanges);
  HistogramUtilities::nullptr_array_1D(embTracksAllCent_EtaPt_VzRanges, numVzRanges);
  HistogramUtilities::nullptr_array_1D(matchTracksAllCent_EtaPt_VzRanges, numVzRanges);
  HistogramUtilities::nullptr_array_1D(embTrackHisto_byRefMult_VzRanges, numVzRanges);
  HistogramUtilities::nullptr_array_1D(matchTrackHisto_byRefMult_VzRanges, numVzRanges);

  for(unsigned int centIndex = 0; centIndex < 16; centIndex++){
    for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
      embeddedTrackHistogram_VzRanges[centIndex][range_index] = nullptr;
      matchTrackHistogram_VzRanges[centIndex][range_index]    = nullptr;
      embeddedTrackHistogram_EtaPt_VzRanges[centIndex][range_index] = nullptr;
      matchTrackHistogram_EtaPt_VzRanges[centIndex][range_index]    = nullptr;
    }
  }

  for(unsigned int rap_index = 0; rap_index < g_EffFit_MaxNumRapBins; rap_index++){
    for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
      for(unsigned int centIndex = 0; centIndex < 16; centIndex++){
        graph_eff_VzRanges[centIndex][rap_index][range_index] = nullptr;
        graph_eff_EtaPt_VzRanges[centIndex][rap_index][range_index] = nullptr;
      }
      graph_eff_AllCent_VzRanges[rap_index][range_index] = nullptr;
      graph_eff_EtaPt_AllCent_VzRanges[rap_index][range_index] = nullptr;
    }
  }


  for(int vZ_bin = 1; vZ_bin <= histo_Vz_structure->GetNbinsX(); vZ_bin++){
    int vZ_index = vZ_bin - 1;
    double vZ_bin_center = histo_Vz_structure->GetBinCenter(vZ_bin);
    int range_bin = histo_Vz_rebin_structure->FindBin(vZ_bin_center);
    if(range_bin < 1 || range_bin > histo_Vz_rebin_structure->GetNbinsX()) continue;
    int range_index = range_bin - 1;

    for(unsigned int file_index = 0; file_index < numFiles; file_index++){
      if(!embTracksAllCent[vZ_index][file_index]) continue;
      if(embTracksAllCent_VzRanges[range_index]){
        embTracksAllCent_VzRanges[range_index]          ->Add( embTracksAllCent[vZ_index][file_index] );
        matchTracksAllCent_VzRanges[range_index]        ->Add( matchTracksAllCent[vZ_index][file_index] );
        embTracksAllCent_EtaPt_VzRanges[range_index]    ->Add( embTracksAllCent_EtaPt[vZ_index][file_index] );
        matchTracksAllCent_EtaPt_VzRanges[range_index]  ->Add( matchTracksAllCent_EtaPt[vZ_index][file_index] );
        if(doRefMultCuts){
          embTrackHisto_byRefMult_VzRanges[range_index]   ->Add( embTrackHisto_byRefMult[vZ_index][file_index] );
          matchTrackHisto_byRefMult_VzRanges[range_index] ->Add( matchTrackHisto_byRefMult[vZ_index][file_index] );
        }
      }else{
        embTracksAllCent_VzRanges[range_index]          = (TH2D*) embTracksAllCent[vZ_index][file_index]->Clone();
        matchTracksAllCent_VzRanges[range_index]        = (TH2D*) matchTracksAllCent[vZ_index][file_index]->Clone();
        embTracksAllCent_EtaPt_VzRanges[range_index]    = (TH2D*) embTracksAllCent_EtaPt[vZ_index][file_index]->Clone();
        matchTracksAllCent_EtaPt_VzRanges[range_index]  = (TH2D*) matchTracksAllCent_EtaPt[vZ_index][file_index]->Clone();
        if(doRefMultCuts){
          embTrackHisto_byRefMult_VzRanges[range_index]   = (TH3D*) embTrackHisto_byRefMult[vZ_index][file_index]->Clone();
          matchTrackHisto_byRefMult_VzRanges[range_index] = (TH3D*) matchTrackHisto_byRefMult[vZ_index][file_index]->Clone();
        }
      }

      for(unsigned int centIndex = 0; centIndex < 16; centIndex++){
        if(!embeddedTrackHistogram[centIndex][vZ_index][file_index]) continue;
        if(embeddedTrackHistogram_VzRanges[centIndex][range_index]){
          embeddedTrackHistogram_VzRanges[centIndex][range_index]       ->Add( embeddedTrackHistogram[centIndex][vZ_index][file_index] );
          matchTrackHistogram_VzRanges[centIndex][range_index]          ->Add( matchTrackHistogram[centIndex][vZ_index][file_index] );
          embeddedTrackHistogram_EtaPt_VzRanges[centIndex][range_index] ->Add( embeddedTrackHistogram_EtaPt[centIndex][vZ_index][file_index] );
          matchTrackHistogram_EtaPt_VzRanges[centIndex][range_index]    ->Add( matchTrackHistogram_EtaPt[centIndex][vZ_index][file_index] );
        }else{
          embeddedTrackHistogram_VzRanges[centIndex][range_index]       = (TH2D*) embeddedTrackHistogram[centIndex][vZ_index][file_index]->Clone();
          matchTrackHistogram_VzRanges[centIndex][range_index]          = (TH2D*) matchTrackHistogram[centIndex][vZ_index][file_index]->Clone();
          embeddedTrackHistogram_EtaPt_VzRanges[centIndex][range_index] = (TH2D*) embeddedTrackHistogram_EtaPt[centIndex][vZ_index][file_index]->Clone();
          matchTrackHistogram_EtaPt_VzRanges[centIndex][range_index]    = (TH2D*) matchTrackHistogram_EtaPt[centIndex][vZ_index][file_index]->Clone();
        }
      }// cent index loop
    } // file index

    if(!embTracksAllCent_AllVzRanges){
      embTracksAllCent_AllVzRanges = (TH2D*) embTracksAllCent_VzRanges[range_index]->Clone();
      matchTracksAllCent_AllVzRanges = (TH2D*) matchTracksAllCent_VzRanges[range_index]->Clone();
      embTracksAllCent_EtaPt_AllVzRanges = (TH2D*) embTracksAllCent_EtaPt_VzRanges[range_index]->Clone();
      matchTracksAllCent_EtaPt_AllVzRanges = (TH2D*) matchTracksAllCent_EtaPt_VzRanges[range_index]->Clone();
    }else{
      embTracksAllCent_AllVzRanges->Add( embTracksAllCent_VzRanges[range_index] );
      matchTracksAllCent_AllVzRanges->Add( matchTracksAllCent_VzRanges[range_index] );
      embTracksAllCent_EtaPt_AllVzRanges->Add( embTracksAllCent_EtaPt_VzRanges[range_index] );
      matchTracksAllCent_EtaPt_AllVzRanges->Add( matchTracksAllCent_EtaPt_VzRanges[range_index] );
    }
  }// vz bin loop



  TH2D* histos_matched[2] = {matchTracksAllCent_AllVzRanges,matchTracksAllCent_EtaPt_AllVzRanges};
  TH2D* histos_embed[2]   = {embTracksAllCent_AllVzRanges,embTracksAllCent_EtaPt_AllVzRanges};
  for(unsigned int histo_index = 0; histo_index < 2; histo_index++){
    for(int rapBin = 1; rapBin <= histos_embed[histo_index]->GetNbinsX(); rapBin++){
      int rapIndex = rapBin - 1;
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        if(histo_index == 0)  cout << "##### Rapidity = [" << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin) << ", " << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
        else cout << "##### Eta = [" << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin) << ", " << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
      #endif

      TH1D* embTrackHisto   = histos_embed[histo_index]->ProjectionY(Form("embTrack_%02d_%03d",histo_index,rapBin),rapBin,rapBin);
      TH1D* matchTrackHisto = histos_matched[histo_index]->ProjectionY(Form("matchTrack_%02d_%03d",histo_index,rapBin),rapBin,rapBin);
      TH1D* diff_histo = (TH1D*) embTrackHisto->Clone();
      for(int bin = 1; bin <= diff_histo->GetNbinsX(); bin++){
        double content = embTrackHisto->GetBinContent(bin) - matchTrackHisto->GetBinContent(bin);
        double error = sqrt(content);
        if(!std::isfinite(error) || error < 0.0) error = 0.0;
        diff_histo->SetBinContent(bin,content);
        diff_histo->SetBinError(bin,error);
      }

      int firstBin = HistogramUtilities::firstBinWithContent(matchTrackHisto);
      int lastBin  = HistogramUtilities::lastBinWithContent(matchTrackHisto);
      double lowEdgeFirstBin = matchTrackHisto->GetBinLowEdge(firstBin);
      double highEdgeLastBin = matchTrackHisto->GetBinLowEdge(lastBin+1);
      if(firstBin > lastBin) continue; // it's empty
      TH1D* truncDiffHisto = HistogramUtilities::truncateHistogram(diff_histo,lowEdgeFirstBin,highEdgeLastBin);
      delete diff_histo;
      TH1D* histo_rebin_structure  = HistogramUtilities::rebinEqualByCDF(truncDiffHisto, m_tpc_eff_no_fit_num_bins,false,false); // rebins so more bins where content is changing, less bins where there is less change
      HistogramUtilities::zeroBins(histo_rebin_structure);
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Diff Histo Truncated mT-m0 = [" << lowEdgeFirstBin << ", " << highEdgeLastBin << "]      struct = [" << HistogramUtilities::getLowEdge(histo_rebin_structure) << ", " << HistogramUtilities::getHighEdge(histo_rebin_structure) << "]";
        cout << "         Match: " << HistogramUtilities::sumOfContent(matchTrackHisto) << "    Embed: " << HistogramUtilities::sumOfContent(embTrackHisto) << endl;
      #endif
      delete truncDiffHisto;
      delete embTrackHisto;
      delete matchTrackHisto;

      for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){

        // IF USING REFMULT ADDITIONS, OVERWRITE PREVIOUS EFF HISTOGRAMS
        if(doRefMultCuts && histo_index == 0){
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "Doing refMult slice algorithm... " << endl;
          #endif
          vector< TH2D* > embTrackHisto_byRefMult_Vector;
          vector< TH2D* > matchTrackHisto_byRefMult_Vector;
          embTrackHisto_byRefMult_Vector   = HistogramUtilities::projeciton2DSlices_Z(embTrackHisto_byRefMult_VzRanges[range_index],  a_refMultCuts);
          matchTrackHisto_byRefMult_Vector = HistogramUtilities::projeciton2DSlices_Z(matchTrackHisto_byRefMult_VzRanges[range_index],a_refMultCuts);
          for(unsigned int centIndex = 0; centIndex < embTrackHisto_byRefMult_Vector.size(); centIndex++){
            HistogramUtilities::ConditionalDelete(embeddedTrackHistogram_VzRanges[centIndex][range_index]);
            HistogramUtilities::ConditionalDelete(matchTrackHistogram_VzRanges[centIndex][range_index]);
            embeddedTrackHistogram_VzRanges[centIndex][range_index] = embTrackHisto_byRefMult_Vector[centIndex];
            matchTrackHistogram_VzRanges[centIndex][range_index]    = matchTrackHisto_byRefMult_Vector[centIndex];
          }// cent loop
        }


        for(int centIndex = 0; centIndex < 16; centIndex++){
          if(histo_index == 0){
            if(!embeddedTrackHistogram_VzRanges[centIndex][range_index]) continue;
            embTrackHisto   = embeddedTrackHistogram_VzRanges[centIndex][range_index]->ProjectionY(Form("embTrack_%d_%d_%d_%d",histo_index,range_index,rapBin,centIndex),rapBin,rapBin);
            matchTrackHisto = matchTrackHistogram_VzRanges[centIndex][range_index]->ProjectionY(Form("matchTrack_%d_%d_%d_%d",histo_index,range_index,rapBin,centIndex),rapBin,rapBin);
          }else{
            if(!embeddedTrackHistogram_EtaPt_VzRanges[centIndex][range_index]) continue;
            embTrackHisto   = embeddedTrackHistogram_EtaPt_VzRanges[centIndex][range_index]->ProjectionY(Form("embTrack_%d_%d_%d_%d",histo_index,range_index,rapBin,centIndex),rapBin,rapBin);
            matchTrackHisto = matchTrackHistogram_EtaPt_VzRanges[centIndex][range_index]->ProjectionY(Form("matchTrack_%d_%d_%d_%d",histo_index,range_index,rapBin,centIndex),rapBin,rapBin);
          }
          if(HistogramUtilities::sumOfContent(matchTrackHisto) <= 0.0) continue;

          TH1D* rebinMatchTrack = (TH1D*) histo_rebin_structure->Clone();
          TH1D* rebinEmbTrack   = (TH1D*) histo_rebin_structure->Clone();
          HistogramUtilities::rebinHistoToAlternateBinStructure(matchTrackHisto,rebinMatchTrack);
          HistogramUtilities::rebinHistoToAlternateBinStructure(embTrackHisto,rebinEmbTrack);
          delete embTrackHisto;
          delete matchTrackHisto;
          HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
          HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);

          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "Cent " << centIndex << " Match Track Histo  mT-m0 = [" << HistogramUtilities::getLowEdge(rebinMatchTrack) << ", " << HistogramUtilities::getHighEdge(rebinMatchTrack) << "]";
            cout << "         Match: " << HistogramUtilities::sumOfContent(rebinMatchTrack) << "    Embed: " << HistogramUtilities::sumOfContent(rebinEmbTrack) << endl;
          #endif

          if(histo_index == 0){
            graph_eff_VzRanges[centIndex][rapIndex][range_index] = new TGraphAsymmErrors();
            graph_eff_VzRanges[centIndex][rapIndex][range_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
          }else{
            graph_eff_EtaPt_VzRanges[centIndex][rapIndex][range_index] = new TGraphAsymmErrors();
            graph_eff_EtaPt_VzRanges[centIndex][rapIndex][range_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
          }
          HistogramUtilities::ConditionalDelete(rebinMatchTrack);
          HistogramUtilities::ConditionalDelete(rebinEmbTrack);
        }// cent loop




        // CENTRALITY INTEGRATED ------------------------------------
        if(histo_index == 0){
          embTrackHisto   = embTracksAllCent_VzRanges[range_index]->ProjectionY(Form("AllCent_embTrack_%d_%d_%d",histo_index,range_index,rapBin),rapBin,rapBin);
          matchTrackHisto = matchTracksAllCent_VzRanges[range_index]->ProjectionY(Form("AllCent_matchTrack_%d_%d_%d",histo_index,range_index,rapBin),rapBin,rapBin);
        }else{
          embTrackHisto   = embTracksAllCent_EtaPt_VzRanges[range_index]->ProjectionY(Form("AllCent_embTrack_%d_%d_%d",histo_index,range_index,rapBin),rapBin,rapBin);
          matchTrackHisto = matchTracksAllCent_EtaPt_VzRanges[range_index]->ProjectionY(Form("AllCent_matchTrack_%d_%d_%d",histo_index,range_index,rapBin),rapBin,rapBin);
        }

        TH1D* rebinMatchTrack = (TH1D*) histo_rebin_structure->Clone();
        TH1D* rebinEmbTrack   = (TH1D*) histo_rebin_structure->Clone();
        HistogramUtilities::rebinHistoToAlternateBinStructure(matchTrackHisto,rebinMatchTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(embTrackHisto,rebinEmbTrack);
        delete embTrackHisto;
        delete matchTrackHisto;
        HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
        HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);

        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " All Cent Match Track Histo  mT-m0 = [" << HistogramUtilities::getLowEdge(rebinMatchTrack) << ", " << HistogramUtilities::getHighEdge(rebinMatchTrack) << "]" << endl;
        #endif

        if(histo_index == 0){
          graph_eff_AllCent_VzRanges[rapIndex][range_index] = new TGraphAsymmErrors();
          graph_eff_AllCent_VzRanges[rapIndex][range_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
        }else{
          graph_eff_EtaPt_AllCent_VzRanges[rapIndex][range_index] = new TGraphAsymmErrors();
          graph_eff_EtaPt_AllCent_VzRanges[rapIndex][range_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
        }
        HistogramUtilities::ConditionalDelete(rebinMatchTrack);
        HistogramUtilities::ConditionalDelete(rebinEmbTrack);

      } // range index loop
    }// rapidity loop
  }// histo loop




  // Calcluate Vz Weights For Weighting Algorithm  ------------------------------------------------
  double weights[numVzRanges];
  double total_weights = 0.0;
  double data_sumOfContent = HistogramUtilities::sumOfContent(histo_Vz_data,-1,false);
  for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
    int range_bin = range_index + 1;
    double lowVzRange = histo_Vz_rebin_structure->GetBinLowEdge(range_bin);
    double highVzRange = histo_Vz_rebin_structure->GetBinLowEdge(range_bin+1);
    weights[range_index] = HistogramUtilities::sumOfContent(histo_Vz_data,-1,false,lowVzRange,highVzRange);
    total_weights += weights[range_index];
  }
  cout << "Vz Weights:" << endl;
  for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
    cout << "     Range Index: " << range_index << "   weight: " << weights[range_index] << "Vz=[" << histo_Vz_rebin_structure->GetBinLowEdge(range_index+1) << ", "
        << histo_Vz_rebin_structure->GetBinLowEdge(range_index+2) << "]" << endl;
    weights[range_index] = weights[range_index]/total_weights;
    cout << "  fractional weight: " << weights[range_index] << endl;

  }




  //##################   REWEIGHT THE DIFFERENT EFFICIENCIES INTO ONE  ###########################
  for(unsigned int graph_index = 0; graph_index < 2; graph_index++){
    int numBins = 0;
    if(graph_index == 0)  numBins =  embTracksAllCent_AllVzRanges->GetNbinsX();
    else                  numBins =  embTracksAllCent_EtaPt_AllVzRanges->GetNbinsX();

    for(int rapBin = 1; rapBin <= numBins; rapBin++){
      int rapIndex = rapBin - 1;
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        if(graph_index == 0)  cout << "##### Rapidity = [" << embTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin) << ", " << embTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
        else cout << "##### Eta = [" << embTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin) << ", " << embTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
      #endif


      vector< double > mTm0_full_array;
      vector< double > mTm0_ErrLow_full_array;
      vector< double > mTm0_ErrHigh_full_array;
      for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
        TGraphAsymmErrors* effGraph_Base = nullptr;
        if(graph_index == 0)  effGraph_Base = graph_eff_AllCent_VzRanges[rapIndex][range_index];
        else                  effGraph_Base = graph_eff_EtaPt_AllCent_VzRanges[rapIndex][range_index];
        if(!effGraph_Base) continue;
        for(int point_index = 0; point_index < effGraph_Base->GetN(); point_index++){
          double mTm0 = effGraph_Base->GetX()[point_index];
          double found = false;
          for(unsigned int index = 0; index < mTm0_full_array.size(); index++){
            if(fabs(mTm0 - mTm0_full_array[index]) < 0.001){
              found = true;
              break;
            }
          }
          if(!found){
            mTm0_full_array.push_back(mTm0);
            mTm0_ErrLow_full_array.push_back(effGraph_Base->GetEXlow()[point_index]);
            mTm0_ErrHigh_full_array.push_back(effGraph_Base->GetEXhigh()[point_index]);
          }
        }
      }

      if(mTm0_full_array.size() == 0){
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "mTm0_full_array is empty: " << mTm0_full_array.size() << endl;
        #endif
        continue;
      }
      for(int centIndex = 0; centIndex < 17; centIndex++){
        TGraphAsymmErrors* graph_eff_combined = nullptr;

        if(centIndex == 16){ // Centrality Integrated
          if(graph_index == 0){
            m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
            m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality Integrated y=[%1.2f,%1.2f]; m_{T}-m_{0}; TPC Efficiency",
                       m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),matchTracksAllCent_AllVzRanges ->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin+1)));
            graph_eff_combined = m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex];
          }else{
            m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
            m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality Integrated #eta=[%1.2f,%1.2f]; p_{T}; TPC Efficiency",
                       m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),matchTracksAllCent_EtaPt_AllVzRanges ->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin+1)));
            graph_eff_combined = m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex];
          }

        }else{ // Centrality Specific
          /*if(!embeddedTrackHistogram[centIndex][0]){
            #ifdef _EFFICIENCY_FITTER_DEBUG_
              cout << "embeddedTrackHistogram graph not made, ptr: " << embeddedTrackHistogram[centIndex][0] << endl;
            #endif
            continue;
          }*/
          if(graph_index == 0){
            //cout << "ptr to delete: " << m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] << endl;
            //HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]);
            m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
            m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality %02d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
             matchTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin+1)));
            graph_eff_combined = m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex];
          }else{
            //HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]);
            m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_Cent%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
            m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality %d Integrated #eta=[%1.2f,%1.2f]; p_{T}; TPC Efficiency",
                 m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,matchTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rapBin+1)));
            graph_eff_combined = m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex];
          }

          #ifdef _EFFICIENCY_FITTER_DEBUG_
            if(graph_index == 0 && rapBin == 12 && centIndex == 0){
              cout << "###### POINT BY POINT CHECK ----------------------" << endl;
              for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
                cout << "EFF GRAPH ---------------------- Range: " << range_index  << endl;
                HistogramUtilities::printGraphInfo(graph_eff_VzRanges[centIndex][rapIndex][range_index]);
              }
            }
          #endif
        }

        int numPoints = mTm0_full_array.size();
        for(unsigned int point_index = 0; point_index < numPoints; point_index++){ //-------------   point index loop  ------------
          double sumOfWeights = 0.0;
          double efficiency_value      = 0.0;
          double efficiency_error_low  = 0.0;
          double efficiency_error_high = 0.0;
          double mTm0 = mTm0_full_array[point_index];

          if(graph_index == 0 && rapBin == 12 && centIndex == 0) cout << "Constructing Point:  mtm0: " << mTm0 << endl;


          //bool atLeastOneFileIsFullCoverage = false;
          for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
            TGraphAsymmErrors* effGraph = nullptr;
            if(centIndex == 16){ //Centrality Integrated
              if(graph_index == 0) effGraph = graph_eff_AllCent_VzRanges[rapIndex][range_index];
              else                 effGraph = graph_eff_EtaPt_AllCent_VzRanges[rapIndex][range_index];
            }else{ // Centrality Specific
              if(graph_index == 0) effGraph = graph_eff_VzRanges[centIndex][rapIndex][range_index];
              else                 effGraph = graph_eff_EtaPt_VzRanges[centIndex][rapIndex][range_index];
            }

            if(!effGraph || effGraph->GetN() <= 0){
              #ifdef _EFFICIENCY_FITTER_DEBUG_
                if(effGraph) cout << "effGraph graph has no points, NumPoints: " << effGraph->GetN() << endl;
                else cout << "effGraph graph not made, ptr: " << effGraph << endl;
              #endif
              continue;
            }
            int point_searched = HistogramUtilities::getPointClosestToVal(effGraph,mTm0);
            if(fabs(effGraph->GetX()[point_searched] - mTm0) > 0.001){
              cout << "WARNING: Didn't find the right corresponding point in Vz Recombination" << endl;
              continue; // didn't find the right point
            }
            //if(isFullVzRange[range_index]) atLeastOneFileIsFullCoverage = true;
            sumOfWeights          += weights[range_index];
            efficiency_value      += weights[range_index]*effGraph->GetY()[point_searched];
            efficiency_error_low  += pow(weights[range_index]*effGraph->GetEYlow()[point_searched],2);
            efficiency_error_high += pow(weights[range_index]*effGraph->GetEYhigh()[point_searched],2);
            if(graph_index == 0 && rapBin == 12 && centIndex == 0) cout << "    Point File " << range_index <<  ":  mtm0: " << effGraph->GetX()[point_searched] << "  eff: " << effGraph->GetY()[point_searched] << endl;

          }// file index loop
          efficiency_value      = efficiency_value/sumOfWeights;
          efficiency_error_low  = sqrt(efficiency_error_low)/sumOfWeights;
          efficiency_error_high = sqrt(efficiency_error_high)/sumOfWeights;

          if(
            //!atLeastOneFileIsFullCoverage ||
            sumOfWeights == 0.0 || !std::isfinite(efficiency_value) || !std::isfinite(efficiency_error_low) || !std::isfinite(efficiency_error_high)){
            cout << "WARNING: Vz Recomb ->  mTm0: " << mTm0 << "  eff:" << efficiency_value << "  el:" << efficiency_error_low << "   eh:" << efficiency_error_high  << endl;
            continue;
          }


          int pointIndexToAdd = graph_eff_combined->GetN();
          graph_eff_combined->SetPoint(pointIndexToAdd,mTm0,efficiency_value);
          graph_eff_combined->SetPointError(pointIndexToAdd,mTm0_ErrLow_full_array[point_index],mTm0_ErrHigh_full_array[point_index],efficiency_error_low,efficiency_error_high);
          if(graph_index == 0 && rapBin == 12 && centIndex == 0) cout << "Adding Point:  mtm0: " << mTm0 << "  eff: " << efficiency_value << endl;
        }// point index loop
        graph_eff_combined->Sort();
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          if(graph_index == 0 && rapBin == 12 && centIndex == 0){
            cout << "RESULTING EFF GRAPH ----------------------" << endl;
            HistogramUtilities::printGraphInfo(m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]);
          }
        #endif


      } // cent index loop
    }// rapidity loop
  }// graph loop


  // Testing outputs to root file for debugging
  TFile* testFile = new TFile(Form("/home/matthewharasty/ResearchData/2024_01_17_7p7GeV_Embedding/test_part%02d_pm%02d.root",a_partIndex,pmIndex),"RECREATE");
  for(unsigned int rap_index = 0; rap_index < g_EffFit_MaxNumRapBins; rap_index++){
    for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
      if(graph_eff_AllCent_VzRanges[rap_index][range_index]){
        graph_eff_AllCent_VzRanges[rap_index][range_index]->SetName(Form("graph_eff_AllCent_VzRanges_Rap%03d_Vz%02d",rap_index,range_index));
        graph_eff_AllCent_VzRanges[rap_index][range_index]->SetTitle(
            Form("graph_eff_AllCent_VzRanges y=[%f,%f] Vz=[%f,%f]",
              embTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+1),
              embTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+2),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+1),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+2)
             ));
      }
      if(graph_eff_EtaPt_AllCent_VzRanges[rap_index][range_index]){
        graph_eff_EtaPt_AllCent_VzRanges[rap_index][range_index]->SetName(Form("graph_eff_EtaPt_AllCent_VzRanges_Rap%03d_Vz%02d",rap_index,range_index));
        graph_eff_EtaPt_AllCent_VzRanges[rap_index][range_index]->SetTitle(
          Form("graph_eff_EtaPt_AllCent_VzRanges y=[%f,%f] Vz=[%f,%f]",
              embTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+1),
              embTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+2),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+1),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+2)
             ));
      }
    }

    for(unsigned int centIndex = 0; centIndex < 16; centIndex++){
      for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
        if(graph_eff_VzRanges[centIndex][rap_index][range_index]){
          graph_eff_VzRanges[centIndex][rap_index][range_index]->SetName(Form("graph_eff_VzRanges_Cent%02d_Rap%03d_Vz%02d",centIndex,rap_index,range_index));
          graph_eff_VzRanges[centIndex][rap_index][range_index]->SetTitle(
            Form("graph_eff_VzRanges Cent %02d y=[%f,%f] Vz=[%f,%f]", centIndex,
              embTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+1),
              embTracksAllCent_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+2),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+1),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+2)
             ));
        }
        if(graph_eff_EtaPt_VzRanges[centIndex][rap_index][range_index]){
          graph_eff_EtaPt_VzRanges[centIndex][rap_index][range_index]->SetName(Form("graph_eff_EtaPt_VzRanges_Cent%02d_Rap%03d_Vz%02d",centIndex,rap_index,range_index));
          graph_eff_EtaPt_VzRanges[centIndex][rap_index][range_index]->SetTitle(
            Form("graph_eff_EtaPt_VzRanges Cent %02d y=[%f,%f] Vz=[%f,%f]", centIndex,
              embTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+1),
              embTracksAllCent_EtaPt_AllVzRanges->GetXaxis()->GetBinLowEdge(rap_index+2),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+1),
              histo_Vz_rebin_structure->GetBinLowEdge(range_index+2)
             ));
        }
      }
    }
  }

  for(unsigned int rap_index = 0; rap_index < g_EffFit_MaxNumRapBins; rap_index++){
    for(unsigned int centIndex = 0; centIndex < 16; centIndex++){
      for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
        HistogramUtilities::ConditionalWrite(graph_eff_VzRanges[centIndex][rap_index][range_index]);
        HistogramUtilities::ConditionalWrite(graph_eff_EtaPt_VzRanges[centIndex][rap_index][range_index]);
      }
    }
    for(unsigned int range_index = 0; range_index < numVzRanges; range_index++){
      HistogramUtilities::ConditionalWrite(graph_eff_AllCent_VzRanges[rap_index][range_index]);
      HistogramUtilities::ConditionalWrite(graph_eff_EtaPt_AllCent_VzRanges[rap_index][range_index]);
    }
  }






  inFile_PicoBinner->Close();

  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "End Call to EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination" << endl;
  #endif
}










/*  working with old incorrect method
void EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination(int a_partIndex, int a_charge, vector<string> a_inFileNames,vector<vector<double>> a_VzRanges, string a_picobinner_file_for_Vz_data_histo, vector<double> a_refMultCuts){
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Call to EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination" << endl;
  #endif


  bool doRefMultCuts = (a_refMultCuts.size() > 0);
  int pmIndex = ((a_charge > 0) ? 0 : 1);

  // delete everything to make sure this function works
  for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
    HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]);
    HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex]);
    m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] = nullptr;
    m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex] = nullptr;
    for(int centIndex = 0; centIndex < 16; centIndex++){
      HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]);
      HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]);
      m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = nullptr;
      m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex] = nullptr;
    }
  }



  TFile* inFile_PicoBinner = new TFile(a_picobinner_file_for_Vz_data_histo.c_str(),"READ");
  if(!inFile_PicoBinner || !inFile_PicoBinner->IsOpen()){
    cout << "ERROR:  EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination.. couldn't find file: " << a_picobinner_file_for_Vz_data_histo << endl;
    return;
  }
  TH1I* histo_Vz_data_Int = (TH1I*) inFile_PicoBinner->Get("zVertex")->Clone();
  if(!histo_Vz_data_Int){
    cout << "Couldn't read Z Vertex Histo! from " << a_picobinner_file_for_Vz_data_histo << endl;
  }

  //inFile_PicoBinner->Close();
  TH1D* histo_Vz_data = HistogramUtilities::convertInt_to_DoubleHisto(histo_Vz_data_Int);



  const unsigned int numFiles = a_inFileNames.size();

  TH2D* embTracksAllCentAllFiles = nullptr;
  TH2D* matchTracksAllCentAllFiles = nullptr;
  TH2D* embTracksAllCent_EtaPtAllFiles = nullptr;
  TH2D* matchTracksAllCent_EtaPtAllFiles = nullptr;

  TH2D* embTracksAllCent[numFiles];
  TH2D* matchTracksAllCent[numFiles];
  TH2D* embTracksAllCent_EtaPt[numFiles];
  TH2D* matchTracksAllCent_EtaPt[numFiles];

  TH2D* embeddedTrackHistogram[16][numFiles];
  TH2D* matchTrackHistogram[16][numFiles];
  TH2D* embeddedTrackHistogram_EtaPt[16][numFiles];
  TH2D* matchTrackHistogram_EtaPt[16][numFiles];

  TH3D* embTrackHisto_byRefMult[numFiles];
  TH3D* matchTrackHisto_byRefMult[numFiles];

  TGraphAsymmErrors* graph_eff_AllCent[g_EffFit_MaxNumRapBins][numFiles];
  TGraphAsymmErrors* graph_eff_EtaPt_AllCent[g_EffFit_MaxNumRapBins][numFiles];
  TGraphAsymmErrors* graph_eff[16][g_EffFit_MaxNumRapBins][numFiles];
  TGraphAsymmErrors* graph_eff_EtaPt[16][g_EffFit_MaxNumRapBins][numFiles];

  unsigned int numFiles_int = a_inFileNames.size();
  unsigned int numRapBins_int = (unsigned int) g_EffFit_MaxNumRapBins;

  HistogramUtilities::nullptr_array_1D(embTracksAllCent,          numFiles_int);
  HistogramUtilities::nullptr_array_1D(matchTracksAllCent,        numFiles_int);
  HistogramUtilities::nullptr_array_1D(embTracksAllCent_EtaPt,    numFiles_int);
  HistogramUtilities::nullptr_array_1D(matchTracksAllCent_EtaPt,  numFiles_int);
  HistogramUtilities::nullptr_array_1D(embTrackHisto_byRefMult,   numFiles_int);
  HistogramUtilities::nullptr_array_1D(matchTrackHisto_byRefMult, numFiles_int);

  for(unsigned int iii = 0; iii < 16; iii++){
    for(unsigned int jjj = 0; jjj < numFiles; jjj++){
      embeddedTrackHistogram[iii][jjj] = nullptr;
      matchTrackHistogram[iii][jjj] = nullptr;
      embeddedTrackHistogram_EtaPt[iii][jjj] = nullptr;
      matchTrackHistogram_EtaPt[iii][jjj] = nullptr;
    }
  }

  for(unsigned int jjj = 0; jjj < g_EffFit_MaxNumRapBins; jjj++){
    for(unsigned int kkk = 0; kkk < numFiles; kkk++){
      for(unsigned int iii = 0; iii < 16; iii++){
        graph_eff[iii][jjj][kkk] = nullptr;
        graph_eff_EtaPt[iii][jjj][kkk] = nullptr;
      }
      graph_eff_AllCent[jjj][kkk] = nullptr;
      graph_eff_EtaPt_AllCent[jjj][kkk] = nullptr;
    }
  }



  //################   READ IN THE DATA  ##################################
  for(unsigned int file_index = 0; file_index < a_inFileNames.size(); file_index++){
    TFile* inFile = new TFile(a_inFileNames[file_index].c_str(),"READ");
    if(!inFile || !inFile->IsOpen()){
      cout << "ERROR:  EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination.. couldn't find file: " << a_inFileNames[file_index] << endl;
      return;
    }

    TH3I* embTrackHisto_byRefMult_Int          = (TH3I*) inFile->Get("EfficiencyHistograms/embTrackHisto_byRefMult")->Clone();
    TH3I* matchTrackHisto_byRefMult_Int        = (TH3I*) inFile->Get("EfficiencyHistograms/matchTrackHisto_byRefMult")->Clone();
    if(!embTrackHisto_byRefMult_Int || !matchTrackHisto_byRefMult_Int){
      cout << "ERROR: Couln't read embTrackHisto_byRefMult in " << a_inFileNames[file_index] << "..." << endl;
      inFile->Close();
      continue;
    }
    embTrackHisto_byRefMult[file_index]        = HistogramUtilities::convertInt_to_DoubleHisto3D(embTrackHisto_byRefMult_Int);
    matchTrackHisto_byRefMult[file_index]      = HistogramUtilities::convertInt_to_DoubleHisto3D(matchTrackHisto_byRefMult_Int);
    delete embTrackHisto_byRefMult_Int;
    delete matchTrackHisto_byRefMult_Int;

    for(int centIndex = 0; centIndex < 16; centIndex++){
      TH2I* test           = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex));
      if(!test){
        cout << "WARNING: couldn't read " << Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex) << " in " << a_inFileNames[file_index] << "..." << endl;
        continue;
      }
      TH2I* embeddedTrackHistogram_Int           = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_Cent%d",centIndex))->Clone();
      TH2I* matchTrackHistogram_Int              = (TH2I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_Cent%d",centIndex))->Clone();
      TH2I* embeddedTrackHistogram_EtaPt_Int     = (TH2I*) inFile->Get(Form("EfficiencyHistograms/embTrackHisto_pTEta_Cent%d",centIndex))->Clone();
      TH2I* matchTrackHistogram_EtaPt_Int        = (TH2I*) inFile->Get(Form("EfficiencyHistograms/matchTrackHisto_pTEta_Cent%d",centIndex))->Clone();



      embeddedTrackHistogram[centIndex][file_index]           = HistogramUtilities::convertInt_to_DoubleHisto2D(embeddedTrackHistogram_Int);
      matchTrackHistogram[centIndex][file_index]              = HistogramUtilities::convertInt_to_DoubleHisto2D(matchTrackHistogram_Int);
      embeddedTrackHistogram_EtaPt[centIndex][file_index]     = HistogramUtilities::convertInt_to_DoubleHisto2D(embeddedTrackHistogram_EtaPt_Int);
      matchTrackHistogram_EtaPt[centIndex][file_index]        = HistogramUtilities::convertInt_to_DoubleHisto2D(matchTrackHistogram_EtaPt_Int);

      delete embeddedTrackHistogram_Int;
      delete matchTrackHistogram_Int;
      delete embeddedTrackHistogram_EtaPt_Int;
      delete matchTrackHistogram_EtaPt_Int;

      if(centIndex == 0){
        embTracksAllCent[file_index]         = (TH2D*) embeddedTrackHistogram[centIndex][file_index]->Clone();
        matchTracksAllCent[file_index]       = (TH2D*) matchTrackHistogram[centIndex][file_index]->Clone();
        embTracksAllCent_EtaPt[file_index]   = (TH2D*) embeddedTrackHistogram_EtaPt[centIndex][file_index]->Clone();
        matchTracksAllCent_EtaPt[file_index] = (TH2D*) matchTrackHistogram_EtaPt[centIndex][file_index]->Clone();
      }else{
        embTracksAllCent[file_index]->Add(embeddedTrackHistogram[centIndex][file_index]);
        matchTracksAllCent[file_index]->Add(matchTrackHistogram[centIndex][file_index]);
        embTracksAllCent_EtaPt[file_index]->Add(embeddedTrackHistogram_EtaPt[centIndex][file_index]);
        matchTracksAllCent_EtaPt[file_index]->Add(matchTrackHistogram_EtaPt[centIndex][file_index]);
      }
    }// cent index

    if(file_index == 0){
      embTracksAllCentAllFiles = (TH2D*) embTracksAllCent[file_index]->Clone();
      matchTracksAllCentAllFiles = (TH2D*) matchTracksAllCent[file_index]->Clone();
      embTracksAllCent_EtaPtAllFiles = (TH2D*) embTracksAllCent_EtaPt[file_index]->Clone();
      matchTracksAllCent_EtaPtAllFiles = (TH2D*) matchTracksAllCent_EtaPt[file_index]->Clone();
    }else{
      embTracksAllCentAllFiles->Add(embTracksAllCent[file_index]);
      matchTracksAllCentAllFiles->Add(matchTracksAllCent[file_index]);
      embTracksAllCent_EtaPtAllFiles->Add(embTracksAllCent_EtaPt[file_index]);
      matchTracksAllCent_EtaPtAllFiles->Add(matchTracksAllCent_EtaPt[file_index]);
    }
  }// file_index loop




  // FROM ALLFILE HISTOS, DETERMINE BEST BINNING, CREATE EFFICIENCY GRAPHS
  TH2D* histos_matched[2] = {matchTracksAllCentAllFiles,matchTracksAllCent_EtaPtAllFiles};
  TH2D* histos_embed[2]   = {embTracksAllCentAllFiles,embTracksAllCent_EtaPtAllFiles};
  for(unsigned int histo_index = 0; histo_index < 2; histo_index++){
    for(int rapBin = 1; rapBin <= histos_embed[histo_index]->GetNbinsX(); rapBin++){
      int rapIndex = rapBin - 1;
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        if(histo_index == 0)  cout << "##### Rapidity = [" << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin) << ", " << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
        else cout << "##### Eta = [" << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin) << ", " << histos_embed[histo_index]->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
      #endif

      TH1D* embTrackHisto   = histos_embed[histo_index]->ProjectionY(Form("embTrack_%02d_%03d",histo_index,rapBin),rapBin,rapBin);
      TH1D* matchTrackHisto = histos_matched[histo_index]->ProjectionY(Form("matchTrack_%02d_%03d",histo_index,rapBin),rapBin,rapBin);
      TH1D* diff_histo = (TH1D*) embTrackHisto->Clone();
      for(int bin = 1; bin <= diff_histo->GetNbinsX(); bin++){
        double content = embTrackHisto->GetBinContent(bin) - matchTrackHisto->GetBinContent(bin);
        double error = sqrt(content);
        if(!std::isfinite(error) || error < 0.0) error = 0.0;
        diff_histo->SetBinContent(bin,content);
        diff_histo->SetBinError(bin,error);
      }

      int firstBin = HistogramUtilities::firstBinWithContent(matchTrackHisto);
      int lastBin  = HistogramUtilities::lastBinWithContent(matchTrackHisto);
      double lowEdgeFirstBin = matchTrackHisto->GetBinLowEdge(firstBin);
      double highEdgeLastBin = matchTrackHisto->GetBinLowEdge(lastBin+1);
      if(firstBin > lastBin) continue; // it's empty
      TH1D* truncDiffHisto = HistogramUtilities::truncateHistogram(diff_histo,lowEdgeFirstBin,highEdgeLastBin);
      delete diff_histo;
      TH1D* histo_rebin_structure  = HistogramUtilities::rebinEqualByCDF(truncDiffHisto, m_tpc_eff_no_fit_num_bins,false,false); // rebins so more bins where content is changing, less bins where there is less change
      HistogramUtilities::zeroBins(histo_rebin_structure);
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Diff Histo Truncated mT-m0 = [" << lowEdgeFirstBin << ", " << highEdgeLastBin << "]      struct = [" << HistogramUtilities::getLowEdge(histo_rebin_structure) << ", " << HistogramUtilities::getHighEdge(histo_rebin_structure) << "]";
        cout << "         Match: " << HistogramUtilities::sumOfContent(matchTrackHisto) << "    Embed: " << HistogramUtilities::sumOfContent(embTrackHisto) << endl;
      #endif
      delete truncDiffHisto;
      delete embTrackHisto;
      delete matchTrackHisto;

      for(unsigned int file_index = 0; file_index < a_inFileNames.size(); file_index++){
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "File Index: " << file_index << endl;
        #endif
        // IF USING REFMULT ADDITIONS, OVERWRITE PREVIOUS EFF HISTOGRAMS
        if(doRefMultCuts && histo_index == 0){
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "Doing refMult slice algorithm... " << endl;
          #endif
          vector< TH2D* > embTrackHisto_byRefMult_Vector;
          vector< TH2D* > matchTrackHisto_byRefMult_Vector;
          embTrackHisto_byRefMult_Vector   = HistogramUtilities::projeciton2DSlices_Z(embTrackHisto_byRefMult[file_index],  a_refMultCuts);
          matchTrackHisto_byRefMult_Vector = HistogramUtilities::projeciton2DSlices_Z(matchTrackHisto_byRefMult[file_index],a_refMultCuts);
          for(unsigned int centIndex = 0; centIndex < embTrackHisto_byRefMult_Vector.size(); centIndex++){
            HistogramUtilities::ConditionalDelete(embeddedTrackHistogram[centIndex][file_index]);
            HistogramUtilities::ConditionalDelete(matchTrackHistogram[centIndex][file_index]);
            embeddedTrackHistogram[centIndex][file_index] = embTrackHisto_byRefMult_Vector[centIndex];
            matchTrackHistogram[centIndex][file_index]    = matchTrackHisto_byRefMult_Vector[centIndex];
          }// cent loop
        }


        for(int centIndex = 0; centIndex < 16; centIndex++){
          if(histo_index == 0){
            if(!embeddedTrackHistogram[centIndex][file_index]) continue;
            embTrackHisto   = embeddedTrackHistogram[centIndex][file_index]->ProjectionY(Form("embTrack_%d_%d_%d_%d",histo_index,file_index,rapBin,centIndex),rapBin,rapBin);
            matchTrackHisto = matchTrackHistogram[centIndex][file_index]->ProjectionY(Form("matchTrack_%d_%d_%d_%d",histo_index,file_index,rapBin,centIndex),rapBin,rapBin);
          }else{
            if(!embeddedTrackHistogram_EtaPt[centIndex][file_index]) continue;
            embTrackHisto   = embeddedTrackHistogram_EtaPt[centIndex][file_index]->ProjectionY(Form("embTrack_%d_%d_%d_%d",histo_index,file_index,rapBin,centIndex),rapBin,rapBin);
            matchTrackHisto = matchTrackHistogram_EtaPt[centIndex][file_index]->ProjectionY(Form("matchTrack_%d_%d_%d_%d",histo_index,file_index,rapBin,centIndex),rapBin,rapBin);
          }
          if(HistogramUtilities::sumOfContent(matchTrackHisto) <= 0.0) continue;

          TH1D* rebinMatchTrack = (TH1D*) histo_rebin_structure->Clone();
          TH1D* rebinEmbTrack   = (TH1D*) histo_rebin_structure->Clone();
          HistogramUtilities::rebinHistoToAlternateBinStructure(matchTrackHisto,rebinMatchTrack);
          HistogramUtilities::rebinHistoToAlternateBinStructure(embTrackHisto,rebinEmbTrack);
          delete embTrackHisto;
          delete matchTrackHisto;
          HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
          HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);

          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "Cent " << centIndex << " Match Track Histo  mT-m0 = [" << HistogramUtilities::getLowEdge(rebinMatchTrack) << ", " << HistogramUtilities::getHighEdge(rebinMatchTrack) << "]";
            cout << "         Match: " << HistogramUtilities::sumOfContent(rebinMatchTrack) << "    Embed: " << HistogramUtilities::sumOfContent(rebinEmbTrack) << endl;
          #endif

          if(histo_index == 0){
            graph_eff[centIndex][rapIndex][file_index] = new TGraphAsymmErrors();
            graph_eff[centIndex][rapIndex][file_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
          }else{
            graph_eff_EtaPt[centIndex][rapIndex][file_index] = new TGraphAsymmErrors();
            graph_eff_EtaPt[centIndex][rapIndex][file_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
          }
          HistogramUtilities::ConditionalDelete(rebinMatchTrack);
          HistogramUtilities::ConditionalDelete(rebinEmbTrack);
        }// cent loop




        // CENTRALITY INTEGRATED ------------------------------------
        if(histo_index == 0){
          embTrackHisto   = embTracksAllCent[file_index]->ProjectionY(Form("AllCent_embTrack_%d_%d_%d",histo_index,file_index,rapBin),rapBin,rapBin);
          matchTrackHisto = matchTracksAllCent[file_index]->ProjectionY(Form("AllCent_matchTrack_%d_%d_%d",histo_index,file_index,rapBin),rapBin,rapBin);
        }else{
          embTrackHisto   = embTracksAllCent_EtaPt[file_index]->ProjectionY(Form("AllCent_embTrack_%d_%d_%d",histo_index,file_index,rapBin),rapBin,rapBin);
          matchTrackHisto = matchTracksAllCent_EtaPt[file_index]->ProjectionY(Form("AllCent_matchTrack_%d_%d_%d",histo_index,file_index,rapBin),rapBin,rapBin);
        }

        TH1D* rebinMatchTrack = (TH1D*) histo_rebin_structure->Clone();
        TH1D* rebinEmbTrack   = (TH1D*) histo_rebin_structure->Clone();
        HistogramUtilities::rebinHistoToAlternateBinStructure(matchTrackHisto,rebinMatchTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(embTrackHisto,rebinEmbTrack);
        delete embTrackHisto;
        delete matchTrackHisto;
        HistogramUtilities::makeSqrtNErrors(rebinEmbTrack);
        HistogramUtilities::makeSqrtNErrors(rebinMatchTrack);

        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " All Cent Match Track Histo  mT-m0 = [" << HistogramUtilities::getLowEdge(rebinMatchTrack) << ", " << HistogramUtilities::getHighEdge(rebinMatchTrack) << "]" << endl;
        #endif

        if(histo_index == 0){
          graph_eff_AllCent[rapIndex][file_index] = new TGraphAsymmErrors();
          graph_eff_AllCent[rapIndex][file_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
        }else{
          graph_eff_EtaPt_AllCent[rapIndex][file_index] = new TGraphAsymmErrors();
          graph_eff_EtaPt_AllCent[rapIndex][file_index]->Divide(rebinMatchTrack,rebinEmbTrack,"cl=0.68268949 b(1,1) mode");
        }
        HistogramUtilities::ConditionalDelete(rebinMatchTrack);
        HistogramUtilities::ConditionalDelete(rebinEmbTrack);


      }// file index loop
    }// rapidity loop
  }// histo loop




  // Calcluate Vz Weights For Weighting Algorithm
  bool   isFullVzRange[numFiles];
  double weights[numFiles];
  double total_weights = 0.0;
  double data_sumOfContent = HistogramUtilities::sumOfContent(histo_Vz_data,-1,false);
  for(unsigned int file_index = 0; file_index < a_inFileNames.size(); file_index++){
    if(a_VzRanges[file_index].size() != 2){
      cout << "ERROR: you need to give the Vz reweighting algorithm a vz range with size two!!!!!!!!!" << endl;
      return;
    }
    double lowVzRange = a_VzRanges[file_index][0];
    double highVzRange = a_VzRanges[file_index][1];
    weights[file_index] = HistogramUtilities::sumOfContent(histo_Vz_data,-1,false,lowVzRange,highVzRange);
    total_weights += weights[file_index];
    isFullVzRange[file_index] = (weights[file_index]/data_sumOfContent > 0.96);
  }
  cout << "Vz Weights:" << endl;
  for(unsigned int file_index = 0; file_index < a_inFileNames.size(); file_index++){
    cout << "     File Index: " << file_index << "   weight: " << weights[file_index];
    weights[file_index] = weights[file_index]/total_weights;
    cout << "  fractional weight: " << weights[file_index] << "    isFullRange: " << (isFullVzRange[file_index] ? "TRUE" : "FALSE") << endl;

  }





  //##################   REWEIGHT THE DIFFERENT EFFICIENCIES INTO ONE  ###########################
  for(unsigned int graph_index = 0; graph_index < 2; graph_index++){
    int numBins = 0;
    if(graph_index == 0)  numBins =  matchTracksAllCentAllFiles->GetNbinsX();
    else                  numBins =  matchTracksAllCent_EtaPtAllFiles->GetNbinsX();

    for(int rapBin = 1; rapBin <= numBins; rapBin++){
      int rapIndex = rapBin - 1;
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        if(graph_index == 0)  cout << "##### Rapidity = [" << matchTracksAllCentAllFiles->GetXaxis()->GetBinLowEdge(rapBin) << ", " << matchTracksAllCentAllFiles->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
        else cout << "##### Eta = [" << matchTracksAllCent_EtaPtAllFiles->GetXaxis()->GetBinLowEdge(rapBin) << ", " << matchTracksAllCent_EtaPtAllFiles->GetXaxis()->GetBinLowEdge(rapBin+1) << "]" << endl;
      #endif


      vector< double > mTm0_full_array;
      vector< double > mTm0_ErrLow_full_array;
      vector< double > mTm0_ErrHigh_full_array;
      for(unsigned int file_index = 0; file_index < numFiles; file_index++){
        TGraphAsymmErrors* effGraph_Base = nullptr;
        if(graph_index == 0)  effGraph_Base = graph_eff_AllCent[rapIndex][file_index];
        else                  effGraph_Base = graph_eff_EtaPt_AllCent[rapIndex][file_index];
        if(!effGraph_Base) continue;
        for(int point_index = 0; point_index < effGraph_Base->GetN(); point_index++){
          double mTm0 = effGraph_Base->GetX()[point_index];
          double found = false;
          for(unsigned int index = 0; index < mTm0_full_array.size(); index++){
            if(fabs(mTm0 - mTm0_full_array[index]) < 0.001) found = true;
          }
          if(!found){
            mTm0_full_array.push_back(mTm0);
            mTm0_ErrLow_full_array.push_back(effGraph_Base->GetEXlow()[point_index]);
            mTm0_ErrHigh_full_array.push_back(effGraph_Base->GetEXhigh()[point_index]);
          }
        }
      }

      if(mTm0_full_array.size() == 0){
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "mTm0_full_array is empty: " << mTm0_full_array.size() << endl;
        #endif
        continue;
      }
      for(int centIndex = 0; centIndex < 17; centIndex++){
        TGraphAsymmErrors* graph_eff_combined = nullptr;

        if(centIndex == 16){
          if(graph_index == 0){
            m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
            m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality Integrated y=[%1.2f,%1.2f]; m_{T}-m_{0}; TPC Efficiency",
                       m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),matchTracksAllCentAllFiles->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCentAllFiles->GetXaxis()->GetBinLowEdge(rapBin+1)));
            graph_eff_combined = m_TPC_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex];
          }else{
            m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
            m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality Integrated #eta=[%1.2f,%1.2f]; p_{T}; TPC Efficiency",
                       m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),matchTracksAllCent_EtaPtAllFiles->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCent_EtaPtAllFiles->GetXaxis()->GetBinLowEdge(rapBin+1)));
            graph_eff_combined = m_TPC_Efficiency_AllCent_EtaPt[a_partIndex][rapIndex][pmIndex];
          }

        }else{
          / *if(!embeddedTrackHistogram[centIndex][0]){
            #ifdef _EFFICIENCY_FITTER_DEBUG_
              cout << "embeddedTrackHistogram graph not made, ptr: " << embeddedTrackHistogram[centIndex][0] << endl;
            #endif
            continue;
          }* /
          if(graph_index == 0){
            //cout << "ptr to delete: " << m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] << endl;
            //HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]);
            m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
            m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality %02d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
             matchTracksAllCentAllFiles->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCentAllFiles->GetXaxis()->GetBinLowEdge(rapBin+1)));
            graph_eff_combined = m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex];
          }else{
            //HistogramUtilities::ConditionalDelete(m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]);
            m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
            m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("tpcEfficiencyGraph_%s_Cent%02d_etaIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
            m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("TPC Effiency %s Centrality %d Integrated #eta=[%1.2f,%1.2f]; p_{T}; TPC Efficiency",
                 m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,matchTracksAllCent_EtaPtAllFiles->GetXaxis()->GetBinLowEdge(rapBin),matchTracksAllCent_EtaPtAllFiles->GetXaxis()->GetBinLowEdge(rapBin+1)));
                 graph_eff_combined = m_TPC_Efficiency_ByCent_EtaPt[a_partIndex][centIndex][rapIndex][pmIndex];
          }

          #ifdef _EFFICIENCY_FITTER_DEBUG_
            if(graph_index == 0 && rapBin == 12 && centIndex == 0){
              cout << "###### POINT BY POINT CHECK ----------------------" << endl;
              for(unsigned int file_index = 0; file_index < numFiles; file_index++){
                cout << "EFF GRAPH ---------------------- File: " << file_index  << "   " << (isFullVzRange[file_index] ? "COVERS VZ" : "PARTIAL VZ") << endl;
                HistogramUtilities::printGraphInfo(graph_eff[centIndex][rapIndex][file_index]);
              }
            }
          #endif
        }

        int numPoints = mTm0_full_array.size();
        for(unsigned int point_index = 0; point_index < numPoints; point_index++){ //-------------   point index loop  ------------
          double sumOfWeights = 0.0;
          double efficiency_value      = 0.0;
          double efficiency_error_low  = 0.0;
          double efficiency_error_high = 0.0;
          double mTm0 = mTm0_full_array[point_index];

          if(graph_index == 0 && rapBin == 12 && centIndex == 0) cout << "Constructing Point:  mtm0: " << mTm0 << endl;


          bool atLeastOneFileIsFullCoverage = false;
          for(unsigned int file_index = 0; file_index < numFiles; file_index++){
            TGraphAsymmErrors* effGraph = nullptr;
            if(centIndex == 16){
              if(graph_index == 0) effGraph = graph_eff_AllCent[rapIndex][file_index];
              else                 effGraph = graph_eff_EtaPt_AllCent[rapIndex][file_index];
            }else{
              if(graph_index == 0) effGraph = graph_eff[centIndex][rapIndex][file_index];
              else                 effGraph = graph_eff_EtaPt[centIndex][rapIndex][file_index];
            }

            if(!effGraph || effGraph->GetN() <= 0){
              #ifdef _EFFICIENCY_FITTER_DEBUG_
                if(effGraph) cout << "effGraph graph has no points, NumPoints: " << effGraph->GetN() << endl;
                else cout << "effGraph graph not made, ptr: " << effGraph << endl;
              #endif
              continue;
            }
            int point_searched = HistogramUtilities::getPointClosestToVal(effGraph,mTm0);
            if(fabs(effGraph->GetX()[point_searched] - mTm0) > 0.001){
              cout << "WARNING: Didn't find the right corresponding point in Vz Recombination" << endl;
              continue; // didn't find the right point
            }
            if(isFullVzRange[file_index]) atLeastOneFileIsFullCoverage = true;
            sumOfWeights          += weights[file_index];
            efficiency_value      += weights[file_index]*effGraph->GetY()[point_searched];
            efficiency_error_low  += pow(weights[file_index]*effGraph->GetEYlow()[point_searched],2);
            efficiency_error_high += pow(weights[file_index]*effGraph->GetEYhigh()[point_searched],2);
            if(graph_index == 0 && rapBin == 12 && centIndex == 0) cout << "    Point File " << file_index <<  ":  mtm0: " << effGraph->GetX()[point_searched] << "  eff: " << effGraph->GetY()[point_searched] << endl;

          }// file index loop
          efficiency_value      = efficiency_value/sumOfWeights;
          efficiency_error_low  = sqrt(efficiency_error_low)/sumOfWeights;
          efficiency_error_high = sqrt(efficiency_error_high)/sumOfWeights;

          if(!atLeastOneFileIsFullCoverage || sumOfWeights == 0.0 || !std::isfinite(efficiency_value) || !std::isfinite(efficiency_error_low) || !std::isfinite(efficiency_error_high)){
            cout << "WARNING: Vz Recomb ->  mTm0: " << mTm0 << "  eff:" << efficiency_value << "  el:" << efficiency_error_low << "   eh:" << efficiency_error_high  << endl;
            continue;
          }


          int pointIndexToAdd = graph_eff_combined->GetN();
          graph_eff_combined->SetPoint(pointIndexToAdd,mTm0,efficiency_value);
          graph_eff_combined->SetPointError(pointIndexToAdd,mTm0_ErrLow_full_array[point_index],mTm0_ErrHigh_full_array[point_index],efficiency_error_low,efficiency_error_high);
          if(graph_index == 0 && rapBin == 12 && centIndex == 0) cout << "Adding Point:  mtm0: " << mTm0 << "  eff: " << efficiency_value << endl;
        }// point index loop
        graph_eff_combined->Sort();
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          if(graph_index == 0 && rapBin == 12 && centIndex == 0){
            cout << "RESULTING EFF GRAPH ----------------------" << endl;
            HistogramUtilities::printGraphInfo(m_TPC_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]);
          }
        #endif


      } // cent index loop
    }// rapidity loop
  }// graph loop

  inFile_PicoBinner->Close();

  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "End Call to EfficiencyFitter::loadEmbeddingFiles_ForVzRecombination" << endl;
  #endif
}*/









//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________




void EfficiencyFitter::loadPicoBinnerDataFile( string a_inFileName, int a_partIndex, int a_charge, int a_numPicoBinnerCentBins, int a_mtm0Rebin, bool a_convert16to9CentBins, int a_knockout_dca_rebin, int a_knockout_mTm0_rebin){
  string partNameCharge = m_partInfo->GetParticleName(a_partIndex,a_charge).Data();
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "\n\nStarting EfficiencyFitter::loadPicoBinnerDataFile(" << a_inFileName << "," << a_partIndex << "," << a_charge<< ","
         << a_numPicoBinnerCentBins<< "," << a_mtm0Rebin<< "," << a_convert16to9CentBins <<")" << endl;
    cout << "  Working on " << partNameCharge << endl;
  #endif

  //Get the 2D and 3D Histograms from the file
  gSystem->mkdir(Form("%s/%s",m_imageDir.c_str(),partNameCharge.c_str()));
  //gSystem->mkdir(Form("%s/%s/BTOF_Eff_AllCent",m_imageDir.c_str(),partNameCharge.c_str()));
  //gSystem->mkdir(Form("%s/%s/BTOF_Eff_ByCent",m_imageDir.c_str(),partNameCharge.c_str()));

  TFile* inFile = new TFile(a_inFileName.c_str(),"READ");
  if(!inFile || !inFile->IsOpen()){
    cout << "Warning in  EfficiencyFitter::loadPicoBinnerDataFile.. couldn't find file: " << a_inFileName << endl;
    return;
  }
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Loaded picoBinner file successfully: " << a_inFileName << endl;
  #endif

  string firstDir = "yields";
  string partName = m_partInfo->GetParticleName(a_partIndex).Data();
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;
  inFile->cd();


  //##########    STORE THE RAPIDITY STRUCTURE #####################
  TH3I* tpc_temp = (TH3I*) inFile->Get(Form("yields/%s_space/zTPCPlus_%s_Cent0",m_partInfo->GetParticleName(a_partIndex,0).Data(),m_partInfo->GetParticleName(a_partIndex,0).Data()));
  vector<double> tpc_rap_bin_edges = HistogramUtilities::getBinEdges_X(tpc_temp);
  m_picoBinner_rapidity_structure[a_partIndex][0] = new TH1D(
          Form("picoBinner_%s_tpc_rap_struct",m_partInfo->GetParticleName(a_partIndex,0).Data()),
          Form("Rapidity Structure from PicoBinner for %s; y_{%s}",m_partInfo->GetParticleSymbol(a_partIndex,0).Data(),m_partInfo->GetParticleSymbol(a_partIndex,0).Data()),
            (int)(tpc_rap_bin_edges.size() - 1),&tpc_rap_bin_edges[0]);
  delete tpc_temp;

  TH3I* btof_temp = (TH3I*) inFile->Get(Form("yields/%s_space/zBTOFPlus_%s_Cent0",m_partInfo->GetParticleName(a_partIndex,0).Data(),m_partInfo->GetParticleName(a_partIndex,0).Data()));
  vector<double> btof_rap_bin_edges = HistogramUtilities::getBinEdges_X(btof_temp);
  m_picoBinner_rapidity_structure[a_partIndex][1] = new TH1D(
          Form("picoBinner_%s_btof_rap_struct",m_partInfo->GetParticleName(a_partIndex,0).Data()),
          Form("BTOF Rapidity Structure from PicoBinner for %s; y_{%s}",m_partInfo->GetParticleSymbol(a_partIndex,0).Data(),m_partInfo->GetParticleSymbol(a_partIndex,0).Data()),
            (int)(btof_rap_bin_edges.size() - 1),&btof_rap_bin_edges[0]);
  delete btof_temp;

  TH3I* etof_temp = (TH3I*) inFile->Get(Form("yields/%s_space/zETOFPlus_%s_Cent0",m_partInfo->GetParticleName(a_partIndex,0).Data(),m_partInfo->GetParticleName(a_partIndex,0).Data()));
  if(etof_temp){
    vector<double> etof_rap_bin_edges = HistogramUtilities::getBinEdges_X(etof_temp);
    m_picoBinner_rapidity_structure[a_partIndex][2] = new TH1D(
            Form("picoBinner_%s_etof_rap_struct",m_partInfo->GetParticleName(a_partIndex,0).Data()),
            Form("ETOF Rapidity Structure from PicoBinner for %s; y_{%s}",m_partInfo->GetParticleSymbol(a_partIndex,0).Data(),m_partInfo->GetParticleSymbol(a_partIndex,0).Data()),
              (int)(etof_rap_bin_edges.size() - 1),&etof_rap_bin_edges[0]);
    delete etof_temp;
  }


  TH1D* centEvents = (TH1D*) inFile->Get(Form("centEvents%s", m_partInfo->GetParticleName(a_partIndex,0).Data()));
  if(centEvents){ //if every particle was done at once in a single root file
    m_centEvents[a_partIndex] = centEvents;
  }else{
    cout << "ERROR: Can't Load " << Form("centEvents%s", m_partInfo->GetParticleName(a_partIndex,0).Data()) << endl;
  }


  if(a_convert16to9CentBins){
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << " Converting 16 Centrality bins to 9 Centrality Bins" << endl;
    #endif
    if(m_centEvents[a_partIndex]->GetNbinsX() != 16){
      cout << "ERROR: You cant use a_convert16to9 if you don't have 16 centrality bins... centbins = " << m_centEvents[a_partIndex]->GetNbinsX() << endl;
    }
    TH1D* newCentEvents = new TH1D(Form("centEvents%s_16to9", m_partInfo->GetParticleName(a_partIndex,0).Data()),
     Form("Centrality Events for %s",m_partInfo->GetParticleName(a_partIndex,0).Data()), 9, -0.5, 8.5);

    // 0-5% and 5-10% don't change
    double newBinContent = m_centEvents[a_partIndex]->GetBinContent(1);
    double newBinError   = sqrt(newBinContent);
    newCentEvents->SetBinContent(1,newBinContent);
    newCentEvents->SetBinError(1,newBinError);
    newBinContent = m_centEvents[a_partIndex]->GetBinContent(2);
    newBinError   = sqrt(newBinContent);
    newCentEvents->SetBinContent(2,newBinContent);
    newCentEvents->SetBinError(2,newBinError);

    //make 10-20% = 10-15% + 15-20%
    for(int centIndex = 2; centIndex < 16; centIndex+=2){
      int newCentIndex = 1+centIndex/2;
      double newBinContent = m_centEvents[a_partIndex]->GetBinContent(centIndex+1) + m_centEvents[a_partIndex]->GetBinContent(centIndex+2);
      double newBinError = sqrt(newBinContent);
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "cent: " <<  newCentIndex << " = " << centIndex << " + " << centIndex+1 << "NewBinContent: " << newBinContent  << " +/- " << newBinError << endl;
      #endif
      newCentEvents->SetBinContent(newCentIndex + 1,newBinContent);
      newCentEvents->SetBinError(newCentIndex + 1,newBinError);
    }
    for(int centIndex = 9; centIndex < 16; centIndex++){
      newCentEvents->SetBinContent(centIndex,0.0);
      newCentEvents->SetBinError(centIndex,0.0);
    }
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "CentPtrBefor: " << newCentEvents;
    #endif
    m_centEvents[a_partIndex] = newCentEvents;
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "  after: " << m_centEvents[a_partIndex] << endl;
    #endif
  }
  cout << "Number of Events Per Centrality Bin Histogram Now Loaded with " << m_centEvents[a_partIndex]->GetNbinsX() << " centrality bins" << endl;


  // LOAD PROTON DCA HISTOGRAMS
  if(a_partIndex == 2 && a_charge > 0){
    //TH3I* dca_plus_3D_AllCent = NULL;
    //TH3I* dca_minus_3D_AllCent = NULL;
    for(int centIndex = 0; centIndex < 16; centIndex++){
      TH3I* dca_plus_3D_Int  = (TH3I*) inFile->Get(Form("ProtonDCA/rap_mTm0_dca_ProtonPlus_Cent_%02d",centIndex));
      TH3I* dca_minus_3D_Int = (TH3I*) inFile->Get(Form("ProtonDCA/rap_mTm0_dca_ProtonMinus_Cent_%02d",centIndex));
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Loading Proton DCA for Knockout. Ptrs: " << dca_plus_3D_Int << "  " << dca_minus_3D_Int << endl;
      #endif
      if(!dca_plus_3D_Int || !dca_minus_3D_Int) continue;
      TH3I* dca_plus_3D  = NULL;
      TH3I* dca_minus_3D = NULL;
      if(a_knockout_dca_rebin > 1 || a_knockout_mTm0_rebin > 1){
        dca_plus_3D = HistogramUtilities::reBinVar3D(dca_plus_3D_Int,1,a_knockout_mTm0_rebin,a_knockout_dca_rebin);
        dca_minus_3D = HistogramUtilities::reBinVar3D(dca_minus_3D_Int,1,a_knockout_mTm0_rebin,a_knockout_dca_rebin);
        delete dca_plus_3D_Int;
        delete dca_minus_3D_Int;
      }else{
        dca_plus_3D  = dca_plus_3D_Int;
        dca_minus_3D = dca_minus_3D_Int;
      }

      vector<double> centEdges      = HistogramUtilities::makeEvenBinEdges(17,-0.5,16.5);
      vector<double> rap_bin_edges  = HistogramUtilities::getBinEdges_X(dca_plus_3D);
      vector<double> mTm0_bin_edges = HistogramUtilities::getBinEdges_Y(dca_plus_3D);

      if(centIndex == 0){
        m_knockout_mTm0_structure = new TH1D(
          "knockout_mTm0_struct", "Structure from PicoBinner for Proton Knockout; m_{T}-m_{0} [GeV/c^{2}]",
            (int)(mTm0_bin_edges.size() - 1),&mTm0_bin_edges[0]);
      }

      if(dca_plus_3D && dca_minus_3D){
        if(centIndex == 0){
          m_proton_knockout_ChiSqrd = new TH2D(
            "ProtonKnockout_ChiSqrd","#Chi^{2}/dof for Knockout Fraction Fits for p; y; Centrality",rap_bin_edges.size()-1,&rap_bin_edges[0],
            centEdges.size()-1,&centEdges[0]);
          for(int iii = 0; iii < 16; iii++){
            m_proton_knockout_ChiSqrd->GetYaxis()->SetBinLabel(iii+1, Form("Cent %02d",iii));
          }
          m_proton_knockout_ChiSqrd->GetYaxis()->SetBinLabel(17, "All Cent");
        }

        m_proton_dca_ChiSqrd[centIndex] = new TH2D(
          Form("ProtonDCA_ChiSqrd_Cent%d",centIndex),Form("#Chi^{2}/dof for DCA Fits for p Cent %02d; y; m_{T}-m_{p} [GeV/c^{2}]",centIndex),
          rap_bin_edges.size()-1,&rap_bin_edges[0],mTm0_bin_edges.size()-1,&mTm0_bin_edges[0]);
        m_proton_knockout_histogram_ByCent[centIndex] = new TH2D(
          Form("Proton_Knockout_Histogram_Cent%d",centIndex),Form("#Chi^{2}/dof for DCA Fits for p Cent %02d; y; m_{T}-m_{p} [GeV/c^{2}]",centIndex),
          rap_bin_edges.size()-1,&rap_bin_edges[0],mTm0_bin_edges.size()-1,&mTm0_bin_edges[0]);


        for(unsigned int rapBin = 1; rapBin <= dca_plus_3D->GetNbinsX(); rapBin++){
          int rapIndex = rapBin - 1;
          for(unsigned int mTm0Bin = 1; mTm0Bin <= dca_plus_3D->GetNbinsY(); mTm0Bin++){
            int mTm0Index = mTm0Bin - 1;
            m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][0] = dca_plus_3D->ProjectionZ(Form("ProtonPlus_DCA_Cent%02d_yIndex_%02d_mTm0Index_%02d",centIndex,rapIndex,mTm0Index),rapBin,rapBin,mTm0Bin,mTm0Bin);
            HistogramUtilities::makeSqrtNErrors(m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][0]);
            m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][0]->SetTitle(Form("Proton DCA y=[%1.3f,%1.3f]  m_{T}-m_{p}=[%1.3f,%1.3f] p_{T}=[%1.3f,%1.3f] Cent %d; DCA [cm]; Number of Tracks",
                              dca_plus_3D->GetXaxis()->GetBinLowEdge(rapBin),
                              dca_plus_3D->GetXaxis()->GetBinLowEdge(rapBin+1),
                              dca_plus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin),
                              dca_plus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin+1),
                              PhysMath::pT_from_mTm0(dca_plus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin),m_partInfo->GetParticleMass(2)),
                              PhysMath::pT_from_mTm0(dca_plus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin+1),m_partInfo->GetParticleMass(2)),
                              centIndex
                              ) );

            m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][1] = dca_minus_3D->ProjectionZ(Form("ProtonMinus_DCA_Cent%02d_yIndex_%02d_mTm0Index_%02d",centIndex,rapIndex,mTm0Index),rapBin,rapBin,mTm0Bin,mTm0Bin);
            HistogramUtilities::makeSqrtNErrors(m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][1]);
            m_proton_dca_ByCent[centIndex][rapIndex][mTm0Index][1]->SetTitle(Form("Anti-Proton DCA y=[%1.3f,%1.3f]  m_{T}-m_{p}=[%1.3f,%1.3f]  p_{T}=[%1.3f,%1.3f] Cent %d; DCA [cm]; Number of Tracks",
                              dca_minus_3D->GetXaxis()->GetBinLowEdge(rapBin),
                              dca_minus_3D->GetXaxis()->GetBinLowEdge(rapBin+1),
                              dca_minus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin),
                              dca_minus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin+1),
                              PhysMath::pT_from_mTm0(dca_plus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin),m_partInfo->GetParticleMass(2)),
                              PhysMath::pT_from_mTm0(dca_plus_3D->GetYaxis()->GetBinLowEdge(mTm0Bin+1),m_partInfo->GetParticleMass(2)),
                              centIndex
                              ) );
          }//mTm0Index


        }//yIndex


      }else{
        cout << "WARNING / ERROR: Could not load proton DCA histograms from the picoBinner output files..." << endl;
        continue;
      }
    }//cent loop
  }//partIndex




  TH2D* tpc_forBTOFAllCent = NULL;
  TH2D* btof_forBTOFAllCent = NULL;
  TH2D* tpc_forBTOFHistogramRebin[16];// cent
  TH2D* btof_forBTOFHistogramRebin[16];

  for (int centIndex = 0; centIndex < a_numPicoBinnerCentBins; centIndex++){
    string tpc_track_name = "";
    string btof_track_name = "";
    if(a_charge > 0){
      tpc_track_name = Form("/%s/%s_space/nTPCTracksPlus_ForBTOF_%s_Cent%d",firstDir.c_str(), partName.c_str(),partName.c_str(),centIndex);
      btof_track_name = Form("/%s/%s_space/nBTOFTracksPlus_%s_Cent%d",firstDir.c_str(), partName.c_str(), partName.c_str(),centIndex);
    }else{
      tpc_track_name = Form("/%s/%s_space/nTPCTracksMinus_ForBTOF_%s_Cent%d",firstDir.c_str(), partName.c_str(),partName.c_str(),centIndex);
      btof_track_name = Form("/%s/%s_space/nBTOFTracksMinus_%s_Cent%d",firstDir.c_str(), partName.c_str(), partName.c_str(),centIndex);
    }
    TH2I* tpc_forBTOF_Histogram_Int = (TH2I*) inFile->Get(tpc_track_name.c_str());
    TH2I* btof_forBTOF_Histogram_Int = (TH2I*) inFile->Get(btof_track_name.c_str());

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << " Grabbing: " << endl;
      cout << "     " << tpc_track_name << "   " << tpc_forBTOF_Histogram_Int << endl;
      cout << "     " << btof_track_name << "   " << btof_forBTOF_Histogram_Int << endl;
    #endif

    TH2D* tpc_forBTOF_Histogram  = HistogramUtilities::convertInt_to_DoubleHisto2D(tpc_forBTOF_Histogram_Int);
    TH2D* btof_forBTOF_Histogram = HistogramUtilities::convertInt_to_DoubleHisto2D(btof_forBTOF_Histogram_Int);

    if(!tpc_forBTOF_Histogram){
        if(centIndex == 0){
          cout << " Something is not loading right... tpc_forBTOF_Histogram is NULL for Cent00" << endl;
        }
        continue;
    }
    gSystem->mkdir(Form("%s/%s/BTOF_DataDriven_Efficiency_ByCent/Cent%02d",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex));

    tpc_forBTOFHistogramRebin[centIndex]   = HistogramUtilities::reBinVar2D(tpc_forBTOF_Histogram,1,a_mtm0Rebin);
    btof_forBTOFHistogramRebin[centIndex]  = HistogramUtilities::reBinVar2D(btof_forBTOF_Histogram,1,a_mtm0Rebin);
    HistogramUtilities::makeSqrtNErrors2D(tpc_forBTOFHistogramRebin[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(btof_forBTOFHistogramRebin[centIndex]);

    if(a_mtm0Rebin != 1){
      delete tpc_forBTOF_Histogram;
      delete btof_forBTOF_Histogram;
    }
    delete tpc_forBTOF_Histogram_Int;
    delete btof_forBTOF_Histogram_Int;

    //############      CREATE CENTRALITY INTEGRATED  HISTOS ##################
    if (centIndex == 0){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  Made centrality integrated histograms..." << endl;
      #endif
      tpc_forBTOFAllCent = (TH2D*) tpc_forBTOFHistogramRebin[0]->Clone();
      tpc_forBTOFAllCent->SetName(Form("tpc_forBTOF_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      btof_forBTOFAllCent = (TH2D*) btof_forBTOFHistogramRebin[0]->Clone();
      btof_forBTOFAllCent->SetName(Form("btof_forBTOF_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Name: " << tpc_forBTOFAllCent->GetName()   << " Ptr: "  << tpc_forBTOFAllCent << endl;
        cout << " Name: " << btof_forBTOFAllCent->GetName() << " Ptr: "  << btof_forBTOFAllCent << endl;
      #endif

    }else{
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Added to centrality histogram" << endl;
      #endif
      tpc_forBTOFAllCent->Add(tpc_forBTOFHistogramRebin[centIndex]);
      btof_forBTOFAllCent->Add(btof_forBTOFHistogramRebin[centIndex]);
    }


    //#### MAKE CHI^2 HISTOGRAMS ####
    if(!m_BTOF_fromData_Efficiency_ChiSqrd[a_partIndex][pmIndex]){
      m_BTOF_fromData_Efficiency_ChiSqrd[a_partIndex][pmIndex] = new TH2D(
        Form("BTOFEff_%s_ChiSqrd",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()),
        Form("#Chi^{2}/dof for Nominal BTOF Efficiency Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()),
        btof_forBTOFAllCent->GetNbinsX(),btof_forBTOFAllCent->GetXaxis()->GetBinLowEdge(1),btof_forBTOFAllCent->GetXaxis()->GetBinLowEdge(btof_forBTOFAllCent->GetNbinsX()+1),
        17,-0.5,16.5);
      for(int iii = 0; iii < 16; iii++){
        m_BTOF_fromData_Efficiency_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(iii+1, Form("Cent %02d",iii));
      }
      m_BTOF_fromData_Efficiency_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(17, "All Cent");
      m_BTOF_fromData_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex] = (TH2D*) m_BTOF_fromData_Efficiency_ChiSqrd[a_partIndex][pmIndex]->Clone();
      m_BTOF_fromData_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetName(Form("BTOFEff_%s_ChiSqrd_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
      m_BTOF_fromData_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetTitle(Form("#Chi^{2}/dof for Systematic BTOF Efficiency Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()));
    }


  } //cent Index





  //##########     COMBINE 16 to 9 Centrality Bins  ################
  if(a_convert16to9CentBins){
    cout << "  Converting 16 to 9 Cent Bins..." << endl;
    if(!tpc_forBTOFHistogramRebin[15]){
      cout << "  ERROR: You cant use the a_convert16to9CentBins unless you actually have 16 bins in the feed down file..." << endl;
      return;
    }
    for(int centIndex = 2; centIndex < a_numPicoBinnerCentBins; centIndex+=2){
      int newCentIndex = 1+centIndex/2;
      tpc_forBTOFHistogramRebin[newCentIndex]->Add(tpc_forBTOFHistogramRebin[centIndex],tpc_forBTOFHistogramRebin[centIndex+1]);
      btof_forBTOFHistogramRebin[newCentIndex]->Add(btof_forBTOFHistogramRebin[centIndex],btof_forBTOFHistogramRebin[centIndex+1]);
    }

    for(int centIndex = 9; centIndex < a_numPicoBinnerCentBins; centIndex++){
      delete tpc_forBTOFHistogramRebin[centIndex];
      delete btof_forBTOFHistogramRebin[centIndex];
    }
  }

  //######################      Make The BTOF Projections       ########################
  cout << "  Making Projections for BTOF Efficiency..." << endl;
  for(int rapBin = 1; rapBin <= tpc_forBTOFAllCent->GetNbinsX(); rapBin++){
    int rapIndex = rapBin - 1;
    if(m_testingMidRapidity && rapBin != tpc_forBTOFAllCent->GetXaxis()->FindBin(m_ycm)) continue;
    TH1D* tpcTrackHisto = tpc_forBTOFAllCent->ProjectionY(Form("proj_tpc_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);
    TH1D* btofTrackHsito = btof_forBTOFAllCent->ProjectionY(Form("proj_btof_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);

    if(m_useChebyshevSpacing){
      int firstBin = HistogramUtilities::firstBinWithContent(btofTrackHsito);
      int lastBin  = HistogramUtilities::lastBinWithContent(btofTrackHsito);
      double lowEdgeFirstBin = btofTrackHsito->GetBinLowEdge(firstBin);
      double highEdgeLastBin = btofTrackHsito->GetBinLowEdge(lastBin+1);
      TH1D* rebinBTOFTrack =  HistogramUtilities::rebinHistoToChebyshevSpacing(btofTrackHsito, lowEdgeFirstBin, highEdgeLastBin,3,20);
      TH1D* rebinTPCTrack = (TH1D*) rebinBTOFTrack->Clone();
      HistogramUtilities::zeroBins(rebinTPCTrack);
      HistogramUtilities::rebinHistoToAlternateBinStructure(tpcTrackHisto,rebinTPCTrack);
      delete tpcTrackHisto;
      delete btofTrackHsito;
      tpcTrackHisto = rebinTPCTrack;
      btofTrackHsito = rebinBTOFTrack;
    }else{
      int firstBin = HistogramUtilities::firstBinWithContent(btofTrackHsito);
      int lastBin  = HistogramUtilities::lastBinWithContent(btofTrackHsito);
      double lowEdgeFirstBin = btofTrackHsito->GetBinLowEdge(firstBin);
      double highEdgeLastBin = btofTrackHsito->GetBinLowEdge(lastBin+1);
      TH1D* truncBTOFHisto = HistogramUtilities::truncateHistogram(btofTrackHsito,lowEdgeFirstBin,highEdgeLastBin);
      delete btofTrackHsito;
      btofTrackHsito = truncBTOFHisto;
      TH1D* rebinBTOFTrack = HistogramUtilities::rebinEqualByCDF(btofTrackHsito, m_btof_eff_no_fit_num_bins,false,true); // rebins so more bins where content is changing, less bins where there is less change
      TH1D* rebinTPCTrack = (TH1D*) rebinBTOFTrack->Clone();
      HistogramUtilities::zeroBins(rebinTPCTrack);
      HistogramUtilities::rebinHistoToAlternateBinStructure(tpcTrackHisto,rebinTPCTrack);
      delete tpcTrackHisto;
      delete btofTrackHsito;
      tpcTrackHisto = rebinTPCTrack;
      btofTrackHsito = rebinBTOFTrack;
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "  ## Rapidity Bin = " << rapBin << endl;
      cout << "  TH1D*'s : " << tpcTrackHisto  << " " << btofTrackHsito << endl;
      cout << "  NumBins : " << tpcTrackHisto->GetNbinsX()  << " " << btofTrackHsito->GetNbinsX() << endl;
    #endif

    m_BTOF_fromData_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
    m_BTOF_fromData_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("btofFromDataEfficiencyGraph_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
    m_BTOF_fromData_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("BTOF Effiency %s Centrality Integrated y=[%1.2f,%1.2f]; m_{T}-m_{0}; BTOF Efficiency",
               m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin),tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
    m_BTOF_fromData_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->Divide(btofTrackHsito,tpcTrackHisto,"cl=0.68268949 b(1,1) mode");

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "    BTOF Eff All Cent Ptr: " << m_BTOF_fromData_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex] << "    " <<  m_BTOF_fromData_Efficiency_AllCent[a_partIndex][rapIndex][pmIndex]->GetName() << endl;
    #endif

    delete tpcTrackHisto;
    delete btofTrackHsito;

    for(int centIndex = 0; centIndex < a_numPicoBinnerCentBins; centIndex++){
      if(!btof_forBTOFHistogramRebin[centIndex] || !tpc_forBTOFHistogramRebin[centIndex]){
        cout << "WARNING: read pointers for numer and denom are null: " 
        << btof_forBTOFHistogramRebin[centIndex] << " " << tpc_forBTOFHistogramRebin[centIndex] << endl;
        m_BTOF_fromData_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = NULL;
        continue;
      }
      if(m_testingOneCent && centIndex != 0) continue;

      TH1D* tpcTrackHistoCent = tpc_forBTOFHistogramRebin[centIndex]->ProjectionY(Form("proj_tpc_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);
      TH1D* btofTrackHistoCent = btof_forBTOFHistogramRebin[centIndex]->ProjectionY(Form("proj_btof_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);
      if(m_useChebyshevSpacing){
        int firstBin = HistogramUtilities::firstBinWithContent(btofTrackHistoCent);
        int lastBin  = HistogramUtilities::lastBinWithContent(btofTrackHistoCent);
        double lowEdgeFirstBin = btofTrackHistoCent->GetBinLowEdge(firstBin);
        double highEdgeLastBin = btofTrackHistoCent->GetBinLowEdge(lastBin+1);
        TH1D* rebinBTOFTrack =  HistogramUtilities::rebinHistoToChebyshevSpacing(btofTrackHistoCent, lowEdgeFirstBin, highEdgeLastBin,3,20);
        TH1D* rebinTPCTrack = (TH1D*) rebinBTOFTrack->Clone();
        HistogramUtilities::zeroBins(rebinTPCTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(tpcTrackHistoCent,rebinTPCTrack);
        delete tpcTrackHistoCent;
        delete btofTrackHistoCent;
        tpcTrackHistoCent  = rebinTPCTrack;
        btofTrackHistoCent = rebinBTOFTrack;
      }else{
        int firstBin = HistogramUtilities::firstBinWithContent(btofTrackHistoCent);
        int lastBin  = HistogramUtilities::lastBinWithContent(btofTrackHistoCent);
        double lowEdgeFirstBin = btofTrackHistoCent->GetBinLowEdge(firstBin);
        double highEdgeLastBin = btofTrackHistoCent->GetBinLowEdge(lastBin+1);
        TH1D* truncBTOFHisto = HistogramUtilities::truncateHistogram(btofTrackHistoCent,lowEdgeFirstBin,highEdgeLastBin);
        delete btofTrackHistoCent;
        btofTrackHistoCent = truncBTOFHisto;
        TH1D* rebinBTOFTrack = HistogramUtilities::rebinEqualByCDF(btofTrackHistoCent, 80,false,true); // rebins so more bins where content is changing, less bins where there is less change
        TH1D* rebinTPCTrack = (TH1D*) rebinBTOFTrack->Clone();
        HistogramUtilities::zeroBins(rebinTPCTrack);
        HistogramUtilities::rebinHistoToAlternateBinStructure(tpcTrackHistoCent,rebinTPCTrack);
        delete tpcTrackHistoCent;
        delete btofTrackHistoCent;
        tpcTrackHistoCent = rebinTPCTrack;
        btofTrackHistoCent = rebinBTOFTrack;
      }


      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  TH1D*'s Cent : " << tpcTrackHistoCent  << " " << btofTrackHistoCent <<  endl;
        cout << "     NumBins :   " << tpcTrackHistoCent->GetNbinsX()  << " " << btofTrackHistoCent->GetNbinsX() << endl;
      #endif

      m_BTOF_fromData_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
      m_BTOF_fromData_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("btofFromDataEfficiencyGraph_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
      m_BTOF_fromData_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("BTOF Effiency %s Centrality %02d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
         tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin),tpc_forBTOFAllCent->GetXaxis()->GetBinLowEdge(rapBin+1)));
      m_BTOF_fromData_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->Divide(btofTrackHistoCent,tpcTrackHistoCent,"cl=0.68268949 b(1,1) mode");

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "    BTOF Eff  Cent " << centIndex << " Ptr: " << m_BTOF_fromData_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] << "  " << m_BTOF_fromData_Efficiency_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->GetName() << endl;
      #endif

      delete tpcTrackHistoCent;
      delete btofTrackHistoCent;
      //delete total_histo_cent;

    } // end cent index loop
  }// end rap loop





















  /*
  string firstDir = "yields";
  inFile->cd();
  const char* partName = m_partInfo->GetParticleName(a_partIndex).Data();


  TH2I* TPCTracks_AllCent_Plus = NULL;
  TH2I* TPCTracks_AllCent_Minus = NULL;
  TH2I* TPCTracks_AllCent_Combined = NULL;
  TH2I* BTOFTracks_AllCent_Plus = NULL;
  TH2I* BTOFTracks_AllCent_Minus = NULL;
  TH2I* BTOFTracks_AllCent_Combined = NULL;


  for(int centIndex = 0; centIndex < a_numPicoBinnerCentBins; centIndex++){
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "   Loading picoBinner info for cent: " << centIndex << endl;
    #endif
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

    if(!DTBP_int){
      cout << " Error in EfficiencyFitter::loadPicoBinnerDataFile: Cannot find the TPC Tracks Histogram in the PicoBinner Data Root File" << endl;
    }


    if(centIndex == 0){
      TPCTracks_AllCent_Plus = dynamic_cast<TH2I*>(DTBP_int->Clone());
      TPCTracks_AllCent_Minus = dynamic_cast<TH2I*>(DTBM_int->Clone());
      TPCTracks_AllCent_Combined = dynamic_cast<TH2I*>(DTBP_int->Clone());
      TPCTracks_AllCent_Combined->Add(DTBM_int);
      BTOFTracks_AllCent_Plus = dynamic_cast<TH2I*>(TTP_int->Clone());
      BTOFTracks_AllCent_Minus = dynamic_cast<TH2I*>(TTM_int->Clone());
      BTOFTracks_AllCent_Combined = dynamic_cast<TH2I*>(TTP_int->Clone());
      BTOFTracks_AllCent_Combined->Add(TTM_int);
    }else{
      TPCTracks_AllCent_Plus->Add(DTBP_int);
      TPCTracks_AllCent_Minus->Add(DTBM_int);
      TPCTracks_AllCent_Combined->Add(DTBP_int);
      TPCTracks_AllCent_Combined->Add(DTBM_int);
      BTOFTracks_AllCent_Plus->Add(TTP_int);
      BTOFTracks_AllCent_Minus->Add(TTM_int);
      BTOFTracks_AllCent_Combined->Add(TTP_int);
      BTOFTracks_AllCent_Combined->Add(TTM_int);
    }



    TH2D* DTBP = HistogramUtilities::convertInt_to_DoubleHisto2D(DTBP_int);
    TH2D* DTBM = HistogramUtilities::convertInt_to_DoubleHisto2D(DTBM_int);
    TH2D* TTP = HistogramUtilities::convertInt_to_DoubleHisto2D(TTP_int);
    TH2D* TTM = HistogramUtilities::convertInt_to_DoubleHisto2D(TTM_int);
    delete DTBP_int;
    delete DTBM_int;
    delete TTP_int;
    delete TTM_int;


    //reduce the TPC Tracks to the BTOF Histo's Space (They should already have the same binning, but just in case)
    TH2D* tpcReducedPlus = (TH2D*) TTP->Clone();
    HistogramUtilities::zeroBins2D(tpcReducedPlus);
    for(int binx = 1; binx <= tpcReducedPlus->GetNbinsX(); binx++ ){
      int origBin = DTBP->GetXaxis()->FindBin(tpcReducedPlus->GetXaxis()->GetBinCenter(binx));
      for(int biny = 1; biny <= tpcReducedPlus->GetNbinsY(); biny++ ){
        tpcReducedPlus->SetBinContent(binx, biny, DTBP->GetBinContent(origBin,biny));
      }
    }

    TH2D* tpcReducedMinus = (TH2D*) TTM->Clone();
    HistogramUtilities::zeroBins2D(tpcReducedMinus);
    for(int binx = 1; binx <= tpcReducedMinus->GetNbinsX(); binx++ ){
      int origBin = DTBM->GetXaxis()->FindBin(tpcReducedMinus->GetXaxis()->GetBinCenter(binx));
      for(int biny = 1; biny <= tpcReducedMinus->GetNbinsY(); biny++ ){
        tpcReducedMinus->SetBinContent(binx, biny, DTBM->GetBinContent(origBin,biny));
      }
    }

    TH2D* tpcHisto_Combined = HistogramUtilities::addTH2D(tpcReducedPlus,tpcReducedMinus);
    TH2D* tofHisto_Combined = HistogramUtilities::addTH2D(TTP,TTM);

    HistogramUtilities::makeSqrtNErrors2D(tpcReducedPlus);
    HistogramUtilities::makeSqrtNErrors2D(TTP);
    HistogramUtilities::makeSqrtNErrors2D(tpcReducedMinus);
    HistogramUtilities::makeSqrtNErrors2D(TTM);
    HistogramUtilities::makeSqrtNErrors2D(tpcHisto_Combined);
    HistogramUtilities::makeSqrtNErrors2D(tofHisto_Combined);


    TH2D* ratioHistoPlus = (TH2D*) tpcReducedPlus->Clone();
    ratioHistoPlus->SetTitle(Form("BTOF Efficiency %s Centrality Index = %d", m_partInfo->GetParticleName(a_partIndex,1).Data(), centIndex));
    ratioHistoPlus->SetName(Form("BTofEff_%s_Cent%02d", m_partInfo->GetParticleName(a_partIndex,1).Data(), centIndex));
    ratioHistoPlus->Divide(TTP,tpcReducedPlus,1.0,1.0,"B");//binomial

    TH2D* ratioHistoMinus = (TH2D*) tpcReducedMinus->Clone();
    ratioHistoMinus->SetTitle(Form("BTOF Efficiency %s Centrality Index = %d", m_partInfo->GetParticleName(a_partIndex,-1).Data(), centIndex));
    ratioHistoMinus->SetName(Form("BTofEff_%s_Cent%02d", m_partInfo->GetParticleName(a_partIndex,-1).Data(), centIndex));
    ratioHistoMinus->Divide(TTM,tpcReducedMinus,1.0,1.0,"B");//binomial

    TH2D* ratioHistoCombined = (TH2D*) tpcHisto_Combined->Clone();
    ratioHistoCombined->SetTitle(Form("BTOF Efficiency %s Centrality Index = %d", m_partInfo->GetParticleName(a_partIndex).Data(), centIndex));
    ratioHistoCombined->SetName(Form("BTofEff_%s_Cent%02d", m_partInfo->GetParticleName(a_partIndex).Data(), centIndex));
    ratioHistoCombined->Divide(tofHisto_Combined,tpcHisto_Combined,1.0,1.0,"B");//binomial


    m_BTOFEff_Data_Plus[a_partIndex][centIndex] = ratioHistoPlus;
    m_BTOFEff_Data_Minus[a_partIndex][centIndex] = ratioHistoMinus;
    m_BTOFEff_Data_Combined[a_partIndex][centIndex] = ratioHistoCombined;



    //make the fit parameter histograms
    TH1D* structureHistogram = TTP->ProjectionX("structHisto",2,2);
    HistogramUtilities::zeroBins(structureHistogram);
    for(int fitIndex = 0; fitIndex < 8; fitIndex++){
      for(int paramIndex = 0; paramIndex < 6; paramIndex++){
        m_BTOFEff_Data_FitParams_Plus[a_partIndex][centIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
        m_BTOFEff_Data_FitParams_Plus[a_partIndex][centIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Data_%s_Cent%02d_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,1).Data(),centIndex,fitIndex,paramIndex));
        m_BTOFEff_Data_FitParams_Plus[a_partIndex][centIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Centrality %02d Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,1).Data(),centIndex,fitIndex,paramIndex));

        m_BTOFEff_Data_FitParams_Minus[a_partIndex][centIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
        m_BTOFEff_Data_FitParams_Minus[a_partIndex][centIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Data_%s_Cent%02d_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,-1).Data(),centIndex,fitIndex,paramIndex));
        m_BTOFEff_Data_FitParams_Minus[a_partIndex][centIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Centrality %02d Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,-1).Data(),centIndex,fitIndex,paramIndex));

        m_BTOFEff_Data_FitParams_Combined[a_partIndex][centIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
        m_BTOFEff_Data_FitParams_Combined[a_partIndex][centIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Data_%s_Cent%02d_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex).Data(),centIndex,fitIndex,paramIndex));
        m_BTOFEff_Data_FitParams_Combined[a_partIndex][centIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Centrality %02d Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,1).Data(),centIndex,fitIndex,paramIndex));

        m_BTOFEff_Emb_FitParams_Plus[a_partIndex][centIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
        m_BTOFEff_Emb_FitParams_Plus[a_partIndex][centIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Emb_%s_Cent%02d_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,1).Data(),centIndex,fitIndex,paramIndex));
        m_BTOFEff_Emb_FitParams_Plus[a_partIndex][centIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Centrality %02d Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,1).Data(),centIndex,fitIndex,paramIndex));

        m_BTOFEff_Emb_FitParams_Minus[a_partIndex][centIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
        m_BTOFEff_Emb_FitParams_Minus[a_partIndex][centIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Emb_%s_Cent%02d_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,-1).Data(),centIndex,fitIndex,paramIndex));
        m_BTOFEff_Emb_FitParams_Minus[a_partIndex][centIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Centrality %02d Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,-1).Data(),centIndex,fitIndex,paramIndex));

        m_BTOFEff_Emb_FitParams_Combined[a_partIndex][centIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
        m_BTOFEff_Emb_FitParams_Combined[a_partIndex][centIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Emb_%s_Cent%02d_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex).Data(),centIndex,fitIndex,paramIndex));
        m_BTOFEff_Emb_FitParams_Combined[a_partIndex][centIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Centrality %02d Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,1).Data(),centIndex,fitIndex,paramIndex));
      }
    } // number of fit parameters

  } // end cent loop

  if(!TPCTracks_AllCent_Plus){
    cout << " ERROR: can't find tpc tracks for btof eff histogram in the picoBinner output files" << endl;
  }

  TH2D* TPCTracks_AllCent_Plus_Double      = HistogramUtilities::convertInt_to_DoubleHisto2D(TPCTracks_AllCent_Plus);
  TH2D* TPCTracks_AllCent_Minus_Double     = HistogramUtilities::convertInt_to_DoubleHisto2D(TPCTracks_AllCent_Minus);
  TH2D* TPCTracks_AllCent_Combined_Double  = HistogramUtilities::convertInt_to_DoubleHisto2D(TPCTracks_AllCent_Combined);
  TH2D* BTOFTracks_AllCent_Plus_Double     = HistogramUtilities::convertInt_to_DoubleHisto2D(BTOFTracks_AllCent_Plus);
  TH2D* BTOFTracks_AllCent_Minus_Double    = HistogramUtilities::convertInt_to_DoubleHisto2D(BTOFTracks_AllCent_Minus);
  TH2D* BTOFTracks_AllCent_Combined_Double = HistogramUtilities::convertInt_to_DoubleHisto2D(BTOFTracks_AllCent_Combined);

  HistogramUtilities::makeSqrtNErrors2D(TPCTracks_AllCent_Plus_Double);
  HistogramUtilities::makeSqrtNErrors2D(TPCTracks_AllCent_Minus_Double);
  HistogramUtilities::makeSqrtNErrors2D(TPCTracks_AllCent_Combined_Double);
  HistogramUtilities::makeSqrtNErrors2D(BTOFTracks_AllCent_Plus_Double);
  HistogramUtilities::makeSqrtNErrors2D(BTOFTracks_AllCent_Minus_Double);
  HistogramUtilities::makeSqrtNErrors2D(BTOFTracks_AllCent_Combined_Double);

  TH2D* ratioHistoPlus = (TH2D*) TPCTracks_AllCent_Plus_Double->Clone();
  ratioHistoPlus->SetTitle(Form("BTOF Efficiency %s All Centralities", m_partInfo->GetParticleName(a_partIndex,1).Data()));
  ratioHistoPlus->SetName(Form("BTofEff_%s_AllCent", m_partInfo->GetParticleName(a_partIndex,1).Data()));
  ratioHistoPlus->Divide(BTOFTracks_AllCent_Plus_Double,TPCTracks_AllCent_Plus_Double,1.0,1.0,"B");//bayesian

  TH2D* ratioHistoMinus = (TH2D*) TPCTracks_AllCent_Minus_Double->Clone();
  ratioHistoMinus->SetTitle(Form("BTOF Efficiency %s All Centralities", m_partInfo->GetParticleName(a_partIndex,-1).Data()));
  ratioHistoMinus->SetName(Form("BTofEff_%s_AllCent", m_partInfo->GetParticleName(a_partIndex,-1).Data()));
  ratioHistoMinus->Divide(BTOFTracks_AllCent_Minus_Double,TPCTracks_AllCent_Minus_Double,1.0,1.0,"B");//bayesian

  TH2D* ratioHistoCombined = (TH2D*) TPCTracks_AllCent_Combined_Double->Clone();
  ratioHistoCombined->SetTitle(Form("BTOF Efficiency %s All Centralities", m_partInfo->GetParticleName(a_partIndex).Data()));
  ratioHistoCombined->SetName(Form("BTofEff_%s_AllCent", m_partInfo->GetParticleName(a_partIndex).Data()));
  ratioHistoCombined->Divide(BTOFTracks_AllCent_Combined_Double,TPCTracks_AllCent_Combined_Double,1.0,1.0,"B");//bayesian


  m_BTOFEff_Data_AllCent_Plus[a_partIndex] = ratioHistoPlus;
  m_BTOFEff_Data_AllCent_Minus[a_partIndex] = ratioHistoMinus;
  m_BTOFEff_Data_AllCent_Combined[a_partIndex] = ratioHistoCombined;

  cout << "BTOF Efficieincy Data Loaded and Calculated" << endl;


  //make the fit parameter histograms
  TH1D* structureHistogram = BTOFTracks_AllCent_Plus->ProjectionX("structHisto",2,2);
  HistogramUtilities::zeroBins(structureHistogram);
  for(int fitIndex = 0; fitIndex < 8; fitIndex++){
    for(int paramIndex = 0; paramIndex < 6; paramIndex++){
      m_BTOFEff_Data_FitParams_AllCent_Plus[a_partIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
      m_BTOFEff_Data_FitParams_AllCent_Plus[a_partIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Data_AllCent_%s_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,1).Data(),fitIndex,paramIndex));
      m_BTOFEff_Data_FitParams_AllCent_Plus[a_partIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,1).Data(),fitIndex,paramIndex));

      m_BTOFEff_Data_FitParams_AllCent_Minus[a_partIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
      m_BTOFEff_Data_FitParams_AllCent_Minus[a_partIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Data_AllCent_%s_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,-1).Data(),fitIndex,paramIndex));
      m_BTOFEff_Data_FitParams_AllCent_Minus[a_partIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,-1).Data(),fitIndex,paramIndex));

      m_BTOFEff_Data_FitParams_AllCent_Combined[a_partIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
      m_BTOFEff_Data_FitParams_AllCent_Combined[a_partIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Data_AllCent_%s_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex).Data(),fitIndex,paramIndex));
      m_BTOFEff_Data_FitParams_AllCent_Combined[a_partIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex).Data(),fitIndex,paramIndex));

      m_BTOFEff_Emb_FitParams_AllCent_Plus[a_partIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
      m_BTOFEff_Emb_FitParams_AllCent_Plus[a_partIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Emb_AllCent_%s_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,1).Data(),fitIndex,paramIndex));
      m_BTOFEff_Emb_FitParams_AllCent_Plus[a_partIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,1).Data(),fitIndex,paramIndex));

      m_BTOFEff_Emb_FitParams_AllCent_Minus[a_partIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
      m_BTOFEff_Emb_FitParams_AllCent_Minus[a_partIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Emb_AllCent_%s_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex,-1).Data(),fitIndex,paramIndex));
      m_BTOFEff_Emb_FitParams_AllCent_Minus[a_partIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex,-1).Data(),fitIndex,paramIndex));

      m_BTOFEff_Emb_FitParams_AllCent_Combined[a_partIndex][fitIndex][paramIndex] = dynamic_cast<TH1D*>(structureHistogram->Clone());
      m_BTOFEff_Emb_FitParams_AllCent_Combined[a_partIndex][fitIndex][paramIndex]->SetName(Form("BTOFEff_Emb_AllCent_%s_Fit%d_Param%d",m_partInfo->GetParticleName(a_partIndex).Data(),fitIndex,paramIndex));
      m_BTOFEff_Emb_FitParams_AllCent_Combined[a_partIndex][fitIndex][paramIndex]->SetTitle(Form("BTOF Efficiency %s Fit Funciton %d Parameter %d; Rapidity; Parameter Value",m_partInfo->GetParticleName(a_partIndex).Data(),fitIndex,paramIndex));
    }
  }
  */

}







//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________








void EfficiencyFitter::loadFeedDownFile( string a_inFileName, int a_partIndex, int a_charge, bool a_convert16to9CentBins){

  cout << "\n\nStarting EfficiencyFitter::loadFeedDownFile" << a_inFileName << "," << a_partIndex << "," << a_charge << "," << a_convert16to9CentBins << ")" << endl;

  unsigned int pmIndex = 0;
  if(a_charge == -1) pmIndex = 1;

  cout << "  Loading TFile..." << endl;
  TFile* inFile = new TFile(a_inFileName.c_str(),"READ");
  if(!inFile){
    cout << "Warning in  EfficiencyFitter::loadFeedDownFile.. couldn't find file: " << a_inFileName << endl;
    return;
  }


   // y vs mtm0 for particles
   TH2I* particleYields_primary[16];
   TH2I* particleYields_feedDown[16];
   TH2I* particleYields_total[16];
   for(int iii = 0; iii < 16; iii++){
     particleYields_primary[iii] = NULL;
     particleYields_feedDown[iii] = NULL;
     particleYields_total[iii] = NULL;
   }
   TH2I* particleYields_primary_AllCent = NULL;
   TH2I* particleYields_feedDown_AllCent = NULL;
   TH2I* particleYields_total_AllCent = NULL;

   for(int centIndex = 0; centIndex < 16; centIndex++){
     string primaryName = Form("/%s/Cent%02d/particleYields%s_primary_Cent%d",
                m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), centIndex, m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), centIndex);
     particleYields_primary[centIndex] = (TH2I*) inFile->Get(primaryName.c_str());
     particleYields_feedDown[centIndex] = (TH2I*) inFile->Get(Form("/%s/Cent%02d/particleYields%s_feedDown_Cent%d",
                m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), centIndex, m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), centIndex));
     particleYields_total[centIndex] = (TH2I*) inFile->Get(Form("/%s/Cent%02d/particleYields%s_total_Cent%d",
                m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), centIndex, m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), centIndex));

     #ifdef _EFFICIENCY_FITTER_DEBUG_
       cout << " Trying to grab " << primaryName.c_str() << endl;
       cout << " Read in FeedDown Ptrs: " << particleYields_primary[centIndex] << " " << particleYields_feedDown[centIndex] << " " << particleYields_total[centIndex] << endl;
     #endif

     if(centIndex == 0){
       if(!particleYields_total[0]){
         cout << "ERROR: in loadFeedDownFile. The particleYields---_total_Cent__ can't be accessed!" << endl;
         return;
       }

       particleYields_primary_AllCent = (TH2I*) particleYields_primary[0]->Clone();
       particleYields_primary_AllCent->SetName(Form("particleYields_primary_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
       particleYields_feedDown_AllCent = (TH2I*) particleYields_feedDown[0]->Clone();
       particleYields_feedDown_AllCent->SetName(Form("particleYields_feedDown_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
       particleYields_total_AllCent = (TH2I*) particleYields_total[0]->Clone();
       particleYields_total_AllCent->SetName(Form("particleYields_total_AllCent_%s",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
     }else{
      if(!particleYields_total[centIndex]) continue;
      particleYields_primary_AllCent->Add(particleYields_primary[centIndex]);
      particleYields_feedDown_AllCent->Add(particleYields_feedDown[centIndex]);
      particleYields_total_AllCent->Add(particleYields_total[centIndex]);
     }
   } // end cent index loop

  int maxCentBins = 16;
  if(a_convert16to9CentBins){
    cout << "  Converting 16 to 9 Cent Bins..." << endl;
    if(!particleYields_total[15]){
      cout << "  ERROR: You cant use the a_convert16to9CentBins unless you actually have 16 bins in the feed down file..." << endl;
      return;
    }
    maxCentBins = 9;
    for(int centIndex = 2; centIndex < 16; centIndex+=2){
      int newCentIndex = 1+centIndex/2;
      particleYields_primary[newCentIndex]->Add(particleYields_primary[centIndex],particleYields_primary[centIndex+1]);
      particleYields_feedDown[newCentIndex]->Add(particleYields_feedDown[centIndex],particleYields_feedDown[centIndex+1]);
      particleYields_total[newCentIndex]->Add(particleYields_total[centIndex],particleYields_total[centIndex+1]);
    }

    for(int centIndex = 9; centIndex < 16; centIndex++){
      delete particleYields_primary[centIndex];
      delete particleYields_feedDown[centIndex];
      delete particleYields_total[centIndex];
      //particleYields_primary[centIndex] = NULL;
      //particleYields_feedDown[centIndex] = NULL;
      //particleYields_total[centIndex] = NULL;
    }

  }


  if(!m_FeedDown_ChiSqrd[a_partIndex][pmIndex]){
    m_FeedDown_ChiSqrd[a_partIndex][pmIndex] = new TH2D(
      Form("feeddown_%s_ChiSqrd",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()),
      Form("#Chi^{2}/dof for Nominal Feed-Down Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()),
      particleYields_primary_AllCent->GetNbinsX(),particleYields_primary_AllCent->GetXaxis()->GetBinLowEdge(1),particleYields_primary_AllCent->GetXaxis()->GetBinLowEdge(particleYields_primary_AllCent->GetNbinsX()+1),
      17,-0.5,16.5);
    for(int iii = 0; iii < 16; iii++){
      m_FeedDown_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(iii+1, Form("Cent %02d",iii));
    }
    m_FeedDown_ChiSqrd[a_partIndex][pmIndex]->GetYaxis()->SetBinLabel(17, "All Cent");
    m_FeedDown_ChiSqrd_Sys[a_partIndex][pmIndex] = (TH2D*) m_FeedDown_ChiSqrd[a_partIndex][pmIndex]->Clone();
    m_FeedDown_ChiSqrd_Sys[a_partIndex][pmIndex]->SetName(Form("feeddown_%s_ChiSqrd_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    m_FeedDown_ChiSqrd_Sys[a_partIndex][pmIndex]->SetTitle(Form("#Chi^{2}/dof for Systematic Feed-Down Fits for %s; y; Centrality",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data()));
  }





  //Convert TH2I to TH2D and make proper square root errors
  cout << "  Converting TH2I to TH2D..." << endl;
  TH2D* particleYields_primary_double[16];
  TH2D* particleYields_feedDown_double[16];
  TH2D* particleYields_total_double[16];
   for(int iii = 0; iii < 16; iii++){
     particleYields_primary_double[iii] = NULL;
     particleYields_feedDown_double[iii] = NULL;
     particleYields_total_double[iii] = NULL;
   }

  TH2D* particleYields_primary_AllCent_double= HistogramUtilities::convertInt_to_DoubleHisto2D(particleYields_primary_AllCent);
  TH2D* particleYields_feedDown_AllCent_double= HistogramUtilities::convertInt_to_DoubleHisto2D(particleYields_feedDown_AllCent);
  TH2D* particleYields_total_AllCent_double= HistogramUtilities::convertInt_to_DoubleHisto2D(particleYields_total_AllCent);
  HistogramUtilities::makeSqrtNErrors2D(particleYields_primary_AllCent_double);
  HistogramUtilities::makeSqrtNErrors2D(particleYields_feedDown_AllCent_double);
  HistogramUtilities::makeSqrtNErrors2D(particleYields_total_AllCent_double);
  for(int centIndex = 0; centIndex < maxCentBins; centIndex++){
    if(!particleYields_primary[centIndex]){
      continue;
    }
    particleYields_primary_double[centIndex]  = HistogramUtilities::convertInt_to_DoubleHisto2D(particleYields_primary[centIndex]);
    particleYields_feedDown_double[centIndex] = HistogramUtilities::convertInt_to_DoubleHisto2D(particleYields_feedDown[centIndex]);
    particleYields_total_double[centIndex]    = HistogramUtilities::convertInt_to_DoubleHisto2D(particleYields_total[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(particleYields_primary_double[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(particleYields_feedDown_double[centIndex]);
    HistogramUtilities::makeSqrtNErrors2D(particleYields_total_double[centIndex]);
  } // end cent index loop





  // #####  MAKE TGraphAsymmErrors for the Primary/Total and FeedDown/Total Feed Down Corrections #######
  cout << "  Creating Primary/Total and FeedDown/Total TGraphAsymmErrors..." << endl;
  for(int rapBin = 1; rapBin <= particleYields_primary_AllCent_double->GetNbinsX(); rapBin++){
    int rapIndex = rapBin - 1;

    TH1D* primary_histo = particleYields_primary_AllCent_double->ProjectionY(Form("proj_primary_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);
    TH1D* feeddown_histo = particleYields_feedDown_AllCent_double->ProjectionY(Form("proj_feeddown_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);
    TH1D* total_histo = particleYields_total_AllCent_double->ProjectionY(Form("proj_total_rap%02d_%02d_%d",rapIndex,a_partIndex,pmIndex),rapBin,rapBin);

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "  ## Rapidity Bin = " << rapBin << endl;
      cout << "  TH1D*'s : " << primary_histo  << " " << feeddown_histo << " " <<  total_histo << endl;
      cout << "  NumBins : " << primary_histo->GetNbinsX()  << " " << feeddown_histo->GetNbinsX() << " " <<  total_histo->GetNbinsX() << endl;

    #endif

    m_FeedDown_primaryOverTotal_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
    m_FeedDown_primaryOverTotal_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("feedDown_primaryOverTotal_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
    m_FeedDown_primaryOverTotal_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("Primary / Total %s Centrality Integrated y=[%1.2f,%1.2f]",
               m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin),particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin+1)));
    m_FeedDown_primaryOverTotal_AllCent[a_partIndex][rapIndex][pmIndex]->Divide(primary_histo,total_histo,"cl=0.68268949 b(1,1) mode");

    m_FeedDown_feeddownOverTotal_AllCent[a_partIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
    m_FeedDown_feeddownOverTotal_AllCent[a_partIndex][rapIndex][pmIndex]->SetName(Form("feedDown_feedDownOverTotal_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),rapIndex));
    m_FeedDown_feeddownOverTotal_AllCent[a_partIndex][rapIndex][pmIndex]->SetTitle(Form("Feed-Down / Total %s Centrality Integrated y=[%1.2f,%1.2f]",
               m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin),particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin+1)));
    m_FeedDown_feeddownOverTotal_AllCent[a_partIndex][rapIndex][pmIndex]->Divide(feeddown_histo,total_histo,"cl=0.68268949 b(1,1) mode");

    delete primary_histo;
    delete feeddown_histo;
    delete total_histo;

    for(int centIndex = 0; centIndex < maxCentBins; centIndex++){
      if(!particleYields_primary[centIndex]){
        //m_FeedDown_primaryOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = NULL;
        //m_FeedDown_feeddownOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = NULL;
        continue;
      }

      TH1D* primary_histo_cent = particleYields_primary_double[centIndex]->ProjectionY(Form("proj_primary_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);
      TH1D* feeddown_histo_cent = particleYields_feedDown_double[centIndex]->ProjectionY(Form("proj_feeddown_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);
      TH1D* total_histo_cent = particleYields_total_double[centIndex]->ProjectionY(Form("proj_total_rap%02d_cent%02d_%02d_%d",rapIndex,centIndex,a_partIndex,pmIndex),rapBin,rapBin);

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "  TH1D*'s Cent : " << primary_histo_cent  << " " << feeddown_histo_cent << " " <<  total_histo_cent << endl;
        cout << "     NumBins :   " << primary_histo_cent->GetNbinsX()  << " " << feeddown_histo_cent->GetNbinsX() << " " <<  total_histo_cent->GetNbinsX() << endl;
      #endif

      m_FeedDown_primaryOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
      m_FeedDown_primaryOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("feedDown_primaryOverTotal_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
      m_FeedDown_primaryOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("Primary / Total %s Centrality %02d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
         particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin),particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin+1)));
      m_FeedDown_primaryOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->Divide(primary_histo_cent,total_histo_cent,"cl=0.68268949 b(1,1) mode");


      m_FeedDown_feeddownOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex] = new TGraphAsymmErrors();
      m_FeedDown_feeddownOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetName(Form("feedDown_feedDownOverTotal_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,rapIndex));
      m_FeedDown_feeddownOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->SetTitle(Form("Feed-Down / Total %s Centrality %02d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
         particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin),particleYields_primary_AllCent_double->GetXaxis()->GetBinLowEdge(rapBin+1)));
      m_FeedDown_feeddownOverTotal_ByCent[a_partIndex][centIndex][rapIndex][pmIndex]->Divide(feeddown_histo_cent,total_histo_cent,"cl=0.68268949 b(1,1) mode");


      delete primary_histo_cent;
      delete feeddown_histo_cent;
      delete total_histo_cent;

    } // end cent index loop



  }



  cout << "  Deleting Unnecessary Histograms..." << endl;



  // delete all of the TH2I & TH2D  histograms
  if(particleYields_primary_AllCent_double) delete particleYields_primary_AllCent_double;
  if(particleYields_feedDown_AllCent_double) delete particleYields_feedDown_AllCent_double;
  if(particleYields_total_AllCent_double) delete particleYields_total_AllCent_double;
  if(particleYields_primary_AllCent) delete particleYields_primary_AllCent;
  if(particleYields_feedDown_AllCent) delete particleYields_feedDown_AllCent;
  if(particleYields_total_AllCent) delete particleYields_total_AllCent;
  for(int centIndex = 0; centIndex < 16; centIndex++){
    if(!particleYields_primary_double[centIndex]) continue;
    if(particleYields_primary_double[centIndex])  delete particleYields_primary_double[centIndex];
    if(particleYields_feedDown_double[centIndex]) delete particleYields_feedDown_double[centIndex];
    if(particleYields_total_double[centIndex])    delete particleYields_total_double[centIndex];
    if(particleYields_primary[centIndex])         delete particleYields_primary[centIndex];
    if(particleYields_feedDown[centIndex])        delete particleYields_feedDown[centIndex];
    if(particleYields_total[centIndex])           delete particleYields_total[centIndex];
  }

  m_particleLoaded_feeddown[a_partIndex][pmIndex] = true;

}



//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void EfficiencyFitter::loadHybridFeedDown(string a_filename, int a_partIndex, int a_charge){

  cout << "EfficiencyFitter::loadHybridFeedDown( " << a_filename << " , " << a_partIndex << " , " << a_charge << " )" << endl;
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile){
    cout << "ERROR: in EfficiencyFitter::loadHybridFeedDown couldn't find file: " << a_filename << endl;
    return;
  }

  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;

  //###### CENTRALITY INDEX LOOP ####################
  for(int centIndex = 0; centIndex < 16; centIndex++){
    TH2D* feedDownSpectra = (TH2D*) inFile->Get(Form("totalFeedDownSpectra_%s_Cent%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex));
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "FeedDown Spectra Ptrs Cent " << centIndex << " : " << feedDownSpectra <<  endl;
    #endif
    if(!feedDownSpectra){
      cout << "WARNING: Did not load Feed-Down Spectra for " << m_partInfo->GetParticleName(a_partIndex,a_charge) << " Cent " << centIndex << endl;
      continue;
    }

    vector< TH1D* > projList           = HistogramUtilities::VerticalProjectionHistos(feedDownSpectra, "y", "y");
    vector< TGraphErrors* > projGraphs = HistogramUtilities::VerticalProjectionGraphs(feedDownSpectra, "y", "y");
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "     Size of projections " << projList.size() <<  endl;
    #endif
    //## USE THE TPC BINNING TO DETERMINE THE CORRECT RAPIDITY BINS TO USE
    //TAxis* rapidity_bin_axis = m_BTOF_fromData_Efficiency_ChiSqrd[a_partIndex][pmIndex]->GetXaxis();
    //int numOfficialYBins = rapidity_bin_axis->GetNbins(); //->FindBin(m_ycm)
    for(int yIndex = 0; yIndex < feedDownSpectra->GetNbinsX(); yIndex++){
      //double binCenter = rapidity_bin_axis->GetBinCenter(yIndex+1);
      //int yIndexHybrid = feedDownSpectra->GetXaxis()->FindBin(binCenter) - 1;
      m_FeedDown_Hybrid_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = (TH1D*) projList[yIndex]->Clone();
      m_FeedDown_Hybrid_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("hybridFeedDown_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex,yIndex));
      m_FeedDown_Hybrid_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("hybridFeedDown Spectra for %s Cent%d y=[%1.2f,%1.2f];m_{T}-m_{0} [GeV/c^{2}];#frac{1}{N_{evt}} #frac{1}{2#pi m_{T}} #frac{d^{2}N}{dm_{T}dy}",
           m_partInfo->GetParticleSymbol(a_partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex,feedDownSpectra->GetXaxis()->GetBinLowEdge(yIndex+1),feedDownSpectra->GetXaxis()->GetBinLowEdge(yIndex+2)));

      m_FeedDown_Hybrid_ByCent_graphs[a_partIndex][centIndex][yIndex][pmIndex] = (TGraphErrors*) projGraphs[yIndex]->Clone();
      m_FeedDown_Hybrid_ByCent_graphs[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("hybridFeedDown_%s_Cent%02d_yIndex%02d_graph",m_partInfo->GetParticleName(a_partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex,yIndex));
      m_FeedDown_Hybrid_ByCent_graphs[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("hybridFeedDown Spectra for %s Cent%d y=[%1.2f,%1.2f];m_{T}-m_{0} [GeV/c^{2}];#frac{1}{N_{evt}} #frac{1}{2#pi m_{T}} #frac{d^{2}N}{dm_{T}dy}",
           m_partInfo->GetParticleSymbol(a_partIndex,pmIndex == 0 ? 1 : -1).Data(),centIndex,feedDownSpectra->GetXaxis()->GetBinLowEdge(yIndex+1),feedDownSpectra->GetXaxis()->GetBinLowEdge(yIndex+2)));


      //#ifdef _EFFICIENCY_FITTER_DEBUG_
      //  cout << "Finding coresponding bin for y=" << binCenter << " found to be yIndex_feeddown=" << yIndexHybrid << " histoPtr: " << m_FeedDown_Hybrid_ByCent[a_partIndex][centIndex][yIndex][pmIndex] <<  endl;
      //#endif

    }//yIndex


    //##########    STORE THE RAPIDITY STRUCTURE #####################
    if(centIndex == 0 && !m_hybridFeedDown_rapidity_structure[a_partIndex]){
      vector<double> rap_bin_edges = HistogramUtilities::getBinEdges_X(feedDownSpectra);
      m_hybridFeedDown_rapidity_structure[a_partIndex] = new TH1D(Form("hybridFeedDown_%s_rap_struct",m_partInfo->GetParticleName(a_partIndex,0).Data()),
        Form("Efficiency Rapidity Structure from Embedding for %s; y_{%s}",
              m_partInfo->GetParticleSymbol(a_partIndex,0).Data(),m_partInfo->GetParticleSymbol(a_partIndex,0).Data()),
              (int)(rap_bin_edges.size() - 1),&rap_bin_edges[0]);
    }


    delete feedDownSpectra;
  }//cent

  //inFile->Close();


}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void EfficiencyFitter::copyAllCentTPCFitsToSpecCent(int a_partIndex, int a_charge){
  int pmIndex = 0;
  if(a_charge < 1) pmIndex = 1;
  TAxis* rapidity_bin_axis = m_embedding_rapidity_structure[a_partIndex]->GetXaxis();
  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int yIndex = 0; yIndex < g_EffFit_MaxNumRapBins; yIndex++){
      TF1* allCentFit                      = m_TPC_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex];
      TMatrixDSym* allCentFitCovMat        = m_TPC_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex];
      TF1* allCentFit_highPt               = m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex];
      TMatrixDSym* allCentFitCovMat_highPt = m_TPC_Eff_CovMat_HighPt_AllCent[a_partIndex][yIndex][pmIndex];
      if(allCentFit){
        m_TPC_Eff_Fits_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = (TF1*) allCentFit->Clone();
        m_TPC_Eff_Fits_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,yIndex));
        m_TPC_Eff_Fits_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent %d y=[%1.2f,%1.2f] (copied from Centrality Integrated)",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      }
      if(allCentFitCovMat){
        m_TPC_Eff_CovMat_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = allCentFitCovMat;
      }
      if(allCentFit_highPt){
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = allCentFit_highPt;
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_HighPt_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,yIndex));
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s High p_{T} Cent %d y=[%1.2f,%1.2f] (copied from Centrality Integrated)",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      }
      if(allCentFitCovMat_highPt){
        m_TPC_Eff_CovMat_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = allCentFitCovMat_highPt;
      }
    }
  }

  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int yIndex = 0; yIndex < g_EffFit_MaxNumRapBins; yIndex++){
      TF1* allCentFit = m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex];
      TMatrixDSym* allCentFitCovMat = m_TPC_Eff_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex];
      if(allCentFit){
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex] = (TF1*) allCentFit->Clone();
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,yIndex));
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent %d y=[%1.2f,%1.2f] (copied from Centrality Integrated)",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      }
      if(allCentFitCovMat){
        m_TPC_Eff_CovMat_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex] = allCentFitCovMat;
      }
    }
  }
}

//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void EfficiencyFitter::copyAllCentBTOFFitsToSpecCent(int a_partIndex, int a_charge){
  int pmIndex = 0;
  if(a_charge < 1) pmIndex = 1;
  TAxis* rapidity_bin_axis = m_picoBinner_rapidity_structure[a_partIndex][1]->GetXaxis();

  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int yIndex = 0; yIndex < g_EffFit_MaxNumRapBins; yIndex++){
      TF1* allCentFit = m_BTOF_fromData_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex];
      TMatrixDSym* allCentFitCovMat = m_BTOF_fromData_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex];
      TF1* allCentFit_highPt = m_BTOF_fromData_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex];
      TMatrixDSym* allCentFitCovMat_highPt = m_BTOF_fromData_Eff_CovMat_HighPt_AllCent[a_partIndex][yIndex][pmIndex];
      if(allCentFit){
        m_BTOF_fromData_Eff_Fits_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = (TF1*) allCentFit->Clone();
          m_BTOF_fromData_Eff_Fits_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,yIndex));
          m_BTOF_fromData_Eff_Fits_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent %d y=[%1.2f,%1.2f] (copied from centrality integrated)",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      }
      if(allCentFitCovMat){
        m_BTOF_fromData_Eff_CovMat_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = allCentFitCovMat;
      }
      if(allCentFit_highPt){
        m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = allCentFit_highPt;
            m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_HighPt_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,yIndex));
            m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s High p_{T}Cent %d y=[%1.2f,%1.2f] (copied from centrality integrated)",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      }
      if(allCentFitCovMat_highPt){
        m_BTOF_fromData_Eff_CovMat_HighPt_ByCent[a_partIndex][centIndex][yIndex][pmIndex] = allCentFitCovMat_highPt;
      }
    }
  }

  for(int centIndex = 0; centIndex < 16; centIndex++){
    for(int yIndex = 0; yIndex < g_EffFit_MaxNumRapBins; yIndex++){
      TF1* allCentFit = m_BTOF_fromData_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex];
      TMatrixDSym* allCentFitCovMat = m_BTOF_fromData_Eff_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex];
      if(allCentFit){
        m_BTOF_fromData_Eff_Fits_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex] = (TF1*) allCentFit->Clone();
          m_BTOF_fromData_Eff_Fits_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,yIndex));
          m_BTOF_fromData_Eff_Fits_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent %d y=[%1.2f,%1.2f] (copied from centrality integrated)",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      }
      if(allCentFitCovMat){
        m_BTOF_fromData_Eff_CovMat_ByCent_Sys[a_partIndex][centIndex][yIndex][pmIndex] = allCentFitCovMat;
      }
    }
  }

}



//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void EfficiencyFitter::write(string a_outFileName){

  //Create the Output Correction File
  cout << "Writing to output" << endl;
  TFile *outFile = new TFile(a_outFileName.c_str(),"RECREATE");
  if(!outFile->IsOpen()) cout << "Correction File is not open!!" << endl;
  for(int pmIndex = 0; pmIndex < 2; pmIndex++){
    int charge = 1;
    if(pmIndex == 1) charge = -1;
    for(int partIndex = 0; partIndex < 9; partIndex++){
      if(!m_particleLoaded_embedding[partIndex]) continue;
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Writing Part Index: " << partIndex << " charge: " << charge << endl;
      #endif
      outFile->mkdir(Form("%s",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s",m_partInfo->GetParticleName(partIndex,charge).Data()));

      HistogramUtilities::ConditionalWrite(m_picoBinner_rapidity_structure[partIndex][0]);
      HistogramUtilities::ConditionalWrite(m_picoBinner_rapidity_structure[partIndex][1]);
      HistogramUtilities::ConditionalWrite(m_picoBinner_rapidity_structure[partIndex][2]);
      HistogramUtilities::ConditionalWrite(m_embedding_rapidity_structure[partIndex]);
      HistogramUtilities::ConditionalWrite(m_embedding_pseudorapidity_structure[partIndex]);
      HistogramUtilities::ConditionalWrite(m_hybridFeedDown_rapidity_structure[partIndex]);
      HistogramUtilities::ConditionalWrite(m_knockout_mTm0_structure);

      HistogramUtilities::ConditionalWrite(m_TPC_Efficiency_ChiSqrd[partIndex][pmIndex]);
      HistogramUtilities::ConditionalWrite(m_TPC_Efficiency_ChiSqrd_Sys[partIndex][pmIndex]);
      HistogramUtilities::ConditionalWrite(m_EnergyLoss_ChiSqrd[partIndex][pmIndex]);
      HistogramUtilities::ConditionalWrite(m_EnergyLoss_ChiSqrd_Sys[partIndex][pmIndex]);
      HistogramUtilities::ConditionalWrite(m_FeedDown_ChiSqrd[partIndex][pmIndex]);
      HistogramUtilities::ConditionalWrite(m_FeedDown_ChiSqrd_Sys[partIndex][pmIndex]);

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Writing Efficiency Information" << endl;
      #endif
      outFile->mkdir(Form("%s/EfficiencyGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EfficiencyGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        //cout << "Trying to write : " << m_TPC_Efficiency_AllCent[partIndex][rapIndex][pmIndex] << "  part: " << partIndex << "  rap: " << rapIndex << "  pmIndex: " << pmIndex << endl;
        HistogramUtilities::ConditionalWrite(m_TPC_Efficiency_AllCent[partIndex][rapIndex][pmIndex]);
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_TPC_Efficiency_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
        }
      }
      outFile->mkdir(Form("%s/EfficiencyFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EfficiencyFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_TPC_Eff_Fits_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_TPC_Eff_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_TPC_Eff_CovMat_AllCent[partIndex][rapIndex][pmIndex],Form("tpcEfficiencyCovMat_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
        HistogramUtilities::ConditionalWrite(m_TPC_Eff_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex],Form("tpcEfficiencyCovMat_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
        HistogramUtilities::ConditionalWrite(m_TPC_Eff_Fits_HighPt_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_TPC_Eff_CovMat_HighPt_AllCent[partIndex][rapIndex][pmIndex],Form("tpcEfficiencyCovMat_HighPt_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_TPC_Eff_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_TPC_Eff_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          //HistogramUtilities::ConditionalWrite(m_TPC_Eff_FitErrGraphs_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          //HistogramUtilities::ConditionalWrite(m_TPC_Eff_FitErrGraphs_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          //cout << "PTR: " << m_TPC_Eff_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex] << endl;
          HistogramUtilities::ConditionalWrite(m_TPC_Eff_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex],Form("tpcEfficiencyCovMat_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_TPC_Eff_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex],Form("tpcEfficiencyCovMat_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_TPC_Eff_Fits_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          cout << "PTR: " << m_TPC_Eff_CovMat_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex] <<  " cent" << centIndex << " rapIndex:" << rapIndex << endl;
          HistogramUtilities::ConditionalWrite(m_TPC_Eff_CovMat_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex],Form("tpcEfficiencyCovMat_HighPt_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
        }
      }

      outFile->mkdir(Form("%s/EfficiencyGraphs_EtaPt",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EfficiencyGraphs_EtaPt",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int etaIndex = 0; etaIndex < g_EffFit_MaxNumRapBins; etaIndex++){
        HistogramUtilities::ConditionalWrite(m_TPC_Efficiency_AllCent_EtaPt[partIndex][etaIndex][pmIndex]);
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int etaIndex = 0; etaIndex < g_EffFit_MaxNumRapBins; etaIndex++){
          HistogramUtilities::ConditionalWrite(m_TPC_Efficiency_ByCent_EtaPt[partIndex][centIndex][etaIndex][pmIndex]);
        }
      }

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Writing BTOF Efficiency Information" << endl;
      #endif
      outFile->mkdir(Form("%s/BTOFEfficiencyGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/BTOFEfficiencyGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_BTOF_Efficiency_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Efficiency_AllCent[partIndex][rapIndex][pmIndex]);
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_BTOF_Efficiency_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Efficiency_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
        }
      }
      outFile->mkdir(Form("%s/BTOFEfficiencyFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/BTOFEfficiencyFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_BTOF_Eff_Fits_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_BTOF_Eff_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_BTOF_Eff_CovMat_AllCent[partIndex][rapIndex][pmIndex],Form("btofEfficiencyCovMat_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
        HistogramUtilities::ConditionalWrite(m_BTOF_Eff_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex],Form("btofEfficiencyCovMat_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));

        HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_Fits_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_CovMat_AllCent[partIndex][rapIndex][pmIndex],    Form("btofEfficiencyCovMat_DataDriven_%s_AllCent_yIndex%d",    m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
        HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex],Form("btofEfficiencyCovMat_DataDriven_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
        HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_Fits_HighPt_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_CovMat_HighPt_AllCent[partIndex][rapIndex][pmIndex],Form("btofEfficiencyCovMat_DataDriven_HighPt_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_BTOF_Eff_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_BTOF_Eff_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_BTOF_Eff_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex]    ,Form("btofEfficiencyCovMat_%s_Cent%d_yIndex%d",    m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_BTOF_Eff_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex],Form("btofEfficiencyCovMat_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));

          HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex],    Form("btofEfficiencyCovMat_DataDriven_%s_Cent%d_yIndex%d",    m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex],Form("btofEfficiencyCovMat_DataDriven_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_Fits_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_BTOF_fromData_Eff_CovMat_HighPt_ByCent[partIndex][centIndex][rapIndex][pmIndex],Form("btofEfficiencyCovMat_DataDriven_HighPt_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
        }
      }

      outFile->mkdir(Form("%s/EfficiencyStudy",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EfficiencyStudy",m_partInfo->GetParticleName(partIndex,charge).Data()));
	    HistogramUtilities::ConditionalWrite( m_TPC_Efficiency_ChiSqrd_Compare[partIndex][pmIndex] );

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Writing Energy Loss Information" << endl;
      #endif
      outFile->mkdir(Form("%s/EnergyLossGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EnergyLossGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_EnergyLoss_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_EnergyLossProfile_AllCent[partIndex][rapIndex][pmIndex]);
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_EnergyLoss_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_EnergyLossProfile_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
        }
      }

      outFile->mkdir(Form("%s/EnergyLossFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EnergyLossFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_EnergyLoss_Fits_AllCent[partIndex][rapIndex][pmIndex]);
        //HistogramUtilities::ConditionalWrite(m_EnergyLoss_FitErrGraphs_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_EnergyLoss_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex]);
        //HistogramUtilities::ConditionalWrite(m_EnergyLoss_FitErrGraphs_AllCent_Sys[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_EnergyLoss_CovMat_AllCent[partIndex][rapIndex][pmIndex],Form("energyLossCovMat_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
        HistogramUtilities::ConditionalWrite(m_EnergyLoss_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex],Form("energyLossCovMat_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_EnergyLoss_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          //HistogramUtilities::ConditionalWrite(m_EnergyLoss_FitErrGraphs_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_EnergyLoss_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          //HistogramUtilities::ConditionalWrite(m_EnergyLoss_FitErrGraphs_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_EnergyLoss_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex],Form("energyLossCovMat_%s_Cent%02d_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_EnergyLoss_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex],Form("energyLossCovMat_%s_Cent%02d_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
        }
      }

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Writing Feed-Down Information" << endl;
      #endif
      outFile->mkdir(Form("%s/FeedDownBackgroundGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/FeedDownBackgroundGraphs",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_FeedDown_feeddownOverTotal_AllCent[partIndex][rapIndex][pmIndex]);
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_FeedDown_feeddownOverTotal_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_FeedDown_Hybrid_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_FeedDown_Hybrid_ByCent_graphs[partIndex][centIndex][rapIndex][pmIndex]);
        }
      }

      outFile->mkdir(Form("%s/FeedDownBackgroundFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/FeedDownBackgroundFits",m_partInfo->GetParticleName(partIndex,charge).Data()));
      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_FeedDown_Fits_AllCent[partIndex][rapIndex][pmIndex]);
        //HistogramUtilities::ConditionalWrite(m_FeedDown_FitErrGraphs_AllCent[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_FeedDown_Fits_AllCent_Sys[partIndex][rapIndex][pmIndex]);
        //HistogramUtilities::ConditionalWrite(m_FeedDown_FitErrGraphs_AllCent_Sys[partIndex][rapIndex][pmIndex]);
        HistogramUtilities::ConditionalWrite(m_FeedDown_CovMat_AllCent[partIndex][rapIndex][pmIndex],Form("FeedDownCovMat_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));
        HistogramUtilities::ConditionalWrite(m_FeedDown_CovMat_AllCent_Sys[partIndex][rapIndex][pmIndex],Form("FeedDownCovMat_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),rapIndex));


      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
          HistogramUtilities::ConditionalWrite(m_FeedDown_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          //HistogramUtilities::ConditionalWrite(m_FeedDown_FitErrGraphs_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_FeedDown_Fits_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          //HistogramUtilities::ConditionalWrite(m_FeedDown_FitErrGraphs_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_FeedDown_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex],Form("FeedDownCovMat_%s_Cent%02d_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_FeedDown_CovMat_ByCent_Sys[partIndex][centIndex][rapIndex][pmIndex],Form("FeedDownCovMat_%s_Cent%02d_yIndex%d_Sys",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_FeedDown_Hybrid_Fits_ByCent[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_FeedDown_Hybrid_CovMat_ByCent[partIndex][centIndex][rapIndex][pmIndex],Form("FeedDownCovMat_Hybrid_%s_Cent%02d_yIndex%d",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));
          HistogramUtilities::ConditionalWrite(m_FeedDown_Hybrid_Fits_ByCent_lowPt[partIndex][centIndex][rapIndex][pmIndex]);
          HistogramUtilities::ConditionalWrite(m_FeedDown_Hybrid_CovMat_ByCent_lowPt[partIndex][centIndex][rapIndex][pmIndex],Form("FeedDownCovMat_Hybrid_%s_Cent%02d_yIndex%d_lowPt",m_partInfo->GetParticleName(partIndex,charge).Data(),centIndex,rapIndex));

        }
      }


      outFile->mkdir(Form("%s/Knockout",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/Knockout",m_partInfo->GetParticleName(partIndex,charge).Data()));
      HistogramUtilities::ConditionalWrite(m_proton_knockout_ChiSqrd);
      for(int centIndex = 0; centIndex < 16; centIndex++){
        HistogramUtilities::ConditionalWrite(m_proton_knockout_histogram_ByCent[centIndex]);
        HistogramUtilities::ConditionalWrite(m_proton_dca_ChiSqrd[centIndex]);
      }

      for(int rapIndex = 0; rapIndex < g_EffFit_MaxNumRapBins; rapIndex++){
        HistogramUtilities::ConditionalWrite(m_proton_knockout_AllCent[rapIndex]);
        HistogramUtilities::ConditionalWrite(m_proton_knockout_Fits_AllCent[rapIndex]);
        for(int centIndex = 0; centIndex < 16; centIndex++){
          HistogramUtilities::ConditionalWrite(m_proton_knockout_ByCent[centIndex][rapIndex]);
          HistogramUtilities::ConditionalWrite(m_proton_knockout_Fits_ByCent[centIndex][rapIndex]);
          HistogramUtilities::ConditionalWrite(m_proton_knockout_CovMat_ByCent[centIndex][rapIndex],Form("KnockoutCovMat_Cent%02d_yIndex%d",centIndex,rapIndex));
          for(int mTm0Index = 0; mTm0Index < g_EffFit_MaxNumMtM0Bins; mTm0Index++){
            HistogramUtilities::ConditionalWrite(m_proton_dca_Fits_ByCent[centIndex][rapIndex][mTm0Index]);
            if(centIndex == 0){
              HistogramUtilities::ConditionalWrite(m_proton_dca_Fits_AllCent[rapIndex][mTm0Index]);
            }
          }
        }
      }





      /*
      outFile->mkdir(Form("%s/EfficiencyGraphs2D",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EfficiencyGraphs2D",m_partInfo->GetParticleName(partIndex,charge).Data()));

      outFile->mkdir(Form("%s/EfficiencyGraphs2D",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/EfficiencyGraphs2D",m_partInfo->GetParticleName(partIndex,charge).Data()));
      */

    }
  }


  /*
    int charge = 1;
    for(int partIndex = 0; partIndex < 9; partIndex++){
      //if(!m_particleLoaded[partIndex]) continue;
      cout << " Part Index: " << partIndex << endl;
      outFile->mkdir(Form("%s/BTOF_Efficiency",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/BTOF_Efficiency",m_partInfo->GetParticleName(partIndex,charge).Data()));
      m_BTOFEff_Data_AllCent_Plus[partIndex]->Write();
      m_BTOFEff_Data_ChiSqrd_AllCent_Plus[partIndex]->Write();; // part rap
      for(int functIndex = 0; functIndex < 8; functIndex++){
        for(int paramIndex = 0; paramIndex < 6; paramIndex++){
          if(!HistogramUtilities::isEmpty(m_BTOFEff_Data_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex])) m_BTOFEff_Data_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->Write();
        }
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        cout << "centIndex: " << centIndex << endl;
        m_BTOFEff_Data_Plus[partIndex][centIndex]->Write();
        for(int functIndex = 0; functIndex < 8; functIndex++){
          for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty(m_BTOFEff_Data_FitParams_Plus[partIndex][centIndex][functIndex][paramIndex])) m_BTOFEff_Data_FitParams_Plus[partIndex][centIndex][functIndex][paramIndex]->Write(); //part cent fit param
          }
        }
      }

      m_BTOFEff_Emb_AllCent_Plus[partIndex]->Write();
      m_BTOFEff_Emb_ChiSqrd_AllCent_Plus[partIndex]->Write();; // part rap
      for(int functIndex = 0; functIndex < 8; functIndex++){
        for(int paramIndex = 0; paramIndex < 6; paramIndex++){
          if(!HistogramUtilities::isEmpty(m_BTOFEff_Emb_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex])) m_BTOFEff_Emb_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->Write();
        }
      }
      for(int centIndex = 0; centIndex < 16; centIndex++){
        cout << "centIndex: " << centIndex << endl;
        m_BTOFEff_Emb_Plus[partIndex][centIndex]->Write();
        for(int functIndex = 0; functIndex < 8; functIndex++){
          for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty(m_BTOFEff_Emb_FitParams_Plus[partIndex][centIndex][functIndex][paramIndex])) m_BTOFEff_Emb_FitParams_Plus[partIndex][centIndex][functIndex][paramIndex]->Write(); //part cent fit param
          }
        }
      }

    }




    charge = 0;
    for(int partIndex = 0; partIndex < 9; partIndex++){
      //if(!m_particleLoaded[partIndex]) continue;
      outFile->mkdir(Form("%s",m_partInfo->GetParticleName(partIndex).Data()));
      outFile->cd(Form("%s",m_partInfo->GetParticleName(partIndex).Data()));
      outFile->mkdir(Form("%s/BTOF_Efficiency",m_partInfo->GetParticleName(partIndex).Data()));
      outFile->cd(Form("%s/BTOF_Efficiency",m_partInfo->GetParticleName(partIndex).Data()));
      m_BTOFEff_Data_AllCent_Combined[partIndex]->Write();
      m_BTOFEff_Data_ChiSqrd_AllCent_Combined[partIndex]->Write();; // part rap

      for(int functIndex = 0; functIndex < 8; functIndex++){
        for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty(m_BTOFEff_Data_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex])) m_BTOFEff_Data_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->Write();
        }
      }

      for(int centIndex = 0; centIndex < 16; centIndex++){
        cout << "centIndex: " << centIndex << endl;
        m_BTOFEff_Data_Combined[partIndex][centIndex]->Write();
        for(int functIndex = 0; functIndex < 8; functIndex++){
          for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty( m_BTOFEff_Data_FitParams_Combined[partIndex][centIndex][functIndex][paramIndex])) m_BTOFEff_Data_FitParams_Combined[partIndex][centIndex][functIndex][paramIndex]->Write();
          }
        }
      }

      m_BTOFEff_Emb_AllCent_Combined[partIndex]->Write();
      m_BTOFEff_Emb_ChiSqrd_AllCent_Combined[partIndex]->Write();; // part rap

      for(int functIndex = 0; functIndex < 8; functIndex++){
        for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty(m_BTOFEff_Emb_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex])) m_BTOFEff_Emb_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->Write();
        }
      }

      for(int centIndex = 0; centIndex < 16; centIndex++){
        cout << "centIndex: " << centIndex << endl;
        m_BTOFEff_Emb_Combined[partIndex][centIndex]->Write();
        for(int functIndex = 0; functIndex < 8; functIndex++){
          for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty( m_BTOFEff_Emb_FitParams_Combined[partIndex][centIndex][functIndex][paramIndex])) m_BTOFEff_Emb_FitParams_Combined[partIndex][centIndex][functIndex][paramIndex]->Write();
          }
        }
      }

    }


    charge = -1;
    for(int partIndex = 0; partIndex < 9; partIndex++){
      //if(!m_particleLoaded[partIndex]) continue;
      outFile->mkdir(Form("%s/BTOF_Efficiency",m_partInfo->GetParticleName(partIndex,charge).Data()));
      outFile->cd(Form("%s/BTOF_Efficiency",m_partInfo->GetParticleName(partIndex,charge).Data()));
      m_BTOFEff_Data_AllCent_Minus[partIndex]->Write();
      m_BTOFEff_Data_ChiSqrd_AllCent_Minus[partIndex]->Write();; // part rap

      for(int functIndex = 0; functIndex < 8; functIndex++){
        for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty( m_BTOFEff_Data_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex])) m_BTOFEff_Data_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->Write();
        }
      }

      for(int centIndex = 0; centIndex < 16; centIndex++){
        cout << "centIndex: " << centIndex << endl;
        m_BTOFEff_Data_Minus[partIndex][centIndex]->Write();
        for(int functIndex = 0; functIndex < 8; functIndex++){
          for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty( m_BTOFEff_Data_FitParams_Minus[partIndex][centIndex][functIndex][paramIndex])) m_BTOFEff_Data_FitParams_Minus[partIndex][centIndex][functIndex][paramIndex]->Write();
          }
        }
      }

      m_BTOFEff_Emb_AllCent_Minus[partIndex]->Write();
      m_BTOFEff_Emb_ChiSqrd_AllCent_Minus[partIndex]->Write();; // part rap

      for(int functIndex = 0; functIndex < 8; functIndex++){
        for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty( m_BTOFEff_Emb_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex])) m_BTOFEff_Emb_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->Write();
        }
      }

      for(int centIndex = 0; centIndex < 16; centIndex++){
        cout << "centIndex: " << centIndex << endl;
        m_BTOFEff_Emb_Minus[partIndex][centIndex]->Write();
        for(int functIndex = 0; functIndex < 8; functIndex++){
          for(int paramIndex = 0; paramIndex < 6; paramIndex++){
            if(!HistogramUtilities::isEmpty( m_BTOFEff_Emb_FitParams_Minus[partIndex][centIndex][functIndex][paramIndex])) m_BTOFEff_Emb_FitParams_Minus[partIndex][centIndex][functIndex][paramIndex]->Write();
          }
        }
      }

    }
    */


  outFile->Close();
  for(unsigned int iii = 0; iii < m_splinesToDelete.size(); iii++){
    delete m_splinesToDelete[iii];
  }

  cout << "####################################################################################" << endl;
  cout << "                Search for \"ERROR:\" to make sure everything ran properly!  \n\n\n" << endl;
  cout << "                Search for \"Future Seeds\" to copy into your macro for better fitting \n\n\n" << endl;

  return;
}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double EfficiencyFitter::legendre_ChiSqr(const double* a_param){
  //TF1* m_functToFitMinimizer;
  //TGraphAsymmErrors* m_graphToFitMinimizer;

  //for()
  cout << "ERROR: NOT IMPLEMENTED: EfficiencyFitter::legendre_ChiSqr(const double* a_param)" << endl;
  return 1.0;


}

//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double EfficiencyFitter::getChiSquared_DCA(const double* a_params){
  double proton_anitproton_ratio = a_params[0];
  double amp_bkgd                = a_params[1];
  double DCA_naught              = a_params[2];
  double alpha                   = a_params[3];
  double chiSqr = 0.0;
  for(int binX = m_currentDCALowerFitRangeBin; binX <= m_currentDCAUpperFitRangeBin; binX++){
    double binCenter       = m_currentProtonDCA->GetBinCenter(binX);
    double binContentPlus  = m_currentProtonDCA->GetBinContent(binX);
    double binErrorPlus    = m_currentProtonDCA->GetBinError(binX);
    double binContentMinus = m_currentAntiProtonDCA->GetBinContent(binX);
    double binErrorMinus   = m_currentAntiProtonDCA->GetBinError(binX);
    double predictValue = proton_anitproton_ratio * binContentMinus + amp_bkgd * pow(1.0 - exp(-binCenter/DCA_naught),1.0/alpha);
    double newError = sqrt(binErrorPlus*binErrorPlus + proton_anitproton_ratio*proton_anitproton_ratio*binErrorMinus*binErrorMinus);
    //double newError = sqrt(binErrorPlus*binErrorPlus + binErrorMinus*binErrorMinus);
    //double newError = binErrorPlus;
    //double newError = 1.0;
    chiSqr += pow((binContentPlus - predictValue)/newError,2);
  }
  return chiSqr;
}





//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void EfficiencyFitter::compareTurnOnFunctions(){
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  for(int partIndex = 0; partIndex < 9; partIndex++){
    TAxis* rapidity_bin_axis = m_embedding_rapidity_structure[partIndex]->GetXaxis();
    int numRapBins = rapidity_bin_axis->GetNbins();
  	TH2D* histoChiPlus  = new TH2D(Form("functCompare_chiSqrd_%s", m_partInfo->GetParticleName(partIndex,1).Data()),
  		                           Form("#chi^{2} for Different Turn-On Functions %s",m_partInfo->GetParticleSymbol(partIndex,1).Data()),
  		                           rapidity_bin_axis->GetNbins(),
  		                           rapidity_bin_axis->GetBinLowEdge(1),
  		                           rapidity_bin_axis->GetBinLowEdge(rapidity_bin_axis->GetNbins() + 1),
  		                           8,0,8);
  	histoChiPlus->GetYaxis()->SetBinLabel(1,"Logistic Function");
  	histoChiPlus->GetYaxis()->SetBinLabel(2,"Generalized Logistic Function");
  	histoChiPlus->GetYaxis()->SetBinLabel(3,"TanH");
  	histoChiPlus->GetYaxis()->SetBinLabel(4,"ArcTan");
  	histoChiPlus->GetYaxis()->SetBinLabel(5,"Gundimannian");
  	histoChiPlus->GetYaxis()->SetBinLabel(6,"Erf");
  	histoChiPlus->GetYaxis()->SetBinLabel(7,"Algebraic 1");
  	histoChiPlus->GetYaxis()->SetBinLabel(8,"Algebraic 2");

  	TH2D* histoChiMinus  = new TH2D(Form("functCompare_chiSqrd_%s", m_partInfo->GetParticleName(partIndex,-1).Data()),
  		         Form("#chi^{2} for Different Turn-On Functions %s",m_partInfo->GetParticleSymbol(partIndex,-1).Data()),
  		         rapidity_bin_axis->GetNbins(),
  		         rapidity_bin_axis->GetBinLowEdge(1),
  		         rapidity_bin_axis->GetBinLowEdge(rapidity_bin_axis->GetNbins() + 1),
  		         8,0,8);
  	histoChiMinus->GetYaxis()->SetBinLabel(1,"Logistic Function");
  	histoChiMinus->GetYaxis()->SetBinLabel(2,"Generalized Logistic Function");
  	histoChiMinus->GetYaxis()->SetBinLabel(3,"TanH");
  	histoChiMinus->GetYaxis()->SetBinLabel(4,"ArcTan");
  	histoChiMinus->GetYaxis()->SetBinLabel(5,"Gundimannian");
  	histoChiMinus->GetYaxis()->SetBinLabel(6,"Erf");
  	histoChiMinus->GetYaxis()->SetBinLabel(7,"Algebraic 1");
  	histoChiMinus->GetYaxis()->SetBinLabel(8,"Algebraic 2");

    m_TPC_Efficiency_ChiSqrd_Compare[partIndex][0] = histoChiPlus;
    m_TPC_Efficiency_ChiSqrd_Compare[partIndex][1] = histoChiMinus;

  	if(numRapBins <= 0) numRapBins = 0;
    //if(!m_particleLoaded[partIndex]) continue;
    for(int pmIndex = 0; pmIndex < 2; pmIndex ++){
      int charge = 0;
      if(pmIndex == 0){
      	charge = 1;
      }else{
      	charge = -1;
      }



      for (int yIndex = 0; yIndex < numRapBins; yIndex++){
        TF1*  effFitFuncts[8];
        effFitFuncts[0] = PhysMath::getLogisticFunction();
        effFitFuncts[1] = PhysMath::getGeneralizedLogisticFunction();
        effFitFuncts[2] = PhysMath::getHyperbolicTangentFunction();
        effFitFuncts[3] = PhysMath::getArcTangentFunction();
        effFitFuncts[4] = PhysMath::getGundimannianFunction();
        effFitFuncts[5] = PhysMath::getErrorFunction();
        effFitFuncts[6] = PhysMath::getAlgebraicTurnOnFunction();
        effFitFuncts[7] = PhysMath::getAlgebraicTurnOnFunction_2();
        TGraphAsymmErrors* effGraph = NULL;

        effGraph = m_TPC_Efficiency_AllCent[partIndex][yIndex][pmIndex];

        if(!effGraph || effGraph->GetN() < 5) continue;

        for(int functIndex = 0; functIndex < 8; functIndex++){
          effFitFuncts[functIndex]->SetParLimits(0,0.1,5.0);
          effFitFuncts[functIndex]->SetParLimits(1,0.0001,20.0);
          effFitFuncts[functIndex]->FixParameter(2,0.0);
          effFitFuncts[functIndex]->FixParameter(3,0.0);
          //effFitFuncts[functIndex]->SetParLimits(2,-1.0,1.0);
          //effFitFuncts[functIndex]->SetParLimits(3,-5.0,5.0);

          TFitResultPtr resultPtr = effGraph->Fit(effFitFuncts[functIndex],"S EX0");
          TF1* fitResult = effGraph->GetFunction(effFitFuncts[functIndex]->GetName());
          double chiSqrd = effGraph->Chisquare(fitResult);
          double dof = effGraph->GetN() - fitResult->GetNpar();
          double chiSqrd_dof = chiSqrd / ((double) dof);


          m_TPC_Efficiency_ChiSqrd_Compare[partIndex][pmIndex]->SetBinContent(yIndex+1, functIndex + 1, chiSqrd_dof);


          fittingCanvas->cd();
          gPad->SetRightMargin(0.45);
          TH1F *frame = fittingCanvas->DrawFrame(0.0,0.0,4.0,1.0);
          frame->SetTitle(Form("%s Spectra y=[%.2f,%.2f] 0-80%% Central",m_partInfo->GetParticleSymbol(partIndex,charge).Data(),
                      rapidity_bin_axis->GetBinLowEdge(yIndex + 1), rapidity_bin_axis->GetBinLowEdge(yIndex + 2)));
          frame->GetXaxis()->SetTitle("m_{T}-m_{0} (GeV/c^{2})");
          frame->GetYaxis()->SetTitle("TPC Efficiency x Acceptance");
          //frame->GetYaxis()->SetRangeUser(0.0,1.0);
          effGraph->Draw("P SAME");
          gPad->Update();
          gSystem->ProcessEvents();

          fittingCanvas->SaveAs(Form("%s/%s/TPC_Eff_AllCent/RapIndex_%02d_%s.png",m_imageDir.c_str(),m_partInfo->GetParticleName(partIndex,charge).Data(),yIndex,fitResult->GetName() ));

        }
  		//for(int iii = 0; iii < 8; iii++) delete effFitFuncts[iii];

      } // end rapidity loop
    } // end plus minus loop
  } // end particle loop

}

//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

double EfficiencyFitter::btofEffOperator(double* a_x, const double* a_params){
  //if(a_x[0] > m_btofLowMtM0Ignore && a_x[0] < m_btofHighMtM0Ignore){
  //  return m_currentBTOF_Eff_ToFit->GetBinContent(m_currentBTOF_Eff_ToFit->FindBin(a_x[0]));
  //}
  return a_params[0]*exp(a_params[1]*pow(a_x[0],-a_params[3])) + a_params[2] + a_params[4]*a_x[0];
}

/*void EfficiencyFitter::makeTSpline3(string a_name, int a_numPoints, vector<double> a_x_knots, double a_low_deriv, double a_high_deriv){
  if(m_spline) delete m_spline;
  vector<double> yLocs;
  for(int iii = 0; iii < a_numPoints; iii++){
    yLocs.push_back(0.0);
  }
  m_spline = new TSpline3(a_name.c_str(),&a_x_knots[0],&yLocs[0],a_numPoints,"b1e1", a_low_deriv, a_high_deriv);
}*/



//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________




void EfficiencyFitter::fitBTOFEfficiency( bool a_useEmbedding, bool a_defineContaminationRegionToSkip,  int a_specificFitFunctionIndex){
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Starting EfficiencyFitter::fitBTOFEfficiency( bool " << a_useEmbedding << " , bool " <<  a_defineContaminationRegionToSkip  << " , int " <<  a_specificFitFunctionIndex << " )" << endl;
  #endif
  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  for(int partIndex = 0; partIndex < 9; partIndex++){
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "  Working on Particle  " << partIndex  << endl;
    #endif

    TAxis* rapidity_bin_axis = NULL;
    if(a_useEmbedding){
      rapidity_bin_axis = m_embedding_rapidity_structure[partIndex]->GetXaxis();
    }else{
      rapidity_bin_axis = m_picoBinner_rapidity_structure[partIndex][1]->GetXaxis();
    }


    int numRapBins = rapidity_bin_axis->GetNbins();
    if(numRapBins <= 0) numRapBins = 0;
    //if(!m_particleLoaded[partIndex]) continue;


    TH2D* histoChiPlus  = new TH2D(Form("functCompare_BTOFEff_Data_chiSqrd_%s", m_partInfo->GetParticleName(partIndex,1).Data()),
                                 Form("#chi^{2} for Different Turn-On Functions %s",m_partInfo->GetParticleSymbol(partIndex,1).Data()),
                                 rapidity_bin_axis->GetNbins(),
                                 rapidity_bin_axis->GetBinLowEdge(1),
                                 rapidity_bin_axis->GetBinLowEdge(rapidity_bin_axis->GetNbins() + 1),
                                 8,0,8);
    histoChiPlus->GetYaxis()->SetBinLabel(1,"Logistic Function");
    histoChiPlus->GetYaxis()->SetBinLabel(2,"Generalized Logistic Function");
    histoChiPlus->GetYaxis()->SetBinLabel(3,"TanH");
    histoChiPlus->GetYaxis()->SetBinLabel(4,"ArcTan");
    histoChiPlus->GetYaxis()->SetBinLabel(5,"Gundimannian");
    histoChiPlus->GetYaxis()->SetBinLabel(6,"Erf");
    histoChiPlus->GetYaxis()->SetBinLabel(7,"Algebraic 1");
    histoChiPlus->GetYaxis()->SetBinLabel(8,"Algebraic 2");


    TH2D* histoChiMinus  = new TH2D(Form("functCompare_BTOFEff_Data_chiSqrd_%s", m_partInfo->GetParticleName(partIndex,-1).Data()),
               Form("#chi^{2} for Different Turn-On Functions %s",m_partInfo->GetParticleSymbol(partIndex,-1).Data()),
               rapidity_bin_axis->GetNbins(),
               rapidity_bin_axis->GetBinLowEdge(1),
               rapidity_bin_axis->GetBinLowEdge(rapidity_bin_axis->GetNbins() + 1),
               8,0,8);
    histoChiMinus->GetYaxis()->SetBinLabel(1,"Logistic Function");
    histoChiMinus->GetYaxis()->SetBinLabel(2,"Generalized Logistic Function");
    histoChiMinus->GetYaxis()->SetBinLabel(3,"TanH");
    histoChiMinus->GetYaxis()->SetBinLabel(4,"ArcTan");
    histoChiMinus->GetYaxis()->SetBinLabel(5,"Gundimannian");
    histoChiMinus->GetYaxis()->SetBinLabel(6,"Erf");
    histoChiMinus->GetYaxis()->SetBinLabel(7,"Algebraic 1");
    histoChiMinus->GetYaxis()->SetBinLabel(8,"Algebraic 2");


    TH2D* histoChiCombined  = new TH2D(Form("functCompare_BTOFEff_Data_chiSqrd_%s", m_partInfo->GetParticleName(partIndex).Data()),
               Form("#chi^{2} for Different Turn-On Functions %s",m_partInfo->GetParticleSymbol(partIndex).Data()),
               rapidity_bin_axis->GetNbins(),
               rapidity_bin_axis->GetBinLowEdge(1),
               rapidity_bin_axis->GetBinLowEdge(rapidity_bin_axis->GetNbins() + 1),
               8,0,8);
    histoChiCombined->GetYaxis()->SetBinLabel(1,"Logistic Function");
    histoChiCombined->GetYaxis()->SetBinLabel(2,"Generalized Logistic Function");
    histoChiCombined->GetYaxis()->SetBinLabel(3,"TanH");
    histoChiCombined->GetYaxis()->SetBinLabel(4,"ArcTan");
    histoChiCombined->GetYaxis()->SetBinLabel(5,"Gundimannian");
    histoChiCombined->GetYaxis()->SetBinLabel(6,"Erf");
    histoChiCombined->GetYaxis()->SetBinLabel(7,"Algebraic 1");
    histoChiCombined->GetYaxis()->SetBinLabel(8,"Algebraic 2");

    m_BTOFEff_Data_ChiSqrd_AllCent_Plus[partIndex] = histoChiPlus;
    m_BTOFEff_Data_ChiSqrd_AllCent_Minus[partIndex] = histoChiMinus;
    m_BTOFEff_Data_ChiSqrd_AllCent_Combined[partIndex] = histoChiCombined;

    m_BTOFEff_Emb_ChiSqrd_AllCent_Plus[partIndex] = (TH2D*) histoChiPlus->Clone();
    m_BTOFEff_Emb_ChiSqrd_AllCent_Minus[partIndex] = (TH2D*) histoChiMinus->Clone();
    m_BTOFEff_Emb_ChiSqrd_AllCent_Combined[partIndex] = (TH2D*) histoChiMinus->Clone();
    m_BTOFEff_Emb_ChiSqrd_AllCent_Plus[partIndex]->SetName(Form("functCompare_BTOFEff_Emb_chiSqrd_%s", m_partInfo->GetParticleName(partIndex,1).Data()));
    m_BTOFEff_Emb_ChiSqrd_AllCent_Minus[partIndex]->SetName(Form("functCompare_BTOFEff_Emb_chiSqrd_%s", m_partInfo->GetParticleName(partIndex,-1).Data()));
    m_BTOFEff_Emb_ChiSqrd_AllCent_Combined[partIndex]->SetName(Form("functCompare_BTOFEff_Emb_chiSqrd_%s_Combined", m_partInfo->GetParticleName(partIndex,0).Data()));

    //############     FIT CENTRALITY INTEGRATED BTOF EFFICIENCIES ######################
    //############              PLUS, COMBINED, AND MINUS ###############################
    for(int charge = -1; charge < 2; charge++){
      for (int yIndex = 0; yIndex < numRapBins; yIndex++){
        double rapidity = m_BTOFEff_Data_AllCent_Plus[partIndex]->GetXaxis()->GetBinCenter(yIndex+1);
        TF1*  effFitFuncts[8];
        effFitFuncts[0] = PhysMath::getLogisticFunction();
        effFitFuncts[1] = PhysMath::getGeneralizedLogisticFunction();
        effFitFuncts[2] = PhysMath::getHyperbolicTangentFunction();
        effFitFuncts[3] = PhysMath::getArcTangentFunction();
        effFitFuncts[4] = PhysMath::getGundimannianFunction();
        effFitFuncts[5] = PhysMath::getErrorFunction();
        effFitFuncts[6] = PhysMath::getAlgebraicTurnOnFunction();
        effFitFuncts[7] = PhysMath::getAlgebraicTurnOnFunction_2();
        TH1D* effHisto = NULL;
        TH2D* chiSqrdHisto = NULL;
        if(charge == 1){
          if(a_useEmbedding){
            effHisto = m_BTOFEff_Emb_AllCent_Plus[partIndex]->ProjectionY("",yIndex+1,yIndex+1);
            chiSqrdHisto = m_BTOFEff_Emb_ChiSqrd_AllCent_Plus[partIndex];
          }else{
            effHisto = m_BTOFEff_Data_AllCent_Plus[partIndex]->ProjectionY("",yIndex+1,yIndex+1);
            chiSqrdHisto = m_BTOFEff_Data_ChiSqrd_AllCent_Plus[partIndex];
          }

        }else if(charge == -1){
          if(a_useEmbedding){
            effHisto = m_BTOFEff_Emb_AllCent_Minus[partIndex]->ProjectionY("",yIndex+1,yIndex+1);
            chiSqrdHisto = m_BTOFEff_Emb_ChiSqrd_AllCent_Minus[partIndex];
          }else{
            effHisto = m_BTOFEff_Data_AllCent_Minus[partIndex]->ProjectionY("",yIndex+1,yIndex+1);
            chiSqrdHisto = m_BTOFEff_Data_ChiSqrd_AllCent_Minus[partIndex];
          }
        }else{
          if(a_useEmbedding){
            effHisto = m_BTOFEff_Emb_AllCent_Combined[partIndex]->ProjectionY("",yIndex+1,yIndex+1);
            chiSqrdHisto = m_BTOFEff_Emb_ChiSqrd_AllCent_Combined[partIndex];
          }else{
            effHisto = m_BTOFEff_Data_AllCent_Combined[partIndex]->ProjectionY("",yIndex+1,yIndex+1);
            chiSqrdHisto = m_BTOFEff_Data_ChiSqrd_AllCent_Combined[partIndex];
          }
        }

        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "  BTOF Eff Fit:   charge = " << charge << "   yIndex = " << yIndex  << endl;
        #endif


        effHisto->SetName(Form("%s_btofEff_%02d",m_partInfo->GetParticleName(partIndex,charge).Data(), yIndex));
        if(effHisto->Integral() < 0.01){
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "  continuing because of low integral " << endl;
          #endif
          continue;
        }

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

        TLine* lowIndexLine = new TLine(lowIndexMtM0Cont,0.0,lowIndexMtM0Cont,1.0);
        TLine* highIndexLine = new TLine(highIndexMtM0Cont,0.0,highIndexMtM0Cont,1.0);
        lowIndexLine->SetLineColor(kViolet);
        highIndexLine->SetLineColor(kBlue);
        if(a_defineContaminationRegionToSkip) lowIndexLine->Draw("same");
        if(a_defineContaminationRegionToSkip) highIndexLine->Draw("same");

        TLine* lowRangeLine = new TLine(m_btofLowMtM0Ignore,0.0,m_btofLowMtM0Ignore,1.0);
        TLine* highRangeLine = new TLine(m_btofHighMtM0Ignore,0.0,m_btofHighMtM0Ignore,1.0);
        lowRangeLine->SetLineColor(kRed);
        highRangeLine->SetLineColor(kRed+1);
        if(a_defineContaminationRegionToSkip) lowRangeLine->Draw("same");
        if(a_defineContaminationRegionToSkip) highRangeLine->Draw("same");


        if(!a_defineContaminationRegionToSkip){
          m_btofHighMtM0Ignore = -1;
          m_btofLowMtM0Ignore = 999;
        }
        //cout << "Current Ignore region, mT-m0 != [" << m_btofLowMtM0Ignore << "," << m_btofHighMtM0Ignore << "]"<< endl;

        for(int functIndex = 0; functIndex < 8; functIndex++){
          if(a_specificFitFunctionIndex >= 0 && functIndex != a_specificFitFunctionIndex) continue;
          effFitFuncts[functIndex]->SetParLimits(0,0.0,5.0);
          effFitFuncts[functIndex]->SetParLimits(1,0.1,20.0);
          effFitFuncts[functIndex]->FixParameter(2,0.0);
          effFitFuncts[functIndex]->FixParameter(3,0.0);
          //effFitFuncts[functIndex]->SetParLimits(2,-1.0,1.0);
          //effFitFuncts[functIndex]->SetParLimits(3,-5.0,5.0);

          TFitResultPtr resultPtr = effHisto->Fit(effFitFuncts[functIndex],"S EX0");
          TF1* fitResult = effHisto->GetFunction(effFitFuncts[functIndex]->GetName());
          if(!fitResult){
            //continue;
            cout << "ERROR: No Fit Result, just using bad parameters" << endl;
            double chiSqrd = effHisto->Chisquare(effFitFuncts[functIndex]);
            double dof =  (double)(HistogramUtilities::numBinsWithContentAboveThreshold(effHisto, 0.0)) - (double)(effFitFuncts[functIndex]->GetNpar());
            double chiSqrd_dof = chiSqrd / ((double) dof);

            chiSqrdHisto->SetBinContent(yIndex+1, functIndex + 1, chiSqrd_dof);

            for(int paramIndex = 0; paramIndex < effFitFuncts[functIndex]->GetNpar(); paramIndex++){
              if(a_useEmbedding){
                if(charge == 1){
                  m_BTOFEff_Emb_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, effFitFuncts[functIndex]->GetParameter(paramIndex));
                  m_BTOFEff_Emb_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, effFitFuncts[functIndex]->GetParError(paramIndex));
                }
                if(charge == 0){
                  m_BTOFEff_Emb_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, effFitFuncts[functIndex]->GetParameter(paramIndex));
                  m_BTOFEff_Emb_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, effFitFuncts[functIndex]->GetParError(paramIndex));
                }
                if(charge == -1){
                  m_BTOFEff_Emb_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, effFitFuncts[functIndex]->GetParameter(paramIndex));
                  m_BTOFEff_Emb_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, effFitFuncts[functIndex]->GetParError(paramIndex));
                }
              }else{
                if(charge == 1){
                  m_BTOFEff_Data_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, effFitFuncts[functIndex]->GetParameter(paramIndex));
                  m_BTOFEff_Data_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, effFitFuncts[functIndex]->GetParError(paramIndex));
                }
                if(charge == 0){
                  m_BTOFEff_Data_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, effFitFuncts[functIndex]->GetParameter(paramIndex));
                  m_BTOFEff_Data_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, effFitFuncts[functIndex]->GetParError(paramIndex));
                }
                if(charge == -1){
                  m_BTOFEff_Data_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, effFitFuncts[functIndex]->GetParameter(paramIndex));
                  m_BTOFEff_Data_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, effFitFuncts[functIndex]->GetParError(paramIndex));
                }
              }

            }
          }else{
            double chiSqrd = effHisto->Chisquare(fitResult);
            double dof =  (double)(HistogramUtilities::numBinsWithContentAboveThreshold(effHisto, 0.0)) - (double)(fitResult->GetNpar());
            double chiSqrd_dof = chiSqrd / ((double) dof);
            chiSqrdHisto->SetBinContent(yIndex+1, functIndex + 1, chiSqrd_dof);

            for(int paramIndex = 0; paramIndex < fitResult->GetNpar(); paramIndex++){
              if(a_useEmbedding){
                if(charge == 1){
                  m_BTOFEff_Emb_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, fitResult->GetParameter(paramIndex));
                  m_BTOFEff_Emb_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, fitResult->GetParError(paramIndex));
                }
                if(charge == 0){
                  m_BTOFEff_Emb_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, fitResult->GetParameter(paramIndex));
                  m_BTOFEff_Emb_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, fitResult->GetParError(paramIndex));
                }
                if(charge == -1){
                  m_BTOFEff_Emb_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, fitResult->GetParameter(paramIndex));
                  m_BTOFEff_Emb_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, fitResult->GetParError(paramIndex));
                }
              }else{
                if(charge == 1){
                  m_BTOFEff_Data_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, fitResult->GetParameter(paramIndex));
                  m_BTOFEff_Data_FitParams_AllCent_Plus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, fitResult->GetParError(paramIndex));
                }
                if(charge == 0){
                  m_BTOFEff_Data_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, fitResult->GetParameter(paramIndex));
                  m_BTOFEff_Data_FitParams_AllCent_Combined[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, fitResult->GetParError(paramIndex));
                }
                if(charge == -1){
                  m_BTOFEff_Data_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinContent(yIndex+1, fitResult->GetParameter(paramIndex));
                  m_BTOFEff_Data_FitParams_AllCent_Minus[partIndex][functIndex][paramIndex]->SetBinError(yIndex+1, fitResult->GetParError(paramIndex));
                }
              }

            }
          }





          fittingCanvas->cd();
          gPad->SetRightMargin(0.45);
          effHisto->SetTitle(Form("%s y=[%.2f,%.2f] 0-80%% Central",m_partInfo->GetParticleSymbol(partIndex,charge).Data(),
                      rapidity_bin_axis->GetBinLowEdge(yIndex+1), rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          effHisto->GetXaxis()->SetTitle("m_{T} - m_{0} (GeV/c^{2})");
          effHisto->GetYaxis()->SetTitle("bTOF Efficiency");
          //frame->GetYaxis()->SetRangeUser(0.0,1.0);
          effHisto->Draw("E");
          fitResult->Draw("same");

          lowIndexLine->Draw("same");
          highIndexLine->Draw("same");
          //lowRangeLine->Draw("same");
          highRangeLine->Draw("same");



          gPad->Update();
          gSystem->ProcessEvents();

          if(a_useEmbedding){
            fittingCanvas->SaveAs(Form("%s/%s/BTOF_Eff_AllCent/EmbFit_RapIndex_%02d_%s.png",m_imageDir.c_str(),m_partInfo->GetParticleName(partIndex,charge).Data(),yIndex,fitResult->GetName() ));
          }else{
            fittingCanvas->SaveAs(Form("%s/%s/BTOF_Eff_AllCent/DataFit_RapIndex_%02d_%s.png",m_imageDir.c_str(),m_partInfo->GetParticleName(partIndex,charge).Data(),yIndex,fitResult->GetName() ));
          }

        }
        //for(int iii = 0; iii < 8; iii++) delete effFitFuncts[iii];
        //for(int centIndex = 0; centIndex < m_numCentBins; centIndex++){
        //}//end centrality loop



      } // end rapidity loop
    } // end charge loop
  } // end particle loop

}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________



/*  This function is legacy, just use the fitEfficiency function because BTOF and TPC are fit the same...
void EfficiencyFitter::fitTPCEff(int a_partIndex, int a_charge, int a_centIndex){//, int* a_numFitPoints){

  bool doSplineFitting = false;
  bool fixXLocationsOfSpline = false;
  bool removeFirstZerosForLegendreFits = true;

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);

  bool allCentMode = false;
  if(a_centIndex < 0) allCentMode = true;
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;

  cout << "############    TPC EFFICIENCY for " << m_partInfo->GetParticleName(a_partIndex,a_charge).Data() << " Cent " << a_centIndex << "   (-1 is centrality integrated)" << endl;

  gSystem->mkdir(Form("%s/%s/TPC_Efficiency",m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  if(allCentMode){
    gSystem->mkdir(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_AllCent",
      m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  }else{
    gSystem->mkdir(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_Cent%02d",
      m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex));
  }

  int numRapBins = m_cuts->getNRapidityBins(a_partIndex,1); // get number of tpc rapidity bins
  for(int yIndex = 0; yIndex < numRapBins; yIndex++){
    int rapBin = yIndex + 1;
    if(m_testingMidRapidity && rapBin != m_cuts->getRapMtM0BinStructure(a_partIndex,1)->FindBin(m_ycm)){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Continuing because yIndex: " << yIndex << " is not midRapidity" << endl;
      #endif
      continue;
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "######  Working on yIndex: " << yIndex <<  endl;
    #endif

    TGraphAsymmErrors* effGraph = NULL;
    if(allCentMode){
      effGraph = m_TPC_Efficiency_AllCent[a_partIndex][yIndex][pmIndex];
    }else{
      effGraph = m_TPC_Efficiency_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
    }
    //cout << " effGraph ptr: " << effGraph << endl;
    //cout << " effGraph name: " << effGraph->GetName() << endl;
    if(!effGraph){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " TPC Eff not fit for part: " << a_partIndex << " cent " << a_centIndex << "   yIndex: " << yIndex << "   ptr: " << effGraph << endl;
      #endif
      continue;
    }
    if(HistogramUtilities::getNumNonZeroPts(effGraph) < 5){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " TPC Eff not fit for part: " << a_partIndex << " cent " << a_centIndex << "   yIndex: " << yIndex << "   ptr: " << effGraph << endl;
        cout << "   Because it has fewer than 5 nonzero points." << endl;
      #endif
      continue;
    }


    //###################### Fit with Nominal Function #########################
    SplineFitter* fitter = NULL;
    TF1* fitFunction = NULL;
    TGraph* splinePoints = NULL;
    TF1* lineFitFunction = NULL;
    double chiSqr = 0.0;
    double ndof = 0.0;
    TMatrixDSym* legendreCovMat = NULL;
    TMatrixDSym* lineCovMat = NULL;
    int legendreStatus = -1;
    if(doSplineFitting){   // #########  fit with spline fitter #############
      fitter = new SplineFitter();
      fitter->setGraphToFit(effGraph);
      if(fixXLocationsOfSpline){
        int firstGoodPoint = HistogramUtilities::getFirstNonZeroPoint(effGraph);
        int lastGoodPoint = HistogramUtilities::getLastNonZeroPoint(effGraph);
        double xOfFirstGoodPoint = effGraph->GetX()[firstGoodPoint];
        int anotherPoint = HistogramUtilities::getPointClosestToVal(effGraph, xOfFirstGoodPoint + 0.5);
        int numCurvePoints = 6;
        if(anotherPoint - firstGoodPoint < 13) numCurvePoints = (anotherPoint - firstGoodPoint)/3;
        double stepForCurve = 0.5/((double) numCurvePoints);
        vector<double> fixedPoints;
        fixedPoints.push_back(xOfFirstGoodPoint);
        for(int iii = 1; iii <= numCurvePoints; iii++){
          fixedPoints.push_back(xOfFirstGoodPoint + ((double)iii)* stepForCurve);
        }
        vector< double > extraFixedPoints = HistogramUtilities::getEdgesByVerticalChange(effGraph, 4, true, fixedPoints[fixedPoints.size()-1],effGraph->GetX()[lastGoodPoint]);
        for(int jjj = 0; jjj < extraFixedPoints.size(); jjj++){
          fixedPoints.push_back(extraFixedPoints[jjj]);
        }
        fitter->analyzeGraph(fixedPoints.size(),fitFunction, splinePoints,chiSqr,ndof,false,false,&fixedPoints[0],false,0.0);
      }else{
        fitter->analyzeGraph(5,fitFunction, splinePoints,chiSqr,ndof,false,false,NULL,false,0.0);
      }
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Fit Funciton Ptr: " << fitFunction << endl;
      #endif
      if(fitFunction) fitFunction->SetNpx(300);
      if(fitFunction) fitFunction->SetLineColor(kViolet);
    }else{          //#########  FIT WITH LEGENDRE POLYNOMIAL AND LINEAR FIT At High mTm0  ##########
      //Legendre Fit
      fitFunction = PhysMath::getLegendreExpansion("legendreFunct", m_polynomialOrderForTPCEff[a_partIndex][pmIndex],false,0.0, 1.5);
      PhysMath::initializeLegendreExpansion(fitFunction, effGraph, 0.0, 1.5, false);
      TFitResultPtr resultPtr = effGraph->Fit(fitFunction,"EX0 S","",0,1.5);
      ndof =  HistogramUtilities::getNPointsInRange(effGraph,0,1.5) - fitFunction->GetNpar() - 2; //-2 becasue window bounds are fit parameters
      chiSqr = effGraph->Chisquare(fitFunction,"R EX0");
      cout << " issues? ptr: " << resultPtr.Get() << endl;
      legendreStatus = resultPtr->Status();
      cout << " status: " << legendreStatus << endl;
      legendreCovMat = (TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();
      fitFunction = (TF1*) effGraph->GetFunction("legendreFunctBtofEff")->Clone();

      //########   FIT HIGH PT WITH A LINE    #############
      int lowPoint = HistogramUtilities::getPointClosestToVal(effGraph, 1.5);
      int highPoint = HistogramUtilities::getPointClosestToVal(effGraph, 2.0);
      if(lowPoint == highPoint) cout << "ERROR: Line Eff Fit Failed, same point for high and low: " << lowPoint << endl;
      double lowY = effGraph->GetY()[lowPoint];
      double lowX = effGraph->GetX()[lowPoint];
      double highY = effGraph->GetY()[highPoint];
      double highX = effGraph->GetX()[highPoint];
      double deltaY = fabs(lowY-highY);
      double slopeGuess = (highY - lowY)/(highX - lowX);
      TF1* lineFunct = PhysMath::linearFunct("linearHighPtTPCEFF", lowY, slopeGuess, lowX);
      lineFunct->SetParLimits(0,lowY - 10.0*deltaY,lowY + 10.0*deltaY);
      if(slopeGuess > 0.0){
        lineFunct->SetParLimits(1,-5.0*slopeGuess,10.0*slopeGuess);
      }else{
        lineFunct->SetParLimits(1,10.0*slopeGuess,-5.0*slopeGuess);
      }
      lineFunct->SetRange(1.2,5.0);
      TFitResultPtr fitPtr = effGraph->Fit(lineFunct,"S EX0","",1.2,5.0);
      lineCovMat = (TMatrixDSym*) fitPtr->GetCovarianceMatrix().Clone();
      lineFitFunction = (TF1*) effGraph->GetFunction("linearHighPtTPCEFF")->Clone();
      lineFitFunction->SetRange(1.2,5.0);
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << " FitFunctPtr: " << fitFunction << endl;
      cout << " SplinePtr: " << splinePoints << endl;
      cout << " Chisquare: " << chiSqr << endl;
      cout << " NDOF: " << ndof << endl;
    #endif

    if(!fitFunction){
      cout << " ERROR: Spline Fitting Failed!!!!!" << endl;
      return;
    }
    fitFunction->SetNpx(300);

    if(allCentMode){
      m_TPC_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunction->Clone();
      m_TPC_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
      m_TPC_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                             m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
      if(doSplineFitting) m_TPC_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = fitter->getCovMatrix();
      else m_TPC_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = legendreCovMat;
      m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiSqr/ndof);

      if(lineFitFunction){
        m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineFitFunction;
        m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_HighPt_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s High p_{T} Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                             m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineCovMat;
      }

    }else{
      m_TPC_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunction->Clone();
      m_TPC_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
      m_TPC_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                             m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
      m_TPC_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = fitter->getCovMatrix();
      m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiSqr/ndof);
      if(doSplineFitting) m_TPC_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = fitter->getCovMatrix();
      else  m_TPC_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = legendreCovMat;

      if(lineFitFunction){
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineFitFunction;
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_HighPt_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s High p_{T} Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                             m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineCovMat;
      }


    }



    //#################     FIT THE GRAPH WITH SYSTEMATIC #####################
    TF1* sysFunct = new TF1("sysFunct","[0]*exp(-[1]*pow(x,-[2]))",0.0,5.0);
    sysFunct->SetParameter(0,0.5);
    sysFunct->SetParameter(1,0.15);
    sysFunct->SetParameter(2,1.0);
    sysFunct->SetParLimits(0,0,1.0);
    sysFunct->SetParLimits(1,0.0,50.0);
    sysFunct->SetParLimits(2,0.0,15);

    sysFunct->SetParNames("y_0","scale","exponent");
    sysFunct->SetNpx(300);
    sysFunct->SetLineColor(kBlue);
    string fitOptions = "E EX0 S";
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      fitOptions = fitOptions + " V ";
    #endif
    TFitResultPtr resultPtr = effGraph->Fit(sysFunct,fitOptions.c_str());
    TF1* fitFunctSys = effGraph->GetFunction("sysFunct");

    if(fitFunctSys){
      double chiPerNDF = resultPtr->Chi2() / resultPtr->Ndf();
      if(allCentMode){
        m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();;
        m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent Integrated y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();
        m_TPC_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiPerNDF);


      }else{
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();;
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent %d y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();
        m_TPC_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiPerNDF);
      }
    }



    fittingCanvas->cd();
    double lowYRange = 0;
    double highYRange = 0;
    HistogramUtilities::verticalDrawRange(effGraph, lowYRange,highYRange);
    TH1F *frame = fittingCanvas->DrawFrame(0.0,lowYRange,2.5,highYRange);
    frame->SetTitle(Form("TPC Efficiency %s Cent %1.2f - %1.2f %% y=[%1.2f,%1.2f]; m_{T}-m_{0} (GeV/c^{2}); TPC Efficiency",
          m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
          m_cuts->getLowCentPercent(a_centIndex),
          m_cuts->getHighCentPercent(a_centIndex),
          m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),
          m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
    frame->GetYaxis()->SetTitle("TPC Efficiency");
    frame->GetXaxis()->SetTitle("m_{T}-m_{0} (GeV/c^{2})");
    effGraph->Draw("PSame");
    fitFunction->Draw("same");
    if(splinePoints) splinePoints->Draw("PSame");
    if(lineFitFunction) lineFitFunction->Draw("PSame");

    gPad->Update();
    gSystem->ProcessEvents();
    if(allCentMode){
      fittingCanvas->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_AllCent/yIndex_%02d.png",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_Cent%02d/yIndex_%02d.png",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
    }

    TGraphErrors* confidenceGraph = NULL;
    if(allCentMode){
      m_TPC_Eff_FitErrGraphs_AllCent[a_partIndex][yIndex][pmIndex] = new TGraphErrors();
      m_TPC_Eff_FitErrGraphs_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_AllCent_yIndex%d_Conf",
            m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
      confidenceGraph = m_TPC_Eff_FitErrGraphs_AllCent[a_partIndex][yIndex][pmIndex];
    }else{
      m_TPC_Eff_FitErrGraphs_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = new TGraphErrors();
      m_TPC_Eff_FitErrGraphs_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d_Conf",
            m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
      confidenceGraph = m_TPC_Eff_FitErrGraphs_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
    }

    for(double iii = 0; iii < 5.0; iii+=0.01){
      int point = confidenceGraph->GetN();
      confidenceGraph->SetPoint(point,iii,fitFunction->Eval(iii));
      double integral = 0.0;
      double integralErr = 0.0;
      fitter->getIntegralAndError(iii-0.05, iii+0.05, integral, integralErr, true);
      confidenceGraph->SetPointError(point,0.0,integralErr);
    }
    m_splinesToDelete.push_back(fitter);
    if(splinePoints) delete splinePoints;
    //delete fitter;
    delete frame;

  } // end rapidity loop

}
*/



//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________




void EfficiencyFitter::fitEfficiency(int a_partIndex, int a_charge, int a_centIndex, bool a_doTPCEff, bool a_doDataDrivenBTOF, bool a_doLegendreResiduleFits){

  bool     doSplineFitting                 = false; // Legacy (this is SplineFitter, not TSpline3)
  bool     removeFirstZerosForLegendreFits = true;
  bool     useAllCentFitsAsSeeds           = true;
  bool     initiallyFitTSplineWithFixedHorLocs = false;
  bool     refitTSplineWithFixedHorLocs    = true;
  bool     doInitialFitWithHistogram       = !m_fixed_target_mode;
  double   whereToAddExtraEffPoints        = (m_fixed_target_mode ? -1.0 : 6.0);  // if < 0 it wont add extra points to constrain high pT end
  double   threshold_for_remove_zeros      = 0.001;





  TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
  fittingCanvas->SetWindowSize(1600,950);
  fittingCanvas->SetCanvasSize(1500,850);
  gPad->SetRightMargin(0.38);
  fittingCanvas->Draw();
  TPad* topPad_base    = new TPad("topPad_Eff_base","topPad_Eff",0.0,0.3,1.0,1.0);
  TPad* bottomPad_base = new TPad("bottomPad_Eff_base","bottomPad_Eff",0.0,0.0,1.0,0.3);
  topPad_base->Draw();
  gPad->SetRightMargin(0.38);
  bottomPad_base->Draw();
  gPad->SetRightMargin(0.38);



  double highTotalMtM0Range = 4.5;
  if(whereToAddExtraEffPoints > highTotalMtM0Range) highTotalMtM0Range = whereToAddExtraEffPoints + 0.3;

  double mass = m_partInfo->GetParticleMass(a_partIndex);

  bool allCentMode = false;
  if(a_centIndex < 0) allCentMode = true;
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;
  bool doBTOFEff = !a_doTPCEff;
  string modeName = "BTOF";
  if(a_doTPCEff) modeName = "TPC";
  TAxis* rapidity_bin_axis = NULL;
  if(doBTOFEff){
    cout << "############    BTOF EFFICIENCY for " << m_partInfo->GetParticleName(a_partIndex,a_charge).Data() << " Cent " << a_centIndex << "   (-1 is centrality integrated)" << endl;
    gSystem->mkdir(Form("%s/%s/BTOF_Efficiency",m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    if(allCentMode){
      gSystem->mkdir(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_AllCent",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    }else{
      gSystem->mkdir(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_Cent%02d",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex));
    }
    if(a_doDataDrivenBTOF){
      if(!m_BTOF_fromData_Efficiency_ChiSqrd[a_partIndex][pmIndex]){
        cout << "ERROR: BTOF (data) NOT PROPERLY LOADED!!!!!!!!!" << endl;
      }
      rapidity_bin_axis = m_picoBinner_rapidity_structure[a_partIndex][1]->GetXaxis();
    }else{
      if(!m_BTOF_Efficiency_ChiSqrd[a_partIndex][pmIndex]){
        cout << "ERROR: BTOF (emb) NOT PROPERLY LOADED!!!!!!!!!" << endl;
      }
      rapidity_bin_axis = m_embedding_rapidity_structure[a_partIndex]->GetXaxis();
    }
  }

  if(a_doTPCEff){
    cout << "############    TPC EFFICIENCY for " << m_partInfo->GetParticleName(a_partIndex,a_charge).Data() << " Cent " << a_centIndex << "   (-1 is centrality integrated)" << endl;
    gSystem->mkdir(Form("%s/%s/TPC_Efficiency",m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    if(allCentMode){
      gSystem->mkdir(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_AllCent",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    }else{
      gSystem->mkdir(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_Cent%02d",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex));
    }
    if(!m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]){
      cout << "ERROR: TPC Eff ChiSqrd NOT PROPERLY LOADED!!!!!!!!!" << endl;
    }
    rapidity_bin_axis = m_embedding_rapidity_structure[a_partIndex]->GetXaxis();

  }

  // partIndex pmIndex rapidity lowMtM0Range highMtM0Range xVal_1 xVal_2 ...
  vector< vector< double > > tspline3_initialization_array = {
    {2, 0,  0.0, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0,  0.1, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0,  0.2, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0,  0.3, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0,  0.4, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0,  0.5, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0,  0.6, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0,  0.7, 0.0, 10.0, 0.05, 0.07, 0.23, 0.5, 2.5},
    {2, 0,  0.8, 0.0, 10.0, 0.06, 0.09, 0.23, 0.5, 2.5},
    {2, 0,  0.9, 0.0, 10.0, 0.06, 0.10, 0.23, 0.5, 2.5},
    {2, 0,  1.0, 0.0, 10.0, 0.10, 0.15, 0.23, 0.5, 2.5},
    {2, 0,  1.1, 0.0, 10.0, 0.15, 0.20, 0.35, 0.6, 2.5},
    {2, 0,  1.2, 0.0, 10.0, 0.15, 0.25, 0.60, 1.5, 2.5},
    {2, 0,  1.3, 0.0, 10.0, 0.20, 0.35, 0.60, 1.5, 2.5},
    {2, 0,  1.4, 0.0, 10.0, 0.30, 0.5, 0.60, 1.5, 2.5},
    {2, 0,  1.5, 0.0, 10.0, 0.40, 0.6, 0.8, 1.5, 2.5},
    {2, 0,  1.6, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 0,  1.7, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 0,  1.8, 0.0, 10.0, 1.0, 1.6, 3.5},

    {2, 0, -0.0, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0, -0.1, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0, -0.2, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0, -0.3, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0, -0.4, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0, -0.5, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0, -0.6, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 0, -0.7, 0.0, 10.0, 0.05, 0.07, 0.23, 0.5, 2.5},
    {2, 0, -0.8, 0.0, 10.0, 0.06, 0.09, 0.23, 0.5, 2.5},
    {2, 0, -0.9, 0.0, 10.0, 0.06, 0.10, 0.23, 0.5, 2.5},
    {2, 0, -1.0, 0.0, 10.0, 0.10, 0.15, 0.23, 0.5, 2.5},
    {2, 0, -1.1, 0.0, 10.0, 0.15, 0.20, 0.35, 0.6, 2.5},
    {2, 0, -1.2, 0.0, 10.0, 0.15, 0.25, 0.60, 1.5, 2.5},
    {2, 0, -1.3, 0.0, 10.0, 0.20, 0.35, 0.60, 1.5, 2.5},
    {2, 0, -1.4, 0.0, 10.0, 0.30, 0.5, 0.60, 1.5, 2.5},
    {2, 0, -1.5, 0.0, 10.0, 0.40, 0.6, 0.8, 1.5, 2.5},
    {2, 0, -1.6, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 0, -1.7, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 0, -1.8, 0.0, 10.0, 1.0, 1.6, 3.5},


    {2, 1,  0.0, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1,  0.1, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1,  0.2, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1,  0.3, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1,  0.4, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1,  0.5, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1,  0.6, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1,  0.7, 0.0, 10.0, 0.05, 0.07, 0.23, 0.5, 2.5},
    {2, 1,  0.8, 0.0, 10.0, 0.06, 0.09, 0.23, 0.5, 2.5},
    {2, 1,  0.9, 0.0, 10.0, 0.06, 0.10, 0.23, 0.5, 2.5},
    {2, 1,  1.0, 0.0, 10.0, 0.10, 0.15, 0.23, 0.5, 2.5},
    {2, 1,  1.1, 0.0, 10.0, 0.15, 0.20, 0.35, 0.6, 2.5},
    {2, 1,  1.2, 0.0, 10.0, 0.15, 0.25, 0.60, 1.5, 2.5},
    {2, 1,  1.3, 0.0, 10.0, 0.20, 0.35, 0.60, 1.5, 2.5},
    {2, 1,  1.4, 0.0, 10.0, 0.30, 0.5, 0.60, 1.5, 2.5},
    {2, 1,  1.5, 0.0, 10.0, 0.40, 0.6, 0.8, 1.5, 2.5},
    {2, 1,  1.6, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 1,  1.7, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 1,  1.8, 0.0, 10.0, 1.0, 1.6, 3.5},

    {2, 1, -0.0, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1, -0.1, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1, -0.2, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1, -0.3, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1, -0.4, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1, -0.5, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1, -0.6, 0.0, 10.0, 0.03, 0.05, 0.14, 0.3, 0.75, 2.5},
    {2, 1, -0.7, 0.0, 10.0, 0.05, 0.07, 0.23, 0.5, 2.5},
    {2, 1, -0.8, 0.0, 10.0, 0.06, 0.09, 0.23, 0.5, 2.5},
    {2, 1, -0.9, 0.0, 10.0, 0.06, 0.10, 0.23, 0.5, 2.5},
    {2, 1, -1.0, 0.0, 10.0, 0.10, 0.15, 0.23, 0.5, 2.5},
    {2, 1, -1.1, 0.0, 10.0, 0.15, 0.20, 0.35, 0.6, 2.5},
    {2, 1, -1.2, 0.0, 10.0, 0.15, 0.25, 0.60, 1.5, 2.5},
    {2, 1, -1.3, 0.0, 10.0, 0.20, 0.35, 0.60, 1.5, 2.5},
    {2, 1, -1.4, 0.0, 10.0, 0.30, 0.5, 0.60, 1.5, 2.5},
    {2, 1, -1.5, 0.0, 10.0, 0.40, 0.6, 0.8, 1.5, 2.5},
    {2, 1, -1.6, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 1, -1.7, 0.0, 10.0, 1.0, 1.6, 3.5},
    {2, 1, -1.8, 0.0, 10.0, 1.0, 1.6, 3.5}

  };


  // partIndex pmIndex rapidity lowMtM0Range highMtM0Range xVal_1 xVal_2 ... (-999 at end means initially fix the x values)
  vector< vector< double > > tspline5_initialization_array = {

    // ______  PION PLUS  _______
    {0, 0,  0.0, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.1, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.2, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.3, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.4, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.5, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.6, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.7, 0.0, 10.0, 0.1, 1.5 },
    {0, 0,  0.8, 0.0, 10.0, 0.1, 1.2 },
    {0, 0,  0.9, 0.0, 10.0, 0.1, 1.2 },
    {0, 0,  1.0, 0.0, 10.0, 0.1, 0.3 },
    {0, 0,  1.1, 0.0, 10.0, 0.1, 0.3 },
    {0, 0,  1.2, 0.0, 10.0, 0.15, 0.5 },
    {0, 0,  1.3, 0.0, 10.0, 0.15, 0.6 },
    {0, 0,  1.4, 0.0, 10.0, 0.2, 0.7 },
    {0, 0,  1.5, 0.0, 10.0, 0.2, 0.8 },
    {0, 0,  1.6, 0.0, 10.0, 0.4, 1.0 },
    {0, 0,  1.7, 0.0, 10.0, 0.5, 1.2 },
    {0, 0,  1.8, 0.0, 10.0, 2, 3},

    {0, 0, -0.0, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.1, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.2, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.3, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.4, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.5, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.6, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.7, 0.0, 10.0, 0.1, 1.5 },
    {0, 0, -0.8, 0.0, 10.0, 0.1, 1.2 },
    {0, 0, -0.9, 0.0, 10.0, 0.1, 1.2 },
    {0, 0, -1.0, 0.0, 10.0, 0.1, 0.3 },
    {0, 0, -1.1, 0.0, 10.0, 0.1, 0.3 },
    {0, 0, -1.2, 0.0, 10.0, 0.15, 0.5 },
    {0, 0, -1.3, 0.0, 10.0, 0.15, 0.6 },
    {0, 0, -1.4, 0.0, 10.0, 0.2, 0.7 },
    {0, 0, -1.5, 0.0, 10.0, 0.2, 0.8 },
    {0, 0, -1.6, 0.0, 10.0, 0.4, 1.0 },
    {0, 0, -1.7, 0.0, 10.0, 0.5, 1.2 },
    {0, 0, -1.8, 0.0, 10.0, 2, 3},

    // ______  PION MINUS  _______
    {0, 1,  0.0, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  0.1, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  0.2, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  0.3, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  0.4, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  0.5, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  0.6, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  0.7, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  0.8, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  0.9, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  1.0, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  1.1, 0.0, 10.0, 0.4, 2.5},
    {0, 1,  1.2, 0.0, 10.0, 0.4, 2.5},
    {0, 1,  1.3, 0.0, 10.0, 0.6, 2.5},
    {0, 1,  1.4, 0.0, 10.0, 1.2, 2.5},
    {0, 1,  1.5, 0.0, 10.0, 1.2, 2.5},
    {0, 1,  1.6, 0.0, 10.0, 2, 3},
    {0, 1,  1.7, 0.0, 10.0, 2, 3},
    {0, 1,  1.8, 0.0, 10.0, 2, 3},

    {0, 1,  -0.0, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  -0.1, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  -0.2, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  -0.3, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  -0.4, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  -0.5, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  -0.6, 0.0, 10.0, 0.1, 2.5},
    {0, 1,  -0.7, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  -0.8, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  -0.9, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  -1.0, 0.0, 10.0, 0.2, 2.5},
    {0, 1,  -1.1, 0.0, 10.0, 0.4, 2.5},
    {0, 1,  -1.2, 0.0, 10.0, 0.4, 2.5},
    {0, 1,  -1.3, 0.0, 10.0, 0.6, 2.5},
    {0, 1,  -1.4, 0.0, 10.0, 1.2, 2.5},
    {0, 1,  -1.5, 0.0, 10.0, 1.2, 2.5},
    {0, 1,  -1.6, 0.0, 10.0, 2, 3},
    {0, 1,  -1.7, 0.0, 10.0, 2, 3},
    {0, 1,  -1.8, 0.0, 10.0, 2, 3},


    // ______  KAON PLUS  _______
    {1, 0,  0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  0.7, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  0.8, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  0.9, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  1.0, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  1.1, 0.0, 10.0, 0.4, 2.5},
    {1, 0,  1.2, 0.0, 10.0, 0.4, 2.5},
    {1, 0,  1.3, 0.0, 10.0, 0.6, 2.5},
    {1, 0,  1.4, 0.0, 10.0, 1.2, 2.5},
    {1, 0,  1.5, 0.0, 10.0, 1.2, 2.5},
    {1, 0,  1.6, 0.0, 10.0, 2, 3},
    {1, 0,  1.7, 0.0, 10.0, 2, 3},
    {1, 0,  1.8, 0.0, 10.0, 2, 3},

    {1, 0,  -0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  -0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  -0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  -0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  -0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  -0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  -0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 0,  -0.7, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  -0.8, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  -0.9, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  -1.0, 0.0, 10.0, 0.2, 2.5},
    {1, 0,  -1.1, 0.0, 10.0, 0.4, 2.5},
    {1, 0,  -1.2, 0.0, 10.0, 0.4, 2.5},
    {1, 0,  -1.3, 0.0, 10.0, 0.6, 2.5},
    {1, 0,  -1.4, 0.0, 10.0, 1.2, 2.5},
    {1, 0,  -1.5, 0.0, 10.0, 1.2, 2.5},
    {1, 0,  -1.6, 0.0, 10.0, 2, 3},
    {1, 0,  -1.7, 0.0, 10.0, 2, 3},
    {1, 0,  -1.8, 0.0, 10.0, 2, 3},

    // ______  KAON MINUS  _______
    {1, 1,  0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  0.7, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  0.8, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  0.9, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  1.0, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  1.1, 0.0, 10.0, 0.4, 2.5},
    {1, 1,  1.2, 0.0, 10.0, 0.4, 2.5},
    {1, 1,  1.3, 0.0, 10.0, 0.6, 2.5},
    {1, 1,  1.4, 0.0, 10.0, 1.2, 2.5},
    {1, 1,  1.5, 0.0, 10.0, 1.2, 2.5},
    {1, 1,  1.6, 0.0, 10.0, 2, 3},
    {1, 1,  1.7, 0.0, 10.0, 2, 3},
    {1, 1,  1.8, 0.0, 10.0, 2, 3},

    {1, 1,  -0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  -0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  -0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  -0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  -0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  -0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  -0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {1, 1,  -0.7, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  -0.8, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  -0.9, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  -1.0, 0.0, 10.0, 0.2, 2.5},
    {1, 1,  -1.1, 0.0, 10.0, 0.4, 2.5},
    {1, 1,  -1.2, 0.0, 10.0, 0.4, 2.5},
    {1, 1,  -1.3, 0.0, 10.0, 0.6, 2.5},
    {1, 1,  -1.4, 0.0, 10.0, 1.2, 2.5},
    {1, 1,  -1.5, 0.0, 10.0, 1.2, 2.5},
    {1, 1,  -1.6, 0.0, 10.0, 2, 3},
    {1, 1,  -1.7, 0.0, 10.0, 2, 3},
    {1, 1,  -1.8, 0.0, 10.0, 2, 3},

    // ______  PROTON PLUS  _______
    {2, 0,  0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  0.7, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  0.8, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  0.9, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  1.0, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  1.1, 0.0, 10.0, 0.4, 2.5},
    {2, 0,  1.2, 0.0, 10.0, 0.4, 2.5},
    {2, 0,  1.3, 0.0, 10.0, 0.6, 2.5},
    {2, 0,  1.4, 0.0, 10.0, 1.2, 2.5},
    {2, 0,  1.5, 0.0, 10.0, 1.2, 2.5},
    {2, 0,  1.6, 0.0, 10.0, 2, 3},
    {2, 0,  1.7, 0.0, 10.0, 2, 3},
    {2, 0,  1.8, 0.0, 10.0, 2, 3},

    {2, 0,  -0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  -0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  -0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  -0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  -0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  -0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  -0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 0,  -0.7, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  -0.8, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  -0.9, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  -1.0, 0.0, 10.0, 0.2, 2.5},
    {2, 0,  -1.1, 0.0, 10.0, 0.4, 2.5},
    {2, 0,  -1.2, 0.0, 10.0, 0.4, 2.5},
    {2, 0,  -1.3, 0.0, 10.0, 0.6, 2.5},
    {2, 0,  -1.4, 0.0, 10.0, 1.2, 2.5},
    {2, 0,  -1.5, 0.0, 10.0, 1.2, 2.5},
    {2, 0,  -1.6, 0.0, 10.0, 2, 3},
    {2, 0,  -1.7, 0.0, 10.0, 2, 3},
    {2, 0,  -1.8, 0.0, 10.0, 2, 3},

    // ______  PROTON MINUS  _______
    {2, 1,  0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  0.7, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  0.8, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  0.9, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  1.0, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  1.1, 0.0, 10.0, 0.4, 2.5},
    {2, 1,  1.2, 0.0, 10.0, 0.4, 2.5},
    {2, 1,  1.3, 0.0, 10.0, 0.6, 2.5},
    {2, 1,  1.4, 0.0, 10.0, 1.2, 2.5},
    {2, 1,  1.5, 0.0, 10.0, 1.2, 2.5},
    {2, 1,  1.6, 0.0, 10.0, 2, 3},
    {2, 1,  1.7, 0.0, 10.0, 2, 3},
    {2, 1,  1.8, 0.0, 10.0, 2, 3},

    {2, 1,  -0.0, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  -0.1, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  -0.2, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  -0.3, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  -0.4, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  -0.5, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  -0.6, 0.0, 10.0, 0.1, 1.0, 2.5},
    {2, 1,  -0.7, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  -0.8, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  -0.9, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  -1.0, 0.0, 10.0, 0.2, 2.5},
    {2, 1,  -1.1, 0.0, 10.0, 0.4, 2.5},
    {2, 1,  -1.2, 0.0, 10.0, 0.4, 2.5},
    {2, 1,  -1.3, 0.0, 10.0, 0.6, 2.5},
    {2, 1,  -1.4, 0.0, 10.0, 1.2, 2.5},
    {2, 1,  -1.5, 0.0, 10.0, 1.2, 2.5},
    {2, 1,  -1.6, 0.0, 10.0, 2, 3},
    {2, 1,  -1.7, 0.0, 10.0, 2, 3},
    {2, 1,  -1.8, 0.0, 10.0, 2, 3}

  };



  vector<vector<double>> paramStorageForFutureSeeds;

  int polynomialOrder            = m_polynomialOrderForBTOFEff[a_partIndex][pmIndex];
  if(a_doTPCEff) polynomialOrder = m_polynomialOrderForTPCEff[a_partIndex][pmIndex];
  int numRapBins = rapidity_bin_axis->GetNbins();
  //int numRapBins = m_cuts->getNRapidityBins(a_partIndex,1); // get number of btof rapidity bins
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << " Number of rapidity bins: " << numRapBins << endl;
  #endif
  //################    RAPIDITY LOOP   ##########################
  for(int yIndex = 0; yIndex < numRapBins; yIndex++){  // ---------------------------------------------   RAPIDITY LOOP
    int    rapBin   = yIndex + 1;
    double rapidity = rapidity_bin_axis->GetBinCenter(rapBin);
    if(m_testingMidRapidity && rapBin != rapidity_bin_axis->FindBin(m_ycm)){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Continuing because yIndex: " << yIndex << " is not midRapidity" << endl;
      #endif
      continue;
    }
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - rapidity) < 0.005){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating rapidity = " << rapidity << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "######  Working on yIndex: " << yIndex <<  endl;
    #endif


    // Determine which initializion vector to use for the tspline
    int tspline3_initialization_array_index = -1;
    for(int iii = 0; iii < tspline3_initialization_array.size(); iii++){
      if(tspline3_initialization_array[iii].size() < 6){
        cout << "ERROR: something is wrong with the tspline3_initialization_array!" << endl;
      }
      int particle_index = (int) tspline3_initialization_array[iii][0];
      int plus_minus_index = (int) tspline3_initialization_array[iii][1];
      int charge = (plus_minus_index == 0) ? 1 : -1;
      double rapidity_init = tspline3_initialization_array[iii][2];
      if(particle_index != a_partIndex || charge != a_charge || fabs(rapidity - rapidity_init) > 0.01){
        continue;
      }else{
        tspline3_initialization_array_index = iii;
        break;
      }
    }
    #ifdef _EFFICIENCY_FITTER_DEBUG_
    if(tspline3_initialization_array_index < 0){
      cout << "WARNING: No TSpline3 initialization for this particle and rapidity" << endl;
    }
    #endif

    int tspline5_initialization_array_index = -1;
    for(int iii = 0; iii < tspline5_initialization_array.size(); iii++){
      if(tspline5_initialization_array[iii].size() < 6){
        cout << "ERROR: something is wrong with the tspline5_initialization_array!" << endl;
      }
      int particle_index   = (int) tspline5_initialization_array[iii][0];
      int plus_minus_index = (int) tspline5_initialization_array[iii][1];
      int charge = (plus_minus_index == 0) ? 1 : -1;
      double rapidity_init = tspline5_initialization_array[iii][2];
      if(particle_index != a_partIndex || charge != a_charge || fabs(rapidity - rapidity_init) > 0.01){
        continue;
      }else{
        tspline5_initialization_array_index = iii;
        break;
      }
    }
    #ifdef _EFFICIENCY_FITTER_DEBUG_
    if(tspline5_initialization_array_index < 0){
      cout << "WARNING: No TSpline3 initialization for this particle and rapidity" << endl;
    }
    #endif





    TGraphAsymmErrors* effGraph = NULL; // THE EFFICIENCY GRAPH
    TH1D*              effHisto = NULL; // Efficiency Histogram (Approximation due to symmetric error bars)
    if(doBTOFEff){
      if(allCentMode){
        if(a_doDataDrivenBTOF) effGraph = m_BTOF_fromData_Efficiency_AllCent[a_partIndex][yIndex][pmIndex];
        else                   effGraph = m_BTOF_Efficiency_AllCent[a_partIndex][yIndex][pmIndex];
      }else{
        if(a_doDataDrivenBTOF) effGraph = m_BTOF_fromData_Efficiency_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
        else                   effGraph = m_BTOF_Efficiency_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
      }
    }
    if(a_doTPCEff){
      if(allCentMode){
        effGraph = m_TPC_Efficiency_AllCent[a_partIndex][yIndex][pmIndex];
      }else{
        effGraph = m_TPC_Efficiency_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
      }
    }
    if(!effGraph || HistogramUtilities::getNumNonZeroPts(effGraph) < 4){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << modeName << " Eff not fit for part: " << a_partIndex << " cent " << a_centIndex << "   yIndex: " << yIndex << "   ptr: " << effGraph << endl;
        if(effGraph) cout << "     getNumNonZeroPts=" <<  HistogramUtilities::getNumNonZeroPts(effGraph) << endl;
      #endif
      continue;
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "######    First point of eff graph: x: " << effGraph->GetX()[0] << " y: " << effGraph->GetY()[0] << "  eyl: " << effGraph->GetEYlow()[0] << " eyh: " << effGraph->GetEYhigh()[0] << endl;
    #endif

    effHisto = HistogramUtilities::makeHistogram_Asymm(effGraph);




    /*TF1* kaonDecayFunction = NULL;
    if(a_partIndex == 1){
      kaonDecayFunction = new TF1("kaonDecayFunct","PhysMath::fractionRemainingAfterDistance([0],x,[1],[2],[3],[4])");
      kaonDecayFunction->SetParNames("rapidity","mass","tau","R","L");
      kaonDecayFunction->SetParameter(0,rapidity);
      kaonDecayFunction->SetParameter(1,mass);
      kaonDecayFunction->SetParameter(2,1.238E-8);
      kaonDecayFunction->SetParameter(3,0.1);
      if(m_fixed_target_mode){
        kaonDecayFunction->SetParameter(4,4.0);
      }else{
        kaonDecayFunction->SetParameter(4,2.0);
      }
      for(int pointIndex = 0; pointIndex < effGraph->GetN(); pointIndex++){
        double decay_scale = kaonDecayFunction->Eval(effGraph->GetX()[pointIndex]);
        effGraph->SetPoint(pointIndex,effGraph->GetX()[pointIndex],effGraph->GetY()[pointIndex]/decay_scale);
        effGraph->SetPointError(pointIndex,effGraph->GetEXlow()[pointIndex],effGraph->GetEXhigh()[pointIndex],
                  effGraph->GetEYlow()[pointIndex]/decay_scale,effGraph->GetEYhigh()[pointIndex]/decay_scale);
      }

    }*/


    //###################### Fit with Spline Fitter #########################
    SplineFitter* sp_fitter              = NULL;
    TF1*          fitFunction         = NULL;
    TF1*          lineFitFunction     = NULL;
    TF1*          initialLegendreSeed = NULL;
    TGraph*       splinePoints        = NULL;
    double chiSqr = 0.0;
    double ndof   = 0.0;
    //TFitResultPtr* fitPtr;
    TMatrixDSym* legendreCovMat = NULL;
    TMatrixDSym* lineCovMat     = NULL;
    int legendreStatus = -1;
    //    we now fit residules, so this shouldn't matter
    if(doSplineFitting){
      sp_fitter = new SplineFitter();
      sp_fitter->setGraphToFit(effGraph);
      sp_fitter->analyzeGraph(4,fitFunction, splinePoints,chiSqr,ndof,false,false,NULL,false,0.0);
      if(fitFunction) fitFunction->SetNpx(300);
      if(fitFunction) fitFunction->SetLineColor(kViolet);

    }else{ //#########  FIT WITH LEGENDRE POLYNOMIAL AND LINEAR FIT At High mTm0  ##########
      //Legendre Fit
      ////if(3*polynomialOrder > effGraph->GetN()*2) polynomialOrder = effGraph->GetN()*2/3;
      //double lowRange = 0.0;
      //if(removeFirstZerosForLegendreFits){
      //  for(int pointIndex = 0; pointIndex < effGraph->GetN(); pointIndex++){
      //    if(effGraph->GetX()[pointIndex] < 1.0 && effGraph->GetY()[pointIndex] < threshold_for_remove_zeros) lowRange = effGraph->GetX()[pointIndex] + 0.5*effGraph->GetEXhigh()[pointIndex];
      //  }
      //}
      double mTm0CutOffLine = m_BTOF_highMtM0ToLine[a_partIndex]->Eval(fabs(rapidity));
      double mTm0LegendreCutOffForFit = mTm0CutOffLine+0.5;
      //if the mTm0 matrix is populated, set the lower limit that way (by hand in macro)
      if(a_doTPCEff && m_TPC_Eff_highMtM0ToLine_Tuple.size() > 0){
        double toleranceForSeedRap = 0.01;
        for(int seedIndex = 0; seedIndex < m_TPC_Eff_highMtM0ToLine_Tuple.size(); seedIndex++){
          if(a_partIndex == m_TPC_Eff_highMtM0ToLine_Tuple[seedIndex][0] && fabs(m_TPC_Eff_highMtM0ToLine_Tuple[seedIndex][1] - rapidity) < toleranceForSeedRap){
            mTm0CutOffLine = m_TPC_Eff_highMtM0ToLine_Tuple[seedIndex][2];
            #ifdef _EFFICIENCY_FITTER_DEBUG_
              cout << "  Using seed matrix : " << m_TPC_Eff_highMtM0ToLine_Tuple[seedIndex][0] << " " << m_TPC_Eff_highMtM0ToLine_Tuple[seedIndex][1] << " " << m_TPC_Eff_highMtM0ToLine_Tuple[seedIndex][2] << endl;
            #endif
            break;
          }
        }
      }
      //fitFunction = PhysMath::getLegendreExpansion("legendreFunctEff", polynomialOrder,false,lowRange, mTm0LegendreCutOffForFit);
      //PhysMath::initializeLegendreExpansion(fitFunction, effGraph, lowRange, mTm0LegendreCutOffForFit, false,2000); // fixXRange? num Iterations
      ////fitFunction->SetParameter(fitFunction->GetNpar()-2,0.0);
      ////fitFunction->SetParameter(fitFunction->GetNpar()-1,mTm0LegendreCutOffForFit);
      //initialLegendreSeed = (TF1*) fitFunction->Clone();
      //initialLegendreSeed->SetName("initialLegendreSeed");
      //initialLegendreSeed->SetLineColor(kYellow);
      //TFitResultPtr fitPtr = effGraph->Fit(fitFunction,"S EX0","",lowRange,mTm0LegendreCutOffForFit);
      //ndof =  HistogramUtilities::getNPointsInRange(effGraph,lowRange,mTm0LegendreCutOffForFit) - fitFunction->GetNpar() - 2; //-2 becasue window bounds are fit parameters
      //chiSqr = effGraph->Chisquare(fitFunction,"R EX0");
      //cout << " issues? ptr: " << fitPtr.Get() << endl;
      //legendreStatus = fitPtr->Status();
      //cout << " status: " << legendreStatus << endl;
      //legendreCovMat = (TMatrixDSym*) fitPtr->GetCovarianceMatrix().Clone();
      //fitFunction = (TF1*) effGraph->GetFunction("legendreFunctEff")->Clone();
      //
      //########   FIT HIGH PT WITH A LINE    #############

      int lowPoint  = HistogramUtilities::getPointClosestToVal(effGraph, mTm0CutOffLine);
      int highPoint = HistogramUtilities::getPointClosestToVal(effGraph, highTotalMtM0Range);
      bool enoughToLineFit = true;
      if(lowPoint == highPoint){
        cout << "ERROR: Line Eff Fit Failed, same point for high and low: " << lowPoint << endl;
        enoughToLineFit = false;
      }else{
        double lowY = effGraph->GetY()[lowPoint];
        double lowX = effGraph->GetX()[lowPoint];
        double highY = effGraph->GetY()[highPoint];
        double highX = effGraph->GetX()[highPoint];
        double deltaY = fabs(lowY-highY);
        double slopeGuess = (highY - lowY)/(highX - lowX);
        TF1* lineFunct = PhysMath::linearFunct("linearHighPtBTOFEFFDD", lowY, slopeGuess, lowX);
        lineFunct->SetParLimits(0,lowY - 10.0*deltaY,lowY + 10.0*deltaY);
        if(slopeGuess > 0.0){
          lineFunct->SetParLimits(1,-5.0*slopeGuess,10.0*slopeGuess);
        }else{
          lineFunct->SetParLimits(1,10.0*slopeGuess,-5.0*slopeGuess);
        }

        lineFunct->SetRange(mTm0CutOffLine,5.0);
        //fitFunction->SetRange(lowRange,mTm0CutOffLine);
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " Linear Fit: Low Point: " << lowPoint << "(" << lowX << "," << lowY << ")" << "  High Point: " << highPoint << "(" << highX << "," << highY << ")" << endl;
          cout << "       slope guess: " << slopeGuess << " big delta y max : " << deltaY << endl;
        #endif
        lineFunct->FixParameter(1,0.0);
        TFitResultPtr lineFitPtr = effGraph->Fit(lineFunct,"S EX0","",mTm0CutOffLine,highTotalMtM0Range);
        lineCovMat = (TMatrixDSym*) lineFitPtr->GetCovarianceMatrix().Clone();
        lineFitFunction = (TF1*) effGraph->GetFunction("linearHighPtBTOFEFFDD")->Clone();
        //lineFitFunction->SetRange(1.2,5.0);
        lineFitFunction->SetLineColor(kRed+1);
      }
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << " FitFunctPtr: " << fitFunction << " SplinePtr: " << splinePoints << "  HighPtLinePtr: " << lineFitFunction << endl;
      cout << " Chisquare: " << chiSqr << " NDOF: " << ndof << endl;
    #endif

    double chiSqrPerDOF = chiSqr/ndof;

    //if(!fitFunction){
    //  cout << " ERROR: Nominal Fitting for BTOF Eff Failed!!!!!" << endl;
    //  return;
    //}
    //fitFunction->SetNpx(300);

    //######################     SAVE THE LINEAR FIT   ########################
    if(doBTOFEff && a_doDataDrivenBTOF && lineFitFunction){
      if(allCentMode){
        m_BTOF_fromData_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineFitFunction;
        m_BTOF_fromData_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_HighPt_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_BTOF_fromData_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s High p_{T} Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                             rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_BTOF_fromData_Eff_CovMat_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineCovMat;
      }else{
        m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineFitFunction;
        m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_HighPt_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s High p_{T}Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                             rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_BTOF_fromData_Eff_CovMat_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineCovMat;
      }
    }
    /*if(doBTOFEff && !a_doDataDrivenBTOF){
      if(allCentMode){
        m_BTOF_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex] = fitFunction;
        m_BTOF_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_BTOF_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
        if(doSplineFitting) m_BTOF_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = fitter->getCovMatrix();
        else  if(legendreStatus == 0) m_BTOF_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = legendreCovMat;
      }else{
        m_BTOF_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = fitFunction;
        m_BTOF_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_BTOF_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+1),m_cuts->getRapMtM0BinStructure(a_partIndex,1)->GetBinLowEdge(yIndex+2)));
        if(doSplineFitting) m_BTOF_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = fitter->getCovMatrix();
        else  if(legendreStatus == 0)  m_BTOF_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = legendreCovMat;
      }
    }*/
    if(a_doTPCEff && lineFitFunction){
      if(allCentMode){
        m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineFitFunction;
        m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_HighPt_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s High p_{T} Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                             rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineCovMat;
      }else{
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineFitFunction;
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_HighPt_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s High p_{T} Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                            rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineCovMat;
      }
    }




    // create a point at high pT to constrain the high end because we know it's pretty flat
    double highestAvgVal = -1.0;
    if(whereToAddExtraEffPoints > 0.0){
      int    highest_point_index = HistogramUtilities::findLargestXValPoint(effGraph);
      double highest_point_mTm0  = effGraph->GetX()[highest_point_index];
      int    low_range_point     = HistogramUtilities::getPointClosestToXVal(effGraph, highest_point_mTm0 - 0.5);
      double integral_range      = (effGraph->GetX()[highest_point_index] + effGraph->GetEXhigh()[highest_point_index]) - (effGraph->GetX()[low_range_point] - effGraph->GetEXlow()[low_range_point]);
      double integral_temp       = 0.0;
      double sum_of_weights      = 0.0;
      for(int point_index = low_range_point; point_index <= highest_point_index; point_index++){
        double weight   = 1.0/pow( 0.5*(effGraph->GetEYlow()[point_index]+effGraph->GetEYhigh()[point_index]),2);
        sum_of_weights += weight;
        integral_temp  += weight * effGraph->GetY()[point_index]; // * (effGraph->GetEXlow()[point_index] + effGraph->GetEXhigh()[point_index]);
      }
      integral_temp = integral_temp/(sum_of_weights); //*integral_range);
      highestAvgVal = integral_temp;
      int point_index_new = effGraph->GetN();
      effGraph->SetPoint(point_index_new, whereToAddExtraEffPoints, integral_temp);
      double half_integral_range = 0.5*integral_range;
      double error_of_the_weighted_mean = sqrt(1.0/sum_of_weights);
      effGraph->SetPointError(point_index_new, half_integral_range, half_integral_range, error_of_the_weighted_mean, error_of_the_weighted_mean);
      cout << "ExtraPointTest: " << highest_point_index << "  " << highest_point_mTm0 << "  " << low_range_point << "  " << integral_range << "  " << integral_temp << "    " << sum_of_weights << endl;

    }






    //#################     FIT THE GRAPH WITH SYSTEMATIC ##############################################################################################
    //bool mTm0BasicFunctMode = false;
    bool doSimpleTurnOn     = false;
    bool doDoubleTurnOn     = false; // this won't do anything if piecewise is true
    bool doPieceWise         = false; //(a_partIndex == 1);
    bool doPieceWiseTwoLines = false;
    bool doBumpFunct         = false;
    bool doGausBumpFunct     = false;
    bool doSamHepFunct       = false;
    bool doSamHepFunctQuad   = false;
    bool doDoubleTurnOnBump  = false;
    bool doTSpline3          = false;
    bool doTSpline5          = false;
    bool fixLineParam_BumpFuncts = false;
    bool fixBumpMeanAtZero   = false;
    bool fixSamHepConstTerm  = true;
    bool fixTSplineHorLocInSpecCentFromAllCent = true;

    bool localTSplineRefitOverride = true;

    //##############################  FIXED TARGET ###################################
    if(m_fixed_target_mode){

      if(a_partIndex > 3){
        doSamHepFunct = true;
      }

      if(a_doTPCEff){
        if(a_charge > 0){
          if(a_partIndex == 0){
            doSimpleTurnOn = true;
            // fixLineParam_BumpFuncts= true;
            // fixBumpMeanAtZero = true;
          }
          if(a_partIndex == 1){
            doSimpleTurnOn = true;
          }
          if(a_partIndex == 2){
            //doSamHepFunct = true;
          doSimpleTurnOn = true;
          }
        }else{ //negative
          if(a_partIndex == 0){
              doSimpleTurnOn = true;
          }
          if(a_partIndex == 1){
            doSimpleTurnOn = true;
          }
          if(a_partIndex == 2){
            doSimpleTurnOn = true;
          }
        }
      }else{ //BTOF Efficiency
        if(a_partIndex == 0){
          doSimpleTurnOn = true;
        }
        if(a_partIndex == 1){
          doSimpleTurnOn = true;
        }
        if(a_partIndex == 2){
            doSimpleTurnOn = true;
          //doGausBumpFunct = true;
          //fixBumpMeanAtZero = true;
          //fixLineParam_BumpFuncts = true;
        }
      }


    //##############################  COLLIDER MODE ###################################

    }else{

      /*
      if(a_doTPCEff){
        if(a_charge > 0){
          if(a_partIndex == 0){
            doSamHepFunctQuad  = true;
            fixSamHepConstTerm = false;
          }
          if(a_partIndex == 1){
            //doSamHepFunct = true;
            doSamHepFunctQuad  = true;
            fixSamHepConstTerm = false;
          }
          if(a_partIndex == 2){
            //doSamHepFunct = true;
            //doGausBumpFunct = true;
            //fixLineParam_BumpFuncts = true;
            //fixBumpMeanAtZero = true;
            doSamHepFunctQuad  = true;
            fixSamHepConstTerm = false;
          }
        }else{ //negative
          if(a_partIndex == 0){
            doSamHepFunctQuad  = true;
            fixSamHepConstTerm = false;
          }
          if(a_partIndex == 1){
            //doSamHepFunct = true;
            doSamHepFunctQuad  = true;
            fixSamHepConstTerm = false;
          }
          if(a_partIndex == 2){
            //doSimpleTurnOn = true;
            doSamHepFunctQuad  = true;
            fixSamHepConstTerm = false;
          }
        }
      }else{ //BTOF Efficiency
        if(a_partIndex == 0){
          //doSamHepFunct = true;
          doSamHepFunctQuad  = true;
          fixSamHepConstTerm = false;
        }
        if(a_partIndex == 1){
          //doSamHepFunct = true;
          doSamHepFunctQuad  = true;
          fixSamHepConstTerm = false;
        }
        if(a_partIndex == 2){
          //doSamHepFunct = true;
          //fixSamHepConstTerm = false;
          doSamHepFunctQuad  = true;
          fixSamHepConstTerm = false;
          //doGausBumpFunct = true;
          //fixBumpMeanAtZero = true;
          //fixLineParam_BumpFuncts = true;
        }
      }*/

      //doTSpline3 = true;
      doTSpline5 = true;
    }


    double effThresholdToFit = 0.01; // the points below 1% will not be fit
    double fracOfMaxToFit = 0.05;
    TF1* sysFunct = NULL;
    if(doSimpleTurnOn){
      sysFunct = new TF1("sysFunct",
      // "[0]*exp(-[1]*pow((pow(x+[3],2.0)-[3]*[3]),-[2]*0.5))",
      "[0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2]))",
      0.0,highTotalMtM0Range); // function for A*exp(-b*pT^-c). x = mT-m0 if you switch back to above
      sysFunct->SetParNames("y_inf","pT_half","exponent","m");
    }
    if(doDoubleTurnOn){
      sysFunct = new TF1("sysFunct",
      //"[0]*exp(-[1]*pow((pow(x+[3],2.0)-[3]*[3]),-[2]*0.5)) + [4]*exp(-[5]*pow((pow(x+[3],2.0)-[3]*[3]),-[6]*0.5))",
      "[0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2])) + [4]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[5],-[6]))",
      0.0,highTotalMtM0Range); // function for A*exp(-b*pT^-c). x = mT-m0 if you switch back to above
      sysFunct->SetParNames("y_inf_1","pT_half_1","exponent_1","m","y_inf_2","pT_half_2","exponent_2");
    }
    if(doSamHepFunct){
      //sysFunct = new TF1("sysFunct","[0]*exp(-[1]*pow((pow(x+[3],2.0)-[3]*[3]),-[2]*0.5)) + [4]*pow((pow(x+[3],2.0)-[3]*[3]),-[2]*0.5) + [5]*pow((pow(x+[3],2.0)-[3]*[3]),-[2])",0.0,highTotalMtM0Range);
      if(a_partIndex == 1){
        sysFunct = new TF1("sysFunct",
          "([0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2])))*PhysMath::fractionRemainingAfterDistance([4],x,[3], [5], [6], [7]) ",
          0.0,highTotalMtM0Range);
        sysFunct->SetParNames("y_inf","pT_half","exponent","m","a","b","c","rapidity","tau","R","L");
      }else{
        sysFunct = new TF1("sysFunct",
          //"[0]*exp(-[1]*pow((pow(x+[3],2.0)-[3]*[3]),-[2]*0.5)) + [4]*x + [5]*x*x",
          "[0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2])) + [4]*x + [5]*x*x + [6]", 0.0,highTotalMtM0Range);
        sysFunct->SetParNames("y_inf","pT_half","exponent","m","a","b","c");
      }

    }
    if(doSamHepFunctQuad){
      sysFunct = new TF1("sysFunct",
        "[0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2])) + [4]*x + [5]*x*x + [6] + [7]*x*x*x + [8]*x*x*x*x", 0.0,highTotalMtM0Range);
      sysFunct->SetParNames("y_inf","pT_half","exponent","m","a","b","c","d","e");
    }

    if(doBumpFunct){
      sysFunct = new TF1("sysFunct",
        //"[0]*exp(-[1]*pow(sqrt(pow(x+[3],2.0)-[3]*[3]),-[2])) + [4]*PhysMath::skewNormal([5],[6],[7],sqrt(pow(x+[3],2.0)-[3]*[3])) + [8]*sqrt(pow(x+[3],2.0)-[3]*[3])",
        "[0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2])) + [4]*PhysMath::skewNormal([5],[6],[7],sqrt(pow(x+[3],2.0)-[3]*[3])) + [8]*sqrt(pow(x+[3],2.0)-[3]*[3])",
        0.0,highTotalMtM0Range);
      sysFunct->SetParNames("y_inf","pT_half","exponent","m","N_skew","#mu","#sigma","#gamma","a");
    }
    if(doGausBumpFunct){
      sysFunct = new TF1("sysFunct",
        //"[0]*exp(-[1]*pow(sqrt(pow(x+[3],2.0)-[3]*[3]),-[2])) + [4]*PhysMath::gaus([5],[6],sqrt(pow(x+[3],2.0)-[3]*[3])) + [7]*sqrt(pow(x+[3],2.0)-[3]*[3])",
        "[0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2])) + [4]*PhysMath::gaus([5],[6],sqrt(pow(x+[3],2.0)-[3]*[3])) + [7]*sqrt(pow(x+[3],2.0)-[3]*[3])",
        0.0,highTotalMtM0Range);
      sysFunct->SetParNames("y_inf","pT_half","exponent","m","N_gaus","#mu","#sigma","a");
    }
    if(doDoubleTurnOnBump){
      sysFunct = new TF1("sysFunct",
        //"[0]*exp(-[1]*pow(sqrt(pow(x+[3],2.0)-[3]*[3]),-[2])) + [4]*PhysMath::skewNormal([5],[6],[7],sqrt(pow(x+[3],2.0)-[3]*[3])) + [8]*exp(-[9]*pow(x,-[10])) + [11]*x + [12]*x*x",
        "[0]*exp(-log(2.0)*pow(sqrt(pow(x+[3],2.0)-[3]*[3])/[1],-[2])) + [4]*PhysMath::skewNormal([5],[6],[7],sqrt(pow(x+[3],2.0)-[3]*[3])) + [8]*exp(-log(2.0)*pow(x/[9],-[10])) + [11]*x + [12]*x*x",
        0.0,highTotalMtM0Range);
      sysFunct->SetParNames("y_inf_1","pT_half_1","exponent_1","m","N_skew","#mu","#sigma","#gamma","y_inf_2","pT_half_2","exponent_2");
      sysFunct->SetParName(11,"a");
      sysFunct->SetParName(12,"b");
    }
    if(doPieceWise){
      sysFunct = PhysMath::piecewiseTurnOnFunct("sysFunct", 1.0, 1.0, 1.0, 10, 0.0, 1.0);
    }
    if(doPieceWiseTwoLines){//proton
      sysFunct = PhysMath::piecewiseTurnOnTwoLinesFunct("sysFunct",1.0,1.0,1.0,0.3,0.5,0.85,0.0,1.0);
    }
    if(doTSpline3){
      //vector<double> knots;
      m_spline_numKnots = 5;
      if(tspline3_initialization_array_index >= 0 ){
        m_spline_numKnots = ((int) tspline3_initialization_array[tspline3_initialization_array_index].size()) - 5;
        if(tspline3_initialization_array[tspline3_initialization_array_index][tspline3_initialization_array[tspline3_initialization_array_index].size()-1] == -999){
          m_spline_numKnots = m_spline_numKnots - 1;
          localTSplineRefitOverride = false;
        }
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "TSpline Initialization Array:  ";
          for(unsigned int iii = 0; iii < tspline3_initialization_array[tspline3_initialization_array_index].size(); iii++){
            cout << tspline3_initialization_array[tspline3_initialization_array_index][iii] << "  ";
          }
          cout << endl;
        #endif
      }

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        else{
          cout << "WARNING: TSpline3 initialization doesnt exist!" << endl;
        }
      #endif
      /*
      double lowestKnot  = HistogramUtilities::findSmallestXVal(effGraph);
      double highestKnot = HistogramUtilities::findLargestXVal(effGraph);
      double width = highestKnot - lowestKnot;
      knots.push_back(lowestKnot + 0.1*width);
      double step = width/((double) (numKnots-1));
      for(int iii = 1; iii < numKnots - 1; iii++){
        knots.push_back(lowestKnot + ((double) iii)*step);
      }
      knots.push_back(highestKnot - 0.1*width);
      makeTSpline3("Sys_spline", numKnots, knots,1.0,0.0);
      */
      sysFunct = new TF1("sysFunct", this, &EfficiencyFitter::spline_operator_3, 0.0, highTotalMtM0Range, (2*m_spline_numKnots) + 2,"EfficiencyFitter","spline_operator_3");
      for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
        sysFunct->SetParName(knotIndex, Form("x_{%d}", knotIndex));
        //if(knotIndex == 0) sysFunct->SetParName(knotIndex, Form("x_{%d}", knotIndex));
        //else               sysFunct->SetParName(knotIndex, Form("#Deltax_{%d}", knotIndex));
        sysFunct->SetParName(m_spline_numKnots + knotIndex, Form("y_{%d}", knotIndex));
      }
      sysFunct->SetParName(2*m_spline_numKnots, Form("y'_{%d}", 0));
      sysFunct->SetParName(2*m_spline_numKnots+1, Form("y'_{%d}", m_spline_numKnots - 1));

    }
    if(doTSpline5){
      //vector<double> knots;
      m_spline_numKnots = 5;
      if(tspline5_initialization_array_index >= 0 ){
        m_spline_numKnots = ((int) tspline5_initialization_array[tspline5_initialization_array_index].size()) - 5;
        if(tspline5_initialization_array[tspline5_initialization_array_index][tspline5_initialization_array[tspline5_initialization_array_index].size()-1] == -999){
          m_spline_numKnots = m_spline_numKnots - 1;
          localTSplineRefitOverride = false;
        }
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "TSpline Initialization Array:  ";
          for(unsigned int iii = 0; iii < tspline5_initialization_array[tspline5_initialization_array_index].size(); iii++){
            cout << tspline5_initialization_array[tspline5_initialization_array_index][iii] << "  ";
          }
          cout << endl;
        #endif
      }
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Init found " << m_spline_numKnots << " knots to use" << endl;
      #endif

      /*
      double lowestKnot  = HistogramUtilities::findSmallestXVal(effGraph);
      double highestKnot = HistogramUtilities::findLargestXVal(effGraph);
      double width = highestKnot - lowestKnot;
      knots.push_back(lowestKnot + 0.1*width);
      double step = width/((double) (numKnots-1));
      for(int iii = 1; iii < numKnots - 1; iii++){
        knots.push_back(lowestKnot + ((double) iii)*step);
      }
      knots.push_back(highestKnot - 0.1*width);
      makeTSpline3("Sys_spline", numKnots, knots,1.0,0.0);
      */
      if(m_spline_numKnots == 2){
        sysFunct = new TF1("sysFunct", this, &EfficiencyFitter::spline_operator_5_two_point, 0.0, highTotalMtM0Range, (2*m_spline_numKnots) + 4,"EfficiencyFitter","spline_operator_5_two_point");
        sysFunct->SetParName(0, Form("x_{%d}", 0));
        sysFunct->SetParName(1, Form("#Deltax_{%d}", 1));
        sysFunct->SetParName(m_spline_numKnots, Form("y_{%d}", 0));
        sysFunct->SetParName(m_spline_numKnots+1, Form("y_{%d}", 1));
      }else{
        sysFunct = new TF1("sysFunct", this, &EfficiencyFitter::spline_operator_5, 0.0, highTotalMtM0Range, (2*m_spline_numKnots) + 4,"EfficiencyFitter","spline_operator_5");
        for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
          sysFunct->SetParName(knotIndex, Form("x_{%d}", knotIndex));
          //if(knotIndex == 0) sysFunct->SetParName(knotIndex, Form("x_{%d}", knotIndex));
          //else               sysFunct->SetParName(knotIndex, Form("#Deltax_{%d}", knotIndex));
          sysFunct->SetParName(m_spline_numKnots + knotIndex, Form("y_{%d}", knotIndex));
        }
      }

      sysFunct->SetParName(2*m_spline_numKnots, Form("y'_{%d}", 0));
      sysFunct->SetParName(2*m_spline_numKnots+1, Form("y'_{%d}", m_spline_numKnots - 1));
      sysFunct->SetParName(2*m_spline_numKnots+2, Form("y''_{%d}", 0));
      sysFunct->SetParName(2*m_spline_numKnots+3, Form("y''_{%d}", m_spline_numKnots - 1));
    }


    TF1* fit_fromAllCent = NULL;
    if(a_centIndex >= 0 && doBTOFEff && a_doDataDrivenBTOF && m_BTOF_fromData_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]){
        fit_fromAllCent = m_BTOF_fromData_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex];
    }
    if(a_centIndex >= 0 && doBTOFEff && !a_doDataDrivenBTOF && m_BTOF_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]){
        fit_fromAllCent = m_BTOF_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex];
    }
    if(a_centIndex >= 0 && a_doTPCEff && m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]){
        fit_fromAllCent = m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex];
    }
    if(fit_fromAllCent){
      fit_fromAllCent->SetLineColor(kGreen);
      fit_fromAllCent->SetLineStyle(10);
      fit_fromAllCent->SetLineWidth(1);
    }


    /*
    double largestYVal = HistogramUtilities::findLargestYVal(effGraph,0,highTotalMtM0Range);
    double y_inf = effGraph->GetY()[HistogramUtilities::getPointClosestToVal(effGraph, highTotalMtM0Range)];
    // Do the seeding procedure 3 times and average
    double y1_hope = y_inf/8.0;
    double y2_hope = 7.0*y1_hope;
    int lowPoint  = HistogramUtilities::getPointClosestToYVal(effGraph, y1_hope,0,highTotalMtM0Range);
    int highPoint = HistogramUtilities::getPointClosestToYVal(effGraph, y2_hope,0,highTotalMtM0Range);
    if(lowPoint == highPoint) highPoint = lowPoint+1;
    double x1 = PhysMath::pT_from_mTm0(effGraph->GetX()[lowPoint],mass);
    double x2 = PhysMath::pT_from_mTm0(effGraph->GetX()[highPoint],mass);
    double y1 = effGraph->GetY()[lowPoint];
    double y2 = effGraph->GetY()[highPoint];
    double z1 = log(y_inf/y1);
    double z2 = log(y_inf/y2);
    double expon = -log(z1/z2)/log(x1/x2);
    double scale = z1*pow(x1,expon);
    //scale = pow(z1,-1.0/expon)*x1;
    */

    //for adjusted formulas
    /*double y_inf = HistogramUtilities::findLargestYVal(effGraph,0,highTotalMtM0Range);
    int lowPoint = HistogramUtilities::getPointClosestToYVal(effGraph, 0.5*y_inf,0,highTotalMtM0Range);
    double x1 = PhysMath::pT_from_mTm0(effGraph->GetX()[lowPoint],mass);
    double y1 = effGraph->GetY()[lowPoint];
    int highPoint = HistogramUtilities::getPointClosestToXVal(effGraph, 0.5,0,highTotalMtM0Range);
    double x2 = PhysMath::pT_from_mTm0(effGraph->GetX()[highPoint],mass);
    double y2 = effGraph->GetY()[highPoint];
    double scale = x1;
    double expon = log(log(y_inf/y2)*log(2.0))/log(scale/x2);
    //recorrect y_inf because it is overshot on purpose for seeding...
    y_inf = effGraph->GetY()[HistogramUtilities::getPointClosestToXVal(effGraph, 2.0,0,highTotalMtM0Range)];*/

    int maxPoint = HistogramUtilities::findLargestYValPoint(effGraph,0,highTotalMtM0Range);
    double y_inf = effGraph->GetY()[maxPoint] + 0.05;
    if(highestAvgVal < 0.0) highestAvgVal = effGraph->GetY()[maxPoint]; // if avg not computed, use the max y point
    if(y_inf > 1.0)  y_inf = 1.0;
    if(y_inf <= 0.0) y_inf = 0.5;
    int lowPoint = HistogramUtilities::getPointClosestToYVal(effGraph, 0.5*y_inf,0, effGraph->GetX()[maxPoint]);
    double x1 = PhysMath::pT_from_mTm0(effGraph->GetX()[lowPoint],mass);
    double y1 = effGraph->GetY()[lowPoint];
    double x2 = PhysMath::pT_from_mTm0(effGraph->GetX()[maxPoint],mass);
    double y2 = effGraph->GetY()[maxPoint];
    double scale = x1;
    double expon = log(log(y_inf/y2)*log(2.0))/log(scale/x2);


    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "Seeding Systematic Funciton: " << endl;
      cout << "  Low Point: (" << effGraph->GetX()[lowPoint] << ", " << y1 << ")" << endl;
      cout << "  High Point: (" << effGraph->GetX()[maxPoint] << ", " << y2 << ")" << endl;
      cout << "  y_inf:       " << y_inf << endl;
      cout << "  pT_half:  " << scale << "   mTm0_half: " << PhysMath::mTm0(mass,scale) << endl;
      cout << "  Exponent: " << expon << endl;
      cout << "  Check: " << y2 << "  =?  " << y_inf*exp(-log(2.0)*pow(x2/scale,-expon)) << endl;
    #endif






    double massVar   = mass;

    //for peicewise
    double xCutVar   = 0.5;
    double x1CutVar  = 0.6;
    double y1CutVar  = effGraph->Eval(x1CutVar);
    if(y1CutVar <=0 || y1CutVar >1.0) y1CutVar = 0.8;
    double ySlopeVar = 0.0;


    double initialParameters[15];
    initialParameters[0] = y_inf;
    initialParameters[1] = scale;
    initialParameters[2] = expon;
    initialParameters[3] = mass;

    if(doSamHepFunct){
      initialParameters[4] = 0.0; // linear term
      initialParameters[5] = 0.0; // quadratic term
      initialParameters[6] = 0.0; // const term
    }
    if(doSamHepFunctQuad){
      initialParameters[4] = 0.0; // linear term
      initialParameters[5] = 0.0; // quadratic term
      initialParameters[6] = 0.0; // const term
      initialParameters[7] = 0.0; // x^3 term
      initialParameters[8] = 0.0; // x^4 term
    }
    if(doBumpFunct){
      int highest_mTm0_Point = HistogramUtilities::findLargestYValPoint(effGraph, 0,1.0);
      double pT_of_highest   = PhysMath::pT_from_mTm0(effGraph->GetX()[highest_mTm0_Point],mass);
      double eff_of_hightest = effGraph->GetY()[highest_mTm0_Point];
      double turnOnFuctVal   = y_inf*exp(-log(2.0)*pow(pT_of_highest/scale,-expon));
      if(eff_of_hightest > turnOnFuctVal){
        if(fixBumpMeanAtZero) initialParameters[5] = 0.0;
        else                  initialParameters[5] = pT_of_highest;
        initialParameters[6] = 0.4;
        initialParameters[7] = 0.0; // gamma
        initialParameters[8] = 0.0; // linear term
        initialParameters[4] = (eff_of_hightest - turnOnFuctVal)/PhysMath::skewNormal(initialParameters[5],initialParameters[6],initialParameters[7],initialParameters[5]);
      }else{
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " WARNING: the highest point is lower than turnOnVal, seeding wont work for bump funciton... " << endl;
        #endif
        if(fixBumpMeanAtZero) initialParameters[5] = 0.0;
        else                  initialParameters[5] = pT_of_highest;
        initialParameters[6] = 0.4;
        initialParameters[7] = 0.0; // gamma
        initialParameters[8] = 0.0; // linear term
        initialParameters[4] = (0.02)/PhysMath::skewNormal(initialParameters[5],initialParameters[6],initialParameters[7],initialParameters[5]);
      }
    }
    if(doGausBumpFunct){
      int highest_mTm0_Point = HistogramUtilities::findLargestYValPoint(effGraph, 0,1.0);
      double pT_of_highest   = PhysMath::pT_from_mTm0(effGraph->GetX()[highest_mTm0_Point],mass);
      double eff_of_hightest = effGraph->GetY()[highest_mTm0_Point];
      double turnOnFuctVal   = y_inf*exp(-log(2.0)*pow(pT_of_highest/scale,-expon));

      if(eff_of_hightest > turnOnFuctVal){
        if(fixBumpMeanAtZero) initialParameters[5] = 0.0;
        else                  initialParameters[5] = pT_of_highest;
        initialParameters[6] = 0.4;
        initialParameters[7] = 0.0;// linear term
        initialParameters[4] = (eff_of_hightest - turnOnFuctVal)/PhysMath::gaus(initialParameters[5],initialParameters[6],initialParameters[5]);
      }else{
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " WARNING: the highest point is lower than turnOnFuctVal, seeding wont work for bump funciton... " << endl;
        #endif
        if(fixBumpMeanAtZero) initialParameters[5] = 0.0;
        else                  initialParameters[5] = pT_of_highest;
        initialParameters[6] = 0.4;
        initialParameters[7] = 0.0; // linear term
        initialParameters[4] = (0.02)/PhysMath::skewNormal(initialParameters[5],initialParameters[6],initialParameters[7],initialParameters[5]);
      }
    }
    if(doDoubleTurnOnBump){
      int highest_mTm0_Point = HistogramUtilities::findLargestYValPoint(effGraph, 0,1.0);
      double pT_of_highest   = PhysMath::pT_from_mTm0(effGraph->GetX()[highest_mTm0_Point],mass);
      double eff_of_hightest = effGraph->GetY()[highest_mTm0_Point];
      double turnOnFuctVal   = y_inf*exp(-log(2.0)*pow(pT_of_highest/scale,-expon));
      if(eff_of_hightest > turnOnFuctVal){
        if(fixBumpMeanAtZero) initialParameters[5] = 0.0;
        else                  initialParameters[5] = pT_of_highest;
        initialParameters[6] = 0.4;
        initialParameters[7] = 0.0; // gamma
        initialParameters[8]  = 0.001; //second yInf
        initialParameters[9]  = scale;
        initialParameters[10] = expon; // gamma
        initialParameters[11] = 0.0; // linear term
        initialParameters[12] = 0.0; // square term
        initialParameters[4] = (eff_of_hightest - turnOnFuctVal)/PhysMath::skewNormal(initialParameters[5],initialParameters[6],initialParameters[7],initialParameters[5]);
      }else{
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " WARNING: the highest point is lower than turnOnFuctVal, seeding wont work for bump funciton... " << endl;
        #endif
        if(fixBumpMeanAtZero) initialParameters[5] = 0.0;
        else                  initialParameters[5] = pT_of_highest;
        initialParameters[6] = 0.4;
        initialParameters[7] = 0.0; // gamma
        initialParameters[8]  = 0.001; // second yInf
        initialParameters[9]  = scale;
        initialParameters[10] = expon; // gamma
        initialParameters[11] = 0.0; // linear term
        initialParameters[4] = (0.02)/PhysMath::skewNormal(initialParameters[5],initialParameters[6],initialParameters[7],initialParameters[5]);
      }
    }
    if(doTSpline3 || doTSpline5){
      //effThresholdToFit = 0.8;
      fracOfMaxToFit = 0.05;
      vector<double> knots;
      if(!fixTSplineHorLocInSpecCentFromAllCent || (fixTSplineHorLocInSpecCentFromAllCent && !fit_fromAllCent )){
        if(doTSpline3){
          if(tspline3_initialization_array_index < 0){
            // if not seeded, just spread evenly (this doesn't do well)
            double lowestKnot  = HistogramUtilities::findSmallestXVal(effGraph);
            double highestKnot = HistogramUtilities::findLargestXVal(effGraph);
            double width = highestKnot - lowestKnot;
            knots.push_back(lowestKnot + 0.1*width);
            double step = width/((double) (m_spline_numKnots-1));
            for(int iii = 1; iii < m_spline_numKnots - 1; iii++){
              knots.push_back(lowestKnot + ((double) iii)*step);
            }
            knots.push_back(highestKnot - 0.1*width);
          }else{
            for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
              knots.push_back(tspline3_initialization_array[tspline3_initialization_array_index][5 + knotIndex]);
            }
          }
        }

        if(doTSpline5){
          if(tspline5_initialization_array_index < 0){
            // if not seeded, just spread evenly (this doesn't do well)
            double lowestKnot  = HistogramUtilities::findSmallestXVal(effGraph);
            double highestKnot = HistogramUtilities::findLargestXVal(effGraph);
            double width = highestKnot - lowestKnot;
            knots.push_back(lowestKnot + 0.1*width);
            double step = width/((double) (m_spline_numKnots-1));
            for(int iii = 1; iii < m_spline_numKnots - 1; iii++){
              knots.push_back(lowestKnot + ((double) iii)*step);
            }
            knots.push_back(highestKnot - 0.1*width);
          }else{
            for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
              knots.push_back(tspline5_initialization_array[tspline5_initialization_array_index][5 + knotIndex]);
            }
          }
        }



      }else{
        for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
          knots.push_back(fit_fromAllCent->GetParameter(knotIndex));
        }
      }



      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Knot Locations: ";
      #endif
      for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
        if(doTSpline5 && m_spline_numKnots == 2 && knotIndex == 1){
          initialParameters[knotIndex] = knots[knotIndex] - knots[knotIndex-1];        // horizontal location difference
        }else{
          initialParameters[knotIndex] = knots[knotIndex];        // horizontal location
        }
        //if(knotIndex == 0) initialParameters[knotIndex] = knots[knotIndex];        // horizontal location
        //else              initialParameters[knotIndex] = knots[knotIndex] - knots[knotIndex-1];        // horizontal location difference
        initialParameters[m_spline_numKnots + knotIndex] = effGraph->Eval(knots[knotIndex]); // vertical location
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "   (" << knots[knotIndex] << ", " << initialParameters[m_spline_numKnots + knotIndex] << ")";
        #endif
      }
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << endl;
      #endif

      if(fit_fromAllCent){
        initialParameters[2*m_spline_numKnots]   = fit_fromAllCent->GetParameter(2*m_spline_numKnots); // derivative at low end
        initialParameters[2*m_spline_numKnots+1] = fit_fromAllCent->GetParameter(2*m_spline_numKnots+1); // derivative at high end
        if(doTSpline5){
          initialParameters[2*m_spline_numKnots+2]   = fit_fromAllCent->GetParameter(2*m_spline_numKnots+2); // 2nd derivative at low end
          initialParameters[2*m_spline_numKnots+3]   = fit_fromAllCent->GetParameter(2*m_spline_numKnots+3); // 2nd derivative at high end
        }
      }else{
        initialParameters[2*m_spline_numKnots]   = 1.0; // derivative at low end
        if(m_spline_numKnots > 1){
          initialParameters[2*m_spline_numKnots] = (initialParameters[m_spline_numKnots+1] - initialParameters[m_spline_numKnots])/(initialParameters[1] - initialParameters[0]);
        }
        if(std::isfinite(initialParameters[2*m_spline_numKnots])) initialParameters[2*m_spline_numKnots] = 1.0;
        initialParameters[2*m_spline_numKnots+1] = 0.0; // derivative at high end
        if(doTSpline5){
          initialParameters[2*m_spline_numKnots+2]   = 0.0; // 2nd derivative at low end
          initialParameters[2*m_spline_numKnots+3]   = 0.0; // 2nd derivative at high end
        }
      }

      /*
      for(int knotIndex = 0; knotIndex < m_spline->GetNp(); knotIndex++){
        double knot_x, knot_y;
        m_spline->GetKnot(knotIndex,knot_x,knot_y);
        knot_y = effGraph->Eval(knot_x);
        initialParameters[knotIndex] = knot_y;
        m_spline->SetPoint(knotIndex,knot_x,knot_y);
      }
      initialParameters[m_spline->GetNp()]   = 1.0;
      initialParameters[m_spline->GetNp()+1] = 0.0;
      */
    }


    //#### Populate if you have a seed matrix
    bool didSeedingSys = false;
    vector<vector<double> > seedMatrix;
    if(a_doTPCEff) seedMatrix = m_tpcEff_Systematic_SeedTuple;
    else           seedMatrix = m_btofEff_Systematic_SeedTuple;

    double toleranceForSeedRap = 0.02; // can't be exact because 0 is 1E-17 or something in the histogram file
    for(int seedIndex = 0; seedIndex < seedMatrix.size(); seedIndex++){
      int    particle_index_of_seed = (int) (seedMatrix[seedIndex][0]+0.1);
      double rapidity_of_seed       = seedMatrix[seedIndex][1];

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Using Seed Matrix for Initialization..." << endl;
      #endif

      if(a_partIndex == particle_index_of_seed && fabs(rapidity_of_seed - rapidity) < toleranceForSeedRap){
        didSeedingSys = true;
        initialParameters[0] = seedMatrix[seedIndex][2];
        initialParameters[1] = seedMatrix[seedIndex][3];
        initialParameters[2] = seedMatrix[seedIndex][4];
        if(doPieceWise && seedMatrix[seedIndex].size() > 7){
          initialParameters[3]   = seedMatrix[seedIndex][5];
          initialParameters[4]   = seedMatrix[seedIndex][6];
          initialParameters[5]   = seedMatrix[seedIndex][7];
        }
        if(doPieceWiseTwoLines && seedMatrix[seedIndex].size() > 9){
          initialParameters[3]   = seedMatrix[seedIndex][5];
          initialParameters[4]   = seedMatrix[seedIndex][6];
          initialParameters[5]   = seedMatrix[seedIndex][7];
          initialParameters[6]   = seedMatrix[seedIndex][8];
          initialParameters[7]   = seedMatrix[seedIndex][9];
        }
        if(doSamHepFunct  && seedMatrix[seedIndex].size() > 7){
          initialParameters[3]   = seedMatrix[seedIndex][5];
          initialParameters[4]   = seedMatrix[seedIndex][6];
          initialParameters[5]   = seedMatrix[seedIndex][7];
        }
        if(doSamHepFunctQuad  && seedMatrix[seedIndex].size() > 9){
          cout << "initializing samHepFunctQuad by seed matrix" << endl;
          initialParameters[3]   = seedMatrix[seedIndex][5];
          initialParameters[4]   = seedMatrix[seedIndex][6];
          initialParameters[5]   = seedMatrix[seedIndex][7];
          initialParameters[6]   = seedMatrix[seedIndex][8];
          initialParameters[7]   = seedMatrix[seedIndex][9];
          initialParameters[8]   = seedMatrix[seedIndex][10];
        }
        if(doBumpFunct && seedMatrix[seedIndex].size() > 10){
          initialParameters[3]   = seedMatrix[seedIndex][5];
          initialParameters[4]   = seedMatrix[seedIndex][6];
          initialParameters[5]   = seedMatrix[seedIndex][7];
          initialParameters[6]   = seedMatrix[seedIndex][8];
          initialParameters[7]   = seedMatrix[seedIndex][9];
          initialParameters[8]   = seedMatrix[seedIndex][10];
        }
        if(doGausBumpFunct && seedMatrix[seedIndex].size() > 9){
          initialParameters[3]   = seedMatrix[seedIndex][5];
          initialParameters[4]   = seedMatrix[seedIndex][6];
          initialParameters[5]   = seedMatrix[seedIndex][7];
          initialParameters[6]   = seedMatrix[seedIndex][8];
          initialParameters[7]   = seedMatrix[seedIndex][9];
        }
        if(doDoubleTurnOnBump && seedMatrix[seedIndex].size() > 11){
          initialParameters[3]   = seedMatrix[seedIndex][5];
          initialParameters[4]   = seedMatrix[seedIndex][6];
          initialParameters[5]   = seedMatrix[seedIndex][7];
          initialParameters[6]   = seedMatrix[seedIndex][8];
          initialParameters[7]   = seedMatrix[seedIndex][9];
          initialParameters[8]   = seedMatrix[seedIndex][10];
          initialParameters[9]   = seedMatrix[seedIndex][11];
          initialParameters[10]  = seedMatrix[seedIndex][12];
          initialParameters[11]  = seedMatrix[seedIndex][13];
          //initialParameters[12]   = seedMatrix[seedIndex][14];
        }
        if(doTSpline3 && seedMatrix[seedIndex].size() - 2 == 2*m_spline_numKnots + 2){
          for(unsigned int iii = 0; iii < 2*m_spline_numKnots + 2; iii++){
            initialParameters[iii]   = seedMatrix[seedIndex][iii + 2];
          }
        }
        if(doTSpline5 && seedMatrix[seedIndex].size() - 2 == 2*m_spline_numKnots + 4){
          for(unsigned int iii = 0; iii < 2*m_spline_numKnots + 4; iii++){
            initialParameters[iii]   = seedMatrix[seedIndex][iii + 2];
          }
        }
        break;
      }
    }


    double lowMtM0FitRange_Array[3][2];
    double highMtM0FitRange_Array[3][2];



    if(m_fixed_target_mode){
      if(a_doTPCEff){
        lowMtM0FitRange_Array[0][0] = 0.0;
        lowMtM0FitRange_Array[0][1] = 0.0;
        lowMtM0FitRange_Array[1][0] = 0.0;
        lowMtM0FitRange_Array[1][1] = 0.0;
        lowMtM0FitRange_Array[2][0] = 0.0;
        lowMtM0FitRange_Array[2][1] = 0.0;
        highMtM0FitRange_Array[0][0] = 8.0;
        highMtM0FitRange_Array[0][1] = 8.0;
        highMtM0FitRange_Array[1][0] = 8.0;
        highMtM0FitRange_Array[1][1] = 8.0;
        highMtM0FitRange_Array[2][0] = 8.0;
        highMtM0FitRange_Array[2][1] = 8.0;
      }else{
        lowMtM0FitRange_Array[0][0] = 0.0;
        lowMtM0FitRange_Array[0][1] = 0.0;
        lowMtM0FitRange_Array[1][0] = 0.0;
        lowMtM0FitRange_Array[1][1] = 0.0;
        lowMtM0FitRange_Array[2][0] = 0.0;
        lowMtM0FitRange_Array[2][1] = 0.0;
        highMtM0FitRange_Array[0][0] = 8.0;
        highMtM0FitRange_Array[0][1] = 8.0;
        highMtM0FitRange_Array[1][0] = 8.0;
        highMtM0FitRange_Array[1][1] = 8.0;
        highMtM0FitRange_Array[2][0] = 8.0;
        highMtM0FitRange_Array[2][1] = 8.0;
      }
    }else{
      if(a_doTPCEff){
        lowMtM0FitRange_Array[0][0] = 0.0;
        lowMtM0FitRange_Array[0][1] = 0.0;
        lowMtM0FitRange_Array[1][0] = 0.0;
        lowMtM0FitRange_Array[1][1] = 0.0;
        lowMtM0FitRange_Array[2][0] = 0.08;
        lowMtM0FitRange_Array[2][1] = 0.08;
        highMtM0FitRange_Array[0][0] = 8.0;
        highMtM0FitRange_Array[0][1] = 8.0;
        highMtM0FitRange_Array[1][0] = 8.0;
        highMtM0FitRange_Array[1][1] = 8.0;
        highMtM0FitRange_Array[2][0] = 8.0;
        highMtM0FitRange_Array[2][1] = 8.0;
      }else{
        lowMtM0FitRange_Array[0][0] = 0.0;
        lowMtM0FitRange_Array[0][1] = 0.0;
        lowMtM0FitRange_Array[1][0] = 0.0;
        lowMtM0FitRange_Array[1][1] = 0.0;
        lowMtM0FitRange_Array[2][0] = 0.0;
        lowMtM0FitRange_Array[2][1] = 0.0;
        highMtM0FitRange_Array[0][0] = 8.0;
        highMtM0FitRange_Array[0][1] = 8.0;
        highMtM0FitRange_Array[1][0] = 8.0;
        highMtM0FitRange_Array[1][1] = 8.0;
        highMtM0FitRange_Array[2][0] = 8.0;
        highMtM0FitRange_Array[2][1] = 8.0;
      }
    }




    double lowFitRange  = lowMtM0FitRange_Array[a_partIndex][pmIndex];
    double highFitRange = highMtM0FitRange_Array[a_partIndex][pmIndex];
    double lowByEffThresh = 0.0;
    for(int pointIndex = 0; effThresholdToFit > 0.0 && pointIndex < effGraph->GetN(); pointIndex++){
      double eff  = effGraph->GetY()[pointIndex];
      double mTm0 = effGraph->GetX()[pointIndex];
      if(eff > effThresholdToFit && mTm0 > lowFitRange){
        lowByEffThresh = mTm0 - effGraph->GetEXlow()[pointIndex];
        break;
      }
    }

    double lowByFracThresh = 0.0;
    for(int pointIndex = 0; fracOfMaxToFit > 0.0 && pointIndex < effGraph->GetN(); pointIndex++){
      double eff  = effGraph->GetY()[pointIndex];
      double mTm0 = effGraph->GetX()[pointIndex];
      if(eff/highestAvgVal > fracOfMaxToFit && mTm0 > lowFitRange){
        lowByFracThresh = mTm0 - effGraph->GetEXlow()[pointIndex];
        break;
      }
    }

    if(lowByEffThresh > lowByFracThresh){
      lowFitRange = lowByEffThresh;
      cout << "Low Range By EffThresh" << endl;
    }else{
      lowFitRange = lowByFracThresh;
      cout << "Low Range By FracThresh" << endl;
    }




    bool seededFromAllCent = false;
    if(useAllCentFitsAsSeeds && fit_fromAllCent){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Using All Cent Fit as a Seed" << endl;
      #endif
      for(int iii = 0; iii < fit_fromAllCent->GetNpar(); iii++){
        initialParameters[iii] = fit_fromAllCent->GetParameter(iii);
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "     Par " << iii << ": " << initialParameters[iii] << endl;
        #endif
      }
      seededFromAllCent = true;
    }

    cout << "seededFromAllCent " << seededFromAllCent << endl;
    cout << "didSeedingSys " << didSeedingSys << endl;
    cout << "scale " << scale << endl;
    cout << "expon " << expon << endl;

    cout << "thing   " << !seededFromAllCent << endl;
    cout << "thing   " <<  !didSeedingSys  << endl;
    cout << "thing   " <<  !std::isfinite(scale)  << endl;
    cout << "thing   " <<   !std::isfinite(expon) << endl;
    cout << "thing   " <<    (expon <= 0.0) << endl;
    cout << "thing   " <<     (scale <= 0.0) << endl;
    cout << "sysFunct   " <<     sysFunct << endl;


    if(!seededFromAllCent && !didSeedingSys && (!std::isfinite(scale) || !std::isfinite(expon) || expon <= 0.0 || scale <= 0.0)){
      cout << "WARNING: Initializion didn't work..." << endl;
      sysFunct->SetParameter(0,initialParameters[0]);
      sysFunct->SetParameter(1,0.15);
      sysFunct->SetParameter(2,1.0);
      sysFunct->SetParameter(3,initialParameters[0]);
      sysFunct->SetParLimits(0,-1.0,2.0);
      sysFunct->SetParLimits(1,-50.0,50.0);
      sysFunct->SetParLimits(2,-5.0,15);
      sysFunct->SetParLimits(3,0.0,100.0);
      if(doDoubleTurnOn){
        sysFunct->SetParameter(0,initialParameters[0]*0.51);
        sysFunct->SetParameter(4,initialParameters[0]*0.49);
        sysFunct->SetParameter(5,0.15);
        sysFunct->SetParameter(6,1.0);
        sysFunct->SetParLimits(4,0,2.0);
        sysFunct->SetParLimits(5,0.0,50.0);
        sysFunct->SetParLimits(6,0.0,15);
      }
    }else{
      sysFunct->SetParameter(0,initialParameters[0]);
      sysFunct->SetParameter(1,initialParameters[1]);
      sysFunct->SetParameter(2,initialParameters[2]);
      sysFunct->SetParameter(3,initialParameters[3]);
      sysFunct->SetParLimits(0,-5.0,20.0); // sometimes, the region you fit is all < 1, but shape likes >1 so I'll allow it
      sysFunct->SetParLimits(1,0.0,200);
      sysFunct->SetParLimits(2,0.0,20);
      sysFunct->SetParLimits(3,0.0,100.0);
      if(doDoubleTurnOn){
        sysFunct->SetParameter(0,initialParameters[0]*0.51);
        sysFunct->SetParameter(4,initialParameters[0]*0.49);
        sysFunct->SetParameter(5,initialParameters[1]);
        sysFunct->SetParameter(6,initialParameters[2]);
        sysFunct->SetParLimits(4,0,2.0); // sometimes, the region you fit is all < 1, but shape likes >1 so I'll allow it
        sysFunct->SetParLimits(5,0.0,200);
        sysFunct->SetParLimits(6,0.0,20);
      }
    }
    if(doPieceWise){
        sysFunct->SetParameter(3,initialParameters[3]);
        sysFunct->SetParameter(4,initialParameters[4]);
        sysFunct->SetParameter(5,initialParameters[5]);
        sysFunct->SetParLimits(1,0.0,200.0);
        sysFunct->SetParLimits(3,0.0,5.0); // x cut param
        sysFunct->SetParLimits(4,-2.0,10.0); // slope param
        //sysFunct->SetParLimits(5,-10.0,1000.0); // allow mass to very? Why not. It helps with the low pT curvature.
        sysFunct->FixParameter(5,mass); // m is tightly correlated with b (par 1)
    }
    if(doPieceWiseTwoLines){
      sysFunct->SetParameter(3,initialParameters[3]);
      sysFunct->SetParameter(4,initialParameters[4]);
      sysFunct->SetParameter(5,initialParameters[5]);
      sysFunct->SetParameter(6,initialParameters[6]);
      sysFunct->SetParameter(7,initialParameters[7]);
      sysFunct->SetParLimits(1,0.0,200.0);
      sysFunct->SetParLimits(3,0.0,5.0); // x0 cut param
      sysFunct->SetParLimits(4,0.0,5.0); // x1 cut param
      sysFunct->SetParLimits(5,0.0,1.0); // y1 cut param
      sysFunct->SetParLimits(6,-2.0,10.0); // slope param
      sysFunct->SetParLimits(7,0.0,200);

      //sysFunct->SetParLimits(5,-10.0,1000.0); // allow mass to very? Why not. It helps with the low pT curvature.
      //sysFunct->FixParameter(5,mass); // m is tightly correlated with b (par 1)
      //if(mTm0BasicFunctMode) sysFunct->FixParameter(5,0.0);
      //sysFunct->FixParameter(3,10.0);
      //sysFunct->FixParameter(4,0.0);
      //sysFunct->SetParLimits(5,-10,1000);
      if(fabs(rapidity) > 0.65){
        sysFunct->FixParameter(3,10.0);
        sysFunct->FixParameter(4,12.0);
        sysFunct->FixParameter(5,1.0);
        sysFunct->FixParameter(6,0.0);
      }
    }
    if(doSamHepFunct){
      sysFunct->FixParameter(3,initialParameters[3]);
      if(a_partIndex != 1){  // not a Kaon
        sysFunct->SetParameter(4,initialParameters[4]);
        sysFunct->SetParameter(5,initialParameters[5]);
        //sysFunct->SetParLimits(4,0.0,5.0); // linear param
        //sysFunct->SetParLimits(5,0.0,1.0); // quadratic param
        if(fixSamHepConstTerm){
          sysFunct->FixParameter(6,0.0);
        }else{
          sysFunct->SetParameter(0,0.01); //rescale to make it a const instead
          sysFunct->SetParameter(6,y_inf);
          sysFunct->SetParLimits(6,-1.0,2.0);
        }
      }else{
        sysFunct->FixParameter(4,0.0);
        sysFunct->FixParameter(5,0.0);
        sysFunct->FixParameter(6,0.0);
        sysFunct->FixParameter(4,rapidity);
        sysFunct->FixParameter(5,1.238E-8); //mean lifetime kaon +/-
        sysFunct->SetParameter(6,1.0);     // Radius of TPC
        sysFunct->SetParLimits(6,0.0,3.0); // Radius of TPC
        sysFunct->FixParameter(7,2.0); // Half length of TPC
      }
    }
    if(doSamHepFunctQuad){
      sysFunct->FixParameter(3,initialParameters[3]);
      sysFunct->SetParameter(4,initialParameters[4]);
      sysFunct->SetParameter(5,initialParameters[5]);
      sysFunct->SetParameter(7,initialParameters[7]);
      sysFunct->SetParameter(8,initialParameters[8]);
      //sysFunct->SetParLimits(4,0.0,5.0); // linear param
      //sysFunct->SetParLimits(5,0.0,1.0); // quadratic param
      if(fixSamHepConstTerm){
        sysFunct->FixParameter(6,0.0);
      }else{
        sysFunct->SetParameter(6,initialParameters[6]);

        // sysFunct->SetParameter(0,0.01); //rescale to make it a const instead
        // sysFunct->SetParameter(6,y_inf);
        // sysFunct->SetParLimits(6,-1.0,2.0);
      }

    }
    if(doBumpFunct){
      sysFunct->FixParameter(3,initialParameters[3]);
      sysFunct->SetParameter(4,initialParameters[4]); // N
      sysFunct->SetParameter(6,initialParameters[6]); // Sig
      sysFunct->SetParameter(7,initialParameters[7]); // gamma
      sysFunct->SetParLimits(4,0.0,15.0);
      sysFunct->SetParLimits(6,0.00001,10.0);
      sysFunct->SetParLimits(7,-1.0,1.0); //hard limits?
      if(fixLineParam_BumpFuncts){
        sysFunct->FixParameter(8,0.0); // linear term
      }else{
        sysFunct->SetParameter(8,initialParameters[8]); // linear term
        sysFunct->SetParLimits(8,-10.0,20.0);
      }
      if(fixBumpMeanAtZero){
        sysFunct->FixParameter(5,0.0); // Mu
      }else{
        sysFunct->SetParameter(5,initialParameters[5]); // Mu
        sysFunct->SetParLimits(5,-15.0,15.0);
      }

    }
    if(doGausBumpFunct){
      sysFunct->FixParameter(3,initialParameters[3]);
      sysFunct->SetParameter(4,initialParameters[4]); // N
      sysFunct->SetParameter(6,initialParameters[6]); // Sig
      sysFunct->SetParLimits(4,-15.0,15.0);
      sysFunct->SetParLimits(6,1E-16,20.0);
      if(fixLineParam_BumpFuncts){
        sysFunct->FixParameter(7,0.0); // linear term
      }else{
        sysFunct->SetParameter(7,initialParameters[7]); // linear term
        sysFunct->SetParLimits(7,-10.0,20.0);
      }
      if(fixBumpMeanAtZero){
        sysFunct->FixParameter(5,0.0); // Mu
      }else{
        sysFunct->SetParameter(5,initialParameters[5]); // Mu
        sysFunct->SetParLimits(5,-15.0,15.0);
      }
    }
    if(doDoubleTurnOnBump){
      sysFunct->FixParameter(3,initialParameters[3]);
      sysFunct->SetParameter(4,initialParameters[4]); // N
      sysFunct->SetParameter(6,initialParameters[6]); // Sig
      sysFunct->SetParameter(7,initialParameters[7]); // gamma
      sysFunct->SetParameter(8,initialParameters[8]); // y_inf 2
      sysFunct->SetParameter(9,initialParameters[9]); // scale 2
      sysFunct->SetParameter(10,initialParameters[10]); // expon 2
      sysFunct->SetParLimits(4,0.0,15.0);
      sysFunct->SetParLimits(6,0.00001,10.0);
      sysFunct->SetParLimits(7,-1.0,1.0); //hard limits?

      sysFunct->SetParLimits(8,-1.0,2.0); // sometimes, the region you fit is all < 1, but shape likes >1 so I'll allow it
      sysFunct->SetParLimits(9,0.0,200);
      sysFunct->SetParLimits(10,0.0,20);
      sysFunct->FixParameter(12,0.0); // square term
      if(fixLineParam_BumpFuncts){
        sysFunct->FixParameter(11,0.0); // linear term
        sysFunct->FixParameter(12,0.0); // square term
      }else{
        sysFunct->SetParameter(11,initialParameters[11]); // linear term
        sysFunct->SetParLimits(11,-10.0,20.0);
        //sysFunct->SetParameter(12,initialParameters[12]); // linear term
        //sysFunct->SetParLimits(12,-10.0,20.0);
      }
      if(fixBumpMeanAtZero){
        sysFunct->FixParameter(5,0.0); // Mu
      }else{
        sysFunct->SetParameter(5,initialParameters[5]); // Mu
        sysFunct->SetParLimits(5,-15.0,15.0);
      }
    } // if didnt seed

    // Scale the function up to match a point
    if(doPieceWise){
      sysFunct->SetParameter(0,1.0);
      sysFunct->SetParameter(0,effGraph->Eval(sysFunct->GetParameter(3))/sysFunct->Eval(sysFunct->GetParameter(3)));
    }
    if(doPieceWiseTwoLines){
      sysFunct->SetParameter(5,effGraph->Eval(sysFunct->GetParameter(4)));
    }

    TGraph* tspline_initial_points = NULL;
    if(doTSpline3 || doTSpline5){
      tspline_initial_points = new TGraph();
      tspline_initial_points->SetMarkerStyle(20);
      tspline_initial_points->SetMarkerColor(kYellow+1);

      double width = highFitRange - lowFitRange;


      for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
        if((a_centIndex >= 0 && fixTSplineHorLocInSpecCentFromAllCent && fit_fromAllCent) || !localTSplineRefitOverride){
          sysFunct->FixParameter(knotIndex,initialParameters[knotIndex]);
        }else{
          sysFunct->SetParameter(knotIndex,initialParameters[knotIndex]);
          if(doTSpline5 && m_spline_numKnots == 2 && knotIndex == 1){
            double temp_max = 0.7*width;
            if(temp_max < 0.7) temp_max = 5.0;
            sysFunct->SetParLimits(knotIndex,0.5,temp_max);
          }else{
            //if(knotIndex == 0) sysFunct->SetParLimits(knotIndex,-0.1,5.5);
            //else               sysFunct->SetParLimits(knotIndex,0.0001,5.5);
            //sysFunct->SetParLimits(knotIndex,lowFitRange + 0.03*width,highFitRange - 0.03*width);
            sysFunct->SetParLimits(knotIndex,-1,highFitRange);

          }
        }

        if(knotIndex == 0){
          sysFunct->SetParameter(m_spline_numKnots + knotIndex,initialParameters[m_spline_numKnots + knotIndex]);
          sysFunct->SetParLimits(m_spline_numKnots + knotIndex,-500,500);
        }else{
          sysFunct->SetParameter(m_spline_numKnots + knotIndex,initialParameters[m_spline_numKnots + knotIndex]);
          sysFunct->SetParLimits(m_spline_numKnots + knotIndex,0.0,1.1);
        }


        if(doTSpline5 && m_spline_numKnots == 2 && knotIndex == 1){
          tspline_initial_points->SetPoint(knotIndex, initialParameters[knotIndex] + tspline_initial_points->GetX()[knotIndex-1], initialParameters[m_spline_numKnots + knotIndex]);
        }else{
          tspline_initial_points->SetPoint(knotIndex, initialParameters[knotIndex],initialParameters[m_spline_numKnots + knotIndex]);
        }
        //if(knotIndex == 0) tspline_initial_points->SetPoint(knotIndex, initialParameters[knotIndex],initialParameters[m_spline_numKnots + knotIndex]);
        //else               tspline_initial_points->SetPoint(knotIndex, initialParameters[knotIndex] + tspline_initial_points->GetX()[knotIndex-1], initialParameters[m_spline_numKnots + knotIndex]);

      }
      sysFunct->SetParameter(2*m_spline_numKnots,initialParameters[2*m_spline_numKnots]);
      sysFunct->SetParLimits(2*m_spline_numKnots,-50,100.0);
      //sysFunct->FixParameter(m_spline_numKnots+1,0.0);
      sysFunct->SetParameter(2*m_spline_numKnots+1,initialParameters[2*m_spline_numKnots+1]);
      sysFunct->SetParLimits(2*m_spline_numKnots+1,-10.0,10.0);

      if(doTSpline5){
        sysFunct->SetParameter(2*m_spline_numKnots+2,initialParameters[2*m_spline_numKnots+2]);
        sysFunct->SetParLimits(2*m_spline_numKnots+2,-300,300.0);
        sysFunct->SetParameter(2*m_spline_numKnots+3,initialParameters[2*m_spline_numKnots+3]);
        sysFunct->SetParLimits(2*m_spline_numKnots+3,-300.0,300.0);
      }
    }

    sysFunct->SetNpx(300);
    sysFunct->SetLineColor(kBlue);

    TF1* initSysFunct = (TF1*) sysFunct->Clone();
    initSysFunct->SetName("initSysFunctEff");
    initSysFunct->SetLineColor(kYellow+1);
    initSysFunct->SetLineWidth(1);



    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << " Starting simple fit range at mT-m0 = " << lowFitRange << endl;
      cout << "\n\n  Systematic funciton initialization: " << endl;
      HistogramUtilities::printTF1Info(sysFunct);
    #endif



    TF1*          fitFunctSys = NULL;
    TMatrixDSym*  sysCovMat   = NULL;
    TGraphErrors* conf_graph  = NULL;
    GraphFitter*  gr_fitter   = NULL;

    double chiPerNDF_Sys = 0.0;
    if(m_useIntegralChiSqr){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Using Integral Method for Chi Squared " << endl;
      #endif
      gr_fitter = new GraphFitter();
      gr_fitter->loadGraphAsym(effGraph);
      gr_fitter->loadTF1(sysFunct, !(doTSpline3 || doTSpline5));
      gr_fitter->setFitRange(lowFitRange,highFitRange);

      if((doTSpline5 || doTSpline3) && initiallyFitTSplineWithFixedHorLocs){
        vector<int> paramsToFix;
        for(int iii = 0; iii < m_spline_numKnots; iii++){
          paramsToFix.push_back(iii);
        }
        gr_fitter->setParametersToBeFixed(paramsToFix);
      }else{
        vector<int> paramsToFix; // empty
        gr_fitter->setParametersToBeFixed(paramsToFix);
      }

      gr_fitter->minimize();
      int status = gr_fitter->getStatus();
      fitFunctSys = gr_fitter->getResultFunction();
      cout << "First Status: " << status << endl;

      if((doTSpline5 || doTSpline3) && refitTSplineWithFixedHorLocs){
        if(fitFunctSys){
          fitFunctSys->SetName("sys_funct_again");
          gr_fitter->loadTF1(fitFunctSys, false);
        }
        //else{
        //  gr_fitter->loadTF1(sysFunct,false);
        //}

        vector<int> paramsToFix;
        for(int iii = 0; iii < m_spline_numKnots; iii++){
          paramsToFix.push_back(iii);
        }
        gr_fitter->setParametersToBeFixed(paramsToFix);
        gr_fitter->minimize();
        status = gr_fitter->getStatus();
        cout << "Second Status: " << status << endl;
      }


      if(fitFunctSys){//status == 0){
        //fitFunctSys = gr_fitter->getResultFunction();
        sysCovMat = (TMatrixDSym*) gr_fitter->getCovMat()->Clone();
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "####  Covariance Matrix: " << endl;
          gr_fitter->getCovMat()->Print();
        #endif
        conf_graph = HistogramUtilities::makeErrorBandGraph(fitFunctSys, sysCovMat, 300);

      }else{
        cout << "WARNING: BAD FIT.... " << status << "  " << m_partInfo->GetParticleName(a_partIndex,a_charge) << " Cent: " << a_centIndex << " y:" << rapidity  << " FAILED! " << endl;
      }

    }else{

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Not Using Integral Method for Chi Squared " << endl;
      #endif
      string fitOptions = "EX0 S"; //E
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        fitOptions = fitOptions + " ";
      #else
        fitOptions = fitOptions + " Q";
      #endif
      TFitResult*  result_ptr;
      if(doInitialFitWithHistogram){
        fitOptions += "I";
        TFitResultPtr resultPtr = effHisto->Fit(sysFunct,fitOptions.c_str(),"",lowFitRange,highFitRange);
        result_ptr = resultPtr.Get();
        fitFunctSys = effHisto->GetFunction("sysFunct");
      }else{
        TFitResultPtr resultPtr = effGraph->Fit(sysFunct,fitOptions.c_str(),"",lowFitRange,highFitRange);
        result_ptr = resultPtr.Get();
        fitFunctSys = effGraph->GetFunction("sysFunct");
      }
      if(!result_ptr){
        cout << "  ERROR: Efficiency fitting for " << m_partInfo->GetParticleName(a_partIndex,a_charge) << " Cent: " << a_centIndex << " y:" << rapidity  << " FAILED! " << endl;
        cout << "           PTR: "<< result_ptr <<" CONTINUING..." << endl;
        continue;
      }
      int status = result_ptr->Status();
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "---- Fit Status: " << status << endl;
      #endif


      if((doTSpline5 || doTSpline3) && refitTSplineWithFixedHorLocs && fitFunctSys && localTSplineRefitOverride){
        TF1* fitFunctSysTwo = NULL;
        if(doTSpline3) fitFunctSysTwo = new TF1("sysFunctTwo", this, &EfficiencyFitter::spline_operator_3, 0.0, highTotalMtM0Range, (2*m_spline_numKnots) + 2,"EfficiencyFitter","spline_operator_3");
        if(doTSpline5){
          if(m_spline_numKnots == 2){
            fitFunctSysTwo = new TF1("sysFunctTwo", this, &EfficiencyFitter::spline_operator_5_two_point, 0.0, highTotalMtM0Range, (2*m_spline_numKnots) + 4,"EfficiencyFitter","spline_operator_5_two_point");
          }else{
            fitFunctSysTwo = new TF1("sysFunctTwo", this, &EfficiencyFitter::spline_operator_5, 0.0, highTotalMtM0Range, (2*m_spline_numKnots) + 4,"EfficiencyFitter","spline_operator_5");
          }

        }
        fitFunctSysTwo->SetLineColor(fitFunctSys->GetLineColor());
        for(int iii = 0; iii < 2*m_spline_numKnots+2; iii++){
          fitFunctSysTwo->SetParName(iii,fitFunctSys->GetParName(iii));
        }
        if(doTSpline5){
          fitFunctSysTwo->SetParName(2*m_spline_numKnots+2,fitFunctSys->GetParName(2*m_spline_numKnots+2));
          fitFunctSysTwo->SetParName(2*m_spline_numKnots+3,fitFunctSys->GetParName(2*m_spline_numKnots+3));
        }


        for(int iii = 0; iii < m_spline_numKnots; iii++){
          fitFunctSysTwo->FixParameter(iii,fitFunctSys->GetParameter(iii));
        }
        TFitResultPtr resultPtrSecond = effGraph->Fit(fitFunctSysTwo,fitOptions.c_str(),"",lowFitRange,highFitRange);
        fitFunctSys = effGraph->GetFunction("sysFunctTwo");
        result_ptr  = resultPtrSecond.Get();
        status      = result_ptr->Status();
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "---- Second Fit Status: " << status << endl;
        #endif
      }

      //fitFunctSys = effGraph->GetFunction("sysFunct");
      if(fitFunctSys) fitFunctSys->SetRange(lowFitRange,highFitRange);
      sysCovMat     = (TMatrixDSym*) result_ptr->GetCovarianceMatrix().Clone();
      chiPerNDF_Sys = result_ptr->Chi2() / result_ptr->Ndf();
      if(fitFunctSys) conf_graph = HistogramUtilities::makeErrorBandGraph(fitFunctSys, sysCovMat, 300);

    }


    if(fitFunctSys && PhysMath::isTMatrixDSymEmpty(*sysCovMat)){
      cout << "WARNING: Covariance Matrix is Empty!!!!!!!!!!!!!!!" << endl;
    }



    //##############    BTOF EFF   #########################
    if(doBTOFEff && a_doDataDrivenBTOF){
       if(allCentMode){
        m_BTOF_fromData_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();
        m_BTOF_fromData_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_BTOF_fromData_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent Integrated y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_BTOF_fromData_Eff_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex] = sysCovMat;
        m_BTOF_fromData_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiPerNDF_Sys);

      }else{
        m_BTOF_fromData_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();;
        m_BTOF_fromData_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_BTOF_fromData_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent %d y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_BTOF_fromData_Eff_CovMat_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = sysCovMat;
        m_BTOF_fromData_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiPerNDF_Sys);

      }
    }
    if(doBTOFEff && !a_doDataDrivenBTOF){  // embedding BTOF Efficiency
      if(allCentMode){
        m_BTOF_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();;
        m_BTOF_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_BTOF_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent Integrated y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_BTOF_Eff_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex] = sysCovMat;
        m_BTOF_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiPerNDF_Sys);

      }else{
        m_BTOF_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();;
        m_BTOF_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_BTOF_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent %d y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_BTOF_Eff_CovMat_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = sysCovMat;
        m_BTOF_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiPerNDF_Sys);

      }
    }

    //##############    TPC EFF   #########################
    if(a_doTPCEff && fitFunctSys){
      if(allCentMode){
        m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();;
        m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_AllCent_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_TPC_Eff_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent Integrated y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex] = sysCovMat;
        m_TPC_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiPerNDF_Sys);
      }else{
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();;
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_TPC_Eff_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent %d y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_TPC_Eff_CovMat_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = sysCovMat;
        m_TPC_Efficiency_ChiSqrd_Sys[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiPerNDF_Sys);
      }
    }

    vector<double> finalPars;
    finalPars.push_back(a_partIndex);
    finalPars.push_back(rapidity);
    for(int parIndex = 0;fitFunctSys && parIndex < fitFunctSys->GetNpar(); parIndex++){
      finalPars.push_back(fitFunctSys->GetParameter(parIndex));
    }
    paramStorageForFutureSeeds.push_back(finalPars);

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      if(fitFunctSys){
        cout << " Systematic Covariance Matrix: " << endl;
        sysCovMat->Print();
      }
    #endif


    TGraphAsymmErrors* pullGraph_Sys = NULL;
    if(fitFunctSys){
      pullGraph_Sys = HistogramUtilities::makeFitPullGraph(effGraph, fitFunctSys, m_useIntegralChiSqr);
      pullGraph_Sys->SetMarkerStyle(20);
      pullGraph_Sys->SetMarkerSize(0.5);
      pullGraph_Sys->SetMarkerColor(kBlack);
      pullGraph_Sys->SetLineColor(kBlack);
    }


    TGraph* tspline_final_points = NULL;
    if( (doTSpline3 || doTSpline5) && fitFunctSys){
      tspline_final_points = new TGraph();
      tspline_final_points->SetMarkerStyle(20);
      tspline_final_points->SetMarkerColor(kBlue);
      for(int knotIndex = 0; knotIndex < m_spline_numKnots; knotIndex++){
        if(doTSpline5 && m_spline_numKnots == 2 && knotIndex == 1){
          tspline_final_points->SetPoint(knotIndex, fitFunctSys->GetParameter(knotIndex) + tspline_final_points->GetX()[knotIndex-1], fitFunctSys->GetParameter(m_spline_numKnots + knotIndex));
        }else{
          tspline_final_points->SetPoint(knotIndex, fitFunctSys->GetParameter(knotIndex),fitFunctSys->GetParameter(m_spline_numKnots + knotIndex));
        }
        //if(knotIndex == 0) tspline_final_points->SetPoint(knotIndex, fitFunctSys->GetParameter(knotIndex),fitFunctSys->GetParameter(m_spline_numKnots + knotIndex));
        //else               tspline_final_points->SetPoint(knotIndex, fitFunctSys->GetParameter(knotIndex) + tspline_final_points->GetX()[knotIndex-1], fitFunctSys->GetParameter(m_spline_numKnots + knotIndex));
      }
    }



    fittingCanvas->cd();
    topPad_base->cd();
    double lowYRange = 0;
    double highYRange = 0;
    HistogramUtilities::verticalDrawRange(effGraph, lowYRange,highYRange);
    TH1F* frame = fittingCanvas->DrawFrame(0.0,lowYRange,highTotalMtM0Range,highYRange);
    frame->SetTitle(Form("%s Efficiency %s Cent %d y=[%1.2f,%1.2f]; m_{T}-m_{0} (GeV/c^{2}); %s Efficiency",
          modeName.c_str(),
          m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
          rapidity_bin_axis->GetBinLowEdge(yIndex+1),
          rapidity_bin_axis->GetBinLowEdge(yIndex+2),modeName.c_str()));
    gPad->SetRightMargin(0.4);
    effGraph->Draw("PSame");
    if(initialLegendreSeed) initialLegendreSeed->Draw("same");
    initSysFunct->Draw("same");
    if(tspline_initial_points) tspline_initial_points->Draw("P");
    if(fit_fromAllCent) fit_fromAllCent->Draw("same");
    if(conf_graph) conf_graph->Draw("3");
    if(fitFunctSys) fitFunctSys->Draw("same");
    if(splinePoints) splinePoints->Draw("P");
    if(tspline_final_points) tspline_final_points->Draw("P");
    if(lineFitFunction) lineFitFunction->Draw("Same");
    gPad->Update();

    bottomPad_base->cd();
    TH1F* frame_pull_sys = fittingCanvas->DrawFrame(0.0,-4.0,highTotalMtM0Range,4.0);
    frame_pull_sys->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Standardized Pulls"));
    gPad->SetRightMargin(0.4);
    if(pullGraph_Sys) pullGraph_Sys->Draw("PSame");
    TLine* oneLineSys = new TLine(0.0,0.0,highTotalMtM0Range,0.0);
    oneLineSys->SetLineStyle(7);
    oneLineSys->SetLineColor(kBlack);
    oneLineSys->Draw("same");
    gSystem->ProcessEvents();
    if(doBTOFEff && a_doDataDrivenBTOF){
      if(allCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_AllCent/DataDriven_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_Cent%02d/DataDriven_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
      }
    }
    if(doBTOFEff && !a_doDataDrivenBTOF){
      if(allCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_AllCent/yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_Cent%02d/yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
      }
    }
    if(a_doTPCEff){
      if(allCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_AllCent/yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_Cent%02d/yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
      }
    }
    delete frame;



    if(!a_doLegendreResiduleFits) continue; ////////////////#############################################################      Legendre fit of residuals below this point



   //############ form the deviation of the efficiency and fit with legendre polynomials
   if(fitFunctSys){
    TGraphAsymmErrors* effGraphResidulals = (TGraphAsymmErrors*) effGraph->Clone();
    string residuleGraphName = effGraph->GetName();
    residuleGraphName = residuleGraphName + "_residules";
    effGraphResidulals->SetName(residuleGraphName.c_str());
    for(int pointIndex = 0; pointIndex < effGraphResidulals->GetN(); pointIndex++){
      double X_ofPoint = effGraphResidulals->GetX()[pointIndex];
      double Y_ofPoint = effGraphResidulals->GetY()[pointIndex];
      double mainFunctEval = fitFunctSys->Eval(X_ofPoint);
      double residule = Y_ofPoint - mainFunctEval;
      //effGraphResidulals->SetPointY(pointIndex,residule);
      effGraphResidulals->SetPoint(pointIndex,X_ofPoint,residule);
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Residule for " << effGraphResidulals->GetX()[pointIndex] << " is " <<  effGraphResidulals->GetY()[pointIndex] << " with shift " << mainFunctEval << endl;
      #endif
    }
    //TF1* fitFunctionRes = NULL;
    //TF1* initialLegendreSeedRes = NULL;
    //TMatrixDSym* legendreCovMat;
    //int legendreStatus = -1;
    //#########  FIT WITH LEGENDRE POLYNOMIAL ##########

    double lowRange = 0.0;
    if(removeFirstZerosForLegendreFits){
      for(int pointIndex = 0; pointIndex < effGraph->GetN(); pointIndex++){
        if(effGraph->GetX()[pointIndex] < 1.0 && effGraph->GetY()[pointIndex] < threshold_for_remove_zeros) lowRange = effGraph->GetX()[pointIndex] + 0.5*effGraph->GetEXhigh()[pointIndex];
      }
    }

    double pT_LowCut = 0.5;
    if(a_doTPCEff) pT_LowCut = 0.15;
    double mTm0_LowCut = PhysMath::mTm0(m_partInfo->GetParticleMass(a_partIndex),pT_LowCut);
    if(lowRange < mTm0_LowCut) lowRange = mTm0_LowCut;
    if(lowRange < effGraph->GetX()[0]) lowRange = effGraph->GetX()[0];
    if(m_TPC_Eff_lowMtM0Legendre_Tuple.size() > 0){
      double toleranceForSeedRap = 0.01;
      for(int seedIndex = 0; seedIndex < m_TPC_Eff_lowMtM0Legendre_Tuple.size(); seedIndex++){
        if(a_partIndex == m_TPC_Eff_lowMtM0Legendre_Tuple[seedIndex][0] && fabs(m_TPC_Eff_lowMtM0Legendre_Tuple[seedIndex][1] - rapidity) < toleranceForSeedRap){
          if(lowRange < m_TPC_Eff_lowMtM0Legendre_Tuple[seedIndex][2]) lowRange = m_TPC_Eff_lowMtM0Legendre_Tuple[seedIndex][2];
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "  Using seed matrix : " << m_TPC_Eff_lowMtM0Legendre_Tuple[seedIndex][0] << " " << m_TPC_Eff_lowMtM0Legendre_Tuple[seedIndex][1] << " " << m_TPC_Eff_lowMtM0Legendre_Tuple[seedIndex][2] << endl;
          #endif
          break;
        }
      }
    }
    double mTm0CutOffLine = m_BTOF_highMtM0ToLine[a_partIndex]->Eval(fabs(rapidity));
    double mTm0LegendreCutOffForFit = mTm0CutOffLine + 0.5;//extra lengendre bit that gets trimmed later
    TF1* fitFunctionRes = PhysMath::getLegendreExpansion("legendreFunctBtofEffRes", polynomialOrder,false,lowRange, mTm0LegendreCutOffForFit);
    PhysMathHistogramUtilities::initializeLegendreExpansion(fitFunctionRes, effGraphResidulals, lowRange, mTm0LegendreCutOffForFit, false);
    TF1* initialLegendreSeedRes = (TF1*) fitFunctionRes->Clone();
    initialLegendreSeedRes->SetName("initialLegendreSeedRes");
    initialLegendreSeedRes->SetLineColor(kYellow);

    /*
    m_currentFitFunction = fitFunctionRes;
    m_currentGraphToFit = effGraphResidulals;
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer();
    ROOT::Math::Functor* chiSqrdFunctor = new ROOT::Math::Functor(this,&EfficiencyFitter::TGraphAsymmErrors_ChiSqr,fitFunctionRes->GetNpar());
    cout << "created chi squared sys functor with NDim " << chiSqrdFunctor->NDim() << "  functor: " << &EfficiencyFitter::TGraphAsymmErrors_ChiSqr  << endl;
    minimizer->SetFunction(*chiSqrdFunctor);
    for(int paramIndex = 0; paramIndex < fitFunctionRes->GetNpar(); paramIndex++){
      double lowLimit = 0;
      double highLimit = 0;
      fitFunctionRes->GetParLimits(paramIndex,lowLimit,highLimit);
      minimizer->SetVariable(paramIndex, fitFunctionRes->GetParName(paramIndex), fitFunctionRes->GetParameter(paramIndex) , 0.1*fitFunctionRes->GetParameter(paramIndex));
      minimizer->SetVariableLimits(paramIndex,lowLimit,highLimit);
    }

    minimizer->SetPrintLevel(5);
    minimizer->SetMaxIterations(3000);
    minimizer->SetMaxFunctionCalls(200000);
    minimizer->SetValidError(true);
    minimizer->Minimize();
    int legendreStatus = minimizer->Status();
    const double* finalParams = minimizer->X();
    if(legendreStatus == 0){
      cout << "Fit Status: " << legendreStatus << endl;
    }else{
      cout << "################################ STATUS OF FIT = " << legendreStatus << "   IS BAD!!!!!  ################################" << endl;
      cout << "Current Parameters: " << endl;
      for(int iii = 0; iii < fitFunctionRes->GetNpar(); iii++){
        cout << "  " << finalParams[iii];
      }
      cout << endl;
      continue;
    }
    int largeSize = fitFunctionRes->GetNpar() + 3;
    TMatrixDSym matrixTemp(largeSize);
    for(int iii = 0; iii < largeSize; iii++){
      for(int jjj = 0; jjj < largeSize; jjj++){
        matrixTemp[iii][jjj] = 0.0;
      }
    }
    for(int iii = 0; iii < fitFunctionRes->GetNpar(); iii++){
      for(int jjj = 0; jjj < fitFunctionRes->GetNpar(); jjj++){
        matrixTemp[iii][jjj] = minimizer->CovMatrix(iii,jjj);
      }
    }
    TMatrixDSym* legendreCovMatRes = (TMatrixDSym*) matrixTemp.Clone();
    //legendreCovMatRes = &matrixTemp;
    TF1* functOrig = fitFunctionRes;
    for(int iii = 0; iii < fitFunctionRes->GetNpar(); iii++){
        functOrig->SetParameter(iii,finalParams[iii]);
    }
    */

    int numZerosToAdd = 4;

    TF1* functOrig      = NULL;
    int  legendreStatus = -1;
    TMatrixDSym* legendreCovMatRes = NULL;
    if(m_useIntegralChiSqr){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Using Integral Method for Chi Squared " << endl;
      #endif
      GraphFitter* gr_fitter = new GraphFitter();
      gr_fitter->loadGraphAsym(effGraphResidulals);
      gr_fitter->loadTF1(fitFunctionRes);
      //void setParametersToBeFixed(vector<int> a_paramsToFix){m_paramsToFix = a_paramsToFix;};
      gr_fitter->minimize();
      legendreStatus = gr_fitter->getStatus();
      if(legendreStatus == 0){
        functOrig = gr_fitter->getResultFunction();
        TMatrixDSym matrixTemp = PhysMath::addZerosToCovMatrix(*gr_fitter->getCovMat(), numZerosToAdd);
        legendreCovMatRes = (TMatrixDSym*) matrixTemp.Clone();
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "####  initial covariance matrix: " << endl;
          gr_fitter->getCovMat()->Print();
          cout << "####  final covariance matrix: " << endl;
          legendreCovMatRes->Print();
        #endif
      }
    }else{
      string fitOptions = "S EX0";
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        fitOptions = fitOptions + " ";
      #else
        fitOptions = fitOptions + " Q";
      #endif
      TFitResultPtr fitPtrRes = effGraphResidulals->Fit(fitFunctionRes,fitOptions.c_str(),"",lowRange,mTm0LegendreCutOffForFit);
      //ndof =  HistogramUtilities::getNPointsInRange(effGraph,lowRange,mTm0LegendreCutOffForFit) - fitFunction->GetNpar() - 2; //-2 becasue window bounds are fit parameters
      //chiSqr = effGraph->Chisquare(fitFunction,"R EX0");
      if(!fitPtrRes.Get()){
        cout << "  ERROR: This is probably a bin you don't care about, but residual efficiency fitting for " << m_partInfo->GetParticleName(a_partIndex,a_charge) << " Cent: " << a_centIndex << " yIndex:" << yIndex  << " FAILED! " << endl;
        cout << "     PTR: "<< fitPtrRes.Get() <<" CONTINUING..." << endl;
        continue;
      }
      legendreStatus = fitPtrRes->Status();
      functOrig = effGraphResidulals->GetFunction("legendreFunctBtofEffRes");
      TMatrixDSym matrixTemp = PhysMath::addZerosToCovMatrix(fitPtrRes->GetCovarianceMatrix(), numZerosToAdd);
      legendreCovMatRes = (TMatrixDSym*) matrixTemp.Clone();
    }

    cout << " status: " << legendreStatus << endl;
    fitFunctionRes = (TF1*) functOrig->Clone();
    fittingCanvas->cd();
    double lowYRange  = 0;
    double highYRange = 0;
    HistogramUtilities::verticalDrawRange(effGraphResidulals, lowYRange,highYRange,0.0,highTotalMtM0Range);
    TH1F* frame2 = fittingCanvas->DrawFrame(0.0,lowYRange,highTotalMtM0Range,highYRange);
    frame2->SetTitle(Form("%s Efficiency Residules %s Cent %d y=[%1.2f,%1.2f]; m_{T}-m_{0} (GeV/c^{2}); %s Efficiency",
          modeName.c_str(),
          m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
          rapidity_bin_axis->GetBinLowEdge(yIndex+1),
          rapidity_bin_axis->GetBinLowEdge(yIndex+2),modeName.c_str())) ;
    effGraphResidulals->Draw("PSame");
    if(initialLegendreSeedRes) initialLegendreSeedRes->Draw("same");
    if(fitFunctionRes) fitFunctionRes->Draw("same");

    gPad->Update();
    gSystem->ProcessEvents();
    if(doBTOFEff && a_doDataDrivenBTOF){
      if(allCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_AllCent/DataDriven_Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_Cent%02d/DataDriven_Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
      }
    }
    if(doBTOFEff && !a_doDataDrivenBTOF){
      if(allCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_AllCent/Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_Cent%02d/Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
      }
    }
    if(a_doTPCEff){
      if(allCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_AllCent/Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_Cent%02d/Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
      }
    }

    delete frame2;

    if(fitFunctionRes){
      string totalFormula = PhysMath::getLegendreExpansionString(polynomialOrder);
      int numParams = fitFunctionRes->GetNpar();
      string additionalFormula = "";
      if(doPieceWise){
        if(doPieceWiseTwoLines){
          additionalFormula = Form(" + PhysMath::piecewiseTurnOnTwoLines([%d],[%d],[%d],[%d],[%d],[%d],[%d],[%d],x)",numParams,numParams+1,numParams+2,numParams+3,numParams+4,numParams+5,numParams+6,numParams+7);
        }else{
          additionalFormula = Form(" + PhysMath::piecewiseTurnOn([%d],[%d],[%d],[%d],[%d],[%d],x)",numParams,numParams+1,numParams+2,numParams+3,numParams+4,numParams+5);
        }
      }else{
        //if(mTm0BasicFunctMode) additionalFormula = Form("+[%d]*exp(-[%d]*pow(x,-[%d]))",numParams,numParams+1,numParams+2);
        //else
        additionalFormula = Form("+[%d]*exp(-[%d]*pow((pow(x+[%d],2.0)-[%d]*[%d]),-[%d]*0.5))",numParams,numParams+1,numParams+3,numParams+3,numParams+3,numParams+2);
      }
      totalFormula = totalFormula + additionalFormula;
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "combined formula: " << totalFormula << endl;
      #endif

      TF1* combinedFunct = new TF1("combinedFunct",totalFormula.c_str(),lowRange,mTm0CutOffLine);
      for(int paramIndex = 0; paramIndex < numParams; paramIndex++){
        combinedFunct->SetParameter(paramIndex,fitFunctionRes->GetParameter(paramIndex));
      }
      if(doPieceWise){
        if(doPieceWiseTwoLines){
          combinedFunct->SetParameter(numParams,  fitFunctSys->GetParameter(0));
          combinedFunct->SetParameter(numParams+1,fitFunctSys->GetParameter(1));
          combinedFunct->SetParameter(numParams+2,fitFunctSys->GetParameter(2));
          combinedFunct->SetParameter(numParams+3,fitFunctSys->GetParameter(3));
          combinedFunct->SetParameter(numParams+4,fitFunctSys->GetParameter(4));
          combinedFunct->SetParameter(numParams+5,fitFunctSys->GetParameter(5));
          combinedFunct->SetParameter(numParams+6,fitFunctSys->GetParameter(6));
          combinedFunct->SetParameter(numParams+7,fitFunctSys->GetParameter(7));
        }else{
          combinedFunct->SetParameter(numParams,  fitFunctSys->GetParameter(0));
          combinedFunct->SetParameter(numParams+1,fitFunctSys->GetParameter(1));
          combinedFunct->SetParameter(numParams+2,fitFunctSys->GetParameter(2));
          combinedFunct->SetParameter(numParams+3,fitFunctSys->GetParameter(3));
          combinedFunct->SetParameter(numParams+4,fitFunctSys->GetParameter(4));
          combinedFunct->SetParameter(numParams+5,fitFunctSys->GetParameter(5));
        }

      }else{
        combinedFunct->SetParameter(numParams,  fitFunctSys->GetParameter(0));
        combinedFunct->SetParameter(numParams+1,fitFunctSys->GetParameter(1));
        combinedFunct->SetParameter(numParams+2,fitFunctSys->GetParameter(2));
        combinedFunct->SetParameter(numParams+3,fitFunctSys->GetParameter(3));
      }


      TGraphAsymmErrors* pullGraph = HistogramUtilities::makeFitPullGraph(effGraph, combinedFunct, false);
      pullGraph->SetMarkerStyle(20);
      pullGraph->SetMarkerSize(0.5);
      pullGraph->SetMarkerColor(kBlack);
      pullGraph->SetLineColor(kBlack);

      TGraphAsymmErrors* pullGraph_line = NULL;
      if(lineFitFunction){
        pullGraph_line = HistogramUtilities::makeFitPullGraph(effGraph, lineFitFunction, false);
        pullGraph_line->SetMarkerStyle(20);
        pullGraph_line->SetMarkerSize(0.5);
        pullGraph_line->SetMarkerColor(kBlack);
        pullGraph_line->SetLineColor(kBlack);
      }





      if(doBTOFEff && a_doDataDrivenBTOF){
        if(allCentMode){
          m_BTOF_fromData_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex] = combinedFunct;
          m_BTOF_fromData_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
          m_BTOF_fromData_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          /*if(lineFitFunction){
            m_BTOF_fromData_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineFitFunction;
            m_BTOF_fromData_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_HighPt_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
            m_BTOF_fromData_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s High p_{T} Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
            m_BTOF_fromData_Eff_CovMat_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineCovMat;
          }*/
          if(doSplineFitting) m_BTOF_fromData_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = sp_fitter->getCovMatrix();
          else                m_BTOF_fromData_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = legendreCovMatRes;

        }else{
          m_BTOF_fromData_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = combinedFunct;
          m_BTOF_fromData_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
          m_BTOF_fromData_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          /*if(lineFitFunction){
            m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineFitFunction;
            m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_DataDriven_%s_HighPt_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
            m_BTOF_fromData_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s High p_{T}Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
            m_BTOF_fromData_Eff_CovMat_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineCovMat;
          }*/
          if(doSplineFitting) m_BTOF_fromData_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = sp_fitter->getCovMatrix();
          else                m_BTOF_fromData_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = legendreCovMatRes;
        }
      }
      if(doBTOFEff && !a_doDataDrivenBTOF){
        if(allCentMode){
          m_BTOF_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex] = combinedFunct;
          m_BTOF_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
          m_BTOF_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          if(doSplineFitting)          m_BTOF_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = sp_fitter->getCovMatrix();
          else if(legendreStatus == 0) m_BTOF_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = legendreCovMatRes;
        }else{
          m_BTOF_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = combinedFunct;
          m_BTOF_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("btofEfficiencyFit_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
          m_BTOF_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("BTOF Efficiency %s Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          if(doSplineFitting) m_BTOF_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = sp_fitter->getCovMatrix();
          else                m_BTOF_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = legendreCovMatRes;
        }
      }
      if(a_doTPCEff){
        if(allCentMode){
          m_TPC_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex] = (TF1*) combinedFunct->Clone();
          m_TPC_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
          m_TPC_Eff_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          if(doSplineFitting) m_TPC_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = sp_fitter->getCovMatrix();
          else                m_TPC_Eff_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = legendreCovMatRes;
          //m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiSqr/ndof);

          /*if(lineFitFunction){
            m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineFitFunction;
            m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_HighPt_AllCent_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
            m_TPC_Eff_Fits_HighPt_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s High p_{T} Cent Integrated y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
            m_TPC_Eff_CovMat_HighPt_AllCent[a_partIndex][yIndex][pmIndex] = lineCovMat;
          }*/

        }else{
          m_TPC_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) combinedFunct->Clone();
          m_TPC_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
          m_TPC_Eff_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          //m_TPC_Efficiency_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiSqr/ndof);
          if(doSplineFitting) m_TPC_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = sp_fitter->getCovMatrix();
          else                m_TPC_Eff_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = legendreCovMatRes;

          /*if(lineFitFunction){
            m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineFitFunction;
            m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("tpcEfficiencyFit_%s_HighPt_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
            m_TPC_Eff_Fits_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("TPC Efficiency %s High p_{T} Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
            m_TPC_Eff_CovMat_HighPt_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = lineCovMat;
          }*/
        }
      }







      //############     DRAW THE FINAL RESULT    ###############################
      TCanvas* fittingCanvas_final = new TCanvas("fittingCanvas_btofFinal","Fitting Canvas BTOFEFF",1300,1000);
      fittingCanvas_final->Draw();
      TPad* topPad    = new TPad("topPad_btofEff","topPad_btofEff",0.0,0.3,1.0,1.0);
      TPad* bottomPad = new TPad("bottomPad_btofEff","bottomPad_btofEff",0.0,0.0,1.0,0.3);
      topPad->Draw();
      gPad->SetRightMargin(0.38);
      bottomPad->Draw();
      gPad->SetRightMargin(0.38);

      topPad->cd();
      //gStyle->SetOptStat(1000000001); // just name of histogram
      //gStyle->SetOptTitle(1); // print name of histogram
      topPad->SetLogy(false);
      double lowYRange = 0;
      double highYRange = 0;
      HistogramUtilities::verticalDrawRange(effGraph, lowYRange,highYRange,0.0,highTotalMtM0Range);
      TH1F* frame3 = fittingCanvas->DrawFrame(0.0,lowYRange,highTotalMtM0Range,highYRange);
      frame3->SetTitle(Form("%s Efficiency  %s Cent %d y=[%1.2f,%1.2f]; m_{T}-m_{0} (GeV/c^{2}); %s Efficiency",modeName.c_str(),
            m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
            rapidity_bin_axis->GetBinLowEdge(yIndex+1),
            rapidity_bin_axis->GetBinLowEdge(yIndex+2),
            modeName.c_str()));
      effGraph->Draw("PSame");
      combinedFunct->Draw("same");
      if(lineFitFunction) lineFitFunction->Draw("same");

      bottomPad->cd();
      //gStyle->SetOptStat(0);
      //gStyle->SetOptTitle(0);
      bottomPad->SetLogy(false);
      TH1F* frame_pull = fittingCanvas->DrawFrame(0.0,-4.0,highTotalMtM0Range,4.0);
      frame_pull->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Standardized Pulls"));
      pullGraph->Draw("PSame");
      if(pullGraph_line) pullGraph_line->Draw("PSame");
      TLine* oneLine = new TLine(0.0,0.0,highTotalMtM0Range,0.0);
      oneLine->SetLineStyle(7);
      oneLine->SetLineColor(kBlack);
      oneLine->Draw("same");

      fittingCanvas_final->Update();
      gSystem->ProcessEvents();

      if(doBTOFEff && a_doDataDrivenBTOF){
        if(allCentMode){
          fittingCanvas_final->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_AllCent/DataDriven_CombinedFinal_yIndex_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
        }else{
          fittingCanvas_final->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_Cent%02d/DataDriven_CombinedFinal_yIndex_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
        }
      }
      if(doBTOFEff && !a_doDataDrivenBTOF){
        if(allCentMode){
          fittingCanvas_final->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_AllCent/CombinedFinal_yIndex_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
        }else{
          fittingCanvas_final->SaveAs(Form("%s/%s/BTOF_Efficiency/BTOF_Efficiency_Cent%02d/CombinedFinal_yIndex_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
        }
      }
      if(a_doTPCEff){
        if(allCentMode){
          fittingCanvas_final->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_AllCent/CombinedFinal_yIndex_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
        }else{
          fittingCanvas_final->SaveAs(Form("%s/%s/TPC_Efficiency/TPC_Efficiency_Cent%02d/CombinedFinal_yIndex_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
        }
      }
      delete frame3;
      delete frame_pull;
      delete oneLine;
      delete fittingCanvas_final;
      delete pullGraph;
    } // if fit funct residual
   }//if fit funct sys



    if(sp_fitter) m_splinesToDelete.push_back(sp_fitter);
    if(splinePoints) delete splinePoints;
    //delete fitter;
    //delete frame;
    if(effHisto) delete effHisto;

  } // end rapidity loop



  // this prints out so it can be copied directly into a macro
  cout << "Parameters For Future Seeds:" << endl;
  for(int iii = 0;iii < paramStorageForFutureSeeds.size(); iii++){
    cout << "{";
    int sm1 = paramStorageForFutureSeeds[iii].size() - 1;
    for(int parIndex = 0; parIndex < sm1; parIndex++){
      cout << paramStorageForFutureSeeds[iii][parIndex] << ",";
    }
    if(sm1 >= 0){
      cout << paramStorageForFutureSeeds[iii][sm1];
    }
    cout << "}," << endl;
  }




}




//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

void EfficiencyFitter::fitFeedDown(int a_partIndex,int a_charge, int a_centIndex){

  bool doSplineFitting = false;

  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Call to EfficiencyFitter::fitFeedDown(int " << a_partIndex << ",int " << a_charge << ", int " << a_centIndex << ") " << endl;
  #endif


  bool allCentMode = false;
  if(a_centIndex < 0) allCentMode = true;
  gSystem->mkdir(Form("%s/%s",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  if(allCentMode){
    gSystem->mkdir(Form("%s/%s/FeedDown_AllCent",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  }else{
    gSystem->mkdir(Form("%s/%s/FeedDown_ByCent",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    gSystem->mkdir(Form("%s/%s/FeedDown_ByCent/Cent%02d",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex));
  }

  TCanvas* fittingCanvas = new TCanvas("fittingCanvas_feeddown","Fitting Canvas");
  fittingCanvas->SetWindowSize(1300,700);
  fittingCanvas->SetCanvasSize(1200,600);
  fittingCanvas->cd();


  unsigned int pmIndex = 0;
  if(a_charge == -1) pmIndex = 1;

  TAxis* rapidity_bin_axis = m_FeedDown_ChiSqrd[a_partIndex][pmIndex]->GetXaxis();
  for(int rapBin = 1; rapBin <= g_EffFit_MaxNumRapBins; rapBin++){
    if(m_testingMidRapidity && rapBin != rapidity_bin_axis->FindBin(m_ycm)) continue;

    int yIndex = rapBin - 1;
    //TGraphAsymmErrors* graph_primary = m_FeedDown_primaryOverTotal_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
    TGraphAsymmErrors* graph_feeddown = NULL;
    if(allCentMode) graph_feeddown = m_FeedDown_feeddownOverTotal_AllCent[a_partIndex][yIndex][pmIndex];
    else graph_feeddown = m_FeedDown_feeddownOverTotal_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];

    if(!graph_feeddown || graph_feeddown->GetN() < 4){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Skipping " << m_partInfo->GetParticleName(a_partIndex,a_charge).Data() << "  Cent: " << a_centIndex << " yIndex:" << yIndex << "  feeddown ptr: " << graph_feeddown << endl;
      #endif
      continue;
    }
    string fitOptions = "E EX0 S";
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      fitOptions = fitOptions + " ";
    #else
      fitOptions = fitOptions + " Q";
    #endif

    //##### FIT THE GRAPH #########
    TF1* powerLaw = new TF1("PowerLawFunct","[0]*pow(x,-[1])",0.0,5.0);
    powerLaw->SetParameter(0,0.1);
    powerLaw->SetParameter(1,1.0);
    powerLaw->SetParLimits(0,0,1.0);
    powerLaw->SetParLimits(1,0.001,50.0);
    powerLaw->SetNpx(300);
    powerLaw->SetLineColor(kBlue);
    TFitResultPtr resultPtr = graph_feeddown->Fit(powerLaw,fitOptions.c_str(),"",0.0,1.5);
    TF1* fitFunctNom = graph_feeddown->GetFunction("PowerLawFunct");
    if(fitFunctNom){
      if(allCentMode){
        m_FeedDown_Fits_AllCent[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunctNom->Clone();
        m_FeedDown_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("feeddownBackground_%s_AllCent_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_FeedDown_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("Feed-Down %s Cent Integrated y=[%1.2f,%1.2f] Nominal",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_FeedDown_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = (TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();
      }else{
        m_FeedDown_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunctNom->Clone();
        m_FeedDown_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("feeddownBackground_%s_Cent%02d_yIndex%02d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_FeedDown_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("Feed-Down %s Cent %d y=[%1.2f,%1.2f] Nominal",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_FeedDown_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = (TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();
      }
    }


    fitOptions = fitOptions + " + ";

    TF1* expFunct = new TF1("expFunct","[0]*exp(-[1]*x)+[2]",0.0,5.0);
    expFunct->SetParameter(0,0.1);
    expFunct->SetParameter(1,1.0);
    expFunct->SetParameter(2,0.0);
    expFunct->SetParLimits(0,0,1.0);
    expFunct->SetParLimits(1,0.0001,50.0);
    expFunct->SetParLimits(2,0.0,1.0);
    expFunct->SetNpx(300);
    expFunct->SetLineColor(kRed);
    TFitResultPtr resultPtrSys = graph_feeddown->Fit(expFunct,fitOptions.c_str(),"",0.0,1.5);
    TF1* fitFunctSys = graph_feeddown->GetFunction("expFunct");

    if(fitFunctSys){
      if(allCentMode){
        m_FeedDown_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();
        m_FeedDown_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetName(Form("feeddownBackground_%s_AllCent_yIndex%02d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_FeedDown_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetTitle(Form("Feed-Down %s Cent Integrated y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_FeedDown_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TMatrixDSym*) resultPtrSys->GetCovarianceMatrix().Clone();
      }else{
        m_FeedDown_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunctSys->Clone();
        m_FeedDown_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("feeddownBackground_%s_Cent%02d_yIndex%02d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_FeedDown_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("Feed-Down %s Cent %d y=[%1.2f,%1.2f] Systematic",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_FeedDown_CovMat_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TMatrixDSym*) resultPtrSys->GetCovarianceMatrix().Clone();
      }
    }



    SplineFitter* fitter = NULL;
    TF1* fitFunction = NULL;
    TGraph* splinePoints = NULL;
    double chiSqr = 0.0;
    double ndof = 0.0;
    if(doSplineFitting){
      fitter = new SplineFitter();
      fitter->setGraphToFit(graph_feeddown);
      //double xLocs[4];
      //xLocs[0] = 0.0; xLocs[1] = 0.25; xLocs[2] = 0.5; xLocs[3] = 2;
      //fitter->analyzeGraph(4,fitFunction, splinePoints,chiSqr,ndof,true,false,xLocs);
      fitter->analyzeGraph(3,fitFunction, splinePoints,chiSqr,ndof,false,false,NULL,false,0.0);
      if(fitFunction) fitFunction->SetNpx(300);
      if(fitFunction) fitFunction->SetLineColor(kViolet);
    }


    fittingCanvas->cd();
    gPad->SetRightMargin(0.45);
    double lowYRange = 0;
    double highYRange = 0;
    HistogramUtilities::verticalDrawRange(graph_feeddown, lowYRange,highYRange);
    TH1F *frame = fittingCanvas->DrawFrame(0.0,lowYRange,1.5,highYRange);
    if(allCentMode){
      frame->SetTitle(Form("%s Feed-Down Fraction y=[%.2f,%.2f] 0-80%% Central",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),
                rapidity_bin_axis->GetBinLowEdge(yIndex+1), rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
    }else{
      frame->SetTitle(Form("%s Feed-Down Fraction y=[%.2f,%.2f] Cent %02d",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),
                rapidity_bin_axis->GetBinLowEdge(yIndex+1), rapidity_bin_axis->GetBinLowEdge(yIndex+2),a_centIndex));
    }
    frame->GetYaxis()->SetTitle("Feed-Down Fraction");
    frame->GetXaxis()->SetTitle("m_{T}-m_{0} (GeV/c^{2})");
    graph_feeddown->Draw("Psame");
    if(fitFunction) fitFunction->Draw("same");
    if(splinePoints) splinePoints->Draw("Psame");
    gPad->Update();
    gSystem->ProcessEvents();
    if(!allCentMode){
      fittingCanvas->SaveAs(Form("%s/%s/FeedDown_ByCent/Cent%02d/feeddown_yIndex%02d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/FeedDown_AllCent/feeddown_yIndex%02d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
    }
  }
  delete fittingCanvas;

  return;
}






//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________





void EfficiencyFitter::fitHybridFeedDown(int a_partIndex,int a_charge, int a_centIndex){

  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Call to EfficiencyFitter::fitHybridFeedDown(int " << a_partIndex << ",int " << a_charge << ", int " << a_centIndex << ") " << endl;
  #endif

  bool allCentMode = false;
  if(a_centIndex < 0){
    allCentMode = true;
    cout << "ALL CENT IS NOT WRITEN YET / Exiting" << endl;
    return;
  }

  gSystem->mkdir(Form("%s/%s",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  if(allCentMode){
    gSystem->mkdir(Form("%s/%s/FeedDown_AllCent",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  }else{
    gSystem->mkdir(Form("%s/%s/HybridFeedDown/",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
    gSystem->mkdir(Form("%s/%s/HybridFeedDown/FeedDown_Cent%02d",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex));
  }


  /*TCanvas* fittingCanvas = new TCanvas("fittingCanvas_feeddownHybrid","Fitting Canvas",1300,900);
  //fittingCanvas->SetWindowSize(1300,1300);
  //fittingCanvas->SetCanvasSize(1200,1200);
  fittingCanvas->Draw();
  TPad* topPad = new TPad("topPad","topPad",0.0,0.0,1.0,0.7);
  TPad* bottomPad = new TPad("bottomPad","bottomPad",0.0,0.7,1.0,1.0);
  topPad->Draw();
  bottomPad->Draw();
  //fittingCanvas->Divide(1,2,0,0);*/


  unsigned int pmIndex = 0;
  if(a_charge == -1) pmIndex = 1;
  TAxis* rapidity_bin_axis = m_hybridFeedDown_rapidity_structure[a_partIndex]->GetXaxis();

  for(int rapBin = 1; rapBin <= g_EffFit_MaxNumRapBins; rapBin++){
    if(m_testingMidRapidity && rapBin != rapidity_bin_axis->FindBin(m_ycm)) continue;
    int yIndex = rapBin - 1;
    double rapidity = rapidity_bin_axis->GetBinCenter(rapBin);
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - rapidity) < 0.005){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating rapidity = " << rapidity << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }






    TH1D* feedDownHisto = m_FeedDown_Hybrid_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
    TGraphErrors* feedDownGraph = m_FeedDown_Hybrid_ByCent_graphs[a_partIndex][a_centIndex][yIndex][pmIndex];
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "\n\n#####################################################################################################################" << endl;
      cout << "Fitting Feed-Down " << m_partInfo->GetParticleName(a_partIndex,a_charge) << " Cent " << a_centIndex << " y: " << Form("%1.2f",rapidity) << "  ptr: " << feedDownHisto <<  endl;
    #endif
    if(!feedDownHisto){
      cout << "\nWARNING: No FeedDown Histo to fit for " << m_partInfo->GetParticleName(a_partIndex,a_charge) << "  Cent: " << a_centIndex << "  yIndex: " << yIndex << "\n" << endl;
      continue;
    }
    feedDownHisto->SetLineWidth(2);
    feedDownHisto->SetLineColor(kBlack);
    feedDownGraph->SetLineColor(kRed+1);
    int numFilledBins = HistogramUtilities::numBinsWithContentAboveThreshold(feedDownHisto, 0.0);
    if(numFilledBins < 12) continue;
    //HistogramUtilities::removeHighCountCorruptData(m_FeedDown_Hybrid_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]);
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      HistogramUtilities::printHistoInfo(feedDownHisto);
    #endif

    int firstGoodBin = HistogramUtilities::firstBinWithContent(feedDownHisto);
    int lastGoodBin  = HistogramUtilities::lastBinWithContent(feedDownHisto);
    double highEdge  = feedDownHisto->GetBinLowEdge(lastGoodBin+1);

    //TH1D* originalFeedDownHisto = feedDownHisto;
    //TH1D* truncHisto = HistogramUtilities::truncateHistogram(feedDownHisto,0.0,highEdge);
    //if(numFilledBins > 25){
      //feedDownHisto = HistogramUtilities::rebinEqualByCDF(truncHisto,25,true,true);
      //delete truncHisto;
    //}else{
      //feedDownHisto = truncHisto;
    //}
    //m_FeedDown_Hybrid_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = feedDownHisto;
    //delete originalFeedDownHisto;

    //cout << " lastGoodBin: " << lastGoodBin << endl;
    //cout << " Bins " << originalFeedDownHisto->GetNbinsX() << " -> " << truncHisto->GetNbinsX() << " -> " << feedDownHisto->GetNbinsX() << endl;


    //######  FIT top 50 with a exponential #####
    double x_high = feedDownHisto->GetBinCenter(lastGoodBin);
    double y_high = feedDownHisto->GetBinContent(lastGoodBin);

    //int halfWayBin = feedDownHisto->GetXaxis()->FindBin(x_high/2.0);
    //double halfWayLowEdge = feedDownHisto->GetBinLowEdge(halfWayBin);
    int halfWayBin = 1;
    double halfWayLowEdge = 0.0;
    double x_mid = feedDownHisto->GetBinCenter(halfWayBin);
    double y_mid = feedDownHisto->GetBinContent(halfWayBin);

    if(x_mid == x_high) cout << " ERROR: x locations are the same. Imposible to seed high pT exp for feed down" << endl;

    double expon = log(y_mid/y_high)/(x_mid-x_high);



    //TF1* topHalfExpFunct = PhysMath::getSimpleExponentialFunction("expHighPtFeedDown", amp, expon);

    /*
    TF1* topHalfExpFunct = new TF1("expHighPtFeedDown", "[0] * exp(-x/[1])/(6.283185307179586*([1]*[1]+[2]*[1]))",0.0,10.0);
    topHalfExpFunct->SetParNames("dN/dy","T","m");
    topHalfExpFunct->SetNpx(300);
    */



    double mass = m_partInfo->GetParticleMass(a_partIndex);
    double temp = -1.0/expon;
    //double initialLowFitRange[3] = {0.0,0.0,0.2};
    int    legendreOrder[3] = {4,2,4};
    //double lowFitRange[3][2] ={
    //  {0.0,0.0},
    //  {0.08,0.08},
    //  {0.2,0.2}
    //};

    double lowFitRange[3][2] ={
      {0.0,0.0},
      {0.0,0.0},
      {0.0,0.0}
    };

    double maxFitRange = highEdge;
    bool refitWithExpoIfFailed = false;
    bool doLowPtPolyFit = false;
    bool doExpo = false;
    bool doDoubleExpo = false;
    bool doTripleExpo = false;
    bool doGaus       = false;
    bool doExpoPoly4  = false;
    bool fixExpoPoly4_highOrders = false;
    if(a_partIndex == 0){
      if(a_charge > 0){
        doDoubleExpo = true;
        //doTripleExpo = true;
        //doExpoPoly4     = true;
      }else{
        doDoubleExpo = true;
        //doTripleExpo = true;
        //doExpoPoly4     = true;
      }
    }else if(a_partIndex == 1){
      if(a_charge > 0){
        //doDoubleExpo = true;
        doExpoPoly4     = true;
        fixExpoPoly4_highOrders = true;
        if(maxFitRange > 0.6) maxFitRange = 0.6;
        refitWithExpoIfFailed = true;
      }else{
        //doDoubleExpo = true;
        doExpoPoly4     = true;
        fixExpoPoly4_highOrders = true;
        if(maxFitRange > 0.6) maxFitRange = 0.6;
        refitWithExpoIfFailed = true;
      }
    }else if(a_partIndex == 2){
      if(a_charge > 0){
        //doDoubleExpo = true;
        doExpoPoly4     = true;
        //doLowPtPolyFit = true;
      }else{
        //doGaus = true;
        doExpoPoly4     = true;
      }
    }

    string fitFunctName = Form("feeddown_fit_%d_%02d_%02d_%s", a_partIndex, a_centIndex, yIndex, a_charge > 0 ? "plus" : "minus" );

    TF1* fitFunct = NULL;

    if(lowFitRange[a_partIndex][pmIndex] > x_high){
      cout << "WARNING: NO CONTENT TO FIT FEED DOWN... TRYING JUST A POLYNOMIAL FIT..." << endl;
    }

    if(doExpo){
      fitFunct = PhysMath::getFitSafeExponentialFunction(fitFunctName.c_str(), 1.0, temp);
      double newAmp = y_mid/fitFunct->Eval(x_mid);
      fitFunct->SetParameter(0,newAmp);
      fitFunct->SetParLimits(0,0.0,1000.0*newAmp);
      if(temp > 0){
        fitFunct->SetParLimits(1,0.001,100.0*temp);
      }else{
        cout << "WARNING: TEMPERATURE SEED IS NEGATIVE" << endl;
        fitFunct->SetParLimits(1,-50.0*temp,0.001);
      }
    }
    if(doDoubleExpo){

      double x_2 = feedDownHisto->GetBinCenter(lastGoodBin);
      double y_2 = feedDownHisto->GetBinContent(lastGoodBin);
      //int halfWayBin = feedDownHisto->GetXaxis()->FindBin(x_high/2.0);
      //double halfWayLowEdge = feedDownHisto->GetBinLowEdge(halfWayBin);
      int halfWayBin = (firstGoodBin + lastGoodBin)/2;
      if(halfWayBin < firstGoodBin) halfWayBin = firstGoodBin + 1;
      double x_1 = feedDownHisto->GetBinCenter(halfWayBin);
      double y_1 = feedDownHisto->GetBinContent(halfWayBin);
      double expon_1 = log(y_1/y_2)/(x_1-x_2);
      double temp_1 = -1.0/expon_1;
      double x_0 = feedDownHisto->GetBinCenter(firstGoodBin);
      double y_0 = feedDownHisto->GetBinContent(firstGoodBin);
      double expon_0 = log(y_0/y_1)/(x_0-x_1);
      double temp_0 = -1.0/expon_0;

      cout << "DoubleExpo Seeding:" << endl;
      cout << "  Bins: " << firstGoodBin << "  " << halfWayBin << "  " << lastGoodBin << endl;
      cout << "  Points: "  << x_0 << " " << y_0 << "     " << x_1 << " " << y_1 << "      " << x_2 << " " << y_2 << endl;
      cout << "  Temps: "  << temp_0 << "     " << temp_1 << endl;




      fitFunct = PhysMath::getFitSafeDoubleExponentialFunction(fitFunctName.c_str(), 1.0, temp_0, 0.0, temp_1);
      double newAmp = y_0/fitFunct->Eval(x_0);
      double valAtTwo = (y_2 - fitFunct->Eval(x_2));
      if(valAtTwo <= 0.0) valAtTwo = y_2;
      fitFunct->SetParameter(0,0.0);
      fitFunct->SetParameter(2,1.0);
      double newAmp2 = valAtTwo/fitFunct->Eval(x_2);
      fitFunct->SetParameter(0,newAmp);
      fitFunct->SetParameter(2,newAmp2);
      fitFunct->SetParLimits(0,-100*newAmp,100.0*newAmp);
      fitFunct->SetParLimits(2,-100*newAmp2,100.0*newAmp2);
      if(temp > 0){
        fitFunct->SetParLimits(1,-2.0*temp_0,100.0*temp_0);
        fitFunct->SetParLimits(3,-2.0*temp_1,100.0*temp_1);
      }else{
        cout << "WARNING: TEMPERATURE SEED IS NEGATIVE" << endl;
        fitFunct->SetParLimits(1,100.0*temp_0,-20.0*temp_0);
        fitFunct->SetParLimits(3,100.0*temp_1,-20.0*temp_1);
      }
    }
    if(doTripleExpo){
      fitFunct = PhysMath::getFitSafeTripleExponentialFunction(fitFunctName.c_str(), 1.0, temp, 0.0, 0.99*temp,0.0, 0.98*temp);
      double newAmp = y_mid/fitFunct->Eval(x_mid);
      fitFunct->SetParameter(0,newAmp/2.97);
      fitFunct->SetParameter(2,newAmp/3.01);
      fitFunct->SetParameter(4,newAmp/3.02);
      fitFunct->SetParLimits(0,-100.0*newAmp,100.0*newAmp);
      fitFunct->SetParLimits(2,-100.0*newAmp,100.0*newAmp);
      fitFunct->SetParLimits(4,-100.0*newAmp,200.0*newAmp);
      if(temp > 0){
        fitFunct->SetParLimits(1,-2.0*temp,100.0*temp);
        fitFunct->SetParLimits(3,-2.0*temp,100.0*temp);
        fitFunct->SetParLimits(5,-2.0*temp,100.0*temp);
      }else{
        cout << "WARNING: TEMPERATURE SEED IS NEGATIVE" << endl;
        fitFunct->SetParLimits(1,100.0*temp,-2.0*temp);
        fitFunct->SetParLimits(3,100.0*temp,-2.0*temp);
        fitFunct->SetParLimits(5,100.0*temp,-2.0*temp);
      }
    }
    if(doGaus){
      fitFunct = PhysMath::getGaussianFunction(fitFunctName.c_str(), 1.0, 0.0, 0.5);
      fitFunct->SetParLimits(0,0.0,100.0);
      fitFunct->SetParLimits(1,-20,20);
      fitFunct->SetParLimits(2,0.0,20.0);
    }
    if(doExpoPoly4){
      fitFunct = PhysMath::getExpoPoly4Function(fitFunctName.c_str(),1.0, 0, expon, 0, 0.0, 0.0);
      double newAmp = y_mid/fitFunct->Eval(x_mid);
      fitFunct->SetParameter(0, newAmp);
      fitFunct->SetParLimits(0,0.0,1000.0);
      if(fixExpoPoly4_highOrders){
        fitFunct->SetParLimits(3,-1000.0,0.0);
        fitFunct->FixParameter(4,0.0);
        fitFunct->FixParameter(5,0.0);
      }


    }


    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "Parameters Initialized to: " << endl;
      for(int paramIndex = 0; paramIndex < fitFunct->GetNpar(); paramIndex++){
        cout << "     " << fitFunct->GetParName(paramIndex) << "   " << fitFunct->GetParameter(paramIndex) << endl;
      }
    #endif

    fitFunct->SetRange(lowFitRange[a_partIndex][pmIndex],15.0);

    TF1* initialExpoFunct = (TF1*) fitFunct->Clone();
    initialExpoFunct->SetName("initialExpoFunct");
    initialExpoFunct->SetLineColor(kYellow);

    string fitOptions = "";
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      fitOptions = fitOptions + " ";
    #else
      fitOptions = fitOptions + " Q";
    #endif

    TFitResult* fitResult = NULL;
    fitOptions += "S";
    TFitResultPtr fitPtr = feedDownGraph->Fit(fitFunct,fitOptions.c_str(),"",lowFitRange[a_partIndex][pmIndex],maxFitRange);
    fitResult = fitPtr.Get();

    TF1* fitResultFunct =  (TF1*) feedDownGraph->GetFunction(fitFunctName.c_str());
    cout << "Eval at 0.6: " << fitResultFunct->Eval(0.6) << endl;
    if((!fitResultFunct || !std::isfinite(fitResultFunct->Eval(0.6))) && refitWithExpoIfFailed){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "WARNING: FIT WAS BAD... Doing Refit with Expon. fit options: " << fitOptions << endl;
      #endif
      delete fitFunct;
      fitFunct = PhysMath::getFitSafeExponentialFunction(fitFunctName.c_str(), 1.0, temp);
      double newAmp = y_mid/fitFunct->Eval(x_mid);
      fitFunct->SetParameter(0,newAmp);
      fitFunct->SetParLimits(0,0.0,1000.0*newAmp);
      if(temp > 0){
        fitFunct->SetParLimits(1,0.000,200.0*temp);
      }else{
        cout << "WARNING: TEMPERATURE SEED IS NEGATIVE" << endl;
        fitFunct->SetParLimits(1,-50.0*temp,0.001);
      }
      TFitResultPtr fitPtr2 = feedDownGraph->Fit(fitFunct,fitOptions.c_str(),"",lowFitRange[a_partIndex][pmIndex],maxFitRange);
      fitResult = fitPtr2.Get();
    }


    //go through the plots and remove things with pulls greater than 10 in range
    bool needsReFit = false;
    for(int bin = feedDownHisto->GetXaxis()->FindBin(lowFitRange[a_partIndex][pmIndex]);fitFunct && bin <= lastGoodBin;bin++){
      double pull = fabs(feedDownHisto->GetBinContent(bin) - fitFunct->Eval(feedDownHisto->GetBinCenter(bin)))/feedDownHisto->GetBinError(bin);
      if(feedDownHisto->GetBinError(bin) > 0.0 && pull > 10.0){
        feedDownHisto->SetBinContent(bin,0.0);
        feedDownHisto->SetBinError(bin,0.0);
        needsReFit = true;
      }
    }


    if(needsReFit){
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "ERROR: FIT WAS REALLY BAD... Doing Refit with fit options: " << fitOptions << endl;
          HistogramUtilities::printHistoInfo(feedDownHisto);
      #endif
      TFitResultPtr fitPtr2 = feedDownGraph->Fit(fitFunct,fitOptions.c_str(),"",lowFitRange[a_partIndex][pmIndex],highEdge);
      fitResult = fitPtr2.Get();
    }
    if(!fitResult){
      cout << "ERROR: Fit Failed... " << endl;
    }else{
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "TFitResult status = " << fitResult->Status()  << endl;
      #endif
    }


    /*if(a_partIndex == 4){ // not done
      TF1* fitFunct_result = (TF1*) feedDownHisto->GetFunction(fitFunctName.c_str());
      double temp_fit = fitFunct_result->GetParameter(1);
      fitFunctName += "_bw";
      fitFunct = PhysMath::blastWaveFunctPlusFitSafeExponentialFunct(fitFunctName.c_str(),1.0, temp_fit, 0.5, 0.5, mass, 0.0, temp_fit);
      //auto sets nice par limits for non amp parameters
      double bw_amp = y_mid/fitFunct->Eval(x_mid);
      fitFunct->SetParameter(0,0.0);
      fitFunct->SetParameter(5,1.0);
      fitFunct->FixParameter(4,mass);

      double exp_amp = y_mid/fitFunct->Eval(x_mid);
      fitFunct->SetParameter(0,bw_amp);
      fitFunct->SetParLimits(0,0.01*bw_amp,1E5*bw_amp);
      //fitFunct->FixParameter(5,0.0);
      //fitFunct->FixParameter(6,1.0);
      //fitFunct->SetParameter(5,0.5*exp_amp);
      //fitFunct->SetParLimits(5,0.00,5.0*exp_amp);
      //fitFunct->SetParLimits(6,0.0,5.0);

      //just a blast wave
      int status = feedDownHisto->Fit(fitFunct,fitOptions.c_str(),"",0.0,0.3);
      TF1* blastWaveFunct = (TF1*) feedDownHisto->GetFunction(fitFunctName.c_str())->Clone();
      //setup exponential correction
      double y_high_lessBW = y_high - blastWaveFunct->Eval(x_high);
      double y_mid_lessBW = y_mid - blastWaveFunct->Eval(x_mid);
      double expon_2 = log(y_mid_lessBW/y_high_lessBW)/(x_mid-x_high);
      double temp_2 = -1.0/expon_2;
      if(temp_2 < 0.0) temp_2 = 0.5;
      bw_amp = blastWaveFunct->GetParameter(0);
      blastWaveFunct->SetParameter(0,0.0);
      blastWaveFunct->SetParameter(5,1.0);
      exp_amp = y_high_lessBW/blastWaveFunct->Eval(x_high);
      blastWaveFunct->SetParameter(0,bw_amp);
      blastWaveFunct->SetParameter(5,exp_amp);
      blastWaveFunct->SetParLimits(5,0,1E5);
      blastWaveFunct->SetParameter(6,temp_2);
      blastWaveFunct->SetParLimits(6,0,10);
      fitOptions += "S";
      TFitResultPtr fitPtr = feedDownHisto->Fit(fitFunct,fitOptions.c_str(),"",0.0,0.5);
      fitResult = fitPtr.Get();
    }*/
    TMatrixDSym* expCovMat = NULL;
    TF1* expFitFunction    = NULL;
    if(fitResult){
      expCovMat = (TMatrixDSym*) fitResult->GetCovarianceMatrix().Clone();
      if(!expCovMat){
        cout << "ERROR: Fitting Failed... Cant get a covariance matrix... " << endl;
        cout << "   Can't continue because of the failure..." << endl;
      }
      expFitFunction    =  (TF1*) feedDownGraph->GetFunction(fitFunctName.c_str())->Clone();
    }

    //do a low pT fit with legendre polynomials
    TF1* optionalLegendreFit = NULL;
    TMatrixDSym* optionalLegendreFit_CovMat = NULL;
    if(doLowPtPolyFit){
      fitFunctName += "_lowPt";
      //make a graph
      TGraphAsymmErrors* feedDownGraph =   HistogramUtilities::makeGraphAsymm(feedDownHisto);
      TF1* legendreFunct = PhysMath::getLegendreExpansion(fitFunctName.c_str(), legendreOrder[a_partIndex],false,0.0, lowFitRange[a_partIndex][pmIndex]);
      PhysMathHistogramUtilities::initializeLegendreExpansion(legendreFunct, feedDownGraph,0.0, lowFitRange[a_partIndex][pmIndex], true);
      delete feedDownGraph;
      //TF1* initialLegendreSeedRes = (TF1*) legendreFunct->Clone();
      //initialLegendreSeedRes->SetName("initialLegendreSeedRes");
      //initialLegendreSeedRes->SetLineColor(kYellow);
      string fitOptions = "S EX0";
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        fitOptions = fitOptions + " ";
      #else
        fitOptions = fitOptions + " Q";
      #endif
      TFitResultPtr fitPtrRes = feedDownHisto->Fit(legendreFunct,fitOptions.c_str(),"",0.0,lowFitRange[a_partIndex][pmIndex]);
      if(!fitPtrRes.Get()){
        cout << "  ERROR: This is probably a bin you don't care about, but feeddown fitting for " << m_partInfo->GetParticleName(a_partIndex,a_charge) << " Cent: " << a_centIndex << " yIndex:" << yIndex  << " FAILED! " << endl;
        cout << "     PTR: "<< fitPtrRes.Get() <<" CONTINUING..." << endl;
      }
      int legendreStatus  = fitPtrRes->Status();
      optionalLegendreFit = feedDownHisto->GetFunction(fitFunctName.c_str());
      TMatrixDSym matrixTemp = PhysMath::addZerosToCovMatrix(fitPtrRes->GetCovarianceMatrix(), 0);
      optionalLegendreFit_CovMat = (TMatrixDSym*) matrixTemp.Clone();

      m_FeedDown_Hybrid_Fits_ByCent_lowPt[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) optionalLegendreFit->Clone();
      m_FeedDown_Hybrid_Fits_ByCent_lowPt[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("FeedDownFit_Hybrid_%s_Cent%d_yIndex%d_lowPt",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
      m_FeedDown_Hybrid_Fits_ByCent_lowPt[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("Feed-Down Spectra for %s Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                   rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      m_FeedDown_Hybrid_CovMat_ByCent_lowPt[a_partIndex][a_centIndex][yIndex][pmIndex] = optionalLegendreFit_CovMat; // names when writing, because not TObject

      //readjust the lowRange of the expo fit when they meet
      if(expFitFunction){
        double mTm0_ofClosestAgreement = lowFitRange[a_partIndex][pmIndex];
        double dev_ofClosestAgreement = expFitFunction->Eval(mTm0_ofClosestAgreement) - optionalLegendreFit->Eval(mTm0_ofClosestAgreement);
        double initialSignOfDeviation = dev_ofClosestAgreement/fabs(dev_ofClosestAgreement);
        for(double mTm0 = 0.0; mTm0 > lowFitRange[a_partIndex][pmIndex]; mTm0 += 0.005){
          double deviation = expFitFunction->Eval(mTm0) - optionalLegendreFit->Eval(mTm0);
          if(fabs(deviation) < fabs(dev_ofClosestAgreement)){
            dev_ofClosestAgreement = deviation;
            mTm0_ofClosestAgreement = mTm0;
          }
          if(deviation * initialSignOfDeviation < 0.0) break;
        }
        for(double mTm0 = mTm0_ofClosestAgreement - 0.005; mTm0 < mTm0_ofClosestAgreement + 0.005; mTm0+= 0.00001){
          double deviation = expFitFunction->Eval(mTm0) - optionalLegendreFit->Eval(mTm0);
          if(fabs(deviation) < fabs(dev_ofClosestAgreement)){
            dev_ofClosestAgreement = deviation;
            mTm0_ofClosestAgreement = mTm0;
          }
          if(deviation * initialSignOfDeviation < 0.0) break;
        }

        expFitFunction->SetRange(mTm0_ofClosestAgreement,15.0);
        m_FeedDown_Hybrid_Fits_ByCent_lowPt[a_partIndex][a_centIndex][yIndex][pmIndex]->SetRange(0.0,mTm0_ofClosestAgreement);
        optionalLegendreFit->SetRange(0.0,mTm0_ofClosestAgreement);
      }

    }



    feedDownHisto->SetTitle(Form("Feed-Down Spectra %s Cent %d y=[%1.2f,%1.2f]; m_{T}-m_{0} (GeV/c^{2}); #frac{1}{N_{Evt}}#times#frac{1}{2#pi m_{T}}#times#frac{d^{2}N_{feed-down}}{dm_{T}dy} (GeV/c^{2})^{-2}",
          m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
          rapidity_bin_axis->GetBinLowEdge(yIndex+1),
          rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
    TH1D* pullHisto = NULL;
    if(expFitFunction){
      pullHisto = HistogramUtilities::makeFitPullHistogram(feedDownHisto, expFitFunction, false);
      pullHisto->SetMarkerStyle(20);
      pullHisto->SetMarkerSize(0.5);
      pullHisto->SetMarkerColor(kBlack);
      pullHisto->SetLineColor(kBlack);
      pullHisto->SetTitle("");
    }

    TCanvas* fittingCanvas = new TCanvas("fittingCanvas_feeddownHybrid","Fitting Canvas",1300,900);
    fittingCanvas->Draw();
    TPad* topPad = new TPad("topPad","topPad",0.0,0.3,1.0,1.0);
    TPad* bottomPad = new TPad("bottomPad","bottomPad",0.0,0.0,1.0,0.3);
    topPad->Draw();
    bottomPad->Draw();

    //##########  DRAW THE FIT ############
    topPad->cd();
    //gStyle->SetOptStat(1000000001); // just name of histogram
    //gStyle->SetOptTitle(1); // print name of histogram
    gPad->SetLogy();
    double lowVertRange;
    double highVertRange;
    HistogramUtilities::verticalDrawRange(feedDownGraph, lowVertRange,highVertRange, -999,-999, true);

    feedDownHisto->Draw();
    feedDownHisto->GetXaxis()->SetRangeUser(0.0,highEdge);
    feedDownHisto->GetYaxis()->SetRangeUser(lowVertRange,highVertRange);
    if(feedDownGraph) feedDownGraph->Draw("P Same");
    if(initialExpoFunct) initialExpoFunct->Draw("same");
    if(expFitFunction) expFitFunction->Draw("same");
    if(optionalLegendreFit) optionalLegendreFit->Draw("same");
    gPad->Update();


    bottomPad->cd();
    //gStyle->SetOptStat(0);
    //gStyle->SetOptTitle(0);
    gPad->SetLogy(false);
    if(pullHisto){
      pullHisto->Draw("E");
      pullHisto->GetXaxis()->SetRangeUser(0.0,highEdge);
      pullHisto->GetYaxis()->SetRangeUser(-15.0,15.0);
    }
    TLine* oneLine = new TLine(0.0,0.0,highEdge,0.0);
    oneLine->SetLineStyle(7);
    oneLine->SetLineColor(kBlack);
    oneLine->Draw("same");
    bottomPad->Update();

    fittingCanvas->Update();
    gSystem->ProcessEvents();

    if(allCentMode){
      fittingCanvas->SaveAs(Form("%s/%s/FeedDown_AllCent/FeedDown_yIndex_%02d.png",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
    }else{
      fittingCanvas->SaveAs(Form("%s/%s/HybridFeedDown/FeedDown_Cent%02d/FeedDown_yIndex_%02d.png",
        m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
    }
    if(expFitFunction){
      m_FeedDown_Hybrid_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) expFitFunction->Clone();
      m_FeedDown_Hybrid_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("FeedDownFit_Hybrid_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
      m_FeedDown_Hybrid_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("Feed-Down Spectra for %s Cent %d y=[%1.2f,%1.2f]",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                   rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
      m_FeedDown_Hybrid_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = expCovMat; // names when writing, because not TObject
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "Histogram after fit: " << endl;
      HistogramUtilities::printHistoInfo(feedDownHisto);
    #endif

    //delete stuff
    delete oneLine;
    if(initialExpoFunct) delete initialExpoFunct;
    if(pullHisto) delete pullHisto;
    delete fittingCanvas;



  }// rapidity lop



}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


void EfficiencyFitter::fitProtonKnockout(int a_centIndex, double a_maxDCA){
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Call to EfficiencyFitter::fitProtonKnockout(int " <<a_centIndex << ", double " << a_maxDCA << ")"<< endl;
  #endif
  bool drawInitializations = false;
  bool fixAnitProtonScalingParameter = false;
  bool fixExponent = true;
  double exponentFixedValue = 0.01;
  double maxDcaForFit = 3.0;


  gSystem->mkdir(Form("%s/%s",m_imageDir.c_str(),m_partInfo->GetParticleName(2,1).Data()));
  gSystem->mkdir(Form("%s/%s/KnockoutProton_ByCent",m_imageDir.c_str(),m_partInfo->GetParticleName(2,1).Data()));
  gSystem->mkdir(Form("%s/%s/KnockoutProton_ByCent/Cent%02d",m_imageDir.c_str(),m_partInfo->GetParticleName(2,1).Data(),a_centIndex));
  gSystem->mkdir(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/DCA_Fits",m_imageDir.c_str(),m_partInfo->GetParticleName(2,1).Data(),a_centIndex));


  //double dcaScaleLowRange  = 0.8;
  //double dcaScaleHighRange = 1.8;


  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  ROOT::Math::Functor chiSquaredFunctor(this,&EfficiencyFitter::getChiSquared_DCA,4);
  minimizer->SetFunction(chiSquaredFunctor);
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "Created Minimizer " << minimizer << "   functor: " << &EfficiencyFitter::getChiSquared_DCA << endl;
    cout << "created chi squared functor with NDim " << chiSquaredFunctor.NDim() << endl;
  #endif




  for(int centIndex = 0; centIndex < 16; centIndex++){
    if(a_centIndex >= 0 && centIndex != a_centIndex) continue;
    m_currentCentIndex = centIndex;

    TH2D* protonKnockoutHisto2D = (TH2D*) m_proton_knockout_histogram_ByCent[centIndex]->Clone();
    HistogramUtilities::zeroBins2D(protonKnockoutHisto2D);
    TH2D* protonTotalHisto2D = (TH2D*) protonKnockoutHisto2D->Clone();
    protonKnockoutHisto2D->SetName(Form("knock2D_%d",centIndex));
    protonTotalHisto2D->SetName(Form("pro2D_%d",centIndex));


    for(unsigned int yIndex = 0; yIndex < g_EffFit_MaxNumRapBins; yIndex++){

      m_currentRapIndex = yIndex;
      double rapidity = m_proton_dca_ChiSqrd[centIndex]->GetXaxis()->GetBinLowEdge(yIndex+1);

      TH1D* protonKnockoutHisto = m_proton_dca_ChiSqrd[centIndex]->ProjectionY(Form("knock_%d_%d",centIndex,yIndex),yIndex,yIndex);
      TH1D* protonTotalHisto    = m_proton_dca_ChiSqrd[centIndex]->ProjectionY(Form("pro_%d_%d",centIndex,yIndex),yIndex,yIndex);
      HistogramUtilities::zeroBins(protonKnockoutHisto);
      HistogramUtilities::zeroBins(protonTotalHisto);
      TH1D* parameterHistos[4];
      string parNames[4] = {"frac","A","DCA_0","exp"};
      for(int iii = 0; iii < 4; iii++){
        parameterHistos[iii] = (TH1D*) protonTotalHisto->Clone();
        parameterHistos[iii]->SetName(Form("knockout_param_%d_%d_%d",iii,centIndex,yIndex));
        parameterHistos[iii]->SetTitle(Form("Knockout %s Cent %d y=[%1.3f,%1.3f]; m_{T}-m_{p} [GeV/c^{2}]; %s",parNames[iii].c_str(),
                         centIndex,m_proton_dca_ChiSqrd[centIndex]->GetXaxis()->GetBinLowEdge(yIndex+1),
                         m_proton_dca_ChiSqrd[centIndex]->GetXaxis()->GetBinLowEdge(yIndex+2),parNames[iii].c_str()));
      }
      TF1* parameterTrendFunctions[4] = {NULL,NULL,NULL,NULL};

      int fittingRound = 0;
      double start_mTm0  = PhysMath::mTm0(m_partInfo->GetParticleMass(2),0.33);
      int    start_index = m_proton_dca_ChiSqrd[0]->GetYaxis()->FindBin(start_mTm0) - 1;
      vector<double> mTm0IndexArray;
      for(int iii = start_index; iii >= 0; iii--){
        mTm0IndexArray.push_back(iii);
      }
      for(int iii = start_index + 1; iii < m_proton_dca_ChiSqrd[0]->GetNbinsY(); iii++){
        mTm0IndexArray.push_back(iii);
      }
      // rerun for second round of fitting
      for(int iii = 0; iii < m_proton_dca_ChiSqrd[0]->GetNbinsY(); iii++){
        mTm0IndexArray.push_back(iii);
      }


      double first_fit_bkgd_params[5]    = {-999,-999,-999,-999,-999};
      double previous_fit_bkgd_params[5] = {-999,-999,-999,-999,-999};
      int numSuccessfulDcaFits = 0;
      for(unsigned int mTm0IndexIndex = 0; mTm0IndexIndex < g_EffFit_MaxNumMtM0Bins; mTm0IndexIndex++){
        unsigned int mTm0Index = mTm0IndexArray[mTm0IndexIndex];
        double mTm0 = m_proton_dca_ChiSqrd[centIndex]->GetYaxis()->GetBinCenter(mTm0Index+1);
        double pT = PhysMath::pT_from_mTm0(mTm0,m_partInfo->GetParticleMass(2));
        if(pT > 0.38) continue;
        if(mTm0Index == 0){
          fittingRound += 1;
          if(fittingRound == 2){
            if(numSuccessfulDcaFits > 4){
              parameterTrendFunctions[0] = fitKnockoutTrendHistogram(0,parameterHistos[0]);
              parameterTrendFunctions[1] = fitKnockoutTrendHistogram(1,parameterHistos[1]);
              parameterTrendFunctions[2] = fitKnockoutTrendHistogram(2,parameterHistos[2]);
              if(!fixExponent) parameterTrendFunctions[3] = fitKnockoutTrendHistogram(3,parameterHistos[3]);
            }else{
              break;
            }
          }
        }
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "######### KNOCKOUT PROTON DCA:  Cent " << centIndex << " y " << yIndex << " mTm0 " << mTm0Index << endl;
        #endif

        m_currentProtonDCA     = m_proton_dca_ByCent[centIndex][yIndex][mTm0Index][0];
        m_currentAntiProtonDCA = m_proton_dca_ByCent[centIndex][yIndex][mTm0Index][1];
        if(!m_currentProtonDCA){
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "Knockout cent: " << centIndex << "  rap: " << yIndex << "  mTm0: " << mTm0Index << " no histos, continuing..." << endl;
          #endif
          continue;
        }


        m_currentProtonDCA->SetLineColor(kBlue);
        m_currentAntiProtonDCA->SetLineColor(kRed);
        m_currentProtonDCA->SetBinContent(m_currentProtonDCA->GetNbinsX(),0.0);
        m_currentProtonDCA->SetBinError(m_currentProtonDCA->GetNbinsX(),0.0);
        m_currentProtonDCA->Scale(1.0,"width");
        m_currentAntiProtonDCA->SetBinContent(m_currentAntiProtonDCA->GetNbinsX(),0.0);
        m_currentAntiProtonDCA->SetBinError(m_currentAntiProtonDCA->GetNbinsX(),0.0);
        m_currentAntiProtonDCA->Scale(1.0,"width");

        minimizer->Clear();

        double totalProtonContent      = HistogramUtilities::sumOfContent(m_currentProtonDCA,-1,true);
        double totalAntiProtonContent  = HistogramUtilities::sumOfContent(m_currentAntiProtonDCA,-1,true);
        int    maxBin_proton       = HistogramUtilities::getMaxBinInRange(m_currentProtonDCA);
        double proton_content      = m_currentProtonDCA->GetBinContent(maxBin_proton);
        double antiproton_content  = m_currentAntiProtonDCA->GetBinContent(maxBin_proton);
        if(totalProtonContent == 0 || totalAntiProtonContent == 0){
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "WARNING: Proton dca fit cant initialize... has counts p:" << proton_content << " pbar:" << antiproton_content <<  " continuing..." << endl;
          #endif
          continue;
        }
        double scale_predict = proton_content/antiproton_content;

        //setup scale better

        double maxDCAForScaleInit = 0.5;
        int    binOfMaxDCAForScaleInit = m_currentProtonDCA->FindBin(maxDCAForScaleInit);
        int    lowProtonBin     = HistogramUtilities::firstBinWithContent(m_currentProtonDCA);
        int    lowAntiProtonBin = HistogramUtilities::firstBinWithContent(m_currentAntiProtonDCA);
        int    lowBinForBoth = max(lowProtonBin,lowAntiProtonBin);
        int    highProtonBin     = HistogramUtilities::lastBinWithContent(m_currentProtonDCA);
        int    highAntiProtonBin = HistogramUtilities::lastBinWithContent(m_currentAntiProtonDCA);
        int    highBinForBoth = min(highProtonBin,highAntiProtonBin);
        if(fittingRound < 2){
          double scale_sum  = 0.0;
          double weight_sum = 0.0;
          for(int bin = 0; bin <= binOfMaxDCAForScaleInit; bin++){
            double pro      = m_currentProtonDCA->GetBinContent(bin);
            double anti     = m_currentAntiProtonDCA->GetBinContent(bin);
            double pro_err  = m_currentProtonDCA->GetBinError(bin);
            double anti_err = m_currentAntiProtonDCA->GetBinError(bin);
            if(pro > 0.0 && anti > 0.0){
              double scale  = pro/anti;
              double weight = 1.0/(scale*scale*(pro_err*pro_err/(pro*pro) + anti_err*anti_err/(anti*anti)));
              scale_sum  += scale*weight;
              weight_sum += weight;
            }
          }
          scale_predict = scale_sum/weight_sum;
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << " scale: " << scale_predict << "  from sum: " << scale_sum << " and weight: " << weight_sum << endl;
          #endif

        }else{
          scale_predict = parameterTrendFunctions[0]->Eval(mTm0);
        }



        double dca_ofTest = 2.5;
        int    bin_ofTest                = m_currentProtonDCA->FindBin(dca_ofTest);
        double proton_content_ofTest     = m_currentProtonDCA->GetBinContent(bin_ofTest);
        double antiproton_content_ofTest = m_currentAntiProtonDCA->GetBinContent(bin_ofTest);
        if(proton_content_ofTest <= 0.0){
          for(int bin = bin_ofTest - 1; bin > 0; bin--){
            bin_ofTest = bin;
            proton_content_ofTest     = m_currentProtonDCA->GetBinContent(bin_ofTest);
            antiproton_content_ofTest = m_currentAntiProtonDCA->GetBinContent(bin_ofTest);
            if(proton_content_ofTest > 0.0) break;
          }
        }
        if(proton_content_ofTest <= 0 ){
          cout << " WARNING: Antiproton histo is above proton histo.. cant initialize background function... continuing..." << endl;
          continue;
        }

        m_currentDCALowerFitRangeBin = 1; //maxBin_proton
        m_currentDCAUpperFitRangeBin = m_currentProtonDCA->FindBin(maxDcaForFit);
        if(m_currentDCAUpperFitRangeBin > highBinForBoth) m_currentDCAUpperFitRangeBin = highBinForBoth;
        double lowFitRange  = m_currentProtonDCA->GetBinLowEdge(m_currentDCALowerFitRangeBin);
        double highFitRange = m_currentProtonDCA->GetBinLowEdge(m_currentDCAUpperFitRangeBin+1);

        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "  Setting Fit Range from DCA=[" <<lowFitRange << "," << highFitRange <<"]" << endl;
        #endif

        double amp_predict   = proton_content_ofTest - scale_predict*antiproton_content_ofTest;
        double dca0_predict  = 0.855;
        double alpha_predict = 0.085;
        if(fittingRound == 2){
          //amp_predict   = parameterTrendFunctions[1]->Eval(mTm0);
          if(parameterTrendFunctions[2]) dca0_predict  = parameterTrendFunctions[2]->Eval(mTm0);
          if(parameterTrendFunctions[3]) alpha_predict = parameterTrendFunctions[3]->Eval(mTm0);
        }else{
          if(previous_fit_bkgd_params[0] > 0){
            dca0_predict  = previous_fit_bkgd_params[1];// (this is not the index of minimizer parameters)
            alpha_predict = previous_fit_bkgd_params[2];
          }
        }

        TF1* initBkgdFunct = new TF1("initBkgdFunct","[0]*pow(1.0 - exp(-x/[1]),1.0/[2])",0.0,3.0);
        initBkgdFunct->SetParNames("Amp","DCA_0","alpha");
        initBkgdFunct->SetNpx(500);
        initBkgdFunct->SetParameter(0,1.0);
        initBkgdFunct->SetParameter(1,dca0_predict);// (this is not the index of minimizer parameters)
        initBkgdFunct->SetParameter(2,alpha_predict);
        initBkgdFunct->SetLineColor(kYellow+1);

        amp_predict = amp_predict/initBkgdFunct->Eval(dca_ofTest);
        if(amp_predict < 0){
          amp_predict = 10;
          if(previous_fit_bkgd_params[0] > 0){
            amp_predict = previous_fit_bkgd_params[0];
          }
        }
        initBkgdFunct->SetParameter(0,amp_predict);
        bool initFunctBad = !(std::isfinite(amp_predict) && std::isfinite(dca0_predict) && std::isfinite(alpha_predict));

        TH1D* combinedDca_init    = (TH1D*) m_currentAntiProtonDCA->Clone();
        combinedDca_init->SetLineColor(kYellow+1);
        HistogramUtilities::zeroBins(combinedDca_init);
        for( int binX = 1; binX <= combinedDca_init->GetNbinsX(); binX++){
          double center  = m_currentAntiProtonDCA->GetBinCenter(binX);
          double content = m_currentAntiProtonDCA->GetBinContent(binX);
          double error   = m_currentAntiProtonDCA->GetBinError(binX);
          double newContent = content*scale_predict;
          double newError   = 0.0;
          if( std::isfinite(newContent) && std::isfinite(newError) ){
            newContent += initBkgdFunct->Eval(center);
            combinedDca_init->SetBinContent(binX,newContent);
          }
        }

        /*double integralPlus  = m_currentProtonDCA->Integral(    m_currentProtonDCA->GetXaxis()->FindBin(dcaScaleLowRange),    m_currentProtonDCA->GetXaxis()->FindBin(dcaScaleHighRange));
        double integralMinus = m_currentAntiProtonDCA->Integral(m_currentAntiProtonDCA->GetXaxis()->FindBin(dcaScaleLowRange),m_currentAntiProtonDCA->GetXaxis()->FindBin(dcaScaleHighRange));
        if(integralPlus == 0 || integralMinus == 0){
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "Neither proton dca histo has counts... continuing..." << endl;
          #endif
          continue;
        }
        double scaleFactor = integralPlus/integralMinus;
        double scaleError = scaleFactor*sqrt(1.0/integralPlus + 1.0/integralMinus);
        */

        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << "Minimizer Seeds: " << endl;
          cout << "  scale_predict: " << scale_predict << endl;
          cout << "  amp_predict: "   << amp_predict   << endl;
          cout << "  dca0_predict: "  << dca0_predict  << endl;
          cout << "  alpha_predict: " << alpha_predict << endl;
        #endif

        //############ START THE FIT ################
        if(fixAnitProtonScalingParameter){
          minimizer->SetFixedVariable(0,"p_pbar_frac", scale_predict);
        }else{
          if(fittingRound == 2){
            minimizer->SetFixedVariable(0,"p_pbar_frac", parameterTrendFunctions[0]->Eval(parameterHistos[0]->GetBinCenter(mTm0Index+1)));
          }else{
            minimizer->SetVariable(0,"p_pbar_frac", scale_predict, 1E-2);
            minimizer->SetVariableLimits(0,0.01*scale_predict,100.0*scale_predict);
          }
        }
        minimizer->SetVariable(1,"A_bkgd", amp_predict, 1E-2);
        minimizer->SetVariableLimits(1,-1000.0*amp_predict,1000.0*amp_predict);
        if(fittingRound == 2){
          minimizer->SetFixedVariable(2,"DCA_naught", dca0_predict);
        }else{
          minimizer->SetVariable(2,"DCA_naught", dca0_predict, 1E-5);
          minimizer->SetVariableLimits(2,0.0,3.0);
        }


        if(fixExponent){
          minimizer->SetFixedVariable(3,"alpha", exponentFixedValue);
        }else{
          minimizer->SetVariable(3,"alpha", alpha_predict, 1E-5);
          minimizer->SetVariableLimits(3,0.0,1.0);
        }


        #ifdef _EFFICIENCY_FITTER_DEBUG_
          minimizer->SetPrintLevel(1);
        #else
          minimizer->SetPrintLevel(0);
        #endif
        minimizer->SetMaxIterations(3000);
        minimizer->SetMaxFunctionCalls(100000);
        minimizer->SetValidError(true);
        minimizer->Minimize();
        int status = minimizer->Status();
        if(status != 0){
          cout << "ERROR: FITTING FAILED..." << endl;
        }else{
          minimizer->Hesse();
        }



        if(status == 0){
          numSuccessfulDcaFits += 1;
          const double* finalParams = minimizer->X();
          double scaleError_Final = sqrt(minimizer->CovMatrix(0,0));
          m_proton_dca_ChiSqrd[centIndex]->SetBinContent(yIndex+1,mTm0Index+1,getChiSquared_DCA(finalParams));

          TF1* backgroundFunct = new TF1(Form("backgroundKnockoutFunction_Cent%02d_yIndex%02d_mTm0Index%02d",centIndex,yIndex,mTm0Index),"[0]*pow(1.0 - exp(-x/[1]),1.0/[2])",0.0,3.0);
          backgroundFunct->SetParNames("Amp","DCA_0","alpha");
          backgroundFunct->SetNpx(500);
          backgroundFunct->SetParameter(0,finalParams[1]);
          backgroundFunct->SetParameter(1,finalParams[2]);//not aligned on purpose
          backgroundFunct->SetParameter(2,finalParams[3]);
          backgroundFunct->SetLineColor(kBlack);
          m_proton_dca_Fits_ByCent[centIndex][yIndex][mTm0Index] = backgroundFunct;

          if(mTm0IndexIndex == 0){
            first_fit_bkgd_params[0] = finalParams[1];
            first_fit_bkgd_params[1] = finalParams[2];
            first_fit_bkgd_params[2] = finalParams[3];
          }
          previous_fit_bkgd_params[0] = finalParams[1];
          previous_fit_bkgd_params[1] = finalParams[2];
          previous_fit_bkgd_params[2] = finalParams[3];
          if(mTm0IndexIndex > 1 && mTm0IndexArray[mTm0IndexIndex - 1] > mTm0Index && mTm0IndexArray[mTm0IndexIndex + 1] > mTm0Index){
            previous_fit_bkgd_params[0] = first_fit_bkgd_params[0];
            previous_fit_bkgd_params[1] = first_fit_bkgd_params[1];
            previous_fit_bkgd_params[2] = first_fit_bkgd_params[2];
          }

          parameterHistos[0]->SetBinContent(mTm0Index+1,finalParams[0]);
          parameterHistos[0]->SetBinError(  mTm0Index+1,sqrt(minimizer->CovMatrix(0,0)));
          parameterHistos[1]->SetBinContent(mTm0Index+1,finalParams[1]/totalProtonContent);
          parameterHistos[1]->SetBinError(  mTm0Index+1,sqrt(minimizer->CovMatrix(1,1))/totalProtonContent);
          parameterHistos[2]->SetBinContent(mTm0Index+1,finalParams[2]);
          parameterHistos[2]->SetBinError(  mTm0Index+1,sqrt(minimizer->CovMatrix(2,2)));
          parameterHistos[3]->SetBinContent(mTm0Index+1,finalParams[3]);
          parameterHistos[3]->SetBinError(  mTm0Index+1,sqrt(minimizer->CovMatrix(3,3)));

          TH1D* scaledDcaMinus_Final = (TH1D*) m_currentAntiProtonDCA->Clone();
          TH1D* combinedDca_Final    = (TH1D*) m_currentAntiProtonDCA->Clone();
          scaledDcaMinus_Final->SetLineColor(kViolet);
          combinedDca_Final->SetLineColor(kGreen+2);
          HistogramUtilities::zeroBins(scaledDcaMinus_Final);
          HistogramUtilities::zeroBins(combinedDca_Final);
          for( int binX = 1; binX <= scaledDcaMinus_Final->GetNbinsX(); binX++){
            double center  = m_currentAntiProtonDCA->GetBinCenter(binX);
            double content = m_currentAntiProtonDCA->GetBinContent(binX);
            double error   = m_currentAntiProtonDCA->GetBinError(binX);
            double newContent = content*finalParams[0];
            double newError   = newContent * sqrt( error*error/(content*content) + scaleError_Final*scaleError_Final/(finalParams[0]*finalParams[0]));
            if( std::isfinite(newContent) && std::isfinite(newError) ){
              scaledDcaMinus_Final->SetBinContent(binX,newContent);
              scaledDcaMinus_Final->SetBinError(binX,newError);
              newContent += backgroundFunct->Eval(center);
              combinedDca_Final->SetBinContent(binX,newContent);
            }
          }

          double minProtonVal          = HistogramUtilities::getMinInRange(m_currentProtonDCA,0,2.9);
          double minAntiProtonVal      = HistogramUtilities::getMinInRange(m_currentAntiProtonDCA,0,2.9);
          double minScaleAntiProtonVal = HistogramUtilities::getMinInRange(scaledDcaMinus_Final,0,2.9);
          double maxProtonVal          = HistogramUtilities::getMaxInRange(m_currentProtonDCA,0.0,2.9);
          double maxAntiProtonVal      = HistogramUtilities::getMaxInRange(m_currentAntiProtonDCA,0.0,2.9);
          double maxScaleAntiProtonVal = HistogramUtilities::getMaxInRange(scaledDcaMinus_Final,0.0,2.9);

          double smallestVal = minProtonVal;
          //if(minAntiProtonVal < smallestVal)      smallestVal = minAntiProtonVal;
          if(minScaleAntiProtonVal < smallestVal) smallestVal = minScaleAntiProtonVal;
          double largestVal = maxProtonVal;
          //if(maxAntiProtonVal > largestVal)      largestVal = maxAntiProtonVal;
          if(maxScaleAntiProtonVal > largestVal) largestVal = maxScaleAntiProtonVal;
          smallestVal = 0.1*smallestVal;
          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << " Y Axis Range = [" << smallestVal << "," << largestVal << "]" << endl;;
            cout << "  from " << minProtonVal << " " << minAntiProtonVal << " " << minScaleAntiProtonVal << endl;
            cout << "  from " << maxProtonVal << " " << maxAntiProtonVal << " " << maxScaleAntiProtonVal << endl;
          #endif
          if(smallestVal <= 0) smallestVal = 0.9;

          TLine* lowFitRangeLine = new TLine(lowFitRange,smallestVal,lowFitRange,largestVal);
          TLine* highFitRangeLine = new TLine(highFitRange,smallestVal,highFitRange,largestVal);
          lowFitRangeLine->SetLineWidth(2);
          highFitRangeLine->SetLineWidth(2);
          lowFitRangeLine->SetLineColor(kRed);
          highFitRangeLine->SetLineColor(kRed);


          TGraphAsymmErrors* graph1 = HistogramUtilities::makeGraphAsymm(m_currentProtonDCA);
          TGraphAsymmErrors* graph2 = HistogramUtilities::makeGraphAsymm(combinedDca_Final);

          TGraphAsymmErrors* pullGraph = HistogramUtilities::makeFitPullGraph(graph1, graph2);
          pullGraph->SetMarkerStyle(20);
          pullGraph->SetMarkerSize(0.5);
          pullGraph->SetMarkerColor(kBlack);
          pullGraph->SetLineColor(kBlack);





          TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
          fittingCanvas->SetWindowSize(1300,700);
          fittingCanvas->SetCanvasSize(1200,600);
          fittingCanvas->Draw();
          TPad* topPad_base = new TPad("topPad_Eff_base","topPad_Eff",0.0,0.3,1.0,1.0);
          TPad* bottomPad_base = new TPad("bottomPad_Eff_base","bottomPad_Eff",0.0,0.0,1.0,0.3);
          topPad_base->Draw();
          topPad_base->SetRightMargin(0.38);
          bottomPad_base->Draw();
          bottomPad_base->SetRightMargin(0.38);
          fittingCanvas->cd();


          topPad_base->cd();
          gPad->SetLogy();
          m_currentProtonDCA->Draw("E");
          m_currentProtonDCA->GetYaxis()->SetRangeUser(0.1*smallestVal,largestVal*2.0);
          if(drawInitializations){
            initBkgdFunct->Draw("same");
            combinedDca_init->Draw("same");
          }
          //m_currentAntiProtonDCA->Draw("E same");
          combinedDca_Final->Draw("same");
          scaledDcaMinus_Final->Draw("E same");
          backgroundFunct->Draw("same");
          lowFitRangeLine->Draw("same");
          highFitRangeLine->Draw("same");

          TPaveText* fitParsTxt = new TPaveText(0.65,0.15,0.95,0.85,"NDC");
          fitParsTxt->AddText(Form("p/#bar{p}  %.2e#pm%.2e",finalParams[0],sqrt(minimizer->CovMatrix(0,0))));
          fitParsTxt->AddText(Form("A           %.2e#pm%.2e",finalParams[1],sqrt(minimizer->CovMatrix(1,1))));
          fitParsTxt->AddText(Form("DCA_{0}     %.2e#pm%.2e",finalParams[2],sqrt(minimizer->CovMatrix(2,2))));
          fitParsTxt->AddText(Form("#alpha      %.2e#pm%.2e",finalParams[3],sqrt(minimizer->CovMatrix(3,3))));

          fitParsTxt->SetFillStyle(0);
          fitParsTxt->SetBorderSize(5);
          fitParsTxt->Draw("same");

          bottomPad_base->cd();
          TH1F *frame_pull_sys = fittingCanvas->DrawFrame(0.0,-10.0,3.0,10.0);
          frame_pull_sys->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Standardized Pulls"));
          if(pullGraph) pullGraph->Draw("PSame");
          TLine* oneLineSys = new TLine(0.0,0.0,3.0,0.0);
          oneLineSys->SetLineStyle(7);
          oneLineSys->SetLineColor(kBlack);
          oneLineSys->Draw("same");
          gSystem->ProcessEvents();

          fittingCanvas->SaveAs(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/DCA_Fits/round%02d_yIndex_%02d_mTm0Index_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(2,1).Data(), a_centIndex,fittingRound, yIndex,mTm0Index));

          delete topPad_base;
          delete bottomPad_base;
          delete fittingCanvas;

          double covMat[16];
          minimizer->GetCovMatrix(covMat);
          double reducedCovMat[9];
          for(int iii = 0; iii < 3; iii++){
            for(int jjj = 0; jjj < 3; jjj++){
              reducedCovMat[(3*iii)+jjj] = covMat[(4*(iii+1))+jjj+1];
            }
          }

          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << "Originial Cov Mat: ";
            for(int iii = 0; iii < 4; iii++){
              cout << endl;
              for(int jjj = 0; jjj < 4; jjj++){
                cout << covMat[(4*iii)+jjj] << " ";
              }
            }
            cout << endl;
            cout << "New Cov Mat: ";
            for(int iii = 0; iii < 3; iii++){
              cout << endl;
              for(int jjj = 0; jjj < 3; jjj++){
                cout << reducedCovMat[(3*iii)+jjj] << " ";
              }
            }
            cout << endl;
            cout << "Integrating background from 0 to " << a_maxDCA << endl;
          #endif


          double background     = backgroundFunct->Integral(0.0,a_maxDCA)/a_maxDCA;
          double backgroundErr  = backgroundFunct->IntegralError(0.0,a_maxDCA, backgroundFunct->GetParameters(), reducedCovMat)/a_maxDCA;
          protonKnockoutHisto->SetBinContent(mTm0Index+1,background);
          protonKnockoutHisto->SetBinError(mTm0Index+1,backgroundErr);
          protonKnockoutHisto2D->SetBinContent(yIndex+1,mTm0Index+1,background);
          protonKnockoutHisto2D->SetBinError(yIndex+1,mTm0Index+1,backgroundErr);

          double countsOfProton = 0;
          int lastBinOfDCA      =  m_currentProtonDCA->FindBin(a_maxDCA);
          for(int binX = 1; binX <= lastBinOfDCA; binX++){
            countsOfProton += m_currentProtonDCA->GetBinContent(binX);
          }
          double countsOfProtonErr = sqrt(countsOfProton);
          protonTotalHisto->SetBinContent(mTm0Index+1,countsOfProton);
          protonTotalHisto->SetBinError(mTm0Index+1,countsOfProtonErr);
          protonTotalHisto2D->SetBinContent(yIndex+1,mTm0Index+1,countsOfProton);
          protonTotalHisto2D->SetBinError(yIndex+1,mTm0Index+1,countsOfProtonErr);

          double fraction          = background/countsOfProton;
          m_proton_knockout_histogram_ByCent[centIndex]->SetBinContent(m_currentRapIndex+1,mTm0Index+1,fraction);

          /*  // basyian errors calculated below
          double fractionErr       = fraction * sqrt( countsOfProtonErr*countsOfProtonErr/(countsOfProton*countsOfProton)
                               + backgroundErr*backgroundErr/(background*background));
          int pointIndex = m_proton_knockout_ByCent[centIndex][yIndex]->GetN();
          double mTm0    = m_knockout_mTm0_structure->GetBinCenter(mTm0Index+1);
          double mTm0HalfWidth = 0.5*m_knockout_mTm0_structure->GetBinWidth(mTm0Index+1);
          if(std::isfinite(fraction)){
            m_proton_knockout_ByCent[centIndex][yIndex]->SetPoint(pointIndex,mTm0,fraction);
            if(std::isfinite(fractionErr)){
              m_proton_knockout_ByCent[centIndex][yIndex]->SetPointError(pointIndex,mTm0HalfWidth,mTm0HalfWidth,fractionErr,fractionErr);
            }
          }*/

          #ifdef _EFFICIENCY_FITTER_DEBUG_
            cout << " Proton Count:      " << countsOfProton << "  +/-  " << countsOfProtonErr  << endl;
            cout << " Knockout Count:    " << background << "  +/-  " << backgroundErr << endl;
            cout << " Fraction Knockout: " << fraction <<  endl;
          #endif

          delete combinedDca_init;
          delete lowFitRangeLine;
          delete highFitRangeLine;
          // status is good for fit
        }else{

          double minProtonVal          = HistogramUtilities::getMinInRange(m_currentProtonDCA);
          double minAntiProtonVal      = HistogramUtilities::getMinInRange(m_currentAntiProtonDCA);
          double maxProtonVal          = HistogramUtilities::getMaxInRange(m_currentProtonDCA);
          double maxAntiProtonVal      = HistogramUtilities::getMaxInRange(m_currentAntiProtonDCA);

          double smallestVal = minProtonVal;
          if(minAntiProtonVal < smallestVal)     smallestVal = minAntiProtonVal;
          double largestVal = maxProtonVal;
          if(maxAntiProtonVal > largestVal)      largestVal = maxAntiProtonVal;

          TLine* lowFitRangeLine = new TLine(lowFitRange,smallestVal,lowFitRange,largestVal);
          TLine* highFitRangeLine = new TLine(highFitRange,smallestVal,highFitRange,largestVal);
          lowFitRangeLine->SetLineWidth(2);
          highFitRangeLine->SetLineWidth(2);
          lowFitRangeLine->SetLineColor(kRed);
          highFitRangeLine->SetLineColor(kRed);

          TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
          fittingCanvas->SetWindowSize(1300,700);
          fittingCanvas->SetCanvasSize(1200,600);
          fittingCanvas->cd();
          gPad->SetRightMargin(0.35);
          gPad->SetLogy(false);

          m_currentProtonDCA->Draw("E");
          m_currentProtonDCA->GetYaxis()->SetRangeUser(0.1*smallestVal,largestVal*2.0);
          initBkgdFunct->Draw("same");
          combinedDca_init->Draw("same");
          m_currentAntiProtonDCA->Draw("E same");
          lowFitRangeLine->Draw("same");
          highFitRangeLine->Draw("same");

          fittingCanvas->SaveAs(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/DCA_Fits/round%02d_yIndex_%02d_mTm0Index_%02d____BAD.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(2,1).Data(), a_centIndex,fittingRound, yIndex,mTm0Index));
          delete lowFitRangeLine;
          delete highFitRangeLine;

          // status is bad for fit
        }
        delete initBkgdFunct;
      }//mTm0Index loop

      if(numSuccessfulDcaFits < 4) continue;



      // ###########################################################
      // FIT KNOCKOUT FRACTION AS A FUNCTION OF mT-m0
      // ###########################################################
      m_proton_knockout_ByCent[centIndex][yIndex] = new TGraphAsymmErrors();
      m_currentGraphToFit = m_proton_knockout_ByCent[centIndex][yIndex];
      m_currentGraphToFit->SetName(Form("knockoutProtonFraction_Cent%02d_yIndex%02d",centIndex,yIndex));
      m_currentGraphToFit->SetTitle(Form("Knockout Proton Fraction Cent %02d y=[%1.3f,%1.3f]; m_{T}-m_{p} [GeV/c^{2}]; N_{knockout}/N_{total}",
                         centIndex,m_proton_dca_ChiSqrd[centIndex]->GetXaxis()->GetBinLowEdge(yIndex+1),m_proton_dca_ChiSqrd[centIndex]->GetXaxis()->GetBinLowEdge(yIndex+2)));
      m_currentGraphToFit->Divide(protonKnockoutHisto,protonTotalHisto,"cl=0.68268949 b(1,1) mode");
      m_currentGraphToFit->Sort();
      for(int pointIndex = 0; pointIndex < m_currentGraphToFit->GetN(); pointIndex++){
        double y_val      = m_currentGraphToFit->GetY()[pointIndex];
        double y_low_err  = m_currentGraphToFit->GetEYlow()[pointIndex];
        double y_high_err = m_currentGraphToFit->GetEYhigh()[pointIndex];
        if(y_low_err > y_val){
          m_currentGraphToFit->SetPointError(pointIndex, m_currentGraphToFit->GetEXlow()[pointIndex],
                                m_currentGraphToFit->GetEXhigh()[pointIndex], y_val, y_high_err);
        }
        if(y_val <= 0.0){
          m_currentGraphToFit->RemovePoint(pointIndex);
          pointIndex--;
        }
      }
      m_currentGraphToFit->SetMarkerStyle(20);
      m_currentGraphToFit->SetMarkerSize(0.5);
      m_currentGraphToFit->SetLineColor(kBlack);
      if(m_currentGraphToFit->GetN() < 4) continue;
      double mTm0_0 = m_currentGraphToFit->GetX()[0];
      double mTm0_1 = m_currentGraphToFit->GetX()[1];
      double frac_0 = m_currentGraphToFit->GetY()[0];
      double frac_1 = m_currentGraphToFit->GetY()[1];
      /*
      double exponent_seed = -log(frac_0/frac_1)/log(mTm0_0/mTm0_1);
      if(!std::isfinite(exponent_seed) || exponent_seed < 0.0){
        exponent_seed = 8.172;
      }*/

      TF1* fitFunction_knockout = new TF1("fitFunction_knockout", "[0]*[4]*([1]-1.0)*[2]*pow([2]*x+[3],-[1])/pow([3],1.0-[1])", 0,3.0); // a power law in pT space // /pow([3],1.0-[1]) if b > 0
      //TF1* fitFunction_knockout = new TF1("fitFunction_knockout", "[0]*[4]*([1]-1.0)*[2]*pow([2]*x+[3],-[1])/pow([3],1.0-[1]) + [5]*exp(-x/[6])/[6]", 0,3.0); // a power law in pT space // /pow([3],1.0-[1]) if b > 0

      //fitFunction_knockout->SetParNames("Scale","Exponent","a","b","A_0","N","T");
      fitFunction_knockout->SetParNames("Scale","Exponent","a","b","A_0");

      fitFunction_knockout->SetNpx(250);
      fitFunction_knockout->SetParameter(0,1.0);
      fitFunction_knockout->SetParameter(1,33.0);
      fitFunction_knockout->SetParameter(2,4.2);
      fitFunction_knockout->SetParameter(3,0.51);
      fitFunction_knockout->SetParameter(4,1.0); //this is a dummy variable to set to A_0 so the fitter doesn't complain about size
      //fitFunction_knockout->SetParameter(5,7E-10); //this is a dummy variable to set to A_0 so the fitter doesn't complain about size
      //fitFunction_knockout->SetParameter(6,15.0); //this is a dummy variable to set to A_0 so the fitter doesn't complain about size


      double xValToUse = m_currentGraphToFit->GetX()[0];
      double ampToUse  = m_currentGraphToFit->GetY()[0];
      if(ampToUse > 0.0){
        double scale = ampToUse/fitFunction_knockout->Eval(xValToUse);
        fitFunction_knockout->FixParameter(4,scale);
        fitFunction_knockout->SetParLimits(0,0.001,1000.0);
      }else{
        fitFunction_knockout->SetParLimits(0,0.0,2.0);
        fitFunction_knockout->FixParameter(4,1E6);
      }
      fitFunction_knockout->SetParLimits(1,0,100.0);
      fitFunction_knockout->SetParLimits(2,0.0,5.0);
      fitFunction_knockout->SetParLimits(3,0.0,5.0);
      //fitFunction_knockout->SetParameter(5,0.01);
      //fitFunction_knockout->SetParLimits(5,0.0,20);
      //fitFunction_knockout->SetParLimits(6,0.0,15.0);

      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Seeding Knockout mTm0 Funct to : " << fitFunction_knockout->GetParameter(0) << " " << fitFunction_knockout->GetParameter(1) << " " << fitFunction_knockout->GetParameter(2)
             << " " << fitFunction_knockout->GetParameter(3) << " " << fitFunction_knockout->GetParameter(4) << endl;
      #endif

      TF1* initialFunct = (TF1*) fitFunction_knockout->Clone();
      initialFunct->SetName("initialFunct_knockout");
      initialFunct->SetLineColor(kYellow);

      int status_knockout = -999;
      string fitOptions = "S EX0";
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        fitOptions += " ";
      #else
        fitOptions += " Q";
      #endif
      TFitResultPtr resultPtrPtr = m_currentGraphToFit->Fit(fitFunction_knockout, "S EX0 V");
      TFitResult*   resultPtr    = resultPtrPtr.Get();
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Result Pointer: " << resultPtr << endl;
      #endif
      if(resultPtr){
        status_knockout = resultPtrPtr->Status();
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " Fit status: " << status_knockout << endl;
        #endif
        TF1* finalfitFunction_knockout = (TF1*) m_currentGraphToFit->GetFunction("fitFunction_knockout")->Clone();
        m_proton_knockout_ChiSqrd->SetBinContent(yIndex+1,centIndex,resultPtr->Chi2()/resultPtr->Ndf());
        m_proton_knockout_Fits_ByCent[centIndex][yIndex] = finalfitFunction_knockout;
        m_proton_knockout_CovMat_ByCent[centIndex][yIndex] =  (TMatrixDSym*) resultPtrPtr->GetCovarianceMatrix().Clone();

        TGraphAsymmErrors* pullGraph_rap = HistogramUtilities::makeFitPullGraph(m_currentGraphToFit, finalfitFunction_knockout);
        pullGraph_rap->SetMarkerStyle(20);
        pullGraph_rap->SetMarkerSize(0.5);
        pullGraph_rap->SetMarkerColor(kBlack);
        pullGraph_rap->SetLineColor(kBlack);

        TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
        fittingCanvas->SetWindowSize(1300,700);
        fittingCanvas->SetCanvasSize(1200,600);
        fittingCanvas->Draw();
        TPad* topPad_base = new TPad("topPad_Eff_base","topPad_Eff",0.0,0.3,1.0,1.0);
        TPad* bottomPad_base = new TPad("bottomPad_Eff_base","bottomPad_Eff",0.0,0.0,1.0,0.3);
        topPad_base->Draw();
        bottomPad_base->Draw();
        fittingCanvas->cd();

        double lowVertRange  = 0;
        double highVertRange = 0;
        HistogramUtilities::verticalDrawRange(m_currentGraphToFit, lowVertRange,highVertRange,-999,-999, true);

        fittingCanvas->cd();
        topPad_base->cd();
        gPad->SetLogy(true);
        TH1F *frame_knockout = gPad->DrawFrame(0.0,lowVertRange,0.15,highVertRange);
        frame_knockout->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Proton Knockout Fraction"));
        m_currentGraphToFit->Draw("P same");
        initialFunct->Draw("same");
        finalfitFunction_knockout->Draw("same");

        bottomPad_base->cd();
        TH1F *frame_pull = gPad->DrawFrame(0.0,-10.0,0.15,10.0);
        frame_pull->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Standardized Pulls"));
        pullGraph_rap->Draw("PSame");
        TLine* oneLineSys = new TLine(0.0,0.0,0.15,0.0);
        oneLineSys->SetLineStyle(7);
        oneLineSys->SetLineColor(kBlack);
        oneLineSys->Draw("same");
        gSystem->ProcessEvents();

        fittingCanvas->SaveAs(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/mTm0Fit_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(2,1).Data(), a_centIndex, yIndex));

        delete  pullGraph_rap;
        delete  oneLineSys;
        delete topPad_base;
        delete bottomPad_base;
        delete fittingCanvas;
      }else{
        cout << " very bad fit, no TFitResultPtr... status: " <<status_knockout<< endl;
        cout << "Check Graph Information | | | | | |" << endl;
        cout << "                        V V V V V V" << endl;
        HistogramUtilities::printGraphInfo(m_currentGraphToFit);
        m_proton_knockout_ChiSqrd->SetBinContent(yIndex+1,centIndex,999);

        TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
        fittingCanvas->SetWindowSize(1300,700);
        fittingCanvas->SetCanvasSize(1200,600);
        fittingCanvas->Draw();
        TPad* topPad_base = new TPad("topPad_Eff_base","topPad_Eff",0.0,0.3,1.0,1.0);
        TPad* bottomPad_base = new TPad("bottomPad_Eff_base","bottomPad_Eff",0.0,0.0,1.0,0.3);
        topPad_base->Draw();
        gPad->SetRightMargin(0.38);
        bottomPad_base->Draw();
        gPad->SetRightMargin(0.38);
        fittingCanvas->cd();

        fittingCanvas->cd();
        topPad_base->cd();
        gPad->SetLogy(false);
        TH1F *frame_knockout = fittingCanvas->DrawFrame(0.0,0,0.15,0.5);
        frame_knockout->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Proton Knockout Fraction"));
        m_currentGraphToFit->Draw("P same");

        bottomPad_base->cd();
        TH1F *frame_pull = fittingCanvas->DrawFrame(0.0,-10.0,0.15,10.0);
        frame_pull->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Standardized Pulls"));
        TLine* oneLineSys = new TLine(0.0,0.0,0.5,0.0);
        oneLineSys->SetLineStyle(7);
        oneLineSys->SetLineColor(kBlack);
        oneLineSys->Draw("same");
        gSystem->ProcessEvents();

        fittingCanvas->SaveAs(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/mTm0Fit_yIndex_%02d___BAD.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(2,1).Data(), a_centIndex, yIndex));

        delete  oneLineSys;
        delete topPad_base;
        delete bottomPad_base;
        delete fittingCanvas;
      }
      delete initialFunct;
    }//END yIndex LOOP




    // ###########################################################
    // FIT KNOCKOUT FRACTION AS A FUNCTION OF Rapidity
    // ###########################################################
      for(unsigned int mTm0Index = 0; mTm0Index < m_proton_dca_ChiSqrd[centIndex]->GetNbinsY(); mTm0Index++){
        double mTm0 = m_proton_dca_ChiSqrd[centIndex]->GetYaxis()->GetBinCenter(mTm0Index+1);
        double pT = PhysMath::pT_from_mTm0(mTm0,m_partInfo->GetParticleMass(2));
        if(pT > 0.38) continue;
        TH1D* protonHisto = protonTotalHisto2D->ProjectionX(Form("prot_cent%02d_mTm0%03d",centIndex,mTm0Index),mTm0Index+1,mTm0Index+1);
        TH1D* knockHisto = protonKnockoutHisto2D->ProjectionX(Form("knock_cent%02d_mTm0%03d",centIndex,mTm0Index),mTm0Index+1,mTm0Index+1);
        HistogramUtilities::zeroBinsBelowThreshold(protonHisto,0.0);
        HistogramUtilities::zeroBinsBelowThreshold(knockHisto,0.0);
        m_currentGraphToFit = new TGraphAsymmErrors();
        m_currentGraphToFit->SetName(Form("knockoutProtonFraction_Cent%02d_mTm0Index%02d",centIndex,mTm0Index));
        m_currentGraphToFit->SetTitle(Form("Knockout Proton Fraction Cent %02d m_{T}-m_{p}=[%1.3f,%1.3f]; y; N_{knockout}/N_{total}",
                           centIndex,m_proton_dca_ChiSqrd[centIndex]->GetYaxis()->GetBinLowEdge(mTm0Index+1),m_proton_dca_ChiSqrd[centIndex]->GetYaxis()->GetBinLowEdge(mTm0Index+2)));
        m_currentGraphToFit->Divide(knockHisto,protonHisto,"v cl=0.68268949 b(1,1) mode");
        m_currentGraphToFit->Sort();
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          HistogramUtilities::printHistoInfo(protonHisto);
          HistogramUtilities::printHistoInfo(knockHisto);
          HistogramUtilities::printGraphInfo(m_currentGraphToFit);
        #endif


        HistogramUtilities::discardNanData(m_currentGraphToFit,false);
        HistogramUtilities::removeBelowThreshold(m_currentGraphToFit,0.0);

        double status = 0;

        if(status == 0){

          TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
          fittingCanvas->SetWindowSize(1300,700);
          fittingCanvas->SetCanvasSize(1200,600);
          fittingCanvas->Draw();
          TPad* topPad_base = new TPad("topPad_Eff_base1","topPad_Eff",0.0,0.3,1.0,1.0);
          TPad* bottomPad_base = new TPad("bottomPad_Eff_base1","bottomPad_Eff",0.0,0.0,1.0,0.3);
          topPad_base->Draw();
          topPad_base->SetRightMargin(0.38);
          bottomPad_base->Draw();
          bottomPad_base->SetRightMargin(0.38);
          fittingCanvas->cd();


          topPad_base->cd();
          gPad->SetLogy();
          double lowVert = 0;
          double highVert = 0;
          HistogramUtilities::verticalDrawRange(m_currentGraphToFit, lowVert,highVert,-999,-999, true);
          double lowestRap  = protonTotalHisto2D->GetXaxis()->GetBinLowEdge(1);
          double highestRap = protonTotalHisto2D->GetXaxis()->GetBinLowEdge(protonTotalHisto2D->GetNbinsX() + 1);
          TH1F *frame_rap = fittingCanvas->DrawFrame(lowestRap,lowVert,highestRap,highVert);
          frame_rap->SetTitle(Form("; y; N_{knockout}/N_{total}"));
          m_currentGraphToFit->Draw("P");

          bottomPad_base->cd();
          TH1F *frame_pull_sys = fittingCanvas->DrawFrame(lowestRap,-10.0,highestRap,10.0);
          frame_pull_sys->SetTitle(Form("; y; Standardized Pulls"));
          TLine* oneLineSys = new TLine(lowestRap,0.0,highestRap,0.0);
          oneLineSys->SetLineStyle(7);
          oneLineSys->SetLineColor(kBlack);
          oneLineSys->Draw("same");
          gSystem->ProcessEvents();

          fittingCanvas->SaveAs(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/rapFit_mTm0Index_%02d.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(2,1).Data(), a_centIndex,mTm0Index));

          delete topPad_base;
          delete bottomPad_base;
          delete fittingCanvas;
          // status is good for fit
        }else{

          TCanvas* fittingCanvas = new TCanvas("fittingCanvas","Fitting Canvas");
          fittingCanvas->SetWindowSize(1300,700);
          fittingCanvas->SetCanvasSize(1200,600);
          fittingCanvas->Draw();
          TPad* topPad_base = new TPad("topPad_Eff_base1","topPad_Eff",0.0,0.3,1.0,1.0);
          TPad* bottomPad_base = new TPad("bottomPad_Eff_base1","bottomPad_Eff",0.0,0.0,1.0,0.3);
          topPad_base->Draw();
          topPad_base->SetRightMargin(0.38);
          bottomPad_base->Draw();
          bottomPad_base->SetRightMargin(0.38);
          fittingCanvas->cd();


          topPad_base->cd();
          gPad->SetLogy();
          double lowVert = 0;
          double highVert = 0;
          HistogramUtilities::verticalDrawRange(m_currentGraphToFit, lowVert,highVert,-999,-999, true);
          double lowestRap  = protonTotalHisto2D->GetXaxis()->GetBinLowEdge(1);
          double highestRap = protonTotalHisto2D->GetXaxis()->GetBinLowEdge(protonTotalHisto2D->GetNbinsX() + 1);
          TH1F *frame_rap = fittingCanvas->DrawFrame(lowestRap,lowVert,highestRap,highVert);
          frame_rap->SetTitle(Form("; y; N_{knockout}/N_{total"));
          m_currentGraphToFit->Draw("P");

          bottomPad_base->cd();
          TH1F *frame_pull_sys = fittingCanvas->DrawFrame(lowestRap,-10.0,highestRap,10.0);
          frame_pull_sys->SetTitle(Form("; y; Standardized Pulls"));
          TLine* oneLineSys = new TLine(lowestRap,0.0,highestRap,0.0);
          oneLineSys->SetLineStyle(7);
          oneLineSys->SetLineColor(kBlack);
          oneLineSys->Draw("same");
          gSystem->ProcessEvents();

          fittingCanvas->SaveAs(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/rapFit_mTm0Index_%02d____BAD.png",
            m_imageDir.c_str(), m_partInfo->GetParticleName(2,1).Data(), a_centIndex,mTm0Index));

          delete topPad_base;
          delete bottomPad_base;
          delete fittingCanvas;

          // status is bad for fit
        }
      }//mTm0Index loop








  }//centIndex
};


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


TF1* EfficiencyFitter::fitKnockoutTrendHistogram(int a_paramIndex,TH1D* a_histo){
  HistogramUtilities::discardNanHistoData(a_histo,true);
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    if(a_histo->GetBinError(binx) > a_histo->GetBinContent(binx)){
      a_histo->SetBinContent(binx,0.0);
      a_histo->SetBinError(binx,0.0);
    }
    if(a_paramIndex == 3 && a_histo->GetBinContent(binx) > 80){
      a_histo->SetBinContent(binx,0.0);
      a_histo->SetBinError(binx,0.0);
    }
    if(a_paramIndex == 2 && a_histo->GetBinContent(binx) > 2.8){
      a_histo->SetBinContent(binx,0.0);
      a_histo->SetBinError(binx,0.0);
    }
  }

  TF1* fitFunct = NULL;
  if(a_paramIndex == 0){ // p/pbar
    fitFunct = PhysMath::linearFunct(Form("par%d_y%d",a_paramIndex,m_currentRapIndex),0.2,-0.2,0.3);
  }else if(a_paramIndex == 1){ // A
    fitFunct = PhysMath::linearFunct(Form("par%d_y%d",a_paramIndex,m_currentRapIndex),0.08,-0.74,0.05);
    //fitFunct = PhysMath::getFitSafeExponentialFunction(Form("par%d_y%d",a_paramIndex,m_currentRapIndex),0.2,2.0);
  }else if(a_paramIndex == 2){ // DCA_0
    fitFunct = PhysMath::linearFunct(Form("par%d_y%d",a_paramIndex,m_currentRapIndex),0.4,2.66,0.05);
    //fitFunct = PhysMath::getFitSafeExponentialFunction(Form("par%d_y%d",a_paramIndex,m_currentRapIndex),0.2,2.0);
  }else if(a_paramIndex == 3){ // alpha
    fitFunct = PhysMath::linearFunct(Form("par%d_y%d",a_paramIndex,m_currentRapIndex),0.2,-0.2,0.3);
  }
  fitFunct->SetRange(0,0.5);
  TF1* initFunct = (TF1*) fitFunct->Clone();
  initFunct->SetName("initFunct");
  initFunct->SetLineColor(kYellow+1);

  int status = a_histo->Fit(fitFunct,"E");
  TF1* fitResult = a_histo->GetFunction(fitFunct->GetName());
  TF1* functToReturn = NULL;
  TH1D* pullGraph_par = NULL;
  if(fitResult){
     functToReturn = (TF1*) fitResult->Clone();
     pullGraph_par = HistogramUtilities::makeFitPullHistogram(a_histo, functToReturn);
     pullGraph_par->SetMarkerStyle(20);
     pullGraph_par->SetMarkerSize(0.5);
     pullGraph_par->SetMarkerColor(kBlack);
     pullGraph_par->SetLineColor(kBlack);
  }


  TCanvas* fittingCanvas_par = new TCanvas("fittingCanvas_par","Fitting Canvas");
  fittingCanvas_par->SetWindowSize(1300,700);
  fittingCanvas_par->SetCanvasSize(1200,600);
  fittingCanvas_par->Draw();
  TPad* topPad_par = new TPad("topPad_Eff_base","topPad_Eff",0.0,0.3,1.0,1.0);
  TPad* bottomPad_par = new TPad("bottomPad_Eff_base","bottomPad_Eff",0.0,0.0,1.0,0.3);
  topPad_par->Draw();
  gPad->SetRightMargin(0.38);
  bottomPad_par->Draw();
  gPad->SetRightMargin(0.38);
  fittingCanvas_par->cd();


  topPad_par->cd();
  a_histo->Draw("E");
  initFunct->Draw("same");
  if(functToReturn) functToReturn->Draw("same");
  double maxVal = HistogramUtilities::getMaxInRange(a_histo,0.0,0.7);
  a_histo->GetXaxis()->SetRangeUser(0.0,0.1);
  a_histo->GetYaxis()->SetRangeUser(0.0,1.2*maxVal);

  bottomPad_par->cd();
  TH1F *frame_pull_par = fittingCanvas_par->DrawFrame(0.0,-10.0,0.1,10.0);
  frame_pull_par->SetTitle(Form("; m_{T}-m_{0} (GeV/c^{2}); Standardized Pulls"));
  if(pullGraph_par) pullGraph_par->Draw("PSame");
  TLine* oneLinePar = new TLine(0.0,0.0,0.1,0.0);
  oneLinePar->SetLineStyle(7);
  oneLinePar->SetLineColor(kBlack);
  oneLinePar->Draw("same");

  gSystem->ProcessEvents();
  fittingCanvas_par->SaveAs(Form("%s/%s/KnockoutProton_ByCent/Cent%02d/param_%02d_yIndex_%02d.png",
  m_imageDir.c_str(), m_partInfo->GetParticleName(2,1).Data(), m_currentCentIndex, a_paramIndex, m_currentRapIndex));
  delete initFunct;
  delete oneLinePar;
  if(pullGraph_par) delete pullGraph_par;
  delete fittingCanvas_par;

  return functToReturn;
}

//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________



void EfficiencyFitter::fitEnergyLoss(int a_partIndex, int a_charge, int a_centIndex, vector<double> a_lowRangeExcludeRegionPt, vector<double> a_highRangeExcludeRegionPt, vector<double> a_lowRangeExcludeRegionEta, vector<double> a_highRangeExcludeRegionEta, bool a_fixFitFunction_y0_toZero){
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    cout << "EfficiencyFitter::fitEnergyLoss(int " << a_partIndex << ", int " << a_charge << ", int " << a_centIndex << ", vector<double> a_lowRangeExcludeRegionPt, vector<double> a_highRangeExcludeRegionPt, vector<double> a_lowRangeExcludeRegionEta, vector<double> a_highRangeExcludeRegionEta, bool " << a_fixFitFunction_y0_toZero << ")" << endl;
  #endif
  bool doSplineFitting = false;
  bool removeFirstZerosForLegendreFits = false;
  double threshold_for_remove_zeros = 0.005;
  double highTotalMtM0Range = 4.5;
  int polynomialOrder = 6;
  bool allCentMode = false;
  if(a_centIndex < 0) allCentMode = true;

  gSystem->mkdir(Form("%s/%s",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  gSystem->mkdir(Form("%s/%s/EnergyLoss_AllCent",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  gSystem->mkdir(Form("%s/%s/EnergyLoss_ByCent",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data()));
  gSystem->mkdir(Form("%s/%s/EnergyLoss_ByCent/Cent%02d",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex));

  if(a_partIndex == 0){
    cout << "WARNING: Pion Energy Loss does not need to be fit. I will not fit that curve. You cant make me..." << endl;
    return;
  }
  int pmIndex = 0;
  if(a_charge < 0) pmIndex = 1;
  TAxis* rapidity_bin_axis = m_embedding_rapidity_structure[a_partIndex]->GetXaxis();

  for(int yIndex = 0; yIndex < g_EffFit_MaxNumRapBins; yIndex++){
    double rapidity = rapidity_bin_axis->GetBinCenter(yIndex+1);
    if(m_rapiditiesToTest.size() > 0){
      bool rapInList = false;
      for(int iii = 0; iii < m_rapiditiesToTest.size(); iii++){
        if(fabs(m_rapiditiesToTest[iii] - rapidity) < 0.005){
          rapInList = true;
          break;
        }
      }
      if(!rapInList){
        cout << "Not calculating rapidity = " << rapidity << " because it isn't in the testing list... continuing..." << endl;
        continue;
      }
    }
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "Energy Loss for " << m_partInfo->GetParticleName(a_partIndex,a_charge) << "  cent " << a_centIndex << " yIndex " << yIndex << endl;
    #endif
    TProfile* mTm0LossProfile        = NULL;
    TH1D* mTm0LossHisto              = NULL;
    TGraphAsymmErrors* mTm0LossGraph = NULL;
    if(m_useEnergyLossTProfiles){
      if(allCentMode) mTm0LossProfile = m_EnergyLossProfile_AllCent[a_partIndex][yIndex][pmIndex];
      else            mTm0LossProfile = m_EnergyLossProfile_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
      mTm0LossHisto = dynamic_cast<TH1D*>(mTm0LossProfile);
      if(!mTm0LossProfile || HistogramUtilities::numBinsWithContentAboveThreshold(mTm0LossHisto,-1) < 5) continue;
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        HistogramUtilities::printHistoInfo(mTm0LossHisto);
      #endif
    }else{
      if(allCentMode) mTm0LossGraph = m_EnergyLoss_AllCent[a_partIndex][yIndex][pmIndex];
      else            mTm0LossGraph = m_EnergyLoss_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex];
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "energy loss graph ptr: " << mTm0LossGraph << endl;
      #endif
      if(!mTm0LossGraph){
        cout << "WARNING: Skipping this energy loss bin becasue the graph doesn't exist " << endl;
        continue;
      }
      if(HistogramUtilities::getNumNonZeroPts(mTm0LossGraph) < 5){
        cout << "WARNING: Skipping this energy loss bin becasue the graph has too few points" << endl;
        continue;
      }
    }

    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << "Fitting energy loss graph ptr: " << mTm0LossGraph << endl;
      cout << "   Name: " << mTm0LossGraph->GetName() << endl;
      cout << "   Title: " << mTm0LossGraph->GetTitle() << endl;
      mTm0LossGraph->Print();
    #endif

    if(a_lowRangeExcludeRegionPt.size() > 0 || a_lowRangeExcludeRegionEta.size() > 0){
      if(m_useEnergyLossTProfiles){
        for(int bin = 1; bin <= mTm0LossHisto->GetNbinsX(); bin++){
          double mTm0 = mTm0LossHisto->GetBinCenter(bin);
          double pT   = PhysMath::pT_from_mTm0(mTm0,m_partInfo->GetParticleMass(a_partIndex));
          if(m_pTLossMode){
            pT = mTm0LossHisto->GetBinCenter(bin);
            mTm0 = PhysMath::mTm0(m_partInfo->GetParticleMass(a_partIndex),pT);
          }

          for(unsigned int excludeIndex = 0; excludeIndex < a_lowRangeExcludeRegionPt.size(); excludeIndex++){
            if( pT > a_lowRangeExcludeRegionPt[excludeIndex] && pT < a_highRangeExcludeRegionPt[excludeIndex] ){ //removing point
              mTm0LossHisto->SetBinContent(bin,0.0);
              mTm0LossHisto->SetBinError(bin,0.0);
            }
          }
          double eta =  PhysMath::eta_from_y_mTm0(rapidity, mTm0, m_partInfo->GetParticleMass(a_partIndex));
          for(unsigned int excludeIndex = 0; excludeIndex < a_lowRangeExcludeRegionEta.size(); excludeIndex++){
            if(eta > a_lowRangeExcludeRegionEta[excludeIndex] && eta < a_highRangeExcludeRegionEta[excludeIndex]){
              mTm0LossHisto->SetBinContent(bin,0.0);
              mTm0LossHisto->SetBinError(bin,0.0);
            }//rap check
          }
        }//point loop
      }else{
        for(int pointIndex = 0; pointIndex < mTm0LossGraph->GetN(); pointIndex++){
          double mTm0 = mTm0LossGraph->GetX()[pointIndex];
          double pT = PhysMath::pT_from_mTm0(mTm0,m_partInfo->GetParticleMass(a_partIndex));
          if(m_pTLossMode){
            pT = mTm0LossGraph->GetX()[pointIndex];
            mTm0 = PhysMath::mTm0(m_partInfo->GetParticleMass(a_partIndex),pT);
          }

          for(unsigned int excludeIndex = 0; excludeIndex < a_lowRangeExcludeRegionPt.size(); excludeIndex++){
            if( pT > a_lowRangeExcludeRegionPt[excludeIndex] && pT < a_highRangeExcludeRegionPt[excludeIndex] ){ //removing point
              mTm0LossGraph->RemovePoint(pointIndex);
              pointIndex--;
              continue; //don't need to do eta if already cut
            }
          }
          double eta =  PhysMath::eta_from_y_mTm0(rapidity, mTm0, m_partInfo->GetParticleMass(a_partIndex));
          for(unsigned int excludeIndex = 0; excludeIndex < a_lowRangeExcludeRegionEta.size(); excludeIndex++){
            if(eta > a_lowRangeExcludeRegionEta[excludeIndex] && eta < a_highRangeExcludeRegionEta[excludeIndex]){
              mTm0LossGraph->RemovePoint(pointIndex);
              pointIndex--;
              continue;
            }//rap check
          }
        }//point loop
      }
    }


    /*
    for(unsigned int excludeIndex = 0; excludeIndex < a_lowRangeExcludeRegionPt.size(); excludeIndex++){
      if(m_useEnergyLossTProfiles){
        for(int bin = 1; bin <= mTm0LossHisto->GetNbinsX(); bin++){
          double mTm0 = mTm0LossHisto->GetBinCenter(bin);
          double pT = PhysMath::pT_from_mTm0(mTm0,m_partInfo->GetParticleMass(a_partIndex));
          if(m_pTLossMode){
            pT = mTm0LossHisto->GetBinCenter(bin);
            mTm0 = PhysMath::mTm0(m_partInfo->GetParticleMass(a_partIndex),pT);
          }
          double eta =  PhysMath::eta_from_y_mTm0(rapidity, mTm0, m_partInfo->GetParticleMass(a_partIndex));
          if(eta > a_lowRangeExcludeRegionEta[excludeIndex] && eta < a_highRangeExcludeRegionEta[excludeIndex]){
            if( pT > a_lowRangeExcludeRegionPt[excludeIndex] && pT < a_highRangeExcludeRegionPt[excludeIndex] ){ //removing point
              mTm0LossHisto->SetBinContent(bin,0.0);
              mTm0LossHisto->SetBinError(bin,0.0);
            }
          }//rap check
        }//point loop
      }else{
        for(int pointIndex = 0; pointIndex < mTm0LossGraph->GetN(); pointIndex++){
          double mTm0 = mTm0LossGraph->GetX()[pointIndex];
          double pT = PhysMath::pT_from_mTm0(mTm0,m_partInfo->GetParticleMass(a_partIndex));
          if(m_pTLossMode){
            pT = mTm0LossGraph->GetX()[pointIndex];
            mTm0 = PhysMath::mTm0(m_partInfo->GetParticleMass(a_partIndex),pT);
          }
          double eta =  PhysMath::eta_from_y_mTm0(rapidity, mTm0, m_partInfo->GetParticleMass(a_partIndex));
          if(eta > a_lowRangeExcludeRegionEta[excludeIndex] && eta < a_highRangeExcludeRegionEta[excludeIndex]){
            if( pT > a_lowRangeExcludeRegionPt[excludeIndex] && pT < a_highRangeExcludeRegionPt[excludeIndex] ){ //removing point
              mTm0LossGraph->RemovePoint(pointIndex);
              pointIndex--;
              continue;
            }
          }//rap check
        }//point loop
      }

    }// region loop*/


    double pTModSeeds[3][5] ={
     {0,0,0,1,0},
     {0.001241,3.484,0.0,1.5,0.2358},
     {0.006105,5.968,0.0,1.022,0.5155}
    };



    string fitOptions = "E S";
    if(!m_useEnergyLossTProfiles) fitOptions = fitOptions + " EX0";
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      fitOptions = fitOptions + " ";
    #else
      fitOptions = fitOptions + " Q";
    #endif
    bool   doPtModification = true;
    double highestMtM0InGraph = -1;
    if(m_useEnergyLossTProfiles) highestMtM0InGraph = mTm0LossHisto->GetBinLowEdge(HistogramUtilities::lastBinWithContent((TH1D*) mTm0LossHisto));
    else highestMtM0InGraph = HistogramUtilities::findLargestXVal(mTm0LossGraph);
    //####   Fit with simple function for systematic  ######
    //TF1* powerLawPt = new TF1("powerLawPt", "[0]*pow( (x+[3])*(x+[3]) + [3]*[3] ,-[1]) + [2]", 0,5.0); // a power law in pT space
    //This is the power law function, but propagated for mT-m0
    //TF1* powerLawPt = new TF1("powerLawPt", "sqrt(pow([0]*pow( (x+[3])*(x+[3]) + [3]*[3] ,-[1]) + [2] + sqrt( (x+[3])*(x+[3]) - [3]*[3]),2.0) + [3]*[3]) - x - [3]", 0,5.0); // a power law in pT space
    //sqrt(pow([0]*pow( (x+[3])*(x+[3]) + [3]*[3] ,-[1]) + [2] + sqrt( (x+[3])*(x+[3]) - [3]*[3]),2.0) + [3]*[3]) - x - [3]
    TF1* powerLawPt = NULL;
    if(!m_pTLossMode) powerLawPt = new TF1("powerLawPt", "sqrt(pow([0]*pow(x*x + 2.0*x*[3],-0.5*[1]) + [2] + sqrt(x*x + 2.0*x*[3]),2.0) + [3]*[3]) - x - [3]", 0,highTotalMtM0Range); // a power law in pT space
    else{
      if(doPtModification){
        powerLawPt = new TF1("powerLawPt", "[0]*pow([3]*x+[4],-[1]) + [2]", 0,highTotalMtM0Range); // a power law in pT space
      }else{
        powerLawPt = new TF1("powerLawPt", "[0]*pow(x,-[1]) + [2]", 0,highTotalMtM0Range); // a power law in pT space
      }
    }
    powerLawPt->SetNpx(250);
    //sqrt(pow([0]*pow(x*x + 2.0*x*[3],-0.5*[1]) + [2] + sqrt(x*x + 2.0*x*[3]),2.0) + [3]*[3]) - x - [3]

    powerLawPt->SetParameter(0,0.0005);
    powerLawPt->SetParLimits(0,0.0,1.0);
    powerLawPt->SetParameter(1,1.5);
    powerLawPt->SetParLimits(1,0,50.0);
    powerLawPt->SetParameter(2,0);
    powerLawPt->SetParLimits(2,-1.0,1.0);
    if(a_fixFitFunction_y0_toZero) powerLawPt->FixParameter(2,0.0);
    if(doPtModification){
      powerLawPt->SetParLimits(1,0.1,15.0);
      powerLawPt->SetParameter(3,1.0);
      powerLawPt->SetParLimits(3,0.001,25.0);
      powerLawPt->SetParameter(4,0.0);
      powerLawPt->SetParLimits(4,0.0,5.0);

      //do the seeding
      powerLawPt->SetParameter(0,1.0);
      powerLawPt->SetParameter(1,pTModSeeds[a_partIndex][1]);
      powerLawPt->SetParameter(3,pTModSeeds[a_partIndex][3]);
      powerLawPt->SetParameter(4,pTModSeeds[a_partIndex][4]);
      if(!m_useEnergyLossTProfiles){

        //use the largest value from [0,1] at the lowest xPoint

        double xValToUse = HistogramUtilities::findSmallestXVal(mTm0LossGraph);
        double yValToUse = HistogramUtilities::findLargestYVal(mTm0LossGraph,0,1);
        //double ampToUse = mTm0LossGraph->Eval(xValToUse);
        if(yValToUse > 0.0){
          powerLawPt->SetParameter(0,yValToUse/powerLawPt->Eval(xValToUse));
        }
      }
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << " Seeding Energy Loss to : " << powerLawPt->GetParameter(0) << " " << powerLawPt->GetParameter(1) << " " << powerLawPt->GetParameter(2) << " " << powerLawPt->GetParameter(3) << " " << powerLawPt->GetParameter(4) << endl;
      #endif
    }
    if(!m_pTLossMode){
      powerLawPt->FixParameter(3,m_partInfo->GetParticleMass(a_partIndex));
      powerLawPt->SetParNames("Scale","Exponent","y_0","m");
    }else{
      if(doPtModification){
        powerLawPt->SetParNames("Scale","Exponent","y_0","a","b");
      }else{
        powerLawPt->SetParNames("Scale","Exponent","y_0");
      }
    }

    TF1* initialFunct = (TF1*) powerLawPt->Clone();
    initialFunct->SetName("initialFunct");
    initialFunct->SetLineColor(kYellow+1);
    initialFunct->SetLineWidth(1);



    TFitResult*     resultPtr        = NULL;
    int             statusEnergyLoss = -999;
    double          highFitRange     = highestMtM0InGraph;
    TF1*            fitFunction_PowerLawPt = NULL;
    TMatrixDSym*    covMat           = NULL;
    TGraphErrors*   conf_graph       = NULL;
    if(m_useEnergyLossTProfiles){
      TFitResultPtr resultPtrPtr = mTm0LossHisto->Fit(powerLawPt, fitOptions.c_str(),"",0.0,highFitRange);
      resultPtr = resultPtrPtr.Get();
      if(resultPtr){
        fitFunction_PowerLawPt = mTm0LossHisto->GetFunction("powerLawPt");
      }
    }else{
      if(!m_fixed_target_mode && m_pTLossMode && doPtModification && highestMtM0InGraph > 1.2) highFitRange = 1.2;
      TFitResultPtr resultPtrPtr = mTm0LossGraph->Fit(powerLawPt, fitOptions.c_str(),"", 0.0, highFitRange);
      resultPtr = resultPtrPtr.Get();
      if(resultPtr){
        fitFunction_PowerLawPt = mTm0LossGraph->GetFunction("powerLawPt");
      }
    }
    if(resultPtr && fitFunction_PowerLawPt){
      statusEnergyLoss = resultPtr->Status();
      covMat = (TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();
      conf_graph = HistogramUtilities::makeErrorBandGraph(fitFunction_PowerLawPt, covMat, 300);
    }
    //int status = mTm0LossGraph->Fit("powerLawPt",fitOptions.c_str());
    if(!resultPtr){
      cout << " very bad fit, no TFitResultPtr... status: " << statusEnergyLoss<< endl;
      continue;
    }
    #ifdef _EFFICIENCY_FITTER_DEBUG_
      cout << " Fit status: " << statusEnergyLoss << endl;
    #endif
    if(!mTm0LossGraph){
      cout << "ERROR: Something went really wrong... no energy loss graph after fit: " << mTm0LossGraph << endl;
    }

    bool badFinalParameters     = false;

    if(!fitFunction_PowerLawPt){
      cout << "ERROR: Energy Loss fit for " << m_partInfo->GetParticleName(a_partIndex,a_charge).Data() << "  Cent: " << a_centIndex << "  yIndex:" << yIndex << "  DID NOT FIT PROPERLY!!!" << endl;
      cout << "      fitFunctionPointer: " << fitFunction_PowerLawPt << endl;
      //return;
    }else{

      /*if(m_pTLossMode && doPtModification){
        for(int iii = 0; iii < 3; iii++){ // b can go a little negative
          if(fitFunction_PowerLawPt->GetParameter(iii) - fitFunction_PowerLawPt->GetParError(iii) < 0.0){
            badFinalParameters = true;
          }
        }
        if(badFinalParameters){
          fitFunction_PowerLawPt->SetParameter(0,0.0);
          fitFunction_PowerLawPt->SetParameter(1,1.0);
          fitFunction_PowerLawPt->SetParameter(2,0.0);
          fitFunction_PowerLawPt->SetParameter(3,1.0);
          fitFunction_PowerLawPt->SetParameter(4,0.0);
        }
      }*/

      fitFunction_PowerLawPt->SetRange(0.0,highestMtM0InGraph);
      fitFunction_PowerLawPt->SetLineColor(kBlue);
      double chiPerNDF = resultPtr->Chi2() / resultPtr->Ndf();
      if(allCentMode){
        m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiPerNDF);
        m_EnergyLoss_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex] = (TF1*) fitFunction_PowerLawPt->Clone();
        m_EnergyLoss_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetName(Form("energyLossFit_%s_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_EnergyLoss_Fits_AllCent_Sys[a_partIndex][yIndex][pmIndex]->SetTitle(Form("Energy Loss %s %s Cent Integrated y=[%1.2f,%1.2f] Systematic", m_pTLossMode ? "pT" : "mT-m0", m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_EnergyLoss_CovMat_AllCent_Sys[a_partIndex][yIndex][pmIndex] = covMat; //(TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();

      }else{
        m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiPerNDF);
        m_EnergyLoss_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = (TF1*) fitFunction_PowerLawPt->Clone();
        m_EnergyLoss_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("energyLossFit_%s_Cent%d_yIndex%d_Sys",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_EnergyLoss_Fits_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("Energy Loss %s %s Cent %d y=[%1.2f,%1.2f] Systematic", m_pTLossMode ? "pT" : "mT-m0", m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_EnergyLoss_CovMat_ByCent_Sys[a_partIndex][a_centIndex][yIndex][pmIndex] = covMat; //(TMatrixDSym*) resultPtr->GetCovarianceMatrix().Clone();
      }
    }


    //####   Fit with Spline Fitter ########
    SplineFitter* fitter       = NULL;
    TF1*          fitFunction  = NULL;
    TGraph*       splinePoints = NULL;
    double        chiSqr       = 0.0;
    double        ndof         = 0.0;
    if(doSplineFitting){
      fitter = new SplineFitter();
      fitter->setGraphToFit(mTm0LossGraph);
      fitter->analyzeGraph(3,fitFunction, splinePoints,chiSqr,ndof,false,false,NULL,false,0.0);
      if(fitFunction) fitFunction->SetNpx(300);
      if(fitFunction) fitFunction->SetLineColor(kViolet);

      if(allCentMode){
        m_EnergyLoss_Fits_AllCent[a_partIndex][yIndex][pmIndex] = fitFunction;
        m_EnergyLoss_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("energyLossFit_%s_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
        m_EnergyLoss_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("Energy Loss %s Cent Integrated y=[%1.2f,%1.2f] Nominal",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_EnergyLoss_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = fitter->getCovMatrix();
        m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiSqr/ndof);

      }else{
        m_EnergyLoss_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = fitFunction;
        m_EnergyLoss_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("energyLossFit_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        m_EnergyLoss_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("Energy Loss %s Cent %d y=[%1.2f,%1.2f] Nominal",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                               rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        m_EnergyLoss_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = fitter->getCovMatrix();
        m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiSqr/ndof);

      }
    } // do Spline fitter

    //############ fit deviation with legendre //############//############//############

    if(fitFunction_PowerLawPt){
      //cloning mTm0LossGraph to make residules
      #ifdef _EFFICIENCY_FITTER_DEBUG_
        cout << "Cloning graph for residules: " << mTm0LossGraph << endl;
      #endif
      if(!mTm0LossGraph){
        cout << "ERROR: Something went really wrong... no energy loss graph: " << mTm0LossGraph << endl;
        cout << "Continuing, but there will be no correction for this rapidity: " << rapidity << endl;
        continue;
      }
      TGraphAsymmErrors* mTm0GraphResidulals = (TGraphAsymmErrors*) mTm0LossGraph->Clone();
      string residuleGraphName = mTm0LossGraph->GetName();
      residuleGraphName = residuleGraphName + "_residules";
      mTm0GraphResidulals->SetName(residuleGraphName.c_str());
      for(int pointIndex = 0; pointIndex < mTm0GraphResidulals->GetN(); pointIndex++){
        double X_ofPoint = mTm0GraphResidulals->GetX()[pointIndex];
        double Y_ofPoint = mTm0GraphResidulals->GetY()[pointIndex];
        double mainFunctEval = fitFunction_PowerLawPt->Eval(X_ofPoint);
        double residule = Y_ofPoint - mainFunctEval;
        //mTm0GraphResidulals->SetPointY(pointIndex,residule);
        mTm0GraphResidulals->SetPoint(pointIndex,X_ofPoint,residule);
        #ifdef _EFFICIENCY_FITTER_DEBUG_
          cout << " Residule for " << mTm0GraphResidulals->GetX()[pointIndex] << " is " <<  mTm0GraphResidulals->GetY()[pointIndex] << " with shift " << mainFunctEval << endl;
        #endif
      }
      double lowRange = 0.0;
      if(removeFirstZerosForLegendreFits){
        for(int pointIndex = 0; pointIndex < mTm0LossGraph->GetN(); pointIndex++){
          if(mTm0LossGraph->GetX()[pointIndex] < 1.0 && mTm0LossGraph->GetY()[pointIndex] < threshold_for_remove_zeros) lowRange = mTm0LossGraph->GetX()[pointIndex] + 0.5*mTm0LossGraph->GetEXhigh()[pointIndex];
        }
      }
      TF1* fitFunctionRes = PhysMath::getLegendreExpansion("legendreFunctMtM0LossRes", polynomialOrder,false,lowRange, highestMtM0InGraph);
      PhysMathHistogramUtilities::initializeLegendreExpansion(fitFunctionRes, mTm0GraphResidulals, lowRange, highestMtM0InGraph, false);
      TF1* initialLegendreSeedRes = (TF1*) fitFunctionRes->Clone();
      initialLegendreSeedRes->SetName("initialLegendreSeedRes");
      initialLegendreSeedRes->SetLineColor(kYellow);
      TFitResultPtr fitPtrRes = mTm0GraphResidulals->Fit(fitFunctionRes,"S EX0","",lowRange,highestMtM0InGraph);
      int legendreStatus = fitPtrRes->Status();
      //TMatrixDSym* legendreCovMatRes = (TMatrixDSym*) fitPtrRes->GetCovarianceMatrix().Clone();
      TMatrixDSym matrixTemp = PhysMath::addZerosToCovMatrix(fitPtrRes->GetCovarianceMatrix(), 4);
      TMatrixDSym* legendreCovMatRes = (TMatrixDSym*) matrixTemp.Clone();

      TF1* functOrig = mTm0GraphResidulals->GetFunction("legendreFunctMtM0LossRes");
      fitFunctionRes = (TF1*) functOrig->Clone();
      cout << " more issues? ptr: " << fitFunctionRes << endl;
      TCanvas* fittingCanvas = new TCanvas("fittingCanvas_EnergyLoss","Fitting Canvas",1300,900);
      fittingCanvas->cd();
      double lowYRange = 0;
      double highYRange = 0;
      HistogramUtilities::verticalDrawRange(mTm0GraphResidulals, lowYRange,highYRange,0.0,highestMtM0InGraph);
      TH1F *frame2 = fittingCanvas->DrawFrame(0.0,lowYRange,highestMtM0InGraph,highYRange);
      if(m_pTLossMode){
        frame2->SetTitle(Form("Energy Loss Residules %s Cent %d y=[%1.2f,%1.2f]; p_{T} (GeV/c); p_{T} Loss Residule",
          m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
          //m_cuts->getLowCentPercent(a_centIndex),
          //m_cuts->getHighCentPercent(a_centIndex),
          rapidity_bin_axis->GetBinLowEdge(yIndex+1),
          rapidity_bin_axis->GetBinLowEdge(yIndex+2))) ;
      }else{
        frame2->SetTitle(Form("Energy Loss Residules %s Cent %d y=[%1.2f,%1.2f]; m_{T}-m_{0} (GeV/c^{2}); MtM0 Loss Residule",
          m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
          //m_cuts->getLowCentPercent(a_centIndex),
          //m_cuts->getHighCentPercent(a_centIndex),
          rapidity_bin_axis->GetBinLowEdge(yIndex+1),
          rapidity_bin_axis->GetBinLowEdge(yIndex+2))) ;
      }

      mTm0GraphResidulals->Draw("PSame");
      if(initialLegendreSeedRes) initialLegendreSeedRes->Draw("same");
      if(fitFunctionRes) fitFunctionRes->Draw("same");

      gPad->Update();
      gSystem->ProcessEvents();

      if(allCentMode){
        fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_AllCent/Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), yIndex));
      }else{
        fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_ByCent/Cent%02d/Residules_yIndex_%02d.png",
          m_imageDir.c_str(), m_partInfo->GetParticleName(a_partIndex,a_charge).Data(), a_centIndex, yIndex));
      }


      delete frame2;
      delete fittingCanvas;



      if(fitFunctionRes){

        string totalFormula = PhysMath::getLegendreExpansionString(polynomialOrder);
        int numParams = fitFunctionRes->GetNpar();
        //sqrt(pow([0]*pow(x*x + 2.0*x*[3],-0.5*[1]) + [2] + sqrt(x*x + 2.0*x*[3]),2.0) + [3]*[3]) - x - [3]
        string additionalFormula = "";
        if(!m_pTLossMode){
          additionalFormula = Form(" + sqrt(pow([%d]*pow(x*x + 2.0*x*[%d],-0.5*[%d]) + [%d] + sqrt(x*x + 2.0*x*[%d]),2.0) + [%d]*[%d]) - x - [%d]",numParams,numParams+3,numParams+1,numParams+2,numParams+3,numParams+3,numParams+3,numParams+3);
        }else{
          if(doPtModification){
            additionalFormula = Form(" + [%d]*pow([%d]*x + [%d],-[%d]) + [%d]",numParams,numParams+3,numParams+4,numParams+1,numParams+2);
          }else{
            additionalFormula = Form(" + [%d]*pow(x,-[%d]) + [%d]",numParams,numParams+1,numParams+2);
          }
        }
        totalFormula = totalFormula + additionalFormula;

        cout << "combined formula: " << totalFormula << endl;

        TF1* combinedFunct = new TF1("combinedFunct",totalFormula.c_str(),lowRange,highestMtM0InGraph);
        for(int paramIndex = 0; paramIndex < numParams; paramIndex++){
          combinedFunct->SetParameter(paramIndex,fitFunctionRes->GetParameter(paramIndex));
        }
        combinedFunct->SetParameter(numParams,  fitFunction_PowerLawPt->GetParameter(0));
        combinedFunct->SetParameter(numParams+1,fitFunction_PowerLawPt->GetParameter(1));
        combinedFunct->SetParameter(numParams+2,fitFunction_PowerLawPt->GetParameter(2));
        if(!m_pTLossMode) combinedFunct->FixParameter(numParams+3,fitFunction_PowerLawPt->GetParameter(3));
        else if(doPtModification){
          combinedFunct->FixParameter(numParams+3,fitFunction_PowerLawPt->GetParameter(3));
          combinedFunct->FixParameter(numParams+4,fitFunction_PowerLawPt->GetParameter(4));
        }
        combinedFunct->SetLineColor(kRed);
        //combinedFunct->SetRange(lowRange,highestMtM0InGraph);

        if(allCentMode){
          m_EnergyLoss_Fits_AllCent[a_partIndex][yIndex][pmIndex] = combinedFunct;
          m_EnergyLoss_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetName(Form("energyLossFit_%s_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
          m_EnergyLoss_Fits_AllCent[a_partIndex][yIndex][pmIndex]->SetTitle(Form("Energy Loss %s Cent Integrated y=[%1.2f,%1.2f] Nominal",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          if(doSplineFitting) m_EnergyLoss_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = fitter->getCovMatrix();
          else m_EnergyLoss_CovMat_AllCent[a_partIndex][yIndex][pmIndex] = legendreCovMatRes;
          //m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,17,chiSqr/ndof);

        }else{
          m_EnergyLoss_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = combinedFunct;
          m_EnergyLoss_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetName(Form("energyLossFit_%s_Cent%d_yIndex%d",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
          m_EnergyLoss_Fits_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex]->SetTitle(Form("Energy Loss %s Cent %d y=[%1.2f,%1.2f] Nominal",m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,
                                 rapidity_bin_axis->GetBinLowEdge(yIndex+1),rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
          m_EnergyLoss_CovMat_ByCent[a_partIndex][a_centIndex][yIndex][pmIndex] = fitter->getCovMatrix();
          //m_EnergyLoss_ChiSqrd[a_partIndex][pmIndex]->SetBinContent(yIndex+1,a_centIndex+1,chiSqr/ndof);

        }

        TGraphAsymmErrors* pullGraph = HistogramUtilities::makeFitPullGraph(mTm0LossGraph, combinedFunct, false);
        pullGraph->SetMarkerStyle(20);
        pullGraph->SetMarkerSize(0.5);
        pullGraph->SetMarkerColor(kBlack);
        pullGraph->SetLineColor(kBlack);


        TCanvas* fittingCanvas = new TCanvas("fittingCanvas_EnergyLoss","Fitting Canvas",1300,900);
        fittingCanvas->Draw();
        TPad* topPad = new TPad("topPad","topPad",0.0,0.3,1.0,1.0);
        TPad* bottomPad = new TPad("bottomPad","bottomPad",0.0,0.0,1.0,0.3);
        topPad->Draw();
        bottomPad->Draw();


        //gPad->SetRightMargin(0.45);
        topPad->cd();
        double lowYRange = 0;
        double highYRange = 0;
        HistogramUtilities::verticalDrawRange(mTm0LossGraph, lowYRange,highYRange);
        TH1F* frame = topPad->DrawFrame(0.0,lowYRange,highestMtM0InGraph+0.01,highYRange);
        frame->SetTitle(Form("%s EnergyLoss y=[%.2f,%.2f] 0-80%% Central",m_partInfo->GetParticleSymbol(a_partIndex,a_charge).Data(),
                    rapidity_bin_axis->GetBinLowEdge(yIndex+1), rapidity_bin_axis->GetBinLowEdge(yIndex+2)));
        if(!m_pTLossMode){
          frame->GetYaxis()->SetTitle("m_{T}^{emb}-m_{T}^{reco} (GeV/c^{2})");
          frame->GetXaxis()->SetTitle("m_{T}-m_{0} (GeV/c^{2})");
        }else{
          frame->GetYaxis()->SetTitle("p_{T}^{emb}-p_{T}^{reco} (GeV/c)");
          frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        }
        mTm0LossGraph->Draw("Psame");
        if(initialFunct) initialFunct->Draw("same");
        if(conf_graph) conf_graph->Draw("3");
        fitFunction_PowerLawPt->Draw("same");
        combinedFunct->Draw("same");

        TLine* highLine = new TLine(highFitRange,lowYRange,highFitRange,highYRange);
        highLine->SetLineColor(kViolet);
        highLine->Draw("SAME");

        //combinedFunct->Draw();
        if(splinePoints) splinePoints->Draw("Psame");
        gPad->Update();
        gSystem->ProcessEvents();

        bottomPad->cd();
        TH1F* frame_3 = bottomPad->DrawFrame(0.0,lowYRange,highestMtM0InGraph+0.01,highYRange);
        frame_3->SetTitle("");
        frame_3->GetYaxis()->SetTitle("Standardized Pulls");
        if(!m_pTLossMode) frame_3->GetXaxis()->SetTitle("m_{T}-m_{0} (GeV/c^{2})");
        else frame_3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        frame_3->GetYaxis()->SetRangeUser(-4,4);
        //gStyle->SetOptStat(0);
        //gStyle->SetOptTitle(0);
        bottomPad->SetLogy(false);
        pullGraph->Draw("PSame");
        TLine* oneLine = new TLine(0.0,0.0,highestMtM0InGraph,0.0);
        oneLine->SetLineStyle(7);
        oneLine->SetLineColor(kBlack);
        oneLine->Draw("same");

        fittingCanvas->Update();
        gSystem->ProcessEvents();
        if(allCentMode){
          if(!m_pTLossMode) fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_AllCent/mTLoss_Fit_yIndex%02d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));
          else             fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_AllCent/pTLoss_Fit_yIndex%02d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),yIndex));

        }else{
          if(!m_pTLossMode) fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_ByCent/Cent%02d/mTLoss_Fit_yIndex%02d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
          else              fittingCanvas->SaveAs(Form("%s/%s/EnergyLoss_ByCent/Cent%02d/pTLoss_Fit_yIndex%02d.png",m_imageDir.c_str(),m_partInfo->GetParticleName(a_partIndex,a_charge).Data(),a_centIndex,yIndex));
        }

        delete oneLine;
        delete highLine;

        if(doSplineFitting){
          m_splinesToDelete.push_back(fitter);
          //delete fitter;
          delete splinePoints;
        }
        delete fittingCanvas;
      }
    }//end sys funct exist
    if(initialFunct) delete initialFunct;
  } // end rap loop
}// end fitEnergyLoss



//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


/*
void EfficiencyFitter::fitTPCEff_Chris(int a_partIndex){

  // *****************************************
  //   TPC TRACKING EFFICIENCY
  // *****************************************

  //----STEP ONE---
  //First Fit the centrality integrated graphs to extract the turn on
  //parameters of the error function. These parameters will be fixed
  //when we fit the efficiency graphs in each centrality bin
  int nRapidityBins = m_cuts->getNRapidityBins(a_partIndex,1);
  double prevPars[3] = {.7,0.005,2};
  for (int yIndex=0; yIndex < nRapidityBins; yIndex++){
    //Compute Kinematics
    double rapidity = m_cuts->getRapidityRangeCenter(yIndex,a_partIndex,1);
    int yBinIndex = m_embTrackHisto2D_AllCent[a_partIndex]->GetXaxis()->FindBin(rapidity);

    //Get the 1D Histograms
    TH1D *embTrackHisto = m_embTrackHisto2D_AllCent[a_partIndex]->ProjectionY("emby",yBinIndex,yBinIndex);
    TH1D *matchTrackHisto = m_matchTrackHisto2D_AllCent[a_partIndex]->ProjectionY("matchy",yBinIndex,yBinIndex);

    if (embTrackHisto->GetEntries() < 10 || matchTrackHisto->GetEntries() < 10){
      delete embTrackHisto;
      delete matchTrackHisto;
      continue;
    }

    //Create the Efficiency Graph
    corrFile->cd();
    TGraphAsymmErrors *efficiencyGraph = new TGraphAsymmErrors();
    efficiencyGraph->BayesDivide(matchTrackHisto,embTrackHisto);

    if (efficiencyGraph->GetN() < 10){
      delete efficiencyGraph;
      continue;
    }

    //Create the Fit
    m__AllCentEfficiencyFits[partIndex][yIndex] = new TF1(Form("tpcEfficiencyFit_%s_AllCent_yIndex%d",
      particleInfo->GetParticleName(a_partIndex,charge).Data(),
      yIndex),
    "[0]*exp(-([1]/pow(x,[2])))",minFitRange,3.0);
    m__AllCentEfficiencyFits[partIndex][yIndex]->SetNpx(10000);
    m__AllCentEfficiencyFits[partIndex][yIndex]->SetLineWidth(3);
    m__AllCentEfficiencyFits[partIndex][yIndex]->SetLineColor(particleInfo->GetParticleColor(a_partIndex));
    m__AllCentEfficiencyFits[partIndex][yIndex]->SetLineStyle(7);

    m__AllCentEfficiencyFits[partIndex][yIndex]->SetParameters(prevPars);
    m__AllCentEfficiencyFits[partIndex][yIndex]->SetParLimits(0,0.005,0.9);
    m__AllCentEfficiencyFits[partIndex][yIndex]->SetParLimits(1,0.000001,.07);
    m__AllCentEfficiencyFits[partIndex][yIndex]->SetParLimits(2,.1,25);


    TFitResultPtr fitResult = efficiencyGraph->Fit(m__AllCentEfficiencyFits[partIndex][yIndex],"RNSEX0");
    int attempt = 0;
    while (fitResult.Get()->Status() != 0 && attempt < 10) {
      fitResult = efficiencyGraph->Fit(m__AllCentEfficiencyFits[partIndex][yIndex],"RNSEX0");
      cout << attempt++ <<endl;;
    }

    prevPars[0] = m__AllCentEfficiencyFits[partIndex][yIndex]->GetParameter(0);
    prevPars[1] = m__AllCentEfficiencyFits[partIndex][yIndex]->GetParameter(1);
    prevPars[2] = m__AllCentEfficiencyFits[partIndex][yIndex]->GetParameter(2);

    //Clean Up
    delete efficiencyGraph;

  }//End Loop Over yIndex to fit centrality integrated efficiencies
  cout << "Done with step 1" << endl;


  //----STEP TWO---
  //Loop Over the centrality bins, for each cent bin construct the 1D efficiency and
  //energy loss graphs for various rapidity bins. Use the fits to the centrality integrated
  //graphs from STEP ONE to fix the turn on parameters.
  for (int centIndex=0; centIndex<nCentBins; centIndex++){

    //Loop Over the Rapidity Bins
    for (int yIndex=0; yIndex < nRapidityBins; yIndex++){

      //Compute Kinematics
      double rapidity = m_cuts->getRapidityRangeCenter(yIndex,a_partIndex,1);
      int yBinIndex = embTrackHisto2D.at(centIndex)->GetXaxis()->FindBin(rapidity);


      //Get the 1D Histograms
      TH1D *embTrackHisto = embTrackHisto2D.at(centIndex)->ProjectionY("emby",yBinIndex,yBinIndex);
      TH1D *matchTrackHisto = matchTrackHisto2D.at(centIndex)->ProjectionY("matchy",yBinIndex,yBinIndex);

      if (embTrackHisto->GetEntries() < 10 || matchTrackHisto->GetEntries() < 10)
       continue;

      //Rebin
     embTrackHisto->Rebin(rebinVal);
     matchTrackHisto->Rebin(rebinVal);

      //Create the Efficiency Graph
     corrFile->cd();

     TGraphAsymmErrors *efficiencyGraph = new TGraphAsymmErrors();
     efficiencyGraph->BayesDivide(matchTrackHisto,embTrackHisto);

     if (efficiencyGraph->GetN() < 10)
       continue;

     efficiencyGraph->SetMarkerStyle(particleInfo->GetParticleMarker(a_partIndex,charge));
     efficiencyGraph->SetMarkerColor(particleInfo->GetParticleColor(a_partIndex));
     efficiencyGraph->SetName(Form("tpcEfficiencyGraph_%s_Cent%d_yIndex%d",
      particleInfo->GetParticleName(a_partIndex,charge).Data(),
      centIndex,yIndex));
     efficiencyGraph->SetTitle(Form("TPC Tracking Efficiency %s | Cent=[%.02g,%.02g]%% | y_{%s}^{Emb}=[%.03g,%.03g];(m_{T}-m_{%s})^{Emb};Efficiency #times Acceptance",
      particleInfo->GetParticleSymbol(a_partIndex,charge).Data(),
				     //				     centIndex == nCentBins-1 ? 0:centralityPercents.at(centIndex+1),
      m_cuts->getLowCentPercent(centIndex),
      m_cuts->getHighCentPercent(centIndex),
      particleInfo->GetParticleSymbol(a_partIndex).Data(),
      m_cuts->getRapidityRangeLow(yIndex,a_partIndex,1),
      m_cuts->getRapidityRangeHigh(yIndex,a_partIndex,1),
      particleInfo->GetParticleSymbol(a_partIndex).Data()));
     efficiencyGraph->SetMarkerStyle(particleInfo->GetParticleMarker(a_partIndex,charge));
     efficiencyGraph->SetMarkerColor(particleInfo->GetParticleColor(a_partIndex));
     efficiencyGraph->SetMarkerSize(a_partIndex == KAON ? 1.5 : 1.2);

     cout << Form("TPC Tracking Efficiency %s | Cent=[%.02g,%.02g]%% | y_{%s}^{Emb}=[%.03g,%.03g];(m_{T}-m_{%s})^{Emb};Efficiency #times Acceptance",
       particleInfo->GetParticleSymbol(a_partIndex,charge).Data(),
       m_cuts->getLowCentPercent(centIndex),
       m_cuts->getHighCentPercent(centIndex),
       particleInfo->GetParticleSymbol(a_partIndex).Data(),
       m_cuts->getRapidityRangeLow(yIndex,a_partIndex,1),
       m_cuts->getRapidityRangeHigh(yIndex,a_partIndex,1),
       particleInfo->GetParticleSymbol(a_partIndex).Data()) << endl;

     cout << "centIndex: " << centIndex << endl;
     cout << "Low percent: " << m_cuts->getLowCentPercent(centIndex) << endl;
     cout << "High percent: " << m_cuts->getHighCentPercent(centIndex) << endl;


      //Create the Fit
     TF1 efficiencyFit(Form("tpcEfficiencyFit_%s_Cent%d_yIndex%d",
      particleInfo->GetParticleName(a_partIndex,charge).Data(),
      centIndex,yIndex),
     "[0]*exp(-([1]/pow(x,[2])))",minFitRange,3.0);
     efficiencyFit.SetNpx(10000);
     efficiencyFit.SetLineWidth(3);
     efficiencyFit.SetLineColor(particleInfo->GetParticleColor(a_partIndex));
     efficiencyFit.SetLineStyle(1);

     efficiencyFit.SetParameter(0,m__AllCentEfficiencyFits.at(yIndex)->GetParameter(0));
     efficiencyFit.SetParameter(1,m__AllCentEfficiencyFits.at(yIndex)->GetParameter(1));
     efficiencyFit.SetParameter(2,m__AllCentEfficiencyFits.at(yIndex)->GetParameter(2));

     efficiencyFit.SetParLimits(0,.2,1);
     efficiencyFit.SetParLimits(1,0.00001,.25);
     efficiencyFit.SetParLimits(2,1.0,6.0);

      //Create the Fit which will be used for the systematic Errors
      //This function is the same form as the nominal fit above,
      //but will not have its parameters fixed.
     TF1 efficiencyFitSys(Form("%s_Sys",efficiencyFit.GetName()),
      Form("%s",efficiencyFit.GetExpFormula().Data()),
      minFitRange,3.0);
     efficiencyFitSys.SetNpx(10000);
     efficiencyFitSys.SetLineWidth(3);
     efficiencyFitSys.SetLineColor(kGray);
     efficiencyFitSys.SetLineStyle(7);

     efficiencyFitSys.SetParameter(0,efficiencyFit.GetParameter(0));
     efficiencyFitSys.FixParameter(1,efficiencyFit.GetParameter(1));
     efficiencyFitSys.SetParameter(2,efficiencyFit.GetParameter(2));

     efficiencyFitSys.SetParLimits(0,efficiencyFitSys.GetParameter(0)*.5,1.0);
     efficiencyFit.SetParLimits(1,0.00001,.25);
      //efficiencyFitSys.SetParLimits(1,TMath::Min(0.000001,efficiencyFitSys.GetParameter(1)),fabs(efficiencyFitSys.GetParameter(1))*5.0);
     efficiencyFitSys.SetParLimits(2,1.0,6.0);


     TFitResultPtr fitResult = efficiencyGraph->Fit(&efficiencyFit,"RS");
     int attempt(0);
     while (fitResult.Get()->Status() != 0 && attempt < 10 ||
      (efficiencyFit.GetParameter(1) - efficiencyFit.GetParError(1)) < 0) {
       cout <<"*********REFITTING!*********" <<endl;

     if (efficiencyFit.GetParameter(1) - efficiencyFit.GetParError(1) < 0){
       efficiencyFit.FixParameter(1,efficiencyFit.GetParameter(1));
     }

     efficiencyFit.SetParameter(0,m__AllCentEfficiencyFits.at(yIndex)->GetParameter(0));
     efficiencyFit.SetParameter(2,m__AllCentEfficiencyFits.at(yIndex)->GetParameter(2));

     fitResult = efficiencyGraph->Fit(&efficiencyFit,"RS");
     attempt++;

     cout <<attempt <<" " <<fitResult.Get()->Status() <<endl;
   }

   prevPars[0] = efficiencyFit.GetParameter(0);
   prevPars[1] = efficiencyFit.GetParameter(1);
   prevPars[2] = efficiencyFit.GetParameter(2);

   efficiencyFit.SetParError(1,m__AllCentEfficiencyFits.at(yIndex)->GetParError(1));
   efficiencyFit.SetParError(2,m__AllCentEfficiencyFits.at(yIndex)->GetParError(2));

   TGraphErrors *efficiencyFitConf = GetConfidenceIntervalOfFit(&efficiencyFit);

    //Fit the Systematic Function
   efficiencyGraph->Fit(&efficiencyFitSys,"REX0+");
   TGraphErrors *efficiencyFitSysConf = GetConfidenceIntervalOfFit(&efficiencyFitSys);


   //Save
   corrFile->cd();
   corrFile->cd(Form("%s",particleInfo->GetParticleName(a_partIndex,charge).Data()));
   gDirectory->cd("EfficiencyGraphs");
   efficiencyGraph->Write(efficiencyGraph->GetName(),TObject::kOverwrite);
   corrFile->cd();
   corrFile->cd(Form("%s",particleInfo->GetParticleName(a_partIndex,charge).Data()));
   gDirectory->cd("EfficiencyFits");
   efficiencyFit->Write(efficiencyFit.GetName(),TObject::kOverwrite);
   efficiencyFitConf->Write(efficiencyFitConf->GetName(),TObject::kOverwrite);
   efficiencyFitSys->Write(efficiencyFitSys.GetName(),TObject::kOverwrite);
   efficiencyFitSysConf->Write(efficiencyFitSysConf->GetName(),TObject::kOverwrite);
   corrFile->cd();


      //Clean Up
   if (efficiencyGraph)
     delete efficiencyGraph;
   if (efficiencyFitConf)
     delete efficiencyFitConf;
   if (efficiencyFitSysConf)
     delete efficiencyFitSysConf;


    }//End Loop Over yIndex

    corrFile->cd();
    corrFile->cd(Form("%s",particleInfo->GetParticleName(a_partIndex,charge).Data()));
    gDirectory->cd("EfficiencyGraphs2D");
    effHisto2D.at(centIndex)->Write(Form("tpcEfficiency2D_%s_Cent%d",particleInfo->GetParticleName(a_partIndex,charge).Data(),centIndex));
    corrFile->cd();
  }//End Loop Over Centrality Bin
}
*/



//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double EfficiencyFitter::TGraphAsymmErrors_ChiSqr(const double* a_param){
  for(int paramIndex = 0; paramIndex < m_currentFitFunction->GetNpar(); paramIndex++){
    m_currentFitFunction->SetParameter(paramIndex,a_param[paramIndex]);
  }
  double chiSquared = 0.0;
  for(int pointIndex = 0; pointIndex < m_currentGraphToFit->GetN(); pointIndex++){
    double lowX  = m_currentGraphToFit->GetX()[pointIndex] - m_currentGraphToFit->GetEXlow()[pointIndex];
    double highX = m_currentGraphToFit->GetX()[pointIndex] + m_currentGraphToFit->GetEXhigh()[pointIndex];
    double integral = m_currentFitFunction->Integral(lowX,highX)/(highX - lowX);
    if(!std::isfinite(integral)) integral = m_currentFitFunction->Eval(m_currentGraphToFit->GetX()[pointIndex]);
    double value = m_currentGraphToFit->GetY()[pointIndex];
    double error = m_currentGraphToFit->GetEYhigh()[pointIndex];
    if(value > integral) error = m_currentGraphToFit->GetEYlow()[pointIndex];
    double addition = (value - integral)*(value - integral)/(error*error);
    if(std::isfinite(addition)) chiSquared += addition;
  }
  return chiSquared;
}

//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________



double EfficiencyFitter::spline_operator_3(double* a_x, double* a_params){
  bool sameParams = m_spline;
  for(int iii = 0; sameParams && iii < (2*m_spline_numKnots)+2; iii++){
    if(m_splineParams[iii] != a_params[iii]){
      sameParams = false;
    }
  }
  if(!sameParams){
    if(m_spline) delete m_spline;
    m_spline = nullptr;

    vector<double> x_loc(m_spline_numKnots,0);
    vector<double> y_loc(m_spline_numKnots,0);

    bool sorted = true;
    x_loc[0] = a_params[0];
    y_loc[0] = a_params[m_spline_numKnots];
    for(int knotIndex = 1; knotIndex < m_spline_numKnots; knotIndex++){
      x_loc[knotIndex] = a_params[knotIndex];
      //x_loc[knotIndex] = x_loc[knotIndex-1] + a_params[knotIndex];
      if(knotIndex > 1 && x_loc[knotIndex] < x_loc[knotIndex-1]) sorted = false;
      y_loc[knotIndex] = a_params[m_spline_numKnots + knotIndex];
    }

    while(!sorted){
      sorted = true;
      for(int knotIndex = 1; knotIndex < m_spline_numKnots; knotIndex++){
        if(x_loc[knotIndex] < x_loc[knotIndex-1]){
          sorted = false;
          double x_temp = x_loc[knotIndex];
          double y_temp = y_loc[knotIndex];
          x_loc[knotIndex] = x_loc[knotIndex-1];
          y_loc[knotIndex] = y_loc[knotIndex-1];
          x_loc[knotIndex-1] = x_temp;
          y_loc[knotIndex-1] = y_temp;
        }
      }
    }


    m_spline = new TSpline3("spline", &x_loc[0], &y_loc[0], m_spline_numKnots, "b1e1", a_params[2*m_spline_numKnots], a_params[2*m_spline_numKnots+1]);
    for(int iii = 0; iii < 2*m_spline_numKnots+2; iii++){
      m_splineParams[iii] = a_params[iii];
    }
  }

  return m_spline->Eval(a_x[0]);


  /*double knot_x;
  double knot_y;
  for(unsigned int pointIndex = 0; pointIndex < m_spline->GetNp(); pointIndex++){
    m_spline->GetKnot(pointIndex, knot_x, knot_y);
    m_spline->SetPoint(pointIndex,knot_x, a_params[pointIndex]);
  }
  double b1,c1,d1,b2,c2,d2;
  m_spline->GetCoeff(0,knot_x,knot_y,b1,c1,d1);
  m_spline->GetCoeff(m_spline->GetNp()-1, knot_x, knot_y,b2,c2,d2);
  m_spline->SetPointCoeff(0,a_params[m_spline->GetNp()],c1,d1);
  m_spline->SetPointCoeff(m_spline->GetNp()-1,a_params[m_spline->GetNp()+1],c2,d2);
  //m_spline->BuildCoeff();
  #ifdef _EFFICIENCY_FITTER_DEBUG_
    m_spline->GetCoeff(0,knot_x,knot_y,b1,c1,d1);
    cout << "Equal? " << a_params[m_spline->GetNp()] << " and " << b1 << endl;
    if(  a_params[m_spline->GetNp()] != b1){
      cout << "-------------------------------------------------------------------------    NOT EQUAL!!!!!!!!!!!!" << endl;
    }
  #endif

  return m_spline->Eval(a_x[0]);*/
}


//_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


double EfficiencyFitter::spline_operator_5(double* a_x, double* a_params){
  bool sameParams = m_spline5;
  for(int iii = 0; sameParams && iii < (2*m_spline_numKnots)+4; iii++){
    if(m_splineParams[iii] != a_params[iii]){
      sameParams = false;
    }
  }
  if(!sameParams){
    if(m_spline5) delete m_spline5;
    m_spline5 = nullptr;

    vector<double> x_loc(m_spline_numKnots,0);
    vector<double> y_loc(m_spline_numKnots,0);

    bool sorted = true;
    x_loc[0] = a_params[0];
    y_loc[0] = a_params[m_spline_numKnots];
    for(int knotIndex = 1; knotIndex < m_spline_numKnots; knotIndex++){
      x_loc[knotIndex] = a_params[knotIndex];
      //x_loc[knotIndex] = x_loc[knotIndex-1] + a_params[knotIndex];
      if(knotIndex > 1 && x_loc[knotIndex] < x_loc[knotIndex-1]) sorted = false;
      y_loc[knotIndex] = a_params[m_spline_numKnots + knotIndex];
    }

    while(!sorted){
      sorted = true;
      for(int knotIndex = 1; knotIndex < m_spline_numKnots; knotIndex++){
        if(x_loc[knotIndex] < x_loc[knotIndex-1]){
          sorted = false;
          double x_temp = x_loc[knotIndex];
          double y_temp = y_loc[knotIndex];
          x_loc[knotIndex] = x_loc[knotIndex-1];
          y_loc[knotIndex] = y_loc[knotIndex-1];
          x_loc[knotIndex-1] = x_temp;
          y_loc[knotIndex-1] = y_temp;
        }
      }
    }

    m_spline5 = new TSpline5("spline5", &x_loc[0], &y_loc[0], m_spline_numKnots, "b1b2e1e2", a_params[2*m_spline_numKnots], a_params[2*m_spline_numKnots+1],a_params[2*m_spline_numKnots+2], a_params[2*m_spline_numKnots+3]);
    for(int iii = 0; iii < 2*m_spline_numKnots+4; iii++){
      m_splineParams[iii] = a_params[iii];
    }
  }

  return m_spline5->Eval(a_x[0]);
}


double EfficiencyFitter::spline_operator_5_two_point(double* a_x, double* a_params){
  //check to see if the parameters are different
  bool sameParams = m_spline5;
  for(int iii = 0; sameParams && iii < 8; iii++){
    if(m_splineParams[iii] != a_params[iii])  sameParams = false;
  }
  if(!sameParams){ // remake the spline if the parameters are different
    if(m_spline5) delete m_spline5;
    m_spline5 = nullptr;
    double x_loc[2] = {a_params[0], a_params[0] + a_params[1]};
    double y_loc[2] = {a_params[2], a_params[3]};
    m_spline5 = new TSpline5("spline5", x_loc, y_loc, 2, "b1b2e1e2", a_params[4], a_params[5],a_params[6], a_params[7]);
    for(int iii = 0; iii < 8; iii++){
      m_splineParams[iii] = a_params[iii];
    }
  }
  return m_spline5->Eval(a_x[0]);
}






#endif
