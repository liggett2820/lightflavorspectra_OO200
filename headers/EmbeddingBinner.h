#ifndef _EMBEDDING_BINNER_
#define _EMBEDDING_BINNER_

#include "TH2D.h"
#include "TH2I.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TFitResult.h"

#include <iostream>
#include <string>

#include "MattMcTrack.h"
#include "MattMcEvent.h"
#include "CutClass.h"
#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"

using namespace std;

static const int g_EmbBin_MaxNumSysCuts = 5;
static const int g_EmbBin_MaxVzBins = 28;


// This class reprocesses the MattMcEvent TTree after MuDstProcessEmbedding
// This is intended to be run once per particle/charge combination
class EmbeddingBinner{
  
public:

  EmbeddingBinner(string a_particleName, CutClass* a_cuts); // pip pim Kp Km pro pbar deu tri hel alp
  ~EmbeddingBinner();
  void loadEmbeddingOutput(string a_filename);
  void loadDataWeightHistograms(string a_filename, int a_mTm0Rebin = 1); // y mtm0 dca histos to compare to embedding and weight them
  void makeHistograms();
  void fillHistograms();
  void fillWeightedHistograms(); // call this before fillHistograms to have weights
  void closeInputFile();
  void write(string a_outFileName);

private:

  TFile*       m_inFile;
  TTree*       m_tree;
  MattMcEvent* m_event;
  MattMcTrack* m_track;

  CutClass* m_cuts;
  ParticleInfo* m_partInfo;

  int    m_nCentBins;
  int    m_partIndex;
  int    m_particleCharge;
  double m_mass;
  string m_pidSymbol;
  int    m_mTm0DivNum;
  bool   m_useCutClassMtM0Binning;

  int    m_nRapidityBins;
  double m_rapidityMin;
  double m_rapidityMax;

  int    m_nEtaBins;
  double m_etaMin;
  double m_etaMax;

  int    m_nMtM0Bins;
  double m_mTm0Min;
  double m_mTm0Max;


  TH1I* m_centEvents;
  TH2I* m_xyHisto;
  TH1I* m_zHisto;
  TH2I* m_tofMatchVsRefMult;
  TH2I* m_tofTrayMultVsRefMult;
  TH2I* m_refMult_byCentIndex;


  // the reconstructed tracks that were embedded
  TH1I* m_track_NHits_reco;
  TH1I* m_track_NHitsFit_reco;
  TH1I* m_track_NHitsPoss_reco;
  TH1I* m_track_NHitsDedx_reco;
  TH1I* m_track_NHitsFitOverPoss_reco;
  TH1I* m_track_pT_reco;
  TH1I* m_track_mTm0_reco;
  TH1I* m_track_eta_reco;
  TH1I* m_track_rapidity_reco;
  TH1I* m_track_phi_reco;
  TH1I* m_track_pTot_reco;
  TH1I* m_track_pZ_reco;
  TH1I* m_track_gDCA_reco;


  //embeded data (the simulated variables, not reconstructed variables)
  TH1I* m_track_pT_emb;
  TH1I* m_track_mTm0_emb;
  TH1I* m_track_eta_emb;
  TH1I* m_track_rapidity_emb;
  TH1I* m_track_phi_emb;
  TH1I* m_track_pTot_emb;
  TH1I* m_track_pZ_emb;

  TH1I* m_RefMult;           
  TH1I* m_RefMultTrue;        
  TH1I* m_numPrimaryDataTracks;

  TH1I* m_vZ_structure;

  //Efficiency Related Histograms
  TH3I* m_embTrackHisto_byRefMult;
  TH2I* m_embTrackHisto[16];
  TH3I* m_matchTrackHisto_byRefMult;
  TH2I* m_matchTrackHisto[16];
  TH2I* m_embTrackHisto_pTEta[16];
  TH2I* m_matchTrackHisto_pTEta[16];

  TH3I* m_embTrackHisto_byRefMult_byVz[g_EmbBin_MaxVzBins];
  TH2I* m_embTrackHisto_byVz[16][g_EmbBin_MaxVzBins];
  TH3I* m_matchTrackHisto_byRefMult_byVz[g_EmbBin_MaxVzBins];
  TH2I* m_matchTrackHisto_byVz[16][g_EmbBin_MaxVzBins];
  TH2I* m_embTrackHisto_pTEta_byVz[16][g_EmbBin_MaxVzBins];
  TH2I* m_matchTrackHisto_pTEta_byVz[16][g_EmbBin_MaxVzBins];

  //TEfficiencies for weighted efficiencies
  TEfficiency* m_unweighted_eff_y_mTm0[16];
  TEfficiency* m_weighted_eff_y_mTm0[16];
  TH3I*        m_emb_weights_y_mTm0_DCA;
  TH3I*        m_data_weights_y_mTm0_DCA; // this is read-in
  TH2D*        m_dca_eff_y_mTm0_data; // this is calculated during loading
  TEfficiency* m_eff_y_mTm0_noDCA[16]; //2D
  TEfficiency* m_weighted_eff_y_mTm0_noDCA[16]; //2D
  TEfficiency* m_dca_eff_y_mTm0[16];   //2D
  TEfficiency* m_eff_y_mTm0_DCA[16];   //3D

  //Energy loss histograms 
  TH3I* m_pTLossHisto[16];
  TH3I* m_mTm0LossHisto[16];
  TH3I* m_rapidityLossHisto[16];

  //mTm0 Response Matrix
  TH3I* m_pTLossReponseMatrixHisto; // Centrality Integrated
  TH3I* m_mTm0ReponseMatrixHisto; // Centrality Integrated
  TH3I* m_mTm0ReponseMatrixHistoBTOF; // Centrality Integrated


  TH2I* m_btof_tracks_emb[16];
  TH2I* m_tpc_tracks_forBTOF_emb[16];
  TH2I* m_btof_tracks_reco[16];
  TH2I* m_tpc_tracks_forBTOF_reco[16];


};





#endif