
//  Include header files.
#include "TFile.h"
#include "StMessMgr.h"
#include "TH1.h"
#include "TH2.h"
#include <map>
#include <assert.h>
#include <algorithm>

#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDebug.h"
#include "StBTofHeader.h"
#include "StMuAnalysisMaker.h"


/*     ##########    INFORMATION ON SOME OF THE EMBEDDING DATA STRUCTURES ###############

  info from: https://drupal.star.bnl.gov/STAR/event/2020/11/10/simulations-embedding-20201109
  IdTruth:
    IdTruth == -1:
    IdTruth >=10000:                 Real Hits
    IdTruth > 0 && IdTruth < 10000:  MC Hits

    'idTruth > 0' (indicating that the detector level track in question was successfully matched to a generated particle).   https://drupal.star.bnl.gov/STAR/blog/dmawxc/update-02232018-run-9-pp-efficiency-efficiency-using-idtruth


     In the event that idTruth is zero, no MC particle was responsible for hits on the track.  https://drupal.star.bnl.gov/STAR/comp/simu/star-simulation-framework

  QATruth:
   QATruth == -1 for “real” hits
   values 1 to 100 would flag a hit w/ MC contribution
   something w/ 65k (ish) would flag a real hit



   'IdTruth' is the index of the generated particle the reconstructed track was associated with by the embedding software, and
   'QaTruth' is the percentage of TPC hits common to the reconstructed track and its associated generated particle.
    https://drupal.star.bnl.gov/STAR/blog/dmawxc/update-06072019-run-9-embeddingmatched-vs-unmatched-track-comparison-including-those-mat


*/



ClassImp(StMuAnalysisMaker)

//__________________________________________________________________________________________________________________________________________________________


// DCA MATH ----------------------------------------------------------------------------------------------------
double crossProduct3D(double* a_vect_1,double* a_vect_2){
  double f_1 = (a_vect_1[1]*a_vect_2[2]) - (a_vect_1[2]*a_vect_2[1]);
  double f_2 = (a_vect_1[2]*a_vect_2[0]) - (a_vect_1[0]*a_vect_2[2]);
  double f_3 = (a_vect_1[0]*a_vect_2[1]) - (a_vect_1[1]*a_vect_2[0]);
  return sqrt( (f_1*f_1) + (f_2*f_2) + (f_3*f_3) );
}

//__________________________________________________________________________________________________________________________________________________________

double dca_from_decay(double a_primary_x, double a_primary_y, double a_primary_z, double a_secondary_x, double a_secondary_y, double a_secondary_z, double daughter_px, double daughter_py, double daughter_pz){
  double displacement[3] = {a_secondary_x - a_primary_x, a_secondary_y - a_primary_y, a_secondary_z - a_primary_z};
  double distance = sqrt(displacement[0]*displacement[0] + displacement[1]*displacement[1] + displacement[2]*displacement[2]);
  double daughter_pTot = sqrt(daughter_px*daughter_px + daughter_py*daughter_py + daughter_pz*daughter_pz);
  double parent_direction_unit[3]   = {displacement[0]/distance,displacement[1]/distance,displacement[2]/distance};
  double daughter_direction_unit[3] = {daughter_px/daughter_pTot,daughter_py/daughter_pTot,daughter_pz/daughter_pTot};
  return fabs(distance * crossProduct3D(daughter_direction_unit,parent_direction_unit));
}

//__________________________________________________________________________________________________________________________________________________________


// magnetic field should be (0.5T)*10^-11
// locations should be in cm (converted to m in code)
// output is given in cm
double dca_from_decay_helix(double a_magnetic_field, int a_charge_sign, double a_primary_x, double a_primary_y, double a_primary_z, double a_secondary_x, double a_secondary_y, double a_secondary_z, double daughter_px, double daughter_py, double daughter_pz, double &a_dcaX,double &a_dcaY, double &a_dcaZ){
  TVector3 momentum(daughter_px,  daughter_py,  daughter_pz);
  TVector3 decay_vertex(0.01*a_secondary_x,  0.01*a_secondary_y,  0.01*a_secondary_z);
  TVector3 event_vertex(0.01*a_primary_x,  0.01*a_primary_y,  0.01*a_primary_z);
  Helix helix(momentum,decay_vertex,a_magnetic_field,a_charge_sign);
  double dca_tot = 0.0;
  helix.dca(event_vertex,dca_tot,a_dcaX,a_dcaY,a_dcaZ);
  a_dcaX *= 100.0;
  a_dcaY *= 100.0;
  a_dcaZ *= 100.0;
  return 100.0*dca_tot;
  //cout << "DCA Debug:    B:" << a_magnetic_field << "   q:" << a_charge_sign << endl;
  //cout << "Mom: " << "(" << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << ")" << endl;
  //cout << "Decay: " << "(" << decay_vertex.X() << ", " << decay_vertex.Y() << ", " << decay_vertex.Z() << ")" << endl;
  //cout << "Event: " << "(" << event_vertex.X() << ", " << event_vertex.Y() << ", " << event_vertex.Z() << ")" << endl;
  //cout << "DCA: " << fabs(helix.geometricSignedDistance(event_vertex)) << endl;

  //return 100.0*fabs(helix.geometricSignedDistance(event_vertex));
}




bool StMuAnalysisMaker::isQualityEmbeddingMatch(StMuTrack* a_track){
    //if(   primaryTrack->type() != 1
    if( a_track->flag() <= 0
       || abs(a_track->charge()) != 1) return false; 
    m_qaTruth->Fill(a_track->qaTruth()); //before cut
    if(a_track->qaTruth() < 50.0) return false; // Has less than half mc hits
    return true;
}






//__________________________________________________________________________________________________________________________________________________________



// Default Constructor
StMuAnalysisMaker::StMuAnalysisMaker(const Char_t *name) : StMaker(name){
   m_EventCounter = 0;
   m_File = 0;
   m_nCentBins = 16;
   m_currentRunId = -999;
   m_currentRefMult = -999;
}

//__________________________________________________________________________________________________________________________________________________________


StMuAnalysisMaker::~StMuAnalysisMaker(){}; // no destructor

//__________________________________________________________________________________________________________________________________________________________


//  Called once at the beginning of each job to Setup Histograms
Int_t StMuAnalysisMaker::Init(){
  #ifdef _MuEmbDebug_
    cout << "Call to StMuAnalysisMaker::Init()" << endl;
  #endif

  //Setup TTree
  m_tree = new TTree("eventTree","Tree of Data Tracks/Events");
  #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
  m_event = new MattMcFDEvent();
  m_track = new MattMcFDTrack();
  m_tree->Branch("m_event","MattMcFDEvent",&m_event);
  #else
  m_event = new MattMcEvent();
  m_track = new MattMcTrack();
  m_tree->Branch("m_event","MattMcEvent",&m_event);
  #endif

  m_numMakeCalls = 0;
  m_numGoodMakeCalls = 0;
  m_numEventsPerDst = new TH1I("numEventsPerDst","Total Events Per MuDst; Number of Events; Number of MuDsts",1000,0,10000);
  m_numGoodEventsPerDst = new TH1I("numGoodEventsPerDst","Good Events Per MuDst; Number of Events; Number of MuDsts",1000,0,10000);

  m_centEvents           = new TH1I("centEvents","Number of Events for Each Centrality Bin; Centrality Bin; Number of Events",16, -0.5, ((double) 16) - 0.5);
  m_numPrimaryDataTracks = new TH1I("numPriamryDataTracks","Number of Real Primary Data Tracks; Number of real primary tracks; Number of events",250,0,2500);
  m_numMcTracks          = new TH1I("numMcTracks","Number of MC Tracks; Number of MC tracks; Number of events",250,0,2500);
  m_RefMult              = new TH1I("refMult_withEmb", "refMult", 1000, 0., 1000.);
  m_RefMultTrue          = new TH1I("refMult_withoutEmb", "refMultTrue", 1000, 0., 1000.);
  m_qaTruth              = new TH1I("qaTruth", "QA Truth", 101, 0., 101.);
 

  m_rejectedEventsEmb = new TH1I("rejectedEventsEmb","Rejected Events by Emb Code",3,-0.5,2.5);
  m_rejectedEventsEmb->GetXaxis()->SetBinLabel(1,"!muDst");
  m_rejectedEventsEmb->GetXaxis()->SetBinLabel(2,"!muEvent");
  m_rejectedEventsEmb->GetXaxis()->SetBinLabel(3,"SelectedVertex");

  m_rejectedTracksEmb = new TH1I("rejectedTracksEmb","Rejected Tracks After Event Cuts by Emb Code",17,-0.5,16.5);
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(1,"!mcTrack");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(2,"IdVtx");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(3,"GeantId");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(4,"!primaryTrack");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(5,"isData");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(6,"idTruth<=0");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(7,"geantId");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(8,"mcIdVsIdTruth");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(9,"IdVtxRecosMc");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(10,"BadRecoCandidate");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(11,"Found");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(12,"!Found");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(13,"!mcTrack");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(14,"IdVtxMc");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(15,"GeantId");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(16,"Duplicates");
  m_rejectedTracksEmb->GetXaxis()->SetBinLabel(17,"idTruth>numMcTracks");


  #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
  m_mother_delta_decay_distance_div_betagamma = new TH1I("mother_delta_decay_distance_div_betagamma","Decay Distance of Mother Divided by #beta#gamma; Decay Distance Location / #beta#gamma [cm]",100,0,100);
  m_mother_delta_decay_distance = new TH1I("mother_delta_decay_distance","Decay Distance of Mother; Decay Distance Location [cm]",300,0,300); 
  m_mother_decay_location = new TH2I("mother_decay_location","Decay Location of Mother;Z Location [cm]; R Location [cm]",300,-300,300,300,0,400);
  m_mother_delta_decay_location = new TH2I("mother_delta_decay_location","Change in Decay Location of Mother; #Delta Z Location [cm]; #Delta R Location [cm]",300,-300,300,300,0,400);
  #endif

  return StMaker::Init();
}

//__________________________________________________________________________________________________________________________________________________________



// This funciton is called every event after Make().
void StMuAnalysisMaker::Clear(Option_t *opt){
   StMaker::Clear();
}



//__________________________________________________________________________________________________________________________________________________________




//  Called once after all events have been processed by Make()
Int_t StMuAnalysisMaker::Finish(){
  #ifdef _MuEmbDebug_
    cout << "StMuAnalysisMaker::Finish()\n";
    cout << "\tProcessed " << m_EventCounter << " events." << endl;
    cout <<  "\tOutput File Name :" << m_FileName.c_str() << endl;
  #endif
  m_File =  new TFile(m_FileName.c_str(), "RECREATE");
  if(!m_File->IsOpen()){
    cout << "Output File is broken... Can't make or can't read TFile " << endl;
  }

  m_numEventsPerDst->Fill(m_numMakeCalls);
  m_numGoodEventsPerDst->Fill(m_numGoodMakeCalls);

  HistogramUtilities::ConditionalWrite(m_rejectedTracksEmb);
  HistogramUtilities::ConditionalWrite(m_rejectedEventsEmb);
  HistogramUtilities::ConditionalWrite(m_numEventsPerDst);
  HistogramUtilities::ConditionalWrite(m_numGoodEventsPerDst);
  HistogramUtilities::ConditionalWrite(m_centEvents);
  HistogramUtilities::ConditionalWrite(m_numPrimaryDataTracks);
  HistogramUtilities::ConditionalWrite(m_numMcTracks);
  HistogramUtilities::ConditionalWrite(m_RefMult);
  HistogramUtilities::ConditionalWrite(m_RefMultTrue);
  HistogramUtilities::ConditionalWrite(m_qaTruth);
  #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
  HistogramUtilities::ConditionalWrite(m_mother_delta_decay_distance_div_betagamma);
  HistogramUtilities::ConditionalWrite(m_mother_delta_decay_distance);
  HistogramUtilities::ConditionalWrite(m_mother_decay_location);
  HistogramUtilities::ConditionalWrite(m_mother_delta_decay_location);
  #endif
  HistogramUtilities::ConditionalWrite(m_tree);

  if(m_File){
    m_File->Write();
    m_File->Close();
  }

  return kStOK;
}














//__________________________________________________________________________________________________________________________________________________________

//  This method is called every event.
Int_t StMuAnalysisMaker::Make(){
  #ifdef _MuEmbDebug_
    cout << "\n#####################################################################\nCall to StMuAnalysisMaker::Make() with Event Index " << m_EventCounter << endl;
  #endif
  m_EventCounter++;  // increase the event counter
  m_numMakeCalls++;
  StMuDst* muDST = (StMuDst*) GetInputDS("MuDst");
  if(!muDST){
     cout << "ERROR: No MuDST! Ptr: " << muDST << endl;
     m_rejectedEventsEmb->Fill(0);
     return kStOK;        // if no event, we're done
  }

  //  Check StMuEvent branch
  StMuEvent* muEvent = (StMuEvent*) muDST->event();
  if(!muEvent){
    cout << "ERROR: No MuEvent!    Ptr: " << muEvent << endl;
    m_rejectedEventsEmb->Fill(1);
    return kStOK;
  }
  m_event->Clear();

  //----------------------------   vertex selection   -------------------------------------------------
  //int const originalVertexId = muDST->currentVertexIndex();
  //muDST->setVertexIndex(0); // choose the default vertex, i.e. the first vertex
  StMuPrimaryVertex* selectedVertex = muDST->primaryVertex();
  // fall back to default vertex if no vertex is selected in the algorithm above.
  if( !selectedVertex ){
    cout << "WARNING: Vertex is not valid" << endl;
    m_rejectedEventsEmb->Fill(2);
    return kStOK;
    //mMuDst->setVertexIndex(originalVertexId);
  }

  // Vertex Selection from 2024_09_30 Embedding Example
  //int const originalVertexId = mMuDst->currentVertexIndex();
  //mMuDst->setVertexIndex(0); // choose the default vertex, i.e. the first vertex
  //StMuPrimaryVertex* selectedVertex = mMuDst->primaryVertex();
  //if ( ! selectedVertex ) return kStOK;





  /*   FROM OLD VERSION, FXT PRIMARY VERTEX SHOULD ALREADY BE IN TARGET, IF NOT, IT IS REJECTED BY EVENT VZ CUT ANYWAY
  int vtxIndex=-1;
  if( m_cuts->getEventConfig() == "ColliderCenter"){
    vtxIndex = 0;
    //muDST->setVertexIndex(0); //probably not necessary because 0 is probably the default
  }else if( m_cuts->getEventConfig() == "FixedTarget"){
    for(unsigned int iVtx = 0; iVtx < muDST->numberOfPrimaryVertices(); iVtx++){
      StMuPrimaryVertex* vtx = muDST->primaryVertex(iVtx);
      if(!vtx) continue;
      if(vtx->position().z() > m_cuts->getLowZVertexCut() && vtx->position().z() < m_cuts->getHighZVertexCut()){
        muDST->setVertexIndex(iVtx);
        vtxIndex = iVtx;
        break;
      }
    }
    if(vtxIndex < 0){
      cout << "ERROR: THERE WAS NOT VERTEX WITHIN THE TARGET WINDOW Vz=[" <<  m_cuts->getLowZVertexCut() << "," << m_cuts->getHighZVertexCut() << "]" << endl;
    }
  }*/




  TClonesArray* mcVertices  = muDST->mcArray(0);
  TClonesArray* mcTracks    = muDST->mcArray(1);
  Int_t numMcVertices       = mcVertices->GetEntriesFast();
  Int_t numMcTracks         = mcTracks->GetEntriesFast();
  TObjArray* tracks         = muDST->primaryTracks();
  TObjArray* global_tracks  = muDST->globalTracks();
  int numPrimaryTracks      = tracks->GetEntriesFast();
  int numGlobalTracks       = global_tracks->GetEntriesFast();
  int  nPrimaryDataTracks = 0;
  int  numRecoTracks      = 0;
  int  refMult_dataOnly   = 0;
  bool passedEventCuts    = true;
  bool isBadRun           = false;

  if(numMcVertices == 0 || numMcTracks == 0) return kStOK;

  vector< StMuTrack* > primaryTracks;
  vector< int > primaryTracksTimesAccessed;
  for(int primaryIndex = 0; primaryIndex < numPrimaryTracks; primaryIndex++){
    StMuTrack* trk = (StMuTrack*) tracks->UncheckedAt(primaryIndex);
    primaryTracks.push_back(trk);
    primaryTracksTimesAccessed.push_back(0);
  }

  vector< StMuTrack* > globalTracks;
  vector< int > globalTracksTimesAccessed;
  for(int index = 0; index < numGlobalTracks; index++){
    StMuTrack* trk = (StMuTrack*) global_tracks->UncheckedAt(index);
    globalTracks.push_back(trk);
    globalTracksTimesAccessed.push_back(0);
  }


  vector< StMuMcTrack* > simTracks;
  vector< int > simTracksTimesAccessed;
  for(int index = 0; index < numMcTracks; index++){
    StMuMcTrack* trk = (StMuMcTrack*) mcTracks->UncheckedAt(index);
    simTracks.push_back(trk);
    simTracksTimesAccessed.push_back(0);
  }


  //#######################################   CALCULATE REFMULT (Ignoring MC Tracks)    #############################################
  for(int primaryIndex = 0; primaryIndex < numPrimaryTracks; primaryIndex++){
    //StMuTrack* primaryTrack = (StMuTrack*) tracks->UncheckedAt(primaryIndex);
    StMuTrack* primaryTrack = primaryTracks[primaryIndex];    
    primaryTracksTimesAccessed[primaryIndex] += 1;
    if(!primaryTrack) continue;
    if(0 < primaryTrack->idTruth() && primaryTrack->idTruth() <= numMcTracks && isQualityEmbeddingMatch(primaryTrack)){ //is a data track (not a reconstructed simulated track)
      numRecoTracks++;
    }else{
      nPrimaryDataTracks++; //FXT Centrality Variable
      //  These are the requirements of the RefMult variable, not your analysis cuts (Do Not Change)
      if(!(primaryTrack->charge() == 0
              || primaryTrack->nHitsFit() < 10
              || fabs(primaryTrack->momentum().mag()) < 1.e-10
              || fabs(primaryTrack->eta()) > 0.5
              || primaryTrack->dca().mag() >= 3.0)){
        refMult_dataOnly++;
      }
    }
  }

  //#############################################    EVENT LEVEL CUTS   ##########################################################
  #ifdef _STREFMULTCORR_
    if(m_cuts->isBadStRefMultCorrRunId(muEvent->runId())){
      #ifdef _MuEmbDebug_
        cout << "Rejecting Run because of RunId by StRefMultCorr: " << muEvent->runId() << endl;
      #endif
      passedEventCuts = false;
      isBadRun = true;
      m_cuts->getEventRejectionHistogram()->Fill(5.5);
      return kStOK; // -------------------  DONT ALLOW BAD RUNS
    }
    #ifdef _MuEmbDebug_
    else cout << "  Event Passed StRefMultCorrRunId Selection" << endl;
    #endif
  #endif



  //#########      DETERMINE CENTRALITY ##################
  #ifdef _STREFMULTCORR_
    #ifdef _MuEmbDebug_
      cout << "   Setting StRefmultCorr : runID=" << muEvent->runId() << "  refMult_dataOnly=" << refMult_dataOnly << "  z=" << muEvent->primaryVertexPosition().z() << "   ZDCRate=" << muEvent->runInfo().zdcCoincidenceRate()  << endl;
    #endif
    m_cuts->setStRefMultCorrRunId(muEvent->runId());
    m_cuts->setStRefMultCorrEvent(refMult_dataOnly,muEvent->primaryVertexPosition().z(), muEvent->runInfo().zdcCoincidenceRate());
    #ifdef _MuEmbDebug_
      cout << "refMult Before: " << refMult_dataOnly << endl;
    #endif
    refMult_dataOnly = m_cuts->getCorrectedRefMult(); // update refmult with corrected value
    #ifdef _MuEmbDebug_
      cout << "refMult After: " << refMult_dataOnly << endl;
    #endif
  #endif

  //set a new variable just incase you're not using StRefMultCorr
  m_currentRefMult = refMult_dataOnly;

  int centralityIndex = -1;
  if(m_cuts->getEventConfig() == "FixedTarget"){
    centralityIndex = m_cuts->centralityIndex(nPrimaryDataTracks); //Primary Data Tracks
  }else{
    centralityIndex = m_cuts->centralityIndex(refMult_dataOnly); // RefMult (without MC tracks)
  }
  if(centralityIndex < 0){
    m_cuts->getEventRejectionHistogram()->Fill(6.5);
    #ifdef _MuEmbDebug_
      cout << "Bad Centrality - CentBin: " << centralityIndex << "    refMult_dataOnly=" << refMult_dataOnly << endl;
    #endif
    passedEventCuts = false;
    return kStOK; // -------------------  DONT ALLOW BAD CENTRALITY
  }else{
    #ifdef _MuEmbDebug_
      cout << "  Event Passed Centrality Selection" << endl;
    #endif
  }

  #ifdef _MuEmbDebug_
    cout << "   Centrality Bin = " << centralityIndex << "   with refMultCorrected = " << refMult_dataOnly << "   and nPrimaryDataTracks = " << nPrimaryDataTracks << endl;
  #endif

  StMuPrimaryVertex* primaryVtx = muDST->primaryVertex();
  int nBTOFMatch = primaryVtx->nBTOFMatch();                         //Should the nBtofMatch and tray multiplicity change if emebedded tracks are added to event?



  if(m_cuts->isPileUp(refMult_dataOnly,nBTOFMatch)){
    #ifdef _MuEmbDebug_
      cout << "  Event rejected becuase of pileup: refMult_dataOnly = " << refMult_dataOnly << "   and bTof Variable: " << nBTOFMatch << endl;
    #endif
    passedEventCuts = false;
    m_cuts->getEventRejectionHistogram()->Fill(7.5);
    return kStOK; // -------------------  DONT ALLOW PILEUP
  }else{
    #ifdef _MuEmbDebug_
      cout << "  Event Passed Simple Pileup Cut" << endl;
    #endif
  }


  if(m_cuts->getUseHuiDT27GevPileup()){
    #ifdef _MuEmbDebug_
      cout << "WARNING: using Hui Deut Tri 27 GeV Pileup Cuts" << endl;
    #endif
    nBTOFMatch = muEvent->btofTrayMultiplicity();
  }

  //##################   GOOD EVENT CUTS ######################
  if(!m_cuts->isGoodEvent(muDST->primaryVertex()->position().x(),
       muDST->primaryVertex()->position().y(),
       muDST->primaryVertex()->position().z(),
       nBTOFMatch,
       muDST->btofHeader()->nTzero())){
    #ifdef _MuEmbDebug_
      cout << " Rejecting Event by isGoodEvent  " << endl;
    #endif
    passedEventCuts = false;
  }else{
    #ifdef _MuEmbDebug_
      cout << "  Event Passed CutClass isGoodEvent Selection" << endl;
    #endif
  }


  //##################   TRIGGER CUTS ######################
  vector<unsigned int> triggers_for_cuts = m_cuts->getTriggers();
  bool hasValidTrigger = false;
  #ifdef _MuEmbDebug_
    cout << " Trigger Check for " << triggers_for_cuts.size() << "  triggers"  << endl;
  #endif
  for(unsigned int triggerIndex = 0; triggerIndex < triggers_for_cuts.size(); triggerIndex++){
     #ifdef _MuEmbDebug_
       cout << "     Trigger " << triggers_for_cuts[triggerIndex] << " good?   = " << muEvent->triggerIdCollection().nominal().isTrigger(triggers_for_cuts[triggerIndex]) << endl;
     #endif
     if(muEvent->triggerIdCollection().nominal().isTrigger(triggers_for_cuts[triggerIndex])){
       hasValidTrigger = true;
       break;
     };
  }
  if(!hasValidTrigger){
    #ifdef _MuEmbDebug_
      cout << " Rejected by Triggers  " << endl;
    #endif
    passedEventCuts = false;
    m_cuts->getEventRejectionHistogram()->Fill(4.5);
    return kStOK; // -------------------  DONT ALLOW BAD TRIGGERS
  }
  #ifdef _MuEmbDebug_
  else  cout << "  Event Passed Trigger Selection" << endl;
  #endif
  
  #ifdef _MuEmbDebug_
    cout << "  Event " << (passedEventCuts ? "Passed" : "Failed") << " Quality Cuts " << endl;
  #endif


  if(passedEventCuts){
    m_centEvents->Fill(centralityIndex);
    m_RefMultTrue->Fill(refMult_dataOnly);
    m_RefMult->Fill(muEvent->refMult());
    m_numPrimaryDataTracks->Fill(nPrimaryDataTracks);
    m_numMcTracks->Fill(numRecoTracks);
  }


  int trackTreeCounter = 0;
 
  #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_ // ------------------------------------------------------------------------  FEED DOWN ---------------------------------------------------------------
  //cout << "DOING _FEEDDOWN_EMBEDDING_TOGGLE_" << endl;
  // Setup Strange Particle Decay MC Information
  vector< int >    daughterMCTrackGeantIds;
  vector< int >    daughterMCTrackIndicies; // from 0 to numMcTracks-1
  vector< int >    daughterMCTrackBestPrimaryIndex;
  vector< double > daughterMCTrackBestPrimaryQATruth;
  vector< int >    daughterMCTrackBestGlobalIndex;
  vector< double > daughterMCTrackBestGlobalQATruth;
  vector< double > daughterMCTrack_eta;
  vector< double > daughterMCTrack_pT;
  vector< double > daughterMCTrack_phi;
  vector< double > daughterMCTrackMCDCA;
  vector< double > daughterMCTrackMCDCA_Straight;
  vector< double > daughterMCTrackDeltaZ;
  vector< double > daughterMCTrackDeltaR;
  vector< bool >   daughterMCTrackValidDecay;
  vector< bool >   daughterMCTrackFound;



  // Loop through the Monte Carlo tracks and find the primary Lambdas.
  // See if their daughters are a proton and a pion and store their indicies for the reco loop.

  for(Int_t mcTrackIndex=0; mcTrackIndex < numMcTracks; mcTrackIndex++){
    //StMuMcTrack* mcTrackDaughterCandidate = (StMuMcTrack*) mcTracks->UncheckedAt(mcTrackIndex);
    StMuMcTrack* mcTrackDaughterCandidate = simTracks[mcTrackIndex];
    simTracksTimesAccessed[mcTrackIndex] = simTracksTimesAccessed[mcTrackIndex] + 1;
    if(!mcTrackDaughterCandidate) continue;
    #ifdef _MuEmbDebug_
      cout << "MC Track " << mcTrackIndex << " with ID: " << mcTrackDaughterCandidate->GePid() <<  "  vtx: " << mcTrackDaughterCandidate->IdVx() << " -> " << mcTrackDaughterCandidate->IdVxEnd() << endl;
    #endif
    if(!isValidFeedDownID(mcTrackDaughterCandidate->GePid())) continue;
    //Check to see if the mcTrack came from the primary vertex (Embedding Example)
    StMuMcVertex* primaryVertex   = nullptr;
    StMuMcVertex* secondaryVertex = nullptr;
    Int_t idMcVx = mcTrackDaughterCandidate->IdVx();
    int motherIndex = -1;
    vector<int> parent_id_decay_chain;
    parent_id_decay_chain.push_back(mcTrackDaughterCandidate->GePid());
    while(idMcVx != 1){
      StMuMcVertex* mcVertex = (StMuMcVertex *) mcVertices->UncheckedAt(idMcVx-1);
      Int_t idMcParentTrack = mcVertex->IdParTrk();
      if(!idMcParentTrack) break;
      //StMuMcTrack* mcTrackParent = (StMuMcTrack *) mcTracks->UncheckedAt(idMcParentTrack-1);
      StMuMcTrack* mcTrackParent = simTracks[idMcParentTrack-1];
      simTracksTimesAccessed[idMcParentTrack-1] = simTracksTimesAccessed[idMcParentTrack-1] + 1;
      idMcVx = mcTrackParent->IdVx();
      parent_id_decay_chain.push_back(mcTrackParent->GePid());
      if(idMcVx == 1 && (mcTrackParent->GePid() == m_mother_id)){ // IS THE EMBEDDED MOTHER
        motherIndex = idMcParentTrack;
        primaryVertex = (StMuMcVertex *) mcVertices->UncheckedAt(idMcVx - 1);
        int idMcVxEnd = mcTrackParent->IdVxEnd();
        if(0 < idMcVxEnd && idMcVxEnd <= numMcVertices){
          secondaryVertex = (StMuMcVertex *) mcVertices->UncheckedAt(mcTrackParent->IdVxEnd() - 1); 
          m_mother_decay_location->Fill(secondaryVertex->XyzV().z(),sqrt(pow(secondaryVertex->XyzV().x(),2) + pow(secondaryVertex->XyzV().y(),2)));
          double dist = sqrt(pow(secondaryVertex->XyzV().x() - primaryVertex->XyzV().x(),2) + pow(secondaryVertex->XyzV().y() - primaryVertex->XyzV().y(),2) + pow(secondaryVertex->XyzV().z() - primaryVertex->XyzV().z(),2));
          m_mother_delta_decay_distance->Fill(dist);
          m_mother_delta_decay_location->Fill(secondaryVertex->XyzV().z() - primaryVertex->XyzV().z(),sqrt(pow(secondaryVertex->XyzV().x() - primaryVertex->XyzV().x(),2) + pow(secondaryVertex->XyzV().y() - primaryVertex->XyzV().y(),2)));
          double pTot = mcTrackParent->Ptot();
          double ETot = mcTrackParent->E();
          double betagamma = pTot/sqrt(pow(ETot,2)-pow(pTot,2));
          m_mother_delta_decay_distance_div_betagamma->Fill(dist/betagamma);
        }
      }
      if(!idMcVx) break;
    }
    if(idMcVx != 1 || motherIndex == -1) continue; //this MC track is not eventually originated from PV
    bool validDecay = isValidDecayChain(parent_id_decay_chain);
    //cout << (validDecay ? "Valid Decay" : "INVALID Decay") << endl;
    if( !validDecay ){
      #ifdef _MuEmbDebug_
      cout << "Invalid Decay Chain: ";
      for(int iii = 0; iii < parent_id_decay_chain.size(); iii++) cout << parent_id_decay_chain[iii] << "  ";
      cout << endl;
      #endif
      continue;
    }

    StMuMcVertex* decayVertex = (StMuMcVertex *) mcVertices->UncheckedAt(mcTrackDaughterCandidate->IdVx() - 1);
    daughterMCTrackGeantIds.push_back(mcTrackDaughterCandidate->GePid());
    daughterMCTrackIndicies.push_back(mcTrackIndex); // from 0 to numMcTracks-1
    daughterMCTrackBestPrimaryIndex.push_back(-1);
    daughterMCTrackBestPrimaryQATruth.push_back(-1.0);
    daughterMCTrackBestGlobalIndex.push_back(-1);
    daughterMCTrackBestGlobalQATruth.push_back(-1.0);
    daughterMCTrack_eta.push_back(mcTrackDaughterCandidate->Pxyz().pseudoRapidity());
    daughterMCTrack_pT.push_back(mcTrackDaughterCandidate->Pxyz().perp());
    daughterMCTrack_phi.push_back(atan2(mcTrackDaughterCandidate->Pxyz().y(), mcTrackDaughterCandidate->Pxyz().x()));
    daughterMCTrackDeltaZ.push_back(decayVertex->XyzV().z() - primaryVertex->XyzV().z());
    daughterMCTrackDeltaR.push_back(sqrt(pow(decayVertex->XyzV().x() - primaryVertex->XyzV().x(),2) + pow(decayVertex->XyzV().y() - primaryVertex->XyzV().y(),2)));
    double mc_dcaX = 0;
    double mc_dcaY = 0;
    double mc_dcaZ = 0;
    double mc_dca_straight = dca_from_decay(primaryVertex->XyzV().x(), primaryVertex->XyzV().y(), primaryVertex->XyzV().z(), decayVertex->XyzV().x(), decayVertex->XyzV().y(), decayVertex->XyzV().z(), mcTrackDaughterCandidate->Pxyz().x(), mcTrackDaughterCandidate->Pxyz().y(), mcTrackDaughterCandidate->Pxyz().z());
    double mc_dca = dca_from_decay_helix(-0.5E-11,mcTrackDaughterCandidate->Charge(),primaryVertex->XyzV().x(), primaryVertex->XyzV().y(), primaryVertex->XyzV().z(), decayVertex->XyzV().x(), decayVertex->XyzV().y(), decayVertex->XyzV().z(), mcTrackDaughterCandidate->Pxyz().x(), mcTrackDaughterCandidate->Pxyz().y(), mcTrackDaughterCandidate->Pxyz().z(),mc_dcaX,mc_dcaY,mc_dcaZ);
    daughterMCTrackMCDCA.push_back(mc_dca);
    daughterMCTrackMCDCA_Straight.push_back(mc_dca_straight);
    daughterMCTrackFound.push_back(false);
  }


  //    ####################   RECO  PRIMARY TRACK LOOP   ####################

  for(int primaryIndex = 0; primaryIndex < numPrimaryTracks; primaryIndex++){
    //primaryTrack = (StMuTrack*) tracks->UncheckedAt(primaryIndex);
    StMuTrack* primaryTrack = primaryTracks[primaryIndex];
    primaryTracksTimesAccessed[primaryIndex] = primaryTracksTimesAccessed[primaryIndex] + 1;
    if(primaryTrack->idTruth() <= 0 || primaryTrack->idTruth() > numMcTracks){
      continue;
    }
    //StMuMcTrack* mcTrack = (StMuMcTrack *) mcTracks->UncheckedAt(primaryTrack->idTruth()-1);
    StMuMcTrack* mcTrack = simTracks[primaryTrack->idTruth()-1];
    simTracksTimesAccessed[primaryTrack->idTruth()-1] += 1;
    if(!mcTrack){
      LOG_WARN << "Inconsistency in mcArray(1), ignored" << endm;
      continue;
    }
    if(mcTrack->Id() != primaryTrack->idTruth()){
      LOG_WARN << "Mismatched idTruth " << primaryTrack->idTruth() << " and mcTrack Id " <<  mcTrack->Id()
               << " this track is ignored" <<  endm;
      cout << "Mismatched idTruth " << primaryTrack->idTruth() << " and mcTrack Id " <<  mcTrack->Id()
               << " this track is ignored" << endl;
      continue;
    }
    //cout << "Primary Track Index: " << primaryIndex << "   ptr: "  << primaryTrack << endl; 
    // Note: Futher cuts will be applied after the TTree is created. These are the widest cuts
    if(!isQualityEmbeddingMatch(primaryTrack)){
      continue;
    }



    int daughter_index = -1;
    for(unsigned int index = 0; index < daughterMCTrackIndicies.size(); index++){
      if(daughterMCTrackIndicies[index] == primaryTrack->idTruth()-1){ // off by one on purpose daughterMCTrackIndicies[] starts at 0, idTruth starts at 1
        daughter_index = index;
        break;
      }
    }
    if(daughter_index != -1){ //THIS IS A RECONSTRUCTED TRACK ASSOCIATED WITH A MC TRACK FROM A FEED DOWN PARENT
      if(daughterMCTrackBestPrimaryQATruth[daughter_index] < primaryTrack->qaTruth()){
        daughterMCTrackBestPrimaryIndex[daughter_index] = primaryIndex;
        daughterMCTrackBestPrimaryQATruth[daughter_index] = primaryTrack->qaTruth();
        daughterMCTrackFound[daughter_index] = true;
      }
    } //end if daugther found
  } // end of primary track loop
  //cout << "Finished Primary Track Loop" << endl;

  for(int index = 0; index < numGlobalTracks; index++){
    //primaryTrack = (StMuTrack*) tracks->UncheckedAt(primaryIndex);
    StMuTrack* globalTrack = globalTracks[index];
    //cout << "Global Track: " << index << " of " << numGlobalTracks << endl;
    globalTracksTimesAccessed[index] += 1;
    if(globalTrack->idTruth() <= 0 || globalTrack->idTruth() > numMcTracks){
      continue;
    }
    //StMuMcTrack* mcTrack = (StMuMcTrack *) mcTracks->UncheckedAt(primaryTrack->idTruth()-1);
    StMuMcTrack* mcTrack = simTracks[globalTrack->idTruth()-1];
    simTracksTimesAccessed[globalTrack->idTruth()-1] = simTracksTimesAccessed[globalTrack->idTruth()-1] + 1;
    if(!mcTrack){
      LOG_WARN << "Inconsistency in mcArray(1), ignored" << endm;
      continue;
    }
    if(mcTrack->Id() != globalTrack->idTruth()){
      LOG_WARN << "Mismatched idTruth " << globalTrack->idTruth() << " and mcTrack Id " <<  mcTrack->Id()
               << " this track is ignored" <<  endm;
      cout << "Mismatched idTruth " << globalTrack->idTruth() << " and mcTrack Id " <<  mcTrack->Id()
               << " this track is ignored" << endl;
      continue;
    }
    //cout << "Primary Track Index: " << primaryIndex << "   ptr: "  << primaryTrack << endl; 
    // Note: Futher cuts will be applied after the TTree is created. These are the widest cuts
    if(!isQualityEmbeddingMatch(globalTrack)){
      //cout << "isn't quality" << endl;
      continue;
    }
    //cout << " is quality " << endl;


    int daughter_index = -1;
    for(unsigned int mcindex = 0; mcindex < daughterMCTrackIndicies.size(); mcindex++){
      //cout << "Finding mc track: " << daughterMCTrackIndicies[mcindex] << " =? " << globalTrack->idTruth()-1 << endl;
      if(daughterMCTrackIndicies[mcindex] == globalTrack->idTruth()-1){ // off by one on purpose daughterMCTrackIndicies[] starts at 0, idTruth starts at 1
        daughter_index = mcindex;
        break;
      }
    }
    //cout << "Daughter Index = " << daughter_index << endl;
    if(daughter_index != -1){ //THIS IS A RECONSTRUCTED TRACK ASSOCIATED WITH A MC TRACK FROM A FEED DOWN PARENT
      if(daughterMCTrackBestGlobalQATruth[daughter_index] < globalTrack->qaTruth()){
        daughterMCTrackBestGlobalIndex[daughter_index] = index;
        daughterMCTrackBestGlobalQATruth[daughter_index] = globalTrack->qaTruth();
      }
    } //end if daugther found

  }



  //#####################    Find the MC Tracks that were not reconstructed  #########################################
  //cout << "MONTE CARLO TRACKS WE CARE ABOUT: " << endl;
  for(unsigned int daughter_index = 0; daughter_index < daughterMCTrackFound.size(); daughter_index++){
      m_track->Clear();
      m_track->set_pT_emb(daughterMCTrack_pT[daughter_index]);
      m_track->set_eta_emb(daughterMCTrack_eta[daughter_index]);
      m_track->set_phi_emb(daughterMCTrack_phi[daughter_index]);
      m_track->set_particle_id(daughterMCTrackGeantIds[daughter_index]);
      m_track->set_dca_emb(daughterMCTrackMCDCA[daughter_index]);
      m_track->set_decay_delta_R(daughterMCTrackDeltaR[daughter_index]);
      m_track->set_decay_delta_Z(daughterMCTrackDeltaZ[daughter_index]);
      //cout << daughterMCTrackGeantIds[daughter_index] << " pT: " << daughterMCTrack_pT[daughter_index] << "  eta: " << daughterMCTrack_eta[daughter_index] << " mcDCA: " << daughterMCTrackMCDCA[daughter_index]; 
      if(daughterMCTrackFound[daughter_index]){
        StMuTrack* primaryTrack = primaryTracks[daughterMCTrackBestPrimaryIndex[daughter_index]];
        primaryTracksTimesAccessed[daughterMCTrackBestPrimaryIndex[daughter_index]] = primaryTracksTimesAccessed[daughterMCTrackBestPrimaryIndex[daughter_index]] + 1; 
        m_track->set_isPrimary(true);
        m_track->set_qaTruth(primaryTrack->qaTruth());
        m_track->set_pT(primaryTrack->pt());
        m_track->set_eta(primaryTrack->eta());
        m_track->set_phi(primaryTrack->phi());
        m_track->set_gDCA(primaryTrack->dcaGlobal().mag());
        m_track->set_nHitsFit(primaryTrack->nHitsFit());
        m_track->set_nHitsDeDx(primaryTrack->nHitsDedx());
        m_track->set_nHitsRatio(((double )primaryTrack->nHitsFit())/((double)primaryTrack->nHitsPoss()));
        if(primaryTrack->btofPidTraits().matchFlag() != 0){
          m_track->set_btofLocalY(primaryTrack->btofPidTraits().yLocal());
          m_track->set_btofLocalZ(primaryTrack->btofPidTraits().zLocal());
        }
        //cout << "  PRIMARY     " << primaryTrack->dcaGlobal().mag();
      }else if(daughterMCTrackBestGlobalQATruth[daughter_index] > 0){ // Look into globals if not in primaries
        StMuTrack* globalTrack = globalTracks[daughterMCTrackBestGlobalIndex[daughter_index]];
        globalTracksTimesAccessed[daughterMCTrackBestGlobalIndex[daughter_index]] += 1;
        m_track->set_isPrimary(false);
        m_track->set_qaTruth(globalTrack->qaTruth());
        m_track->set_pT(globalTrack->pt());
        m_track->set_eta(globalTrack->eta());
        m_track->set_phi(globalTrack->phi());
        m_track->set_gDCA(globalTrack->dcaGlobal().mag());
        m_track->set_nHitsFit(globalTrack->nHitsFit());
        m_track->set_nHitsDeDx(globalTrack->nHitsDedx());
        m_track->set_nHitsRatio(((double )globalTrack->nHitsFit())/((double)globalTrack->nHitsPoss()));
        if(globalTrack->btofPidTraits().matchFlag() != 0){
          m_track->set_btofLocalY(globalTrack->btofPidTraits().yLocal());
          m_track->set_btofLocalZ(globalTrack->btofPidTraits().zLocal());
        }
        //cout << "  GLOBAL    " << globalTrack->dcaGlobal().mag();
      }
      cout << endl;
      m_event->addTrack(trackTreeCounter,m_track);
      trackTreeCounter++;
  }




  #ifdef _MuEmbDebug_
  cout << "ACCESSED COUNT FOR MC TRACKS" << endl;
  for(int index = 0; index < numMcTracks; index++){
    if(simTracksTimesAccessed[index] >= 1){
      cout << simTracks[index]->GePid() << " accessed " << simTracksTimesAccessed[index] << " times" << endl;
    }
  }

  for(int primaryIndex = 0; primaryIndex < numPrimaryTracks; primaryIndex++){
    if(primaryTracksTimesAccessed[primaryIndex] > 1){
      cout << primaryIndex << " accessed " << primaryTracksTimesAccessed[primaryIndex] << " times" << endl;
    }
  }

  #endif




  #else    // ------------------------------------------------------------------------ NOT FEED DOWN ---------------------------------------------------------------



  vector<int> foundMcTracks;
  //##############################    SIMULATED PARTICLE LOOP  ###############################################
  for(Int_t mcTrackIndex = 0; mcTrackIndex < numMcTracks; mcTrackIndex++){
    StMuMcTrack* mcTrack = (StMuMcTrack *) mcTracks->UncheckedAt(mcTrackIndex);
    if(!mcTrack){
      m_rejectedTracksEmb->Fill(0);
      continue;
    }

    Int_t IdVtx = mcTrack->IdVx();
    if(IdVtx != 1){
      m_rejectedTracksEmb->Fill(1);
      continue; // Select only Triggered Mc Vertex, i.e. the MC track should originate from PV (IdVx=1)
    }
    const int geantId = mcTrack->GePid();
    if(!isParticleOfInterest(geantId)){
      m_rejectedTracksEmb->Fill(2);
      continue;
    }
  }



  //#######################################      Data/Reco TRACK LOOP BEFORE EVENT CUTS   #############################################

  for(int primaryIndex = 0; primaryIndex < numPrimaryTracks; primaryIndex++){
    StMuTrack* primaryTrack = (StMuTrack*) tracks->UncheckedAt(primaryIndex);
    if(!primaryTrack){
      m_rejectedTracksEmb->Fill(3);
      continue;
    }

    #ifdef _MuEmbDebug_
      cout << " Processing muTrack Primary...   # " << primaryIndex << " / " << numPrimaryTracks  << endl;
    #endif

    bool isDataTrack = (primaryTrack->idTruth() > numMcTracks);
    bool passedTrackCuts = m_cuts->isGoodTrack((primaryTrack->type() == 1),primaryTrack->nHitsFit(), primaryTrack->nHitsPoss(), primaryTrack->nHitsDedx(),primaryTrack->dcaGlobal().mag(), primaryTrack->phi());
    if(isDataTrack){
      // this is a primary track that is not a simulated track (this is the data originally in the file before embedding)
      m_rejectedTracksEmb->Fill(4);
      continue; //Below is only for MC Reconstructed Tracks
    }

    //#######################################   Determine the associated simulation track   #############################################
    if(primaryTrack->idTruth() <= 0){
      //cout << "WARNING: Invalid idTruth for primaryTrack" << endl;
      m_rejectedTracksEmb->Fill(5);
      continue;
    }
    if(primaryTrack->idTruth() > numMcTracks){
      //cout << "WARNING: Invalid idTruth for primaryTrack" << endl;
      m_rejectedTracksEmb->Fill(16);
      continue;
    }



    StMuMcTrack* mcTrack = (StMuMcTrack *) mcTracks->UncheckedAt(primaryTrack->idTruth()-1);
    const int geantId = mcTrack->GePid();
    if(!isParticleOfInterest(geantId)){
      #ifdef _MuEmbDebug_
        cout << "WARNING: The MC Track Isn't the particle you care about" << endl;
      #endif
      m_rejectedTracksEmb->Fill(6);
      continue;
    }
    if(mcTrack->Id() != primaryTrack->idTruth()){
      cout << "WARNING: Mismatched idTruth " << primaryTrack->idTruth() << " and mcTrack Id " <<  mcTrack->Id()  << " this track is ignored" <<  endl;
      m_rejectedTracksEmb->Fill(7);
      //continue; // this continue wasn't in the example script
    }
    // VERTEX LOGIC FROM EMBEDDING EXAMPLE
    Int_t idMcVx = mcTrack->IdVx();
    /*while(idMcVx != 1){
      StMuMcVertex* mcVertex = (StMuMcVertex *) mcVertices->UncheckedAt(idMcVx-1);
      Int_t idMcTrack = mcVertex->IdParTrk();
      if (!idMcTrack) break;
      StMuMcTrack*  mcTrack  = (StMuMcTrack *) mcTracks->UncheckedAt(idMcTrack-1);
      idMcVx = mcTrack->IdVx();
      if(!idMcVx) break;
    }*/
    if(idMcVx != 1){
      cout << "WARNING: Simulated Track did not originate with the primary vertex" << endl;
      m_rejectedTracksEmb->Fill(8);
      continue; //this MC track is not eventually originated from PV
    }
    //if(mcTrack->IdVx() != 1){
       //cout<<"WARNING: mc track may not directly originate from PV!"<<endl;
       //There was no continue in the embedding example
    //}
    //Selection Criteria from Example Code
    bool isBadRecoCandidate = false;
    if(primaryTrack->qaTruth() < 50.0 || primaryTrack->flag() <= 0 ) isBadRecoCandidate = true; //|| abs(primaryTrack->charge())!=1
    if(abs(m_particleCharge) == 1 && m_particleCharge != primaryTrack->charge()) isBadRecoCandidate = true; // protection for deuterons and charge 2 or 3 particles
    if(isBadRecoCandidate){
      m_rejectedTracksEmb->Fill(9);
      continue;
    }
    m_rejectedTracksEmb->Fill(17);

    //Fill the embedding track info even if reco is bad
    /*if(isBadRecoCandidate){
      if(!isBadRun){
        m_track->Clear();
        m_track->set_pT_emb(mc_pT);
        m_track->set_eta_emb(mc_eta);
        m_track->set_phi_emb(mc_phi);
        m_event->addTrack(trackTreeCounter,m_track);
        trackTreeCounter++;
        foundMcTracks.push_back(mcTrack->Id());
      }
      m_rejectedTracksEmb->Fill(9);
      continue;
    }
    m_rejectedTracksEmb->Fill(17);
    */

    #ifdef _MuEmbDebug_
      cout << " Filling muTrack Reco Information... " << endl;
    #endif

    m_track->Clear();
    m_track->set_pT_emb(mcTrack->pT());
    m_track->set_eta_emb(mcTrack->Eta());
    m_track->set_phi_emb(atan2(mcTrack->Pxyz().y(), mcTrack->Pxyz().x()));
    m_track->set_pT(primaryTrack->pt());
    m_track->set_eta(primaryTrack->eta());
    m_track->set_phi(primaryTrack->phi());
    m_track->set_gDCA(primaryTrack->dcaGlobal().mag());
    m_track->set_nHitsFit(primaryTrack->nHitsFit());
    m_track->set_nHitsDeDx(primaryTrack->nHitsDedx());
    m_track->set_nHitsRatio(((double )primaryTrack->nHitsFit())/((double)primaryTrack->nHitsPoss()));
    if(primaryTrack->btofPidTraits().matchFlag() != 0){
      m_track->set_btofLocalY(primaryTrack->btofPidTraits().yLocal());
      m_track->set_btofLocalZ(primaryTrack->btofPidTraits().zLocal());
    }
    m_event->addTrack(trackTreeCounter,m_track);
    trackTreeCounter++;
    foundMcTracks.push_back(mcTrack->Id());

  }


  //#####################    Find the MC Tracks that were not reconstructed  #########################################
  for(Int_t mcTrackIndex = 0; !isBadRun && mcTrackIndex < numMcTracks; mcTrackIndex++){
    StMuMcTrack* mcTrack = (StMuMcTrack *) mcTracks->UncheckedAt(mcTrackIndex);
    if(!mcTrack) continue;
    bool alreadyFound = false;

    for(int iii = 0; iii < foundMcTracks.size(); iii++){
      if(mcTrack->Id() == foundMcTracks[iii]){
        alreadyFound = true;
        break;
      }
    }

    #ifdef _MuEmbDebug_
      cout << " Processing mcTrack ... #: " << mcTrackIndex << "   Found?  " << (alreadyFound ? "YES" : "NO") << endl;
    #endif

    if(alreadyFound){
      m_rejectedTracksEmb->Fill(10);
      continue;
    }else{
      m_rejectedTracksEmb->Fill(11);
    }
    Int_t IdVtx = mcTrack->IdVx();
    if(IdVtx != 1){
      m_rejectedTracksEmb->Fill(13);
      continue; // Select only Triggered Mc Vertex, i.e. the MC track should originate from PV (IdVx=1)
    }
    const int geantId = mcTrack->GePid();
    if(!isParticleOfInterest(geantId)){
      m_rejectedTracksEmb->Fill(14);
      continue;
    }

    m_track->Clear();
    m_track->set_pT_emb(mcTrack->pT());
    m_track->set_eta_emb(mcTrack->Eta());
    m_track->set_phi_emb(atan2(mcTrack->Pxyz().y(), mcTrack->Pxyz().x()));
    m_event->addTrack(trackTreeCounter,m_track);
    trackTreeCounter++;
    foundMcTracks.push_back(mcTrack->Id());
  }


  //Test for duplicates
  #ifdef _MuEmbDebug_
    cout << "Testing for duplicates:" << endl;
  #endif
  for(unsigned int iii = 0; iii < foundMcTracks.size(); iii++){
    for(unsigned int jjj = iii+1; jjj < foundMcTracks.size(); jjj++){
      if(foundMcTracks[iii] == foundMcTracks[jjj]){
        m_rejectedTracksEmb->Fill(15);
      }
    }
  }


  #endif // ------   END NORMAL MODE (NO FEEDDOWN) ----------------------------

  m_event->set_refMult(muEvent->refMult());
  m_event->set_refMult_data_only(refMult_dataOnly);
  m_event->set_numPrimaryTracks(nPrimaryDataTracks);
  m_event->set_vertex_X(muEvent->primaryVertexPosition().x());
  m_event->set_vertex_Y(muEvent->primaryVertexPosition().y());
  m_event->set_vertex_Z(muEvent->primaryVertexPosition().z());
  m_event->set_numTofMatch(nBTOFMatch);
  m_event->set_numTofT0(muDST->btofHeader()->nTzero());
  if(centralityIndex >= 0) m_event->set_centIndex(centralityIndex); // the storage cant be negative (default 250)

  #ifdef _MuEmbDebug_
    cout << "Filling Tree: " << m_event << endl;
    cout << "Testing Event: ";
    m_event->print();
  #endif

  m_tree->Fill();
  if(passedEventCuts) m_numGoodMakeCalls++;





  #ifdef _MuEmbDebug_
    cout << "  Done Processing the MAKE funciton" << endl;
  #endif


  return kStOK;
}



#ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
bool StMuAnalysisMaker::isValidFeedDownID(int a_geantId){
  return (a_geantId == 8 
       || a_geantId == 9 
       || a_geantId == 11 
       || a_geantId == 12 
       || a_geantId == 14 
       || a_geantId == 15);
}

bool StMuAnalysisMaker::isValidDecayChain(vector<int> a_decay_chain){
  //cout << "DECAY CHAIN: ";
  //for(int iii = 0; iii < a_decay_chain.size(); iii++){
  //  cout << "  " << a_decay_chain[iii];
  //}
  //cout << endl;

  // Note: there are decays like DECAY CHAIN:   14  14  98  40005, that are thrown out by this selection (probably scatters?)
  if(a_decay_chain.size() <= 1) return false;
  if(a_decay_chain[a_decay_chain.size()-1] == 10151) return true; //phi is difficult, this is a bypass
  if(a_decay_chain.size() == 2){
    if(a_decay_chain[0] == 14 && a_decay_chain[1] == 10018) return true; // p L
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 10018) return true; // pi- L
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 40003) return true; // pi- Xi-
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 40002) return true; // pi- Omega-
    if(a_decay_chain[0] == 12 && a_decay_chain[1] == 40002) return true; // K- Omega-
    //Antiparticles
    if(a_decay_chain[0] == 15 && a_decay_chain[1] == 10026) return true; // pbar Lbar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 10026) return true; // pi+ Lbar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 40004) return true; // pi+ XiBar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 40001) return true; // pi+ OmegaBar
    if(a_decay_chain[0] == 11 && a_decay_chain[1] == 40001) return true; // K+ OmegaBar

    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 707) return true; // pi- K0S
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 707) return true; // pi+ K0S
  }
  else if(a_decay_chain.size() == 3){
    if(a_decay_chain[0] == 14 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40003) return true; // p L Xi-
    if(a_decay_chain[0] == 14 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40005) return true; // p L Xi0
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40003) return true; // pi- L Xi-
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40005) return true; // pi- L Xi0
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40002) return true; // pi- Xi- Omega-
    //Antiparticles
    if(a_decay_chain[0] == 15 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40004) return true; // pbar Lbar XiBar
    if(a_decay_chain[0] == 15 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40006) return true; // pbar Lbar Xi0Bar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40004) return true; // pi+ Lbar XiBar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40006) return true; // pi+ Lbar Xi0Bar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40001) return true; // pi+ XiBar OmegaBar

    if(a_decay_chain[0] == 14 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40002) return true; // p L Xi- Omega-
    if(a_decay_chain[0] == 14 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40002) return true; // p L Xi0 Omega-
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40002) return true; // pi- L Xi- Omega-
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 98 && a_decay_chain[2] == 40002) return true; // pi- L Xi0 Omega-
    //Antiparticles
    if(a_decay_chain[0] == 15 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40001) return true; // pbar Lbar XiBar OmegaBar
    if(a_decay_chain[0] == 15 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40001) return true; // pbar Lbar Xi0Bar OmegaBar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40001) return true; // pi+ Lbar XiBar OmegaBar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 97 && a_decay_chain[2] == 40001) return true; // pi+ Lbar Xi0Bar OmegaBar


  }
  /*else if(a_decay_chain.size() == 4){
    if(a_decay_chain[0] == 14 && a_decay_chain[1] == 10018 && a_decay_chain[2] == 40003 && a_decay_chain[3] == 40002) return true; // p L Xi- Omega-
    if(a_decay_chain[0] == 14 && a_decay_chain[1] == 10018 && a_decay_chain[2] == 40005 && a_decay_chain[3] == 40002) return true; // p L Xi0 Omega-
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 10018 && a_decay_chain[2] == 40003 && a_decay_chain[3] == 40002) return true; // pi- L Xi- Omega-
    if(a_decay_chain[0] == 9 && a_decay_chain[1] == 10018 && a_decay_chain[2] == 40005 && a_decay_chain[3] == 40002) return true; // pi- L Xi0 Omega-
    //Antiparticles
    if(a_decay_chain[0] == 15 && a_decay_chain[1] == 10026 && a_decay_chain[2] == 40004 && a_decay_chain[3] == 40001) return true; // pbar Lbar XiBar OmegaBar
    if(a_decay_chain[0] == 15 && a_decay_chain[1] == 10026 && a_decay_chain[2] == 40006 && a_decay_chain[3] == 40001) return true; // pbar Lbar Xi0Bar OmegaBar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 10026 && a_decay_chain[2] == 40004 && a_decay_chain[3] == 40001) return true; // pi+ Lbar XiBar OmegaBar
    if(a_decay_chain[0] == 8 && a_decay_chain[1] == 10026 && a_decay_chain[2] == 40006 && a_decay_chain[3] == 40001) return true; // pi+ Lbar Xi0Bar OmegaBar
  }*/
  return false;
}




#else

bool StMuAnalysisMaker::isParticleOfInterest(int a_geantID){
  int partIndex = -1;
  int charge = 0;
  if     (a_geantID == 8 ) {partIndex = 0; charge =  1;} // pion
  else if(a_geantID == 9 ) {partIndex = 0; charge = -1;}// pion
  else if(a_geantID == 11) {partIndex = 1; charge =  1;} // kaon
  else if(a_geantID == 12) {partIndex = 1; charge = -1;}// kaon
  else if(a_geantID == 14) {partIndex = 2; charge =  1;}// proton
  else if(a_geantID == 15) {partIndex = 2; charge = -1;}// proton
  else if(a_geantID == 45) {partIndex = 4; charge =  2;} // deuteron
  else if(a_geantID == 46) {partIndex = 5; charge =  3;} // triton
  else if(a_geantID == 49) {partIndex = 6; charge =  2;} // helion
  else if(a_geantID == 47) {partIndex = 7; charge =  2;} // alpha
  if(partIndex != m_partIndex || charge != m_particleCharge){
    return false;
  }
  return true;
}


void StMuAnalysisMaker::setPart(string partCStr){
  TString partStr = Form("%s",partCStr.c_str());
  if( partStr.EqualTo("pip",TString::kIgnoreCase) ){
    m_partIndex = 0;
    m_particleCharge = 1;
    m_mass = 0.13957;
  }else if( partStr.EqualTo("pim",TString::kIgnoreCase) ){
    m_partIndex = 0;
    m_particleCharge = -1;
    m_mass = 0.13957;
  }else if( partStr.EqualTo("Kp",TString::kIgnoreCase) ){
    m_partIndex = 1;
    m_particleCharge = 1;
    m_mass = 0.49367;
  }else if( partStr.EqualTo("Km",TString::kIgnoreCase) ){
    m_partIndex = 1;
    m_particleCharge = -1;
    m_mass = 0.49367;
  }else if( partStr.EqualTo("pro",TString::kIgnoreCase) ){
    m_partIndex = 2;
    m_particleCharge = 1;
    m_mass = 0.93827;
  }else if( partStr.EqualTo("pbar",TString::kIgnoreCase) ){
    m_partIndex = 2;
    m_particleCharge = -1;
    m_mass = 0.93827;
  }else if( partStr.EqualTo("deu",TString::kIgnoreCase) ){
    m_partIndex = 4;
    m_particleCharge = 2;
    m_mass = 1.87561;
  }else if( partStr.EqualTo("tri",TString::kIgnoreCase) ){
    m_partIndex = 5;
    m_particleCharge = 3;
    m_mass = 2.80925;
  }else if( partStr.EqualTo("hel",TString::kIgnoreCase) ){
    m_partIndex = 6;
    m_particleCharge = 2;
    m_mass = 2.8094;
  }else if( partStr.EqualTo("alp",TString::kIgnoreCase) ){
    m_partIndex = 7;
    m_particleCharge = 2;
    m_mass = 3.72742;
  }
  return;
}

#endif
