/*!
 * \class  StMuAnalysisMaker
 * \brief  A typical Analysis Class for MuDst
 * \author Wei-Ming Zhang, KSU, Mar 2004
 *
 * This is an example of a maker to perform analysis using MuDst.
 *
 * $Id: StMuAnalysisMaker.h,v 1.2 2014/08/06 11:43:31 jeromel Exp $
 *
 * -------------------------------------------------------------------------
 * $Log: StMuAnalysisMaker.h,v $
 * Revision 1.2  2014/08/06 11:43:31  jeromel
 * Suffix on literals need to be space (later gcc compiler makes it an error) - first wave of fixes
 *
 * Revision 1.1  2004/08/10 16:09:11  perev
 * new GridCollector stuff
 *
 * -------------------------------------------------------------------------
 */
#ifndef StMuAnalysisMaker_hh
#define StMuAnalysisMaker_hh
//
//  Include files
#include "StMaker.h"
#include <string>
#include "TH3F.h"
#include "TTree.h"
#include "../../headers/CutClass.h"
#include "../../headers/makefile_toggles.h"
#include "../../headers/PhysMath.h"
#include "../../headers/HistogramUtilities.h"
#include "../../headers/Helix.h"
#ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
#include "../../headers/MattMcFDEvent.h"
#include "../../headers/MattMcFDTrack.h"
#else
#include "../../headers/MattMcEvent.h"
#include "../../headers/MattMcTrack.h"
#endif
//
//  Forward declarations
class StMuTrack;
class TFile;
class TH1D;
class TH2F;
class TH3F;
class StMuDst;
class StMuEvent;
class StMuTrack;
class StMuMcTrack;

#ifndef ST_NO_NAMESPACES
using std::string;
#endif








//
//  The class declaration. It innherits from StMaker.
class StMuAnalysisMaker : public StMaker{
public:
   StMuAnalysisMaker(const Char_t *name = "muDstEmbeddingAnalysis"); // constructor
   ~StMuAnalysisMaker();                                 // destructor

   void   Clear(Option_t *option = ""); // called after every event to cleanup
   Int_t  Init();                   // called once at the beginning of your job
   Int_t  Make();                   // invoked for every event
   Int_t  Finish();                 // called once at the end

   #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
   void set_mother_id(int a_id){m_mother_id = a_id;};
   bool isValidFeedDownID(int a_geantId);
   bool isValidDecayChain(vector<int> a_decay_chain); // note that phi is not implimented, if phi is in chain, just defaults to true
   #else
   bool isParticleOfInterest(int a_geantId);
   void setPart(string partCStr); 
   #endif
   bool isQualityEmbeddingMatch(StMuTrack* a_track); // things like QA Truth and flags

   void setOutputName(string a_name){m_FileName = a_name;};
   void setCutClass(CutClass* a_Cuts){ m_cuts = a_Cuts; };

   virtual const char *GetCVS() const{
      static const char cvs[] = "Tag $Name:  $ $Id: StMuAnalysisMaker.h,v 1.2 2014/08/06 11:43:31 jeromel Exp $ built " __DATE__ " " __TIME__ ;
      return cvs;
   }

private:

  CutClass* m_cuts;

  // data member
  int        m_EventCounter;  //!
  string     m_FileName;      //!
  TFile*     m_File;         //!

  TTree*       m_tree;
  #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
  int m_mother_id;
  MattMcFDEvent* m_event;
  MattMcFDTrack* m_track;
  #else
  int m_partIndex;
  int m_particleCharge;
  double m_mass;
  MattMcEvent* m_event;
  MattMcTrack* m_track;
  #endif

  double m_currentRefMult;
  int    m_currentRunId;
  int    m_nCentBins;


  int   m_numMakeCalls;
  int   m_numGoodMakeCalls;
  TH1I* m_numEventsPerDst;
  TH1I* m_numGoodEventsPerDst;

  TH1I* m_centEvents;
  TH1I* m_rejectedEventsEmb;
  TH1I* m_rejectedTracksEmb;
  TH1I* m_RefMult;
  TH1I* m_RefMultTrue;
  TH1I* m_numPrimaryDataTracks;
  TH1I* m_numMcTracks;
  TH1I* m_qaTruth;


  #ifdef _FEEDDOWN_EMBEDDING_TOGGLE_
  //QA Plots for Feed down Embedding
  TH2I* m_mother_decay_location;
  TH1I* m_mother_delta_decay_distance;
  TH1I* m_mother_delta_decay_distance_div_betagamma;
  TH2I* m_mother_delta_decay_location;

  #endif

  ClassDef(StMuAnalysisMaker, 0)
};
#endif
