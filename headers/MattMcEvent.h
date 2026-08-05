#ifndef _MATT_MC_EVENT_
#define _MATT_MC_EVENT_

#include "MattMcTrack.h"
#include <iostream>
using namespace std;

class MattMcEvent: public TObject {
public:

  MattMcEvent();
  //~MattMcEvent();
  void Clear(Option_t* a_option = "");
  //void setMaxTracks(int a_numTracks);
  void print();

  //###########  GETTERS ##############
  float refMult()          {return m_refMult;}; // includes the simulation tracks
  float refMult_data_only(){return m_refMult;}; // does not include the simulation tracks
  int   numPrimaryTracks() {return ((int) m_numPrimaryTracks);};
  float vertex_X()         {return ((float) m_vertex_X_1000)/1000.0;};
  float vertex_Y()         {return ((float) m_vertex_Y_1000)/1000.0;};
  float vertex_Z()         {return ((float) m_vertex_Z_100)/100.0;};
  int   numTofMatch()      {return ((int) m_numTofMatch);};
  int   numTofT0()         {return ((int) m_numTofT0);};
  int   centIndex()        {return ((int) m_centIndex);};
  int   numMcTracks()      {return m_tracks->GetEntriesFast();};

  //##########  SETTERS ##############
  void  set_refMult          (float a_val){m_refMult           = a_val;};
  void  set_refMult_data_only(float a_val){m_refMult_data_only = a_val;};
  void  set_numPrimaryTracks (int   a_val){if(a_val < 65534)        m_numPrimaryTracks  = (unsigned short int) (a_val);};
  void  set_vertex_X         (float a_val){if(fabs(a_val) < 32.766) m_vertex_X_1000     = (signed short int) (a_val*1000.0);}; //-32 -> 32 cm
  void  set_vertex_Y         (float a_val){if(fabs(a_val) < 32.766) m_vertex_Y_1000     = (signed short int) (a_val*1000.0);}; //-32 -> 32 cm
  void  set_vertex_Z         (float a_val){if(fabs(a_val) < 327.66) m_vertex_Z_100      = (signed short int) (a_val*100.0);};  //-320 -> 320 cm
  void  set_numTofMatch      (int   a_val){if(0 <= a_val && a_val <= 254) m_numTofMatch = (unsigned char) (a_val);}; //0-255
  void  set_numTofT0         (int   a_val){if(0 <= a_val && a_val <= 254) m_numTofT0    = (unsigned char) (a_val);}; //0-255
  void  set_centIndex        (int   a_val){if(0 <= a_val && a_val <= 254) m_centIndex   = (unsigned char) (a_val);}; //0-255
  void  addTrack(int a_trackIndex, MattMcTrack* a_track);

  MattMcTrack* getTrack(int a_index){return (MattMcTrack*) m_tracks->UncheckedAt(a_index);}
private:

  float              m_refMult;
  float              m_refMult_data_only;
  unsigned short int m_numPrimaryTracks;
  signed short int   m_vertex_X_1000; //-32 -> 32 cm
  signed short int   m_vertex_Y_1000; //-32 -> 32 cm
  signed short int   m_vertex_Z_100;  //-320 -> 320 cm
  unsigned char      m_numTofMatch;   //0-255
  unsigned char      m_numTofT0;      //0-255
  unsigned char      m_centIndex;     //0-255
  TClonesArray*      m_tracks;        //-> 

  ClassDef(MattMcEvent,1);
};


#endif
