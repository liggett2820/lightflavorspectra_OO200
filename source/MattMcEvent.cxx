#ifndef _MATT_MC_EVENT_CXX_
#define _MATT_MC_EVENT_CXX_

#include "../headers/MattMcEvent.h"


ClassImp(MattMcEvent);

MattMcEvent::MattMcEvent(){
  m_tracks = new TClonesArray("MattMcTrack");
  Clear();
}

/*void MattMcEvent::setMaxTracks(int a_numTracks){
  m_tracks = new TClonesArray("MattMcTrack", a_numTracks);
}*/

void MattMcEvent::Clear(Option_t* a_option){
  m_refMult = -999;
  m_refMult_data_only = -999;
  m_numPrimaryTracks = 0;
  m_vertex_X_1000 = -32766; //-32 -> 32 cm
  m_vertex_Y_1000 = -32766; //-32 -> 32 cm
  m_vertex_Z_100 = -32766;  //-320 -> 320 cm
  m_numTofMatch = 0;   //0-255
  m_numTofT0 = 0;      //0-255 
  m_centIndex = 250;   //0-255
  if(m_tracks) m_tracks->Clear("C");
}

void MattMcEvent::print(){
  cout << "MattMcEvent::print()" << endl;
  cout << "     refMult\t" << refMult() << endl;          
  cout << "     refMult_data_only\t" << refMult_data_only() << endl;
  cout << "     numPrimaryTracks\t" << numPrimaryTracks() << endl; 
  cout << "     vertex_X\t" << vertex_X() << endl;         
  cout << "     vertex_Y\t" << vertex_Y() << endl;         
  cout << "     vertex_Z\t" << vertex_Z() << endl;         
  cout << "     numTofMatch\t" << numTofMatch() << endl;      
  cout << "     numTofT0\t" << numTofT0() << endl;     
  cout << "     centIndex\t" << centIndex() << endl;
  cout << "     numMuTracks\t" << m_tracks->GetEntriesFast() << endl;
}


void  MattMcEvent::addTrack(int a_trackIndex, MattMcTrack* a_track){
  if(!a_track){
    cout << "ERROR: MattMcEvent::addTrack  you supplied a trackIndex = " << a_trackIndex << " and ptr= " << a_track << endl;
    return;
  }
  MattMcTrack* track = (MattMcTrack*) m_tracks->ConstructedAt(a_trackIndex);
  track->copy(a_track); 
}







#endif
