#ifndef _MATT_MC_TRACK_CXX_
#define _MATT_MC_TRACK_CXX_

#include "../headers/MattMcTrack.h"

ClassImp(MattMcTrack);

MattMcTrack::MattMcTrack(){
  Clear();
}

void MattMcTrack::Clear(Option_t* a_option){
  m_pT_emb  = -9999;
  m_eta_emb = -999;
  m_phi_emb = -999;
  m_pT      = -999;
  m_eta     = -999;
  m_phi     = -999;
  m_gDCA_20000 = 65534;
  m_nHitsFit   = 0;
  m_nHitsDeDx  = 0;
  m_nHitsRatio_50000 = 0;
  m_btofLocalY_10000 = 65534;
  m_btofLocalZ_10000 = 65534;
}


void MattMcTrack::print(){
  cout << "MattMcTrack::print()" << endl;
  cout << "     pT_emb\t" << pT_emb() << endl;
  cout << "     eta_emb\t" << eta_emb() << endl;
  cout << "     phi_emb\t" << phi_emb() << endl;
  cout << "     pT  \t" << pT() << endl;
  cout << "     eta \t" << eta() << endl;
  cout << "     phi \t" << phi() << endl;
  cout << "     gDCA\t" << gDCA() << endl;
  cout << "     nHitsFit\t" << nHitsFit() << endl;
  cout << "     nHitsDeDx\t" << nHitsDeDx() << endl;
  cout << "     nHitsRatio\t" << nHitsRatio() << endl;  
  cout << "     btofLocalY\t" << btofLocalY() << endl;  
  cout << "     btofLocalZ\t" << btofLocalZ() << endl;

}

#endif
