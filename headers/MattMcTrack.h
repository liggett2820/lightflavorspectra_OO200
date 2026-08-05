#ifndef _MATT_MC_TRACK_
#define _MATT_MC_TRACK_

#include "TClonesArray.h"
#include "PhysMath.h"
#include <iostream>
using namespace std;

// This class contains both simulation (MC) information and reconstructed information (Reco) from STAR embedding data
class MattMcTrack: public TObject {
public:

  MattMcTrack();
  //~MattMcTrack();
  void Clear(Option_t* a_option = "");

  void copy(MattMcTrack* a_track){
    m_pT_emb  = a_track->m_pT_emb;
    m_eta_emb = a_track->m_eta_emb;
    m_phi_emb = a_track->m_phi_emb;
    m_pT = a_track->m_pT;
    m_eta = a_track->m_eta;
    m_phi = a_track->m_phi;
    m_gDCA_20000 = a_track->m_gDCA_20000;
    m_nHitsFit = a_track->m_nHitsFit;
    m_nHitsDeDx = a_track->m_nHitsDeDx;
    m_nHitsRatio_50000 = a_track->m_nHitsRatio_50000;
    m_btofLocalY_10000 = a_track->m_btofLocalY_10000;
    m_btofLocalZ_10000 = a_track->m_btofLocalZ_10000;
  }

  void print();
  bool wasReconstructed(){return (pT() > 0.0);};


  //##########  GETTERS  ###################
  //Simulated Variables
  float pT_emb()    {return m_pT_emb;};
  float eta_emb()   {return m_eta_emb;};
  float phi_emb()   {return m_phi_emb;};

  //Reconstructed Variables
  float pT()        {return m_pT;};
  float eta()       {return m_eta;};
  float phi()       {return m_phi;};
  float gDCA()      {return ((float) m_gDCA_20000)/20000.0;};
  int   nHitsFit()  {return ((int) m_nHitsFit);};
  int   nHitsDeDx() {return ((int) m_nHitsDeDx);};
  float nHitsRatio(){return ((float) m_nHitsRatio_50000)/50000.0;};
  float btofLocalY(){return ((float) m_btofLocalY_10000)/10000.0;};
  float btofLocalZ(){return ((float) m_btofLocalZ_10000)/10000.0;};


  // Calculated Variables
  float pTot_emb()                {return sqrt(pow(pT_emb(),2) + pow(pZ_emb(),2));};
  float pX_emb()                  {return pT_emb()*cos(phi_emb());};
  float pY_emb()                  {return pT_emb()*sin(phi_emb());};
  float pZ_emb()                  {return PhysMath::pZ_from_eta_pT(eta_emb(),pT_emb());};
  float mTm0_emb(float a_mass)    {return PhysMath::mTm0(a_mass,pT_emb());};
  float rapidity_emb(float a_mass){return PhysMath::rapFromEtaPt(a_mass,eta_emb(),pT_emb());};

  float pTot()                    {return sqrt(pow(pT(),2) + pow(pZ(),2));};
  float pX()                      {return pT()*cos(phi());};
  float pY()                      {return pT()*sin(phi());};
  float pZ()                      {return PhysMath::pZ_from_eta_pT(eta(),pT());};
  float mTm0(float a_mass)        {return PhysMath::mTm0(a_mass,pT());};
  float rapidity(float a_mass)    {return PhysMath::rapFromEtaPt(a_mass,eta(),pT());};



  //########   SETTERS  ####################
  void set_pT_emb    (float a_val){m_pT_emb = a_val;};
  void set_eta_emb   (float a_val){m_eta_emb = a_val;};
  void set_phi_emb   (float a_val){m_phi_emb = a_val;};
  void set_pT        (float a_val){m_pT = a_val;};

  void set_eta       (float a_val){m_eta = a_val;};
  void set_phi       (float a_val){m_phi = a_val;};
  void set_gDCA      (float a_val){if(0.0 < a_val && a_val < 3.2766) m_gDCA_20000 = (unsigned short int) (a_val*20000.0);};
  void set_nHitsFit  (int   a_val){if(0 <=a_val && a_val <= 254)     m_nHitsFit  = (unsigned char) (a_val);};
  void set_nHitsDeDx (int   a_val){if(0 <=a_val && a_val <= 254)     m_nHitsDeDx = (unsigned char) (a_val);};
  void set_nHitsRatio(float a_val){if(0.0 <= a_val && a_val < 1.01)  m_nHitsRatio_50000 = (unsigned short int) (a_val*50000.0);};
  void set_btofLocalY(float a_val){if(0.0 < a_val && a_val < 6.5534) m_btofLocalY_10000 = (unsigned short int) (a_val*10000.0);};
  void set_btofLocalZ(float a_val){if(0.0 < a_val && a_val < 6.5534) m_btofLocalZ_10000 = (unsigned short int) (a_val*10000.0);};

//private:
  float              m_pT_emb;
  float              m_eta_emb;
  float              m_phi_emb;
  float              m_pT;
  float              m_eta;
  float              m_phi;
  unsigned short int m_gDCA_20000;
  unsigned char      m_nHitsFit;
  unsigned char      m_nHitsDeDx;
  unsigned short int m_nHitsRatio_50000;
  unsigned short int m_btofLocalY_10000;
  unsigned short int m_btofLocalZ_10000;

  ClassDef(MattMcTrack,1); // last number is version control
};

#endif
