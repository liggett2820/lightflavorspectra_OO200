#ifndef _CALIBRATION_CLASS_CXX_
#define _CALIBRATION_CLASS_CXX_

#include "../headers/CalibrationClass.h"

CalibrationClass::CalibrationClass(){
  #ifdef _CALIBRATIONCLASS_DEBUG_
    cout << "Call to CalibrationClass::CalibrationClass()" << endl;
  #endif
  m_partInfo           = NULL;
  m_doSigmaAdjustment  = false;
  m_usePionCalibForKaon = false;
  m_particleOfInterest = -1;

  m_masses[0] = 0.13957;    //Pion
  m_masses[1] = 0.49367;    //Kaon
  m_masses[2] = 0.93827;    //Proton
  m_masses[3] = 0.00051099; //Electron
  m_masses[4] = 1.87561;    //Deuteron
  m_masses[5] = 2.80925;    //Triton
  m_masses[6] = 2.8094;     //Helion
  m_masses[7] = 3.72742;    //Alpha
  m_masses[8] = 0.105658;   //Muon

  m_particle_loaded[0] = false;
  m_particle_loaded[1] = false;
  m_particle_loaded[2] = false;
  m_particle_loaded[3] = false;
  m_particle_loaded[4] = false;
  m_particle_loaded[5] = false;
  m_particle_loaded[6] = false;
  m_particle_loaded[7] = false;
  m_particle_loaded[8] = false;

  m_dEdx_invBetaGammaSqr_profile  = NULL;
  m_dEdx_invBetaGammaSqr_parabola = NULL;


  m_vZ_calib_struct  = NULL;
  m_eta_calib_struct = NULL;
  for(int partIndex = 0; partIndex < 9; partIndex++){
    m_minMaxVzEtaIndex[partIndex][0] = 999;
    m_minMaxVzEtaIndex[partIndex][1] = -999;
    m_minMaxVzEtaIndex[partIndex][2] = 999;
    m_minMaxVzEtaIndex[partIndex][3] = -999;
    for(int vZ_index = 0; vZ_index < g_CalibClass_MaxNumVzBins; vZ_index++){
      for(int eta_index = 0; eta_index < g_CalibClass_MaxNumEtaBins; eta_index++){
        m_mean_calibrationFunctions[partIndex][vZ_index][eta_index] = NULL;
        m_meanPlusSigma_calibrationFunctions[partIndex][vZ_index][eta_index] = NULL;
        m_sigma_calibrationFunctions[partIndex][vZ_index][eta_index]  = NULL;
      }
    }
  }

  #ifdef _CALIBRATIONCLASS_DEBUG_
    cout << "End of CalibrationClass::CalibrationClass()" << endl;
  #endif

}


void CalibrationClass::loadDeDxCalibration(string a_filename, bool a_byParamHistos){
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  m_vZ_calib_struct  = (TH1I*) inFile->Get("vZ_calib_struct");
  m_eta_calib_struct = (TH1I*) inFile->Get("eta_calib_struct");

  for(int partIndex = 0; partIndex < 9; partIndex++){
    for(int vZ_index = 0; vZ_index < g_CalibClass_MaxNumVzBins; vZ_index++){
      for(int eta_index = 0; eta_index < g_CalibClass_MaxNumEtaBins; eta_index++){

        if(a_byParamHistos){
          TH1D* mean_param_histo = (TH1D*) inFile->Get(Form(
                "calibration_functions/mean_calibration_fit_%s_vZIndex%02d_etaIndex%02d_params",m_partInfo->GetParticleName(partIndex).Data(),vZ_index,eta_index));
          if(!mean_param_histo) continue;
          m_mean_calibrationFunctions[partIndex][vZ_index][eta_index] = PhysMath::getLnBichselFunctionBetaGamma(
                Form("mean_calibration_fit_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(partIndex).Data(),vZ_index,eta_index));
          m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]->SetParameter(0,mean_param_histo->GetBinContent(1));
          m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]->SetParameter(1,mean_param_histo->GetBinContent(2));
          m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]->SetParameter(2,mean_param_histo->GetBinContent(3));
          m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]->SetParameter(3,mean_param_histo->GetBinContent(4));
          m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]->SetParameter(4,mean_param_histo->GetBinContent(5));

        }else{
          m_mean_calibrationFunctions[partIndex][vZ_index][eta_index] = (TF1*) inFile->Get(Form(
        	"calibration_functions/mean_calibration_fit_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(partIndex).Data(),vZ_index,eta_index));
          m_sigma_calibrationFunctions[partIndex][vZ_index][eta_index] = (TF1*) inFile->Get(Form(
        	"calibration_functions/sigma_calibration_fit_%s_vZIndex%02d_etaIndex%02d",m_partInfo->GetParticleName(partIndex).Data(),vZ_index,eta_index));
        }
        if(m_mean_calibrationFunctions[partIndex][vZ_index][eta_index]){
          m_particle_loaded[partIndex] = true;
          if(vZ_index  < m_minMaxVzEtaIndex[partIndex][0]) m_minMaxVzEtaIndex[partIndex][0] = vZ_index;
          if(vZ_index  > m_minMaxVzEtaIndex[partIndex][1]) m_minMaxVzEtaIndex[partIndex][1] = vZ_index;
          if(eta_index < m_minMaxVzEtaIndex[partIndex][2]) m_minMaxVzEtaIndex[partIndex][2] = eta_index;
          if(eta_index > m_minMaxVzEtaIndex[partIndex][3]) m_minMaxVzEtaIndex[partIndex][3] = eta_index;          
        }
      }
    }
    #ifdef  _CALIBRATIONCLASS_DEBUG_
      cout << "PartIndex " << partIndex << endl;
      cout << "     Vz Index = [" <<   m_minMaxVzEtaIndex[partIndex][0] << ", " << m_minMaxVzEtaIndex[partIndex][1] << "]" << endl;
      cout << "     Eta Index = [" <<  m_minMaxVzEtaIndex[partIndex][2] << ", " << m_minMaxVzEtaIndex[partIndex][3] << "]" << endl;
    #endif

  }
}


void CalibrationClass::loadDeDxCalibrationProfile(string a_filename){
  TFile* inFile = new TFile(a_filename.c_str(),"READ");
  if(!inFile->IsOpen()){
    cout << "ERROR: CalibrationClass::loadDeDxCalibrationProfile    Can't open file " << a_filename << endl;
  }
  m_dEdx_invBetaGammaSqr_profile  = (TProfile*) inFile->Get("calibProfile_invBetaGammaSqr");
  m_dEdx_invBetaGammaSqr_parabola = (TF1*) inFile->Get("parab2");

}




void CalibrationClass::setParticleInfo(ParticleInfo* a_info){
	m_partInfo = a_info;
}
void CalibrationClass::setDoSigmaAdjustment(bool a_toggle){
	m_doSigmaAdjustment = a_toggle;
}
void CalibrationClass::setUsePionCalibForKaon(bool a_toggle){
  m_usePionCalibForKaon = a_toggle;
}

void CalibrationClass::setParticleOfInterest(int a_partIndex){
  if(!m_particle_loaded[a_partIndex]){
    cout << "WARNING: CalibrationClass::setParticleOfInterest(" << a_partIndex << ")   Particle Not Loaded. Using Pion Calibration Instead" << endl;
    m_particleOfInterest = 0;
    return;
  }

  m_particleOfInterest = a_partIndex;
  if(a_partIndex == 1 && m_usePionCalibForKaon){
    m_particleOfInterest = 0;
  }
}

double CalibrationClass::getCalibratedDeDx(double a_Vz, double a_eta, double a_momentum){
  #ifdef  _CALIBRATIONCLASS_DEBUG_
    cout << "CalibrationClass::getCalibratedDeDx(" << a_Vz << "," << a_eta << "," << a_momentum << ")" << endl;
  #endif


  int vZ_bin_1 = m_vZ_calib_struct->FindBin(a_Vz);
  //if(vZ_bin_1 < 1) vZ_bin_1 = 1;
  //else if(vZ_bin_1 > m_vZ_calib_struct->GetNbinsX()) vZ_bin_1 = m_vZ_calib_struct->GetNbinsX();
  if(vZ_bin_1 - 1 < m_minMaxVzEtaIndex[m_particleOfInterest][0]) vZ_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][0] + 1;
  if(vZ_bin_1 - 1 > m_minMaxVzEtaIndex[m_particleOfInterest][1]) vZ_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][1] + 1;

  int vZ_bin_2 = -1;
  double vZ_1 = m_vZ_calib_struct->GetBinCenter(vZ_bin_1);
  double vZ_2 = 0;
  if(a_Vz < m_vZ_calib_struct->GetBinCenter(vZ_bin_1)){
  	vZ_bin_2 = vZ_bin_1 - 1;
  	vZ_2 = m_vZ_calib_struct->GetBinCenter(vZ_bin_2);
  }else{
  	vZ_bin_2 = vZ_bin_1 + 1;
        vZ_2 = m_vZ_calib_struct->GetBinCenter(vZ_bin_2);
  }

  int eta_bin_1 = m_eta_calib_struct->FindBin(a_eta);
  if(eta_bin_1 - 1 < m_minMaxVzEtaIndex[m_particleOfInterest][2]) eta_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][2] + 1;
  if(eta_bin_1 - 1 > m_minMaxVzEtaIndex[m_particleOfInterest][3]) eta_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][3] + 1;

  int eta_bin_2 = -1;
  double eta_1  = m_eta_calib_struct->GetBinCenter(eta_bin_1);
  double eta_2  = 0;
  if(a_eta < m_eta_calib_struct->GetBinCenter(eta_bin_1)){
  	eta_bin_2 = eta_bin_1 - 1;
  	eta_2     = m_eta_calib_struct->GetBinCenter(eta_bin_2);
  }else{
  	eta_bin_2 = eta_bin_1 + 1;
        eta_2     = m_eta_calib_struct->GetBinCenter(eta_bin_2);
  }


  //Protection for edges of histograms
  if(vZ_bin_2 < 1 || vZ_bin_2 > m_vZ_calib_struct->GetNbinsX()){
  	vZ_bin_2 = vZ_bin_1;
  	vZ_2     = vZ_1;
  }
  if(eta_bin_2 < 1 || eta_bin_2 > m_eta_calib_struct->GetNbinsX()){
  	eta_bin_2 = eta_bin_1;
  	eta_2     = eta_1;
  }

  #ifdef  _CALIBRATIONCLASS_DEBUG_
    cout << "    Trying Vz Bins: " << vZ_bin_1 << "  " << vZ_bin_2 << "  min Index: " << m_minMaxVzEtaIndex[m_particleOfInterest][0] << "   max Index " << m_minMaxVzEtaIndex[m_particleOfInterest][1]<< endl;
    cout << "    Trying Eta Bins:" << eta_bin_1 << "  " << eta_bin_2 << "  minIndex: " << m_minMaxVzEtaIndex[m_particleOfInterest][2] << "   maxIndex " << m_minMaxVzEtaIndex[m_particleOfInterest][3] << endl;
  #endif


  double betaGamma = a_momentum/m_masses[m_particleOfInterest];
  if(!(m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_1 - 1]
     && m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_1 - 1]
     && m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_2 - 1]
     && m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_2 - 1])) return -999;
  double d1 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  double d2 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  double d3 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  double d4 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  #ifdef  _CALIBRATIONCLASS_DEBUG_
    cout << "      Particle Of Interest: " << m_particleOfInterest << endl;
    cout << "      BetaGamma: " << betaGamma << endl;
    cout << "      Points for 2D Interpolation: " << d1 << "  " << d2 << "  " << d3 << "  " << d4 << endl;
    cout << "      Name Check: " << endl;
    cout << "         " << m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_1 - 1]->GetName() << endl;
    cout << "         " << m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_1 - 1]->GetName() << endl;
    cout << "         " << m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_2 - 1]->GetName() << endl;
    cout << "         " << m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_2 - 1]->GetName() << endl;

  #endif

  return PhysMath::interpolateSquare2D(a_Vz,a_eta, vZ_1, eta_1, vZ_2,eta_2, d1, d2, d3, d4);


}


double CalibrationClass::getCalibratedNSigma(double a_Vz, double a_eta, double a_momentum, double a_dEdx){
  int vZ_bin_1 = m_vZ_calib_struct->FindBin(a_Vz);
  if(vZ_bin_1 - 1 < m_minMaxVzEtaIndex[m_particleOfInterest][0]) vZ_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][0] + 1;
  if(vZ_bin_1 - 1 > m_minMaxVzEtaIndex[m_particleOfInterest][1]) vZ_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][1] + 1;


  int vZ_bin_2 = -1;
  double vZ_1 = m_vZ_calib_struct->GetBinCenter(vZ_bin_1);
  double vZ_2 = 0;
  if(a_Vz < m_vZ_calib_struct->GetBinCenter(vZ_bin_1)){
  	vZ_bin_2 = vZ_bin_1 - 1;
  	vZ_2 = m_vZ_calib_struct->GetBinCenter(vZ_bin_2);
  }else{
  	vZ_bin_2 = vZ_bin_1 + 1;
    vZ_2 = m_vZ_calib_struct->GetBinCenter(vZ_bin_2);
  }

  int eta_bin_1 = m_eta_calib_struct->FindBin(a_eta);
  if(eta_bin_1 - 1 < m_minMaxVzEtaIndex[m_particleOfInterest][2]) eta_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][2] + 1;
  if(eta_bin_1 - 1 > m_minMaxVzEtaIndex[m_particleOfInterest][3]) eta_bin_1 = m_minMaxVzEtaIndex[m_particleOfInterest][3] + 1;
  int eta_bin_2 = -1;
  double eta_1  = m_eta_calib_struct->GetBinCenter(eta_bin_1);
  double eta_2  = 0;
  if(a_eta < m_eta_calib_struct->GetBinCenter(eta_bin_1)){
  	eta_bin_2 = eta_bin_1 - 1;
  	eta_2     = m_eta_calib_struct->GetBinCenter(eta_bin_2);
  }else{
  	eta_bin_2 = eta_bin_1 + 1;
    eta_2     = m_eta_calib_struct->GetBinCenter(eta_bin_2);
  }


  //Protection for edges of histograms
  if(vZ_bin_2 < 1 || vZ_bin_2 > m_vZ_calib_struct->GetNbinsX()){
  	vZ_bin_2 = vZ_bin_1;
  	vZ_2     = vZ_1;
  }
  if(eta_bin_2 < 1 || eta_bin_2 > m_eta_calib_struct->GetNbinsX()){
  	eta_bin_2 = eta_bin_1;
  	eta_2     = eta_1;
  }

  double betaGamma = a_momentum/m_masses[m_particleOfInterest];
  double d1 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  double d2 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  double d3 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  double d4 = m_mean_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  double ln_dE_dx = PhysMath::interpolateSquare2D(a_Vz,a_eta, vZ_1, eta_1, vZ_2,eta_2, d1, d2, d3, d4);

  double t1 = m_meanPlusSigma_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  double t2 = m_meanPlusSigma_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  double t3 = m_meanPlusSigma_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  double t4 = m_meanPlusSigma_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  double ln_dE_dx_Plus_Sigma = PhysMath::interpolateSquare2D(a_Vz,a_eta, vZ_1, eta_1, vZ_2,eta_2, t1, t2, t3, t4);

   

  return (log(a_dEdx) - ln_dE_dx)/fabs(ln_dE_dx_Plus_Sigma - ln_dE_dx);

  //double s1 = m_sigma_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  //double s2 = m_sigma_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_1 - 1]->Eval(betaGamma);
  //double s3 = m_sigma_calibrationFunctions[m_particleOfInterest][vZ_bin_1 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  //double s4 = m_sigma_calibrationFunctions[m_particleOfInterest][vZ_bin_2 - 1][eta_bin_2 - 1]->Eval(betaGamma);
  //double sigma = PhysMath::interpolateSquare2D(a_Vz,a_eta, vZ_1, eta_1, vZ_2,eta_2, s1, s2, s3, s4);
  //return log(a_dEdx/dE_dx)/sigma;

}

double CalibrationClass::getCalibratedZTPC(double a_Vz, double a_eta, double a_momentum, double a_dEdx){
  #ifdef  _CALIBRATIONCLASS_DEBUG_
    cout << "CalibrationClass::getCalibratedZTPC(" << a_Vz << "," << a_eta << "," << a_momentum << "," << a_dEdx << ")" << endl;
    cout << "   without Sigma Adjustment: " << log(a_dEdx) - getCalibratedDeDx(a_Vz,a_eta,a_momentum) << endl;
  #endif


  if(m_doSigmaAdjustment){
    return getCalibratedNSigma(a_Vz,a_eta,a_momentum,a_dEdx)*0.067658648; // number is to ln(1.07) to convert nSigma to dE/dx 
  }else{
    return log(a_dEdx) - getCalibratedDeDx(a_Vz,a_eta,a_momentum);
  }
}

double CalibrationClass::getDeDxFromProfileParab(double a_momentum){
  double rangeLow  = 0;
  double rangeHigh = 0;
  m_dEdx_invBetaGammaSqr_parabola->GetRange(rangeLow,rangeHigh);
  double invBetaGammaSqr = pow(m_masses[m_particleOfInterest]/a_momentum,2.0);
  if(rangeLow < invBetaGammaSqr){
    return m_dEdx_invBetaGammaSqr_parabola->Eval(invBetaGammaSqr);
  }else{
    int lowBin = m_dEdx_invBetaGammaSqr_profile->FindBin(invBetaGammaSqr);
    return PhysMath::interpolate(invBetaGammaSqr, m_dEdx_invBetaGammaSqr_profile->GetBinCenter(lowBin),   m_dEdx_invBetaGammaSqr_profile->GetBinContent(lowBin),
                                           m_dEdx_invBetaGammaSqr_profile->GetBinCenter(lowBin+1), m_dEdx_invBetaGammaSqr_profile->GetBinContent(lowBin+1));
  }
}





#endif
