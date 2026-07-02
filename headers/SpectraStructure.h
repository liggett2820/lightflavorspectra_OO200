

#ifndef _SPECTRA_STRUCTURE_
#define _SPECTRA_STRUCTURE_


#include <iostream>
#include <utility>

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TString.h>
#include <TMath.h>
#include <TFile.h>
#include <TMatrixDSym.h>

#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"
//#include "../headers/SpectraStructure.h"



/*
class TGraphErrors;
class TMultiGraph;
class TDirectory;
class TF1;
class TString;
class TFile;
*/
//class TMatrixDSym; // not a class, but a type def


using namespace std;


class SpectraStructure : public TObject {


 public:
  //#########################     CONSTRUCTOR     #############################
  void Construct();
  void Construct(int a_partIndex, int a_charge, int a_rapdidityIndex, double a_rapdidityValue, int a_centIndex, int a_lowCentralityPercent, int a_highCentralityPercent);
  //virtual ~SpectraStructure();

  //#########################     SPECIAL FUCTIONS     #############################
  void makeSpectraPDF(string a_imagePathAndName);
  void makeSpectraTXT(string a_txtFileName);

  //#########################     SETTERS     #############################

  void set_picoBinnerFileName(string a_name){m_picoBinnerFileName = a_name;};
  void set_picoBinnerRootName(string a_name){m_picoBinnerRootName = a_name;};
  void set_ZFitterFileName(string a_name){m_ZFitterFileName = a_name;};
  void set_ZFitterRootName(string a_name){m_ZFitterRootName = a_name;};
  void set_rawSpectraFileName(string a_name){m_rawSpectraFileName = a_name;};
  void set_rawSpectraName(string a_name){m_rawSpectraName = a_name;}; 
  void set_tpcEfficiencyFileName(string a_name){m_tpcEfficiencyFileName = a_name;};
  void set_tpcEfficiencyRootName(string a_name){m_tpcEfficiencyRootName = a_name;};
  void set_bToFEfficiencyFileName(string a_name){m_bToFEfficiencyFileName = a_name;};
  void set_bToFEfficiencyRootName(string a_name){m_bToFEfficiencyRootName = a_name;};
  void set_energyLossFileName(string a_name){m_energyLossFileName = a_name;};
  void set_energyLossRootName(string a_name){m_energyLossRootName = a_name;};
  void set_feedDownFileName(string a_name){m_feedDownFileName = a_name;};
  void set_feedDownRootName(string a_name){m_feedDownRootName = a_name;};
  void set_knockoutFileName(string a_name){m_knockoutFileName = a_name;};
  void set_knockoutRootName(string a_name){m_knockoutRootName = a_name;};

  /*
  // Detector index is -1 = total / 0=tpc / 1 = btof  / 2 = etof
  void addRawSpectraMtM0LowBinEdge(double a_edge, int a_detectorIndex = -1);
  void addRawSpectraMtM0HighBinEdge(double a_edge, int a_detectorIndex = -1);
  void addEnergyLossCorrectedMtM0LowBinEdge(double a_edge, int a_detectorIndex = -1);
  void addEnergyLossCorrectedMtM0HighBinEdge(double a_edge, int a_detectorIndex = -1);
  */

  void addTPCEffPoint(double a_mTm0, double a_efficiency, double a_effErr, bool a_integratedValue = true); // put false for point evaluation of function
  void addBTOFEffPoint(double a_mTm0, double a_efficiency, double a_effErr, bool a_integratedValue = true); 
  void addETOFEffPoint(double a_mTm0, double a_efficiency, double a_effErr, bool a_integratedValue = true);
  void addEnergyLossCorrectionPoint(double a_mTm0, double a_energyLoss, double a_error, bool a_integratedValue = true);
  void addFeedDownCorrectionPoint(double a_mTm0, double a_energyLoss, double a_error, bool a_integratedValue = true);


  void set_dNdy_Nom(double a_val){m_dNdy_Nom = a_val;};
  void set_dNdy_Nom_TotalErrLow(double a_val){m_dNdy_Nom_TotalErrLow = a_val;};
  void set_dNdy_Nom_TotalErrHigh(double a_val){m_dNdy_Nom_TotalErrHigh = a_val;};
  void set_dNdy_Nom_StatErr(double a_val){m_dNdy_Nom_StatErr = a_val;};
  void set_dNdy_Nom_SysErrHigh(double a_val){m_dNdy_Nom_SysErrHigh = a_val;};
  void set_dNdy_Nom_SysErrLow(double a_val){m_dNdy_Nom_SysErrLow = a_val;};

  void set_dNdy_Sys(double a_val){m_dNdy_Sys = a_val;};
  void set_dNdy_Sys_TotalErrLow(double a_val){m_dNdy_Sys_TotalErrLow = a_val;};
  void set_dNdy_Sys_TotalErrHigh(double a_val){m_dNdy_Sys_TotalErrHigh = a_val;};
  void set_dNdy_Sys_StatErr(double a_val){m_dNdy_Sys_StatErr = a_val;};
  void set_dNdy_Sys_SysErrHigh(double a_val){m_dNdy_Sys_SysErrHigh = a_val;};
  void set_dNdy_Sys_SysErrLow(double a_val){m_dNdy_Sys_SysErrLow = a_val;};
  
  // Clones the passed graphs
  void setRawTPCSpectra(TGraphErrors* a_spectra){m_raw_tpcSpectrum = *((TGraphErrors*) a_spectra->Clone());};
  void setRawBTOFSpectra(TGraphErrors* a_spectra){m_raw_btofSpectrum = *((TGraphErrors*) a_spectra->Clone());};
  void setRawETOFSpectra(TGraphErrors* a_spectra){m_raw_etofSpectrum = *((TGraphErrors*) a_spectra->Clone());};

  void setCorrectedTPCSpectra(TGraphErrors* a_spectra){m_tpcSpectrum = *((TGraphErrors*) a_spectra->Clone());};
  void setCorrectedBTOFSpectra(TGraphErrors* a_spectra){m_btofSpectrum = *((TGraphErrors*) a_spectra->Clone());};
  void setCorrectedETOFSpectra(TGraphErrors* a_spectra){m_etofSpectrum = *((TGraphErrors*) a_spectra->Clone());};
  void setCorrectedTotalSpectra(TGraphErrors* a_spectra){m_totalSpectrum = *((TGraphErrors*) a_spectra->Clone());};

  void set_tpcEff_originalGraph(TGraphAsymmErrors* a_graph){     m_tpcEff_originalGraph = *((TGraphAsymmErrors*) a_graph->Clone());};
  void set_energyLoss_originalGraph(TGraphAsymmErrors* a_graph){ m_energyLoss_originalGraph = *((TGraphAsymmErrors*) a_graph->Clone());};
  void set_bTOFEff_originalGraph(TGraphAsymmErrors* a_graph){    m_bTOFEff_originalGraph = *((TGraphAsymmErrors*) a_graph->Clone());};
  void set_eTOFEff_originalGraph(TGraphAsymmErrors* a_graph){    m_eTOFEff_originalGraph = *((TGraphAsymmErrors*) a_graph->Clone());};
  void set_feedDown_originalGraph(TGraphAsymmErrors* a_graph){   m_feedDown_originalGraph = *((TGraphAsymmErrors*) a_graph->Clone());};


  void set_spectra_nominalFit(TF1* a_funct){m_spectra_nominalFit = *((TF1*) a_funct->Clone());};
  void set_spectra_systematicFit(TF1* a_funct){m_spectra_systematicFit = *((TF1*) a_funct->Clone());};
  void set_tpcEff_nominalFit(TF1* a_funct){m_tpcEff_nominalFit = *((TF1*) a_funct->Clone());};
  void set_tpcEff_systematicFit(TF1* a_funct){m_tpcEff_systematicFit = *((TF1*) a_funct->Clone());};
  void set_bTOFEff_nominalFit(TF1* a_funct){m_bTOFEff_nominalFit = *((TF1*) a_funct->Clone());};
  void set_bTOFEff_systematicFit(TF1* a_funct){m_bTOFEff_systematicFit = *((TF1*) a_funct->Clone());};
  void set_energyLoss_nominalFit(TF1* a_funct){m_energyLoss_nominalFit = *((TF1*) a_funct->Clone());};
  void set_energyLoss_systematicFit(TF1* a_funct){m_energyLoss_systematicFit = *((TF1*) a_funct->Clone());};
  void set_feedDown_nominalFit(TF1* a_funct){m_feedDown_nominalFit = *((TF1*) a_funct->Clone());};
  void set_feedDown_systematicFit(TF1* a_funct){m_feedDown_systematicFit = *((TF1*) a_funct->Clone());};

  void set_spectra_nominalFit_CovMat(TMatrixDSym* a_matrix){m_spectra_nominalFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_spectra_systematicFit_CovMat(TMatrixDSym* a_matrix){m_spectra_systematicFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_tpcEff_nominalFit_CovMat(TMatrixDSym* a_matrix){m_tpcEff_nominalFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_tpcEff_systematicFit_CovMat(TMatrixDSym* a_matrix){m_tpcEff_systematicFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_bTOFEff_nominalFit_CovMat(TMatrixDSym* a_matrix){m_bTOFEff_nominalFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_bTOFEff_systematicFit_CovMat(TMatrixDSym* a_matrix){m_bTOFEff_systematicFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_energyLoss_nominalFit_CovMat(TMatrixDSym* a_matrix){m_energyLoss_nominalFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_energyLoss_systematicFit_CovMat(TMatrixDSym* a_matrix){m_energyLoss_systematicFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_feedDown_nominalFit_CovMat(TMatrixDSym* a_matrix){m_feedDown_nominalFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};
  void set_feedDown_systematicFit_CovMat(TMatrixDSym* a_matrix){m_feedDown_systematicFit_CovMat = *((TMatrixDSym*) a_matrix->Clone());};



  //#########################     GETTERS     #############################
  int getParticleIndex(){return m_partIndex;};
  int getParticleCharge(){return m_charge;};
  int getRapidityIndex(){return m_rapidityIndex;};
  double getRapidityValue(){return m_rapidityIndex;};
  int getCentralityIndex(){return m_centIndex;};



private:
  int m_partIndex;
  int m_charge;
  int m_rapidityIndex;
  double m_rapidityValue;
  int m_centIndex;
  int m_lowCentralityPercentageEdge; //example if 0-5%, this would be the 0
  int m_highCentralityPercentageEdge; //example if 0-5%, this would be the 5


  //#### Information on where things came from #####
  // FileName: full path to the .root file                          Example:  /Users/matthewharasty/RootStuff/thisisthefile.root
  // RootName: inside the root file, the full path and name of thing     Example:   /PionPlus/theSpecificGraph
  string m_picoBinnerFileName;
  string m_picoBinnerRootName;
  string m_ZFitterFileName;
  string m_ZFitterRootName;
  string m_rawSpectraFileName;
  string m_rawSpectraName; 
  string m_tpcEfficiencyFileName;
  string m_tpcEfficiencyRootName;
  string m_bToFEfficiencyFileName;
  string m_bToFEfficiencyRootName;
  string m_energyLossFileName;
  string m_energyLossRootName;
  string m_feedDownFileName;
  string m_feedDownRootName;
  string m_knockoutFileName;
  string m_knockoutRootName;

  vector<double> m_rawSpectraMtM0LowBinEdge; // the order is not the same index as the point, use the lookup funcition above
  vector<double> m_rawSpectraMtM0HighBinEdge; 
  vector<double> m_energyLossCorrectedMtM0LowBinEdge;  // the energy loss correction shifts the widths of each mT-m0 bin
  vector<double> m_energyLossCorrectedMtM0HighBinEdge;  // the energy loss correction shifts the widths of each mT-m0 bin

  //########################     dN/dy rapdidity density      ########################
  double m_dNdy_Nom;        
  double m_dNdy_Nom_TotalErrLow;
  double m_dNdy_Nom_TotalErrHigh;          
  double m_dNdy_Nom_StatErr; 
  double m_dNdy_Nom_SysErrHigh;
  double m_dNdy_Nom_SysErrLow;

  double m_dNdy_Sys;        
  double m_dNdy_Sys_TotalErrLow;
  double m_dNdy_Sys_TotalErrHigh;          
  double m_dNdy_Sys_StatErr; 
  double m_dNdy_Sys_SysErrHigh;
  double m_dNdy_Sys_SysErrLow;

  //########################     RAW SPECTRA      ########################
  TGraphErrors m_raw_tpcSpectrum;
  TGraphErrors m_raw_btofSpectrum;
  TGraphErrors m_raw_etofSpectrum;

  //########################     FULLY CORRECTED SPECTRA     ########################
  TGraphErrors m_totalSpectrum;
  TGraphErrors m_tpcSpectrum;
  TGraphErrors m_btofSpectrum;
  TGraphErrors m_etofSpectrum;

  //########################     FULLY CORRECTED YIELDS (dN/dmtm0*dy)     ########################
  TGraphErrors m_totalNonInvariantSpectrum;
  TGraphErrors m_tpcNonInvariantSpectrum;
  TGraphErrors m_btofNonInvariantSpectrum;
  TGraphErrors m_etofNonInvariantSpectrum;  


  //########################     CORRECTIONS     ########################
  TGraphAsymmErrors m_tpcEff_originalGraph;
  TGraphErrors m_tpcEff_nominal_integrated; // one point for every spectra point, ei: 2 if overlaped tpc and btof
  TGraphErrors m_tpcEff_systematic_integrated;
  TGraphErrors m_tpcEff_nominal_pointEval;
  TGraphErrors m_tpcEff_systematic_pointEval;

  TGraphAsymmErrors m_energyLoss_originalGraph;
  TGraphErrors m_energyLoss_nominal_integrated;
  TGraphErrors m_energyLoss_systematic_integrated;
  TGraphErrors m_energyLoss_nominal_pointEval;
  TGraphErrors m_energyLoss_systematic_pointEval;

  TGraphAsymmErrors m_bTOFEff_originalGraph;
  TGraphErrors m_bTOFEff_nominal_integrated;
  TGraphErrors m_bTOFEff_systematic_integrated;
  TGraphErrors m_bTOFEff_nominal_pointEval;
  TGraphErrors m_bTOFEff_systematic_pointEval;

  TGraphAsymmErrors m_eTOFEff_originalGraph;
  TGraphErrors m_eTOFEff_nominal_integrated;
  TGraphErrors m_eTOFEff_systematic_integrated;
  TGraphErrors m_eTOFEff_nominal_pointEval;
  TGraphErrors m_eTOFEff_systematic_pointEval;

  TGraphAsymmErrors m_feedDown_originalGraph;
  TGraphErrors m_feedDown_nominal_integrated;
  TGraphErrors m_feedDown_systematic_integrated;
  TGraphErrors m_feedDown_nominal_pointEval;
  TGraphErrors m_feedDown_systematic_pointEval;

  //FITS TO CORRECTIONS
  TF1 m_spectra_nominalFit;
  TF1 m_spectra_systematicFit;
  TF1 m_tpcEff_nominalFit;
  TF1 m_tpcEff_systematicFit;
  TF1 m_bTOFEff_nominalFit;
  TF1 m_bTOFEff_systematicFit;
  TF1 m_energyLoss_nominalFit;
  TF1 m_energyLoss_systematicFit;
  TF1 m_feedDown_nominalFit;
  TF1 m_feedDown_systematicFit;

  TMatrixDSym m_spectra_nominalFit_CovMat;
  TMatrixDSym m_spectra_systematicFit_CovMat;
  TMatrixDSym m_tpcEff_nominalFit_CovMat;
  TMatrixDSym m_tpcEff_systematicFit_CovMat;
  TMatrixDSym m_bTOFEff_nominalFit_CovMat;
  TMatrixDSym m_bTOFEff_systematicFit_CovMat;
  TMatrixDSym m_energyLoss_nominalFit_CovMat;
  TMatrixDSym m_energyLoss_systematicFit_CovMat;
  TMatrixDSym m_feedDown_nominalFit_CovMat;
  TMatrixDSym m_feedDown_systematicFit_CovMat;



  ClassDef(SpectraStructure,3); //number is version control. If something is changed in this file, increase the number by one
};


//ClassImp(SpectraStructure);

#endif
