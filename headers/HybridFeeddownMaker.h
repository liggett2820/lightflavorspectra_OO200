#ifndef _HYBRID_FEEDDOWN_MAKER_
#define _HYBRID_FEEDDOWN_MAKER_

#include "TH2D.h"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TVirtualFitter.h"
#include "TGraphAsymmErrors.h"
//#include "TCanvas.h"
#include "TNtuple.h"
#include "TString.h"
#include "TPaveText.h"
#include "TRandom3.h"

#include <string>
#include <cstring>
#include <iostream>

#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"
#include "PhysMath.h"
#include "HistogramUtilities.h"
#include "Helix.h"

#ifdef _MULTI_THREADING_
  #include <mutex>
  #include <thread>
#endif

using namespace std;

//This class will take the mid-rapidity measurements lambdas/Sigma/Omega... and the urqmd spectra file to produce feed-down correction graphs
//The pT spectra of the UrQMD simulation are stretched to fit the data, so that the spectra away from mid-rapidity can be used as well.


class HybridFeeddownMaker{
public:

  // FEED DOWN ID NAMES: 0:Lambda 1:Antilambda 2:SigmaPlus 3:AntiSigmaMinus 4:XiMinus 5:AntiXiPlus 6:XiZero 7:OmegaMinus 8:AntiOmegaPlus  9: SigmaZero 10:ANTISIGMA_0  11: KAON_0_S 12: ETA

  //####  Constructor to initialize member variables
  HybridFeeddownMaker();

  //##### Load the output of ./macros/RunPostProcessingUrQMD (contains dN/dy and invariant spectra by centrality)
  void loadUrQMDSpectraFile(string a_filename, bool a_oldNames = false);

  //#### Load the mid-rapidity invariant pT spectra from measurments (this will be fit later)
  void loadDataSpectra(int a_feeddownParticleId,  // 0:Lambda 1:Antilambda 2:SigmaPlus 3:AntiSigmaMinus 4:XiMinus 5:AntiXiPlus 6:XiZero 7:OmegaMinus 8:AntiOmegaPlus  9: SigmaZero 10:ANTISIGMA_0  11: KAON_0_S 12: ETA
                       int a_centIndex, // Centrality Index: 0 is most central
                       TGraphAsymmErrors* a_graph, // an invariant spectra from measurments
                       bool a_mTm0Space = false // false if the spectra is in pT space, true if the spectra is in mT-m0 space
                       );

  //#### Load the scaling and stretching parameters from a different file ( specifically for 14.5 with 19.6 GeV parameters)
  void loadScaleAndStretchParameters(string a_filename);

  //#### Scales and Stretches the UrQMD simulation invariant spectra to the measured data. Lambda measurments include sigma zero contributions, so they must be included (urqmd doesn't decay them)
  void findSimToDataScaling(int a_feeddownParticleId, // 0:Lambda 1:Antilambda 2:SigmaPlus 3:AntiSigmaMinus 4:XiMinus 5:AntiXiPlus 6:XiZero 7:OmegaMinus 8:AntiOmegaPlus  9: SigmaZero
                            int a_centIndex, //  Centrality Index: 0 is most central
                            vector< int > a_centralityBinsToUse, // 0,5,10,20,30,40,50,60,70,80 are normal cuts. If the measurment is a wider centrality bin, add more indecies to this vector
                            double a_lowRapidityOfData, // Low Rapdiity Window Cut for Meaurment Spectra
                            double a_highRapidityOfData, // High Rapidity WIndow Cut for Measurement Spectra
                            vector<int> a_extraFeedDownIdsToInclude, // If the measurement is a contribution of several particles, add them here (example: Lambda + Sigma Zero);
                            bool a_isMostPeripheralCentralityMeasured = false); // put this to true a the highest centrality to propagate scaling to other centralities not measured

  //#### Fit the measurment with a blast wave and use this fit to resample decays later. If you use findSimToDataScaling, you don't need this. It is just another option
  #ifndef _BEFORE_ROOT_6p24_
  void fitMeasurement(int a_feeddownParticleId,  // 0:Lambda 1:Antilambda 2:SigmaPlus 3:AntiSigmaMinus 4:XiMinus 5:AntiXiPlus 6:XiZero 7:OmegaMinus 8:AntiOmegaPlus  9: SigmaZero
                      int a_centIndex //  Centrality Index: 0 is most central
                      );
  #endif

  //#### Makes the data structure to fill : range and bins for y and mT-m0 response functions (where the protons go after lambda decay)
  void setResponseStructure(int a_numRapBins,
                            double a_lowRapRange,
                            double a_highRapRange,
                            double a_numMtM0Bins,
                            double a_lowMtM0Range,
                            double a_highMtM0Range);

  void setDoNTuple(bool a_toggle = true){m_doNTuple = a_toggle;};
  void setImageDir(string a_imageDir){m_imageDir = a_imageDir;};
  void setDCACut(bool a_use_cut = true, double a_cut = 1.0){m_use_DCA_cut = a_use_cut; m_DCA_Cut = a_cut;};

  //#### After you scale or fit the UrQMD Simulation, use that spectra to feed-down into proton, pions, kaons... This will give you the expected spectra of the daughters
  void generateResponseFunctions(int a_feeddownParticleId, // 0:Lambda 1:Antilambda 2:SigmaPlus 3:AntiSigmaMinus 4:XiMinus 5:AntiXiPlus 6:XiZero 7:OmegaMinus 8:AntiOmegaPlus  9: SigmaZero
                                 int a_centIndex, //  Centrality Index: 0 is most central
                                 int a_numTrials = 100000,
                                 double a_feeddownSpectraIndexToUse = -1, // use a different particles spectra as a proxy for this particle (transforms mT-m0 between each particle)
                                 double a_extraFactor = 1.0, // if you want to use a different particle, you probably want to scale that spectra by a factor (sigma/lambda ratio)
                                 bool a_useFitFunct = false // if you called fitMeasurement, you can use the fit function instead of the histogram to resample decays
                                 );

  //#### After generating the response functions (the spectra of the daughters), combine all of the differetent decays into one histogram
  void combineFeedDownSpectra();

  //#### After combining the feed down spectra, fit them with a blast wave function
  //  This is now handled in EfficiencyFitter, so this is legacy
  void fitTotalFeedDownSpectra(int a_partIndex, int a_charge, int a_centIndex,string a_imageDir = "");

  double getSigmaPlusOverLambdaRatio(vector<int> a_centIndicies);
  double getAntiSigmaMinusOverAntiLambdaRatio(vector<int> a_centIndicies);
  double getSigmaPlusOverLambdaRatio(int a_centIndex);
  double getAntiSigmaMinusOverAntiLambdaRatio(int a_centIndex);
  //#### Write the caluclated spectra to a file (the whole internal data of this class)
  void   write(string a_outFileName);

  //#### Do Not Use : must be public for minimizer to access. This is the chi squared function for matching the data to the simulation
  double simToDataMatchingChiSqr(const double* a_param);
  double simToDataMatchingChiSqr_BinSplitting(const double* a_param);


private:

  ParticleInfo* m_partInfo;
  bool m_doNTuple;
  double m_feeddown_particle_masses[20];
  string m_imageDir;
  string m_feeddown_particle_names[20];
  string m_feeddown_particle_symbols[20];
  int m_currentFeedDownIndex;
  TGraphAsymmErrors* m_currentDataGraphToFit;
  TGraphAsymmErrors* m_currentSimGraphToFit;
  vector<double>     m_currentDataHistoToFit_BinEdges;
  vector<double>     m_currentSimHistoToFit_BinEdges;
  TH1D*              m_currentDataHistoToFit;
  TH1D*              m_currentSimHistoToFit;
  TH1D*              m_currentSimHistoToFit_inDataBinning;
  double m_matchingParameters[20][16][4]; // FeedDownParticleIndex centrality scale/stretch/scaleErr/stretchErr
  TH3D*  m_matchingParametersHisto;
  bool m_use_DCA_cut;
  double m_DCA_Cut;
  bool m_symmeterizeOutputs;


  //#### URQMD ####
  TH2D* m_urqmd_spectra[20][16]; // FeedDownParticleIndex centrality (centrality index of simulation)
  TH2D* m_urqmd_yields[20][16]; // FeedDownParticleIndex centrality : not invariant, just counts in bin (makes it easy)
  TH2D* m_urqmd_dNdy[20]; // FeedDownParticleIndex
  TH2D* m_urqmd_dNdy_scaled[20]; // FeedDownParticleIndex (centrality index of simulation)
  TH1D* m_urqmd_dNdy_scaled_byCent[20][16]; // FeedDownParticleIndex (centrality index of simulation)

  //#### DATA ####
  TGraphAsymmErrors* m_measuredSpectra[20][16]; // FeedDownParticleIndex in mTm0 space (centrality index of measurment)
  TGraphAsymmErrors* m_mappedUrQMDSpectra[20][16]; // stretch and scale (centrality index of measurment)
  TGraphAsymmErrors* m_unmappedUrQMDSpectra[20][16]; // no stretch and scale (centrality index of measurment)
  TH2D* m_mappedUrQMD_dNdy[20]; // FeedDownParticleIndex // this is the dN/dy after scaling (centrality index of measurment)
  TH1D* m_mappedUrQMD_dNdyByCent[20][16]; // FeedDownParticleIndex // this is the dN/dy after scaling (centrality index of measurment)

  TF1* m_measuredSpectra_Fits[20][16];  //(centrality index of measurment)


  //##########  Toy Monte Carlo Parent and Daughter Histograms
  TH2D* m_response_structure;
  TH1D* m_response_structure_1D;
  TH1D* m_produced_feeddown_dNdy[9][20][16][2]; // given the simulation, these are the feed down dN/dy. For example, proton by y and mTm0 from lambda decays (given the spectra of lambda)
  TH2D* m_produced_feeddown_yields[9][20][16][2]; // given the simulation, these are the feed down yields. For example, proton by y and mTm0 from lambda decays (given the spectra of lambda)
  TH2D* m_produced_feeddown_spectra[9][20][16][2]; // given the simulation, these are the feed down invariant spectra. For example, proton by y and mTm0 from lambda decays (given the spectra of lambda)
  TH1D* m_produced_feeddown_DCA[9][20][16][2]; // given the simulation, these are the feed down invariant spectra. For example, proton by y and mTm0 from lambda decays (given the spectra of lambda)
  TH2D* m_produced_feeddown_dNdy_DCA[9][20][16][2]; // given the simulation, these are the feed down dN/dy by DCA. For example, proton by y and mTm0 from lambda decays (given the spectra of lambda)


  TH2D* m_produced_total_feeddown_yields[9][16][2]; // sum of all the spectra from all parents
  TH2D* m_produced_total_feeddown_spectra[9][16][2]; // sum of all the spectra from all parents
  TF1*  m_produced_total_feeddown_spectraFits[9][16][2][50]; // sum of all the spectra from all parents  // part cent +/- rapIndex
  TMatrixDSym* m_produced_total_feeddown_spectraFitsCovMat[9][16][2][50]; // cov matrix of the fit of sum of all the spectra from all parents  // part cent +/- rapIndex
  TH1D* m_produced_total_feeddown_dNdy[9][16][2]; // sum of all the dN/dy from all parents
  TH2D* m_produced_total_feeddown_dNdy_DCA[9][16][2]; // sum of all the dN/dy from all parents with DCA as another axis

  TH1D* m_resampled_feeddown_dNdy[20][16];
  TH2D* m_resampled_feeddown_yields[20][16]; // given the simulation, these are the resampled parent y and mT-m0. Just a check to make sure the resampling actually follows the parent spectra
  TH2D* m_resampled_feeddown_spectra[20][16]; // given the simulation, these are the resampled parent y and mT-m0. Just a check to make sure the resampling actually follows the parent spectra

  // contains all simulated decays
  TNtuple* m_response_nTuple; // ID_feed y_feed mtm0_feed ID_d1 y_d1 mtm0_d1 ID_d2 y_d2 mTm0_d2 ID_d3 y_d3 mTm0_d3



  double calcScalingFactor(double a_par0, double a_par1, double a_mTm0, double a_mTm0_width, double a_mass);
  void symmeterizeTotalHistograms();
  double dca_from_decay_helix(ROOT::Math::XYZVector* a_primaryVertex, ROOT::Math::XYZVector* a_secondaryVertex, ROOT::Math::PxPyPzEVector* a_daughter_vect, int a_daughter_charge);




  //midrapidity data spectra




};






#endif
