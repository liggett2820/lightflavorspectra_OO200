#ifndef _HISTOGRAM_UTILITIES_
#define _HISTOGRAM_UTILITIES_

#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TPaveText.h"

#include "../makefile_toggles.h"
#include "../headers/PhysMath.h"

using namespace std;


//unsigned long g_histoUtil_uniqueIdentifier;


namespace HistogramUtilities{

  //void setUniqueNum(unsigned long a_num);

  //##################    HISTOGRAM OPERATIONS ###########################

  void     setBinValues(TH1D* a_histo, int a_bin, double a_value, double a_error);
  void     setBinValues(TH2D* a_histo, int a_binX, int a_binY, double a_value, double a_error);

  int      firstBinWithContent(TH1D* a_histo); // will return negative bin if empty;
  int      lastBinWithContent(TH1D* a_histo);
  int      firstYBinWithContent(TH2D* a_histo, int a_binX);
  int      lastYBinWithContent(TH2D* a_histo, int a_binX);

  int      numBinsWithContentAboveThreshold(TH1D* a_histo, double a_threshold = -1);
  double   getLowEdge(TH1D* a_histo);
  double   getHighEdge(TH1D* a_histo);
  double   getBinContentAvg(TH1D* a_histo, int a_bin, int a_numBinsToAvg);
  int      getMaxBinInRange(TH1D* a_histo, double a_lowRange = -999, double a_highRange = -999);
  double   getMaxInRange(TH1D* a_histo, double a_lowRange = -999, double a_highRange = -999);
  int      getMinBinInRange(TH1D* a_histo, double a_lowRange = -999, double a_highRange = -999);
  double   getMinInRange(TH1D* a_histo, double a_lowRange = -999, double a_highRange = -999);
  double   getSmallestBinWidth(TH1D* a_histo, int a_lowBin = -1, int a_highBin = -1); // leave as -1 if you want range 1 - nBins
  void     printHistoInfo(TH1D* a_histo);
  void     printHistoInfo(TH2D* a_histo);

  //Takes an initial bin and steps left and right to find when the content is a percentage of the initial bin
  void getPeakRange(TH1D* a_histo, int a_startBin, int &a_lowBinReturn, int &a_highBinReturn, double a_percentageRequired = 0.05);
  #ifdef _CPP11_
    vector< vector< double > > findPeaks(TH1D* a_histo);
  #endif

  double sumOfContent(TH1D* a_histo, double a_threshold= -1, bool a_isDifferential = false, double a_lowRange = -999, double a_highRange = -999); // is inclusive, does not slice bin if range is in middle of bin.
  //double sumOfContent(TH2D* a_histo, double a_threshold= -1, bool a_isDifferential = false);
  double sumOfContent2D(TH2D* a_histo, double a_threshold= -1, bool a_isDifferential = false);

  double integral(TH1D* a_histo,double a_lowRange = -999, double a_highRange = -999, double a_integral_error = -999); // this slices bins if range is in middle of bin
  double meanOfYVals(TH1D* a_histo,double a_lowRange = -999, double a_highRange = -999); // uses the y error bars to calculate the mean
  double meanOfZValsFromProj(TH2D* a_histo, int a_xBin, double a_lowYRange = -999, double a_highYRange = -999); // uses the z error bars to calculate the mean
  bool   isEmpty(TH1D* a_histo);
  bool   isEmpty2D(TH2D* a_histo);
  bool   isYSliceEmpty(TH2D* a_histo, int a_binX);
  bool   outsideHistogramRange(TH1D* a_histo,double a_value);
  TH1*   makeBinStructure(TH1* a_hist);
  TH1D*  makeBinStructure(TH1D* a_hist);

  vector< double > getBinEdges(TH1* a_histo); // lowbin of 1 to low bin edge of nBins+1
  vector< double > getBinEdges_X(TH2* a_histo);
  vector< double > getBinEdges_Y(TH2* a_histo);
  vector< double > getBinEdges_X(TH3* a_histo);
  vector< double > getBinEdges_Y(TH3* a_histo);
  vector< double > getBinEdges_Z(TH3* a_histo);

  vector< double > getSmallerBinEdges(TH1* a_histo, int a_divNum = 2); // divides every bin in pieces
  vector< double > getBinCenters(TH1* a_histo); // {underflow, bin1, bin2, ... lastBin, overflow} size nBins+2
  vector< double > getBinContents(TH1D* a_histo); // {underflow, bin1, bin2, ... lastBin, overflow} size nBins+2
  vector< double > getBinErrors(TH1D* a_histo); // {underflow, bin1, bin2, ... lastBin, overflow} size nBins+2
  vector< double > makeEvenBinEdges(int a_numBins, double a_lowBound, double a_highBound);
  vector<double>   makeLogBinEdges(int a_numBins, double a_lowBound, double a_highBound);
  vector<double>   makeAdaptiveBinEdgesEqualCounts(TH1D* a_histo,double a_minCounts, bool a_roundToMillionsPlace = false); // returns bin edges so that each bin (after a rebin) would have about a_minCounts

  double normalize(TH1D* a_histo, bool a_isDifferential = false, double a_lowRange = -999, double a_highRange = -999); // returns the normalization factor   new = N * old    | N is returned
  void normalizeEachVerticalSlice(TH2D* a_histo, bool a_isDifferential = false);
  void makeDifferential(TH1D* a_histo, double a_extraScale = 1.0);
  void makeDifferential(TH2D* a_histo, double a_extraScale = 1.0);
  void makeNotDifferential(TH1D* a_histo, double a_extraScale = 1.0);
  void makeNotDifferential(TH2D* a_histo, double a_extraScale = 1.0);
  //mostly for the spine fitting fixed points
  vector< double > getEdgesByVerticalChange(TH1D* a_histo, int a_numBins, bool a_weightByVariance = true, double a_lowRange = -999,double a_highRange = -999); // rebins so more bins where content is changing, less bins where there is less change
  vector< double > getEdgesByVerticalChange(TGraphAsymmErrors* a_graph, int a_numBins, bool a_weightByVariance = true, double a_lowRange = -999,double a_highRange = -999); // rebins so more bins where content is changing, less bins where there is less change


  TH2D* addTH2D(TH2D* a_histo1, TH2D* a_histo2);
  void zeroBins(TH1D* a_histo);
  void zeroBins2D(TH2D* a_histo);
  void zeroBins3D(TH3D* a_histo);
  void zeroErrorBars2D(TH2D* a_histo);
  void zeroBins(TH1I* a_histo);
  void zeroBins2D(TH2I* a_histo);
  void zeroBins3D(TH3I* a_histo);
  void zeroInfBins(TH1D* a_histo);
  void zeroBinsBelowThreshold(TH1D* a_histo,double a_threshold);
  void discardNanHistoData(TH1D* a_histo, bool a_coupleDataAndError = true);
  void discardNanHistoData(TH2D* a_histo, bool a_coupleDataAndError = true);
  void removeHighCountCorruptData(TH1D* a_histo, double a_thresholdFactor = 100.0); // some corruption happens in fitting where one or two bins are 10^4 higher than the rest. This removes that point
  void symmeterizeHistoX(TH1D* a_histo);
  void symmeterizeHistoX(TH2D* a_histo);
  void setYAxisRanges_collectionOfHistos(vector<TH1D*> a_histo_vect, bool a_logScale = false); //if you're drawing a bunch of histos, but don't know which histo has the biggest range
  void setAxisRanges(TH1D* a_histo, double a_lowRangeX, double a_highRangeX, bool a_logScale = false, bool a_includeErrorBars = true);

  TH1D* shiftHisto(TH1D* a_histogram, double a_shift);
  TH2D* makeScaleAndStretchHisto(TH2D* a_histo, double a_scale, double a_stretch); // scale in Z and stretch in Y

  void chi_squared(TH1D* a_histo_1, TH1D* a_histo_2,double &a_chi_sqr, double &a_ndf);

  void  interpolate(TH1D* a_histo, double a_x, double &a_y, double &a_yErr);
  void  interpolate2D(TH2D* a_histo, int a_binX, double a_y, double &a_z, double &a_zError);
  TH1D* convertFromLogBins(TH1D* a_histo); // if log spaced bins -> even bins in log space ie. bin edge -> ln(bin edge)
  TH1D* convertInt_to_DoubleHisto(TH1I* a_histo);
  TH2D* convertInt_to_DoubleHisto2D(TH2I* a_histo);
  TH3D* convertInt_to_DoubleHisto3D(TH3I* a_histo);

  TGraphErrors*      makeGraph(TH1D* a_histo, bool a_removeZeros = true,double a_lowRange = -999, double a_highRange = -999, bool a_removeZeroErrors = false);
  TGraphAsymmErrors* makeGraphAsymm(TH1D* a_histo, bool a_removeZeros = true,double a_lowRange = -999, double a_highRange = -999, bool a_removeZeroErrors = false);
  TGraphErrors*      makeGraph(TH2D* a_histo, int a_binX);
  TGraph2DErrors*    makeGraph2D(TH2D* a_histo, TGraph2DErrors* a_graph = NULL); // makes new one or you can give it a pointer to fill
  TGraph2D*          makeErrorGraph2D(TH2D* a_histo); // puts the content errors into a 2D Graph to find error on interpolations

  TH1D* makeFitRatioHistogram(TH1D* a_histo, TF1* a_funct, bool a_doIntegral = false); // makes a histo/funct ratio plot
  TH1D* makeFitPullHistogram(TH1D* a_histo, TF1* a_funct, bool a_doIntegral = false); // makes a histo/funct ratio plot


  //returns counts, mean, and std dev in range
  //if use errors is true, the errors are used for the weights of the mean
  //NOTE: This must be counts or counts/bin width (NO SCALING bc of N/N-1 varaince bias)
  void getGausData(double* a_returnedData, TH1D* a_histo, bool a_isDifferential = false, bool a_isntCounts = false, double a_lowBound = -999, double a_highBound = -999);


  #ifdef _MAC_OSX_
  vector<double> calcGausData(TH1D* a_histo, double a_guessMean, int a_numBinsToScan = 5, int a_binScanWidth = 3);
  #endif


  //this adds a value, but spreads it over a range and cuts it up according to the bins
  void  addToHistogram_binSplitting(TH1D* a_histo, TH1D* a_histoToAdd);
  void  addToHistogram(TH1D* a_histo, double a_value, double a_value_err, double a_lowRange, double a_highRange);
  TH1D* truncateHistogram(TH1D* a_histo, double a_lowCut, double a_highCut);
  TH2D* truncateHistogram2D(TH2D* a_histo, double a_lowCut_X, double a_highCut_X, double a_lowCut_Y, double a_highCut_Y);

  //returns xlow,xhigh, ylow, yhigh
  //it is used for making a smaller projection, cutting white space
  void binBounds2D(TH3D* a_histo, int a_zBin, int* a_returnVars, double a_minThreshold = 1);
  void binBounds2D(TH3I* a_histo, int a_zBin, int* a_returnVars, double a_minThreshold = 1);

  void makeSqrtNErrors(TH1D* a_histo);
  void makeSqrtNErrors2D(TH2D* a_histo);
  void makeSqrtNErrors3D(TH3D* a_histo);

  TH1D* projectionY(TH2D* a_histo, double a_lowVal = -999, double a_highVal = -999);
  TH1I* projectionY(TH2I* a_histo, double a_lowVal = -999, double a_highVal = -999);
  //   binLow <= bin <= binHigh
  // only even binned histograms
  // NOTE: reduction only works for single bin at the moment
  TH2D* projection2D(TH3D* a_histo, int a_binZLow, int a_binZHigh, bool a_reduceBounds = false, double a_boundThreshold = 1);
  TH2I* projection2D(TH3I* a_histo, int a_binZLow, int a_binZHigh, bool a_reduceBounds = false, double a_boundThreshold = 1);

  TH2D* projection2D_XY(TH3D* a_histo, int a_binZ);
  TH2D* projection2D_XZ(TH3D* a_histo, int a_binY);
  TH2D* projection2D_YZ(TH3D* a_histo, int a_binX);


  //This forces new creation of histogram, and makes a unique name if the given histo has a unique name
  TH1D* projectionZ(TH3D* a_histo, int a_binX, int a_binY, string a_uniqueIdentifier = "");
  vector< vector<TH1D*> > projectionsZ(TH3D* a_histo, bool a_makeDifferential = false, bool a_normalize = false);

  vector< TH2I* > projeciton2DSlices_X(TH3I* a_histo); // steps through x producing all 2D slices
  vector< TH2D* > projeciton2DSlices_X(TH3D* a_histo); // steps through x producing all 2D slices
  vector< TH2I* > projeciton2DSlices_Z(TH3I* a_histo, vector<double> a_zValueCuts); // steps through z producing all 2D slices, forming 2D projections within the Z cuts
  vector< TH2D* > projeciton2DSlices_Z(TH3D* a_histo, vector<double> a_zValueCuts); // steps through z producing all 2D slices, forming 2D projections within the Z cuts


  TH2D* crop2DHisto( TH2D* a_histoToCrop,double a_lowX, double a_highX, double a_lowY, double a_highY);

  ///calibrates the histogram passed to the function using a linear calibration of the form a*X +b
  void applyCalibration(TH1*   a_histogram,///pointer to the histogram to apply the calibration to
                        double a_a,///<linear factor in the calibration
                        double a_b///<constant to change the x intercept
                        );



  //must already be in order
  TH1D* stitchHistos(int a_numHistos, TH1D** a_histos);

  TH1D* logifyHisto(TH1D* a_histo, double a_minContent = 0.01);
  TH1D* logifyHistoHorizontalAxis(TH1D* a_histo);

  TH1D* make1DHistLogBins(string a_name, string a_title, int a_numBins, double a_start, double a_stop);
  TH2D* make2DHistLogBins(string a_name, string a_title,
                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY);
  TH3D* make3DHistLogBins(string a_name, string a_title,
                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY,
                          bool a_makeLogZ, int a_numBinsZ, double a_startZ, double a_stopZ);

  TH1I* make1DHistLogBinsInt(string a_name, string a_title, int a_numBins, double a_start, double a_stop);
  TH2I* make2DHistLogBinsInt(string a_name, string a_title,
                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY);
  TH3I* make3DHistLogBinsInt(string a_name, string a_title,
                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY,
                          bool a_makeLogZ, int a_numBinsZ, double a_startZ, double a_stopZ);


  // makes a bunch of TGraphErrors for each y slice of your histogram (index 0 of return is bin 1)
  vector<TH1D* >         HorizontalProjectionHistos(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol);
  vector<TH1D* >         VerticalProjectionHistos(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol);
  vector<TGraphErrors* > VerticalProjectionGraphs(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol, bool a_removeZeros = true);
  vector<TGraph* >       VerticalProjectionGraphsNoErrors(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol, bool a_removeZeros = true);



  double scaledEntropy(TH1D* a_histo);
  TH1D* convolution(TH1D* a_histo_1,TH1D* a_histo_2); // must have equal bin widths. note overall scaling not correct
  TH1D* selfConvolution(TH1D* a_histo); // must have equal bin widths. note overall scaling not correct
  TH1D* rebinMinCount(TH1D* a_histo, double a_minCount); // merges bins if the counts are lower than minCounts. //must be a histogram containing counts (so sqrt(N) error)
  TH1D* rebinEqualCounts(TH1D* a_histo, int a_numBins); // calculates the sum of counts and from the number of bins, tries to make each bin count equal. //must be a histogram containing counts (so sqrt(N) error)
  TH1D* rebinEqualByVerticalChange(TH1D* a_histo, int a_numBins); // rebins so more bins where content is changing, less bins where there is less change
  TH1D* rebinNoZeroBins(TH1D* a_histo, bool a_isDifferential); // from first bin content to last bin content, it widens bins so there are no zero bins. if differential, will return differential
  TH1D* rebinEqualByCDF(TH1D* a_histo, int a_numBins, bool a_isDifferential = false, bool a_useLnWeight = false, bool a_makeNoZeroBins = true); // rebins so more bins where content is changing, less bins where there is less change
  //The stock rebin function for ROOT does not allow for variable bin widts, these functions do what you would expect a rebin to do
  TH3D* reBinVar3D(TH3D* a_histo, int a_rebinX, int a_rebinY, int a_rebinZ);
  TH2D* reBinVar2D(TH2D* a_histo, int a_rebinX, int a_rebinY);
  TH3I* reBinVar3D(TH3I* a_histo, int a_rebinX, int a_rebinY, int a_rebinZ);
  TH2I* reBinVar2D(TH2I* a_histo, int a_rebinX, int a_rebinY);

  // dont do ranges that overlap...
  TH2I* reBinBar2D_specificY(TH2I* a_histo, vector< vector<double> > a_rebinInfo); // rebinInteger lowYRange highYRange
  TH2D* reBinBar2D_specificY(TH2D* a_histo, vector< vector<double> > a_rebinInfo);
  TH3I* reBinBar3D_specificY(TH3I* a_histo, vector< vector<double> > a_rebinInfo); // rebinInteger lowYRange highYRange
  TH3D* reBinBar3D_specificY(TH3D* a_histo, vector< vector<double> > a_rebinInfo);

  //This protects against Runge effect in polynomial fitting.
  TH1D* rebinHistoToChebyshevSpacing(TH1D* a_histo,double a_lowRange, double a_highRange, int a_reBinFactor = 3, int a_numPoints = -999); // numPoints overrides the rebin factor
  TH2D* rebin2DHistoToVertChebyshevSpacing(TH2D* a_histo,double a_lowRange, double a_highRange, int a_reBinFactor = 3, int a_numPoints = -999); // numPoints overrides the rebin factor
  //This rebin will split and add bins to adapt to another bin structure. Must be a simple count histogram so uncertianty makes sense
  TH1D* rebinHistoToAlternateBinStructure(TH1D* a_histo, vector<double> a_binEdges); // this creates a histo with the bin edges supplied
  void  rebinHistoToAlternateBinStructure(TH1D* a_histo, TH1D* a_histoToPopulate); //this doesn't zero bins of populate histo
  TH2D* rebin2DHistoToAlternateVertBinStructure(TH2D* a_histo, vector<double> a_binEdges); // this creates a histo with the bin edges supplied
  void  rebin2DHistoToAlternateVertBinStructure(TH2D* a_histo, TH2D* a_histoToPopulate); //this doesn't zero bins of populate histo

  TH1D* convertFromYieldToSpectra(TH1D* a_mTm0YieldHisto, double a_centEvents, double a_rapWidth, double a_mass, double a_efficiency = 1.0, double a_efficiencyErr = 0.0);
  TH1D* convertFromSpectraToYield(TH1D* a_mTm0SpectraHisto, double a_centEvents, double a_rapWidth, double a_mass);
  TH2D* convertFromYieldToSpectra2D(TH2D* a_mTm0YieldHisto, double a_centEvents,  double a_mass, double a_efficiency = 1.0, double a_efficiencyErr = 0.0);
  TH1D* convertZTOFToMassSqr(TH1D* a_invBetaHisto, double a_p, double a_m_calib, bool a_ZTOF_isDifferential = false, bool a_massSqr_isDifferential = false);
  TH1D* convertMassSqrToInvBeta(TH1D* a_massSqrHisto,double a_p, bool a_massSqr_isDifferential = false, bool a_ZTOF_isDifferential = false, double a_m_calib = -999); // if a_m_calib mass is definined, it will do 1/beta -1/beta_predict = ZTOF
  TH1D* convertSpectraHisto_mTm0_to_pT(TH1D* a_histo, double a_mass);
  TH1D* convertSpectraHisto_pT_to_mTm0(TH1D* a_histo, double a_mass);







  //##################    HISTOGRAM OPERATIONS TESTERS ###########################
  void testRebinEqualByCDF();
  void testRebinNoZeros();
  void test_functionEvalWithError();
  void test_rebinHistoToAlternateBinStructure();
  void test_rebinHistoToChebyshevSpacing();
  void test_rebin2DHistoToVertChebyshevSpacing();
  void test_histoAddWithBinSplitting();
  void test_m2_invBeta_conversion();
  void test_self_convolution();
  void test_convolution();








  //##################    TF1 OPERATIONS ###########################
  void functionEvalWithError(TF1* a_funct,double a_x, double &a_eval, double &a_evalErr, TMatrixDSym* a_covMatrix = NULL);
  void printTF1Info(TF1* a_funct);
  bool parametersAreFinite(TF1* a_funct);
  TPaveText* makeTF1TextBox(TF1* a_funct, double a_x1 = 0.65, double a_y1 = 0.6, double a_x2 = 0.95, double a_y2 = 0.97);

  //##################    TGRAPH OPERATIONS ###########################
  void  printGraphInfo(TGraphAsymmErrors* a_graph);
  void  printGraphInfo(TGraphErrors* a_graph);
  void  printReconstrucitonCode_TGraphErrors(TGraphErrors* a_graph);
  TGraphErrors* makeGraphFromDataTheif(string a_name, int a_numPoints, double* a_xVal, double* a_yVal, double* a_yValPlusErr, double* a_xValPlusErr = NULL);
  TH1D* makeHistogram(TGraphErrors* a_graph); // only if x error bars correspond to bin widths
  TH1D* makeHistogram_Asymm(TGraphAsymmErrors* a_graph);
  int   getPointClosestToVal(TGraph* a_graph, double a_x);
  int   getPointClosestToVal(TGraphErrors* a_graph, double a_x);
  int   getPointClosestToVal(TGraphAsymmErrors* a_graph, double a_x);
  void  evalWithError(TGraphErrors* a_graph, double a_evalX, double &a_yVal, double &a_yErr);
  void  evalWithError(TGraphAsymmErrors* a_graph, double a_evalX, double &a_yVal, double &a_yErr);
  void  evalWithAsymmError(TGraphAsymmErrors* a_graph, double a_evalX, double &a_yVal, double &a_yErrLow, double &a_yErrHigh);
  void  truncateTGraph(TGraph* a_graph, double a_lowCut, double a_highCut);
  void  scaleTGraph(TGraphErrors* a_graph, double a_scale);
  void  scaleTGraph(TGraphAsymmErrors* a_graph, double a_scale);
  void  shiftTGraph(TGraphErrors* a_graph, double a_shift); // x = x + a_shift
  void  zeroTGraph(TGraphErrors* a_graph); // keeps x locations and xErr but zeros y and yErr
  void  removeBelowThreshold(TGraphErrors* a_graph,double a_threshold);
  void  removeBelowThreshold(TGraphAsymmErrors* a_graph,double a_threshold);
  void  discardNanData(TGraphErrors* a_graph, bool a_coupleDataAndError = false);
  void  discardNanData(TGraphAsymmErrors* a_graph, bool a_coupleDataAndError = false);
  TGraphErrors*      addTGraphErrors(TGraphErrors* a_graph1, TGraphErrors* a_graph2 = NULL, TGraphErrors* a_graph3 = NULL);
  TGraphErrors*      convertAsymToNoAsym(TGraphAsymmErrors* a_graph);
  TGraphAsymmErrors* convertNoAsymToAsym(TGraphErrors* a_graph);
  double findSmallestYVal(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);
  double findLargestYVal(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);
  double findSmallestXVal(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);
  double findLargestXVal(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);
  double findSmallestYVal(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);
  double findLargestYVal(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);
  double findSmallestXVal(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);
  double findLargestXVal(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_includeErrorBar = false);

  int findSmallestYValPoint(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);
  int findLargestYValPoint(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);
  int findSmallestXValPoint(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);
  int findLargestXValPoint(TGraphErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);

  int findSmallestYValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);
  int findLargestYValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);
  int findSmallestXValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);
  int findLargestXValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange = -999, double a_highXRange = -999);

  //Uses Trapazoidal rule to find the integral with error. This gets more wrong if it is outside the range of the TGraphErrors
  void integrateGraph(TGraphErrors* a_graph, double a_lowRange, double a_highRange, double &a_integral, double &a_integralError, bool a_divideByRange = false);

  void writeTgraphErrorsToFile(TGraphErrors* a_graph, TString a_outFile);
  void writeTgraphErrorsWithSysToFile(TGraphErrors* a_graphNom, TGraphErrors* a_graphSys, TString a_outFile);

  TGraphErrors*      makeDerivativeGraph(TGraphAsymmErrors* a_graph);
  TGraphErrors*      makeDerivativeGraph(TGraphErrors* a_graph);
  TGraphAsymmErrors* makeFitPullGraph(TGraphAsymmErrors* a_graph, TF1* a_funct, bool a_doIntegral = false);
  TGraphAsymmErrors* makeFitPullGraph(TGraphErrors* a_graph, TF1* a_funct, bool a_doIntegral = false);
  TGraphAsymmErrors* makeFitPullGraph(TGraphAsymmErrors* a_graph, TGraphAsymmErrors* a_graph2);
  TGraphErrors*      makeErrorBandGraph(TF1* a_funct, TMatrixDSym* a_covMatrix, int a_num_sample_points = 300); // draw with "lf,same" option

  // 1/mT * dN/dmT -> 1/pT * dN/dpT    1/mT dN/dmT is input in mTm0 space to 1/pT dN/dpT in pT space
  TGraphAsymmErrors* convertSpectra_mTm0_to_pT(TGraphErrors* a_graph, double a_mass);
  TGraphAsymmErrors* convertSpectra_mTm0_to_pT(TGraphAsymmErrors* a_graph, double a_mass);

  // 1/pT * dN/dpT -> 1/mT * dN/dmT    1/pT dN/dpT is input in pT space to 1/mT dN/dmT in mTm0 space
  TGraphAsymmErrors* convertSpectra_pT_to_mTm0(TGraphErrors* a_graph, double a_mass);
  TGraphAsymmErrors* convertSpectra_pT_to_mTm0(TGraphAsymmErrors* a_graph, double a_mass);

  // dmT/dpT * dN/dmT = dN/dpT    dN/dmT is input in mT space    dN/dpT is output in pT space
  TGraphAsymmErrors* convertPlot_mT_to_pT(TGraphErrors* a_graph, double a_mass);

  //Used for TH1D->TGraphErrors->Invaraint (dividing by rap width & mTm0 width & numEvents...)
  TGraphErrors*      convertFromYieldToSpectra(TGraphErrors* a_mTm0YieldGraph,double a_centEvents, double a_rapWidth, double a_mass, double a_efficiency = 1.0, double a_efficiencyErr = 0.0);
  //Multiply by 2*pi*mT (for dNdy calculations / plotting)
  TGraphErrors*      convertFrom_Invariant_To_NotInvariant_Spectra(TGraphErrors* a_mTm0SpectraGraph, double a_mass);
  TGraphAsymmErrors* convertFrom_Invariant_To_NotInvariant_Spectra(TGraphAsymmErrors* a_mTm0SpectraGraph, double a_mass);

  TGraphAsymmErrors* convertFrom_pT_To_betaGamma(TGraphAsymmErrors* a_graph, double a_mass);


  //Note: bin errors must be full widths, not energy-loss propagated error. (errors in x must touch between points)
  void calculateFiducialDnDy(TGraphErrors*      a_NotInvariantSpectra, double &a_dNdy, double &a_dNdyErr); //use convertFrom_Invariant_To_NotInvariant_Spectra before this if needed
  void calculateFiducialDnDy(TGraphAsymmErrors* a_NotInvariantSpectra, double &a_dNdy, double &a_dNdyErr); //use convertFrom_Invariant_To_NotInvariant_Spectra before this if needed

  // if there are gaps in the middle of the spectra, it will return {{xLow,xHigh},{xLow,xHigh},...}
  vector< vector<double> > findSpectraHoles(TGraphErrors* a_spectra, double a_tolerance = 0.00001);
  vector< vector<double> > findSpectraHoles(TGraphAsymmErrors* a_spectra, double a_tolerance = 0.00001);







  //##############   #NOTE: TEMPLATE FUNCTIONS MUST BE DEFINED AND DECLARED AT THE SAME TIME FOR COMPLIATON!
  template <class T>
  void verticalDrawRange(T* a_graph, double &a_lowRange,double &a_highRange, double a_lowXRange = -999, double a_highXRange = -999, bool a_isLogY = false, bool a_includeErrorBars = true){
    double lowVal  = findSmallestYVal(a_graph,a_lowXRange,a_highXRange,a_includeErrorBars);
    double highVal = findLargestYVal(a_graph,a_lowXRange,a_highXRange,a_includeErrorBars);
    double width   = highVal - lowVal;
    if(a_isLogY){
      a_lowRange  = 0.5*lowVal;
      a_highRange = 2.0*highVal;
      if(a_lowRange < 0.0 || a_highRange < 0.0){
        cout << "ERROR: HistogramUtilities::verticalDrawRange... you are using a negative for a log scale" << endl;
      }
    }else{
      a_lowRange  = lowVal  - 0.1*width;
      a_highRange = highVal + 0.1*width;
    }
  };


  template <class T>
  void horizontalDrawRange(T* a_graph, double &a_lowRange,double &a_highRange, double a_lowXRange = -999, double a_highXRange = -999, bool a_isLogX = false, bool a_includeErrorBars = true){
    int lowPoint  = findSmallestXValPoint(a_graph,a_lowXRange,a_highXRange);
    int highPoint = findLargestXValPoint(a_graph,a_lowXRange,a_highXRange);
    if(lowPoint < 0 || highPoint < 0){
      cout << "WARNING: HistogramUtilities::horizontalDrawRange can't find any points..." << endl;
      return;
    }
    double lowVal  = a_graph->GetX()[lowPoint];
    double highVal = a_graph->GetX()[highPoint];
    if(a_includeErrorBars){
      lowVal  = findSmallestXVal(a_graph,a_lowXRange,a_highXRange,a_includeErrorBars);
      highVal = findLargestXVal(a_graph,a_lowXRange,a_highXRange,a_includeErrorBars);
    }


    double width   = highVal - lowVal;
    if(a_isLogX){
      a_lowRange  = 0.5*lowVal;
      a_highRange = 2.0*highVal;
      if(a_lowRange < 0.0 || a_highRange < 0.0){
        cout << "ERROR: HistogramUtilities::horizontalDrawRange... you are using a negative for a log scale" << endl;
      }
    }else{
      a_lowRange  = lowVal  - 0.1*width;
      a_highRange = highVal + 0.1*width;
    }
  };

  template <class T>
  TH1D* generateGraphFrame(T* a_graph, double a_lowXRange = -999, double a_highXRange = -999, bool a_isLogX = false, bool a_isLogY = false, string a_XAxisLabel = "", string a_YAxisLabel = "",bool a_includeErrorBars = true){
    double lowRangeX = 0;
    double highRangeX = 0;
    double lowRangeY = 0;
    double highRangeY = 0;
    horizontalDrawRange(a_graph, lowRangeX,highRangeX, a_lowXRange, a_highXRange, a_isLogX, a_includeErrorBars);
    verticalDrawRange(a_graph, lowRangeY,highRangeY, a_lowXRange, a_highXRange, a_isLogY, a_includeErrorBars);
    string name = a_graph->GetName();
    name += "_frame";
    string title = Form("%s;%s;%s",a_graph->GetTitle(),a_XAxisLabel.c_str(),a_YAxisLabel.c_str());
    TH1D* frame = new TH1D(name.c_str(),title.c_str(),50,lowRangeX,highRangeX);
    for(int bin = 1; bin < 51; bin++){
      frame->SetBinContent(bin,lowRangeY-1.0);
    }
    frame->GetXaxis()->SetRangeUser(lowRangeX,highRangeX);
    frame->GetYaxis()->SetRangeUser(lowRangeY,highRangeY);
    return frame;
  };


  template <class T>
  void ConditionalWrite(T* a_object,string a_name = ""){
    if(a_object){
      if(a_name != "") a_object->Write(a_name.c_str());
      else             a_object->Write();
    }
  };

  template <class T>
  void ConditionalDraw(T* a_object,string a_option = ""){
    if(a_object)  a_object->Draw(a_option.c_str());
  };

  template <class T>
  void ConditionalDelete(T* a_object){
    if(a_object){
      delete a_object;
      //a_object = nullptr;
    }
  };

  TObject* FileGet(TFile* a_file, string a_name); // will print an error if object not found

  template <class T>
  int getNumNonZeroPts(T* a_object){
    int numGoodPts = 0;
    for(unsigned int iii = 0; iii < a_object->GetN(); iii++) if(a_object->GetY()[iii] != 0.0) numGoodPts++;
    return numGoodPts;
  };

  template <class T>
  int getNPointsInRange(T* a_object, double a_lowRange, double a_highRange){
    int numPoints = 0;
    for(unsigned int pointIndex = 0; pointIndex < a_object->GetN(); pointIndex++) if(a_object->GetX()[pointIndex] > a_lowRange && a_object->GetX()[pointIndex] < a_highRange) numPoints++;
    return numPoints;
  };

  template <class T>
  int getFirstNonZeroPoint(T* a_object){
    for(unsigned int pointIndex = 0; pointIndex < a_object->GetN(); pointIndex++) if(a_object->GetY()[pointIndex] != 0.0) return pointIndex;
    return -1;
  };

  template <class T>
  int getLastNonZeroPoint(T* a_object){
    for(int pointIndex = a_object->GetN() - 1; pointIndex >= 0; pointIndex--) if(a_object->GetY()[pointIndex] != 0.0) return pointIndex;
    return -1;
  };

  template <class T>
  int getPointClosestToYVal(T* a_object, double a_yVal, double a_lowRange = -999, double a_highRange = -999){
    int pointIndexClosest = -1;
    double absDifForClosest = 1e200;
    bool useRange = false;
    if(a_lowRange!=a_highRange) useRange = true;
    for(unsigned int pointIndex = 0; pointIndex < a_object->GetN(); pointIndex++){
      if(!useRange ||(useRange && a_object->GetX()[pointIndex] > a_lowRange && a_object->GetX()[pointIndex] < a_highRange)){
        double diff = fabs(a_object->GetY()[pointIndex] - a_yVal);
        if(diff < absDifForClosest){
          absDifForClosest = diff;
          pointIndexClosest = pointIndex;
        }
      }
    }
    return pointIndexClosest;
  };

  template <class T>
  int getPointClosestToXVal(T* a_object, double a_xVal, double a_lowRange = -999, double a_highRange = -999){
    int pointIndexClosest = -1;
    double absDifForClosest = 1e200;
    bool useRange = false;
    if(a_lowRange!=a_highRange) useRange = true;
    for(unsigned int pointIndex = 0; pointIndex < a_object->GetN(); pointIndex++){
      if(!useRange ||(useRange && a_object->GetX()[pointIndex] > a_lowRange && a_object->GetX()[pointIndex] < a_highRange)){
        double diff = fabs(a_object->GetX()[pointIndex] - a_xVal);
        if(diff < absDifForClosest){
          absDifForClosest = diff;
          pointIndexClosest = pointIndex;
        }
      }
    }
    return pointIndexClosest;
  };


  // --------------     set all array values to nullptr;


  template <class T>
  void nullptr_array_8D(T a_object, unsigned int a_dim_1, unsigned int a_dim_2, unsigned int a_dim_3, unsigned int a_dim_4, unsigned int a_dim_5, unsigned int a_dim_6, unsigned int a_dim_7, unsigned int a_dim_8){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
    for(unsigned int index_2 = 0; index_2 < a_dim_2; index_2++){
    for(unsigned int index_3 = 0; index_3 < a_dim_3; index_3++){
    for(unsigned int index_4 = 0; index_4 < a_dim_4; index_4++){
    for(unsigned int index_5 = 0; index_5 < a_dim_5; index_5++){
    for(unsigned int index_6 = 0; index_6 < a_dim_6; index_6++){
    for(unsigned int index_7 = 0; index_7 < a_dim_7; index_7++){
    for(unsigned int index_8 = 0; index_8 < a_dim_8; index_8++){
      a_object[index_1][index_2][index_3][index_4][index_5][index_6][index_7][index_8] = nullptr;
    }}}}}}}}
  };


  template <class T>
  void nullptr_array_7D(T a_object, unsigned int a_dim_1, unsigned int a_dim_2, unsigned int a_dim_3, unsigned int a_dim_4, unsigned int a_dim_5, unsigned int a_dim_6, unsigned int a_dim_7){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
    for(unsigned int index_2 = 0; index_2 < a_dim_2; index_2++){
    for(unsigned int index_3 = 0; index_3 < a_dim_3; index_3++){
    for(unsigned int index_4 = 0; index_4 < a_dim_4; index_4++){
    for(unsigned int index_5 = 0; index_5 < a_dim_5; index_5++){
    for(unsigned int index_6 = 0; index_6 < a_dim_6; index_6++){
    for(unsigned int index_7 = 0; index_7 < a_dim_7; index_7++){
      a_object[index_1][index_2][index_3][index_4][index_5][index_6][index_7] = nullptr;
    }}}}}}}
  };

  template <class T>
  void nullptr_array_6D(T a_object, unsigned int a_dim_1, unsigned int a_dim_2, unsigned int a_dim_3, unsigned int a_dim_4, unsigned int a_dim_5, unsigned int a_dim_6){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
    for(unsigned int index_2 = 0; index_2 < a_dim_2; index_2++){
    for(unsigned int index_3 = 0; index_3 < a_dim_3; index_3++){
    for(unsigned int index_4 = 0; index_4 < a_dim_4; index_4++){
    for(unsigned int index_5 = 0; index_5 < a_dim_5; index_5++){
    for(unsigned int index_6 = 0; index_6 < a_dim_6; index_6++){
      a_object[index_1][index_2][index_3][index_4][index_5][index_6] = nullptr;
    }}}}}}
  };

  template <class T>
  void nullptr_array_5D(T a_object, unsigned int a_dim_1, unsigned int a_dim_2, unsigned int a_dim_3, unsigned int a_dim_4, unsigned int a_dim_5){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
    for(unsigned int index_2 = 0; index_2 < a_dim_2; index_2++){
    for(unsigned int index_3 = 0; index_3 < a_dim_3; index_3++){
    for(unsigned int index_4 = 0; index_4 < a_dim_4; index_4++){
    for(unsigned int index_5 = 0; index_5 < a_dim_5; index_5++){
      a_object[index_1][index_2][index_3][index_4][index_5] = nullptr;
    }}}}}
  };

  template <class T>
  void nullptr_array_4D(T a_object, unsigned int a_dim_1, unsigned int a_dim_2, unsigned int a_dim_3, unsigned int a_dim_4){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
    for(unsigned int index_2 = 0; index_2 < a_dim_2; index_2++){
    for(unsigned int index_3 = 0; index_3 < a_dim_3; index_3++){
    for(unsigned int index_4 = 0; index_4 < a_dim_4; index_4++){
      a_object[index_1][index_2][index_3][index_4] = nullptr;
    }}}}
  };

  template <class T>
  void nullptr_array_3D(T a_object, unsigned int a_dim_1, unsigned int a_dim_2, unsigned int a_dim_3){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
    for(unsigned int index_2 = 0; index_2 < a_dim_2; index_2++){
    for(unsigned int index_3 = 0; index_3 < a_dim_3; index_3++){
      a_object[index_1][index_2][index_3] = nullptr;
    }}}
  };

  template <class T>
  void nullptr_array_2D(T a_object, unsigned int a_dim_1, unsigned int a_dim_2){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
    for(unsigned int index_2 = 0; index_2 < a_dim_2; index_2++){
      a_object[index_1][index_2] = nullptr;
    }}
  };

  template <class T>
  void nullptr_array_1D(T a_object, unsigned int a_dim_1){
    for(unsigned int index_1 = 0; index_1 < a_dim_1; index_1++){
      a_object[index_1] = nullptr;
    }
  };

  //###################################################################################################################


};



#endif
