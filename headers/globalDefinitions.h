//Definitions used throughout the analysis
#ifndef __GLOBAL_DEFINITIONS__
#define __GLOBAL_DEFINITIONS__

#include <TMath.h>

//THE ETA BINS MUST BE IDENTICAL FOR TPC / bTOF / eTOF



//#####################################################
//THIS IS LEGACY!!!!!!! USE THE BINNING IN CUTCLASS!!!!
//#####################################################

#define etaMinTPC_CC -2.05
#define etaMaxTPC_CC 2.05
#define nEtaTPCBins_CC 41

#define nEtaBTOFBins_CC 25
#define etaMinBTOF_CC -1.25
#define etaMaxBTOF_CC 1.25

#define nEtaETOFBins_CC 6
#define etaMinETOF_CC -1.65
#define etaMaxETOF_CC -1.05

#define pTLogMin_CC 0.05
#define pTLogMax_CC 2.5
#define npTLogBins_CC 80


#define etaMinTPC_FXT -2.45
#define etaMaxTPC_FXT 0.05
#define nEtaTPCBins_FXT 25

#define nEtaBTOFBins_FXT 15
#define etaMinBTOF_FXT -1.45
#define etaMaxBTOF_FXT 0.05

#define nEtaETOFBins_FXT 9
#define etaMinETOF_FXT -2.35
#define etaMaxETOF_FXT -1.45

#define pTLogMin_FXT 0.05
#define pTLogMax_FXT 2.5
#define npTLogBins_FXT 80

#define rapidityMin_CC -1.55
#define rapidityMax_CC 1.55
#define rapidityMin_FXT -2.39
#define rapidityMax_FXT 0.21
#define rapidityBinWidth 0.1

#define mTm0Min 0.0
#define mTm0Max 2.0
#define mTm0BinWidth 0.025

//WARNING: If you change the min, max, or width values
//make sure you change these too!
#define nRapidityBins_CC 31
#define nRapidityBins_FXT 26
#define nmTm0Bins 80

//#####################################################
//THIS IS LEGACY!!!!!!! USE THE BINNING IN CUTCLASS!!!!
//#####################################################


#endif

