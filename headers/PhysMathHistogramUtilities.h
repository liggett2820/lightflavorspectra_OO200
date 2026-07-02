#ifndef _PHYSMATHHISTOGRAMUTILITIES_H_
#define _PHYSMATHHISTOGRAMUTILITIES_H_

#include "../makefile_toggles.h"
#include "../headers/PhysMath.h"
#include "../headers/HistogramUtilities.h"
//using namespace PhysMath;
//using namespace HistogramUtilities;

namespace PhysMathHistogramUtilities{

  void initializeLegendreExpansion(TF1* a_funct, TGraphAsymmErrors* a_graph, double a_lowRange, double a_highRange, bool a_fixRange = false, int a_numIntegrationSteps = -1);
  TF1* determineGoodLegendreOrder(TGraphAsymmErrors* a_graph, double a_lowRange = -999.999, double a_highRange = -999.999, int a_lowOrderRange = 0, int a_highOrderRange = 30, int a_leaveOutEvery = 6); // this will not include Xerror in chi squared, returns best option for legendre function. LeaveOutEvery, deletes about every nth point




  #ifdef _MAC_OSX_
  void testLambdaDecay(); // by rapidity mtm0
  void testNewLambdaDecay(); //by 4D vectors  // WITH DCA INFO
  void testXiMinusDecay();
  void testNewXiMinusDecay(); //by 4D vectors
  void testOmegaMinusDecay();
  void testNewOmegaMinusDecay(); //by 4D vectors


  #endif

  void testLegendreExpansionInitialization(int a_order = 13, bool a_freeLegendreBoundsInFit = true);
  void testDetermineGoodLegendreOrder();

}

#endif
