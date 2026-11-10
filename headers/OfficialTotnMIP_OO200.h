#ifndef __OFFICIAL_TOTNMIP_OO200_H__
#define __OFFICIAL_TOTNMIP_OO200_H__
//______________________________________________________________________________
// Official StRefMultCorr("totnmip") correction for O+O 200 GeV, Run 21 (calibrated on the min-bias
// triggers 860001, 860002, 860011, 860012), as the centrality group asked on Nov 10, 2026: use the
// official min-bias totnMIP centrality classes for the 860003 analysis.
//
// Copied from star-sw StRoot/StRefMultCorr, commit 83175f6ade (Mar 12 2025, "Update centrality
// definition (#735)"): Param.h mParamStr_ref7 (parameter index 0) and oo200_run21_vzCorr_totnMIP;
// StRefMultCorr.cxx initEvent / getRefMultCorr / sampleRefMult / vzCorrection /
// getVzWindowForVzDepCentDef / isZvertexOk / isCentralityOk. It reproduces
//   StRefMultCorr* c = new StRefMultCorr("totnmip"); c->init(runId);
//   c->initEvent((UShort_t) totnMIP, vz, zdcx); c->getRefMultCorr();
// without linking the library (this build carries no StRefMultCorr):
//   corr = ((UShort_t) totnMIP + U(0,1)) * c_vz,  c_vz in 6 cm Vz bins over -30 <= Vz < 30;
//   no luminosity correction (official parameters 0, 0);
//   returned uncorrected (no sampling) for a run outside 22130029-22144006 or Vz outside (-30, 30).
// Raw totnMIP is the same definition as headers/EpdCentrality_OO200.h (EpdCentOO200::truncMip).
//
// Official class edges (lower edges, 80% ... 5%): 37, 46, 55, 66, 77, 91, 105, 122, 141, 162, 186,
// 213, 242, 276, 315, 367. The six analysis classes use 367 / 315 / 242 / 141 / 37 (set in
// macros/SetCutClass.C); below 37 is 80-100%, which the official calibration does not cover.
// Written for ROOT 5 (C++98).
//______________________________________________________________________________
#include "TRandom.h"

namespace OfficialTotnMIPOO200 {
  const int    kRunStart = 22130029;
  const int    kRunStop  = 22144006;
  const double kVzMin    = -30.0;
  const double kVzMax    =  30.0;
  const int    kNVz      = 10;
  const double kVzLow[kNVz]  = {-30, -24, -18, -12, -6, 0, 6, 12, 18, 24};   // each bin [low, low + 6)
  const double kVzCorr[kNVz] = {1.002648, 1.000983, 1.000385, 1.000295, 1.000004,
                                1.000097, 0.999605, 0.999488, 0.999480, 1.000674};
  const int    kNEdges = 16;
  const double kEdges[kNEdges] = {37, 46, 55, 66, 77, 91, 105, 122, 141, 162, 186, 213, 242, 276, 315, 367};

  inline double vzCorr(double a_vz){
    for(int i = 0; i < kNVz; i++) if(kVzLow[i] <= a_vz && a_vz < kVzLow[i] + 6.0) return kVzCorr[i];
    return 1.0;
  }
  // a_rawTotnMIP: sum of EpdCentOO200::truncMip(nMIP) over the event's EPD hits
  inline double correct(double a_rawTotnMIP, double a_vz, int a_runId){
    const unsigned short t = (unsigned short) a_rawTotnMIP;          // initEvent takes UShort_t
    if(a_runId < kRunStart || a_runId > kRunStop) return (double) t;  // isIndexOk() false
    if(!(a_vz > kVzMin && a_vz < kVzMax)) return (double) t;          // isZvertexOk() false
    return ((double) t + gRandom->Rndm()) * vzCorr(a_vz);            // sampleRefMult x vzCorrection
  }
}

#endif
