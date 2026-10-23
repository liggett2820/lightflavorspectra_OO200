#ifndef __EPD_CENTRALITY_OO200_H__
#define __EPD_CENTRALITY_OO200_H__
//______________________________________________________________________________
// Corrected EPD totnMIP for O+O 200 GeV (Run 21, trigger 860003): the centrality variable of the
// NBD-Glauber fit (centrality/star_glauber, fine scan Oct 2026).
//
// Raw totnMIP, as macros/MakeEpdCentInputs_OO200.C (= StRefMultCorr "totnmip" definition):
//   sum over all EPD hits of nMIP, with nMIP < 0.3 skipped and nMIP > 6 counted as 6.
// Corrections, as macros/FitEpdCent_OO200.C pass 1 (OO200_EpdCent_params.txt, Oct 2026):
//   c_vz   = <totnMIP>(|Vz| <= 2 cm) / <totnMIP>(Vz), mean method, 1 cm bins over |Vz| <= 30 cm
//            (1 outside), largest |c_vz - 1| = 0.35%;
//   c_lumi = p0 / (p0 + p1 * ZDCx[kHz]), p0 = 228.772, p1 = -0.762178 per kHz (-0.33% per kHz).
// FitEpdCent_OO200.C evaluated c_lumi at the centre of 0.2 kHz ZDCx bins; here it is evaluated at
// the event's ZDCx (difference below 0.04%).
// Written for ROOT 5 (C++98).
//______________________________________________________________________________
#include <cmath>

namespace EpdCentOO200 {

  const double kMipMin = 0.3;   // tiles below this are skipped
  const double kMipMax = 6.0;   // tiles above this count as 6

  const int    kNVz    = 60;    // 1 cm bins, Vz from -30 to 30 cm
  const double kVzMin  = -30.0;
  const double kVzCorr[kNVz] = {
    0.99813, 0.99847, 0.99750, 0.99902, 0.99881, 0.99750, 0.99848, 0.99934, 0.99930, 0.99885, // -30 .. -20
    0.99938, 0.99918, 0.99965, 0.99960, 0.99975, 0.99984, 0.99952, 0.99970, 0.99996, 0.99990, // -20 .. -10
    1.00024, 1.00017, 1.00026, 0.99985, 0.99981, 1.00014, 1.00018, 1.00006, 1.00011, 0.99989, // -10 ..   0
    0.99994, 1.00006, 0.99994, 1.00010, 0.99984, 0.99983, 0.99954, 0.99960, 0.99997, 0.99916, //   0 ..  10
    0.99985, 0.99932, 0.99913, 0.99946, 0.99962, 0.99953, 0.99847, 0.99921, 0.99968, 0.99888, //  10 ..  20
    0.99907, 0.99919, 0.99852, 0.99859, 0.99928, 0.99802, 0.99912, 0.99691, 0.99782, 0.99647  //  20 ..  30
  };

  const double kLumiP0 = 228.772;
  const double kLumiP1 = -0.762178; // per kHz

  inline double truncMip(double m){ if(m < kMipMin) return 0.0; return (m > kMipMax) ? kMipMax : m; }

  inline double vzCorr(double vz){
    const int i = (int) std::floor(vz - kVzMin);
    if(i < 0 || i >= kNVz) return 1.0;
    return kVzCorr[i];
  }

  // a_zdcxHz: StPicoEvent::ZDCx(), in Hz
  inline double lumiCorr(double a_zdcxHz){
    const double d = kLumiP0 + kLumiP1 * (a_zdcxHz / 1000.0);
    return (d > 0) ? kLumiP0 / d : 1.0;
  }

  inline double correct(double a_rawTotnMIP, double a_vz, double a_zdcxHz){
    return a_rawTotnMIP * vzCorr(a_vz) * lumiCorr(a_zdcxHz);
  }
}

#endif
