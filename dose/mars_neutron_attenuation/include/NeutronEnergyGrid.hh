// NeutronEnergyGrid.hh
//
// Same 50-point, log-spaced, 10 MeV - 200 GeV energy grid used throughout
// this session's other Mars GCR calculations (MarsGCRPenetrationDepth.cpp,
// GCRSpectrum.hh), so this plot lines up point-for-point with those. For a
// neutron (A=1) "per nucleon" and "total" kinetic energy are the same thing,
// so this is just a plain neutron kinetic energy grid.

#ifndef NeutronEnergyGrid_h
#define NeutronEnergyGrid_h 1

#include "globals.hh"
#include <vector>
#include <cmath>

namespace MarsNeutron {

inline std::vector<G4double> EnergyGrid_MeV(){
  std::vector<G4double> grid;
  const int N = 50;
  const double eLo = 10.0, eHi = 200000.0;
  for(int i = 0; i < N; i++){
    double frac = (double) i / (N - 1);
    grid.push_back(eLo * std::pow(eHi / eLo, frac));
  }
  return grid;
}

} // namespace MarsNeutron

#endif
