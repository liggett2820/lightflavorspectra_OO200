// GCRSpectrum.hh
//
// Shared GCR reference-ion table for the Geant4 dose campaign. Same 8
// species and same 50-point, 10 MeV/u - 200 GeV/u log-spaced energy grid as
// MarsGCRPenetrationDepth.cpp (dose/mars_gcr_penetration/), so the Bethe-Bloch
// CSDA calculation and this Geant4 dose calculation are directly comparable
// point-for-point.
//
// RELATIVE ABUNDANCE WEIGHTS -- SOURCING AND APPROXIMATION, READ BEFORE USE:
//   Group-level values are from the Particle Data Group Cosmic Rays review,
//   Table 30.1 (https://pdg.lbl.gov/2022/reviews/rpp2022-rev-cosmic-rays.pdf),
//   relative abundances at 10.6 GeV/nucleon, normalized so oxygen's own flux
//   = 1.00:
//     H = 550, He = 34, "C-O" group = 2.20, "Ne-F" group = 0.30,
//     "Na-Mg" group = 0.22, "Al-Si" group = 0.19, "Fe-Ni" group = 0.12
//   These are REAL, cited numbers. But four of our eight species (C, Ne, Mg,
//   Si, Fe) share a table row with a neighboring element the table doesn't
//   resolve individually. Splitting each group is NOT sourced from the PDG
//   table -- it's my own approximation, done as follows:
//     - C vs O: the table's own O=1.00 normalization is used as the anchor,
//       so C = (group total) - 1.00 = 2.20 - 1.00 = 1.20. This is the one
//       split that follows directly from the table's stated normalization
//       rather than an assumed ratio.
//     - Ne vs F, Mg vs Na, Si vs Al, Fe vs Ni: all four pairs are an
//       even-Z/odd-Z pair. Odd-Z cosmic-ray nuclei are well known to be
//       suppressed relative to their even-Z neighbors (nuclear stability --
//       the same effect visible in solar-system elemental abundances), so
//       I assign the even-Z (more abundant) member 85-95% of the group's
//       flux, heavier weighting for the pairs with a bigger stability gap:
//         Ne = 0.90 x 0.30 = 0.27   (F is a pure spallation secondary, no
//                                    primary source -- strongest suppression)
//         Mg = 0.85 x 0.22 = 0.187
//         Si = 0.85 x 0.19 = 0.1615
//         Fe = 0.95 x 0.12 = 0.114  (real Fe/Ni particle ratios run ~20-25,
//                                    i.e. Ni is a small fraction of the group)
//       These split fractions (0.90/0.85/0.85/0.95) are MY judgment calls,
//       not measured numbers -- flagged here so they're easy to replace with
//       better values later without hunting back through this file's logic.
//
// ENERGY-SPECTRUM WEIGHTING -- DELIBERATELY NOT DONE HERE:
//   This table only weights species RELATIVE TO EACH OTHER. It does *not*
//   weight the 50 energy grid points relative to each other -- every energy
//   point for a given species gets equal Monte Carlo statistics here. A real
//   GCR differential energy spectrum is strongly peaked around a few hundred
//   MeV/u to a few GeV/u and falls off at both ends, so an equal-weight
//   average over energy is a real simplification, not just a statistics
//   choice -- it will not reproduce a true GCR-weighted dose rate by itself.
//   This is intentional: energy-spectrum weighting is applied in the ROOT
//   post-processing macro (AnalyzeMarsGCRDose.C), not baked into the Monte
//   Carlo, specifically so it can be revised later without re-running the
//   (expensive) Geant4 campaign. See that macro's header for the default
//   weighting scheme used and how to change it.

#ifndef GCRSpectrum_h
#define GCRSpectrum_h 1

#include "globals.hh"
#include <vector>
#include <cmath>

namespace MarsGCR {

struct GCRIon {
  G4String name;
  G4int    Z;
  G4int    A;          // integer mass number, for G4IonTable::GetIon(Z,A,0.0)
  G4double relAbundance; // see file header -- O-normalized, NOT a probability
};

inline std::vector<GCRIon> GCRReferenceIons(){
  return {
    {"H",  1,  1,  550.0},
    {"He", 2,  4,   34.0},
    {"C",  6,  12,   1.20},
    {"O",  8,  16,   1.00},
    {"Ne", 10, 20,   0.27},
    {"Mg", 12, 24,   0.187},
    {"Si", 14, 28,   0.1615},
    {"Fe", 26, 56,   0.114},
  };
}

// Same grid as MarsGCRPenetrationDepth.cpp: 10 MeV/nucleon to 200 GeV/nucleon,
// 50 log-spaced points.
inline std::vector<G4double> GCREnergyGrid_MeVPerNucleon(){
  std::vector<G4double> grid;
  const int N = 50;
  const double eLo = 10.0, eHi = 200000.0;
  for(int i = 0; i < N; i++){
    double frac = (double) i / (N - 1);
    grid.push_back(eLo * std::pow(eHi / eLo, frac));
  }
  return grid;
}

} // namespace MarsGCR

#endif
