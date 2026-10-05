// DetectorParameters.hh
//
// A single thick Mars-regolith slab (not a stack of scoring layers like the
// dose-vs-depth app -- this app doesn't need per-layer scoring, since the
// quantity of interest (depth of the primary neutron's FIRST interaction)
// is read directly from the step's global position, see SteppingAction.cc).
//
// SLAB THICKNESS: set generously (400 g/cm^2) to comfortably exceed the
// expected neutron mean free path in rock-like material across the whole
// 10 MeV - 200 GeV range (order 30-60 g/cm^2 for the nuclear interaction
// length alone, and elastic scattering -- which also counts as "the first
// interaction" here -- typically has a similar or shorter mean free path,
// see SteppingAction.cc). RunAction reports the fraction of primaries that
// transit the ENTIRE slab without interacting at each energy; if that
// fraction isn't small/negligible, the slab wasn't thick enough at that
// energy and the reported mean depth is a biased (truncated, i.e. too
// short) underestimate of the true mean free path -- check that report
// before trusting the plot, don't just take the numbers on faith.

#ifndef DetectorParameters_h
#define DetectorParameters_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

namespace MarsNeutron {

const G4double kRegolithDensity = 1.52 * g / cm3; // same nominal density as the other Mars regolith apps this session

const G4double kSlabArealThickness = 400.0 * g / cm2; // see header note above
const G4double kSlabHalfWidth = 100.0 * cm;           // transverse half-width (broad slab, minimal edge leakage)

inline G4double SlabThicknessCm(){
  return (kSlabArealThickness / kRegolithDensity) / cm;
}

const char* const kSlabVolumeName = "RegolithSlabLV";

} // namespace MarsNeutron

#endif
