// DetectorParameters.hh
//
// Shared geometry/material constants used by both DetectorConstruction and
// RunAction (RunAction needs the layer mass to convert accumulated energy
// deposit into Gy/Sv, without re-querying the geometry at run time).
//
// STATUS: none of this has been compiled or run -- see the top-level README
// message delivered with this application for the full verification-status
// summary. The one piece that WAS checked against the user's actual
// installed Geant4 v11.3.2 headers/source is the ion-physics-list swap in
// main.cc (see comments there); everything else follows standard, stable
// Geant4 example patterns (G4Material mixtures, G4Accumulable, sensitive
// stepping-action dose scoring) that have been essentially unchanged for
// many major Geant4 versions.

#ifndef DetectorParameters_h
#define DetectorParameters_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

namespace MarsGCR {

// Regolith depth scan: 0 to 100 g/cm^2, matching the x-axis range of the
// borrowed Norbury & Slaba (2019) aluminum-shielding figure on Andrew's
// motivation slide, so the new regolith plot is a direct visual swap-in.
const G4int    kNLayers            = 50;
const G4double kLayerArealDensity  = 2.0 * g / cm2;  // areal density per layer
const G4double kTotalArealDensity  = kNLayers * kLayerArealDensity; // 100 g/cm^2

// Mars regolith bulk density -- Pathfinder nominal (see MarsGCRPenetrationDepth.cpp
// header for the sourced density range; this Geant4 app uses the single nominal
// value only, unlike the CSDA calc's low/nominal/high band, to keep the already
// large Monte Carlo campaign (8 species x 50 energies x 2 physics lists) bounded).
const G4double kRegolithDensity = 1.52 * g / cm3;

// Transverse extent of the regolith slab and of the broad parallel beam's
// transverse sampling disk. The slab is made much wider than the beam disk
// so secondary-particle showers (mainly spread by nuclear cascade neutrons)
// stay well inside the slab and don't leak out the sides before being
// scored -- a standard broad-beam depth-dose approximation.
const G4double kSlabHalfWidth   = 50.0 * cm;
const G4double kBeamDiskRadius  = 20.0 * cm;

// Derived: physical thickness of one layer and of the full stack, from the
// areal density and the chosen bulk density.
inline G4double LayerThicknessCm(){
  return (kLayerArealDensity / kRegolithDensity) / cm; // returns a plain number in cm
}
inline G4double TotalThicknessCm(){
  return kNLayers * LayerThicknessCm();
}

// Layer mass, needed by RunAction to convert accumulated MeV of energy
// deposit into Gy (J/kg). Same for every layer since all layers are
// identical in size/material.
inline G4double LayerMassKg(){
  G4double volume_cm3 = (2.0*kSlabHalfWidth/cm) * (2.0*kSlabHalfWidth/cm) * LayerThicknessCm();
  G4double mass_g = volume_cm3 * (kRegolithDensity / (g/cm3));
  return mass_g / 1000.0; // kg
}

// Volume naming convention used by both DetectorConstruction (to name the
// placed physical volumes) and SteppingAction (to recognize which layer a
// step occurred in via the placement copy number -- see SteppingAction.cc).
const char* const kLayerNamePrefix = "RegolithLayer_";

} // namespace MarsGCR

#endif
