// DetectorConstruction.hh
//
// World volume (vacuum) containing a stack of MarsGCR::kNLayers thin Mars
// regolith slabs of identical material, stacked along +Z. Scoring dose per
// layer (in SteppingAction) as a beam of GCR ions passes straight through
// the whole stack gives the full depth-dose profile from ONE run per
// (species, energy) point, rather than needing a separate run per depth --
// standard technique for depth-dose Monte Carlo (the same approach used by,
// e.g., HZETRN and published Geant4 GCR-shielding depth-dose studies).

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4Material;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override = default;

  G4VPhysicalVolume* Construct() override;

private:
  G4Material* BuildMarsRegolithMaterial();
};

#endif
