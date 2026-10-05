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
  // Identical to mars_gcr_dose_sim/src/DetectorConstruction.cc's
  // BuildMarsRegolithMaterial() -- same source data, same construction via
  // Geant4's own G4Material compound API. Duplicated here (not shared via a
  // library) since this is a separate small standalone app; if you change
  // the regolith composition in one, change it in the other.
  G4Material* BuildMarsRegolithMaterial();
};

#endif
