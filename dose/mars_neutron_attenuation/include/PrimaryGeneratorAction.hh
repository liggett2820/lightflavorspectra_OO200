#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

// Fires mono-energetic neutrons straight down +Z into the regolith slab's
// front face (z=0). Energy is set from main.cc via SetKineticEnergy()
// before each beamOn() call -- one Run per energy point, matching the
// pattern established in mars_gcr_dose_sim.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* anEvent) override;
  void SetKineticEnergy(G4double energy);

private:
  G4ParticleGun* m_particleGun;
};

#endif
