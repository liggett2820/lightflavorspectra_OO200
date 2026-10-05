// PrimaryGeneratorAction.hh
//
// Fires a broad, parallel, mono-energetic beam of one GCR ion species
// straight into the regolith stack's front face. Species/energy are set
// from main.cc via SetIon()/SetKineticEnergyPerNucleon() before each
// beamOn() call in the species x energy campaign loop -- NOT sampled
// randomly per-event, since each Geant4 Run in this campaign already
// corresponds to exactly one (species, energy) point (see main.cc).

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* anEvent) override;

  // Z = atomic number, A = integer mass number (matches MarsGCR::GCRIon).
  void SetIon(G4int Z, G4int A);
  // Kinetic energy PER NUCLEON; total gun energy = E_per_nucleon * A.
  void SetKineticEnergyPerNucleon(G4double energyPerNucleon);

private:
  G4ParticleGun* m_particleGun;
  G4int m_A; // cached mass number, needed to convert per-nucleon -> total KE
};

#endif
