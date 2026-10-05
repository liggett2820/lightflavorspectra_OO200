#include "PrimaryGeneratorAction.hh"
#include "DetectorParameters.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), m_particleGun(nullptr), m_A(1)
{
  m_particleGun = new G4ParticleGun(1); // 1 primary per event
  // Beam travels in +Z, starting just in front of the regolith stack's
  // front face (stack front face is at z = -TotalThicknessCm/2).
  m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  G4double zStart = -(MarsGCR::TotalThicknessCm() * cm) / 2.0 - 10.0 * cm;
  m_particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, zStart));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  delete m_particleGun;
}

void PrimaryGeneratorAction::SetIon(G4int Z, G4int A){
  m_A = A;
  // Ground-state ion (excitation energy 0). Requires the physics list's
  // ConstructParticle() to have already registered G4GenericIon / enabled
  // dynamic ion lookup -- true here since G4IonPhysics/G4IonQMDPhysics/
  // G4IonINCLXXPhysics (whichever is active) does this as part of building
  // the physics list in main.cc, which runs before RunManager::Initialize().
  G4ParticleDefinition* ionDef = G4IonTable::GetIonTable()->GetIon(Z, A, 0.0);
  m_particleGun->SetParticleDefinition(ionDef);
}

void PrimaryGeneratorAction::SetKineticEnergyPerNucleon(G4double energyPerNucleon){
  m_particleGun->SetParticleEnergy(energyPerNucleon * m_A);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  // Broad parallel beam: randomize the transverse (x,y) start position
  // uniformly over a disk of radius kBeamDiskRadius, well inside the slab's
  // transverse extent (kSlabHalfWidth) so shower spread stays contained.
  G4double r = MarsGCR::kBeamDiskRadius * std::sqrt(G4UniformRand());
  G4double phi = CLHEP::twopi * G4UniformRand();
  G4double x = r * std::cos(phi);
  G4double y = r * std::sin(phi);

  G4double zStart = -(MarsGCR::TotalThicknessCm() * cm) / 2.0 - 10.0 * cm;
  m_particleGun->SetParticlePosition(G4ThreeVector(x, y, zStart));
  m_particleGun->GeneratePrimaryVertex(anEvent);
}
