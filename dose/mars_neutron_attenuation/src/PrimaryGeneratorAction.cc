#include "PrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Neutron.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), m_particleGun(nullptr)
{
  m_particleGun = new G4ParticleGun(1);
  m_particleGun->SetParticleDefinition(G4Neutron::Neutron());
  m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  m_particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -1.0 * cm));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  delete m_particleGun;
}

void PrimaryGeneratorAction::SetKineticEnergy(G4double energy){
  m_particleGun->SetParticleEnergy(energy);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  m_particleGun->GeneratePrimaryVertex(anEvent);
}
