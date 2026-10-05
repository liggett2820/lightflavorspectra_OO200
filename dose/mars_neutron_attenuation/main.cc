// main.cc -- MarsNeutronAttenuation
//
// For each of the 50 log-spaced neutron energies (10 MeV - 200 GeV, same
// grid as the other Mars GCR calculations this session), fires
// nEventsPerPoint mono-energetic neutrons into a 400 g/cm^2 Mars regolith
// slab and records the depth of each primary's FIRST real interaction
// (elastic, inelastic, or capture -- see SteppingAction.cc). The mean of
// those depths is the Monte Carlo estimate of the neutron's attenuation
// length / mean free path at that energy -- the physically correct
// "penetration depth" quantity for a neutral particle (see the AskUserQuestion
// exchange this was scoped from: a neutron doesn't have a CSDA range the way
// a charged particle does).
//
// USAGE:
//   ./MarsNeutronAttenuation [outputFile.root] [nEventsPerPoint]
// Example (defaults: MarsNeutronAttenuation.root, 3000 events/point):
//   ./MarsNeutronAttenuation
//   ./MarsNeutronAttenuation MarsNeutronAttenuation.root 5000
//
// NO ion-physics variant choice here (unlike mars_gcr_dose_sim's QMD vs.
// INCLXX): G4IonQMDPhysics/G4IonINCLXXPhysics only affect ION-nucleus
// (nucleus-nucleus) interactions, i.e. when the PROJECTILE is an ion.
// Neutron-nucleus interactions are handled by FTFP_BERT_HP's own neutron/
// hadron physics (G4HadronPhysicsFTFP_BERT_HP + the _HP high-precision
// neutron data library below 20 MeV), which this app uses unmodified --
// no ReplacePhysics() call needed or appropriate here.
//
// COMPUTE COST: much lighter than mars_gcr_dose_sim -- only 50 runs total
// (not 400), and each primary is killed at its first interaction rather
// than tracked to completion, so this should run considerably faster. Still
// worth a small first test (a few hundred events/point) before committing
// to the full statistics run.

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "NeutronEnergyGrid.hh"

#include "G4RunManagerFactory.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include "FTFP_BERT_HP.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char** argv){
  std::string outputFile = (argc >= 2) ? argv[1] : "MarsNeutronAttenuation.root";
  G4int nEventsPerPoint = (argc >= 3) ? std::atoi(argv[2]) : 3000;

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new FTFP_BERT_HP(1)); // unmodified -- see file header

  auto* actionInit = new ActionInitialization();
  runManager->SetUserInitialization(actionInit);

  runManager->Initialize();

  PrimaryGeneratorAction* primaryGen = actionInit->GetPrimaryGeneratorAction();
  RunAction* runAction = actionInit->GetRunAction();
  if(!primaryGen || !runAction){
    std::cerr << "FATAL: ActionInitialization did not produce the expected user actions.\n";
    return 1;
  }

  auto* analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(outputFile);
  RunAction::BookNtuple();

  std::vector<G4double> energyGrid = MarsNeutron::EnergyGrid_MeV();

  G4int nDone = 0;
  G4int nTotal = (G4int) energyGrid.size();
  for(G4double energyMeV : energyGrid){
    primaryGen->SetKineticEnergy(energyMeV * MeV);
    runAction->SetCurrentEnergy(energyMeV);

    runManager->BeamOn(nEventsPerPoint);

    nDone++;
    G4cout << ">>> Progress: " << nDone << "/" << nTotal
           << "  (" << energyMeV << " MeV)" << G4endl;
  }

  analysisManager->Write();
  analysisManager->CloseFile();

  delete runManager;
  std::cout << "Wrote " << outputFile << " (" << nTotal << " energy points x "
            << nEventsPerPoint << " primaries)" << std::endl;
  return 0;
}
