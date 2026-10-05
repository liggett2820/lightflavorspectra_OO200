// main.cc -- MarsGCRDoseSim
//
// Runs the full GCR dose-vs-regolith-depth campaign: for each of the 8
// reference-ion species x 50 log-spaced energies (MarsGCR::GCRReferenceIons(),
// MarsGCR::GCREnergyGrid_MeVPerNucleon() -- see GCRSpectrum.hh), fires
// nPrimariesPerPoint primaries through the 100 g/cm^2 regolith stack and
// records the per-layer dose/dose-equivalent to a ROOT ntuple. Run once with
// each physics-list variant (QMD, INCLXX) to reproduce the same two-model
// comparison as the Geant4 curves on the existing motivation slide, now for
// Mars regolith instead of aluminum.
//
// USAGE:
//   ./MarsGCRDoseSim <QMD|INCLXX> [outputFile.root] [nPrimariesPerPoint]
//
// Example (defaults: MarsGCRDose_QMD.root, 1000 primaries/point):
//   ./MarsGCRDoseSim QMD
//   ./MarsGCRDoseSim INCLXX MarsGCRDose_INCLXX.root 2000
//
// COMPUTE COST WARNING: this is 8 species x 50 energies = 400 separate Runs
// per physics-list variant, i.e. 400 * nPrimariesPerPoint total primaries,
// each individually tracked (with all its secondaries) through up to
// 100 g/cm^2 of regolith with a high-precision (_HP) hadronic physics list.
// This is meaningfully slower than the standalone Bethe-Bloch calculation --
// expect this to take a while even at modest nPrimariesPerPoint. Start with
// a small value (e.g. 200-500) to confirm the build/run works end to end and
// check the output looks sane before committing to a long, high-statistics
// run. This app is SERIAL (not multithreaded) -- see RunAction.hh for what
// switching to MT would require.
//
// PHYSICS-LIST CONSTRUCTION -- THE ONE PIECE VERIFIED AGAINST THE ACTUAL
// INSTALLED GEANT4 v11.3.2 HEADERS/SOURCE (not just general Geant4
// knowledge): FTFP_BERT_HP's constructor registers "new G4IonPhysics(ver)"
// (physics TYPE bIons) as its ion-ion inelastic physics constructor.
// G4IonQMDPhysics and G4IonINCLXXPhysics both also set physics type bIons
// (confirmed by reading
// source/physics_lists/constructors/ions/src/G4Ion{QMD,INCLXX}Physics.cc),
// and G4VModularPhysicsList::ReplacePhysics() matches by physics TYPE, not
// by name (confirmed by reading source/run/src/G4VModularPhysicsList.cc) --
// so even though G4IonQMDPhysics's own default name ("IonQMD") differs from
// G4IonPhysics's ("ionInelasticFTFP_BIC"), ReplacePhysics() still correctly
// swaps it in. This is why the pattern below is
// physicsList->ReplacePhysics(new G4IonQMDPhysics()) rather than
// RegisterPhysics (which would just add a second, conflicting ion-physics
// constructor alongside the default one).

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "GCRSpectrum.hh"

#include "G4RunManagerFactory.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include "FTFP_BERT_HP.hh"
#include "G4IonQMDPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4VModularPhysicsList.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char** argv){
  if(argc < 2){
    std::cerr << "Usage: " << argv[0] << " <QMD|INCLXX> [outputFile.root] [nPrimariesPerPoint]\n";
    return 1;
  }
  std::string variantArg = argv[1];
  G4int physicsVariantCode; // 0 = QMD, 1 = INCLXX -- matches RunAction ntuple column
  if(variantArg == "QMD"){
    physicsVariantCode = 0;
  } else if(variantArg == "INCLXX"){
    physicsVariantCode = 1;
  } else {
    std::cerr << "First argument must be exactly \"QMD\" or \"INCLXX\", got \""
              << variantArg << "\"\n";
    return 1;
  }

  std::string outputFile = (argc >= 3) ? argv[2] : ("MarsGCRDose_" + variantArg + ".root");
  G4int nPrimariesPerPoint = (argc >= 4) ? std::atoi(argv[3]) : 1000;

  // --- Run manager: serial only, see RunAction.hh header comment ---
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  runManager->SetUserInitialization(new DetectorConstruction());

  // --- Physics list: FTFP_BERT_HP with the ion-ion physics swapped to the
  //     requested variant. See file header for the verification note. ---
  auto* physicsList = new FTFP_BERT_HP(1); // verbose=1
  if(physicsVariantCode == 0){
    physicsList->ReplacePhysics(new G4IonQMDPhysics());
    G4cout << ">>> Ion physics: G4IonQMDPhysics (QMD)" << G4endl;
  } else {
    physicsList->ReplacePhysics(new G4IonINCLXXPhysics());
    G4cout << ">>> Ion physics: G4IonINCLXXPhysics (INCLXX)" << G4endl;
  }
  runManager->SetUserInitialization(physicsList);

  auto* actionInit = new ActionInitialization();
  runManager->SetUserInitialization(actionInit);

  runManager->Initialize(); // triggers ActionInitialization::Build()

  PrimaryGeneratorAction* primaryGen = actionInit->GetPrimaryGeneratorAction();
  RunAction* runAction = actionInit->GetRunAction();
  if(!primaryGen || !runAction){
    std::cerr << "FATAL: ActionInitialization did not produce the expected user actions.\n";
    return 1;
  }

  // --- Analysis manager: one ROOT file for the whole campaign (all 400
  //     species x energy Runs write into the same ntuple, tagged by
  //     Z/A/energy/physicsVariantCode -- see RunAction::EndOfRunAction). ---
  auto* analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(outputFile);
  RunAction::BookNtuple();

  std::vector<MarsGCR::GCRIon> ions = MarsGCR::GCRReferenceIons();
  std::vector<G4double> energyGrid = MarsGCR::GCREnergyGrid_MeVPerNucleon();

  G4int nDone = 0;
  G4int nTotal = (G4int)(ions.size() * energyGrid.size());
  for(const auto& ion : ions){
    for(G4double energyMeVPerU : energyGrid){
      primaryGen->SetIon(ion.Z, ion.A);
      primaryGen->SetKineticEnergyPerNucleon(energyMeVPerU * MeV);
      runAction->SetCurrentTag(ion.Z, ion.A, energyMeVPerU, physicsVariantCode);

      runManager->BeamOn(nPrimariesPerPoint);

      nDone++;
      G4cout << ">>> Campaign progress: " << nDone << "/" << nTotal
             << "  (" << ion.name << ", " << energyMeVPerU << " MeV/u)" << G4endl;
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();

  delete runManager;
  std::cout << "Wrote " << outputFile << " (" << nTotal << " runs x "
            << nPrimariesPerPoint << " primaries)" << std::endl;
  return 0;
}
