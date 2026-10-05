// RunAction.hh
//
// Accumulates per-layer energy deposit (for absorbed dose, Gy) and
// LET-weighted energy deposit (for dose equivalent, Sv, via the ICRP 60/103
// Q(L) quality factor -- see SteppingAction.cc) across one Run's worth of
// primaries, then writes one ROOT ntuple row per layer at EndOfRunAction,
// tagged with which (species, energy, physics-list variant) this Run was.
//
// SIMPLIFIED FOR SERIAL EXECUTION ONLY: this uses plain G4double arrays,
// reset at BeginOfRunAction and read at EndOfRunAction, which is only
// correct with a SERIAL G4RunManager (see main.cc -- G4RunManagerType::Serial
// is the default here). If you switch to a multithreaded run manager, these
// arrays are NOT thread-safe as written and would need to become
// G4Accumulable<G4double> registered with G4AccumulableManager instead
// (Reset() in BeginOfRunAction, Merge() in EndOfRunAction) -- flagged here
// rather than attempted, to keep this first version simple and more likely
// to build correctly without a compiler available to check it against.

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

class RunAction : public G4UserRunAction
{
public:
  RunAction();
  ~RunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  // Called from SteppingAction for every step depositing energy in a
  // regolith layer.
  void AddStep(G4int layerIndex, G4double edep_MeV, G4double doseEqWeighted_MeV);

  // Called from main.cc before each beamOn() so EndOfRunAction knows what
  // to tag this run's ntuple rows with. physicsVariantCode: 0 = QMD ion
  // physics, 1 = INCLXX ion physics (see main.cc) -- an integer code rather
  // than a string so the ntuple only needs I/D columns (species is likewise
  // identified by (Z,A) rather than a name string) -- avoids any dependency
  // on G4AnalysisManager's string-ntuple-column support, which is newer and
  // less certain to be present/working than the basic I/D columns.
  void SetCurrentTag(G4int Z, G4int A, G4double energyMeVPerNucleon,
                      G4int physicsVariantCode);

  // Must be called once, after G4AnalysisManager::Instance()->OpenFile(...),
  // before any beamOn() -- creates the ntuple columns.
  static void BookNtuple();

private:
  std::vector<G4double> m_edep_MeV;          // per layer, this run
  std::vector<G4double> m_doseEqWeighted_MeV; // per layer, this run

  G4int    m_Z = 0;
  G4int    m_A = 0;
  G4double m_energyMeVPerNucleon = 0.0;
  G4int    m_physicsVariantCode = 0;
};

#endif
