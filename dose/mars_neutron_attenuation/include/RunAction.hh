// RunAction.hh
//
// One Run = one incident neutron energy. Accumulates the first-interaction
// depth (cm) of every primary that interacted within the slab (reported
// by EventAction at the end of each event), plus a count of primaries that
// transited the whole slab without interacting (see DetectorParameters.hh
// for why that count matters -- it's the "was the slab thick enough" check).
// At EndOfRunAction, writes ONE ntuple row summarizing this energy point:
// mean first-interaction depth (= the Monte Carlo estimate of the
// attenuation length/mean free path), its standard error, and the
// transmitted (non-interacting) fraction.
//
// SERIAL EXECUTION ONLY -- same caveat as mars_gcr_dose_sim/RunAction.hh:
// plain (non-thread-safe) accumulation, correct only with
// G4RunManagerType::Serial (see main.cc).

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

  // Called once per event, from EventAction::EndOfEventAction.
  void RecordInteraction(G4double depth_cm);
  void RecordTransmitted(); // primary reached the back of the slab without interacting

  void SetCurrentEnergy(G4double energyMeV);

  static void BookNtuple();

private:
  std::vector<G4double> m_interactionDepths_cm;
  G4int m_nTransmitted = 0;
  G4double m_energyMeV = 0.0;
};

#endif
