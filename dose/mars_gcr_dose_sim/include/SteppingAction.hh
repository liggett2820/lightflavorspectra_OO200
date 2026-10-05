// SteppingAction.hh
//
// For every step that deposits energy inside a "RegolithLayer_*" volume,
// computes a local LET estimate (energy deposit / step length within that
// layer) and applies the ICRP 60/103 quality factor Q(L) to get the
// dose-equivalent contribution, then adds both the raw energy deposit and
// the Q-weighted deposit to the current Run's per-layer accumulators (via
// RunAction::AddStep). Runs for every track -- primary ions AND all
// secondaries (delta rays, spallation fragments, neutrons, ...) -- so the
// scored dose includes the full nuclear-interaction shower, not just direct
// primary ionization.

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class RunAction;
class G4Step;

class SteppingAction : public G4UserSteppingAction
{
public:
  explicit SteppingAction(RunAction* runAction);
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

private:
  RunAction* m_runAction;

  // ICRP 60/103 quality factor Q(L), L in keV/um (unrestricted LET in
  // water -- approximated here as LET in the regolith itself, not
  // re-derived in water, see .cc for the caveat this introduces).
  static G4double QualityFactor(G4double LET_keV_per_um);
};

#endif
