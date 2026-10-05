#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;
class G4Event;

// Per-event state for "has the primary neutron's first interaction already
// been recorded for this event". SteppingAction sets the depth via
// SetInteractionDepth() as soon as it sees the primary's first real
// (non-Transportation) physics step; EndOfEventAction reports the result
// (interacted at depth X, or transmitted without interacting) to RunAction.
class EventAction : public G4UserEventAction
{
public:
  explicit EventAction(RunAction* runAction);
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void SetInteractionDepth(G4double depth_cm);

private:
  RunAction* m_runAction;
  G4bool m_hasInteracted;
  G4double m_interactionDepth_cm;
};

#endif
