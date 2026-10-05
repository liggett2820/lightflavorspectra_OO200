#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class G4Step;

// Watches ONLY the primary neutron's track (trackID==1). The first step
// inside the regolith slab whose post-step process is anything other than
// "Transportation" (i.e. a real physics interaction actually fired -- could
// be elastic scattering, inelastic/nuclear reaction, or capture) is that
// primary's first interaction: its depth is reported to EventAction and the
// primary track is killed immediately (we only need the FIRST interaction
// depth -- no reason to keep transporting the primary, or to score anything
// from its secondaries, for this measurement).
class SteppingAction : public G4UserSteppingAction
{
public:
  explicit SteppingAction(EventAction* eventAction);
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

private:
  EventAction* m_eventAction;
};

#endif
