#include "EventAction.hh"
#include "RunAction.hh"

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(), m_runAction(runAction),
    m_hasInteracted(false), m_interactionDepth_cm(0.0) {}

void EventAction::BeginOfEventAction(const G4Event*){
  m_hasInteracted = false;
  m_interactionDepth_cm = 0.0;
}

void EventAction::SetInteractionDepth(G4double depth_cm){
  if(m_hasInteracted) return; // only the FIRST interaction counts
  m_hasInteracted = true;
  m_interactionDepth_cm = depth_cm;
}

void EventAction::EndOfEventAction(const G4Event*){
  if(m_hasInteracted){
    m_runAction->RecordInteraction(m_interactionDepth_cm);
  } else {
    m_runAction->RecordTransmitted();
  }
}
