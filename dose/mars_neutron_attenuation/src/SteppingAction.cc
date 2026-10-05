#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorParameters.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
  : G4UserSteppingAction(), m_eventAction(eventAction) {}

void SteppingAction::UserSteppingAction(const G4Step* step){
  G4Track* track = step->GetTrack();
  if(track->GetTrackID() != 1) return; // only the primary neutron matters here

  const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume* pv = touchable->GetVolume();
  if(!pv || pv->GetLogicalVolume()->GetName() != MarsNeutron::kSlabVolumeName) return;

  const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  if(!process) return;
  if(process->GetProcessName() == "Transportation") return; // just a geometry step, not a real interaction

  // Real physics interaction -- this is the primary's first interaction
  // (EventAction::SetInteractionDepth ignores any subsequent calls, but we
  // also kill the track right here so there ISN'T a subsequent call).
  G4double depth_cm = step->GetPostStepPoint()->GetPosition().z() / cm;

  m_eventAction->SetInteractionDepth(depth_cm);
  track->SetTrackStatus(fStopAndKill);
}
