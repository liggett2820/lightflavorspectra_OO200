#include "RunAction.hh"

#include "G4Run.hh"
#include "G4AnalysisManager.hh"

#include <cmath>

RunAction::RunAction() : G4UserRunAction() {}

void RunAction::SetCurrentEnergy(G4double energyMeV){
  m_energyMeV = energyMeV;
}

void RunAction::BeginOfRunAction(const G4Run*){
  m_interactionDepths_cm.clear();
  m_nTransmitted = 0;
}

void RunAction::RecordInteraction(G4double depth_cm){
  m_interactionDepths_cm.push_back(depth_cm);
}

void RunAction::RecordTransmitted(){
  m_nTransmitted++;
}

void RunAction::BookNtuple(){
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->CreateNtuple("NeutronAttenuation",
      "Mean first-interaction depth (attenuation length estimate) of a "
      "primary neutron in Mars regolith, per incident energy");
  analysisManager->CreateNtupleDColumn("energyMeV");
  analysisManager->CreateNtupleIColumn("nEvents");
  analysisManager->CreateNtupleIColumn("nInteracted");
  analysisManager->CreateNtupleIColumn("nTransmitted");
  analysisManager->CreateNtupleDColumn("meanFirstInteractionDepth_cm");
  analysisManager->CreateNtupleDColumn("stdErrOfMean_cm");
  analysisManager->CreateNtupleDColumn("transmittedFraction");
  analysisManager->FinishNtuple();
}

void RunAction::EndOfRunAction(const G4Run* run){
  G4int nEvents = run->GetNumberOfEvent();
  if(nEvents <= 0) return;

  G4int nInteracted = (G4int) m_interactionDepths_cm.size();

  G4double mean = 0.0;
  for(G4double d : m_interactionDepths_cm) mean += d;
  if(nInteracted > 0) mean /= nInteracted;

  G4double variance = 0.0;
  for(G4double d : m_interactionDepths_cm) variance += (d - mean) * (d - mean);
  G4double stdErr = 0.0;
  if(nInteracted > 1){
    variance /= (nInteracted - 1);
    stdErr = std::sqrt(variance / nInteracted);
  }

  G4double transmittedFraction = (G4double) m_nTransmitted / (G4double) nEvents;

  if(transmittedFraction > 0.02){
    G4cout << ">>> WARNING: at E=" << m_energyMeV << " MeV, "
           << (transmittedFraction * 100.0) << "% of primaries transited the "
           << "whole slab without interacting -- the slab may not be thick "
           << "enough at this energy; the mean depth below is likely a "
           << "truncated UNDERESTIMATE. See DetectorParameters.hh." << G4endl;
  }

  auto analysisManager = G4AnalysisManager::Instance();
  G4int col = 0;
  analysisManager->FillNtupleDColumn(col++, m_energyMeV);
  analysisManager->FillNtupleIColumn(col++, nEvents);
  analysisManager->FillNtupleIColumn(col++, nInteracted);
  analysisManager->FillNtupleIColumn(col++, m_nTransmitted);
  analysisManager->FillNtupleDColumn(col++, mean);
  analysisManager->FillNtupleDColumn(col++, stdErr);
  analysisManager->FillNtupleDColumn(col++, transmittedFraction);
  analysisManager->AddNtupleRow();
}
