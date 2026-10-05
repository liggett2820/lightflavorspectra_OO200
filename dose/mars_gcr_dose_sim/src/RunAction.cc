#include "RunAction.hh"
#include "DetectorParameters.hh"

#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>

namespace {
  const G4double kMeVtoJoule = 1.602176634e-13; // exact, CODATA MeV->J
}

RunAction::RunAction() : G4UserRunAction() {
  m_edep_MeV.assign(MarsGCR::kNLayers, 0.0);
  m_doseEqWeighted_MeV.assign(MarsGCR::kNLayers, 0.0);
}

void RunAction::SetCurrentTag(G4int Z, G4int A, G4double energyMeVPerNucleon,
                               G4int physicsVariantCode){
  m_Z = Z;
  m_A = A;
  m_energyMeVPerNucleon = energyMeVPerNucleon;
  m_physicsVariantCode = physicsVariantCode;
}

void RunAction::BeginOfRunAction(const G4Run*){
  std::fill(m_edep_MeV.begin(), m_edep_MeV.end(), 0.0);
  std::fill(m_doseEqWeighted_MeV.begin(), m_doseEqWeighted_MeV.end(), 0.0);
}

void RunAction::AddStep(G4int layerIndex, G4double edep_MeV, G4double doseEqWeighted_MeV){
  if(layerIndex < 0 || layerIndex >= MarsGCR::kNLayers) return; // defensive
  m_edep_MeV[layerIndex] += edep_MeV;
  m_doseEqWeighted_MeV[layerIndex] += doseEqWeighted_MeV;
}

void RunAction::BookNtuple(){
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->CreateNtuple("MarsGCRDose",
      "Per-layer absorbed dose and dose equivalent per primary, "
      "by species (Z,A) / energy / ion-physics variant");
  analysisManager->CreateNtupleIColumn("Z");
  analysisManager->CreateNtupleIColumn("A");
  analysisManager->CreateNtupleDColumn("energyMeVPerNucleon");
  analysisManager->CreateNtupleIColumn("physicsVariantCode"); // 0=QMD, 1=INCLXX
  analysisManager->CreateNtupleIColumn("layerIndex");
  analysisManager->CreateNtupleDColumn("layerDepthLo_cm");
  analysisManager->CreateNtupleDColumn("layerDepthHi_cm");
  analysisManager->CreateNtupleIColumn("nPrimaries");
  analysisManager->CreateNtupleDColumn("dose_Gy_per_primary");
  analysisManager->CreateNtupleDColumn("doseEq_Sv_per_primary");
  analysisManager->FinishNtuple();
}

void RunAction::EndOfRunAction(const G4Run* run){
  G4int nEvents = run->GetNumberOfEvent();
  if(nEvents <= 0) return; // nothing to write (shouldn't happen)

  G4double layerMassKg = MarsGCR::LayerMassKg();
  auto analysisManager = G4AnalysisManager::Instance();

  for(G4int i = 0; i < MarsGCR::kNLayers; i++){
    G4double dose_Gy_total    = (m_edep_MeV[i] * kMeVtoJoule) / layerMassKg;
    G4double doseEq_Sv_total  = (m_doseEqWeighted_MeV[i] * kMeVtoJoule) / layerMassKg;
    G4double dose_Gy_perPrim   = dose_Gy_total   / (G4double) nEvents;
    G4double doseEq_Sv_perPrim = doseEq_Sv_total / (G4double) nEvents;

    // Linear depth in cm (not areal density) -- layer thickness is derived
    // from the areal-density geometry parameter + regolith bulk density,
    // see DetectorParameters.hh::LayerThicknessCm().
    G4double layerThicknessCm = MarsGCR::LayerThicknessCm();
    G4double depthLo = i * layerThicknessCm;
    G4double depthHi = (i + 1) * layerThicknessCm;

    G4int col = 0;
    analysisManager->FillNtupleIColumn(col++, m_Z);
    analysisManager->FillNtupleIColumn(col++, m_A);
    analysisManager->FillNtupleDColumn(col++, m_energyMeVPerNucleon);
    analysisManager->FillNtupleIColumn(col++, m_physicsVariantCode);
    analysisManager->FillNtupleIColumn(col++, i);
    analysisManager->FillNtupleDColumn(col++, depthLo);
    analysisManager->FillNtupleDColumn(col++, depthHi);
    analysisManager->FillNtupleIColumn(col++, nEvents);
    analysisManager->FillNtupleDColumn(col++, dose_Gy_perPrim);
    analysisManager->FillNtupleDColumn(col++, doseEq_Sv_perPrim);
    analysisManager->AddNtupleRow();
  }
}
