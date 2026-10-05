#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorParameters.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

SteppingAction::SteppingAction(RunAction* runAction)
  : G4UserSteppingAction(), m_runAction(runAction) {}

G4double SteppingAction::QualityFactor(G4double LET_keV_per_um){
  // ICRP 60 (1991), retained unchanged in ICRP 103 (2007) / 127 (2014) --
  // exact constants verified against icrpaedia.org's "Quality factor" page:
  //   Q(L) = 1               for L < 10 keV/um
  //   Q(L) = 0.32*L - 2.2    for 10 <= L <= 100 keV/um
  //   Q(L) = 300 / sqrt(L)   for L > 100 keV/um
  //
  // CAVEAT: the ICRP definition is for *unrestricted LET in water*. Here L
  // is computed from the energy deposit and step length IN THE REGOLITH
  // ITSELF, not re-expressed as an equivalent water LET. Using local-medium
  // LET as a stand-in for water LET is a common practical shortcut in quick
  // Monte Carlo dosimetry, not the rigorous ICRP prescription -- flagged
  // here rather than silently applied, since regolith's stopping power
  // differs from water's (see MarsGCRPenetrationDepth.cpp's Bragg-mixture
  // Z/A and I comparisons, e.g. Z/A=0.493 vs water's ~0.555). The effect on
  // Q is second-order (Q depends on L only logarithmically/as a fractional
  // power in the relevant range) but it is a real, unquantified systematic,
  // not zero.
  if(LET_keV_per_um < 10.0) return 1.0;
  if(LET_keV_per_um <= 100.0) return 0.32 * LET_keV_per_um - 2.2;
  return 300.0 / std::sqrt(LET_keV_per_um);
}

void SteppingAction::UserSteppingAction(const G4Step* step){
  G4double edep = step->GetTotalEnergyDeposit(); // Geant4 internal units
  if(edep <= 0.0) return;

  const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume* pv = touchable->GetVolume();
  if(!pv) return;
  G4LogicalVolume* lv = pv->GetLogicalVolume();
  if(lv->GetName() != "RegolithLayerLV") return; // only score inside the regolith stack

  G4int layerIndex = touchable->GetCopyNumber();

  G4double stepLength = step->GetStepLength(); // Geant4 internal units
  G4double edep_MeV = edep / CLHEP::MeV;
  G4double LET_keV_per_um = 0.0;
  if(stepLength > 0.0){
    // MeV/mm and keV/um are numerically identical (both are "1000 keV per
    // 1000 um" scaled the same way), so once edep and stepLength are each
    // expressed in Geant4's own MeV/mm internal units, their ratio IS the
    // LET in keV/um directly -- no extra conversion factor needed.
    G4double stepLength_mm = stepLength / CLHEP::mm;
    LET_keV_per_um = edep_MeV / stepLength_mm;
  }

  G4double Q = QualityFactor(LET_keV_per_um);
  G4double doseEqWeighted_MeV = edep_MeV * Q;

  m_runAction->AddStep(layerIndex, edep_MeV, doseEqWeighted_MeV);
}
