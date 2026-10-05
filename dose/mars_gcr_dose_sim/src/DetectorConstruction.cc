#include "DetectorConstruction.hh"
#include "DetectorParameters.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <string>

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction() {}

//============================================================================
// Mars regolith material, built the same way as the Bragg's-rule mixture in
// MarsGCRPenetrationDepth.cpp -- same source data (NASA NTRS Gusev Crater
// "Panda" average oxide wt%, sums to 98.37%, renormalized to 100% here just
// like there), but built via Geant4's own G4Material compound/mixture API
// instead of a hand-rolled elemental decomposition, so Geant4 handles the
// molar-mass arithmetic for each oxide directly from NIST element data.
//============================================================================
G4Material* DetectorConstruction::BuildMarsRegolithMaterial(){
  G4NistManager* nist = G4NistManager::Instance();

  G4Element* elSi = nist->FindOrBuildElement("Si");
  G4Element* elAl = nist->FindOrBuildElement("Al");
  G4Element* elFe = nist->FindOrBuildElement("Fe");
  G4Element* elMg = nist->FindOrBuildElement("Mg");
  G4Element* elCa = nist->FindOrBuildElement("Ca");
  G4Element* elNa = nist->FindOrBuildElement("Na");
  G4Element* elK  = nist->FindOrBuildElement("K");
  G4Element* elTi = nist->FindOrBuildElement("Ti");
  G4Element* elS  = nist->FindOrBuildElement("S");
  G4Element* elO  = nist->FindOrBuildElement("O");
  G4Element* elCl = nist->FindOrBuildElement("Cl");

  // Each oxide as its own G4Material. The density given here is a nominal
  // placeholder (1.0 g/cm3) and has NO effect on the final regolith mixture's
  // physics -- G4Material::AddMaterial() below only uses each oxide's
  // elemental mass-fraction breakdown (derived from atom counts + NIST
  // atomic weights), not its own bulk density. Only the final "Mars_Regolith"
  // material's density (set explicitly below to the real regolith density)
  // matters for stopping power / dose.
  auto MakeOxide = [&](const char* name, G4Element* el, G4int nEl, G4int nO){
    G4Material* m = new G4Material(name, 1.0*g/cm3, 2);
    m->AddElement(el, nEl);
    m->AddElement(elO, nO);
    return m;
  };

  G4Material* sio2  = MakeOxide("SiO2_component",  elSi, 1, 2);
  G4Material* al2o3 = MakeOxide("Al2O3_component", elAl, 2, 3);
  G4Material* feo   = MakeOxide("FeO_component",   elFe, 1, 1);
  G4Material* fe2o3 = MakeOxide("Fe2O3_component", elFe, 2, 3);
  G4Material* mgo   = MakeOxide("MgO_component",   elMg, 1, 1);
  G4Material* cao   = MakeOxide("CaO_component",   elCa, 1, 1);
  G4Material* na2o  = MakeOxide("Na2O_component",  elNa, 2, 1);
  G4Material* k2o   = MakeOxide("K2O_component",   elK,  2, 1);
  G4Material* tio2  = MakeOxide("TiO2_component",  elTi, 1, 2);
  G4Material* so3   = MakeOxide("SO3_component",   elS,  1, 3);

  // Source wt% (NASA NTRS 20170005414, Gusev Crater "Panda" average) -- same
  // numbers as MarsGCRPenetrationDepth.cpp's MarsRegolithOxides(). Sums to
  // 98.37%; Cl (0.61%) is reported elemental, not as an oxide.
  const G4double wtSiO2  = 46.52, wtAl2O3 = 10.46, wtFeO = 12.18, wtFe2O3 = 4.20,
                 wtMgO   = 8.93,  wtCaO   = 6.27,  wtNa2O = 3.02, wtK2O   = 0.41,
                 wtTiO2  = 0.87,  wtSO3   = 4.90,  wtClElemental = 0.61;
  const G4double rawSum = wtSiO2 + wtAl2O3 + wtFeO + wtFe2O3 + wtMgO + wtCaO
                         + wtNa2O + wtK2O + wtTiO2 + wtSO3 + wtClElemental; // 98.37

  G4Material* regolith = new G4Material("Mars_Regolith", MarsGCR::kRegolithDensity, 11);
  regolith->AddMaterial(sio2,  wtSiO2  / rawSum);
  regolith->AddMaterial(al2o3, wtAl2O3 / rawSum);
  regolith->AddMaterial(feo,   wtFeO   / rawSum);
  regolith->AddMaterial(fe2o3, wtFe2O3 / rawSum);
  regolith->AddMaterial(mgo,   wtMgO   / rawSum);
  regolith->AddMaterial(cao,   wtCaO   / rawSum);
  regolith->AddMaterial(na2o,  wtNa2O  / rawSum);
  regolith->AddMaterial(k2o,   wtK2O   / rawSum);
  regolith->AddMaterial(tio2,  wtTiO2  / rawSum);
  regolith->AddMaterial(so3,   wtSO3   / rawSum);
  regolith->AddElement(elCl,   wtClElemental / rawSum);

  return regolith;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* regolith = BuildMarsRegolithMaterial();

  G4double layerThickness = MarsGCR::LayerThicknessCm() * cm;
  G4double totalThickness = MarsGCR::TotalThicknessCm() * cm;
  G4double worldHalfZ = totalThickness / 2.0 + 20.0 * cm; // buffer for the beam start point
  G4double worldHalfXY = MarsGCR::kSlabHalfWidth + 10.0 * cm;

  // World
  G4Box* worldBox = new G4Box("World", worldHalfXY, worldHalfXY, worldHalfZ);
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldBox, vacuum, "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VPhysicalVolume* worldPV = new G4PVPlacement(
      nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);

  // Regolith stack: kNLayers identical slabs stacked along +Z, starting at
  // z = -totalThickness/2 (front face) so the stack is centered in World.
  G4Box* layerBox = new G4Box("RegolithLayerSolid",
                               MarsGCR::kSlabHalfWidth, MarsGCR::kSlabHalfWidth,
                               layerThickness / 2.0);
  G4LogicalVolume* layerLV = new G4LogicalVolume(layerBox, regolith, "RegolithLayerLV");
  G4VisAttributes* regolithVis = new G4VisAttributes(G4Colour(0.72, 0.45, 0.20, 0.3));
  layerLV->SetVisAttributes(regolithVis);

  G4double zFront = -totalThickness / 2.0;
  for(G4int i = 0; i < MarsGCR::kNLayers; i++){
    G4double zCenter = zFront + (i + 0.5) * layerThickness;
    std::string pvName = std::string(MarsGCR::kLayerNamePrefix) + std::to_string(i);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCenter), layerLV, pvName,
                       worldLV, false, i /* copy number == layer index */, true);
  }

  return worldPV;
}
