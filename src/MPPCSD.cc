#include "MPPCSD.hh"
#include "MPPCHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

MPPCSD::MPPCSD(G4String name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("MppcCollection");
}

MPPCSD::~MPPCSD()
{}

void MPPCSD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  MppcCollection = new MPPCHitsCollection(GetName(), collectionName[0]);

  if (HCID<0) HCID = GetCollectionID(0);

  HCTE->AddHitsCollection(HCID, MppcCollection);


  

}

G4bool MPPCSD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
{
  
  const G4StepPoint* preStepPoint = astep-> GetPreStepPoint();
  G4Track* atrack = astep->GetTrack();
  G4int pid = astep->GetTrack()->GetDefinition()-> GetPDGEncoding();
  G4ThreeVector worldPos =preStepPoint->GetPosition();
  G4ThreeVector pos = preStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  G4double energy = atrack->GetTotalEnergy();

  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;
  G4double wavelength = ((h*c)/(energy*1.6e-13))*(1e+9); //nm

  

  
  auto hitTime = astep->GetPreStepPoint()->GetGlobalTime();

  atrack->SetTrackStatus(fStopAndKill);

  
  MPPCHit* ahit = new MPPCHit(pos,hitTime,pid,wavelength);
  MppcCollection->insert(ahit);
  return true;
      

}


void MPPCSD::EndOfEvent(G4HCofThisEvent *HCTE)
{
  MppcCollection->PrintAllHits();

}
  
