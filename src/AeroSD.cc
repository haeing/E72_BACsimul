#include "AeroSD.hh"
#include "AeroHit.hh"
#include "typeinfo"

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
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4StepStatus.hh"

AeroSD::AeroSD(G4String name)
  : G4VSensitiveDetector(name)//, fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("AeroCollection"); //have to modify
}

AeroSD::~AeroSD()
{}

void AeroSD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  AeroCollection = new AeroHitsCollection(SensitiveDetectorName, collectionName[0]);

  if(HCID<0) HCID = GetCollectionID(0);
  
  HCTE->AddHitsCollection(HCID, AeroCollection);
  
}

G4bool AeroSD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
{
  const G4StepPoint* preStepPoint = astep-> GetPreStepPoint();
  G4Track* atrack = astep->GetTrack();
  G4int pid = astep->GetTrack()->GetDefinition()->GetPDGEncoding();
  G4ThreeVector pos =preStepPoint->GetPosition();
      

  
  G4double tof = astep->GetPreStepPoint()->GetGlobalTime();


  AeroHit* ahit = new AeroHit(pos,tof,pid);

  AeroCollection->insert(ahit);

  return true;  
}

void AeroSD::EndOfEvent(G4HCofThisEvent* HCTE)
{

  AeroCollection->PrintAllHits();


}

