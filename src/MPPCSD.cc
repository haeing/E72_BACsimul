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
  : G4VSensitiveDetector(name)//, fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("MppcCollection");
}

MPPCSD::~MPPCSD()
{}

void MPPCSD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  MppcCollection = new MPPCHitsCollection(GetName(), collectionName[0]);
  //if (fHCID<0) fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  if (HCID<0) HCID = GetCollectionID(0);

  HCTE->AddHitsCollection(HCID, MppcCollection);
  /*
  auto ahit = new MPPCHit();  //fill hits with zero energy deposition
  //ahit->SetEdep(0);
  ahit->SetTime(0);
  ahit->ClearPhotonCount();
  fHitsCollection ->insert(ahit);
  aa = 0.0;
  std::cout<<"MPPCSD start"<<std::endl;
  */


  

}

G4bool MPPCSD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
{
  
  const G4StepPoint* preStepPoint = astep-> GetPreStepPoint();
  G4Track* atrack = astep->GetTrack();
  //G4double time = track->GetLocalTime();
  //std::cout<<"MPPCSD localtime "<<time<<std::endl;
  //G4String particleName = aStep->GetTrack()->GetDefinition->GetParticleName();
  G4int pid = astep->GetTrack()->GetDefinition()-> GetPDGEncoding();
  //auto edep = step->GetTotalEnergyDeposit();
  G4ThreeVector worldPos =preStepPoint->GetPosition();
  G4ThreeVector pos = preStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);
  //G4double energy = preStepPoint->GetTotalEnergy();
  //G4double energy = preStepPoint->GetKineticEnergy();
  G4double energy = atrack->GetTotalEnergy();

  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;
  G4double wavelength = ((h*c)/(energy*pow(10,6)*1.6e-13))/(1e-9); //nm
  //G4double wavelength = energy;
  //to localPos = step->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  //auto edep = track->GetKineticEnergy();
  //G4ThreeVector hitmom = track->GetMomentum()/CLHEP::eV;
  //G4double E_p = hitmom.mag();
  
  auto hitTime = astep->GetPreStepPoint()->GetGlobalTime();
  //if (edep==0.) return true;
  //if (hitTime==0.) return true;

  MPPCHit* ahit = new MPPCHit(pos,wavelength,pid,hitTime);
  MppcCollection->insert(ahit);
  return true;
  /*
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto physical = touchable->GetVolume();
  auto copyNo = physical ->GetCopyNo();   //volume ID
  auto hit = (MPPCHit*) fHitsCollection->GetHit(0);
  hit->SetPos(localPos);

  

  //hit-> AddEdep(edep);
  hit-> IncPhotonCount();
  
  if(hitTime>aa){
    //hit->SetTime(hitTime);
    aa=hitTime;
    hit->SetTime(aa);
  
  }
  */

  
  //hit->SetTrackLen(tracklen);
  //hit->SetStepNum(numstep);

  //if(particleName == "opticalphoton"){
  //const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
  //if(secondaries->size()>0){
     
      //}
      //}
      //fHitsCollection->insert(hit);

  
  //auto ix=-1;
  //check if this finger already has a hit
  /*
  for(std::size_t i=0;i<fHitsCollection->entries();++i){
    if((*fHitsCollection)[i]->GetID() == copyNo){
      ix = i;
      break;
    }
  }

  if(ix>=0){
    //if it has, then take the earlier time
    if ((*fHitsCollection)[ix]->GetTime()>hitTime){
      (*fHitsCollection)[ix]->SetTime(hitTime);
    }
  }
  else{
    hit->SetTime(hitTime);
  }
  */
  /*
  else{
    //if not, create a new hit and set it to the collection
    auto hit1 = new MPPCHit(copyNo, hitTime);
    hit1->SetTime(hitTime);
  }
  */
  

    
  

  
  
      

}


void MPPCSD::EndOfEvent(G4HCofThisEvent *HCTE)
{
  MppcCollection->PrintAllHits();

}


//void MPPCSD::PrintAll() const
  
