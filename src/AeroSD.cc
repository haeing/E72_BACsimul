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
  : G4VSensitiveDetector(name), fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("aeroColl"); //have to modify
}

AeroSD::~AeroSD()
{}

void AeroSD::Initialize(G4HCofThisEvent *hce)
{

  fHitsCollection = new AeroHitsCollection(GetName(), collectionName[0]);
  if (fHCID<0) fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  hce->AddHitsCollection(fHCID, fHitsCollection);
  auto ahit = new AeroHit();  //fill hits with zero energy deposition
  ahit->SetEdep(0);
  ahit->SetTime(0);
  fHitsCollection ->insert(ahit);
  if (fHCID<0){
    id1 = 0;
    id3 = 0;
  }


  

}

G4bool AeroSD::ProcessHits(G4Step *step, G4TouchableHistory *)
{

  G4Track* track = step->GetTrack();
  G4String particleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  auto edep = step->GetTotalEnergyDeposit();

  /*
  G4int id = track->GetTrackID();
  //if(id == id1){
    id3 = id;
  std::cout<<"trackID"<<id<<std::endl;

  //G4double time = track->GetGlobalTime();
  G4double time = step->GetDeltaTime();
  std::cout<<"deltatime"<<time<<std::endl;
  if(time<0)track->SetLocalTime(0);


  G4double pid = track->GetParentID();
  std::cout<<"parent id"<<pid<<std::endl;

  G4double pretime = step->GetPreStepPoint()->GetLocalTime();
  G4double posttime = step->GetPostStepPoint()->GetLocalTime();
  std::cout<<"pre "<<pretime<<std::endl;
  std::cout<<"post "<<posttime<<std::endl;
  if(posttime<pretime)step->GetPostStepPoint()->AddLocalTime((pretime-posttime)*2);

  const G4ThreeVector& vector = track->GetPosition();
  std::cout<<"x "<<vector[0]<<std::endl;
  std::cout<<"y "<<vector[1]<<std::endl;
  std::cout<<"z "<<vector[2]<<std::endl;

  G4double pretime1 = step->GetPreStepPoint()->GetGlobalTime();
  G4double posttime1 = step->GetPostStepPoint()->GetGlobalTime();
  std::cout<<"pre global "<<pretime1<<std::endl;
  std::cout<<"post global"<<posttime1<<std::endl;

  G4int vol = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
  std::cout<<"volume"<<vol<<std::endl;


  if(track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
    if(track->GetParentID()>0){
      G4ProcessType type =  track->GetCreatorProcess()->GetProcessType();
      std::cout<<G4VProcess::GetProcessTypeName(type)<<std::endl;
    }
  }

  std::cout<<"-----------------------------------------------"<<std::endl;
  */
  //id1+=1;
  //}
  //if(id != id3)id1 = track->GetTrackID();

  //auto edep = track->GetKineticEnergy();
  //G4ThreeVector hitmom = track->GetMomentum()/CLHEP::eV;
  //G4double E_p = hitmom.mag();
  
  //if (edep==0.) return true;
  //if (hitTime==0.) return true;

  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto physical = touchable->GetVolume();
  auto copyNo = physical ->GetCopyNo();   //volume ID
  auto hit = (AeroHit*) fHitsCollection->GetHit(0);


  auto hitTime = step->GetPreStepPoint()->GetGlobalTime();
  hit-> AddEdep(edep);
  hit->SetTime(hitTime);
  //hit->SetTrackLen(tracklen);
  //hit->SetStepNum(numstep);

  //if(particleName == "opticalphoton"){
  //const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
  //if(secondaries->size()>0){
  //hit-> IncPhotonCount();
      //}
      //}

  
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
    auto hit1 = new AeroHit(copyNo, hitTime);
    hit1->SetTime(hitTime);
  }
  */
  

    
  

  
  

  return true;
}
