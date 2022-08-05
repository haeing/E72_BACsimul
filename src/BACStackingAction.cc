#include "BACStackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ClassificationOfNewTrack.hh"

int gCerenkovCounter;


BACStackingAction::BACStackingAction()
  : G4UserStackingAction(),
    fScintillationCounter(0), fCerenkovCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BACStackingAction::~BACStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
BACStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{

  G4ClassificationOfNewTrack classification = fWaiting;
  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;
  
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
        fScintillationCounter++;
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
	if(((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)>320&&((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)<900){
	  fCerenkovCounter++;
	  gCerenkovCounter++;
	}
	else return fKill;
      }
    }
  }
  /*
    G4Track* tr = (G4Track*) aTrack;
    G4String volume;
  
  
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  
  if(particle==G4PionMinus::PionMinus()){
    volume = aTrack->GetNextVolume()->GetName();
    //if(volume == "Part2LW")
    //tr->SetTrackStatus(fStopButAlive);
    std::cout<<" "<<particle<<std::endl;
  }
  */

  
  //auto *Change = new BACChange();
  //if(aTrack == nullptr)std::cout<<"error"<<std::endl;
  //Change->Analysis(*aTrack);
  //Change->ShowInfo();
		   
  return fUrgent;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BACStackingAction::NewStage()
{
  //std::cout << "Number of Scintillation photons produced in this event : "
  //	    << fScintillationCounter << std::endl;
  std::cout << "Number of Cerenkov photons produced in this event : "
	    << fCerenkovCounter << std::endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BACStackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
  gCerenkovCounter = 0;

}
