#include "BACEventAction.hh"
#incldue "BACAnlaysisManager.hh"

#include "G4Event.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "g4analysis.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4Step.hh"
#include "G4PrimaryParticle.hh"
#include <TMath.h>

#include "BACPrimaryGeneratorAction.hh"



namespace{
G4VHitsCollection* GetHC(const G4Event* event, G4int collId){
  auto hce = event->GetHCofThisEvent();
  if(!hce){
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl;
    G4Exception("BACEventAction::EndOfEventAction()","BACCode001",JustWarning,msg);
    return nullptr;
    }

  auto hc = hce->GetHC(collId);
  if(! hc){
    G4ExceptionDescription msg;
    msg << "Hits collection" << collId << " of this event not found" <<G4endl;
    G4Exception("BACEventAction::EndOfEventAction()","BACCode001",JustWarning, msg);
    }
  return hc;
}
}


BACEventAction::BACEventAction()
  : G4UserEventAction(),fAeroHCID(-1),fMPPCHCID(-1)
{
}

BACEventAction::~BACEventAction()
{}

void BACEventAction::BeginOfEventAction(const G4Event* evt)
{

  G4int event_id = evt -> GetEventID();

  if (event_id == 0)
    DefineTree();

  if (fAeroHCID == -1&&fMPPCHCID==-1){
    auto sdManager = G4SDManager::GetSDMpointer();


    
    G4String mppcHCName = "mppcSD/mppcColl";


    fMPPCHCID = sdManager->GetCollectionID(mppcHCName);

    G4String aeroHCName = "aeroSD/aeroColl";
    fAeroHCID = sdManager->GetCollectionID(aeroHCName);




  }
  //phit=0.0;
 
}


  
void BACEventAction::EndOfEventAction(const G4Event* evt/*, const G4Step* step*/)
{

  G4HCofThisEvent *HCE = evt -> GetHCofThisEvent();
  //if(!HCE)return;
  

  

  auto hc = GetHC(evt,fAeroHCID);
  auto hc1 = GetHC(evt,fMPPCHCID);


  //if (!hc) return;
  //auto nhit = hc->GetSize();
  //auto analysisManager = G4AnalysisManager::Instance();
  
  G4double phit = 0.0;

  G4double t;

  
  //for(unsigned int i=0;i<hc->GetSize();++i){
  

  for (int i=0;i<HCE->GetCapacity();i++){
    auto hit = (MPPCHit*)(HCE->GetHC(i)->GetHit(0));
    phit += hit->GetPhotonCount();
    std::cout<<"photon"<<phit<<std::endl;
    t = hit->GetTime();

  }

  /*
  for(int i=0;i<hc1->GetSize();++i){
    auto hit = static_cast<MPPCHit*>(hc1->GetHit(i));
    phit = phit+hit->GetPhotonCount();
  }
  */

  //Change->ShowInfoEvent();
  

  
  G4int event_id = evt->GetEventID();


  G4ThreeVector position = evt->GetPrimaryVertex(0)->GetPosition();
  G4int pdg = evt->GetPrimaryVertex(0)->GetPrimary(0)->GetPDGcode();
  xPrm = position.x();
  yPrm = position.y();
  zPrm = position.z();


  auto primary = evt->GetPrimaryVertex(0)->GetPrimary(0);
  double energy = primary->GetTotalEnergy();
  //double energy = primary->GetKineticEnergy();
  //double momentum = primary-> GetTotal

  numPho = phit;


  eventID = event_id;
  particleID = pdg;
  init_energy = energy;
  time = t;
  tree->Fill();

  

}

void BACEventAction::DefineTree(){
  //tree = dynamic_cast<TTree*>(gFile->Get("tree"));
  tree = (TTree*)(gFile->Get("tree"));
  tree->Branch("eventID",&eventID,"eventID/I");
  tree->Branch("particleID",&particleID,"particleID/I");
  tree->Branch("init_energy",&init_energy,"init_energy/D");
  tree->Branch("numPho",&numPho,"numPho/D");
  //tree->Branch("edep",&edep,"edep/D");
  //tree->Branch("time",&time,"time/D");
  tree->Branch("xPrm",&xPrm,"xPrm/D");
  tree->Branch("yPrm",&yPrm,"yPrm/D");
  tree->Branch("zPrm",&zPrm,"zPrm/D");

	       
}


