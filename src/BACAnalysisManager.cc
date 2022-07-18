#include "BACAnalysisManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "AeroHit.hh"
#include "MPPCHit.hh"



#include "Randomize.hh"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

#include <string>
#include <sstream>

extern int gCerenkovCounter;

BACAnalysisManager::BACAnalysisManager(const G4String & histname)
  :outfile(histname), fActive_(true)
{}

BACAnalysisManager::~BACAnalysisManager()
{
  SaveFile();
}

void BACAnalysisManager::SaveFile(void) const
{
  if (fActive_)
    hfile->Write();
}

void BACAnalysisManager::Terminate(void) const
{
  if(fActive_)
    {
      hfile->Write();
      hfile->Close();
    }
}

void BACAnalysisManager::BeginOfRun(const G4Run*)
{

  

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();

  hfile = new TFile(outfile, "RECREATE");
  tree = new TTree("tree","EvtGen tree");

  tree->Branch("event",&event, "event/I");
  tree->Branch("nEvt",&nEvt, "nEvt/I");
  //tree->Branch("evtid",evtid,"evtid[nEvt]/I");
  tree->Branch("evtpid",&evtpid,"evtpid/I");
  tree->Branch("evtposx",&evtposx,"evtposx/D");
  tree->Branch("evtposy",&evtposy,"evtposy/D");
  tree->Branch("evtposz",&evtposz,"evtposz/D");
  //tree->Branch("evtpid",evtpid,"evtpid[nEvt]/I");
  tree->Branch("evtnumce",&evtnumce,"evtnumce/I");


  //Aerogel
  tree->Branch("nhAero",&nhAero,"nhAero/I");
  tree->Branch("aeropid",aeropid,"aeropid[nhAero]/I");
  tree->Branch("aeroposx",aeroposx,"aeroposx[nhAero]/D");
  tree->Branch("aeroposy",aeroposy,"aeroposy[nhAero]/D");
  tree->Branch("aeroposz",aeroposz,"aeroposz[nhAero]/D");
  tree->Branch("aerotime",aerotime,"aerotime[nhAero]/D");

  
  
  //MPPC
  tree->Branch("nhMppc",&nhMppc,"nhMppc/I");
  //tree->Branch("mppcmulti"&mppcmulti,"mppcmulti/I");
  tree->Branch("mppcpid",mppcpid,"mppcpid[nhMppc]/I");
  tree->Branch("mppcposx",mppcposx,"mppcposx[nhMppc]/D");
  tree->Branch("mppcposy",mppcposy,"mppcposy[nhMppc]/D");
  tree->Branch("mppcposz",mppcposz,"mppcposz[nhMppc]/D");
  tree->Branch("mppctime",mppctime,"mppctime[nhMppc]/D");
  tree->Branch("mppcwavelength",mppcwavelength,"mppcwavelength[nhMppc]/D");

  event = 0;
  nEvt = 0;


}

void BACAnalysisManager::EndOfRun(const G4Run*)
{
  tree->Write();
  hfile->Write();
  hfile->Close();
}

void BACAnalysisManager::BeginOfEvent(const G4Event* anEvent)
{

}

void BACAnalysisManager::EndOfEvent(const G4Event* anEvent)
{

  G4HCofThisEvent* HCTE = anEvent->GetHCofThisEvent();
  if(!HCTE) return;
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhmppc = 0;
  G4int nhaero = 0;
  G4int multiplicity = 0;

  G4int pdg = anEvent->GetPrimaryVertex(0)->GetPrimary(0)->GetPDGcode();
  evtpid = pdg;
  
  G4ThreeVector pripos = anEvent->GetPrimaryVertex(0)->GetPosition();
  G4double prix = pripos.x();
  G4double priy = pripos.y();
  G4double priz = pripos.z();
  evtposx = prix;
  evtposy = priy;
  evtposz = priz;

  evtnumce = gCerenkovCounter;
  

  
  MPPCHitsCollection *MPPCHC = 0;
  G4int ColIdMPPC = SDMan->GetCollectionID("MppcCollection");
  if(ColIdMPPC>=0)
    {
      MPPCHC=dynamic_cast<MPPCHitsCollection *>(HCTE->GetHC( ColIdMPPC ));
      if(MPPCHC)
	{
	  nhmppc = MPPCHC->entries();
	}
    }

  for(int i=0;i<nhmppc;i++)
    {
      MPPCHit* aHit = (*MPPCHC)[i];
      mppctime[i] = aHit->GetTOF();
      mppcposx[i] = aHit->GetPosition().x();
      mppcposy[i] = aHit->GetPosition().y();
      if(aHit->GetWorldPosition().x()<=0)mppcposz[i]=-1.0;
      else if(aHit->GetWorldPosition().x()>0)mppcposz[i]=1.0;
      //mppcposz[i] = aHit->GetPosition().z();
      //mppcposz[i] = aHit->GetWorldPosition().y();
      mppcpid[i] = aHit->GetParticleID();
      mppcwavelength[i] = aHit->GetWavelength();
      
    }
  nhMppc = nhmppc;

  AeroHitsCollection *AEROHC = 0;
  G4int ColIdAERO = SDMan->GetCollectionID("AeroCollection");
  if(ColIdAERO>=0)
    {
      AEROHC=dynamic_cast<AeroHitsCollection *>(HCTE->GetHC(ColIdAERO));
      if(AEROHC)
	{
	  nhaero = AEROHC->entries();
	}
    }

  for(int i=0;i<nhaero;i++)
    {
      AeroHit* aHit = (*AEROHC)[i];
      aerotime[i] = aHit->GetTOF();
      aeroposx[i] = aHit->GetPosition().x();
      aeroposy[i] = aHit->GetPosition().y();
      aeroposz[i] = aHit->GetPosition().z();
      aeropid[i] = aHit->GetParticleID();

    }
  nhAero = nhaero;

  tree->Fill();
  event++;

  nEvt=0;
  nhaero = 0;
  nhmppc = 0;


  
}

void BACAnalysisManager::SetEvtGen(int j,int partnum)
{
  if(nEvt<j) nEvt=j;
  //evtid[j-1]=j;
  //evtpid[j-1]=partnum;
}


  
      



