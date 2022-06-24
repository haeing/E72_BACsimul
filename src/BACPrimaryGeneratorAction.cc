#include "BACDetectorConstruction.hh"
#include "BACPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


#include "TFile.h"
#include "TTree.h"
#include "Randomize.hh"

BACPrimaryGeneratorAction::BACPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  
  fParticleGun = new G4ParticleGun(n_particle);

}

BACPrimaryGeneratorAction::~BACPrimaryGeneratorAction()
{
  
  delete fParticleGun;
}

void BACPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4int event_num = anEvent->GetEventID()+1;
  if(event_num% 1000 == 0)
    {
      G4cout<<"Event# "<<event_num<<G4endl;
    }



  G4ThreeVector D(0*mm,0*mm,-100*mm);
  G4ThreeVector P(0,0,1);

  //GenerateBeamKaonMBr(anEvent,D,P,particle);


 

  
  //homogeneous test------------------------------------
  
  G4double momentum = 0.735;
  G4double beam_size = 125; 

  if(particle=="kaon"){
    energy = (sqrt(mass_kaonm*mass_kaonm+momentum*momentum) - mass_kaonm )*GeV;
    fParticleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  }
  if(particle=="pion"){
    energy = (sqrt(mass_pim*mass_pim+momentum*momentum) - mass_pim )*GeV;
    fParticleGun->SetParticleDefinition (particleTable -> FindParticle("pi-"));
  }
  
  fParticleGun->SetParticleMomentumDirection ( G4ThreeVector(0,0,1) );

  G4double x = beam_size*0.5-G4UniformRand()*beam_size*mm;
  //original 220621
  //G4double y = beam_size*0.5-G4UniformRand()*beam_size*mm;
  G4double y = beam_size*0.5-G4UniformRand()*beam_size*mm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x,y,-10*cm) );
  //fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-10*cm) );
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);


  


}

void BACPrimaryGeneratorAction::GenerateBeamKaonMBr(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P,G4String particle)
{
  G4ThreeVector X_file, P_file;
  ReadBeamProfile(X_file, P_file);

  G4ThreeVector beamx(X_file.x()*mm + D.x(), X_file.y()*mm + D.y(), D.z());

  double beam_file_mag = P_file.mag();
  double beam_mag_set = (P.mag()/0.9)*beam_file_mag; //for 0.9 GeV/c beam




  G4ThreeVector p_dir(P_file.x()/beam_file_mag, P_file.y()/beam_file_mag, P_file.z()/beam_file_mag) ; //from beam profile file


  G4ThreeVector beamp (beam_mag_set*p_dir.x(), beam_mag_set*p_dir.y(), beam_mag_set*p_dir.z());

  //beam rotate angle
  G4double rotate_angle = 0.0*degree; //1.8 GeV/c && field = 0.0 case

  G4ThreeVector beampu =  beamp/beamp.mag();

  

  if(particle=="kaon"){
    energy = (sqrt(mass_kaonm*mass_kaonm+beamp.mag2()) - mass_kaonm )*GeV;
    fParticleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  }

  if(particle=="pion"){
    energy = (sqrt(mass_pim*mass_pim+beamp.mag2()) - mass_pim )*GeV;
    fParticleGun -> SetParticleDefinition (particleTable -> FindParticle("pi-"));
  }
  
  fParticleGun->SetParticleMomentumDirection ( beampu );

  //fParticleGun->SetParticleTime ( 0.0 );
  fParticleGun->SetParticlePosition( beamx );

  fParticleGun->SetParticleEnergy( energy );

  fParticleGun->GeneratePrimaryVertex( anEvent);

  

}

void BACPrimaryGeneratorAction::ReadBeamProfile( G4ThreeVector & X, G4ThreeVector & P )
{

  G4String beamfilename = "../param/beam/beam.k.run69_0130.root";

  TFile *beam_file = new TFile(beamfilename, "read");

  
  TTree *beam_tree = (TTree*)beam_file->Get("tr");

  int ntK18;
  double pointInx[1];
  double pointIny[1];
  double pointInz[1];
  double pInx[1];
  double pIny[1];
  double pInz[1];


  beam_tree->SetBranchAddress("ntK18",&ntK18);
  beam_tree->SetBranchAddress("pointInx",  &pointInx);
  beam_tree->SetBranchAddress("pointIny",  &pointIny);
  beam_tree->SetBranchAddress("pointInz",  &pointInz);


  beam_tree->SetBranchAddress("pInx",  &pInx);
  beam_tree->SetBranchAddress("pIny",  &pIny);
  beam_tree->SetBranchAddress("pInz",  &pInz);

  bp_file_ndata = beam_tree->GetEntries();


  if(bp_file_ndata == bp_nAccess) bp_nAccess = 0;
  
  beam_tree->GetEntry(bp_nAccess);
 
  G4ThreeVector TVp(pInx[0], pIny[0], pInz[0]);
  G4ThreeVector TVx(pointInx[0], pointIny[0], pointInz[0]);


  X=TVx;
  P=TVp;
  bp_nAccess++;

  beam_file->Close();

}



    
  
