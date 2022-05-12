#include "BACRunAction.hh"
#include "BACPrimaryGeneratorAction.hh"
#include "BACDetectorConstruction.hh"
#include "BACActionInitialization.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <stdlib.h>

BACRunAction::BACRunAction()
  :G4UserRunAction(), RfileName("test.root")
{
  
}
	     

BACRunAction::~BACRunAction()
{
}

void BACRunAction::BeginOfRunAction(const G4Run*)
{

  file = new TFile(RfileName, "recreate");
  TTree *tree = new TTree("tree","simulation");
  Tree = (TTree*)file->Get("tree");

  
}

void BACRunAction::EndOfRunAction(const G4Run*)
{

  Tree->Write();
  file->Close();
}


