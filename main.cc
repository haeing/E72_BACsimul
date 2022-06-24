#include "BACDetectorConstruction.hh"
#include "BACPrimaryGeneratorAction.hh"
#include "BACRunAction.hh"
#include "BACEventAction.hh"
#include "BACStackingAction.hh"
#include "BACAnalysisManager.hh"
#include "AeroSD.hh"
#include "AeroHit.hh"
#include "MPPCSD.hh"
#include "MPPCHit.hh"

//#ifdef G4Multithreded
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"
//#include "G4RunManagerFactory.hh"

#include "G4VisExecutive.hh"
#include "G4ScoringManager.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  CLHEP::HepRandom::setTheSeed((unsigned)time(NULL));

  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  G4String histname;
  /*
  if(argc>=3){
    histname = argv[2];
  }
  else
    {
      histname = "geant4_test.root";
    }
  */


  //new test for kekcc sicmulation with input data(detector size)
  G4String num_aerogel;
  G4String reflect_part_length;
  G4String light_guide_length;
  if(argc>=3){
    histname = argv[2];
    num_aerogel = argv[3];
    reflect_part_length = argv[4];
    light_guide_length = argv[5];
  }
  else
    {
      histname = "geant4_test.root";
      num_aerogel = "3";
      reflect_part_length = "110.0";
      light_guide_length= "110.0";
    }
  

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;
  
  
  //    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new BACDetectorConstruction(num_aerogel,reflect_part_length,light_guide_length));
  

  // Physics list
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  //opticalPhysics->SetScintillationYieldFactor(1.);
  //opticalPhysics->SetScintillationExcitationRatio(0.);
  //opticalPhysics->SetWLSTimeProfile("delta");
  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  //opticalPhysics->Configure(kCerenkov,false);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  //opticalPhysics->Configure(kScintillation,false);
  physicsList->RegisterPhysics(opticalPhysics);

  runManager->SetUserInitialization(physicsList);

  BACAnalysisManager *anaMan = new BACAnalysisManager(histname);
  BACPrimaryGeneratorAction *priGen = new BACPrimaryGeneratorAction();
  BACRunAction *runAction = new BACRunAction(anaMan);
  BACEventAction *eventAction = new BACEventAction(anaMan);
  BACStackingAction *stackAction = new BACStackingAction();

  runManager->SetUserAction(priGen);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(stackAction);
    
  // User action initialization
  //runManager->SetUserInitialization(new BACActionInitialization());
  runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);

  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;

  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !


  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
