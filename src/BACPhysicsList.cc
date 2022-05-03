#include "globals.hh"
#include "BACPhysicsList.hh"
//#include "BACPhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleWithCuts.hh"

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ProcessManager.hh"
//#include "G4OpticalPhysics.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4Threading.hh"

using namespace CLHEP;

G4ThreadLocal G4int BACPhysicsList::fVerboseLevel = 0;
G4ThreadLocal G4int BACPhysicsList::fMaxNumPhotonStep = 0;
G4ThreadLocal G4Cerenkov* BACPhysicsList::fCerenkovProcess = 0;
G4ThreadLocal G4Scintillation* BACPhysicsList::fScintillationProcess = 0;
G4ThreadLocal G4OpAbsorption* BACPhysicsList::fAbsorptionProcess = 0;
G4ThreadLocal G4OpBoundaryProcess* BACPhysicsList::fBoundaryProcess = 0;
G4ThreadLocal G4OpRayleigh* BACPhysicsList::fRayleighScatteringProcess = 0;
G4ThreadLocal G4OpMieHG* BACPhysicsList::fMieHGScatteringProcess = 0;

BACPhysicsList::BACPhysicsList()
  : G4VUserPhysicsList()
{}

BACPhysicsList::~BACPhysicsList() {}

void BACPhysicsList::ConstructParticle()
{
  G4BosonConstructor bConstructor;
  bConstructor.ConstructParticle();

  G4LeptonConstructor lConstructor;
  lConstructor.ConstructParticle();

  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  G4BaryonConstructor rConstructor;
  rConstructor.ConstructParticle();

  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void BACPhysicsList::ConstructProcess()
{
  
  G4StepLimiter* StepLimit = new G4StepLimiter();
  G4UserSpecialCuts* UserCuts = new G4UserSpecialCuts();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while((*particleIterator)() ) {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager->AddDiscreteProcess(StepLimit);
    pmanager->AddDiscreteProcess(UserCuts);
  }
  
  AddTransportation();
  ConstructDecay();
  ConstructEM();
  ConstructOp();
}
/*
#include "BACTransportation.hh"

void BACPhysicsList::AddTransportation() {

  BACTransportation* theTransportationProcess = new BACTransportation();
  
  G4VUserPhysicsList::AddTransportation();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)()){
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (!particle->IsShortLived() ) {
	if (pmanager == 0) {
	  G4cout << "Warning" << G4endl;
	}
	else {
     	pmanager->AddProcess(theTransportationProcess);
	pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
	pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
	}
      }	
      else {
      }
  }
}
*/

#include "G4Decay.hh"

void BACPhysicsList::ConstructDecay()
{
  G4Decay* theDecayProcess = new G4Decay();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ( (*particleIterator)())
    {
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      if (theDecayProcess->IsApplicable(*particle))
	{
	  pmanager->AddProcess(theDecayProcess);
	  pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
	  pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
	}
    }
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4LivermoreIonisationModel.hh"

#include "G4UAtomicDeexcitation.hh"

void BACPhysicsList::ConstructEM()
{

  G4EmParameters* param = G4EmParameters::Instance();
  param->SetMaxEnergy(100*GeV);
  param->SetNumberOfBinsPerDecade(20);
  param->SetMscStepLimitType(fMinimal);
  param->SetFluo(true);
  param->SetPixe(true);
  param->SetAuger(true);
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  if (!ad) {
    man->SetAtomDeexcitation(new G4UAtomicDeexcitation());
  }

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while((*particleIterator)() )
    {
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "gamma")
	{
	  pmanager->AddDiscreteProcess(new G4GammaConversion());
	  pmanager->AddDiscreteProcess(new G4ComptonScattering());
	  pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
	}
      else if (particleName == "e-")
	{
	  pmanager = G4Electron::Electron()->GetProcessManager();
	  
	  pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
	  pmanager->AddProcess(new G4eBremsstrahlung(),-1, -1, 3);
	  G4eIonisation* eIonisation = new G4eIonisation();
	  //	  eIonisation->SetEmModel(new G4LivermoreIonisationModel());
	  eIonisation->SetStepFunction(0.2, 100*um);
	  pmanager->AddProcess(eIonisation, -1, 2, 2);
	  
	}
      else if (particleName == "e+")
	{
	   pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
	  pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
	  pmanager->AddProcess(new G4eBremsstrahlung(), -1, 3, 3);
	  pmanager->AddProcess(new G4eplusAnnihilation(), 0 ,-1, 4);
	}
      else if (particleName == "mu-" || particleName == "mu+" || particleName == "mu=")
	{
	  pmanager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
	  pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
	  pmanager->AddProcess(new G4MuBremsstrahlung(), -1, 3, 3);
	  pmanager->AddProcess(new G4MuPairProduction(), -1, 4, 4);
	}
      else {
	if ( (particle->GetPDGCharge() != 0.0) && (particle->GetParticleName() != "chargedgeantino") && !particle->IsShortLived())
	  {
	    pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
	    pmanager->AddProcess(new G4hIonisation(), -1, 2, 2);
	  }
      }
    }
}


void BACPhysicsList::ConstructOp()
{
  fCerenkovProcess = new G4Cerenkov("Cerenkov");
  fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
  fCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  fCerenkovProcess->SetTrackSecondariesFirst(true);

  fScintillationProcess = new G4Scintillation("Scintillation");
  fScintillationProcess->SetScintillationYieldFactor(1.0);
  fScintillationProcess->SetTrackSecondariesFirst(true);
  
  fAbsorptionProcess = new G4OpAbsorption();
  fBoundaryProcess = new G4OpBoundaryProcess();
  fRayleighScatteringProcess = new G4OpRayleigh();
  fMieHGScatteringProcess = new G4OpMieHG();
  
  fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
  fScintillationProcess->SetVerboseLevel(fVerboseLevel);
  fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
  fBoundaryProcess->SetVerboseLevel(fVerboseLevel);
  fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
  fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);

  if(G4Threading::IsMasterThread())
    {
      G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
      fScintillationProcess->AddSaturation(emSaturation);
    }
  
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
    {
      G4ParticleDefinition *particle = particleIterator->value();
      G4ProcessManager * pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      
      if (fCerenkovProcess->IsApplicable(*particle))
	{
	  pmanager->AddProcess(fCerenkovProcess);
	  pmanager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
	}
      
      if (fScintillationProcess->IsApplicable(*particle))
	{
	  pmanager->AddProcess(fScintillationProcess);
	  pmanager->SetProcessOrderingToLast(fScintillationProcess,idxAtRest);
	  pmanager->SetProcessOrderingToLast(fScintillationProcess,idxPostStep);
	}
      
      if (particleName == "opticalphoton")
	{
	  G4cout << "Add Discrete Process for Optical Photon" << G4endl;
	  pmanager->AddDiscreteProcess(fAbsorptionProcess);
	  pmanager->AddDiscreteProcess(fBoundaryProcess);
	  pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
	  pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
	}
    }
}

void BACPhysicsList::SetVerbose(G4int verbose)
{
  fVerboseLevel = verbose;
  fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
  fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
  fBoundaryProcess->SetVerboseLevel(fVerboseLevel);
  fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
  fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
  
}  
void BACPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
  fMaxNumPhotonStep = MaxNumber;
  fCerenkovProcess -> SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
}

void BACPhysicsList::SetCuts()
{
  SetCutsWithDefault();
  G4int temp = GetVerboseLevel();

  SetVerboseLevel(0);

  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  G4double lowerbound = 250*CLHEP::eV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowerbound,10*CLHEP::GeV);

  auto particleIterator=GetParticleIterator(); //added for ver.10.03                                                                                  
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition *particle=particleIterator->value();
    particle->SetApplyCutsFlag( true );
    //////////////////////////////////////////////////////////////////////                                                                               
    //    G4cout << particle->GetParticleName() << " ==> ApplyCutFlag = "                                                                                
    //	   << particle->GetApplyCutsFlag() << G4endl;                                                                                                    
    //////////////////////////////////////////////////////////////////////                                                                               
  }

  
  //  if (verboseLevel > 0) DumpCutValuesTable();
}
	  
				       
	  
			    
