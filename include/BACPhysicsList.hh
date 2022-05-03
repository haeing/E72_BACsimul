#ifndef BACPhysicsList_h
#define BACPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//class BACPhysicsListMessenger;

class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpBoundaryProcess;
class G4OpRayleigh;
class G4OpMieHG;


class BACPhysicsList : public G4VUserPhysicsList
{
public:
  BACPhysicsList();
  virtual ~BACPhysicsList();

public:

  virtual void ConstructParticle();
  virtual void ConstructProcess();

  virtual void SetCuts();
  
  void ConstructDecay();
  void ConstructEM();
  void ConstructOp();

  void SetVerbose(G4int);
  void SetNbOfPhotonsCerenkov(G4int);

private:

  // BACPhysicsListMessenger* fMessenger;

  static G4ThreadLocal G4int fVerboseLevel;
  static G4ThreadLocal G4int fMaxNumPhotonStep;

  static G4ThreadLocal G4Cerenkov* fCerenkovProcess;
  static G4ThreadLocal G4Scintillation* fScintillationProcess;
  static G4ThreadLocal G4OpAbsorption* fAbsorptionProcess;
  static G4ThreadLocal G4OpBoundaryProcess* fBoundaryProcess;
  static G4ThreadLocal G4OpRayleigh* fRayleighScatteringProcess;
  static G4ThreadLocal G4OpMieHG* fMieHGScatteringProcess;
};

#endif


  
