#include "BACActionInitialization.hh"
#include "BACPrimaryGeneratorAction.hh"
#include "BACRunAction.hh"
#include "BACEventAction.hh"
#include "BACAnalysisManager.hh"
#include "BACStackingAction.hh"


BACActionInitialization::BACActionInitialization()
  : G4VUserActionInitialization()
{}

BACActionInitialization::~BACActionInitialization()
{}

void BACActionInitialization::BuildForMaster() const
{
  //SetUserAction(new BACRunAction());
}

void BACActionInitialization::Build() const
{
  SetUserAction(new BACPrimaryGeneratorAction());
  SetUserAction(new BACRunAction());
  //SetUserAction(new BACAnalysisManager());
  SetUserAction(new BACEventAction());
  SetUserAction(new BACStackingAction());

}
