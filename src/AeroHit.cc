#include "AeroHit.hh"
#include "BACDetectorConstruction.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AttDefStore.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<AeroHit>* AeroHitAllocator;

AeroHit::AeroHit()
  : G4VHit(),fEdep(0.)
{}

AeroHit::AeroHit(G4int id, G4double time)
  : G4VHit(),fEdep(0.),fTime(0.0), fId(id)//,fPhotons(0.0)
{}

AeroHit::~AeroHit()
{}

AeroHit::AeroHit(const AeroHit &right)
  : G4VHit(),fId(right.fId),fEdep(right.fEdep),fTime(right.fTime)
{}

const AeroHit& AeroHit::operator=(const AeroHit &right)
{
  fEdep = right.fEdep;
  fTime = right.fTime;
  fId = right.fId;
  return *this;
}



const std::map<G4String,G4AttDef>* AeroHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("AeroHit",isNew);

  if (isNew){
    //(*store)["Energy"] = G4AttDef("Energy","Energy Deposited","Physics","G4BestUnit","G4double");
    (*store)["Energy"] = G4AttDef("Energy","Energy Deposited","Physics","MeV","G4double");
    (*store)["ID"] = G4AttDef("ID","ID","Physics","","G4int");
    (*store)["Time"] = G4AttDef("Time","Time","Physics","ns","G4double");
  }
  return store;
}

std::vector<G4AttValue>* AeroHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  values ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
  values ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
  values ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fId),""));
  return values;
}




