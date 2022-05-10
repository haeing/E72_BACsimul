#ifndef MPPCSD_h
#define MPPCSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MPPCHit.hh"
//#include "TGraph.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class MPPCSD : public G4VSensitiveDetector
{
public:
  MPPCSD(G4String name);
  //MPPCSD(G4String name, ParamMan*);
  virtual ~MPPCSD();

  virtual void Initialize( G4HCofThisEvent *HCE );
  virtual G4bool ProcessHits(G4Step* astep, G4TouchableHistory *ROhist);
  virtual void EndOfEvent(G4HCofThisEvent *);
  
  //void EndOfEvent( G4HCofThisEvent *HCE );
  //G4bool ProcessHits_constStep(const G4Step*, G4TouchableHistory* );
  //G4int QEFlag;
  //TGraph* QETable;
  //void SetQETable(G4int);

  //void PrintAll() const;

private:
  //int EMFlag;
  MPPCHitsCollection *fHitsCollection;
  G4int fHCID;
  G4double aa;
  //G4double fEdep;
  //G4double fTime;
};

#endif
