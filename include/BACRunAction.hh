#ifndef BACRunAction_hh
#define BACRunAction_hh

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include "g4root.hh"


#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>

class G4Run;


class BACRunAction : public G4UserRunAction
{
public:
  BACRunAction();
  virtual ~BACRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

private:
  G4String RfileName; 
  TFile *file;
  TTree *Tree;
};

#endif
