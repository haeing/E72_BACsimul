#ifndef BACActionInitialization_hh
#define BACActionInitialization_hh

#include "G4VUserActionInitialization.hh"

class BACActionInitialization : public G4VUserActionInitialization
{
public:
  BACActionInitialization();
  virtual ~BACActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
};

#endif
