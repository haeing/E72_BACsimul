#ifndef BACDetectorConstruction_hh
#define BACDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4VisAttributes;

class BACDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  BACDetectorConstruction();
  virtual ~BACDetectorConstruction();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:
  G4LogicalVolume* AeroLW;
  G4LogicalVolume* BlackLW;
  G4LogicalVolume* trdworldLW;
  G4LogicalVolume* TrdLW;
  G4LogicalVolume* UpReflLW;
  G4LogicalVolume* DownReflLW;

  std::vector<G4VisAttributes*> fVisAttributes;

  
};

#endif
