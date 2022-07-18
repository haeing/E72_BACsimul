#ifndef BACDetectorConstruction_hh
#define BACDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;

class BACDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  BACDetectorConstruction(const G4String &num_aerogel, const G4String &th1, const G4String &th2, const G4String &th3);
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
  G4LogicalVolume* MPPCLW;
  G4LogicalVolume* WinstonLW;
  G4LogicalVolume* CCPCLW;


  G4LogicalVolume* Part1LW;
  G4LogicalVolume* Part2LW;
  G4LogicalVolume* CheckLW;
  G4LogicalVolume* ReflectLW;
  G4LogicalVolume* ReflectBLW;
  G4LogicalVolume* DetectLW;


  std::vector<G4VisAttributes*> fVisAttributes;

  const int version = 2;

  G4String num_aero;
  G4String theta1_put;
  G4String theta2_put;
  G4String theta3_put;


  
};

#endif
