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
  BACDetectorConstruction(const G4String &num_aerogel, const G4String &reflec, const G4String &light);
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


  G4LogicalVolume* Part1LW;
  G4LogicalVolume* Part2LW;
  G4LogicalVolume* CheckLW;


  std::vector<G4VisAttributes*> fVisAttributes;

  const int version = 3;

  G4String num_aero;
  G4String reflect_part_length;
  G4String light_guide_length;

  
};

#endif
