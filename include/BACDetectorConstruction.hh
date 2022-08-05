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
  BACDetectorConstruction(const G4String &version_put, const G4String &num_aerogel, const G4String &parameter1, const G4String &parameter2, const G4String &parameter3);
  virtual ~BACDetectorConstruction();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:

  G4double reflect_part_length_d;
  G4double light_guide_length_d;
  G4double middle_length_d;

  G4double theta1;
  G4double theta2;
  G4double theta3;

  G4double p;
  G4double ref_z;
  G4double ref_theta;

  G4double mppc_theta;
  G4double Dpartz;

  

  
  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physDetect;
  G4LogicalVolume* AeroLW;
  G4LogicalVolume* BlackLW;
  G4LogicalVolume* trdworldLW;
  G4LogicalVolume* TrdLW;
  G4LogicalVolume* UpReflLW;
  G4LogicalVolume* DownReflLW;
  G4LogicalVolume* MPPCLW;
  G4LogicalVolume* WinstonLW;
  G4LogicalVolume* CCPCLW;
  G4LogicalVolume* UpLW;
  G4LogicalVolume* BottomLW;
  G4LogicalVolume* DetLW;


  G4LogicalVolume* Part1LW;
  G4LogicalVolume* Part2LW;
  G4LogicalVolume* CheckLW;
  G4LogicalVolume* ReflectLW;
  G4LogicalVolume* ReflectBLW;
  G4LogicalVolume* DetectLW;
  G4LogicalVolume* SideLW;



  std::vector<G4VisAttributes*> fVisAttributes;

  G4int version;
  //const int version = stoi(version_put);

  //G4String version_put;
  G4String version_in;
  G4String num_aero;
  G4String parameter1;
  G4String parameter2;
  G4String parameter3;
  
  //BACLens *Lens;


  
};

#endif
