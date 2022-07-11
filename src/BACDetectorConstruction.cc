#include "BACDetectorConstruction.hh"
#include "AeroSD.hh"
#include "MPPCSD.hh"

#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnionSolid.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trap.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "TString.h"
#include "TMath.h"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4RotationMatrix.hh"

double theta = 25*TMath::Pi()/180;
G4double r_out = 10*mm;

double sint = TMath::Sin(theta);
double cost = TMath::Cos(theta);

double f(G4double r, G4double z){
  return pow(r*cost+z*sint,2)+2*r_out*pow(1+sint,2)*r-2*r_out*cost*pow(2+sint,2)*z-pow(r_out,2)*(1+sint)*(3+sint);
}

BACDetectorConstruction::BACDetectorConstruction(const G4String &num_aerogel, const G4String &reflect,const G4String &light,const G4String &middle)
  : G4VUserDetectorConstruction(),num_aero(num_aerogel),reflect_part_length(reflect), light_guide_length(light) ,middle_length(middle)
{
}

BACDetectorConstruction::~BACDetectorConstruction()
{
  for (auto visAttributes: fVisAttributes){
    delete visAttributes;
  }
}


G4VPhysicalVolume* BACDetectorConstruction::Construct()
{
  G4int num_aero_i = stoi(num_aero);
  G4double reflect_part_length_d = stod(reflect_part_length);
  G4double light_guide_length_d = stod(light_guide_length);
  G4double middle_length_d = stod(middle_length);
  
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;


  //Material---------------------------------------------------------------

  G4Material* world_mat = nist-> FindOrBuildMaterial("G4_AIR");
  G4Material* Al = nist-> FindOrBuildMaterial("G4_Al");

  G4Element* C = new G4Element("Carbon","C",6.,12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* O = new G4Element("Oxygen","O",15.,15.9994*g/mole);
  //G4Element* Al = new G4Element("Aluminum","Al",13.,26.981538*g/mole);
  G4Element* B = new G4Element("Boron","B",5.,10.811*g/mole);
  G4Element* Na = new G4Element("Na","Na",11.,23.0*g/mole);
  G4Element* Si = new G4Element("Silicon","Si",14.,28.0844*g/mole);
  G4Element* Cl = new G4Element("Chlorine","Cl",17.,35.453*g/mole);
  G4Element* K = new G4Element("Potassium","K",19.,39.093*g/mole);
  

  G4Material *Aerogel = new G4Material("Aerogel",0.2000*g/cm3,2);
  Aerogel->AddElement(Si,1);
  Aerogel->AddElement(O,2);

  G4Material* blacksheet = new G4Material("blacksheet",0.95*g/cm3,2);
  blacksheet->AddElement(C,1);
  blacksheet->AddElement(H,2);

  G4Material *Mylar = new G4Material("Mylar", 1.39*g/cm3, 3);
  Mylar->AddElement(C, 5);
  Mylar->AddElement(H, 4);
  Mylar->AddElement(O, 2);

  G4Material* Tyvek = new G4Material("Tyvek",0.38*g/cm3,2);
  Tyvek->AddElement(C,1);
  Tyvek->AddElement(H,2);

  G4Material* Epoxi = new G4Material("Epoxi",1.1*g/cm3,4);
  Epoxi->AddElement(C, 21);
  Epoxi->AddElement(H, 25);
  Epoxi->AddElement(O,  5);
  Epoxi->AddElement(Cl, 1);



  //Property--------------------------------------------------------------

  //air property---------------------------------------------------------

   
  G4MaterialPropertiesTable* prop_air = new G4MaterialPropertiesTable();
  G4double air_ep[] = {1.3*eV,7.*eV};
  G4double air_rindex[] = {1.0,1.0};
  prop_air->AddProperty("RINDEX",air_ep,air_rindex,2)->SetSpline(true);
  world_mat->SetMaterialPropertiesTable(prop_air);


  //Aerogel Property------------------------------------------------------
  G4MaterialPropertiesTable* prop_aerogel = new G4MaterialPropertiesTable();

  

  G4double aerogel_ep[] = {1.3*eV,7.*eV};
  G4double aerogel_abs[] = {15*mm,15*mm};
  //G4double aerogel_abs[] = {50*cm,50*cm};
  G4double aerogel_rindex[]={1.10,1.10};
  //G4double aerogel_rindex[]={1.03,1.03};
  G4double aerogel_ray[] = {6.16*pow(10,10),6.16*pow(10,10)};

  assert(sizeof(aerogel_ep_abs)==sizeof(aerogel_abs));
  //const G4int num_aerogel = sizeof(aerogel_ep_abs)/sizeof(G4double);
  
  prop_aerogel->AddProperty("RINDEX",aerogel_ep,aerogel_rindex,2)->SetSpline(true);
  prop_aerogel->AddProperty("ABSLENGTH",aerogel_ep,aerogel_abs,2)->SetSpline(true);
  prop_aerogel->AddProperty("RAYLEIGH",aerogel_ep,aerogel_ray,2)->SetSpline(true);
  //prop_aerogel->AddProperty("FASTCOMPONENT",scin_ep1,scin_fast,numentries_scin1)->SetSpline(true);
  //prop_aerogel->AddProperty("SLOWCOMPONENT",scin_ep1,scin_fast,numentries_scin1)->SetSpline(true);
  //prop_aerogel->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);

  prop_aerogel->AddConstProperty("RESOLUTIONSCALE",1.0);
  prop_aerogel->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  prop_aerogel->AddConstProperty("SLOWTIMECONSTANT",2.8*ns);
  prop_aerogel->AddConstProperty("YIELDRATIO",0.8);
  Aerogel->SetMaterialPropertiesTable(prop_aerogel);


  //Black sheet property-------------------------------------------------
  G4MaterialPropertiesTable* prop_bs = new G4MaterialPropertiesTable();

  G4double bs_ep[] = {1.3*eV,7.*eV};
  G4double bs_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4double bs_rindex[]={1.6,1.6};
  
  prop_bs->AddProperty("RINDEX",bs_ep,bs_rindex,2)->SetSpline(true);
  prop_bs->AddProperty("ABSLENGTH",bs_ep,bs_abs,2)->SetSpline(true);
  blacksheet->SetMaterialPropertiesTable(prop_bs);


  //Mylar property
  G4double mylar_ep[]={1.4*eV,1.48*eV,1.52*eV,1.56*eV,1.6*eV,1.7*eV,1.8*eV,1.9*eV,2*eV,2.2*eV,2.4*eV,2.6*eV,2.8*eV,3*eV,3.4*eV,3.8*eV,4*eV,5*eV,6*eV,7*eV};

  G4double mylar_real[]={2.2802,2.6945,2.7668,2.7675,2.6154,2.1606,1.8301,1.5724,1.366,1.0728,0.8734,0.7278,0.6079,0.52135,0.39877,0.31474,0.28003,0.18137,0.12677,0.094236};
  G4double mylar_ima[]={8.1134,8.1878,8.2573,8.3866,8.4914,8.3565,8.0601,7.7354,7.4052,6.7839,6.2418,5.7781,5.3676,5.0008,4.3957,3.9165,3.7081,2.9029,2.3563,1.9519};

  G4double mylar_ep1[] = {1.3*eV,7.*eV};

  assert (sizeof(mylar_ep) == sizeof(mylar_real));
  assert (sizeof(mylar_ep) == sizeof(mylar_ima));
  
  const G4int numentries_mylar = sizeof(mylar_ep)/sizeof(G4double);
  G4double mylar_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4MaterialPropertiesTable* prop_mylar = new G4MaterialPropertiesTable();
  prop_mylar->AddProperty("REALRINDEX",mylar_ep,mylar_real,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("IMAGINARYRINDEX",mylar_ep,mylar_ima,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("ABSLENGTH",mylar_ep1,mylar_abs,2)->SetSpline(true);
  Mylar->SetMaterialPropertiesTable(prop_mylar);

  //Tyvek property
  G4MaterialPropertiesTable* prop_tyvek = new G4MaterialPropertiesTable();
  
  G4double tyvek_rindex[2]={1.5,1.5};
  G4double tyvek_ep[] = {1.6*eV,7.*eV};
  G4double tyvek_abs[] = {1.0e-9*cm,1.0e-9*cm};

  prop_tyvek->AddProperty("RINDEX",tyvek_ep,tyvek_rindex,2)->SetSpline(true);
  prop_tyvek->AddProperty("ABSLENGTH",tyvek_ep,tyvek_abs,2)->SetSpline(true);
  Tyvek->SetMaterialPropertiesTable(prop_tyvek);

  //MPPC property
  G4MaterialPropertiesTable* prop_mppc = new G4MaterialPropertiesTable();
  
  G4double mppc_rindex[2]={1.,1.};
  G4double mppc_ep[] = {1.6*eV,7.*eV};
  G4double mppc_abs[] = {1.0*cm,1.0*cm};
  //G4double mppc_abs[] = {1.0*cm,1.0*cm};

  prop_mppc->AddProperty("RINDEX",mppc_ep,mppc_rindex,2)->SetSpline(true);
  prop_mppc->AddProperty("ABSLENGTH",mppc_ep,mppc_abs,2)->SetSpline(true);
  Al->SetMaterialPropertiesTable(prop_mppc);

  //MPPC surface (Epoxi) property
  G4MaterialPropertiesTable* prop_epoxi = new G4MaterialPropertiesTable();
 
  G4double epoxi_rindex[2]={1.5,1.5};
  G4double epoxi_ep[] = {1.3*eV,7.*eV};
  G4double epoxi_abs[] = {1.0*cm,1.0*cm};


  prop_epoxi->AddProperty("RINDEX",epoxi_ep,epoxi_rindex,2)->SetSpline(true);
  prop_epoxi->AddProperty("ABSLENGTH",epoxi_ep,epoxi_abs,2)->SetSpline(true);
  Epoxi->SetMaterialPropertiesTable(prop_epoxi);
  

  
  //Geometry---------------------------------------------------------------

  //World------------------------------------------------------------------
  G4double world_size = 1*m;
  G4Box* solidWorld = new G4Box("World",world_size, world_size, world_size); 
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(), logicWorld, "World",0,false,0,checkOverlaps);


  //Size-----------------------------------------------------------------------
  G4double extra = 10*mm;
  G4double Aerox = 125.0 *mm+extra;
  G4double Aeroy = 125.0 *mm+extra;
  G4double Aeroz = 12.0 *mm*num_aero_i;

  G4double Aerox_real = 125.0 *mm;
  G4double Aeroy_real = 125.0 *mm;
  G4double Aeroz_real = 12.0 *mm*num_aero_i;

  G4double reflect_thick = 0.3*mm;
  G4double mppc_thick = 1*mm;
  G4double air_thin = 0.2*mm;

  //G4double empty_part1_z = 12.0*num_aero_i*mm;
  G4double empty_part1_z = 0;
  G4double empty_part2_z = reflect_part_length_d*mm;
  //G4double empty_part2_z = Aeroy+air_thin*2;

  //original
  G4double trd_dxa = 2.4*cm;    //-z position x length
  //G4double trd_dxa = 4.8*cm;    //-z position x length
  G4double trd_dxb = Aerox*0.5+air_thin;
  //G4double trd_dxb = Aerox+air_thin*2;
  G4double trd_dya = 2.4*cm;
  G4double trd_dyb = empty_part2_z;                        
  G4double trd_dz  = light_guide_length_d*mm;


  G4int numRZ = 33;


  G4double extra_trap = -10*mm;
  G4double top_trap =12.0*mm*num_aero_i+extra_trap;
  G4double bottom_trap = empty_part2_z;
  G4double middle_trap = middle_length_d*mm;
  G4double height_trap = Aeroy+air_thin*2;
  G4double length_trap = Aerox+air_thin*2;

    



  //Part1-----------------------------------------------------------------------------
  G4Box* part1_cover = new G4Box("part1_cover",Aerox/2+air_thin+reflect_thick,Aeroy/2+air_thin+reflect_thick,(empty_part1_z+air_thin+reflect_thick)/2);
  G4Box* part1_hole = new G4Box("part1_hole",Aerox/2+air_thin,Aeroy/2+air_thin,(empty_part1_z+air_thin)/2);
  G4SubtractionSolid* Part1 = new G4SubtractionSolid("Part1",part1_cover,part1_hole,0,G4ThreeVector(0,0,+reflect_thick/2));

  Part1LW = new G4LogicalVolume(Part1,Mylar,"Part1");
  new G4PVPlacement(0,G4ThreeVector(0,0,-(empty_part1_z+air_thin+reflect_thick)/2),Part1LW,"Part1",logicWorld,false,0,checkOverlaps);
    

  //Part2 - Aerogel-------------------------------------------------------------------
  G4Box* Aero = new G4Box("Aero",Aerox_real/2,Aeroy_real/2,Aeroz_real/2);
  AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
  new G4PVPlacement(0,G4ThreeVector(0,0*mm,38*mm-Aeroz_real/2),AeroLW,"Aero",logicWorld,false,0,checkOverlaps);



  //Part2---------------------------------------------------------------------------

    

  //new
  G4Box* part2_cover = new G4Box("part2_cover",Aerox/2+air_thin+reflect_thick,Aeroy/2+air_thin+reflect_thick+trd_dz,(empty_part2_z+reflect_thick)/2);


  std::vector<G4TwoVector> polygon(5);
  polygon[0].set(-0.5*bottom_trap,0.5*height_trap);
  polygon[1].set(-0.5*bottom_trap,-height_trap*0.5);
  polygon[2].set(0.5*bottom_trap,-height_trap*0.5);
  polygon[3].set(0.5*bottom_trap-middle_trap,0);
  polygon[4].set(0.5*bottom_trap,0.5*height_trap);


  G4TwoVector offsetA(0.,0.), offsetB(0.,0.);
  G4double scaleA=1., scaleB=1.;
  G4ExtrudedSolid* part2_hole1 = new G4ExtrudedSolid("part2_hole1",polygon,length_trap/2,offsetA,scaleA, offsetB, scaleB);




  G4double anx = 0;
  G4double any = -90*degree;
  G4double anz = 180*degree;
    
  G4double sina = TMath::Sin(anx);
  G4double cosa = TMath::Cos(anx);
  G4double sinb = TMath::Sin(any);
  G4double cosb = TMath::Cos(any);
  G4double sinc = TMath::Sin(anz);
  G4double cosc = TMath::Cos(anz);
  G4RotationMatrix *rot = new G4RotationMatrix(G4ThreeVector(cosb*cosc,sina*sinb*cosc-cosa*sinc,cosa*sinb*cosc+sina*sinc),
					       G4ThreeVector(cosb*sinc,sina*sinb*sinc+cosa*cosc,cosa*sinb*sinc-sina*cosc),
					       G4ThreeVector(-sinb,sina*sinb,cosa*cosb));

  G4SubtractionSolid* part2_cover_second = new G4SubtractionSolid("part2_cover_second",part2_cover,part2_hole1,rot,G4ThreeVector(0,0,-reflect_thick*0.5));
  G4RotationMatrix *rotY = new G4RotationMatrix();
  rotY->rotateY(90*degree);


  G4Trd* trd_hole =   new G4Trd("trd_hole",0.5*trd_dxa, 0.5*trd_dxb,0.5*trd_dya, 0.5*trd_dyb, 0.5*(trd_dz+reflect_thick));
  G4RotationMatrix *rotX = new G4RotationMatrix();
  rotX->rotateX(270*degree);
  G4RotationMatrix *rotX90 = new G4RotationMatrix();
  rotX90->rotateX(90*degree);


  //Two MPPC per each side
  
  //light guide bottom part
  G4SubtractionSolid* Part2_1 = new G4SubtractionSolid("Part2",part2_cover_second,trd_hole,rotX90,G4ThreeVector(-Aerox*0.25-air_thin*0.5,-(Aeroy/2+air_thin+reflect_thick/2+trd_dz/2),-reflect_thick*0.5));
  G4SubtractionSolid* Part2_2 = new G4SubtractionSolid("Part2",Part2_1,trd_hole,rotX90,G4ThreeVector(Aerox*0.25+air_thin*0.5,-(Aeroy/2+air_thin+reflect_thick/2+trd_dz/2),-reflect_thick*0.5));
  
  //light guide top part
  G4SubtractionSolid* Part2_3 = new G4SubtractionSolid("Part2",Part2_2,trd_hole,rotX,G4ThreeVector(-Aerox*0.25-air_thin*0.5,Aeroy/2+air_thin+reflect_thick/2+trd_dz/2,-reflect_thick*0.5));
  G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",Part2_3,trd_hole,rotX,G4ThreeVector(Aerox*0.25+air_thin*0.5,Aeroy/2+air_thin+reflect_thick/2+trd_dz/2,-reflect_thick*0.5));

  
  //One MPPC per each side
  /*  
    G4SubtractionSolid* Part2_1 = new G4SubtractionSolid("Part2",part2_cover_second,trd_hole,rotX90,G4ThreeVector(0,-(Aeroy/2+air_thin+reflect_thick/2+trd_dz/2),-reflect_thick*0.5));
    G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",Part2_1,trd_hole,rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick/2+trd_dz/2,-reflect_thick*0.5));
  */



  G4Box* Block = new G4Box("Block",trd_dxb,trd_dz*0.5,1*mm);
  G4LogicalVolume* BlockLW =new G4LogicalVolume(Block,Mylar,"Block");

  

  Part2LW = new G4LogicalVolume(Part2,Mylar,"Part2");

  new G4PVPlacement(0,G4ThreeVector(0,0,(empty_part2_z+reflect_thick)/2),Part2LW,"Part2",logicWorld,false,0,checkOverlaps);
    

  //MPPC---------------------------------------------------------------------------
  G4Box* MPPC = new G4Box("MPPC",8*cm,8*cm,mppc_thick/2);
  //MPPCLW = new G4LogicalVolume(MPPC,Al,"MPPC");
  MPPCLW = new G4LogicalVolume(MPPC,Epoxi,"MPPC");
  //new G4PVPlacement(rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);

  new G4PVPlacement(rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(0,-(Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2),empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);



    




    

  



  //visattributes------------------------------------------------
  auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  visAttributes -> SetVisibility(false);
  logicWorld->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  
  visAttributes = new G4VisAttributes(G4Color::Blue()); 
  AeroLW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  fVisAttributes.push_back(visAttributes);





  //Surface--------------------------------------------------------------
  //black sheet surface------------------
  
  G4OpticalSurface* surface_bs = new G4OpticalSurface("surface_bs");
  surface_bs->SetType(dielectric_dielectric);
  surface_bs->SetModel(unified);
  surface_bs->SetFinish(polishedtyvekair);

  G4MaterialPropertiesTable* sp_bs = new G4MaterialPropertiesTable();

  G4double bs_effi[] = {0.0,0.0};
  G4double bs_reflec[]= {0.05,0.05};
  G4double bs_specularLobe[] = {0.3,0.3};
  G4double bs_specularSpike[] = {0.2,0.2};
  G4double bs_backScatter[] = {0,0};

  sp_bs->AddProperty("EFFICIENCY",bs_ep,bs_effi,2)->SetSpline(true);
  sp_bs->AddProperty("REFLECTIVITY",bs_ep,bs_reflec,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARLOBECONSTANT",bs_ep,bs_specularLobe,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARSPIKECONSTANT",bs_ep,bs_specularSpike,2)->SetSpline(true);
  sp_bs->AddProperty("BACKSCATTERCONSTANT",bs_ep,bs_backScatter,2)->SetSpline(true);
  surface_bs->SetMaterialPropertiesTable(sp_bs);

  //new G4LogicalSkinSurface("bs_surface",BlackLW,surface_bs);

  //mylar_al surface-------------------
  G4OpticalSurface* surface_mylar = new G4OpticalSurface("surface_mylar");
  surface_mylar->SetType(dielectric_metal);
  surface_mylar->SetFinish(polished);    
  surface_mylar->SetModel(unified);
  
  G4MaterialPropertiesTable* sp_mylar = new G4MaterialPropertiesTable();
  G4double mylar_reflec[] = {1.0,1.0};  //for metal, reflectivity is calculated using rindex, they use polarization, angle, energy
  G4double mylar_effi[] = {0.0,0.0};
  G4double mylar_specularLobe[] = {0.85,0.85};
  G4double mylar_specularSpike[]={0.87,0.87};
  G4double mylar_backScatter[] = {0,0};
  sp_mylar->AddProperty("EFFICIENCY",mylar_ep1,mylar_effi,2)->SetSpline(true);
  //sp_mylar->AddProperty("REFLECTIVITY",air_ep,mylar_reflec,2)->SetSpline(true);
  sp_mylar->AddProperty("SPECULARLOBECONSTANT",mylar_ep1,mylar_specularLobe,2)->SetSpline(true);
  sp_mylar->AddProperty("SPECULARSPIKECONSTANT",mylar_ep1,mylar_specularSpike,2)->SetSpline(true);
  sp_mylar->AddProperty("BACKSCATTERCONSTANT",mylar_ep1,mylar_backScatter,2)->SetSpline(true);
  surface_mylar->SetMaterialPropertiesTable(sp_mylar);
  
  //new G4LogicalSkinSurface("mylar_surface",TrdLW,surface_mylar);

  new G4LogicalSkinSurface("mylar_surface",Part1LW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",Part2LW,surface_mylar);



  
  //tyvek surface--------------------------------------------------
  G4OpticalSurface* surface_tyvek = new G4OpticalSurface("surface_tyvek");
  surface_tyvek->SetType(dielectric_dielectric);
  surface_tyvek->SetModel(unified);
  surface_tyvek->SetFinish(groundfrontpainted);



  //MPPC surface--------------------------------------------------
  //G4OpticalSurface* surface_mppc = new G4OpticalSurface("surface_mppc",glisur, ground, dielectric_metal, polished);


  /*
  G4OpticalSurface* surface_mppc = new G4OpticalSurface("surface_mppc");
  //surface_mppc->SetType(dielectric_metal);
  surface_mppc->SetType(dielectric_dielectric);
  surface_mppc->SetFinish(polished);    
  surface_mppc->SetModel(unified);
  G4MaterialPropertiesTable* sp_mppc = new G4MaterialPropertiesTable();
  G4double mppc_reflec[]={0.0,0.0};
  G4double mppc_effi[]={1.0,1.0};
  
  //G4double mppc_ep1[] = {1.38*eV,1.43*eV,1.47*eV,1.51*eV,1.56*eV,1.61*eV,1.66*eV,1.7*eV,1.74*eV,1.79*eV,1.84*eV,1.88*eV,1.93*eV,1.97*eV,2*eV,2.06*eV,2.1*eV,2.15*eV,2.19*eV,2.24*eV,2.3*eV,2.33*eV,2.4*eV,2.47*eV,2.57*eV,2.7*eV,2.85*eV,2.96*eV,3.05*eV,3.15*eV,3.22*eV,3.28*eV,3.35*eV,3.41*eV,3.5*eV,3.57*eV,3.65*eV,3.67*eV,3.71*eV,3.73*eV,3.77*eV,3.79*eV,3.83*eV,3.9*eV};
  //G4double mppc_effi[] = {0.035,0.048,0.057,0.07,0.085,0.098,0.113,0.126,0.142,0.158,0.172,0.191,0.206,0.226,0.243,0.258,0.276,0.294,0.308,0.326,0.344,0.356,0.371,0.385,0.395,0.399,0.392,0.376,0.360,0.342,0.321,0.300,0.278,0.251,0.228,0.201,0.175,0.141,0.120,0.098,0.079,0.059,0.039,0.021};



  //assert (sizeof(mppc_ep1) == sizeof(mppc_effi));
  //const G4int numentries_mppc = sizeof(mppc_ep1)/sizeof(G4double);
    
  sp_mppc->AddProperty("REFLECTIVITY", mppc_ep, mppc_reflec,2);
  sp_mppc->AddProperty("EFFICIENCY", mppc_ep, mppc_effi,2);
  surface_mppc->SetMaterialPropertiesTable(sp_mppc);

  new G4LogicalSkinSurface("mppc_surface",MPPCLW,surface_mppc);
  */  


 


  return physWorld;
}		    


void BACDetectorConstruction::ConstructSDandField()
{
  
  auto aeroSD = new AeroSD("aeroSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(aeroSD);
  AeroLW->SetSensitiveDetector(aeroSD);

  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  MPPCLW->SetSensitiveDetector(mppcSD);

}



  
