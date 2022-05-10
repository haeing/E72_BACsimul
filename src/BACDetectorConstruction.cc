#include "BACDetectorConstruction.hh"
#include "AeroSD.hh"
#include "MPPCSD.hh"

#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnionSolid.hh"
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

BACDetectorConstruction::BACDetectorConstruction()
  : G4VUserDetectorConstruction()
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



  //Property--------------------------------------------------------------

  //air property---------------------------------------------------------

   
  G4MaterialPropertiesTable* prop_air = new G4MaterialPropertiesTable();
  G4double air_ep[] = {1.6*eV,7.*eV};
  G4double air_rindex[] = {1.0,1.0};
  prop_air->AddProperty("RINDEX",air_ep,air_rindex,2)->SetSpline(true);
  world_mat->SetMaterialPropertiesTable(prop_air);


  //Aerogel Property------------------------------------------------------
  G4MaterialPropertiesTable* prop_aerogel = new G4MaterialPropertiesTable();


  G4double aerogel_ep[] = {1.6*eV,7.*eV};
  G4double aerogel_abs[] = {100*mm,100*mm};
  G4double aerogel_rindex[]={1.12,1.12};
  
  prop_aerogel->AddProperty("RINDEX",aerogel_ep,aerogel_rindex,2)->SetSpline(true);
  prop_aerogel->AddProperty("ABSLENGTH",aerogel_ep,aerogel_abs,2)->SetSpline(true);
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

  G4double bs_ep[] = {1.6*eV,7.*eV};
  G4double bs_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4double bs_rindex[]={1.6,1.6};
  
  prop_bs->AddProperty("RINDEX",bs_ep,bs_rindex,2)->SetSpline(true);
  prop_bs->AddProperty("ABSLENGTH",bs_ep,bs_abs,2)->SetSpline(true);
  blacksheet->SetMaterialPropertiesTable(prop_bs);


  //Mylar property
  G4double mylar_ep[]={1.4*eV,1.48*eV,1.52*eV,1.56*eV,1.6*eV,1.7*eV,1.8*eV,1.9*eV,2*eV,2.2*eV,2.4*eV,2.6*eV,2.8*eV,3*eV,3.4*eV,3.8*eV,4*eV,5*eV,6*eV,7*eV};

  G4double mylar_real[]={2.2802,2.6945,2.7668,2.7675,2.6154,2.1606,1.8301,1.5724,1.366,1.0728,0.8734,0.7278,0.6079,0.52135,0.39877,0.31474,0.28003,0.18137,0.12677,0.094236};
  G4double mylar_ima[]={8.1134,8.1878,8.2573,8.3866,8.4914,8.3565,8.0601,7.7354,7.4052,6.7839,6.2418,5.7781,5.3676,5.0008,4.3957,3.9165,3.7081,2.9029,2.3563,1.9519};

  G4double mylar_ep1[] = {1.6*eV,7.*eV};

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
  //G4double mppc_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4double mppc_abs[] = {1.0*cm,1.0*cm};

  prop_mppc->AddProperty("RINDEX",mppc_ep,mppc_rindex,2)->SetSpline(true);
  prop_mppc->AddProperty("ABSLENGTH",mppc_ep,mppc_abs,2)->SetSpline(true);
  Al->SetMaterialPropertiesTable(prop_mppc);

  

  
  //Geometry---------------------------------------------------------------

  //World------------------------------------------------------------------
  G4double world_size = 1*m;
  G4Box* solidWorld = new G4Box("World",world_size, world_size, world_size); 
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(), logicWorld, "World",0,false,0,checkOverlaps);

  
  //Size-------------------------------------------------------------------
  if(version==1){

    G4double Aerox = 150.0 *mm;
    G4double Aeroy = 150.0 *mm;
    G4double Aeroz = 20.0 *mm;
    G4double empty_length = 50.0 *mm;
    G4double reflect_thick = 0.3*mm;

    G4double Black_thick = 0.5 *mm;

    G4double trd_dxa = 5*cm;    //-z position x length
    G4double trd_dxb = Aeroy;          //+z position x length 
    G4double trd_dya = 3*cm;                        
    G4double trd_dyb = Aeroz+empty_length;                        
    G4double trd_dz  = 5*cm;

    G4double mppc_thick = 0.5*cm;

    //Aerogel---------------------------------------------------------------
    G4Box* Aero = new G4Box("Aero",Aerox/2,Aeroy/2,Aeroz/2);
    AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
    new G4PVPlacement(0,G4ThreeVector(0,0,-empty_length/2+reflect_thick), AeroLW, "Aero",logicWorld,false,0,checkOverlaps);

    //Black sheet----------------------------------------------------------
    G4Box* Black_hole = new G4Box("black_hole",Aerox/2,Aeroy/2+Black_thick,(Aeroz+empty_length)/2);
    G4Box* Black_cover = new G4Box("black_cover",Aerox/2+Black_thick,Aeroy/2,(Aeroz+empty_length)/2+Black_thick);

    G4SubtractionSolid* BlackSL = new G4SubtractionSolid("Black",Black_cover,Black_hole,0,G4ThreeVector(0,2*Black_thick,0));
    BlackLW = new G4LogicalVolume(BlackSL,Tyvek,"Black");
    new G4PVPlacement(0,G4ThreeVector(0,-Black_thick,0), BlackLW, "Black", logicWorld,false,0,checkOverlaps);

  
  


    //Reflector----------------------------------------------------------
    G4Box* trdworld = new G4Box("trdworld",trd_dxb/2+reflect_thick,trd_dyb/2+reflect_thick,trd_dz/2);
    trdworldLW = new G4LogicalVolume(trdworld, world_mat, "trdworld");
    G4RotationMatrix *rotX = new G4RotationMatrix();
    rotX->rotateX(270*degree);
    G4VPhysicalVolume* trdworldPW = new G4PVPlacement(rotX,G4ThreeVector(0,(trd_dz+Aerox)/2,0), trdworldLW, "trdworld",logicWorld,false,0,checkOverlaps);
  
    G4Trd* trd_hole =   new G4Trd("trd_hole",0.5*trd_dxa, 0.5*trd_dxb,0.5*trd_dya, 0.5*trd_dyb, 0.5*trd_dz+0.1*mm);
    G4Trd* trd_cover =   new G4Trd("trd_hole",0.5*trd_dxa+reflect_thick, 0.5*trd_dxb+reflect_thick,0.5*trd_dya+reflect_thick, 0.5*trd_dyb+reflect_thick, 0.5*trd_dz);

    G4SubtractionSolid* TrdSL = new G4SubtractionSolid("Trd",trd_cover,trd_hole,0,G4ThreeVector());

    TrdLW = new G4LogicalVolume(TrdSL,Mylar,"Trd");
    //TrdLW = new G4LogicalVolume(TrdSL,Tyvek,"Trd");
    new G4PVPlacement(0,G4ThreeVector(), TrdLW, "Trd", trdworldLW,false,0,checkOverlaps);


    G4Box* UpRefl = new G4Box("UpRefl",Aerox/2,Aeroy/2,reflect_thick/2);
    //UpReflLW = new G4LogicalVolume(UpRefl,Mylar,"UpRefl");
    UpReflLW = new G4LogicalVolume(UpRefl,Tyvek,"UpRefl");
    new G4PVPlacement(0,G4ThreeVector(0,0,-(Aeroz+empty_length-reflect_thick)/2), UpReflLW, "UpRefl", logicWorld,false,0,checkOverlaps);

    G4Box* DownRefl = new G4Box("DownRefl",Aerox/2,Aeroy/2,reflect_thick/2);
    //DownReflLW = new G4LogicalVolume(DownRefl,Mylar,"DownRefl");
    DownReflLW = new G4LogicalVolume(DownRefl,Tyvek,"DownRefl");
    new G4PVPlacement(0,G4ThreeVector(0,0,(Aeroz+empty_length-reflect_thick)/2), DownReflLW, "DownRefl", logicWorld,false,0,checkOverlaps);

    //MPPC----------------------------------------------------------------
    G4Box* MPPC = new G4Box("MPPC",trd_dxa/2,trd_dya/2,mppc_thick/2);
    MPPCLW = new G4LogicalVolume(MPPC,Al,"MPPC");
  new G4PVPlacement(rotX,G4ThreeVector(0,Aeroy/2+trd_dz+mppc_thick/2,0),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);
  }

  if(version==2){
    //Size-----------------------------------------------------------------------
    G4double Aerox = 150.0 *mm;
    G4double Aeroy = 130.0 *mm;
    G4double Aeroz = 20.0 *mm;
    G4double reflect_thick = 0.3*mm;
    //G4double reflect_thick = 1*cm;
    G4double mppc_thick = 0.5*cm;
    G4double air_thin = 0.2*mm;
    //G4double air_thin = 1*cm;
    G4double empty_part1_z = 5*cm;
    G4double empty_part2_z = 10*cm;

    G4double trd_dxa = 5.0*cm;    //-z position x length
    G4double trd_dxb = Aerox+air_thin*2;
    G4double trd_dya = 3*cm;                        
    G4double trd_dyb = empty_part2_z;                        
    G4double trd_dz  = 5*cm;

    G4double top_trap =3.0*cm;
    G4double bottom_trap = empty_part2_z;
    G4double height_trap = Aeroy+air_thin*2;
    G4double length_trap = Aerox+air_thin*2;


    //Trapezoid position correction
    G4double alpha = 0.5*(top_trap-bottom_trap)/height_trap;
    G4double move = 0.5*height_trap*TMath::Tan(alpha);
    

    //Part1-----------------------------------------------------------------------------
    G4Box* part1_cover = new G4Box("part1_cover",Aerox/2+air_thin+reflect_thick,Aeroy/2+air_thin+reflect_thick,(Aeroz+empty_part1_z+air_thin+reflect_thick)/2);
    G4Box* part1_hole = new G4Box("part1_hole",Aerox/2+air_thin,Aeroy/2+air_thin,(Aeroz+empty_part1_z+air_thin)/2);
    G4SubtractionSolid* Part1 = new G4SubtractionSolid("Part1",part1_cover,part1_hole,0,G4ThreeVector(0,0,+reflect_thick/2));
    Part1LW = new G4LogicalVolume(Part1,Mylar,"Part1");
    new G4PVPlacement(0,G4ThreeVector(0,0,-(Aeroz+empty_part1_z+air_thin+reflect_thick)/2),Part1LW,"Part1",logicWorld,false,0,checkOverlaps);
    

    //Part1 - Aerogel-------------------------------------------------------------------
    G4Box* Aero = new G4Box("Aero",Aerox/2,Aeroy/2,Aeroz/2);
    AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
    new G4PVPlacement(0,G4ThreeVector(0,0,-empty_part1_z-Aeroz/2),AeroLW,"Aero",logicWorld,false,0,checkOverlaps);


    //Part2---------------------------------------------------------------------------
    G4Box* part2_cover = new G4Box("part2_cover",Aerox/2+air_thin+reflect_thick,Aeroy/2+air_thin+reflect_thick+trd_dz/2,(empty_part2_z+reflect_thick)/2);
    G4Trap *part2_hole1 = new G4Trap("part2_hole1",length_trap,height_trap,bottom_trap,top_trap);


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
    G4SubtractionSolid* part2_cover_second = new G4SubtractionSolid("part2_cover_second",part2_cover,part2_hole1,rot,G4ThreeVector(0,-trd_dz/2,-reflect_thick/2+move));
    

    G4Trd* trd_hole =   new G4Trd("trd_hole",0.5*trd_dxa, 0.5*trd_dxb,0.5*trd_dya, 0.5*trd_dyb, 0.5*(trd_dz+reflect_thick));
    G4RotationMatrix *rotX = new G4RotationMatrix();
    rotX->rotateX(270*degree);
    G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",part2_cover_second,trd_hole,rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick/2,-reflect_thick*0.5));
    Part2LW = new G4LogicalVolume(Part2,Mylar,"Part2");
    new G4PVPlacement(0,G4ThreeVector(0,trd_dz/2,(empty_part2_z+reflect_thick)/2),Part2LW,"Part2",logicWorld,false,0,checkOverlaps);
    }
  
  



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
  //new G4LogicalSkinSurface("mylar_surface",UpReflLW,surface_mylar);
  //new G4LogicalSkinSurface("mylar_surface",DownReflLW,surface_mylar);
  //new G4LogicalSkinSurface("bs_surface",BlackLW,surface_mylar);

  
  //tyvek surface--------------------------------------------------
  G4OpticalSurface* surface_tyvek = new G4OpticalSurface("surface_tyvek");
  surface_tyvek->SetType(dielectric_dielectric);
  surface_tyvek->SetModel(unified);
  surface_tyvek->SetFinish(groundfrontpainted);
  //new G4LogicalSkinSurface("tyvek_surface",TrdLW,surface_tyvek);
  //new G4LogicalSkinSurface("tyvek_surface",UpReflLW,surface_tyvek);
  //new G4LogicalSkinSurface("tyvek_surface",DownReflLW,surface_tyvek);
  //new G4LogicalSkinSurface("bs_surface",BlackLW,surface_tyvek);


  //MPPC surface--------------------------------------------------
  //G4OpticalSurface* surface_mppc = new G4OpticalSurface("surface_mppc",glisur, ground, dielectric_metal, polished);
  /*
  G4OpticalSurface* surface_mppc = new G4OpticalSurface("surface_mppc");
  G4MaterialPropertiesTable* sp_mppc = new G4MaterialPropertiesTable();
  G4double mppc_reflec[]={0.0,0.0};
  //G4double mppc_ep[] = {1.6*eV,7.*eV};
  G4double mppc_effi[] = {1.,1.};

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
  //AeroLW->SetSensitiveDetector(aeroSD);

  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  //MPPCLW->SetSensitiveDetector(mppcSD);

}



  
  
