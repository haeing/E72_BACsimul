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
  G4double aerogel_abs[] = {40*mm,40*mm};
  //G4double aerogel_abs[] = {50*cm,50*cm};
  G4double aerogel_rindex[]={1.10,1.10};
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
  G4double mppc_abs[] = {1.0e-9*cm,1.0e-9*cm};
  //G4double mppc_abs[] = {1.0*cm,1.0*cm};

  prop_mppc->AddProperty("RINDEX",mppc_ep,mppc_rindex,2)->SetSpline(true);
  prop_mppc->AddProperty("ABSLENGTH",mppc_ep,mppc_abs,2)->SetSpline(true);
  Al->SetMaterialPropertiesTable(prop_mppc);

  

  
  //Geometry---------------------------------------------------------------

  //World------------------------------------------------------------------
  G4double world_size = 1*m;
  G4Box* solidWorld = new G4Box("World",world_size, world_size, world_size); 
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(), logicWorld, "World",0,false,0,checkOverlaps);

  


  //Size-----------------------------------------------------------------------
  G4double Aerox = 125.0 *mm;
  G4double Aeroy = 125.0 *mm;
  G4double Aeroz = 24.0 *mm;
  G4double reflect_thick = 0.3*mm;
  //G4double reflect_thick = 0.5*cm;
  G4double mppc_thick = 1*mm;
  G4double air_thin = 0.2*mm;
  //G4double air_thin = 1*cm;
  G4double empty_part1_z = 0*cm;
  G4double empty_part2_z = 125*mm;
  //G4double empty_part2_z = Aeroy+air_thin*2;

  //original
  /*
    G4double trd_dxa = 2.0*cm;    //-z position x length
    G4double trd_dxb = Aerox+air_thin*2;
    G4double trd_dya = 2*cm;                        
    G4double trd_dyb = empty_part2_z;                        
    G4double trd_dz  = 5*cm;
  */

  G4int numRZ = 125;
  G4double trd_dxa = 1.2*cm;    //-z position x length
  G4double trd_dxb = Aerox*0.5+air_thin;
  G4double trd_dya = 1.2*cm;                        
  G4double trd_dyb = empty_part2_z;                        
  //G4double trd_dz  = (numRZ-2)*mm-reflect_thick;
  G4double trd_dz  =2*cm;

    
  G4double top_trap =1.0*cm;
  G4double bottom_trap = empty_part2_z;
  G4double height_trap = Aeroy+air_thin*2;
  G4double length_trap = Aerox+air_thin*2;

    


  //Trapezoid position correction
  //G4double alpha = 0.5*(top_trap-bottom_trap)/height_trap;
  //G4double move = 0.5*height_trap*TMath::Tan(alpha);
    
    

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
  //G4Trap *part2_hole1 = new G4Trap("part2_hole1",length_trap,height_trap,bottom_trap,top_trap);

  G4int number = 20;
  G4double step_height[number];
  G4double step_bottom[number];

  G4double a2 = bottom_trap/pow(height_trap,2);
  for(int i=0;i<number;i++){
    step_height[i] = height_trap/(number+1)*(i+1);
    step_bottom[i] = -pow(step_height[i],2)*a2;

  }

  G4double r[numRZ];
  G4double z[numRZ];
  G4double r_in[numRZ];

  std::cout<<"ongoing"<<std::endl;
  for(int i=0;i<numRZ;i++){
    z[i] = i*mm;
    //r_in[i] = 0*mm;
  }
  r[0] = r_out;
  for(int i=1;i<numRZ;i++){
    G4double x = r[i-1];
    while(fabs(f(x,z[i]))>0.1){
      x+=0.0001*mm;
      std::cout<<f(x,z[i])<<std::endl;
    }
    r[i] = x;
    std::cout<<r[i]<<std::endl;
  }


  G4double a = 1*cm;
  G4double b = 5*cm;
  /* 
     std::vector<G4TwoVector> polygon(7);
     polygon[0].set(0,Aeroy*0.5+air_thin);
     polygon[1].set(Aeroy*0.5+air_thin-0.5*b,0.5*b);
     polygon[2].set(Aeroy*0.5+air_thin-0.5*b,Aeroy*0.5+air_thin);
     polygon[3].set(Aeroy*0.5+air_thin-0.5*b+a,Aeroy*0.5+air_thin);
     polygon[4].set(Aeroy*0.5+air_thin-0.5*b+a,-0.5*b);
     polygon[5].set(Aeroy*0.5-0.5*b,-0.5*b);
     polygon[6].set(0,-Aeroy*0.5-air_thin);
  */
    

  std::vector<G4TwoVector> polygon(numRZ+3);
  polygon[0].set(-0.5*bottom_trap,0.5*height_trap);
  polygon[1].set(-0.5*bottom_trap,-height_trap*0.5);
  polygon[2].set(0.5*bottom_trap,-height_trap*0.5);
  for(int i=0;i<numRZ;i++){
    //polygon[i+3].set(0.5*bottom_trap+step_bottom[i],-0.5*height_trap+step_height[i]);
    polygon[i+3].set(0.5*bottom_trap-z[i],-0.5*height_trap+r[i]);
  }

  G4TwoVector offsetA(0.,0.), offsetB(0.,0.);
  G4double scaleA=1., scaleB=1.;
  //G4ExtrudedSolid* test = new G4ExtrudedSolid("test",polygon,length_trap/2,offsetA,scaleA, offsetB, scaleB);
  G4ExtrudedSolid* part2_hole1 = new G4ExtrudedSolid("part2_hole1",polygon,length_trap/2,offsetA,scaleA, offsetB, scaleB);
  //auto testLW = new G4LogicalVolume(test,Aerogel,"test");
  //new G4PVPlacement(0,G4ThreeVector(20*cm, 20*cm, 20*cm),testLW,"test",logicWorld,false,0,checkOverlaps);



  G4double anx1 = 0*degree;
  G4double any1 = -90*degree;
  G4double anz1 = 0*degree;
    
  G4double sina1 = TMath::Sin(anx1);
  G4double cosa1 = TMath::Cos(anx1);
  G4double sinb1 = TMath::Sin(any1);
  G4double cosb1 = TMath::Cos(any1);
  G4double sinc1 = TMath::Sin(anz1);
  G4double cosc1 = TMath::Cos(anz1);

  G4RotationMatrix *rot1 = new G4RotationMatrix(G4ThreeVector(cosb1*cosc1,sina1*sinb1*cosc1-cosa1*sinc1,cosa1*sinb1*cosc1+sina1*sinc1),
						G4ThreeVector(cosb1*sinc1,sina1*sinb1*sinc1+cosa1*cosc1,cosa1*sinb1*sinc1-sina1*cosc1),
						G4ThreeVector(-sinb1,sina1*sinb1,cosa1*cosb1));

    
    
  G4double anx2 = 0*degree;
  G4double any2 = -90*degree;
  G4double anz2 = 180*degree;
    
  G4double sina2 = TMath::Sin(anx2);
  G4double cosa2 = TMath::Cos(anx2);
  G4double sinb2 = TMath::Sin(any2);
  G4double cosb2 = TMath::Cos(any2);
  G4double sinc2 = TMath::Sin(anz2);
  G4double cosc2 = TMath::Cos(anz2);
  G4RotationMatrix *rot2 = new G4RotationMatrix(G4ThreeVector(cosb2*cosc2,sina2*sinb2*cosc2-cosa2*sinc2,cosa2*sinb2*cosc2+sina2*sinc2),
						G4ThreeVector(cosb2*sinc2,sina2*sinb2*sinc2+cosa2*cosc2,cosa2*sinb2*sinc2-sina2*cosc2),
						G4ThreeVector(-sinb2,sina2*sinb2,cosa2*cosb2));


  
  G4SubtractionSolid* part2_cover_first = new G4SubtractionSolid("part2_cover_first",part2_cover,part2_hole1,rot1,G4ThreeVector(0,-trd_dz/2,-reflect_thick*0.5));
  G4SubtractionSolid* part2_cover_second = new G4SubtractionSolid("part2_cover_second",part2_cover_first,part2_hole1,rot2,G4ThreeVector(0,-trd_dz/2,-reflect_thick*0.5));
  
  G4RotationMatrix *rotY = new G4RotationMatrix();
  rotY->rotateY(90*degree);
  //G4SubtractionSolid* part2_cover_second = new G4SubtractionSolid("part2_cover_second",part2_cover,part2_hole1,rotY,G4ThreeVector(0,-trd_dz/2,-(empty_part2_z+reflect_thick)/2));

  //test winstone


  /*
  G4double r[numRZ];
  G4double z[numRZ];
  G4double r_in[numRZ];

  for(int i=0;i<numRZ;i++){
    z[i] = i*mm;
    r_in[i] = 0*mm;
  }
  r[0] = 6*mm;
  for(int i=0;i<numRZ;i++){
    G4double x = r[i-1];
    while(fabs(f(x,z[i]))>0.1){
      x+=0.0001;
    }
    r[i] = x;
  }
  */


  //G4VSolid* poly = new G4Polycone("poly",0.*deg, 360.*deg, numRZ,z, r_in,r);
  //auto* polyLW = new G4LogicalVolume(poly,world_mat,"poly");
    
  //new G4PVPlacement(rotX,G4ThreeVector(20*cm,20*cm,20*cm),polyLW,"poly",logicWorld,false,0,checkOverlaps);
    
    

  G4Trd* trd_hole =   new G4Trd("trd_hole",0.5*trd_dxa, 0.5*trd_dxb,0.5*trd_dya, 0.5*trd_dyb, 0.5*(trd_dz+reflect_thick));
  G4RotationMatrix *rotX = new G4RotationMatrix();
  rotX->rotateX(270*degree);
  //G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",part2_cover_second,trd_hole,rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick/2,-reflect_thick*0.5));

  /*
    G4SubtractionSolid* Part2_1 = new G4SubtractionSolid("Part2_1",part2_cover_second,poly,rotX,G4ThreeVector(-Aerox/3-air_thin*4/3,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,-(empty_part2_z+reflect_thick)*0.25));
    G4SubtractionSolid* Part2_2 = new G4SubtractionSolid("Part2_2",Part2_1,poly,rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,-(empty_part2_z+reflect_thick)*0.25));
    G4SubtractionSolid* Part2_3 = new G4SubtractionSolid("Part2_3",Part2_2,poly,rotX,G4ThreeVector(Aerox/3+air_thin*4/3,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,-(reflect_thick+empty_part2_z)*0.25));
    G4SubtractionSolid* Part2_4 = new G4SubtractionSolid("Part2_4",Part2_3,poly,rotX,G4ThreeVector(-Aerox/3-air_thin*4/3,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,(reflect_thick+empty_part2_z)*0.25));
    G4SubtractionSolid* Part2_5 = new G4SubtractionSolid("Part2_5",Part2_4,poly,rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,(empty_part2_z+reflect_thick)*0.25));
    G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",Part2_5,poly,rotX,G4ThreeVector(Aerox/3+air_thin*4/3,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,(empty_part2_z+reflect_thick)*0.25));
  */
  //G4SubtractionSolid* Part2_1 = new G4SubtractionSolid("Part2_1",part2_cover_second,poly,rotX,G4ThreeVector(-Aerox*0.5/3-air_thin/3,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,-reflect_thick*0.25));

  G4SubtractionSolid* Part2_1 = new G4SubtractionSolid("Part2",part2_cover_second,trd_hole,rotX,G4ThreeVector(-Aerox*0.25-air_thin*0.5,Aeroy/2+air_thin+reflect_thick/2,-reflect_thick*0.5));
  G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",Part2_1,trd_hole,rotX,G4ThreeVector(Aerox*0.25+air_thin*0.5,Aeroy/2+air_thin+reflect_thick/2,-reflect_thick*0.5));
  //G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",Part2_1,poly,rotX,G4ThreeVector(Aerox*0.25+air_thin*0.5,Aeroy/2+air_thin+reflect_thick/2+numRZ*mm*0.5,-reflect_thick*0.5));
  Part2LW = new G4LogicalVolume(Part2,Mylar,"Part2");
  new G4PVPlacement(0,G4ThreeVector(0,trd_dz/2,(empty_part2_z+reflect_thick)/2),Part2LW,"Part2",logicWorld,false,0,checkOverlaps);
    

  //MPPC---------------------------------------------------------------------------
  //G4Box* MPPC = new G4Box("MPPC",trd_dxa/2,trd_dya/2,mppc_thick/2);
  G4Box* MPPC = new G4Box("MPPC",8*cm,8*cm,mppc_thick/2);
  MPPCLW = new G4LogicalVolume(MPPC,Al,"MPPC");
  new G4PVPlacement(rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);


  //Virtual plane
  G4Box* Check = new G4Box("Check",Aerox/2+air_thin-1*mm,0.1*mm,empty_part2_z/2-1*mm);
  CheckLW = new G4LogicalVolume(Check,world_mat,"Check");
  new G4PVPlacement(0,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick/2-5*mm,(empty_part2_z)/2),CheckLW,"Check",logicWorld,false,0,checkOverlaps);


    

  



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
  AeroLW->SetSensitiveDetector(aeroSD);

  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  CheckLW->SetSensitiveDetector(mppcSD);

}



  
