// Based on template provided by G4 example B2a.
// Adapted and developed by Vitor Jose Shen
// at department of physics, Tsinghua University.
// This simulation program "nDa" stands for:
// Simulation of "neutron (liquid scintillation) Detector array"
// 06 Apr 2021
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file nDaDetectorConstruction.cc
/// \brief Implementation of the nDaDetectorConstruction class
 
#include "nDaDetectorConstruction.hh"
#include "nDaDetectorMessenger.hh"
#include "nDaTrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include <mutex>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* nDaDetectorConstruction::fMagFieldMessenger = 0;

nDaDetectorConstruction::nDaDetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfChambers(40),
 fNbOfRows(1),
 fNbOfLevels(40),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new nDaDetectorMessenger(this);

  size_t fullNChambers = fNbOfLevels*fNbOfRows*fNbOfChambers;
  fLogicChamber = new G4LogicalVolume*[fullNChambers];
  fPhysChamber = new G4VPhysicalVolume*[fullNChambers]; //!
  fLockTarget = new std::mutex;                         //!
  //fLockTargetShell = new std::mutex;                    //!

  //  fNbOfChambers = 1;									// Altering num#########
  //fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
nDaDetectorConstruction::~nDaDetectorConstruction()
{
  delete [] fPhysChamber;
  delete [] fLogicChamber; 
  delete fStepLimit;
  delete fMessenger;
  delete fLockTarget;
  //delete fLockTargetShell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* nDaDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaDetectorConstruction::DefineMaterials()
{
  // Material definition 
  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager 
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Lead defined using NIST Manager for Target Material as illustration
  //fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");

  // Xenon gas defined using NIST Manager for Chamber Material as illustration
  //fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Xe");

   // Polyethylene using NIST Manager
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  
  //Adding what we needed #########################
  // Hydrogen defined using NIST Manager
  /*
  G4Material* H = nistManager->FindOrBuildMaterial("G4_H");
  // Carbon defined using NIST Manager
  G4Material* C = nistManager->FindOrBuildMaterial("G4_C");
  // Nitrogen
  G4Material* N = nistManager->FindOrBuildMaterial("G4_N");
  // Oxygen
  G4Material* O = nistManager->FindOrBuildMaterial("G4_O");
  // Aluminum
  G4Material* Al = nistManager->FindOrBuildMaterial("G4_Al"); 
  // Concrete
  G4Material* Concrete = nistManager->FindOrBuildMaterial("G4_CONCRETE"); 
  //Iron
  G4Material* Fe = nistManager->FindOrBuildMaterial("G4_Fe");
  */
  //#####
//***********************************************************
  //how it works!
/*  G4NistManager* man = G4NistManager::Instance();
  G4Material* H2O = man->FindOrBuildMaterial("G4_WATER"); 
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");   */
  //G4Material* air  = G4Material::GetMaterial("G4_AIR");
//************************************************************
  //Elements
  //G4double z, a;
  G4double z, a;
  G4String name, symbol;
  a = 1.0079*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen",symbol="H", z=1, a);
  a=12.0107*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6, a);
//  a=14.0067*g/mole;
//  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7, a);
//  a=15.9994*g/mole;
//  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8, a);

  
  //G4double z, a;
  G4double fractionmass, density;
  G4int ncomponents;
  
  //NE213
  density = 0.874*g/cm3;
  G4Material* NE213 = new G4Material(name="NE213",density,ncomponents=2, kStateLiquid);
  NE213->AddElement(elH,fractionmass=54.812*perCent);
  NE213->AddElement(elC, fractionmass=45.188*perCent);
  
  //"Vacuum"
  // define a vacuum with a restgas pressure  typical for accelerators
  G4double const Torr  = atmosphere/760.;         // 1 Torr
  G4double pressure = 10e-9*Torr, temperature = 296.150*kelvin;    // 23  Celsius
  G4Material* Vacuum = new G4Material("Vacuum", z=7., a=14.01*g/mole, density= 1.516784e-11*kg/m3,
                 kStateGas, temperature, pressure);
  
  // Liquid Deutron at 20.15 K
  // ref: https://www.aqua-calc.com/page/density-table/substance/heavy-blank-hydrogen-coma-and-blank-liquid
  auto D2L = new G4Material("D2L", 1.69*g/cm3, 1, kStateLiquid, 20.15*kelvin, CLHEP::STP_Pressure);
  auto elD = new G4Element("D", "D", 1);
  auto isoD = new G4Isotope("D", 1, 2);
  elD->AddIsotope(isoD, 1.);
  D2L->AddElement(elD, 2);
  
  // U238
  auto U238 = new G4Material("U238", 19.1*g/cm3, 1, kStateSolid);
  auto elU238 = new G4Element("U238", "U238", 1);
  auto isoU238 = new G4Isotope("U238", 92, 238);
  elU238->AddIsotope(isoU238, 1.);
  U238->AddElement(elU238, 1);

//***********************************************************
//Setting the target material
//fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Al");
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Fe");
//fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");
//fTargetMaterial  = G4Material::GetMaterial("U238");

//Setting the chamber material
  fChamberMaterial = G4Material::GetMaterial("NE213");

//Setting the tracker material
  fTrackerMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
  //fTrackerMaterial = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");

//***********************************************************


  // Testing for Hydrogen and Cabron
  nistManager->FindOrBuildMaterial("G4_H");
  nistManager->FindOrBuildMaterial("G4_C");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* nDaDetectorConstruction::DefineVolumes()
{
  G4Material* Air  = G4Material::GetMaterial("G4_AIR");
 // G4Material* Al  = G4Material::GetMaterial("G4_Al");
 // G4Material* Concrete  = G4Material::GetMaterial("G4_CONCRETE");
  
  // Previous Setting and altering
//***************************************************
  //G4double chamberSpacing = 80*cm; // from chamber center to center!
  //G4double chamberWidth = 20.0*cm; // width of the chambers  ######
  //G4double targetLength =  5.0*cm; // full length of Target  ######
  //G4double trackerLength = (fNbOfChambers)*chamberSpacing;	//fNbOfChamber+1 #############
  //G4double worldLength = 1.2 * (2*targetLength + trackerLength);   //####
  //G4double worldLength = 1080*cm;
//***************************************************
// 3 sizes of world length setting (Previous version for nD1)
// corresponding to 3 length of tracker to set the source distance
  	//G4double worldLength = 2*106.35*cm;
  	//G4double worldLength = 2*112.7*cm;
  	//G4double worldLength = 2*125.4*cm;
// **************************************************

// Starting here!
// Sizes of the principal geometrical components (solids)
// For nDa the distance between source and array is 10m, and with material is 1m
     G4double worldLength = 5000*cm;    //50m
     G4double chamberDiam = 5.*cm;     //50mm
     G4double chamberLength = 10.*cm;  //100mm
     
  G4double chamberSpacing = chamberDiam + 0. *cm; // from chamber center to center!
  G4double chamberRowSpacing = chamberLength + 0.*cm;
  G4double chamberLevelSpacing = chamberDiam + 0.*cm;
  
//*************************************************
  // Target material ball's radius equals 1 to 7 cm
  // setting Radius of material ball
  // setting radius of Tub shape Hole

  //G4double targetRadius = 1.*cm;
  G4double targetRadius = 2.*cm;
  //G4double targetRadius = 3.*cm;
  //G4double targetRadius = 4.*cm;
  //G4double targetRadius = 5.*cm;
  //G4double targetRadius = 6.*cm;
  //G4double targetRadius = 7.*cm;

  // Half length og targetBox
  G4double targetBoxHalfXYZ = 6. *cm;
//**************************************************

  G4double targetDiam = 2*targetRadius;
  // Target Placement  9m + 5cm =905 cm
  G4double targetPosZ = -905.*cm;

  //G4double chamberWidth = 20.0*cm; // width of the chambers  ######
  //G4double chamberSpacing = 80*cm; // from chamber center to center!
  //G4double targetLength =  5.0*cm; // full length of Target  ######
//**********************************************
//	G4double InnerL = 4.6*cm;
//	G4double OutterL = 6.6*cm;

  //G4double targetRadius  = 0.5*targetLength;   // Radius of Target		   ##########
  //targetLength = 0.5*targetLength;             // Half length of the Target  ##########
  //G4double trackerSize   = 0.5*trackerLength;  // Half length of the Tracker  which is 200cm=2m
  
  // for single unit (previous version) NE213: length=diameter || 2*length=diameter
  //********************************************
  //G4double trackerDiam = 12.7*cm;
  //G4double trackerDiam = 25.4*cm;
  //********************************************

  // For neutron Detector array
  // the actual detector would be chamber
  // the tracker just for testing
  // the diameter of chamber would be 50mm = 5cm exactly
  // the length of chamber would be 100mm = 10cm
  // tracker would envople those chambers
  // tracker can be chosen as Box
  // x = y = 2m 
  // z = 10cm

  G4double trackerX = 200*cm;   //2m
  G4double trackerY = 200*cm;   //2m
  G4double trackerZ = 10*cm;    //100mm

  //********************************************
  //G4double trackerLength = 12.7*cm;
  //G4double trackerLength = 25.4*cm;
  //G4double trackerLength = 50.8*cm;
  //********************************************
  //G4double trackerRadius = 0.5*trackerDiam;
  //G4double trackerSize = 0.5*trackerLength;
  //*********************************************


  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;												//#####

 
  G4Box* worldS
    = new G4Box("world",                                    //its name
                worldLength/2,worldLength/2,worldLength/2); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 Air,      //its material
                 "World"); //its name
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps           

  // Previous version Target
 //******************************************************************************* 
 /* G4ThreeVector positionTarget = G4ThreeVector(0,0,-(targetLength+trackerSize));

  G4Tubs* targetS
    = new G4Tubs("target",0.,targetRadius,targetLength,0.*deg,360.*deg);
  fLogicTarget
    = new G4LogicalVolume(targetS, fTargetMaterial,"Target",0,0,0);
  new G4PVPlacement(0,               // no rotation
                    positionTarget,  // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << 2*targetLength/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;							*///#########
  //*******************************************************************************
  //Target updated
  G4ThreeVector targetPosition = G4ThreeVector(0,0,targetPosZ);
  //Define Box and Spheres and Tubs
 
  G4Box* TargetS = new G4Box("targetBox", targetBoxHalfXYZ, targetBoxHalfXYZ, targetBoxHalfXYZ);  //double tubs series1
  //Half xyz same-> 6cm

  //G4Box* TargetS = new G4Box("targetBox", 6.5*cm, targetBoxHalfXYZ, targetBoxHalfXYZ); //THU hole serie
  //Construction of THU holes below
  /*
  G4Box* T1 = new G4Box("targetHoleT1", 3.*cm/2, 2.*cm/2, targetBoxHalfXYZ);
  G4Box* T2 = new G4Box("targetHoleT2", 1.*cm/2, 6.*cm/2, targetBoxHalfXYZ);
  G4Box* H1 = new G4Box("targetHoleH1", 1.*cm/2, 8.*cm/2, targetBoxHalfXYZ);
  G4Box* H2 = new G4Box("targetHoleH2", 1.*cm/2, 2.*cm/2, targetBoxHalfXYZ);
  G4Box* H3 = new G4Box("targetHoleH3", 1.*cm/2, 8.*cm/2, targetBoxHalfXYZ);
  G4Box* U1 = new G4Box("targetHoleU1", 1.*cm/2, 8.*cm/2, targetBoxHalfXYZ);
  G4Box* U2 = new G4Box("targetHoleU2", 1.*cm/2, 2.*cm/2, targetBoxHalfXYZ);
  G4Box* U3 = new G4Box("targetHoleU3", 1.*cm/2, 8.*cm/2, targetBoxHalfXYZ);
  */

  //Construct 2 Tubs for 2 holes sim
  G4Tubs* TargetS1 = new G4Tubs("TargetHole1", 0, targetRadius, 6.*cm, 0.*deg, 360.*deg);
  G4Tubs* TargetS2 = new G4Tubs("TargetHole2", 0, targetRadius, 6.*cm, 0.*deg, 360.*deg);
  

  //Notice: abandon construction below, changed to tubs as shown above
  //Construction of double spheres below
  //G4Sphere* TargetS1 = new G4Sphere("targetSp1", 0., targetRadius,  0.*deg, 360.*deg, 0.*deg, 180.*deg);
  //G4Sphere* TargetS2 = new G4Sphere("targetSp2", 0., targetRadius,  0.*deg, 360.*deg, 0.*deg, 180.*deg);
  // targetRadius = 1cm default 

  //Define y-Rotation and Z-Transformation Action for Solid 2
  G4RotationMatrix* yRot = new G4RotationMatrix;
  yRot->rotateY(0. *rad);
  G4ThreeVector yTrans1(0, 3.*cm,0);
  G4ThreeVector yTrans2(0,-3.*cm,0);
  //Use Subtraction
  G4SubtractionSolid* subtraction = new G4SubtractionSolid("Box-Tub", TargetS, TargetS1, yRot, yTrans1);
  G4SubtractionSolid* FinalTargetS= new G4SubtractionSolid("Target", subtraction, TargetS2, yRot, yTrans2); 
  //Box-Tub*2
  /*G4Tubs* targetS
    = new G4Tubs("Target",0.,1.*cm,1.*cm,0.*deg,360.*deg);*/
  //Udated THU hole shape
  /*
  G4ThreeVector T1Trans(-4.*cm, 3.*cm,0);
  G4ThreeVector T2Trans(-4.*cm,-1.*cm,0);
  G4ThreeVector H1Trans(-1.*cm,0,0);
  G4ThreeVector H2Trans(0,0,0);
  G4ThreeVector H3Trans(1.*cm,0,0);
  G4ThreeVector U1Trans(3.*cm,0,0);
  G4ThreeVector U2Trans(4.*cm,-3.*cm,0);
  G4ThreeVector U3Trans(5.*cm,0,0);
  G4SubtractionSolid* BoxT1 = new G4SubtractionSolid("Box-T1", TargetS, T1, yRot, T1Trans);
  G4SubtractionSolid* BoxT2 = new G4SubtractionSolid("Box-T2", BoxT1, T2, yRot, T2Trans);
  G4SubtractionSolid* BoxH1 = new G4SubtractionSolid("Box-H1", BoxT2, H1, yRot, H1Trans);
  G4SubtractionSolid* BoxH2 = new G4SubtractionSolid("Box-H2", BoxH1, H2, yRot, H2Trans);
  G4SubtractionSolid* BoxH3 = new G4SubtractionSolid("Box-H3", BoxH2, H3, yRot, H3Trans);
  G4SubtractionSolid* BoxU1 = new G4SubtractionSolid("Box-U1", BoxH3, U1, yRot, U1Trans);
  G4SubtractionSolid* BoxU2 = new G4SubtractionSolid("Box-U2", BoxU1, U2, yRot, U2Trans);
  G4SubtractionSolid* BoxTHU = new G4SubtractionSolid("THU"  , BoxU2, U3, yRot, U3Trans);
  */

  fSolidTarget
    //=BoxTHU;
    = FinalTargetS;
    //=subtraction;
    //=FinalTargetS;

  fLogicTarget  
    = new G4LogicalVolume(fSolidTarget, fTargetMaterial,"Target",0,0,0);
  
  fPhysTarget = new G4PVPlacement( NULL,
                    targetPosition,
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    true,           // no (or with) boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << targetDiam/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;						//#########

  /*
  fSolidTarget
    =new G4Sphere("target", 0., targetRadius, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
  fLogicTarget
    = new G4LogicalVolume(fSolidTarget, fTargetMaterial,"Target",0,0,0);
  fPhysTarget = new G4PVPlacement( NULL,
                    targetPosition,
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  G4cout << "Target is " << targetDiam/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;*/							//#########


  // Tracker
 																		//what about only tracker
  double x=0;
  double y=0;
  double z=0;
//*******************************
// Moving set
  //double x = 1 *cm;         // R1
  //double x = 2 *cm;         // R2
  //double x = -1 *cm;        // L1
  //double x = -2 *cm;        // L2

  //double y= 1 *cm;          // U1
  //double y= 2 *cm;          // U2
  //double y= -1 *cm;         // D1
  //double y= -2 *cm;         // D2

// And also with their combination RU, RD, LU, LD

//End of moving set
//********************************
  G4ThreeVector positionTracker = G4ThreeVector(x,y,z);
//  G4ThreeVector positionTracker = G4ThreeVector(fNbOfChambers%2?0:chamberSpacing/2.,fNbOfLevels%2?0:chamberLevelSpacing/2.,0);
// above?

// previous version
//***********************************************************
/*  G4Tubs* trackerS
    = new G4Tubs("tracker",			//pName: it's name	
    			 0,					//pRmin: inner radius
    			 trackerRadius,		//pRmax: outer radius
    			 trackerSize,		//pDz: half length in z
    			 0.*deg, 			//pSphi: phi angle in radius
    			 360.*deg);			//pDphi: angle of segment of radians */
 //***********************************************************

 // update tracker to Box
  G4VSolid* trackerS
    = new G4Box("tracker", trackerX/2., trackerY/2., trackerZ/2.);

  G4LogicalVolume* trackerLV
    = new G4LogicalVolume(trackerS, fTrackerMaterial, "Tracker",0,0,0);  
  //= new G4LogicalVolume(trackerS, G4Material::GetMaterial("G4_WATER"), "Tracker",0,0,0);  
  
  fPhysTracker = new G4PVPlacement(0,               // no rotation
                    positionTracker, // at (x,y,z)
                    trackerLV,       // its logical volume
                    "Tracker",       // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));  //#####

  worldLV      ->SetVisAttributes(boxVisAtt);
  fLogicTarget ->SetVisAttributes(boxVisAtt);  
  trackerLV    ->SetVisAttributes(boxVisAtt);

  // Tracker segments

  G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
         << G4endl
         << "The chambers are " << chamberDiam/cm << " cm of "
         << fChamberMaterial->GetName() << G4endl
         << "The distance between chamber is " << chamberSpacing/cm << " cm" 
         << G4endl;																 //#####
  
 /* G4double firstPosition = -trackerSize + chamberSpacing;
  G4double firstLength   = trackerLength/10;
  G4double lastLength    = trackerLength;
  G4double halfWidth = 0.5*chamberWidth;
  G4double rmaxFirst = 0.5 * firstLength;
  G4double rmaxIncr = 0.0;*/

  // chamberDiam --> chamberWidth or chamberHeight --> X or Y
  // chamberLength --> chamberDepth  --> Z

  G4double firstPositionX = -trackerX/2. + chamberDiam/2.;
  G4double firstPositionY = -trackerY/2. + chamberDiam/2.;
  G4double firstPositionZ = -trackerZ/2. + chamberLength/2.;



  if( fNbOfChambers > 0 ){
  //  rmaxIncr =  0.5 * (lastLength-firstLength)/(fNbOfChambers-1);
    if (chamberSpacing  < chamberDiam) {
       G4Exception("nDaDetectorConstruction::DefineVolumes()",
                   "InvalidSetup", FatalException,
                   "Diameter>Spacing");
    }
  }																			 //########

  for (G4int copyNo=0; copyNo<fNbOfLevels*fNbOfRows*fNbOfChambers; copyNo++) {			//#####################
// copyNo = nbLevel * fNbOfRows * fNbOfChambers
	//        + nbRow * fNbOfChambers
	//        + nbInRow
      G4int nbLevel = copyNo/(fNbOfRows*fNbOfChambers);
      G4int nbInLevel = copyNo-nbLevel*fNbOfRows*fNbOfChambers;
      G4int nbRow = nbInLevel/fNbOfChambers;
      G4int nbInRow = nbInLevel%fNbOfChambers;
      G4int isOdd = nbRow%2;
      if(isOdd && (nbInRow==fNbOfChambers-1)) 
      continue; // skip the last one for odd Rows
      G4double Xposition = firstPositionX + (nbInRow + isOdd/2.) * chamberSpacing; //###################
      G4double Yposition = firstPositionY + nbLevel * chamberLevelSpacing;
      G4double Zposition = firstPositionZ + (nbRow) * chamberRowSpacing;

      G4Tubs* chamberS
    = new G4Tubs("Chamber",	//pName: it's name	
    			 0,					      //pRmin: inner radius
    			 chamberDiam/2,		//pRmax: outer radius
    			 chamberLength/2,		//pDz: half length in z
    			 0.*deg, 			    //pSphi: phi angle in radius
    			 360.*deg);		  	//pDphi: angle of segment of radians */

//*****************************************************************************
      /*G4VSolid* chamberS
	=new G4Box("Chamber_solid", chamberDiam/2., chamberDiam/2., chamberLength/2.);*/
//*****************************************************************************

      fLogicChamber[copyNo] =											//#################
              new G4LogicalVolume(chamberS,fChamberMaterial,"Chamber_LV",0,0,0);

      fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);			//#################

      fPhysChamber[copyNo]=new G4PVPlacement(0,       // no rotation
                        G4ThreeVector(Xposition,Yposition,Zposition), // at (x,y,z)
                        fLogicChamber[copyNo],        // its logical volume
                        "Chamber_PV",                 // its name
                        trackerLV,                    // its mother  volume
                        false,                        // no boolean operations
                        copyNo,                       // copy number   ###################
                        fCheckOverlaps);              // checking overlaps 
  }			
      //G4double Zposition = firstPosition + copyNo * chamberSpacing; //###################
      //G4double Zposition = firstPosition + chamberSpacing;
      //G4double rmax =  rmaxFirst + copyNo * rmaxIncr;				  //###################
//      G4double rmax =  rmaxFirst + rmaxIncr;

//      G4Tubs* chamberS
//        = new G4Tubs("Chamber_solid", 0, rmax, halfWidth, 0.*deg, 360.*deg);

      //fLogicChamber[copyNo] =											//#################
      //        new G4LogicalVolume(chamberS,fChamberMaterial,"Chamber_LV",0,0,0);

    //  fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);			//#################

//      new G4PVPlacement(0,                            // no rotation
//                        G4ThreeVector(0,0,Zposition), // at (x,y,z)
//                        fLogicChamber
                        //fLogicChamber[copyNo],        // its logical volume
//                        "Chamber_PV",                 // its name
//                        trackerLV,                    // its mother  volume
//                        false,                        // no boolean operations
//                        0,
                        //copyNo,                       // copy number   ###################
//                        fCheckOverlaps);              // checking overlaps 

//  }																	#####################

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  //G4double maxStep = 0.5*chamberWidth;	//###
  //fStepLimit = new G4UserLimits(maxStep);	//###
  //trackerLV->SetUserLimits(fStepLimit);	//###
 
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void nDaDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "nDa/TrackerChamberSD";
  nDaTrackerSD* aTrackerSD = new nDaTrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
//  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);		//##############################
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  //G4ThreeVector fieldValue = G4ThreeVector();
  //fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  //fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  //G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void nDaDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout 
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
G4String nDaDetectorConstruction::GetTargetMaterial() const
{
   if(fTargetMaterial) return fTargetMaterial->GetName();
   else return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaDetectorConstruction::SetTargetRadius(G4double radius)
{
   auto geoMgr = G4GeometryManager::GetInstance();
   auto targetS = (G4Sphere*)(fSolidTarget);
   //auto targetShellS = (G4Sphere*)(fSolidTargetShell);
   if(!fPhysTarget) { G4cerr << "Physical Volume for Target Not Found!" << G4endl; return; }
   if(!targetS) { G4cerr << "Solid for Target Not Found!" << G4endl; return; }
   //if(!fPhysTargetShell) { G4cerr << "Physical Volume for Shell Not Found!" << G4endl; return; }
   //if(!targetShellS) { G4cerr << "Solid for Shell Not Found!" << G4endl; return; }
   //if(fLockTarget) std::lock_guard<std::mutex> targetRadiusLockGuard(*fLockTarget);
   //if(fLockTargetShell) std::lock_guard<std::mutex> targetRadiusLockGuard(*fLockTargetShell);
   geoMgr->OpenGeometry(fPhysTarget);
   //geoMgr->OpenGeometry(fPhysTargetShell);
   targetS->SetOuterRadius(radius);
   //targetShellS->SetInnerRadius(depth);
   //geoMgr->CloseGeometry(true, false, fPhysTargetShell);
   geoMgr->CloseGeometry(true, false, fPhysTarget);
}

G4double nDaDetectorConstruction::GetTargetRadius() const
{
   if(!fPhysTarget) return 0.;
   auto targetS = (G4Sphere*)(fSolidTarget);
   return targetS->GetOuterRadius();
}


//...ooooOOOO0000....
void nDaDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                                               SetMaterial(fChamberMaterial);
        }
        G4cout 
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

G4String nDaDetectorConstruction::GetChamberMaterial() const
{
   return fChamberMaterial?fChamberMaterial->GetName():"";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}
G4double nDaDetectorConstruction::GetStepMax() const
{
  return 0.; //TODO: FIXME
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
