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
/// \file nDaDetectorConstruction.hh
/// \brief Definition of the nDaDetectorConstruction class

#ifndef nDaDetectorConstruction_h
#define nDaDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class nDaDetectorMessenger;

namespace std { class mutex; }

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class nDaDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    nDaDetectorConstruction();
    virtual ~nDaDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // Set methods
    void SetTargetMaterial (G4String );
    void SetChamberMaterial(G4String );

    //void SetTargetShellMaterial (G4String ); //!
    //void SetTargetDepth (G4double);          //!
    //void SetTargetShellDepth(G4double);      //!

    void SetTargetRadius (G4double);
    void SetMaxStep (G4double );
    void SetCheckOverlaps(G4bool );
    void SetTrackerPos (G4ThreeVector);
    void SetTrackerPosX (G4double);
    void SetTrackerPosY (G4double);
    void SetTrackerPosZ (G4double);

    // Get methods                           //!
    //G4double GetTargetDepth() const;
    //G4double GetTargetShellDepth() const;
    G4String GetTargetMaterial() const;
    //G4String GetTargetShellMaterial() const;
    G4String GetChamberMaterial() const;
    G4double GetTargetRadius() const;
    G4double GetStepMax() const;
    G4ThreeVector GetTrackerPos() const;
    //G4double GetTrackerPosY() const;
    //G4double GetTrackerPosZ() const;

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
  
    // data members
    G4int fNbOfChambers;
    G4int fNbOfRows;
    G4int fNbOfLevels;

    G4double trackerPosX;
    G4double trackerPosY;
    G4double trackerPosZ;
    G4ThreeVector     fTracPos;        //the vector of tracker position

    G4VPhysicalVolume* fPhysTarget;   //!// pointer to the physical Target
    G4LogicalVolume*   fLogicTarget;     // pointer to the logical Target
    G4VSolid*          fSolidTarget;
    std::mutex*        fLockTarget;

    //G4VPhysicalVolume* fPhysTargetShell;
    //G4LogicalVolume*   fLogicTargetShell;
    //G4VSolid*          fSolidTargetShell;   

    G4VPhysicalVolume** fPhysChamber;   //!// pointer to the physical Target
    G4LogicalVolume**  fLogicChamber;    // pointer to the logical Chamber

    G4LogicalVolume*   fLogicTracker;
    G4VPhysicalVolume* fPhysTracker;  //!// pointer to the physical Tracker  

    G4Material*        fTargetMaterial;  // pointer to the target  material
    //G4Material*        fTargetShellMaterial;   //!
    G4Material*        fChamberMaterial; // pointer to the chamber material
    G4Material*        fTrackerMaterial; //pointer to the tracker material

    G4UserLimits* fStepLimit;            // pointer to user step limits

    nDaDetectorMessenger*  fMessenger;   // messenger

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                         // magnetic field messenger
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 

  public:
    inline G4int GetNbOfChambers() const noexcept { return fNbOfChambers; }
    inline G4int GetNbOfRows() const noexcept { return fNbOfRows; }
    inline G4int GetNbOfLevels() const noexcept { return fNbOfLevels; }
    static std::array<G4int,3> ResolveCopyNb(G4int new_copyNb)
    {
      // calculate row number
      // return as: {nbInRow, nbRow, nbLevel}
      auto runMgr = G4RunManager::GetRunManager();
      const auto det = (nDaDetectorConstruction*)(runMgr->GetUserDetectorConstruction());
      G4int fNbOfRows = det->GetNbOfRows(), fNbOfChambers = det->GetNbOfChambers(); //fNbOfLevels = det->GetNbOfLevels();
      G4int nbLevel = new_copyNb/(fNbOfRows*fNbOfChambers);
      G4int nbInLevel = new_copyNb-nbLevel*fNbOfRows*fNbOfChambers;
      G4int nbRow = nbInLevel/fNbOfChambers;
      G4int nbInRow = nbInLevel%fNbOfChambers;
      return {nbInRow, nbRow, nbLevel};
    }
    static G4ThreeVector GetChamberPosition(G4int copyNb)
    {
      // calculate row number
      // return as: {nbInRow, nbRow, nbLevel}
      auto runMgr = G4RunManager::GetRunManager();
      const auto det = (nDaDetectorConstruction*)(runMgr->GetUserDetectorConstruction());
      auto rela_pos = det->fPhysChamber[copyNb]->GetTranslation();
      auto mother_pos = det->fPhysTracker->GetTranslation();
      return rela_pos + mother_pos;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
