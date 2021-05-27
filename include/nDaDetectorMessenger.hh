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
/// \file nDaDetectorMessenger.hh
/// \brief Definition of the nDaDetectorMessenger class

#ifndef nDaDetectorMessenger_h
#define nDaDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class nDaDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for nDaDetectorConstruction.
///
/// It implements commands:
/// - /nDa/det/setTargetMaterial name
/// - /nDa/det/setChamberMaterial name
/// - /nDa/det/stepMax value unit

class nDaDetectorMessenger: public G4UImessenger
{
  public:
    nDaDetectorMessenger(nDaDetectorConstruction* );
    virtual ~nDaDetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    //virtual G4String GetCurrentValue(G4UIcommand* command);
    
  private:
    nDaDetectorConstruction*  fDetectorConstruction;

    G4UIdirectory*           fnDaDirectory;
    G4UIdirectory*           fDetDirectory;

    G4UIcmdWithAString*      fTargMatCmd;
    G4UIcmdWithAString*      fChamMatCmd;

    G4UIcmdWithADoubleAndUnit* fStepMaxCmd;
    G4UIcmdWithADoubleAndUnit* fTargRadCmd;

    G4UIcmdWithADoubleAndUnit* fTracPosXCmd;
    G4UIcmdWithADoubleAndUnit* fTracPosYCmd;
    G4UIcmdWithADoubleAndUnit* fTracPosZCmd;

    G4UIcmdWithAString* fTracPosCmd;

    /*enum cmd_name_t {
      CMD_OTHER=0,
      CMD_Target_Material,
      //CMD_Target_Shell_Material,
      CMD_Target_Radius,
      //CMD_Target_Depth,
      //CMD_Target_Shell_Depth,
      CMD_Chamber_Material,
      CMD_Step_Max //,
      //CMD_Mode_Change
    };
    class cmd_map_t;
    cmd_map_t* cmd_map;
    G4UIcommand* GetCommand(cmd_name_t name) const;
    inline G4String GetCurrentValue(cmd_name_t name) { return GetCurrentValue(GetCommand(name)); }*/

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
