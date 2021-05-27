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
/// \file nDaDetectorMessenger.cc
/// \brief Implementation of the nDaDetectorMessenger class

#include "nDaDetectorMessenger.hh"
#include "nDaDetectorConstruction.hh"
#include "nDaPrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include <unordered_map>
/*class nDaDetectorMessenger::cmd_map_t :
  public std::unordered_map<G4UIcommand*, cmd_name_t> {};

G4UIcommand* nDaDetectorMessenger::GetCommand(cmd_name_t name) const
{
  return std::find_if(cmd_map->begin(), cmd_map->end(), [&](cmd_map_t::value_type& i){return i.second==name;})->first;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaDetectorMessenger::nDaDetectorMessenger(nDaDetectorConstruction* Det)
 : G4UImessenger(),
   fDetectorConstruction(Det)
{
  fnDaDirectory = new G4UIdirectory("/nDa/");
  fnDaDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/nDa/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  fTargMatCmd = new G4UIcmdWithAString("/nDa/det/setTargetMaterial",this);
  fTargMatCmd->SetGuidance("Select Material of the Target.");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //cmd_map->insert({fTargMatCmd, CMD_Target_Material});

  fChamMatCmd = new G4UIcmdWithAString("/nDa/det/setChamberMaterial",this);
  fChamMatCmd->SetGuidance("Select Material of the Chamber.");
  fChamMatCmd->SetParameterName("choice",false);
  fChamMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //cmd_map->insert({fChamMatCmd, CMD_Chamber_Material});

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/nDa/det/stepMax",this);
  fStepMaxCmd->SetGuidance("Define a step max");
  fStepMaxCmd->SetParameterName("stepMax",false);
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_Idle);

  fTargRadCmd = new G4UIcmdWithADoubleAndUnit("/nDa/det/setTargetRadius",this);
  fTargRadCmd->SetGuidance("Set the radius of the Target Ball. [radius] [unit]");
  fTargRadCmd->SetParameterName("targetRadius", true);
  fStepMaxCmd->SetUnitCategory("Length");
  fTargRadCmd->SetDefaultValue(1.);
  fTargRadCmd->SetDefaultUnit("cm");
  fTargRadCmd->SetUnitCandidates("um mm cm dm m");
  fTargRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //cmd_map->insert({fTargRadCmd, CMD_Target_Radius});
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaDetectorMessenger::~nDaDetectorMessenger()
{
  delete fTargMatCmd;
  delete fChamMatCmd;
  delete fStepMaxCmd;
  delete fTargRadCmd;
  delete fnDaDirectory;
  delete fDetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fTargMatCmd )
   { fDetectorConstruction->SetTargetMaterial(newValue);}

  if( command == fChamMatCmd )
   { fDetectorConstruction->SetChamberMaterial(newValue);}

  if( command == fStepMaxCmd ) {
    fDetectorConstruction
      ->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));
  }   
  if( command == fTargRadCmd )
   { fDetectorConstruction->SetTargetRadius(fTargRadCmd->GetNewDoubleValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void nDaDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  auto cmd_name = (*cmd_map)[command];
  switch(cmd_name)
  {
    case CMD_Target_Material:
      fDetectorConstruction->SetTargetMaterial(newValue);
      break;

    case CMD_Chamber_Material:
      fDetectorConstruction->SetChamberMaterial(newValue);
      break;
    case CMD_Step_Max:
    case CMD_Target_Radius:
      {
        auto newRawValue = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue);
        switch(cmd_name)
        {
          case CMD_Step_Max:
            fDetectorConstruction->SetMaxStep(newRawValue);
            break;
          case CMD_Target_Radius:
            fDetectorConstruction->SetTargetRadius(newRawValue);
            break;
          default: break;
        }
      } break;
  
    default: break;
  }
}*/
//G4String nDaDetectorMessenger::GetCurrentValue(G4UIcommand* command)
/*{
  auto cmd_name = (*cmd_map)[command];
  switch (cmd_name)
  {
  case CMD_Target_Radius:
    return static_cast<G4UIcmdWithADoubleAndUnit*>(command)->ConvertToStringWithBestUnit(fDetectorConstruction->GetTargetRadius());
  case CMD_Chamber_Material:
    return fDetectorConstruction->GetChamberMaterial();
  case CMD_Target_Material:
    return fDetectorConstruction->GetTargetMaterial();
  case CMD_Step_Max:
    return static_cast<G4UIcmdWithADoubleAndUnit*>(command)->ConvertToStringWithBestUnit(fDetectorConstruction->GetStepMax());
  default:
    return "";
  }
}*/