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
/// \file nDaRunAction.hh
/// \brief Definition of the nDaRunAction class

#ifndef nDaRunAction_h
#define nDaRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <vector>
#include "G4ThreeVector.hh"

#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

/// Run action class

class nDaRunAction : public G4UserRunAction
{
  public:
    nDaRunAction();
    virtual ~nDaRunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void   EndOfRunAction(const G4Run* run);

private:
    std::vector<G4double> Edep;
    std::vector<G4double> PosX, PosY, PosZ;
    std::vector<G4int> copyNb, InRowNb, LevelNb, RowNb;
    std::vector<G4double> GlobalTime;
    std::vector<G4int> det_copyNb, det_InRowNb, det_LevelNb, det_RowNb;
    std::vector<G4double> det_Edep, det_GlobalTime;

  public:
    G4long Add_Hit(G4double new_Edep, G4ThreeVector new_Pos, G4int new_copyNb, G4double new_GlobalTime);
    void Clear_Hit() noexcept;
    G4bool Calculate_Detector(); // Calculate all Detector branches

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
