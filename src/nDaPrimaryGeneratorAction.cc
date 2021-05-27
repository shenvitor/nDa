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
/// \file nDaPrimaryGeneratorAction.cc
/// \brief Implementation of the nDaPrimaryGeneratorAction class

#include "nDaPrimaryGeneratorAction.hh"
#include "nDaAnalysisManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include "g4root.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//double GunEnergy =  5.*MeV;
//double GunEnergy = 10.*MeV;
//double GunEnergy = 14.*MeV;
//double GunEnergy = 15.*MeV;
//double GunEnergy = 26.*MeV;
//double GunEnergy = 50.*MeV;
//double GunEnergy = 75.*MeV;
//double GunEnergy = 100.*MeV;

nDaPrimaryGeneratorAction::nDaPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction()
{
  //G4int nofParticles = 1;
  //fParticleGun = new G4GeneralParticleSource(nofParticles);
  fGPS = new G4GeneralParticleSource();
  // default particle kinematic

  G4ParticleDefinition* particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("neutron");

  fGPS->SetParticleDefinition(particleDefinition);
  //fGPS->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //fGPS->SetParticleEnergy(GunEnergy);
  //fParticleGun->SetParticleEnergy(100.*MeV);
  //fParticleGun->SetParticleEnergy(50.*MeV);
  //fParticleGun->SetParticleEnergy(10. *MeV);
  //fParticleGun->SetParticleEnergy(14.*MeV);

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4double worldZHalfLength = 0;
  G4LogicalVolume* worldLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = NULL;
  if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if ( worldBox ) worldZHalfLength = worldBox->GetZHalfLength();
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  // Note that this particular case of starting a primary particle on the world boundary
  // requires shooting in a direction towards inside the world.
  //fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));
  
  //Setting to be in a distance of 5cm+10m = 1005cm
  fGPS->SetParticlePosition(G4ThreeVector(0., 0., -1005.*cm));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaPrimaryGeneratorAction::~nDaPrimaryGeneratorAction()
{
  delete fGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void nDaPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

 
  fParticleGun->GeneratePrimaryVertex(anEvent);
}*/

void nDaPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  fGPS->GeneratePrimaryVertex(anEvent);
  // Record the GPS status
  auto anaMgr = G4AnalysisManager::Instance();
  anaMgr->FillNtupleIColumn(0, anEvent->GetEventID());
  anaMgr->FillNtupleDColumn(1, fGPS->GetParticleEnergy());
  auto direction = fGPS->GetParticleMomentumDirection();
  anaMgr->FillNtupleDColumn(2, direction.x());
  anaMgr->FillNtupleDColumn(3, direction.y());
  anaMgr->FillNtupleDColumn(4, direction.z());
  anaMgr->FillNtupleDColumn(5, direction.r());
  anaMgr->FillNtupleDColumn(6, direction.theta());
  anaMgr->FillNtupleDColumn(7, direction.phi());
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// how to get info from particle gun?
//auto anaMgr=G4AnalysisManager::Instance();

