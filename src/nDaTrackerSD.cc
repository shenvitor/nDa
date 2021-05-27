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
/// \file nDaTrackerSD.cc
/// \brief Implementation of the nDaTrackerSD class

#include "nDaTrackerSD.hh"
#include "nDaDetectorConstruction.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaTrackerSD::nDaTrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaTrackerSD::~nDaTrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new nDaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool nDaTrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory* aTH)
{  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  //G4double eesc = aStep->GetDeltaEnergy(); 
  //G4double eesc = aStep->GetPostStepPoint()->GetTotalEnergy();

  if (edep==0.) return false;

  nDaTrackerHit* newHit = new nDaTrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
                                               ->GetCopyNumber());
  newHit->SetEdep(edep);
  //newHit->SetEesc(aStep->GetPostStepPoint()->GetTotalEnergy());
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetDirection(aStep->GetPostStepPoint()->GetMomentumDirection());
  newHit->SetIsBoundary(aStep->IsLastStepInVolume());
  newHit->SetGlobalTime(aStep->GetPostStepPoint()->GetGlobalTime());

  // Beware of how to get KE and momentum!
  // G4Track* track = step -> GetTrack();
  // To get Etot, Ek, momentum, should it be track or point!?
  //newHit->SetEkin(edep); 
  newHit->SetEkin(aStep->GetPostStepPoint()->GetKineticEnergy()); 
  newHit->SetMomentum(aStep->GetPostStepPoint()->GetMomentum());  

  //G4Track* aTrack = aStep->GetTrack();
  //newHit->SetEkin(aTrack->GetKineticEnergy());
  //newHit->SetMomentum(aTrack->GetMomentum());

  auto runMgr = G4RunManager::GetRunManager();
  nDaDetectorConstruction* det = (nDaDetectorConstruction*)(runMgr->GetUserDetectorConstruction());
  det->GetNbOfRows();

  fHitsCollection->insert( newHit );

  //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
