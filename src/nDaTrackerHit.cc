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
/// \file nDaTrackerHit.cc
/// \brief Implementation of the nDaTrackerHit class

#include "nDaTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<nDaTrackerHit>* nDaTrackerHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaTrackerHit::nDaTrackerHit()
 : G4VHit(),
   fTrackID(-1),
   fChamberNb(-1),
   fEdep(0.),
   fEkin(0.),
   fPos(G4ThreeVector()),
   fIsBoundary(false),
   fDirection(0,0,0),
   fMomentum(G4ThreeVector())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaTrackerHit::~nDaTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaTrackerHit::nDaTrackerHit(const nDaTrackerHit& right)
  : G4VHit()
{
  //(*this)=right;
  this->operator=(right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const nDaTrackerHit& nDaTrackerHit::operator=(const nDaTrackerHit& right)
{
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fEkin      = right.fEkin;
  fPos       = right.fPos;
  fIsBoundary= right.fIsBoundary;
  fDirection = right.fDirection;
  fMomentum  = right.fMomentum;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool nDaTrackerHit::operator==(const nDaTrackerHit& right) const
{
  return ( this == &right ) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaTrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaTrackerHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID << " chamberNb: " << fChamberNb
     << "Edep: "
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " Position: "
     << std::setw(7) << G4BestUnit( fPos,"Length")
     << "Direction: "
     << std::setw(7) << G4BestUnit(fDirection, "Degree")
     << "Eese: "
     << std::setw(7) << G4BestUnit(fEkin,"Energy")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
