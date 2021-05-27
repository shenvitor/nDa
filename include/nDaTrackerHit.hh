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
/// \file nDaTrackerHit.hh
/// \brief Definition of the nDaTrackerHit class

#ifndef nDaTrackerHit_h
#define nDaTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class nDaTrackerHit : public G4VHit
{
  public:
    nDaTrackerHit();
    nDaTrackerHit(const nDaTrackerHit&);
    virtual ~nDaTrackerHit();

    // operators
    const nDaTrackerHit& operator=(const nDaTrackerHit&);
    G4bool operator==(const nDaTrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    //Can be continued
    void SetTrackID    (G4int track)      { fTrackID = track; };
    void SetChamberNb  (G4int chamb)      { fChamberNb = chamb; };
    void SetEdep       (G4double de)      { fEdep = de; };
    void SetEkin       (G4double ki)      { fEkin = ki; };
    void SetPos        (G4ThreeVector xyz){ fPos = xyz; };
    void SetIsBoundary (G4bool is)       { fIsBoundary = is; };
    void SetDirection  (G4ThreeVector xyz){ fDirection = xyz;};
    void SetMomentum   (G4ThreeVector xyz){ fMomentum = xyz; };
    void SetGlobalTime (G4double time)    { fGlobalTime = time; };

    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4int GetChamberNb() const   { return fChamberNb; };
    G4double GetEdep() const     { return fEdep; };
    G4double GetEkin() const     { return fEkin; }
    G4ThreeVector GetPos() const { return fPos; };
    G4bool GetIsBoundary() const { return fIsBoundary; };
    G4ThreeVector GetDirection() const {return fDirection; };
    G4ThreeVector GetMomentum() const {return fMomentum; };
    G4double GetGlobalTime() const {return fGlobalTime;};

  private:

      G4int         fTrackID;
      G4int         fChamberNb;
      G4double      fEdep;
      G4double      fEkin;
      G4ThreeVector fPos;
      G4bool        fIsBoundary;
      G4ThreeVector fDirection;
      G4ThreeVector fMomentum;
      G4double      fGlobalTime;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<nDaTrackerHit> nDaTrackerHitsCollection;

extern G4ThreadLocal G4Allocator<nDaTrackerHit>* nDaTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* nDaTrackerHit::operator new(size_t)
{
  if(!nDaTrackerHitAllocator)
      nDaTrackerHitAllocator = new G4Allocator<nDaTrackerHit>;
  return (void *) nDaTrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void nDaTrackerHit::operator delete(void *hit)
{
  nDaTrackerHitAllocator->FreeSingle((nDaTrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
