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
/// \file nDaRunAction.cc
/// \brief Implementation of the nDaRunAction class

#include "nDaRunAction.hh"
#include "nDaDetectorConstruction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

#include <g4root.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaRunAction::nDaRunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);  
  // Create analysis manager
	auto anaMgr = G4AnalysisManager::Instance();
	anaMgr->SetFileName("neutron");
	anaMgr->SetVerboseLevel(1);
	anaMgr->SetFirstHistoId(1);
	anaMgr->CreateH1("nEdet", "the Neutron Energy that Detected [MeV]", 100, 0., 120.);		//1
	//anaMgr->CreateH1("Edep", "Total Energy Deposit of Neutron [MeV]", 100, 0., 120.);		//2
	// Ntuple for data collection
  anaMgr->SetNtupleMerging(true);
	anaMgr->CreateNtuple("neutron", "Detected Neutron");
	// simulation data of GPS status
	anaMgr->CreateNtupleIColumn("EventID");				      	  //0
	anaMgr->CreateNtupleDColumn("IncidentEnergy");	  		  //1
	anaMgr->CreateNtupleDColumn("IncidentDirectionX");	  	//2
	anaMgr->CreateNtupleDColumn("IncidentDirectionY");		  //3
	anaMgr->CreateNtupleDColumn("IncidentDirectionZ");	  	//4
	anaMgr->CreateNtupleDColumn("IncidentDirectionR");	  	//5
	anaMgr->CreateNtupleDColumn("IncidentDirectionTheta");	//6
	anaMgr->CreateNtupleDColumn("IncidentDirectionPhi");	  //7
	// simulation data of neutron Energy detected
	anaMgr->CreateNtupleDColumn("nEdet");					          //8
	anaMgr->CreateNtupleDColumn("nEdet_PosX");		       	  //9
	anaMgr->CreateNtupleDColumn("nEdet_PosY");				      //10
	anaMgr->CreateNtupleDColumn("nEdet_PosZ");			       	//11
	// simulation data of detector & detector hit
	anaMgr->CreateNtupleDColumn("Hit_Edep", Edep);				
	anaMgr->CreateNtupleDColumn("Hit_PosX", PosX);			
	anaMgr->CreateNtupleDColumn("Hit_PosY", PosY);				
	anaMgr->CreateNtupleDColumn("Hit_PosZ", PosZ);				
	anaMgr->CreateNtupleIColumn("Hit_copyNb", copyNb);				
	anaMgr->CreateNtupleIColumn("Hit_InRowNb", InRowNb);
	anaMgr->CreateNtupleIColumn("Hit_RowNb", RowNb);
	anaMgr->CreateNtupleIColumn("Hit_LevelNb", LevelNb);
	anaMgr->CreateNtupleDColumn("Hit_GlobalTime", GlobalTime);
	anaMgr->CreateNtupleIColumn("Detector_copyNb", det_copyNb);
	anaMgr->CreateNtupleIColumn("Detector_InRowNb", det_InRowNb);
	anaMgr->CreateNtupleIColumn("Detector_Level_Nb", det_LevelNb);
	anaMgr->CreateNtupleIColumn("Detector_Row_Nb", det_RowNb);
	anaMgr->CreateNtupleDColumn("Detector_Edep", det_Edep);
	anaMgr->CreateNtupleDColumn("Detector_GlobalTime", det_GlobalTime);

	//anaMgr->CreateNtupleIColumn("Escaped");//5
	//anaMgr->CreateNtupleDColumn("EscapedTheta");//6
	//anaMgr->CreateNtupleDColumn("EscapedPhi");//7
	//anaMgr->CreateNtupleDColumn("EscapedPosX");//8
	//anaMgr->CreateNtupleDColumn("EscapedPosY");//9
	//anaMgr->CreateNtupleDColumn("EscapedPosZ");//10
	//anaMgr->CreateNtupleDColumn("Eesc");//11
	//anaMgr->CreateNtupleDColumn("EscMomentumX");//12
    //anaMgr->CreateNtupleDColumn("EscMomentumY");//13
	//anaMgr->CreateNtupleDColumn("EscMomentumZ");//14

	anaMgr->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDaRunAction::~nDaRunAction()
{
	Clear_Hit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4long nDaRunAction::Add_Hit(G4double new_Edep, G4ThreeVector new_Pos, G4int new_copyNb, G4double new_GlobalTime)
{
  Edep.push_back(new_Edep);
  PosX.push_back(new_Pos.x());
  PosY.push_back(new_Pos.y());
  PosZ.push_back(new_Pos.z());
  copyNb.push_back(new_copyNb);
  auto resolved_copyNb = nDaDetectorConstruction::ResolveCopyNb(new_copyNb);
  InRowNb.push_back(resolved_copyNb[0]);
  RowNb.push_back(resolved_copyNb[1]);
  LevelNb.push_back(resolved_copyNb[2]);
  GlobalTime.push_back(new_GlobalTime);
  return Edep.size();
}

void nDaRunAction::Clear_Hit() noexcept
{
  Edep.clear();
  PosX.clear();
  PosY.clear();
  PosZ.clear();
  copyNb.clear();
  InRowNb.clear();
  LevelNb.clear();
  RowNb.clear();
  GlobalTime.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaRunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  
  // Initialize Analysis Manager
  auto anaMgr = G4AnalysisManager::Instance();
  anaMgr->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nDaRunAction::EndOfRunAction(const G4Run* )
{
	auto anaMgr = G4AnalysisManager::Instance();
	anaMgr->Write();
	anaMgr->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Calculate all detector branches
G4bool nDaRunAction::Calculate_Detector()
    {
      constexpr G4double THRESHOLD = 1 * CLHEP::MeV;
      G4double nEdet=0.;
      G4ThreeVector nEdet_Pos;

      // Initialize Detectors
      auto runMgr = G4RunManager::GetRunManager();
      auto det = (nDaDetectorConstruction*)(runMgr->GetUserDetectorConstruction());
      G4int fNbOfRows = det->GetNbOfRows(), fNbOfChambers = det->GetNbOfChambers(), fNbOfLevels = det->GetNbOfLevels();
      G4int NCopyNbs = fNbOfRows * fNbOfChambers * fNbOfLevels;
      det_copyNb.resize(NCopyNbs);
      det_InRowNb.resize(NCopyNbs);
      det_LevelNb.resize(NCopyNbs);
      det_RowNb.resize(NCopyNbs);
      det_Edep.resize(NCopyNbs);
      det_GlobalTime.resize(NCopyNbs);

      // populate row and level id
      for(G4int this_copyNb=0; this_copyNb<NCopyNbs; ++this_copyNb)
      {
        auto resolved_copyNb = nDaDetectorConstruction::ResolveCopyNb(this_copyNb);
        det_copyNb[this_copyNb]=this_copyNb;
        det_InRowNb[this_copyNb]=resolved_copyNb[0];
        det_RowNb[this_copyNb]=resolved_copyNb[1];
        det_LevelNb[this_copyNb]=resolved_copyNb[2];
        det_Edep[this_copyNb]=0.;
        det_GlobalTime[this_copyNb]=-1.;
      }

      // Calculate Detectors
      for(size_t i=0; i<Edep.size(); ++i)
      {
        auto this_copyNb = copyNb[i];
        det_Edep[this_copyNb]+=Edep[i];
        if(det_Edep[this_copyNb]>=THRESHOLD && det_GlobalTime[i]<0) det_GlobalTime[this_copyNb] = GlobalTime[i];
      }

      // Calculate nEdet and Center of Mass
      for(G4int this_copyNb=0; this_copyNb<NCopyNbs; ++this_copyNb)
      {
        auto this_nEdet = det_Edep[this_copyNb];
        auto position = det->GetChamberPosition(this_copyNb);
        nEdet+=this_nEdet;
        nEdet_Pos+=this_nEdet*position;
      }
      nEdet_Pos/=nEdet;

      auto anaMgr=G4AnalysisManager::Instance();
      anaMgr->FillH1(1, nEdet);
      //anaMgr->FillH1(2, Edep);
      anaMgr->FillNtupleDColumn(8, nEdet);
      anaMgr->FillNtupleDColumn(9, nEdet_Pos.x());
      anaMgr->FillNtupleDColumn(10, nEdet_Pos.y());
      anaMgr->FillNtupleDColumn(11, nEdet_Pos.z());
      return anaMgr->AddNtupleRow();
    }