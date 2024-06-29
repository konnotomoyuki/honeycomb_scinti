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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  nistManager->FindOrBuildMaterial("G4_Al");

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double calorSizeXY  = 4 * cm;
  G4double calorThickness  = 2. * cm;
  G4double detSize = 40 * cm;
  /*
  G4double calorSizeXY2 = 1. * cm;
  G4double calorThickness2 = 1. * cm;
  */
  
  auto worldSizeXY = 10 * m;
  auto worldSizeZ  = 10 * m;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto calorMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto boardMaterial = G4Material::GetMaterial("G4_Al");

  if ( ! defaultMaterial || ! calorMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Calorimeter
  //
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
  auto vetoS
    = new G4Box("veto",     // its name
                 50 * cm / 2., 50 * cm / 2., 20 * cm / 2.); // its size

  auto veto2S
    = new G4Box("veto2",     // its name
                 20 * cm / 2., 50 * cm / 2., 50 * cm / 2.); // its size

  auto veto3S
    = new G4Box("veto3",     // its name
                 50 * cm / 2., 20 * cm / 2., 50 * cm / 2.); // its size
  auto boardS
    = new G4Box("board",     // its name
                 40 * cm / 2., 40 * cm / 2., 2 * mm / 2.); // its size

  auto board2S
    = new G4Box("board2",     // its name
                 20 * mm / 2., 40 * cm / 2., 40 * cm / 2.); // its size

  auto board3S
    = new G4Box("board3",     // its name
                 40 * cm / 2., 2 * mm / 2., 40 * cm / 2.); // its size

  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 calorMaterial,  // its material
                 "Calorimeter");   // its name
  auto vetoLV
    = new G4LogicalVolume(
                 vetoS,     // its solid
                 calorMaterial,  // its material
                 "Veto");   // its name
  auto veto2LV
    = new G4LogicalVolume(
                 veto2S,     // its solid
                 calorMaterial,  // its material
                 "Veto2");   // its name
  
  auto veto3LV
    = new G4LogicalVolume(
                 veto3S,     // its solid
                 calorMaterial,  // its material
                 "Veto3");   // its name
  
  auto boardLV
    = new G4LogicalVolume(
                 boardS,     // its solid
                 boardMaterial,  // its material
                 "Board");   // its name
  auto board2LV
    = new G4LogicalVolume(
                 board2S,     // its solid
                 boardMaterial,  // its material
                 "Board2");   // its name
  
  auto board3LV
    = new G4LogicalVolume(
                 board3S,     // its solid
                 boardMaterial,  // its material
                 "Board3");   // its name
  
  G4Colour calor_color = G4Colour::Magenta();
  G4VisAttributes *calor_attributes = new G4VisAttributes(calor_color);
  calor_attributes->SetVisibility(true);
  calor_attributes->SetForceWireframe(true);
  calorLV->SetVisAttributes(calor_attributes);
  vetoLV->SetVisAttributes(calor_attributes);
  veto2LV->SetVisAttributes(calor_attributes);
  veto3LV->SetVisAttributes(calor_attributes);
  
  G4Colour veto_color(0., 1., 0., 0.5);    // red, green, blue, alpha
  G4VisAttributes *veto_attributes = new G4VisAttributes(veto_color);
  veto_attributes->SetVisibility(true);
  veto_attributes->SetForceWireframe(true);
  boardLV->SetVisAttributes(veto_attributes);
  board2LV->SetVisAttributes(veto_attributes);
  board3LV->SetVisAttributes(veto_attributes);
  
  int nx = (int)(detSize / calorSizeXY);
  int ny = (int)(detSize / calorSizeXY);
  int nz = (int)(detSize / calorThickness);
  double xmax = -100 * m;
  double xmin = 100 * m;
  double ymax = -100 * m;
  double ymin = 100 * m;
  double zmax = -100 * m;
  double zmin = 100 * m;
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
	double x = calorSizeXY * (ix);
	double y = calorSizeXY * (iy);
	double z = calorThickness * (iz);
	if (x < xmin) xmin = x;
	if (x > xmax) xmax = x;
	if (y < ymin) ymin = y;
	if (y > ymax) ymax = y;
	if (z < zmin) zmin = z;
	if (z > zmax) zmax = z;
	new G4PVPlacement(
			  0,                // no rotation
			  G4ThreeVector(x, y, z),  // at (0,0,0)
			  calorLV,          // its logical volume
			  "Calorimeter",    // its name
			  worldLV,          // its mother  volume
			  false,            // no boolean operation
			  (ix + 1) + iy * 100 + iz * 10000,                // copy number
			  fCheckOverlaps);  // checking overlaps
	if (true 
	    /*iz - nz/2 == 0 &&
	    iy - ny/2 == 0 &&
	    ix - nx/2 == 0*/) {
	  std::cout << "detector center " << x << " " << y << " " << z << std::endl;
	}
      }
    }
  }
  double xmean = (xmin + xmax) / 2.0;
  double ymean = (ymin + ymax) / 2.0;
  double zmean = (zmin + zmax) / 2.0;
  //std::cout << "detector center " << xmean << " " << ymean << " " << zmean << std::endl;
  int i = 1;
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean,
				  zmean + (detSize/2 + 15 * cm)),  // at (0,0,0)
		    vetoLV,          // its logical volume
			  "Calorimeter",    // its name
			  worldLV,          // its mother  volume
			  false,            // no boolean operation
			  -i,                // copy number
			  fCheckOverlaps);  // checking overlaps
  i++;
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean,
				  zmean - (detSize/2 + 15 * cm)),  // at (0,0,0)
		    vetoLV,          // its logical volume
			  "Calorimeter",    // its name
			  worldLV,          // its mother  volume
			  false,            // no boolean operation
			  -i,                // copy number
			  fCheckOverlaps);  // checking overlaps
  i++;
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean + (detSize/2 + 15 * cm),
				  ymean,
				  zmean),  // at (0,0,0)
		    veto2LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    -i,                // copy number
		    fCheckOverlaps);  // checking overlaps
  i++;
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean - (detSize/2 + 15 * cm),
				  ymean,
				  zmean),  // at (0,0,0)
		    veto2LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    -i,                // copy number
		    fCheckOverlaps);  // checking overlaps
  i++;
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean + (detSize/2 + 15 * cm),
				  zmean),  // at (0,0,0)
		    veto3LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    -i,                // copy number
		    fCheckOverlaps);  // checking overlaps
  i++;
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean - (detSize/2 + 15 * cm),
				  zmean),  // at (0,0,0)
		    veto3LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    -i,                // copy number
		    fCheckOverlaps);  // checking overlaps
  i++;

  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean,
				  zmean + (detSize/2 + 1 * cm)),  // at (0,0,0)
		    boardLV,          // its logical volume
			  "Calorimeter",    // its name
			  worldLV,          // its mother  volume
			  false,            // no boolean operation
			  0,                // copy number
			  fCheckOverlaps);  // checking overlaps
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean,
				  zmean - (detSize/2 + 1 * cm)),  // at (0,0,0)
		    boardLV,          // its logical volume
			  "Calorimeter",    // its name
			  worldLV,          // its mother  volume
			  false,            // no boolean operation
			  0,                // copy number
			  fCheckOverlaps);  // checking overlaps
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean + (detSize/2 + 1 * cm),
				  ymean,
				  zmean),  // at (0,0,0)
		    board2LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean - (detSize/2 + 1 * cm),
				  ymean,
				  zmean),  // at (0,0,0)
		    board2LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean + (detSize/2 + 1 * cm),
				  zmean),  // at (0,0,0)
		    board3LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps
  new G4PVPlacement(
		    0,                // no rotation
		    G4ThreeVector(xmean,
				  ymean - (detSize/2 + 1 * cm),
				  zmean),  // at (0,0,0)
		    board3LV,          // its logical volume
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps
  
  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

