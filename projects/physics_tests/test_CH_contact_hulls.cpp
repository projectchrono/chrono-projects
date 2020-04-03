// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Test program for collision with contact hulls
//
// =============================================================================

#include "chrono/ChConfig.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/collision/ChCollisionUtils.h"

#include "chrono_irrlicht/ChIrrApp.h"

using namespace chrono;
using namespace chrono::irrlicht;

void AddWallBox(std::shared_ptr<ChBody> body, std::shared_ptr<ChMaterialSurface> mat, const ChVector<>& dim, const ChVector<>& loc) {
    body->GetCollisionModel()->AddBox(mat, dim.x(), dim.y(), dim.z(), loc);

    auto box = chrono_types::make_shared<ChBoxShape>();
    box->GetBoxGeometry().Size = dim;
    box->GetBoxGeometry().Pos = loc;
    body->AddAsset(box);

    body->AddAsset(chrono_types::make_shared<ChColorAsset>(0.5f, 0.5f, 0.0f));
}

void AddWallMesh(std::shared_ptr<ChBody> body,
                 std::shared_ptr<ChMaterialSurface> mat,
                 const ChVector<>& dim,
                 const ChVector<>& loc) {
    auto trimesh = chrono_types::make_shared<geometry::ChTriangleMeshConnected>();

    std::vector<ChVector<> >& vertices = trimesh->getCoordsVertices();
    std::vector<ChVector<> >& normals = trimesh->getCoordsNormals();
    std::vector<ChVector<int> >& idx_vertices = trimesh->getIndicesVertexes();
    std::vector<ChVector<int> >& idx_normals = trimesh->getIndicesNormals();

    int num_vert = 8;
    int num_faces = 12;

    vertices.resize(num_vert);
    normals.resize(num_vert);
    trimesh->getCoordsUV().resize(num_vert);
    trimesh->getCoordsColors().resize(num_vert);

    idx_vertices.resize(num_faces);
    idx_normals.resize(num_faces);

    vertices[0] = ChVector<>(-dim.x(), -dim.y(), -dim.z()) + loc;
    vertices[1] = ChVector<>(-dim.x(), +dim.y(), -dim.z()) + loc;
    vertices[2] = ChVector<>(+dim.x(), +dim.y(), -dim.z()) + loc;
    vertices[3] = ChVector<>(+dim.x(), -dim.y(), -dim.z()) + loc;
    
    normals[0] = ChVector<>(-1, -1, -1).Normalize();
    normals[1] = ChVector<>(-1, +1, -1).Normalize();
    normals[2] = ChVector<>(+1, +1, -1).Normalize();
    normals[3] = ChVector<>(+1, -1, -1).Normalize();

    vertices[4] = ChVector<>(-dim.x(), -dim.y(), +dim.z()) + loc;
    vertices[5] = ChVector<>(-dim.x(), +dim.y(), +dim.z()) + loc;
    vertices[6] = ChVector<>(+dim.x(), +dim.y(), +dim.z()) + loc;
    vertices[7] = ChVector<>(+dim.x(), -dim.y(), +dim.z()) + loc;

    normals[4] = ChVector<>(-1, -1, +1).Normalize();
    normals[5] = ChVector<>(-1, +1, +1).Normalize();
    normals[6] = ChVector<>(+1, +1, +1).Normalize();
    normals[7] = ChVector<>(+1, -1, +1).Normalize();

    for (int i = 0; i < num_vert; i++) {
        trimesh->getCoordsColors()[i] = ChVector<float>(0, 0, 1);
    }

    idx_vertices[0] = ChVector<int>(0, 1, 3);
    idx_vertices[1] = ChVector<int>(1, 2, 3);
    idx_vertices[2] = ChVector<int>(4, 7, 5);
    idx_vertices[3] = ChVector<int>(7, 6, 5);

    idx_vertices[4] = ChVector<int>(0, 3, 4);
    idx_vertices[5] = ChVector<int>(3, 7, 4);
    idx_vertices[6] = ChVector<int>(1, 5, 2);
    idx_vertices[7] = ChVector<int>(2, 5, 6);

    idx_vertices[8] = ChVector<int>(3, 6, 7);
    idx_vertices[9] = ChVector<int>(3, 2, 6);
    idx_vertices[10] = ChVector<int>(0, 4, 5);
    idx_vertices[11] = ChVector<int>(0, 5, 1);

    for (int i = 0; i < num_faces; i++)
      idx_normals[i] = idx_vertices[i];

    body->GetCollisionModel()->AddTriangleMesh(mat, trimesh, true, true, loc);

    auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    trimesh_shape->SetMesh(trimesh);
    body->AddAsset(trimesh_shape);

    body->AddAsset(chrono_types::make_shared<ChColorAsset>(0.0f, 0.0f, 0.5f));
}

void AddWallHull(std::shared_ptr<ChBody> body,
                 std::shared_ptr<ChMaterialSurface> mat,
                 const ChVector<>& dim,
                 const ChVector<>& loc) {
    std::vector<ChVector<>> points;

    points.push_back(ChVector<>(-dim.x(), -dim.y(), -dim.z()) + loc);
    points.push_back(ChVector<>(-dim.x(), +dim.y(), -dim.z()) + loc);
    points.push_back(ChVector<>(+dim.x(), +dim.y(), -dim.z()) + loc);
    points.push_back(ChVector<>(+dim.x(), -dim.y(), -dim.z()) + loc);
    points.push_back(ChVector<>(-dim.x(), -dim.y(), +dim.z()) + loc);
    points.push_back(ChVector<>(-dim.x(), +dim.y(), +dim.z()) + loc);
    points.push_back(ChVector<>(+dim.x(), +dim.y(), +dim.z()) + loc);
    points.push_back(ChVector<>(+dim.x(), -dim.y(), +dim.z()) + loc);

    body->GetCollisionModel()->AddConvexHull(mat, points);

    auto shape = chrono_types::make_shared<ChTriangleMeshShape>();
    collision::ChConvexHullLibraryWrapper lh;
    lh.ComputeHull(points, *shape->GetMesh());
    body->AddAsset(shape);

    body->AddAsset(chrono_types::make_shared<ChColorAsset>(0.5f, 0.0f, 0.0f));
}

void BuildContainerBoxes(std::shared_ptr<ChBody> body,
                         std::shared_ptr<ChMaterialSurface> mat,
                         double hdimX,
                         double hdimY,
                         double hdimZ,
                         double hthick) {
  std::cout << "Using boxes for container" << std::endl;

  body->GetCollisionModel()->ClearModel();

  // Bottom box
  AddWallBox(body, mat, ChVector<>(hdimX, hthick, hdimY), ChVector<>(0, 0, 0));

  // Side box
  AddWallBox(body, mat, ChVector<>(hthick, hdimZ, hdimY), ChVector<>(hdimX - hthick, hdimZ, 0));

  body->GetCollisionModel()->BuildModel();
}

void BuildContainerMeshes(std::shared_ptr<ChBody> body,
                          std::shared_ptr<ChMaterialSurface> mat,
                          double hdimX,
                          double hdimY,
                          double hdimZ,
                          double hthick) {
  std::cout << "Using meshes for container" << std::endl;

  body->GetCollisionModel()->ClearModel();

  // Bottom mesh
  AddWallMesh(body, mat, ChVector<>(hdimX, hthick, hdimY), ChVector<>(0, 0, 0));

  // Side mesh
  AddWallMesh(body, mat, ChVector<>(hthick, hdimZ, hdimY), ChVector<>(hdimX - hthick, hdimZ, 0));

  body->GetCollisionModel()->BuildModel();
}

void BuildContainerHulls(std::shared_ptr<ChBody> body,
                         std::shared_ptr<ChMaterialSurface> mat,
                         double hdimX,
                         double hdimY,
                         double hdimZ,
                         double hthick) {
  std::cout << "Using convex hulls for container" << std::endl;

  body->GetCollisionModel()->ClearModel();

  // Bottom hull
  AddWallHull(body, mat, ChVector<>(hdimX, hthick, hdimY), ChVector<>(0, 0, 0));

  // Side hull
  AddWallHull(body, mat, ChVector<>(hthick, hdimZ, hdimY), ChVector<>(hdimX - hthick, hdimZ, 0));

  body->GetCollisionModel()->BuildModel();
}


int main(int argc, char* argv[]) {
  SetChronoDataPath(CHRONO_DATA_DIR);
  
  // Simulation parameters
    double gravity = -9.81;
    double time_step = 1e-4;

    // Parameters for the falling ball
    int ballId = 100;
    double radius = 0.5;
    double mass = 1000;
    ChVector<> pos(0.2, 2, 0.4);
    ChQuaternion<> rot(1, 0, 0, 0);
    ChVector<> init_vel(1, 0, 0);

    // Parameters for the containing bin
    int binId = 200;
    double hdimX = 3;
    double hdimY = 2;
    double hdimZ = 1;
    double hthick = 0.2;

    // Create the system (Y up)
    ChSystemSMC msystem;
    msystem.Set_G_acc(ChVector<>(0, gravity, 0));

    msystem.SetContactForceModel(ChSystemSMC::ContactForceModel::Hertz);
    msystem.SetAdhesionForceModel(ChSystemSMC::AdhesionForceModel::Constant);

    // Create a material (will be used by both objects)
    auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    material->SetYoungModulus(1.0e7f);
    material->SetRestitution(0.1f);
    material->SetFriction(0.4f);

    // Create the falling ball
    auto ball = chrono_types::make_shared<ChBody>();

    ball->SetIdentifier(ballId);
    ball->SetMass(mass);
    ball->SetInertiaXX(0.4 * mass * radius * radius * ChVector<>(1, 1, 1));
    ball->SetPos(pos);
    ball->SetRot(rot);
    ball->SetPos_dt(init_vel);
    // ball->SetWvel_par(ChVector<>(0,0,3));
    ball->SetBodyFixed(false);

    ball->SetCollide(true);
    ball->GetCollisionModel()->ClearModel();
    ball->GetCollisionModel()->AddSphere(material, radius);
    ball->GetCollisionModel()->BuildModel();

    auto sphere = chrono_types::make_shared<ChSphereShape>();
    sphere->GetSphereGeometry().rad = radius;
    ball->AddAsset(sphere);

    auto mtexture = chrono_types::make_shared<ChTexture>();
    mtexture->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
    ball->AddAsset(mtexture);

    msystem.AddBody(ball);

    // Create container
    auto bin = chrono_types::make_shared<ChBody>();

    bin->SetIdentifier(binId);
    bin->SetMass(1);
    bin->SetPos(ChVector<>(0, 0, 0));
    bin->SetRot(ChQuaternion<>(1, 0, 0, 0));
    bin->SetCollide(true);
    bin->SetBodyFixed(true);

    //BuildContainerBoxes(bin, material, hdimX, hdimY, hdimZ, hthick);
    //BuildContainerMeshes(bin, material, hdimX, hdimY, hdimZ, hthick);
    BuildContainerHulls(bin, material, hdimX, hdimY, hdimZ, hthick);

    msystem.AddBody(bin);

    // Create the Irrlicht visualization
    ChIrrApp application(&msystem, L"Collision test", irr::core::dimension2d<irr::u32>(800, 600), false, true);

    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(irr::core::vector3df(0, 3, -6));

    // Enable contact forces visualization in Irrlicht application
    application.SetSymbolscale(1e-4);
    application.SetContactsDrawMode(ChIrrTools::eCh_ContactsDrawMode::CONTACT_FORCES);

    // Complete asset construction
    application.AssetBindAll();
    application.AssetUpdateAll();

    // Simulation loop
    while (application.GetDevice()->run()) {
        application.BeginScene();
        application.DrawAll();
        application.EndScene();

        msystem.DoStepDynamics(time_step);
    }

    return 0;
}
