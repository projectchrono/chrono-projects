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
#include "chrono/collision/bullet/ChCollisionUtilsBullet.h"
#include "chrono/physics/ChSystemSMC.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

using namespace chrono;
using namespace chrono::irrlicht;

void AddWallBox(std::shared_ptr<ChBody> body,
                std::shared_ptr<ChContactMaterial> mat,
                const ChVector3d& dim,
                const ChVector3d& loc) {
    auto ct_shape = chrono_types::make_shared<ChCollisionShapeBox>(mat, dim);
    body->AddCollisionShape(ct_shape, ChFrame<>(loc, QUNIT));

    auto box = chrono_types::make_shared<ChVisualShapeBox>(dim);
    box->SetColor(ChColor(0.5f, 0.5f, 0.0f));
    body->AddVisualShape(box, ChFrame<>(loc));
}

void AddWallMesh(std::shared_ptr<ChBody> body,
                 std::shared_ptr<ChContactMaterial> mat,
                 const ChVector3d& dim,
                 const ChVector3d& loc) {
    auto trimesh = chrono_types::make_shared<geometry::ChTriangleMeshConnected>();

    std::vector<ChVector3d>& vertices = trimesh->getCoordsVertices();
    std::vector<ChVector3d>& normals = trimesh->getCoordsNormals();
    std::vector<ChVector3i>& idx_vertices = trimesh->getIndicesVertexes();
    std::vector<ChVector3i>& idx_normals = trimesh->getIndicesNormals();

    int num_vert = 8;
    int num_faces = 12;

    vertices.resize(num_vert);
    normals.resize(num_vert);
    trimesh->getCoordsUV().resize(num_vert);
    trimesh->getCoordsColors().resize(num_vert);

    idx_vertices.resize(num_faces);
    idx_normals.resize(num_faces);

    vertices[0] = ChVector3d(-dim.x(), -dim.y(), -dim.z()) + loc;
    vertices[1] = ChVector3d(-dim.x(), +dim.y(), -dim.z()) + loc;
    vertices[2] = ChVector3d(+dim.x(), +dim.y(), -dim.z()) + loc;
    vertices[3] = ChVector3d(+dim.x(), -dim.y(), -dim.z()) + loc;

    normals[0] = ChVector3d(-1, -1, -1).Normalize();
    normals[1] = ChVector3d(-1, +1, -1).Normalize();
    normals[2] = ChVector3d(+1, +1, -1).Normalize();
    normals[3] = ChVector3d(+1, -1, -1).Normalize();

    vertices[4] = ChVector3d(-dim.x(), -dim.y(), +dim.z()) + loc;
    vertices[5] = ChVector3d(-dim.x(), +dim.y(), +dim.z()) + loc;
    vertices[6] = ChVector3d(+dim.x(), +dim.y(), +dim.z()) + loc;
    vertices[7] = ChVector3d(+dim.x(), -dim.y(), +dim.z()) + loc;

    normals[4] = ChVector3d(-1, -1, +1).Normalize();
    normals[5] = ChVector3d(-1, +1, +1).Normalize();
    normals[6] = ChVector3d(+1, +1, +1).Normalize();
    normals[7] = ChVector3d(+1, -1, +1).Normalize();

    for (int i = 0; i < num_vert; i++) {
        trimesh->getCoordsColors()[i] = ChColor(0, 0, 1);
    }

    idx_vertices[0] = ChVector3i(0, 1, 3);
    idx_vertices[1] = ChVector3i(1, 2, 3);
    idx_vertices[2] = ChVector3i(4, 7, 5);
    idx_vertices[3] = ChVector3i(7, 6, 5);

    idx_vertices[4] = ChVector3i(0, 3, 4);
    idx_vertices[5] = ChVector3i(3, 7, 4);
    idx_vertices[6] = ChVector3i(1, 5, 2);
    idx_vertices[7] = ChVector3i(2, 5, 6);

    idx_vertices[8] = ChVector3i(3, 6, 7);
    idx_vertices[9] = ChVector3i(3, 2, 6);
    idx_vertices[10] = ChVector3i(0, 4, 5);
    idx_vertices[11] = ChVector3i(0, 5, 1);

    for (int i = 0; i < num_faces; i++)
        idx_normals[i] = idx_vertices[i];

    auto body_ct_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mat, trimesh, true, true, 0);
    body->AddCollisionShape(body_ct_shape, ChFrame<>(loc, QUNIT));

    auto trimesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    trimesh_shape->SetMesh(trimesh);
    trimesh_shape->SetColor(ChColor(0.0f, 0.0f, 0.5f));
    body->AddVisualShape(trimesh_shape);
}

void AddWallHull(std::shared_ptr<ChBody> body,
                 std::shared_ptr<ChContactMaterial> mat,
                 const ChVector3d& dim,
                 const ChVector3d& loc) {
    auto hdim = dim / 2;
    std::vector<ChVector3d> points;

    points.push_back(ChVector3d(-hdim.x(), -hdim.y(), -hdim.z()) + loc);
    points.push_back(ChVector3d(-hdim.x(), +hdim.y(), -hdim.z()) + loc);
    points.push_back(ChVector3d(+hdim.x(), +hdim.y(), -hdim.z()) + loc);
    points.push_back(ChVector3d(+hdim.x(), -hdim.y(), -hdim.z()) + loc);
    points.push_back(ChVector3d(-hdim.x(), -hdim.y(), +hdim.z()) + loc);
    points.push_back(ChVector3d(-hdim.x(), +hdim.y(), +hdim.z()) + loc);
    points.push_back(ChVector3d(+hdim.x(), +hdim.y(), +hdim.z()) + loc);
    points.push_back(ChVector3d(+hdim.x(), -hdim.y(), +hdim.z()) + loc);

    auto body_ct_shape = chrono_types::make_shared<ChCollisionShapeConvexHull>(mat, points);
    body->AddCollisionShape(body_ct_shape);

    auto shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    bt_utils::ChConvexHullLibraryWrapper lh;
    lh.ComputeHull(points, *shape->GetMesh());
    shape->SetColor(ChColor(0.5f, 0.0f, 0.0f));
    body->AddVisualShape(shape);
}

void BuildContainerBoxes(std::shared_ptr<ChBody> body,
                         std::shared_ptr<ChContactMaterial> mat,
                         double dimX,
                         double dimY,
                         double dimZ,
                         double thick) {
    std::cout << "Using boxes for container" << std::endl;

    // Bottom box
    AddWallBox(body, mat, ChVector3d(dimX, thick, dimY), ChVector3d(0, 0, 0));

    // Side box
    AddWallBox(body, mat, ChVector3d(thick, dimZ, dimY), ChVector3d(dimX / 2 - thick / 2, dimZ / 2, 0));
}

void BuildContainerMeshes(std::shared_ptr<ChBody> body,
                          std::shared_ptr<ChContactMaterial> mat,
                          double hdimX,
                          double hdimY,
                          double hdimZ,
                          double hthick) {
    std::cout << "Using meshes for container" << std::endl;

    // Bottom mesh
    AddWallMesh(body, mat, ChVector3d(hdimX, hthick, hdimY), ChVector3d(0, 0, 0));

    // Side mesh
    AddWallMesh(body, mat, ChVector3d(hthick, hdimZ, hdimY), ChVector3d(hdimX - hthick, hdimZ, 0));
}

void BuildContainerHulls(std::shared_ptr<ChBody> body,
                         std::shared_ptr<ChContactMaterial> mat,
                         double dimX,
                         double dimY,
                         double dimZ,
                         double thick) {
    std::cout << "Using convex hulls for container" << std::endl;

    // Bottom hull
    AddWallHull(body, mat, ChVector3d(dimX, thick, dimY), ChVector3d(0, 0, 0));

    // Side hull
    AddWallHull(body, mat, ChVector3d(thick, dimZ, dimY), ChVector3d(dimX / 2 - thick / 2, dimZ / 2, 0));
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
    ChVector3d pos(0.2, 2, 0.4);
    ChQuaternion<> rot(1, 0, 0, 0);
    ChVector3d init_vel(1, 0, 0);

    // Parameters for the containing bin
    int binId = 200;
    double dimX = 6;
    double dimY = 4;
    double dimZ = 2;
    double thick = 0.4;

    // Create the system (Y up)
    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.Set_G_acc(ChVector3d(0, gravity, 0));

    sys.SetContactForceModel(ChSystemSMC::ContactForceModel::Hertz);
    sys.SetAdhesionForceModel(ChSystemSMC::AdhesionForceModel::Constant);

    // Create a material (will be used by both objects)
    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetYoungModulus(1.0e7f);
    material->SetRestitution(0.1f);
    material->SetFriction(0.4f);

    // Create the falling ball
    auto ball = chrono_types::make_shared<ChBody>();

    ball->SetIdentifier(ballId);
    ball->SetMass(mass);
    ball->SetInertiaXX(0.4 * mass * radius * radius * ChVector3d(1, 1, 1));
    ball->SetPos(pos);
    ball->SetRot(rot);
    ball->SetPos_dt(init_vel);
    // ball->SetWvel_par(ChVector3d(0,0,3));
    ball->SetBodyFixed(false);
    ball->SetCollide(true);

    auto ball_ct_shape = chrono_types::make_shared<ChCollisionShapeSphere>(material, radius);
    ball->AddCollisionShape(ball_ct_shape);

    auto sphere = chrono_types::make_shared<ChVisualShapeSphere>(radius);
    sphere->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
    ball->AddVisualShape(sphere);

    sys.AddBody(ball);

    // Create container
    auto bin = chrono_types::make_shared<ChBody>();

    bin->SetIdentifier(binId);
    bin->SetMass(1);
    bin->SetPos(ChVector3d(0, 0, 0));
    bin->SetRot(ChQuaternion<>(1, 0, 0, 0));
    bin->SetCollide(true);
    bin->SetBodyFixed(true);

    BuildContainerHulls(bin, material, dimX, dimY, dimZ, thick);

    sys.AddBody(bin);

    // Create the Irrlicht visualization
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Collision test");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector3d(0, 3, -6));
    vis->AddTypicalLights();
    vis->SetSymbolScale(1e-4);
    vis->EnableContactDrawing(ContactsDrawMode::CONTACT_FORCES);
    vis->AttachSystem(&sys);

    // Simulation loop
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        sys.DoStepDynamics(time_step);
    }

    return 0;
}
