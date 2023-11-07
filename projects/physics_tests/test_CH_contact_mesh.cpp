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
// Test mesh collision
//
// =============================================================================

#include "chrono/ChConfig.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

using namespace chrono;
using namespace chrono::collision;
using namespace chrono::irrlicht;

// ====================================================================================

int main(int argc, char* argv[]) {
    // ---------------------------------
    // Set path to Chrono data directory
    // ---------------------------------
    SetChronoDataPath(CHRONO_DATA_DIR);

    // ---------------------
    // Simulation parameters
    // ---------------------

    double gravity = 9.81;    // gravitational acceleration
    double time_step = 1e-4;  // integration step size

    // ---------------------------
    // Contact material properties
    // ---------------------------
    ChContactMethod contact_method = ChContactMethod::SMC;
    bool use_mat_properties = true;

    float object_friction = 0.9f;
    float object_restitution = 0.1f;
    float object_young_modulus = 2e7f;
    float object_poisson_ratio = 0.3f;
    float object_adhesion = 0.0f;
    float object_kn = 2e5;
    float object_gn = 40;
    float object_kt = 2e5;
    float object_gt = 20;

    float ground_friction = 0.9f;
    float ground_restitution = 0.01f;
    float ground_young_modulus = 2e7f;
    float ground_poisson_ratio = 0.3f;
    float ground_adhesion = 0.0f;
    float ground_kn = 2e5;
    float ground_gn = 40;
    float ground_kt = 2e5;
    float ground_gt = 20;

    // ---------------------------------
    // Parameters for the falling object
    // ---------------------------------

    int objectId = 100;
    double mass = 1000;
    ChVector<> pos(0, 0, 1);
    ChQuaternion<> rot(1, 0, 0, 0);
    ChVector<> init_vel(0, 0, 0);
    ChVector<> init_omg(0, 0, 0);

    // ---------------------------------
    // Parameters for the containing bin
    // ---------------------------------

    int groundId = 200;
    double width = 4;
    double length = 2;
    double thickness = 0.2;

    // -----------------
    // Create the system
    // -----------------

    ChSystem* system;

    switch (contact_method) {
        case ChContactMethod::NSC:
            system = new ChSystemNSC();
            break;
        case ChContactMethod::SMC:
            system = new ChSystemSMC(use_mat_properties);
            break;
    }

    system->Set_G_acc(ChVector<>(0, 0, -gravity));

    // Create the falling object
    auto object = std::shared_ptr<ChBody>(system->NewBody());
    system->AddBody(object);

    object->SetIdentifier(objectId);
    object->SetMass(mass);
    object->SetInertiaXX(400.0 * ChVector<>(1, 1, 1));
    object->SetPos(pos);
    object->SetRot(rot);
    object->SetPos_dt(init_vel);
    object->SetWvel_par(init_omg);
    object->SetCollide(true);
    object->SetBodyFixed(false);

    std::shared_ptr<ChMaterialSurface> object_mat;
    switch (contact_method) {
        case ChContactMethod::NSC: {
            auto matNSC = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            matNSC->SetFriction(object_friction);
            matNSC->SetRestitution(object_restitution);
            object_mat = matNSC;
            break;
        }
        case ChContactMethod::SMC: {
            auto matSMC = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            matSMC->SetFriction(object_friction);
            matSMC->SetRestitution(object_restitution);
            matSMC->SetYoungModulus(object_young_modulus);
            matSMC->SetPoissonRatio(object_poisson_ratio);
            matSMC->SetKn(object_kn);
            matSMC->SetGn(object_gn);
            matSMC->SetKt(object_kt);
            matSMC->SetGt(object_gt);
            object_mat = matSMC;
            break;
        }
    }

    auto trimesh = chrono_types::make_shared<geometry::ChTriangleMeshConnected>();
    trimesh->LoadWavefrontMesh(GetChronoDataFile("vehicle/hmmwv/hmmwv_tire.obj"), true, false);

    object->GetCollisionModel()->Clear();
    auto object_ct_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(object_mat, trimesh, false, false, 0.01);
    object->GetCollisionModel()->AddShape(object_ct_shape);
    object->GetCollisionModel()->Build();

    auto trimesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    trimesh_shape->SetMesh(trimesh);
    trimesh_shape->SetColor(ChColor(0.3f, 0.3f, 0.3f));
    object->AddVisualShape(trimesh_shape);

    // Create ground body
    auto ground = std::shared_ptr<ChBody>(system->NewBody());
    system->AddBody(ground);

    ground->SetIdentifier(groundId);
    ground->SetMass(1);
    ground->SetPos(ChVector<>(0, 0, 0));
    ground->SetRot(ChQuaternion<>(1, 0, 0, 0));
    ground->SetCollide(true);
    ground->SetBodyFixed(true);

    std::shared_ptr<ChMaterialSurface> ground_mat;
    switch (contact_method) {
        case ChContactMethod::NSC: {
            auto matNSC = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            matNSC->SetFriction(ground_friction);
            matNSC->SetRestitution(ground_restitution);
            ground_mat = matNSC;
            break;
        }
        case ChContactMethod::SMC: {
            auto matSMC = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            matSMC->SetFriction(ground_friction);
            matSMC->SetRestitution(ground_restitution);
            matSMC->SetYoungModulus(ground_young_modulus);
            matSMC->SetPoissonRatio(ground_poisson_ratio);
            matSMC->SetKn(ground_kn);
            matSMC->SetGn(ground_gn);
            matSMC->SetKt(ground_kt);
            matSMC->SetGt(ground_gt);
            ground_mat = matSMC;
            break;
        }
    }

    auto ground_ct_shape = chrono_types::make_shared<ChCollisionShapeBox>(ground_mat, width, length, thickness);
    ground->GetCollisionModel()->AddShape(ground_ct_shape, ChFrame<>(ChVector<>(0, 0, -thickness/2), QUNIT));
    ground->GetCollisionModel()->Build();

    auto box = chrono_types::make_shared<ChVisualShapeBox>(width, length, thickness);
    ground->AddVisualShape(box, ChFrame<>(ChVector<>(0, 0, -thickness/2)));

    // Create the Irrlicht visualization
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Mesh collision test");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0, 4, -0.2));
    vis->AddTypicalLights();
    vis->SetSymbolScale(1e-4);
    vis->EnableContactDrawing(ContactsDrawMode::CONTACT_FORCES);
    vis->AttachSystem(system);

    // ---------------
    // Simulation loop
    // ---------------
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        system->DoStepDynamics(time_step);
    }

    return 0;
}
