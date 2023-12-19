// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora
// =============================================================================
//
// Demo code about loading a .chulls file, with xyz points of cluters of convex
// hulls that define a complicate concave shape. The shape is a wheel for
// tractors, with large knobs, that has been decomposed using demo_decomposition
// from .obj shape to a .chull.
//
// =============================================================================

#include "chrono/core/ChRealtimeStep.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include <irrlicht.h>

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::irrlicht;

std::shared_ptr<ChBody> create_wheel(ChVector<> mposition, ChSystem& sys) {
    ChCollisionModel::SetDefaultSuggestedEnvelope(0.005);
    ChCollisionModel::SetDefaultSuggestedMargin(0.004);

    // create a basic rigid body, it comes with no visualization or collision shapes
    auto mrigidBody = chrono_types::make_shared<ChBody>();
    sys.Add(mrigidBody);
    mrigidBody->SetMass(50);
    mrigidBody->SetInertiaXX(ChVector<>(10, 10, 10));
    mrigidBody->SetPos(mposition);

    // now attach a visualization shape, as a mesh from disk
    auto tireMesh = chrono_types::make_shared<ChVisualShapeModelFile>();
    tireMesh->SetFilename(GetChronoDataFile("models/tractor_wheel/tractor_wheel.obj").c_str());
    mrigidBody->AddVisualShape(tireMesh);

    // contact material shared by all collision shapes of the wheel
    auto mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat->SetFriction(0.5f);

    // now attach collision shape, as a compound of convex hulls (for each thread pair):
    std::string knobs_filename(GetChronoDataFile("models/tractor_wheel/tractor_wheel_knobs.chulls"));
    std::string slice_filename(GetChronoDataFile("models/tractor_wheel/tractor_wheel_slice.chulls"));
    auto knobs_shapes = ChCollisionShapeConvexHull::Read(mat, knobs_filename);
    auto slice_shapes = ChCollisionShapeConvexHull::Read(mat, slice_filename);
    for (double mangle = 0; mangle < 360.; mangle += (360. / 15.)) {
        auto q = Q_from_AngX(mangle * CH_C_DEG_TO_RAD);
        for (const auto& s : knobs_shapes)
            mrigidBody->AddCollisionShape(s, ChFrame<>(VNULL, q));
        for (const auto& s : slice_shapes)
            mrigidBody->AddCollisionShape(s, ChFrame<>(VNULL, q));
    }
    mrigidBody->SetCollide(true);

    return mrigidBody;
}

void create_some_falling_items(ChSystemNSC& sys) {
    ChCollisionModel::SetDefaultSuggestedEnvelope(0.003);
    ChCollisionModel::SetDefaultSuggestedMargin(0.002);

    ChQuaternion<> rot;
    rot.Q_from_AngAxis(ChRandom() * CH_C_2PI, VECT_Y);

    auto pebble_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    pebble_mat->SetFriction(0.4f);

    double bed_x = 0.6;
    double bed_z = 1;

    int n_pebbles = 30;
    for (int bi = 0; bi < n_pebbles; bi++) {
        double sphrad = 0.02 + 0.02 * ChRandom();
        double sphdens = 1;
        ChQuaternion<> randrot(ChRandom(), ChRandom(), ChRandom(), ChRandom());
        randrot.Normalize();

        auto mrigidBody = chrono_types::make_shared<ChBodyEasySphere>(sphrad, sphdens, true, true, pebble_mat);
        sys.Add(mrigidBody);
        mrigidBody->SetRot(randrot);
        mrigidBody->SetPos(ChVector<>(-0.5 * bed_x + ChRandom() * bed_x, 0.01 + 0.04 * ((double)bi / (double)n_pebbles),
                                      -0.5 * bed_z + ChRandom() * bed_z));
    }

    // Create the a plane using body of 'box' type:
    auto ground_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    ground_mat->SetFriction(0.5f);

    auto mrigidBodyB = chrono_types::make_shared<ChBodyEasyBox>(10, 1, 10, 1000, true, true, ground_mat);
    sys.Add(mrigidBodyB);
    mrigidBodyB->SetBodyFixed(true);
    mrigidBodyB->SetPos(ChVector<>(0, -0.5, 0));
}

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Create a Chrono physical system
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Create some debris
    create_some_falling_items(sys);

    // Create the wheel
    std::shared_ptr<ChBody> mwheelBody = create_wheel(ChVector<>(0, 1, 0), sys);

    // Create the Irrlicht visualization sys
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Convex decomposed wheel");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(3.5, 2.5, -2.4));
    vis->AddTypicalLights();

    // Simulation loop
    double timestep = 0.01;
    ChRealtimeStepTimer realtime_timer;
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(timestep);
        realtime_timer.Spin(timestep);
    }

    return 0;
}
