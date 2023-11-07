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
// Authors:
// =============================================================================
//
// Test using 4 flexible tires (Abaqus) on a simple vehicle framework.
//
// =============================================================================

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoaderUV.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

void MakeWheel(ChSystemSMC& sys,
               const ChVector<> tire_center,
               const ChQuaternion<> tire_alignment,
               const double tire_scale_R,
               const double tire_scale_W,
               const double tire_w0,
               const double tire_vel_z0,
               std::shared_ptr<ChMaterialSurface> mysurfmaterial,
               std::shared_ptr<ChContinuumElastic> mmaterial,
               std::shared_ptr<ChBody>& mrim) {
    ChMatrix33<> malign(tire_alignment);
    ChMatrix33<> mscale;
    mscale(0, 0) = tire_scale_W;
    mscale(1, 1) = tire_scale_R;
    mscale(2, 2) = tire_scale_R;

    // Create a mesh, that is a container for groups
    // of FEA elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Load an ABAQUS .INP tetahedron mesh file from disk, defining a tetahedron mesh.
    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets;

    try {
        ChMeshFileLoader::FromAbaqusFile(my_mesh, GetChronoDataFile("fea/tractor_wheel_coarse.INP").c_str(), mmaterial,
                                         node_sets, tire_center, mscale * malign);
    } catch (ChException myerr) {
        GetLog() << myerr.what();
        return;
    }

    // Create the contact surface(s).
    // In this case it is a ChContactSurfaceNodeCloud, so just pass
    // all nodes to it.
    auto mcontactsurf = chrono_types::make_shared<ChContactSurfaceNodeCloud>(mysurfmaterial);
    my_mesh->AddContactSurface(mcontactsurf);

    mcontactsurf->AddAllNodes();

    // Apply initial speed and angular speed
    double speed_x0 = 0.5;
    for (unsigned int i = 0; i < my_mesh->GetNnodes(); ++i) {
        auto node_pos = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->GetPos();
        ChVector<> tang_vel = Vcross(ChVector<>(tire_w0, 0, 0), node_pos - tire_center);
        std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))
            ->SetPos_dt(ChVector<>(0, 0, tire_vel_z0) + tang_vel);
    }

    // Remember to add the mesh to the system!
    sys.Add(my_mesh);

    // Add a rim
    auto mwheel_rim = chrono_types::make_shared<ChBody>();
    mwheel_rim->SetMass(80);
    mwheel_rim->SetInertiaXX(ChVector<>(60, 60, 60));
    mwheel_rim->SetPos(tire_center);
    mwheel_rim->SetRot(tire_alignment);
    mwheel_rim->SetPos_dt(ChVector<>(0, 0, tire_vel_z0));
    mwheel_rim->SetWvel_par(ChVector<>(tire_w0, 0, 0));
    sys.Add(mwheel_rim);

    auto mobjmesh = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    mobjmesh->GetMesh()->LoadWavefrontMesh(GetChronoDataFile("fea/tractor_wheel_rim.obj"));
    mobjmesh->GetMesh()->Transform(VNULL, mscale);
    mwheel_rim->AddVisualShape(mobjmesh);

    mrim = mwheel_rim;

    // Conect rim and tire using constraints.
    // the BC_RIMTIRE nodeset, in the Abaqus INP file, lists the nodes involved
    auto nodeset_sel = "BC_RIMTIRE";
    for (auto i = 0; i < node_sets.at(nodeset_sel).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets[nodeset_sel][i]), mwheel_rim);
        sys.Add(mlink);
    }

    /// Create a mesh surface, for applying loads:
    auto mmeshsurf = chrono_types::make_shared<ChMeshSurface>();
    my_mesh->AddMeshSurface(mmeshsurf);

    // In the .INP file there are two additional NSET nodesets, the 1st is used to mark load surface:
    nodeset_sel = "BC_SURF";
    mmeshsurf->AddFacesFromNodeSet(node_sets[nodeset_sel]);

    /// Apply load to all surfaces in the mesh surface
    auto mloadcontainer = chrono_types::make_shared<ChLoadContainer>();
    sys.Add(mloadcontainer);

    for (int i = 0; i < mmeshsurf->GetFacesList().size(); ++i) {
        auto aface = std::shared_ptr<ChLoadableUV>(mmeshsurf->GetFacesList()[i]);
        auto faceload = chrono_types::make_shared<ChLoad<ChLoaderPressure>>(aface);
        faceload->loader.SetPressure(10000);  // low pressure... the tire has no ply!
        mloadcontainer->Add(faceload);
    }
    // ==Asset== attach a visualization of the FEM mesh.
    auto mvisualizemesh = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    mvisualizemesh->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
    mvisualizemesh->SetColorscaleMinMax(0.0, 10);
    mvisualizemesh->SetSmoothFaces(true);
    my_mesh->AddVisualShapeFEA(mvisualizemesh);
}

int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Global parameter for tire:
    double tire_rad = 0.8;
    double tire_vel_z0 = -3;
    ChVector<> vehicle_center(0, 0.9, 0.5);
    ChVector<> rel_tire_center_BL(1, -0.2, 1.7);
    ChVector<> tire_center_BL = vehicle_center + rel_tire_center_BL;
    ChVector<> rel_tire_center_BR(-1, -0.2, 1.7);
    ChVector<> tire_center_BR = vehicle_center + rel_tire_center_BR;
    ChVector<> rel_tire_center_FL(1, -0.2, -1.7);
    ChVector<> tire_center_FL = vehicle_center + rel_tire_center_FL;
    ChVector<> rel_tire_center_FR(-1, -0.2, -1.7);
    ChVector<> tire_center_FR = vehicle_center + rel_tire_center_FR;

    ChQuaternion<> tire_alignment = Q_from_AngAxis(CH_C_PI, VECT_Y);  // create rotated 180� on y
    double tire_scaleR = 0.85;
    double tire_scaleW = 0.95;

    double tire_w0 = tire_vel_z0 / tire_rad;

    // Create a Chrono::Engine physical system
    ChSystemSMC sys;


    //
    // CREATE THE PHYSICAL SYSTEM
    //

    // Create the surface material, containing information
    // about friction etc.
    auto mysurfmaterial = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mysurfmaterial->SetYoungModulus(10e4);
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0.2f);
    mysurfmaterial->SetAdhesion(0);
    auto mysurfmaterial2 = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mysurfmaterial->SetYoungModulus(26e4);
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0.2f);
    mysurfmaterial->SetAdhesion(0);

    // RIGID BODIES
    // Create some rigid bodies, for instance a floor:
    auto mfloor = chrono_types::make_shared<ChBodyEasyBox>(15, 0.2, 15, 2700, true, true, mysurfmaterial);
    mfloor->SetBodyFixed(true);
    mfloor->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    sys.Add(mfloor);

    // Create the car truss
    auto mtruss = chrono_types::make_shared<ChBodyAuxRef>();
    mtruss->SetPos(vehicle_center);
    mtruss->SetPos_dt(ChVector<>(0, 0, tire_vel_z0));
    mtruss->SetBodyFixed(false);
    mtruss->SetMass(100);
    mtruss->SetInertiaXX(ChVector<>(100, 100, 100));
    sys.Add(mtruss);

    auto mtrussmesh = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    mtrussmesh->GetMesh()->LoadWavefrontMesh(GetChronoDataFile("vehicle/hmmwv/hmmwv_chassis_simple.obj"));
    mtruss->AddVisualShape(mtrussmesh,
                           ChFrame<>(VNULL, Q_from_AngAxis(CH_C_PI_2, VECT_Z) % Q_from_AngAxis(CH_C_PI_2, VECT_Y)));

    // Create a step
    if (true) {
        auto mfloor_step = chrono_types::make_shared<ChBodyEasyBox>(3, 0.2, 0.5, 2700, true, true, mysurfmaterial);
        mfloor_step->SetPos(ChVector<>(2, 0.1, -1.8));
        mfloor_step->SetBodyFixed(true);
        sys.Add(mfloor_step);
    }

    // Create some bent rectangular fixed slabs
    if (false) {
        for (int i = 0; i < 50; ++i) {
            auto mcube = chrono_types::make_shared<ChBodyEasyBox>(0.25, 0.2, 0.25, 2700, true, true, mysurfmaterial);
            ChQuaternion<> vrot;
            vrot.Q_from_AngAxis(ChRandom() * CH_C_2PI, VECT_Y);
            mcube->Move(ChCoordsys<>(VNULL, vrot));
            vrot.Q_from_AngAxis((ChRandom() - 0.5) * 2 * CH_C_DEG_TO_RAD * 20,
                                ChVector<>(ChRandom() - 0.5, 0, ChRandom() - 0.5).Normalize());
            mcube->Move(ChCoordsys<>(VNULL, vrot));
            mcube->SetPos(ChVector<>((ChRandom() - 0.5) * 2.8, ChRandom() * 0.1, -ChRandom() * 3.2 - 1.1));
            mcube->SetBodyFixed(true);
            mcube->GetVisualShape(0)->SetColor(ChColor(0.3f, 0.3f, 0.3f));
            sys.Add(mcube);
        }
    }

    // Create some stones / obstacles on the ground
    if (false) {
        for (int i = 0; i < 150; ++i) {
            auto mcube = chrono_types::make_shared<ChBodyEasyBox>(0.18, 0.04, 0.18, 2700, true, true, mysurfmaterial2);
            ChQuaternion<> vrot;
            vrot.Q_from_AngAxis(ChRandom() * CH_C_2PI, VECT_Y);
            mcube->Move(ChCoordsys<>(VNULL, vrot));
            mcube->SetPos(ChVector<>((ChRandom() - 0.5) * 1.4, ChRandom() * 0.2 + 0.05, -ChRandom() * 2.6 + 0.2));
            mcube->GetVisualShape(0)->SetColor(ChColor(0.3f, 0.3f, 0.3f));
            sys.Add(mcube);
        }
    }

    //
    // THE WHEELS
    //

    // Create a material, that must be assigned to rubber of tires
    auto mtirematerial = chrono_types::make_shared<ChContinuumElastic>();
    mtirematerial->Set_E(0.016e9);  // rubber 0.01e9, steel 200e9
    mtirematerial->Set_v(0.4);
    mtirematerial->Set_RayleighDampingK(0.004);
    mtirematerial->Set_density(1000);

    // Make a wheel and connect it to truss:
    std::shared_ptr<ChBody> mrim_BL;
    MakeWheel(sys, tire_center_BL, tire_alignment, tire_scaleR, tire_scaleW, tire_w0, tire_vel_z0, mysurfmaterial,
              mtirematerial, mrim_BL);

    auto mrevolute_BL = chrono_types::make_shared<ChLinkLockRevolute>();
    sys.Add(mrevolute_BL);
    mrevolute_BL->Initialize(mtruss, mrim_BL, ChCoordsys<>(tire_center_BL, Q_from_AngAxis(CH_C_PI_2, VECT_Y)));

    // Make a wheel and connect it to truss:
    std::shared_ptr<ChBody> mrim_FL;
    MakeWheel(sys, tire_center_FL, tire_alignment, tire_scaleR, tire_scaleW, tire_w0, tire_vel_z0, mysurfmaterial,
              mtirematerial, mrim_FL);
    auto mrevolute_FL = chrono_types::make_shared<ChLinkLockRevolute>();
    sys.Add(mrevolute_FL);
    mrevolute_FL->Initialize(mtruss, mrim_FL, ChCoordsys<>(tire_center_FL, Q_from_AngAxis(CH_C_PI_2, VECT_Y)));

    // Make a wheel and connect it to truss:
    std::shared_ptr<ChBody> mrim_BR;
    MakeWheel(sys, tire_center_BR, tire_alignment, tire_scaleR, tire_scaleW, tire_w0, tire_vel_z0, mysurfmaterial,
              mtirematerial, mrim_BR);
    auto mrevolute_BR = chrono_types::make_shared<ChLinkLockRevolute>();
    sys.Add(mrevolute_BR);
    mrevolute_BR->Initialize(mtruss, mrim_BR, ChCoordsys<>(tire_center_BR, Q_from_AngAxis(CH_C_PI_2, VECT_Y)));

    // Make a wheel and connect it to truss:
    std::shared_ptr<ChBody> mrim_FR;
    MakeWheel(sys, tire_center_FR, tire_alignment, tire_scaleR, tire_scaleW, tire_w0, tire_vel_z0, mysurfmaterial,
              mtirematerial, mrim_FR);

    auto mrevolute_FR = chrono_types::make_shared<ChLinkLockRevolute>();
    sys.Add(mrevolute_FR);
    mrevolute_FR->Initialize(mtruss, mrim_FR, ChCoordsys<>(tire_center_FR, Q_from_AngAxis(CH_C_PI_2, VECT_Y)));

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(1280, 720);
    vis->SetWindowTitle("FEA contacts");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(3, 3, -4), ChVector<>(0, tire_rad, 0));
    vis->AddTypicalLights();
    vis->AttachSystem(&sys);

    //
    // THE SOFT-REAL-TIME CYCLE
    //

    /*
        // Change solver to embedded MINRES
        auto solver = chrono_types::make_shared<ChSolverMINRES>();
        solver->EnableWarmStart(true);
        solver->SetMaxIterations(90);
        solver->SetVerbose(false);
        sys.SetSolver(solver);
        sys.SetSolverForceTolerance(1e-10);
    */
    // Change solver to pluggable PardisoMKL
    auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
    sys.SetSolver(mkl_solver);
    mkl_solver->LockSparsityPattern(true);
    sys.Update();

    // Change type of integrator:
    sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);  // fast, less precise
    // sys.SetTimestepperType(chrono::ChTimestepper::Type::HHT);  // precise,slower, might iterate each step

    // if later you want to change integrator settings:
    if (auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(sys.GetTimestepper())) {
        mystepper->SetAlpha(-0.2);
        mystepper->SetMaxiters(2);
        mystepper->SetAbsTolerances(1e-6);
    }

    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(1e-3);
    }

    return 0;
}
