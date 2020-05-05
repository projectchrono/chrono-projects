//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

#include "chrono/ChConfig.h"
#include "chrono/assets/ChColorAsset.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/assets/ChObjShapeFile.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoaderUV.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemSMC.h"

#include "chrono/fea/ChElementTetra_4.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"
#include "chrono/fea/ChLinkPointFrame.h"

#ifdef CHRONO_IRRLICHT
#include "chrono_irrlicht/ChIrrApp.h"
#endif

#ifdef CHRONO_MKL
#include "chrono_mkl/ChSolverMKL.h"
#endif

#ifdef CHRONO_OPENMP_ENABLED
#include <omp.h>
#endif

using namespace chrono;
using namespace chrono::fea;

#ifdef CHRONO_IRRLICHT
using namespace chrono::irrlicht;
using namespace irr;
#endif

int main(int argc, char* argv[]) {
    // ---------------------------------
    // Set path to Chrono data directory
    // ---------------------------------
    SetChronoDataPath(CHRONO_DATA_DIR);

    // --------------------------
    // Simulation parameters
    // --------------------------

    int num_threads = 2;

    enum SolverType { MINRES, MKL };
    SolverType solver_type = MINRES;

    enum IntegratorType { EULER, HHT };
    IntegratorType integrator_type = HHT;

    double step_size = 5e-4;
    int num_steps = 20;

    bool visualization = false;

    // --------------------------
    // Model parameters
    // --------------------------

    double tire_vel_z0 = 0;//// -3;

    bool include_wheel_body = true;
    bool include_tire_contact = false;
    bool include_tire_pressure = false;

    bool include_obstacles = false;

    // --------------------------
    // Set number of threads
    // --------------------------

#ifdef CHRONO_OPENMP_ENABLED
    int max_threads = CHOMPfunctions::GetNumProcs();

    if (num_threads > max_threads)
        num_threads = max_threads;

    CHOMPfunctions::SetNumThreads(num_threads);
    GetLog() << "Using " << num_threads << " thread(s)\n";
#else
    GetLog() << "No OpenMP\n";
#endif

    // --------------------------
    // Create the physical system
    // --------------------------

    ChSystemSMC my_system;

    // Global parameter for tire
    double tire_rad = 0.8;
    double tire_w0 = tire_vel_z0 / tire_rad;
    ChVector<> tire_center(0, 0.02 + tire_rad, 0.5);
    ChMatrix33<> tire_alignment(Q_from_AngAxis(CH_C_PI, VECT_Y));

    // Contact material
    auto mysurfmaterial = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mysurfmaterial->SetYoungModulus(30e4);
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0.2f);
    mysurfmaterial->SetAdhesion(0);

    // ---------------
    // Create FEA mesh
    // ---------------

    // Mesh material
    auto mesh_material = chrono_types::make_shared<ChContinuumElastic>();
    mesh_material->Set_E(0.016e9);  // rubber 0.01e9, steel 200e9
    mesh_material->Set_v(0.4);
    mesh_material->Set_RayleighDampingK(0.004);
    mesh_material->Set_density(1000);

    // Create tire mesh from ABAQUS input file
    auto my_mesh = chrono_types::make_shared<ChMesh>();
    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets;

    try {
        ChMeshFileLoader::FromAbaqusFile(my_mesh, GetChronoDataFile("fea/tractor_wheel_coarse.INP").c_str(),
                                         mesh_material, node_sets, tire_center, tire_alignment);
    } catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    // Mesh visualization
    if (visualization) {
        auto mvisualizemesh = chrono_types::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
        mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
        mvisualizemesh->SetColorscaleMinMax(0.0, 10);
        mvisualizemesh->SetSmoothFaces(true);
        my_mesh->AddAsset(mvisualizemesh);
    }

    // Apply initial speed and angular speed
    for (unsigned int i = 0; i< my_mesh->GetNnodes(); ++i) {
        ChVector<> node_pos = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->GetPos();
        ChVector<> tang_vel = Vcross(ChVector<>(tire_w0, 0, 0), node_pos - tire_center);
        std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->SetPos_dt(ChVector<>(0, 0, tire_vel_z0) + tang_vel);
    }

    // Create the contact surface
    if (include_tire_contact) {
        // In this case it is a ChContactSurfaceNodeCloud, so just pass  all nodes to it.
        auto mcontactsurf = chrono_types::make_shared<ChContactSurfaceNodeCloud>(mysurfmaterial);
        my_mesh->AddContactSurface(mcontactsurf);
        mcontactsurf->AddAllNodes();
    }

    // Create tire pressure load
    if (include_tire_pressure) {
        // Create a mesh surface, for applying loads:
        auto mmeshsurf = chrono_types::make_shared<ChMeshSurface>();
        my_mesh->AddMeshSurface(mmeshsurf);

        // In the .INP file there are two additional NSET nodesets, the 1st is used to mark load surface:
        auto nodeset_sel = "BC_SURF";
        mmeshsurf->AddFacesFromNodeSet(node_sets[nodeset_sel]);

        // Apply load to all surfaces in the mesh surface
        auto mloadcontainer = chrono_types::make_shared<ChLoadContainer>();
        my_system.Add(mloadcontainer);

        for (int i = 0; i < mmeshsurf->GetFacesList().size(); ++i) {
            auto aface = std::shared_ptr<ChLoadableUV>(mmeshsurf->GetFacesList()[i]);
            auto faceload = chrono_types::make_shared<ChLoad< ChLoaderPressure > >(aface);
            faceload->loader.SetPressure(10000); // low pressure... the tire has no ply!
            mloadcontainer->Add(faceload);
        }
    }

    my_system.Add(my_mesh);

    // -------------------
    // Create rigid bodies
    // -------------------

    // Create the rim body
    if (include_wheel_body) {
        auto wheel = chrono_types::make_shared<ChBody>();
        wheel->SetMass(80);
        wheel->SetInertiaXX(ChVector<>(60, 60, 60));
        wheel->SetPos(tire_center);
        wheel->SetRot(tire_alignment);
        wheel->SetPos_dt(ChVector<>(0, 0, tire_vel_z0));
        wheel->SetWvel_par(ChVector<>(tire_w0, 0, 0));
        my_system.Add(wheel);

        if (visualization) {
            auto mobjmesh = chrono_types::make_shared<ChObjShapeFile>();
            mobjmesh->SetFilename(GetChronoDataFile("fea/tractor_wheel_rim.obj"));
            wheel->AddAsset(mobjmesh);
        }

        // Conect rim and tire using constraints.
        // the BC_RIMTIRE nodeset, in the Abaqus INP file, lists the nodes involved
        auto nodeset_sel = "BC_RIMTIRE";
        for (auto i = 0; i < node_sets.at(nodeset_sel).size(); ++i) {
            auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
            mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets[nodeset_sel][i]), wheel);
            my_system.Add(mlink);
        }
    }

    // Create ground
    auto mfloor = chrono_types::make_shared<ChBodyEasyBox>(2, 0.2, 6, 2700, true, true, mysurfmaterial);
    mfloor->SetBodyFixed(true);
    my_system.Add(mfloor);

    if (visualization) {
        auto mtexture = chrono_types::make_shared<ChTexture>();
        mtexture->SetTextureFilename(GetChronoDataFile("concrete.jpg"));
        mfloor->AddAsset(mtexture);
    }

    // Create obstacles
    if (include_obstacles) {
        for (int i = 0; i < 150; ++i) {
            auto mcube = chrono_types::make_shared<ChBodyEasyBox>(0.18, 0.04, 0.18, 2700, true, true, mysurfmaterial);
            ChQuaternion<> vrot;
            vrot.Q_from_AngAxis(ChRandom() * CH_C_2PI, VECT_Y);
            mcube->Move(ChCoordsys<>(VNULL, vrot));
            mcube->SetPos(ChVector<>((ChRandom() - 0.5) * 1.4, ChRandom() * 0.2 + 0.05, -ChRandom() * 2.6 + 0.2));
            my_system.Add(mcube);

            if (visualization) {
                auto mcubecol = chrono_types::make_shared<ChColorAsset>();
                mcubecol->SetColor(ChColor(0.3f, 0.3f, 0.3f));
                mcube->AddAsset(mcubecol);
            }
        }
    }

    // ----------------------------
    // Set up solver and integrator
    // ----------------------------

    // Set up solver
#ifndef CHRONO_MKL
    solver_type = MINRES;
#endif

    switch (solver_type) {
        case MINRES: {
            GetLog() << "Using MINRES solver\n";
            auto solver = chrono_types::make_shared<ChSolverMINRES>();
            solver->EnableWarmStart(true);
            solver->SetMaxIterations(40);
            solver->SetVerbose(false);
            my_system.SetSolver(solver);
            my_system.SetSolverForceTolerance(1e-10);
            break;
        }
        case MKL: {
#ifdef CHRONO_MKL
            GetLog() << "Using MKL solver\n";
            auto mkl_solver = chrono_types::make_shared<ChSolverMKL>();
            mkl_solver->LockSparsityPattern(true);
            my_system.SetSolver(mkl_solver);
#endif
            break;
        }
    }

    // Set up integrator
    switch (integrator_type) {
        case EULER:
            GetLog() << "Using EULER_IMPLICIT_LINEARIZED integrator\n";
            my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
            break;
        case HHT: {
            GetLog() << "Using HHT integrator\n";
            my_system.SetTimestepperType(ChTimestepper::Type::HHT);
            auto integrator = std::static_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
            integrator->SetAlpha(-0.2);
            integrator->SetMaxiters(10);
            integrator->SetAbsTolerances(1e-3, 1e-2);
            integrator->SetVerbose(true);
            break;
        }
    }

    // ------------------
    // Perform simulation
    // ------------------
#ifndef CHRONO_IRRLICHT
    visualization = false;
#endif

    if (visualization) {
#ifdef CHRONO_IRRLICHT
        // Create Irrlicht visualization application
        ChIrrApp application(&my_system, L"ABAQUS tire demo", core::dimension2d<u32>(1280, 720), false, true);

        application.AddTypicalLogo();
        application.AddTypicalSky();
        application.AddTypicalLights();
        application.AddTypicalCamera(core::vector3df(1.0f, 1.4f, -1.2f), core::vector3df(0, (f32)tire_rad, 0));
        application.AddLightWithShadow(core::vector3df(1.5f, 5.5f, -2.5f), core::vector3df(0, 0, 0), 3, 2.2, 7.2, 40, 512,
                                       video::SColorf(0.8f, 0.8f, 1.0f));
        application.AddShadowAll();
        ////application.SetContactsDrawMode(ChIrrTools::CONTACT_DISTANCES);

        // Bind assets
        application.AssetBindAll();
        application.AssetUpdateAll();

        // Simulation loop
        application.SetTimestep(step_size);

        while (application.GetDevice()->run()) {
            application.BeginScene();
            application.DrawAll();
            application.DoStep();
            application.EndScene();
        }
#endif
    } else {
        // Simulation loop
        ChTimer<> timer;
        timer.start();
        for (int istep = 0; istep < num_steps; istep++) {
            my_system.DoStepDynamics(step_size);
        }
        timer.stop();

        // Report run time.
        GetLog() << "Simulation time:  " << timer() << "\n";
        GetLog() << "Internal forces (" << my_mesh->GetNumCallsInternalForces()
                 << "):  " << my_mesh->GetTimeInternalForces() << "\n";
        GetLog() << "Jacobian (" << my_mesh->GetNumCallsJacobianLoad() << "):  " << my_mesh->GetTimeJacobianLoad()
                 << "\n";
        GetLog() << "Extra time:  " << timer() - my_mesh->GetTimeInternalForces() - my_mesh->GetTimeJacobianLoad()
                 << "\n";
    }

    return 0;
}
