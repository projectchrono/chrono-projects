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

#include <iomanip>

#include "chrono/ChConfig.h"
#include "chrono/physics/ChContactContainerSMC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#ifdef CHRONO_PARDISO_MKL
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#endif

#include <irrlicht.h>

using namespace chrono;
using namespace chrono::irrlicht;

// Custom contact container -- get access to the contact lists in the base class.
class MyContactContainer : public ChContactContainerSMC {
  public:
    MyContactContainer() {}

    // Traverse the list contactlist_6_6
    void ScanContacts(std::shared_ptr<ChBody> ball) {
        auto iter = contactlist_6_6.begin();
        while (iter != contactlist_6_6.end()) {
            ChContactable* objA = (*iter)->GetObjA();
            ChContactable* objB = (*iter)->GetObjB();
            ChVector3d pA = (*iter)->GetContactP1();
            ChVector3d pB = (*iter)->GetContactP2();

            std::cout << std::scientific << std::setw(12) << std::setprecision(3);

            int iball;
            if (objA == ball.get()) {
                iball = 0;
                std::cout << "pA = " << pA << "\n";
            } else if (objB == ball.get()) {
                iball = 1;
                std::cout << "pB = " << pB << "\n";
            }

            const ChKblockGeneric* KRM = (*iter)->GetJacobianKRM();
            const ChMatrixDynamic<double>* K = (*iter)->GetJacobianK();
            const ChMatrixDynamic<double>* R = (*iter)->GetJacobianR();

            if (KRM) {
                std::cout << "K = " << *K << "\n";
                std::cout << "R = " << *R << "\n";

                ChMatrixNM<double, 6, 6> Kball;
                ChMatrixNM<double, 6, 6> Rball;
                Kball = K->block(iball * 6, iball * 6, 6, 6);
                Rball = R->block(iball * 6, iball * 6, 6, 6);

                std::cout << "Kball = " << Kball << "\n";
                std::cout << "Rball = " << Rball << "\n";
            }

            ++iter;
        }
    }
};

// ====================================================================================

int main(int argc, char* argv[]) {
    // ---------------------------------
    // Set path to Chrono data directory
    // ---------------------------------
    SetChronoDataPath(CHRONO_DATA_DIR);

    // ---------------------
    // Simulation parameters
    // ---------------------

    double gravity = -9.81;   // gravitational acceleration
    double time_step = 1e-4;  // integration step size

    ChSolver::Type solver_type = ChSolver::Type::PARDISO_MKL;

    bool stiff_contact = true;

    // ---------------------------
    // Contact material properties
    // ---------------------------

    bool use_mat_properties = false;
    ChSystemSMC::ContactForceModel force_model = ChSystemSMC::Hooke;
    ChSystemSMC::TangentialDisplacementModel tdispl_model = ChSystemSMC::None;

    float young_modulus = 2e9f;
    float friction = 0.4f;
    float restitution = 0.1f;
    float adhesion = 0.0f;

    float kn = 1e8;
    float gn = 0;
    float kt = 0;
    float gt = 0;

    // -------------------------------
    // Parameters for the falling ball
    // -------------------------------

    int ballId = 100;
    double radius = 1;
    double mass = 1000;
    ChVector3d pos(0, 2, 0);
    ChQuaternion<> rot(1, 0, 0, 0);
    ChVector3d init_vel(0, 0, 0);
    ChVector3d init_omg(0, 0, 0);

    // ---------------------------------
    // Parameters for the containing bin
    // ---------------------------------

    int binId = 200;
    double width = 4;
    double length = 4;
    double thickness = 0.2;

    // -----------------
    // Create the system
    // -----------------

    ChSystemSMC system(use_mat_properties);
    system.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Set the SMC contact force model
    system.SetContactForceModel(force_model);

    // Set tangential displacement model
    system.SetTangentialDisplacementModel(tdispl_model);

    // Set contact forces as stiff (to force Jacobian computation) or non-stiff
    system.SetStiffContact(stiff_contact);

    system.Set_G_acc(ChVector3d(0, gravity, 0));

    // Create a material (will be used by both objects)
    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetYoungModulus(young_modulus);
    material->SetRestitution(restitution);
    material->SetFriction(friction);
    material->SetAdhesion(adhesion);
    material->SetKn(kn);
    material->SetGn(gn);
    material->SetKt(kt);
    material->SetGt(gt);

    // Create the falling ball
    auto ball = chrono_types::make_shared<ChBody>();

    ball->SetIdentifier(ballId);
    ball->SetMass(mass);
    ball->SetInertiaXX(0.4 * mass * radius * radius * ChVector3d(1, 1, 1));
    ball->SetPos(pos);
    ball->SetRot(rot);
    ball->SetPos_dt(init_vel);
    ball->SetWvel_par(init_omg);
    ball->SetCollide(true);
    ball->SetBodyFixed(false);

    auto ball_ct_shape = chrono_types::make_shared<ChCollisionShapeSphere>(material, radius);
    ball->AddCollisionShape(ball_ct_shape);

    auto sphere = chrono_types::make_shared<ChVisualShapeSphere>(radius);
    sphere->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
    ball->AddVisualShape(sphere);

    system.AddBody(ball);

    // Create ground
    auto ground = chrono_types::make_shared<ChBody>();

    ground->SetIdentifier(binId);
    ground->SetMass(1);
    ground->SetPos(ChVector3d(0, 0, 0));
    ground->SetRot(ChQuaternion<>(1, 0, 0, 0));
    ground->SetCollide(true);
    ground->SetBodyFixed(true);

    auto ground_ct_shape = chrono_types::make_shared<ChCollisionShapeBox>(material, width, thickness, length);
    ground->AddCollisionShape(ground_ct_shape, ChFrame<>(ChVector3d(0, -thickness / 2, 0), QUNIT));

    auto box = chrono_types::make_shared<ChVisualShapeBox>(width, thickness, length);
    ground->AddVisualShape(box, ChFrame<>(ChVector3d(0, -thickness / 2, 0)));

    system.AddBody(ground);

    // Create the Irrlicht visualization
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("SMC demo");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector3d(0, 3, -6));
    vis->AddTypicalLights();
    vis->SetSymbolScale(1e-4);
    vis->EnableContactDrawing(ContactsDrawMode::CONTACT_FORCES);
    vis->AttachSystem(&system);

    // ----------------------------
    // Use custom contact container
    // ----------------------------

    auto container = chrono_types::make_shared<MyContactContainer>();
    system.SetContactContainer(container);

    // -------------------
    // Setup linear solver
    // -------------------

    // Note that not all solvers support stiffness matrices (that includes the default SolverSMC).
#ifndef CHRONO_PARDISO_MKL
    if (solver_type == ChSolver::Type::PARDISO_MKL) {
        std::cout << "PardisoMKL support not enabled.  Solver reset to default.\n";
        solver_type = ChSolver::Type::MINRES;
    }
#endif

    switch (solver_type) {
        case ChSolver::Type::MINRES: {
            std::cout << "Using MINRES solver.\n";
            auto minres_solver = chrono_types::make_shared<ChSolverMINRES>();
            minres_solver->EnableDiagonalPreconditioner(true);
            minres_solver->SetMaxIterations(100);
            system.SetSolver(minres_solver);
            system.SetSolverForceTolerance(1e-6);
            break;
        }
        case ChSolver::Type::PARDISO_MKL: {
#ifdef CHRONO_PARDISO_MKL
            std::cout << "Using PardisoMKL solver.\n";
            auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
            mkl_solver->LockSparsityPattern(true);
            system.SetSolver(mkl_solver);
#endif
            break;
        }
        default: {
            std::cout << "Using DEFAULT solver.\n";
            system.SetSolverMaxIterations(100);
            system.SetSolverForceTolerance(1e-6);
            break;
        }
    }

    // ----------------
    // Setup integrator
    // ----------------

    system.SetTimestepperType(ChTimestepper::Type::HHT);
    auto integrator = std::static_pointer_cast<ChTimestepperHHT>(system.GetTimestepper());
    integrator->SetAlpha(0.0);
    integrator->SetMaxiters(100);
    integrator->SetAbsTolerances(1e-08);
    ////integrator->SetStepControl(false);
    ////integrator->SetVerbose(true);

    // ---------------
    // Simulation loop
    // ---------------
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        system.DoStepDynamics(time_step);
        if (ball->GetPos().y() <= radius) {
            container->ScanContacts(ball);
            std::cout << "t = " << system.GetChTime() << "  NR iters. = " << integrator->GetNumIterations() << "\n";
        }
    }

    return 0;
}
