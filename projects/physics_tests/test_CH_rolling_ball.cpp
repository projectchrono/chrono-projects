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
// Chrono demo for slinding-rolling transition.
//
// The model simulated here is a rough uniform sphere of unit radius and mass in
// contact with a fixed horizontal plane in a uniform gravitational field.
//
// The ball is initially in contact with the plane and has an initial velocity
// in the x direction of 2 m/s (no angular velocity).
//
// The coefficient of friction is 0.2
//
// The global reference frame has Z up and gravitational acceleration is g=9.81
// =============================================================================

#include <cstdio>
#include <cmath>

#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/assets/ChVisualShapeSphere.h"
#include "chrono/physics/ChSystemNSC.h"

#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

using namespace chrono;

// -----------------------------------------------------------------------------
// Callback functor for contact reporting
// -----------------------------------------------------------------------------
class ContactManager : public ChContactContainer::ReportContactCallback {
public:
    const ChVector3d& GetForce() const { return m_force; }

  private:
    virtual bool OnReportContact(const ChVector3d& pA,
                                 const ChVector3d& pB,
                                 const ChMatrix33<>& plane_coord,
                                 const double& distance,
                                 const double& eff_radius,
                                 const ChVector3d& cforce,
                                 const ChVector3d& ctorque,
                                 ChContactable* modA,
                                 ChContactable* modB) override {
        m_force = cforce;
        return true;
    }

    ChVector3d m_force;
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Output(double time, std::shared_ptr<ChBody> ball, int ncontacts, const ChVector3d& cforce) {
    auto pos = ball->GetPos();
    auto vel = ball->GetLinVel();
    auto omg = ball->GetAngVelLocal();

    std::cout << time << "  ";
    std::cout << ncontacts << "     ";
    std::cout << pos.x() << "  " << pos.y() << "  " << pos.z() << "     ";
    std::cout << vel.x() << " " << vel.y() << "  " << vel.z() << "     ";
    std::cout << omg.x() << "  " << omg.y() << "  " << omg.z() << "     ";
    std::cout << cforce.x() << "  " << cforce.y() << "  " << cforce.z();
    std::cout << std::endl;
}

// -----------------------------------------------------------------------------
// Create the mechanical system, define bodies and joints, and perform
// simulation.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Create the physical system and set gravity along negative Z
    ChSystemNSC system;
    system.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    system.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));

    system.SetSolverType(ChSolver::Type::APGD);
    system.GetSolver()->AsIterative()->SetMaxIterations(1000);
    system.GetSolver()->AsIterative()->SetTolerance(1e-12);
    system.SetMaxPenetrationRecoverySpeed(0);

    // Create a material (will be used by both objects)
    auto material = chrono_types::make_shared<ChContactMaterialNSC>();
    material->SetRestitution(0);
    material->SetFriction(0.2f);

    // Ground
    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetFixed(true);
    ground->EnableCollision(true);

    auto ground_ct_shape = chrono_types::make_shared<ChCollisionShapeBox>(material, 4, 2, 1);
    ground->AddCollisionShape(ground_ct_shape, ChFrame<>(ChVector3d(0, 0, -1), QUNIT));

    auto box = chrono_types::make_shared<ChVisualShapeBox>(8, 4, 2);
    box->SetColor(ChColor(1, 0, 0));
    ground->AddVisualShape(box, ChFrame<>(ChVector3d(0, 0, -0.5)));

    // Crank
    auto ball = chrono_types::make_shared<ChBody>();
    system.AddBody(ball);
    ball->SetMass(1);
    ball->SetInertiaXX(ChVector3d(0.4, 0.4, 0.4));
    ball->SetPos(ChVector3d(0, 0, 1));
    ball->SetLinVel(ChVector3d(2, 0, 0));
    ball->EnableCollision(true);

    auto ball_ct_shape = chrono_types::make_shared<ChCollisionShapeSphere>(material, 1);
    ball->AddCollisionShape(ball_ct_shape);

    auto sphere = chrono_types::make_shared<ChVisualShapeSphere>(1);
    ball->AddVisualShape(sphere);

    // Functor class for contact reporting
    auto cmanager = chrono_types::make_shared<ContactManager>();

    // Perform the simulation.
    double time = 0;
    double time_step = 1e-3;

    Output(0.0, ball, 0, ChVector3d(0, 0, 0));

    while (time < 0.6) {
        // Advance system state for one step.
        system.DoStepDynamics(time_step);
        time += time_step;

        // Process contacts
        system.GetContactContainer()->ReportAllContacts(cmanager);

        // Print ball states
        Output(time, ball, system.GetNumContacts(), cmanager->GetForce());        
    }

    return 0;
}


