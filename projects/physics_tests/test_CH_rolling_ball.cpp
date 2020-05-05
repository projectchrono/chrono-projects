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

#include "chrono/core/ChStream.h"
#include "chrono/assets/ChBoxShape.h"
#include "chrono/assets/ChSphereShape.h"
#include "chrono/physics/ChSystemNSC.h"

#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

using namespace chrono;

// -----------------------------------------------------------------------------
// Callback functor for contact reporting
// -----------------------------------------------------------------------------
class ContactManager : public ChContactContainer::ReportContactCallback {
public:
    const ChVector<>& GetForce() const { return m_force; }

  private:
    virtual bool OnReportContact(const ChVector<>& pA,
                                 const ChVector<>& pB,
                                 const ChMatrix33<>& plane_coord,
                                 const double& distance,
                                 const double& eff_radius,
                                 const ChVector<>& cforce,
                                 const ChVector<>& ctorque,
                                 ChContactable* modA,
                                 ChContactable* modB) override {
        m_force = cforce;
        return true;
    }

    ChVector<> m_force;
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Output(double time, std::shared_ptr<ChBody> ball, int ncontacts, const ChVector<>& cforce) {
    auto pos = ball->GetPos();
    auto vel = ball->GetPos_dt();
    auto omg = ball->GetWvel_loc();

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
    system.Set_G_acc(ChVector<>(0, 0, -9.81));

    system.SetSolverType(ChSolver::Type::APGD);
    system.SetSolverMaxIterations(1000);
    system.SetSolverTolerance(1e-9);
    system.SetMaxPenetrationRecoverySpeed(0);

    // Create a material (will be used by both objects)
    auto material = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    material->SetRestitution(0);
    material->SetFriction(0.2f);

    // Ground
    auto ground = chrono_types::make_shared<ChBody>();
    system.AddBody(ground);
    ground->SetBodyFixed(true);

    ground->SetCollide(true);
    ground->GetCollisionModel()->ClearModel();
    ground->GetCollisionModel()->AddBox(material, 4, 2, 1, ChVector<>(0, 0, -1));
    ground->GetCollisionModel()->BuildModel();

    auto box = chrono_types::make_shared<ChBoxShape>();
    box->GetBoxGeometry().Size = ChVector<>(4, 2, 1);
    box->GetBoxGeometry().Pos = ChVector<>(0, 0, -0.5);
    box->SetColor(ChColor(1, 0, 0));
    box->SetFading(0.6f);
    ground->AddAsset(box);

    // Crank
    auto ball = chrono_types::make_shared<ChBody>();
    system.AddBody(ball);
    ball->SetMass(1);
    ball->SetInertiaXX(ChVector<>(0.4, 0.4, 0.4));
    ball->SetPos(ChVector<>(0, 0, 1));
    ball->SetPos_dt(ChVector<>(2, 0, 0));

    ball->SetCollide(true);
    ball->GetCollisionModel()->ClearModel();
    ball->GetCollisionModel()->AddSphere(material, 1);
    ball->GetCollisionModel()->BuildModel();

    auto sphere = chrono_types::make_shared<ChSphereShape>();
    sphere->GetSphereGeometry().rad = 1;
    ball->AddAsset(sphere);

    // Functor class for contact reporting
    auto cmanager = chrono_types::make_shared<ContactManager>();

    // Perform the simulation.
    double time = 0;
    double time_step = 1e-3;

    Output(0.0, ball, 0, ChVector<>(0, 0, 0));

    while (time < 0.6) {
        // Advance system state for one step.
        system.DoStepDynamics(time_step);
        time += time_step;

        // Process contacts
        system.GetContactContainer()->ReportAllContacts(cmanager);

        // Print ball states
        Output(time, ball, system.GetNcontacts(), cmanager->GetForce());        
    }

    return 0;
}


