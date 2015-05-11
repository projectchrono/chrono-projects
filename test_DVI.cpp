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
// Author: Radu Serban
// =============================================================================
//
// Basic test program for DVI solver.  The test consist of two plates, the first
// fixed to ground and the second connected through a translational joint to the
// first one.  A single ball is placed between the two plates.
//
// The global reference frame has Y up.
//
// =============================================================================

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "demo_utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// System type parallel?

bool parallel = true;
SOLVERTYPE solver_type = APGD;

// Simulation parameters

double UNIT_SCALE = 1;     // 1 for m, 1000 for mm

double gravity = 10 * UNIT_SCALE;
double normal_pressure = 1e3 / UNIT_SCALE;

bool shear = true;   // move bottom plate?
double shear_speed = 0.05 * UNIT_SCALE;

// Solver settings

double time_step = 1e-4;
double tolerance = 0.01 * UNIT_SCALE;
double contact_recovery_speed = 1e30;
int max_iteration_bilateral = 0;
int max_iteration_sliding = 10000;

// Parameters for the ball

double radius = 0.0025 * UNIT_SCALE;
double density = 2500 / (UNIT_SCALE * UNIT_SCALE * UNIT_SCALE);
double mass = density * (4.0 / 3) * CH_C_PI * radius * radius * radius;
float COR = 0.9f;
float mu = 0.5f;

// Parameters for the lower and upper plates

double w = 0.06 * UNIT_SCALE;
double l = 0.06 * UNIT_SCALE;
double h = 0.06 * UNIT_SCALE;
double thick = 0.01 * UNIT_SCALE;


// =============================================================================
// Utility for adding (visible or invisible) walls
// =============================================================================
void AddWall(ChSharedPtr<ChBody> body, const ChVector<>& dim, const ChVector<>& loc, bool visible) {
  body->GetCollisionModel()->AddBox(dim.x, dim.y, dim.z, loc);

  if (visible) {
    ChSharedPtr<ChBoxShape> box(new ChBoxShape);
    box->GetBoxGeometry().Size = dim;
    //    box->GetBoxGeometry().Pos = loc;		// for Chrono
    box->Pos = loc;  // for Chrono-Parallel
    box->SetColor(ChColor(1, 0, 0));
    box->SetFading(0.6f);
    body->AddAsset(box);
  }
}

// =============================================================================
// Utility for adding a visualization cylinder
// =============================================================================
void AddVisCylinder(ChSharedPtr<ChBody> body, double radius, double hlen, const ChVector<>& loc) {
  ChSharedPtr<ChCylinderShape> cyl(new ChCylinderShape);
  cyl->GetCylinderGeometry().rad = radius;
  cyl->GetCylinderGeometry().p1 = ChVector<>(0, -hlen, 0);
  cyl->GetCylinderGeometry().p2 = ChVector<>(0,  hlen, 0);
  cyl->Pos = loc;
  cyl->SetColor(ChColor(1, 0, 0));
  body->AddAsset(cyl);
}

// =============================================================================
// Main program
// =============================================================================
int main(int argc, char* argv[]) {

#ifndef CHRONO_PARALLEL_HAS_OPENGL
  // Do nothing if OpenGL not available
  return 1;
#endif

  // -----------------
  // Create the system
  // -----------------

  ChSystem* system;

  if (parallel) {
    ChSystemParallelDVI* systemP = new ChSystemParallelDVI();

    systemP->GetSettings()->solver.tolerance = tolerance;
    systemP->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    systemP->GetSettings()->solver.clamp_bilaterals = false;
    systemP->GetSettings()->solver.bilateral_clamp_speed = 0.1;

    systemP->GetSettings()->solver.solver_mode = SLIDING;
    systemP->GetSettings()->solver.max_iteration_normal = 0;
    systemP->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    systemP->GetSettings()->solver.max_iteration_spinning = 0;
    systemP->GetSettings()->solver.alpha = 0;
    systemP->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    systemP->ChangeSolverType(solver_type);

    systemP->GetSettings()->collision.collision_envelope = 0.05 * radius;
    systemP->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

    system = systemP;

  } else {
    ChSystem* systemS = new ChSystem();

    systemS->SetTolForce(tolerance);
    systemS->SetIterLCPmaxItersSpeed(max_iteration_sliding);
    systemS->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
    systemS->SetUseSleeping(false);
    systemS->SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
    //systemS->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SYMMSOR);
    systemS->SetIntegrationType(ChSystem::INT_ANITESCU);
    systemS->SetIterLCPwarmStarting(true);

    system = systemS;
  }

  system->Set_G_acc(ChVector<>(0, -gravity, 0));

  // --------------------
  // Create the mechanism
  // --------------------

  ChSharedPtr<ChBody> ground;
  ChSharedPtr<ChBody> bin;
  ChSharedPtr<ChBody> plate;
  ChSharedPtr<ChBody> ball;

  // Create a material (will be used by all objects)

  ChSharedPtr<ChMaterialSurface> material;
  material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  material->SetRestitution(COR);
  material->SetFriction(mu);

  // Create ground body

  ground = parallel ?
    ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel)) :
    ChSharedPtr<ChBody>(new ChBody());

  ground->SetMass(1);
  ground->SetPos(ChVector<>(0, 0, 0));
  ground->SetBodyFixed(true);
  ground->SetCollide(false);

  system->AddBody(ground);
  unsigned int id_ground = ground->GetId();
  cout << "Create ground body, Id = " << id_ground << endl;

  // Create lower bin

  bin = parallel ?
    ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel)) :
    ChSharedPtr<ChBody>(new ChBody());

  bin->SetMass(1);
  bin->SetPos(ChVector<>(0, -thick / 2, 0));
  bin->SetBodyFixed(true);
  bin->SetCollide(true);
  bin->SetMaterialSurface(material);

  bin->GetCollisionModel()->ClearModel();
  AddWall(bin, ChVector<>(w / 2, thick / 2, l / 2), ChVector<>(0, 0, 0), true);
  bin->GetCollisionModel()->SetFamily(1);
  bin->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);
  bin->GetCollisionModel()->BuildModel();

  AddVisCylinder(bin, radius / 4, radius / 2, ChVector<>(w / 4, thick / 2, l / 4));

  system->AddBody(bin);
  unsigned int id_bin = bin->GetId();
  cout << "Create bin body, Id = " << id_bin << endl;

  // Create upper load plate

  plate = parallel ?
    ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel)) :
    ChSharedPtr<ChBody>(new ChBody());

  double plate_mass = normal_pressure * (w * l) / gravity;

  plate->SetMass(plate_mass);
  //plate->SetInertia(utils::CalcBoxGyration(plate_hdim) * plate_mass);
  plate->SetPos(ChVector<>(0, 2 * radius + thick / 2, 0));
  plate->SetBodyFixed(false);
  plate->SetCollide(true);
  plate->SetMaterialSurface(material);

  plate->GetCollisionModel()->ClearModel();
  AddWall(plate, ChVector<>(w / 2, thick / 2, l / 2), ChVector<>(0, 0, 0), true);
  plate->GetCollisionModel()->SetFamily(2);
  plate->GetCollisionModel()->BuildModel();

  AddVisCylinder(plate, radius / 4, radius / 2, ChVector<>(w / 4, -thick / 2, l / 4));

  system->AddBody(plate);
  unsigned int id_plate = plate->GetId();
  cout << "Create plate body, Id = " << id_plate << endl;
  cout << "                   mass = " << plate_mass << endl;
  cout << "                   weight = " << plate_mass * gravity << endl;

  // Create ball

  ball = parallel ?
    ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel)) :
    ChSharedPtr<ChBody>(new ChBody());

  ball->SetMass(mass);
  ball->SetInertiaXX((2.0 / 5.0) * mass * radius * radius * ChVector<>(1, 1, 1));
  ball->SetPos(ChVector<>(w / 4, radius, l / 4));
  ball->SetBodyFixed(false);
  ball->SetCollide(true);

  ball->SetMaterialSurface(material);

  ball->GetCollisionModel()->ClearModel();
  ball->GetCollisionModel()->AddSphere(radius);
  ball->GetCollisionModel()->BuildModel();

  ChSharedPtr<ChSphereShape> sphere(new ChSphereShape);
  sphere->GetSphereGeometry().rad = radius;
  sphere->SetColor(ChColor(1, 0, 1));
  ball->AddAsset(sphere);

  system->AddBody(ball);
  unsigned int id_ball = ball->GetId();
  cout << "Create ball body, Id = " << id_ball << endl;
  cout << "                  mass = " << mass << endl;
  cout << "                  weight = " << mass * gravity << endl;

  // Create prismatic (translational) joint between ground and upper plate.
  // The translational axis of a prismatic joint is along the Z axis of the
  // specified joint coordinate system.  Here, we apply the 'z2y' rotation to
  // align it with the Y axis of the global reference frame.

  ChQuaternion<> z2y;
  z2y.Q_from_AngX(-CH_C_PI / 2);

  ChSharedPtr<ChLinkLockPrismatic> prismatic(new ChLinkLockPrismatic);
  prismatic->Initialize(plate, ground, ChCoordsys<>(ChVector<>(0, 0, 0), z2y));

  system->AddLink(prismatic);

  // -------------------------------
  // Create the OpenGL visualization
  // -------------------------------

  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "DVI test", system);
  gl_window.SetCamera(ChVector<>(1.5 * w, 0, 0), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0), radius, radius);
  gl_window.SetRenderMode(opengl::WIREFRAME);

  // ---------------
  // Simulation loop
  // ---------------

  int frame = 0;

  while (gl_window.Active()) {

    if (shear) {
      double displ = system->GetChTime() * shear_speed;
      bin->SetPos(ChVector<>(0, -thick / 2, displ));
      bin->SetPos_dt(ChVector<>(0, 0, shear_speed));
      bin->SetRot(QUNIT);
    }

    gl_window.DoStepDynamics(time_step);
    gl_window.Render();

    if (system->GetNcontacts() < 2)
      break;

    TimingOutput(system);

    if (parallel && frame % 10 == 0) {
      ChSystemParallelDVI* systemP = dynamic_cast<ChSystemParallelDVI*>(system);
      real3 force;
      systemP->CalculateContactForces();
      cout << endl;
      force = systemP->GetBodyContactForce(id_bin);
      cout << force.x << "  " << force.y << "  " << force.z << endl;
      force = systemP->GetBodyContactForce(id_plate);
      cout << force.x << "  " << force.y << "  " << force.z << endl;
      force = systemP->GetBodyContactForce(id_ball);
      cout << force.x << "  " << force.y << "  " << force.z << endl;
      cout << endl;
    }

    frame++;
  }

  return 0;
}
