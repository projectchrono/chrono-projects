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
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

// Control use of OpenGL run-time rendering
//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// Simulation parameters

double UNIT_SCALE = 1000;     // 1 for m, 1000 for mm

double gravity = 9.81 * UNIT_SCALE;
double time_step = 1e-3;
double settling_time = 0.01;   // 0.2;
double normal_pressure = 1e3 / UNIT_SCALE;

// Parameters for the ball

int ballId = 1;
double radius = 0.0025 * UNIT_SCALE;
double density = 2500 / (UNIT_SCALE * UNIT_SCALE * UNIT_SCALE);
double mass = density * (4.0 / 3) * CH_C_PI * radius * radius * radius;
double COR = 0.9;
double mu = 0.5;
ChVector<> ball_pos(2 * radius, 0, 2 * radius);
//ChVector<> ball_pos(0, 0, 0);

// Parameters for the lower and upper plates

int binId = -1;
int plateId = -3;
double w = 0.06 * UNIT_SCALE;
double l = 0.06 * UNIT_SCALE;
double h = 0.06 * UNIT_SCALE;
double thick = 0.01 * UNIT_SCALE;


// =============================================================================
// Utility for adding (visible or invisible) walls
// =============================================================================
void AddWall(ChSharedBodyPtr& body, const ChVector<>& dim, const ChVector<>& loc, bool visible) {
  body->GetCollisionModel()->AddBox(dim.x, dim.y, dim.z, loc);

  if (visible == true) {
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
// Create mechanism bodies and links
// =============================================================================
void CreateMechanism(ChSystem* system) {

  // Get the type of system.
  utils::SystemType sys_type = utils::GetSystemType(system);

  // Create a material (will be used by all objects)

  ChSharedPtr<ChMaterialSurface> material;
  material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  material->SetRestitution(COR);
  material->SetFriction(mu);

  // Create lower bin

  ChSharedPtr<ChBody> bin;
  switch (sys_type) {
    case utils::PARALLEL_DVI:
      bin = ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel));
      break;
    case utils::SEQUENTIAL_DVI:
      bin = ChSharedPtr<ChBody>(new ChBody());
      break;
  }

  bin->SetIdentifier(binId);
  bin->SetMass(1);
  bin->SetPos(ChVector<>(0, -h / 2, 0));
  bin->SetBodyFixed(true);
  bin->SetCollide(true);
  bin->SetMaterialSurface(material);

  bin->GetCollisionModel()->ClearModel();
  AddWall(bin, ChVector<>(w / 2, thick / 2, l / 2), ChVector<>(0, 0, 0), true);
  bin->GetCollisionModel()->SetFamily(1);
  bin->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);
  bin->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
  bin->GetCollisionModel()->SetFamilyMaskDoCollisionWithFamily(4);
  bin->GetCollisionModel()->BuildModel();

  system->AddBody(bin);

  // Create upper load plate

  ChSharedPtr<ChBody> plate;
  switch (sys_type) {
  case utils::PARALLEL_DVI:
    plate = ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel));
    break;
  case utils::SEQUENTIAL_DVI:
    plate = ChSharedPtr<ChBody>(new ChBody());
    break;
  }

  double plate_mass = normal_pressure * (w * l) / gravity;
  ChVector<> plate_hdim(w / 2, thick / 2, l / 2);

  plate->SetIdentifier(plateId);
  plate->SetMass(plate_mass);
  //plate->SetInertia(utils::CalcBoxGyration(plate_hdim) * plate_mass);
  plate->SetPos(ChVector<>(0, h, 0));
  plate->SetBodyFixed(false);
  plate->SetCollide(true);
  plate->SetMaterialSurface(material);

  plate->GetCollisionModel()->ClearModel();
  AddWall(plate, plate_hdim, ChVector<>(0, 0, 0), true);
  plate->GetCollisionModel()->SetFamily(2);
  plate->GetCollisionModel()->BuildModel();

  system->AddBody(plate);

  // Create ball

  ChSharedPtr<ChBody> ball;
  switch (sys_type) {
  case utils::PARALLEL_DVI:
    ball = ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel));
    break;
  case utils::SEQUENTIAL_DVI:
    ball = ChSharedPtr<ChBody>(new ChBody());
    break;
  }

  ball->SetIdentifier(ballId);
  ball->SetMass(mass);
  ball->SetInertiaXX((2.0 / 5.0) * mass * radius * radius * ChVector<>(1, 1, 1));
  ball->SetPos(ball_pos);
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

  // Create prismatic (translational) joint between plates.
  // The translational axis of a prismatic joint is along the Z axis of the
  // specified joint coordinate system.  Here, we apply the 'z2y' rotation to
  // align it with the Y axis of the global reference frame.

  ChQuaternion<> z2y;
  z2y.Q_from_AngX(-CH_C_PI / 2);

  ChSharedPtr<ChLinkLockPrismatic> prismatic(new ChLinkLockPrismatic);
  prismatic->Initialize(plate, bin, ChCoordsys<>(ChVector<>(0, 0, 0), z2y));

  system->AddLink(prismatic);
}

// =============================================================================
// Main program
// =============================================================================
int main(int argc, char* argv[]) {

  // Solver parameters

  double tolerance = 0.01;
  double contact_recovery_speed = 0.1;
  int max_iteration_bilateral = 100;
  int max_iteration_sliding = 1000;

  // Create the serial system

  ChSystem* systemS = new ChSystem();
  systemS->Set_G_acc(ChVector<>(0, -gravity, 0));
  systemS->SetTolForce(tolerance);
  systemS->SetIterLCPmaxItersSpeed(max_iteration_sliding);
  systemS->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  systemS->SetUseSleeping(false);
  //systemS->SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
  systemS->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SYMMSOR);
  systemS->SetIntegrationType(ChSystem::INT_ANITESCU);
  systemS->SetIterLCPwarmStarting(true);

  // Create the parallel system

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
  systemP->ChangeSolverType(APGDREF);

  systemP->GetSettings()->collision.collision_envelope = 0.05 * radius;
  systemP->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

  systemP->Set_G_acc(ChVector<>(0, -gravity, 0));

  // Select which system is simulated

  ChSystem* system = systemP;

  // Create the mechanism

  CreateMechanism(system);

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  // Create the OpenGL visualization

  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "DVI test", system);
  gl_window.SetCamera(ChVector<>(3 * w, 0, 0), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0), radius, radius);
  gl_window.SetRenderMode(opengl::WIREFRAME);

  // Simulation loop

  while (gl_window.Active()) {
    gl_window.DoStepDynamics(time_step);
    gl_window.Render();
  }
#else
  systemP->DoStepDynamics(time_step);
#endif

  return 0;
}
