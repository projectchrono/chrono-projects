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

#include "chrono_opengl/ChOpenGLWindow.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// Simulation parameters

double gravity = 9.81;
double time_step = 1e-3;
double settling_time = 0.01;// 0.2;
double normal_pressure = 1e3;  // Pa

// Parameters for the ball

int ballId = 1;
double radius = 0.0025;  // m
double density = 2500;   // kg/m^3
double mass = density * (4.0 / 3) * CH_C_PI * radius * radius * radius;
double COR = 0.9;
double mu = 0.5;
ChVector<> ball_pos(2 * radius, 0, 2 * radius);
//ChVector<> ball_pos(0, 0, 0);

// Parameters for the lower and upper plates

int groundId = 0;
int binId = -1;
int plateId = -3;
double w = 0.06;
double l = 0.06;
double h = 0.06;
double thick = 0.01;


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

  // Create a material (will be used by all objects)

  ChSharedPtr<ChMaterialSurface> material;
  material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  material->SetRestitution(COR);
  material->SetFriction(mu);

  // Create lower bin

  ChSharedPtr<ChBody> bin(new ChBody(new ChCollisionModelParallel));
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

  double shear_Area = w * l;

  ChSharedPtr<ChBody> plate(new ChBody(new ChCollisionModelParallel));
  plate->SetIdentifier(binId);
  plate->SetMass(normal_pressure * shear_Area / gravity);
  plate->SetPos(ChVector<>(0, h, 0));
  plate->SetBodyFixed(false);
  plate->SetCollide(true);
  plate->SetMaterialSurface(material);

  plate->GetCollisionModel()->ClearModel();
  AddWall(plate, ChVector<>(w / 2, thick / 2, l / 2), ChVector<>(0, 0, 0), true);
  plate->GetCollisionModel()->SetFamily(2);
  plate->GetCollisionModel()->BuildModel();

  system->AddBody(plate);

  // Create ball

  ChSharedPtr<ChBody> ball(new ChBody(new ChCollisionModelParallel));
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

  // Create the system

  ChSystemParallelDVI* my_system = new ChSystemParallelDVI();
  my_system->GetSettings()->solver.tolerance = 0.01;
  my_system->GetSettings()->solver.max_iteration_bilateral = 100;
  my_system->GetSettings()->solver.clamp_bilaterals = false;
  my_system->GetSettings()->solver.bilateral_clamp_speed = 0.1;

  my_system->GetSettings()->solver.solver_mode = SLIDING;
  my_system->GetSettings()->solver.max_iteration_normal = 0;
  my_system->GetSettings()->solver.max_iteration_sliding = 1000;
  my_system->GetSettings()->solver.max_iteration_spinning = 0;
  my_system->GetSettings()->solver.alpha = 0;
  my_system->GetSettings()->solver.contact_recovery_speed = 0.1;
  my_system->ChangeSolverType(APGDREF);

  my_system->GetSettings()->collision.collision_envelope = 0.05 * radius;
  my_system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

  my_system->Set_G_acc(ChVector<>(0, -gravity, 0));

  // Create the mechanism

  CreateMechanism(my_system);

  // Create the OpenGL visualization

  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "DVI test", my_system);
  gl_window.SetCamera(ChVector<>(3 * w, 0, 0), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0), radius, radius);
  gl_window.SetRenderMode(opengl::WIREFRAME);

  // Simulation loop

  while (gl_window.Active()) {
    gl_window.DoStepDynamics(time_step);
    gl_window.Render();
  }

  return 0;
}
