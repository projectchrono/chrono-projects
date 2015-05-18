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
// Author: Daniel Melanz
// =============================================================================
//
// ChronoParallel program for testing the DVI contact, there are 4 cases that
// are being tested (coefficient of friction is 1.0):
//   1) Ball sitting on plane, no gravity
//   2) Ball sitting on plane, with gravity - no applied force in lateral dir
//   3) Ball sitting on plane, with gravity - small applied force in lateral dir
//   4) Ball sitting on plane, with gravity - large applied force in lateral dir
//
// The contact force is measured in each case, the results should be as follows:
//   1) Normal Force: Zero, Tangent Force: Zero
//   2) Normal Force: Weight, Tangent Force: Zero
//   3) Normal Force: Weight, Tangent Force: Equal to applied force (mu = 1)
//   4) Normal Force: Weight, Tangent Force: Weight (mu = 1)
//
// The global reference frame has Z up.
// All units SI.
//
// =============================================================================

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

// Control use of OpenGL run-time rendering
//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "demo_utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

ChSharedPtr<ChBody> createBall(ChSystem* mphysicalSystem, double radius, ChVector<> position, ChQuaternion<> rotation)
{
  // create the body
  ChSharedPtr<ChBody> ball(new ChBody(new ChCollisionModelParallel));
  ball->SetPos(position);
  ball->SetRot(rotation);
  mphysicalSystem->Add(ball);

  // Define the collision shape
  ball->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball.get_ptr(), radius, ChVector<>(0,0,0));
  ball->GetCollisionModel()->BuildModel();
  ball->SetCollide(true);

  return ball;
}

ChSharedPtr<ChBody> createBox(ChSystem* mphysicalSystem, ChVector<> size, ChVector<> position, ChQuaternion<> rotation)
{
  double L = size.x;
  double H = size.y;
  double W = size.z;

  // create the body
  ChSharedPtr<ChBody> box(new ChBody(new ChCollisionModelParallel));
  box->SetPos(position);
  box->SetRot(rotation);
  mphysicalSystem->Add(box);

  // Define the collision shape
  box->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(box.get_ptr(), ChVector<>(L*.5,H*.5,W*.5), ChVector<>(0,0,0));
  box->GetCollisionModel()->BuildModel();
  box->SetCollide(true);

  return box;
}
// =============================================================================

int main(int argc, char* argv[])
{
  // working in [m, kg, s]
  ChVector<> gravity = ChVector<>(0,0,-9.81);
  ChVector<> appliedForce = ChVector<>(100,0,0);
  double hh = 0.001;
  double tolerance = 0.001;
  double envelope = 0.005;
  double mass = 1.0;
  double mu_sliding = 1.0;
  int testCase = 4;

  // Choose a test case
  if (argc > 1)
  {
    testCase = atoi(argv[1]);
    switch (testCase) {
      case 1:
    	// 1) Ball sitting on plane, no gravity
    	gravity = ChVector<>(0,0,0);
    	appliedForce = ChVector<>(0,0,0);
        break;
      case 2:
    	// 2) Ball sitting on plane, with gravity - no applied force in lateral dir
      	gravity = ChVector<>(0,0,-9.81);
      	appliedForce = ChVector<>(0,0,0);
        break;
      case 3:
    	// 3) Ball sitting on plane, with gravity - small applied force in lateral dir
      	gravity = ChVector<>(0,0,-9.81);
      	appliedForce = ChVector<>(1,0,0);
        break;
      case 4:
    	// 4) Ball sitting on plane, with gravity - large applied force in lateral dir
      	gravity = ChVector<>(0,0,-9.81);
      	appliedForce = ChVector<>(100,0,0);
        break;
      }
  }

  // -------------
  // Create system
  // -------------

  cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
  msystem->Set_G_acc(gravity);

  // -------------
  // Create bodies
  // -------------

  // Create a material
  ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
  mmaterial->SetFriction(mu_sliding); // Friction coefficient of steel

  // Create the ground
  ChSharedPtr<ChBody> ground = createBox(msystem, ChVector<>(2,2,2), ChVector<>(0,0,-1), QUNIT);
  ground->SetBodyFixed(true);
  ground->SetMaterialSurface(mmaterial);

  // Create a unit ball
  ChSharedPtr<ChBody> ball = createBall(msystem, 1.0, ChVector<>(0,0,1), QUNIT);
  ball->SetMaterialSurface(mmaterial);
  ball->SetMass(mass);
  ball->Empty_forces_accumulators();
  ball->Accumulate_force(appliedForce, ball->GetPos(), false);

  // Set number of threads.
  int max_threads = msystem->GetParallelThreadNumber();
  if (threads > max_threads) threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;
  msystem->GetSettings()->max_threads = threads;
  msystem->GetSettings()->perform_thread_tuning = thread_tuning;

  // Edit system settings
  msystem->GetSettings()->solver.tolerance = tolerance;
  msystem->GetSettings()->solver.max_iteration_bilateral = 0;
  msystem->GetSettings()->solver.clamp_bilaterals = false;
  msystem->GetSettings()->solver.bilateral_clamp_speed = 10e30;
  msystem->GetSettings()->solver.solver_mode = SLIDING;
  msystem->GetSettings()->solver.max_iteration_normal = 0;
  msystem->GetSettings()->solver.max_iteration_sliding = 50000;
  msystem->GetSettings()->solver.max_iteration_spinning = 0;
  msystem->GetSettings()->solver.alpha = 0;
  msystem->GetSettings()->solver.contact_recovery_speed = 10e30;
  msystem->SetMaxPenetrationRecoverySpeed(10e30);
  msystem->ChangeSolverType(APGD);

  msystem->GetSettings()->collision.collision_envelope = envelope;
  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

  // ----------------------
  // Perform the simulation
  // ----------------------

  // Initialize counters
  int currentStage = 1;

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow &gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "Contact Test", msystem);
  gl_window.SetCamera(ChVector<>(0,10,0), ChVector<>(0,0,0),ChVector<>(0,0,1));
  gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

  // Loop until reaching the end time...
  while (currentStage<=3) {
    // Advance simulation by one step
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(hh);
      gl_window.Render();
    }
#else
    msystem->DoStepDynamics(hh);
#endif

    std::vector<double> history = ((ChLcpIterativeSolver*) (msystem->GetLcpSolverSpeed()))->GetViolationHistory();
    int numIters = history.size();
    printf("Iterations: %d\n",numIters);

    msystem->CalculateContactForces();
    real3 force(0, 0, 0);
    force = msystem->GetBodyContactForce(1);
    printf("Contact force on ball: (%f, %f, %f)\n",force.x,force.y,force.z);
    switch (testCase) {
      case 1:
    	// 1) Ball sitting on plane, no gravity
    	printf("The contact force should be (0, 0, 0)\n");
        break;
      case 2:
    	// 2) Ball sitting on plane, with gravity - no applied force in lateral dir
    	printf("The contact force should be (0, 0, %f)\n",-mass*gravity.z);
        break;
      case 3:
    	// 3) Ball sitting on plane, with gravity - small applied force in lateral dir
    	printf("The contact force should be (%f, 0, %f)\n",-appliedForce.x,-mass*gravity.z);
        break;
      case 4:
    	// 4) Ball sitting on plane, with gravity - large applied force in lateral dir
    	printf("The contact force should be (%f, 0, %f)\n",mass*gravity.z*mu_sliding,-mass*gravity.z);
        break;
      }

    std::cin.get();   
  }

  return 0;
}
