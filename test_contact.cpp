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
// ChronoParallel program for testing the DVI contact
//
// The global reference frame has Z up.
// All units SI (CGS, i.e., centimeter - gram - second)
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

// Save PovRay post-processing data?
bool write_povray_data = true;

// Output
std::string out_dir = "../TEST_SHEAR";

std::string pov_dir = out_dir + "/POVRAY";
std::string fill_file = out_dir + "/filling.dat";
std::string stats_file = out_dir + "/stats.dat";
std::string settled_ckpnt_file = out_dir + "/settled.dat";

// Frequency for visualization output
int out_fps = 60;

// Frequency for writing results to output file
int write_fps = 1000;

// Simulation frame at which detailed timing information is printed
int timing_frame = -1;

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
  double gravity = -9.81;
  double hh = 0.001;
  double tolerance = 1;
  double envelope = 0.05;
  if (argc > 1)
  {
    tolerance      = atof(argv[1]);
    envelope       = atof(argv[2]);   // [cm]
  }

  // -------------
  // Create system
  // -------------

  cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
  msystem->Set_G_acc(ChVector<>(0, gravity, 0));

  // -------------
  // Create bodies
  // -------------

  // Create a material
  ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
  mmaterial->SetFriction(1.0); // Friction coefficient of steel

  // Create the ground
  ChSharedPtr<ChBody> ground = createBox(msystem, ChVector<>(10,2,2), ChVector<>(0,0,0), QUNIT);
  ground->SetBodyFixed(true);
  ground->SetMaterialSurface(mmaterial);

  // Create a stationary ball on top
  ChSharedPtr<ChBody> ball1 = createBall(msystem, 1, ChVector<>(-4,2,0), QUNIT);
  ball1->SetMaterialSurface(mmaterial);
  ball1->SetMass(1.0);

  // Create a stationary ball on bottom
  ChSharedPtr<ChBody> ball2 = createBall(msystem, 1, ChVector<>(-2,-2,0), QUNIT);
  ball2->SetMaterialSurface(mmaterial);
  ball2->SetMass(1.0);

  // Create a ball with a small force applied in lateral direction
  ChSharedPtr<ChBody> ball3 = createBall(msystem, 1, ChVector<>(0,2,0), QUNIT);
  ball3->SetMaterialSurface(mmaterial);
  ball3->SetMass(1.0);
  ball3->Empty_forces_accumulators();
  ball3->Accumulate_force(ChVector<>(1, 0, 0), ball3->GetPos(), false);

  // Create a ball with a large force applied in lateral direction
  ChSharedPtr<ChBody> ball4 = createBall(msystem, 1, ChVector<>(4,2,0), QUNIT);
  ball4->SetMaterialSurface(mmaterial);
  ball4->SetMass(1.0);
  ball4->Empty_forces_accumulators();
  ball4->Accumulate_force(ChVector<>(100, 0, 0), ball4->GetPos(), false);

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
  gl_window.SetCamera(ChVector<>(0,3,20), ChVector<>(0,3,0),ChVector<>(0,1,0));
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
    printf("Ball1: %f %f %f\n",force.x,force.y,force.z);
    force = msystem->GetBodyContactForce(2);
    printf("Ball2: %f %f %f\n",force.x,force.y,force.z);
    force = msystem->GetBodyContactForce(3);
    printf("Ball3: %f %f %f\n",force.x,force.y,force.z);
    force = msystem->GetBodyContactForce(4);
    printf("Ball4: %f %f %f\n",force.x,force.y,force.z);
    std::cin.get();   
  }

  return 0;
}
