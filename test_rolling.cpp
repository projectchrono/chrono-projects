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
// ChronoParallel program for testing the DVI contact, this case tests a ball
// that is initially given a translational velocity in the lateral direction.
// The transition from sliding to rolling is tested in this example.
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

ChSharedPtr<ChBody> createBallSerial(ChSystem* mphysicalSystem, double radius, ChVector<> position, ChQuaternion<> rotation)
{
  // create the body
  ChSharedPtr<ChBody> ball(new ChBody());
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

ChSharedPtr<ChBody> createBoxSerial(ChSystem* mphysicalSystem, ChVector<> size, ChVector<> position, ChQuaternion<> rotation)
{
  double L = size.x;
  double H = size.y;
  double W = size.z;

  // create the body
  ChSharedPtr<ChBody> box(new ChBody());
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
  double hh = 0.001;
  double tolerance = 0.001;
  double envelope = 0.05;
  double mass = 1.0;
  double mu_sliding = 0.2;
  double radius = 1.0;
  double initialVel = 2.0;
  double simTime = 0.6;

  // -------------
  // Create DVI Parallel system
  // -------------

  // Set up the serial system
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
  msystem->Set_G_acc(gravity);

  // Set number of threads
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

  // Create a material
  ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
  mmaterial->SetFriction(mu_sliding); // Friction coefficient of steel

  // Create the ground
  ChSharedPtr<ChBody> ground = createBox(msystem, ChVector<>(10,2,2), ChVector<>(0,0,-1), QUNIT);
  ground->SetBodyFixed(true);
  ground->SetMaterialSurface(mmaterial);

  // Create a unit ball
  ChSharedPtr<ChBody> ball = createBall(msystem, radius, ChVector<>(0,0,radius), QUNIT);
  ball->SetMaterialSurface(mmaterial);
  ball->SetMass(mass);
  ball->SetInertiaXX(0.4*radius*radius*ChVector<>(mass,mass,mass));
  ball->SetPos_dt(ChVector<>(initialVel,0,0));

  // -------------
  // Create DVI Serial system
  // -------------

  // Set up the serial system
  ChSystem* serialSystem = new ChSystem();
  serialSystem->Set_G_acc(gravity);
  serialSystem->SetTolForce(tolerance);
  serialSystem->SetIterLCPmaxItersSpeed(50000);
  serialSystem->SetMaxPenetrationRecoverySpeed(10e30);  // used by Anitescu stepper only
  serialSystem->SetUseSleeping(false);
  serialSystem->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
  serialSystem->SetIntegrationType(ChSystem::INT_ANITESCU);
  serialSystem->SetIterLCPwarmStarting(false);

  // Create the ground for serial system
  ChSharedPtr<ChBody> groundS = createBoxSerial(serialSystem, ChVector<>(10,2,2), ChVector<>(0,0,-1), QUNIT);
  groundS->SetBodyFixed(true);
  groundS->SetMaterialSurface(mmaterial);

  // Create a unit ball for serial system
  ChSharedPtr<ChBody> ballS = createBallSerial(serialSystem, radius, ChVector<>(0,0,radius), QUNIT);
  ballS->SetMaterialSurface(mmaterial);
  ballS->SetMass(mass);
  ballS->SetInertiaXX(0.4*radius*radius*ChVector<>(mass,mass,mass));
  ballS->SetPos_dt(ChVector<>(initialVel,0,0));

  // ----------------------
  // Perform the simulation
  // ----------------------

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow &gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "Rolling Test", msystem);
  gl_window.SetCamera(ChVector<>(0,10,0), ChVector<>(0,0,0),ChVector<>(0,0,1));
  gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

  ChStreamOutAsciiFile statStream("test_rolling.csv");

  // Loop until reaching the end time...
  while (msystem->GetChTime()<=simTime) {

    // Advance simulation by one step
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(hh);
      gl_window.Render();
    }
#else
    msystem->DoStepDynamics(hh);
#endif
    serialSystem->DoStepDynamics(hh);

    std::vector<double> history = ((ChLcpIterativeSolver*) (msystem->GetLcpSolverSpeed()))->GetViolationHistory();
    int numIters = history.size();
    printf("Time: %f, Iterations: %d\n",msystem->GetChTime(),numIters);

    msystem->CalculateContactForces();
    real3 force(0, 0, 0);
    force = msystem->GetBodyContactForce(1);
    printf("  Contact force on ball: (%f, %f, %f)\n",force.x,force.y,force.z);
    printf("  Trans. velocity of ball: (%f, %f, %f)\n",ball->GetPos_dt().x,ball->GetPos_dt().y,ball->GetPos_dt().z);
    printf("  Rot. velocity of ball  : (%f, %f, %f)\n",ball->GetWvel_par().x,ball->GetWvel_par().y,ball->GetWvel_par().z);
    statStream << msystem->GetChTime() << ", " << numIters << ", " << force.x << ", " << force.y << ", " << force.z << ", " << ball->GetPos_dt().x << ", " << ball->GetPos_dt().y << ", " << ball->GetPos_dt().z << ", " << ball->GetWvel_par().x << ", " << ball->GetWvel_par().y << ", " << ball->GetWvel_par().z << ", " << ballS->GetPos_dt().x << ", " << ballS->GetPos_dt().y << ", " << ballS->GetPos_dt().z << ", " << ballS->GetWvel_par().x << ", " << ballS->GetWvel_par().y << ", " << ballS->GetWvel_par().z << ", \n"; 
  }

  return 0;
}
