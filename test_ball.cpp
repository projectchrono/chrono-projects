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
// Authors: Jonathan Fleischmann & Radu Serban
// =============================================================================
//
// Chrono-Parallel test program to determine the *statically indeterminate*
// normal and tangential (static friction) contact forces on a single sphere
// supported by two walls in the shape of a "V"
//
// The global reference frame has Z up.
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

#include "demo_utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// Comment the following line to use DVI contact
#define USE_DEM

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Solver settings
#ifdef USE_DEM
double time_step = 1e-4;
double tolerance = 0.1;
int max_iteration_bilateral = 100;
CONTACTFORCEMODEL contact_force_model = HOOKE;
TANGENTIALDISPLACEMENTMODE tangential_displ_mode = MULTI_STEP;
#else
double time_step = 1e-4;
double tolerance = 0.1;
int max_iteration_normal = 0;
int max_iteration_sliding = 10000;
int max_iteration_spinning = 0;
int max_iteration_bilateral = 100;
double contact_recovery_speed = 10e30;
#endif

bool clamp_bilaterals = false;
double bilateral_clamp_speed = 0.1;

// Tilt angle (about global Y axis) of the container.
double tilt_angle = 45.0 * CH_C_PI / 180;

// Simulation duration
double end_simulation_time = 10;

// Ball properties
double mass = 1;
double radius = 0.5;

// Material properties (same on walls and ball)
float Y = 2e6f;
float nu = 0.3f;
float cr = 0.1f;

float mu = 10.0f;

// Output
#ifdef USE_DEM
const std::string out_dir = "../BALL_DEM";
#else
const std::string out_dir = "../BALL_DVI";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string force_file = out_dir + "/forces.dat";
const std::string stats_file = out_dir + "/stats.dat";

bool write_povray_data = false;

double data_out_step = 1e-3;       // time interval between data outputs
double visual_out_step = 1e-1;     // time interval between PovRay outputs


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  // Create output directories

  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

   // Simulation parameters

   double gravity = 9.81; // m/s^2

   // Create system

 #ifdef USE_DEM
   cout << "Create DEM system" << endl;
   const std::string title = "soft-sphere (DEM) single sphere test";
   ChBody::ContactMethod contact_method = ChBody::DEM;
   ChSystemParallelDEM* my_system = new ChSystemParallelDEM();
 #else
   cout << "Create DVI system" << endl;
   const std::string title = "hard-sphere (DVI) single sphere test";
   ChBody::ContactMethod contact_method = ChBody::DVI;
   ChSystemParallelDVI* my_system = new ChSystemParallelDVI();
 #endif

   my_system->Set_G_acc(ChVector<>(0, 0, -gravity));

   // Set number of threads

   int max_threads = omp_get_num_procs();
   if (threads > max_threads) threads = max_threads;
   my_system->SetParallelThreadNumber(threads);
   omp_set_num_threads(threads);

   my_system->GetSettings()->max_threads = threads;
   my_system->GetSettings()->perform_thread_tuning = thread_tuning;

   // Edit system settings

   my_system->GetSettings()->solver.use_full_inertia_tensor = false;
   my_system->GetSettings()->solver.tolerance = tolerance;
   my_system->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
   my_system->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
   my_system->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

 #ifdef USE_DEM
   my_system->GetSettings()->solver.contact_force_model = contact_force_model;
   my_system->GetSettings()->solver.tangential_displ_mode = tangential_displ_mode;
 #else
   my_system->GetSettings()->solver.solver_mode = SLIDING;
   my_system->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
   my_system->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
   my_system->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
   my_system->GetSettings()->solver.alpha = 0;
   my_system->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
   my_system->ChangeSolverType(APGD);

   my_system->GetSettings()->collision.collision_envelope = 0.05 * radius;
 #endif

   my_system->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
   my_system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

   // Create common material

#ifdef USE_DEM
  ChSharedPtr<ChMaterialSurfaceDEM> material;
  material = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  material->SetYoungModulus(Y);
  material->SetPoissonRatio(nu);
  material->SetRestitution(cr);
  material->SetFriction(mu);
#else
  ChSharedPtr<ChMaterialSurface> material;
  material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  material->SetRestitution(cr);
  material->SetFriction(mu);
#endif

  // Create walls

  int plate1_Id = -100;
  int plate2_Id = -200;

  ChVector<> hdim(2, 2, 0.5);
  double hthick = 0.1;

  ChSharedPtr<ChBody> plate1(new ChBody(new ChCollisionModelParallel, contact_method));
  plate1->SetMaterialSurface(material);
  plate1->SetIdentifier(plate1_Id);
  plate1->SetMass(1);
  plate1->SetPos(ChVector<>(0, 0, 0));
  plate1->SetRot(Q_from_AngY(tilt_angle));
  plate1->SetCollide(true);
  plate1->SetBodyFixed(true);

  plate1->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(plate1.get_ptr(), ChVector<>(hdim.x, hdim.y, hthick), ChVector<>(0, 0, -hthick));
  plate1->GetCollisionModel()->SetFamily(1);
  plate1->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);
  plate1->GetCollisionModel()->BuildModel();

  my_system->AddBody(plate1);

  ChSharedPtr<ChBody> plate2(new ChBody(new ChCollisionModelParallel, contact_method));
  plate2->SetMaterialSurface(material);
  plate2->SetIdentifier(plate2_Id);
  plate2->SetMass(1);
  plate2->SetPos(ChVector<>(0, 0, 0));
  plate2->SetRot(Q_from_AngY(-tilt_angle));
  plate2->SetCollide(true);
  plate2->SetBodyFixed(true);

  plate2->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(plate2.get_ptr(), ChVector<>(hdim.x, hdim.y, hthick), ChVector<>(0, 0, -hthick));
  plate2->GetCollisionModel()->SetFamily(2);
  plate2->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
  plate2->GetCollisionModel()->BuildModel();

  my_system->AddBody(plate2);

  // Create ball

  int ballId = 0;

  ChVector<> inertia = (2.0 / 5.0) * mass * radius * radius * ChVector<>(1, 1, 1);
  ChVector<> pos(0, 0, radius/cos(tilt_angle));

  ChSharedPtr<ChBody> ball(new ChBody(new ChCollisionModelParallel, contact_method));
  ball->SetMaterialSurface(material);

  ball->SetIdentifier(ballId);
  ball->SetMass(mass);
  ball->SetInertiaXX(inertia);
  ball->SetPos(pos);
  ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
  ball->SetBodyFixed(false);
  ball->SetCollide(true);

  ball->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball.get_ptr(), radius);
  ball->GetCollisionModel()->SetFamily(3);
  ball->GetCollisionModel()->BuildModel();

  my_system->AddBody(ball);

  // Setup output

  ChStreamOutAsciiFile forceStream(force_file.c_str());
  ChStreamOutAsciiFile statsStream(stats_file.c_str());
  forceStream.SetNumFormat("%16.4e");

  // Create the OpenGL visualization window

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(800, 600, title.c_str(), my_system);
  gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
  gl_window.SetRenderMode(opengl::SOLID);
#endif

  // Begin simulation

  int data_out_frame = 0;
  int visual_out_frame = 0;

  while (my_system->GetChTime() < end_simulation_time) {

    //  Do time step

#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(time_step);
      gl_window.Render();
    }
    else
      break;
#else
    my_system->DoStepDynamics(time_step);
#endif

//    TimingOutput(my_system, &statsStream);

    //  Output to screen and file

    if (my_system->GetChTime() >= data_out_frame * data_out_step) {
#ifndef USE_DEM
      my_system->CalculateContactForces();
#endif
      real3 frc0 = my_system->GetBodyContactForce(0);
      real3 frc1 = my_system->GetBodyContactForce(1);
      real F = frc1.z*sin(tilt_angle) + frc1.x*cos(tilt_angle);
      real N = frc1.z*cos(tilt_angle) - frc1.x*sin(tilt_angle);
      cout << frc0.x << "  " << frc0.y << "  " << frc0.z << "  " << frc1.x << "  " << frc1.y << "  " << frc1.z
           << "  F = " << -F << "  N = " << -N << "\n";
      forceStream << frc0.x << "  " << frc0.y << "  " << frc0.z << "  " << frc1.x << "  " << frc1.y << "  " << frc1.z
                  << "  F = " << -F << "  N = " << -N << "\n";

      data_out_frame++;
    }

    //  Output to POV-Ray

    if (write_povray_data && my_system->GetChTime() >= visual_out_frame * visual_out_step) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), visual_out_frame + 1);
      utils::WriteShapesPovray(my_system, filename, false);

      visual_out_frame++;
    }
  }

  return 0;
}
