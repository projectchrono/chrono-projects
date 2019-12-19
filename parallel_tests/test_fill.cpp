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
// ChronoParallel demo program for testing the filling process.
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

// Comment the following line to use parallel collision detection
//#define BULLET

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

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Save PovRay post-processing data?
bool write_povray_data = true;

// Load the bodies from a checkpoint file?
bool loadCheckPointFile = false;

// Simulation times
double time_settling_min = 0.1;
double time_settling_max = 3.0;

// Stopping criteria for settling (fraction of particle radius)
double settling_tol = 0.2;

// Solver settings
double time_step = 1e-3;
int max_iteration_normal = 0;
int max_iteration_sliding = 1000;
int max_iteration_spinning = 0;
int max_iteration_bilateral = 1000;
float contact_recovery_speed = 10e30;
bool clamp_bilaterals = false;
double bilateral_clamp_speed = 10e30;
double tolerance = 1e-2;

// Output
std::string out_dir = "../TEST_FILL";

std::string pov_dir = out_dir + "/POVRAY";
std::string fill_file = out_dir + "/filling.dat";
std::string stats_file = out_dir + "/stats.dat";
std::string settled_ckpnt_file = out_dir + "/settled.dat";

// Frequency for visualization output
int out_fps_settling = 120;

// Frequency for writing results to output file
int write_fps = 1000;

// Simulation frame at which detailed timing information is printed
int timing_frame = -1;

// Gravitational acceleration [m/s^2]
double gravity = 9.81;

// Parameters for the mechanism
int Id_container = 0;  // body ID for the containing bin
int Id_ground = 1;     // body ID for the ground

double hdimX = 2.0 / 2;   // [m] bin half-length in x direction
double hdimY = 2.0 / 2;   // [m] bin half-depth in y direction
double hdimZ = 2.0 / 2;   // [m] bin half-height in z direction
double hthick = 0.1 / 2;  // [m] bin half-thickness of the walls
float mu_walls = 0.3f;

// Parameters for the granular material
int Id_g = 2;         // start body ID for particles
double r_g = 0.1;     // [m] radius of granular sphers
double rho_g = 1000;  // [kg/m^3] density of granules
float mu_g = 0.5f;

// =============================================================================
// Create the containing bin (the ground) and the load plate.
//
// The load plate is constructed fixed to the ground.
// No joints between bodies are defined at this time.
// =============================================================================

void CreateMechanismBodies(ChSystemParallel* system) {
  // -------------------------------
  // Create a material for the walls
  // -------------------------------

  ChSharedPtr<ChMaterialSurface> mat_walls(new ChMaterialSurface);
  mat_walls->SetFriction(mu_walls);

  // ----------------------
  // Create the container body -- always FIRST body in system
  // ----------------------

  ChSharedPtr<ChBody> container(new ChBody(
#ifndef BULLET
      new ChCollisionModelParallel
#endif
      ));
  container->SetMaterialSurface(mat_walls);
  container->SetIdentifier(Id_container);
  container->SetBodyFixed(false);
  container->SetCollide(true);
  container->SetMass(10000.0);

  // Attach geometry of the containing bin
  container->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(container.get_ptr(), ChVector<>(hdimX, hdimY, hthick), ChVector<>(0, 0, -hthick));
  utils::AddBoxGeometry(container.get_ptr(), ChVector<>(hthick, hdimY, hdimZ), ChVector<>(-hdimX - hthick, 0, hdimZ));
  utils::AddBoxGeometry(container.get_ptr(), ChVector<>(hthick, hdimY, hdimZ), ChVector<>(hdimX + hthick, 0, hdimZ));
  utils::AddBoxGeometry(container.get_ptr(), ChVector<>(hdimX, hthick, hdimZ), ChVector<>(0, -hdimY - hthick, hdimZ));
  utils::AddBoxGeometry(container.get_ptr(), ChVector<>(hdimX, hthick, hdimZ), ChVector<>(0, hdimY + hthick, hdimZ));
  container->GetCollisionModel()->BuildModel();

  system->AddBody(container);

  // ----------------------
  // Create the ground body -- always SECOND body in system
  // ----------------------

  ChSharedPtr<ChBody> ground(new ChBody(
#ifndef BULLET
      new ChCollisionModelParallel
#endif
      ));
  ground->SetMaterialSurface(mat_walls);
  ground->SetIdentifier(Id_ground);
  ground->SetBodyFixed(true);
  ground->SetCollide(true);
  ground->SetMass(1.0);

  system->AddBody(ground);
}

// =============================================================================
// Create the granular material
//
// Granular material consisting of identical spheres with specified radius and
// material properties; the spheres are generated in a number of vertical
// layers with locations within each layer obtained using Poisson Disk sampling,
// thus ensuring that no two spheres are closer than twice the radius.
// =============================================================================

int CreateGranularMaterial(ChSystemParallel* system) {
  // -------------------------------------------
  // Create a material for the granular material
  // -------------------------------------------

  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
  mat_g->SetFriction(mu_g);

  // ---------------------------------------------
  // Create a mixture entirely made out of spheres
  // ---------------------------------------------

  // Create the particle generator with a mixture of 100% spheres
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
  m1->setDefaultMaterialDVI(mat_g);
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  // Ensure that all generated particle bodies will have positive IDs.
  gen.setBodyIdentifier(Id_g);

  // ----------------------
  // Generate the particles
  // ----------------------

  double r = 1.01 * r_g;
  ChVector<> hdims(hdimX - r, hdimY - r, 0);
  ChVector<> center(0, 0, 2 * r);

  while (center.z < 2 * hdimZ) {
    gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims);
    center.z += 2 * r;
  }

  // Return the number of generated particles.
  return gen.getTotalNumBodies();
}

double calculateContactForceOnBody(ChSystemParallel* system, int bodyIndex) {
  uint num_contacts = system->data_manager->num_rigid_contacts;
  uint num_unilaterals = system->data_manager->num_unilaterals;
  uint num_bilaterals = system->data_manager->num_bilaterals;

  DynamicVector<real> gamma = system->data_manager->host_data.gamma;
  blaze::DenseSubvector<DynamicVector<real> > gamma_b = blaze::subvector(gamma, num_unilaterals, num_bilaterals);
  blaze::DenseSubvector<DynamicVector<real> > gamma_n = blaze::subvector(gamma, 0, num_contacts);
  blaze::DenseSubvector<DynamicVector<real> > gamma_t = blaze::subvector(gamma, num_contacts, num_contacts * 2);

  CompressedMatrix<real> D_n = system->data_manager->host_data.D_n;
  CompressedMatrix<real> D_t = system->data_manager->host_data.D_t;
  CompressedMatrix<real> D_b = system->data_manager->host_data.D_b;

  // F_contact = D*gamma/h;
  blaze::DynamicVector<real> contactForces;
  contactForces.resize(system->data_manager->host_data.gamma.size());
  contactForces = (D_n * gamma_n + D_t * gamma_t) / time_step;  // Don't include the bilateral

  // NOTE: contactForces is now a vector of length 6*numBodies that contains the contact force/torque on each body
  //       unfortunately, the D matrix does not include entries for inactive bodies meaning that there will be 0's
  //       for any inactive body! (As a hack, you can go into ChConstraintRigidRigid.cpp and force the solver to
  //       fill in the entries for the fixed bodies in D).

  double force_z = contactForces[6 * bodyIndex + 2];

  return force_z;
}

// =============================================================================

int main(int argc, char* argv[]) {
  if (argc > 1) {
    tolerance = atof(argv[1]);
    hdimZ = atof(argv[2]);
    out_dir = argv[3];
    pov_dir = out_dir + "/POVRAY";
    fill_file = out_dir + "/filling.dat";
    stats_file = out_dir + "/stats.dat";
    settled_ckpnt_file = out_dir + "/settled.dat";
  }

  // Create output directories.
  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    // cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    // cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

  // -------------
  // Create system
  // -------------

  // cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
  msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

#ifdef BULLET
  msystem->ChangeCollisionSystem(COLLSYS_BULLET_PARALLEL);
#endif

  // Set number of threads.
  int max_threads = omp_get_num_procs();
  if (threads > max_threads)
    threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  // cout << "Using " << threads << " threads" << endl;
  msystem->GetSettings()->max_threads = threads;
  msystem->GetSettings()->perform_thread_tuning = thread_tuning;

  // Edit system settings
  msystem->GetSettings()->solver.tolerance = tolerance;
  msystem->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
  msystem->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
  msystem->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;
  msystem->GetSettings()->solver.solver_mode = SLIDING;
  msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
  msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
  msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
  msystem->GetSettings()->solver.alpha = 0;
  msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  msystem->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  msystem->ChangeSolverType(APGD);
  msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);

  // --------------
  // Problem set up
  // --------------

  // Depending on problem type:
  // - Select end simulation time
  // - Select output FPS
  // - Create / load objects

  double time_min = 0;
  double time_end;
  int out_fps;

  time_min = time_settling_min;
  time_end = time_settling_max;
  out_fps = out_fps_settling;

  int num_particles = 0;
  if (loadCheckPointFile) {
    // cout << "Read checkpoint data from " << settled_ckpnt_file;
    utils::ReadCheckpoint(msystem, settled_ckpnt_file);
    // cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;
  } else {
    // Create the mechanism bodies (all fixed).
    CreateMechanismBodies(msystem);

    // Create granular material.
    num_particles = CreateGranularMaterial(msystem);
    // cout << "Granular material:  " << num_particles << " particles" << endl;
  }

  // Lock the container to the ground
  ChSharedPtr<ChBody> container = ChSharedPtr<ChBody>(msystem->Get_bodylist()->at(0));
  ChSharedPtr<ChBody> ground = ChSharedPtr<ChBody>(msystem->Get_bodylist()->at(1));
  msystem->Get_bodylist()->at(0)->AddRef();
  msystem->Get_bodylist()->at(1)->AddRef();

  ChSharedPtr<ChLinkLockLock> lock(new ChLinkLockLock);
  lock->Initialize(
      container, ground, false, ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT), ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));
  msystem->AddLink(lock);

  // ----------------------
  // Perform the simulation
  // ----------------------

  // Set number of simulation steps and steps between successive output
  int num_steps = (int)std::ceil(time_end / time_step);
  int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);
  int write_steps = (int)std::ceil((1.0 / time_step) / write_fps);

  // Initialize counters
  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  int next_out_frame = 0;
  double exec_time = 0;
  int num_contacts = 0;

  // Circular buffer with highest particle location
  // (only used for SETTLING or PRESSING)
  int buffer_size = std::ceil(time_min / time_step);
  std::valarray<double> hdata(0.0, buffer_size);

  // Create output files
  ChStreamOutAsciiFile statsStream(stats_file.c_str());
  ChStreamOutAsciiFile fillStream(fill_file.c_str());
  fillStream.SetNumFormat("%16.4e");
  statsStream.SetNumFormat("%16.4e");

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "Filling Test", msystem);
  gl_window.SetCamera(ChVector<>(0, -10 * hdimY, hdimZ), ChVector<>(0, 0, hdimZ), ChVector<>(0, 0, 1));
  gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

  // Loop until reaching the end time...
  while (time < time_end) {
    // If at an output frame, write PovRay file and print info
    if (sim_frame == next_out_frame) {
      cout << "------------ Output frame:   " << out_frame + 1 << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Execution time: " << exec_time << endl;

      // Save PovRay post-processing data.
      if (write_povray_data) {
        char filename[100];
        sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
        utils::WriteShapesPovray(msystem, filename, false);
      }

      // Create a checkpoint from the current state.
      cout << "             Write checkpoint data " << flush;
      utils::WriteCheckpoint(msystem, settled_ckpnt_file);
      cout << msystem->Get_bodylist()->size() << " bodies" << endl;

      // Increment counters
      out_frame++;
      next_out_frame += out_steps;
    }

// Advance simulation by one step
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(time_step);
      gl_window.Render();
    }
#else
    msystem->DoStepDynamics(time_step);
#endif

    // Record stats about the simulation
    if (sim_frame % write_steps == 0) {
      // Compute contact force on container (container should be the first body)
      double force = 0;
      if (msystem->GetNcontacts()) {
        force = -calculateContactForceOnBody(msystem, 0);
      }
      double actualWeight =
          (msystem->Get_bodylist()->size() - 1) * (4.0 / 3.0) * CH_C_PI * pow(r_g, 3.0) * rho_g * gravity;

      // get maximum body velocity
      double maxVelocity = 0;
      for (int i = 0; i < msystem->Get_bodylist()->size(); ++i) {
        ChBody* body = (ChBody*)msystem->Get_bodylist()->at(i);
        double vel =
            sqrt(pow(body->GetPos_dt().x, 2.0) + pow(body->GetPos_dt().y, 2.0) + pow(body->GetPos_dt().z, 2.0));
        if (vel > maxVelocity)
          maxVelocity = vel;
      }

      // write fill info
      fillStream << time << ", " << -(lock->Get_react_force().z + gravity * container->GetMass()) << ", " << force
                 << ", " << actualWeight << ", " << maxVelocity << ", \n";
      fillStream.GetFstream().flush();

      // write stat info
      int numIters = msystem->data_manager->measures.solver.maxd_hist.size();
      double residual = 0;
      if (numIters)
        residual = msystem->data_manager->measures.solver.residual;
      statsStream << time << ", " << exec_time << ", " << num_contacts / write_steps << ", " << numIters << ", "
                  << residual << ", \n";
      statsStream.GetFstream().flush();

      num_contacts = 0;
    }

    // Increment counters
    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
    num_contacts += msystem->GetNcontacts();

    // If requested, output detailed timing information for this step
    if (sim_frame == timing_frame)
      msystem->PrintStepStats();
  }

  // ----------------
  // Final processing
  // ----------------

  // Create a checkpoint from the last state
  cout << "             Write checkpoint data " << flush;
  utils::WriteCheckpoint(msystem, settled_ckpnt_file);
  cout << msystem->Get_bodylist()->size() << " bodies" << endl;

  // Final stats
  cout << "==================================" << endl;
  cout << "Number of bodies:  " << msystem->Get_bodylist()->size() << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}
