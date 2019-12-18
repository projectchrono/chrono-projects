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
// ChronoParallel test for comparing APGD solver to Chrono.
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

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// =============================================================================
// Create the containing bin (the ground) and the load plate.
//
// The load plate is constructed fixed to the ground.
// No joints between bodies are defined at this time.
// =============================================================================

void CreateMechanismBodies(ChSystem* system, double hdimX, double hdimY, double hdimZ, double hthick, double mu_walls) {
  // Parameters for the mechanism
  int Id_container = 0;  // body ID for the containing bin

  // -------------------------------
  // Create a material for the walls
  // -------------------------------

  ChSharedPtr<ChMaterialSurface> mat_walls(new ChMaterialSurface);
  mat_walls->SetFriction(mu_walls);

  // ----------------------
  // Create the container body -- always FIRST body in system
  // ----------------------

  ChSharedPtr<ChBody> container(new ChBody());
  container->SetMaterialSurface(mat_walls);
  container->SetIdentifier(Id_container);
  container->SetBodyFixed(true);
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
}

// =============================================================================
// Create the granular material
//
// Granular material consisting of identical spheres with specified radius and
// material properties; the spheres are generated in a number of vertical
// layers with locations within each layer obtained using Poisson Disk sampling,
// thus ensuring that no two spheres are closer than twice the radius.
// =============================================================================

int CreateGranularMaterial(ChSystem* system,
                           double hdimX,
                           double hdimY,
                           double hdimZ,
                           double r_g,
                           double rho_g,
                           double mu_g) {
  int Id_g = 1;  // start body ID for particles

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

// =============================================================================

int generateSettledFile(ChSystem* msystem, double timeStep, double endTime, std::string settled_ckpnt_file) {
  // Parameters
  double r_g = 0.8;     // [m] radius of granular sphers
  double rho_g = 1000;  // [kg/m^3] density of granules
  double mu_g = 0.5f;
  double hdimX = 2.0 / 2;   // [m] bin half-length in x direction
  double hdimY = 2.0 / 2;   // [m] bin half-depth in y direction
  double hdimZ = 2.0 / 2;   // [m] bin half-height in z direction
  double hthick = 0.1 / 2;  // [m] bin half-thickness of the walls
  double mu_walls = 0.3f;

  // Create the bodies
  CreateMechanismBodies(msystem, hdimX, hdimY, hdimZ, hthick, mu_walls);

  // Create the granular material
  CreateGranularMaterial(msystem, hdimX, hdimY, hdimZ, r_g, rho_g, mu_g);

  while (msystem->GetChTime() < endTime) {
    msystem->DoStepDynamics(timeStep);
    printf("  Time: %f\n", msystem->GetChTime());
  }

  utils::WriteCheckpoint(msystem, settled_ckpnt_file);

  return 0;
}

int main(int argc, char* argv[]) {
  // Step 0: Set up the simulation parameters

  // Solver settings
  double time_step = 1e-3;
  int max_iteration_normal = 0;
  int max_iteration_sliding = 1000;
  int max_iteration_spinning = 0;
  int max_iteration_bilateral = 0;
  float contact_recovery_speed = 10e30;
  bool clamp_bilaterals = false;
  double bilateral_clamp_speed = 10e30;
  double tolerance = 1e-2;

  // Desired number of OpenMP threads (will be clamped to maximum available)
  int threads = 1;

  // Perform dynamic tuning of number of threads?
  bool thread_tuning = false;

  // Gravitational acceleration [m/s^2]
  double gravity = 9.81;

  // Output
  std::string out_dir = "../TEST_SERIAL";

  std::string pov_dir = out_dir + "/POVRAY";
  std::string stats_file = out_dir + "/stats.dat";
  std::string settled_ckpnt_file = out_dir + "/settled.dat";
  std::string parallel_out_dir = out_dir + "/parallel/";
  std::string serial_out_dir = out_dir + "/serial/";

  // Create output directories.
  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(parallel_out_dir.c_str()) < 0) {
    cout << "Error creating directory " << parallel_out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(serial_out_dir.c_str()) < 0) {
    cout << "Error creating directory " << serial_out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }
  // End Step 0

  // Step 1: Set up the serial system
  ChSystem* serialSystem = new ChSystem();
  serialSystem->Set_G_acc(ChVector<>(0, 0, -gravity));
  serialSystem->SetTolForce(tolerance);
  serialSystem->SetIterLCPmaxItersSpeed(max_iteration_sliding);
  serialSystem->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);  // used by Anitescu stepper only
  serialSystem->SetUseSleeping(false);
  serialSystem->SetLcpSolverType(ChSystem::LCP_ITERATIVE_APGD);
  serialSystem->SetIntegrationType(ChSystem::INT_ANITESCU);
  serialSystem->SetIterLCPwarmStarting(true);
  cout << "Serial system created!" << endl;
  // End Step 1

  // Step 2: Set up the parallel system
  ChSystemParallelDVI* parallelSystem = new ChSystemParallelDVI();
  parallelSystem->Set_G_acc(ChVector<>(0, 0, -gravity));

  // Make sure to use consistent collision detection
  parallelSystem->ChangeCollisionSystem(COLLSYS_BULLET_PARALLEL);

  int max_threads = omp_get_num_procs();
  if (threads > max_threads)
    threads = max_threads;
  parallelSystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  parallelSystem->GetSettings()->max_threads = threads;
  parallelSystem->GetSettings()->perform_thread_tuning = thread_tuning;

  parallelSystem->GetSettings()->solver.tolerance = tolerance;
  parallelSystem->GetSettings()->solver.max_iteration_bilateral = 0;
  parallelSystem->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
  parallelSystem->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;
  parallelSystem->GetSettings()->solver.solver_mode = SLIDING;
  parallelSystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
  parallelSystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
  parallelSystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
  parallelSystem->GetSettings()->solver.alpha = 0;
  parallelSystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  parallelSystem->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  parallelSystem->ChangeSolverType(APGDREF);
  parallelSystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
  cout << "Parallel system created!" << endl;
  // End Step 2

  // Step 3: Generate the settled file with serial solver, load into parallel solver
  generateSettledFile(serialSystem, time_step, 1.0, settled_ckpnt_file);
  utils::ReadCheckpoint(parallelSystem, settled_ckpnt_file);
  cout << "Generated settled file!" << endl << endl;
  // End Step 3

  // Step 4: Verify that the systems match each other
  cout << "Verifying system is the same before step ... " << endl;
  // Check position and velocity
  ChBody* bodySerial;
  ChBody* bodyParallel;
  double posError = 0;
  double velError = 0;
  for (int i = 0; i < serialSystem->Get_bodylist()->size(); i++) {
    bodySerial = (ChBody*)serialSystem->Get_bodylist()->at(i);
    bodyParallel = (ChBody*)parallelSystem->Get_bodylist()->at(i);

    posError += fabs(bodySerial->GetPos().x - bodyParallel->GetPos().x);
    posError += fabs(bodySerial->GetPos().y - bodyParallel->GetPos().y);
    posError += fabs(bodySerial->GetPos().z - bodyParallel->GetPos().z);

    velError += fabs(bodySerial->GetPos_dt().x - bodyParallel->GetPos_dt().x);
    velError += fabs(bodySerial->GetPos_dt().y - bodyParallel->GetPos_dt().y);
    velError += fabs(bodySerial->GetPos_dt().z - bodyParallel->GetPos_dt().z);
  }
  cout << "  Position error: " << posError << endl;
  cout << "  Velocity error: " << velError << endl;

  cout << "System verified!" << endl << endl;
  // End Step 4

  // Step 5: Run a single APGD iteration with serial system
  cout << "Running 1 APGD iteration (serial) ... ";
  serialSystem->SetIterLCPmaxItersSpeed(1);
  serialSystem->DoStepDynamics(time_step);
  serialSystem->GetLcpSystemDescriptor()->DumpLastMatrices(serial_out_dir.c_str());
  cout << "Complete!" << endl;
  // End Step 5

  // Step 6: Run a single APGD iteration with parallel system
  cout << "Running 1 APGD iteration (parallel) ... ";
  parallelSystem->GetSettings()->solver.max_iteration_sliding = 1;
  parallelSystem->DoStepDynamics(time_step);
  parallelSystem->data_manager->ExportCurrentSystem(parallel_out_dir);
  cout << "Complete!" << endl << endl;
  // End Step 6

  // Step 7: Compare results
  cout << "Verifying system is the same after step ... " << endl;
  // Check position and velocity
  posError = 0;
  velError = 0;
  for (int i = 0; i < serialSystem->Get_bodylist()->size(); i++) {
    bodySerial = (ChBody*)serialSystem->Get_bodylist()->at(i);
    bodyParallel = (ChBody*)parallelSystem->Get_bodylist()->at(i);

    posError += fabs(bodySerial->GetPos().x - bodyParallel->GetPos().x);
    posError += fabs(bodySerial->GetPos().y - bodyParallel->GetPos().y);
    posError += fabs(bodySerial->GetPos().z - bodyParallel->GetPos().z);

    velError += fabs(bodySerial->GetPos_dt().x - bodyParallel->GetPos_dt().x);
    velError += fabs(bodySerial->GetPos_dt().y - bodyParallel->GetPos_dt().y);
    velError += fabs(bodySerial->GetPos_dt().z - bodyParallel->GetPos_dt().z);
  }
  cout << "  Position error: " << posError << endl;
  cout << "  Velocity error: " << velError << endl;
  cout << "System verified!" << endl;
  // End Step 7

  return 0;
}
