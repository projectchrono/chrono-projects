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
// Authors: Daniel Melanz, Radu Serban
// =============================================================================
//
// ChronoParallel demo program for mass flow rate studies.
//
// The model simulated here consists of a granular material that flows out of a 
// container and the mass of the collected material is measured over time, using 
// either penalty or complementarity method for frictional contact.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <stdio.h>
#include <vector>
#include <cmath>

#include "chrono_parallel/ChSystemParallel.h"
#include "chrono_parallel/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

using namespace chrono;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use DVI contact
////#define DEM

enum ProblemType {
  SETTLING,
  DROPPING
};

ProblemType problem = SETTLING;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 100;

// Perform dynamic tuning of number of threads?
bool thread_tuning = true;

// Simulation parameters
double gravity = 9.81;

double time_settling_min = 0.1;
double time_settling_max = 1.0;
double time_dropping_max = 6.0;

#ifdef DEM
double time_step = 1e-5;
int max_iteration = 20;
#else
double time_step = 1e-4;
int max_iteration_normal = 30;
int max_iteration_sliding = 20;
int max_iteration_spinning = 0;
float contact_recovery_speed = 0.1;
#endif

int max_iteration_bilateral = 0;

// Output
#ifdef DEM
const char* out_folder = "../MASSFLOW_DEM/POVRAY";
const char* flow_file = "../MASSFLOW_DEM/flow.dat";
const char* stats_file = "../MASSFLOW_DEM/stats.dat";
const char* checkpoint_file = "../MASSFLOW_DEM/settled.dat";
#else
const char* out_folder = "../MASSFLOW_DVI/POVRAY";
const char* flow_file = "../MASSFLOW_DVI/flow.dat";
const char* stats_file = "../MASSFLOW_DVI/stats.dat";
const char* checkpoint_file = "../MASSFLOW_DVI/settled.dat";
#endif

int out_fps_settling = 200;
int out_fps_dropping = 200;

int timing_frame = -1;   // output detailed step timing at this frame

// Parameters for the granular material
double r_g = 0.25e-3;
double rho_g = 2500.0;

float  Y_g = 1e7;
float  alpha_g = 0.6;
float  mu_g = 0.3;

// Desired number of particles and X-Y dimensions of the sampling volume for
// granular material.
int    desired_num_particles = 400;

// Parameters for the mechanism material
float  Y_c = 2e6;
float  alpha_c = 0.6;
float  mu_c = 0.4;

// Dimensions of mechanism
double height = 6.0e-2;       // height of the cavity
////double width = 0.9525e-2;     // width of the cavity
double width = 0.9398e-2;     // width of the cavity
double thickness = 0.25e-2;   // thickness of walls

double height_insert = sqrt(2*height*height);
double delta = sqrt(thickness*thickness/8);

////double speed = 1.5e-3;       // speed of the angled insert
double speed = 1.0e-3;       // speed of the angled insert
double gap = 2e-3;           // size of gap

double time_opening = gap / speed;

// Dimensions of collector
double pos_collector = 4.0e-2;    // position below measuring line
double size_collector = 8.0e-2;   // width and length of collector bin
double height_collector = 1.0e-2; // height of collector walls


// -----------------------------------------------------------------------------
// Create mechanism
// -----------------------------------------------------------------------------
ChBody* CreateMechanism(ChSystemParallel* system)
{

  // Create the common material
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_b;
  mat_b = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_b->SetYoungModulus(Y_c);
  mat_b->SetFriction(mu_c);
  mat_b->SetDissipationFactor(alpha_c);
#else
  ChSharedPtr<ChMaterialSurface> mat_b(new ChMaterialSurface);
  mat_b->SetFriction(mu_c);
#endif

  // Angled insert
#ifdef DEM
  ChSharedBodyDEMPtr insert(new ChBodyDEM(new ChCollisionModelParallel));
  insert->SetMaterialSurfaceDEM(mat_b);
#else
  ChSharedBodyPtr insert(new ChBody(new ChCollisionModelParallel));
  insert->SetMaterialSurface(mat_b);
#endif

  insert->SetIdentifier(0);
  insert->SetMass(1);
  insert->SetInertiaXX(ChVector<>(1,1,1));
  insert->SetPos(ChVector<>(-0.5*height-delta, 0, 0.5*height-delta));
  insert->SetRot(chrono::Q_from_AngAxis(-CH_C_PI/4, VECT_Y));
  insert->SetCollide(true);
  insert->SetBodyFixed(true);

  insert->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(insert.get_ptr(), ChVector<>(thickness*0.5,width*0.5,height_insert*0.5));
  insert->GetCollisionModel()->BuildModel();

  system->AddBody(insert);

  // Static slot (back wall)
#ifdef DEM
  ChSharedBodyDEMPtr slot(new ChBodyDEM(new ChCollisionModelParallel));
  slot->SetMaterialSurfaceDEM(mat_b);
#else
  ChSharedBodyPtr slot(new ChBody(new ChCollisionModelParallel));
  slot->SetMaterialSurface(mat_b);
#endif

  slot->SetIdentifier(-1);
  slot->SetMass(1);
  slot->SetInertiaXX(ChVector<>(1,1,1));
  slot->SetPos(ChVector<>(0.5*thickness, 0, 0.5*height));
  slot->SetRot(ChQuaternion<>(1, 0, 0, 0));
  slot->SetCollide(true);
  slot->SetBodyFixed(true);

  slot->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(slot.get_ptr(), ChVector<>(thickness/2, width/2, height/2), ChVector<>(0, 0, 0));
  slot->GetCollisionModel()->BuildModel();

  system->AddBody(slot);

  // Lateral walls
#ifdef DEM
  ChSharedBodyDEMPtr wall(new ChBodyDEM(new ChCollisionModelParallel));
  wall->SetMaterialSurfaceDEM(mat_b);
#else
  ChSharedBodyPtr wall(new ChBody(new ChCollisionModelParallel));
  wall->SetMaterialSurface(mat_b);
#endif

  wall->SetIdentifier(-2);
  wall->SetMass(1);
  wall->SetInertiaXX(ChVector<>(1, 1, 1));
  wall->SetPos(ChVector<>(0, 0, 0));
  wall->SetRot(ChQuaternion<>(1, 0, 0, 0));
  wall->SetCollide(true);
  wall->SetBodyFixed(true);

  wall->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(wall.get_ptr(), ChVector<>(3*height/2, thickness/2, height), ChVector<>(0,  width/2+thickness/2, height/2));
  utils::AddBoxGeometry(wall.get_ptr(), ChVector<>(3*height/2, thickness/2, height), ChVector<>(0, -width/2-thickness/2, height/2));
  wall->GetCollisionModel()->BuildModel();

  system->AddBody(wall);

  // Containing bin
#ifdef DEM
  utils::CreateBoxContainerDEM(system, -3, mat_b, ChVector<>(size_collector/2, size_collector/2, height_collector/2), thickness/2, ChVector<>(0, 0, -pos_collector));
#else
  utils::CreateBoxContainerDVI(system, -3, mat_b, ChVector<>(size_collector/2, size_collector/2, height_collector/2), thickness/2, ChVector<>(0, 0, -pos_collector));
#endif

  // Return the angled insert body
  return insert.get_ptr();
}


// -----------------------------------------------------------------------------
// Create granular material
// -----------------------------------------------------------------------------
void CreateParticles(ChSystemParallel* system)
{
  // Create a material for the granular material
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_g;
  mat_g = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_g->SetYoungModulus(Y_g);
  mat_g->SetFriction(mu_g);
  mat_g->SetDissipationFactor(alpha_g);
#else
  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
  mat_g->SetFriction(mu_g);
#endif

  // Create a mixture entirely made out of spheres
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
#ifdef DEM
  m1->setDefaultMaterialDEM(mat_g);
#else
  m1->setDefaultMaterialDVI(mat_g);
#endif
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  gen.setBodyIdentifier(1);

  ChVector<> hdims(0.3 * height, 0.3 * width, 0);
  ChVector<> center(-0.4 * height, 0, 0.8 * height);
  ChVector<> vel(0, 0, 0);
  double r = 1.01 * r_g;

  while (gen.getTotalNumBodies() < desired_num_particles) {
    gen.createObjectsBox(utils::POISSON_DISK, 2 * r, center, hdims, vel);
    center.z += 2 * r;
  }

  std::cout << "Number of particles: " << gen.getTotalNumBodies() << std::endl;
}


// -----------------------------------------------------------------------------
// Find and return the body with specified identifier.
// -----------------------------------------------------------------------------
ChBody* FindBodyById(ChSystemParallel* sys, int id)
{
  for (int i = 0; i < sys->GetNumBodies(); ++i) {
    ChBody* body = (ChBody*) sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() == id)
      return body;
  }

  return NULL;
}


// -----------------------------------------------------------------------------
// Find the number of particles whose height is below and above, respectively,
// the specified value.
// -----------------------------------------------------------------------------
int GetNumParticlesBelowHeight(ChSystemParallel* sys, double value)
{
  int count = 0;
  for (int i = 0; i < sys->GetNumBodies(); ++i) {
    ChBody* body = (ChBody*) sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z < value)
      count++;
  }
  return count;
}

int GetNumParticlesAboveHeight(ChSystemParallel* sys, double value)
{
  int count = 0;
  for (int i = 0; i < sys->GetNumBodies(); ++i) {
    ChBody* body = (ChBody*) sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z > value)
      count++;
  }
  return count;
}


// -----------------------------------------------------------------------------
// Return true if all bodies in the granular mix have a linear velocity whose
// magnitude is below the specified value.
// -----------------------------------------------------------------------------
bool CheckSettled(ChSystemParallel* sys, double threshold)
{
  double t2 = threshold * threshold;

  for (int i = 0; i < sys->GetNumBodies(); ++i) {
    ChBody* body = (ChBody*) sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0) {
      double vel2 = body->GetPos_dt().Length2();
      if (vel2 > t2)
        return false;
    }
  }

  return true;
}


// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  // Create system
#ifdef DEM
  cout << "Create DEM system" << endl;
  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();
#else
  cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
#endif

  // Set number of threads.
  int max_threads = msystem->GetParallelThreadNumber();
  if (threads > max_threads)
    threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;

  msystem->DoThreadTuning(thread_tuning);

  // Set gravitational acceleration
  msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

  // Edit system settings
  msystem->SetTol(1e-3);
  msystem->SetTolSpeeds(1e-3);
  msystem->SetStep(time_step);

#ifdef DEM
  ((ChLcpSolverParallelDEM*) msystem->GetLcpSolverSpeed())->SetMaxIterationBilateral(max_iteration_bilateral);
  ((ChLcpSolverParallelDEM*) msystem->GetLcpSolverSpeed())->SetTolerance(1e-3);

  ((ChCollisionSystemParallel*) msystem->GetCollisionSystem())->ChangeNarrowphase(new ChCNarrowphaseR);
#else
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetMaxIterationNormal(max_iteration_normal);
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetMaxIterationSliding(max_iteration_sliding);
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetMaxIterationSpinning(max_iteration_spinning);
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetMaxIterationBilateral(max_iteration_bilateral);
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetTolerance(1e-3);
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetCompliance(0);
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetContactRecoverySpeed(contact_recovery_speed);
  ((ChLcpSolverParallelDVI*) msystem->GetLcpSolverSpeed())->SetSolverType(APGDRS);

  ((ChCollisionSystemParallel*) msystem->GetCollisionSystem())->SetCollisionEnvelope(0.05 * r_g);
#endif

  ((ChCollisionSystemParallel*) msystem->GetCollisionSystem())->setBinsPerAxis(I3(10, 10, 10));

  // Set simulation duration and create bodies (depending on problem type).
  double time_end;
  int out_fps;
  ChBody* insert;

  switch (problem) {
  case SETTLING:
    time_end = time_settling_max;
    out_fps = out_fps_settling;
    insert = CreateMechanism(msystem);
    CreateParticles(msystem);
    break;

  case DROPPING:
    time_end = time_dropping_max;
    out_fps = out_fps_dropping;
    utils::ReadCheckpoint(msystem, checkpoint_file);
    insert = FindBodyById(msystem, 0);
    break;
  }

  // Number of steps
  int num_steps = std::ceil(time_end / time_step);
  int out_steps = std::ceil((1.0 / time_step) / out_fps);

  // Zero velocity level for settling check
  double zero_v = 2 * r_g;

  // Perform the simulation
  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  int next_out_frame = 0;
  double exec_time = 0;
  int num_contacts = 0;
  ChStreamOutAsciiFile sfile(stats_file);
  ChStreamOutAsciiFile ffile(flow_file);

  while (time < time_end) {
    // Output data
    if (sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", out_folder, out_frame + 1);
      utils::WriteShapesPovray(msystem, filename);

      cout << "------------ Output frame:   " << out_frame << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
      cout << "             Execution time: " << exec_time << endl;

      double opening = insert->GetPos().x + 0.5 * height + delta;
      int count = GetNumParticlesBelowHeight(msystem, 0);

      cout << "             Gap:            " << -opening << endl;
      cout << "             Flow:           " << count << endl;

      sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";

      switch (problem) {
      case SETTLING:
        // Create a checkpoint from the current state.
        utils::WriteCheckpoint(msystem, checkpoint_file);
        cout << "             Checkpoint:     " << msystem->Get_bodylist()->size() << " bodies" << endl;
        break;
      case DROPPING:
        // Save current gap opening and number of dropped particles.
        ffile << time << "  " << -opening << "  " << count << "\n";
        break;
      }

      out_frame++;
      next_out_frame += out_steps;
      num_contacts = 0;
    }

    // Check for early termination of settling phase.
    if (problem == SETTLING && time > time_settling_min && CheckSettled(msystem, zero_v)) {
      cout << "Granular material settled...  time = " << time << endl;
      break;
    }

    // Check for early termination of dropping phase.
    if (problem == DROPPING && time > time_opening && GetNumParticlesAboveHeight(msystem, -pos_collector/2) == 0) {
      cout << "Granular material exhausted... time = " << time << endl;
      break;
    }

    // Advance system state by one step.
    msystem->DoStepDynamics(time_step);

    // Open the gate until it reaches the specified gap distance.
    if (problem == DROPPING && time <= time_opening) {
        insert->SetPos(ChVector<>(-0.5 * height - delta - time * speed, 0, 0.5 * height - delta));
        insert->SetPos_dt(ChVector<>(-speed, 0, 0));
    }

    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
    num_contacts += msystem->GetNcontacts();

    // If requested, output detailed timing information for this step
    if (sim_frame == timing_frame)
      msystem->PrintStepStats();
  }

  // Create a checkpoint from the last state
  if (problem == SETTLING) {
    cout << "Write checkpoint data to " << checkpoint_file;
    utils::WriteCheckpoint(msystem, checkpoint_file);
    cout << "  done.  Wrote " << msystem->Get_bodylist()->size() << " bodies." << endl;
  }

  // Final stats
  cout << "==================================" << endl;
  cout << "Number of bodies:  " << msystem->GetNumBodies() << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}
