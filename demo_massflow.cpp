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
bool thread_tuning = false;

// Simulation parameters
double gravity = 981; // gravity, cm/s^2

double time_settling_min = 0.1;
double time_settling_max = 1.0;

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
const char* height_file = "../MASSFLOW_DEM/height.dat";
const char* stats_file = "../MASSFLOW_DEM/stats.dat";
const char* checkpoint_file = "../MASSFLOW_DEM/settled.dat";
#else
const char* out_folder = "../MASSFLOW_DVI/POVRAY";
const char* height_file = "../MASSFLOW_DVI/height.dat";
const char* stats_file = "../MASSFLOW_DVI/stats.dat";
const char* checkpoint_file = "../MASSFLOW_DVI/settled.dat";
#endif

int out_fps_settling = 120;
int out_fps_dropping = 1200;

int timing_frame = -1;   // output detailed step timing at this frame

// Parameters for the granular material
double r_g = 0.025;        // cm
double rho_g = 2.5;        // g/cm^3

float  Y_g = 1e6;          // g/cm/s^2
float  alpha_g = 0.1;      // s/cm
float  mu_g = 0.3;

// Desired number of particles and X-Y dimensions of the sampling volume for
// granular material.
int    desired_num_particles = 400;

// Parameters for the mechanism material
float  Y_c = 1e6;          // g/cm/s^2
float  alpha_c = 0.1;      // s/cm
float  mu_c = 0.4;

// Dimensions of mechanism
double height = 3.175;     // height of the cavity, cm
double width = 1;          // width of the cavity, cm
double thickness = 0.25;   // thickness of walls, cm
//double angle = 45;         // angle of the cavity, degrees

double height_insert = sqrt(2*height*height);
double delta = sqrt(thickness*thickness/8);

double speed = 0.15;       // speed of the angled insert, cm/s
double gap = 0.2;          // size of gap, cm


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

  // Static slot
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
  utils::AddBoxGeometry(slot.get_ptr(), ChVector<>(thickness*0.5,width*0.5,height*0.5), ChVector<>(0,0,0));
  utils::AddBoxGeometry(slot.get_ptr(), ChVector<>(3*height*0.5,width*0.5,height),      ChVector<>(0,width,0.5*height));
  utils::AddBoxGeometry(slot.get_ptr(), ChVector<>(3*height*0.5,width*0.5,height),      ChVector<>(0,-width,0.5*height));
  slot->GetCollisionModel()->BuildModel();

  system->AddBody(slot);

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
  utils::AddBoxGeometry(insert.get_ptr(), ChVector<>(thickness*0.5,width*0.5,height_insert*0.5),  ChVector<>(0,0,0));
  insert->GetCollisionModel()->BuildModel();

  system->AddBody(insert);

  // Containing bin
#ifdef DEM
  utils::CreateBoxContainerDEM(system, -3, mat_b, ChVector<>(0.5*height, 0.5*height, 0.1*height), thickness, ChVector<>(0, 0, -1.5*height));
#else
  utils::CreateBoxContainerDVI(system, -3, mat_b, ChVector<>(0.5*height, 0.5*height, 0.1*height), thickness, ChVector<>(0, 0, -1.5*height));
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

  double r = 1.01 * r_g;
  //gen.createObjectsBox(utils::POISSON_DISK, 2 * r, ChVector<>(-0.25*height, 0, 1.8*height), ChVector<>(0.25*height, 0.5*width, 0.8*height), ChVector<>(0, 0, -1));
  //gen.createObjectsBox(utils::POISSON_DISK, 2 * r, ChVector<>(-0.25*height, 0, 1.01*height), ChVector<>(0.25*height, 0.5*width, 0.01*height), ChVector<>(0, 0, -1));

  ChVector<> hdims(0.45 * height, 0.45 * width, 0);
  ChVector<> center(-0.5 * height, 0, height);
  ChVector<> vel(0, 0, 0);

  while (gen.getTotalNumBodies() < desired_num_particles) {
    gen.createObjectsBox(utils::POISSON_DISK, 2 * r, center, hdims, vel);
    center.z += 2 * r;
  }

  std::cout << "Number of particles: " << gen.getTotalNumBodies() << std::endl;
}


// -----------------------------------------------------------------------------
// Find the number of particles whose height is below a specified value.
// -----------------------------------------------------------------------------
int getNumParticlesBelowValue(ChSystemParallel* sys, double value)
{
  int numBelowZero = 0;
  for (int i = 0; i < sys->GetNumBodies(); ++i) {
    ChBody* body = (ChBody*) sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z < value)
      numBelowZero++;
  }
  return numBelowZero;
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

  ////((ChCollisionSystemParallel*) msystem->GetCollisionSystem())->setBinsPerAxis(I3(50, 50, 50));
  ((ChCollisionSystemParallel*) msystem->GetCollisionSystem())->setBinsPerAxis(I3(10, 10, 10));

  // Create bodies
  ChBody* insert = CreateMechanism(msystem);
  CreateParticles(msystem);

  // Simulation length
  double time_end = 20.0;
  int out_fps = 60;

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
  ChStreamOutAsciiFile hfile(height_file);
  int currentStage = 1;

  while (time < time_end && currentStage <= 3) {

    // Switch between the direct shear stages
    switch(currentStage) {
      case 1:
        // Stage 1: Fill
        if((time > time_settling_min && CheckSettled(msystem, zero_v)) || time > time_settling_max) 
          currentStage = 2;
        break;
      case 2:
        // Stage 2: Open
        insert->SetPos(insert->GetPos()-ChVector<>(speed*time_step,0,0));
        if(insert->GetPos().x <= -0.5*height-delta-gap) currentStage = 3;
        break;
      case 3:
        // Stage 3: Pour
        if(CheckSettled(msystem, zero_v)) currentStage = 4;
        break;
    }

    // output data
    if (sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", out_folder, out_frame + 1);
      utils::WriteShapesPovray(msystem, filename);

      cout << "------------ Output frame:   " << out_frame << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Stage:          " << currentStage << endl;
      cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
      cout << "             Execution time: " << exec_time << endl;

      sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";
/*
      // Create a checkpoint from the current state.
      if (problem == SETTLING) {
        cout << "             Write checkpoint data " << flush;
        utils::WriteCheckpoint(msystem, checkpoint_file);
        cout << msystem->Get_bodylist()->size() << " bodies" << endl;
      }
*/
      // Save current projectile height.
      if (currentStage != 1) {
        int numBelow = getNumParticlesBelowValue(msystem, 0);
        hfile << time << "  " << numBelow << "\n";
        cout << "             Flow:           " << numBelow << endl;
      }

      out_frame++;
      next_out_frame += out_steps;
      num_contacts = 0;
    }
/*
    if (problem == SETTLING && time > time_settling_min && CheckSettled(msystem, zero_v)) {
      cout << "Granular material settled...  time = " << time << endl;
      break;
    }
*/
    msystem->DoStepDynamics(time_step);

    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
    num_contacts += msystem->GetNcontacts();

    // If requested, output detailed timing information for this step
    if (sim_frame == timing_frame)
      msystem->PrintStepStats();

  }
/*
  // Create a checkpoint from the last state
  if (problem == SETTLING) {
    cout << "Write checkpoint data to " << checkpoint_file;
    utils::WriteCheckpoint(msystem, checkpoint_file);
    cout << "  done.  Wrote " << msystem->Get_bodylist()->size() << " bodies." << endl;
  }
*/
  // Final stats
  cout << "==================================" << endl;
  cout << "Number of bodies:  " << msystem->GetNumBodies() << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}
