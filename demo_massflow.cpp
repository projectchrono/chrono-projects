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
double time_settling_max = 0.8;
double time_dropping = 0.06;

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

int timing_frame = 10;   // output detailed step timing at this frame

// Parameters for the granular material
int        Id_g = 1;
double     r_g = 1e-3 / 2;
double     rho_g = 2500;
double     vol_g = (4.0/3) * CH_C_PI * r_g * r_g * r_g;
double     mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1,1,1);

float      Y_g = 1e8;
float      mu_g = 0.3;
float      alpha_g = 0.1;

// Parameters for the falling ball
int        Id_b = 0;
double     R_b = 2.54e-2 / 2;
double     rho_b = 700;
double     vol_b = (4.0/3) * CH_C_PI * R_b * R_b * R_b;
double     mass_b = rho_b * vol_b;
ChVector<> inertia_b = 0.4 * mass_b * R_b * R_b * ChVector<>(1,1,1);

float      Y_b = 1e8;
float      mu_b = 0.3;
float      alpha_b = 0.1;

// Parameters for the containing bin
int        binId = -200;
double     hDimX = 4e-2;            // length in x direction
double     hDimY = 4e-2;            // depth in y direction
double     hDimZ = 7.5e-2;          // height in z direction
double     hThickness = 0.5e-2;     // wall thickness

float      Y_c = 2e6;
float      mu_c = 0.3;
float      alpha_c = 0.6;

// Number of layers and height of one layer for generator domain
////int        numLayers = 10;
////double     layerHeight = 1e-2;
int        numLayers = 1;
double     layerHeight = 3e-3;

// Drop height (above surface of settled granular material)
double h = 10e-2;

// -----------------------------------------------------------------------------
// Create the dynamic objects:
// - granular material consisting of identical spheres with specified radius and
//   material properties; the spheres are generated in a number of vertical
//   layers with locations within each layer obtained using Poisson Disk
//   sampling (thus ensuring that no two spheres are closer than twice the
//   radius)
// - a containing bin consisting of five boxes (no top)
// -----------------------------------------------------------------------------
int CreateObjects(ChSystemParallel* system)
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

  gen.setBodyIdentifier(Id_g);

  double r = 1.01 * r_g;

  for (int i = 0; i < numLayers; i++) {
    double center = r + layerHeight / 2 + i * (2 * r + layerHeight);
    gen.createObjectsBox(utils::POISSON_DISK,
                         2 * r,
                         ChVector<>(0, 0, center),
                         ChVector<>(hDimX - r, hDimY - r, layerHeight/2));
    cout << "Layer " << i << "  total bodies: " << gen.getTotalNumBodies() << endl;
  }

  // Create the containing bin
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_c;
  mat_c = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_c->SetYoungModulus(Y_c);
  mat_c->SetFriction(mu_c);
  mat_c->SetDissipationFactor(alpha_c);

  utils::CreateBoxContainerDEM(system, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);
#else
  ChSharedPtr<ChMaterialSurface> mat_c(new ChMaterialSurface);
  mat_c->SetFriction(mu_c);

  utils::CreateBoxContainerDVI(system, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);
#endif

  return gen.getTotalNumBodies();
}


// -----------------------------------------------------------------------------
// Create the falling ball such that its bottom point is at the specified height
// and its downward initial velocity has the specified magnitude.
// -----------------------------------------------------------------------------
ChBody* CreateFallingBall(ChSystemParallel* system, double z, double vz)
{
  // Create a material for the falling ball
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_b;
  mat_b = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_b->SetYoungModulus(1e8f);
  mat_b->SetFriction(0.4f);
  mat_b->SetDissipationFactor(0.1f);
#else
  ChSharedPtr<ChMaterialSurface> mat_b(new ChMaterialSurface);
  mat_b->SetFriction(mu_c);
#endif

  // Create the falling ball
#ifdef DEM
  ChSharedBodyDEMPtr ball(new ChBodyDEM(new ChCollisionModelParallel));
  ball->SetMaterialSurfaceDEM(mat_b);
#else
  ChSharedBodyPtr ball(new ChBody(new ChCollisionModelParallel));
  ball->SetMaterialSurface(mat_b);
#endif

  ball->SetIdentifier(Id_b);
  ball->SetMass(mass_b);
  ball->SetInertiaXX(inertia_b);
  ball->SetPos(ChVector<>(0, 0, z + r_g + R_b));
  ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
  ball->SetPos_dt(ChVector<>(0, 0, -vz));
  ball->SetCollide(true);
  ball->SetBodyFixed(false);

  ball->GetCollisionModel()->ClearModel();
  ball->GetCollisionModel()->AddSphere(R_b);
  ball->GetCollisionModel()->BuildModel();

  ChSharedPtr<ChSphereShape> ball_shape = ChSharedPtr<ChAsset>(new ChSphereShape);
  ball_shape->SetColor(ChColor(0, 0, 1));
  ball_shape->GetSphereGeometry().rad = R_b;
  ball_shape->Pos = ChVector<>(0, 0, 0);
  ball_shape->Rot = ChQuaternion<>(1, 0, 0, 0);
  ball->GetAssets().push_back(ball_shape);

  system->AddBody(ball);

  return ball.get_ptr();
}


// -----------------------------------------------------------------------------
// Find the height of the highest and lowest, respectively, sphere in the
// granular mix, respectively.  We only look at bodies whith stricty positive
// identifiers (to exclude the containing bin).
// -----------------------------------------------------------------------------
double FindHighest(ChSystem* sys)
{
  double highest = 0;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*) sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z > highest)
      highest = body->GetPos().z;
  }
  return highest;
}

int getNumParticlesBelowValue(ChSystem* sys, double value)
{
  int numBelowZero = 0;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
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
bool CheckSettled(ChSystem* sys, double threshold)
{
  double t2 = threshold * threshold;

  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
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

  // Create the mass flow meter (TODO: move outside main to reflect the demo_crater example)
  double height = 3.175; // height of the cavity, cm
  double width = 1; // width of the cavity, cm
  double thickness = 1; // thickness of walls, cm	
  double angle = 45; // angle of the cavity, degrees
  double particleRadius = 0.025; // radius of the particles, cm
  double particleDensity = 2.5; // density of the particles, g/cm^3
  double speed = 0.15; // speed of the angled insert, cm/s
  double gap = 0.2; // size of gap, cm

  double time_end = 20.0;
  int out_fps = 60;
  

  // Step 1: The static slot portion
  #ifdef DEM
    ChSharedPtr<ChMaterialSurfaceDEM> mat_b;
    mat_b = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
    mat_b->SetYoungModulus(1e8f);
    mat_b->SetFriction(0.4f);
    mat_b->SetDissipationFactor(0.1f);
  #else
    ChSharedPtr<ChMaterialSurface> mat_b(new ChMaterialSurface);
    mat_b->SetFriction(mu_c);
  #endif

  #ifdef DEM
    ChSharedBodyDEMPtr slot(new ChBodyDEM(new ChCollisionModelParallel));
    slot->SetMaterialSurfaceDEM(mat_b);
  #else
    ChSharedBodyPtr slot(new ChBody(new ChCollisionModelParallel));
    slot->SetMaterialSurface(mat_b);
  #endif

  slot->SetIdentifier(1);
  slot->SetMass(1);
  slot->SetInertiaXX(ChVector<>(1,1,1));
  slot->SetPos(ChVector<>(0.5*thickness, 0, 0.5*height));
  slot->SetRot(ChQuaternion<>(1, 0, 0, 0));
  slot->SetCollide(true);
  slot->SetBodyFixed(true);

  slot->GetCollisionModel()->ClearModel();
  slot->GetCollisionModel()->AddBox(thickness*0.5,width*0.5,height*0.5, ChVector<>(0,0,0));
  slot->GetCollisionModel()->AddBox(3*height*0.5,width*0.5,height, ChVector<>(0,width,0.5*height));
  slot->GetCollisionModel()->AddBox(3*height*0.5,width*0.5,height, ChVector<>(0,-width,0.5*height));
  // additional collision geometry added to prevent particles from falling out
  slot->GetCollisionModel()->AddBox(thickness*0.5,width*0.5,height*0.5, ChVector<>(0,0,height));
  slot->GetCollisionModel()->AddBox(3*height*0.5,width*0.5,height, ChVector<>(0,width,1.5*height));
  slot->GetCollisionModel()->AddBox(3*height*0.5,width*0.5,height, ChVector<>(0,-width,1.5*height));
  slot->GetCollisionModel()->AddBox(thickness*0.5,width*0.5,height*0.5, ChVector<>(0,0,2*height));
  slot->GetCollisionModel()->BuildModel();

  ChSharedPtr<ChBoxShape> slot_shape = ChSharedPtr<ChAsset>(new ChBoxShape);
  slot_shape->SetColor(ChColor(0, 0, 1));
  slot_shape->GetBoxGeometry().Size = ChVector<>(thickness*0.5,width*0.5,height*0.5);
  slot_shape->Pos = ChVector<>(0, 0, 0);
  slot_shape->Rot = ChQuaternion<>(1, 0, 0, 0);
  slot->GetAssets().push_back(slot_shape);

  msystem->AddBody(slot);

  // Step 2: The angled insert portion
  double height_insert = sqrt(2*height*height);
  #ifdef DEM
    ChSharedBodyDEMPtr insert(new ChBodyDEM(new ChCollisionModelParallel));
    insert->SetMaterialSurfaceDEM(mat_b);
  #else
    ChSharedBodyPtr insert(new ChBody(new ChCollisionModelParallel));
    insert->SetMaterialSurface(mat_b);
  #endif

  insert->SetIdentifier(2);
  insert->SetMass(1);
  insert->SetInertiaXX(ChVector<>(1,1,1));
  double delta = sqrt(thickness*thickness/8);
  insert->SetPos(ChVector<>(-0.5*height-delta, 0, 0.5*height-delta));
  insert->SetRot(chrono::Q_from_AngAxis(-CH_C_PI/4, VECT_Y));
  insert->SetCollide(true);
  insert->SetBodyFixed(true);

  insert->GetCollisionModel()->ClearModel();
  insert->GetCollisionModel()->AddBox(thickness*0.5,width*0.5,height_insert*0.5, ChVector<>(0,0,0));
  // additional collision geometry added to prevent particles from falling out
  insert->GetCollisionModel()->AddBox(thickness*0.5,width*0.5,height_insert*0.5, ChVector<>(0,0,height_insert));
  insert->GetCollisionModel()->BuildModel();

  ChSharedPtr<ChBoxShape> insert_shape = ChSharedPtr<ChAsset>(new ChBoxShape);
  insert_shape->SetColor(ChColor(0, 0, 1));
  insert_shape->GetBoxGeometry().Size = ChVector<>(thickness*0.5,width*0.5,height_insert*0.5);
  insert_shape->Pos = ChVector<>(0, 0, 0);
  insert_shape->Rot = ChQuaternion<>(1, 0, 0, 0);
  insert->GetAssets().push_back(insert_shape);

  msystem->AddBody(insert);

  // Step 3: Create the containing bin
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_c;
  mat_c = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_c->SetYoungModulus(Y_c);
  mat_c->SetFriction(mu_c);
  mat_c->SetDissipationFactor(alpha_c);

  utils::CreateBoxContainerDEM(msystem, 3, mat_c, ChVector<>(0.5*height, 0.5*height, 0.1*height), thickness, ChVector<>(0, 0, -1.5*height));
#else
  ChSharedPtr<ChMaterialSurface> mat_c(new ChMaterialSurface);
  mat_c->SetFriction(mu_c);

  utils::CreateBoxContainerDVI(msystem, 3, mat_c, ChVector<>(0.5*height, 0.5*height, 0.1*height), thickness, ChVector<>(0, 0, -1.5*height));
#endif

  // Step 4: Add granular material
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
  utils::Generator gen(msystem);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
  #ifdef DEM
    m1->setDefaultMaterialDEM(mat_g);
  #else
    m1->setDefaultMaterialDVI(mat_g);
  #endif
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(particleRadius);

  gen.setBodyIdentifier(Id_g);

  double r = 1.01 * particleRadius;	
  //gen.createObjectsBox(utils::POISSON_DISK, 2 * r, ChVector<>(-0.25*height, 0, 1.8*height), ChVector<>(0.25*height, 0.5*width, 0.8*height), ChVector<>(0, 0, -1));
  gen.createObjectsBox(utils::POISSON_DISK, 2 * r, ChVector<>(-0.25*height, 0, 1.01*height), ChVector<>(0.25*height, 0.5*width, 0.01*height), ChVector<>(0, 0, -1));

  // Number of steps
  int num_steps = std::ceil(time_end / time_step);
  int out_steps = std::ceil((1.0 / time_step) / out_fps);

  // Zero velocity level for settling check
  // (fraction of a grain radius per second)
  double zero_v = 0.1 * r_g;

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
    switch(currentStage)
    {
      case 1:
        // Stage 1: Fill
        if(CheckSettled(msystem, zero_v)) currentStage = 2;
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
  cout << "Number of bodies:  " << msystem->Get_bodylist()->size() << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}
