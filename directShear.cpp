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
// Author: Daniel Melanz, Radu Serban
// =============================================================================
//
// ChronoParallel demo program for shearing studies.
//
// The system contains a shearing box composed of three bodies: (1) a containing
// bin, (2) a shearing plate, and (3) a load body. Granular material sits inside
// of the containing bin and shearing plate and is compressed from the top by
// the load body. During the shearing mode, the shear plate is translated in the
// x-direction at a specified velocity.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <stdio.h>
#include <vector>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "chrono_parallel/ChSystemParallel.h"
#include "chrono_parallel/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use DVI contact
#define DEM

// Comment the following line to use parallel collision detection
//#define BULLET

enum ProblemType {
  SETTLING,
  SHEARING
};

ProblemType problem = SETTLING;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Perform shearing action via constraints?
bool kinematic_toggle = false;

// Simulation parameters
double gravity = 981;

double time_settling_min = 0.1;
double time_settling_max = 0.2;
double time_shearing = 10;//0.06;

#ifdef DEM
double time_step = 1e-5;
int max_iteration = 200;
#else
double time_step = 1e-4;
int max_iteration_normal = 0;
int max_iteration_sliding = 200;
int max_iteration_spinning = 0;
float contact_recovery_speed = 10e30;
#endif

int max_iteration_bilateral = 0;
double tolerance = 1e-4;

// Output
#ifdef DEM
const std::string out_dir = "../DIRECTSHEAR_DEM/";
#else
const std::string out_dir = "../DIRECTSHEAR_DVI/";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string shear_file = out_dir + "/shear.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string checkpoint_file = out_dir + "/settled.dat";

int out_fps_settling = 120;
int out_fps_shearing = 120;

int timing_frame = 10;   // output detailed step timing at this frame

// Parameters for the granular material
int        Id_g = 1;
double     r_g = .3;//0.2;//4e-3 / 2;
double     rho_g = 2.500;
double     vol_g = (4.0/3) * CH_C_PI * r_g * r_g * r_g;
double     mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1,1,1);
double     desiredBulkDensity = 1.3894; // 1.3894 g/cm^3

float      Y_g = 1e8;
float      mu_g = 0.5;
float      alpha_g = 0;

// Parameters for the shear box
int        Id_b = 0;
double     hDimX_plate = 6*0.5;                     // length in x direction
double     hDimY_plate = 6*0.5;                     // depth in y direction
double     hDimZ_plate = 3*0.5;                     // height in z direction
double     area = hDimX_plate*2*hDimY_plate*2;
double     rho_b = .700;
double     normalPressure = 168881;                 // apply normal force to cieling// 16,888.1 Pa // 44,127.0 Pa// 71,365.9 Pa
double     mass_b = normalPressure*area / gravity;    // mass of the top body
ChVector<> inertia_b = ChVector<>(1,1,1);           // inertia of shear box
double     desiredVelocity = 0.066;                 // desired shearing velocity

float      Y_b = 1e8;
float      mu_b = 0;
float      alpha_b = 0.1f;

// Parameters for the containing bin
int        binId = -200;
double     hDimX = hDimX_plate;            // length in x direction
double     hDimY = hDimY_plate;            // depth in y direction
double     hDimZ = 3*0.5;                  // height in z direction
double     hThickness = 3;                 // wall thickness

float      Y_c = 2e7;
float      mu_c = 0;//0.3;
float      alpha_c = 0.1f;

// Number of layers and height of one layer for generator domain
int        numLayers = 1;
int        scaleCon = 8;
double     layerHeight = scaleCon*hDimZ/numLayers-r_g;//3e-3;

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

  utils::CreateBoxContainerDEM(system, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness/2);
#else
  ChSharedPtr<ChMaterialSurface> mat_c(new ChMaterialSurface);
  mat_c->SetFriction(mu_c);

  utils::CreateBoxContainerDVI(system, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness/2);
#endif

  return gen.getTotalNumBodies();
}


// -----------------------------------------------------------------------------
// Create the load body such that its bottom point is at the specified height
// and its downward initial velocity has the specified magnitude.
// -----------------------------------------------------------------------------
ChBody* CreateLoadBody(ChSystemParallel* system, double z)
{
  // Create a material for the body
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_b;
  mat_b = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_b->SetYoungModulus(Y_b);
  mat_b->SetFriction(mu_b);
  mat_b->SetDissipationFactor(alpha_b);
#else
  ChSharedPtr<ChMaterialSurface> mat_b(new ChMaterialSurface);
  mat_b->SetFriction(mu_b);
#endif

  // Create the load body
#ifdef DEM
#ifdef BULLET
  ChSharedBodyDEMPtr loadBody(new ChBodyDEM);
#else
  ChSharedBodyDEMPtr loadBody(new ChBodyDEM(new ChCollisionModelParallel));
#endif
  loadBody->SetMaterialSurfaceDEM(mat_b);
#else
#ifdef BULLET
  ChSharedBodyPtr loadBody(new ChBody);
#else
  ChSharedBodyPtr loadBody(new ChBody(new ChCollisionModelParallel));
#endif
  loadBody->SetMaterialSurface(mat_b);
#endif

  loadBody->SetIdentifier(Id_b);
  loadBody->SetMass(mass_b);
  loadBody->SetInertiaXX(inertia_b);
  loadBody->SetPos(ChVector<>(0, 0, z + 2*r_g + hDimZ_plate));
  loadBody->SetRot(ChQuaternion<>(1, 0, 0, 0));
  loadBody->SetCollide(true);
  loadBody->SetBodyFixed(false);

  loadBody->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(loadBody.get_ptr(), ChVector<>(hDimX_plate, hDimY_plate, hDimZ_plate));
  loadBody->GetCollisionModel()->BuildModel();

  system->AddBody(loadBody);

  return loadBody.get_ptr();
}

void AddWall(ChBody*              body,
             const ChVector<>&    loc,
             const ChVector<>&    hdim)
{
  // Append to collision geometry
  body->GetCollisionModel()->AddBox(hdim.x, hdim.y, hdim.z, loc);

  // Append to assets
  ChSharedPtr<ChBoxShape> box_shape(new ChBoxShape);
  box_shape->Pos = loc;
  box_shape->Rot = ChQuaternion<>(1,0,0,0);
  box_shape->GetBoxGeometry().Size = hdim;

  body->GetAssets().push_back(box_shape);
}

ChSharedPtr<ChBody> CreateShearPlate(ChSystemParallel* system)
{
  // Create a material for the plate
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_b;
  mat_b = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_b->SetYoungModulus(1e8f);
  mat_b->SetFriction(mu_c);
  mat_b->SetDissipationFactor(0.1f);
#else
  ChSharedPtr<ChMaterialSurface> mat_b(new ChMaterialSurface);
  mat_b->SetFriction(mu_c);
#endif

  // Create the body
#ifdef DEM
#ifdef BULLET
  ChSharedBodyDEMPtr plate(new ChBodyDEM);
#else
  ChSharedBodyDEMPtr plate(new ChBodyDEM(new ChCollisionModelParallel));
#endif
  plate->SetMaterialSurfaceDEM(mat_b);
#else
#ifdef BULLET
  ChSharedBodyPtr plate(new ChBody);
#else
  ChSharedBodyPtr plate(new ChBody(new ChCollisionModelParallel));
#endif
  plate->SetMaterialSurface(mat_b);
#endif

  plate->SetIdentifier(-11);
  plate->SetMass(1);
  plate->SetInertiaXX(inertia_b);
  plate->SetPos(ChVector<>(0, 0, 2*hDimZ));
  plate->SetRot(ChQuaternion<>(1, 0, 0, 0));
  plate->SetPos_dt(ChVector<>(0, 0, 0));
  plate->SetCollide(true);
  plate->SetBodyFixed(true);

  plate->GetCollisionModel()->ClearModel();
  AddWall(plate.get_ptr(), ChVector<>(-hDimX-hThickness, 0, 5*hDimZ), ChVector<>(hThickness, hDimY, 5*hDimZ));
  AddWall(plate.get_ptr(), ChVector<>( hDimX+hThickness, 0, 5*hDimZ), ChVector<>(hThickness, hDimY, 5*hDimZ));
  AddWall(plate.get_ptr(), ChVector<>(0, -hDimY-hThickness, 5*hDimZ), ChVector<>(hDimX, hThickness, 5*hDimZ));
  AddWall(plate.get_ptr(), ChVector<>(0,  hDimY+hThickness, 5*hDimZ), ChVector<>(hDimX, hThickness, 5*hDimZ));

  plate->GetCollisionModel()->BuildModel();

  system->AddBody(plate);

  return plate;
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
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z > highest)
      if ((body->GetPos().x <= hDimX_plate || body->GetPos().x >= -hDimX_plate) && (body->GetPos().y <= hDimY_plate || body->GetPos().y >= -hDimY_plate))
      {
      highest = body->GetPos().z;
      }
  }
  return highest;
}

double FindLowest(ChSystem* sys)
{
  double lowest = 1000;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*) sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z < lowest)
      lowest = body->GetPos().z;
  }
  return lowest;
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

int setBulkDensity(ChSystem* sys, double bulkDensity)
{
  double normalPlateHeight = sys->Get_bodylist()->at(1)->GetPos().z - hDimZ_plate;
  double bottomHeight = 0;
  int numBodies = sys->Get_bodylist()->size();
  double boxVolume = hDimX_plate * 2 * hDimX_plate * 2 * (normalPlateHeight - bottomHeight);
  double granularVolume = (numBodies - 3)*vol_g;
  double reqDensity = bulkDensity*boxVolume / granularVolume;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 1) {
      body->SetMass(reqDensity*vol_g);
    }
  }
  cout << "N Bodies: " << numBodies << endl;
  cout << "Box Volume: " << boxVolume << endl;
  cout << "Granular Volume: " << granularVolume << endl;
  cout << "Desired bulk density = " << bulkDensity << ", Required Body Density = " << reqDensity << endl;

  return 0;
}


// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  // Create output directories.
  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if(ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

  // Create system
#ifdef DEM
  cout << "Create DEM system" << endl;
  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();
#else
  cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
#endif

  msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

  // Set number of threads.
  int max_threads = msystem->GetParallelThreadNumber();
  if (threads > max_threads) threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;

  msystem->GetSettings()->max_threads = threads;
  msystem->GetSettings()->perform_thread_tuning = thread_tuning;

  // Edit system settings
  msystem->SetTol(tolerance);
  msystem->SetTolSpeeds(tolerance);
  msystem->SetStep(time_step);

#ifdef BULLET
  ChCollisionSystemBulletParallel * bullet_coll = new ChCollisionSystemBulletParallel();
  bullet_coll->Clear();
  msystem->ChangeCollisionSystem(bullet_coll);
#endif

  msystem->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
  msystem->GetSettings()->solver.tolerance = tolerance;

#ifdef DEM
  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;
#else
  msystem->GetSettings()->solver.solver_mode = SLIDING;
  msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
  msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
  msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
  msystem->GetSettings()->solver.alpha = 0;
  msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  msystem->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  msystem->ChangeSolverType(APGDBLAZE);

  msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
  msystem->GetSettings()->collision.min_body_per_bin = 50;
  msystem->GetSettings()->collision.max_body_per_bin = 100;

  // Depending on problem type:
  // - Select end simulation time
  // - Select output FPS
  // - Create granular material and container
  // - Create falling ball
  double time_end;
  int out_fps;
  ChBody* loadBody;
  ChSharedPtr<ChBody> shearPlate;
  ChLinkLockLock* translational;

  if (problem == SETTLING) {
    time_end = time_settling_max;
    out_fps = out_fps_settling;

    // Create the shear plate
    shearPlate = CreateShearPlate(msystem);

    // Create the load body just above the granular material
    double z = scaleCon * hDimZ_plate;
    cout << "Create load body with center at  " << z + 2*r_g + hDimZ_plate << endl;
    loadBody = CreateLoadBody(msystem, z);

    // Create granular material and containing bin
    cout << "Create granular material" << endl;
    CreateObjects(msystem);

  } else {
    time_end = time_shearing;
    out_fps = out_fps_shearing;

    // Create the granular material and the container from the checkpoint file.
    cout << "Read checkpoint data from " << checkpoint_file;
    utils::ReadCheckpoint(msystem, checkpoint_file);
    cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;

    shearPlate = (ChSharedPtr<ChBody>) msystem->Get_bodylist()->at(0);
    loadBody = (ChBody*) msystem->Get_bodylist()->at(1);
    loadBody->SetBodyFixed(false);
    setBulkDensity(msystem, desiredBulkDensity);

    if(kinematic_toggle)
    {
#ifdef BULLET
      ChSharedPtr<ChBody> ground(new ChBody());
#else
      ChSharedPtr<ChBody> ground(new ChBody(new ChCollisionModelParallel));
#endif
      ground->SetBodyFixed(true);
      msystem->AddBody(ground);

      ChVector<> topCM = shearPlate->GetPos();
      shearPlate->SetBodyFixed(false);
      translational = new ChLinkLockLock();
      translational->Initialize(ground, shearPlate, ChCoordsys<>(topCM, QUNIT));
      msystem->AddLink(translational);

      // apply motion
      ChFunction_Ramp* motionFunc = new ChFunction_Ramp(0, desiredVelocity);
      translational->SetMotion_X(motionFunc);
    }
  }

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
  ChStreamOutAsciiFile sfile(stats_file.c_str());
  ChStreamOutAsciiFile shearStream(shear_file.c_str());

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow &gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "Direct Shear Test", msystem);
  gl_window.SetCamera(ChVector<>(0,-10,0), ChVector<>(0,0,0),ChVector<>(0,0,1));
#endif


  while (time < time_end) {
    if (sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
      utils::WriteShapesPovray(msystem, filename, false);

      cout << "------------ Output frame:   " << out_frame << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Lowest point:   " << FindLowest(msystem) << endl;
      cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
      cout << "             Execution time: " << exec_time << endl;

      sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";

      // Create a checkpoint from the current state.
      if (problem == SETTLING) {
        cout << "             Write checkpoint data " << flush;
        utils::WriteCheckpoint(msystem, checkpoint_file);
        cout << msystem->Get_bodylist()->size() << " bodies" << endl;
      }

      out_frame++;
      next_out_frame += out_steps;
      num_contacts = 0;
    }

    if (problem == SETTLING && time > time_settling_min && CheckSettled(msystem, zero_v)) {
      cout << "Granular material settled...  time = " << time << endl;
      break;
    }

    double vx_old = 0;
    double vx_new = 0;
    double px_old = 0;
    double px_new = 0;
    if (problem == SHEARING) {
      vx_old = shearPlate->GetPos_dt().x;
      px_old = shearPlate->GetPos().x;
      shearPlate->SetBodyFixed(false);
    }

#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(time_step);
      gl_window.Render();
    }
#else
    msystem->DoStepDynamics(time_step);
#endif

    if (problem == SHEARING) {
      vx_new = shearPlate->GetPos_dt().x;
      px_new = shearPlate->GetPos().x;

      double FconX = (vx_new - vx_old) / time_step;

      if (kinematic_toggle)
      {
        FconX = translational->Get_react_force().x;
      }
      else
      {
        shearPlate->SetBodyFixed(true);
        vx_new = desiredVelocity;
        px_new = px_old + vx_new*time_step;
        shearPlate->SetPos(ChVector<>(px_new, 0, 2 * hDimZ));
        shearPlate->SetPos_dt(ChVector<>(vx_new, 0, 0));
        shearPlate->SetRot(ChQuaternion<>(1, 0, 0, 0));
        shearPlate->SetRot_dt(ChQuaternion<>(1, 0, 0, 0));
      }
      shearStream << time << ", " << px_new << ", " << FconX << ", " << loadBody->GetPos().z << ", \n";
      std::cout << "Pos: " << px_new << " FconX: " << FconX << endl;
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
  cout << "Number of bodies:  " << msystem->Get_bodylist()->size() << endl;
  cout << "Lowest position:   " << FindLowest(msystem) << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}

