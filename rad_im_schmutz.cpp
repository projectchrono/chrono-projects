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
// Authors: Radu Serban
// =============================================================================
//
// =============================================================================

#include <stdio.h>
#include <vector>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_parallel/ChSystemParallel.h"
#include "chrono_parallel/ChLcpSystemDescriptorParallel.h"

#include "chrono_parallel/collision/ChCNarrowphaseRUtils.h"

#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::endl;
using std::flush;

// -----------------------------------------------------------------------------
// Problem setup
// -----------------------------------------------------------------------------

// Comment the following line to use DVI contact
////#define DEM

// Simulation phase
enum ProblemType {
  SETTLING,
  PUSHING,
  TESTING
};

ProblemType problem = TESTING;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 100;

// Perform dynamic tuning of number of threads?
bool thread_tuning = true;

// Simulation duration.
double time_settling = 5;
double time_pushing = 2;

// Solver parameters
int max_iteration_bilateral = 50;
#ifdef DEM
double time_step = 5e-4;
#else
double time_step = 5e-4;
int max_iteration_normal = 30;
int max_iteration_sliding = 20;
int max_iteration_spinning = 0;
float contact_recovery_speed = 0.1;
#endif

// Output
#ifdef DEM
const std::string out_dir = "../SCHMUTZ_DEM";
#else
const std::string out_dir = "../SCHMUTZ_DVI";
#endif
const std::string pov_dir = out_dir + "/POVRAY";
const std::string checkpoint_file = out_dir + "/settled.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string results_file = out_dir + "/results.dat";

int out_fps_settling = 30;
int out_fps_pushing = 60;

// -----------------------------------------------------------------------------
// Parameters for the granular material (identical spheres)
// -----------------------------------------------------------------------------
double r_g = 0.01;
double rho_g = 2700;
int    desired_num_particles = 10000;

// -----------------------------------------------------------------------------
// Parameters for the test rig
// -----------------------------------------------------------------------------
double mass1 = 895;
double mass_wheel = 661;

ChVector<> inertia_sled(1, 1, 1);
ChVector<> inertia_wheel(1, 1, 1);

double a = 0.987;
double b = 0.112;
double c = 0.5;

double e = 1;

double r_w = 0.3;
double w_w = 0.1;
double s_w = 0.03;

double L = 4;
double W = 0.8;
double H = 0.4;

double d = L / 2 - 2 * w_w;

double init_vel = 5;
double init_angle = (CH_C_PI / 180) * 4;


// =============================================================================
// This class encapsulates the rig's mechanism
// =============================================================================
class Mechanism {
public:
  Mechanism(ChSystemParallel* system, double h);

  void WriteReactionForce(ChStreamOutAsciiFile& f, double time);

private:
  ChVector<> calcLocationWheel(double h)
  {
    double ca = std::cos(init_angle);
    double sa = std::sin(init_angle);

    return ChVector<>(-d - ca * (c + w_w / 2) + sa * (b + r_w),
                                 0,
                      h + sa * (c + w_w / 2) + ca * (b + r_w));
  }

  ChVector<> calcLocationRevolute(double h)
  {
    double ca = std::cos(init_angle);
    double sa = std::sin(init_angle);

    return ChVector<>(-d - ca * (a + c + w_w / 2) + sa * r_w,
                                 0,
                      h + sa * (a + c + w_w / 2) + ca * r_w);
  }

  ChSharedPtr<ChBody> m_ground;
  ChSharedPtr<ChBody> m_sled;
  ChSharedPtr<ChBody> m_wheel;

  ChSharedPtr<ChLinkLockPrismatic> m_prismatic;
  ChSharedPtr<ChLinkLockRevolute> m_revolute;
};

Mechanism::Mechanism(ChSystemParallel* system, double h)
{
  // Calculate hardpoint locations at initial configuration (expressed in the
  // global frame)
  ChVector<> loc_revolute = calcLocationRevolute(h);
  ChVector<> loc_wheel = calcLocationWheel(h);
  ChVector<> loc_sled = loc_revolute - ChVector<>(e, 0, 0);
  ChVector<> loc_prismatic = loc_sled - ChVector<>(0, 0, e / 4);

  // Create the ground body
#ifdef DEM
  m_ground = ChSharedPtr<ChBodyDEM>(new ChBodyDEM(new ChCollisionModelParallel));
#else
  m_ground = ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel));
#endif
  m_ground->SetIdentifier(-1);
  m_ground->SetBodyFixed(true);
  m_ground->SetCollide(false);

  system->AddBody(m_ground);

  // Create the sled body
#ifdef DEM
  m_sled = ChSharedPtr<ChBodyDEM>(new ChBodyDEM(new ChCollisionModelParallel));
#else
  m_sled = ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel));
#endif
  m_sled->SetIdentifier(1);
  m_sled->SetMass(mass1);
  m_sled->SetInertiaXX(inertia_sled);
  m_sled->SetPos(loc_sled);
  m_sled->SetPos_dt(ChVector<>(init_vel, 0, 0));
  m_sled->SetBodyFixed(false);
  m_sled->SetCollide(false);

  ChSharedPtr<ChBoxShape> box_sled(new ChBoxShape);
  box_sled->GetBoxGeometry().Size = ChVector<>(e, e / 3, e / 3);
  box_sled->Pos = ChVector<>(0, 0, 0);
  box_sled->Rot = ChQuaternion<>(1, 0, 0, 0);
  m_sled->AddAsset(box_sled);

  system->AddBody(m_sled);

  // Create the wheel body
#ifdef DEM
  m_wheel = ChSharedPtr<ChBodyDEM>(new ChBodyDEM(new ChCollisionModelParallel));
#else
  m_wheel = ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel));
#endif
  m_wheel->SetIdentifier(2);
  m_wheel->SetMass(mass_wheel);
  m_wheel->SetInertiaXX(inertia_wheel);
  m_wheel->SetPos(loc_wheel);
  m_wheel->SetRot(Q_from_AngY(init_angle));
  m_wheel->SetBodyFixed(false);
  m_wheel->SetCollide(true);

  m_wheel->GetCollisionModel()->ClearModel();
  utils::AddRoundedCylinderGeometry(m_wheel.get_ptr(), r_w, w_w / 2, s_w, ChVector<>(c, 0, -b), Q_from_AngZ(CH_C_PI_2));
  m_wheel->GetCollisionModel()->BuildModel();

  ChSharedPtr<ChCapsuleShape> cap_wheel(new ChCapsuleShape);
  cap_wheel->GetCapsuleGeometry().hlen = (a + c) / 2 - w_w / 2;
  cap_wheel->GetCapsuleGeometry().rad = w_w / 2;
  cap_wheel->Pos = ChVector<>((c - a) / 2, 0, -b);
  cap_wheel->Rot = Q_from_AngZ(CH_C_PI_2);
  m_wheel->AddAsset(cap_wheel);

  system->AddBody(m_wheel);

  // Create and initialize translational joint ground - sled
  m_prismatic = ChSharedPtr<ChLinkLockPrismatic>(new ChLinkLockPrismatic);
  m_prismatic->Initialize(m_sled, m_ground, ChCoordsys<>(loc_prismatic, Q_from_AngY(CH_C_PI_2)));
  system->AddLink(m_prismatic);

  // Create and initialize revolute joint sled - wheel
  m_revolute = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
  m_revolute->Initialize(m_wheel, m_sled, ChCoordsys<>(loc_revolute, Q_from_AngX(CH_C_PI_2)));
  system->AddLink(m_revolute);
}

void Mechanism::WriteReactionForce(ChStreamOutAsciiFile& f, double time)
{
  // Coordinate system of the revolute joint (relative to the frame of body2, in
  // this case the sled)
  ChCoordsys<> revCoordsys = m_revolute->GetLinkRelativeCoords();

  ChVector<> force_jointsys = m_revolute->Get_react_force();
  ChVector<> force_bodysys = revCoordsys.TransformDirectionLocalToParent(force_jointsys);
  ChVector<> force_abssys = m_sled->GetCoord().TransformDirectionLocalToParent(force_bodysys);

  f << time << "  "
    << force_jointsys.x << "  " << force_jointsys.y << "  " << force_jointsys.z << "  "
    << force_bodysys.x << "  " << force_bodysys.y << "  " << force_bodysys.z << "  "
    << force_abssys.x << "  " << force_abssys.y << "  " << force_abssys.z << "\n";
}


// =============================================================================
// Create container bin.
// =============================================================================
void CreateContainer(ChSystemParallel* system)
{
  int    id_c = -200;
  double thickness = 0.2;

#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_c;
  mat_c = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_c->SetYoungModulus(2e6f);
  mat_c->SetFriction(0.4f);
  mat_c->SetDissipationFactor(0.6f);

  utils::CreateBoxContainerDEM(system, id_c, mat_c, ChVector<>(L/2, W/2, H/2), thickness/2);
#else
  ChSharedPtr<ChMaterialSurface> mat_c(new ChMaterialSurface);
  mat_c->SetFriction(0.4f);

  utils::CreateBoxContainerDVI(system, id_c, mat_c, ChVector<>(L/2, W/2, H/2), thickness/2);

#endif
}


// =============================================================================
// Create granular material.
// =============================================================================
void CreateParticles(ChSystemParallel* system)
{
  // Create a material for the ball mixture.
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_g;
  mat_g = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_g->SetYoungModulus(1e8f);
  mat_g->SetFriction(0.4f);
  mat_g->SetDissipationFactor(0.1f);
#else
  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
  mat_g->SetFriction(0.4f);
#endif

  // Create a mixture entirely made out of spheres.
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
#ifdef DEM
  m1->setDefaultMaterialDEM(mat_g);
#else
  m1->setDefaultMaterialDVI(mat_g);
#endif
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  // Create particles, one layer at a time, until the desired number is reached.
  gen.setBodyIdentifier(100);

  double     r = 1.01 * r_g;
  ChVector<> hdims(L/2 - r, W/2 - r, 0);
  ChVector<> center(0, 0, 2 * r);

  int layer = 1;
  while(gen.getTotalNumBodies() < desired_num_particles) {
    gen.createObjectsBox(utils::POISSON_DISK, 2 * r, center, hdims);
    cout << "layer " << layer << "    total particles: " << gen.getTotalNumBodies() << endl;
    center.z += 2 * r;
    layer++;
  }
}


// =============================================================================
// Find the height of the highest and lowest sphere in the granular mix.
// We only look at bodies whith identifiers larger than 100 (to exclude all
// other bodies).
// =============================================================================
void FindRange(ChSystem* sys, double& lowest, double& highest)
{
  highest = -1000;
  lowest = 1000;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() < 100)
      continue;
    double h = body->GetPos().z;
    if (h < lowest) lowest = h;
    if (h > highest) highest = h;
  }
}

// =============================================================================
// =============================================================================
int main(int argc, char* argv[])
{
  // --------------------------
  // Create output directories.
  // --------------------------

  if(ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if(ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

  // --------------
  // Create system.
  // --------------

#ifdef DEM
  cout << "Create DEM system" << endl;
  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();
#else
  cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
#endif

  msystem->Set_G_acc(ChVector<>(0, 0, -9.81));

  // ----------------------
  // Set number of threads.
  // ----------------------

  int max_threads = msystem->GetParallelThreadNumber();
  if (threads > max_threads)
    threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;

  // ---------------------
  // Edit system settings.
  // ---------------------

  msystem->SetTol(1e-3);
  msystem->SetTolSpeeds(1e-3);
  msystem->SetStep(time_step);

  msystem->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
  msystem->GetSettings()->solver.tolerance = 1e-3;

#ifdef DEM
  msystem->ChangeCollisionNarrowphase(NARROWPHASE_R);
#else
  msystem->GetSettings()->solver.solver_mode = SLIDING;
  msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
  msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
  msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
  msystem->GetSettings()->solver.alpha = 0;
  msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  msystem->ChangeSolverType(APGDRS);

  msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);

  // ----------------------------------------
  // Depending on problem type:
  // - Select end simulation time
  // - Create granular material and container
  // - Create mechanism
  // ----------------------------------------

  double time_end;
  int out_fps;
  Mechanism* mech = NULL;

  switch (problem) {
  case SETTLING:
    time_end = time_settling;
    out_fps = out_fps_settling;

    cout << "Create granular material" << endl;
    CreateContainer(msystem);
    CreateParticles(msystem);

    break;

  case PUSHING:
  {
    time_end = time_pushing;
    out_fps = out_fps_pushing;

    // Create the granular material and the container from the checkpoint file.
    cout << "Read checkpoint data from " << checkpoint_file;
    utils::ReadCheckpoint(msystem, checkpoint_file);
    cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;

    // Create the mechanism with the wheel just above the granular material.
    double lowest, highest;
    FindRange(msystem, lowest, highest);
    cout << "Create mechanism above height" << highest + r_g << endl;
    mech = new Mechanism(msystem, highest + r_g);
  }

    break;

  case TESTING:
    time_end = time_pushing;
    out_fps = out_fps_pushing;

    mech = new Mechanism(msystem, 0.9 * H);

    break;
  }

  // Number of steps.
  int num_steps = std::ceil(time_end / time_step);
  int out_steps = std::ceil((1.0 / time_step) / out_fps);

  // -----------------------
  // Perform the simulation.
  // -----------------------
  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  int next_out_frame = 0;
  double exec_time = 0;
  int num_contacts = 0;
  ChStreamOutAsciiFile sfile(stats_file.c_str());
  ChStreamOutAsciiFile rfile(results_file.c_str());

  while (time < time_end) {
    if (sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
      utils::WriteShapesPovray(msystem, filename, problem == TESTING);
      double highest, lowest;
      FindRange(msystem, lowest, highest);
      cout << "------------ Output frame:   " << out_frame << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Lowest point:   " << lowest << endl;
      cout << "             Highest point:  " << highest << endl;
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

    // Advance dynamics.
    msystem->DoStepDynamics(time_step);

    // Increment counters
    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
    num_contacts += msystem->GetNcontacts();

    // Save results
    if (problem != SETTLING) {
      mech->WriteReactionForce(rfile, time);
    }
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
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}

