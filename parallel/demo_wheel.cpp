#include <stdio.h>
#include <vector>
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
using std::endl;

// =======================================================================

enum ProblemType { SETTLING, SIMULATION };

ProblemType problem = SETTLING;

// =======================================================================
// Global problem definitions

int threads = 8;

// Simulation parameters
double gravity = 9.81;
double time_step_settling = 4e-5;
double time_step_simulation = 1e-5;
double time_settling_min = 0.1;
double time_settling_max = 5;
double time_simulation = 10;

int max_iteration = 20;

// Output
const std::string out_dir = "../WHEEL";
const std::string pov_dir = out_dir + "/POVRAY";
const std::string checkpoint_file = out_dir + "/settled.dat";
double out_fps = 60;

// Parameters for the granular material
int Id_g = 100;
double r_g = 0.02;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

float Y_g = 2e8;
float mu_g = 0.5;
float cr_g = 0.1;
float cohesion_g = 20;

// Parameters for the wheel
const std::string obj_mesh_file("../WHEEL/wheel.obj");
const std::string mesh_name("wheel");
const std::string pov_mesh_file("../WHEEL/wheel.inc");

int Id_w = 0;
double mass_w = 600;  ////60;
ChVector<> inertia_w = ChVector<>(1.85, 1.85, 3.675);

float Y_w = 1e8;
float mu_w = 1.0;
float cr_w = 0.1;
float cohesion_w = 20;

// Parameters for the containing bin
int binId = -200;
double hDimX = 4.0;        // length in x direction
double hDimY = 1.0;        // width in y direction
double hDimZ = 0.5;        // height in z direction
double hThickness = 0.04;  // wall thickness

float Y_c = 2e6;
float mu_c = 1.0;
float cr_c = 0.1;
float cohesion_c = 0;

// Height of layer for generator domain
double layerHeight = 1.0;

// =======================================================================

int CreateObjects(ChSystemParallel* system) {
  // Create a material for the granular material
  ChSharedPtr<ChMaterialSurfaceDEM> mat_g;
  mat_g = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_g->SetYoungModulus(Y_g);
  mat_g->SetFriction(mu_g);
  mat_g->SetRestitution(cr_g);
  mat_g->SetCohesion(cohesion_g);

  // Create a material for the container
  ChSharedPtr<ChMaterialSurfaceDEM> mat_c;
  mat_c = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_c->SetYoungModulus(Y_c);
  mat_c->SetFriction(mu_c);
  mat_c->SetRestitution(cr_c);
  mat_c->SetCohesion(cohesion_c);

  // Create a mixture entirely made out of spheres
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
  m1->setDefaultMaterialDEM(mat_g);
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  gen.setBodyIdentifier(Id_g);

  double r = 1.01 * r_g;
  gen.createObjectsBox(utils::POISSON_DISK,
                       2 * r,
                       ChVector<>(0, 0, r + layerHeight / 2),
                       ChVector<>(hDimX - r, hDimY - r, layerHeight / 2));
  cout << "total granules: " << gen.getTotalNumBodies() << endl;

  // Create the containing bin
  utils::CreateBoxContainer(system, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);

  return gen.getTotalNumBodies();
}

// =======================================================================
// Create the wheel body at the specified height.

ChSharedPtr<ChBody> CreateWheel(ChSystemParallel* system, double z) {
  // Create a material for the wheel
  ChSharedPtr<ChMaterialSurfaceDEM> mat_w;
  mat_w = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_w->SetYoungModulus(Y_w);
  mat_w->SetFriction(mu_w);
  mat_w->SetRestitution(cr_w);
  mat_w->SetCohesion(cohesion_w);

  // Create the wheel body
  ChSharedPtr<ChBody> wheel(new ChBody(new ChCollisionModelParallel, ChBody::DEM));

  wheel->SetMaterialSurface(mat_w);

  wheel->SetIdentifier(Id_w);
  wheel->SetMass(mass_w);
  wheel->SetInertiaXX(inertia_w);
  wheel->SetPos(ChVector<>(0, 0, z));
  wheel->SetRot(ChQuaternion<>(1, 0, 0, 0));
  wheel->SetCollide(true);
  wheel->SetBodyFixed(false);

  wheel->GetCollisionModel()->ClearModel();
  utils::AddTriangleMeshGeometry(wheel.get_ptr(), obj_mesh_file, mesh_name);
  wheel->GetCollisionModel()->BuildModel();

  wheel->SetInertiaXX(inertia_w);

  system->AddBody(wheel);

  // Write POV-Ray mesh model.
  utils::WriteMeshPovray(obj_mesh_file, mesh_name, pov_mesh_file);

  return wheel;
}

// ========================================================================
// This utility function returns true if all bodies in the granular mix
// have a linear velocity whose magnitude is below the specified value.

bool CheckSettled(ChSystem* sys, double threshold) {
  double t2 = threshold * threshold;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() >= Id_g) {
      double vel2 = body->GetPos_dt().Length2();
      if (vel2 > t2)
        return false;
    }
  }

  return true;
}

// ========================================================================
// These utility functions find the height of the highest or lowest sphere
// in the granular mix, respectively.  We only look at bodies whith
// identifiers larger than Id_g.

double FindHighest(ChSystem* sys) {
  double highest = 0;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() >= Id_g && body->GetPos().z > highest)
      highest = body->GetPos().z;
  }
  return highest;
}

double FindLowest(ChSystem* sys) {
  double lowest = DBL_MAX;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() >= Id_g && body->GetPos().z < lowest)
      lowest = body->GetPos().z;
  }
  return lowest;
}

// ========================================================================
int main(int argc, char* argv[]) {
  // Create system
  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();

  // Set number of threads.
  int max_threads = omp_get_num_procs();
  if (threads > max_threads)
    threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);

  // Set gravitational acceleration
  msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

  // Edit system settings
  msystem->GetSettings()->solver.tolerance = 1e-3;

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);

  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;

  // Create output directories.
  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

  // Depending on problem type:
  // - Select time step
  // - Select end simulation times
  // - Create granular material and container
  // - Create wheel
  double time_step;
  double time_end;
  ChSharedPtr<ChBody> wheel;

  switch (problem) {
    case SETTLING:
      time_step = time_step_settling;
      time_end = time_settling_max;

      // Create containing bin and the granular material at randomized initial positions
      CreateObjects(msystem);

      break;

    case SIMULATION:
      time_step = time_step_simulation;
      time_end = time_simulation;

      // Create the granular material bodies and the container from the checkpoint file.
      cout << "Read checkpoint data from " << checkpoint_file;
      utils::ReadCheckpoint(msystem, checkpoint_file);
      cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;

      // Create the wheel.
      double z = FindHighest(msystem);
      wheel = CreateWheel(msystem, z + r_g + 0.4);

      break;
  }

  msystem->SetStep(time_step);

  // Number of steps
  int num_steps = std::ceil(time_end / time_step);
  int out_steps = std::ceil((1 / time_step) / out_fps);

  // Zero velocity level for settling check (fraction of a grain radius per second)
  double zero_v = 0.9 * r_g;

  // Perform the simulation
  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  int next_out_frame = 0;
  double exec_time = 0;

  while (time < time_end) {
    if (sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame);
      utils::WriteShapesPovray(msystem, filename);
      utils::WriteCheckpoint(msystem, checkpoint_file);

      cout << " --------------------------------- Output frame:   " << out_frame << endl;
      cout << "                                   Sim frame:      " << sim_frame << endl;
      cout << "                                   Time:           " << time << endl;
      cout << "                                   Execution time: " << exec_time << endl;

      // Check if already settled.
      if (problem == SETTLING && time > time_settling_min && CheckSettled(msystem, zero_v)) {
        cout << "Granular material settled...  time = " << time << endl;
        break;
      }

      out_frame++;
      next_out_frame += out_steps;
    }

    msystem->DoStepDynamics(time_step);

    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
  }

  // Create a checkpoint from the last state
  if (problem == SETTLING)
    utils::WriteCheckpoint(msystem, checkpoint_file);

  // Final stats
  cout << "==================================" << endl;
  cout << "Number of bodies: " << msystem->Get_bodylist()->size() << endl;
  cout << "Simulation time: " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}
