#include <stdio.h>
#include <vector>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"
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
#define USE_DEM

// Simulation phase
enum ProblemType {
  SETTLING,
  DROPPING,
};
ProblemType problem = DROPPING;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 100;

// Perform dynamic tuning of number of threads?
bool thread_tuning = true;

// Simulation duration.
double time_settling = 5;
double time_dropping = 2;

// Solver parameters
#ifdef USE_DEM
double time_step = 1e-4;
int max_iteration = 20;
#else
double time_step = 1e-4;
int max_iteration_normal = 30;
int max_iteration_sliding = 20;
int max_iteration_spinning = 0;
float contact_recovery_speed = 0.1;
#endif

// Output
#ifdef USE_DEM
const std::string out_dir = "../SOILBIN_DEM";
#else
const std::string out_dir = "../SOILBIN_DVI";
#endif
const std::string pov_dir = out_dir + "/POVRAY";
const std::string checkpoint_file = out_dir + "/settled.dat";
const std::string stats_file = out_dir + "/stats.dat";

int out_fps_settling = 30;
int out_fps_dropping = 60;

// -----------------------------------------------------------------------------
// Parameters for the granular material (identical spheres)
// -----------------------------------------------------------------------------
double r_g = 0.1;
double rho_g = 2000;
int desired_num_particles = 1000;

// -----------------------------------------------------------------------------
// Parameters for the falling object
// -----------------------------------------------------------------------------
// Shape of dropped object
collision::ShapeType shape_o = collision::ROUNDEDCYL;

ChQuaternion<> initRot(1.0, 0.0, 0.0, 0.0);
ChVector<> initLinVel(0.0, 0.0, 0.0);
ChVector<> initAngVel(0.0, 0.0, 0.0);

// -----------------------------------------------------------------------------
// Half-dimensions of the container bin
// -----------------------------------------------------------------------------
double hDimX = 5;
double hDimY = 2;
double hDimZ = 2;

// =============================================================================
// Create container bin.
// =============================================================================
void CreateContainer(ChSystemParallel* system) {
  int id_c = -200;
  double hThickness = 0.1;

#ifdef USE_DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_c;
  mat_c = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_c->SetYoungModulus(2e6f);
  mat_c->SetFriction(0.4f);
  mat_c->SetRestitution(0.1f);

  utils::CreateBoxContainer(system, id_c, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);
#else
  ChSharedPtr<ChMaterialSurface> mat_c(new ChMaterialSurface);
  mat_c->SetFriction(0.4f);

  utils::CreateBoxContainer(system, id_c, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);

#endif
}

// =============================================================================
// Create granular material.
// =============================================================================
void CreateParticles(ChSystemParallel* system) {
// Create a material for the ball mixture.
#ifdef USE_DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_g;
  mat_g = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_g->SetYoungModulus(1e8f);
  mat_g->SetFriction(0.4f);
  mat_g->SetRestitution(0.1f);
#else
  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
  mat_g->SetFriction(0.4f);
#endif

  // Create a mixture entirely made out of spheres.
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
#ifdef USE_DEM
  m1->setDefaultMaterialDEM(mat_g);
#else
  m1->setDefaultMaterialDVI(mat_g);
#endif
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  // Create particles, one layer at a time, until the desired number is reached.
  gen.setBodyIdentifier(1);

  double r = 1.01 * r_g;
  ChVector<> hdims(hDimX - r, hDimY - r, 0);
  ChVector<> center(0, 0, 2 * r);

  while (gen.getTotalNumBodies() < desired_num_particles) {
    gen.createObjectsBox(utils::POISSON_DISK, 2 * r, center, hdims);
    center.z += 2 * r;
  }

  cout << "Number of particles: " << gen.getTotalNumBodies() << endl;
}

// =============================================================================
// Create falling object.
// =============================================================================
void CreateObject(ChSystemParallel* system, double z) {
  double rho_o = 2000.0;

// -----------------------------------------
// Create a material for the falling object.
// -----------------------------------------

#ifdef USE_DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_o;
  mat_o = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_o->SetYoungModulus(1e8f);
  mat_o->SetFriction(0.4f);
  mat_o->SetRestitution(0.1f);
#else
  ChSharedPtr<ChMaterialSurface> mat_o(new ChMaterialSurface);
  mat_o->SetFriction(0.4f);
#endif

// --------------------------
// Create the falling object.
// --------------------------

#ifdef USE_DEM
  ChSharedPtr<ChBody> obj(new ChBody(new ChCollisionModelParallel, ChBody::DEM));
  obj->SetMaterialSurface(mat_o);
#else
  ChSharedBodyPtr obj(new ChBody(new ChCollisionModelParallel));
  obj->SetMaterialSurface(mat_o);
#endif

  obj->SetIdentifier(0);
  obj->SetCollide(true);
  obj->SetBodyFixed(false);

  // ----------------------------------------------------
  // Depending on the shape of the falling object,
  //    - Calculate bounding radius, volume, and gyration
  //    - Set contact and visualization shape
  // ----------------------------------------------------

  double rb;
  double vol;
  ChMatrix33<> J;

  obj->GetCollisionModel()->ClearModel();

  switch (shape_o) {
    case collision::SPHERE: {
      double radius = 0.5;
      rb = utils::CalcSphereBradius(radius);
      vol = utils::CalcSphereVolume(radius);
      J = utils::CalcSphereGyration(radius);
      utils::AddSphereGeometry(obj.get_ptr(), radius);
    } break;
    case collision::BOX: {
      ChVector<> hdims(0.5, 0.75, 1.0);
      rb = utils::CalcBoxBradius(hdims);
      vol = utils::CalcBoxVolume(hdims);
      J = utils::CalcBoxGyration(hdims);
      utils::AddBoxGeometry(obj.get_ptr(), hdims);
    } break;
    case collision::CAPSULE: {
      double radius = 0.25;
      double hlen = 0.5;
      rb = utils::CalcCapsuleBradius(radius, hlen);
      vol = utils::CalcCapsuleVolume(radius, hlen);
      J = utils::CalcCapsuleGyration(radius, hlen);
      utils::AddCapsuleGeometry(obj.get_ptr(), radius, hlen);
    } break;
    case collision::CYLINDER: {
      double radius = 0.25;
      double hlen = 0.5;
      rb = utils::CalcCylinderBradius(radius, hlen);
      vol = utils::CalcCylinderVolume(radius, hlen);
      J = utils::CalcCylinderGyration(radius, hlen);
      utils::AddCylinderGeometry(obj.get_ptr(), radius, hlen);
    } break;
    case collision::ROUNDEDCYL: {
      double radius = 0.25;
      double hlen = 0.1;
      double srad = 0.1;
      rb = utils::CalcRoundedCylinderBradius(radius, hlen, srad);
      vol = utils::CalcRoundedCylinderVolume(radius, hlen, srad);
      J = utils::CalcRoundedCylinderGyration(radius, hlen, srad);
      utils::AddRoundedCylinderGeometry(obj.get_ptr(), radius, hlen, srad);
    } break;
  }

  obj->GetCollisionModel()->BuildModel();

  // ---------------------
  // Set mass and inertia.
  // ---------------------

  double mass = rho_o * vol;
  obj->SetMass(mass);
  obj->SetInertia(J * mass);

  // ------------------
  // Set initial state.
  // ------------------

  obj->SetPos(ChVector<>(0, 0, z + rb));
  obj->SetRot(initRot);
  obj->SetPos_dt(initLinVel);
  obj->SetWvel_loc(initAngVel);

  // ---------------------
  // Add object to system.
  // ---------------------
  system->AddBody(obj);
}

// =============================================================================
// Find the height of the highest and lowest, respectively, sphere in the
// granular mix, respectively.  We only look at bodies whith stricty positive
// identifiers (to exclude the containing bin).
// =============================================================================
double FindHighest(ChSystem* sys) {
  double highest = 0;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z > highest)
      highest = body->GetPos().z;
  }
  return highest;
}

double FindLowest(ChSystem* sys) {
  double lowest = 1000;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 0 && body->GetPos().z < lowest)
      lowest = body->GetPos().z;
  }
  return lowest;
}

// =============================================================================
// =============================================================================
int main(int argc, char* argv[]) {
  // --------------------------
  // Create output directories.
  // --------------------------

  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

// --------------
// Create system.
// --------------

#ifdef USE_DEM
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

  int max_threads = omp_get_num_procs();
  if (threads > max_threads)
    threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;

  // ---------------------
  // Edit system settings.
  // ---------------------

  msystem->GetSettings()->solver.tolerance = 1e-3;

#ifdef USE_DEM
  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;
#else
  msystem->GetSettings()->solver.solver_mode = SLIDING;
  msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
  msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
  msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
  msystem->GetSettings()->solver.alpha = 0;
  msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  msystem->ChangeSolverType(APGDREF);

  msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);

  // ----------------------------------------
  // Depending on problem type:
  // - Select end simulation time
  // - Create granular material and container
  // - Create falling object
  // ----------------------------------------

  double time_end;
  int out_fps;
  ChBody* ball;

  if (problem == SETTLING) {
    time_end = time_settling;
    out_fps = out_fps_settling;

    cout << "Create granular material" << endl;
    CreateContainer(msystem);
    CreateParticles(msystem);
  } else {
    time_end = time_dropping;
    out_fps = out_fps_dropping;

    // Create the granular material and the container from the checkpoint file.
    cout << "Read checkpoint data from " << checkpoint_file;
    utils::ReadCheckpoint(msystem, checkpoint_file);
    cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;

    // Create the falling object just above the granular material.
    double z = FindHighest(msystem);
    cout << "Create falling object above height" << z + r_g << endl;
    CreateObject(msystem, z + r_g);
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

  while (time < time_end) {
    if (sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
      utils::WriteShapesPovray(msystem, filename);

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

    // Advance dynamics.
    msystem->DoStepDynamics(time_step);

    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
    num_contacts += msystem->GetNcontacts();
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
