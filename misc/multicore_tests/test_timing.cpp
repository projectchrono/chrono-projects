#include <stdio.h>
#include <vector>
#include <cmath>

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
// Global problem definitions

int threads = 8;

// Simulation parameters
double gravity = -9.81;
double time_step = 1e-4;
int num_steps_settle = std::ceil(3 / time_step);
int num_steps_measure = 1000;

int max_iteration = 50;

// Parameters for the mixture balls
double radius = 0.1;
double density = 2000;
double vol = (4.0 / 3) * CH_C_PI * radius * radius * radius;
double mass = density * vol;
ChVector<> inertia = 0.4 * mass * radius * radius * ChVector<>(1, 1, 1);

// Parameters for the containing bin
int binId = -200;
double hDimX = 2;         // length in x direction
double hDimY = 2;         // depth in y direction
double hDimZ = 5;         // height in z direction
double hThickness = 0.1;  // wall thickness

// =======================================================================

enum ProblemType { SETTLE, MEASURE };

ProblemType problem = SETTLE;

const char* checkpoint_file = "../TEST/settled_ckp.dat";
const char* povray_file = "../TEST/settled_povray.dat";

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CreateObjects(ChSystemParallel* system) {
  // Create a material for the ball mixture
  ChSharedPtr<ChMaterialSurfaceDEM> ballMixMat;
  ballMixMat = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  ballMixMat->SetYoungModulus(1e8f);
  ballMixMat->SetFriction(0.4f);
  ballMixMat->SetRestitution(0.4f);

  // Create a mixture entirely made out of spheres
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
  m1->setDefaultMaterialDEM(ballMixMat);
  m1->setDefaultDensity(density);
  m1->setDefaultSize(radius);

  // Generate the objects
  gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * radius, ChVector<>(0, 0, 2.5),
                       ChVector<>(0.8 * hDimX, 0.8 * hDimY, 1));

  // Create the containing bin
  ChSharedPtr<ChMaterialSurfaceDEM> binMat;
  binMat = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  binMat->SetYoungModulus(2e6f);
  binMat->SetFriction(0.4f);
  binMat->SetRestitution(0.1f);

  utils::CreateBoxContainer(system, binId, binMat, ChVector<>(hDimX, hDimY, hDimZ), hThickness);

  return gen.getTotalNumBodies();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  // Create system
  // -------------

  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();

  // Set number of threads.
  int max_threads = omp_get_num_procs();
  if (threads > max_threads)
    threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);

  msystem->GetSettings()->perform_thread_tuning = false;

  // Set gravitational acceleration
  msystem->Set_G_acc(ChVector<>(0, 0, gravity));

  // Edit system settings
  msystem->GetSettings()->solver.tolerance = 1e-3;

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);

  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;

  // Create objects and perform simulation
  // -------------------------------------

  switch (problem) {
    case SETTLE: {
      cout << "Create objects  ";
      int num_obj = CreateObjects(msystem);
      cout << num_obj << endl;

      for (int i = 0; i < num_steps_settle; i++)
        msystem->DoStepDynamics(time_step);

      utils::WriteCheckpoint(msystem, checkpoint_file);
      utils::WriteShapesPovray(msystem, povray_file);
    }

    break;
    case MEASURE: {
      cout << "Read checkpoint data from " << checkpoint_file;
      utils::ReadCheckpoint(msystem, checkpoint_file);
      cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;

      double run_time = 0;
      double contact_time = 0;
      double collision_time = 0;
      double collision_broad_time = 0;
      double collision_narrow_time = 0;
      double update_time = 0;
      double solver_time = 0;
      int num_contacts = 0;

      for (int i = 0; i < num_steps_measure; i++) {
        msystem->DoStepDynamics(time_step);
        run_time += msystem->GetTimerStep();
        update_time += msystem->GetTimerUpdate();
        solver_time += msystem->GetTimerLcp();
        contact_time += msystem->GetTimerProcessContact();
        collision_time += msystem->GetTimerCollision();
        collision_broad_time += msystem->GetTimerCollisionBroad();
        collision_narrow_time += msystem->GetTimerCollisionNarrow();
        num_contacts += msystem->GetNcontacts();
      }

      cout << "==================================" << endl;
      cout << "Number of threads:     " << threads << endl;
      cout << "Number of bodies:      " << msystem->Get_bodylist()->size() << endl;
      cout << "Number of steps:       " << num_steps_measure << endl;
      cout << endl;
      cout << "Average num. contacts: " << (1.0 * num_contacts) / num_steps_measure << endl;
      cout << endl;
      cout << "Simulation time:       " << run_time << endl;
      cout << "   Collision detection:       " << collision_time << endl;
      cout << "      Broad phase:                 " << collision_broad_time << endl;
      cout << "      Narrow phase:                " << collision_narrow_time << endl;
      cout << "   Update:                    " << update_time << endl;
      cout << "   Solver:                    " << solver_time << endl;
      cout << "      Contact force calculation:   " << contact_time << endl;
    }

    break;
  }

  return 0;
}
