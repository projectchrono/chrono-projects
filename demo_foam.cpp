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

//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "demo_utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::endl;

// =======================================================================
// Global problem definitions

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 8;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Simulation parameters
double gravity = 9.81;
double time_step = 1e-4;
double time_end = 5;

int max_iteration = 20;

// Output
const std::string out_dir = "../FOAM";
const std::string pov_dir = out_dir + "/POVRAY";
const std::string out_file = out_dir + "/timing.dat";
double out_fps = 50;

// Parameters for the granular material
int Id_g = 1;
double r_g = 0.01;
double rho_g = 2000;

float Y_g = 2e7;
float mu_g = 0.2;
float cr_g = 0.1;
float cohesion_g = 300;

// Parameters for the containing bin
int binId = -200;
double hDimX = 10;        // length in x direction
double hDimY = 10;        // depth in y direction
double hDimZ = 1;         // height in z direction
double hThickness = 0.4;  // wall thickness

float Y_c = 2e6;
float mu_c = 0.3;
float cr_c = 0.1;
float cohesion_c = 5;

// Particle generator
utils::Generator* gen;

double initVel = 5;  // initial particle velocity in negative X direction

int maxNumParticles = 100000;

// =======================================================================

int SpawnParticles() {
  double dist = 2 * 0.99 * r_g;

  ////gen->createObjectsBox(utils::POISSON_DISK,
  ////                     dist,
  ////                     ChVector<>(9, 0, 3),
  ////                     ChVector<>(0, 1, 0.5),
  ////                     ChVector<>(-initVel, 0, 0));
  gen->createObjectsCylinderX(utils::POISSON_DISK, dist, ChVector<>(9, 0, 3), 0.2, 0, ChVector<>(-initVel, 0, 0));
  cout << "  total bodies: " << gen->getTotalNumBodies() << endl;

  return gen->getTotalNumBodies();
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

  msystem->GetSettings()->perform_thread_tuning = thread_tuning;

  // Set gravitational acceleration
  msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

  // Edit system settings
  msystem->GetSettings()->solver.tolerance = 1e-4;

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);

  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;

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

  // Create the containing bin
  utils::CreateBoxContainer(msystem, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);

  // Create a mixture entirely made out of spheres
  double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
  double mass_g = rho_g * vol_g;
  ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

  gen = new utils::Generator(msystem);

  utils::MixtureIngredientPtr& m1 = gen->AddMixtureIngredient(utils::SPHERE, 1.0);
  m1->setDefaultMaterialDEM(mat_g);
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  gen->setBodyIdentifier(Id_g);

  // Number of steps
  int num_steps = std::ceil(time_end / time_step);
  int out_steps = std::ceil((1 / time_step) / out_fps);
  int gen_steps = std::ceil(3 * r_g / initVel / time_step);

  // Create output directories.
  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  // Initialize OpenGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "Foam demo", msystem);
  gl_window.SetCamera(ChVector<>(4, -5, 4), ChVector<>(9, 0, 3), ChVector<>(0, 0, 1));
  gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

  // Perform the simulation
  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  double exec_time = 0;
  ChStreamOutAsciiFile ofile(out_file.c_str());

  while (time < time_end) {
    int numParticles = msystem->Get_bodylist()->size() - 1;

    if (numParticles < maxNumParticles && sim_frame % gen_steps == 0) {
      SpawnParticles();
    }

    if (sim_frame % out_steps == 0) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
      utils::WriteShapesPovray(msystem, filename);

      cout << " --------------------------------- Output frame:   " << out_frame + 1 << endl;
      cout << "                                   Sim frame:      " << sim_frame << endl;
      cout << "                                   Time:           " << time << endl;
      cout << "                                   Execution time: " << exec_time << endl;
      cout << "                                   Num. bodies:    " << numParticles << endl;

      ofile << sim_frame << "  " << time << "  " << exec_time << "  " << numParticles << "  " << msystem->GetNcontacts()
            << "\n";

      out_frame++;
    }

    // Advance dynamics.
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(time_step);
      gl_window.Render();
    } else
      break;
#else
    msystem->DoStepDynamics(time_step);
#endif

    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
  }

  // Final stats
  cout << "==================================" << endl;
  cout << "Number of bodies: " << msystem->Get_bodylist()->size() << endl;
  cout << "Simulation time: " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}
