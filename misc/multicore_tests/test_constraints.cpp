#include "physics/ChSystem.h"

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsInputOutput.h"

using namespace chrono;
using namespace chrono::collision;

// Define this to save the data when using the OpenGL code
//#define SAVE_DATA

#define SLIDER_CRANK

int main(int argc, char* argv[]) {
  int threads = 1;

  // Simulation parameters
  double gravity = 9.81;
  double time_step = 0.01;
  double time_end = 4.0;
  double out_fps = 50;

  uint max_iteration = 50;
  real tolerance = real(1e-8);

  const std::string out_dir = "../TEST_CONSTRAINTS";
  const std::string pov_dir = out_dir + "/POVRAY";

  // Create system.
  // --------------

  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();

  // Set gravitational acceleration
  msystem->Set_G_acc(ChVector<>(0, -gravity, 0));

  // Set number of threads.
  int max_threads = omp_get_num_procs();
  if (threads > max_threads)
    threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);

  // Specify narrowphase collision detection
  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;

  // Set tolerance and maximum number of iterations for bilateral constraint solver
  msystem->GetSettings()->solver.max_iteration_bilateral = max_iteration;
  msystem->GetSettings()->solver.tolerance = tolerance;

  // Create the rigid bodies of the slider-crank mechanical system.
  // --------------------------------------------------------------

  // Rotation of -90 degrees around z
  ChVector<> zero(0, 0, 0);
  ChQuaternion<> y2x(1, 0, 0, 0);
  y2x.Q_from_AngAxis(-CH_PI / 2, ChVector<>(0, 0, 1));

  // ground
  ChSharedPtr<ChBody> ground(new ChBody(new ChCollisionModelParallel, ChBody::DEM));
  msystem->AddBody(ground);
  ground->SetIdentifier(0);
  ground->SetBodyFixed(true);
  ground->SetCollide(false);

  // crank
  ChSharedPtr<ChBody> crank(new ChBody(new ChCollisionModelParallel, ChBody::DEM));
  msystem->AddBody(crank);
  crank->SetIdentifier(1);
  crank->SetPos(ChVector<>(1, 0, 0));
  crank->SetCollide(false);
  crank->GetCollisionModel()->ClearModel();
  utils::AddCapsuleGeometry(crank.get_ptr(), 0.1, 1, zero, y2x);
  crank->GetCollisionModel()->BuildModel();

#ifdef SLIDER_CRANK
  // rod
  ChSharedPtr<ChBody> rod(new ChBody(new ChCollisionModelParallel, ChBody::DEM));
  msystem->AddBody(rod);
  rod->SetIdentifier(2);
  rod->SetPos(ChVector<>(4, 0, 0));
  rod->SetCollide(false);
  rod->GetCollisionModel()->ClearModel();
  utils::AddCapsuleGeometry(rod.get_ptr(), 0.1, 2, zero, y2x);
  rod->GetCollisionModel()->BuildModel();
#endif

  // HACK!!!!
  ChSharedBodyPtr ground_(ground);
  ChSharedBodyPtr crank_(crank);
#ifdef SLIDER_CRANK
  ChSharedBodyPtr rod_(rod);
#endif

// Create joint constraints. Joint locations are specified in global frame.
// ------------------------------------------------------------------------

#ifdef SLIDER_CRANK
  // Revolute joint between crank and rod
  ChSharedPtr<ChLinkLockRevolute> rev_crank_rod(new ChLinkLockRevolute);
  rev_crank_rod->Initialize(crank_, rod_, ChCoordsys<>(ChVector<>(2, 0, 0)));
  msystem->AddLink(rev_crank_rod);

  // Slider (point on line) joint between rod and ground
  ChSharedPtr<ChLinkLockPointLine> slider_rod_ground(new ChLinkLockPointLine);
  slider_rod_ground->Initialize(rod_, ground_, ChCoordsys<>(ChVector<>(6, 0, 0)));
  msystem->AddLink(slider_rod_ground);
#endif

  // Engine between ground and crank (also acts as a revolute joint)
  ////ChSharedPtr<ChLinkEngine> engine_ground_crank(new ChLinkEngine);
  ////engine_ground_crank->Initialize(ground_, crank_, ChCoordsys<>(ChVector<>(0,0,0)));
  ////engine_ground_crank->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
  ////if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(engine_ground_crank->Get_spe_funct()))
  ////	mfun->Set_yconst(CH_PI); // speed w=3.145 rad/sec
  ////msystem->AddLink(engine_ground_crank);

  // Revolute between ground and crank
  ChSharedPtr<ChLinkLockRevolute> rev_ground_crank(new ChLinkLockRevolute);
  rev_ground_crank->Initialize(ground_, crank_, ChCoordsys<>(ChVector<>(0, 0, 0)));
  msystem->AddLink(rev_ground_crank);

  // Create output directories
  // -------------------------

  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    std::cout << "Error creating directory " << out_dir << std::endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    std::cout << "Error creating directory " << pov_dir << std::endl;
    return 1;
  }

  // Perform the simulation.
  // -----------------------

  int num_steps = (int)std::ceil(time_end / time_step);
  int out_steps = (int)std::ceil((1 / time_step) / out_fps);

  double time = 0;
  int out_frame = 0;
  char filename[100];

  for (int i = 0; i < num_steps; i++) {
    if (i % out_steps == 0) {
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame);
        utils::WriteVisualizationAssets(msystem, filename);

      out_frame++;
    }

    msystem->DoStepDynamics(time_step);
    time += time_step;
  }

  return 0;
}
