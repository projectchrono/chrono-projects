#include <stdio.h>
#include <vector>
#include <cmath>

#include "chrono_parallel/ChSystemParallel.h"
#include "chrono_parallel/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_utils/opengl/ChOpenGL.h"
#endif

// Define this to save the data when using the OpenGL code
//#define SAVE_DATA

// Define macro DEM to use penalty-based contact. Otherwise, DVI
#define DEM


// Global variables (for callback)
int out_steps;

#ifdef DEM
  ChStreamOutAsciiFile sph_file("../TEST_BALL_DEM/sphere_pos.txt");
  const char* data_folder = "../TEST_BALL_DEM/POVRAY";
#else
  ChStreamOutAsciiFile sph_file("../TEST_BALL_DVI/sphere_pos.txt");
  const char* data_folder = "../TEST_BALL_DVI/POVRAY";
#endif


// -----------------------------------------------------------------------------
// Callback function invoked at every simulation frame. Generates output with a
// frequency based on the specified FPS rate.
// -----------------------------------------------------------------------------
template<class T>
void SimFrameCallback(T* mSys, const int frame)
{
  if (frame % out_steps != 0)
    return;

  chrono::Vector bodyAngs;

  sph_file << frame << "     ";
  std::cout << frame << "     ";

  for (int i = 0; i < mSys->Get_bodylist()->size(); ++i) {
    ChBody* abody = (ChBody*) mSys->Get_bodylist()->at(i);
#ifdef DEM
    assert(typeid(*abody) == typeid(ChBodyDEM));
#endif
    const ChVector<>& bodypos = abody->GetPos();
    bodyAngs = abody->GetRot().Q_to_NasaAngles();
    sph_file << bodypos.x  << "  " << bodypos.y  << "  " << bodypos.z  << "  ";
    sph_file << bodyAngs.x << "  " << bodyAngs.y << "  " << bodyAngs.z << "       ";
    std::cout << bodypos.x << "  " << bodypos.y << "  " << bodypos.z << "   |   ";
  }

  sph_file << "\n";
  std::cout << std::endl;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int threads = 8;

  // Simulation parameters
  double gravity = -9.81;
  double time_step = 1e-4;
  double time_end = 1;
  int    num_steps = std::ceil(time_end / time_step);

  int max_iteration = 20;

  // Output
  double out_fps = 1000;

  // Parameters for the falling ball
  int             ballId = 100;
  double          radius = 0.1;
  double          density = 2000;
  double          volume = (4.0/3) * CH_C_PI * radius * radius * radius;
  double          mass = density * volume;
  ChVector<>      inertia = 0.4 * mass * radius * radius * ChVector<>(1,1,1);
  ChVector<>      initPos(1, -1, 0.5);
  ChVector<>      initPos2(0.629447, -1.45809, 0.5);
  ChQuaternion<>  initRot(1,0,0,0);//(0.89, 0.4, -0.2, 0.088889);
  ChVector<>      initVel(0,0,0);//(0.75, 0.2, 0);

  // Parameters for the containing bin
  int    binId = -200;
  double hDimX = 5;          // length in x direction
  double hDimY = 2;          // depth in y direction
  double hDimZ = 0.5;        // height in z direction
  double hThickness = 0.1;   // wall thickness

  // Create system
#ifdef DEM
  ChSystemParallelDEM msystem;
#else
  ChSystemParallelDVI msystem;
#endif

  msystem.Set_G_acc(ChVector<>(0, 0, gravity));

  // Edit system settings
  msystem.SetParallelThreadNumber(threads);
  msystem.SetMaxiter(max_iteration);
  msystem.SetIterLCPmaxItersSpeed(max_iteration);
  msystem.SetTol(1e-3);
  msystem.SetTolSpeeds(1e-3);
  msystem.SetStep(time_step);

#ifdef DEM
  ((ChCollisionSystemParallel*) msystem.GetCollisionSystem())->setBinsPerAxis(I3(10, 10, 10));
  ((ChCollisionSystemParallel*) msystem.GetCollisionSystem())->setBodyPerBin(100, 50);

  ((ChCollisionSystemParallel*) msystem.GetCollisionSystem())->ChangeNarrowphase(new ChCNarrowphaseR);
#else
  ((ChLcpSolverParallelDVI*) msystem.GetLcpSolverSpeed())->SetMaxIteration(max_iteration);
  ((ChLcpSolverParallelDVI*) msystem.GetLcpSolverSpeed())->SetTolerance(0);
  ((ChLcpSolverParallelDVI*) msystem.GetLcpSolverSpeed())->SetCompliance(0);
  ((ChLcpSolverParallelDVI*) msystem.GetLcpSolverSpeed())->SetContactRecoverySpeed(1);
  ((ChLcpSolverParallelDVI*) msystem.GetLcpSolverSpeed())->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);

  ((ChCollisionSystemParallel*) msystem.GetCollisionSystem())->SetCollisionEnvelope(radius * 0.05);
  ((ChCollisionSystemParallel*) msystem.GetCollisionSystem())->setBinsPerAxis(I3(10, 10, 10));
  ((ChCollisionSystemParallel*) msystem.GetCollisionSystem())->setBodyPerBin(100, 50);
#endif

  omp_set_num_threads(threads);

  // Create a material for the balls
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> ballMat(new ChMaterialSurfaceDEM);
  ballMat->SetYoungModulus(2e5f);
  ballMat->SetFriction(0.4f);
  ballMat->SetDissipationFactor(0.6f);
#else
  ChSharedPtr<ChMaterialSurface> ballMat(new ChMaterialSurface);
  ballMat->SetFriction(0.4f);
#endif

  // Create a material for the bin
#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> binMat(new ChMaterialSurfaceDEM);
  binMat->SetYoungModulus(2e5f);
  binMat->SetFriction(0.4f);
  binMat->SetDissipationFactor(0.6f);
#else
  ChSharedPtr<ChMaterialSurface> binMat(new ChMaterialSurface);
  binMat->SetFriction(0.4f);
#endif

  // Create the falling ball
#ifdef DEM
  ChSharedBodyDEMPtr ball(new ChBodyDEM(new ChCollisionModelParallel));
  ball->SetMaterialSurfaceDEM(ballMat);
#else
  ChSharedBodyPtr ball(new ChBody(new ChCollisionModelParallel));
  ball->SetMaterialSurface(ballMat);
#endif

  ball->SetIdentifier(ballId);
  ball->SetMass(mass);
  ball->SetInertiaXX(inertia);
  ball->SetPos(initPos);
  ball->SetRot(initRot);
  ball->SetPos_dt(initVel);
  ball->SetBodyFixed(false);
  ball->SetCollide(true);

  ball->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball.get_ptr(), radius);
  ball->GetCollisionModel()->BuildModel();

  msystem.AddBody(ball);

  // Create a second ball
#ifdef DEM
  ChSharedBodyDEMPtr ball2(new ChBodyDEM(new ChCollisionModelParallel));
  ball2->SetMaterialSurfaceDEM(ballMat);
#else
  ChSharedBodyPtr ball2(new ChBody(new ChCollisionModelParallel));
  ball2->SetMaterialSurface(ballMat);
#endif

  ball2->SetIdentifier(ballId);
  ball2->SetMass(mass);
  ball2->SetInertiaXX(inertia);
  ball2->SetPos(initPos2);
  ball2->SetRot(initRot);
  ball2->SetPos_dt(initVel);
  ball2->SetBodyFixed(false);
  ball2->SetCollide(true);

  ball2->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball2.get_ptr(), radius);
  ball2->GetCollisionModel()->BuildModel();

  msystem.AddBody(ball2);


  // Create the containing bin
#ifdef DEM
  utils::CreateBoxContainerDEM(&msystem, binId, binMat, ChVector<>(hDimX, hDimY, hDimZ), hThickness);
#else
  utils::CreateBoxContainerDVI(&msystem, binId, binMat, ChVector<>(hDimX, hDimY, hDimZ), hThickness);
#endif

  // Perform the simulation
  // ----------------------

  out_steps = std::ceil((1 / time_step) / out_fps);
  int out_frame = 0;

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  // The OpenGL manager will automatically run the simulation
  utils::ChOpenGLManager * window_manager = new  utils::ChOpenGLManager();
  utils::ChOpenGL openGLView(window_manager, &msystem, 800, 600, 0, 0, "mixerDEM");
  openGLView.render_camera->camera_position = glm::vec3(0, 5, 10);
  openGLView.render_camera->camera_look_at = glm::vec3(0, 0, 0);
  openGLView.render_camera->camera_scale = 1;
#ifdef SAVE_DATA
  openGLView.SetCustomCallback(SimFrameCallback);
#endif
  openGLView.StartSpinning(window_manager);
  window_manager->CallGlutMainLoop();
#else
  // Run simulation for specified time
  double time = 0;

  for (int i = 0; i < num_steps; i++) {
    SimFrameCallback(&msystem, i);
    msystem.DoStepDynamics(time_step);
    time += time_step;
  }
#endif

  return 0;
}

