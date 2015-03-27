#include <stdio.h>
#include <vector>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"
#include "chrono_parallel/collision/ChCNarrowphaseRUtils.h"

#include "chrono_utils/ChUtilsVehicle.h"
#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::endl;

// =============================================================================

// -----------------------------------------------------------------------------
// Specification of the vehicle model
// -----------------------------------------------------------------------------

enum WheelType {CYLINDRICAL, LUGGED};

// Type of wheel/tire (controls both contact and visualization)
WheelType wheel_type = LUGGED;

// JSON files for vehicle model (using different wheel visualization meshes)
std::string vehicle_file_cyl("hmmwv/vehicle/HMMWV_Vehicle_simple.json");
std::string vehicle_file_lug("hmmwv/vehicle/HMMWV_Vehicle_simple_lugged.json");

// JSON files for powertrain (simple)
std::string simplepowertrain_file("hmmwv/powertrain/HMMWV_SimplePowertrain.json");

// Initial vehicle position and orientation
ChVector<> initLoc(0, 0, 1.0);
ChQuaternion<> initRot(1, 0, 0, 0);

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 100;

// Perform dynamic tuning of number of threads?
bool thread_tuning = true;

// Simulation duration.
double time_end = 5;

// Solver parameters
double time_step = 1e-3;

double tolerance = 1e-1;

int max_iteration_normal = 0;
int max_iteration_sliding = 200;
int max_iteration_spinning = 0;

float contact_recovery_speed = 0.1;

// Output
bool povray_output = false;

const std::string out_dir = "../HMMWV";
const std::string pov_dir = out_dir + "/POVRAY";

int out_fps = 60;

// Continuous loop (only if OpenGL available)
bool loop = false;

// =============================================================================

// Callback class for providing driver inputs.
class MyDriverInputs : public utils::DriverInputsCallback {
 public:
  virtual void onCallback(double time, double& throttle, double& steering, double& braking) {
    throttle = 0;
    steering = 0;
    braking = 0;

    if (time > 0.5)
      throttle = 1.0;
    else if (time > 0.25)
      throttle = 4 * (time - 0.25);
  }
};

// Callback class for specifying rigid tire contact model.
// This version uses cylindrical contact shapes.
class MyCylindricalTire : public utils::TireContactCallback {
 public:
  virtual void onCallback(ChSharedPtr<ChBody> wheelBody, double radius, double width) {
    wheelBody->GetCollisionModel()->ClearModel();
    wheelBody->GetCollisionModel()->AddCylinder(0.46, 0.46, width / 2);
    wheelBody->GetCollisionModel()->BuildModel();

    wheelBody->GetMaterialSurface()->SetFriction(0.6f);
  }
};

// Callback class for specifying rigid tire contact model.
// This version uses a collection of convex contact shapes (meshes).
class MyLuggedTire : public utils::TireContactCallback {
 public:
  MyLuggedTire() {
    std::string lugged_file("hmmwv/lugged_wheel_section.obj");
    geometry::ChTriangleMeshConnected lugged_mesh;
    utils::LoadConvexMesh(vehicle::GetDataFile(lugged_file), lugged_mesh, lugged_convex);
    num_hulls = lugged_convex.GetHullCount();
  }

  virtual void onCallback(ChSharedPtr<ChBody> wheelBody, double radius, double width) {
    ChCollisionModelParallel* coll_model = (ChCollisionModelParallel*) wheelBody->GetCollisionModel();
    coll_model->ClearModel();

    // Assemble the tire contact from 15 segments, properly offset.
    // Each segment is further decomposed in convex hulls.
    for (int iseg = 0; iseg < 15; iseg++) {
      ChQuaternion<> rot = Q_from_AngAxis(iseg * 24 * CH_C_DEG_TO_RAD, VECT_Y);
      for (int ihull = 0; ihull < num_hulls; ihull++) {
        std::vector<ChVector<> > convexhull;
        lugged_convex.GetConvexHullResult(ihull, convexhull);
        coll_model->AddConvexHull(convexhull, VNULL, rot);
      }
    }

    // Add a cylinder to represent the wheel hub.
    coll_model->AddCylinder(0.223, 0.223, 0.126);

    coll_model->BuildModel();

    wheelBody->GetMaterialSurface()->SetFriction(0.8f);
  }

 private:
  ChConvexDecompositionHACDv2 lugged_convex;
  int num_hulls;
};

// Callback class for specifying rigid tire contact model.
// This version uses a collection of convex contact shapes (meshes).
// In addition, this version overrides the visualization assets of the provided
// wheel body with the collision meshes.
class MyLuggedTire_vis : public utils::TireContactCallback {
public:
  MyLuggedTire_vis() {
    std::string lugged_file("hmmwv/lugged_wheel_section.obj");
    utils::LoadConvexMesh(vehicle::GetDataFile(lugged_file), lugged_mesh, lugged_convex);
  }

  virtual void onCallback(ChSharedPtr<ChBody> wheelBody, double radius, double width) {
    // Clear any existing assets (will be overriden)
    wheelBody->GetAssets().clear();

    wheelBody->GetCollisionModel()->ClearModel();
    for (int j = 0; j < 15; j++) {
      utils::AddConvexCollisionModel(wheelBody, lugged_mesh, lugged_convex, VNULL,
        Q_from_AngAxis(j * 24 * CH_C_DEG_TO_RAD, VECT_Y), false);
    }
    // This cylinder acts like the rims
    utils::AddCylinderGeometry(wheelBody.get_ptr(), 0.223, 0.126);
    wheelBody->GetCollisionModel()->BuildModel();

    wheelBody->GetMaterialSurface()->SetFriction(0.8f);
  }

private:
  ChConvexDecompositionHACDv2 lugged_convex;
  geometry::ChTriangleMeshConnected lugged_mesh;
};


// =============================================================================

int main(int argc, char* argv[]) {
  // Set path to ChronoVehicle data files
  vehicle::SetDataPath(CHRONO_VEHICLE_DATA_DIR);

  // --------------------------
  // Create output directories.
  // --------------------------

  if (povray_output) {
    if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
      cout << "Error creating directory " << out_dir << endl;
      return 1;
    }
    if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
      cout << "Error creating directory " << pov_dir << endl;
      return 1;
    }
  }

  // --------------
  // Create system.
  // --------------

  ChSystemParallelDVI* system = new ChSystemParallelDVI();

  system->Set_G_acc(ChVector<>(0, 0, -9.81));

  // ----------------------
  // Set number of threads.
  // ----------------------

  int max_threads = system->GetParallelThreadNumber();
  if (threads > max_threads)
    threads = max_threads;
  system->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;

  // ---------------------
  // Edit system settings.
  // ---------------------

  system->GetSettings()->solver.tolerance = tolerance;

  system->GetSettings()->solver.solver_mode = SLIDING;
  system->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
  system->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
  system->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
  system->GetSettings()->solver.alpha = 0;
  system->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  system->ChangeSolverType(APGD);

  system->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);

  // -----------------------------------------
  // Create and initialize the vehicle system.
  // -----------------------------------------

  utils::VehicleSystem* vehicle;
  utils::TireContactCallback* tire_cb;

  switch (wheel_type) {
    case CYLINDRICAL: {
      vehicle = new utils::VehicleSystem(system, vehicle_file_cyl, simplepowertrain_file);
      tire_cb = new MyCylindricalTire();
    } break;
    case LUGGED: {
      vehicle = new utils::VehicleSystem(system, vehicle_file_lug, simplepowertrain_file);
      tire_cb = new MyLuggedTire_vis();
    } break;
  }

  vehicle->SetTireContactCallback(tire_cb);

  MyDriverInputs driver_cb;
  vehicle->SetDriverInputsCallback(&driver_cb);

  vehicle->Initialize(initLoc, initRot);

  // ------------------------------------------------
  // Create the ground body and set contact geometry.
  // ------------------------------------------------
  double hdimX = 100;
  double hdimY = 100;
  double hdimZ = 5;

  ChSharedPtr<ChBody> ground = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
  ground->SetIdentifier(-1);
  ground->SetBodyFixed(true);
  ground->SetCollide(true);
  ground->SetPos(ChVector<>(0, 0, -hdimZ));

  ground->GetMaterialSurface()->SetFriction(0.8f);

  ground->GetCollisionModel()->ClearModel();
  ground->GetCollisionModel()->AddBox(hdimX, hdimY, hdimZ);
  ground->GetCollisionModel()->BuildModel();

  ChSharedPtr<ChBoxShape> box_ground(new ChBoxShape);
  box_ground->GetBoxGeometry().Size = ChVector<>(hdimX, hdimY, hdimZ);
  ground->AddAsset(box_ground);

  system->AddBody(ground);

// -----------------------
// Perform the simulation.
// -----------------------

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  // Initialize OpenGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "HMMWV", system);
  gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));

  // Let the OpenGL manager run the simulation until interrupted.
  // NOTE: we need to add a user callback to the OpenGL library to give back
  // control to the user before advancing the simulation by one step (in this
  // case to call the vehicle update function).
  if (loop) {
    gl_window.StartDrawLoop(time_step);
    return 0;
  }
#endif

  // Run simulation for specified time.
  int out_steps = std::ceil((1.0 / time_step) / out_fps);

  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  int next_out_frame = 0;
  double exec_time = 0;
  int num_contacts = 0;

  while (time < time_end) {
    if (povray_output && sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
      utils::WriteShapesPovray(system, filename);

      cout << "------------ Output frame:   " << out_frame + 1 << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
      cout << "             Execution time: " << exec_time << endl;

      out_frame++;
      next_out_frame += out_steps;
      num_contacts = 0;
    }

    // Update vehicle
    vehicle->Update(time);

// Advance dynamics.
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(time_step);
      gl_window.Render();
    } else
      break;
#else
    system->DoStepDynamics(time_step);
#endif

    // Update counters.
    time += time_step;
    sim_frame++;
    exec_time += system->GetTimerStep();
    num_contacts += system->GetNcontacts();
  }

  // Final stats
  cout << "==================================" << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  delete vehicle;
  delete tire_cb;

  return 0;
}
