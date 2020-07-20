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
//
// =============================================================================

#include <iostream>
#include <memory>

// Chrono::Engine header files
#include "chrono/ChConfig.h"
#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsInputOutput.h"

// Chrono::Parallel header files
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"
#include "chrono_parallel/collision/ChNarrowphaseRUtils.h"

// Chrono::Parallel OpenGL header files
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

// Chrono utility header files
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

// Chrono vehicle header files
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"

// M113 model header files
#include "chrono_models/vehicle/m113/M113_Vehicle.h"
#include "chrono_models/vehicle/m113/M113_SimplePowertrain.h"

#include "chrono_thirdparty/filesystem/path.h"

// Utilities
#include "../../utils.h"

using namespace chrono;
using namespace chrono::collision;
using namespace chrono::vehicle;
using namespace chrono::vehicle::m113;

// =============================================================================
// USER SETTINGS
// =============================================================================

// -----------------------------------------------------------------------------
// Specification of the terrain
// -----------------------------------------------------------------------------

enum TerrainType { RIGID_TERRAIN, GRANULAR_TERRAIN };

// Type of terrain
TerrainType terrain_type = GRANULAR_TERRAIN;

// Control visibility of containing bin walls
bool visible_walls = true;

// Dimensions
double hdimX = 4.5;
double hdimY = 1.75;
double hdimZ = 0.5;
double hthick = 0.25;
double numLayers = 8;

// Parameters for granular material
int Id_g = 100;
double r_g = 18e-3;
double rho_g = 2500;
double coh_pressure = 3e4;
float mu_g = 0.9f;

double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);
double coh_force = CH_C_PI * r_g * r_g * coh_pressure;

// -----------------------------------------------------------------------------
// Specification of the vehicle model
// -----------------------------------------------------------------------------

// Initial vehicle position and orientation
ChVector<> initLoc(-hdimX + 4.5, 0, 1.0);
ChQuaternion<> initRot(1, 0, 0, 0);

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Input file names for the path-follower driver model
std::string steering_controller_file("generic/driver/SteeringController.json");
std::string speed_controller_file("generic/driver/SpeedController.json");
std::string path_file("paths/straight10km.txt");

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Total simulation duration.
double time_end = 7;

// Duration of the "hold time" (vehicle chassis fixed and no driver inputs).
// This can be used to allow the granular material to settle.
double time_hold = 0.2;

// Solver parameters
double time_step = 1e-3;
double tolerance = 1e-5;

int max_iteration_bilateral = 1000;
int max_iteration_normal = 0;
int max_iteration_sliding = 100;
int max_iteration_spinning = 0;

float contact_recovery_speed = 12;

// Periodically monitor maximum bilateral constraint violation
bool monitor_bilaterals = false;
int bilateral_frame_interval = 100;

// Output directories
bool povray_output = true;

const std::string out_dir = "../M113_PARALLEL_NSC";
const std::string pov_dir = out_dir + "/POVRAY";

int out_fps = 60;

// =============================================================================

class MyDriver : public ChDriver {
  public:
    MyDriver(ChVehicle& vehicle, double delay) : ChDriver(vehicle), m_delay(delay) {}
    ~MyDriver() {}

    virtual void Synchronize(double time) override {
        m_throttle = 0;
        m_steering = 0;
        m_braking = 0;

        double eff_time = time - m_delay;

        // Do not generate any driver inputs for a duration equal to m_delay.
        if (eff_time < 0)
            return;

        if (eff_time > 4.2) {
            m_throttle = 0;
            m_braking = 1;
        } else if (eff_time > 0.2) {
            m_throttle = 0.8;
        } else {
            m_throttle = 4 * eff_time;
        }
    }

  private:
    double m_delay;
};

// =============================================================================

double CreateParticles(ChSystem* system) {
    // Create a material
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(0.0f);
    mat_g->SetCohesion(static_cast<float>(coh_force));

    // Create a particle generator and a mixture entirely made out of spheres
    utils::Generator gen(system);
    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    // Set starting value for body identifiers
    gen.setBodyIdentifier(Id_g);

    // Create particles in layers until reaching the desired number of particles
    double r = 1.01 * r_g;
    ChVector<> hdims(hdimX - r, hdimY - r, 0);
    ChVector<> center(0, 0, 2 * r);

    double layerCount = 0;
    while (layerCount < numLayers) {
        gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims);
        center.z() += 2 * r;
        layerCount++;
    }

    std::cout << "Created " << gen.getTotalNumBodies() << " particles." << std::endl;

    return center.z();
}

// =============================================================================
int main(int argc, char* argv[]) {
    // -------------------------------------------------------
    // Set path to Chrono and Chrono::Vehicle data directories
    // -------------------------------------------------------

    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_VEHICLE_DATA_DIR);

    // ---------------------------------
    // Create system and modify settings
    // ---------------------------------

    std::cout << "Create Parallel NSC system" << std::endl;
    ChSystemParallelNSC system;

    system.Set_G_acc(ChVector<>(0, 0, -9.80665));

    // Set number of threads
    int max_threads = CHOMPfunctions::GetNumProcs();
    if (threads > max_threads)
        threads = max_threads;
    CHOMPfunctions::SetNumThreads(threads);
    std::cout << "Using " << threads << " threads" << std::endl;

    system.GetSettings()->perform_thread_tuning = thread_tuning;

    // Set solver parameters
    system.GetSettings()->solver.use_full_inertia_tensor = false;
    system.GetSettings()->solver.tolerance = tolerance;
    system.GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    system.GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    system.GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    system.GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    system.GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    system.GetSettings()->solver.alpha = 0;
    system.GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    system.GetSettings()->solver.bilateral_clamp_speed = 1e8;
    system.ChangeSolverType(SolverType::BB);

    system.GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
    system.GetSettings()->collision.collision_envelope = 0.001;

    system.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    int factor = 2;
    int binsX = (int)std::ceil(hdimX / r_g) / factor;
    int binsY = (int)std::ceil(hdimY / r_g) / factor;
    int binsZ = 1;
    system.GetSettings()->collision.bins_per_axis = vec3(binsX, binsY, binsZ);
    std::cout << "broad-phase bins: " << binsX << " x " << binsY << " x " << binsZ << std::endl;

    // ------------------
    // Create the terrain
    // ------------------

    // Contact material
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);

    // Ground body
    auto ground = std::shared_ptr<ChBody>(system.NewBody());
    ground->SetIdentifier(-1);
    ground->SetBodyFixed(true);
    ground->SetCollide(true);

    ground->GetCollisionModel()->ClearModel();

    // Bottom box
    utils::AddBoxGeometry(ground.get(), mat_g, ChVector<>(hdimX, hdimY, hthick), ChVector<>(0, 0, -hthick),
                          ChQuaternion<>(1, 0, 0, 0), true);
    if (terrain_type == GRANULAR_TERRAIN) {
        // Front box
        utils::AddBoxGeometry(ground.get(), mat_g, ChVector<>(hthick, hdimY, hdimZ + hthick),
                              ChVector<>(hdimX + hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);
        // Rear box
        utils::AddBoxGeometry(ground.get(), mat_g, ChVector<>(hthick, hdimY, hdimZ + hthick),
                              ChVector<>(-hdimX - hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0),
                              visible_walls);
        // Left box
        utils::AddBoxGeometry(ground.get(), mat_g, ChVector<>(hdimX, hthick, hdimZ + hthick),
                              ChVector<>(0, hdimY + hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);
        // Right box
        utils::AddBoxGeometry(ground.get(), mat_g, ChVector<>(hdimX, hthick, hdimZ + hthick),
                              ChVector<>(0, -hdimY - hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0),
                              visible_walls);
    }

    ground->GetCollisionModel()->BuildModel();

    system.AddBody(ground);

    // Create the granular material.
    double vertical_offset = 0;

    if (terrain_type == GRANULAR_TERRAIN) {
        vertical_offset = CreateParticles(&system);
    }

    // --------------------------
    // Construct the M113 vehicle
    // --------------------------

    // Create and initialize vehicle systems
    auto vehicle = chrono_types::make_shared<M113_Vehicle>(true, TrackShoeType::SINGLE_PIN, BrakeType::SIMPLE, &system);
    auto powertrain = chrono_types::make_shared<M113_SimplePowertrain>("Powertrain");

    vehicle->Initialize(ChCoordsys<>(initLoc + ChVector<>(0.0, 0.0, vertical_offset), initRot));

    // Set visualization type for subsystems
    vehicle->SetChassisVisualizationType(VisualizationType::PRIMITIVES);
    vehicle->SetSprocketVisualizationType(VisualizationType::MESH);
    vehicle->SetIdlerVisualizationType(VisualizationType::MESH);
    vehicle->SetRoadWheelAssemblyVisualizationType(VisualizationType::PRIMITIVES);
    vehicle->SetRoadWheelVisualizationType(VisualizationType::MESH);
    vehicle->SetTrackShoeVisualizationType(VisualizationType::MESH);

    ////vehicle->SetCollide(TrackCollide::NONE);
    ////vehicle->SetCollide(TrackCollide::WHEELS_LEFT | TrackCollide::WHEELS_RIGHT);
    ////vehicle->SetCollide(TrackCollide::ALL & (~TrackCollide::SPROCKET_LEFT) & (~TrackCollide::SPROCKET_RIGHT));

    // Initialize the powertrain system
    vehicle->InitializePowertrain(powertrain);

    // Create the driver system (temporarily 1 for steering 1 for steering & brakes)
    MyDriver driver_speed(*vehicle, 0.5);
    driver_speed.Initialize();

	auto path = ChBezierCurve::read(vehicle::GetDataFile(path_file));
    ChPathFollowerDriver driver_steering(*vehicle, vehicle::GetDataFile(steering_controller_file),
                                vehicle::GetDataFile(speed_controller_file), path, "my_path", 0.0);
    driver_steering.Initialize();
	
    // ------------------------------------
    // Prepare output directories and files
    // ------------------------------------

    // Create output directories
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    if (povray_output) {
        if (!filesystem::create_directory(filesystem::path(pov_dir))) {
            std::cout << "Error creating directory " << pov_dir << std::endl;
            return 1;
        }
    }

    utils::CSV_writer csv("\t");
    csv.stream().setf(std::ios::scientific | std::ios::showpos);
    csv.stream().precision(6);

    // ---------------
    // Simulation loop
    // ---------------

    // Number of simulation steps between two 3D view render frames
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

#ifdef CHRONO_OPENGL
    // Initialize OpenGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "M113", &system);
    gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

    // Run simulation for specified time.
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;

    // Inter-module communication data
    BodyStates shoe_states_left(vehicle->GetNumTrackShoes(LEFT));
    BodyStates shoe_states_right(vehicle->GetNumTrackShoes(RIGHT));
    TerrainForces shoe_forces_left(vehicle->GetNumTrackShoes(LEFT));
    TerrainForces shoe_forces_right(vehicle->GetNumTrackShoes(RIGHT));

    while (time < time_end) {
        // Collect output data from modules
        ChDriver::Inputs driver_inputs;
        driver_inputs.m_throttle = driver_speed.GetThrottle();
        driver_inputs.m_braking = driver_speed.GetBraking();
        driver_inputs.m_steering = driver_steering.GetSteering();
        vehicle->GetTrackShoeStates(LEFT, shoe_states_left);
        vehicle->GetTrackShoeStates(RIGHT, shoe_states_right);

        const ChVector<>& pos_CG = vehicle->GetChassis()->GetPos();
        ChVector<> vel_CG = vehicle->GetChassisBody()->GetPos_dt();
        vel_CG = vehicle->GetChassisBody()->GetCoord().TransformDirectionParentToLocal(vel_CG);

        // Vehicle and Control Values
        csv << time << driver_inputs.m_steering << driver_inputs.m_throttle << driver_inputs.m_braking;
        csv << vehicle->GetTrackAssembly(LEFT)->GetSprocket()->GetAxleSpeed()
            << vehicle->GetTrackAssembly(RIGHT)->GetSprocket()->GetAxleSpeed();
        csv << powertrain->GetMotorSpeed() << powertrain->GetMotorTorque();
        csv << powertrain->GetOutputTorque() << vehicle->GetDriveshaftSpeed();
        // Chassis Position & Velocity
        csv << pos_CG.x() << pos_CG.y() << pos_CG.z();
        csv << vel_CG.x() << vel_CG.y() << vel_CG.z();
        csv << std::endl;

        // Output
        if (sim_frame == next_out_frame) {
            std::cout << std::endl;
            std::cout << "---- Frame:          " << out_frame + 1 << std::endl;
            std::cout << "     Sim frame:      " << sim_frame << std::endl;
            std::cout << "     Time:           " << time << std::endl;
            std::cout << "     Avg. contacts:  " << num_contacts / out_steps << std::endl;
            std::cout << "     Throttle input: " << driver_inputs.m_throttle << std::endl;
            std::cout << "     Braking input:  " << driver_inputs.m_braking << std::endl;
            std::cout << "     Steering input: " << driver_inputs.m_steering << std::endl;
            std::cout << "     Execution time: " << exec_time << std::endl;

            if (povray_output) {
                char filename[100];
                sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
                utils::WriteShapesPovray(&system, filename);
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;

            csv.write_to_file(out_dir + "/output.dat");
        }

        // Release the vehicle chassis at the end of the hold time.
        if (vehicle->GetChassis()->IsFixed() && time > time_hold) {
            std::cout << std::endl << "Release vehicle t = " << time << std::endl;
            vehicle->GetChassisBody()->SetBodyFixed(false);
        }

        // Update modules (process inputs from other modules)
        driver_speed.Synchronize(time);
		driver_steering.Synchronize(time);
        vehicle->Synchronize(time, driver_inputs, shoe_forces_left, shoe_forces_right);

        // Advance simulation for one timestep for all modules
        driver_speed.Advance(time_step);
		driver_steering.Advance(time_step);
        vehicle->Advance(time_step);

#ifdef CHRONO_OPENGL
        if (gl_window.Active())
            gl_window.Render();
        else
            break;
#endif

        progressbar(out_steps + sim_frame - next_out_frame + 1, out_steps);

        // Periodically display maximum constraint violation
        if (monitor_bilaterals && sim_frame % bilateral_frame_interval == 0) {
            vehicle->LogConstraintViolations();
        }

        // Update counters.
        time += time_step;
        sim_frame++;
        exec_time += system.GetTimerStep();
        num_contacts += system.GetNcontacts();
    }

    // Final stats
    std::cout << "==================================" << std::endl;
    std::cout << "Simulation time:   " << exec_time << std::endl;
    std::cout << "Number of threads: " << threads << std::endl;

    csv.write_to_file(out_dir + "/output.dat");

    return 0;
}
