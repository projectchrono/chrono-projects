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

// Chrono::Engine header files
#include "chrono/ChConfig.h"
#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsInputOutput.h"

// Chrono::Parallel header files
#include "chrono_parallel/collision/ChNarrowphaseRUtils.h"
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"

// Chrono::Parallel OpenGL header files
#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

// Chrono utility header files
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsInputOutput.h"

// Chrono vehicle header files
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/utils/ChSpeedController.h"
// M113 model header files

#include "fording_setup.h"
#include "input_output.h"
#include "m113_custom.h"
using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::endl;

// =============================================================================
// USER SETTINGS
// =============================================================================

// -----------------------------------------------------------------------------
// Specification of the terrain
// -----------------------------------------------------------------------------

enum TerrainType { RIGID_TERRAIN, GRANULAR_TERRAIN };

// Type of terrain
TerrainType terrain_type = RIGID_TERRAIN;

// Control visibility of containing bin walls
bool visible_walls = false;

// Dimensions
double hdimX = 5.5;  //// 2.5;
double hdimY = 2.5;
double hdimZ = 0.5;
double hthick = 0.25;

// Parameters for granular material
int Id_g = 100;
double r_g = 0.02;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

float mu_g = 0.8f;

unsigned int num_particles = 100;  //// 40000;
std::vector<real3> forces;
std::vector<real3> torques;
// -----------------------------------------------------------------------------
// Specification of the vehicle model
// -----------------------------------------------------------------------------
// Initial vehicle position and orientation

std::string data_output_path = "m113_fording/";
ChSpeedController m_speedPIDVehicle;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Total simulation duration.
double time_end = 20;

// Duration of the "hold time" (vehicle chassis fixed and no driver inputs).
// This can be used to allow the granular material to settle.
double time_hold = 0.0;

// Solver parameters
double time_step = 1e-3;  // 2e-4;

double tolerance = 0.00001;

int max_iteration_bilateral = 1000;  // 1000;
int max_iteration_normal = 0;
int max_iteration_sliding = 100;  // 2000;
int max_iteration_spinning = 0;

float contact_recovery_speed = 12;

// Periodically monitor maximum bilateral constraint violation
bool monitor_bilaterals = false;
int bilateral_frame_interval = 100;

int out_fps = 60;

double target_speed = 2;
// =============================================================================

void static WriteTrackedVehicleData(M113_Vehicle_Custom& vehicle,
                                    ChDriver::Inputs driver_inputs,
                                    std::vector<real3> forces,
                                    std::vector<real3> torques,
                                    std::string filename) {
    CSVGen csv_output;
    csv_output.OpenFile(filename.c_str(), false);

    auto m_driveline = vehicle.GetDriveline();
    auto m_powertrain = vehicle.GetPowertrain();

    csv_output << vehicle.GetChassisBody()->GetPos();
    csv_output << vehicle.GetVehicleSpeed();
    csv_output << m_driveline->GetDriveshaftSpeed();
    csv_output << m_powertrain->GetMotorTorque();
    csv_output << m_powertrain->GetMotorSpeed();
    csv_output << m_powertrain->GetOutputTorque();

    csv_output << driver_inputs.m_throttle;
    csv_output << driver_inputs.m_braking;

    double total_force_len = 0;
    double total_torque_len = 0;
    for (int i = 0; i < forces.size(); i++) {
        total_force_len += Length(forces[i]);
        total_torque_len += Length(torques[i]);
    }
    csv_output << Length(forces[0]);
    csv_output << Length(torques[0]);

    csv_output << total_force_len;
    csv_output << total_torque_len;

    csv_output << m_driveline->GetSprocketTorque(LEFT);
    csv_output << m_driveline->GetSprocketTorque(RIGHT);

    csv_output.endline();
    csv_output.CloseFile();
}

double CreateParticles(ChSystemParallelNSC* system) {
    // Create a material
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);

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

    while (gen.getTotalNumBodies() < num_particles) {
        gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims);
        center.z() += 2 * r;
    }

    std::cout << "Created " << gen.getTotalNumBodies() << " particles." << std::endl;

    return center.z();
}

// =============================================================================
// Utility function for displaying an ASCII progress bar for the quantity x
// which must be a value between 0 and n. The width 'w' represents the number
// of '=' characters corresponding to 100%.

void progressbar(unsigned int x, unsigned int n, unsigned int w = 50) {
    if ((x != n) && (x % (n / 100 + 1) != 0))
        return;

    float ratio = x / (float)n;
    unsigned int c = (unsigned int)(ratio * w);

    std::cout << std::setw(3) << (int)(ratio * 100) << "% [";
    for (unsigned int x = 0; x < c; x++)
        std::cout << "=";
    for (unsigned int x = c; x < w; x++)
        std::cout << " ";
    std::cout << "]\r" << std::flush;
}

// =============================================================================
int main(int argc, char* argv[]) {
    // --------------
    // Create system.
    // --------------
    // ----  Parallel
    std::cout << "Create Parallel DVI system" << std::endl;
    ChSystemParallelNSC* system = new ChSystemParallelNSC();

    system->Set_G_acc(ChVector<>(0, 0, -9.81));

    // ---------------------
    // Edit system settings.
    // ---------------------
    // Set solver parameters
    system->GetSettings()->solver.tolerance = tolerance;
    system->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    system->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    system->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    system->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    system->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;  // make 1000, should be about 220
    system->GetSettings()->solver.compute_N = false;
    system->GetSettings()->solver.alpha = 0;
    system->GetSettings()->solver.cache_step_length = true;
    system->GetSettings()->solver.use_full_inertia_tensor = false;
    system->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    system->GetSettings()->solver.bilateral_clamp_speed = 1e8;
    system->GetSettings()->min_threads = threads;
    system->ChangeSolverType(SolverType::BB);
    system->SetLoggingLevel(LoggingLevel::LOG_INFO);
    system->SetLoggingLevel(LoggingLevel::LOG_TRACE);

    system->GetSettings()->collision.collision_envelope = 0.1 * r_g;

    system->GetSettings()->collision.bins_per_axis = vec3(100, 20, 25);
    system->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
    system->GetSettings()->collision.fixed_bins = true;

    // -------------------
    // Create the terrain.
    // -------------------

    CreateContainer(system);

    CreateFluid(system);

    // --------------------------
    // Construct the M113 vehicle
    // --------------------------

    m_speedPIDVehicle.SetGains(4.0, 1.0, 0.00);

    // Create and initialize vehicle system
    M113_Vehicle_Custom vehicle(true, TrackShoeType::SINGLE_PIN, system);
    ////vehicle.SetStepsize(0.0001);

    vehicle.Initialize(ChCoordsys<>(initLoc, initRot));

    // vehicle.SetChassisVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetSprocketVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetRoadWheelVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetIdlerVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetRoadWheelAssemblyVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetTrackShoeVisualizationType(VisualizationType::PRIMITIVES);

    // Create the powertrain system
    auto powertrain = chrono_types::make_shared<M113_SimplePowertrain>("powertrain");
    vehicle.InitializePowertrain(powertrain);

    // ---------------
    // Simulation loop
    // ---------------

#ifdef CHRONO_OPENGL
    // Initialize OpenGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "M113", system);
    gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

    // Number of simulation steps between two 3D view render frames
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

    // Run simulation for specified time.
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;

    // Inter-module communication data
    BodyStates shoe_states_left(vehicle.GetNumTrackShoes(LEFT));
    BodyStates shoe_states_right(vehicle.GetNumTrackShoes(RIGHT));
    TerrainForces shoe_forces_left(vehicle.GetNumTrackShoes(LEFT));
    TerrainForces shoe_forces_right(vehicle.GetNumTrackShoes(RIGHT));

    bool set_time = true;

    while (time < time_end) {
        // Driver inputs
        ChDriver::Inputs driver_inputs = {0, 0, 0};
        {
            double out_speed = m_speedPIDVehicle.Advance(vehicle, target_speed, time_step);

            if (time > .3) {
                ChClampValue(out_speed, -2.0, 2.0);
            } else {
                out_speed = 0;
            }

            if (out_speed > 0) {
                driver_inputs.m_braking = 0;
                driver_inputs.m_throttle = out_speed;
            } else {
                driver_inputs.m_braking = -out_speed;
                driver_inputs.m_throttle = 0;
            }

            if (vehicle.GetChassis()->GetPos().x() > dist_end) {
                driver_inputs.m_braking = 1;
                driver_inputs.m_throttle = 0;
                if (set_time) {
                    real time_end_temp = time + 1;
                    time_end = Max(time_end_temp, time_end);
                    set_time = false;
                }
            }
        }

        double powertrain_torque = powertrain->GetOutputTorque();
        double driveshaft_speed = vehicle.GetDriveshaftSpeed();
        vehicle.GetTrackShoeStates(LEFT, shoe_states_left);
        vehicle.GetTrackShoeStates(RIGHT, shoe_states_right);

        // Output
        if (sim_frame == next_out_frame) {
            std::cout << "write: " << out_frame << std::endl;
            forces.resize(0);
            torques.resize(0);

            system->CalculateContactForces();

            // get total force and torques on tracks and chassis
            forces.push_back(system->GetBodyContactForce(vehicle.GetChassisBody()->GetId()));
            torques.push_back(system->GetBodyContactTorque(vehicle.GetChassisBody()->GetId()));

            for (int i = 0; i < vehicle.GetNumTrackShoes(LEFT); i++) {
                forces.push_back(system->GetBodyContactForce(vehicle.GetTrackShoe(LEFT, i)->GetShoeBody()->GetId()));
                torques.push_back(system->GetBodyContactTorque(vehicle.GetTrackShoe(LEFT, i)->GetShoeBody()->GetId()));
            }
            for (int i = 0; i < vehicle.GetNumTrackShoes(RIGHT); i++) {
                forces.push_back(system->GetBodyContactForce(vehicle.GetTrackShoe(RIGHT, i)->GetShoeBody()->GetId()));
                torques.push_back(system->GetBodyContactTorque(vehicle.GetTrackShoe(RIGHT, i)->GetShoeBody()->GetId()));
            }

            DumpFluidData(system, data_output_path + "data_" + std::to_string(out_frame) + ".dat", true);
            DumpAllObjectsWithGeometryPovray(system, data_output_path + "vehicle_" + std::to_string(out_frame) + ".dat",
                                             true);
            WriteTrackedVehicleData(vehicle, driver_inputs, forces, torques,
                                    data_output_path + "stats_" + std::to_string(out_frame) + ".dat");

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

        // Release the vehicle chassis at the end of the hold time.
        if (vehicle.GetChassisBody()->GetBodyFixed() && time > time_hold) {
            vehicle.GetChassisBody()->SetBodyFixed(false);
        }

        vehicle.Synchronize(time, driver_inputs, shoe_forces_left, shoe_forces_right);
        vehicle.Advance(time_step);

#ifdef CHRONO_OPENGL
        // gl_window.Pause();

        if (gl_window.Active()) {
            if (gl_window.DoStepDynamics(time_step)) {
                // Update counters.
                time += time_step;
                sim_frame++;
                exec_time += system->GetTimerStep();
                num_contacts += system->GetNcontacts();
            }
            gl_window.Render();
        } else {
            break;
        }
#else
        system->DoStepDynamics(time_step);
        time += time_step;
        sim_frame++;
        exec_time += system->GetTimerStep();
        num_contacts += system->GetNcontacts();

#endif
    }

    // Final stats
    std::cout << "==================================" << std::endl;
    std::cout << "Simulation time:   " << exec_time << std::endl;

    return 0;
}
