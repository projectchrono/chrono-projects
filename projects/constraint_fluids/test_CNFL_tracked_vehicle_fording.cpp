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

#include <algorithm>
#include <iostream>

// Chrono::Engine header files
#include "chrono/ChConfig.h"
#include "chrono/core/ChDataPath.h"

// Chrono::Multicore header files
#include "chrono_multicore/physics/ChSystemMulticore.h"

// Control use of VSG run-time rendering
// Note: CHRONO_VSG is defined in ChConfig.h
#undef CHRONO_VSG

#ifdef CHRONO_VSG
    #include "chrono_vsg/ChVisualSystemVSG.h"
#endif

// Chrono utility header files
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"

// Chrono vehicle header files
#include "chrono_vehicle/ChVehicleDataPath.h"
#include "chrono_vehicle/utils/ChSpeedController.h"

// M113 model header files
#include "chrono_models/vehicle/m113/M113.h"

#include "fording_setup.h"
#include "input_output.h"

using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::m113;

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
double vol_g = (4.0 / 3) * CH_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector3d inertia_g = 0.4 * mass_g * r_g * r_g * ChVector3d(1, 1, 1);

float mu_g = 0.8f;

unsigned int num_particles = 100;  //// 40000;
std::vector<real3> forces;
std::vector<real3> torques;
// -----------------------------------------------------------------------------
// Specification of the vehicle model
// -----------------------------------------------------------------------------
// Initial vehicle position and orientation

std::string data_output_path = "m113_fording/";
ChSpeedControllerPID m_speedPIDVehicle;

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

void static WriteTrackedVehicleData(M113& m113, DriverInputs driver_inputs, std::vector<real3> forces, std::vector<real3> torques, std::string filename) {
    CSVGen csv_output;
    csv_output.OpenFile(filename.c_str(), false);

    auto& vehicle = m113.GetVehicle();
    auto m_driveline = vehicle.GetDriveline();
    auto m_engine = vehicle.GetEngine();
    auto m_transmission = vehicle.GetTransmission();

    csv_output << vehicle.GetChassisBody()->GetPos();
    csv_output << vehicle.GetSpeed();
    csv_output << m_driveline->GetOutputDriveshaftSpeed();
    csv_output << m_engine->GetOutputMotorshaftTorque();
    csv_output << m_engine->GetMotorSpeed();
    csv_output << m_transmission->GetOutputDriveshaftTorque();

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

double CreateParticles(ChSystemMulticoreNSC* system) {
    // Create a material
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(mu_g);

    // Create a particle generator and a mixture entirely made out of spheres
    double r = 1.01 * r_g;
    utils::ChPDSampler<double> sampler(2 * r);
    utils::ChGenerator gen(system);
    std::shared_ptr<utils::ChMixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->SetDefaultMaterial(mat_g);
    m1->SetDefaultDensity(rho_g);
    m1->SetDefaultSize(ChVector3d(r_g));

    // Set starting value for body tags
    gen.SetStartTag(Id_g);

    // Create particles in layers until reaching the desired number of particles
    ChVector3d hdims(hdimX - r, hdimY - r, 0);
    ChVector3d center(0, 0, 2 * r);

    while (gen.GetTotalNumBodies() < num_particles) {
        gen.CreateObjectsBox(sampler, center, hdims);
        center.z() += 2 * r;
    }

    std::cout << "Created " << gen.GetTotalNumBodies() << " particles." << std::endl;

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
    SetChronoDataPath(CHRONO_DATA_DIR);
    SetVehicleDataPath(CHRONO_VEHICLE_DATA_DIR);

    if (!CreateOutputDirectory(data_output_path)) {
        std::cout << "Error creating directory " << data_output_path << std::endl;
        return 1;
    }

    // --------------
    // Create system.
    // --------------
    std::cout << "Create multicore NSC system" << std::endl;
    ChSystemMulticoreNSC* system = new ChSystemMulticoreNSC();
    system->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);
    system->SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));
    system->SetNumThreads(threads);

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
    system->ChangeSolverType(SolverType::BB);

    system->GetSettings()->collision.collision_envelope = 0.1 * r_g;

    system->GetSettings()->collision.bins_per_axis = vec3(100, 20, 25);
    system->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    // -------------------
    // Create the terrain.
    // -------------------

    CreateContainer(system);

    CreateFluid(system);

    // --------------------------
    // Construct the M113 vehicle
    // --------------------------

    m_speedPIDVehicle.SetGains(4.0, 1.0, 0.00);

    // Create and initialize vehicle system.
    // The chassis collides through its convex hulls, but not with the track assembly components.
    M113 m113(system);
    m113.SetChassisFixed(true);
    m113.SetChassisCollisionType(CollisionType::HULLS);
    m113.SetTrackShoeType(TrackShoeType::SINGLE_PIN);
    m113.SetBrakeType(BrakeType::SIMPLE);
    m113.SetDrivelineType(DrivelineTypeTV::SIMPLE);
    m113.SetEngineType(EngineModelType::SIMPLE);
    m113.SetInitPosition(ChCoordsys<>(initLoc, initRot));
    m113.Initialize();

    m113.SetChassisVisualizationType(VisualizationType::NONE);
    m113.SetSprocketVisualizationType(VisualizationType::PRIMITIVES);
    m113.SetRoadWheelVisualizationType(VisualizationType::PRIMITIVES);
    m113.SetIdlerVisualizationType(VisualizationType::PRIMITIVES);
    m113.SetSuspensionVisualizationType(VisualizationType::PRIMITIVES);
    m113.SetTrackShoeVisualizationType(VisualizationType::PRIMITIVES);

    auto& vehicle = m113.GetVehicle();

    // ---------------
    // Simulation loop
    // ---------------

#ifdef CHRONO_VSG
    // Create the run-time visualization system
    auto vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
    vis->AttachSystem(system);
    vis->SetWindowTitle("M113");
    vis->SetWindowSize(1280, 720);
    vis->AddCamera(ChVector3d(0, -10, 0), ChVector3d(0, 0, 0));
    vis->SetCameraVertical(CameraVerticalDir::Z);
    vis->Initialize();
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

    bool set_time = true;

    while (time < time_end) {
        // Driver inputs
        DriverInputs driver_inputs = {0, 0, 0, 0};
        {
            double out_speed = m_speedPIDVehicle.Advance(vehicle.GetRefFrame(), target_speed, time, time_step);

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
                    time_end = std::max(time + 1, time_end);
                    set_time = false;
                }
            }
        }

        // Output
        if (sim_frame == next_out_frame) {
            std::cout << "write: " << out_frame << std::endl;
            forces.resize(0);
            torques.resize(0);

            system->CalculateContactForces();

            // get total force and torques on tracks and chassis
            forces.push_back(system->GetBodyContactForce(vehicle.GetChassisBody()));
            torques.push_back(system->GetBodyContactTorque(vehicle.GetChassisBody()));

            for (size_t i = 0; i < vehicle.GetNumTrackShoes(LEFT); i++) {
                forces.push_back(system->GetBodyContactForce(vehicle.GetTrackShoe(LEFT, i)->GetShoeBody()));
                torques.push_back(system->GetBodyContactTorque(vehicle.GetTrackShoe(LEFT, i)->GetShoeBody()));
            }
            for (size_t i = 0; i < vehicle.GetNumTrackShoes(RIGHT); i++) {
                forces.push_back(system->GetBodyContactForce(vehicle.GetTrackShoe(RIGHT, i)->GetShoeBody()));
                torques.push_back(system->GetBodyContactTorque(vehicle.GetTrackShoe(RIGHT, i)->GetShoeBody()));
            }

            DumpFluidData(system, data_output_path + "data_" + std::to_string(out_frame) + ".dat", true);
            DumpAllObjectsWithGeometryPovray(system, data_output_path + "vehicle_" + std::to_string(out_frame) + ".dat",
                                             true);
            WriteTrackedVehicleData(m113, driver_inputs, forces, torques, data_output_path + "stats_" + std::to_string(out_frame) + ".dat");

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

        // Release the vehicle chassis at the end of the hold time.
        if (vehicle.GetChassisBody()->IsFixed() && time > time_hold) {
            vehicle.GetChassisBody()->SetFixed(false);
        }

        m113.Synchronize(time, driver_inputs);
        m113.Advance(time_step);

#ifdef CHRONO_VSG
        if (!vis->Run())
            break;
        vis->Render();
#endif

        system->DoStepDynamics(time_step);
        time += time_step;
        sim_frame++;
        exec_time += system->GetTimerStep();
        num_contacts += system->GetNumContacts();
    }

    FinishFluidDataOutput();

    // Final stats
    std::cout << "==================================" << std::endl;
    std::cout << "Simulation time:   " << exec_time << std::endl;

    return 0;
}
