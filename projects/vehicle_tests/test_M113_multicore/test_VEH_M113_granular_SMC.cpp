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
#include "chrono/input_output/ChUtilsInputOutput.h"

// Chrono::Multicore header files
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChSystemDescriptorMulticore.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

// Chrono utility header files
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/input_output/ChUtilsInputOutput.h"

// Chrono vehicle header files
#include "chrono_vehicle/ChVehicleDataPath.h"
#include "chrono_vehicle/driver/ChDataDriver.h"

// M113 model header files
#include "chrono_models/vehicle/m113/powertrain/M113_EngineSimpleMap.h"
#include "chrono_models/vehicle/m113/powertrain/M113_AutomaticTransmissionSimpleMap.h"
#include "chrono_models/vehicle/m113/M113_Vehicle.h"

#include "chrono_thirdparty/filesystem/path.h"

// Utilities
#include "../../utils.h"

using namespace chrono;
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
int tag_particles = 0;
double r_g = 18e-3;
double rho_g = 2500;
double coh_pressure = 3e4;
float mu_g = 0.9f;

double vol_g = (4.0 / 3) * CH_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector3d inertia_g = 0.4 * mass_g * r_g * r_g * ChVector3d(1, 1, 1);
double coh_force = CH_PI * r_g * r_g * coh_pressure;

// -----------------------------------------------------------------------------
// Specification of the vehicle model
// -----------------------------------------------------------------------------

// Initial vehicle position and orientation
ChVector3d initLoc(-hdimX + 4.5, 0, 1.0);
ChQuaternion<> initRot(1, 0, 0, 0);

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Total simulation duration.
double time_end = 7;

// Duration of the "hold time" (vehicle chassis fixed and no driver inputs).
// This can be used to allow the granular material to settle.
double time_hold = 0.2;

// Solver parameters
double time_step = 5e-5;
double tolerance = 1e-5;

int max_iteration_bilateral = 1000;

// Periodically monitor maximum bilateral constraint violation
bool monitor_bilaterals = false;
int bilateral_frame_interval = 100;

// Output directories
bool povray_output = true;

const std::string out_dir = "../M113_MULTICORE_SMC";
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
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(0.0f);
    mat_g->SetYoungModulus(8e5f);
    mat_g->SetPoissonRatio(0.3f);
    mat_g->SetAdhesion(static_cast<float>(coh_force));
    mat_g->SetKn(1.0e6f);
    mat_g->SetGn(6.0e1f);
    mat_g->SetKt(4.0e5f);
    mat_g->SetGt(4.0e1f);

    // Create a particle generator and a mixture entirely made out of spheres
    double r = 1.01 * r_g;
    chrono::utils::ChPDSampler<double> sampler(2 * r);
    chrono::utils::ChGenerator gen(system);
    std::shared_ptr<chrono::utils::ChMixtureIngredient> m1 =
        gen.AddMixtureIngredient(chrono::utils::MixtureType::SPHERE, 1.0);
    m1->SetDefaultMaterial(mat_g);
    m1->SetDefaultDensity(rho_g);
    m1->SetDefaultSize(r_g);

    // Set starting value for body identifiers
    gen.SetStartTag(tag_particles);

    // Create particles in layers until reaching the desired number of particles
    ChVector3d hdims(hdimX - r, hdimY - r, 0);
    ChVector3d center(0, 0, 2 * r);

    double layerCount = 0;
    while (layerCount < numLayers) {
        gen.CreateObjectsBox(sampler, center, hdims);
        center.z() += 2 * r;
        layerCount++;
    }

    std::cout << "Created " << gen.GetTotalNumBodies() << " particles." << std::endl;

    return center.z();
}

// =============================================================================
int main(int argc, char* argv[]) {
    // -------------------------------------------------------
    // Set path to Chrono and Chrono::Vehicle data directories
    // -------------------------------------------------------

    SetChronoDataPath(CHRONO_DATA_DIR);
    SetVehicleDataPath(CHRONO_VEHICLE_DATA_DIR);

    // ---------------------------------
    // Create system and modify settings
    // ---------------------------------

    std::cout << "Create multicore SMC system" << std::endl;
    ChSystemMulticoreSMC system;
    system.SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);
    system.SetGravitationalAcceleration(ChVector3d(0, 0, -9.80665));

    // Set number of threads
    system.SetNumThreads(std::min(threads, ChOMP::GetNumProcs()));

    // Set solver parameters
    system.GetSettings()->solver.use_full_inertia_tensor = false;
    system.GetSettings()->solver.tolerance = tolerance;
    system.GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    system.GetSettings()->solver.contact_force_model = ChSystemSMC::Hertz;
    system.GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::OneStep;
    system.GetSettings()->solver.use_material_properties = true;

    system.GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

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
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(1e8f);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(0.4f);

    // Ground body
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    ground->EnableCollision(true);

    // Bottom box
    chrono::utils::AddBoxGeometry(ground.get(), mat_g, ChVector3d(hdimX, hdimY, hthick), ChVector3d(0, 0, -hthick),
                                  ChQuaternion<>(1, 0, 0, 0), true);
    if (terrain_type == GRANULAR_TERRAIN) {
        // Front box
        chrono::utils::AddBoxGeometry(ground.get(), mat_g, ChVector3d(hthick, hdimY, hdimZ + hthick),
                                      ChVector3d(hdimX + hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0),
                                      visible_walls);
        // Rear box
        chrono::utils::AddBoxGeometry(ground.get(), mat_g, ChVector3d(hthick, hdimY, hdimZ + hthick),
                                      ChVector3d(-hdimX - hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0),
                                      visible_walls);
        // Left box
        chrono::utils::AddBoxGeometry(ground.get(), mat_g, ChVector3d(hdimX, hthick, hdimZ + hthick),
                                      ChVector3d(0, hdimY + hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0),
                                      visible_walls);
        // Right box
        chrono::utils::AddBoxGeometry(ground.get(), mat_g, ChVector3d(hdimX, hthick, hdimZ + hthick),
                                      ChVector3d(0, -hdimY - hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0),
                                      visible_walls);
    }

    system.AddBody(ground);

    // Create the granular material.
    double vertical_offset = 0;

    if (terrain_type == GRANULAR_TERRAIN) {
        vertical_offset = CreateParticles(&system);
    }

    // --------------------------
    // Construct the M113 vehicle
    // --------------------------

    auto vehicle = chrono_types::make_shared<M113_Vehicle_SinglePin>(true, DrivelineTypeTV::SIMPLE, BrakeType::SIMPLE,
                                                                     false, false, false, &system);
    auto engine = chrono_types::make_shared<M113_EngineSimpleMap>("Engine");
    auto transmission = chrono_types::make_shared<M113_AutomaticTransmissionSimpleMap>("Transmission");
    auto powertrain = chrono_types::make_shared<ChPowertrainAssembly>(engine, transmission);

    vehicle->Initialize(ChCoordsys<>(initLoc + ChVector3d(0.0, 0.0, vertical_offset), initRot));

    // Set visualization type for subsystems
    vehicle->SetChassisVisualizationType(VisualizationType::PRIMITIVES);
    vehicle->SetSprocketVisualizationType(VisualizationType::MESH);
    vehicle->SetIdlerVisualizationType(VisualizationType::MESH);
    vehicle->SetSuspensionVisualizationType(VisualizationType::PRIMITIVES);
    vehicle->SetRoadWheelVisualizationType(VisualizationType::MESH);
    vehicle->SetTrackShoeVisualizationType(VisualizationType::MESH);

    ////vehicle->EnableCollision(TrackCollide::NONE);
    ////vehicle->EnableCollision(TrackCollide::WHEELS_LEFT | TrackCollide::WHEELS_RIGHT);
    ////vehicle->EnableCollision(TrackCollide::ALL & (~TrackCollide::SPROCKET_LEFT) & (~TrackCollide::SPROCKET_RIGHT));

    // Initialize the powertrain system
    vehicle->InitializePowertrain(powertrain);

    // Create the driver system
    MyDriver driver(*vehicle, 0.5);
    driver.Initialize();

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

    chrono::ChWriterCSV csv("\t");
    csv.Stream().setf(std::ios::scientific | std::ios::showpos);
    csv.Stream().precision(6);

    // ---------------
    // Simulation loop
    // ---------------

    // Number of simulation steps between two 3D view render frames
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

#ifdef CHRONO_OPENGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(&system);
    vis.SetWindowTitle("M113");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(0, -10, 0), ChVector3d(0, 0, 0));
    vis.SetCameraVertical(CameraVerticalDir::Z);
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
        DriverInputs driver_inputs = driver.GetInputs();
        vehicle->GetTrackShoeStates(LEFT, shoe_states_left);
        vehicle->GetTrackShoeStates(RIGHT, shoe_states_right);

        const ChVector3d& pos_CG = vehicle->GetChassis()->GetPos();
        ChVector3d vel_CG = vehicle->GetChassisBody()->GetLinVel();
        vel_CG = vehicle->GetChassisBody()->GetCoordsys().TransformDirectionParentToLocal(vel_CG);

        // Vehicle and Control Values
        csv << time << driver_inputs.m_steering << driver_inputs.m_throttle << driver_inputs.m_braking;
        csv << vehicle->GetTrackAssembly(LEFT)->GetSprocket()->GetAxleSpeed()
            << vehicle->GetTrackAssembly(RIGHT)->GetSprocket()->GetAxleSpeed();
        csv << engine->GetMotorSpeed() << engine->GetOutputMotorshaftTorque();
        csv << powertrain->GetOutputTorque() << vehicle->GetDriveline()->GetOutputDriveshaftSpeed();
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
                chrono::utils::WriteVisualizationAssets(&system, filename);
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;

            csv.WriteToFile(out_dir + "/output.dat");
        }

        // Release the vehicle chassis at the end of the hold time.
        if (vehicle->GetChassis()->IsFixed() && time > time_hold) {
            std::cout << std::endl << "Release vehicle t = " << time << std::endl;
            vehicle->GetChassisBody()->SetFixed(false);
        }

        // Update modules (process inputs from other modules)
        driver.Synchronize(time);
        vehicle->Synchronize(time, driver_inputs, shoe_forces_left, shoe_forces_right);

        // Advance simulation for one timestep for all modules
        driver.Advance(time_step);
        vehicle->Advance(time_step);

#ifdef CHRONO_OPENGL
        if (vis.Run())
            vis.Render();
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
        num_contacts += system.GetNumContacts();
    }

    // Final stats
    std::cout << "==================================" << std::endl;
    std::cout << "Simulation time:   " << exec_time << std::endl;
    std::cout << "Number of threads: " << threads << std::endl;

    csv.WriteToFile(out_dir + "/output.dat");

    return 0;
}
