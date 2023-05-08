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
// Black-box program for using an external optimization program for tuning
// parameters of a PID steering controller.
//
// =============================================================================

#include <vector>
#include <valarray>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "chrono/geometry/ChLineBezier.h"
#include "chrono/assets/ChLineShape.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"
#include "chrono_vehicle/wheeled_vehicle/tire/RigidTire.h"
#include "chrono_vehicle/wheeled_vehicle/tire/FialaTire.h"
#include "chrono_vehicle/powertrain/EngineSimple.h"
#include "chrono_vehicle/powertrain/AutomaticTransmissionSimpleMap.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"

#include "chrono_vehicle/driver/ChPathFollowerDriver.h"

using namespace chrono;
using namespace chrono::vehicle;
using namespace geometry;

// =============================================================================
// Global definitions

typedef std::valarray<double> DataArray;

struct Data {
    Data(int n) {
        time.resize(n);
        err_x.resize(n);
        err_y.resize(n);
        err_z.resize(n);
        err_speed.resize(n);
    }

    DataArray time;       // current time
    DataArray err_x;      // x component of vehicle location error
    DataArray err_y;      // y component of vehicle location error
    DataArray err_z;      // z component of vehicle location error
    DataArray err_speed;  // vehicle speed error
};

// Type of tire model
TireModelType tire_model = TireModelType::RIGID;

// Input file names for the path-follower driver model
std::string steering_controller_file("generic/driver/SteeringController.json");
std::string speed_controller_file("generic/driver/SpeedController.json");
std::string path_file("paths/curve.txt");

// Output file name
std::string out_file("results.out");

// JSON file names for vehicle model, tire models, (simple) powertrain, and (rigid) terrain
std::string vehicle_file("generic/vehicle/Vehicle_DoubleWishbones.json");
std::string rigidtire_file("generic/tire/RigidTire.json");
std::string fialatire_file("generic/tire/FialaTire.json");
std::string engine_file("generic/powertrain/EngineSimple.json");
std::string transmission_file("generic/powertrain/AutomaticTransmissionSimpleMap.json");
std::string rigidterrain_file("terrain/RigidPlane.json");

// Initial vehicle position and orientation
ChVector<> initLoc(-125, -125, 0.6);
ChQuaternion<> initRot(1, 0, 0, 0);

// Desired vehicle speed (m/s)
double target_speed = 10;

// Rigid terrain dimensions
double terrainHeight = 0;
double terrainLength = 300.0;  // size in X direction
double terrainWidth = 300.0;   // size in Y direction

// Simulation step size and simulation length
double step_size = 2e-3;        // integration step size
int num_steps_settling = 3000;  // number of steps for settling
int num_steps = 5000;           // number of steps for data colection

// =============================================================================
// Forward declarations

void processData(const utils::CSV_writer& csv, const Data& data);

// =============================================================================
// Main driver program

int main(int argc, char* argv[]) {
    // Set path to Chrono and Chrono::Vehicle data directories
    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_VEHICLE_DATA_DIR);

    // Create and initialize the vehicle system
    WheeledVehicle vehicle(vehicle::GetDataFile(vehicle_file));
    vehicle.Initialize(ChCoordsys<>(initLoc, initRot));

    // Create the terrain
    RigidTerrain terrain(vehicle.GetSystem(), vehicle::GetDataFile(rigidterrain_file));

    // Create and initialize the powertrain system
    auto engine = chrono_types::make_shared<EngineSimple>(vehicle::GetDataFile(engine_file));
    auto transmission = chrono_types::make_shared<AutomaticTransmissionSimpleMap>(vehicle::GetDataFile(transmission_file));
    auto powertrain = chrono_types::make_shared<ChPowertrainAssembly>(engine, transmission);
    vehicle.InitializePowertrain(powertrain);

    // Create and initialize the tires
    for (auto& axle : vehicle.GetAxles()) {
        switch (tire_model) {
            default:
            case TireModelType::RIGID: {
                auto tireL = chrono_types::make_shared<RigidTire>(vehicle::GetDataFile(rigidtire_file));
                auto tireR = chrono_types::make_shared<RigidTire>(vehicle::GetDataFile(rigidtire_file));
                vehicle.InitializeTire(tireL, axle->m_wheels[0], VisualizationType::MESH);
                vehicle.InitializeTire(tireR, axle->m_wheels[1], VisualizationType::MESH);
                break;
            }
            case TireModelType::FIALA: {
                auto tireL = chrono_types::make_shared<FialaTire>(vehicle::GetDataFile(fialatire_file));
                auto tireR = chrono_types::make_shared<FialaTire>(vehicle::GetDataFile(fialatire_file));
                vehicle.InitializeTire(tireL, axle->m_wheels[0], VisualizationType::MESH);
                vehicle.InitializeTire(tireR, axle->m_wheels[1], VisualizationType::MESH);
                break;
            }
        }
    }

    // Create the driver system
    auto path = ChBezierCurve::read(vehicle::GetDataFile(path_file));
    ChPathFollowerDriver driver(vehicle, vehicle::GetDataFile(steering_controller_file),
                                vehicle::GetDataFile(speed_controller_file), path, "my_path", target_speed);

    // Create a path tracker to keep track of the error in vehicle location.
    ChBezierCurveTracker tracker(path);

    // ---------------
    // Simulation loop
    // ---------------

    // Initialize data collectors
    utils::CSV_writer csv("\t");
    csv.stream().setf(std::ios::scientific | std::ios::showpos);
    csv.stream().precision(6);

    Data data(num_steps);

    std::cout << "Total number of steps:  " << num_steps_settling + num_steps << std::endl;
    for (int it = 0; it < num_steps_settling + num_steps; it++) {
        bool settling = (it < num_steps_settling);

        // Collect data
        if (!settling) {
            const ChVector<> sentinel = driver.GetSteeringController().GetSentinelLocation();
            const ChVector<> target = driver.GetSteeringController().GetTargetLocation();
            const ChVector<> vehicle_location = vehicle.GetPos();
            ChVector<> vehicle_target;
            tracker.calcClosestPoint(vehicle_location, vehicle_target);
            ChVector<> vehicle_err = vehicle_target - vehicle_location;
            double speed_err = target_speed - vehicle.GetSpeed();

            csv << vehicle.GetChTime() << vehicle_location << vehicle_target << vehicle_err << speed_err << std::endl;

            int id = it - num_steps_settling;
            data.time[id] = vehicle.GetChTime();
            data.err_x[id] = vehicle_err.x();
            data.err_y[id] = vehicle_err.y();
            data.err_z[id] = vehicle_err.z();
            data.err_speed = speed_err;
        }

        // Collect output data from modules (for inter-module communication)
        DriverInputs driver_inputs = driver.GetInputs();
        if (settling) {
            driver_inputs.m_throttle = 0;
            driver_inputs.m_steering = 0;
            driver_inputs.m_braking = 0;
        }

        // Update modules (process inputs from other modules)
        double time = vehicle.GetChTime();
        driver.Synchronize(time);
        vehicle.Synchronize(time, driver_inputs, terrain);
        terrain.Synchronize(time);

        // Advance simulation for one timestep for all modules
        driver.Advance(step_size);
        vehicle.Advance(step_size);
        terrain.Advance(step_size);

        std::cout << '\r' << std::fixed << std::setprecision(6) << time << "  (" << it << ")" << std::flush;
    }

    processData(csv, data);

    return 0;
}

// =============================================================================
// Simulation data post-processing

void processData(const utils::CSV_writer& csv, const Data& data) {
    // Optionally, write simulation results to file for external post-processing
    csv.write_to_file(out_file);

    // Alternatively, post-process simulation results here and write out results
    DataArray loc_err_norm2 = data.err_x * data.err_x + data.err_y * data.err_y + data.err_z * data.err_z;
    double loc_L2_norm = std::sqrt(loc_err_norm2.sum());
    double loc_RMS_norm = std::sqrt(loc_err_norm2.sum() / num_steps);
    double loc_INF_norm = std::sqrt(loc_err_norm2.max());

    std::cout << "|location err|_L2 =  " << loc_L2_norm << std::endl;
    std::cout << "|location err|_RMS = " << loc_RMS_norm << std::endl;
    std::cout << "|location err|_INF = " << loc_INF_norm << std::endl;

    ////std::ofstream ofile(out_file.c_str());
    ////ofile << loc_L2_norm << std::endl;
    ////ofile.close();

    double speed_L2_norm = std::sqrt((data.err_speed * data.err_speed).sum());
    double speed_RMS_norm = std::sqrt((data.err_speed * data.err_speed).sum() / num_steps);
    double speed_INF_norm = std::abs(data.err_speed).max();

    std::cout << "|speed err|_L2 =  " << speed_L2_norm << std::endl;
    std::cout << "|speed err|_RMS = " << speed_RMS_norm << std::endl;
    std::cout << "|speed err|_INF = " << speed_INF_norm << std::endl;

    ////std::ofstream ofile(out_file.c_str());
    ////ofile << speed_L2_norm << std::endl;
    ////ofile.close();
}
