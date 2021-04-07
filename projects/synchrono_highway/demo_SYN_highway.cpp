// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Aaron Young
// =============================================================================
//
// Basic demonstration of multiple wheeled vehicles in a single simulation using
// the SynChrono wrapper
//
// =============================================================================

#include <chrono>

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChIrrGuiDriver.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/utils/ChUtilsJSON.h"

#include "chrono_vehicle/wheeled_vehicle/utils/ChWheeledVehicleIrrApp.h"
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"

#include "chrono_synchrono/SynChronoManager.h"
#include "chrono_synchrono/SynConfig.h"
#include "chrono_synchrono/agent/SynWheeledVehicleAgent.h"
#include "chrono_synchrono/communication/mpi/SynMPICommunicator.h"
#include "chrono_synchrono/utils/SynDataLoader.h"
#include "chrono_synchrono/utils/SynLog.h"

#include "chrono_sensor/ChCameraSensor.h"
#include "chrono_sensor/ChSensorManager.h"
#include "chrono_sensor/filters/ChFilterSave.h"
#include "chrono_sensor/filters/ChFilterVisualize.h"

#include "chrono_thirdparty/cxxopts/ChCLI.h"

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::irrlicht;
using namespace chrono::synchrono;
using namespace chrono::vehicle;
using namespace chrono::sensor;

// =============================================================================

// Initial vehicle location and orientation
ChVector<> initLoc(907, -230, -64.8);
ChQuaternion<> initRot(1, 0, 0, 0);

// Visualization type for vehicle parts (PRIMITIVES, MESH, or NONE)
VisualizationType chassis_vis_type = VisualizationType::MESH;
VisualizationType suspension_vis_type = VisualizationType::PRIMITIVES;
VisualizationType steering_vis_type = VisualizationType::PRIMITIVES;
VisualizationType wheel_vis_type = VisualizationType::MESH;
VisualizationType tire_vis_type = VisualizationType::MESH;

// Type of tire model
TireModelType tire_model = TireModelType::TMEASY;

// Type of vehicle
enum VehicleType { SEDAN, AUDI, TRUCK, VAN, SUV, CITYBUS };

// Point on chassis tracked by the camera
ChVector<> trackPoint(0.0, 0.0, 1.75);

// Contact method
ChContactMethod contact_method = ChContactMethod::NSC;

// Simulation step sizes
double step_size = 2e-3;

// Simulation end time
double end_time = 1000;

// How often SynChrono state messages are interchanged
double heartbeat = 1e-2;  // 100[Hz]

// Time interval between two render frames
double render_step_size = 1.0 / 50;  // FPS = 50

bool save = false;

std::string demo_data_path = std::string(STRINGIFY(HIGHWAY_DATA_DIR));

// =============================================================================

// Forward declares for straight forward helper functions
void LogCopyright(bool show);
void AddCommandLineOptions(ChCLI& cli);
void GetVehicleModelFiles(VehicleType type,
                          std::string& vehicle,
                          std::string& powertrain,
                          std::string& tire,
                          std::string& zombie,
                          double& cam_distance);

void AddSceneMeshes(ChSystem* chsystem, RigidTerrain* terrain);

class IrrAppWrapper {
  public:
    IrrAppWrapper(std::shared_ptr<ChWheeledVehicleIrrApp> app = nullptr) : app(app) {}

    void Synchronize(const std::string& msg, const ChDriver::Inputs& driver_inputs) {
        if (app)
            app->Synchronize(msg, driver_inputs);
    }

    void Advance(double step) {
        if (app)
            app->Advance(step);
    }

    void Render() {
        if (app) {
            app->BeginScene(true, true, irr::video::SColor(255, 140, 161, 192));
            app->DrawAll();
            app->EndScene();
        }
    }

    void Set(std::shared_ptr<ChWheeledVehicleIrrApp> app) { this->app = app; }
    bool IsOk() { return app ? app->GetDevice()->run() : true; }

    std::shared_ptr<ChWheeledVehicleIrrApp> app;
};

class DriverWrapper : public ChDriver {
  public:
    DriverWrapper(ChVehicle& vehicle) : ChDriver(vehicle) {}

    /// Update the state of this driver system at the specified time.
    virtual void Synchronize(double time) override {
        if (irr_driver) {
            irr_driver->Synchronize(time);
            m_throttle = irr_driver->GetThrottle();
            m_steering = irr_driver->GetSteering();
            m_braking = irr_driver->GetBraking();
        }
    }

    /// Advance the state of this driver system by the specified time step.
    virtual void Advance(double step) override {
        if (irr_driver)
            irr_driver->Advance(step);
    }

    void Set(std::shared_ptr<ChIrrGuiDriver> irr_driver) { this->irr_driver = irr_driver; }

    std::shared_ptr<ChIrrGuiDriver> irr_driver;
};

// =============================================================================

int main(int argc, char* argv[]) {
    // -----------------------
    // Create SynChronoManager
    // -----------------------
    auto communicator = chrono_types::make_shared<SynMPICommunicator>(argc, argv);
    int node_id = communicator->GetRank();
    int num_nodes = communicator->GetNumRanks();
    SynChronoManager syn_manager(node_id, num_nodes, communicator);

    // SetChronoDataPath(CHRONO_DATA_DIR);
    // vehicle::SetDataPath(CHRONO_DATA_DIR + std::string("vehicle/"));
    // synchrono::SetDataPath(CHRONO_DATA_DIR + std::string("synchrono/"));

    // all the demo data will be in user-specified location
    SetChronoDataPath(demo_data_path);
    vehicle::SetDataPath(demo_data_path + std::string("/vehicles/"));
    synchrono::SetDataPath(demo_data_path + std::string("/synchrono/"));

    // Copyright
    LogCopyright(node_id == 0);

    // -----------------------------------------------------
    // CLI SETUP - Get most parameters from the command line
    // -----------------------------------------------------

    ChCLI cli(argv[0]);

    AddCommandLineOptions(cli);
    if (!cli.Parse(argc, argv, node_id == 0))
        return 0;

    // Normal simulation options
    step_size = cli.GetAsType<double>("step_size");
    end_time = cli.GetAsType<double>("end_time");
    heartbeat = cli.GetAsType<double>("heartbeat");
    save = cli.GetAsType<bool>("save");

    // Change SynChronoManager settings
    syn_manager.SetHeartbeat(heartbeat);

    // --------------
    // Create systems
    // --------------

    // Adjust position of each vehicle so they aren't on top of each other
    initLoc.x() += node_id * 5;

    // Get the vehicle JSON filenames
    double cam_distance;
    std::string vehicle_filename, powertrain_filename, tire_filename, zombie_filename;
    GetVehicleModelFiles((VehicleType)cli.GetAsType<int>("vehicle"), vehicle_filename, powertrain_filename,
                         tire_filename, zombie_filename, cam_distance);

    // Create the vehicle, set parameters, and initialize
    WheeledVehicle vehicle(vehicle_filename, contact_method);
    vehicle.Initialize(ChCoordsys<>(initLoc, initRot));
    vehicle.GetChassis()->SetFixed(false);
    vehicle.SetChassisVisualizationType(chassis_vis_type);
    vehicle.SetSuspensionVisualizationType(suspension_vis_type);
    vehicle.SetSteeringVisualizationType(steering_vis_type);
    vehicle.SetWheelVisualizationType(wheel_vis_type);

    // Create and initialize the powertrain system
    auto powertrain = ReadPowertrainJSON(powertrain_filename);
    vehicle.InitializePowertrain(powertrain);

    // Create and initialize the tires
    for (auto& axle : vehicle.GetAxles()) {
        for (auto& wheel : axle->GetWheels()) {
            auto tire = ReadTireJSON(tire_filename);
            vehicle.InitializeTire(tire, wheel, tire_vis_type);
        }
    }

    // Add vehicle as an agent and initialize SynChronoManager
    auto agent = chrono_types::make_shared<SynWheeledVehicleAgent>(&vehicle, zombie_filename);
    syn_manager.AddAgent(agent);
    syn_manager.Initialize(vehicle.GetSystem());

    RigidTerrain terrain(vehicle.GetSystem());
    AddSceneMeshes(vehicle.GetSystem(), &terrain);

    MaterialInfo minfo;  // values from RigidPlane.json
    minfo.mu = 0.9;      // coefficient of friction
    minfo.cr = 0.01;     // coefficient of restitution
    minfo.Y = 2e7;       // Young's modulus
    minfo.nu = 0.3;      // Poisson ratio
    minfo.kn = 2e5;      // normal stiffness
    minfo.gn = 40.0;     // normal viscous damping
    minfo.kt = 2e5;      // tangential stiffness
    minfo.gt = 20.0;     // tangential viscous damping
    auto patch_mat = minfo.CreateMaterial(contact_method);
    auto patch = terrain.AddPatch(patch_mat, ChVector<>({0, 0, -65.554}), ChVector<>({0, 0, 1}), 10000.0, 10000.0, 2,
                                  false, 1, false);
    terrain.Initialize();

    // Create the vehicle Irrlicht interface
    IrrAppWrapper app;
    DriverWrapper driver(vehicle);
    if (cli.HasValueInVector<int>("irr", node_id)) {
        auto temp_app = chrono_types::make_shared<ChWheeledVehicleIrrApp>(&vehicle, L"SynChrono Wheeled Vehicle Demo");
        temp_app->SetSkyBox();
        temp_app->AddTypicalLights(irr::core::vector3df(30.f, -30.f, 100.f), irr::core::vector3df(30.f, 50.f, 100.f),
                                   250, 130);
        temp_app->SetChaseCamera(trackPoint, cam_distance, 0.5);
        temp_app->SetTimestep(step_size);
        temp_app->AssetBindAll();
        temp_app->AssetUpdateAll();

        // Create the interactive driver system
        auto irr_driver = chrono_types::make_shared<ChIrrGuiDriver>(*temp_app);

        // optionally force the gui driver to use keyboard rather than joystick
        if (cli.HasValueInVector<int>("keyboard", node_id))
            irr_driver->SetInputMode(ChIrrGuiDriver::KEYBOARD);

        // Set the time response for steering and throttle keyboard inputs.
        // double steering_time = 1.0;  // time to go from 0 to +1 (or from 0 to -1)
        // double throttle_time = 1.0;  // time to go from 0 to +1
        // double braking_time = 0.3;   // time to go from 0 to +1
        // irr_driver->SetSteeringDelta(render_step_size / steering_time);
        // irr_driver->SetThrottleDelta(render_step_size / throttle_time);
        // irr_driver->SetBrakingDelta(render_step_size / braking_time);
        irr_driver->SetSteeringDelta(0.02);
        irr_driver->SetThrottleDelta(0.02);
        irr_driver->SetBrakingDelta(0.06);
        irr_driver->Initialize();

        app.Set(temp_app);
        driver.Set(irr_driver);
    }

    // add a sensor manager
    auto manager = chrono_types::make_shared<ChSensorManager>(vehicle.GetSystem());
    float brightness = 0.5f;
    manager->scene->AddPointLight({1000, 1000, 1000}, {brightness, brightness, brightness}, 10000);
    manager->scene->AddPointLight({-1000, 1000, 1000}, {brightness, brightness, brightness}, 10000);
    manager->scene->AddPointLight({1000, -1000, 1000}, {brightness, brightness, brightness}, 10000);
    manager->scene->AddPointLight({-1000, -1000, 1000}, {brightness, brightness, brightness}, 10000);

    // add a camera to the vehicle
    auto camera = chrono_types::make_shared<ChCameraSensor>(
        vehicle.GetChassisBody(),                                            // body camera is attached to
        60.f,                                                                // update rate in Hz
        chrono::ChFrame<double>({-12, 0, 3}, Q_from_AngAxis(0, {0, 1, 0})),  // offset pose
        1920,                                                                // image width
        1080,                                                                // image height
        3.14 / 2, 1);
    // camera->SetLag(1 / 30.f);
    // camera->SetName("Camera Sensor");
    camera->PushFilter(chrono_types::make_shared<ChFilterVisualize>(1280, 720));
    if (save)
        camera->PushFilter(chrono_types::make_shared<ChFilterSave>("DEMO_OUTPUT/cam/"));
    manager->AddSensor(camera);

    // ---------------
    // Simulation loop
    // ---------------
    // Number of simulation steps between miscellaneous events
    int render_steps = (int)std::ceil(render_step_size / step_size);

    // Initialize simulation frame counters
    int step_number = 0;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    float orbit_radius = 8.f;
    float orbit_rate = 1;

    while (app.IsOk() && syn_manager.IsOk()) {
        double time = vehicle.GetSystem()->GetChTime();

        // std::cout << "t=" << time << std::endl;

        // End simulation
        if (time >= end_time)
            break;

        // Render scene
        if (step_number % render_steps == 0)
            app.Render();

        // Get driver inputs
        ChDriver::Inputs driver_inputs = driver.GetInputs();

        // Update modules (process inputs from other modules)
        syn_manager.Synchronize(time);  // Synchronize between nodes
        driver.Synchronize(time);
        vehicle.Synchronize(time, driver_inputs, terrain);
        terrain.Synchronize(time);
        app.Synchronize("", driver_inputs);

        // Advance simulation for one timestep for all modules
        driver.Advance(step_size);
        vehicle.Advance(step_size);
        terrain.Advance(step_size);
        app.Advance(step_size);

        // update sensors
        // camera->SetOffsetPose(
        //     chrono::ChFrame<double>({-orbit_radius * cos(time * orbit_rate), -orbit_radius * sin(time * orbit_rate), 3},
        //                             Q_from_AngAxis(time * orbit_rate, {0, 0, 1})));
        manager->Update();

        // Increment frame number
        step_number++;

        // Log clock time
        if (step_number % 500 == 0 && node_id == 0) {
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            auto time_span = std::chrono::duration_cast<std::chrono::seconds>(end - start);

            SynLog() << (time_span.count()) / time << "\n";
        }
    }
    // Properly shuts down other ranks when one rank ends early
    syn_manager.QuitSimulation();

    return 0;
}

void LogCopyright(bool show) {
    if (!show)
        return;

    SynLog() << "Copyright (c) 2020 projectchrono.org\n";
    SynLog() << "Chrono version: " << CHRONO_VERSION << "\n\n";
}

void AddCommandLineOptions(ChCLI& cli) {
    // Standard demo options
    cli.AddOption<double>("Simulation", "s,step_size", "Step size", std::to_string(step_size));
    cli.AddOption<double>("Simulation", "e,end_time", "End time", std::to_string(end_time));
    cli.AddOption<double>("Simulation", "b,heartbeat", "Heartbeat", std::to_string(heartbeat));
    cli.AddOption<bool>("Simulation", "save", "save", std::to_string(save));

    // Irrlicht options
    cli.AddOption<std::vector<int>>("Irrlicht", "i,irr", "Nodes for irrlicht usage", "-1");
    cli.AddOption<std::vector<int>>("Keyboard", "k,keyboard", "Force irrlicht driver into keyboard control", "-1");

    // Other options
    cli.AddOption<int>("Demo", "v,vehicle", "Vehicle Options [0-4]: Sedan, HMMWV, UAZ, CityBus, MAN", "0");
}

void GetVehicleModelFiles(VehicleType type,
                          std::string& vehicle,
                          std::string& powertrain,
                          std::string& tire,
                          std::string& zombie,
                          double& cam_distance) {
    switch (type) {
        case VehicleType::SEDAN:
            vehicle = vehicle::GetDataFile("sedan/vehicle/Sedan_Vehicle.json");
            powertrain = vehicle::GetDataFile("sedan/powertrain/Sedan_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("sedan/tire/Sedan_TMeasyTire.json");
            zombie = vehicle::GetDataFile("sedan/Sedan.json");
            cam_distance = 6.0;
            break;
        case VehicleType::AUDI:
            vehicle = vehicle::GetDataFile("audi/json/audi_Vehicle.json");
            powertrain = vehicle::GetDataFile("audi/json/audi_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("audi/json/audi_TMeasyTire.json");
            zombie = vehicle::GetDataFile("audi/json/audi.json");
            cam_distance = 6.0;
            break;
        case VehicleType::TRUCK:
            // vehicle::SetDataPath("/mnt/data/code/chrono/data/vehicle/");
            // synchrono::SetDataPath("/mnt/data/code/chrono/data/synchrono/");

            // vehicle = vehicle::GetDataFile("MAN_Kat1/vehicle/MAN_7t_Vehicle_6WD.json");
            // powertrain = vehicle::GetDataFile("MAN_Kat1/powertrain/MAN_7t_SimpleCVTPowertrain.json");
            // tire = vehicle::GetDataFile("MAN_Kat1/tire/MAN_5t_TMeasyTire.json");
            // zombie = synchrono::GetDataFile("vehicle/MAN_8WD.json");

            vehicle = vehicle::GetDataFile("truck/json/truck_Vehicle.json");
            powertrain = vehicle::GetDataFile("truck/json/truck_SimpleCVTPowertrain.json");
            tire = vehicle::GetDataFile("truck/json/truck_TMeasyTire.json");
            zombie = vehicle::GetDataFile("truck/json/truck.json");
            cam_distance = 14.0;
            break;
        case VehicleType::VAN:
            vehicle = vehicle::GetDataFile("van/json/van_Vehicle.json");
            powertrain = vehicle::GetDataFile("van/json/van_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("van/json/van_TMeasyTire.json");
            zombie = vehicle::GetDataFile("van/json/van.json");
            cam_distance = 12.0;
            break;
        case VehicleType::SUV:
            vehicle = vehicle::GetDataFile("suv/json/suv_Vehicle.json");
            powertrain = vehicle::GetDataFile("suv/json/suv_ShaftsPowertrain.json");
            tire = vehicle::GetDataFile("suv/json/suv_TMeasyTire.json");
            zombie = vehicle::GetDataFile("suv/json/suv.json");
            cam_distance = 6.0;
            break;
        case VehicleType::CITYBUS:
            vehicle = vehicle::GetDataFile("citybus/vehicle/CityBus_Vehicle.json");
            powertrain = vehicle::GetDataFile("citybus/powertrain/CityBus_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("citybus/tire/CityBus_TMeasyTire.json");
            zombie = vehicle::GetDataFile("citybus/CityBus.json");
            cam_distance = 14.0;
            break;
    }
}

void AddSceneMeshes(ChSystem* chsystem, RigidTerrain* terrain) {
    // load all meshes in input file, using instancing where possible
    std::string base_path = GetChronoDataFile("/Environments/SanFrancisco/components/");
    std::string input_file = base_path + "instance_map_02.csv";
    // std::string input_file = base_path + "instance_map_roads_only.csv";

    std::ifstream infile(input_file);
    if (!infile.is_open())
        throw std::runtime_error("Could not open file " + input_file);
    std::string line, col;
    std::vector<std::string> result;

    std::unordered_map<std::string, std::shared_ptr<ChTriangleMeshConnected>> mesh_map;

    int mesh_offset = 0;
    int num_meshes = 20000;
    if (infile.good()) {
        int mesh_count = 0;
        int mesh_limit = mesh_offset + num_meshes;
        while (std::getline(infile, line) && mesh_count < mesh_limit) {
            if (mesh_count < mesh_offset) {
                mesh_count++;
            } else {
                mesh_count++;
                result.clear();
                std::stringstream ss(line);
                while (std::getline(ss, col, ',')) {
                    result.push_back(col);
                }
                // std::cout << "Name: " << result[0] << ", mesh: " << result[1] << std::endl;
                std::string mesh_name = result[0];
                std::string mesh_obj = base_path + result[1] + ".obj";

                // std::cout << mesh_name << std::endl;

                // check if mesh is in map
                bool instance_found = false;
                std::shared_ptr<ChTriangleMeshConnected> mmesh;
                if (mesh_map.find(mesh_obj) != mesh_map.end()) {
                    mmesh = mesh_map[mesh_obj];
                    instance_found = true;
                } else {
                    mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
                    mmesh->LoadWavefrontMesh(mesh_obj, false, true);
                    mesh_map[mesh_obj] = mmesh;
                }

                ChVector<double> pos = {std::stod(result[2]), std::stod(result[3]), std::stod(result[4])};
                ChQuaternion<double> rot = {std::stod(result[5]), std::stod(result[6]), std::stod(result[7]),
                                            std::stod(result[8])};
                ChVector<double> scale = {std::stod(result[9]), std::stod(result[10]), std::stod(result[11])};

                // if not road, only add visualization with new pos,rot,scale
                auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
                trimesh_shape->SetMesh(mmesh);
                trimesh_shape->SetName(mesh_name);
                trimesh_shape->SetStatic(true);
                trimesh_shape->SetScale(scale);

                auto mesh_body = chrono_types::make_shared<ChBody>();
                mesh_body->SetPos(pos);
                mesh_body->SetRot(rot);
                mesh_body->AddAsset(trimesh_shape);
                mesh_body->SetBodyFixed(true);
                chsystem->Add(mesh_body);
            }
        }
        ChVector<> pos = {0, 0, 100};
        std::cout << "Terrrain height at <" << pos.x() << "," << pos.y() << "," << pos.z()
                  << ">: " << terrain->GetHeight(pos) << std::endl;
        std::cout << "Total meshes: " << mesh_count - mesh_offset << " | Unique meshes: " << mesh_map.size()
                  << std::endl;
    }
}