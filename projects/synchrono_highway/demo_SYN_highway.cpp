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
// =============================================================================
//
// Basic demonstration of multiple wheeled vehicles in a single simulation using
// the SynChrono wrapper
//
// =============================================================================

#include <chrono>
#include <functional>

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
#ifdef USE_FAST_DDS
    #include "chrono_synchrono/communication/dds/SynDDSCommunicator.h"
#endif
#include "chrono_synchrono/communication/mpi/SynMPICommunicator.h"
#include "chrono_synchrono/controller/SynControllerFunctions.h"
#include "chrono_synchrono/controller/driver/SynMultiPathDriver.h"
#include "chrono_synchrono/utils/SynDataLoader.h"
#include "chrono_synchrono/utils/SynLog.h"
#include "chrono_vehicle/driver/ChPathFollowerACCDriver.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"

#include "chrono_sensor/ChCameraSensor.h"
#include "chrono_sensor/ChLidarSensor.h"
#include "chrono_sensor/ChSensorManager.h"
#include "chrono_sensor/filters/ChFilterAccess.h"
#include "chrono_sensor/filters/ChFilterLidarNoise.h"
#include "chrono_sensor/filters/ChFilterLidarReduce.h"
#include "chrono_sensor/filters/ChFilterPCfromDepth.h"
#include "chrono_sensor/filters/ChFilterSave.h"
#include "chrono_sensor/filters/ChFilterSavePtCloud.h"
#include "chrono_sensor/filters/ChFilterVisualize.h"
#include "chrono_sensor/filters/ChFilterVisualizePointCloud.h"

#include "chrono_thirdparty/cxxopts/ChCLI.h"

#include "extras/driver/ChCSLDriver.h"
#include "extras/driver/ChLidarWaypointDriver.h"
#include "extras/filters/ChFilterFullScreenVisualize.h"

// =============================================================================

#ifdef USE_FAST_DDS

    // FastDDS
    #undef ALIVE
    #include <fastdds/dds/domain/qos/DomainParticipantQos.hpp>
    #include <fastdds/rtps/transport/UDPv4TransportDescriptor.h>

using namespace eprosima::fastdds::dds;
using namespace eprosima::fastdds::rtps;
using namespace eprosima::fastrtps::rtps;

#endif

// =============================================================================

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::irrlicht;
using namespace chrono::synchrono;
using namespace chrono::vehicle;
using namespace chrono::sensor;

// =============================================================================

// Visualization type for vehicle parts (PRIMITIVES, MESH, or NONE)
VisualizationType chassis_vis_type = VisualizationType::MESH;
VisualizationType suspension_vis_type = VisualizationType::PRIMITIVES;
VisualizationType steering_vis_type = VisualizationType::PRIMITIVES;
VisualizationType wheel_vis_type = VisualizationType::MESH;
VisualizationType tire_vis_type = VisualizationType::MESH;

// Type of tire model
TireModelType tire_model = TireModelType::TMEASY;

// Type of vehicle
enum VehicleType { SEDAN, AUDI, SUV, VAN, TRUCK, CITYBUS };

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

unsigned int leader;
bool save = false;
bool use_fullscreen = false;
bool no_sensing = false;

// Resolution of the CSL 3-monitor setup
const int FS_WIDTH = 3840;
const int FS_HEIGHT = 720;

std::string demo_data_path = std::string(STRINGIFY(HIGHWAY_DATA_DIR));

struct PathVehicleSetup {
    VehicleType vehicle_type;
    ChVector<double> pos;
    ChQuaternion<double> rot;
    std::string path_file;
    double lookahead;
    double speed_gain_p;
};

// starting locations and paths
std::vector<PathVehicleSetup> demo_config = {
    {AUDI, {925.434, -164.87, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/2.txt", 10.0, 0.1},

    {AUDI, {845.534, -131.97, -64.8}, Q_from_AngZ(3.14), "/paths/4.txt", 8.0, 0.1},
    {VAN, {763.334, -131.37, -64.8}, Q_from_AngZ(3.14), "/paths/4.txt", 8.0, 1.0},
    {SUV, {727.834, -158.07, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/4.txt", 8.0, 0.1},
    {SUV, {727.834, -203.57, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/4.txt", 8.0, 0.1},
    {AUDI, {759.734, -225.07, -64.8}, Q_from_AngZ(0), "/paths/4.txt", 8.0, 0.1},
    {SUV, {897.934, -223.27, -64.8}, Q_from_AngZ(0), "/paths/4.txt", 8.0, 0.1},
    {AUDI, {925.434, -199.77, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/4.txt", 8.0, 0.1},
    {AUDI, {897.434, -132.07, -64.8}, Q_from_AngZ(3.14), "/paths/4.txt", 8.0, 0.1},

    {SUV, {925.434, -53.47, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {AUDI, {925.434, 0.47, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {VAN, {925.434, 50.47, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/2.txt", 10.0, 1.0},
    {SUV, {925.434, 75.47, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {AUDI, {903.134, 149.13, -64.8}, Q_from_AngZ(3.14), "/paths/2.txt", 10.0, 0.1},
    {AUDI, {825.134, 149.13, -64.8}, Q_from_AngZ(3.14), "/paths/2.txt", 10.0, 0.1},
    {SUV, {751.234, 148.93, -64.8}, Q_from_AngZ(3.14), "/paths/2.txt", 10.0, 0.1},
    {CITYBUS, {727.834, 124.13, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/2.txt", 10.0, 1.0},
    {SUV, {727.834, 85.13, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {AUDI, {727.834, 40.13, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {SUV, {727.834, -34.27, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {AUDI, {727.834, -100.27, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {AUDI, {727.834, -212.97, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/2.txt", 10.0, 0.1},
    {VAN, {748.234, -225.07, -64.8}, Q_from_AngZ(0), "/paths/2.txt", 10.0, 1.0},
    {AUDI, {855.934, -222.77, -64.8}, Q_from_AngZ(0), "/paths/2.txt", 10.0, 0.1},
    {CITYBUS, {925.634, -214.17, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/2.txt", 10.0, 1.0},

    {AUDI, {867.634, 140.83, -64.8}, Q_from_AngZ(0), "/paths/3.txt", 7.0, 0.1},
    {AUDI, {847.634, 140.83, -64.8}, Q_from_AngZ(0), "/paths/3.txt", 7.0, 0.1},
    {AUDI, {917.234, 116.63, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/3.txt", 7.0, 0.1},
    {SUV, {917.234, 60.63, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/3.txt", 5.0, 0.1},
    {SUV, {917.234, -10.63, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/3.txt", 5.0, 0.1},
    {AUDI, {917.334, -95.67, -64.8}, Q_from_AngZ(-3.14 / 2), "/paths/3.txt", 7.0, 0.1},
    {SUV, {892.334, -120.17, -64.8}, Q_from_AngZ(3.14), "/paths/3.txt", 5.0, 0.1},
    {SUV, {850.334, -120.17, -64.8}, Q_from_AngZ(3.14), "/paths/3.txt", 5.0, 0.1},
    {AUDI, {752.934, -119.47, -64.8}, Q_from_AngZ(3.14), "/paths/3.txt", 7.0, 0.1},
    {SUV, {735.734, -102.97, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/3.txt", 5.0, 0.1},
    {AUDI, {735.734, -75.97, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/3.txt", 5.0, 0.1},
    {AUDI, {735.734, 1.43, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/3.txt", 7.0, 0.1},
    {AUDI, {735.734, 123.63, -64.8}, Q_from_AngZ(3.14 / 2), "/paths/3.txt", 7.0, 0.1},
    {SUV, {755.634, 140.93, -64.8}, Q_from_AngZ(0), "/paths/3.txt", 7.0, 0.1},
    {SUV, {785.634, 140.93, -64.8}, Q_from_AngZ(0), "/paths/3.txt", 7.0, 0.1}};

// =============================================================================

// Forward declares for straight forward helper functions
void LogCopyright(bool show);
void AddCommandLineOptions(ChCLI& cli);
void GetVehicleModelFiles(VehicleType type,
                          std::string& vehicle,
                          std::string& powertrain,
                          std::string& tire,
                          std::string& zombie,
                          ChVector<>& lidar_pos,
                          double& cam_distance);

void AddSceneMeshes(ChSystem* chsystem, RigidTerrain* terrain);

void VehicleProcessMessageCallback(std::shared_ptr<SynMessage> message,
                                   WheeledVehicle& vehicle,
                                   std::shared_ptr<SynWheeledVehicleAgent> agent,
                                   std::shared_ptr<ChPathFollowerACCDriver> driver);

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

// =============================================================================

int main(int argc, char* argv[]) {
    // -----------------------------------------------------
    // CLI SETUP - Get most parameters from the command line
    // -----------------------------------------------------

    ChCLI cli(argv[0]);

    AddCommandLineOptions(cli);
    if (!cli.Parse(argc, argv, true))
        return 0;

        // -----------------------
        // Create SynChronoManager
        // -----------------------
#ifdef USE_FAST_DDS
    int node_id, num_nodes;
    std::shared_ptr<SynCommunicator> communicator;
    if (cli.GetAsType<bool>("dds")) {
        node_id = cli.GetAsType<int>("node_id");
        num_nodes = cli.GetAsType<int>("num_nodes");

        // Set up the qos
        DomainParticipantQos qos;
        qos.name("/syn/node/" + AgentKey(node_id, 0).GetKeyString());

        // Use UDP by default
        auto udp_transport = chrono_types::make_shared<UDPv4TransportDescriptor>();
        udp_transport->maxInitialPeersRange = num_nodes;
        qos.transport().user_transports.push_back(udp_transport);
        qos.transport().use_builtin_transports = false;

        // Set up the initial peers list
        std::vector<std::string> ip_list = {"10.8.0.2", "127.0.0.1"};
        for (auto& ip : ip_list) {
            Locator_t peer;
            IPLocator::setIPv4(peer, ip);
            peer.port = 0;
            qos.wire_protocol().builtin.initialPeersList.push_back(peer);
        }

        auto dds_communicator = chrono_types::make_shared<SynDDSCommunicator>(node_id);
        communicator = dds_communicator;
    } else {
        auto mpi_communicator = chrono_types::make_shared<SynMPICommunicator>(argc, argv);

        node_id = mpi_communicator->GetRank();
        num_nodes = mpi_communicator->GetNumRanks();
        communicator = mpi_communicator;
    }
    SynChronoManager syn_manager(node_id, num_nodes, communicator);
#else
    auto communicator = chrono_types::make_shared<SynMPICommunicator>(argc, argv);

    int node_id = communicator->GetRank();
    int num_nodes = communicator->GetNumRanks();
    SynChronoManager syn_manager(node_id, num_nodes, communicator);
#endif

    // all the demo data will be in user-specified location
    SetChronoDataPath(demo_data_path);
    vehicle::SetDataPath(demo_data_path + std::string("/vehicles/"));
    synchrono::SetDataPath(demo_data_path + std::string("/synchrono/"));

    // Normal simulation options
    step_size = cli.GetAsType<double>("step_size");
    end_time = cli.GetAsType<double>("end_time");
    heartbeat = cli.GetAsType<double>("heartbeat");
    leader = cli.GetAsType<unsigned int>("leader");
    save = cli.GetAsType<bool>("save");
    use_fullscreen = cli.GetAsType<bool>("fullscreen");
    no_sensing = cli.GetAsType<bool>("nosensing");

    // Change SynChronoManager settings
    syn_manager.SetHeartbeat(heartbeat);

    // Copyright
    LogCopyright(node_id == leader);

    // Sanity check
    assert(leader < num_nodes);

    // --------------
    // Create systems
    // --------------

    // Adjust position of each vehicle so they aren't on top of each other

    // Get the vehicle JSON filenames
    double cam_distance;
    std::string vehicle_filename, powertrain_filename, tire_filename, zombie_filename;
    ChVector<> lidar_pos;
    GetVehicleModelFiles(demo_config[node_id].vehicle_type, vehicle_filename, powertrain_filename, tire_filename,
                         zombie_filename, lidar_pos, cam_distance);

    // Create the vehicle, set parameters, and initialize
    WheeledVehicle vehicle(vehicle_filename, contact_method);
    if (node_id < demo_config.size()) {
        vehicle.Initialize(ChCoordsys<>(demo_config[node_id].pos, demo_config[node_id].rot));
    } else {
        vehicle.Initialize(ChCoordsys<>(demo_config[0].pos, demo_config[0].rot));
    }

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

    std::shared_ptr<ChLidarSensor> lidar;
    std::shared_ptr<ChCameraSensor> camera;
    std::shared_ptr<ChSensorManager> manager;
    if (node_id == leader || !no_sensing) {
        // add a sensor manager
        manager = chrono_types::make_shared<ChSensorManager>(vehicle.GetSystem());
        manager->SetRayRecursions(11);
        Background b;
        b.mode = BackgroundMode::ENVIRONMENT_MAP;  // GRADIENT
        b.color_zenith = {.5f, .6f, .7f};
        b.color_horizon = {.9f, .8f, .7f};
        b.env_tex = GetChronoDataFile("/Environments/sky_2_4k.hdr");
        manager->scene->SetBackground(b);
        float brightness = 1.5f;
        manager->scene->AddPointLight({0, 0, 10000}, {brightness, brightness, brightness}, 100000);

        const int image_width = use_fullscreen ? FS_WIDTH : 1280;
        const int image_height = use_fullscreen ? FS_HEIGHT : 720;
        if (node_id == leader) {
            camera = chrono_types::make_shared<ChCameraSensor>(
                vehicle.GetChassisBody(),  // body camera is attached to
                30.f,                      // update rate in Hz
                chrono::ChFrame<double>({-cam_distance, 0, .45 * cam_distance},
                                        Q_from_AngAxis(0, {0, 1, 0})),  // offset pose
                image_width,                                            // image width
                image_height,                                           // image height
                3.14 / 4,                                               // fov
                1);

            camera->PushFilter(chrono_types::make_shared<ChFilterFullScreenVisualize>(
                image_width, image_height, "Camera 1, Super Sampled", use_fullscreen));
            if (save)
                camera->PushFilter(chrono_types::make_shared<ChFilterSave>("DEMO_OUTPUT/cam/"));
            manager->AddSensor(camera);

            // auto camera2 = chrono_types::make_shared<ChCameraSensor>(
            //     patch->GetGroundBody(),  // body camera is attached to
            //     30.f,                    // update rate in Hz
            //     chrono::ChFrame<double>({830, -40.87, 200.0}, Q_from_AngAxis(3.14 / 2, {0, 1, 0})),  // offset
            //     pose image_width,                                                                         //
            //     image width image_height, // image height 3.14 / 2, // fov 1);

            // camera2->PushFilter(
            //     chrono_types::make_shared<ChFilterFullScreenVisualize>(image_width, image_height, "Camera 2",
            //     false));
            // if (save)
            //     camera2->PushFilter(chrono_types::make_shared<ChFilterSave>("DEMO_OUTPUT/cam2/"));
            // manager->AddSensor(camera2);
        }

        lidar = chrono_types::make_shared<ChLidarSensor>(
            vehicle.GetChassisBody(),                                          // body lidar is attached to
            20.f,                                                              // scanning rate in Hz
            chrono::ChFrame<double>(lidar_pos, Q_from_AngAxis(0, {0, 1, 0})),  // offset pose
            900,                                                               // number of horizontal samples
            16,                                                                // number of vertical channels
            6.28318530718,                                                     // horizontal field of view
            0.261799,
            -0.261799,                         // vertical field of view
            100.f,                             // max distance
            LidarBeamShape::RECTANGULAR,       // beam shape
            2,                                 // sample radius
            0.003,                             // vertical divergence angle
            0.003,                             // horizontal divergence angle
            LidarReturnMode::STRONGEST_RETURN  // return mode for the lidar
        );
        lidar->SetName("Lidar Sensor 1");
        lidar->SetLag(0.01);
        lidar->SetCollectionWindow(.05);
        lidar->PushFilter(chrono_types::make_shared<ChFilterPCfromDepth>());
        lidar->PushFilter(chrono_types::make_shared<ChFilterLidarNoiseXYZI>(0.01f, 0.001f, 0.001f, 0.01f));
        lidar->PushFilter(chrono_types::make_shared<ChFilterVisualizePointCloud>(640, 480, 2, "Lidar Point Cloud"));
        lidar->PushFilter(chrono_types::make_shared<ChFilterXYZIAccess>());
        if (save)
            lidar->PushFilter(chrono_types::make_shared<ChFilterSavePtCloud>("DEMO_OUTPUT/lidar/"));
        manager->AddSensor(lidar);
    }

    // Create the vehicle Irrlicht interface
    IrrAppWrapper app;
    std::shared_ptr<ChDriver> driver;
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
        double steering_time = 1.0;  // time to go from 0 to +1 (or from 0 to -1)
        double throttle_time = 1.0;  // time to go from 0 to +1
        double braking_time = 0.3;   // time to go from 0 to +1
        irr_driver->SetSteeringDelta(render_step_size / steering_time);
        irr_driver->SetThrottleDelta(render_step_size / throttle_time);
        irr_driver->SetBrakingDelta(render_step_size / braking_time);
        irr_driver->Initialize();

        app.Set(temp_app);
        // driver = irr_driver;
    } else if (cli.HasValueInVector<int>("console", node_id)) {
        // Use custom CSL driver instead of irr driver
        auto csl_driver = chrono_types::make_shared<ChCSLDriver>(vehicle);
        driver = csl_driver;
    } else {
        auto path = ChBezierCurve::read(GetChronoDataFile(demo_config[node_id].path_file));
        double target_speed = 11.2;
        bool isPathClosed = true;
        double following_time = 4.0;
        double following_distance = 10;
        double current_distance = 100;

        if (no_sensing) {
            auto path_driver = chrono_types::make_shared<ChPathFollowerACCDriver>(
                vehicle, path, "NSF", target_speed, following_time, following_distance, current_distance, isPathClosed);
            path_driver->GetSpeedController().SetGains(demo_config[node_id].speed_gain_p, 0.01, 0.0);
            path_driver->GetSteeringController().SetGains(0.5, 0.0, 0.0);
            path_driver->GetSteeringController().SetLookAheadDistance(demo_config[node_id].lookahead);
            path_driver->Initialize();

            // Set the callback so that we can check the state of other vehicles
            auto callback =
                std::bind(&VehicleProcessMessageCallback, std::placeholders::_1, std::ref(vehicle), agent, path_driver);
            agent->SetProcessMessageCallback(callback);

            driver = path_driver;
        } else {
            auto path_driver = chrono_types::make_shared<ChLidarWaypointDriver>(
                vehicle, lidar, path, "NSF", target_speed, following_time, following_distance, current_distance,
                isPathClosed);
            path_driver->SetGains(demo_config[node_id].lookahead, 0.5, 0.0, 0.0, demo_config[node_id].speed_gain_p,
                                  0.01, 0.0);
            path_driver->Initialize();

            driver = path_driver;
        }
    }

    // ---------------
    // Simulation loop
    // ---------------
    // Number of simulation steps between miscellaneous events
    int render_steps = (int)std::ceil(render_step_size / step_size);

    // Initialize simulation frame counters
    int step_number = 0;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    double last_time = 0;

    float orbit_radius = 10.f;
    float orbit_rate = .25;
    double time = 0;
    while (app.IsOk() && syn_manager.IsOk() && time < end_time) {
        time = vehicle.GetSystem()->GetChTime();

        // std::cout << "t=" << time << std::endl;

        // End simulation
        if (time >= end_time)
            break;

        // Render scene
        if (step_number % render_steps == 0)
            app.Render();

        // Get driver inputs
        ChDriver::Inputs driver_inputs = driver->GetInputs();

        // Update modules (process inputs from other modules)
        syn_manager.Synchronize(time);  // Synchronize between nodes
        driver->Synchronize(time);
        vehicle.Synchronize(time, driver_inputs, terrain);
        terrain.Synchronize(time);
        app.Synchronize("", driver_inputs);

        // Advance simulation for one timestep for all modules
        driver->Advance(step_size);
        vehicle.Advance(step_size);
        terrain.Advance(step_size);
        app.Advance(step_size);

        if (manager) {
            // if (camera) {
            //     camera->SetOffsetPose(chrono::ChFrame<double>(
            //         {-orbit_radius * cos(time * orbit_rate), -orbit_radius * sin(time * orbit_rate), 1},
            //         Q_from_AngAxis(time * orbit_rate, {0, 0, 1})));
            // }
            manager->Update();
        }

        // Increment frame number
        step_number++;

        // Log clock time
        if (step_number % 500 == 0 && node_id == leader) {
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            // auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()*1000;
            std::chrono::duration<double> wall_time =
                std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

            SynLog() << (wall_time.count()) / (time - last_time) << "\n";
            last_time = time;
            start = std::chrono::high_resolution_clock::now();
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
    cli.AddOption<unsigned int>("Simulation", "l,leader", "The leader rank/node", "0");

    // Irrlicht options
    cli.AddOption<std::vector<int>>("Irrlicht", "i,irr", "Nodes for irrlicht usage", "-1");

    // Controller Options
    cli.AddOption<std::vector<int>>("Keyboard", "k,keyboard", "Force irrlicht driver into keyboard control", "-1");
    cli.AddOption<std::vector<int>>("Simulation", "console", "Nodes to drive with the steering wheel and pedals", "-1");

    // Sensor options
    cli.AddOption<bool>("Simulation", "fullscreen", "Use full screen camera display", std::to_string(use_fullscreen));
    cli.AddOption<bool>("Simulation", "save", "save", std::to_string(save));

    // disable sensing for additional vehicles (not rank 0)
    cli.AddOption<bool>("Simulation", "nosensing", "Disable sensing on non-human vehicles", std::to_string(no_sensing));

// SynChrono/DDS options
#ifdef USE_FAST_DDS
    cli.AddOption<bool>("DDS", "dds", "Use DDS as the communication mechanism", "false");
    cli.AddOption<int>("DDS", "d,node_id", "ID for this Node", "1");
    cli.AddOption<int>("DDS", "n,num_nodes", "Number of Nodes", "2");
#endif

    // Other options
    cli.AddOption<int>("Demo", "v,vehicle", "Vehicle Options [0-4]: Sedan, Audi, SUV, Van, Truck, CityBus", "0");
}

void GetVehicleModelFiles(VehicleType type,
                          std::string& vehicle,
                          std::string& powertrain,
                          std::string& tire,
                          std::string& zombie,
                          ChVector<>& lidar_pos,
                          double& cam_distance) {
    switch (type) {
        case VehicleType::SEDAN:
            vehicle = vehicle::GetDataFile("sedan/vehicle/Sedan_Vehicle.json");
            powertrain = vehicle::GetDataFile("sedan/powertrain/Sedan_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("sedan/tire/Sedan_TMeasyTire.json");
            zombie = vehicle::GetDataFile("sedan/Sedan.json");
            lidar_pos = {1.0, 0, 0.25};
            cam_distance = 6.0;
            break;
        case VehicleType::AUDI:
            vehicle = vehicle::GetDataFile("audi/json/audi_Vehicle.json");
            powertrain = vehicle::GetDataFile("audi/json/audi_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("audi/json/audi_TMeasyTire.json");
            zombie = vehicle::GetDataFile("audi/json/audi.json");
            lidar_pos = {2.3, 0, .4};
            cam_distance = 6.0;
            break;
        case VehicleType::TRUCK:
            vehicle = vehicle::GetDataFile("truck/json/truck_Vehicle.json");
            powertrain = vehicle::GetDataFile("truck/json/truck_SimpleCVTPowertrain.json");
            tire = vehicle::GetDataFile("truck/json/truck_TMeasyTire.json");
            zombie = vehicle::GetDataFile("truck/json/truck.json");
            lidar_pos = {1.92, 0, 0.88};
            cam_distance = 14.0;
            break;
        case VehicleType::VAN:
            vehicle = vehicle::GetDataFile("van/json/van_Vehicle.json");
            powertrain = vehicle::GetDataFile("van/json/van_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("van/json/van_TMeasyTire.json");
            zombie = vehicle::GetDataFile("van/json/van.json");
            lidar_pos = {1.1, 0, 0.5};
            cam_distance = 5.0;
            break;
        case VehicleType::SUV:
            vehicle = vehicle::GetDataFile("suv/json/suv_Vehicle.json");
            powertrain = vehicle::GetDataFile("suv/json/suv_ShaftsPowertrain.json");
            tire = vehicle::GetDataFile("suv/json/suv_TMeasyTire.json");
            zombie = vehicle::GetDataFile("suv/json/suv.json");
            lidar_pos = {.95, 0, 0.45};
            cam_distance = 6.0;
            break;
        case VehicleType::CITYBUS:
            vehicle = vehicle::GetDataFile("citybus/vehicle/CityBus_Vehicle.json");
            powertrain = vehicle::GetDataFile("citybus/powertrain/CityBus_SimpleMapPowertrain.json");
            tire = vehicle::GetDataFile("citybus/tire/CityBus_TMeasyTire.json");
            zombie = vehicle::GetDataFile("citybus/CityBus.json");
            lidar_pos = {2.32, 0, 0.5};
            cam_distance = 14.0;
            break;
    }
}

void AddSceneMeshes(ChSystem* chsystem, RigidTerrain* terrain) {
    // load all meshes in input file, using instancing where possible
    std::string base_path = GetChronoDataFile("/Environments/SanFrancisco/components_new/");
    std::string input_file = base_path + "instance_map_03.csv";
    // std::string input_file = base_path + "instance_map_roads_only.csv";

    std::ifstream infile(input_file);
    if (!infile.is_open())
        throw std::runtime_error("Could not open file " + input_file);
    std::string line, col;
    std::vector<std::string> result;

    std::unordered_map<std::string, std::shared_ptr<ChTriangleMeshConnected>> mesh_map;

    auto mesh_body = chrono_types::make_shared<ChBody>();
    mesh_body->SetBodyFixed(true);
    mesh_body->SetCollide(false);
    chsystem->Add(mesh_body);

    int meshes_added = 0;
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
                if (mesh_name.find("EmissionOn") == std::string::npos &&
                    mesh_name.find("Road") != std::string::npos) {  // exlude items with
                                                                    // emission on
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
                    trimesh_shape->Pos = pos;
                    trimesh_shape->Rot = ChMatrix33<>(rot);

                    mesh_body->AddAsset(trimesh_shape);

                    meshes_added++;
                }
            }
        }
        ChVector<> pos = {0, 0, 100};
        std::cout << "Terrain height at <" << pos.x() << "," << pos.y() << "," << pos.z()
                  << ">: " << terrain->GetHeight(pos) << std::endl;
        std::cout << "Total meshes: " << meshes_added << " | Unique meshes: " << mesh_map.size() << std::endl;
    }
}

void VehicleProcessMessageCallback(std::shared_ptr<SynMessage> message,
                                   WheeledVehicle& vehicle,
                                   std::shared_ptr<SynWheeledVehicleAgent> agent,
                                   std::shared_ptr<ChPathFollowerACCDriver> driver) {
    if (auto vehicle_message = std::dynamic_pointer_cast<SynWheeledVehicleStateMessage>(message)) {
        // The IsInsideBox function will determine whether the a passsed point is inside a box defined by a
        // front position, back position and the width of the box. Rotation of the vectors are taken into account.
        // The positions must be in the same reference, i.e. local OR global, not both.

        // First calculate the box
        // We'll do everything in the local frame
        constexpr double width = 5;
        ChVector<> back(0, 0, 0);
        ChVector<> front(0, 10, 0);

        // Get the zombies position relative to this vehicle
        auto zombie_pos = vehicle_message->chassis.GetFrame().GetPos() - vehicle.GetVehicleCOMPos();

        // Check if the zombie is in our way
        if (IsInsideBox(zombie_pos, front, back, width)) {
            // If yes, slow down
            driver->SetCurrentDistance(zombie_pos.Length());
            // std::cout << message->GetSourceKey().GetNodeID() << " is in " << agent->GetKey().GetNodeID() << "'s way."
            //           << std::endl;
        } else {
            // Otherwise, restore the old distance value so the vehicle continues to drive
            constexpr double current_distance = 20;
            driver->SetCurrentDistance(current_distance);
        }
    }
}
