// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban, Asher Elmquist
// =============================================================================
//
// Chrono demonstration of the sensor module
// Attach multiple sensors to a hmmwv full vehicle model
//
// =============================================================================

#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/driver/ChIrrGuiDriver.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/wheeled_vehicle/utils/ChWheeledVehicleIrrApp.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "chrono_sensor/sensors/ChCameraSensor.h"
#include "chrono_sensor/ChSensorManager.h"
#include "chrono_sensor/filters/ChFilterVisualize.h"
#include "chrono_sensor/utils/ChVisualMaterialUtils.h"
#include "chrono_thirdparty/cxxopts/ChCLI.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"
#include "chrono_vehicle/utils/ChUtilsJSON.h"
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"
#include "extras/driver/ChCSLDriver.h"

#ifdef CHRONO_IRRKLANG
#include "extras/ChCSLSoundEngine.h"
#endif

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::vehicle;
using namespace chrono::sensor;
using namespace chrono::synchrono;

// -----------------------------------------------------------------------------
// Vehicle parameters
// -----------------------------------------------------------------------------

// Initial vehicle location and orientation
// ChVector<> initLoc(-788, -195, -1);
// ChVector<> initLoc(3990.1, -1200.3124, .75); //middle-ish of oval highway
ChVector<> initLoc(3982.1, -12837, .75);  // extreme-ish of oval highway
// ChQuaternion<> initRot(1, 0, 0, 0);
ChQuaternion<> initRot = Q_from_AngZ(-CH_C_PI_2);

enum DriverMode { DEFAULT, RECORD, PLAYBACK };
DriverMode driver_mode = DEFAULT;

// Visualization type for vehicle parts (PRIMITIVES, MESH, or NONE)
VisualizationType chassis_vis_type = VisualizationType::MESH;
VisualizationType suspension_vis_type = VisualizationType::PRIMITIVES;
VisualizationType steering_vis_type = VisualizationType::PRIMITIVES;
VisualizationType wheel_vis_type = VisualizationType::MESH;
VisualizationType tire_vis_type = VisualizationType::MESH;

// Collision type for chassis (PRIMITIVES, HULLS, or NONE)
CollisionType chassis_collision_type = CollisionType::NONE;

// Rigid terrain
RigidTerrain::PatchType terrain_model = RigidTerrain::PatchType::BOX;
double terrainHeight = 0;        // terrain height (FLAT terrain only)
double terrainLength = 40000.0;  // size in X direction
double terrainWidth = 40000.0;   // size in Y direction

// Point on chassis tracked by the camera
ChVector<> trackPoint(0.0, 0.0, 1.75);

// Contact method
ChContactMethod contact_method = ChContactMethod::SMC;

// -----------------------------------------------------------------------------
// Sensor parameters
// -----------------------------------------------------------------------------

// camera parameters
float frame_rate = 30.0;
int super_samples = 1;
unsigned int image_width = 1920;   // / 2;
unsigned int image_height = 1080;  // / 2;
unsigned int fullscreen_image_width = 3840;
unsigned int fullscreen_image_height = 720;
// float cam_fov = 1.608f;
float cam_fov = .524;
bool use_fullscreen = false;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Simulation step sizes
double step_size = 2e-3;

// Simulation end time
double t_end = 10000;

// Save sensor data
bool sensor_save = false;

// Visualize sensor data
bool sensor_vis = true;

std::string demo_data_path = std::string(STRINGIFY(HIGHWAY_DATA_DIR));

using namespace std::chrono;

// =============================================================================

// button callback placeholder
void customButtonCallback(){
    std::cout << "I AM CALLING THE CALLBACK FROM THE JOYSTICK BUTTON, SEE?? \n";
}

void AddCommandLineOptions(ChCLI& cli) {
    cli.AddOption<double>("Simulation", "s,step_size", "Step size", std::to_string(step_size));
    cli.AddOption<double>("Simulation", "e,end_time", "End time", std::to_string(t_end));

    // options for human driver
    cli.AddOption<bool>("Simulation", "nojoystick", "Turn off joystick control", "false");
    cli.AddOption<bool>("Simulation", "fullscreen", "Use full screen camera display", std::to_string(use_fullscreen));
    cli.AddOption<bool>("Simulation", "record", "Record human driver inputs to file", "false");
    cli.AddOption<bool>("Simulation", "replay", "Replay human driver inputs from file", "false");
}

int main(int argc, char* argv[]) {
    ChCLI cli(argv[0]);
    AddCommandLineOptions(cli);

    if (!cli.Parse(argc, argv, true))
        return 0;

    // parse command line inputs
    step_size = cli.GetAsType<double>("step_size");
    t_end = cli.GetAsType<double>("end_time");
    use_fullscreen = cli.GetAsType<bool>("fullscreen");
    if (use_fullscreen) {
        image_width = fullscreen_image_width;
        image_height = fullscreen_image_height;
	cam_fov = 1.608f;
    }
    bool disable_joystick = cli.GetAsType<bool>("nojoystick");
    std::cout << "disable_joystick=" << disable_joystick << std::endl;
    //  = cli.GetAsType<bool>("record");
    //  = cli.GetAsType<bool>("replay");

    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_DATA_DIR + std::string("vehicle/"));

    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    std::string vehicle_file = vehicle::GetDataFile("audi/json/audi_Vehicle.json");
    std::string powertrain_file = vehicle::GetDataFile("audi/json/audi_SimpleMapPowertrain.json");
    std::string tire_file = vehicle::GetDataFile("audi/json/audi_TMeasyTire.json");

    WheeledVehicle vehicle(vehicle_file, ChContactMethod::NSC);
    vehicle.Initialize(ChCoordsys<>(initLoc, initRot));
    vehicle.GetChassis()->SetFixed(false);
    vehicle.SetChassisVisualizationType(chassis_vis_type);
    vehicle.SetSuspensionVisualizationType(suspension_vis_type);
    vehicle.SetSteeringVisualizationType(steering_vis_type);
    vehicle.SetWheelVisualizationType(wheel_vis_type);
    auto powertrain = ReadPowertrainJSON(powertrain_file);
    vehicle.InitializePowertrain(powertrain);
    // Create and initialize the tires
    for (auto& axle : vehicle.GetAxles()) {
        for (auto& wheel : axle->GetWheels()) {
            auto tire = ReadTireJSON(tire_file);
            vehicle.InitializeTire(tire, wheel, tire_vis_type);
        }
    }

    // Create the terrain
    RigidTerrain terrain(vehicle.GetSystem());

    MaterialInfo minfo;
    minfo.mu = 0.9f;
    minfo.cr = 0.01f;
    minfo.Y = 2e7f;
    auto patch_mat = minfo.CreateMaterial(contact_method);

    std::shared_ptr<RigidTerrain::Patch> patch;
    switch (terrain_model) {
        case RigidTerrain::PatchType::BOX:
            patch = terrain.AddPatch(patch_mat, ChVector<>(0, 0, 0), ChVector<>(0, 0, 1), terrainLength, terrainWidth,
                                     2, false, 1, true);
            // patch->SetTexture(vehicle::GetDataFile("terrain/textures/tile4.jpg"), 200, 200);
            break;
        case RigidTerrain::PatchType::HEIGHT_MAP:
            patch = terrain.AddPatch(patch_mat, CSYSNORM, vehicle::GetDataFile("terrain/height_maps/test64.bmp"), 128,
                                     128, 0, 4);
            patch->SetTexture(vehicle::GetDataFile("terrain/textures/grass.jpg"), 16, 16);
            break;
        case RigidTerrain::PatchType::MESH:
            // patch = terrain.AddPatch(patch_mat, CSYSNORM, vehicle::GetDataFile("terrain/meshes/test.obj"));
            patch = terrain.AddPatch(patch_mat, CSYSNORM, demo_data_path + "/Environments/Iowa/road.obj");
            // patch->SetTexture(vehicle::GetDataFile("terrain/textures/grass.jpg"), 100, 100);
            break;
    }

    // auto ground_body = patch->GetGroundBody();
    // auto visual_asset = std::dynamic_pointer_cast<ChVisualization>(ground_body->GetAssets()[0]);
    // auto vis_mat = chrono_types::make_shared<ChVisualMaterial>();
    // vis_mat->SetKdTexture(vehicle::GetDataFile("terrain/textures/grass.jpg"));
    // vis_mat->SetSpecularColor({.0f, .0f, .0f});
    // vis_mat->SetRoughness(1.f);
    // vis_mat->SetUseSpecularWorkflow(false);
    // visual_asset->material_list.push_back(vis_mat);

    terrain.Initialize();

    // // add terrain as visual only
    auto mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    mmesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/oval_highway.obj", false, true);
    mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    trimesh_shape->SetMesh(mmesh);
    trimesh_shape->SetName("highway");

    auto mesh_body = chrono_types::make_shared<ChBody>();
    mesh_body->SetPos({0, 0, .01});
    mesh_body->AddAsset(trimesh_shape);
    mesh_body->SetBodyFixed(true);
    vehicle.GetSystem()->Add(mesh_body);

    // -----------------
    // Initialize output
    // -----------------

    // Set up vehicle output
    vehicle.SetChassisOutput(true);
    vehicle.SetSuspensionOutput(0, true);
    vehicle.SetSteeringOutput(0, true);

    // ------------------------
    // Create the driver system
    // ------------------------

    ChWheeledVehicleIrrApp app(&vehicle, L"Highway Demo");
    ChRealtimeStepTimer realtime_timer;
#ifdef CHRONO_IRRKLANG
    GetLog() << "USING IRRKLANG" << "\n\n";
    ChCSLSoundEngine soundEng(&vehicle);
#endif
    // Create the interactive driver system
    std::shared_ptr<ChDriver> driver;
    if (!disable_joystick) {
        // ChCSLDriver driver(vehicle);
        // driver = chrono_types::make_shared<ChCSLDriver>(vehicle);
        auto IGdriver = chrono_types::make_shared<ChIrrGuiDriver>(app);
        IGdriver->SetButtonCallback(7, &customButtonCallback);
        driver = IGdriver;
    } else {
        double mph_to_ms = 0.44704;
        std::string path_file = demo_data_path + "/Environments/Iowa/oval_highway_path.csv";
        auto path = ChBezierCurve::read(path_file);
        std::string steering_controller_file("hmmwv/SteeringController.json");
        std::string speed_controller_file("hmmwv/SpeedController.json");
        driver = chrono_types::make_shared<ChPathFollowerDriver>(
            vehicle, vehicle::GetDataFile(steering_controller_file), vehicle::GetDataFile(speed_controller_file), path,
            "road", 65 * mph_to_ms, true);
        std::cout << "Using path follower driver\n";
    }

    driver->Initialize();

    // ---------------
    // Simulation loop
    // ---------------

    // output vehicle mass
    std::cout << "VEHICLE MASS: " << vehicle.GetVehicleMass() << std::endl;

    // Initialize simulation frame counter and simulation time
    int step_number = 0;
    int render_frame = 0;
    double time = 0;

    // ---------------------------------------------
    // Create a sensor manager and add a point light
    // ---------------------------------------------
    auto manager = chrono_types::make_shared<ChSensorManager>(vehicle.GetSystem());
    float intensity = 2.0;
    manager->scene->AddPointLight({0, 0, 1e8}, {intensity, intensity, intensity}, 1e12);
    manager->scene->SetAmbientLight({.1, .1, .1});
    manager->scene->SetSceneEpsilon(.01);

    // Set environment map
    Background b;
    b.mode = BackgroundMode::ENVIRONMENT_MAP;
    b.env_tex = GetChronoDataFile("sensor/textures/sunflowers_4k.hdr");
    manager->scene->SetBackground(b);

    // ------------------------------------------------
    // Create a camera and add it to the sensor manager
    // ------------------------------------------------
    auto cam = chrono_types::make_shared<ChCameraSensor>(
        vehicle.GetChassisBody(),                                                    // body camera is attached to
        frame_rate,                                                                  // update rate in Hz
        chrono::ChFrame<double>({0, 0, 500}, Q_from_AngAxis(CH_C_PI_2, {0, 1, 0})),  // offset pose
        1280,                                                                        // image width
        720,                                                                         // image height
        CH_C_PI_4,
        super_samples);  // fov, lag, exposure
    cam->SetName("Camera Sensor");
    if (sensor_vis)
        cam->PushFilter(chrono_types::make_shared<ChFilterVisualize>(1280, 720));

    // add sensor to the manager
    // manager->AddSensor(cam);

    // -------------------------------------------------------
    // Create a second camera and add it to the sensor manager
    // -------------------------------------------------------
    auto cam2 = chrono_types::make_shared<ChCameraSensor>(
        vehicle.GetChassisBody(),                                               // body camera is attached to
        frame_rate,                                                             // update rate in Hz
        chrono::ChFrame<double>({-.2, .4, .95}, Q_from_AngAxis(0, {1, 0, 0})),  // offset pose
        image_width,                                                            // image width
        image_height,                                                           // image height
        cam_fov,
        super_samples);  // fov, lag, exposure
    cam2->SetName("Camera Sensor");

    if (sensor_vis)
        cam2->PushFilter(
            chrono_types::make_shared<ChFilterVisualize>(image_width, image_height, "Driver View", use_fullscreen));

    // add sensor to the manager
    manager->AddSensor(cam2);

    // ---------------
    // Simulate system
    // ---------------
    float orbit_radius = 15.f;
    float orbit_rate = .5;

    auto t0 = high_resolution_clock::now();

    double extra_time = 0.0;
    double last_sim_sync = 0;

    while (app.GetDevice()->run()) {
        time = vehicle.GetSystem()->GetChTime();

        // End simulation
        if (time >= t_end)
            break;

        cam->SetOffsetPose(
            chrono::ChFrame<double>({-orbit_radius * cos(time * orbit_rate), -orbit_radius * sin(time * orbit_rate), 3},
                                    Q_from_AngAxis(time * orbit_rate, {0, 0, 1})));

        // Collect output data from modules (for inter-module communication)
        ChDriver::Inputs driver_inputs = driver->GetInputs();
        // printf("Driver inputs: %f,%f,%f\n", driver_inputs.m_throttle, driver_inputs.m_braking,
        //        driver_inputs.m_steering);
        driver_inputs.m_steering *= -1;
        // driver_inputs.m_throttle = 0;
        if (step_number % int(1 / step_size) == 0) {
            auto speed = vehicle.GetVehicleSpeed();
            auto wall_time = high_resolution_clock::now();
            printf("Sim Time=%f, \tWall Time=%f, \tExtra Time=%f, \tSpeed=%f\n", time,
                   duration_cast<duration<double>>(wall_time - t0).count(), extra_time, speed);
            extra_time = 0.0;
        }

        // Update modules (process inputs from other modules)
        driver->Synchronize(time);
        terrain.Synchronize(time);
        vehicle.Synchronize(time, driver_inputs, terrain);
        app.Synchronize("", driver_inputs);
#ifdef CHRONO_IRRKLANG
        soundEng.Synchronize(time);
#endif
        // Advance simulation for one timestep for all modules
        double step = step_size;
        driver->Advance(step);
        terrain.Advance(step);
        vehicle.Advance(step);
        app.Advance(step_size);

        // Update the sensor manager
        manager->Update();

        // Increment frame number
        step_number++;

        if (step_number % (int)(2.0 / frame_rate / step_size) == 0) {
            double since_last_sync = time - last_sim_sync;
            last_sim_sync = time;
            auto tt0 = high_resolution_clock::now();
            realtime_timer.Spin(since_last_sync);
            auto tt1 = high_resolution_clock::now();
            extra_time += duration_cast<duration<double>>(tt1 - tt0).count();
        }
    }

    auto t1 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t1 - t0);
    std::cout << "Simulation Time: " << t_end << ", Wall Time: " << time_span.count() << std::endl;

    return 0;
}
