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

#include "chrono/utils/ChUtilsInputOutput.h"

#include <irrlicht.h>
#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/driver/ChInteractiveDriverIRR.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemIrrlicht.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "chrono_sensor/ChSensorManager.h"
#include "chrono_sensor/filters/ChFilterVisualize.h"
#include "chrono_sensor/sensors/ChCameraSensor.h"
#include "chrono_thirdparty/cxxopts/ChCLI.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"
#include "chrono_vehicle/utils/ChUtilsJSON.h"
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"
#include "extras/driver/ChCSLDriver.h"

#ifdef CHRONO_IRRKLANG
#include "extras/ChCSLSoundEngine.h"
#endif

using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::sensor;
using namespace chrono::synchrono;

// -----------------------------------------------------------------------------
// Vehicle parameters
// -----------------------------------------------------------------------------

// Initial vehicle location and orientation
// ChVector3d initLoc(-788, -195, -1);
// ChVector3d initLoc(3990.1, -1200.3124, .75); //middle-ish of oval highway
// ChVector3d initLoc(3991.5, 0.0, .75);
ChVector3d initLoc(4011.5, -345, .75);  // near mile marker
// ChQuaternion<> initRot(1, 0, 0, 0);
// ChQuaternion<> initRot = QuatFromAngleZ(-CH_C_PI_2);
ChQuaternion<> initRot = QuatFromAngleZ(CH_C_PI_2);

// ChVector3d driver_eye(-.2, .4, .95);
ChVector3d driver_eye(-.3, .4, .98);
// ChVector3d driver_eye(1.0, .4, .95);
ChQuaternion<> driver_view_direction = QuatFromAngleX(0);

enum DriverMode { HUMAN, AUTONOMOUS };
DriverMode driver_mode = AUTONOMOUS;

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
ChVector3d trackPoint(0.0, 0.0, 1.75);

// Contact method
ChContactMethod contact_method = ChContactMethod::SMC;

// -----------------------------------------------------------------------------
// Sensor parameters
// -----------------------------------------------------------------------------

// camera parameters
float frame_rate = 30.0;
int super_samples = 1;
unsigned int image_width = 3840 / 2;  // 1920;   // / 2;
unsigned int image_height = 720 / 2;  // / 2;
unsigned int fullscreen_image_width = 3840;
unsigned int fullscreen_image_height = 720;
float cam_fov = 1.608f;
// float cam_fov = .524;
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

/// Speedometer image data
irr::core::position2d<irr::s32> sm_center(300, 200);
double sm_needle = 140;

std::string demo_data_path = std::string(STRINGIFY(HIGHWAY_DATA_DIR));

using namespace std::chrono;

// =============================================================================

// button callback placeholder
void customButtonCallback();

void AddCommandLineOptions(ChCLI& cli) {
    cli.AddOption<double>("Simulation", "s,step_size", "Step size", std::to_string(step_size));
    cli.AddOption<double>("Simulation", "e,end_time", "End time", std::to_string(t_end));

    // options for human driver
    cli.AddOption<bool>("Simulation", "nojoystick", "Turn off joystick control", "false");
    cli.AddOption<bool>("Simulation", "lbj", "Switch Joystick axes as used by LBJ", "false");
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
    bool lbj_joystick = cli.GetAsType<bool>("lbj");
    std::cout << "disable_joystick=" << disable_joystick << std::endl;
    //  = cli.GetAsType<bool>("record");
    //  = cli.GetAsType<bool>("replay");

    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_DATA_DIR + std::string("vehicle/"));

    std::cout << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    std::string vehicle_file = vehicle::GetDataFile("audi/json/audi_Vehicle.json");
    std::string engine_file = vehicle::GetDataFile("audi/json/audi_EngineSimpleMap.json");
    std::string transmission_file = vehicle::GetDataFile("audi/json/audi_AutomaticTransmissionSimpleMap.json");
    std::string tire_file = vehicle::GetDataFile("audi/json/audi_TMeasyTire.json");

    WheeledVehicle vehicle(vehicle_file, ChContactMethod::NSC);
    vehicle.Initialize(ChCoordsys<>(initLoc, initRot));
    vehicle.GetChassis()->SetFixed(false);
    vehicle.SetChassisVisualizationType(chassis_vis_type);
    vehicle.SetSuspensionVisualizationType(suspension_vis_type);
    vehicle.SetSteeringVisualizationType(steering_vis_type);
    vehicle.SetWheelVisualizationType(wheel_vis_type);

    auto engine = ReadEngineJSON(engine_file);
    auto transmission = ReadTransmissionJSON(transmission_file);
    auto powertrain = chrono_types::make_shared<ChPowertrainAssembly>(engine, transmission);
    vehicle.InitializePowertrain(powertrain);

    // Create and initialize the tires
    for (auto& axle : vehicle.GetAxles()) {
        for (auto& wheel : axle->GetWheels()) {
            auto tire = ReadTireJSON(tire_file);
            vehicle.InitializeTire(tire, wheel, tire_vis_type);
        }
    }

    // change the ego vehicle vis out for windowless audi

    auto audi_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    audi_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/vehicles/audi_chassis_windowless.obj", false,
                                 true);
    audi_mesh->Transform(ChVector3d(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
    auto audi_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    audi_shape->SetMesh(audi_mesh);
    audi_shape->SetName("Windowless Audi");
    audi_shape->SetMutable(false);
    vehicle.GetChassisBody()->AddVisualShape(audi_shape);

    // add rearview mirror

    auto rvw_mirror_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    rvw_mirror_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/vehicles/audi_rearview_mirror.obj", false,
                                       true);
    rvw_mirror_mesh->Transform(ChVector3d(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
    auto rvw_mirror_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    rvw_mirror_shape->SetMesh(rvw_mirror_mesh);
    rvw_mirror_shape->SetName("Windowless Audi");
    rvw_mirror_shape->SetMutable(false);

    auto mirror_mat = chrono_types::make_shared<ChVisualMaterial>();
    mirror_mat->SetDiffuseColor({0.2f, 0.2f, 0.2f});
    mirror_mat->SetRoughness(0.f);
    mirror_mat->SetMetallic(1.0f);
    mirror_mat->SetUseSpecularWorkflow(false);
    rvw_mirror_shape->AddMaterial(mirror_mat);

    vehicle.GetChassisBody()->AddVisualShape(
        rvw_mirror_shape, ChFrame<>(ChVector3d(0.442, 0.0, 1.096), QuatFromAngleY(-.08) * QuatFromAngleZ(-.25)));

    // Add a leader vehicle
    WheeledVehicle lead_vehicle(vehicle.GetSystem(), vehicle_file);
    lead_vehicle.Initialize(ChCoordsys<>(initLoc + initRot.Rotate(ChVector3d(20, 0, 0)), initRot));
    lead_vehicle.GetChassis()->SetFixed(false);
    lead_vehicle.SetChassisVisualizationType(chassis_vis_type);
    lead_vehicle.SetSuspensionVisualizationType(suspension_vis_type);
    lead_vehicle.SetSteeringVisualizationType(steering_vis_type);
    lead_vehicle.SetWheelVisualizationType(wheel_vis_type);
    lead_vehicle.InitializePowertrain(powertrain);

    // Create and initialize the tires
    for (auto& axle : lead_vehicle.GetAxles()) {
        for (auto& wheel : axle->GetWheels()) {
            auto tire = ReadTireJSON(tire_file);
            lead_vehicle.InitializeTire(tire, wheel, tire_vis_type);
        }
    }

    // Create the terrain
    RigidTerrain terrain(vehicle.GetSystem());

    ChContactMaterialData minfo;
    minfo.mu = 0.9f;
    minfo.cr = 0.01f;
    minfo.Y = 2e7f;
    auto patch_mat = minfo.CreateMaterial(contact_method);

    std::shared_ptr<RigidTerrain::Patch> patch;
    switch (terrain_model) {
        case RigidTerrain::PatchType::BOX:
            patch = terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector3d(0, 0, 0)), terrainLength, terrainWidth, 2,
                                     false, 1, false);
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

    terrain.Initialize();

    // add terrain with weighted textures
    auto terrain_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    terrain_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/terrain/terrain.obj", false, true);
    terrain_mesh->Transform(ChVector3d(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
    auto terrain_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    terrain_shape->SetMesh(terrain_mesh);
    terrain_shape->SetName("terrain");
    terrain_shape->SetMutable(false);

    auto vis_mat2 = chrono_types::make_shared<ChVisualMaterial>();
    vis_mat2->SetKdTexture(demo_data_path + "/Environments/Iowa/terrain/Gravel/GroundGravel017_COL_1K.jpg");
    vis_mat2->SetRoughnessTexture(demo_data_path + "/Environments/Iowa/terrain/Gravel/GroundGravel017_ROUGH_1K.png");
    vis_mat2->SetNormalMapTexture(demo_data_path + "/Environments/Iowa/terrain/Gravel/GroundGravel017_NRM_1K.jpg");
    vis_mat2->SetSpecularColor({.0f, .0f, .0f});
    vis_mat2->SetRoughness(1.f);
    vis_mat2->SetUseSpecularWorkflow(false);
    terrain_shape->AddMaterial(vis_mat2);

    auto terrain_body = chrono_types::make_shared<ChBody>();
    terrain_body->SetPos({0, 0, -.01});
    terrain_body->AddVisualShape(terrain_shape);
    terrain_body->SetBodyFixed(true);
    vehicle.GetSystem()->Add(terrain_body);

    std::vector<std::string> environment_meshes = {"/Environments/Iowa/signs/mile_markers_inner.obj",
                                                   "/Environments/Iowa/signs/mile_markers_outer.obj",
                                                   "/Environments/Iowa/terrain/oval_highway.obj"};
    std::vector<ChVector3d> offsets = {{0, 0, -128.22}, {0, 0, 0.0}, {0, 0, 0.01}};

    for (int i = 0; i < environment_meshes.size(); i++) {  // auto file_name : environment_meshes) {
        // additional environment assets
        auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        trimesh->LoadWavefrontMesh(demo_data_path + environment_meshes[i], false, true);
        trimesh->Transform(ChVector3d(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
        auto trimesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
        trimesh_shape->SetMesh(trimesh);
        trimesh_shape->SetName(environment_meshes[i]);
        trimesh_shape->SetMutable(false);
        auto mesh_body = chrono_types::make_shared<ChBody>();
        mesh_body->SetPos(offsets[i]);
        mesh_body->AddVisualShape(trimesh_shape);
        mesh_body->SetBodyFixed(true);
        vehicle.GetSystem()->Add(mesh_body);
    }

    // add in corn for testing

    auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    trimesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/trees/tree_01.obj", false, true);
    trimesh->Transform(ChVector3d(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
    auto trimesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    trimesh_shape->SetMesh(trimesh);
    trimesh_shape->SetName("Tree");
    trimesh_shape->SetMutable(false);
    trimesh_shape->SetScale({1, 1, 1});

    double x_step = 100;
    double y_step = 200;

    double x_start = initLoc.x() - 60;
    double y_start = initLoc.y() - 7500;

    int x_count = 2;
    int y_count = 75;

    for (int i = 0; i < x_count; i++) {
        for (int j = 0; j < y_count; j++) {
            auto mesh_body = chrono_types::make_shared<ChBody>();
            mesh_body->SetPos({i * x_step + x_start, j * y_step + y_start + .2 * (i % 3), 0.0});
            // mesh_body->SetRot(QuatFromAngleZ(CH_C_PI_2));
            mesh_body->AddVisualShape(trimesh_shape);
            mesh_body->SetBodyFixed(true);
            vehicle.GetSystem()->Add(mesh_body);
        }
    }

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

    ChWheeledVehicleVisualSystemIrrlicht app;
    app.SetWindowTitle("Highway Demo");
    app.SetWindowSize(600, 400);
    app.AttachVehicle(&vehicle);

    ChRealtimeStepTimer realtime_timer;
    /*
    SPEEDOMETER: we want to use the irrlicht app to display the speedometer, but calling endscene would update the
    entire (massive) scenario. In order to do so, we first have to clean app the Irrlichr app. Once we delete the node,
    we remove all cached meshes and textures. The order is important, otherwise meshes are re-cached!!
    */

    //// TODO
    ////
    ////irr::scene::ISceneNode* mnode = app.GetContainer();
    ////mnode->remove();
    ////irr::scene::IMeshCache* cache = app.GetDevice()->getSceneManager()->getMeshCache();
    ////cache->clear();
    ////app.GetVideoDriver()->removeAllTextures();

    app.Initialize();

#ifdef CHRONO_IRRKLANG
    std::cout << "USING IRRKLANG" << std::endl;
    ChCSLSoundEngine soundEng(&vehicle);
#endif

    // Create the interactive driver system
    // ChCSLDriver driver(vehicle);
    // driver = chrono_types::make_shared<ChCSLDriver>(vehicle);
    auto IGdriver = chrono_types::make_shared<ChInteractiveDriverIRR>(app);
    IGdriver->SetButtonCallback(6, &customButtonCallback);
    //// TODO
    ////
    ////IGdriver->SetJoystickAxes(ChIrrGuiDriver::JoystickAxes::AXIS_Z, ChIrrGuiDriver::JoystickAxes::AXIS_R,
    ////                          ChIrrGuiDriver::JoystickAxes::AXIS_X, ChIrrGuiDriver::JoystickAxes::NONE);

    IGdriver->Initialize();

    double mph_to_ms = 0.44704;
    std::string path_file = demo_data_path + "/Environments/Iowa/terrain/oval_highway_path.csv";
    auto path = ChBezierCurve::read(path_file);
    std::string steering_controller_file("hmmwv/SteeringController.json");
    std::string speed_controller_file("hmmwv/SpeedController.json");
    auto PFdriver = chrono_types::make_shared<ChPathFollowerDriver>(
        vehicle, vehicle::GetDataFile(steering_controller_file), vehicle::GetDataFile(speed_controller_file), path,
        "road", 65 * mph_to_ms);
    PFdriver->Initialize();

    if (!disable_joystick) {
        driver_mode = HUMAN;
    } else {
        driver_mode = AUTONOMOUS;
        std::cout << "Using path follower driver\n";
    }

    // Leader Driver
    auto lead_PFdriver = chrono_types::make_shared<ChPathFollowerDriver>(
        vehicle, vehicle::GetDataFile(steering_controller_file), vehicle::GetDataFile(speed_controller_file), path,
        "road", 65 * mph_to_ms);
    lead_PFdriver->Initialize();

    // ---------------
    // Simulation loop
    // ---------------

    // output vehicle mass
    std::cout << "VEHICLE MASS: " << vehicle.GetMass() << std::endl;

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
        lead_vehicle.GetChassisBody(),                                   // body camera is attached to
        frame_rate,                                                      // update rate in Hz
        chrono::ChFrame<double>({0, 0, 50}, QuatFromAngleY(CH_C_PI_2)),  // offset pose
        1280,                                                            // image width
        720,                                                             // image height
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
        vehicle.GetChassisBody(),                                    // body camera is attached to
        frame_rate,                                                  // update rate in Hz
        chrono::ChFrame<double>(driver_eye, driver_view_direction),  // offset pose
        image_width,                                                 // image width
        image_height,                                                // image height
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
    float orbit_radius = 1000.f;
    float orbit_rate = .5;

    auto t0 = high_resolution_clock::now();

    double extra_time = 0.0;
    double last_sim_sync = 0;

    while (app.GetDevice()->run()) {
        time = vehicle.GetSystem()->GetChTime();

        // End simulation
        if (time >= t_end)
            break;

        cam->SetOffsetPose(chrono::ChFrame<double>(
            {-orbit_radius * cos(time * orbit_rate), -orbit_radius * sin(time * orbit_rate), orbit_radius / 5.0},
            QuatFromAngleZ(time * orbit_rate)));

        // Collect output data from modules (for inter-module communication)
        DriverInputs driver_inputs;
        if (driver_mode == AUTONOMOUS)
            driver_inputs = PFdriver->GetInputs();
        else
            driver_inputs = IGdriver->GetInputs();

        DriverInputs lead_driver_inputs = lead_PFdriver->GetInputs();
        // printf("Driver inputs: %f,%f,%f\n", driver_inputs.m_throttle, driver_inputs.m_braking,
        //        driver_inputs.m_steering);
        driver_inputs.m_steering *= -1;
        // driver_inputs.m_throttle = 0;
        if (step_number % int(1 / step_size) == 0) {
            auto speed = vehicle.GetSpeed();
            auto wall_time = high_resolution_clock::now();
            printf("Sim Time=%f, \tWall Time=%f, \tExtra Time=%f, \tSpeed=%f\n", time,
                   duration_cast<duration<double>>(wall_time - t0).count(), extra_time, speed);
            extra_time = 0.0;
        }

        // Update modules (process inputs from other modules)
        if (driver_mode == AUTONOMOUS)
            PFdriver->Synchronize(time);
        else
            IGdriver->Synchronize(time);
        lead_PFdriver->Synchronize(time);
        terrain.Synchronize(time);
        vehicle.Synchronize(time, driver_inputs, terrain);
        lead_vehicle.Synchronize(time, lead_driver_inputs, terrain);
        app.Synchronize(time, driver_inputs);
#ifdef CHRONO_IRRKLANG
        soundEng.Synchronize(time);
#endif
        // Advance simulation for one timestep for all modules
        double step = step_size;

        if (driver_mode == AUTONOMOUS)
            PFdriver->Advance(step);
        else
            IGdriver->Advance(step);
        lead_PFdriver->Advance(step);
        terrain.Advance(step);
        lead_vehicle.Advance(step);
        vehicle.Advance(step);
        app.Advance(step_size);
        if (step_number % int(1 / (60 * step_size)) == 0) {
            /// irrlicht::tools::drawSegment(app.GetVideoDriver(), v1, v2, video::SColor(255, 80, 0, 0), false);
            app.GetDevice()->getVideoDriver()->draw2DImage(
                app.GetDevice()->getVideoDriver()->getTexture(
                    (demo_data_path + "/miscellaneous/Speedometer.png").c_str()),
                irr::core::position2d<irr::s32>(0, 0));
            /*app.GetDevice()->getVideoDriver()->draw2DImage(app.GetDevice()->getVideoDriver()->getTexture((demo_data_path
               + "/miscellaneous/Needle.png").c_str()), irr::core::position2d<irr::s32>(200, 200));*/
            double speed_mph = vehicle.GetSpeed() * 2.23694;
            double theta = ((270 / 140) * speed_mph) * (CH_C_PI / 180);
            app.GetDevice()->getVideoDriver()->draw2DLine(
                sm_center + irr::core::position2d<irr::s32>(-sm_needle * sin(theta), sm_needle * cos(theta)), sm_center,
                irr::video::SColor(255, 255, 0, 0));
            app.GetDevice()->getVideoDriver()->endScene();
        }

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

void customButtonCallback() {
    std::cout << "Button Callback Invoked";
    driver_mode = HUMAN;
}
