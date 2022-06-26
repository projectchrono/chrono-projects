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
// Authors: Radu Serban, Asher Elmquist, Simone Benatti, Jason Zhou
// =============================================================================
//
// Iowa Highway Simulation
// This demo includes IG interactive vehicle, autonomous dynamic vehicles and
// autonomous dummy vehicles
//
// =============================================================================

#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include <irrlicht.h>
#include <time.h>
#include <limits>
#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/driver/ChIrrGuiDriver.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/wheeled_vehicle/utils/ChWheeledVehicleIrrApp.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "chrono_sensor/ChSensorManager.h"
#include "chrono_sensor/filters/ChFilterVisualize.h"
#include "chrono_sensor/sensors/ChCameraSensor.h"
#include "chrono_sensor/utils/ChVisualMaterialUtils.h"
#include "chrono_thirdparty/cxxopts/ChCLI.h"
//#include "chrono_vehicle/driver/ChPathFollowerDriver.h"

#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"
#include "extras/driver/ChCSLDriver.h"
#include "extras/driver/ChNSF_Drivers.h"

#include "chrono/utils/ChFilters.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef CHRONO_IRRKLANG
#include "extras/ChCSLSoundEngine.h"
#endif

#include "chrono_vehicle/utils/ChUtilsJSON.h"

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::vehicle;
using namespace chrono::sensor;
using namespace chrono::synchrono;

// -----------------------------------------------------------------------------
// rad to RPM conversion parameters
// -----------------------------------------------------------------------------
const double rads2rpm = 30 / CH_C_PI;

// -----------------------------------------------------------------------------
// Vehicle parameters
// -----------------------------------------------------------------------------

// Initial vehicle location and orientation
ChVector<> initLoc(5011.5, -445, 100.75);  // near mile marker

ChQuaternion<> initRot = Q_from_AngZ(-CH_C_PI_2);

ChVector<> driver_eye(-.3, .4, .98);

ChQuaternion<> driver_view_direction = Q_from_AngAxis(0, {1, 0, 0});

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
ChVector<> trackPoint(0.0, 0.0, 1.75);

// Contact method
ChContactMethod contact_method = ChContactMethod::SMC;

// Global Variable to track miles the car traveled
utils::ChRunningAverage IG_speed_avg(100);

// -----------------------------------------------------------------------------
// Steering wheel parameters
// -----------------------------------------------------------------------------
// steering wheel button settings
// change mapping if on a different steering wheel
int r_3 = 19;
int r_1 = 6;

// camera parameters
float frame_rate = 30.0;
int super_samples = 2;
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
irr::core::position2d<irr::s32> sm_center(219, 215);
irr::core::position2d<irr::s32> rpm_center(609, 215);
irr::core::position2d<irr::s32> gr_center(388, 23);
irr::core::position2d<irr::s32> auto_center(358, 338);
irr::core::position2d<irr::s32> cur_left(950, 70);
irr::core::position2d<irr::s32> end_left(950, 190);
irr::core::position2d<irr::s32> eta_left(950, 310);
double sm_needle = 140;

std::string demo_data_path = std::string(STRINGIFY(HIGHWAY_DATA_DIR));

// Driver parameters
std::vector<double> followerParam;
std::string scenario_parameters = "scenario_parameters.json";
std::string simulation_parameters = "simulation_parameters.json";
std::string lead_parameters = demo_data_path + "lead_parameters_1.json";
// cruise speed [mph]
double cruise_speed = 45;

int correct_ratio = 200;

// Data saving
bool save_driver = true;
// time interval between data savings
double tsave = 1e-2;
// path where the output is saved
std::string output_file_path = "./output.csv";
std::string dummy_button_path = "./buttoninfo.csv";
std::stringstream buffer;
std::stringstream button_buffer;
std::ofstream filestream;
std::ofstream buttonstream;

// Fog parameters
bool fog_enabled = false;
ChVector<float> fog_color = {.8, .8, .8};
float fog_distance = 2000.0;

// mirrors position and rotations
ChVector<> mirror_rearview_pos = {0.0, 0.0, 0.0};
ChQuaternion<> mirror_rearview_rot = {1.0, 0.0, 0.0, 0.0};
ChVector<> mirror_wingleft_pos = {0.0, 0.0, 0.0};
ChQuaternion<> mirror_wingleft_rot = {1.0, 0.0, 0.0, 0.0};
ChVector<> mirror_wingright_pos = {0.0, 0.0, 0.0};
ChQuaternion<> mirror_wingright_rot = {1.0, 0.0, 0.0, 0.0};

ChVector<> arrived_sign_pos = {0.0, 0.0, 0.0};
ChQuaternion<> arrived_sign_rot = {1.0, 0.0, 0.0, 0.0};

// benchmarking
bool benchmark = false;

// Dummy vehicle offset
float dummy_audi_z_offset = 0.25;
float dummy_patrol_z_offset = 0.5;
using namespace std::chrono;

// Define all lead vehicles
int num_dummy = 0;    // number of dummy vehicles
int num_dynamic = 0;  // number of dynamic vehicles
int lead_count = 0;
// Dummies
std::vector<ChVector<>> dummy_pos;
std::vector<float> dummy_cruise_speed;
std::vector<int> dummy_lane;  // 0 for inner, 1 for outer
// 0 for no control(constant speed), 1 for distance control, 2 for time control
std::vector<int> dummy_control;
// time control variables
std::vector<std::vector<float>> dummy_control_x;
std::vector<std::vector<float>> dummy_control_y;
std::vector<int> dummy_time_mode;  // 0 if dist mode or this field is N/A, 1 means using sim time, 2 means using wall
                                   // time (real time)
std::vector<double>
    dummy_dist;  // since dummy vehicle do not have a tracker, we need to track the distance of the dummy vehicle
std::vector<ChVector<>> dummy_prev_pos;

// Dynamic Leaders
std::vector<ChVector<>> dynamic_pos;
std::vector<float> dynamic_cruise_speed;
std::vector<int> dynamic_lane;  // 0 for inner, 1 for outer
// 0 for no control(constant speed), 1 for distance control, 2 for time control
std::vector<int> dynamic_control;
// time control variables
std::vector<std::vector<float>> dynamic_control_x;
std::vector<std::vector<float>> dynamic_control_y;
std::vector<int> dynamic_time_mode;  // 0 if dist mode or this field is N/A, 1 means using sim time, 2 means using wall
                                     // time (real time)

std::vector<std::vector<std::vector<double>>> leaderParam;

// Comment section in csv output
std::string csv_comments;
std::string filename;

// distance variable
float IG_dist = 0;
ChVector<> IG_prev_pos;

// driver global variables
// TODO: maybe there is a better way to handle this
std::shared_ptr<ChNSFFollowererDriver> PF_driver_ptr;
float cur_follower_speed;

// Experiment parameters
int meet_time = 4;      // this meet time is defined in full minutes
double eta_dist = 3.6;  // eta counter distance

// Start Stamp
float start_sim_time;
float start_wall_time;
// =============================================================================

// button callback placeholder
void customButtonCallback();
void dummyButtonCallback_r_3();

// distribution of trees near the road
void AddTrees(ChSystem* chsystem);
// terrain and texture system
void AddTerrain(ChSystem* chsystem);
// signs along the road and the road itself
void AddRoadway(ChSystem* chsystem);
// buildings in the environment
void AddBuildings(ChSystem* chsystem);

// helper function to update dummy vehicles
void updateDummy(std::shared_ptr<ChBodyAuxRef> dummy_vehicle,
                 std::shared_ptr<ChBezierCurve> curve,
                 float dummy_speed,
                 float step_size,
                 float z_offset,
                 ChBezierCurveTracker tracker,
                 double& dummy_dist,
                 ChVector<>& dummy_prev_pos);

// compute desired speed for dynamic vehicles
// based on json defined, time-dependent, piecewise speed data
float controlFindSpeed_x_y(std::vector<float> time_vec, std::vector<float> speed_vec, float time, float default_speed);

void ReadParameterFiles() {
    {  // Simulation parameter file
        rapidjson::Document d;
        vehicle::ReadFileJSON(simulation_parameters, d);

        if (d.HasMember("TimeStep")) {
            step_size = d["TimeStep"].GetDouble();
        }

        if (d.HasMember("Camera")) {
            const rapidjson::Value& camera_params = d["Camera"];
            if (camera_params.HasMember("DriverEye")) {
                driver_eye = vehicle::ReadVectorJSON(camera_params["DriverEye"]);
            }
            if (camera_params.HasMember("FrameRate")) {
                frame_rate = camera_params["FrameRate"].GetFloat();
            }
            if (camera_params.HasMember("SuperSamples")) {
                super_samples = camera_params["SuperSamples"].GetInt();
            }
            if (camera_params.HasMember("FieldOfView")) {
                cam_fov = camera_params["FieldOfView"].GetFloat() * CH_C_DEG_TO_RAD;
            }
        }
        if (d.HasMember("Fog")) {
            const rapidjson::Value& fog_params = d["Fog"];
            if (fog_params.HasMember("Enabled")) {
                fog_enabled = fog_params["Enabled"].GetBool();
            }
            if (fog_params.HasMember("Color")) {
                fog_color = vehicle::ReadVectorJSON(fog_params["Color"]);
            }
            if (fog_params.HasMember("Distance")) {
                fog_distance = fog_params["Distance"].GetFloat();
            }
        }

        if (d.HasMember("Mirrors")) {
            const rapidjson::Value& mirror_params = d["Mirrors"];
            if (mirror_params.HasMember("Rearview")) {
                const rapidjson::Value& rearview_params = mirror_params["Rearview"];
                if (rearview_params.HasMember("Position")) {
                    mirror_rearview_pos = vehicle::ReadVectorJSON(rearview_params["Position"]);
                }
                if (rearview_params.HasMember("Rotation")) {
                    mirror_rearview_rot =
                        Q_from_Euler123(CH_C_DEG_TO_RAD * vehicle::ReadVectorJSON(rearview_params["Rotation"]));
                }
            }
            if (mirror_params.HasMember("WingLeft")) {
                const rapidjson::Value& wingleft_params = mirror_params["WingLeft"];
                if (wingleft_params.HasMember("Position")) {
                    mirror_wingleft_pos = vehicle::ReadVectorJSON(wingleft_params["Position"]);
                }
                if (wingleft_params.HasMember("Rotation")) {
                    mirror_wingleft_rot =
                        Q_from_Euler123(CH_C_DEG_TO_RAD * vehicle::ReadVectorJSON(wingleft_params["Rotation"]));
                }
            }
            if (mirror_params.HasMember("WingRight")) {
                const rapidjson::Value& wingright_params = mirror_params["WingRight"];
                if (wingright_params.HasMember("Position")) {
                    mirror_wingright_pos = vehicle::ReadVectorJSON(wingright_params["Position"]);
                }
                if (wingright_params.HasMember("Rotation")) {
                    mirror_wingright_rot =
                        Q_from_Euler123(CH_C_DEG_TO_RAD * vehicle::ReadVectorJSON(wingright_params["Rotation"]));
                }
            }
        }
    }

    {  // Scenario parameter file
        rapidjson::Document d;
        vehicle::ReadFileJSON(scenario_parameters, d);

        if (d.HasMember("StartLocation")) {
            initLoc = vehicle::ReadVectorJSON(d["StartLocation"]);
        }
        if (d.HasMember("EndTime")) {
            t_end = d["EndTime"].GetDouble();
        }

        if (d.HasMember("MeetTime")) {
            meet_time = d["MeetTime"].GetInt();
        }

        if (d.HasMember("ETADist")) {
            eta_dist = d["ETADist"].GetDouble();
        }
        // TODO: figure out what is happening there, not sure necessary
        if (d.HasMember("FollowerDriverParam")) {
            auto marr = d["FollowerDriverParam"].GetArray();
            int msize = marr.Size();
            assert(msize == 6);
            followerParam.resize(msize);
            for (int j = 0; j < marr.Size(); j++) {
                followerParam[j] = marr[j].GetDouble();
            }
        } else {
            followerParam.resize(6);
            followerParam = {30, 1.5, 2.0, 5.0, 3.0, 4.0};
        }

        if (d.HasMember("CruiseSpeed")) {
            cruise_speed = d["CruiseSpeed"].GetDouble();
        }
        if (d.HasMember("ArrivedSign")) {
            const rapidjson::Value& arrived_sign_params = d["ArrivedSign"];
            if (arrived_sign_params.HasMember("Position")) {
                arrived_sign_pos = vehicle::ReadVectorJSON(arrived_sign_params["Position"]);
            }
            if (arrived_sign_params.HasMember("Rotation")) {
                arrived_sign_rot =
                    Q_from_Euler123(CH_C_DEG_TO_RAD * vehicle::ReadVectorJSON(arrived_sign_params["Rotation"]));
            }
        }
    }

    {
        // Leading Vehicles Scenario parameter file
        // Scenario parameter file
        rapidjson::Document d;
        vehicle::ReadFileJSON(lead_parameters, d);

        while (true) {
            std::string entry_name = "lead_" + std::to_string(lead_count);

            std::cout << entry_name << std::endl;

            if (d.HasMember(entry_name.c_str())) {
                std::cout << "lead_count:" << lead_count << std::endl;
                if (d[entry_name.c_str()]["is_dummy"].GetInt() == 0) {
                    dynamic_pos.push_back(vehicle::ReadVectorJSON(d[entry_name.c_str()]["initial_pos"]));
                    dynamic_cruise_speed.push_back(d[entry_name.c_str()]["cruise_speed"].GetDouble());
                    dynamic_lane.push_back(d[entry_name.c_str()]["lane"].GetInt());
                    num_dynamic++;

                    // read speed control
                    if (d[entry_name.c_str()].HasMember("time_speed_control")) {
                        dynamic_control.push_back(2);
                        std::vector<float> temp_time;
                        std::vector<float> temp_speed;

                        auto marr = d[entry_name.c_str()]["time_speed_control"].GetArray();
                        for (int i = 0; i < marr[0].Size(); i++) {
                            for (int j = 0; j < marr[1].Size(); j++) {
                                if (i == 0) {
                                    temp_time.push_back(marr[i][j].GetDouble());
                                } else if (i == 1) {
                                    temp_speed.push_back(marr[i][j].GetDouble());
                                }
                            }
                        }
                        dynamic_control_x.push_back(temp_time);
                        dynamic_control_y.push_back(temp_speed);

                        // update the dummy_time_mode
                        if (d[entry_name.c_str()].HasMember("time_mode")) {
                            dynamic_time_mode.push_back(d[entry_name.c_str()]["time_mode"].GetInt());
                        } else {
                            dynamic_time_mode.push_back(1);
                        }
                    } else if ((d[entry_name.c_str()].HasMember("dist_speed_control"))) {
                        dynamic_control.push_back(1);
                        std::vector<float> temp_dist;
                        std::vector<float> temp_speed;

                        auto marr = d[entry_name.c_str()]["dist_speed_control"].GetArray();
                        for (int i = 0; i < marr[0].Size(); i++) {
                            for (int j = 0; j < marr[1].Size(); j++) {
                                if (i == 0) {
                                    temp_dist.push_back(marr[i][j].GetDouble());
                                } else if (i == 1) {
                                    temp_speed.push_back(marr[i][j].GetDouble());
                                }
                            }
                        }
                        dynamic_control_x.push_back(temp_dist);
                        dynamic_control_y.push_back(temp_speed);

                        // update the dummy_time_mode
                        dynamic_time_mode.push_back(0);
                    } else {
                        std::vector<float> empty_vec_x;
                        std::vector<float> empty_vec_y;
                        dynamic_control_x.push_back(empty_vec_x);
                        dynamic_control_y.push_back(empty_vec_y);
                        dynamic_control.push_back(0);

                        dynamic_time_mode.push_back(0);
                    }

                    if (d[entry_name.c_str()].HasMember("LeaderDriverParam")) {
                        std::vector<std::vector<double>> temp_leaderParam;
                        auto marr = d[entry_name.c_str()]["LeaderDriverParam"].GetArray();
                        int msize0 = marr.Size();
                        int msize1 = marr[0].Size();
                        assert(msize1 == 6);
                        temp_leaderParam.resize(msize0);
                        // printf("ARRAY DIM = %i \n", msize);
                        for (auto it = marr.begin(); it != marr.end(); ++it) {
                            auto i = std::distance(marr.begin(), it);
                            temp_leaderParam[i].resize(msize1);
                            for (int j = 0; j < marr[i].Size(); j++) {
                                temp_leaderParam[i][j] = marr[i][j].GetDouble();
                                // std::cout<< "param :" << i << "," << j << ":" << marr[i][j].GetDouble() << "\n";
                            }
                        }
                        leaderParam.push_back(temp_leaderParam);
                    } else {
                        std::vector<std::vector<double>> temp_leaderParam;
                        temp_leaderParam.resize(1);
                        temp_leaderParam[0] = {0.5, 1.5, 55.0, 5.0, 628.3, 0.0};
                        leaderParam.push_back(temp_leaderParam);
                    }

                } else {
                    dummy_pos.push_back(vehicle::ReadVectorJSON(d[entry_name.c_str()]["initial_pos"]));
                    dummy_cruise_speed.push_back(d[entry_name.c_str()]["cruise_speed"].GetDouble());
                    dummy_lane.push_back(d[entry_name.c_str()]["lane"].GetInt());
                    num_dummy++;

                    // read speed control
                    if (d[entry_name.c_str()].HasMember("time_speed_control")) {
                        dummy_control.push_back(2);
                        std::vector<float> temp_time;
                        std::vector<float> temp_speed;

                        auto marr = d[entry_name.c_str()]["time_speed_control"].GetArray();
                        for (int i = 0; i < marr[0].Size(); i++) {
                            for (int j = 0; j < marr[1].Size(); j++) {
                                if (i == 0) {
                                    temp_time.push_back(marr[i][j].GetDouble());
                                } else if (i == 1) {
                                    temp_speed.push_back(marr[i][j].GetDouble());
                                }
                            }
                        }
                        dummy_control_x.push_back(temp_time);
                        dummy_control_y.push_back(temp_speed);
                        dummy_dist.push_back(0.0);
                        dummy_prev_pos.push_back(dummy_pos[dummy_pos.size() - 1]);

                        // update the dummy_time_mode
                        if (d[entry_name.c_str()].HasMember("time_mode")) {
                            dummy_time_mode.push_back(d[entry_name.c_str()]["time_mode"].GetInt());
                        } else {
                            dummy_time_mode.push_back(1);
                        }
                    } else if ((d[entry_name.c_str()].HasMember("dist_speed_control"))) {
                        dummy_control.push_back(1);
                        std::vector<float> temp_dist;
                        std::vector<float> temp_speed;

                        auto marr = d[entry_name.c_str()]["dist_speed_control"].GetArray();
                        for (int i = 0; i < marr[0].Size(); i++) {
                            for (int j = 0; j < marr[1].Size(); j++) {
                                if (i == 0) {
                                    temp_dist.push_back(marr[i][j].GetDouble());
                                } else if (i == 1) {
                                    temp_speed.push_back(marr[i][j].GetDouble());
                                }
                            }
                        }
                        dummy_control_x.push_back(temp_dist);
                        dummy_control_y.push_back(temp_speed);
                        dummy_dist.push_back(0.0);
                        dummy_prev_pos.push_back(dummy_pos[dummy_pos.size() - 1]);

                        dummy_time_mode.push_back(0);

                    } else {
                        std::vector<float> empty_vec_x;
                        std::vector<float> empty_vec_y;
                        dummy_control_x.push_back(empty_vec_x);
                        dummy_control_y.push_back(empty_vec_y);
                        dummy_control.push_back(0);
                        dummy_dist.push_back(0.0);
                        dummy_prev_pos.push_back(dummy_pos[dummy_pos.size() - 1]);

                        dummy_time_mode.push_back(0);
                    }
                }
                lead_count++;
            } else {
                break;
            }
        }
    }
}

void AddCommandLineOptions(ChCLI& cli) {
    cli.AddOption<double>("Simulation", "s,step_size", "Step size", std::to_string(step_size));
    cli.AddOption<double>("Simulation", "e,end_time", "End time", std::to_string(t_end));

    // options for human driver
    cli.AddOption<bool>("Simulation", "nojoystick", "Turn off joystick control", "false");
    cli.AddOption<bool>("Simulation", "lbj", "Switch Joystick axes as used by LBJ", "false");
    cli.AddOption<bool>("Simulation", "fullscreen", "Use full screen camera display", std::to_string(use_fullscreen));
    cli.AddOption<bool>("Simulation", "record", "Record human driver inputs to file", "false");
    cli.AddOption<bool>("Simulation", "replay", "Replay human driver inputs from file", "false");

    cli.AddOption<bool>("Simulation", "birdseye", "Enable birds eye camera", "false");
    cli.AddOption<bool>("Simulation", "benchmark", "Benchmark the simulation", "false");

    cli.AddOption<std::string>("Simulation", "scenario_params", "Path to scenario configuration file",
                               scenario_parameters);
    cli.AddOption<std::string>("Simulation", "sim_params", "Path to simulation configuration file",
                               simulation_parameters);
    cli.AddOption<std::string>("Simulation", "lead_params", "Path to lead configuration file", lead_parameters);

    cli.AddOption<std::string>("Simulation", "filename", "Output Filenames", filename);
    cli.AddOption<std::string>("Simulation", "csv_comments", "CSV output comments", csv_comments);
}

int main(int argc, char* argv[]) {
    // create cli tool
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

    scenario_parameters =
        demo_data_path + "/Environments/Iowa/parameters/" + cli.GetAsType<std::string>("scenario_params");
    simulation_parameters =
        demo_data_path + "/Environments/Iowa/parameters/" + cli.GetAsType<std::string>("sim_params");
    lead_parameters = demo_data_path + "/Environments/Iowa/parameters/" + cli.GetAsType<std::string>("lead_params");

    std::cout << "Scen params: " << scenario_parameters << "\n";
    std::cout << "Sim params:" << simulation_parameters << "\n";
    std::cout << "Lead params:" << lead_parameters << "\n";

    ReadParameterFiles();

    benchmark = cli.GetAsType<bool>("benchmark");
    if (benchmark) {
        disable_joystick = true;
        lbj_joystick = false;
        t_end = cli.GetAsType<double>("end_time");
    }

    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_DATA_DIR + std::string("vehicle/"));

    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    std::string vehicle_file = vehicle::GetDataFile("audi/json/audi_Vehicle.json");
    std::string powertrain_file = vehicle::GetDataFile("audi/json/audi_SimpleMapPowertrain.json");
    std::string tire_file = vehicle::GetDataFile("audi/json/audi_TMeasyTire.json");

    // std::string path_file = demo_data_path + "/Environments/Iowa/terrain/oval_highway_path.csv";
    std::string outer_path_file = demo_data_path + "/Environments/Iowa/Driver/OnOuterLane.txt";
    auto outer_path = ChBezierCurve::read(outer_path_file);
    auto inner_path_file = demo_data_path + "/Environments/Iowa/Driver/OnInnerLane.txt";
    auto inner_path = ChBezierCurve::read(inner_path_file);

    // IG vehicle lane number tracker
    // lane 1 - inner lane; lane 2 - outer lanes
    auto lane_0_path_file = demo_data_path + "/Environments/Iowa/Driver/OnInnerLane.txt";
    auto lane_0_path = ChBezierCurve::read(lane_0_path_file);
    auto lane_1_path_file = demo_data_path + "/Environments/Iowa/Driver/OnOuterLane.txt";
    auto lane_1_path = ChBezierCurve::read(lane_1_path_file);

    WheeledVehicle vehicle(vehicle_file, ChContactMethod::SMC);
    auto ego_chassis = vehicle.GetChassis();
    vehicle.Initialize(ChCoordsys<>(initLoc, initRot));
    vehicle.GetChassis()->SetFixed(false);
    vehicle.SetChassisVisualizationType(chassis_vis_type);
    vehicle.SetSuspensionVisualizationType(suspension_vis_type);
    vehicle.SetSteeringVisualizationType(steering_vis_type);
    vehicle.SetWheelVisualizationType(wheel_vis_type);
    auto powertrain = ReadPowertrainJSON(powertrain_file);
    vehicle.InitializePowertrain(powertrain);

    output_file_path = "./output_" + cli.GetAsType<std::string>("filename") + ".csv ";
    dummy_button_path = "./buttoninfo_" + cli.GetAsType<std::string>("filename") + ".csv ";
    filestream = std::ofstream(output_file_path);
    buttonstream = std::ofstream(dummy_button_path);

    // Create and initialize the tires
    for (auto& axle : vehicle.GetAxles()) {
        for (auto& wheel : axle->GetWheels()) {
            auto tire = ReadTireJSON(tire_file);
            vehicle.InitializeTire(tire, wheel, tire_vis_type);
        }
    }

    // change the ego vehicle vis out for windowless audi
    vehicle.GetChassisBody()->GetAssets().clear();

    auto audi_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    audi_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/vehicles/audi_chassis_windowless_2.obj", false,
                                 true);
    audi_mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
    auto audi_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    audi_shape->SetMesh(audi_mesh);
    audi_shape->SetName("Windowless Audi");
    audi_shape->SetStatic(true);
    vehicle.GetChassisBody()->AddAsset(audi_shape);

    // add rearview mirror
    auto mirror_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    mirror_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/vehicles/audi_rearview_mirror.obj", false,
                                   true);
    mirror_mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    auto mirror_mat = chrono_types::make_shared<ChVisualMaterial>();
    mirror_mat->SetDiffuseColor({0.2f, 0.2f, 0.2f});
    mirror_mat->SetRoughness(0.f);
    mirror_mat->SetMetallic(1.0f);
    mirror_mat->SetUseSpecularWorkflow(false);

    auto rvw_mirror_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    rvw_mirror_shape->SetMesh(mirror_mesh);
    rvw_mirror_shape->SetName("Windowless Audi");
    rvw_mirror_shape->SetStatic(true);
    rvw_mirror_shape->SetScale({1, 1.8, 1.2});
    rvw_mirror_shape->Pos = mirror_rearview_pos;
    rvw_mirror_shape->Rot = mirror_rearview_rot;  // Q_from_AngY(-.08) * Q_from_AngZ(-.25);
    rvw_mirror_shape->material_list.push_back(mirror_mat);
    vehicle.GetChassisBody()->AddAsset(rvw_mirror_shape);

    // add left wing mirror
    auto lwm_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    lwm_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/vehicles/audi_left_wing_mirror.obj", false, true);
    lwm_mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    auto lwm_mirror_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    lwm_mirror_shape->SetMesh(lwm_mesh);
    lwm_mirror_shape->SetName("Windowless Audi");
    lwm_mirror_shape->SetStatic(true);
    lwm_mirror_shape->SetScale({1, .95, .95});
    lwm_mirror_shape->Pos = mirror_wingleft_pos;
    lwm_mirror_shape->Rot = mirror_wingleft_rot;  // Q_from_AngY(-.08) * Q_from_AngZ(-.25);
    lwm_mirror_shape->material_list.push_back(mirror_mat);
    vehicle.GetChassisBody()->AddAsset(lwm_mirror_shape);

    // add left wing mirror
    auto rwm_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    rwm_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/vehicles/audi_right_wing_mirror.obj", false, true);
    rwm_mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    auto rwm_mirror_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    rwm_mirror_shape->SetMesh(rwm_mesh);
    rwm_mirror_shape->SetName("Windowless Audi");
    rwm_mirror_shape->SetStatic(true);
    rwm_mirror_shape->SetScale({1, .98, .98});
    rwm_mirror_shape->Pos = mirror_wingright_pos;
    rwm_mirror_shape->Rot = mirror_wingright_rot;  // Q_from_AngY(-.08) * Q_from_AngZ(-.25);
    rwm_mirror_shape->material_list.push_back(mirror_mat);
    vehicle.GetChassisBody()->AddAsset(rwm_mirror_shape);

    // Add leader vehicles
    std::vector<std::shared_ptr<WheeledVehicle>> lead_vehicles;
    for (int i = 0; i < num_dynamic; i++) {
        auto lead_vehicle = chrono_types::make_shared<WheeledVehicle>(vehicle.GetSystem(), vehicle_file);
        auto lead_powertrain = ReadPowertrainJSON(powertrain_file);
        // lead_vehicle->Initialize(ChCoordsys<>(dynamic_pos[i] + initRot.Rotate(ChVector<>(lead_heading * (i + 1), 0,
        // 0)), initRot));
        lead_vehicle->Initialize(ChCoordsys<>(dynamic_pos[i], initRot));
        lead_vehicle->GetChassis()->SetFixed(false);
        lead_vehicle->SetChassisVisualizationType(chassis_vis_type);
        lead_vehicle->SetSuspensionVisualizationType(suspension_vis_type);
        lead_vehicle->SetSteeringVisualizationType(steering_vis_type);
        lead_vehicle->SetWheelVisualizationType(wheel_vis_type);
        lead_vehicle->InitializePowertrain(lead_powertrain);

        // Create and initialize the tires
        for (auto& axle : lead_vehicle->GetAxles()) {
            for (auto& wheel : axle->GetWheels()) {
                auto tire = ReadTireJSON(tire_file);
                lead_vehicle->InitializeTire(tire, wheel, tire_vis_type);
            }
        }
        lead_vehicles.push_back(lead_vehicle);
    }

    // Create the terrain
    RigidTerrain terrain(vehicle.GetSystem());

    MaterialInfo minfo;
    minfo.mu = 0.9f;
    minfo.cr = 0.01f;
    minfo.Y = 2e7f;
    minfo.gt = 3000;
    auto patch_mat = minfo.CreateMaterial(contact_method);

    std::shared_ptr<RigidTerrain::Patch> patch;
    switch (terrain_model) {
        case RigidTerrain::PatchType::BOX:
            patch = terrain.AddPatch(patch_mat, ChVector<>(0, 0, 0), ChVector<>(0, 0, 1), terrainLength, terrainWidth,
                                     2, false, 1, false);
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

    // add terrain with field and grass textures
    AddTerrain(vehicle.GetSystem());

    // add signs along the road
    AddRoadway(vehicle.GetSystem());

    // Add trees near the road
    AddTrees(vehicle.GetSystem());

    // Add buildings away from the road to break up the horizon
    AddBuildings(vehicle.GetSystem());

    // ======================================================================
    // Add dummy vehicles
    // Note that dummy vhicles are placed on the inner lane of the outer loop
    // ======================================================================
    std::vector<ChBezierCurveTracker> tracker_vec;  // bezier curve tracker for dummy's path following functionality
    std::vector<ChVector<>> dummy_start;

    // tracker objects initialization
    for (int i = 0; i < num_dummy; i++) {
        if (dummy_lane[i] == 0) {
            ChBezierCurveTracker tracker(inner_path);
            tracker.setIsClosedPath(true);
            tracker_vec.push_back(tracker);
        } else {
            ChBezierCurveTracker tracker(outer_path);
            tracker.setIsClosedPath(true);
            tracker_vec.push_back(tracker);
        }
    }

    // start location initialization
    for (int i = 0; i < num_dummy; i++) {
        dummy_start.push_back(ChVector<>(0, 0, 0));
        tracker_vec[i].calcClosestPoint(dummy_pos[i], dummy_start[i]);
    }

    // vector which stores dummy vehicles
    std::vector<std::shared_ptr<ChBodyAuxRef>> dummies;

    // declare universal dummy mesh for multiple dummy objects
    // full nissan patrol mesh
    std::string suv_mesh_name = "/vehicles/Nissan_Patrol/FullPatrol.obj";
    auto suv_mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    suv_mmesh->LoadWavefrontMesh(demo_data_path + suv_mesh_name, false, true);
    suv_mmesh->RepairDuplicateVertexes(1e-9);

    auto suv_trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    suv_trimesh_shape->SetMesh(suv_mmesh);

    // full audi mesh
    std::string audi_mesh_name = "/vehicles/audi/Full_Audi.obj";
    auto audi_mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    audi_mmesh->LoadWavefrontMesh(demo_data_path + audi_mesh_name, false, true);
    audi_mmesh->RepairDuplicateVertexes(1e-9);

    auto audi_trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    audi_trimesh_shape->SetMesh(audi_mmesh);

    for (int i = 0; i < num_dummy; i++) {
        std::string mesh_name;
        float dummy_z_offset;
        if (i % 2 == 0) {
            // if the index%2 == 0, we initialize the dummy as a Nissan Patrol
            auto dummy = chrono_types::make_shared<ChBodyAuxRef>();
            dummy->SetCollide(false);
            dummy->SetPos(dummy_start[i] + ChVector<>(0, 0, dummy_patrol_z_offset));
            dummy->SetBodyFixed(true);
            dummy->AddAsset(suv_trimesh_shape);
            vehicle.GetSystem()->AddBody(dummy);
            dummies.push_back(dummy);
        } else if (i % 2 == 1) {
            // if the index%2 == 1, we initialize the dummy as an audi
            auto dummy = chrono_types::make_shared<ChBodyAuxRef>();
            dummy->SetCollide(false);
            dummy->SetPos(dummy_start[i] + ChVector<>(0, 0, dummy_audi_z_offset));
            dummy->SetBodyFixed(true);
            dummy->AddAsset(audi_trimesh_shape);
            vehicle.GetSystem()->AddBody(dummy);
            dummies.push_back(dummy);
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

    ChWheeledVehicleIrrApp app(&vehicle, L"  ", irr::core::dimension2d<irr::u32>(1360, 420));
    ChRealtimeStepTimer realtime_timer;
    /*
    SPEEDOMETER: we want to use the irrlicht app to display the speedometer, but calling endscene would update the
    entire (massive) scenario. In order to do so, we first have to clean app the Irrlichr app. Once we delete the
    node, we remove all cached meshes and textures. The order is important, otherwise meshes are re-cached!!
    */
    irr::scene::ISceneNode* mnode = app.GetContainer();
    mnode->remove();
    irr::scene::IMeshCache* cache = app.GetDevice()->getSceneManager()->getMeshCache();
    cache->clear();
    app.GetVideoDriver()->removeAllTextures();

#ifdef CHRONO_IRRKLANG
    GetLog() << "USING IRRKLANG"
             << "\n\n";
    ChCSLSoundEngine soundEng(&vehicle);
#endif

    // Create the interactive driver system
    // ChCSLDriver driver(vehicle);
    // driver = chrono_types::make_shared<ChCSLDriver>(vehicle);
    auto IGdriver = chrono_types::make_shared<ChIrrGuiDriver>(app);
    IGdriver->SetButtonCallback(r_1, &customButtonCallback);
    IGdriver->SetButtonCallback(r_3, &dummyButtonCallback_r_3);

    // for (int a = 0; a < 23; a++) {
    //    IGdriver->SetButtonCallback(a, &dummyButtonCallback_test);
    //}

    IGdriver->SetJoystickAxes(ChIrrGuiDriver::JoystickAxes::AXIS_Z, ChIrrGuiDriver::JoystickAxes::AXIS_R,
                              ChIrrGuiDriver::JoystickAxes::AXIS_X, ChIrrGuiDriver::JoystickAxes::NONE);
    IGdriver->Initialize();

    std::string steering_controller_file_IG = demo_data_path + "/Environments/Iowa/Driver/SteeringController_IG.json";
    std::string steering_controller_file_IG_nl =
        demo_data_path + "/Environments/Iowa/Driver/SteeringController_IG_nl.json";
    std::string steering_controller_file_LD = demo_data_path + "/Environments/Iowa/Driver/SteeringController_LD.json";

    std::string speed_controller_file_IG = demo_data_path + "/Environments/Iowa/Driver/SpeedController_IG.json";
    std::string speed_controller_file_IG_nl = demo_data_path + "/Environments/Iowa/Driver/SpeedController_IG_nl.json";
    std::string speed_controller_file_LD = demo_data_path + "/Environments/Iowa/Driver/SpeedController_LD.json";

    std::cout << "lead_count: " << lead_count << std::endl;
    // lead_count
    std::shared_ptr<ChNSFFollowererDriver> PFdriver;
    if (lead_count == 0) {
        PFdriver = chrono_types::make_shared<ChNSFFollowererDriver>(vehicle, steering_controller_file_IG_nl,
                                                                    speed_controller_file_IG_nl, outer_path, "road",
                                                                    cruise_speed * MPH_TO_MS, followerParam, true);
    } else {
        PFdriver = chrono_types::make_shared<ChNSFFollowererDriver>(
            vehicle, steering_controller_file_IG, speed_controller_file_IG, outer_path, "road",
            cruise_speed * MPH_TO_MS, lead_vehicles[0], followerParam, true);
    }

    PFdriver->Initialize();

    PF_driver_ptr = PFdriver;

    // we call the callback explicitly to start the timer
    customButtonCallback();
    dummyButtonCallback_r_3();

    if (!disable_joystick) {
        driver_mode = HUMAN;
    } else {
        driver_mode = AUTONOMOUS;
        std::cout << "Using path follower driver\n";
    }

    // Leader Driver
    std::vector<std::shared_ptr<ChNSFLeaderDriver>> lead_PFdrivers;
    for (int i = 0; i < num_dynamic; i++) {
        if (dynamic_lane[i] == 0) {
            auto lead_PFdriver = chrono_types::make_shared<ChNSFLeaderDriver>(
                *lead_vehicles[i], steering_controller_file_LD, speed_controller_file_LD, inner_path, "road",
                dynamic_cruise_speed[i] * MPH_TO_MS, leaderParam[i], true);
            lead_PFdriver->Initialize();
            lead_PFdrivers.push_back(lead_PFdriver);
        } else {
            auto lead_PFdriver = chrono_types::make_shared<ChNSFLeaderDriver>(
                *lead_vehicles[i], steering_controller_file_LD, speed_controller_file_LD, outer_path, "road",
                dynamic_cruise_speed[i] * MPH_TO_MS, leaderParam[i], true);
            lead_PFdriver->Initialize();
            lead_PFdrivers.push_back(lead_PFdriver);
        }
    }

    if (save_driver) {
        filestream << "csv comments: " << cli.GetAsType<std::string>("csv_comments") << " \n";

        if (lead_count != 0) {
            filestream << "tstamp,time,wallTime,isManual,Steering,Throttle,Braking,x[m],y[m],speed[mph],"
                          "acceleration[m/s^2],"
                          "dist[m],dist_projected[m],IG_mile[mile],IG_lane[0-inner/1-outer/-1-invalid],LD_x[m],"
                          "LD_y[m],LD_speed[mph],"
                          "LD_acc[m/s^2],LD_mile[mile]\n";
        } else {
            filestream << "tstamp,time,wallTime,isManual,Steering,Throttle,Braking,x[m],y[m],speed[mph],"
                          "acceleration[m/s^2],"
                          "IG_mile[mile],IG_lane[0-inner/1-outer/-1-invalid]\n";
        }

        buttonstream << "start recording buttons pressed \n";
    }

    // ---------------
    // Simulation loop
    // ---------------

    // output vehicle mass
    std::cout << "VEHICLE MASS: " << vehicle.GetVehicleMass() << std::endl;

    // Initialize simulation frame counter and simulation time
    int step_number = 0;
    int render_frame = 0;
    double sim_time = 0;

    // ---------------------------------------------
    // Create a sensor manager and add a point light
    // ---------------------------------------------
    auto manager = chrono_types::make_shared<ChSensorManager>(vehicle.GetSystem());
    float intensity = 2.0;
    manager->scene->AddPointLight({0, 0, 1e8}, {intensity, intensity, intensity}, 1e12);
    manager->scene->SetAmbientLight({.1, .1, .1});
    manager->scene->SetSceneEpsilon(1e-3);
    manager->scene->EnableDynamicOrigin(true);
    manager->scene->SetOriginOffsetThreshold(500.f);

    // Set environment map
    Background b;
    b.mode = BackgroundMode::ENVIRONMENT_MAP;
    b.env_tex = GetChronoDataFile("sensor/textures/sunflowers_4k.hdr");
    manager->scene->SetBackground(b);
    if (fog_enabled) {
        manager->scene->SetFogScatteringFromDistance(fog_distance);
        manager->scene->SetFogColor(fog_color);
    }

    // ------------------------------------------------
    // Create a camera and add it to the sensor manager
    // ------------------------------------------------
    auto cam = chrono_types::make_shared<ChCameraSensor>(
        vehicle.GetChassisBody(),                                                     // body camera is attached to
        10,                                                                           // update rate in Hz
        chrono::ChFrame<double>({0, 0, 3000}, Q_from_AngAxis(CH_C_PI_2, {0, 1, 0})),  // offset pose
        1920,                                                                         // image width
        1080,                                                                         // image height
        CH_C_PI_4,
        super_samples);  // fov, lag, exposure
    cam->SetName("Camera Sensor");
    if (sensor_vis)
        cam->PushFilter(chrono_types::make_shared<ChFilterVisualize>(1280, 720));

    // add sensor to the manager
    if (cli.GetAsType<bool>("birdseye"))
        manager->AddSensor(cam);

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

    ChVector<> prev_IG_pos;
    bool IG_started_driving = false;

    while (app.GetDevice()->run()) {
        sim_time = vehicle.GetSystem()->GetChTime();

        // End simulation
        if (sim_time >= t_end)
            break;

        // cam->SetOffsetPose(
        //     chrono::ChFrame<double>({-orbit_radius * cos(time * orbit_rate), -orbit_radius * sin(time *
        //     orbit_rate), orbit_radius/5.0},
        //                             Q_from_AngAxis(time * orbit_rate, {0, 0,1})));

        // Collect output data from modules (for inter-module communication)
        ChDriver::Inputs driver_inputs;
        if (driver_mode == AUTONOMOUS)
            driver_inputs = PFdriver->GetInputs();
        else {
            driver_inputs = IGdriver->GetInputs();
            driver_inputs.m_steering *= -1;
        }

        // printf("Driver inputs: %f,%f,%f\n", driver_inputs.m_throttle, driver_inputs.m_braking,
        //        driver_inputs.m_steering);
        // driver_inputs.m_throttle = 0;
        // driver_inputs.m_steering *= -1;
        if (step_number % int(1 / step_size) == 0 && !benchmark) {
            std::cout << "dynamic dist: " << lead_PFdrivers[0]->Get_Dist() << std::endl;

            if (lead_count != 0) {
                auto ld_speed = lead_vehicles[0]->GetVehicleSpeed() * MS_TO_MPH;
                auto ig_speed = vehicle.GetVehicleSpeed() * MS_TO_MPH;
                auto wall_time = high_resolution_clock::now();
                printf("Sim Time=%f, \tWall Time=%f, \tExtra Time=%f, \tLD_Speed mph=%f, \tIG_Speed mph=%f\n", sim_time,
                       duration_cast<duration<double>>(wall_time - t0).count(), extra_time, ld_speed, ig_speed);
                // std::cout << "Current Gear: " << vehicle.GetPowertrain()->GetCurrentTransmissionGear() << std::endl;
                extra_time = 0.0;
            } else {
                auto ig_speed = vehicle.GetVehicleSpeed() * MS_TO_MPH;
                auto wall_time = high_resolution_clock::now();
                printf("Sim Time=%f, \tWall Time=%f, \tExtra Time=%f, \tIG_Speed mph=%f\n", sim_time,
                       duration_cast<duration<double>>(wall_time - t0).count(), extra_time, ig_speed);
                // std::cout << "Current Gear: " << vehicle.GetPowertrain()->GetCurrentTransmissionGear() << std::endl;
                extra_time = 0.0;
            }
        }

        // update current vehicle speed
        cur_follower_speed = vehicle.GetVehicleSpeed();

        // Update modules (process inputs from other modules)
        if (driver_mode == AUTONOMOUS)
            PFdriver->Synchronize(sim_time, step_size);
        else
            IGdriver->Synchronize(sim_time);

        terrain.Synchronize(sim_time);
        vehicle.Synchronize(sim_time, driver_inputs, terrain);

        app.Synchronize("", driver_inputs);
#ifdef CHRONO_IRRKLANG
        soundEng.Synchronize(sim_time);
#endif
        // Advance simulation for one timestep for all modules
        double step = step_size;

        if (driver_mode == AUTONOMOUS)
            PFdriver->Advance(step);
        else
            IGdriver->Advance(step);

        terrain.Advance(step);
        vehicle.Advance(step);

        app.Advance(step_size);

        auto t2 = high_resolution_clock::now();
        float wall_time = duration_cast<duration<double>>(t2 - t0).count();

        // dummy update
        for (int i = 0; i < num_dummy; i++) {
            float temp_z_offset;
            if (i % 2 == 0) {
                temp_z_offset = dummy_patrol_z_offset;
            } else if (i % 2 == 1) {
                temp_z_offset = dummy_audi_z_offset;
            }

            if (dummy_lane[i] == 0) {
                // when in inner lane
                if (dummy_control[i] == 0) {
                    // if dummy doesn't have any control, we just use cruise speed setting always
                    updateDummy(dummies[i], inner_path, dummy_cruise_speed[i], step_size, temp_z_offset, tracker_vec[i],
                                dummy_dist[i], dummy_prev_pos[i]);
                } else if (dummy_control[i] == 1) {
                    float target_speed;
                    target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i], dummy_dist[i],
                                                        dummy_cruise_speed[i]);
                    target_speed = target_speed * MPH_TO_MS;
                    updateDummy(dummies[i], inner_path, target_speed, step_size, temp_z_offset, tracker_vec[i],
                                dummy_dist[i], dummy_prev_pos[i]);
                } else {
                    float target_speed;
                    if (dummy_time_mode[i] == 1) {
                        target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i], sim_time,
                                                            dummy_cruise_speed[i]);
                    } else if (dummy_time_mode[i] == 2) {
                        target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i], wall_time,
                                                            dummy_cruise_speed[i]);
                    } else if (dummy_time_mode[i] == 3) {
                        target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i],
                                                            sim_time - start_sim_time, dummy_cruise_speed[i]);
                    } else if (dummy_time_mode[i] == 4) {
                        target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i],
                                                            wall_time - start_wall_time, dummy_cruise_speed[i]);
                    }

                    target_speed = target_speed * MPH_TO_MS;
                    // if dummy has a control parameter setting, we adjust speed based on given data
                    updateDummy(dummies[i], inner_path, target_speed, step_size, temp_z_offset, tracker_vec[i],
                                dummy_dist[i], dummy_prev_pos[i]);
                }

            } else {
                // when in outer lane
                if (dummy_control[i] == 0) {
                    // if dummy doesn't have any control, we just use cruise speed setting always
                    updateDummy(dummies[i], outer_path, dummy_cruise_speed[i], step_size, temp_z_offset, tracker_vec[i],
                                dummy_dist[i], dummy_prev_pos[i]);
                } else if (dummy_control[i] == 1) {
                    float target_speed;
                    target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i], dummy_dist[i],
                                                        dummy_cruise_speed[i]);
                    target_speed = target_speed * MPH_TO_MS;
                    updateDummy(dummies[i], outer_path, target_speed, step_size, temp_z_offset, tracker_vec[i],
                                dummy_dist[i], dummy_prev_pos[i]);
                } else {
                    float target_speed;
                    if (dummy_time_mode[i] == 1) {
                        target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i], sim_time,
                                                            dummy_cruise_speed[i]);
                    } else if (dummy_time_mode[i] == 2) {
                        target_speed = controlFindSpeed_x_y(dummy_control_x[i], dummy_control_y[i], wall_time,
                                                            dummy_cruise_speed[i]);
                    }
                    target_speed = target_speed * MPH_TO_MS;
                    updateDummy(dummies[i], outer_path, target_speed, step_size, temp_z_offset, tracker_vec[i],
                                dummy_dist[i], dummy_prev_pos[i]);
                }
            }
        }

        for (int i = 0; i < num_dynamic; i++) {
            // set speed control
            if (dynamic_control[i] == 2) {
                float target_speed;
                if (dynamic_time_mode[i] == 1) {
                    target_speed = controlFindSpeed_x_y(dynamic_control_x[i], dynamic_control_y[i], sim_time,
                                                        dynamic_cruise_speed[i]);
                } else if (dynamic_time_mode[i] == 2) {
                    target_speed = controlFindSpeed_x_y(dynamic_control_x[i], dynamic_control_y[i], wall_time,
                                                        dynamic_cruise_speed[i]);
                } else if (dynamic_time_mode[i] == 3) {
                    target_speed = controlFindSpeed_x_y(dynamic_control_x[i], dynamic_control_y[i],
                                                        sim_time - start_sim_time, dynamic_cruise_speed[i]);
                } else if (dynamic_time_mode[i] == 4) {
                    target_speed = controlFindSpeed_x_y(dynamic_control_x[i], dynamic_control_y[i],
                                                        wall_time - start_wall_time, dynamic_cruise_speed[i]);
                }
                lead_PFdrivers[i]->SetCruiseSpeed(target_speed * MPH_TO_MS);
            } else if (dynamic_control[i] == 1) {
                float target_speed;
                target_speed = controlFindSpeed_x_y(dynamic_control_x[i], dynamic_control_y[i],
                                                    lead_PFdrivers[i]->Get_Dist(), dynamic_cruise_speed[i]);

                lead_PFdrivers[i]->SetCruiseSpeed(target_speed * MPH_TO_MS);
            }
            ChDriver::Inputs lead_driver_inputs = lead_PFdrivers[i]->GetInputs();
            lead_PFdrivers[i]->Synchronize(sim_time);
            lead_vehicles[i]->Synchronize(sim_time, lead_driver_inputs, terrain);
            lead_PFdrivers[i]->Advance(step);
            lead_vehicles[i]->Advance(step);
        }

        // std::cout << "sim_time: " << sim_time
        //           << "  lead_speed: " << lead_vehicles[0]->GetChassis()->GetSpeed() * MS_TO_MPH << std::endl;

        if (step_number % int(1 / (60 * step_size)) == 0) {
            /// irrlicht::tools::drawSegment(app.GetVideoDriver(), v1, v2, video::SColor(255, 80, 0, 0), false);
            app.GetDevice()->getVideoDriver()->draw2DImage(
                app.GetDevice()->getVideoDriver()->getTexture((demo_data_path + "/miscellaneous/dash_4_4.jpg").c_str()),
                irr::core::position2d<irr::s32>(0, 0));

            int curr_gear = vehicle.GetPowertrain()->GetCurrentTransmissionGear();

            switch (curr_gear) {
                case 1:
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/GEAR1.png").c_str()),
                        gr_center);
                    break;

                case 2:
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/GEAR2.png").c_str()),
                        gr_center);
                    break;

                case 3:
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/GEAR3.png").c_str()),
                        gr_center);
                    break;

                case 4:
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/GEAR4.png").c_str()),
                        gr_center);
                    break;

                case 5:
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/GEAR5.png").c_str()),
                        gr_center);
                    break;

                case 6:
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/GEAR6.png").c_str()),
                        gr_center);
                    break;

                default:
                    break;
            }

            double speed_mph = vehicle.GetVehicleSpeedCOM() * MS_TO_MPH;
            double theta = ((265 / 130) * speed_mph) * (CH_C_PI / 180);
            app.GetDevice()->getVideoDriver()->draw2DLine(
                sm_center + irr::core::position2d<irr::s32>(-sm_needle * sin(theta), sm_needle * cos(theta)), sm_center,
                irr::video::SColor(255, 255, 0, 0));

            double engine_rpm = powertrain->GetMotorSpeed() * rads2rpm;
            double alpha = ((265.0 / 6500.0) * engine_rpm) * (CH_C_PI / 180);
            app.GetDevice()->getVideoDriver()->draw2DLine(
                rpm_center + irr::core::position2d<irr::s32>(-sm_needle * sin(alpha), sm_needle * cos(alpha)),
                rpm_center, irr::video::SColor(255, 255, 0, 0));

            if (driver_mode == AUTONOMOUS) {
                app.GetDevice()->getVideoDriver()->draw2DImage(
                    app.GetDevice()->getVideoDriver()->getTexture((demo_data_path + "/miscellaneous/auto.png").c_str()),
                    auto_center);
            }

            // display current time
            time_t curr_time = time(NULL);
            struct tm* tmp = localtime(&curr_time);

            int h = tmp->tm_hour;
            int m = tmp->tm_min;
            int s = tmp->tm_sec;

            std::vector<int> display_cur_int;
            display_cur_int.push_back(h / 10);
            display_cur_int.push_back(h % 10);
            display_cur_int.push_back(m / 10);
            display_cur_int.push_back(m % 10);
            display_cur_int.push_back(s / 10);
            display_cur_int.push_back(s % 10);

            for (int i = 0; i < display_cur_int.size() + 2; i++) {
                irr::core::position2d<irr::s32> offset(50, 0);
                irr::core::position2d<irr::s32> colon_offset1(25, 0);
                irr::core::position2d<irr::s32> colon_offset2(25, 0);

                if (i < 2) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/numerical/" + std::to_string(display_cur_int[i]) + ".png")
                                .c_str()),
                        cur_left + offset * i);
                } else if (i == 2) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/numerical/colon.png").c_str()),
                        cur_left + offset * i);
                } else if (i < 5) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture((demo_data_path + "/miscellaneous/numerical/" +
                                                                       std::to_string(display_cur_int[i - 1]) + ".png")
                                                                          .c_str()),
                        cur_left + offset * (i - 1) + colon_offset1);
                } else if (i == 5) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/numerical/colon.png").c_str()),
                        cur_left + offset * (i - 1) + colon_offset1);
                } else if (i <= 7) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture((demo_data_path + "/miscellaneous/numerical/" +
                                                                       std::to_string(display_cur_int[i - 2]) + ".png")
                                                                          .c_str()),
                        cur_left + offset * (i - 2) + colon_offset1 + colon_offset2);
                }
            }

            static int end_s;
            static int end_m;
            static int end_h;

            if (vehicle.GetVehicleSpeedCOM() >= 0.5 && IG_started_driving == false) {
                start_sim_time = sim_time;
                auto t1_temp = high_resolution_clock::now();
                start_wall_time += duration_cast<duration<double>>(t1_temp - t0).count();

                time_t curr_time = time(NULL);
                struct tm* tmp = localtime(&curr_time);

                int temp_h = tmp->tm_hour;
                int temp_m = tmp->tm_min;
                int temp_s = tmp->tm_sec;

                end_s = s;
                end_m = m + meet_time;
                end_h = h;

                if (end_s >= 60) {
                    end_s = end_s - 60;
                    end_m = end_m + 1;
                }
                if (end_m >= 60) {
                    end_m = end_m - 60;
                    end_h = end_h + 1;
                }
                if (end_h >= 23) {
                    end_h = end_h % 24;
                }

                IG_started_driving = true;
            }

            if (IG_started_driving == true) {
                std::vector<int> display_end_int;
                display_end_int.push_back(end_h / 10);
                display_end_int.push_back(end_h % 10);
                display_end_int.push_back(end_m / 10);
                display_end_int.push_back(end_m % 10);
                display_end_int.push_back(end_s / 10);
                display_end_int.push_back(end_s % 10);

                for (int i = 0; i < display_end_int.size() + 2; i++) {
                    irr::core::position2d<irr::s32> offset(50, 0);
                    irr::core::position2d<irr::s32> colon_offset1(25, 0);
                    irr::core::position2d<irr::s32> colon_offset2(25, 0);

                    if (i < 2) {
                        app.GetDevice()->getVideoDriver()->draw2DImage(
                            app.GetDevice()->getVideoDriver()->getTexture((demo_data_path +
                                                                           "/miscellaneous/numerical/" +
                                                                           std::to_string(display_end_int[i]) + ".png")
                                                                              .c_str()),
                            end_left + offset * i);
                    } else if (i == 2) {
                        app.GetDevice()->getVideoDriver()->draw2DImage(
                            app.GetDevice()->getVideoDriver()->getTexture(
                                (demo_data_path + "/miscellaneous/numerical/colon.png").c_str()),
                            end_left + offset * i);
                    } else if (i < 5) {
                        app.GetDevice()->getVideoDriver()->draw2DImage(
                            app.GetDevice()->getVideoDriver()->getTexture(
                                (demo_data_path + "/miscellaneous/numerical/" + std::to_string(display_end_int[i - 1]) +
                                 ".png")
                                    .c_str()),
                            end_left + offset * (i - 1) + colon_offset1);
                    } else if (i == 5) {
                        app.GetDevice()->getVideoDriver()->draw2DImage(
                            app.GetDevice()->getVideoDriver()->getTexture(
                                (demo_data_path + "/miscellaneous/numerical/colon.png").c_str()),
                            end_left + offset * (i - 1) + colon_offset1);
                    } else if (i <= 7) {
                        app.GetDevice()->getVideoDriver()->draw2DImage(
                            app.GetDevice()->getVideoDriver()->getTexture(
                                (demo_data_path + "/miscellaneous/numerical/" + std::to_string(display_end_int[i - 2]) +
                                 ".png")
                                    .c_str()),
                            end_left + offset * (i - 2) + colon_offset1 + colon_offset2);
                    }
                }
            }

            // compute and display ETA

            static int sec_remaining = 0;

            if (step_number == 0) {
                IG_prev_pos = ego_chassis->GetPos();
            }

            IG_dist = IG_dist + (ego_chassis->GetPos() - IG_prev_pos).Length();
            IG_prev_pos = ego_chassis->GetPos();

            if (step_number % 50 == 0) {
                float remaining = eta_dist * MILE_TO_M - IG_dist;
                float avg_speed = IG_speed_avg.Add(ego_chassis->GetSpeed());
                sec_remaining = remaining / avg_speed;
            }

            // panic algorithm
            // if below 0, set to 0
            // if above max, set to 0

            if (sec_remaining < 0 || sec_remaining > 356518) {
                sec_remaining = 0;
            }

            int eta_h = sec_remaining / 3600;
            int eta_m = (sec_remaining % 3600) / 60;
            int eta_s = sec_remaining % 60;

            std::vector<int> display_eta_int;
            display_eta_int.push_back(eta_h / 10);
            display_eta_int.push_back(eta_h % 10);
            display_eta_int.push_back(eta_m / 10);
            display_eta_int.push_back(eta_m % 10);
            display_eta_int.push_back(eta_s / 10);
            display_eta_int.push_back(eta_s % 10);

            for (int i = 0; i < display_eta_int.size() + 2; i++) {
                irr::core::position2d<irr::s32> offset(50, 0);
                irr::core::position2d<irr::s32> colon_offset1(25, 0);
                irr::core::position2d<irr::s32> colon_offset2(25, 0);

                if (i < 2) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/numerical/" + std::to_string(display_eta_int[i]) + ".png")
                                .c_str()),
                        eta_left + offset * i);
                } else if (i == 2) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/numerical/colon.png").c_str()),
                        eta_left + offset * i);
                } else if (i < 5) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture((demo_data_path + "/miscellaneous/numerical/" +
                                                                       std::to_string(display_eta_int[i - 1]) + ".png")
                                                                          .c_str()),
                        eta_left + offset * (i - 1) + colon_offset1);
                } else if (i == 5) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture(
                            (demo_data_path + "/miscellaneous/numerical/colon.png").c_str()),
                        eta_left + offset * (i - 1) + colon_offset1);
                } else if (i <= 7) {
                    app.GetDevice()->getVideoDriver()->draw2DImage(
                        app.GetDevice()->getVideoDriver()->getTexture((demo_data_path + "/miscellaneous/numerical/" +
                                                                       std::to_string(display_eta_int[i - 2]) + ".png")
                                                                          .c_str()),
                        eta_left + offset * (i - 2) + colon_offset1 + colon_offset2);
                }
            }

            app.GetDevice()->getVideoDriver()->endScene();
        }

        // Update the sensor manager
        manager->Update();

        // Increment frame number
        step_number++;

        if (step_number % (int)(2.0 / frame_rate / step_size) == 0 && !benchmark) {
            double since_last_sync = sim_time - last_sim_sync;
            last_sim_sync = sim_time;
            auto tt0 = high_resolution_clock::now();
            realtime_timer.Spin(since_last_sync);
            auto tt1 = high_resolution_clock::now();
            extra_time += duration_cast<duration<double>>(tt1 - tt0).count();
        }

        ChBezierCurveTracker lane_0_tracker(lane_0_path);
        ChBezierCurveTracker lane_1_tracker(lane_1_path);

        lane_0_tracker.setIsClosedPath(true);
        lane_1_tracker.setIsClosedPath(true);

        if (save_driver) {
            if (step_number % int(tsave / step_size) == 0) {
                buffer << std::fixed << std::setprecision(3);

                time_t my_time = time(NULL);
                buffer << strtok(ctime(&my_time), "\n");
                buffer << ",";
                buffer << std::to_string(sim_time) + ",";
                buffer << std::to_string(wall_time) << ",";
                ChDriver* currDriver;
                bool isManual;
                if (driver_mode == HUMAN) {
                    currDriver = IGdriver.get();
                    isManual = true;
                } else {
                    currDriver = PFdriver.get();
                    isManual = false;
                }

                buffer << isManual << ",";
                buffer << currDriver->GetSteering() << ",";
                buffer << currDriver->GetThrottle() << ",";
                buffer << currDriver->GetBraking() << ",";
                buffer << ego_chassis->GetPos().x() << ",";
                buffer << ego_chassis->GetPos().y() << ",";
                buffer << ego_chassis->GetSpeed() * MS_TO_MPH << ",";
                buffer << ego_chassis->GetBody()->GetFrame_REF_to_abs().GetPos_dtdt().Length() << ",";

                if (lead_count != 0) {
                    // Obtain lead vehicle chassis
                    auto lead_chassis = lead_vehicles[lead_vehicles.size() - 1]->GetChassis();
                    double dist =
                        (ego_chassis->GetPos() - lead_vehicles[0]->GetChassis()->GetPos()).Length() - AUDI_LENGTH;
                    buffer << dist << ",";  // Distance bumper-to-bumber
                    ChVector<> dist_v = lead_vehicles[0]->GetChassis()->GetPos() - ego_chassis->GetPos();
                    ChVector<> car_xaxis = ChMatrix33<>(ego_chassis->GetRot()).Get_A_Xaxis();
                    double proj_dist = (dist_v ^ car_xaxis) - AUDI_LENGTH;
                    buffer << proj_dist << ",";  // Projected distance bumper-to-bumber
                }

                // output mile marker
                buffer << IG_dist * M_TO_MILE << ",";

                ChVector<> lane_0_target;
                ChVector<> lane_1_target;

                lane_0_tracker.calcClosestPoint(ego_chassis->GetPos(), lane_0_target);
                lane_1_tracker.calcClosestPoint(ego_chassis->GetPos(), lane_1_target);

                ChVector<> chassis_pos = ego_chassis->GetPos();
                float dist_0 = (lane_0_target.x() - chassis_pos.x()) * (lane_0_target.x() - chassis_pos.x()) +
                               (lane_0_target.y() - chassis_pos.y()) * (lane_0_target.y() - chassis_pos.y());

                float dist_1 = (lane_1_target.x() - chassis_pos.x()) * (lane_1_target.x() - chassis_pos.x()) +
                               (lane_1_target.y() - chassis_pos.y()) * (lane_1_target.y() - chassis_pos.y());

                int lane_num = -2;
                float min_dist;
                if (dist_0 < dist_1) {
                    lane_num = 0;
                    min_dist = dist_0;
                } else {
                    lane_num = 1;
                    min_dist = dist_1;
                }

                if (min_dist > 5) {
                    lane_num = -1;
                }

                // output IG vehicle lane number
                buffer << lane_num << ",";

                // the last lead vehicle data
                if (lead_count != 0) {
                    // Obtain lead vehicle chassis
                    auto lead_chassis = lead_vehicles[lead_vehicles.size() - 1]->GetChassis();
                    buffer << lead_chassis->GetPos().x() << ",";
                    buffer << lead_chassis->GetPos().y() << ",";
                    buffer << lead_chassis->GetSpeed() * MS_TO_MPH << ",";
                    buffer << lead_chassis->GetBody()->GetFrame_REF_to_abs().GetPos_dtdt().Length()
                           << ",";  // output mile marker
                    buffer << lead_PFdrivers[lead_PFdrivers.size() - 1]->Get_Dist();
                }

                buffer << "\n";
                if ((step_number % int(20 / step_size) == 0)) {
                    printf("Writing to output file...=%i", buffer.tellp());
                    filestream << buffer.rdbuf();
                    buffer.str("");

                    printf("Writing to button file...=%i", button_buffer.tellp());
                    buttonstream << button_buffer.rdbuf();
                    buttonstream.flush();
                    button_buffer.str("");
                }
            }
        }
    }

    auto t1 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t1 - t0);
    std::cout << "Simulation Time: " << t_end << ", Wall Time: " << time_span.count() << std::endl;

    return 0;
}

void AddTrees(ChSystem* chsystem) {
    // add in tress next to road
    auto tree_mesh_0 = chrono_types::make_shared<ChTriangleMeshConnected>();
    tree_mesh_0->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/trees/tree_01.obj", false, true);
    tree_mesh_0->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    auto tree_mesh_1 = chrono_types::make_shared<ChTriangleMeshConnected>();
    tree_mesh_1->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/foliage/trees/oaktree_01.obj", false, true);
    tree_mesh_1->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    auto tree_mesh_2 = chrono_types::make_shared<ChTriangleMeshConnected>();
    tree_mesh_2->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/foliage/trees/oaktree_02.obj", false, true);
    tree_mesh_2->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    auto tree_mesh_3 = chrono_types::make_shared<ChTriangleMeshConnected>();
    tree_mesh_3->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/foliage/trees/oaktree_03.obj", false, true);
    tree_mesh_3->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size

    std::vector<std::shared_ptr<ChTriangleMeshConnected>> tree_meshes = {tree_mesh_0, tree_mesh_1, tree_mesh_2,
                                                                         tree_mesh_3};

    ChSetRandomSeed(4);

    // tree placement parameters
    double x_step = 90;
    double y_step = 400;
    int x_count = 2;
    int y_count = 50;
    float y_variation = 300.f;
    float scale_nominal = .5;
    float scale_variation = .4f;
    float x_variation = .1f;

    double x_start = 4000.0 - 45;  // initLoc.x()-38;
    double y_start = -10000.0;     // initLoc.y()-7500;

    for (int i = 0; i < x_count; i++) {
        for (int j = 0; j < y_count; j++) {
            auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
            trimesh_shape->SetMesh(tree_meshes[int(ChRandom() * tree_meshes.size() - .001)]);
            trimesh_shape->SetName("Tree");
            trimesh_shape->SetStatic(true);
            float scale = scale_nominal + scale_variation * (ChRandom() - .5);
            trimesh_shape->SetScale({scale, scale, scale});

            auto mesh_body = chrono_types::make_shared<ChBody>();
            mesh_body->SetPos({i * x_step + x_start + x_variation * (ChRandom() - .5),
                               j * y_step + y_start + y_variation * (ChRandom() - .5), 0.0});
            mesh_body->SetRot(Q_from_AngZ(CH_C_PI_2 * ChRandom()));
            mesh_body->AddAsset(trimesh_shape);
            mesh_body->SetBodyFixed(true);
            chsystem->Add(mesh_body);
        }
    }

    x_start = -4052.0 - 45;  // initLoc.x()-38;
    y_start = -10000.0;      // initLoc.y()-7500;

    for (int i = 0; i < x_count; i++) {
        for (int j = 0; j < y_count; j++) {
            auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
            trimesh_shape->SetMesh(tree_meshes[int(ChRandom() * tree_meshes.size() - .001)]);
            trimesh_shape->SetName("Tree");
            trimesh_shape->SetStatic(true);
            float scale = scale_nominal + scale_variation * (ChRandom() - .5);
            trimesh_shape->SetScale({scale, scale, scale});

            auto mesh_body = chrono_types::make_shared<ChBody>();
            mesh_body->SetPos({i * x_step + x_start + x_variation * (ChRandom() - .5),
                               j * y_step + y_start + y_variation * (ChRandom() - .5), 0.0});
            mesh_body->SetRot(Q_from_AngZ(CH_C_PI_2 * ChRandom()));
            mesh_body->AddAsset(trimesh_shape);
            mesh_body->SetBodyFixed(true);
            chsystem->Add(mesh_body);
        }
    }
}

void AddRoadway(ChSystem* chsystem) {
    std::vector<std::string> environment_meshes = {"/Environments/Iowa/signs/mile_markers_inner.obj",
                                                   "/Environments/Iowa/signs/mile_markers_outer.obj",
                                                   "/Environments/Iowa/terrain/oval_highway.obj"};
    std::vector<ChVector<>> offsets = {{0, 0, -128.22}, {0, 0, 0.0}, {0, 0, 0.01}};

    for (int i = 0; i < environment_meshes.size(); i++) {  // auto file_name : environment_meshes) {
        // additional environment assets
        auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        trimesh->LoadWavefrontMesh(demo_data_path + environment_meshes[i], false, true);
        trimesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
        auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
        trimesh_shape->SetMesh(trimesh);
        trimesh_shape->SetName(environment_meshes[i]);
        trimesh_shape->SetStatic(true);
        auto mesh_body = chrono_types::make_shared<ChBody>();
        mesh_body->SetPos(offsets[i]);
        mesh_body->AddAsset(trimesh_shape);
        mesh_body->SetBodyFixed(true);
        chsystem->Add(mesh_body);
    }

    // Add arrived sign
    {
        std::string meshname = "/Environments/Iowa/signs/arrived.obj";
        auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        trimesh->LoadWavefrontMesh(demo_data_path + meshname, false, true);
        trimesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
        auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
        trimesh_shape->SetMesh(trimesh);
        trimesh_shape->SetName(meshname);
        trimesh_shape->SetStatic(true);
        auto mesh_body = chrono_types::make_shared<ChBody>();
        mesh_body->SetPos(arrived_sign_pos);
        mesh_body->SetRot(arrived_sign_rot);
        mesh_body->AddAsset(trimesh_shape);
        mesh_body->SetBodyFixed(true);
        chsystem->Add(mesh_body);
    }
}

void AddBuildings(ChSystem* chsystem) {
    std::vector<std::string> environment_meshes = {"/Environments/Iowa/buildings/farm_01.obj",  //
                                                   "/Environments/Iowa/buildings/farm_02.obj",  //
                                                   "/Environments/Iowa/buildings/farm_04.obj",  //
                                                   "/Environments/Iowa/buildings/farm_03.obj",  //
                                                   "/Environments/Iowa/buildings/farm_02.obj",  //
                                                   "/Environments/Iowa/buildings/farm_04.obj",  //
                                                   "/Environments/Iowa/buildings/farm_01.obj",  //
                                                   "/Environments/Iowa/buildings/farm_03.obj",  //
                                                   "/Environments/Iowa/buildings/farm_02.obj",  //
                                                   "/Environments/Iowa/buildings/farm_04.obj",  //
                                                   "/Environments/Iowa/buildings/farm_03.obj",  //
                                                   "/Environments/Iowa/buildings/farm_02.obj",  //
                                                   "/Environments/Iowa/buildings/farm_01.obj",  //
                                                   "/Environments/Iowa/buildings/farm_04.obj",  //
                                                   "/Environments/Iowa/buildings/farm_03.obj",  //
                                                   "/Environments/Iowa/buildings/farm_04.obj",  //
                                                   "/Environments/Iowa/buildings/farm_02.obj",  //
                                                   "/Environments/Iowa/buildings/farm_03.obj",  //
                                                   "/Environments/Iowa/buildings/farm_01.obj",  //
                                                   "/Environments/Iowa/buildings/farm_04.obj",  //
                                                   "/Environments/Iowa/buildings/farm_02.obj",  //
                                                   "/Environments/Iowa/buildings/farm_01.obj",  //
                                                   "/Environments/Iowa/buildings/farm_03.obj",  //

                                                   "/Environments/Iowa/buildings/radio_tower.obj",  //
                                                   "/Environments/Iowa/buildings/radio_tower.obj",  //
                                                   "/Environments/Iowa/buildings/radio_tower.obj",  //
                                                   "/Environments/Iowa/buildings/radio_tower.obj",  //
                                                   "/Environments/Iowa/buildings/water_tower.obj",  //
                                                   "/Environments/Iowa/buildings/water_tower.obj"};
    std::vector<ChVector<>> offsets = {
        {4150, 12776, 0},            // farm 4000
        {3800, 10209, 0},            // farm
        {3920, 8372, 0},             // farm
        {4200, 5251, 0},             // farm
        {3850, 500, 0.0},            // farm
        {4120, -1787, 0},            // farm
        {3900, -3623, 0},            // farm
        {3860, -5458, 0},            // farm
        {4180, -8000, 0},            // farm
        {3510, -14004, 0},           // farm
        {-1832, -15665, 0},          // farm
        {-4050 + 200, -10654, 0.0},  // farm -4050 200
        {-4050 + 180, -8683, 0.0},   // farm 180
        {-4050 - 120, -6634, 0.0},   // farm -120
        {-4050 + 150, -2990, 0.0},   // farm 150
        {-4050 - 120, -1040, 0.0},   // farm -120
        {-4050 - 180, -797, 0.0},    // farm -180
        {-4050 + 160, 2626, 0.0},    // farm 160
        {-4050 - 110, 4461, 0.0},    // farm -110
        {-4050 + 130, 6292, 0.0},    // farm 130
        {-4050 + 100, 8730, 0.0},    // farm 100
        {-2602, 15320, 0.0},         // farm
        {-1312, 15459, 0.0},         // farm
        {-4167, 7087, 0.0},          // radio tower
        {-4167, 10630, 0.0},         // radio tower
        {4100, -10630, 0.0},         // radio tower
        {4100, 3543, 0.0},           // radio tower
        {2657, 14488, 0.0},          // water tower
        {-2922, -14611, 0.0}         // water tower
    };

    if (offsets.size() != environment_meshes.size()) {
        std::cout << "ERROR: incorrect number of offsets for building meshes\n";
        return;
    }
    for (int i = 0; i < environment_meshes.size(); i++) {
        // additional environment assets
        auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        trimesh->LoadWavefrontMesh(demo_data_path + environment_meshes[i], false, true);
        trimesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
        auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
        trimesh_shape->SetMesh(trimesh);
        trimesh_shape->SetName(environment_meshes[i]);
        trimesh_shape->SetStatic(true);
        trimesh_shape->SetScale({3, 3, 3});
        auto mesh_body = chrono_types::make_shared<ChBody>();
        mesh_body->SetPos(offsets[i]);
        mesh_body->SetRot(Q_from_AngZ(CH_C_2PI * ChRandom()));
        mesh_body->AddAsset(trimesh_shape);
        mesh_body->SetBodyFixed(true);
        chsystem->Add(mesh_body);
    }
}

void AddTerrain(ChSystem* chsystem) {
    // add terrain with weighted textures
    auto terrain_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    terrain_mesh->LoadWavefrontMesh(demo_data_path + "/Environments/Iowa/terrain/terrain_v2.obj", false, true);
    terrain_mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));  // scale to a different size
    auto terrain_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    terrain_shape->SetMesh(terrain_mesh);
    terrain_shape->SetName("terrain");
    terrain_shape->SetStatic(true);

    auto gravel_tex = chrono_types::make_shared<ChVisualMaterial>();
    gravel_tex->SetKdTexture(demo_data_path + "/Environments/Iowa/terrain/Grass/GroundMudTireTracks001_COL_500.jpg");
    gravel_tex->SetRoughnessTexture(demo_data_path +
                                    "/Environments/Iowa/terrain/Grass/GroundMudTireTracks001_ROUGH_500.png");
    gravel_tex->SetNormalMapTexture(demo_data_path +
                                    "/Environments/Iowa/terrain/Grass/GroundMudTireTracks001_NRM_500.jpg");
    gravel_tex->SetWeightTexture(demo_data_path + "/Environments/Iowa/terrain/Terrain_Weightmap_Gravel_v3.png");
    gravel_tex->SetSpecularColor({.0f, .0f, .0f});
    gravel_tex->SetTextureScale({1000.0, 1000.0, 1.0});
    gravel_tex->SetRoughness(1.f);
    gravel_tex->SetUseSpecularWorkflow(false);
    terrain_shape->material_list.push_back(gravel_tex);

    auto grass_tex_1 = chrono_types::make_shared<ChVisualMaterial>();
    grass_tex_1->SetKdTexture(demo_data_path + "/Environments/Iowa/terrain/Grass/GroundGrassGreen001_COL_500.jpg");
    grass_tex_1->SetRoughnessTexture(demo_data_path +
                                     "/Environments/Iowa/terrain/Grass/GroundGrassGreen001_ROUGH_500.jpg");
    grass_tex_1->SetNormalMapTexture(demo_data_path +
                                     "/Environments/Iowa/terrain/Grass/GroundGrassGreen001_NRM_500.jpg");
    grass_tex_1->SetWeightTexture(demo_data_path + "/Environments/Iowa/terrain/Terrain_Weightmap_Grass_A_v3.png");
    grass_tex_1->SetTextureScale({1000.0, 1000.0, 1.0});
    grass_tex_1->SetSpecularColor({.0f, .0f, .0f});
    grass_tex_1->SetRoughness(1.f);
    grass_tex_1->SetUseSpecularWorkflow(false);
    terrain_shape->material_list.push_back(grass_tex_1);

    auto grass_tex_2 = chrono_types::make_shared<ChVisualMaterial>();
    grass_tex_2->SetKdTexture(demo_data_path +
                              "/Environments/Iowa/terrain/Grass/GroundGrassGreenPatchy002_COL_500.jpg");
    grass_tex_2->SetRoughnessTexture(demo_data_path +
                                     "/Environments/Iowa/terrain/Grass/GroundGrassGreenPatchy002_ROUGH_500.png");
    grass_tex_2->SetNormalMapTexture(demo_data_path +
                                     "/Environments/Iowa/terrain/Grass/GroundGrassGreenPatchy002_NRM_500.jpg");
    grass_tex_2->SetWeightTexture(demo_data_path + "/Environments/Iowa/terrain/Terrain_Weightmap_Grass_B_v3.png");
    grass_tex_2->SetSpecularColor({.0f, .0f, .0f});
    grass_tex_2->SetTextureScale({1000.0, 1000.0, 1.0});
    grass_tex_2->SetRoughness(1.f);
    grass_tex_2->SetUseSpecularWorkflow(false);
    terrain_shape->material_list.push_back(grass_tex_2);

    // auto field_tex = chrono_types::make_shared<ChVisualMaterial>();
    // field_tex->SetKdTexture(demo_data_path + "/Environments/Iowa/terrain/Grass/GroundMudCracked006_COL_500.jpg");
    // field_tex->SetRoughnessTexture(demo_data_path +
    //                                "/Environments/Iowa/terrain/Grass/GroundMudCracked006_ROUGH_500.png");
    // field_tex->SetNormalMapTexture(demo_data_path +
    // "/Environments/Iowa/terrain/Grass/GroundMudCracked006_NRM_500.jpg");
    // field_tex->SetWeightTexture(demo_data_path +
    // "/Environments/Iowa/terrain/Terrain_Weightmap_DirtFields_v2.png"); field_tex->SetSpecularColor({.0f, .0f,
    // .0f}); field_tex->SetTextureScale({1000.0, 1000.0, 1.0}); field_tex->SetRoughness(1.f);
    // field_tex->SetUseSpecularWorkflow(false);
    // terrain_shape->material_list.push_back(field_tex);

    auto terrain_body = chrono_types::make_shared<ChBody>();
    terrain_body->SetPos({0, 0, -.01});
    terrain_body->AddAsset(terrain_shape);
    terrain_body->SetBodyFixed(true);
    chsystem->Add(terrain_body);
}

// Wheel button callback to switch between driving modes
void customButtonCallback() {
    // We use an "anti bounce":
    static auto last_invoked = std::chrono::system_clock::now().time_since_epoch();
    auto current_invoke = std::chrono::system_clock::now().time_since_epoch();
    if (std::chrono::duration_cast<std::chrono::seconds>(current_invoke - last_invoked).count() > 3.0) {
        std::cout << "Button Callback Invoked \n";
        if (driver_mode == HUMAN) {
            driver_mode = AUTONOMOUS;
            PF_driver_ptr->Set_TheroSpeed(cur_follower_speed);
        }

        else
            driver_mode = HUMAN;
        // last, update the last call
        last_invoked = current_invoke;
    } else
        std::cout << "Callback Not Invoked, call was too close to the last one \n";
}

// Update dummy vehicles based on bazier curve and speed
void updateDummy(std::shared_ptr<ChBodyAuxRef> dummy_vehicle,
                 std::shared_ptr<ChBezierCurve> curve,
                 float dummy_speed,
                 float step_size,
                 float z_offset,
                 ChBezierCurveTracker tracker,
                 double& dummy_dist,
                 ChVector<>& dummy_prev_pos) {
    // sentinel point fo the current dummy vehicle position
    ChVector<> sen = dummy_vehicle->GetPos();
    ChVector<> prev_pos = sen;

    ChVector<> target;
    ChFrame<> frame;
    double curv;

    // obtain the tangent direction on the cloest bezier curve
    tracker.calcClosestPoint(sen, frame, curv);
    ChVector<> vel_dir = -frame.TransformDirectionLocalToParent(ChVector<>(1, 0, 0));
    // normalize velocity vector
    vel_dir.Normalize();

    // proceed vehicle by time_step*vel_dir, calculate the closest target point on the bezier curve
    tracker.calcClosestPoint(sen + (vel_dir * dummy_speed * step_size), target);

    target = target + ChVector<>(0, 0, z_offset);

    // compute angle
    float angle = atan2(vel_dir[1], vel_dir[0]);

    // finally update dummy vehicle position and rotated direction
    dummy_vehicle->SetRot(Q_from_Euler123(ChVector<>(0, 0, angle)));
    dummy_vehicle->SetPos(target);

    // update distance and previous position
    dummy_dist = dummy_dist + ((dummy_vehicle->GetPos() - dummy_prev_pos).Length()) * M_TO_MILE;
    dummy_prev_pos = dummy_vehicle->GetPos();
}

float controlFindSpeed_x_y(std::vector<float> time_vec, std::vector<float> speed_vec, float time, float default_speed) {
    int n = time_vec.size();
    int target_idx = n;
    // find out which range the 'time' argument falls into
    for (int i = 0; i < n; i++) {
        if (time <= time_vec[i]) {
            target_idx = i;
            break;
        }
    }
    if (target_idx == 0) {
        return default_speed;
    } else {
        return speed_vec[target_idx - 1];
    }
}

// dummy button call function
// this section of the code should be optimized !
// Wheel botton callback to record button press without any functionalities
void dummyButtonCallback_r_3() {
    static auto last_invoked_dummy_1 = std::chrono::system_clock::now().time_since_epoch();
    auto current_invoke_dummy_1 = std::chrono::system_clock::now().time_since_epoch();

    if (std::chrono::duration_cast<std::chrono::seconds>(current_invoke_dummy_1 - last_invoked_dummy_1).count() > 1.0) {
        std::cout << "Button dummy r_3 Callback Invoked: " << std::endl;
        time_t my_time = time(NULL);
        button_buffer << "button dummy r_3 pressed; time: ";
        button_buffer << ctime(&my_time);
    }

    last_invoked_dummy_1 = current_invoke_dummy_1;
}
