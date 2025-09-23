// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2023 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Rainer Gericke
// =============================================================================
//
// FEDA wall to wall turning test & calculation of some driver model parameters
//
// The vehicle reference frame has Z up, X towards the front of the vehicle, and
// Y pointing to the left.
//
// =============================================================================

#include "chrono/utils/ChFilters.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/output/ChOutputASCII.h"

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/driver/ChDataDriver.h"

#include "chrono_models/vehicle/feda/FEDA.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "chrono_vehicle/driver/ChPathFollowerDriver.h"
#include "chrono_vehicle/utils/ChVehiclePath.h"

#ifdef CHRONO_IRRLICHT
#include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemIrrlicht.h"
using namespace chrono::irrlicht;
#endif

#ifdef CHRONO_VSG
#include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemVSG.h"
using namespace chrono::vsg3d;
#endif

using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::feda;

// =============================================================================

// Run-time visualization system (IRRLICHT or VSG)
ChVisualSystem::Type vis_type = ChVisualSystem::Type::VSG;

// Initial vehicle location
ChVector3d initLoc(0, 0, 1.6);

// Visualization type for vehicle parts (PRIMITIVES, MESH, or NONE)
VisualizationType chassis_vis_type = VisualizationType::MESH;
VisualizationType suspension_vis_type = VisualizationType::PRIMITIVES;
VisualizationType steering_vis_type = VisualizationType::PRIMITIVES;
VisualizationType wheel_vis_type = VisualizationType::MESH;
VisualizationType tire_vis_type = VisualizationType::MESH;

// Collision type for chassis (PRIMITIVES, MESH, or NONE)
CollisionType chassis_collision_type = CollisionType::NONE;

// Type of engine model (SHAFTS, SIMPLE, SIMPLE_MAP)
EngineModelType engine_model = EngineModelType::SIMPLE_MAP;

// Type of transmission model (SHAFTS, SIMPLE_MAP)
TransmissionModelType transmission_model = TransmissionModelType::AUTOMATIC_SIMPLE_MAP;

// Drive type (FWD, RWD, or AWD)
DrivelineTypeWV drive_type = DrivelineTypeWV::AWD;

// Steering type (PITMAN_ARM or PITMAN_ARM_SHAFTS)
SteeringTypeWV steering_type = SteeringTypeWV::PITMAN_ARM;

// Brake type (SIMPLE or SHAFTS)
BrakeType brake_type = BrakeType::SHAFTS;

// Model tierods as bodies (true) or as distance constraints (false)
bool use_tierod_bodies = true;

// Type of tire model (RIGID, RIGID_MESH, TMEASY, FIALA, PAC89, PAC02, TMSIMPLE)
TireModelType tire_model = TireModelType::TMSIMPLE;

// Contact method
ChContactMethod contact_method = ChContactMethod::SMC;

// Simulation step sizes
double step_size = 1e-3;
double tire_step_size = 1e-3;

// Time interval between two render frames
double render_step_size = 1.0 / 50;  // FPS = 50

// Debug logging
bool debug_output = false;
double debug_step_size = 1.0 / 1;  // FPS = 1

// Simulation end time
double t_end = 0;

// --------- calculate turn radius from circular 2d-vertices ----------
double CalculateTurnDiameter(std::vector<std::pair<double, double>>& v, ChVector2d& c) {
    double d = 0.0;
    double x_min = 1.0e15;
    double x_max = -1.0e15;
    double y_min = x_min;
    double y_max = y_max;
    for (auto it = v.begin(); it != v.end(); it++) {
        if (it->first < x_min)
            x_min = it->first;
        if (it->first > x_max)
            x_max = it->first;
        if (it->second < y_min)
            y_min = it->second;
        if (it->second > y_max)
            y_max = it->second;
    }
    d = ((x_max - x_min) + (x_max - x_min)) / 2.0;
    c.x() = (x_max + x_min) / 2.0;
    c.y() = (y_max + y_min) / 2.0;
    return d;
}
// =============================================================================

int main(int argc, char* argv[]) {
    // Set path to Chrono and Chrono::Vehicle data directories
    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_VEHICLE_DATA_DIR);

    const double ft2m = 0.3048;
    const double in2m = 0.0254;

    // KRC Test Results from document 'Wall to Wall Turning Diameter_180520.ppt'
    const double leftWall2WallDistanceRef{50.8 * ft2m};
    const double rightWall2WallDistanceRef{51.1 * ft2m};

    ChFunctionInterp steeringGear;
    steeringGear.AddPoint(-1.0, -648.0);
    steeringGear.AddPoint(0.0, 0.0);
    steeringGear.AddPoint(1.0, 648.0);

    // --------------
    // Create systems
    // --------------

    // Create the FEDA vehicle, set parameters, and initialize
    std::string vehicleName = "FEDA";
    FEDA feda;
    feda.SetContactMethod(contact_method);
    feda.SetChassisCollisionType(chassis_collision_type);
    feda.SetChassisFixed(false);
    feda.SetInitPosition(ChCoordsys<>(initLoc, QUNIT));
    feda.SetEngineType(engine_model);
    feda.SetTransmissionType(transmission_model);
    feda.SetBrakeType(brake_type);
    feda.SetTireType(tire_model);
    feda.SetTireStepSize(tire_step_size);
    feda.SetRideHeight_OnRoad();
    feda.Initialize();

    if (tire_model == TireModelType::RIGID_MESH)
        tire_vis_type = VisualizationType::MESH;

    feda.SetChassisVisualizationType(chassis_vis_type);
    feda.SetSuspensionVisualizationType(suspension_vis_type);
    feda.SetSteeringVisualizationType(steering_vis_type);
    feda.SetWheelVisualizationType(wheel_vis_type);
    feda.SetTireVisualizationType(tire_vis_type);

    // Associate a collision system
    feda.GetSystem()->SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Create the terrain
    RigidTerrain terrain(feda.GetSystem());

    ChContactMaterialData minfo;
    minfo.mu = 0.8f;
    minfo.cr = 0.01f;
    minfo.Y = 2e7f;
    auto patch_mat = minfo.CreateMaterial(contact_method);

    auto patch = terrain.AddPatch(patch_mat, CSYSNORM, 200.0, 200.0);
    patch->SetTexture(vehicle::GetDataFile("terrain/textures/dirt.jpg"), 200, 200);
    patch->SetColor(ChColor(0.8f, 0.8f, 0.5f));

    terrain.Initialize();

    double vehicleWith = 90.1 * in2m;  // inch -> m
    double wheelBase = 130.0 * in2m;   // inch -> m
    double tireRadius = 0.3665;

    // for calculation of the wall to wall distance we need the front corner positions
    // of the car body in the vehicle reference system
    ChVector3d leftCornerPt{1.3 * tireRadius, vehicleWith / 2.0, 0.5};
    ChVector3d rightCornerPt{1.3 * tireRadius, -vehicleWith / 2.0, 0.5};

    // -----------------
    // Initialize output
    // -----------------

    const std::string out_dir = GetChronoOutputPath() + vehicleName;
    const std::string pov_dir = out_dir + "/POVRAY";

    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    // Initialize output file for driver inputs
    std::string driver_file = out_dir + "/driver_inputs.txt";
    utils::ChWriterCSV driver_csv(" ");

    std::string w2w_file_left = out_dir + "/w2w_left_";
    std::string w2w_file_right = out_dir + "/w2w_right_";
    if (tire_model == TireModelType::PAC02) {
        w2w_file_left.append("pac02");
        w2w_file_right.append("pac02");
    }
    if (tire_model == TireModelType::TMSIMPLE) {
        w2w_file_left.append("tmsimple");
        w2w_file_right.append("tmsimple");
    }
    if (tire_model == TireModelType::TMEASY) {
        w2w_file_left.append("tmeasy");
        w2w_file_right.append("tmeasy");
    }
    w2w_file_left.append(".txt");
    w2w_file_right.append(".txt");
    utils::ChWriterCSV w2w_left_csv(" ");
    utils::ChWriterCSV w2w_right_csv(" ");

    // Set up vehicle output
    feda.GetVehicle().SetChassisOutput(true);
    feda.GetVehicle().SetSuspensionOutput(0, true);
    feda.GetVehicle().SetSteeringOutput(0, true);
    feda.GetVehicle().SetOutput(ChOutput::Type::ASCII, ChOutput::Mode::FRAMES, out_dir, "output", 0.1);

    // Generate JSON information with available output channels
    feda.GetVehicle().ExportComponentList(out_dir + "/component_list.json");

    double setThrottle = 0.15;
    double tStart = 10.0;  // settle the vehicle
    double tSteer = 2.0;   // time of changing the steering wheel
    double tHold = 60.0;
    t_end = tStart + 2 * tSteer + 2 * tHold;
    std::vector<std::pair<double, double>> leftTrace, leftCornerTrace;
    std::vector<std::pair<double, double>> rightTrace, rightCornerTrace;
    std::vector<ChDataDriver::Entry> drSignal{{0, 0, 0, 0},
                                              {tStart, 0, setThrottle, 0},
                                              {tStart + tSteer, 1, setThrottle, 0},
                                              {tStart + tSteer + tHold, 1, setThrottle, 0},
                                              {tStart + 2 * tSteer + tHold, -1, setThrottle, 0},
                                              {tStart + 2 * tSteer + 2 * tHold, -1, setThrottle, 0}};

    // ------------------------------------------------------------------------------
    // Create the vehicle run-time visualization interface and the interactive driver
    // ------------------------------------------------------------------------------

#ifndef CHRONO_IRRLICHT
    if (vis_type == ChVisualSystem::Type::IRRLICHT)
        vis_type = ChVisualSystem::Type::VSG;
#endif
#ifndef CHRONO_VSG
    if (vis_type == ChVisualSystem::Type::VSG)
        vis_type = ChVisualSystem::Type::IRRLICHT;
#endif

    ChDataDriver driver(feda.GetVehicle(), drSignal);

    std::shared_ptr<ChVehicleVisualSystem> vis;
    switch (vis_type) {
        case ChVisualSystem::Type::IRRLICHT: {
#ifdef CHRONO_IRRLICHT
            // Create the vehicle Irrlicht interface
            auto vis_irr = chrono_types::make_shared<ChWheeledVehicleVisualSystemIrrlicht>();
            vis_irr->SetWindowTitle("FEDA Steady State Cornering Demo");
            vis_irr->SetChaseCamera(ChVector3d(0.0, 0.0, 1.75), 7.0, 0.5);
            vis_irr->Initialize();
            vis_irr->AddLightDirectional();
            vis_irr->AddSkyBox();
            vis_irr->AddLogo();
            vis_irr->AttachVehicle(&feda.GetVehicle());
            vis = vis_irr;
#endif
            break;
        }
        default:
        case ChVisualSystem::Type::VSG: {
#ifdef CHRONO_VSG
            // Create the vehicle VSG interface
            auto vis_vsg = chrono_types::make_shared<ChWheeledVehicleVisualSystemVSG>();
            vis_vsg->SetWindowTitle(vehicleName + " Wall To Wall Turning Test");
            vis_vsg->AttachVehicle(&feda.GetVehicle());
            vis_vsg->SetChaseCamera(ChVector3d(0.0, 0.0, 1.75), 9.0, 0.5);
            vis_vsg->SetWindowSize(ChVector2i(1200, 800));
            vis_vsg->SetWindowPosition(ChVector2i(100, 300));
            vis_vsg->EnableSkyBox(true);
            vis_vsg->SetCameraAngleDeg(40);
            vis_vsg->SetLightIntensity(1.0f);
            vis_vsg->SetLightDirection(1.5 * CH_PI_2, CH_PI_4);
            vis_vsg->EnableShadows(true);
            vis = vis_vsg;
#endif
            break;
        }
    }

    int sentinelID = 0;
    int targetID = 0;
    if (vis) {
        vis->Initialize();
    }

    // ---------------
    // Simulation loop
    // ---------------

    feda.GetVehicle().LogSubsystemTypes();

    if (debug_output) {
        std::cout << "\n\n============ System Configuration ============\n";
        feda.LogHardpointLocations();
    }

    // Number of simulation steps between miscellaneous events
    int render_steps = (int)std::ceil(render_step_size / step_size);
    int debug_steps = (int)std::ceil(debug_step_size / step_size);

    // Initialize simulation frame counters
    int step_number = 0;
    int render_frame = 0;

    feda.GetVehicle().EnableRealtime(true);

    double real_speed = feda.GetVehicle().GetChassis()->GetSpeed();

    while (true) {
        double time = feda.GetSystem()->GetChTime();
        real_speed = feda.GetVehicle().GetSpeed();
        double real_throttle = driver.GetThrottle();

        double x = feda.GetVehicle().GetPos().x();
        double y = feda.GetVehicle().GetPos().y();
        if (time >= tStart + tSteer && time <= tStart + tSteer + tHold) {
            leftTrace.push_back({x, y});
            // at left turn the right front vehicle corner is the outmost position
            double xl = feda.GetChassis()->GetPointLocation(rightCornerPt).x();
            double yl = feda.GetChassis()->GetPointLocation(rightCornerPt).y();
            leftCornerTrace.push_back({xl, yl});
            w2w_left_csv << x << y << xl << yl << std::endl;
        }
        if (time >= tStart + 2 * tSteer + tHold && time <= t_end) {
            rightTrace.push_back({x, y});
            // at left turn the right front vehicle corner is the outmost position
            double xr = feda.GetChassis()->GetPointLocation(leftCornerPt).x();
            double yr = feda.GetChassis()->GetPointLocation(leftCornerPt).y();
            rightCornerTrace.push_back({xr, yr});
            w2w_right_csv << x << y << xr << yr << std::endl;
        }
        // w2w_left_csv << time << steeringGear.GetVal(driver.GetSteering()) << real_speed << real_accy2 << std::endl;

        // End simulation
        if (time >= t_end) {
            std::cout << "Manoever ended because max. time " << t_end << " s is reached.\n";
            break;
        }

        // Render scene and output POV-Ray data
        if (vis) {
            if (!vis->Run())
                break;

            if (step_number % render_steps == 0) {
                vis->BeginScene();
                vis->Render();
                vis->EndScene();

                render_frame++;
            }
        }

        // Debug logging
        if (debug_output && step_number % debug_steps == 0) {
            std::cout << "\n\n============ System Information ============\n";
            std::cout << "Time = " << time << "\n\n";
            feda.DebugLog(OUT_SPRINGS | OUT_SHOCKS | OUT_CONSTRAINTS);

            auto marker_driver = feda.GetChassis()->GetMarkers()[0]->GetAbsCoordsys().pos;
            auto marker_com = feda.GetChassis()->GetMarkers()[1]->GetAbsCoordsys().pos;
            std::cout << "Markers\n";
            std::cout << "  Driver loc:      " << marker_driver.x() << " " << marker_driver.y() << " "
                      << marker_driver.z() << std::endl;
            std::cout << "  Chassis COM loc: " << marker_com.x() << " " << marker_com.y() << " " << marker_com.z()
                      << std::endl;
        }

        // Driver inputs
        DriverInputs driver_inputs = driver.GetInputs();

        // Update modules (process inputs from other modules)
        driver.Synchronize(time);
        terrain.Synchronize(time);
        feda.Synchronize(time, driver_inputs, terrain);
        if (vis)
            vis->Synchronize(time, driver_inputs);

        // Advance simulation for one timestep for all modules
        driver.Advance(step_size);
        terrain.Advance(step_size);
        feda.Advance(step_size);
        if (vis)
            vis->Advance(step_size);

        // Increment frame number
        step_number++;
    }

    // w2w_left_csv.WriteToFile(w2w_file_left);
    // w2w_right_csv.WriteToFile(w2w_file_right);

    ChVector2d leftCenter, rightCenter, dum1, dum2;
    double leftCenterTurnDiameter = CalculateTurnDiameter(leftTrace, leftCenter);
    double rightCenterTurnDiameter = CalculateTurnDiameter(rightTrace, rightCenter);

    double leftWall2WallDistance = CalculateTurnDiameter(leftCornerTrace, dum1);
    double rightWall2WallDistance = CalculateTurnDiameter(rightCornerTrace, dum2);

    // compare to reference values
    double errLeft = 100.0 * (leftWall2WallDistanceRef - leftWall2WallDistance) / leftWall2WallDistanceRef;
    double errRight = 100.0 * (rightWall2WallDistanceRef - rightWall2WallDistance) / rightWall2WallDistanceRef;

    double wheelAngleLeft = asin(2.0 * wheelBase / leftCenterTurnDiameter);
    double wheelAngleRight = asin(2.0 * wheelBase / rightCenterTurnDiameter);
    double wheelAngleAvg = (wheelAngleLeft + wheelAngleRight) / 2.0;
    std::cout << "Wall-to-Wall Test Results for FEDA:" << std::endl;
    std::cout << "Centerline:" << std::endl;
    std::cout << "  Left Turn Diameter  = " << leftCenterTurnDiameter << "m" << std::endl;
    std::cout << "  Avg. Wheel Turn Angle  = " << wheelAngleLeft * CH_RAD_TO_DEG << "deg (" << wheelAngleLeft << "rad)"
              << std::endl;
    std::cout << "  Right Turn Diameter = " << rightCenterTurnDiameter << "m" << std::endl;
    std::cout << "  Avg. Wheel Turn Angle  = " << wheelAngleRight * CH_RAD_TO_DEG << "deg (" << wheelAngleRight
              << "rad)" << std::endl;
    std::cout << "Vehicle Front Corner Trace:" << std::endl;
    std::cout << "  Left Turn Wall to Wall Distance  = " << leftWall2WallDistance << "m  ("
              << leftWall2WallDistance / ft2m << "ft)" << std::endl;
    std::cout << "  Right Turn Wall to Wall Distance  = " << rightWall2WallDistance << "m  ("
              << rightWall2WallDistance / ft2m << "ft)" << std::endl;
    std::cout << "Vehicle Parameters to use with driver models:" << std::endl;
    std::cout << "  Vehicle wheel base = " << wheelBase << "m" << std::endl;
    std::cout << "  Vehicle minimum turn radius = " << std::max(leftCenterTurnDiameter, rightCenterTurnDiameter) / 2.0
              << "m" << std::endl;
    std::cout << "  Average maximum wheel turn angle = " << wheelAngleAvg << "rad  (" << wheelAngleAvg * CH_RAD_TO_DEG
              << "deg)" << std::endl;
    std::string pltName = out_dir + "/" + vehicleName + "_wall2wall_plot.gpl";
    std::ofstream plt(pltName);
    plt << "$LeftCorner << EOD" << std::endl;
    for (auto it = leftCornerTrace.begin(); it != leftCornerTrace.end(); it++) {
        plt << it->first << "  " << it->second << std::endl;
    }
    plt << "EOD" << std::endl;
    plt << "$RightCorner << EOD" << std::endl;
    for (auto it = rightCornerTrace.begin(); it != rightCornerTrace.end(); it++) {
        plt << it->first << "  " << it->second << std::endl;
    }
    plt << "EOD" << std::endl;
    plt << "set title '" << vehicleName << " Wall To Wall Test Results" << std::endl;
    plt << "set xlabel 'X (m)'" << std::endl;
    plt << "set ylabel 'Y (m)'" << std::endl;
    char lDist[30], rDist[30];
    snprintf(lDist, 29, "Left Turn D = %.1fm", leftWall2WallDistance);
    snprintf(rDist, 29, "Right Turn D = %.1fm", rightWall2WallDistance);
    plt << "set label '" << lDist << "' at " << leftCenter.x() << "," << leftCenter.y() << " center tc \"#00A000\""
        << std::endl;
    plt << "set label '" << rDist << "' at " << rightCenter.x() << "," << rightCenter.y() << " center tc \"#0000A0\""
        << std::endl;
    plt << "plot $LeftCorner t 'Left Turn Corner Trace' with lines lc \"#00A000\", \\" << std::endl;
    plt << "$RightCorner t 'Right Turn Corner Trace' with lines lc \"#0000A0\"" << std::endl;
    plt << "pause -1" << std::endl;

    plt.close();

    std::cout << "Deviation from KRC test results:" << std::endl;
    std::cout << "  Left turn distance error  = " << errLeft << "%" << std::endl;
    std::cout << "  Right turn distance error = " << errRight << "%" << std::endl;

    std::string sysCmd = "gnuplot " + pltName;
    system(sysCmd.c_str());
    return 0;
}
