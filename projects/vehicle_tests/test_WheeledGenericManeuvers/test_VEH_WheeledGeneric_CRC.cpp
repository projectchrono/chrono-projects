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
// Authors: Radu Serban, Mike Taylor
// =============================================================================
//
// Test program for the generic vehicle running a constant radius turn
//
// The vehicle reference frame has Z up, X towards the front of the vehicle, and
// Y pointing to the left.
//
// =============================================================================

#include <vector>

#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChFilters.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"

#include "chrono_models/vehicle/generic/Generic_Vehicle.h"
#include "chrono_models/vehicle/generic/powertrain/Generic_EngineSimpleMap.h"
#include "chrono_models/vehicle/generic/powertrain/Generic_AutomaticTransmissionSimpleMap.h"
#include "chrono_models/vehicle/generic/tire/Generic_FialaTire.h"
#include "chrono_models/vehicle/generic/driveline/Generic_Driveline2WD.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"

#include "chrono_thirdparty/filesystem/path.h"

// Uncomment the following line to unconditionally disable Irrlicht support
//#undef CHRONO_IRRLICHT

#ifdef CHRONO_IRRLICHT
#include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemIrrlicht.h"
#endif

// DEBUGGING:  Uncomment the following line to print shock data
// #define DEBUG_LOG

using namespace chrono;
using namespace chrono::irrlicht;
using namespace chrono::vehicle;
using namespace chrono::vehicle::generic;

// =============================================================================

// Input file names for the path-follower driver model
std::string steering_controller_file("generic/driver/SteeringController.json");
std::string speed_controller_file("generic/driver/SpeedController.json");

// Rigid terrain dimensions
double terrainHeight = 0;
double terrainLength = 500.0;  // size in X direction
double terrainWidth = 500.0;   // size in Y direction

// Simulation step size
double step_size = 1e-4;

// Time interval between two render frames
double render_step_size = 1.0 / 50;  // FPS = 50

// Time interval between two output frames
double output_step_size = 1.0 / 1;  // once a second

// Point on chassis tracked by the camera (Irrlicht only)
ChVector3d trackPoint(0.0, 0.0, 1.75);

// Simulation length (set to a negative value to disable for Irrlicht)
double tend = 30.0;

// Output directories
const std::string out_dir = "../GENERIC_VEHICLE_CRC";
const std::string pov_dir = out_dir + "/POVRAY";

// POV-Ray output
bool povray_output = false;

// Vehicle state output (forced to true if povray output enabled)
bool state_output = true;
int filter_window_size = 20;

// =============================================================================

void CalcControlPoints(double run,
    double radius,
    int nturns,
    std::vector<ChVector3d>& points,
    std::vector<ChVector3d>& inCV,
    std::vector<ChVector3d>& outCV) {
    // Height of path
    double z = 0.1;

    // Approximate circular path using 4 points
    double direction = radius > 0 ? 1 : -1;
    radius = std::abs(radius);
    double factor = radius * 0.55191502449;


    ChVector3d P1(0, direction*radius, z);
    ChVector3d P1_in = P1 - ChVector3d(factor, 0, 0);
    ChVector3d P1_out = P1 + ChVector3d(factor, 0, 0);

    ChVector3d P2(radius, 0, z);
    ChVector3d P2_in = P2 + ChVector3d(0, direction*factor, 0);
    ChVector3d P2_out = P2 - ChVector3d(0, direction*factor, 0);

    ChVector3d P3(0, -direction*radius, z);
    ChVector3d P3_in = P3 + ChVector3d(factor, 0, 0);
    ChVector3d P3_out = P3 - ChVector3d(factor, 0, 0);

    ChVector3d P4(-radius, 0, z);
    ChVector3d P4_in = P4 - ChVector3d(0, direction*factor, 0);
    ChVector3d P4_out = P4 + ChVector3d(0, direction*factor, 0);

    // Start point
    ChVector3d P0(-run, direction*radius, z);
    ChVector3d P0_in = P0 - ChVector3d(run/2., 0, 0);
    ChVector3d P0_out = P0 + ChVector3d(run/2., 0, 0);

    points.push_back(P0);
    inCV.push_back(P0_in);
    outCV.push_back(P0_out);

    for (int i = 0; i < nturns; i++) {
        points.push_back(P1);
        inCV.push_back(P1_in);
        outCV.push_back(P1_out);

        points.push_back(P2);
        inCV.push_back(P2_in);
        outCV.push_back(P2_out);

        points.push_back(P3);
        inCV.push_back(P3_in);
        outCV.push_back(P3_out);

        points.push_back(P4);
        inCV.push_back(P4_in);
        outCV.push_back(P4_out);
    }

    points.push_back(P1);
    inCV.push_back(P1_in);
    outCV.push_back(P1_out);
}


// =============================================================================

int main(int argc, char* argv[]) {
    // Set path to Chrono and Chrono::Vehicle data directories
    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_VEHICLE_DATA_DIR);

    double initFwdSpd = 30.0 / 3.6;  // kph to m/s
    double finalFwdSpd = 100.0 / 3.6;  // kph to m/s
    double cornerRadius = 200;

    // Check for input arguments for running this test in batch
    // First argument is the initial vehicle speed in m/s
    // Second argument is the target final speed in m/s
    // Third argument is the radius of the turn in m
    if (argc > 1)
        initFwdSpd = std::atof(argv[1]);
    if (argc > 2)
        finalFwdSpd = std::atof(argv[2]);
    if (argc > 3)
        cornerRadius = std::atof(argv[4]);

    // ------------------------------------
    // Parameters for the Bezier curve path
    // ------------------------------------

    double run = 10;
    int nturns = 1 + int(std::ceil(((finalFwdSpd + initFwdSpd)/2*tend) / (cornerRadius * CH_C_2PI)));

    // Initial vehicle location
    ChVector3d initLoc(- run -5, cornerRadius, 0.6);


    // --------------------------
    // Create the various modules
    // --------------------------

    // Create the vehicle: specify if chassis is fixed, the suspension type
    // and the inital forward speed
    Generic_Vehicle vehicle(false, SuspensionTypeWV::DOUBLE_WISHBONE, SuspensionTypeWV::DOUBLE_WISHBONE,
                            SteeringTypeWV::PITMAN_ARM, DrivelineTypeWV::AWD, BrakeType::SHAFTS);
    vehicle.Initialize(ChCoordsys<>(initLoc), initFwdSpd);
    vehicle.SetChassisVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetSuspensionVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetSteeringVisualizationType(VisualizationType::PRIMITIVES);
    vehicle.SetWheelVisualizationType(VisualizationType::NONE);

    // Create the ground
    auto patch_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    patch_mat->SetFriction(0.9f);
    patch_mat->SetRestitution(0.01f);
    RigidTerrain terrain(vehicle.GetSystem());
    auto patch =
        terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector3d(0, 0, terrainHeight), QUNIT), terrainLength, terrainWidth);
    patch->SetColor(ChColor(0.5f, 0.8f, 0.5f));
    patch->SetTexture(vehicle::GetDataFile("terrain/textures/tile4.jpg"), 600, 600);
    terrain.Initialize();

    // Create and initialize the powertrain system
    auto engine = chrono_types::make_shared<Generic_EngineSimpleMap>("Engine");
    auto transmission = chrono_types::make_shared<Generic_AutomaticTransmissionSimpleMap>("Transmission");
    auto powertrain = chrono_types::make_shared<ChPowertrainAssembly>(engine, transmission);
    vehicle.InitializePowertrain(powertrain);

    // Create the tires
    auto tire_FL = chrono_types::make_shared<Generic_FialaTire>("FL");
    auto tire_FR = chrono_types::make_shared<Generic_FialaTire>("FR");
    auto tire_RL = chrono_types::make_shared<Generic_FialaTire>("RL");
    auto tire_RR = chrono_types::make_shared<Generic_FialaTire>("RR");

    vehicle.InitializeTire(tire_FL, vehicle.GetWheel(0, LEFT), VisualizationType::PRIMITIVES);
    vehicle.InitializeTire(tire_FR, vehicle.GetWheel(0, RIGHT), VisualizationType::PRIMITIVES);
    vehicle.InitializeTire(tire_RL, vehicle.GetWheel(1, LEFT), VisualizationType::PRIMITIVES);
    vehicle.InitializeTire(tire_RR, vehicle.GetWheel(1, RIGHT), VisualizationType::PRIMITIVES);

    // -------------------------------------
    // Create the path and the driver system
    // -------------------------------------

    std::vector<ChVector3d> points;
    std::vector<ChVector3d> inCV;
    std::vector<ChVector3d> outCV;
    CalcControlPoints(run, cornerRadius, nturns, points, inCV, outCV);
    auto path = chrono_types::make_shared<ChBezierCurve>(points, inCV, outCV);

    ChPathFollowerDriver driver(vehicle, vehicle::GetDataFile(steering_controller_file),
                                vehicle::GetDataFile(speed_controller_file), path, "my_path", initFwdSpd);
    driver.Initialize();

    // Report out the mass of the entire vehicle to the screen
    std::cout << "Vehicle Mass: " << vehicle.GetMass() << std::endl;

#ifdef CHRONO_IRRLICHT

    // ---------------------------------------
    // Create the vehicle Irrlicht application
    // ---------------------------------------

    auto vis = chrono_types::make_shared<ChWheeledVehicleVisualSystemIrrlicht>();
    vis->SetWindowTitle("Generic Wheeled Vehicle Constant Radius Cornering Test");
    vis->SetChaseCamera(trackPoint, 6.0, 0.5);
    vis->Initialize();
    vis->AddTypicalLights();
    vis->AddSkyBox();
    vis->AddLogo();
    vis->AttachVehicle(&vehicle);

    // Visualization of controller points (sentinel & target)
    irr::scene::IMeshSceneNode* ballS = vis->GetSceneManager()->addSphereSceneNode(0.1f);
    irr::scene::IMeshSceneNode* ballT = vis->GetSceneManager()->addSphereSceneNode(0.1f);
    ballS->getMaterial(0).EmissiveColor = irr::video::SColor(0, 255, 0, 0);
    ballT->getMaterial(0).EmissiveColor = irr::video::SColor(0, 0, 255, 0);

#endif

    // ------------------------------------
    // Prepare output directories and files
    // ------------------------------------

    state_output = state_output || povray_output;

    // Create output directories
    if (state_output) {
        if (!filesystem::create_directory(filesystem::path(out_dir))) {
            std::cout << "Error creating directory " << out_dir << std::endl;
            return 1;
        }
    }
    if (povray_output) {
        if (!filesystem::create_directory(filesystem::path(pov_dir))) {
            std::cout << "Error creating directory " << pov_dir << std::endl;
            return 1;
        }
        driver.ExportPathPovray(out_dir);
    }

    utils::CSV_writer csv("\t");
    csv.stream().setf(std::ios::scientific | std::ios::showpos);
    csv.stream().precision(6);

    utils::ChRunningAverage fwd_acc_GC_filter(filter_window_size);
    utils::ChRunningAverage lat_acc_GC_filter(filter_window_size);
    utils::ChRunningAverage vert_acc_GC_filter(filter_window_size);

    utils::ChRunningAverage fwd_acc_driver_filter(filter_window_size);
    utils::ChRunningAverage lat_acc_driver_filter(filter_window_size);
    utils::ChRunningAverage vert_acc_driver_filter(filter_window_size);

    // Driver location in vehicle local frame
    ChVector3d driver_pos = vehicle.GetChassis()->GetLocalDriverCoordsys().pos;

// ---------------
// Simulation loop
// ---------------

#ifdef DEBUG_LOG
    std::cout << "\n\n============ System Configuration ============\n";
    vehicle.LogHardpointLocations();
#endif

    // Number of simulation steps between two 3D view render frames
    int render_steps = (int)std::ceil(render_step_size / step_size);

    // Number of simulation steps between two output frames
    int output_steps = (int)std::ceil(output_step_size / step_size);

    // Initialize simulation frame counter and simulation time
    double time = 0;
    int step_number = 0;
    int render_frame = 0;

#ifdef CHRONO_IRRLICHT

    while (vis->Run()) {
        time = vehicle.GetChTime();

        // End simulation
        if ((time > tend) && (tend > 0))
            break;
#else

    while (time <= tend) {
        time = vehicle.GetChTime();

#endif

        // std::cout << vehicle.GetSystem()->GetSolverCallsCount() << std::endl;
        // Extract accelerations to add to the filter
        ChVector3d acc_CG = vehicle.GetChassisBody()->GetPos_dtdt();
        acc_CG = vehicle.GetChassisBody()->GetCoord().TransformDirectionParentToLocal(acc_CG);
        ChVector3d acc_driver = vehicle.GetPointAcceleration(driver_pos);
        double fwd_acc_CG = fwd_acc_GC_filter.Add(acc_CG.x());
        double lat_acc_CG = lat_acc_GC_filter.Add(acc_CG.y());
        double vert_acc_CG = vert_acc_GC_filter.Add(acc_CG.z());
        double fwd_acc_driver = fwd_acc_driver_filter.Add(acc_driver.x());
        double lat_acc_driver = lat_acc_driver_filter.Add(acc_driver.y());
        double vert_acc_driver = vert_acc_driver_filter.Add(acc_driver.z());

#ifdef CHRONO_IRRLICHT
        // Update sentinel and target location markers for the path-follower controller.
        // Note that we do this whether or not we are currently using the path-follower driver.
        const ChVector3d& pS = driver.GetSteeringController().GetSentinelLocation();
        const ChVector3d& pT = driver.GetSteeringController().GetTargetLocation();
        ballS->setPosition(irr::core::vector3df((irr::f32)pS.x(), (irr::f32)pS.y(), (irr::f32)pS.z()));
        ballT->setPosition(irr::core::vector3df((irr::f32)pT.x(), (irr::f32)pT.y(), (irr::f32)pT.z()));
#endif

        // Driver inputs
        DriverInputs driver_inputs = driver.GetInputs();

        // Render scene
        if (step_number % render_steps == 0) {
#ifdef CHRONO_IRRLICHT
            vis->BeginScene();
            vis->Render();
            vis->EndScene();
#endif

#ifdef DEBUG_LOG
            std::cout << "\n\n============ System Information ============\n";
            std::cout << "Time = " << time << "\n\n";
            // vehicle.DebugLog(DBG_SPRINGS | DBG_SHOCKS | DBG_CONSTRAINTS);
            vehicle.DebugLog(OUT_SPRINGS | OUT_SHOCKS | OUT_CONSTRAINTS);
#endif

            if (povray_output) {
                char filename[100];
                sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), render_frame + 1);
                utils::WriteVisualizationAssets(vehicle.GetSystem(), filename);
            }

            if (state_output) {
                ChVector3d vel_CG = vehicle.GetChassisBody()->GetPos_dt();
                vel_CG = vehicle.GetChassisBody()->GetCoord().TransformDirectionParentToLocal(vel_CG);

                ChVector3d vel_driver_abs =
                    vehicle.GetChassisBody()->GetFrame_REF_to_abs().PointSpeedLocalToParent(driver_pos);
                ChVector3d vel_driver_local =
                    vehicle.GetChassisBody()->GetFrame_REF_to_abs().TransformDirectionParentToLocal(vel_driver_abs);

                int axle = vehicle.GetDriveline()->GetDrivenAxleIndexes()[0];

                // Vehicle and Control Values
                csv << time << driver_inputs.m_steering << driver_inputs.m_throttle << driver_inputs.m_braking;
                csv << engine->GetMotorSpeed() << engine->GetOutputMotorshaftTorque();
                // Chassis Position, Velocity, & Acceleration (Unfiltered and Filtered)
                csv << vehicle.GetChassis()->GetPos().x() << vehicle.GetChassis()->GetPos().y()
                    << vehicle.GetChassis()->GetPos().z();
                csv << vel_CG.x() << vel_CG.y() << vel_CG.z();
                csv << acc_CG.x() << acc_CG.y() << acc_CG.z();
                csv << fwd_acc_CG << lat_acc_CG << vert_acc_CG;
                // Driver Position, Velocity, & Acceleration (Unfiltered and Filtered)
                csv << vehicle.GetDriverPos().x() << vehicle.GetDriverPos().y() << vehicle.GetDriverPos().z();
                csv << vel_driver_local.x() << vel_driver_local.y() << vel_driver_local.z();
                csv << acc_driver.x() << acc_driver.y() << acc_driver.z();         // Chassis CSYS
                csv << fwd_acc_driver << lat_acc_driver << vert_acc_driver;  // filtered Chassis CSYS
                // Torque to the rear wheels
                csv << vehicle.GetDriveline()->GetSpindleTorque(axle, LEFT);
                csv << vehicle.GetDriveline()->GetSpindleTorque(axle, RIGHT);
                // Tire Slip Angles
                csv << tire_FL->GetSlipAngle() << tire_FL->GetLongitudinalSlip() << tire_FL->GetCamberAngle();
                csv << tire_FR->GetSlipAngle() << tire_FR->GetLongitudinalSlip() << tire_FR->GetCamberAngle();
                csv << tire_RL->GetSlipAngle() << tire_RL->GetLongitudinalSlip() << tire_RL->GetCamberAngle();
                csv << tire_RR->GetSlipAngle() << tire_RR->GetLongitudinalSlip() << tire_RR->GetCamberAngle();
                // Suspension Lengths
                for (const auto& susp_force : vehicle.GetSuspension(0)->ReportSuspensionForce(LEFT))
                    csv << susp_force.length;
                for (const auto& susp_force : vehicle.GetSuspension(0)->ReportSuspensionForce(RIGHT))
                    csv << susp_force.length;
                for (const auto& susp_force : vehicle.GetSuspension(1)->ReportSuspensionForce(LEFT))
                    csv << susp_force.length;
                for (const auto& susp_force : vehicle.GetSuspension(1)->ReportSuspensionForce(RIGHT))
                    csv << susp_force.length;
                // tire normal forces
                csv << tire_FL->ReportTireForce(&terrain).force;
                csv << tire_FR->ReportTireForce(&terrain).force;
                csv << tire_RL->ReportTireForce(&terrain).force;
                csv << tire_RR->ReportTireForce(&terrain).force;
                //
                csv << vehicle.GetChassis()->GetRot();
                csv << std::endl;
            }

            render_frame++;
        }

        driver.Synchronize(time);
        terrain.Synchronize(time);
        vehicle.Synchronize(time, driver_inputs, terrain);
#ifdef CHRONO_IRRLICHT
        vis->Synchronize(time, driver_inputs);
#endif

        // Update for the new target vehicle speed
        driver.SetDesiredSpeed((finalFwdSpd - initFwdSpd) / tend * time +initFwdSpd);

        // Advance simulation for one timestep for all modules
        driver.Advance(step_size);
        terrain.Advance(step_size);
        vehicle.Advance(step_size);
#ifdef CHRONO_IRRLICHT
        vis->Advance(step_size);
#endif

        // Increment frame number
        step_number++;
    }
    if (state_output) {
        char filename[100];
        if (cornerRadius>0)
            sprintf(filename, "%s/output_%dmps_to_%dmps_CW_Rad%dm.dat", out_dir.c_str(), int(std::round(initFwdSpd)), int(std::round(finalFwdSpd)), int(std::round(std::abs(cornerRadius))));
        else
            sprintf(filename, "%s/output_%dmps_to_%dmps_CCW_Rad%dm.dat", out_dir.c_str(), int(std::round(initFwdSpd)), int(std::round(finalFwdSpd)), int(std::round(std::abs(cornerRadius))));
        csv.write_to_file(filename);
    }
    return 0;
}
