// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2016 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Hammad Mazhar
// =============================================================================

#include <cmath>
#include <filesystem>

#include "chrono/ChConfig.h"
#include "chrono/core/ChDataPath.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"

#include "chrono_models/vehicle/hmmwv/HMMWV.h"

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleDataPath.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"

#include "chrono_vehicle/wheeled_vehicle/suspension/ChDoubleWishbone.h"

// Chrono::Multicore header files
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"
#include "input_output.h"

// Control use of VSG run-time rendering
// Note: CHRONO_VSG is defined in ChConfig.h
#undef CHRONO_VSG

#ifdef CHRONO_VSG
    #include "chrono_vsg/ChVisualSystemVSG.h"
#endif

using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::hmmwv;
using namespace chrono::utils;

// =============================================================================
enum SimModes { TURN, LANE_CHANGE, GRAVEL };

SimModes simulation_mode = GRAVEL;
// Initial vehicle location and orientation
ChVector3d initLoc(0, 0, .6);
ChQuaterniond initRot(1, 0, 0, 0);
// ChQuaterniond initRot(0.866025, 0, 0, 0.5);
// ChQuaterniond initRot(0.7071068, 0, 0, 0.7071068);
// ChQuaterniond initRot(0.25882, 0, 0, 0.965926);
// ChQuaterniond initRot(0, 0, 0, 1);

double mph_to_m_s = 0.44704;
double target_speed = 35 * mph_to_m_s;  // 35 mph

// Visualization type for vehicle parts (PRIMITIVES, MESH, or NONE)
VisualizationType chassis_vis_type = VisualizationType::NONE;
VisualizationType suspension_vis_type = VisualizationType::PRIMITIVES;
VisualizationType steering_vis_type = VisualizationType::PRIMITIVES;
VisualizationType wheel_vis_type = VisualizationType::PRIMITIVES;
VisualizationType tire_vis_type = VisualizationType::PRIMITIVES;
// Type of engine model (SHAFTS, SIMPLE, SIMPLE_MAP)
EngineModelType engine_model = EngineModelType::SHAFTS;

// Type of transmission model (SHAFTS, SIMPLE_MAP)
TransmissionModelType transmission_model = TransmissionModelType::AUTOMATIC_SHAFTS;

// Drive type (FWD, RWD, or AWD)
DrivelineTypeWV drive_type = DrivelineTypeWV::AWD;

// Type of tire model (RIGID, RIGID_MESH, FIALA)
TireModelType tire_model = TireModelType::RIGID;

// Point on chassis tracked by the camera
ChVector3d trackPoint(0.0, 0.0, 1.75);

// Contact method
ChContactMethod contact_method = ChContactMethod::NSC;
bool contact_vis = false;

// Simulation step sizes
double time_step = 1e-3;
double tire_step_size = time_step;

double tolerance = 0.00001;

int max_iteration_bilateral = 1000;  // 1000;
int max_iteration_normal = 0;
int max_iteration_sliding = 50;  // 2000;
int max_iteration_spinning = 0;

float contact_recovery_speed = 40;

int threads = 20;

// Simulation end time
double time_end = 20;
int out_fps = 60;

// Output directories
const std::string out_dir = "../HMMWV";

// fluid settings

double rho = 1000;
double fluid_r = 0.016;

real container_friction = 1.0;

// Driver paths and controllers (the circle path is not part of the Chrono vehicle data)
std::string circle_file("humvee_input/circlepath_100ft_radius.txt");
std::string lane_change_file("paths/ISO_double_lane_change.txt");
std::string gravel_driver_file("paths/straightOrigin.txt");
std::string steering_controller_file("hmmwv/SteeringController.json");
std::string speed_controller_file("hmmwv/SpeedController.json");

std::string data_output_path = "";

double hdimX = 4;
double hdimY = 1.5;
double hdimZ = 0.1;
double hthick = 0.1;

int Id_g = 100;
double r_g = 0.02;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector3d inertia_g = 0.4 * mass_g * r_g * r_g * ChVector3d(1, 1, 1);

float mu_g = 1.0f;
std::vector<real3> forces;
std::vector<real3> torques;

void static WriteVehicleData(hmmwv::HMMWV_Full& my_hmmwv, DriverInputs driver_inputs, std::vector<real3> forces, std::vector<real3> torques, std::string filename) {
    CSVGen csv_output;
    csv_output.OpenFile(filename.c_str(), false);

    auto m_driveline = my_hmmwv.GetVehicle().GetDriveline();
    auto m_suspension_front = std::dynamic_pointer_cast<ChDoubleWishbone>(my_hmmwv.GetVehicle().GetSuspension(0));
    auto m_suspension_back = std::dynamic_pointer_cast<ChDoubleWishbone>(my_hmmwv.GetVehicle().GetSuspension(1));

    csv_output << my_hmmwv.GetVehicle().GetChassis()->GetPos();
    csv_output << my_hmmwv.GetVehicle().GetSpeed();
    csv_output << m_driveline->GetOutputDriveshaftSpeed();
    csv_output << my_hmmwv.GetVehicle().GetEngine()->GetOutputMotorshaftTorque();
    csv_output << my_hmmwv.GetVehicle().GetEngine()->GetMotorSpeed();
    csv_output << my_hmmwv.GetVehicle().GetTransmission()->GetOutputDriveshaftTorque();

    csv_output << driver_inputs.m_throttle;
    csv_output << driver_inputs.m_braking;
    csv_output << Length(forces[0]) + Length(forces[1]) + Length(forces[2]) + Length(forces[3]) + Length(forces[4]);
    csv_output << Length(torques[0]) + Length(torques[1]) + Length(torques[2]) + Length(torques[3]) +
                      Length(torques[4]);
    csv_output << Length(forces[5]) + Length(forces[6]) + Length(forces[7]) + Length(forces[8]) + Length(forces[9]);
    csv_output << Length(torques[5]) + Length(torques[6]) + Length(torques[7]) + Length(torques[8]) +
                      Length(torques[9]);

    csv_output << m_driveline->GetSpindleTorque(0, LEFT);
    csv_output << m_driveline->GetSpindleTorque(0, RIGHT);
    csv_output << m_driveline->GetSpindleTorque(1, LEFT);
    csv_output << m_driveline->GetSpindleTorque(1, RIGHT);

    csv_output << my_hmmwv.GetVehicle().GetSpindleLinVel(0, LEFT);
    csv_output << my_hmmwv.GetVehicle().GetSpindleLinVel(0, RIGHT);
    csv_output << my_hmmwv.GetVehicle().GetSpindleLinVel(1, LEFT);
    csv_output << my_hmmwv.GetVehicle().GetSpindleLinVel(1, RIGHT);

    csv_output << my_hmmwv.GetVehicle().GetSpindleAngVel(0, LEFT);
    csv_output << my_hmmwv.GetVehicle().GetSpindleAngVel(0, RIGHT);
    csv_output << my_hmmwv.GetVehicle().GetSpindleAngVel(1, LEFT);
    csv_output << my_hmmwv.GetVehicle().GetSpindleAngVel(1, RIGHT);

    csv_output << m_suspension_front->GetSpringDeformation(LEFT);
    csv_output << m_suspension_front->GetSpringDeformation(RIGHT);
    csv_output << m_suspension_back->GetSpringDeformation(LEFT);
    csv_output << m_suspension_back->GetSpringDeformation(RIGHT);
    csv_output << m_suspension_front->GetShockLength(LEFT);
    csv_output << m_suspension_front->GetShockLength(RIGHT);
    csv_output << m_suspension_back->GetShockLength(LEFT);
    csv_output << m_suspension_back->GetShockLength(RIGHT);

    csv_output << forces[0];
    csv_output << forces[1];
    csv_output << forces[2];
    csv_output << forces[3];
    csv_output << forces[4];
    csv_output << torques[0];
    csv_output << torques[1];
    csv_output << torques[2];
    csv_output << torques[3];
    csv_output << torques[4];
    csv_output << forces[5];
    csv_output << forces[6];
    csv_output << forces[7];
    csv_output << forces[8];
    csv_output << forces[9];
    csv_output << torques[5];
    csv_output << torques[6];
    csv_output << torques[7];
    csv_output << torques[8];
    csv_output << torques[9];
    csv_output.endline();
    csv_output.CloseFile();
}

void CreateFluid(ChSystemMulticoreNSC* system) {
    auto fluid_container = chrono_types::make_shared<ChFluidContainer>();
    system->Add3DOFContainer(fluid_container);

    fluid_container->kernel_radius = fluid_r;
    fluid_container->collision_envelope = 0;
    fluid_container->contact_recovery_speed = 20;
    fluid_container->max_velocity = 20;
    fluid_container->contact_cohesion = 0;
    fluid_container->contact_mu = 0;

    real alpha = 0.001;
    fluid_container->alpha = 0.001;
    fluid_container->mass = 1;
    fluid_container->contact_compliance = 1e-9;
    fluid_container->epsilon = 1e-8;
    fluid_container->rho = 1000;
    fluid_container->tau = time_step * 2;
    fluid_container->viscosity = 1;  // Viscosity of water.
    fluid_container->enable_viscosity = false;
    fluid_container->artificial_pressure = false;
    fluid_container->artificial_pressure_k = .01;
    fluid_container->artificial_pressure_dq = .2 * fluid_container->kernel_radius;
    fluid_container->artificial_pressure_n = 4;

    real dist = fluid_container->kernel_radius * .9;

    utils::ChHCPSampler<> sampler(dist);
    utils::ChGenerator::PointVector points = sampler.SampleCylinderY(initLoc + ChVector3d(-1.9306 + .24, .025, 1.2), .4, .5);
    // sampler.SampleSphere(initLoc + ChVector<>(-1.93344,0 , 1.29964), .4);
    std::vector<real3> pos_fluid;
    std::vector<real3> vel_fluid;
    for (int i = 0; i < points.size(); i++) {
        pos_fluid.push_back(real3(points[i].x(), points[i].y(), points[i].z()));
        vel_fluid.push_back(real3(0));
    }
    fluid_container->AddBodies(pos_fluid, vel_fluid);

    real vol = dist * dist * dist * .8;
    real mass = 1000 * vol;
    fluid_container->mass = mass;
}
std::shared_ptr<ChBody> bottom_plate;
// =============================================================================
void CreateBase(ChSystemMulticoreNSC* system) {
    ChVector3d c_pos = ChVector3d(0, 0, 0);

    bottom_plate = chrono_types::make_shared<ChBody>();
    bottom_plate->SetMass(1);
    bottom_plate->SetPos(ChVector3d(0, 0, 0) + c_pos);
    bottom_plate->SetRot(QUNIT);
    bottom_plate->EnableCollision(true);
    bottom_plate->SetFixed(true);

    auto material = chrono_types::make_shared<ChContactMaterialNSC>();
    material->SetFriction(container_friction);
    material->SetCompliance(1e-9);
    material->SetCohesion(0);

    // Box dimensions below are half-lengths
    if (simulation_mode == TURN) {
        AddBoxGeometry(bottom_plate.get(), material, 2.0 * ChVector3d(4, 4, .1), ChVector3d(0, 0, 0), QUNIT);
    } else if (simulation_mode == LANE_CHANGE) {
        AddBoxGeometry(bottom_plate.get(), material, 2.0 * ChVector3d(4, 3, .1), ChVector3d(0, 0, 0), QUNIT);
    } else if (simulation_mode == GRAVEL) {
        AddBoxGeometry(bottom_plate.get(), material, 2.0 * ChVector3d(4, 2, .1), ChVector3d(0, 0, 0), QUNIT);

        AddBoxGeometry(bottom_plate.get(), material, 2.0 * ChVector3d(10, 2, .1), ChVector3d(4 + hdimX * 2 + 10.0, 0, 0), QUNIT);

        // extra side plates
        AddBoxGeometry(bottom_plate.get(), material, 2.0 * ChVector3d(12.5, .1, 2), ChVector3d(4 + hdimX, -hdimY - .1, 2 - .1), QUNIT);
        AddBoxGeometry(bottom_plate.get(), material, 2.0 * ChVector3d(12.5, .1, 2), ChVector3d(4 + hdimX, hdimY + .1, 2 - .1), QUNIT);
        // Top
        AddBoxGeometry(bottom_plate.get(), material, 2.0 * ChVector3d(12.5, 2 + .1 * 2, .1), ChVector3d(4 + hdimX, 0, 4 - .1), QUNIT);
    }

    bottom_plate->GetCollisionModel()->SetFamily(2);
    bottom_plate->GetCollisionModel()->DisallowCollisionsWith(2);
    system->AddBody(bottom_plate);

    if (simulation_mode == GRAVEL) {
        CreateBoxContainer(system, material, 2.0 * ChVector3d(hdimX, hdimY, hdimZ), 2.0 * hthick, ChVector3d(4 + hdimX, 0, -hdimZ * 2 - hthick));
    }
}

// Volume of a bisphere (two overlapping spheres with centers at distance c_dist along the Y axis).
double CalcBiSphereVolume(double radius, double c_dist) {
    double delta = 2 * radius - c_dist;
    double cos_theta = (radius - 0.5 * delta) / radius;
    return (4.0 / 3.0) * CH_PI * radius * radius * radius * (1 + cos_theta);
}

// Principal radii of gyration (squared) of a bisphere.
ChVector3d CalcBiSphereGyration(double radius, double c_dist) {
    double delta = 2 * radius - c_dist;
    double cos_theta = (radius - 0.5 * delta) / radius;
    double z_prim = radius - 0.5 * delta;

    double comp1 = 0.4 * radius * radius * (1 + cos_theta);
    double comp2 = -0.2 * radius * radius * (1. / 3. * (-cos_theta * cos_theta * cos_theta - 1) + (1 + cos_theta));
    double comp3 = 2. / 3. * z_prim * z_prim * (1 + cos_theta);
    double comp4 = 0.5 * radius * z_prim * std::sqrt(1 - cos_theta * cos_theta);
    double numerator = 2 * (comp1 + comp2 + comp3 + comp4);
    double denominator = 4. / 3. * (1 + cos_theta);
    double Jxx = numerator / denominator;
    double Jyy = 0.6 * radius * radius * (1. / 3. * (-cos_theta * cos_theta * cos_theta - 1) + (1 + cos_theta)) / (1 + cos_theta);

    return ChVector3d(Jxx, Jyy, Jxx);
}

// Create bisphere gravel particles (two spheres of radius r_g, with centers at a distance r_g).
void CreateGravel(ChSystemMulticoreNSC* system) {
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(mu_g);

    double mass = rho_g * CalcBiSphereVolume(r_g, r_g);
    ChVector3d inertia = mass * CalcBiSphereGyration(r_g, r_g);

    double r = 1.01 * r_g;
    utils::ChPDSampler<double> sampler(3 * r);

    // Create particles in layers until reaching the desired number of particles
    ChVector3d hdims(hdimX - r * 3, hdimY - r * 3, 0);
    ChVector3d center(4 + hdimX, 0, -.2);

    int num_particles = 0;
    while (num_particles < 100000) {
        for (const auto& pos : sampler.SampleBox(center, hdims)) {
            auto body = chrono_types::make_shared<ChBody>();
            body->SetTag(Id_g + num_particles++);
            body->SetPos(pos);
            body->SetMass(mass);
            body->SetInertiaXX(inertia);
            utils::AddBiSphereGeometry(body.get(), mat_g, r_g, r_g);
            body->EnableCollision(true);
            system->AddBody(body);
        }
        center.z() += 2 * r;
    }

    std::cout << "Created " << num_particles << " particles." << std::endl;
}

int main(int argc, char* argv[]) {
    SetChronoDataPath(CHRONO_DATA_DIR);
    SetVehicleDataPath(CHRONO_VEHICLE_DATA_DIR);

    int target_speed_mph = 35;
    if (argc == 3) {
        simulation_mode = SimModes(atoi(argv[1]));
        target_speed_mph = atoi(argv[2]);
    }
    target_speed = target_speed_mph * mph_to_m_s;

    // --------------
    // Create systems
    // --------------
    ChSystemMulticoreNSC* system = new ChSystemMulticoreNSC();
    system->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);
    system->SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));
    system->SetNumThreads(threads);

    // ---------------------
    // Edit system settings.
    // ---------------------
    system->GetSettings()->solver.tolerance = tolerance;
    system->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    system->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    system->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    system->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    system->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    system->GetSettings()->solver.compute_N = false;
    system->GetSettings()->solver.alpha = 0;
    system->GetSettings()->solver.cache_step_length = true;
    system->GetSettings()->solver.use_full_inertia_tensor = false;
    system->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    system->GetSettings()->solver.bilateral_clamp_speed = 1e8;
    system->ChangeSolverType(SolverType::BB);

    system->GetSettings()->collision.collision_envelope = 0.1 * fluid_r;

    system->GetSettings()->collision.bins_per_axis = vec3(300, 40, 50);
    system->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    std::shared_ptr<ChBezierCurve> path;

    if (simulation_mode == TURN) {
        if (!std::filesystem::exists(circle_file)) {
            std::cout << "Missing path file " << circle_file << std::endl;
            return 1;
        }
        initLoc = ChVector3d(0, 30.48, .6);
        path = ChBezierCurve::Read(circle_file);
        data_output_path = "humvee_tank_turn_" + std::to_string(target_speed_mph) + "/";
        system->GetSettings()->collision.bins_per_axis = vec3(100, 40, 50);
    } else if (simulation_mode == LANE_CHANGE) {
        initLoc = ChVector3d(-50, -125, .6);
        path = ChBezierCurve::Read(GetVehicleDataFile(lane_change_file));
        data_output_path = "humvee_tank_dlc_" + std::to_string(target_speed_mph) + "/";
        system->GetSettings()->collision.bins_per_axis = vec3(100, 40, 50);
    } else if (simulation_mode == GRAVEL) {
        //		target_speed = 15 * mph_to_m_s;
        path = ChBezierCurve::Read(GetVehicleDataFile(gravel_driver_file));
        data_output_path = "humvee_tank_gravel_" + std::to_string(target_speed_mph) + "/";
        system->GetSettings()->collision.bins_per_axis = vec3(300, 40, 50);
    }

    std::cout << "writing data to: " << data_output_path << std::endl;
    if (!CreateOutputDirectory(data_output_path)) {
        std::cout << "Error creating directory " << data_output_path << std::endl;
        return 1;
    }
    CreateBase(system);
    CreateFluid(system);

    // Create the HMMWV vehicle, set parameters, and initialize
    HMMWV_Full my_hmmwv(system);
    my_hmmwv.SetChassisFixed(false);
    my_hmmwv.SetInitPosition(ChCoordsys<>(initLoc, initRot));
    my_hmmwv.SetEngineType(engine_model);
    my_hmmwv.SetTransmissionType(transmission_model);
    my_hmmwv.SetDriveType(drive_type);
    my_hmmwv.SetTireType(tire_model);
    my_hmmwv.SetTireStepSize(tire_step_size);

    my_hmmwv.Initialize();

    my_hmmwv.SetChassisVisualizationType(chassis_vis_type);
    my_hmmwv.SetSuspensionVisualizationType(suspension_vis_type);
    my_hmmwv.SetSteeringVisualizationType(steering_vis_type);
    my_hmmwv.SetWheelVisualizationType(wheel_vis_type);
    my_hmmwv.SetTireVisualizationType(tire_vis_type);

    // The chassis has no collision shapes of its own (chassis collision type NONE). Attach the water tank to the
    // chassis. Collision shapes added before the first step are included when the collision system is initialized.
    std::shared_ptr<ChBodyAuxRef> m_chassis = my_hmmwv.GetChassisBody();

    auto material = chrono_types::make_shared<ChContactMaterialNSC>();

    // add sphere to center to mark chassis pos
    utils::AddSphereGeometry(m_chassis.get(), material, 0.115, ChVector3d(0, 0, 0), QUNIT);

    ChVector3d hdim(.1, .6, .05);

    for (int j = 0; j < 30; j++) {
        ChQuaterniond qq = QuatFromAngleY(j * 12 * CH_DEG_TO_RAD);
        ChVector3d ppp = ChVector3d(sin(12 * CH_DEG_TO_RAD * j), 0, cos(12 * CH_DEG_TO_RAD * j)) * .6;

        utils::AddBoxGeometry(m_chassis.get(), material, 2.0 * hdim, ppp + ChVector3d(-1.9306 + .24, 0, 1.3011), qq);
    }
    // End caps: cylinders of half-length .05 along the Y axis
    utils::AddCylinderGeometry(m_chassis.get(), material, .6, ChVector3d(-1.9306 + .24, -.65, 1.3011), ChVector3d(-1.9306 + .24, -.55, 1.3011));
    utils::AddCylinderGeometry(m_chassis.get(), material, .6, ChVector3d(-1.9306 + .24, .55, 1.3011), ChVector3d(-1.9306 + .24, .65, 1.3011));

    m_chassis->GetCollisionModel()->SetFamily(VehicleCollisionFamily::CHASSIS_FAMILY);
    m_chassis->EnableCollision(true);

    RigidTerrain terrain(my_hmmwv.GetSystem());

    if (simulation_mode == GRAVEL) {
        CreateGravel(system);
    }

    // ------------------------
    // Create the driver system
    // ------------------------

    ChPathFollowerDriver driver(my_hmmwv.GetVehicle(), GetVehicleDataFile(steering_controller_file), GetVehicleDataFile(speed_controller_file), path, "my_path", target_speed);
    driver.Initialize();

    // ---------------
    // Simulation loop
    // ---------------

    double time = 0;
    double exec_time = 0;

    int sim_frame = 0, out_frame = 0, next_out_frame = 0;

    int out_steps = std::ceil((1.0 / time_step) / out_fps);

#ifdef CHRONO_VSG
    // Create the run-time visualization system
    auto vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
    vis->AttachSystem(system);
    vis->SetWindowTitle("Humvee");
    vis->SetWindowSize(1280, 720);
    vis->AddCamera(ChVector3d(0, -10, 0) + initLoc, initLoc);
    vis->SetCameraVertical(CameraVerticalDir::Z);
    vis->Initialize();
#endif

    ChVector3d driver_pos = my_hmmwv.GetVehicle().GetChassis()->GetLocalDriverCoordsys().pos;

    while (time < time_end) {
        if (simulation_mode != GRAVEL) {
            bottom_plate->SetPos(ChVector3d(my_hmmwv.GetChassis()->GetBody()->GetPos().x(), my_hmmwv.GetChassis()->GetBody()->GetPos().y(), 0));
        }

        DriverInputs driver_inputs = driver.GetInputs();

        ChVector3d com_pos = my_hmmwv.GetChassisBody()->GetPos();
        printf("Chassis: %f %f %f\n", com_pos.x(), com_pos.y(), com_pos.z());
        if (simulation_mode == LANE_CHANGE) {
            if (com_pos.x() > 30) {
                driver_inputs.m_throttle = 0;
                driver_inputs.m_braking = 1.0;
                driver_inputs.m_steering = 0;
            }
        } else if (simulation_mode == GRAVEL) {
            if (com_pos.x() > 4 + 4 * 2 + 3) {
                driver_inputs.m_throttle = 0;
                driver_inputs.m_braking = 1.0;
                driver_inputs.m_steering = 0;
            }
        }

        driver.Synchronize(time);
        my_hmmwv.Synchronize(time, driver_inputs, terrain);

        driver.Advance(time_step);

#ifdef CHRONO_VSG
        if (!vis->Run())
            break;
        vis->Render();
#endif

        system->DoStepDynamics(time_step);

        forces.resize(10);
        torques.resize(10);
        std::fill(forces.begin(), forces.end(), real3(0));
        std::fill(torques.begin(), torques.end(), real3(0));

        if (sim_frame == next_out_frame) {
            std::cout << "write: " << out_frame << std::endl;
            DumpFluidData(system, data_output_path + "data_" + std::to_string(out_frame) + ".dat", true);

            DumpAllObjectsWithGeometryPovray(system, data_output_path + "vehicle_" + std::to_string(out_frame) + ".dat",
                                             true);
            WriteVehicleData(my_hmmwv, driver_inputs, forces, torques,
                             data_output_path + "stats_" + std::to_string(out_frame) + ".dat");

            out_frame++;
            next_out_frame += out_steps;
        }

        time += time_step;
        sim_frame++;
        exec_time += system->GetTimerStep();
    }

    FinishFluidDataOutput();

    cout << "==================================" << endl;
    cout << "Simulation time:   " << exec_time << endl;
    return 0;
}
