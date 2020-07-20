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

#include "chrono/core/ChRealtimeStep.h"
#include "chrono/core/ChStream.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_models/vehicle/hmmwv/HMMWV.h"

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"

#include "chrono_vehicle/wheeled_vehicle/suspension/ChDoubleWishbone.h"

// Chrono::Parallel header files
#include "chrono_parallel/collision/ChCollisionSystemParallel.h"
#include "chrono_parallel/collision/ChNarrowphaseRUtils.h"
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"
#include "input_output.h"

#undef CHRONO_OPENGL
#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;
using namespace chrono::vehicle;
using namespace chrono::vehicle::hmmwv;
using namespace chrono::utils;

// =============================================================================
enum SimModes { TURN, LANE_CHANGE, GRAVEL };

SimModes simulation_mode = GRAVEL;
// Initial vehicle location and orientation
ChVector<> initLoc(0, 0, .6);
ChQuaternion<> initRot(1, 0, 0, 0);
// ChQuaternion<> initRot(0.866025, 0, 0, 0.5);
// ChQuaternion<> initRot(0.7071068, 0, 0, 0.7071068);
// ChQuaternion<> initRot(0.25882, 0, 0, 0.965926);
// ChQuaternion<> initRot(0, 0, 0, 1);

double mph_to_m_s = 0.44704;
double target_speed = 35 * mph_to_m_s;  // 35 mph

// Visualization type for vehicle parts (PRIMITIVES, MESH, or NONE)
VisualizationType chassis_vis_type = VisualizationType::NONE;
VisualizationType suspension_vis_type = VisualizationType::PRIMITIVES;
VisualizationType steering_vis_type = VisualizationType::PRIMITIVES;
VisualizationType wheel_vis_type = VisualizationType::PRIMITIVES;
VisualizationType tire_vis_type = VisualizationType::PRIMITIVES;
// Type of powertrain model (SHAFTS, SIMPLE)
PowertrainModelType powertrain_model = PowertrainModelType::SHAFTS;

// Drive type (FWD, RWD, or AWD)
DrivelineType drive_type = DrivelineType::AWD;

// Type of tire model (RIGID, RIGID_MESH, PACEJKA, LUGRE, FIALA)
TireModelType tire_model = TireModelType::RIGID;

// Point on chassis tracked by the camera
ChVector<> trackPoint(0.0, 0.0, 1.75);

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

std::string circle_file("humvee_input/circlepath_100ft_radius.txt");
std::string lane_change_file("humvee_input/ISO_double_lane_change.txt");
std::string gravel_driver_file("humvee_input/straightOrigin.txt");
std::string steering_controller_file("humvee_input/SteeringController.json");
std::string speed_controller_file("humvee_input/SpeedController.json");

std::string data_output_path = "";

#define ERASE_MACRO(x, y) x.erase(x.begin() + y);
#define ERASE_MACRO_LEN(x, y, z) x.erase(x.begin() + y, x.begin() + y + z);

double hdimX = 4;
double hdimY = 1.5;
double hdimZ = 0.1;
double hthick = 0.1;

int Id_g = 100;
double r_g = 0.02;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

float mu_g = 1.0f;
std::vector<real3> forces;
std::vector<real3> torques;

void static WriteVehicleData(hmmwv::HMMWV_Full& my_hmmwv,
                             ChDriver::Inputs driver_inputs,
                             std::vector<real3> forces,
                             std::vector<real3> torques,
                             std::string filename) {
    CSVGen csv_output;
    csv_output.OpenFile(filename.c_str(), false);

    auto m_driveline = my_hmmwv.GetVehicle().GetDriveline();
    auto m_suspension_front = std::dynamic_pointer_cast<ChDoubleWishbone>(my_hmmwv.GetVehicle().GetSuspension(0));
    auto m_suspension_back = std::dynamic_pointer_cast<ChDoubleWishbone>(my_hmmwv.GetVehicle().GetSuspension(1));

    csv_output << my_hmmwv.GetVehicle().GetChassis()->GetPos();
    csv_output << my_hmmwv.GetVehicle().GetVehicleSpeed();
    csv_output << m_driveline->GetDriveshaftSpeed();
    csv_output << my_hmmwv.GetPowertrain()->GetMotorTorque();
    csv_output << my_hmmwv.GetPowertrain()->GetMotorSpeed();
    csv_output << my_hmmwv.GetPowertrain()->GetOutputTorque();

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

void RemoveCollisionModel(ChSystemParallelNSC* system, ChCollisionModel* model) {
    ChParallelDataManager* data_manager = system->data_manager;
#if 1
    ChCollisionModelParallel* pmodel = static_cast<ChCollisionModelParallel*>(model);
    int body_id = pmodel->GetBody()->GetId();
    // loop over the models we nned to remove
    // std::cout << "removing: " << pmodel->GetNObjects() << " objects" << std::endl;
    for (int j = 0; j < pmodel->GetNumShapes(); j++) {
        // find a model to remove
        bool removed = false;
        int index = -1;
        for (int i = 0; i < data_manager->shape_data.id_rigid.size(); i++) {
            if (data_manager->shape_data.id_rigid[i] == body_id) {
                index = i;
                data_manager->num_rigid_shapes--;

                int start = data_manager->shape_data.start_rigid[index];
                int length = data_manager->shape_data.length_rigid[index];
                int type = data_manager->shape_data.typ_rigid[index];

                // std::cout << "removing: type " << type << " " << start<< " " <<j << std::endl;

                switch (type) {
                    case ChCollisionShape::Type::SPHERE:
                        ERASE_MACRO_LEN(data_manager->shape_data.sphere_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::ELLIPSOID:
                        ERASE_MACRO_LEN(data_manager->shape_data.box_like_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::BOX:
                        ERASE_MACRO_LEN(data_manager->shape_data.box_like_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::CYLINDER:
                        ERASE_MACRO_LEN(data_manager->shape_data.box_like_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::CONE:
                        ERASE_MACRO_LEN(data_manager->shape_data.box_like_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::CAPSULE:
                        ERASE_MACRO_LEN(data_manager->shape_data.capsule_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::ROUNDEDBOX:
                        ERASE_MACRO_LEN(data_manager->shape_data.rbox_like_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::ROUNDEDCYL:
                        ERASE_MACRO_LEN(data_manager->shape_data.rbox_like_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::ROUNDEDCONE:
                        ERASE_MACRO_LEN(data_manager->shape_data.rbox_like_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::CONVEX:
                        ERASE_MACRO_LEN(data_manager->shape_data.convex_rigid, start, length);
                        break;
                    case ChCollisionShape::Type::TRIANGLE:
                        ERASE_MACRO_LEN(data_manager->shape_data.convex_rigid, start, 3);
                        break;
                }

                ERASE_MACRO(data_manager->shape_data.ObA_rigid, index);
                ERASE_MACRO(data_manager->shape_data.ObR_rigid, index);
                ERASE_MACRO(data_manager->shape_data.start_rigid, index);
                ERASE_MACRO(data_manager->shape_data.length_rigid, index);

                ERASE_MACRO(data_manager->shape_data.fam_rigid, index);
                ERASE_MACRO(data_manager->shape_data.typ_rigid, index);
                ERASE_MACRO(data_manager->shape_data.id_rigid, index);
                removed = true;
                break;
            }
        }
        // std::cout << "decrement start "<< std::endl;
        if (removed) {
            // we removed a model, all of the starts are off by one for indices past the index removed, decrement all
            // starts before removing a second model
            for (int i = index; i < data_manager->shape_data.start_rigid.size(); i++) {
                if (data_manager->shape_data.start_rigid[i] != 0) {
                    data_manager->shape_data.start_rigid[i] -= 1;
                }
            }
        }
    }
#endif
}

ChFluidContainer* fluid_container;

void CreateFluid(ChSystemParallelNSC* system) {
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

    utils::HCPSampler<> sampler(dist);
    utils::Generator::PointVector points =
        sampler.SampleCylinderY(initLoc + ChVector<>(-1.9306 + .24, .025, 1.2), .4, .5);
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
void CreateBase(ChSystemParallelNSC* system) {
    Vector c_pos = Vector(0, 0, 0);

    bottom_plate = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
    bottom_plate->SetMass(1);
    bottom_plate->SetPos(Vector(0, 0, 0) + c_pos);
    bottom_plate->SetRot(Quaternion(1, 0, 0, 0));
    bottom_plate->SetCollide(true);
    bottom_plate->SetBodyFixed(true);

    auto material = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    material->SetFriction(container_friction);
    material->SetCompliance(1e-9);
    material->SetCohesion(0);

    bottom_plate->GetCollisionModel()->ClearModel();
    bottom_plate->GetCollisionModel()->SetFamily(2);
    bottom_plate->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);

    if (simulation_mode == TURN) {
        AddBoxGeometry(bottom_plate.get(), material, Vector(4, 4, .1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
    } else if (simulation_mode == LANE_CHANGE) {
        AddBoxGeometry(bottom_plate.get(), material, Vector(4, 3, .1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
    } else if (simulation_mode == GRAVEL) {

        AddBoxGeometry(bottom_plate.get(), material, Vector(4, 2, .1), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

        AddBoxGeometry(bottom_plate.get(), material, Vector(10, 2, .1), Vector(4 + hdimX * 2 + 10.0, 0, 0),
                       Quaternion(1, 0, 0, 0));

        // extra side plates
        AddBoxGeometry(bottom_plate.get(), material, Vector(12.5, .1, 2), Vector(4 + hdimX, -hdimY - .1, 2 - .1),
                       Quaternion(1, 0, 0, 0));
        AddBoxGeometry(bottom_plate.get(), material, Vector(12.5, .1, 2), Vector(4 + hdimX, hdimY + .1, 2 - .1),
                       Quaternion(1, 0, 0, 0));
        // Top
        AddBoxGeometry(bottom_plate.get(), material, Vector(12.5, 2 + .1 * 2, .1), Vector(4 + hdimX, 0, 4 - .1),
                       Quaternion(1, 0, 0, 0));
    }

    bottom_plate->GetCollisionModel()->BuildModel();
    system->AddBody(bottom_plate);

    if (simulation_mode == GRAVEL) {
        CreateBoxContainer(system, 1, material, ChVector<>(hdimX, hdimY, hdimZ), hthick,
                           ChVector<>(4 + hdimX, 0, -hdimZ * 2 - hthick));
    }
}

void CreateGravel(ChSystemParallelNSC* system) {
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);

    utils::Generator gen(system);
    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::BISPHERE, 1.0);
    // std::shared_ptr<utils::MixtureIngredient> m2 = gen.AddMixtureIngredient(utils::MixtureType::CYLINDER, 1.0);
    // std::shared_ptr<utils::MixtureIngredient> m3 = gen.AddMixtureIngredient(utils::MixtureType::BOX, 1.0);

    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);
    // m1->setDistributionSize(r_g, r_g, r_g, 1.5 * r_g);
    // Set starting value for body identifiers
    gen.setBodyIdentifier(Id_g);

    // Create particles in layers until reaching the desired number of particles
    double r = 1.01 * r_g;
    ChVector<> hdims(hdimX - r * 3, hdimY - r * 3, 0);
    ChVector<> center(4 + hdimX, 0, -.2);

    while (gen.getTotalNumBodies() < 100000) {
        gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 3 * r, center, hdims);
        center.z() += 2 * r;
    }

    std::cout << "Created " << gen.getTotalNumBodies() << " particles." << std::endl;

    // return center.z();
}

int main(int argc, char* argv[]) {
    if (argc == 3) {
        simulation_mode = SimModes(atoi(argv[1]));
        target_speed = atoi(argv[2]) * mph_to_m_s;
    }

    // --------------
    // Create systems
    // --------------
    ChSystemParallelNSC* system = new ChSystemParallelNSC();
    system->Set_G_acc(ChVector<>(0, 0, -9.81));

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
    system->GetSettings()->min_threads = threads;
    system->ChangeSolverType(SolverType::BB);
    system->SetLoggingLevel(LoggingLevel::LOG_INFO);
    system->SetLoggingLevel(LoggingLevel::LOG_TRACE);

    system->GetSettings()->collision.collision_envelope = 0.1 * fluid_r;

    system->GetSettings()->collision.bins_per_axis = vec3(300, 40, 50);
    system->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
    system->GetSettings()->collision.fixed_bins = true;

    std::shared_ptr<ChBezierCurve> path;

    if (simulation_mode == TURN) {
        initLoc = ChVector<>(0, 30.48, .6);
        path = ChBezierCurve::read(circle_file);
        data_output_path = "humvee_tank_turn_" + std::to_string(atoi(argv[2])) + "/";
        system->GetSettings()->collision.bins_per_axis = vec3(100, 40, 50);
    } else if (simulation_mode == LANE_CHANGE) {
        initLoc = ChVector<>(-50, -125, .6);
        path = ChBezierCurve::read(lane_change_file);
        data_output_path = "humvee_tank_dlc_" + std::to_string(atoi(argv[2])) + "/";
        system->GetSettings()->collision.bins_per_axis = vec3(100, 40, 50);
    } else if (simulation_mode == GRAVEL) {
        //		target_speed = 15 * mph_to_m_s;
        path = ChBezierCurve::read(gravel_driver_file);
        data_output_path = "humvee_tank_gravel_" + std::to_string(atoi(argv[2])) + "/";
        system->GetSettings()->collision.bins_per_axis = vec3(300, 40, 50);
    }

    std::cout << "writing data to: " << data_output_path << std::endl;
    CreateBase(system);
    CreateFluid(system);

    // Create the HMMWV vehicle, set parameters, and initialize
    HMMWV_Full my_hmmwv(system);
    my_hmmwv.SetContactMethod(contact_method);
    my_hmmwv.SetChassisFixed(false);
    my_hmmwv.SetInitPosition(ChCoordsys<>(initLoc, initRot));
    my_hmmwv.SetPowertrainType(powertrain_model);
    my_hmmwv.SetDriveType(drive_type);
    my_hmmwv.SetTireType(tire_model);
    my_hmmwv.SetTireStepSize(tire_step_size);

    my_hmmwv.Initialize();

    my_hmmwv.SetChassisVisualizationType(chassis_vis_type);
    my_hmmwv.SetSuspensionVisualizationType(suspension_vis_type);
    my_hmmwv.SetSteeringVisualizationType(steering_vis_type);
    my_hmmwv.SetWheelVisualizationType(wheel_vis_type);
    my_hmmwv.SetTireVisualizationType(tire_vis_type);

    RemoveCollisionModel(system, my_hmmwv.GetChassisBody()->GetCollisionModel().get());
    my_hmmwv.GetChassisBody()->GetCollisionModel()->ClearModel();
    my_hmmwv.GetChassisBody()->GetAssets().clear();

    std::shared_ptr<ChBodyAuxRef> m_chassis = my_hmmwv.GetChassisBody();

    auto material = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    // add sphere to center to mark chassis pos
    utils::AddSphereGeometry(m_chassis.get(), material, 0.115, ChVector<>(0, 0, 0), QUNIT);

    ChVector<> hdim(.1, .6, .05);

    for (int i = 0; i < m_chassis->GetAssets().size(); i++) {
        ((ChVisualization*)m_chassis->GetAssets().at(i).get())->Pos = ChVector<>(0);
    }

    for (int j = 0; j < 30; j++) {
        ChQuaternion<> qq = Q_from_AngAxis(j * 12 * CH_C_DEG_TO_RAD, VECT_Y);
        ChVector<> ppp = ChVector<>(sin(12 * CH_C_DEG_TO_RAD * j), 0, cos(12 * CH_C_DEG_TO_RAD * j)) * .6;

        utils::AddBoxGeometry(m_chassis.get(), material, hdim, ppp + ChVector<>(-1.9306 + .24, 0, 1.3011), qq);
    }
    utils::AddCylinderGeometry(m_chassis.get(), material, .6, .05, ChVector<>(-1.9306 + .24, -.6, 1.3011), QUNIT);
    utils::AddCylinderGeometry(m_chassis.get(), material, .6, .05, ChVector<>(-1.9306 + .24, .6, 1.3011), QUNIT);

    my_hmmwv.GetChassisBody()->SetCollide(true);
    my_hmmwv.GetChassisBody()->GetCollisionModel()->BuildModel();
    auto collision_system = std::static_pointer_cast<ChCollisionSystemParallel>(system->GetCollisionSystem());
    collision_system->Add(my_hmmwv.GetChassisBody()->GetCollisionModel().get());

    RigidTerrain terrain(my_hmmwv.GetSystem());

    if (simulation_mode == GRAVEL) {
        CreateGravel(system);
    }

    // ------------------------
    // Create the driver system
    // ------------------------

    ChPathFollowerDriver driver(my_hmmwv.GetVehicle(), steering_controller_file, speed_controller_file, path, "my_path",
                                target_speed, true);
    driver.Initialize();

    // ---------------
    // Simulation loop
    // ---------------

    double time = 0;
    double exec_time = 0;

    int sim_frame = 0, out_frame = 0, next_out_frame = 0;

    int out_steps = std::ceil((1.0 / time_step) / out_fps);

#ifdef CHRONO_OPENGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Humvee", system);

    gl_window.SetCamera(ChVector<>(0, -10, 0) + initLoc, initLoc, ChVector<>(0, 0, 1));

    gl_window.SetRenderMode(opengl::WIREFRAME);

    gl_window.Pause();
#endif

    ChVector<> driver_pos = my_hmmwv.GetVehicle().GetChassis()->GetLocalDriverCoordsys().pos;

    while (time < time_end) {
        if (simulation_mode != GRAVEL) {
            bottom_plate->SetPos(ChVector<>(my_hmmwv.GetChassis()->GetBody()->GetPos().x(),
                                            my_hmmwv.GetChassis()->GetBody()->GetPos().y(), 0));
        }

        ChDriver::Inputs driver_inputs = driver.GetInputs();

        ChVector<> com_pos = my_hmmwv.GetChassis()->GetCOMPos();
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

#ifdef CHRONO_OPENGL
        if (gl_window.Active()) {
            if (gl_window.DoStepDynamics(time_step)) {
                // Update counters.
                time += time_step;
                sim_frame++;
                exec_time += system->GetTimerStep();
            }
            gl_window.Render();
        } else {
            break;
        }
#else
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

#endif
    }
    cout << "==================================" << endl;
    cout << "Simulation time:   " << exec_time << endl;
    return 0;
}
