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
// Author: Daniel Melanz, Radu Serban
// =============================================================================
//
// Chrono::Multicore demo program for single wheel studies.
//
// The system contains a wheel test composed of three bodies: (1) a containing
// bin, (2) a wheel, (3) a chassis, and (4) an axle. Granular material sits in
// the containing bin. The wheel body is settled from the top. During the
// rolling mode, the wheel is translated and rotated in the x-direction at a
// specified slip.
//
// The global reference frame has Z up.
// All units SI (CGS, i.e., centimeter - gram - second)
//
// =============================================================================

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChSystemDescriptorMulticore.h"

#include "chrono_thirdparty/filesystem/path.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

using namespace chrono;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use NSC contact
//#define USE_SMC

enum ProblemType { SETTLING, PRESSING, ROLLING, TESTING };

ProblemType problem = SETTLING;

// -----------------------------------------------------------------------------
// Conversion factors
// -----------------------------------------------------------------------------

// Conversion for weight of wheel
// [Y] = N = kg * m / s^2
double N2cgs = 100000;

// Conversion for Young's modulus and pressure
// [Y] = Pa = N / m^2 = kg / m / s^2
double Pa2cgs = 10;

// Conversion for dissipation factor alpha
// [alpha] = s / m
double alpha2cgs = 0.01;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Save PovRay post-processing data?
bool write_povray_data = true;

// Simulation times
double time_settling_min = 0.1;
double time_settling_max = 1.0;

double time_pressing_min = 0.1;
double time_pressing_max = 1.0;

double time_rolling = 10;

double time_testing = 2;

// Stopping criteria for settling (as fraction of particle radius)
double settling_tol = 0.2;

// Solver settings
#ifdef USE_SMC
double time_step = 1e-5;
int max_iteration_bilateral = 100;
#else
double time_step = 1e-3;
int max_iteration_normal = 0;
int max_iteration_sliding = 10000;
int max_iteration_spinning = 0;
int max_iteration_bilateral = 0;
double contact_recovery_speed = 10e30;
#endif

bool clamp_bilaterals = false;
double bilateral_clamp_speed = 10e30;
double tolerance = 1;

// Output
#ifdef USE_SMC
const std::string out_dir = "../SINGLEWHEEL_SMC";
#else
const std::string out_dir = "../SINGLEWHEEL_NSC";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string roll_file = out_dir + "/roll.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string settled_ckpnt_file = out_dir + "/settled.dat";
const std::string pressed_ckpnt_file = out_dir + "/pressed.dat";

// Frequency for visualization output
int out_fps_settling = 120;
int out_fps_pressing = 120;
int out_fps_rolling = 60;
int out_fps_testing = 60;

// Frequency for writing results to output file (ROLLING only)
int write_fps = 1000;

// Simulation frame at which detailed timing information is printed
int timing_frame = -1;

// Gravitational acceleration [cm/s^2]
double gravity = 981;

// Parameters for the mechanism
double hdimX = 100.0 / 2;  // [cm] bin half-length in x direction
double hdimY = 60.0 / 2;   // [cm] bin half-depth in y direction
double hdimZ = 32.0 / 2;   // [cm] bin half-height in z direction
double hthick = 1.0 / 2;   // [cm] bin half-thickness of the walls

double wheelRadius = 13;               // [cm] radius of the wheel
double wheelWidth = 16;                // [cm] width of the wheel
double wheelSlip = 0.3;                // [-]  ratio of angular and translation velocity of the wheel
double angVel = 17 * CH_PI / 180.0;  // [rad/s] angular velocity of the wheel
double wheelWeight = 80;               // * N2cgs;   // Normal load of the wheel
double velocity = angVel * wheelRadius * (1.0 - wheelSlip);

float Y_walls = (float)(Pa2cgs * 2e6);
float mu_walls = 0.3f;

// Parameters for the granular material
int tag_particles = 0; // start tag ID for particles
double r_g = 2.0;      // [cm] radius of granular sphers
double rho_g = 2.500;  // [g/cm^3] density of granules

double desiredBulkDensity = 1.3894;  // [g/cm^3] desired bulk density

float Y_g = (float)(Pa2cgs * 5e7);
float mu_g = 0.5f;

// Parameters of the testing ball
int Id_ball = -4;
double mass_ball = 200;                  // [g] mass of testing ball
double radius_ball = 0.1 * wheelRadius;  // [cm] radius of testing ball

// =============================================================================
// Create the containing bin (the ground), wheel, chassis, and axle.
//
// The wheel, chassis, and axle are constructed fixed to the ground.
// No joints between bodies are defined at this time.
// =============================================================================

void CreateMechanismBodies(ChSystemMulticore* system) {
    // -------------------------------
    // Create a material for the walls
    // -------------------------------

#ifdef USE_SMC
    auto mat_walls = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_walls->SetYoungModulus(Y_walls);
    mat_walls->SetFriction(mu_walls);
#else
    auto mat_walls = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_walls->SetFriction(mu_walls);
#endif

    // ----------------------
    // Create the ground body -- always FIRST body in system
    // ----------------------

    auto ground = chrono_types::make_shared<ChBody>();

    ground->SetFixed(true);
    ground->EnableCollision(true);

    // Attach geometry of the containing bin.  Disable contact ground-shearBox
    // and ground-loadPlate.
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector3d(hdimX, hdimY, hthick), ChVector3d(0, 0, -hthick));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector3d(hthick, hdimY, hdimZ), ChVector3d(-hdimX - hthick, 0, hdimZ));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector3d(hthick, hdimY, hdimZ), ChVector3d(hdimX + hthick, 0, hdimZ));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector3d(hdimX, hthick, hdimZ), ChVector3d(0, -hdimY - hthick, hdimZ));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector3d(hdimX, hthick, hdimZ), ChVector3d(0, hdimY + hthick, hdimZ));

    system->AddBody(ground);

    // --------------------
    // Create the wheel -- always SECOND body in system
    // --------------------

    // Initially, the wheel is fixed to ground.

    auto wheel = chrono_types::make_shared<ChBody>();

    // Set wheel starting pos, near back of the containing bin
    double hdimX_w = -hdimX + 1.01 * wheelRadius;

    wheel->SetPos(ChVector3d(hdimX_w, 0, 2 * hdimZ + wheelRadius));
    wheel->SetMass(1.0);
    wheel->EnableCollision(true);
    wheel->SetFixed(true);

    // Add geometry of the wheel.
    utils::AddCylinderGeometry(wheel.get(), mat_walls, wheelRadius, wheelWidth / 2, ChVector3d(0, 0, 0), QUNIT);

    system->AddBody(wheel);

    // ---------------------
    // Create the chassis -- always THIRD body in the system
    // ---------------------

    // Initially, the chassis is fixed to ground.
    // It is released after the settling phase.

    auto chassis = chrono_types::make_shared<ChBody>();

    chassis->SetMass(1.0);
    chassis->SetPos(wheel->GetPos());
    chassis->EnableCollision(false);
    chassis->SetFixed(true);

    // Add geometry of the wheel.
    utils::AddBoxGeometry(chassis.get(), mat_walls, ChVector3d(.1 * wheelRadius, .1 * wheelRadius, .1 * wheelRadius),
                          ChVector3d(0, 0, 0));

    system->AddBody(chassis);

    // ---------------------
    // Create the axle -- always FOURTH body in the system
    // ---------------------

    // Initially, the axle is fixed to ground.
    // It is released after the settling phase.

    auto axle = chrono_types::make_shared<ChBody>();

    axle->SetMass(wheelWeight / gravity);
    axle->SetPos(wheel->GetPos());
    axle->EnableCollision(false);
    axle->SetFixed(true);

    // Add geometry of the wheel.
    utils::AddSphereGeometry(axle.get(), mat_walls, 0.1 * wheelRadius, ChVector3d(0, 0, 0));

    system->AddBody(axle);
}

// =============================================================================
// Connect the chassis to the containing bin (ground) through a translational
// joint and create a linear actuator.
// =============================================================================

void ConnectChassisToGround(ChSystemMulticore* system, std::shared_ptr<ChBody> ground, std::shared_ptr<ChBody> chassis) {
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(ground, chassis, ChFrame<>(chassis->GetPos(), QuatFromAngleY(CH_PI_2)));
    prismatic->SetName("prismatic_chassis_ground");
    system->AddLink(prismatic);

    velocity = angVel * wheelRadius * (1.0 - wheelSlip);
    auto actuator_fun = chrono_types::make_shared<ChFunctionRamp>(0.0, velocity);

    auto actuator = chrono_types::make_shared<ChLinkLockLinActuator>();
    actuator->Initialize(ground, chassis, false, ChFrame<>(chassis->GetPos(), QUNIT),
                         ChFrame<>(chassis->GetPos() + ChVector3d(1, 0, 0), QUNIT));
    actuator->SetName("actuator");
    actuator->SetDistanceOffset(1);
    actuator->SetActuatorFunction(actuator_fun);
    system->AddLink(actuator);
}

// =============================================================================
// Connect the axle to the chassis through a vertical translational
// joint.
// =============================================================================

void ConnectChassisToAxle(ChSystemMulticore* system, std::shared_ptr<ChBody> chassis, std::shared_ptr<ChBody> axle) {
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(chassis, axle, ChFrame<>(chassis->GetPos(), QUNIT));
    prismatic->SetName("prismatic_axle_chassis");
    system->AddLink(prismatic);
}

// =============================================================================
// Connect the wheel to the axle through a engine joint.
// =============================================================================

void ConnectWheelToAxle(ChSystemMulticore* system, std::shared_ptr<ChBody> wheel, std::shared_ptr<ChBody> axle) {
    auto motor = chrono_types::make_shared<ChLinkMotorRotationAngle>();
    motor->SetName("engine_wheel_axle");
    motor->Initialize(wheel, axle, ChFrame<>(wheel->GetPos(), chrono::QuatFromAngleAxis(CH_PI / 2.0, VECT_X)));
    motor->SetAngleFunction(chrono_types::make_shared<ChFunctionRamp>(0, -angVel));
    system->AddLink(motor);
}

// =============================================================================
// Create the granular material
//
// Granular material consisting of identical spheres with specified radius and
// material properties; the spheres are generated in a number of vertical
// layers with locations within each layer obtained using Poisson Disk sampling,
// thus ensuring that no two spheres are closer than twice the radius.
// =============================================================================

int CreateGranularMaterial(ChSystemMulticore* system) {
// -------------------------------------------
// Create a material for the granular material
// -------------------------------------------

#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
#else
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // ---------------------------------------------
    // Create a mixture entirely made out of spheres
    // ---------------------------------------------

    // Create the particle generator with a mixture of 100% spheres
    double r = 1.01 * r_g;
    utils::ChPDSampler<double> sampler(2 * r);
    utils::ChGenerator gen(system);

    std::shared_ptr<utils::ChMixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->SetDefaultMaterial(mat_g);
    m1->SetDefaultDensity(rho_g);
    m1->SetDefaultSize(r_g);

    // Ensure that all generated particle bodies will have positive IDs.
    gen.SetStartTag(tag_particles);

    // ----------------------
    // Generate the particles
    // ----------------------

    ChVector3d hdims(hdimX - r, hdimY - r, 0);
    ChVector3d center(0, 0, 2 * r);

    while (center.z() < 2 * hdimZ) {
        gen.CreateObjectsBox(sampler, center, hdims);
        center.z() += 2 * r;
    }

    // Return the number of generated particles.
    return gen.GetTotalNumBodies();
}

// =============================================================================
// Create a single large sphere (for use in TESTING)
// =============================================================================

void CreateBall(ChSystemMulticore* system) {
    // ------------------------------
    // Create a material for the ball
    // ------------------------------

#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
#else
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // ---------------
    // Create the ball
    // ---------------

    auto ball = chrono_types::make_shared<ChBody>();

    ball->SetMass(mass_ball);
    ball->SetPos(ChVector3d(0, 0, 1.01 * radius_ball));
    ball->EnableCollision(true);
    ball->SetFixed(false);

    utils::AddSphereGeometry(ball.get(), mat_g, radius_ball);

    system->AddBody(ball);
}

// =============================================================================
// Find the height of the highest and lowest sphere in the granular mix. We only
// look at bodies with positive identifiers (to exclude all other bodies).
// =============================================================================

void FindHeightRange(ChSystemMulticore* sys, double& lowest, double& highest) {
    highest = -1000;
    lowest = 1000;
    for (auto body : sys->GetBodies()) {
        if (body->GetTag() < tag_particles)
            continue;
        double h = body->GetPos().z();
        if (h < lowest)
            lowest = h;
        else if (h > highest)
            highest = h;
    }
}

// =============================================================================
//
//// TODO:  cannot do this with SMC!!!!!
//
// =============================================================================

void setBulkDensity(ChSystem* sys, double bulkDensity) {
    double vol_g = (4.0 / 3) * CH_PI * r_g * r_g * r_g;

    double normalPlateHeight = sys->GetBodies().at(1)->GetPos().z() - hdimZ;
    double bottomHeight = 0;
    double boxVolume = hdimX * 2 * hdimX * 2 * (normalPlateHeight - bottomHeight);
    double granularVolume = (sys->GetBodies().size() - 3) * vol_g;
    double reqDensity = bulkDensity * boxVolume / granularVolume;
    for (auto body : sys->GetBodies()) {
        if (body->GetTag() >= tag_particles) {
            body->SetMass(reqDensity * vol_g);
        }
    }

    cout << "N Bodies: " << sys->GetBodies().size() << endl;
    cout << "Box Volume: " << boxVolume << endl;
    cout << "Granular Volume: " << granularVolume << endl;
    cout << "Desired bulk density = " << bulkDensity << ", Required Body Density = " << reqDensity << endl;
}

// =============================================================================

int main(int argc, char* argv[]) {
    // Create output directories.
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

// -------------
// Create system
// -------------

#ifdef USE_SMC
    cout << "Create SMC system" << endl;
    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemMulticoreNSC* sys = new ChSystemMulticoreNSC();
#endif

    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    sys->SetGravitationalAcceleration(ChVector3d(0, 0, -gravity));

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    sys->SetNumThreads(threads);
    cout << "Using " << threads << " threads" << endl;

    // Edit system settings
    sys->GetSettings()->solver.tolerance = tolerance;
    sys->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    sys->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
    sys->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

#ifdef USE_SMC
    sys->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_R;
#else
    sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    sys->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    sys->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    sys->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    sys->GetSettings()->solver.alpha = 0;
    sys->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    sys->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
    sys->ChangeSolverType(SolverType::APGDREF);

    sys->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // --------------
    // Problem set up
    // --------------

    // Depending on problem type:
    // - Select end simulation time
    // - Select output FPS
    // - Create / load objects

    double time_min = 0;
    double time_end;
    int out_fps;
    std::shared_ptr<ChBody> ground;
    std::shared_ptr<ChBody> wheel;
    std::shared_ptr<ChBody> chassis;
    std::shared_ptr<ChBody> axle;
    std::shared_ptr<ChLinkLockPrismatic> prismatic_chassis_ground;
    std::shared_ptr<ChLinkLockPrismatic> prismatic_axle_chassis;
    std::shared_ptr<ChLinkLockLinActuator> actuator;
    std::shared_ptr<ChLinkMotorRotationAngle> engine_wheel_axle;

    switch (problem) {
        case SETTLING: {
            time_min = time_settling_min;
            time_end = time_settling_max;
            out_fps = out_fps_settling;

            // Create the mechanism bodies (all fixed).
            CreateMechanismBodies(sys);

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = sys->GetBodies().at(0);
            wheel = sys->GetBodies().at(1);
            chassis = sys->GetBodies().at(2);
            axle = sys->GetBodies().at(3);

            // Create granular material.
            int num_particles = CreateGranularMaterial(sys);
            cout << "Granular material:  " << num_particles << " particles" << endl;

            break;
        }

        case PRESSING: {
            time_min = time_pressing_min;
            time_end = time_pressing_max;
            out_fps = out_fps_pressing;

            // Create bodies from checkpoint file.
            cout << "Read checkpoint data from " << settled_ckpnt_file;
            utils::ReadCheckpoint(sys, settled_ckpnt_file);
            cout << "  done.  Read " << sys->GetBodies().size() << " bodies." << endl;

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = sys->GetBodies().at(0);
            wheel = sys->GetBodies().at(1);
            chassis = sys->GetBodies().at(2);
            axle = sys->GetBodies().at(3);

            // Move the load plate just above the granular material.
            double highest, lowest;
            FindHeightRange(sys, lowest, highest);
            ChVector3d pos = wheel->GetPos();
            double z_new = highest + 1.01 * r_g + wheelRadius;
            wheel->SetPos(ChVector3d(pos.x(), pos.y(), z_new));
            chassis->SetPos(wheel->GetPos());
            axle->SetPos(wheel->GetPos());

            // Connect the chassis to the axle.
            ConnectChassisToAxle(sys, chassis, axle);
            prismatic_axle_chassis = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_axle_chassis"));

            // Release the axle.
            axle->SetFixed(false);

            // Set axle mass from desired applied normal load
            axle->SetMass(wheelWeight / gravity);

            // Connect the axle to the chassis.
            angVel = 0.0;  // Don't rotate the wheel in this stage
            ConnectWheelToAxle(sys, wheel, axle);
            engine_wheel_axle = std::static_pointer_cast<ChLinkMotorRotationAngle>(sys->SearchLink("engine_wheel_axle"));

            // Release the axle.
            wheel->SetFixed(false);

            break;
        }

        case ROLLING: {
            time_end = time_rolling;
            out_fps = out_fps_rolling;

            // Create bodies from checkpoint file.
            cout << "Read checkpoint data from " << pressed_ckpnt_file;
            utils::ReadCheckpoint(sys, pressed_ckpnt_file);
            cout << "  done.  Read " << sys->GetBodies().size() << " bodies." << endl;

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = sys->GetBodies().at(0);
            wheel = sys->GetBodies().at(1);
            chassis = sys->GetBodies().at(2);
            axle = sys->GetBodies().at(3);

            // Connect the chassis and get a handle to the actuator.
            ConnectChassisToGround(sys, ground, chassis);
            prismatic_chassis_ground = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_chassis_ground"));
            actuator = std::static_pointer_cast<ChLinkLockLinActuator>(sys->SearchLink("actuator"));

            // Release the load plate.
            chassis->SetFixed(false);

            // Connect the load plate to the shear box.
            ConnectChassisToAxle(sys, chassis, axle);
            prismatic_axle_chassis = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_axle_chassis"));

            // Release the axle.
            axle->SetFixed(false);

            // Set axle mass from desired applied normal load
            axle->SetMass(wheelWeight / gravity);

            // Connect the axle to the chassis.
            ConnectWheelToAxle(sys, wheel, axle);
            engine_wheel_axle = std::static_pointer_cast<ChLinkMotorRotationAngle>(sys->SearchLink("engine_wheel_axle"));

            // Release the axle.
            wheel->SetFixed(false);

            // setBulkDensity(sys, desiredBulkDensity);

            break;
        }

        case TESTING: {
            time_end = time_testing;
            out_fps = out_fps_testing;
            angVel = 10 * CH_PI;

            // Create the mechanism bodies (all fixed).
            CreateMechanismBodies(sys);

            // Create the test ball.
            CreateBall(sys);

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = sys->GetBodies().at(0);
            wheel = sys->GetBodies().at(1);
            chassis = sys->GetBodies().at(2);
            axle = sys->GetBodies().at(3);

            // Move the wheel just above the ground.
            ChVector3d pos = wheel->GetPos();
            double z_new = wheelRadius;
            axle->SetPos(ChVector3d(pos.x(), pos.y(), z_new));
            wheel->SetPos(ChVector3d(pos.x(), pos.y(), z_new));

            // Connect the chassis and get a handle to the actuator.
            ConnectChassisToGround(sys, ground, chassis);
            prismatic_chassis_ground = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_chassis_ground"));
            actuator = std::static_pointer_cast<ChLinkLockLinActuator>(sys->SearchLink("actuator"));

            chassis->SetFixed(false);

            // Connect the axle to the chassis.
            ConnectChassisToAxle(sys, chassis, axle);
            prismatic_axle_chassis = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_axle_chassis"));

            // Release the axle.
            axle->SetFixed(false);

            // Set axle mass from desired applied normal load
            axle->SetMass(wheelWeight / gravity);

            // Connect the axle to the chassis.
            ConnectWheelToAxle(sys, wheel, axle);
            engine_wheel_axle = std::static_pointer_cast<ChLinkMotorRotationAngle>(sys->SearchLink("engine_wheel_axle"));

            // Release the axle.
            wheel->SetFixed(false);

            break;
        }
    }

    // ----------------------
    // Perform the simulation
    // ----------------------

    // Set number of simulation steps and steps between successive output
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);
    int write_steps = (int)std::ceil((1.0 / time_step) / write_fps);

    // Initialize counters
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;
    double max_cnstr_viol[3] = {0, 0, 0};

    // Circular buffer with highest particle location
    // (only used for SETTLING or PRESSING)
    int buffer_size = (int)std::ceil(time_min / time_step);
    std::valarray<double> hdata(0.0, buffer_size);

    // Create output files
    std::ofstream statsStream(stats_file.c_str());
    std::ofstream rollStream(roll_file.c_str());

#ifdef CHRONO_OPENGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(sys);
    vis.SetWindowTitle("Single Wheel Test");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(0, -10 * hdimY, hdimZ), ChVector3d(0, 0, hdimZ));
    vis.SetCameraVertical(CameraVerticalDir::Z);
#endif

    // Loop until reaching the end time...
    while (time < time_end) {
        // Current position and velocity of the wheel
        ChVector3d pos_old = wheel->GetPos();
        ChVector3d vel_old = wheel->GetLinVel();

        // Calculate minimum and maximum particle heights
        double highest, lowest;
        FindHeightRange(sys, lowest, highest);

        // If at an output frame, write PovRay file and print info
        if (sim_frame == next_out_frame) {
            cout << "------------ Output frame:     " << out_frame + 1 << endl;
            cout << "             Sim frame:        " << sim_frame << endl;
            cout << "             Time:             " << time << endl;
            cout << "             Wheel pos:        " << pos_old.x() << endl;
            cout << "                   vel:        " << vel_old.x() << endl;
            cout << "             Particle lowest:  " << lowest << endl;
            cout << "                      highest: " << highest << endl;
            cout << "             Execution time:   " << exec_time << endl;

            // Save PovRay post-processing data.
            if (write_povray_data) {
                char filename[100];
                sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
                utils::WriteVisualizationAssets(sys, filename, false);
            }

            // Create a checkpoint from the current state.
            if (problem == SETTLING || problem == PRESSING) {
                cout << "             Write checkpoint data " << flush;
                if (problem == SETTLING)
                    utils::WriteCheckpoint(sys, settled_ckpnt_file);
                else
                    utils::WriteCheckpoint(sys, pressed_ckpnt_file);
                cout << sys->GetBodies().size() << " bodies" << endl;
            }

            // Increment counters
            out_frame++;
            next_out_frame += out_steps;
        }

        // Check for early termination of a settling phase.
        if (problem == SETTLING || problem == PRESSING) {
            // Store maximum particle height in circular buffer
            hdata[sim_frame % buffer_size] = highest;

            // Check variance of data in circular buffer
            if (time > time_min) {
                double mean_height = hdata.sum() / buffer_size;
                std::valarray<double> x = hdata - mean_height;
                double var = std::sqrt((x * x).sum() / buffer_size);

                // Consider the material settled when the variance is below the
                // specified fraction of a particle radius
                if (var < settling_tol * r_g) {
                    cout << "Granular material settled...  time = " << time << endl;
                    break;
                }
            }
        }

// Advance simulation by one step
#ifdef CHRONO_OPENGL
        if (vis.Run()) {
            sys->DoStepDynamics(time_step);
            vis.Render();
        } else
            break;
#else
        sys->DoStepDynamics(time_step);
#endif

        // Record stats about the simulation
        if (sim_frame % write_steps == 0) {
            // write stat info
            size_t numIters = sys->data_manager->measures.solver.maxd_hist.size();
            double residual = 0;
            if (numIters != 0)
                residual = sys->data_manager->measures.solver.residual;
            statsStream << time << ", " << exec_time << ", " << num_contacts / write_steps << ", " << numIters << ", "
                        << residual << ", \n";
            statsStream.flush();

            num_contacts = 0;
        }

        if (problem == ROLLING || problem == TESTING) {
            // Get the current reaction force or impose shear box position
            ChVector3d rforce_chassis(0, 0, 0);
            ChVector3d rtorque_chassis(0, 0, 0);

            ChVector3d rforce_actuator(0, 0, 0);
            ChVector3d rtorque_actuator(0, 0, 0);

            ChVector3d rforce_wheel(0, 0, 0);
            ChVector3d rtorque_wheel(0, 0, 0);

            auto reaction_chassis = prismatic_chassis_ground->GetReaction2();
            rforce_chassis = reaction_chassis.force;
            rtorque_chassis = reaction_chassis.torque;

            auto reaction_actuator = actuator->GetReaction2();
            rforce_actuator = reaction_actuator.force;
            rtorque_actuator = reaction_actuator.torque;

            auto reaction_wheel = engine_wheel_axle->GetReaction2();
            rforce_wheel = reaction_wheel.force;
            rtorque_wheel = reaction_wheel.torque;

            if (sim_frame % write_steps == 0) {
                ////cout << "X pos: " << xpos_new << " X react: " << cnstr_force << endl;
                rollStream << time << ", " << wheel->GetPos().z() << ", ";

                rollStream << rforce_chassis.x() << ", " << rforce_chassis.y() << ", " << rforce_chassis.z() << ", ";
                rollStream << rforce_actuator.x() << ", " << rforce_actuator.y() << ", " << rforce_actuator.z() << ", ";
                rollStream << rforce_wheel.x() << ", " << rforce_wheel.y() << ", " << rforce_wheel.z() << ", ";

                rollStream << rtorque_chassis.x() << ", " << rtorque_chassis.y() << ", " << rtorque_chassis.z() << ", ";
                rollStream << rtorque_actuator.x() << ", " << rtorque_actuator.y() << ", " << rtorque_actuator.z() << ", ";
                rollStream << rtorque_wheel.x() << ", " << rtorque_wheel.y() << ", " << rtorque_wheel.z() << ", \n";

                rollStream.flush();
            }
        }

        // Increment counters
        time += time_step;
        sim_frame++;
        exec_time += sys->GetTimerStep();
        num_contacts += sys->GetNumContacts();

        // If requested, output detailed timing information for this step
        if (sim_frame == timing_frame)
            sys->PrintStepStats();
    }

    // ----------------
    // Final processing
    // ----------------

    // Create a checkpoint from the last state
    if (problem == SETTLING || problem == PRESSING) {
        cout << "             Write checkpoint data " << flush;
        if (problem == SETTLING)
            utils::WriteCheckpoint(sys, settled_ckpnt_file);
        else
            utils::WriteCheckpoint(sys, pressed_ckpnt_file);
        cout << sys->GetBodies().size() << " bodies" << endl;
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies:  " << sys->GetBodies().size() << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}
