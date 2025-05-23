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
// Chrono::Multicore demo program for shearing studies.
//
// The system contains a shearing box composed of three bodies: (1) a containing
// bin, (2) a shearing plate, and (3) a load body. Granular material sits inside
// of the containing bin and shearing plate and is compressed from the top by
// the load body. During the shearing mode, the shear plate is translated in the
// x-direction at a specified velocity.
//
// The global reference frame has Z up.
// All units SI (CGS, i.e., centimeter - gram - second)
// =============================================================================

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChSystemDescriptorMulticore.h"

#include "chrono_thirdparty/filesystem/path.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

#include "../utils.h"

using namespace chrono;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use NSC contact
//#define USE_SMC

enum ProblemType { SETTLING, PRESSING, SHEARING, TESTING };

ProblemType problem = TESTING;

// -----------------------------------------------------------------------------
// Conversion factors
// -----------------------------------------------------------------------------

// Conversion for Young's modulus and pressure
// [Y] = Pa = N / m^2 = kg / m / s^2
double Pa2cgs = 10;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform shearing action via a linear actuator or kinematically?
bool use_actuator = true;

// Save PovRay post-processing data?
bool write_povray_data = true;

// Simulation times
double time_settling_min = 0.1;
double time_settling_max = 1.0;

double time_pressing_min = 0.1;
double time_pressing_max = 1.0;

double time_shearing = 5;

double time_testing = 2;

// Stopping criteria for settling (as fraction of particle radius)
double settling_tol = 0.2;

// Solver settings
#ifdef USE_SMC
double time_step = 1e-5;
int max_iteration_bilateral = 100;
#else
double time_step = 1e-4;
int max_iteration_normal = 0;
int max_iteration_sliding = 10000;
int max_iteration_spinning = 0;
int max_iteration_bilateral = 100;
double contact_recovery_speed = 10e30;
#endif

bool clamp_bilaterals = false;
double bilateral_clamp_speed = 10e30;
double tolerance = 1;

// Contact force model
#ifdef USE_SMC
ChSystemSMC::ContactForceModel contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
ChSystemSMC::TangentialDisplacementModel tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::MultiStep;
#endif

// Output
#ifdef USE_SMC
const std::string out_dir = "../DIRECTSHEAR_SMC";
#else
const std::string out_dir = "../DIRECTSHEAR_NSC";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string shear_file = out_dir + "/shear.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string settled_ckpnt_file = out_dir + "/settled.dat";
const std::string pressed_ckpnt_file = out_dir + "/pressed.dat";

// Frequency for visualization output
int out_fps_settling = 120;
int out_fps_pressing = 120;
int out_fps_shearing = 60;
int out_fps_testing = 60;

// Frequency for writing results to output file (SHEARING only)
int write_fps = 1000;

// Simulation frame at which detailed timing information is printed
int timing_frame = -1;

// Gravitational acceleration [cm/s^2]
double gravity = 981;

// Parameters for the mechanism
double hdimX = 12.0 / 2;  // [cm] bin half-length in x direction
double hdimY = 12.0 / 2;  // [cm] bin half-depth in y direction
double hdimZ = 3.0 / 2;   // [cm] bin half-height in z direction
double hthick = 1.0 / 2;  // [cm] bin half-thickness of the walls

double h_scaling = 5;  // ratio of shear box height to bin height

float Y_walls = (float)(Pa2cgs * 1e7);
float cr_walls = 0.6f;
float nu_walls = 0.3f;
float mu_walls = 0.08f;

int ground_coll_fam = 1;  // collision family for bin contact shapes
int box_coll_fam = 2;     // collision family for shear box contact shapes
int plate_coll_fam = 3;   // collision family for load plate contact shapes

// Applied normal pressure
double normalPressure = Pa2cgs * 3.1e3;  // 3.1 kPa // 6.4 kPa // 12.5 kPa // 24.2 kPa

// Desired shearing velocity [cm/s]
double desiredVelocity = 0.166;  // 10 cm/min (about 100 times faster than experiment)

// Parameters for the granular material
int tag_particles = 0; // start body tag for particles
double r_g = 0.3;      // [cm] radius of granular spheres
double rho_g = 2.550;  // [g/cm^3] density of granules

double desiredBulkDensity = 1.5;  // [g/cm^3] desired bulk density

float Y_g = (float)(Pa2cgs * 4e7);  // (1,000 times softer than experiment on glass beads)
float cr_g = 0.87f;
float nu_g = 0.22f;
float mu_g = 0.18f;

// Parameters of the testing ball
int Id_ball = -4;
double mass_ball = 200;            // [g] mass of testing ball
double radius_ball = 0.9 * hdimX;  // [cm] radius of testing ball

// =============================================================================
// Create the containing bin (the ground), the shear box, and the load plate.
//
// Both the shear box and load plate are constructed fixed to the ground.
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
    mat_walls->SetRestitution(cr_walls);
    mat_walls->SetPoissonRatio(nu_walls);
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
    ground->GetCollisionModel()->SetFamily(ground_coll_fam);
    ground->GetCollisionModel()->DisallowCollisionsWith(box_coll_fam);
    ground->GetCollisionModel()->DisallowCollisionsWith(plate_coll_fam);

    system->AddBody(ground);

// --------------------
// Create the shear box -- always SECOND body in system
// --------------------

// Initially, the shear box is fixed to ground.
// During the shearing phase it may be released (if using an actuator)

    auto box = chrono_types::make_shared<ChBody>();

    box->SetPos(ChVector3d(0, 0, 2 * hdimZ + r_g));
    box->EnableCollision(true);
    box->SetFixed(true);

    // Add geometry of the shear box.  Disable contact with the load plate.
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector3d(hthick, hdimY, h_scaling * hdimZ),
                          ChVector3d(-hdimX - hthick, 0, h_scaling * hdimZ));
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector3d(hthick, hdimY, h_scaling * hdimZ),
                          ChVector3d(hdimX + hthick, 0, h_scaling * hdimZ));
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector3d(hdimX, hthick, h_scaling * hdimZ),
                          ChVector3d(0, -hdimY - hthick, h_scaling * hdimZ));
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector3d(hdimX, hthick, h_scaling * hdimZ),
                          ChVector3d(0, hdimY + hthick, h_scaling * hdimZ));
    box->GetCollisionModel()->SetFamily(box_coll_fam);
    box->GetCollisionModel()->DisallowCollisionsWith(plate_coll_fam);

    system->AddBody(box);

    // ---------------------
    // Create the plate body -- always THIRD body in the system
    // ---------------------

    // Initially, the load plate is fixed to ground.
    // It is released after the settling phase.

    // Set plate dimensions, increasing the X dimension to accommodate the
    // shearing phase (use 3 times as much as needed)
    double hdimX_p = hdimX + 3 * time_shearing * desiredVelocity;

    // Estimate plate mass from desired applied normal pressure
    double area = 4 * hdimX * hdimY;
    double mass = normalPressure * area / gravity;

    auto plate = chrono_types::make_shared<ChBody>();

    plate->SetMass(mass);
    plate->SetPos(ChVector3d(0, 0, (1 + 2 * h_scaling) * hdimZ));
    plate->EnableCollision(true);
    plate->SetFixed(true);

    // Add geometry of the load plate.
    utils::AddBoxGeometry(plate.get(), mat_walls, ChVector3d(hdimX_p, hdimY, hdimZ), ChVector3d(0, 0, hdimZ));
    plate->GetCollisionModel()->SetFamily(plate_coll_fam);

    system->AddBody(plate);
}

// =============================================================================
// Connect the shear box to the containing bin (ground) through a translational
// joint and create a linear actuator.
// =============================================================================

void ConnectShearBox(ChSystemMulticore* system, std::shared_ptr<ChBody> ground, std::shared_ptr<ChBody> box) {
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(ground, box, ChFrame<>(ChVector3d(0, 0, 2 * hdimZ), QuatFromAngleY(CH_PI_2)));
    prismatic->SetName("prismatic_box_ground");
    system->AddLink(prismatic);

    auto actuator_fun = chrono_types::make_shared<ChFunctionRamp>(0.0, desiredVelocity);

    auto actuator = chrono_types::make_shared<ChLinkLockLinActuator>();
    ChVector3d pt1(0, 0, 2 * hdimZ);
    ChVector3d pt2(1, 0, 2 * hdimZ);
    actuator->Initialize(ground, box, false, ChFrame<>(pt1, QUNIT), ChFrame<>(pt2, QUNIT));
    actuator->SetName("actuator");
    actuator->SetDistanceOffset(1);
    actuator->SetActuatorFunction(actuator_fun);
    system->AddLink(actuator);
}

// =============================================================================
// Connect the load plate to the bin (ground) through a vertical translational
// joint.
// =============================================================================

void ConnectLoadPlate(ChSystemMulticore* system, std::shared_ptr<ChBody> ground, std::shared_ptr<ChBody> plate) {
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(ground, plate, ChFrame<>(ChVector3d(0, 0, 2 * hdimZ), QUNIT));
    prismatic->SetName("prismatic_plate_ground");
    system->AddLink(prismatic);
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
    mat_g->SetRestitution(cr_g);
    mat_g->SetPoissonRatio(nu_g);
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

    while (center.z() < 2 * h_scaling * hdimZ) {
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
    mat_g->SetRestitution(cr_g);
    mat_g->SetPoissonRatio(nu_g);
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
// look at bodies whith positive identifiers (to exclude all other bodies).
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
    sys->GetSettings()->solver.use_full_inertia_tensor = false;
    sys->GetSettings()->solver.tolerance = tolerance;
    sys->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    sys->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
    sys->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

#ifdef USE_SMC
    sys->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_R;
    sys->GetSettings()->solver.contact_force_model = contact_force_model;
    sys->GetSettings()->solver.tangential_displ_mode = tangential_displ_mode;
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
    std::shared_ptr<ChBody> shearBox;
    std::shared_ptr<ChBody> loadPlate;
    std::shared_ptr<ChLinkLockPrismatic> prismatic_box_ground;
    std::shared_ptr<ChLinkLockPrismatic> prismatic_plate_ground;
    std::shared_ptr<ChLinkLockLinActuator> actuator;

    switch (problem) {
        case SETTLING: {
            time_min = time_settling_min;
            time_end = time_settling_max;
            out_fps = out_fps_settling;

            // Create the mechanism bodies (all fixed).
            CreateMechanismBodies(sys);

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = sys->GetBodies().at(0);
            shearBox = sys->GetBodies().at(1);
            loadPlate = sys->GetBodies().at(2);

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
            shearBox = sys->GetBodies().at(1);
            loadPlate = sys->GetBodies().at(2);

            // Move the load plate just above the granular material.
            double highest, lowest;
            FindHeightRange(sys, lowest, highest);
            ChVector3d pos = loadPlate->GetPos();
            double z_new = highest + 2 * r_g;
            loadPlate->SetPos(ChVector3d(pos.x(), pos.y(), z_new));

            // Connect the load plate to the shear box.
            ConnectLoadPlate(sys, ground, loadPlate);
            prismatic_plate_ground = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_plate_ground"));

            // Release the load plate.
            loadPlate->SetFixed(false);

            // Set plate mass from desired applied normal pressure
            double area = 4 * hdimX * hdimY;
            double mass = normalPressure * area / gravity;
            loadPlate->SetMass(mass);

            break;
        }

        case SHEARING: {
            time_end = time_shearing;
            out_fps = out_fps_shearing;

            // Create bodies from checkpoint file.
            cout << "Read checkpoint data from " << pressed_ckpnt_file;
            utils::ReadCheckpoint(sys, pressed_ckpnt_file);
            cout << "  done.  Read " << sys->GetBodies().size() << " bodies." << endl;

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = sys->GetBodies().at(0);
            shearBox = sys->GetBodies().at(1);
            loadPlate = sys->GetBodies().at(2);

            // If using an actuator, connect the shear box and get a handle to the actuator.
            if (use_actuator) {
                ConnectShearBox(sys, ground, shearBox);
                prismatic_box_ground = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_box_ground"));
                actuator = std::static_pointer_cast<ChLinkLockLinActuator>(sys->SearchLink("actuator"));
            }

            // Release the shear box when using an actuator.
            shearBox->SetFixed(!use_actuator);

            // Connect the load plate to the shear box.
            ConnectLoadPlate(sys, ground, loadPlate);
            prismatic_plate_ground = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_plate_ground"));

            // Release the load plate.
            loadPlate->SetFixed(false);

            // setBulkDensity(sys, desiredBulkDensity);

            // Set plate mass from desired applied normal pressure
            double area = 4 * hdimX * hdimY;
            double mass = normalPressure * area / gravity;
            loadPlate->SetMass(mass);

            break;
        }

        case TESTING: {
            time_end = time_testing;
            out_fps = out_fps_testing;

            // For TESTING only, increase shearing velocity.
            desiredVelocity = 0.5;

            // Create the mechanism bodies (all fixed).
            CreateMechanismBodies(sys);

            // Create the test ball.
            CreateBall(sys);

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = sys->GetBodies().at(0);
            shearBox = sys->GetBodies().at(1);
            loadPlate = sys->GetBodies().at(2);

            // Move the load plate just above the test ball.
            ChVector3d pos = loadPlate->GetPos();
            double z_new = 2.1 * radius_ball;
            loadPlate->SetPos(ChVector3d(pos.x(), pos.y(), z_new));

            // If using an actuator, connect the shear box and get a handle to the actuator.
            if (use_actuator) {
                ConnectShearBox(sys, ground, shearBox);
                prismatic_box_ground = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_box_ground"));
                actuator = std::static_pointer_cast<ChLinkLockLinActuator>(sys->SearchLink("actuator"));
            }

            // Release the shear box when using an actuator.
            shearBox->SetFixed(!use_actuator);

            // Connect the load plate to the shear box.
            ConnectLoadPlate(sys, ground, loadPlate);
            prismatic_plate_ground = std::static_pointer_cast<ChLinkLockPrismatic>(sys->SearchLink("prismatic_plate_ground"));

            // Release the load plate.
            loadPlate->SetFixed(false);

            // Set plate mass from desired applied normal pressure
            double area = 4 * hdimX * hdimY;
            double mass = normalPressure * area / gravity;
            loadPlate->SetMass(mass);

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
    std::ofstream shearStream(shear_file.c_str());

#ifdef CHRONO_OPENGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(sys);
    vis.SetWindowTitle("Direct Shear Test");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(0, -10 * hdimY, hdimZ), ChVector3d(0, 0, hdimZ));
    vis.SetCameraVertical(CameraVerticalDir::Z);
#endif

    // Loop until reaching the end time...
    while (time < time_end) {
        // Current position and velocity of the shear box
        ChVector3d pos_old = shearBox->GetPos();
        ChVector3d vel_old = shearBox->GetLinVel();

        // Calculate minimum and maximum particle heights
        double highest, lowest;
        FindHeightRange(sys, lowest, highest);

        // If at an output frame, write PovRay file and print info
        if (sim_frame == next_out_frame) {
            cout << "------------ Output frame:     " << out_frame + 1 << endl;
            cout << "             Sim frame:        " << sim_frame << endl;
            cout << "             Time:             " << time << endl;
            cout << "             Shear box pos:    " << pos_old.x() << endl;
            cout << "                       vel:    " << vel_old.x() << endl;
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

        ////progressbar(out_steps + sim_frame - next_out_frame + 1, out_steps);
        TimingOutput(sys);

        // Record stats about the simulation
        if (sim_frame % write_steps == 0) {
            // write stat info
            size_t numIters = sys->data_manager->measures.solver.maxd_hist.size();
            double residual = 0;
            if (numIters != 0)
                residual = sys->data_manager->measures.solver.residual;
            statsStream << time << ", " << exec_time << ", " << num_contacts / write_steps << ", " << numIters << ", "
                        << residual << ", " << max_cnstr_viol[0] << ", " << max_cnstr_viol[1] << ", "
                        << max_cnstr_viol[2] << ", \n";
            statsStream.flush();

            num_contacts = 0;
            max_cnstr_viol[0] = 0;
            max_cnstr_viol[1] = 0;
            max_cnstr_viol[2] = 0;
        }

        if (problem == SHEARING || problem == TESTING) {
            // Get the current reaction force or impose shear box position
            ChVector3d rforcePbg(0, 0, 0);
            ChVector3d rtorquePbg(0, 0, 0);

            ChVector3d rforcePpb(0, 0, 0);
            ChVector3d rtorquePpb(0, 0, 0);

            ChVector3d rforceA(0, 0, 0);
            ChVector3d rtorqueA(0, 0, 0);

            if (use_actuator) {
                auto reactionPbg = prismatic_box_ground->GetReaction2();
                rforcePbg = reactionPbg.force;
                rtorquePbg = reactionPbg.torque;

                auto reactionPpb = prismatic_plate_ground->GetReaction2();
                rforcePpb = reactionPpb.force;
                rtorquePpb = reactionPpb.torque;

                auto reactionA = actuator->GetReaction2();
                rforceA = reactionA.force;
                rtorqueA = reactionA.torque;
            } else {
                double xpos_new = pos_old.x() + desiredVelocity * time_step;
                shearBox->SetPos(ChVector3d(xpos_new, pos_old.y(), pos_old.z()));
                shearBox->SetLinVel(ChVector3d(desiredVelocity, 0, 0));
            }

            if (sim_frame % write_steps == 0) {
                ////cout << "X pos: " << xpos_new << " X react: " << cnstr_force << endl;
                shearStream << time << "  " << shearBox->GetPos().x() << "     ";
                shearStream << rforceA.x() << "  " << rforceA.y() << "  " << rforceA.z() << "     ";
                shearStream << rtorqueA.x() << "  " << rtorqueA.y() << "  " << rtorqueA.z() << "\n";
                shearStream.flush();
            }
        }

        // Find maximum constraint violation
        if (prismatic_box_ground) {
            ChVectorDynamic<> C = prismatic_box_ground->GetConstraintViolation();
            for (int i = 0; i < 5; i++)
                max_cnstr_viol[0] = std::max(max_cnstr_viol[0], std::abs(C(i)));
        }
        if (prismatic_plate_ground) {
            ChVectorDynamic<> C = prismatic_plate_ground->GetConstraintViolation();
            for (int i = 0; i < 5; i++)
                max_cnstr_viol[1] = std::max(max_cnstr_viol[1], std::abs(C(i)));
        }
        if (actuator) {
            ChVectorDynamic<> C = actuator->GetConstraintViolation();
            max_cnstr_viol[2] = std::max(max_cnstr_viol[2], std::abs(C(0)));
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
