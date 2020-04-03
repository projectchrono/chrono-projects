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
// ChronoParallel demo program for shearing studies.
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
#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"

#include "chrono_thirdparty/filesystem/path.h"

// Control use of OpenGL run-time rendering
// Note: CHRONO_OPENGL is defined in ChConfig.h
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "../utils.h"

using namespace chrono;
using namespace chrono::collision;

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

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

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
int Id_ground = -1;  // body ID for the ground (containing bin)
int Id_box = -2;     // body ID for the shear box
int Id_plate = -3;   // body ID for the load plate

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
int Id_g = 1;          // start body ID for particles
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

void CreateMechanismBodies(ChSystemParallel* system) {
// -------------------------------
// Create a material for the walls
// -------------------------------

#ifdef USE_SMC
    auto mat_walls = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_walls->SetYoungModulus(Y_walls);
    mat_walls->SetFriction(mu_walls);
    mat_walls->SetRestitution(cr_walls);
    mat_walls->SetPoissonRatio(nu_walls);
#else
    auto mat_walls = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_walls->SetFriction(mu_walls);
#endif

    // ----------------------
    // Create the ground body -- always FIRST body in system
    // ----------------------

    auto ground = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());

    ground->SetIdentifier(Id_ground);
    ground->SetBodyFixed(true);
    ground->SetCollide(true);

    // Attach geometry of the containing bin.  Disable contact ground-shearBox
    // and ground-loadPlate.
    ground->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector<>(hdimX, hdimY, hthick), ChVector<>(0, 0, -hthick));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector<>(hthick, hdimY, hdimZ), ChVector<>(-hdimX - hthick, 0, hdimZ));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector<>(hthick, hdimY, hdimZ), ChVector<>(hdimX + hthick, 0, hdimZ));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector<>(hdimX, hthick, hdimZ), ChVector<>(0, -hdimY - hthick, hdimZ));
    utils::AddBoxGeometry(ground.get(), mat_walls, ChVector<>(hdimX, hthick, hdimZ), ChVector<>(0, hdimY + hthick, hdimZ));
    ground->GetCollisionModel()->SetFamily(ground_coll_fam);
    ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(box_coll_fam);
    ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(plate_coll_fam);
    ground->GetCollisionModel()->BuildModel();

    system->AddBody(ground);

// --------------------
// Create the shear box -- always SECOND body in system
// --------------------

// Initially, the shear box is fixed to ground.
// During the shearing phase it may be released (if using an actuator)

    auto box = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());

    box->SetIdentifier(Id_box);
    box->SetPos(ChVector<>(0, 0, 2 * hdimZ + r_g));
    box->SetCollide(true);
    box->SetBodyFixed(true);

    // Add geometry of the shear box.  Disable contact with the load plate.
    box->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector<>(hthick, hdimY, h_scaling * hdimZ),
                          ChVector<>(-hdimX - hthick, 0, h_scaling * hdimZ));
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector<>(hthick, hdimY, h_scaling * hdimZ),
                          ChVector<>(hdimX + hthick, 0, h_scaling * hdimZ));
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector<>(hdimX, hthick, h_scaling * hdimZ),
                          ChVector<>(0, -hdimY - hthick, h_scaling * hdimZ));
    utils::AddBoxGeometry(box.get(), mat_walls, ChVector<>(hdimX, hthick, h_scaling * hdimZ),
                          ChVector<>(0, hdimY + hthick, h_scaling * hdimZ));
    box->GetCollisionModel()->SetFamily(box_coll_fam);
    box->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(plate_coll_fam);
    box->GetCollisionModel()->BuildModel();

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

    auto plate = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());

    plate->SetIdentifier(Id_plate);
    plate->SetMass(mass);
    plate->SetPos(ChVector<>(0, 0, (1 + 2 * h_scaling) * hdimZ));
    plate->SetCollide(true);
    plate->SetBodyFixed(true);

    // Add geometry of the load plate.
    plate->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(plate.get(), mat_walls, ChVector<>(hdimX_p, hdimY, hdimZ), ChVector<>(0, 0, hdimZ));
    plate->GetCollisionModel()->SetFamily(plate_coll_fam);
    plate->GetCollisionModel()->BuildModel();

    system->AddBody(plate);
}

// =============================================================================
// Connect the shear box to the containing bin (ground) through a translational
// joint and create a linear actuator.
// =============================================================================

void ConnectShearBox(ChSystemParallel* system, std::shared_ptr<ChBody> ground, std::shared_ptr<ChBody> box) {
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(ground, box, ChCoordsys<>(ChVector<>(0, 0, 2 * hdimZ), Q_from_AngY(CH_C_PI_2)));
    prismatic->SetName("prismatic_box_ground");
    system->AddLink(prismatic);

    auto actuator_fun = chrono_types::make_shared<ChFunction_Ramp>(0.0, desiredVelocity);

    auto actuator = chrono_types::make_shared<ChLinkLinActuator>();
    ChVector<> pt1(0, 0, 2 * hdimZ);
    ChVector<> pt2(1, 0, 2 * hdimZ);
    actuator->Initialize(ground, box, false, ChCoordsys<>(pt1, QUNIT), ChCoordsys<>(pt2, QUNIT));
    actuator->SetName("actuator");
    actuator->Set_lin_offset(1);
    actuator->Set_dist_funct(actuator_fun);
    system->AddLink(actuator);
}

// =============================================================================
// Connect the load plate to the bin (ground) through a vertical translational
// joint.
// =============================================================================

void ConnectLoadPlate(ChSystemParallel* system, std::shared_ptr<ChBody> ground, std::shared_ptr<ChBody> plate) {
    auto prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic->Initialize(ground, plate, ChCoordsys<>(ChVector<>(0, 0, 2 * hdimZ), QUNIT));
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

int CreateGranularMaterial(ChSystemParallel* system) {
// -------------------------------------------
// Create a material for the granular material
// -------------------------------------------

#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
    mat_g->SetPoissonRatio(nu_g);
#else
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // ---------------------------------------------
    // Create a mixture entirely made out of spheres
    // ---------------------------------------------

    // Create the particle generator with a mixture of 100% spheres
    utils::Generator gen(system);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    // Ensure that all generated particle bodies will have positive IDs.
    gen.setBodyIdentifier(Id_g);

    // ----------------------
    // Generate the particles
    // ----------------------

    double r = 1.01 * r_g;
    ChVector<> hdims(hdimX - r, hdimY - r, 0);
    ChVector<> center(0, 0, 2 * r);

    while (center.z() < 2 * h_scaling * hdimZ) {
        gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims);
        center.z() += 2 * r;
    }

    // Return the number of generated particles.
    return gen.getTotalNumBodies();
}

// =============================================================================
// Create a single large sphere (for use in TESTING)
// =============================================================================

void CreateBall(ChSystemParallel* system) {
// ------------------------------
// Create a material for the ball
// ------------------------------

#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
    mat_g->SetPoissonRatio(nu_g);
#else
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // ---------------
    // Create the ball
    // ---------------

    auto ball = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());

    ball->SetIdentifier(Id_ball);
    ball->SetMass(mass_ball);
    ball->SetPos(ChVector<>(0, 0, 1.01 * radius_ball));
    ball->SetCollide(true);
    ball->SetBodyFixed(false);

    ball->GetCollisionModel()->ClearModel();
    utils::AddSphereGeometry(ball.get(), mat_g, radius_ball);
    ball->GetCollisionModel()->BuildModel();

    system->AddBody(ball);
}

// =============================================================================
// Find the height of the highest and lowest sphere in the granular mix. We only
// look at bodies whith positive identifiers (to exclude all other bodies).
// =============================================================================

void FindHeightRange(ChSystemParallel* sys, double& lowest, double& highest) {
    highest = -1000;
    lowest = 1000;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() <= 0)
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
    double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;

    double normalPlateHeight = sys->Get_bodylist().at(1)->GetPos().z() - hdimZ;
    double bottomHeight = 0;
    double boxVolume = hdimX * 2 * hdimX * 2 * (normalPlateHeight - bottomHeight);
    double granularVolume = (sys->Get_bodylist().size() - 3) * vol_g;
    double reqDensity = bulkDensity * boxVolume / granularVolume;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() > 1) {
            body->SetMass(reqDensity * vol_g);
        }
    }

    cout << "N Bodies: " << sys->Get_bodylist().size() << endl;
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
    ChSystemParallelSMC* msystem = new ChSystemParallelSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemParallelNSC* msystem = new ChSystemParallelNSC();
#endif

    msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);
    cout << "Using " << threads << " threads" << endl;

    msystem->GetSettings()->max_threads = threads;
    msystem->GetSettings()->perform_thread_tuning = thread_tuning;

    // Edit system settings
    msystem->GetSettings()->solver.use_full_inertia_tensor = false;
    msystem->GetSettings()->solver.tolerance = tolerance;
    msystem->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    msystem->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
    msystem->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

#ifdef USE_SMC
    msystem->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_R;
    msystem->GetSettings()->solver.contact_force_model = contact_force_model;
    msystem->GetSettings()->solver.tangential_displ_mode = tangential_displ_mode;
#else
    msystem->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    msystem->GetSettings()->solver.alpha = 0;
    msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    msystem->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
    msystem->ChangeSolverType(SolverType::APGDREF);

    msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

    msystem->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

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
    std::shared_ptr<ChLinkLinActuator> actuator;

    switch (problem) {
        case SETTLING: {
            time_min = time_settling_min;
            time_end = time_settling_max;
            out_fps = out_fps_settling;

            // Create the mechanism bodies (all fixed).
            CreateMechanismBodies(msystem);

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = msystem->Get_bodylist().at(0);
            shearBox = msystem->Get_bodylist().at(1);
            loadPlate = msystem->Get_bodylist().at(2);

            // Create granular material.
            int num_particles = CreateGranularMaterial(msystem);
            cout << "Granular material:  " << num_particles << " particles" << endl;

            break;
        }

        case PRESSING: {
            time_min = time_pressing_min;
            time_end = time_pressing_max;
            out_fps = out_fps_pressing;

            // Create bodies from checkpoint file.
            cout << "Read checkpoint data from " << settled_ckpnt_file;
            utils::ReadCheckpoint(msystem, settled_ckpnt_file);
            cout << "  done.  Read " << msystem->Get_bodylist().size() << " bodies." << endl;

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = msystem->Get_bodylist().at(0);
            shearBox = msystem->Get_bodylist().at(1);
            loadPlate = msystem->Get_bodylist().at(2);

            // Move the load plate just above the granular material.
            double highest, lowest;
            FindHeightRange(msystem, lowest, highest);
            ChVector<> pos = loadPlate->GetPos();
            double z_new = highest + 2 * r_g;
            loadPlate->SetPos(ChVector<>(pos.x(), pos.y(), z_new));

            // Connect the load plate to the shear box.
            ConnectLoadPlate(msystem, ground, loadPlate);
            prismatic_plate_ground = std::static_pointer_cast<ChLinkLockPrismatic>(msystem->SearchLink("prismatic_plate_ground"));

            // Release the load plate.
            loadPlate->SetBodyFixed(false);

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
            utils::ReadCheckpoint(msystem, pressed_ckpnt_file);
            cout << "  done.  Read " << msystem->Get_bodylist().size() << " bodies." << endl;

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = msystem->Get_bodylist().at(0);
            shearBox = msystem->Get_bodylist().at(1);
            loadPlate = msystem->Get_bodylist().at(2);

            // If using an actuator, connect the shear box and get a handle to the actuator.
            if (use_actuator) {
                ConnectShearBox(msystem, ground, shearBox);
                prismatic_box_ground = std::static_pointer_cast<ChLinkLockPrismatic>(msystem->SearchLink("prismatic_box_ground"));
                actuator = std::static_pointer_cast<ChLinkLinActuator>(msystem->SearchLink("actuator"));
            }

            // Release the shear box when using an actuator.
            shearBox->SetBodyFixed(!use_actuator);

            // Connect the load plate to the shear box.
            ConnectLoadPlate(msystem, ground, loadPlate);
            prismatic_plate_ground = std::static_pointer_cast<ChLinkLockPrismatic>(msystem->SearchLink("prismatic_plate_ground"));

            // Release the load plate.
            loadPlate->SetBodyFixed(false);

            // setBulkDensity(msystem, desiredBulkDensity);

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
            CreateMechanismBodies(msystem);

            // Create the test ball.
            CreateBall(msystem);

            // Grab handles to mechanism bodies (must increase ref counts)
            ground = msystem->Get_bodylist().at(0);
            shearBox = msystem->Get_bodylist().at(1);
            loadPlate = msystem->Get_bodylist().at(2);

            // Move the load plate just above the test ball.
            ChVector<> pos = loadPlate->GetPos();
            double z_new = 2.1 * radius_ball;
            loadPlate->SetPos(ChVector<>(pos.x(), pos.y(), z_new));

            // If using an actuator, connect the shear box and get a handle to the actuator.
            if (use_actuator) {
                ConnectShearBox(msystem, ground, shearBox);
                prismatic_box_ground = std::static_pointer_cast<ChLinkLockPrismatic>(msystem->SearchLink("prismatic_box_ground"));
                actuator = std::static_pointer_cast<ChLinkLinActuator>(msystem->SearchLink("actuator"));
            }

            // Release the shear box when using an actuator.
            shearBox->SetBodyFixed(!use_actuator);

            // Connect the load plate to the shear box.
            ConnectLoadPlate(msystem, ground, loadPlate);
            prismatic_plate_ground = std::static_pointer_cast<ChLinkLockPrismatic>(msystem->SearchLink("prismatic_plate_ground"));

            // Release the load plate.
            loadPlate->SetBodyFixed(false);

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
    ChStreamOutAsciiFile statsStream(stats_file.c_str());
    ChStreamOutAsciiFile shearStream(shear_file.c_str());

    shearStream.SetNumFormat("%16.4e");

#ifdef CHRONO_OPENGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Direct Shear Test", msystem);
    gl_window.SetCamera(ChVector<>(0, -10 * hdimY, hdimZ), ChVector<>(0, 0, hdimZ), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

    // Loop until reaching the end time...
    while (time < time_end) {
        // Current position and velocity of the shear box
        ChVector<> pos_old = shearBox->GetPos();
        ChVector<> vel_old = shearBox->GetPos_dt();

        // Calculate minimum and maximum particle heights
        double highest, lowest;
        FindHeightRange(msystem, lowest, highest);

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
                utils::WriteShapesPovray(msystem, filename, false);
            }

            // Create a checkpoint from the current state.
            if (problem == SETTLING || problem == PRESSING) {
                cout << "             Write checkpoint data " << flush;
                if (problem == SETTLING)
                    utils::WriteCheckpoint(msystem, settled_ckpnt_file);
                else
                    utils::WriteCheckpoint(msystem, pressed_ckpnt_file);
                cout << msystem->Get_bodylist().size() << " bodies" << endl;
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
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else
            break;
#else
        msystem->DoStepDynamics(time_step);
#endif

        ////progressbar(out_steps + sim_frame - next_out_frame + 1, out_steps);
        TimingOutput(msystem);

        // Record stats about the simulation
        if (sim_frame % write_steps == 0) {
            // write stat info
            size_t numIters = msystem->data_manager->measures.solver.maxd_hist.size();
            double residual = 0;
            if (numIters != 0)
                residual = msystem->data_manager->measures.solver.residual;
            statsStream << time << ", " << exec_time << ", " << num_contacts / write_steps << ", " << numIters << ", "
                        << residual << ", " << max_cnstr_viol[0] << ", " << max_cnstr_viol[1] << ", "
                        << max_cnstr_viol[2] << ", \n";
            statsStream.GetFstream().flush();

            num_contacts = 0;
            max_cnstr_viol[0] = 0;
            max_cnstr_viol[1] = 0;
            max_cnstr_viol[2] = 0;
        }

        if (problem == SHEARING || problem == TESTING) {
            // Get the current reaction force or impose shear box position
            ChVector<> rforcePbg(0, 0, 0);
            ChVector<> rtorquePbg(0, 0, 0);

            ChVector<> rforcePpb(0, 0, 0);
            ChVector<> rtorquePpb(0, 0, 0);

            ChVector<> rforceA(0, 0, 0);
            ChVector<> rtorqueA(0, 0, 0);

            if (use_actuator) {
                rforcePbg = prismatic_box_ground->Get_react_force();
                rtorquePbg = prismatic_box_ground->Get_react_torque();

                rforcePpb = prismatic_plate_ground->Get_react_force();
                rtorquePpb = prismatic_plate_ground->Get_react_torque();

                rforceA = actuator->Get_react_force();
                rtorqueA = actuator->Get_react_torque();
            } else {
                double xpos_new = pos_old.x() + desiredVelocity * time_step;
                shearBox->SetPos(ChVector<>(xpos_new, pos_old.y(), pos_old.z()));
                shearBox->SetPos_dt(ChVector<>(desiredVelocity, 0, 0));
            }

            if (sim_frame % write_steps == 0) {
                ////cout << "X pos: " << xpos_new << " X react: " << cnstr_force << endl;
                shearStream << time << "  " << shearBox->GetPos().x() << "     ";
                shearStream << rforceA.x() << "  " << rforceA.y() << "  " << rforceA.z() << "     ";
                shearStream << rtorqueA.x() << "  " << rtorqueA.y() << "  " << rtorqueA.z() << "\n";
                shearStream.GetFstream().flush();
            }
        }

        // Find maximum constraint violation
        if (prismatic_box_ground) {
            ChVectorDynamic<> C = prismatic_box_ground->GetC();
            for (int i = 0; i < 5; i++)
                max_cnstr_viol[0] = std::max(max_cnstr_viol[0], std::abs(C(i)));
        }
        if (prismatic_plate_ground) {
            ChVectorDynamic<> C = prismatic_plate_ground->GetC();
            for (int i = 0; i < 5; i++)
                max_cnstr_viol[1] = std::max(max_cnstr_viol[1], std::abs(C(i)));
        }
        if (actuator) {
            ChVectorDynamic<> C = actuator->GetC();
            max_cnstr_viol[2] = std::max(max_cnstr_viol[2], std::abs(C(0)));
        }

        // Increment counters
        time += time_step;
        sim_frame++;
        exec_time += msystem->GetTimerStep();
        num_contacts += msystem->GetNcontacts();

        // If requested, output detailed timing information for this step
        if (sim_frame == timing_frame)
            msystem->PrintStepStats();
    }

    // ----------------
    // Final processing
    // ----------------

    // Create a checkpoint from the last state
    if (problem == SETTLING || problem == PRESSING) {
        cout << "             Write checkpoint data " << flush;
        if (problem == SETTLING)
            utils::WriteCheckpoint(msystem, settled_ckpnt_file);
        else
            utils::WriteCheckpoint(msystem, pressed_ckpnt_file);
        cout << msystem->Get_bodylist().size() << " bodies" << endl;
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies:  " << msystem->Get_bodylist().size() << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}
